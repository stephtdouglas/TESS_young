import os, glob
import multiprocessing as mp

import numpy as np
import matplotlib.pyplot as plt
import astropy.io.ascii as at
from astropy.table import Table

import get_colors
norm, mapper, cmap2, colors, shapes = get_colors.get_colors()
plt.style.use('./paper.mplstyle')
PAPER_DIR = os.path.expanduser("~/my_papers/TESS_young/")

# from plot_periods import plot_periodcolor_histogram
from tau_sq import PeriodMassDistribution, SpinModel




clusters = ["IC_2391","Collinder_135","NGC_2451A","NGC_2547","IC_2602"]
dates = ["2021-06-22","2021-06-18","2021-06-21","2021-06-21","2021-07-02"]


display_names = {"UpSco_Mattea2015":"Matt+15; UpSco initialization",
                 "UpSco_Mattea2022":"Matt+in prep; UpSco initialization",
                 "UpSco_ZeroTorque":"Zero Torque; UpSco initialization",
                 "WideHat8Myr_Mattea2015":"Matt+15; uniform initialization",
                 "WideHat8Myr_Mattea2022":"Matt+in prep; uniform initialization",
                 "WideHat8Myr_ZeroTorque":"Zero Torque; uniform initialization"}
model_names = ["UpSco_Mattea2015","UpSco_Mattea2022","UpSco_ZeroTorque",
               "WideHat8Myr_Mattea2015","WideHat8Myr_Mattea2022","WideHat8Myr_ZeroTorque"]
nmod = 6
nage = 118

def run_one_model(age,pmd,model,period_scale,init_type):
    sm = SpinModel(model,age,period_scale,init_type=init_type)

    # Normalize the model and calculate tau-squared
    if init_type!="kde":
        sm.normalize()
    sm.add_mask()
    sm.calc_tau_sq(pmd)

    return sm.tau_sq

def run_one_model_binned(age,pmd,model,period_scale,init_type):
    sm = SpinModel(model,age,period_scale,init_type=init_type)

    # Normalize the model and calculate tau-squared
    if init_type!="kde":
        sm.normalize()
    sm.add_mask()
    sm.calc_tau_sq_binned(pmd)

    return sm.tau_sq


def run_all_models(config_file):
        # max_q=0,include_blends=True,include_lit=False,
        #            period_scale="linear",output_filebase="tausq_ZAMS_Compare",
        #            models_to_plot=model_names,zoom_ymax=None,
        #            mass_limits=None,pmd=None,to_plot=True):
    
    # parse config file
    config_file = os.path.abspath(os.path.expanduser(config_file))
    with open(config_file, 'r') as f:
        config = yaml.load(f.read())
        config['config_file'] = config_file

    print(config)

    max_q=config["max_q"]
    models_to_plot=config["models"]
    output_filebase=config["output_filebase"]
    zoom_ymax=config["zoom_ymax"]
    mass_limits=config["mass_limits"]
    include_blends=config["include_blends"]
    include_lit=config["include_lit"]
    to_plot=config["to_plot"]
    period_scale=config["period_scale"]
    
    nmod_l = len(models_to_plot)

    pmd = PeriodMassDistribution(max_q,include_blends,include_lit,
                                    mass_limits=mass_limits)

    # Check for the matching output csv and skip straight to plotting if found
    outfilename = f"{output_filebase}_{pmd.param_string}"
    if (config["overwrite"]==False) and (os.path.exists(f"tables/{outfilename}.csv")):
        print("computation already completed")
        run_fits = False
        if to_plot:
            ttab = at.read(f"tables/{outfilename}.csv")
        else:
            return
    else:
        run_fits = True
        ttab = Table(np.zeros(nmod_l*nage).reshape(nage,nmod_l),names=models_to_plot)

    # Set up figure
    if to_plot:
        fig = plt.figure()
        fig.patch.set_facecolor('w')
        fig.patch.set_alpha(1.0)
        ax = plt.subplot(111)

    # Determine cpu count for parallelization
    try:
        # The number of actually available CPUs, important for HPC
        cpu_count = len(os.sched_getaffinity(0))
    except AttributeError:
        # On a regular computer, just use the total CPU count
        cpu_count = mp.cpu_count()-2


    # Run comparison to data for every model
    if run_fits:
        for j,model in enumerate(models_to_plot):

            print(model)

            if "WideHat" in model:
                init_type="kde"
            else:
                init_type="cluster"

            models = glob.glob(os.path.expanduser(f"~/Dropbox/Models/{model}/{model}*Myr.txt"))
            # print(models)

            model_ages = np.sort([int(mod.split("_")[-1][:5]) for mod in models])
            # print(model_ages)

            model_ages = model_ages#[model_ages<=300]
            age_colname = f"Age_{model}"
            ttab[age_colname] = model_ages

            all_tau_sq = np.zeros(len(model_ages))

            pool = mp.Pool(cpu_count)
            chunksize = len(model_ages) // cpu_count

            print(f"{cpu_count} CPUs, chunk size {chunksize}")
            print("starting multiprocessing run",model)

            tau_sq_args = [[age,pmd,model,period_scale,init_type] for
                            age in model_ages]

            all_tau_sq = pool.starmap(run_one_model,tau_sq_args,
                                          chunksize=chunksize)

            if "UpSco" in model:
                ls = "--"
            else:
                ls = "-"

            if to_plot:
                ax.plot(model_ages,all_tau_sq,ls,label=display_names[model],
                        color=mapper.to_rgba((j % 3)+1),alpha=0.75)
            ttab[model] = all_tau_sq

            if j==0:
                ttab["Age(Myr)"] = model_ages

    # If the comparison was already run, just re-plot
    else:
        for j,model in enumerate(models_to_plot):
            age_colname = f"Age_{model}"
            if "UpSco" in model:
                ls = "--"
            else:
                ls = "-"

            if to_plot:
                ax.plot(ttab[age_colname],ttab[model],ls,
                        label=display_names[model],
                        color=mapper.to_rgba((j % 3)+1),alpha=0.75)
        # Plot the models from the saved file

    if to_plot:
        ax.legend(loc=2)
        ax.set_xlabel("Model age (Myr)",fontsize=16)
        ax.set_ylabel("tau squared",fontsize=16)

        ax.tick_params(labelsize=12)
        ax.set_xticks(np.arange(0,300,25),minor=True)

        fig.savefig(f"plots/{outfilename}.png",bbox_inches="tight",dpi=600)

        ax.set_xlim(0,300)
        ylims = ax.get_ylim()
        if zoom_ymax is None:
            ymax = max(ttab[models_to_plot[-1]][ttab["Age(Myr)"]<350])
        else:
            ymax = zoom_ymax
        ax.set_ylim(ylims[0],ymax)
        fig.savefig(f"plots/{outfilename}_zoom.png",bbox_inches="tight",dpi=600)

    if run_fits:
        ttab.write(f"tables/{outfilename}.csv",delimiter=",",overwrite=True)



def run_model_binned(model_name,max_q=0,include_blends=True,
                     include_lit=False,period_scale = "linear",
                     output_filebase="tausq_ZAMS_Compare",
                     zoom_ymax=None):
    pmd = PeriodMassDistribution(max_q,include_blends,include_lit)

    output_filebase = f"{output_filebase}_{model_name}"
    model = model_name

    mass_bins = np.arange(0.05,1.4,0.1)
    mass_labels = [f"{mass:.2f} Msun" for mass in mass_bins]
    nbins = len(mass_labels)-1

    # # Check for the matching output csv and skip straight to plotting if found
    outfilename = f"{output_filebase}_binned_{pmd.param_string}"
    # if os.path.exists(f"tables/{outfilename}.csv"):
    #     print("computation already completed")
    #     run_fits = False
    #     ttab = at.read(f"tables/{outfilename}.csv")
    # else:
    run_fits = True
    colnames = np.append(["Age(Myr)"],mass_labels[:-1])
    ttab = Table(np.zeros((nbins+1)*nage).reshape(nage,(nbins+1)),
                 names=colnames)
    print(ttab.dtype)


    # Determine cpu count for parallelization
    try:
        # The number of actually available CPUs, important for HPC
        cpu_count = len(os.sched_getaffinity(0))
    except AttributeError:
        # On a regular computer, just use the total CPU count
        cpu_count = mp.cpu_count()-2


    # Run comparison to data for every model
    if run_fits:

        if "WideHat" in model:
            init_type="kde"
        else:
            init_type="cluster"

        models = glob.glob(os.path.expanduser(f"~/Dropbox/Models/{model}/{model}*Myr.txt"))
        # print(models)

        model_ages = np.sort([int(mod.split("_")[-1][:5]) for mod in models])
        # print(model_ages)

        model_ages = model_ages

        pool = mp.Pool(cpu_count)
        chunksize = len(model_ages) // cpu_count

        print(f"{cpu_count} CPUs, {len(model_ages)} models, chunk size {chunksize}")
        print("starting multiprocessing run",model)

        tau_sq_args = [[age,pmd,model,period_scale,init_type] for
                        age in model_ages]

        all_tau_sq0 = pool.starmap(run_one_model_binned,tau_sq_args,
                                      chunksize=chunksize)

        all_tau_sq = np.asarray(all_tau_sq0).T
        print(np.shape(all_tau_sq))

        for i, mass_label in enumerate(mass_labels[:-1]):
            ttab[mass_label] = all_tau_sq[i]

        # ttab[model] = all_tau_sq
        #
        ttab["Age(Myr)"] = model_ages

    # If the comparison was already run, just re-plot
    else:
        for j,model in enumerate(models_to_plot):
            age_colname = f"Age_{model}"
            if "UpSco" in model:
                ls = "--"
            else:
                ls = "-"
            ax.plot(ttab[age_colname],ttab[model],ls,
                    label=display_names[model],
                    color=mapper.to_rgba((j % 3)+1),alpha=0.75)
        # Plot the models from the saved file

    if "UpSco" in model:
        ls = "--"
    else:
        ls = "-"

    # print(nbins,mass_labels)

    # Unnormalized tau_sq
    fig = plt.figure()
    fig.patch.set_facecolor('w')
    fig.patch.set_alpha(1.0)
    ax = plt.subplot(111)

    for i,mass in enumerate(mass_bins[:-1]):
        ax.plot(model_ages,all_tau_sq[i],ls,label=mass_labels[i],
                color=mapper2.to_rgba((i % nbins)+1),alpha=0.75)

    ax.legend(loc=2,ncol=3)
    ax.set_xlabel("Model age (Myr)",fontsize=16)
    ax.set_ylabel("tau squared",fontsize=16)

    ax.tick_params(labelsize=12)
    ax.set_xticks(np.arange(0,300,25),minor=True)

    fig.savefig(f"plots/{outfilename}.png",bbox_inches="tight",dpi=600)

    # Normalized tau_sq
    fig = plt.figure()
    fig.patch.set_facecolor('w')
    fig.patch.set_alpha(1.0)
    ax = plt.subplot(111)

    for i,mass in enumerate(mass_bins[:-1]):
        ax.plot(model_ages,all_tau_sq[i]/min(all_tau_sq[i]),ls,
                label=mass_labels[i],
                color=mapper2.to_rgba((i % nbins)+1),alpha=0.75)
        min_i = np.argmin(all_tau_sq[i])
        ax.plot(model_ages[min_i],0.99,"^",color=mapper2.to_rgba((i % nbins)+1))
        print(mass,model_ages[min_i])

    ax.legend(loc=2,ncol=3)
    ax.set_xlabel("Model age (Myr)",fontsize=16)
    ax.set_ylabel("tau squared (minimum = 1)",fontsize=16)

    ax.tick_params(labelsize=12)
    ax.set_xticks(np.arange(0,300,25),minor=True)

    fig.savefig(f"plots/{outfilename}_normalized.png",bbox_inches="tight",dpi=600)


    ax.set_xlim(0,300)
    if zoom_ymax is not None:
        ymax = zoom_ymax
        ax.set_ylim(0.99,ymax)
    fig.savefig(f"plots/{outfilename}_zoom.png",bbox_inches="tight",dpi=600)
    #
    #
    ttab.write(f"tables/{outfilename}.csv",delimiter=",",overwrite=True)


if __name__=="__main__":
    from argparse import ArgumentParser
    import yaml
    # import logging

    # Define parser object
    parser = ArgumentParser(description="")
    parser.add_argument("-c", "--config", dest="config_file", required=True,
                        type=str, help="Path to config file that specifies the "
                                       "parameters for this run.")

    args = parser.parse_args()

    run_all_models(args.config_file)
