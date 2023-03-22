import os, glob, pathlib
import multiprocessing as mp

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.cm as cm
import astropy.io.ascii as at
from astropy.table import Table

import tess_young
from tess_young.get_const import *
_DIR = pathlib.Path(tess_young.__file__).resolve().parent.parent
plt.style.use(os.path.join(_DIR,'paper.mplstyle'))

from periodmass import PeriodMassDistribution
from spinmodel import SpinModel

# mapper and normalization for binning by mass
norm2 = mpl.colors.Normalize(vmin=0, vmax=14)
mapper2 = cm.ScalarMappable(norm=norm2, cmap=cm.viridis)

def run_one_model(age,pmd,model,period_scale,init_type):
    """
    Generate a SpinModel and calculate tau-squared from the input PeriodMassDistribution
    """
    sm = SpinModel(model,age,period_scale,init_type=init_type)

    # Normalize the model and calculate tau-squared
    if init_type!="kde":
        sm.normalize()
    sm.add_mask()
    sm.calc_tau_sq(pmd)

    return sm.tau_sq

def run_one_model_binned(age,pmd,model,period_scale,init_type):
    """
    Generate a SpinModel and calculate tau-squared for each individual mass bin
    """
    sm = SpinModel(model,age,period_scale,init_type=init_type)

    # Normalize the model and calculate tau-squared
    if init_type!="kde":
        sm.normalize()
    sm.add_mask()
    sm.calc_tau_sq_binned(pmd)

    return sm.tau_sq

def run_all_models_yaml(config_file):
    """

    """
    # parse config file
    config_file = os.path.abspath(os.path.expanduser(config_file))
    with open(config_file, 'r') as f:
        config = yaml.load(f.read())
        config['config_file'] = config_file

    print(config)

    run_all_models(max_q=config["max_q"],
                   include_blends=config["include_blends"],
                   include_lit=config["include_lit"],
                   period_scale="linear",
                   output_filebase=config["output_filebase"],
                   models_to_plot=config["models"],
                   zoom_ymax=config["zoom_ymax"],
                   mass_limits=config["mass_limits"],
                   pmd=None,
                   to_plot=True,
                   overwrite=config["overwrite"],
                   init_types=config["init_types"])


def run_all_models(max_q=0,include_blends=True,include_lit=False,
                   period_scale="linear",output_filebase="tausq_ZAMS_Compare",
                   models_to_plot=model_names,zoom_ymax=None,
                   mass_limits=None,pmd=None,to_plot=True,overwrite=False, 
                   init_types=None):
    """
    Compare a single period-mass distribution to multiple sets of models.

    Inputs
    ------
    max_q: integer, maximum quality flag to include (should be 0 or 1)
    include_blends: boolean, whether or not to include potentially blended targets
    include_lit: boolean, whether or not to include literature values
    period_scale: (string) "log" or "linear"
    output_filebase: (string) identifier for this run, will be included in output filenames
    models_to_plot: (list of strings) identifiers for the torque laws/models to be used. Options:
                         "UpSco_Mattea2015", "UpSco_Mattea2022",
                         "UpSco_"ZeroTorque", "WideHat8Myr_Mattea2015",
                         "WideHat8Myr_Mattea2022", "WideHat8Myr_ZeroTorque"
    zoom_ymax: (number, optional) if included, set the maximum for the tau-squared axis
    mass_limits: (tuple or list, optional) minimum and maximum masses to include
    pmd: (PeriodMass object, optional) if included, will override max_q, include_blends, and include_lit
    to_plot: (boolean), whether to generate output plots. Default=True
    overwrite: (boolean), whether to overwrite existing table files. Default=False. 

    """
    

    nmod_l = len(models_to_plot)

    if pmd is None:
        pmd = PeriodMassDistribution(max_q,include_blends,include_lit,
                                     mass_limits=mass_limits)
    else:
        print("WARNING: Using input period-mass distribution.")
        print("Ignoring input q, include_*, and scale match.")

 
    # output_filebase = f"{output_filebase}_{model_name}_{init_type}"

    # Check for the matching output csv and skip straight to plotting if found
    outfilename = f"{output_filebase}_{pmd.param_string}"
    outfilepath = os.path.join(_DIR,f"tables/{outfilename}.csv")
    outplotpath = os.path.join(_DIR,f"plots/{outfilename}.png")

    if (init_types is None):
        init_types = np.zeros(nmod_l,"U8")
        for j, model in enumerate(models_to_plot):
            if("WideHat" in model):
                init_types[j] = "tophat"
            elif ("UpSco" in model):
                init_type="cluster"
            else:
                print("ERROR: Unknown model, ", model)
                print("requires init_type to be specified")
                return

    if (overwrite==False) and (os.path.exists(outfilepath)):
        print("computation already completed")
        run_fits = False
        if to_plot:
            ttab = at.read(outfilepath)
        else:
            return
    else:
        run_fits = True
        colnames = [f"{models_to_plot[i]}_{init_types[i]}" for i in range(nmod_l)]
        print("Cols",colnames)
        ttab = Table(np.zeros(nmod_l*nage).reshape(nage,nmod_l),names=colnames)

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
            init_type = init_types[j]

            print(model, init_type)

            models = glob.glob(os.path.join(MODEL_DIR,f"{model}/{model}*Myr.txt"))
            # print(models)

            model_ages = np.sort([int(mod.split("_")[-1][:5]) for mod in models])
            # print(model_ages)

            model_ages = model_ages#[model_ages<=300]
            age_colname = f"Age_{model}_{init_type}"
            ttab[age_colname] = model_ages

            all_tau_sq = np.zeros(len(model_ages))

            pool = mp.Pool(cpu_count)
            chunksize = len(model_ages) // cpu_count

            print(f"{cpu_count} CPUs, chunk size {chunksize}")
            print("starting multiprocessing run",model, init_type)

            tau_sq_args = [[age,pmd,model,period_scale,init_type] for
                            age in model_ages]

            all_tau_sq = pool.starmap(run_one_model,tau_sq_args,
                                          chunksize=chunksize)

            if init_type=="kde":
                ls = "--"
            else:
                ls = "-"

            if to_plot:
                ax.plot(model_ages,all_tau_sq,ls,label=display_names[model]+init_names[init_type],
                        color=mapper.to_rgba((j % 3)+1),alpha=0.75)
            colname = f"{model}_{init_type}"
            ttab[colname] = all_tau_sq

            if j==0:
                ttab["Age(Myr)"] = model_ages

    # If the comparison was already run, just re-plot
    else:
        for j,model in enumerate(models_to_plot):
            init_type = init_types[j]

            print(model, init_type)

            age_colname = f"Age_{model}_{init_type}"
            colname = f"{model}_{init_type}"
            if init_type=="kde":
                ls = "--"
            else:
                ls = "-"

            if to_plot:
                ax.plot(ttab[age_colname],ttab[colname],ls,
                        label=display_names[model]+init_names[init_type],
                        color=mapper.to_rgba((j % 3)+1),alpha=0.75)
        # Plot the models from the saved file

    if to_plot:
        ax.legend(loc=2)
        ax.set_xlabel("Model age (Myr)",fontsize=16)
        ax.set_ylabel("tau squared",fontsize=16)

        ax.tick_params(labelsize=12)
        ax.set_xticks(np.arange(0,300,25),minor=True)
        ax.set_title(outfilename)

        fig.savefig(outplotpath,bbox_inches="tight",dpi=600)

        ax.set_xlim(0,300)
        ylims = ax.get_ylim()
        # colname = f"{models_to_plot[-1]}_{init_types[-1]}"
        if zoom_ymax is not None:
            ymax = zoom_ymax
            ax.set_ylim(ylims[0],ymax)
        ax.set_title(outfilename)
        fig.savefig(outplotpath.replace(".png","_zoom.png"),bbox_inches="tight",dpi=600)

    if run_fits:
        ttab.write(outfilepath,delimiter=",",overwrite=True)



def run_model_binned(model_name,max_q=0,include_blends=True,
                     include_lit=False,period_scale = "linear",
                     output_filebase="tausq_ZAMS_Compare",
                     zoom_ymax=None, init_type=None):
    pmd = PeriodMassDistribution(max_q,include_blends,include_lit)
    
    model = model_name
    if init_type is None:
            if("WideHat" in model):
                init_type = "tophat"
            elif ("UpSco" in model):
                init_type="cluster"
            else:
                print("ERROR: Unknown model, ", model)
                print("requires init_type to be specified")
                return


    output_filebase = f"{output_filebase}_{model_name}_{init_type}"

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


        models = glob.glob(os.path.join(MODEL_DIR,f"{model}/{model}*Myr.txt"))
        # print(models)

        model_ages = np.sort([int(mod.split("_")[-1][:5]) for mod in models])
        # print(model_ages)

        model_ages = model_ages

        pool = mp.Pool(cpu_count)
        chunksize = len(model_ages) // cpu_count

        print(f"{cpu_count} CPUs, {len(model_ages)} models, chunk size {chunksize}")
        print("starting multiprocessing run",model, init_type)

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

    # # If the comparison was already run, just re-plot
    # else:
    #     for j,model in enumerate(models_to_plot):
    #         age_colname = f"Age_{model}"
    #         if "UpSco" in model:
    #             ls = "--"
    #         else:
    #             ls = "-"
    #         ax.plot(ttab[age_colname],ttab[model],ls,
    #                 label=display_names[model],
    #                 color=mapper.to_rgba((j % 3)+1),alpha=0.75)
    #     # Plot the models from the saved file

    # if init_type=="kde":
    #     ls = "--"
    # else:
    #     ls = "-"

    # print(nbins,mass_labels)

    # # Unnormalized tau_sq
    # fig = plt.figure()
    # fig.patch.set_facecolor('w')
    # fig.patch.set_alpha(1.0)
    # ax = plt.subplot(111)

    # for i,mass in enumerate(mass_bins[:-1]):
    #     ax.plot(model_ages,all_tau_sq[i],ls,label=mass_labels[i],
    #             color=mapper2.to_rgba((i % nbins)+1),alpha=0.75)

    # ax.legend(loc=2,ncol=3)
    # ax.set_xlabel("Model age (Myr)",fontsize=16)
    # ax.set_ylabel("tau squared",fontsize=16)

    # ax.tick_params(labelsize=12)
    # ax.set_xticks(np.arange(0,300,25),minor=True)

    # fig.savefig(os.path.join(_DIR,f"plots/{outfilename}.png"),bbox_inches="tight",dpi=600)

    # # Normalized tau_sq
    # fig = plt.figure()
    # fig.patch.set_facecolor('w')
    # fig.patch.set_alpha(1.0)
    # ax = plt.subplot(111)

    # for i,mass in enumerate(mass_bins[:-1]):
    #     ax.plot(model_ages,all_tau_sq[i]/min(all_tau_sq[i]),ls,
    #             label=mass_labels[i],
    #             color=mapper2.to_rgba((i % nbins)+1),alpha=0.75)
    #     min_i = np.argmin(all_tau_sq[i])
    #     ax.plot(model_ages[min_i],0.99,"^",color=mapper2.to_rgba((i % nbins)+1))
    #     print(mass,model_ages[min_i])

    # ax.legend(loc=2,ncol=3)
    # ax.set_xlabel("Model age (Myr)",fontsize=16)
    # ax.set_ylabel("tau squared (minimum = 1)",fontsize=16)

    # ax.tick_params(labelsize=12)
    # ax.set_xticks(np.arange(0,300,25),minor=True)

    # fig.savefig(os.path.join(_DIR,f"plots/{outfilename}_normalized.png"),bbox_inches="tight",dpi=600)


    # ax.set_xlim(0,300)
    # if zoom_ymax is not None:
    #     ymax = zoom_ymax
    #     ax.set_ylim(0.99,ymax)
    # fig.savefig(os.path.join(_DIR,f"plots/{outfilename}_zoom.png"),bbox_inches="tight",dpi=600)
    # #
    # #
    ttab.write(os.path.join(_DIR,f"tables/{outfilename}.csv"),delimiter=",",overwrite=True)


if __name__=="__main__":
    from argparse import ArgumentParser
    import yaml
    # import logging

    # Define parser object
    parser = ArgumentParser(description="")
    c_group = parser.add_mutually_exclusive_group(required=True)
    c_group.add_argument("-c", "--config", dest="config_file", #required=True,
                        type=str, help="Path to config file that specifies the "
                                       "parameters for this run.")
    c_group.add_argument("-b", "--binned", dest="binned_config", #required=True,
                        type=str, help="Path to config file that specifies the "
                                       "run parameters (binned by mass)")
    # parser.add_argument("-m", "--style", dest="style_file", required=False,
    #                     type=str, help="Path to matplotlib style file")

    args = parser.parse_args()

    # if os.path.exists(args.style_file):
    #     plt.style.use(args.style_file)

    if args.config_file is not None:
        run_all_models_yaml(args.config_file)
    else:
        # parse config file
        config_file = os.path.abspath(os.path.expanduser(args.binned_config))
        with open(config_file, 'r') as f:
            config = yaml.load(f.read())
            config['config_file'] = config_file

        print(config)

        nmod_l = len(config["models"])
        for i in range(nmod_l):
            model_name = config["models"][i]
            init_type = config["init_types"][i]
            run_model_binned(model_name,max_q=config["max_q"],
                           include_blends=config["include_blends"],
                           include_lit=config["include_lit"],
                           period_scale="linear",
                           init_type=init_type,
                           output_filebase=config["output_filebase"]
                           )

