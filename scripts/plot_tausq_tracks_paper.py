import os, glob, pathlib
import multiprocessing as mp

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.cm as cm
from matplotlib.patches import Rectangle
import astropy.io.ascii as at
from astropy.table import Table
import yaml

import tess_young
from tess_young.get_const import *
_DIR = pathlib.Path(tess_young.__file__).resolve().parent.parent
plt.style.use(os.path.join(_DIR,'paper.mplstyle'))

from tess_young.tau_sq.periodmass import PeriodMassDistribution

def plot_all_tracks(pmd=None,
                   period_scale="linear",output_filebase="tausq_ZAMS_Compare",
                   models_to_plot=model_names,zoom_ymax=None,
                   mass_limits=None,to_plot=True,overwrite=False, 
                   init_types=None,cluster="all",
                   max_q=0,include_blends=True,include_lit=False,
                   config_num=None, ax=None):
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
                                     mass_limits=mass_limits,cluster=cluster)
    else:
        print("WARNING: Using input period-mass distribution.")
        print("Ignoring input q, include_*, scale, mass_limits, and cluster.")

 
    # output_filebase = f"{output_filebase}_{model_name}_{init_type}"

    # Check for the matching output csv and skip straight to plotting if found
    outfilename = f"{output_filebase}_{pmd.param_string}"
    outfilepath = os.path.join(_DIR,f"tables/MINESweeper_v7/{outfilename}.csv")
    outplotpath = os.path.join(_DIR,f"plots/MINESweeper_v7/{outfilename}.png")

    if (init_types is None):
        init_types = np.zeros(nmod_l,"U8")
        for j, model in enumerate(models_to_plot):
            if("WideHat" in model):
                init_types[j] = "tophat"
            elif ("UpSco" in model):
                init_types[j]="cluster"
            else:
                print("ERROR: Unknown model, ", model)
                print("requires init_type to be specified")
                return

    if (os.path.exists(outfilepath)) is False:
        print("computation not completed; please run tau_sq_run.py\n",outfilepath)
        return 

    ttab = at.read(outfilepath)

    # Set up figure
    if ax is None:
        fig = plt.figure()
        fig.patch.set_facecolor('w')
        fig.patch.set_alpha(1.0)
        ax = plt.subplot(111)
        ax.set_xlabel("Model age (Myr)")
        ax.set_ylabel(r"$\tau^2$",fontsize=20)

        ax.set_xticks(np.arange(0,300,25),minor=True)
        ax.set_title(outfilename)
        ax.set_xlim(0,300)
        new_plot = True
    else:
        new_plot=False


    for j,model in enumerate(models_to_plot):

        init_type = init_types[j]
        if init_type.lower()!="kde":
            continue

        print(model, init_type)

        age_colname = f"Age_{model}_{init_type}"
        colname = f"{model}_{init_type}"
        
        if "Zero" in model:
            ls = ":"
        elif "2015" in model:
            ls = "-."
        else:
            ls = "-"

        ax.plot(ttab[age_colname],ttab[colname],ls,
                label=display_names[model],
                color=mapper.to_rgba((j % 3)+1),alpha=0.75)

        best_i = np.argmin(ttab[colname])
        best_age = ttab[age_colname][best_i]
        ax.plot(best_age,ttab[colname][best_i],'^',color=mapper.to_rgba((j % 3)+1))
        ax.text(best_age+5,ttab[colname][best_i],f"{best_age:.0f} Myr",
                horizontalalignment="left",verticalalignment="top",
                color=mapper.to_rgba((j % 3)+1))
    # Plot the models from the saved file

    if new_plot:
        ax.legend(loc=2)

        ylims = ax.get_ylim()
        # colname = f"{models_to_plot[-1]}_{init_types[-1]}"
        if zoom_ymax is None:
            ymax = ylims[1]*0.75
        else:
            ymax = zoom_ymax
        ax.set_ylim(ymin,ymax)
        ax.set_title(outfilename)

    search_string = f"tables/MINESweeper_v7/best_ages_*SYN{config_num}{period_scale}_{cluster}*kde.csv"
    age_files = glob.glob(os.path.join(_DIR,search_string))

    ylims = ax.get_ylim()
    ydiff = ylims[1]-ylims[0]
    print(ylims,ydiff)
    for filename in age_files:
        print(filename)
        dat = at.read(filename)
        colname = dat.dtype.names[0]
        print(colname)
        print(dat)

        if "Zero" in colname:
            print("Zero")
            ls = "-" #":"
            yvals = np.ones_like(dat[colname])*ylims[0]+10#+ydiff*0.9
            j=2
        elif "2015" in colname:
            print("2015")
            ls = "-" #"-."
            yvals = np.ones_like(dat[colname])*ylims[0]+10#+ydiff*0.85
            j=0
        else:
            print("2022")
            ls = "-"
            yvals = np.ones_like(dat[colname])*ylims[0]+10#+ydiff*0.8
            j=1



        # min_age = min(dat[colname])
        # age_diff = max(dat[colname])-min_age
        # ylims = ax.get_ylim()
        # ydiff = ylims[1]-ylims[0]
        # print(min_age,age_diff,ylims)
        # rect = Rectangle([min_age,ylims[0]],width=age_diff,height=ydiff,
        #                  color=mapper.to_rgba((j % 3)+1),alpha=0.75,
        #                  zorder=-1*j)

        # ax.add_patch(rect)

        ax.plot(dat[colname],yvals,ls,lw=3,color=mapper.to_rgba((j % 3)+1))

    # fig.savefig(outplotpath.replace(".png","_unc.png"),bbox_inches="tight",dpi=600)

def make_paper_plots():

    # Plotting the results of the blendsFalse, litTrue fits
    lin_config_file = os.path.join(_DIR,"config/tau_sq4.yml")
    log_config_file = os.path.join(_DIR,"config/tau_sq4_log.yml")

    fig, axes = plt.subplots(nrows=2,ncols=1,sharex=True,sharey=True)

    plt.subplots_adjust(hspace=0)
    axes[1].set_xlabel("Model age (Myr)")
    axes[1].set_ylabel(r"$\tau^2$",fontsize=20)
    axes[0].set_ylabel(r"$\tau^2$",fontsize=20)
    axes[0].set_xticks(np.arange(0,300,25),minor=True)
    axes[0].set_xlim(0,300)
    axes[0].set_ylim(1275,1650)

    axes[0].text(5,1625,"Linear fit",fontsize="large")
    axes[1].text(5,1625,"Log fit",fontsize="large")

    for j, config_file in enumerate([lin_config_file,log_config_file]):
        with open(config_file, 'r') as f:
            config = yaml.load(f.read())

        print(config)
        name = config.pop("name")

        cfilename = config_file.split("/")[-1]
        config["config_num"] = int(cfilename[6])

        _ = config.setdefault("cluster","all")

        plot_all_tracks(pmd=None,ax=axes[j],**config)

    axes[0].legend(loc=4)

    # plt.show()
    # fig.savefig(os.path.join(PAPER_DIR,"fig_tausq_tracks.pdf"),bbox_inches="tight")
    fig.savefig(os.path.join(_DIR,"plots/MINESweeper_v7/fig_tausq_tracks.png"),bbox_inches="tight",dpi=600)

if __name__=="__main__":

    make_paper_plots()
    # from argparse import ArgumentParser
    # import yaml
    # # import logging

    # # Define parser object
    # parser = ArgumentParser(description="")
    # c_group = parser.add_mutually_exclusive_group(required=True)
    # c_group.add_argument("-c", "--config", dest="config_file", #required=True,
    #                     type=str, help="Path to config file that specifies the "
    #                                    "parameters for this run.")

    # parser.add_argument("-g", dest="cluster", default="all", required=False,
    #                     help="Which group/cluster to fit (default is all)")

    # args = parser.parse_args()

    # if args.config_file is not None:
    #     # Run a regular fit
    #     config_file = os.path.abspath(os.path.expanduser(args.config_file))
    #     with open(config_file, 'r') as f:
    #         config = yaml.load(f.read())

    #     print(config)
    #     name = config.pop("name")

    #     cfilename = config_file.split("/")[-1]
    #     config["config_num"] = int(cfilename[6])

    #     _ = config.setdefault("cluster",args.cluster)

    #     if args.cluster=="UpSco":
    #         from upsco_periodmass import UpSco
    #         pmd = UpSco(mass_limits=config["mass_limits"])
    #     else:
    #         pmd = None

    #     plot_all_tracks(pmd=pmd,**config)
