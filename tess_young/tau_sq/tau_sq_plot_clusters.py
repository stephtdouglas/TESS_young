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



def plot_all_clusters(max_q=0,include_blends=True,include_lit=False,
                       period_scale="linear",output_filebase="tausq_ZAMS_Compare",
                       models_to_plot=model_names,zoom_ymax=None,
                       mass_limits=None,pmd=None,to_plot=True, 
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


    """
    

    fig, axes = plt.subplots(nrows=len(clusters),figsize=(8,10), sharex=True)
    fig.patch.set_facecolor('w')
    fig.patch.set_alpha(1.0)

    plt.subplots_adjust(hspace=0)


    param_string1 = f"Qmax{max_q}_blends{include_blends}_lit{include_lit}"
    outfilename1 = f"all_clusters_{output_filebase}_{param_string1}"
    outplotpath = os.path.join(_DIR,f"plots/{outfilename1}.png")


    for i, cluster in enumerate(clusters):
        pmd = PeriodMassDistribution(max_q,include_blends,include_lit,
                                     mass_limits=mass_limits,cluster=cluster)

        param_string = f"Qmax{max_q}_blends{include_blends}_lit{include_lit}_{cluster}"

        outfilename = f"{output_filebase}_{param_string}"
        outfilepath = os.path.join(_DIR,f"tables/{outfilename}.csv")

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

        if os.path.exists(outfilepath) is False:
            print("run computation first")
            return
        else:
            ttab = at.read(outfilepath)


        ax = axes[i]
        for j,model in enumerate(models_to_plot):
            init_type = init_types[j]
            if init_type.lower()!="kde":
                continue

            # print(model, init_type)

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
                    color=colors[cluster],alpha=0.9)
                # Plot the models from the saved file

        ylims = ax.get_ylim()
        ax.set_ylim(top=ylims[1]*0.7)

        ax.text(0.05,0.9, cluster, horizontalalignment='left',
                verticalalignment='center', transform=ax.transAxes,
                color=colors[cluster])


    axes[0].legend(loc=1)
    axes[0].set_xlim(0,150)

    axes[0].set_xlabel("Model age (Myr)",fontsize=16)
    axes[2].set_ylabel("tau squared",fontsize=16)

    axes[0].tick_params(labelsize=12)
    axes[0].set_xticks(np.arange(0,151,15),minor=True)
    axes[0].set_title(param_string1)

    fig.savefig(outplotpath,bbox_inches="tight",dpi=600)




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

    args = parser.parse_args()

    if args.config_file is not None:
        # Run a regular fit
        config_file = os.path.abspath(os.path.expanduser(args.config_file))
        with open(config_file, 'r') as f:
            config = yaml.load(f.read())

        # print(config)
        for argkey in ["name","overwrite","to_plot"]:
            _ = config.pop(argkey)

        plot_all_clusters(pmd=None,**config)
