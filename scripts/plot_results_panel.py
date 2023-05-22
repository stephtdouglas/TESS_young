import os, glob, pathlib

import yaml
import numpy as np
import astropy.io.ascii as at
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import LogNorm

import tess_young
from tess_young.get_const import *
_DIR = pathlib.Path(tess_young.__file__).resolve().parent.parent
# plt.style.use(os.path.join(_DIR,'paper.mplstyle'))

from tess_young.tau_sq.periodmass import PeriodMassDistribution
from tess_young.tau_sq.spinmodel import SpinModel

def plot_results_panel(lin_config_file,log_config_file):

    # Set up axes: 2 rows 3 cols and a colorbar
    fig, axes = plt.subplots(nrows=2,ncols=3,sharex=True,sharey="row",
                             figsize=(12,7))
    fig.subplots_adjust(left=0.1,right=0.9,hspace=0,wspace=0)

    # left bottom width height
    cbar_ax = fig.add_axes([0.9, 0.11, 0.02, 0.77])

    # Remove white space between axes and add labels
    axes[0,0].set_ylim(0.07,18)
    axes[0,0].set_xlim(1.3,0.1)
    axes[0,0].set_ylabel("Period (d)", fontsize=16)

    axes[1,0].set_ylim(0.07,18)
    axes[1,0].set_xlim(1.3,0.1)
    axes[1,0].set_ylabel("Period (d)", fontsize=16)
    axes[1,0].set_xlabel(r"Mass (M$_\odot$)", fontsize=16)
    axes[1,1].set_xlabel(r"Mass (M$_\odot$)", fontsize=16)
    axes[1,2].set_xlabel(r"Mass (M$_\odot$)", fontsize=16)


    # Define the min and max values for the colorbar
    vmin, vmax = 1e-4, 4e-1

    # in the top row: plot three models with the linear fit
    # in the bottom row: plot three models with the log fit
    # and overplot the data on both

    for j,config_file in enumerate([lin_config_file,log_config_file]):
        with open(config_file, 'r') as f:
            config = yaml.load(f.read())

        keep_keys = ["max_q","include_blends","include_lit"]
        pmd_config = {}
        for key in keep_keys:
            pmd_config[key] = config[key]

        pmd = PeriodMassDistribution(**pmd_config)

        output_filebase = config["output_filebase"]
        outfilename = f"{output_filebase}_{pmd.param_string}"
        outfilepath = os.path.join(_DIR,f"tables/{outfilename}.csv")
        dat = at.read(outfilepath)

        model = "WideHat8Myr_Mattea2022"
        init_type = "kde"

        colname = f"{model}_{init_type}"
        age_colname = "Age_"+colname
        # identify the best-fitting model age
        best_i = np.argmin(dat[colname])
        best_age = dat[age_colname][best_i]

        # Plot the initial model for reference
        ilow = np.where(dat[age_colname]==20)[0]

        # Then find a model roughly the same separation above the best-fitting age
        age_diff = dat[age_colname][best_i] - dat[age_colname][ilow]
        up_age = dat[age_colname][best_i] + age_diff
        up_i = np.argmin(np.abs(dat[age_colname]-up_age))

        for ct,i in enumerate([ilow,best_i,up_i]):

            age = int(dat[age_colname][i])
            tau2 = float(dat[colname][i])


            # Create the appropriate spinmodel
            print(model,age,config["period_scale"],init_type)
            sm = SpinModel(model,age,config["period_scale"],
                           init_type=init_type)

            ax = axes[j,ct]

            sm.add_mask()
            sm.plot_hist(ax=ax, vmin=vmin, vmax=vmax)
            print(np.percentile(sm.img[sm.mask==False],[0,5,10,99,99.5,100]))

            ax.plot(pmd.mass,pmd.prot,'w.',mec='k',zorder=15,ms=4)

            if j==1:
                va="top"
            else:
                va="bottom"
            ax.text(1.2,16,f"age = {age} Myr\n"+r"$\tau^2$"+f"= {tau2:.1f}",
                    verticalalignment=va,fontsize=8)


    # Set up the colorbar
    norm = LogNorm(vmin=vmin,vmax=vmax)
    # cbar_ax.tick_params(labelbottom=False,labelright=True)
    plt.colorbar(cm.ScalarMappable(norm=norm, cmap="Greens"),cax=cbar_ax,
                 label=r"$\rho_{number}$")

    plt.savefig(os.path.join(_DIR,"plots/panel_tausq.png"),bbox_inches="tight")
    plt.savefig(os.path.join(PAPER_DIR,"fig_panel_tausq.pdf"),bbox_inches="tight")


if __name__=="__main__":

    # from argparse import ArgumentParser
    # import yaml
    # # import logging

    # # Define parser object
    # parser = ArgumentParser(description="")
    # parser.add_argument("-c", "--config", dest="config_file", required=True,
    #                     type=str, help="Path to config file that specifies the "
    #                                    "parameters for this run.")
    # parser.add_argument("-f", "--file", dest="output_file", required=True,
    #                     type=str,help="Path to the output .csv file from this run.")
    # # TODO: need to get two output files for this - linear and log fits

    # args = parser.parse_args()


    # config_file = os.path.abspath(os.path.expanduser(args.config_file))
    # with open(config_file, 'r') as f:
    #     config = yaml.load(f.read())

    # max_q = config["max_q"]
    # include_blends = config["include_blends"]
    # include_lit = config["include_lit"]
    # output_filebase = config["output_filebase"]

    # param_string_wild = f"{output_filebase}*Qmax{max_q}_blends{include_blends}_lit{include_lit}*"

    # filename = os.path.abspath(os.path.expanduser(args.output_file))

    # dat = at.read(filename)
    # ncol = len(config["models_to_plot"])
    # colnames = dat.dtype.names[:ncol]
    # for colname in colnames:
    #     if "tophat" in colname:
    #         continue
    #     print("\t",colname,dat["Age_"+colname][np.argmin(dat[colname])])


    lin_config_file = os.path.join(_DIR,"config/tau_sq4.yml")
    log_config_file = os.path.join(_DIR,"config/tau_sq4_log.yml")
    # lin_filename = os.path.join(_DIR,"tables/")
    # log_filename = os.path.join(_DIR,"tables/")


    plot_results_panel(lin_config_file,log_config_file)
