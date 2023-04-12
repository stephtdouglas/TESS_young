import os, sys, glob, pathlib

import numpy as np
import astropy.io.ascii as at
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

import tess_young
from tess_young.get_const import *
_DIR = pathlib.Path(tess_young.__file__).resolve().parent.parent
plt.style.use(os.path.join(_DIR,'paper.mplstyle'))

from periodmass import PeriodMassDistribution
from spinmodel import SpinModel


def plot_all_models_yaml(config_file, multipanel=False):
    """

    """
    # parse config file
    config_file = os.path.abspath(os.path.expanduser(config_file))
    with open(config_file, 'r') as f:
        config = yaml.load(f.read())
        # config['config_file'] = config_file

    print(config)
    print(multipanel)
    name = config.pop("name")
    _ = config.pop("zoom_ymax")
    _ = config.pop("to_plot")
    _ = config.pop("overwrite")
    
    if multipanel:
        plot_multipanel(**config
                       )

    else:
        plot_all_models(**config
                       )


def plot_all_models(max_q=0,include_blends=True,include_lit=False,
                    period_scale = "linear",models_to_plot=None,
                    init_types=None,
                    output_filebase="tausq_ZAMS_Compare",mass_limits=None):
    pmd = PeriodMassDistribution(max_q,include_blends,include_lit,mass_limits)

    ymin, ymax = 0,14

    if mass_limits is not None:
        pmd_all = PeriodMassDistribution(max_q,include_blends,include_lit)

    # outfilename = f"{output_filebase}_{pmd.param_string}"
    # outfilepath = os.path.join(_DIR,f"tables/{outfilename}.csv")
    # if os.path.exists(outfilepath)==False:
    #     print("No matching output file found.")
    #     print("Please run tau-squared fits first")
    #     sys.exit()
    # else:
    #     results = at.read(outfilepath)

    plt.figure()
    for j,model in enumerate(models_to_plot):
        init_type = init_types[j]
        print(model,init_type)

        age_col = f"Age_{model}_{init_type}"
        colname = f"{model}_{init_type}"
        if init_type=="kde":
            ls = "--"
        else:
            ls = "-"

        models = glob.glob(os.path.join(MODEL_DIR,f"{model}/{model}*Myr.txt"))

        model_ages = np.sort([int(mod.split("_")[-1][:5]) for mod in models])
        # print(model_ages)

        model_ages = model_ages[(model_ages<=250) & (model_ages>=0)]

        for i, age in enumerate(model_ages):
            # if results["age_col"][i]!=age:
            #     print("ERROR: table and file ages do not match")
            #     print(results["age_col"][i],age)
            #     sys.exit()
            sm = SpinModel(model,age,period_scale,init_type=init_type)
            # if j==0 and i==0:
            #     pmd.select_obs(sm)

            # Normalize the model and calculate tau-squared
            if init_type!="kde":
                sm.normalize()
            sm.add_mask()
            # sm.calc_tau_sq(pmd)


            # Plot
            ax = sm.plot_hist()
            ax.set_ylim(ymin,ymax)

            if mass_limits is None:
                pmd.plot_obs(ax)
            else:
                pmd_all.plot_obs(ax)
                if mass_limits[0]>0.05:
                    rect_low = Rectangle((0.05,ymin),width=mass_limits[0]-0.05,height=ymax-ymin,
                                         alpha=1,color="w",zorder=10)
                    ax.add_patch(rect_low)
                if mass_limits[1]<1.3:
                    rect_high = Rectangle((mass_limits[1],ymin),width=1.3-mass_limits[1],
                                          height=ymax-ymin,
                                          alpha=1,color="w",zorder=10)
                    ax.add_patch(rect_high)

            plt.savefig(os.path.join(_DIR,f"plots/model_frames/tausq_{output_filebase}_{colname}_{pmd.param_string}_{period_scale}_{age:05d}Myr_ZAMS.png"),bbox_inches="tight",dpi=600)
            plt.close()


def plot_multipanel(max_q=0,include_blends=True,include_lit=False,
                    period_scale = "linear",output_filebase="tausq_panel_ZAMS_Compare",
                    models_to_plot=model_names,init_types=None,
                    mass_limits=None):
    """
    Create a two-panel plot showing the model/data and the point on the tau-squared vs. age plot
    """

    pmd = PeriodMassDistribution(max_q,include_blends,include_lit,mass_limits)

    ymin, ymax = 0,14

    if mass_limits is not None:
        pmd_all = PeriodMassDistribution(max_q,include_blends,include_lit)

    outfilename = f"{output_filebase}_{pmd.param_string}"
    outfilepath = os.path.join(_DIR,f"tables/{outfilename}.csv")
    if os.path.exists(outfilepath)==False:
        print("No matching output file found.")
        print("Please run tau-squared fits first")
        sys.exit()
    else:
        results = at.read(outfilepath)

    plot_dir = os.path.join(_DIR,"plots/model_frames/")
    if os.path.exists(plot_dir) is False:
        os.mkdir(plot_dir)

    for j,model in enumerate(models_to_plot):

        init_type = init_types[j]
        print(model,init_type)

        fig, axes = plt.subplots(ncols=2,figsize=(10,5))
        fig.patch.set_facecolor('w')
        fig.patch.set_alpha(1.0)
        ax1 = axes[0]
        ax2 = axes[1]


        age_col = f"Age_{model}_{init_type}"
        colname = f"{model}_{init_type}"
        if init_type=="kde":
            ls = "--"
        else:
            ls = "-"


        models = glob.glob(os.path.join(MODEL_DIR,f"{model}/{model}*Myr.txt"))

        model_ages = np.sort([int(mod.split("_")[-1][:5]) for mod in models])
        # print(model_ages)

        max_age = 250
        model_ages = model_ages[(model_ages<=max_age) & (model_ages>=0)]

        cj = np.where(model==model_names)[0][0]
        mcolor=mapper.to_rgba((cj % 3)+1)


        ax1.plot(results[age_col],results[colname],'-',color=mcolor,alpha=0.5)
        ax1.set_xlabel("Model age (Myr)",fontsize=16)
        ax1.set_ylabel("tau squared",fontsize=16)
        ax1.set_xlim(0,max_age)
        ax1.set_xticks(np.arange(0,max_age,25),minor=True)


        for i, age in enumerate(model_ages):
        #     print("\n",age)
            if results[age_col][i]!=age:
                print("ERROR: table and file ages do not match")
                print(results["age_col"][i],age)
                sys.exit()

            sm = SpinModel(model,age,period_scale,init_type=init_type)
            if j==0 and i==0:
                pmd.select_obs(sm)

            # Normalize the model and calculate tau-squared
            if init_type!="kde":
                sm.normalize()
            sm.add_mask()
            # sm.calc_tau_sq(pmd)

    #         print(model,age,np.max(sm.img))
    #         print(sm.img)

            # Add this age's dot to the tausquared histogram
            ax1.plot(results[age_col][i],results[colname][i],'.',color=mcolor,alpha=1)


            # Plot
            _ = sm.plot_hist(ax2)
            ax2.set_xlim(1.3,0.1)
            ax2.set_xlabel(r"Mass (M$_\odot$)")
            ax2.set_ylabel("Period (d)")
            ax2.patch.set_facecolor('w')
            ax2.patch.set_alpha(1.0)
            if mass_limits is None:
                pmd.plot_obs(ax2)
            else:
                pmd_all.plot_obs(ax2)
                if mass_limits[0]>0.05:
                    rect_low = Rectangle((0.05,ymin),width=mass_limits[0]-0.05,height=ymax-ymin,
                                         alpha=1,color="w",zorder=10)
                    ax2.add_patch(rect_low)
                if mass_limits[1]<1.3:
                    rect_high = Rectangle((mass_limits[1],ymin),width=1.3-mass_limits[1],
                                          height=ymax-ymin,
                                          alpha=1,color="w",zorder=10)
                    ax2.add_patch(rect_high)

            if period_scale=="log":
                ax2.set_ylim(0.1,20)
            else:
                ax2.set_ylim(0,14)
            
            plt.savefig(os.path.join(_DIR,f"plots/model_frames/tausq_panel_{model}_{init_type}_{pmd.param_string}_{period_scale}_{age:05d}Myr_ZAMS.png"),bbox_inches="tight",dpi=600)
            ax2.cla()

        #     if i>5:
        #         break
        # break


def plot_tausq_tracks(ttab,models_to_plot=None,ax=None,
                      output_filebase="tausq_tracks"):

    outfilename = output_filebase
    if ax is None:
        new_fig = True
        fig = plt.figure()
        fig.patch.set_facecolor('w')
        fig.patch.set_alpha(1.0)
        ax = plt.subplot(111)
    else:
        new_fig = False

    if models_to_plot is None:
        models_to_plot = []
        for colname in ttab.columns:
            if ("Age" in colname)==False:
                models_to_plot.append(colname)

    for j, model in enumerate(models_to_plot):
        # print(j,model)
        age_colname = f"Age_{model}"
        if "UpSco" in model:
            ls = "--"
        else:
            ls = "-"

        cj = np.where(model==model_names)[0][0]
        ax.plot(ttab[age_colname],ttab[model],ls,
                label=display_names[model],
                color=mapper.to_rgba((cj % 3)+1),alpha=0.75)

    if new_fig:
        ax.legend(loc=2)
        ax.set_xlabel("Model age (Myr)",fontsize=16)
        ax.set_ylabel("tau squared",fontsize=16)

        ax.tick_params(labelsize=12)
        ax.set_xticks(np.arange(0,300,25),minor=True)


if __name__=="__main__":
    from argparse import ArgumentParser
    import yaml
    # import logging

    # Define parser object
    parser = ArgumentParser(description="")
    parser.add_argument("-c", "--config", dest="config_file", required=True,
                        type=str, help="Path to config file that specifies the "
                                       "parameters for this run.")
    parser.add_argument("-m", "--multi", dest="to_plot_multipanel", required=False,
                        action="store_true")

    args = parser.parse_args()

    print(args)

    plot_all_models_yaml(args.config_file,args.to_plot_multipanel)

