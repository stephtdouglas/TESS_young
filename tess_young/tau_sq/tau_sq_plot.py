import os, sys, glob, pathlib

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

import tess_young
from tess_young.get_const import *
_DIR = pathlib.Path(tess_young.__file__).resolve().parent.parent
plt.style.use(os.path.join(_DIR,'paper.mplstyle'))

from periodmass import PeriodMassDistribution
from spinmodel import SpinModel

def plot_all_models_yaml(config_file):
    """

    """
    # parse config file
    config_file = os.path.abspath(os.path.expanduser(config_file))
    with open(config_file, 'r') as f:
        config = yaml.load(f.read())
        config['config_file'] = config_file

    print(config)

    plot_all_models(max_q=config["max_q"],
                   include_blends=config["include_blends"],
                   include_lit=config["include_lit"],
                   period_scale="linear",
                   output_filebase=config["output_filebase"],
                   models_to_plot=config["models"],
                   mass_limits=config["mass_limits"]
                   )


def plot_all_models(max_q=0,include_blends=True,include_lit=False,
                    period_scale = "linear",models_to_plot=model_names,
                    output_filebase="tausq_ZAMS_Compare",mass_limits=None):
    pmd = PeriodMassDistribution(max_q,include_blends,include_lit,mass_limits)

    ymin, ymax = 0,14

    if mass_limits is not None:
        pmd_all = PeriodMassDistribution(max_q,include_blends,include_lit)
        if mass_limits[0]>0.05:
            rect_low = Rectangle((0.05,ymin),width=mass_limits[0]-0.05,height=ymax-ymin,
                                 alpha=0.5,color="w",zorder=10)
        else:
            rect_low = None
        if mass_limits[1]<1.3:
            rect_high = Rectangle((mass_limits[1],ymin),width=1.3-mass_limits[1],
                                  height=ymax-ymin,
                                  alpha=0.5,color="w",zorder=10)
        else:
            rect_high = None
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
        print(model)

        if "WideHat" in model:
            init_type="kde"
        else:
            init_type="cluster"

        models = glob.glob(os.path.join(MODEL_DIR,f"{model}/{model}*Myr.txt"))

        model_ages = np.sort([int(mod.split("_")[-1][:5]) for mod in models])
        # print(model_ages)

        model_ages = model_ages[(model_ages<=150) & (model_ages>=0)]

        for i, age in enumerate(model_ages):
            # if results["age_col"][i]!=age:
            #     print("ERROR: table and file ages do not match")
            #     print(results["age_col"][i],age)
            #     sys.exit()
            sm = SpinModel(model,age,period_scale,init_type=init_type)
            # if j==0 and i==0:
            #     pmd.select_obs(sm)

            # Normalize the model and calculate tau-squared
            if ("WideHat" in model)==False:
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

            plt.savefig(os.path.join(_DIR,f"plots/model_frames/tausq_{output_filebase}_{model}_{pmd.param_string}_{period_scale}_{age:05d}Myr_ZAMS.png"),bbox_inches="tight",dpi=600)
            plt.close()


# def plot_multipanel(max_q=0,include_blends=True,include_lit=False,
#                     period_scale = "linear",models_to_plot=model_names,
#                     mass_limits=None):
#     """
#     Create a two-panel plot showing the model/data and the point on the tau-squared vs. age plot
#     """

#     pmd = PeriodMassDistribution(max_q,include_blends,include_lit,mass_limits)
#     outfilename = f"{output_filebase}_{pmd.param_string}"
#     outfilepath = os.path.join(_DIR,f"tables/{outfilename}.csv")
#     if os.path.exists(outfilepath)==False:
#         print("No matching output file found.")
#         print("Please run tau-squared fits first")
#         sys.exit()
#     else:
#         results = at.read(outfilepath)


#     fix, axes = plt.subplots(ncols=2,figsize=(10,5))

#     for j,model in enumerate(models_to_plot):
#         print(model)

#         if "WideHat" in model:
#             init_type="kde"
#         else:
#             init_type="cluster"

#         models = glob.glob(os.path.join(MODEL_DIR,f"{model}/{model}*Myr.txt"))

#         model_ages = np.sort([int(mod.split("_")[-1][:5]) for mod in models])
#         # print(model_ages)

#         model_ages = model_ages[(model_ages<=150) & (model_ages>=0)]

#         age_col = f"Age_{model}"

#         for i, age in enumerate(model_ages):
#         #     print("\n",age)
#             if results["age_col"][i]!=age:
#                 print("ERROR: table and file ages do not match")
#                 print(results["age_col"][i],age)
#                 sys.exit()

#             sm = SpinModel(model,age,period_scale,init_type=init_type)
#             if j==0 and i==0:
#                 pmd.select_obs(sm)

#             # Normalize the model and calculate tau-squared
#             if ("WideHat" in model)==False:
#                 sm.normalize()
#             sm.add_mask()
#             sm.calc_tau_sq(pmd)

#     #         print(model,age,np.max(sm.img))
#     #         print(sm.img)


#             # Plot
#             ax = sm.plot_hist()
#             pmd.plot_obs(ax)
#             ax.set_ylim(0,14)
#             plt.savefig(os.path.join(_DIR,f"plots/model_frames/tausq2_{model}_{pmd.param_string}_{period_scale}_{age:05d}Myr_ZAMS.png"),bbox_inches="tight",dpi=600)
#             plt.close()


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
                color=mapper.to_rgba((cj % 3)+1),alpha=0.5)

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
    # parser.add_argument("-m", "--style", dest="style_file", required=False,
    #                     type=str, help="Path to matplotlib style file")

    args = parser.parse_args()

    # if os.path.exists(args.style_file):
    #     plt.style.use(args.style_file)

    plot_all_models

    plot_all_models_yaml(args.config_file)

