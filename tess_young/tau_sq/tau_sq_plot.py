import os, glob

import numpy as np
import matplotlib.pyplot as plt

from tess_young.get_const import *
# plt.style.use('../../paper.mplstyle')

from periodmass import PeriodMassDistribution
from spinmodel import SpinModel

def plot_all_models(max_q=0,include_blends=True,include_lit=False,
                    period_scale = "linear",models_to_plot=model_names):
    pmd = PeriodMassDistribution(max_q,include_blends,include_lit)

    plt.figure()
#     for j,model in enumerate(model_names):
    for j,model in enumerate(models_to_plot):
        print(model)

        if "WideHat" in model:
            init_type="kde"
        else:
            init_type="cluster"

        models = glob.glob(os.path.expanduser(f"~/Dropbox/Models/{model}/{model}*Myr.txt"))

        model_ages = np.sort([int(mod.split("_")[-1][:5]) for mod in models])
        # print(model_ages)

        model_ages = model_ages[(model_ages<=150) & (model_ages>=0)]

        for i, age in enumerate(model_ages):
        #     print("\n",age)
            sm = SpinModel(model,age,period_scale,init_type=init_type)
            if j==0 and i==0:
                pmd.select_obs(sm)

            # Normalize the model and calculate tau-squared
            if ("WideHat" in model)==False:
                sm.normalize()
            sm.add_mask()
            sm.calc_tau_sq(pmd)

    #         print(model,age,np.max(sm.img))
    #         print(sm.img)


            # Plot
            ax = sm.plot_hist()
            pmd.plot_obs(ax)
            ax.set_ylim(0,14)
            plt.savefig(f"plots/model_frames/tausq_{model}_{pmd.param_string}_{period_scale}_{age:05d}Myr_ZAMS.png",bbox_inches="tight",dpi=600)
            plt.close()

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

    #### Plot model comparisons with data
    # Original run
    plot_all_models(max_q=0,models_to_plot=["WideHat8Myr_Mattea2015","WideHat8Myr_Mattea2022"])

    # Replace blends with literature
    # Only q=0
    plot_all_models(max_q=0,include_blends=False,include_lit=True,
                    period_scale = "linear",models_to_plot=["WideHat8Myr_Mattea2015","WideHat8Myr_Mattea2022"])
    # allow q=1
    plot_all_models(max_q=1,include_blends=False,include_lit=True,
                    period_scale = "linear",models_to_plot=["WideHat8Myr_Mattea2015","WideHat8Myr_Mattea2022"])


    # # Demo PMD/model plots for the CS poster

    pmd = PeriodMassDistribution(max_q=0,include_blends=False,include_lit=True)
    pmd.calc_mass_percentiles(mass_bins=np.linspace(0.05,1.35,14))
    ax = pmd.plot_obs(plot_errors=True)
    _ = pmd.plot_period_perc(ax)
    _ = ax.set_title(pmd.param_string)
    ax.set_yscale("log")

    sm = SpinModel("WideHat8Myr_Mattea2015",130,"linear",
                   init_type="kde")
    pmd.select_obs(sm)

    # Normalize the model and calculate tau-squared
    # sm.normalize()
    sm.add_mask()
    sm.calc_tau_sq(pmd)

    # Plot
    ax = sm.plot_hist()
    pmd.plot_obs(ax)
    ax.set_ylim(0,12)

    ax.set_title(f"Matt+2015; 130 Myr; tau2={sm.tau_sq:.0f}")
    plt.savefig(f"plots/tausq_WideHat8Myr_Mattea2015_{pmd.param_string}_linear_00130Myr_ZAMS.png",bbox_inches="tight",dpi=600)

    like = np.exp(-0.5*sm.tau_sq)
    print(f"likelihood {like}")

    sm = SpinModel("WideHat8Myr_Mattea2022",80,"linear",
                   init_type="kde")
    pmd.select_obs(sm)

    # Normalize the model and calculate tau-squared
    # sm.normalize()
    sm.add_mask()
    sm.calc_tau_sq(pmd)

    # Plot
    ax = sm.plot_hist()
    pmd.plot_obs(ax)
    ax.set_ylim(0,12)
    ax.set_title(f"Matt+ in prep; 80 Myr; tau2={sm.tau_sq:.0f}")
    plt.savefig(f"plots/tausq_WideHat8Myr_Mattea2022_{pmd.param_string}_linear_00080Myr_ZAMS.png",bbox_inches="tight",dpi=600)

    like = np.exp(-0.5*sm.tau_sq)
    print(f"likelihood {like}")
