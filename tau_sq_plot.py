import os, sys, glob, time
import itertools
import multiprocessing as mp

import numpy as np
from numpy.random import default_rng
import matplotlib.pyplot as plt
import astropy.io.ascii as at
from astropy.io import fits
import astropy.units as u
from astropy import table
from astropy.table import join,vstack,Table
from astropy.coordinates import SkyCoord

import matplotlib as mpl
import matplotlib.cm as cm

norm = mpl.colors.Normalize(vmin=0, vmax=5)
mapper = cm.ScalarMappable(norm=norm, cmap=cm.viridis)

norm2 = mpl.colors.Normalize(vmin=0, vmax=14)
mapper2 = cm.ScalarMappable(norm=norm2, cmap=cm.viridis)

from analyze_cluster_output import colors, shapes
from plot_periods import plot_periodcolor_histogram
from tau_sq import PeriodMassDistribution, SpinModel


# In[2]:


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
# In[3]:

def plot_all_upsco():

    model_names = ["UpSco_Mattea2015","UpSco_Mattea2022","UpSco_ZeroTorque"]
    period_scale = "linear"

    pmd = PeriodMassDistribution()

    fig = plt.figure()
    fig.patch.set_facecolor('w')
    fig.patch.set_alpha(1.0)
    for j,model in enumerate(model_names):
        # print(model)
        models = glob.glob(os.path.expanduser(f"~/Dropbox/Models/{model}/{model}*Myr.txt"))
        # print(models)

        model_ages = np.sort([int(mod.split("_")[-1][:5]) for mod in models])
        # print(model_ages)

        model_ages = model_ages[model_ages<=300]

        all_tau_sq = np.zeros(len(model_ages))
        for i, age in enumerate(model_ages):
            sm = SpinModel(model,age,period_scale)

            # Normalize the model and calculate tau-squared
            sm.normalize()
            sm.add_mask()
            sm.calc_tau_sq(pmd)

            all_tau_sq[i] = sm.tau_sq

        plt.plot(model_ages,all_tau_sq,'--',label=model,color=mapper.to_rgba(j+1),alpha=0.75)
        plt.legend(loc=2)
        plt.xlabel("Model age (Myr)",fontsize=16)
        plt.ylabel("tau squared",fontsize=16)
    #     plt.title(model,fontsize=14)

        ax = plt.gca()
        ax.tick_params(labelsize=12)
        ax.set_xticks(np.arange(0,300,25),minor=True)
    plt.savefig(f"plots/tausq_ZAMS_UpSco.png",bbox_inches="tight",dpi=600)

def plot_all_widehat():

    model_names = ["WideHat8Myr_Mattea2015","WideHat8Myr_Mattea2022","WideHat8Myr_ZeroTorque"]
    period_scale = "linear"

    pmd = PeriodMassDistribution()

    fig = plt.figure()
    fig.patch.set_facecolor('w')
    fig.patch.set_alpha(1.0)
    for j,model in enumerate(model_names):
        print(model)
        models = glob.glob(os.path.expanduser(f"~/Dropbox/Models/{model}/{model}*Myr.txt"))
        # print(models)

        model_ages = np.sort([int(mod.split("_")[-1][:5]) for mod in models])
        # print(model_ages)

        model_ages = model_ages[model_ages<=300]

        all_tau_sq = np.zeros(len(model_ages))
        for i, age in enumerate(model_ages):
            sm = SpinModel(model,age,period_scale)

            # Normalize the model and calculate tau-squared
            sm.normalize()
            sm.add_mask()
            sm.calc_tau_sq(pmd)

            all_tau_sq[i] = sm.tau_sq

        plt.plot(model_ages,all_tau_sq,'-',label=model,color=mapper.to_rgba(j+1),alpha=0.75)
        plt.legend(loc=2)
        plt.xlabel("Model age (Myr)",fontsize=16)
        plt.ylabel("tau squared",fontsize=16)
    #     plt.title(model,fontsize=14)

        ax = plt.gca()
        ax.tick_params(labelsize=12)
        ax.set_xticks(np.arange(0,300,25),minor=True)
    plt.savefig(f"plots/tausq_ZAMS_WideHat8Myr.png",bbox_inches="tight",dpi=600)



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