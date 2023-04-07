import os, sys, glob, pathlib

import numpy as np
import astropy.io.ascii as at
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.cm as cm
from matplotlib.patches import Rectangle

import tess_young
from tess_young.get_const import *
_DIR = pathlib.Path(tess_young.__file__).resolve().parent.parent
plt.style.use(os.path.join(_DIR,'paper.mplstyle'))

from periodmass import PeriodMassDistribution
from spinmodel import SpinModel


# mapper and normalization for binning by mass
norm2 = mpl.colors.Normalize(vmin=0, vmax=14)
mapper2 = cm.ScalarMappable(norm=norm2, cmap=cm.viridis)

def plot_tracks_binned(ttab_filename, model_name, init_type, axes = None,
                       mass_limits=[0.05,1.4]):

    outfilename = ttab_filename.split("/")[-1].split(".")[0]
    ttab_filepath = os.path.join(_DIR,ttab_filename)
    model = model_name
    

    if os.path.exists(ttab_filepath) is False:
        print("ERROR: run binned tau-sq fit first")
        print(ttab_filepath)
        return

    ttab = at.read(ttab_filepath)
    # # The first column is the age, so remove that
    # mass_bins = ttab.dtype.names[1:]
    model_ages = ttab["Age(Myr)"]


    
    mass_bins = np.arange(mass_limits[0],mass_limits[1],0.1)
    mass_labels = [f"{mass:.2f} Msun" for mass in mass_bins]
    nbins = len(mass_labels)-1


    if axes is None:
        fig, axes = plt.subplots(nrows=nbins,figsize=(8,10), sharex=True)
        fig.patch.set_facecolor('w')
        fig.patch.set_alpha(1.0)

        plt.subplots_adjust(hspace=0)

    if len(axes)<nbins:
        print("ERROR: not enough axes for selected mass bins")
        return

    if "Zero" in model:
        ls = ":"
    elif "2015" in model:
        ls = "-."
    else:
        ls = "-"

    for i,mass in enumerate(mass_bins[:-1]):
        mcolor = mapper2.to_rgba((i % nbins)+1)
        axes[i].plot(model_ages,ttab[mass_labels[i]],ls,label=display_names[model_name],
                color=mcolor,alpha=0.75)
        if i==nbins:
            axes[i].legend()

        # min_i = np.argmin(ttab[mass][i])
        # ax.plot(model_ages[min_i],0.99,"^",color=mapper2.to_rgba((i % nbins)+1))

        if axes[i].get_ylim()[1]>30:
            max_ts = np.max(ttab[mass_labels[i]][model_ages<150])
            new_ymax = max_ts+10
            axes[i].set_ylim(top=new_ymax)

        textx = 0.01 if i>0 else 0.4
        axes[i].text(textx,0.9,mass_labels[i],horizontalalignment="left",
                     verticalalignment="top",transform=axes[i].transAxes,
                     color=mcolor)

        if i==(nbins // 2):
            axes[i].set_ylabel("tau squared")

    axes[0].legend(prop={'size': 10})
    axes[-1].set_xlabel("Model age (Myr)")
    axes[-1].set_xlim(0,150)

    return axes

if __name__=="__main__":

    # These give the basic information about the selected run, regardless of model
    fname1 = "tables/tausq_linear_massall_"
    fname3 = "_binned_Qmax0_blendsFalse_litTrue.csv"
    plot_string = "Qmax0_blendsFalse_litTrue"


    # Plot the higher mass bins for the tophat models
    models = model_names[3:]

    axes = None
    for j,model in enumerate(models):
        fname = fname1 + model + "_tophat" + fname3
        # print(fname)
        axes = plot_tracks_binned(fname,model,"tophat",axes,
                                  mass_limits=[0.65,1.3])

    axes[0].set_title(plot_string+"_WideHat")
    plt.savefig(os.path.join(_DIR,f"plots/test_tausq_tracks_masshigh_widehat.png"),bbox_inches="tight",dpi=600)

    # Plot the lower mass bins for the tophat models

    axes = None
    for j,model in enumerate(models):
        fname = fname1 + model + "_tophat"  + fname3
        axes = plot_tracks_binned(fname,model,"tophat",axes,
                                  mass_limits=[0.25,0.8])

    axes[0].set_title(plot_string+"_WideHat")
    plt.savefig(os.path.join(_DIR,f"plots/test_tausq_tracks_masslow_widehat.png"),bbox_inches="tight",dpi=600)


    # Plot the higher mass bins for the USco KDE models

    axes = None
    for j,model in enumerate(models):
        fname = fname1 + model + "_kde"  + fname3
        axes = plot_tracks_binned(fname,model,"kde",axes,
                                  mass_limits=[0.65,1.3])

    axes[0].set_title(plot_string+"_USco")
    plt.savefig(os.path.join(_DIR,f"plots/test_tausq_tracks_masshigh_usco.png"),bbox_inches="tight",dpi=600)


    # Plot the lower mass bins for the USco KDE models
    axes = None
    for j,model in enumerate(models):
        fname = fname1 + model + "_kde" + fname3
        axes = plot_tracks_binned(fname,model,"kde",axes,
                                  mass_limits=[0.25,0.8])


    axes[0].set_title(plot_string+"_USco")
    plt.savefig(os.path.join(_DIR,f"plots/test_tausq_tracks_masslow_usco.png"),bbox_inches="tight",dpi=600)

