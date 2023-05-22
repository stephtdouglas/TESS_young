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

from tess_young.tau_sq.periodmass import PeriodMassDistribution
from tess_young.tau_sq.spinmodel import SpinModel
from tess_young.tau_sq.upsco_periodmass import UpSco

if __name__=="__main__":

    # Plot the upper sco stars on the KDE itself


    mass_bins = np.arange(0.05,1.4,0.1)
    nbins = 30
    # period_bins = np.linspace(0,40,nbins)
    # period_scale="linear"
    period_bins = np.logspace(np.log10(0.08),np.log10(40),nbins)
    period_scale="log"

    kde_pfile = os.path.expanduser("~/Dropbox/Models/UpSco_KDE_init/all_KDE_corrected.csv")
    kde_prob = at.read(kde_pfile)

    mod = kde_prob


    img_raw, xedges, yedges = np.histogram2d(mod["Mass"],mod["prot"],
                                               weights=mod["Prob"],density=True,
                                               bins=[mass_bins,period_bins])


    X, Y = np.meshgrid(xedges, yedges)

    # Determine the cell area
    xdiff = np.diff(xedges) # mass
    ydiff = np.diff(yedges) #period

    Xdiff, Ydiff = np.meshgrid(xdiff,ydiff)
    cell_area = Xdiff*Ydiff

    # Now calculate the number density in each cell, for better visualization
    img = img_raw.T * cell_area

    # Mask out any regions where the KDE does not exist

    img_nomask = np.copy(img)

    # mask the image so it doesn't show cells outside the kde
    mass_array = mod["Mass"]
    prot_array = mod["prot"]
    model_exists = np.ones(np.shape(img_nomask),bool)
    mask = np.zeros(np.shape(img_nomask),bool)
    for i in range(len(mass_bins)-1):
        mass_loc = ((mass_array>=mass_bins[i]) &
                    (mass_array<mass_bins[i+1]))
        # Calculate whether there are periods in each individual bins
        for j in range(len(period_bins)-1):
            per_loc = ((prot_array>=period_bins[j]) &
                       (prot_array<period_bins[j+1]))
            in_this_bin = np.where(mass_loc & per_loc)[0]
    #             print(mass_bins[i],period_bins[j],in_this_bin)
            if len(in_this_bin)==0:
                model_exists[j,i] = False
        # Now, for this mass range, define the mask to only exclude bins
        # beyond the range of the model
    #         print(model_exists[:,i])
        mod_min_j = min(np.where(model_exists[:,i]==True)[0])
        mod_max_j = max(np.where(model_exists[:,i]==True)[0])
        mask[:mod_min_j,i] = True
        mask[mod_max_j+1:,i] = True

    img = np.ma.masked_array(img_nomask,mask=mask)


    fig, ax = plt.subplots(1,1)
    norm = LogNorm(vmin=2e-4,vmax=1e-2)
    pcm = ax.pcolormesh(X, Y, img,cmap="Greens",norm=norm)
    print(np.nanmin(img.flatten()),np.nanmax(img.flatten()))
    plt.xlim(1.25,0.05)
    if period_scale=="linear":
        plt.ylim(0,14)
    else:
        plt.ylim(0.1,30)
        plt.yscale(period_scale)


    usco = UpSco(None)
    ax.plot(usco.mass,usco.prot,'w.',mec="k")
    ax.set_xlabel("Mass")
    ax.set_ylabel("Period")
    plt.colorbar(pcm,label=r"$\rho_{number}$")

    plt.savefig(os.path.join(PAPER_DIR,"fig_upsco_kde.pdf"),
                bbox_inches="tight")
        
