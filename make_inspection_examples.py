import sys,os

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.gridspec as gridspec
from lightkurve import search_lightcurve
import astropy.io.ascii as at

from k2spin import prot

from make_inspection_plots_tess import plot_phased, read_tess

cmap = plt.cm.viridis_r

def plot_lcs(TIC, row, data_dir):
    """
    Plot light curves for a given target
    """

    t, f = read_tess(row, data_dir)

    ylims = np.percentile(f[np.isfinite(f)],[0.05,99.05])

    color1 = cmap(0.3)#plt.cm.inferno(0.5)
    color2 = cmap(0.9)
    label_fontsize = 12

    ### Periodogram
    ls_out = prot.run_ls(t,f,np.ones_like(f),0.1,prot_lims=[0.1,70],run_bootstrap=False)
    periods, pgram = ls_out[2], ls_out[3]


    plt.figure(figsize=(5,3))
    ax0 = plt.subplot(211)
    ### Full light curve
    ax0.plot(t,f,'k.',ms=3)
    ax0.set_ylim(ylims)
    ax0.set_xlim(t[0],t[-1])
    # ax0.set_ylabel("Light curve",fontsize=label_fontsize)
    ax0.set_yticklabels([])
    # ax0.tick_params(labelbottom=False)
    # ax0.set_xlabel("Time (d)",fontsize=label_fontsize)

    ax2 = plt.subplot(223)
    ### Phase folded 1
    if row["sig_periods"]>0:
        plot_phased(ax2,t,f,row["sig_periods"],row["sig_powers"],color=color1,ms=3)
        ax2.set_ylim(ylims)
        ax2.set_ylabel(f"P1={row['sig_periods']:.2f} d",fontsize=label_fontsize)
        ax2.get_legend().remove()
        ax2.tick_params(labelleft=False,labelright=False,labelbottom=False)

        if row["sig_periods"]>=2:
            repeats = np.arange(t[0],t[-1],row["sig_periods"])
            for r in repeats:
                ax0.axvline(r,ls="--",color=color1,lw=2)
        # axes[4].tick_params(labelbottom=False)
        ax2.set_xlabel("Phase",fontsize=label_fontsize)
    else:
        ax2.set_axis_off()


    ax3 = plt.subplot(224)
    ### Phase folded 2
    if row["sec_periods"]>0:
        plot_phased(ax3,t,f,row["sec_periods"],row["sec_powers"],color=color2,ms=3)
        ax3.set_ylim(ylims)
        ax3.set_ylabel(f"P2={row['sec_periods']:.2f} d",fontsize=label_fontsize)
        ax3.yaxis.set_ticks_position("right")
        ax3.yaxis.set_label_position("right")
        ax3.tick_params(labelleft=False,labelright=False,labelbottom=False)
        ax3.get_legend().remove()

        # print(row["sec_periods"])
        if row["sec_periods"]>=2:
            repeats = np.arange(t[0],t[-1],row["sec_periods"])
            # print(repeats)
            for r in repeats:
                # print(r)
                ax0.axvline(r,ls=":",color=color2,lw=2)
        ax3.set_xlabel("Phase",fontsize=label_fontsize)
    else:
        # axes[5].set_yticks([])
        # axes[5].set_yticklabels([])
        ax3.set_axis_off()

    plt.subplots_adjust(hspace=0,wspace=0)

    return [ax0,ax2,ax3]

if __name__=="__main__":

    ### Peaks don't match light curve
    cluster = "Collinder_135"
    date = "2021-06-18"
    TIC, lc_type, flux_col, sector = 23160350,"CDIPS","PCA1",7

    # Read in results for this cluster
    resfile = f"tables/{cluster}_{date}_results_raw.csv"
    res = at.read(resfile,delimiter=",")

    # Look for the lightcurve file; download it if needed
    search = search_lightcurve(f"TIC {TIC}",author=lc_type,sector=sector)
    lc_dir = os.path.expanduser("~/data/.lightkurve-cache/")
    lc = search.download_all(download_dir=lc_dir)

    # send it to the plotting function
    loc = np.where((res["target_name"]==TIC) & (res["provenance_name"]==lc_type)
                   & (res["sequence_number"]==sector) & (res["flux_cols"]==flux_col))[0]
    if len(loc)==1:
        axes = plot_lcs(TIC, res[loc][0], os.path.join(lc_dir,"mastDownload/HLSP"))
    else:
        print("uh oh",loc)

    plt.suptitle(f"TIC {TIC}; Peaks don't match light curve",y=0.93)
    plt.savefig("plots/inspect_example_1.png")

    ### Double dip, selected longer period
    cluster = "Collinder_135"
    date = "2021-06-18"
    TIC, lc_type, flux_col, sector = 97796539,"CDIPS","PCA1",7

    # Read in results for this cluster
    resfile = f"tables/{cluster}_{date}_results_raw.csv"
    res = at.read(resfile,delimiter=",")

    # Look for the lightcurve file; download it if needed
    search = search_lightcurve(f"TIC {TIC}",author=lc_type,sector=sector)
    lc_dir = os.path.expanduser("~/data/.lightkurve-cache/")
    lc = search.download_all(download_dir=lc_dir)

    # send it to the plotting function
    loc = np.where((res["target_name"]==TIC) & (res["provenance_name"]==lc_type)
                   & (res["sequence_number"]==sector) & (res["flux_cols"]==flux_col))[0]
    if len(loc)==1:
        axes = plot_lcs(TIC, res[loc][0], os.path.join(lc_dir,"mastDownload/HLSP"))
    else:
        print("uh oh",loc)

    plt.suptitle(f"TIC {TIC}; Double dip",y=0.93)
    plt.savefig("plots/inspect_example_2.png")


    ### Spurious high power
    cluster = "Collinder_135"
    date = "2021-06-18"
    TIC, lc_type, flux_col, sector = 284492865,"QLP","sap_flux",7

    # Read in results for this cluster
    resfile = f"tables/{cluster}_{date}_results_raw.csv"
    res = at.read(resfile,delimiter=",")

    # Look for the lightcurve file; download it if needed
    search = search_lightcurve(f"TIC {TIC}",author=lc_type,sector=sector)
    lc_dir = os.path.expanduser("~/data/.lightkurve-cache/")
    lc = search.download_all(download_dir=lc_dir)

    # send it to the plotting function
    loc = np.where((res["target_name"]==TIC) & (res["provenance_name"]==lc_type)
                   & (res["sequence_number"]==sector) & (res["flux_cols"]==flux_col))[0]
    if len(loc)==1:
        axes = plot_lcs(TIC, res[loc][0], os.path.join(lc_dir,"mastDownload/HLSP"))
    else:
        print("uh oh",loc)

    plt.suptitle(f"TIC {TIC}; Spurious, high power",y=0.93)
    plt.savefig("plots/inspect_example_3.png")


    ### Two periods
    cluster = "Collinder_135"
    date = "2021-06-18"
    TIC, lc_type, flux_col, sector = 173134217,"CDIPS","PCA1",7

    # Read in results for this cluster
    resfile = f"tables/{cluster}_{date}_results_raw.csv"
    res = at.read(resfile,delimiter=",")

    # Look for the lightcurve file; download it if needed
    search = search_lightcurve(f"TIC {TIC}",author=lc_type,sector=sector)
    lc_dir = os.path.expanduser("~/data/.lightkurve-cache/")
    lc = search.download_all(download_dir=lc_dir)

    # send it to the plotting function
    loc = np.where((res["target_name"]==TIC) & (res["provenance_name"]==lc_type)
                   & (res["sequence_number"]==sector) & (res["flux_cols"]==flux_col))[0]
    if len(loc)==1:
        axes = plot_lcs(TIC, res[loc][0], os.path.join(lc_dir,"mastDownload/HLSP"))
    else:
        print("uh oh",loc)

    plt.suptitle(f"TIC {TIC}; Two periods",y=0.93)
    plt.savefig("plots/inspect_example_4.png")
