import os, sys
import glob
import itertools

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.cm as cm
import astropy.io.ascii as at
from astropy.table import join,vstack,Table
from scipy import stats
from scipy.interpolate import interp1d

from analyze_cluster_output import read_cluster_visual

norm = mpl.colors.LogNorm(vmin=0.1, vmax=30)
mapper = cm.ScalarMappable(norm=norm, cmap=cm.viridis)

# Read in and merge the outputs from k2spin

# colors = {"IC_2391": "C0",
#          "IC_2602": "C4",
#          "NGC_2547": "C3",
#          "NGC_2451A": "C2",
#          "Collinder_135": "C1"}

cmap2 = cm.get_cmap("viridis",7)
colors = {"IC_2391": cmap2(0),
         "IC_2602": cmap2(4),
         "NGC_2547": cmap2(3),
         "NGC_2451A": cmap2(2),
         "Collinder_135": cmap2(1)}


shapes= {"IC_2391": "o",
         "IC_2602": "d",
         "NGC_2547": "v",
         "NGC_2451A": "^",
         "Collinder_135": "s"}


def plot_all(clean_limit=10,to_plot_indiv=True):
    bp_rp_IC_2391, prot_IC_2391 = read_cluster_visual("IC_2391","2021-06-22",clean_limit,to_plot=to_plot_indiv)
    bp_rp_Collinder_135, prot_Collinder_135 = read_cluster_visual("Collinder_135","2021-06-18",clean_limit,to_plot=to_plot_indiv)
    bp_rp_NGC_2451A, prot_NGC_2451A = read_cluster_visual("NGC_2451A","2021-06-21",clean_limit,to_plot=to_plot_indiv)
    bp_rp_NGC_2547, prot_NGC_2547 = read_cluster_visual("NGC_2547","2021-06-21",clean_limit,to_plot=to_plot_indiv)
    bp_rp_IC_2602, prot_IC_2602 = read_cluster_visual("IC_2602","2021-07-02",clean_limit,to_plot=to_plot_indiv)


    plt.figure()
    plt.plot(bp_rp_IC_2391, prot_IC_2391,"o",label="IC_2391",ms=5,color=colors["IC_2391"])
    plt.plot(bp_rp_Collinder_135, prot_Collinder_135,"s",label="Collinder_135",ms=5,color=colors["Collinder_135"])
    plt.plot(bp_rp_NGC_2451A, prot_NGC_2451A,"^",label="NGC_2451A",ms=5,color=colors["NGC_2451A"])
    plt.plot(bp_rp_NGC_2547, prot_NGC_2547,"v",label="NGC_2547",ms=5,color=colors["NGC_2547"])
    plt.plot(bp_rp_IC_2602, prot_IC_2602,"d",label="IC_2602",ms=5,color=colors["IC_2602"])
    plt.legend(loc=1)

    plt.ylim(0.1,50)
    plt.xlim(0.5,3.5)
    plt.yscale("log")

    plt.xlabel(r"G$_{BP}$ - G$_{RP}$")
    plt.ylabel("Period (d)")

    ax = plt.gca()
    ax.axhline(12,linestyle="--",color="tab:grey")

    plt.savefig(f"plots/periodmass_all_clean{clean_limit}.png")

def plot_periodcolor(clean_limit=10,to_plot_indiv=False):
    bp_rp_IC_2391, prot_IC_2391 = read_cluster_visual("IC_2391","2021-06-22",clean_limit,to_plot=to_plot_indiv)
    bp_rp_Collinder_135, prot_Collinder_135 = read_cluster_visual("Collinder_135","2021-06-18",clean_limit,to_plot=to_plot_indiv)
    bp_rp_NGC_2451A, prot_NGC_2451A = read_cluster_visual("NGC_2451A","2021-06-21",clean_limit,to_plot=to_plot_indiv)
    bp_rp_NGC_2547, prot_NGC_2547 = read_cluster_visual("NGC_2547","2021-06-21",clean_limit,to_plot=to_plot_indiv)
    bp_rp_IC_2602, prot_IC_2602 = read_cluster_visual("IC_2602","2021-07-02",clean_limit,to_plot=to_plot_indiv)


    plt.figure()
    plt.plot(bp_rp_IC_2391, prot_IC_2391,"o",label="IC_2391",ms=5,color="k")#colors["IC_2391"])
    plt.plot(bp_rp_Collinder_135, prot_Collinder_135,"s",label="Collinder_135",ms=5,color="k")#colors["Collinder_135"])
    plt.plot(bp_rp_NGC_2451A, prot_NGC_2451A,"^",label="NGC_2451A",ms=5,color="k")#colors["NGC_2451A"])
    plt.plot(bp_rp_NGC_2547, prot_NGC_2547,"v",label="NGC_2547",ms=5,color="k")#colors["NGC_2547"])
    plt.plot(bp_rp_IC_2602, prot_IC_2602,"d",label="IC_2602",ms=5,color="k")#colors["IC_2602"])
    plt.legend(loc=1)

    plt.ylim(0.1,50)
    plt.xlim(0.5,3.5)
    plt.yscale("log")

    plt.xlabel(r"G$_{BP}$ - G$_{RP}$")
    plt.ylabel("Period (d)")

    ax = plt.gca()
    ax.axhline(12,linestyle="--",color="tab:grey")

    return ax


def plot_periodcolor_ruwe(clean_limit=0,to_plot_indiv=False):
    cat_IC_2391 = read_cluster_visual("IC_2391","2021-06-22",clean_limit,to_plot=to_plot_indiv,return_periodcolor=False)
    cat_Collinder_135 = read_cluster_visual("Collinder_135","2021-06-18",clean_limit,to_plot=to_plot_indiv,return_periodcolor=False)
    cat_NGC_2451A = read_cluster_visual("NGC_2451A","2021-06-21",clean_limit,to_plot=to_plot_indiv,return_periodcolor=False)
    cat_NGC_2547 = read_cluster_visual("NGC_2547","2021-06-21",clean_limit,to_plot=to_plot_indiv,return_periodcolor=False)
    cat_IC_2602 = read_cluster_visual("IC_2602","2021-07-02",clean_limit,to_plot=to_plot_indiv,return_periodcolor=False)
    clusters = ["IC_2391","Collinder_135","NGC_2451A","NGC_2547","IC_2602"]
    cats = [cat_IC_2391, cat_Collinder_135, cat_NGC_2451A, cat_NGC_2547, cat_IC_2602]

    fig1 = plt.figure()
    ax = plt.subplot(111)


    for i, match in enumerate(cats):
        # print(match.dtype)
        bp_rp = match["GAIAEDR3_BP"] - match["GAIAEDR3_RP"]

        fig2 = plt.figure()
        ax2 = plt.subplot(111)


        if clean_limit is not None:
            if clean_limit==0:
                all_possible = (match["final_Q"]==0)
            else:
                all_possible = (match["final_Q"]==0) | (match["final_Q"]==1)
        else:
            all_possible = (match["final_Q"]==0) | (match["final_Q"]==1)

        sc = ax.scatter(bp_rp[all_possible],match["final_period"][all_possible],
                    c=match["GAIAEDR3_RUWE"][all_possible],vmin=0.6,vmax=2,alpha=0.75)
        sc = ax2.scatter(bp_rp[all_possible],match["final_period"][all_possible],
                    c=match["GAIAEDR3_RUWE"][all_possible],vmin=0.6,vmax=2,alpha=0.75)
        plt.colorbar(sc,label="Gaia RUWE - EDR3",ax=ax2)

        ax2.set_ylim(0.1,50)
        ax2.set_xlim(0.5,3.5)
        ax2.set_yscale("log")

        ax2.set_xlabel(r"G$_{BP}$ - G$_{RP}$")
        ax2.set_ylabel("Period (d)")
        ax2.axhline(12,linestyle="--",color="tab:grey")

        ax2.set_title(clusters[i])
        fig2.savefig(f"plots/periodcolor_{clusters[i]}_ruwe.png",bbox_inches="tight",dpi=300)
        plt.close(fig2)


    # TODO: add gemini observations and detections

    plt.colorbar(sc,label="Gaia RUWE - EDR3",ax=ax)

    ax.set_ylim(0.1,50)
    ax.set_xlim(0.5,3.5)
    ax.set_yscale("log")

    ax.set_xlabel(r"G$_{BP}$ - G$_{RP}$")
    ax.set_ylabel("Period (d)")

    ax.axhline(12,linestyle="--",color="tab:grey")

    return ax

def plot_periodcolor_histogram(clean_limit=10,to_plot_indiv=False,ax=None):
    bp_rp_IC_2391, prot_IC_2391 = read_cluster_visual("IC_2391","2021-06-22",clean_limit,to_plot=to_plot_indiv)
    bp_rp_Collinder_135, prot_Collinder_135 = read_cluster_visual("Collinder_135","2021-06-18",clean_limit,to_plot=to_plot_indiv)
    bp_rp_NGC_2451A, prot_NGC_2451A = read_cluster_visual("NGC_2451A","2021-06-21",clean_limit,to_plot=to_plot_indiv)
    bp_rp_NGC_2547, prot_NGC_2547 = read_cluster_visual("NGC_2547","2021-06-21",clean_limit,to_plot=to_plot_indiv)
    bp_rp_IC_2602, prot_IC_2602 = read_cluster_visual("IC_2602","2021-07-02",clean_limit,to_plot=to_plot_indiv)

    bp_rp = np.concatenate([bp_rp_IC_2391,bp_rp_Collinder_135,bp_rp_NGC_2451A,bp_rp_NGC_2547,bp_rp_IC_2602])
    prot = np.concatenate([prot_IC_2391,prot_Collinder_135,prot_NGC_2451A,prot_NGC_2547,prot_IC_2602])

    if ax is None:
        plt.figure()
        ax = plt.subplot(111)
    counts, xedges, yedges, image = ax.hist2d(bp_rp,prot,bins=[np.linspace(0.5,3.5,22),
                                          np.logspace(-1,np.log10(50),20)],cmap="Greys")
    xx = xedges[:-1] + (xedges[1:] - xedges[:-1])/2
    yy = yedges[:-1] + (yedges[1:] - yedges[:-1])/2
    cs = ax.contour(xx,yy,counts.transpose(),linewidths=2,levels=[1,3,10,13,20],colors="k",antialiased=True)
    # plt.clabel(cs,fmt="%1.0f",inline=False)

    ax.set_ylim(0.1,50)
    ax.set_xlim(0.5,3.5)
    ax.set_yscale("log")

    ax.set_xlabel(r"G$_{BP}$ - G$_{RP}$")
    ax.set_ylabel("Period (d)")

    ax.axhline(12,linestyle="--",color="tab:grey")

    return ax

if __name__=="__main__":
    plot_periodcolor_ruwe()
    plt.savefig("plots/periodcolor_all_ruwe.png",bbox_inches="tight",dpi=300)
    # plt.show()
