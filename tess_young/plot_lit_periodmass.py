"""
Script to plot literature periods for motivation
"""

import os, sys, pathlib

import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.cm as cm
import numpy as np
import astropy.io.ascii as at
from astropy.table import join, Table, vstack
from astropy import units as u


from tess_young.get_const import *
from tess_young import plot_spt_axis
import tess_young
_DIR = pathlib.Path(tess_young.__file__).resolve().parent.parent
plt.style.use(os.path.join(_DIR,'paper.mplstyle'))
plt.style.use('./paper.mplstyle')

lit_clusters = ["IC_2391","NGC_2547","IC_2602"]

PAPER_DIR = os.path.expanduser("~/my_papers/TESS_young/")


def plot_literature_periodcolor(ax=None):

    if ax is None:
        plt.figure()
        ax = plt.subplot(111)
    # lit_color = "#595959"

    tab = at.read("tab_all_stars.csv")
    ltab = tab[tab["LitPeriod"].mask==False]
    bp_rp = ltab["GAIAEDR3_BP"] - ltab["GAIAEDR3_RP"]

    for i,cluster in enumerate(lit_clusters):
        loc = (ltab["Cluster"]==cluster) & (ltab["to_plot"]==1)

        ax.plot(bp_rp[loc],ltab["LitPeriod"][loc],shapes[cluster],
                color=colors[cluster],label=cluster.replace("_"," "))

    ax.legend(loc=1)

    for i,cluster in enumerate(lit_clusters):
        loc = (ltab["Cluster"]==cluster) & (ltab["to_plot"]==0)

        ax.plot(bp_rp[loc],ltab["LitPeriod"][loc],shapes[cluster],
                mec=colors[cluster],mfc="none",zorder=-5)

    ax.set_ylim(0.1,50)
    ax.set_xlim(0.5,3.5)
    ax.set_yscale("log")

    ax.set_xlabel(r"G$_{BP}$ - G$_{RP}$ (EDR3)")
    ax.set_ylabel("Period (d)")

    plot_spt_axis.add_spt_gaia(ax)

    return ax

def add_literature_periodcolor(ax=None):
    """
    Add the literature Prot to a periodcolor plot in grey
    """
    if ax is None:
        plt.figure()
        ax = plt.subplot(111)
    lit_color = "#595959"

    tab = at.read("tab_all_stars.csv")
    ltab = tab[tab["LitPeriod"].mask==False]
    bp_rp = ltab["GAIAEDR3_BP"] - ltab["GAIAEDR3_RP"]

    for i,cluster in enumerate(lit_clusters):
        loc = (ltab["Cluster"]==cluster) & (ltab["to_plot"]==1)

        ax.plot(bp_rp[loc],ltab["LitPeriod"][loc],shapes[cluster],
                color=lit_colors,label=cluster.replace("_"," ")+" Lit",
                zorder=-10)


if __name__=="__main__":
    _ = plot_literature_periodcolor()
    plt.savefig(f"plots/periodcolor_literature.png",bbox_inches="tight")
    plt.savefig(os.path.join(PAPER_DIR,"fig_periodcolor_literature.pdf"),
                bbox_inches="tight")
    plt.close("all")
