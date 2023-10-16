"""
Script to compare literature periods to TESS periods
"""
import os, sys, pathlib
from datetime import date

import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.cm as cm
import numpy as np
import astropy.io.ascii as at

import tess_young
from tess_young.get_const import *
from tess_young import plot_spt_axis
_DIR = pathlib.Path(tess_young.__file__).resolve().parent.parent
plt.style.use(os.path.join(_DIR,'paper.mplstyle'))

lit_clusters = ["IC_2391","IC_2602","NGC_2547"]

lit_papers = {"IC_2391":"Patten & Simon (1996),\nMessina+ (2011)",
              "IC_2602":"Patten+ (1996), Barnes+ (1999)",
              "NGC_2547":"Irwin+ (2008)"
              }

def compare_literature():

    tab = at.read("tab_all_stars.csv")
    ltab = tab[tab["LitPeriod"].mask==False]

    fig, axes = plt.subplots(1,3,figsize=(16,5.5),sharex=True,sharey=True)


    for i, cluster in enumerate(lit_clusters):

        ax = axes[i]

        q0 = (ltab["Q1"]==0) & (ltab["Cluster"]==cluster)
        q1 = (ltab["Q1"]==1) & (ltab["Cluster"]==cluster)
        ax.plot(ltab["LitPeriod"][q0],ltab["Prot1"][q0],shapes[cluster],
                 color=colors[cluster],label=lit_papers[cluster])
        ax.plot(ltab["LitPeriod"][q1],ltab["Prot1"][q1],shapes[cluster],
                color=colors[cluster],mfc="none",label="TESS low-quality")
        ax.legend(loc=2)

        multi = np.where((q0 | q1) & (ltab["Q2"]<=1))[0]
        for i in multi:
            ax.plot([ltab["LitPeriod"][i],ltab["LitPeriod"][i]],
                    [ltab["Prot1"][i],ltab["Prot2"][i]],"-",color=colors[cluster],
                    zorder=-5)
            if ltab["Q2"][i]==0:
                ax.plot(ltab["LitPeriod"][i],ltab["Prot2"][i],shapes[cluster],
                        color=colors[cluster],label=lit_papers[cluster])
            else:
                ax.plot(ltab["LitPeriod"][i],ltab["Prot2"][i],shapes[cluster],
                        color=colors[cluster],mfc="none",label="TESS low-quality")



        # Identify significant discrepancies
        dbl = ltab["Prot1"]*2.1
        half = ltab["Prot1"]/2.1
        sig_diff_i = ((ltab["Q1"]<=1) & (ltab["Cluster"]==cluster) &
                      ((ltab["LitPeriod"]<half) | (ltab["LitPeriod"]>dbl)))
        print(ltab[sig_diff_i])
        at.write(ltab[sig_diff_i],f"{cluster}_significant_differences.csv",delimiter=",")

        x = np.linspace(0.1,50,20)
        ax.plot(x,x,"-",zorder=-5,color="grey")
        ax.plot(x,x/2,"--",zorder=-5,color="grey")
        ax.plot(x,x*2,"--",zorder=-5,color="grey")

        if cluster=="NGC_2547":
            x = np.logspace(-1,2,200)
            for p in [1,2,3,4,5,6,7]:
                beat_period = 1/(1/x+1/p)
                ax.plot(x,beat_period,":",zorder=-5,color="lightgrey")
                ax.plot(beat_period,x,":",zorder=-5,color="lightgrey")

        if cluster=="IC_2391":
            # Connect the components of Cl* IC 2391 SHJM 4/5
            ax.plot([0.21,0.93],[0.11,0.11],'-',color=colors[cluster],zorder=-5)

        ax.set_xlim(0.1,50)
        ax.set_ylim(0.1,50)
        ax.set_xscale("log")
        ax.set_yscale("log")

        ax.tick_params(labelleft=True,labelbottom=True)
        ax.set_xlabel("Literature Period (d)")
        ax.set_title(cluster.replace("_"," "))
    axes[0].set_ylabel("TESS Period (d)")

if __name__=="__main__":

    compare_literature()

    plt.savefig(f"plots/literature_comparison.png",
                bbox_inches="tight")
    plt.savefig(os.path.join(PAPER_DIR,"fig_literature_comparison.pdf"),
                bbox_inches="tight")
    plt.close("all")
