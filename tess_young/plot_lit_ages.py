
import os, sys, pathlib

import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.cm as cm
import numpy as np
import astropy.io.ascii as at
from astropy.table import join, Table, vstack


from tess_young.get_const import *
import tess_young
_DIR = pathlib.Path(tess_young.__file__).resolve().parent.parent
plt.style.use(os.path.join(_DIR,'paper.mplstyle'))

clusters = ["IC_2391","Collinder_135","NGC_2451A","NGC_2547","IC_2602"]

if __name__=="__main__":

    filename = os.path.join(_DIR,"tables/TESS_Cluster_parameters_Comparison_Ages.csv")
    dat = at.read(filename,delimiter=",")
    print(dat.dtype)

    dat.remove_row(np.where(dat["Source"]=="Brutus")[0][0])

    # Add year, then sort by year and author name
    dat["year"] = np.ones(len(dat),int)*9999
    for j, source in enumerate(dat["Source"]):
        if "Brutus" in source:
            continue
        else:
            dat["year"][j] = int(source[-4:])
    dat.sort(["year","Source"])

    # Plot each source of ages on a separate line
    dat["yvals"] = np.arange(len(dat))

    # Make figure
    fig = plt.figure(figsize=(5,8))
    ax = plt.subplot(111)
    
    # Plot the ages for each cluster on the appropriate level
    for i,cluster in enumerate(clusters):
        good = dat[f"age_{cluster}"].mask==False
        age = dat[f"age_{cluster}"][good]
        asym_err = np.asarray([dat[f"e_dn_age_{cluster}"][good],
                               dat[f"e_up_age_{cluster}"][good]])
        yvals = dat["yvals"][good]+(i/25)

        ax.errorbar(age,yvals,xerr=asym_err,marker=shapes[cluster],
                    color=colors[cluster],linewidth=0,elinewidth=1.5,
                    label=cluster.replace("_"," "))

    ax.legend(loc=2,ncol=2)

    # Add the source names (really the latex citekeys)
    # And plot two vertical lines encapsulating most of the ages
    textx = 1.1
    right_line = 55
    left_line = 25
    for j, source in enumerate(dat["Source"]):
        if "cummings" in source:
            ax.text(right_line*1.1,j,dat["Source_formatted"][j].replace("; ",";\n"),
                    fontsize=8,ha="left",verticalalignment="center")
        else:
            ax.text(textx,j,dat["Source_formatted"][j],fontsize=10)
    ax.axvline(left_line,color="Gray",linestyle="--",zorder=-10)
    ax.axvline(right_line,color="Gray",linestyle="--",zorder=-10)

    # Final plot adjustments
    ax.set_xlim(1,300)
    ax.set_xscale("log")

    ax.set_ylim(-0.5,max(dat["yvals"])+5)
    ax.tick_params(labelleft=False)
    ax.set_yticks([])
    
    ax.set_xlabel("Age [Myr]")
    plt.savefig(os.path.join(_DIR,"plots/literature_ages.png"),
                bbox_inches="tight")
    plt.savefig(os.path.join(PAPER_DIR,"fig_literature_ages.pdf"),
                bbox_inches="tight")
    plt.close("all")
