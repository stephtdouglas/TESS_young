
import os, sys

import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.cm as cm
import numpy as np
import astropy.io.ascii as at


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

if __name__=="__main__":

    clusters = ["IC_2391","Collinder_135","NGC_2451A","NGC_2547","IC_2602"]

    plt.figure(figsize=(4,8))
    ax = plt.subplot(111)
    for i, cluster in enumerate(clusters):
        cf = cluster.replace("_"," ")
        filename = f"tables/TESS Cluster parameters - {cf}.csv"
        dat = at.read(filename,delimiter=",")
        dat = dat[(dat["Year"]>=2010) & (dat["Age (Myr)"].mask==False)]
        dat.sort("Year")
        dat["yvals"] = np.arange(len(dat))

        print(cluster,dat.dtype)

        # TODO: I should fix those inputs so they have an upper and lower
        # uncertainty for everything
        # TODO: also need to figure out something for those range or ~ ages
        # because right now that's making the columns into strings which
        # is causing some issues
        ax.plot(dat["Age (Myr)"],dat["yvals"],marker=shapes[cluster],
                color=colors[cluster],linewidth=0)


    ax.set_xlim(1,300)
    ax.set_xscale("log")
    ax.tick_params(labelleft=False)

    ax.axvline(25,color="k",linestyle="--",zorder=-10)
    ax.axvline(55,color="k",linestyle="--",zorder=-10)
    plt.show()
