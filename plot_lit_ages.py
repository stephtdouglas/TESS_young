
import os, sys

import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.cm as cm
import numpy as np
import astropy.io.ascii as at
from astropy.table import join, Table, vstack


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

    all_dat = []

    # First, stack all the values together
    for i, cluster in enumerate(clusters):
        cf = cluster.replace("_"," ")
        filename = f"tables/TESS Cluster parameters - {cf}.csv"
        dat = at.read(filename,delimiter=",")
        dat = dat[(dat["Year"]>=2010) & (dat["Age (Myr)"].mask==False)]
        dat["Cluster"] = np.empty(len(dat),"U14")
        dat["Cluster"][:] = cluster
        # These columns are strings in only some tables, just remove them
        dat.remove_column("e_Age (Myr)")
        dat.remove_column("e_Distance (pc)")
        all_dat.append(dat)

    dat = vstack(all_dat)
    dat.sort(["Year","Source"])
    print(dat["Source"])

    dat_grouped = dat.group_by(["Year","Source"])
    print(dat_grouped.groups.keys)

    # for i in range(len(dat_grouped.groups.keys)):

    yval = 0
    plt.figure(figsize=(4,8))
    ax = plt.subplot(111)   # TODO: plot each cluster individually, and use text to label the source
    for i, keys in enumerate(dat_grouped.groups.keys):
        sub = dat_grouped.groups[i]
        source = keys["Source"]


        sub["yvals"] = np.ones(len(sub))*yval

        # print(cluster,dat.dtype)


        asym_err = np.asarray([sub["e_dn_Age"],sub["e_up_Age"]])
        # print(asym_err)

        # ax.plot(dat["Age (Myr)"],dat["yvals"],marker=shapes[cluster],
        #         color=colors[cluster],linewidth=0)
        ax.errorbar(sub["Age (Myr)"],sub["yvals"],
                    xerr=asym_err,marker=shapes[cluster],
                    color=colors[cluster],linewidth=0,elinewidth=1.5)
    #
    #     # break
    #
    ax.set_xlim(1,300)
    ax.set_xscale("log")
    ax.tick_params(labelleft=False)
    
    ax.axvline(25,color="k",linestyle="--",zorder=-10)
    ax.axvline(55,color="k",linestyle="--",zorder=-10)
    plt.show()
