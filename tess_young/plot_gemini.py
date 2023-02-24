import os, sys

import matplotlib.pyplot as plt
import numpy as np
import astropy.io.ascii as at
import matplotlib as mpl
import matplotlib.cm as cm
# mpl.rcParams['image.cmap'] = 'viridis_r'

cmap2 = cm.get_cmap("viridis",14)
# colors = {"IC_2391": cmap2(0),
#          "IC_2602": cmap2(4),
#          "NGC_2547": cmap2(3),
#          "NGC_2451A": cmap2(2),
#          "Collinder_135": cmap2(1)}


# C466 = cmap2(0.5)
# C562 = cmap2(1.5)
# C716 = cmap2(2.5)
# C832 = cmap2(3.5)
C466 = cmap2(9)
C562 = cmap2(7)
C716 = cmap2(3)
C832 = cmap2(1)

def plot_contrast():

    clist = at.read(os.path.expanduser("~/Dropbox/data/Speckle/IC2391_contrast_files"))
    spath = os.path.expanduser("~/Dropbox/data/Speckle/")

    plt.figure(figsize=(8,6))

    plt.plot([],[],color=C466,alpha=0.75,label="Limit, 466nm filter")
    plt.plot([],[],color=C562,alpha=0.75,label="Limit, 562nm filter")
    plt.plot([],[],color=C716,alpha=0.75,label="Limit, 716nm filter")
    plt.plot([],[],color=C832,alpha=0.75,label="Limit, 832nm filter")
    plt.plot([],[],"^",color=C716,label="Detection, 716nm filter")
    plt.plot([],[],"v",color=C832,label="Detection, 832nm filter")
    plt.legend(loc=3)

    for filename in clist["filename"]:
        # print(filename)
        # Set up the appropriate arrays
        separation = []
        mag_limit = []

        # read in the contrast curves
        f = open(os.path.join(spath,filename),"r")
        l = f.readline()
        while l.startswith("# fit")==False:
            l = f.readline()
        l = f.readline()
        while l!="":
            lsplit = l.split()
            separation.append(lsplit[0])
            mag_limit.append(lsplit[1])
            l = f.readline()
        separation = np.asarray(separation,"float32")*1000
        mag_limit = np.asarray(mag_limit,"float32")
        if "_466" in filename:
            plt.plot(separation,mag_limit,color=C466,alpha=0.2)
        elif "_562" in filename:
            plt.plot(separation,mag_limit,color=C562,alpha=0.2)
        elif "_716" in filename:
            plt.plot(separation,mag_limit,color=C716,alpha=0.2)
        elif "_832" in filename:
            plt.plot(separation,mag_limit,color=C832,alpha=0.2)
        else:
            print(filename,"band not found")
            plt.plot(separation,mag_limit,color="grey",alpha=0.2)
        f.close()

    # add the detections themselves
    gem_dir = os.path.expanduser("~/proposals/Gemini/2021A/")
    det19 = at.read(os.path.join(gem_dir,"Douglas_Binaries_2019AB.txt"),delimiter="\s")
    for row in det19:
        # if (("IC2391" in row["ID"])==False) or (row["rho"]*1000==0.922):
        #     # Skip the IC2602 star, and the duplicate detection
        #     continue
        # else:
        if row["wav"]==832:
            plt.plot(row["rho"]*1000,row["delm"],"v",color=C832)
        elif row["wav"]==716:
            plt.plot(row["rho"]*1000,row["delm"],"^",color=C716)
        else:
            plt.plot(row["rho"],row["delm"],"o",color="grey")


    plt.xscale("log")
    plt.ylim(10,-1)
    plt.xlim(20,4000)
    ax = plt.gca()
    ax.tick_params(labelsize=13)
    ax.set_xlabel("Projected Separation (mas)",fontsize=14)
    ax.set_ylabel(r"Contrast (mag)",fontsize=14)
    plt.plot(np.array([2000,2000,700,550]),[12,5,2,0],linestyle="--",color="k")
    plt.text(800,2,"Gaia Catalog\nResolved (Approx.)",color="k")

    # plt.title("IC 2391",fontsize=16)
    plt.savefig("plots/Gemini_detections_limits.png",bbox_inches="tight",dpi=600)

    plt.show()

if __name__=="__main__":

    plot_contrast()
