import os, sys, pathlib
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

from tess_young.get_const import *
from tess_young import plot_spt_axis
import tess_young
_DIR = pathlib.Path(tess_young.__file__).resolve().parent.parent
plt.style.use(os.path.join(_DIR,'paper.mplstyle'))

mamajek_file = os.path.expanduser("~/Dropbox/Models/mamajek_colors.dat")
mamajek = at.read(mamajek_file,fill_values=[("...","NaN")])
mamajek_vmi_to_mass = interp1d(mamajek["V-Ic"],mamajek["Msun"],bounds_error=False)
mamajek_vmk_to_mass = interp1d(mamajek["V-Ks"],mamajek["Msun"],bounds_error=False)

def vmi_to_mass(v_minus_i):
    return mamajek_vmi_to_mass(v_minus_i)

def vmk_to_mass(v_minus_k):
    return mamajek_vmk_to_mass(v_minus_k)


def setup_axes(ax):
    ax.set_xlim(1.3,0.1)
    ax.set_ylim(0.1,50)
    ax.set_yscale("log")
    ax.set_xlabel(r"Mass (M$_{\odot}$)",fontsize=16)
    ax.set_ylabel("Period (d)",fontsize=16)
    ax.tick_params(labelsize=14)


if __name__=="__main__":
    clusters = ["IC_2391","Collinder_135","NGC_2451A","NGC_2547","IC_2602"]
    dates = ["2021-06-22","2021-06-18","2021-06-21","2021-06-21","2021-07-02"]

    ## h Per ##############################################
    print("\nhPer")

    cat_file = os.path.expanduser("~/Dropbox/data/catalogs/hper_rotation_moraux2013.tsv")
    cat = at.read(cat_file,delimiter="|",data_start=3)
    cat = Table(cat, masked=True, copy=False)

    plt.figure()
    ax = plt.subplot(111)
    setup_axes(ax)
    ax.plot(cat["Mass"],cat["Per"],'.',color="Grey",label="h Per (13 Myr; Moraux+2013)")
    ax.legend(loc=2)

    plt.savefig("plots/periodmass_intro_hper.png",bbox_inches="tight",dpi=600)
    plt.close()

    ## ZAMS ##############################################
    print("\nZAMS")

    files = {"patten1996":"IC2391_rotation_patten1996_simbad.csv",
             "barnes1999":"IC2602_rotation_barnes1999_simbad.csv",
             "prosser":"IC2602_rotation_prosserstauffer_simbad.csv",
             "irwin2008":"NGC2547_rotation_irwin2008_simbad.csv",
             "messina2011":"IC2391_rotation_messina2011_simbad.csv"}

    labels = {"patten1996":"IC 2391 (Patten & Simon 1996)",
              "barnes1999":"IC 2602 (Barnes+ 1999)",
              "prosser":"IC 2602 (Patten+ 1996)",
              "irwin2008":"NGC 2547 (Irwin+ 2008)",
             "messina2011":"IC 2391 (Messina+ 2011)"}
    clusters = {"patten1996":"IC_2391",
              "barnes1999":"IC_2602",
              "prosser":"IC_2602",
              "irwin2008":"NGC_2547",
             "messina2011":"IC_2391"}

    plt.figure()
    ax = plt.subplot(111)
    setup_axes(ax)

    for cluster in ["IC_2391","NGC_2547","IC_2602"]:
        ax.plot([],[],shapes[cluster],color=colors[cluster],label=cluster.replace("_"," "))
    ax.legend(loc=2)

    for i,source in enumerate(labels.keys()):
        fname = files[source]
        cluster = clusters[source]
        zorder = 5-i

        dat = at.read(fname)

        if source=="patten1996":
            mass = mamajek_vmi_to_mass(dat["V-I"])
            period = dat["Period"]
        elif source=="barnes1999":
            mass = mamajek_vmi_to_mass(dat["(V-I)e"])
            period = dat["Prot"]
        elif source=="prosser":
            use = dat["REF"]==1
            mass = mamajek_vmi_to_mass(dat["V-Ic"][use])
            period = dat["P(days)"][use]
        elif source=="messina2011":
            mass = mamajek_vmi_to_mass(dat["(V-I)"])
            period = dat["Per"]
        else:
            mass = dat["Mass"]
            period = dat["Per"]

        # if source=="prosser":
        #     ax.plot(mass,period,'D',
        #              mec="#174f34",mfc="none",mew=1.5,zorder=zorder,
        #              label=labels[source])
        # else:
        ax.plot(mass,period,shapes[cluster],
                  color=colors[cluster],label=labels[source],zorder=zorder)
    # ax.legend(loc=2)
    plt.savefig("plots/periodmass_intro_zams.png",bbox_inches="tight",dpi=600)
    plt.close()


    ## Pleiades ##############################################
    print("\nPleiades")

    cat_file = os.path.expanduser("~/Dropbox/data/catalogs/pleiades_rotation_rebull2016.csv")
    cat = at.read(cat_file)
    cat = Table(cat, masked=False, copy=False)

    plt.figure()
    ax = plt.subplot(111)
    setup_axes(ax)
    mass = vmk_to_mass(cat["(V-K)0"])
    ax.plot(mass,cat["Per1"],'.',color="Grey",label="Pleiades (125 Myr; Rebull+2016)")
    ax.legend(loc=2)

    plt.savefig("plots/periodmass_intro_pleiades.png",bbox_inches="tight",dpi=600)
    plt.close()
