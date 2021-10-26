"""
Script to plot literature periods for motivation
"""

import os, sys

import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.cm as cm
import numpy as np
import astropy.io.ascii as at
# import astropy.io.fits as fits
# from astropy.coordinates import SkyCoord
# from astroquery.simbad import Simbad
# from astroquery.vizier import Vizier
from astroquery.gaia import Gaia
from astropy.table import join, Table, vstack
# from astroquery.xmatch import XMatch
from astropy import units as u

# from analyze_cluster_output import process_cluster, read_cluster_visual


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


def plot_literature_periodcolor(ax=None):

    basic_query = ("SELECT TOP 5 gaia_source.source_id,gaia_source.ra,"
                  "gaia_source.ra_error,gaia_source.dec,gaia_source.dec_error,"
                  "gaia_source.parallax,gaia_source.parallax_error,"
                  "gaia_source.phot_g_mean_mag,gaia_source.bp_rp "
                  "FROM gaiadr2.gaia_source ")
    if ax is None:
        plt.figure()
        ax = plt.subplot(111)
    lit_color = "#595959"

    ############################################################################
    ############################################################################
    # IC 2391
    cluster = "IC_2391"
    print(cluster,"\n-------")
    simbadfile = "IC2391_rotation_patten1996_simbad.csv"
    simbad = at.read(simbadfile, delimiter=",")
    period_col="Period"

    query_results = []

    for row in simbad:
        print(row["GaiaDR2"])
        job = Gaia.launch_job(f"{basic_query} WHERE (gaiadr2.gaia_source.source_id={row['GaiaDR2']}) ")
        r = job.get_results()
        print(r)
        query_results.append(r)

    gaia = vstack(query_results)
    gaia.rename_column("source_id","GaiaDR2")

    match = join(simbad,gaia,join_type="left",keys=["GaiaDR2"],
                 table_names=["lit","Gaia"])

    bp_rp = match["bp_rp"]
    ax.plot(bp_rp,match[period_col],shapes[cluster],color=lit_color,ms=6,
             label="IC 2391: Patten & Simon (1996)")

    ############################################################################
    ############################################################################
    # NGC 2547
    cluster = "NGC_2547"
    print(cluster,"\n-------")
    simbadfile = "NGC2547_rotation_irwin2008_simbad.csv"
    simbad = at.read(simbadfile, delimiter=",")

    query_results = []

    for row in simbad[simbad["GaiaDR2"].mask==False]:
        print(row["GaiaDR2"])
        job = Gaia.launch_job(f"{basic_query} WHERE (gaiadr2.gaia_source.source_id={row['GaiaDR2']}) ")
        r = job.get_results()
        print(r)
        query_results.append(r)

    gaia = vstack(query_results)
    gaia.rename_column("source_id","GaiaDR2")

    match = join(simbad[simbad["GaiaDR2"].mask==False],gaia,join_type="left",keys=["GaiaDR2"],
                 table_names=["lit","Gaia"])

    bp_rp = match["bp_rp"]
    ax.plot(bp_rp,match["Per"],shapes[cluster],color=lit_color,ms=6,
             label="NGC 2547: Irwin+ (2008)")

    ############################################################################
    ############################################################################
    # IC 2602
    cluster = "IC_2602"
    print(cluster,"\n-------")

    simbadfile1 = "IC2602_rotation_barnes1999_simbad.csv"
    simbad1= at.read(simbadfile1, delimiter=",")
    # good1 = (simbad1["GaiaDR2"].mask==False) & (simbad1["TIC"].mask==False)

    simbadfile2 = "IC2602_rotation_tschape2001_simbad.csv"
    simbad2= at.read(simbadfile2, delimiter=",")
    simbad = join(simbad1,simbad2,join_type="outer",keys=["GaiaDR2","TIC"],
                  table_names=["barnes","tschape"])
    print(simbad)

    print(np.where((simbad["Prot_barnes"].mask==False) & (simbad["Prot_tschape"].mask==False))[0])

    query_results = []

    for row in simbad:
        print(row["GaiaDR2"])
        job = Gaia.launch_job(f"{basic_query} WHERE (gaiadr2.gaia_source.source_id={row['GaiaDR2']}) ")
        r = job.get_results()
        print(r)
        query_results.append(r)

    gaia = vstack(query_results)
    gaia.rename_column("source_id","GaiaDR2")

    match = join(simbad,gaia,join_type="left",keys=["GaiaDR2"],
                 table_names=["lit","Gaia"])

    bp_rp = match["bp_rp"]
    ax.plot(bp_rp,match["Prot_barnes"],shapes[cluster],color=lit_color,ms=6,
             label="IC 2602: Barnes+ (1999)")
    ax.plot(bp_rp,match["Prot_barnes"],"D",color=lit_color,ms=6,
             label="IC 2602: Tschape & Rudiger (2001)")


    ax.legend(loc=2)

    ax.set_ylim(0.1,50)
    ax.set_xlim(0.5,3.5)
    ax.set_yscale("log")

    ax.set_xlabel(r"G$_{BP}$ - G$_{RP}$ (DR2)")
    ax.set_ylabel("Period (d)")

    return ax

def plot_literature_periodcolor_comparisons(ax=None):
    if ax is None:
        plt.figure()
        ax = plt.subplot(111)


if __name__=="__main__":
    _ = plot_literature_periodcolor()
    plt.savefig(f"plots/periodcolor_literature.png")
    plt.close("all")
