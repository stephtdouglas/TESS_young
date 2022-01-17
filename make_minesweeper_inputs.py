import os, sys
import glob
import itertools

import numpy as np
import astropy.io.fits as fits
import matplotlib.pyplot as plt
from astropy import table
from astropy.table import join,vstack,Table

from analyze_cluster_output import read_cluster_visual


def make_cluster_table(cluster,date):
    print("\n",cluster,date)

    cat1 = read_cluster_visual(clusters[i],dates[i],return_periodcolor=False)
    cat1.remove_columns(["GAIAEDR3_RA","GAIAEDR3_DEC",
                        "GAIAEDR3_PMRA","GAIAEDR3_PMDEC",#"GAIAEDR3_G",
                        # "GAIAEDR3_BP","GAIAEDR3_RP",
                        "GAIAEDR3_RUWE",
                        "GAIAEDR3_G_CORRECTED","MemBool"])
    if cluster!="IC_2602":
        cat2 = read_cluster_visual(clusters[i],dates[i],return_periodcolor=False,
                                   which=2)

        cat2.remove_columns(["GAIAEDR3_RA","GAIAEDR3_DEC",
                            "GAIAEDR3_PMRA","GAIAEDR3_PMDEC",#"GAIAEDR3_G",
                            "GAIAEDR3_BP","GAIAEDR3_RP",
                            "GAIAEDR3_RUWE",
                            "GAIAEDR3_G_CORRECTED","MemBool"])
        allcat0 = join(cat1,cat2,keys=["TIC"])
    else:
        cat1b = cat1.copy()
        cat1b.remove_columns(["GAIAEDR3_BP","GAIAEDR3_RP"])
        allcat0 = join(cat1,cat1b,keys=["TIC"])
    allcat = table.unique(allcat0,silent=True)

    allcat.rename_column("GAIAEDR3_ID_1","GAIAEDR3_ID")
    allcat.rename_column("GAIAEDR3_G_1","GAIAEDR3_G")
    allcat.remove_column("GAIAEDR3_ID_2")
    allcat.remove_column("GAIAEDR3_G_2")

    any_good = (allcat["final_Q_1"]<=1) | (allcat["final_Q_2"]<=1)
    # print(len(np.where(any_good)[0]))

    consistently_good = (allcat["final_Q_1"]<=1) & (allcat["final_Q_2"]<=1)
    # print(len(np.where(consistently_good)[0]))

    consistently_best = (allcat["final_Q_1"]==0) & (allcat["final_Q_2"]==0)
    # print(len(np.where(consistently_best)[0]))

    # Now match back to my initial catalog files
    cat_init_file = f"tables/{cluster}_crossmatch_allcolumns.fits"
    with fits.open(cat_init_file) as hdu:
        cat_init = Table(hdu[1].data)

    if cat_init.masked == False:
        cat_init = Table(cat_init, masked=True, copy=False)

    cat = join(cat_init,allcat[allcat["GAIAEDR3_ID"].mask==False],
               join_type="left",keys=["GAIAEDR3_ID"])

    cat.rename_column("GAIAEDR3_G_1","GAIAEDR3_G")
    cat.rename_column("GAIAEDR3_BP_1","GAIAEDR3_BP")
    cat.rename_column("GAIAEDR3_RP_1","GAIAEDR3_RP")

    any_good2 = (((cat["final_Q_1"]<=1) | (cat["final_Q_2"]<=1)) &
                 (cat["final_Q_1"].mask==False) & (cat["final_Q_1"].mask==False))

    bp_rp = cat["GAIAEDR3_BP"]-cat["GAIAEDR3_RP"]
    g = cat["GAIAEDR3_G"]
    gmax = max(cat["GAIAEDR3_G"][any_good2])
    plt.plot(bp_rp,g,'.')
    plt.plot(bp_rp[any_good2],g[any_good2],'.')
    plt.ylim(22,0)
    plt.show()

    bp_rp_max = max(bp_rp[any_good2])

    print(gmax)
    print(bp_rp_max)
    print(len(cat))
    print(len(cat[(bp_rp<bp_rp_max)]))
    print(len(cat[g<gmax]))
    print(len(cat[g<(gmax*1.05)]))

    if cluster=="NGC_2547":
        # There's a non-ascii character in there somewhere
        # Easiest just to remove those columns for now
        cat.remove_columns(["Notes_1","Notes_2"])

    cat[g<(gmax*1.05)].write(f"tables/{cluster}_for_minesweeper.fits",
                             overwrite=False)

if __name__=="__main__":
    clusters = ["IC_2391","Collinder_135","NGC_2451A","NGC_2547","IC_2602"]
    dates = ["2021-06-22","2021-06-18","2021-06-21","2021-06-21","2021-07-02"]

    for i in range(5):
        make_cluster_table(clusters[i],dates[i])
        # break
