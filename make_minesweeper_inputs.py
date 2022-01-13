import os, sys
import glob
import itertools

import numpy as np
import astropy.io.fits as fits
from astropy import table
from astropy.table import join,vstack,Table

from analyze_cluster_output import read_cluster_visual


def make_cluster_table(cluster,date):
    print("\n",cluster,date)

    cat1 = read_cluster_visual(clusters[i],dates[i],return_periodcolor=False)
    cat1.remove_columns(["GAIAEDR3_RA","GAIAEDR3_DEC",
                        "GAIAEDR3_PMRA","GAIAEDR3_PMDEC",#"GAIAEDR3_G",
                        "GAIAEDR3_BP","GAIAEDR3_RP","GAIAEDR3_RUWE",
                        "GAIAEDR3_G_CORRECTED","MemBool"])
    if cluster!="IC_2602":
        cat2 = read_cluster_visual(clusters[i],dates[i],return_periodcolor=False,
                                   which=2)

        cat2.remove_columns(["GAIAEDR3_RA","GAIAEDR3_DEC",
                            "GAIAEDR3_PMRA","GAIAEDR3_PMDEC",#"GAIAEDR3_G",
                            # "GAIAEDR3_BP","GAIAEDR3_RP",
                            "GAIAEDR3_RUWE",
                            "GAIAEDR3_G_CORRECTED","MemBool"])
        allcat0 = join(cat1,cat2,keys=["TIC"])
    else:
        allcat0 = join(cat1,cat1,keys=["TIC"])
    allcat = table.unique(allcat0,silent=True)

    allcat.rename_column("GAIAEDR3_ID_1","GAIAEDR3_ID")
    allcat.rename_column("GAIAEDR3_G_1","GAIAEDR3_G")
    allcat.remove_column("GAIAEDR3_ID_2")
    allcat.remove_column("GAIAEDR3_G_2")

    any_good = (allcat["final_Q_1"]<=1) | (allcat["final_Q_2"]<=1)
    print(len(np.where(any_good)[0]))

    consistently_good = (allcat["final_Q_1"]<=1) & (allcat["final_Q_2"]<=1)
    print(len(np.where(consistently_good)[0]))

    consistently_best = (allcat["final_Q_1"]==0) & (allcat["final_Q_2"]==0)
    print(len(np.where(consistently_best)[0]))

    print(len(np.where(allcat["GAIAEDR3_ID"].mask==True)[0]))

    # print(allcat["final_period_1","final_period_2"][allcat["GAIAEDR3_ID"].mask==True])


    bp_rp = allcat["GAIAEDR3_BP"]-allcat["GAIAEDR3_RP"]
    print(max(bp_rp[any_good]))
    print(max(allcat["GAIAEDR3_G"][any_good]))

    # Now match back to the HDBScan files, for those stars that have Gaia IDs
    hdbscanfile = os.path.expanduser(f"~/Dropbox/EDR3/scats/{cluster}.fits")
    with fits.open(hdbscanfile) as hdu:
        hdbscan = Table(hdu[1].data)
        # hdbscan_memb = hdbscan[hdbscan["MemBool"]==1]

    if hdbscan.masked == False:
        hdbscan = Table(hdbscan, masked=True, copy=False)
    hdbscan.rename_column("GAIADER3_ID","GAIAEDR3_ID")
    # print(len(hdbscan))
    # print(hdbscan.dtype)

    cat = join(hdbscan,allcat[allcat["GAIAEDR3_ID"].mask==False],
               join_type="left",keys=["GAIAEDR3_ID"])
    # print(cat.dtype)

    # print(cat["TIC"])
    bp_rp = cat["GAIAEDR3_BP_1"]-cat["GAIAEDR3_RP_1"]
    to_out = (cat["TIC"].mask==False)
    print(len(cat[to_out]))
    print(len(cat[to_out & (bp_rp<3.5)]))


if __name__=="__main__":
    clusters = ["IC_2391","Collinder_135","NGC_2451A","NGC_2547","IC_2602"]
    dates = ["2021-06-22","2021-06-18","2021-06-21","2021-06-21","2021-07-02"]

    for i in range(4):
        make_cluster_table(clusters[i],dates[i])
        # break




    # vis_file1 = f"tables/{cluster}_{date}_results_comments.csv"
    # vis_file2 = f"tables/{cluster}_{date}_results_comments2.csv"
    # vis1 = at.read(vis_file1,delimiter=",")
    # vis2 = at.read(vis_file2,delimiter=",")
    # vis1.rename_column("\ufefftarget_name","TIC")
    # vis2.rename_column("\ufefftarget_name","TIC")
    #
    # for vis in [vis1,vis2]
    #     good2 = np.where(vis["Select"].mask==False)[0]
    #     vis["final_period"] = np.copy(vis["sig_periods"])
    #     vis["final_Q"] = np.copy(vis["Q"])
    #     replace2 = (vis["Q"]==2) & ((vis["Q2"]==1) | (vis["Q2"]==0))
    #     replace3 = (vis["Q"]==2) & ((vis["Q3"]==1) | (vis["Q3"]==0))
    #     vis["final_period"][replace2] = vis["sec_periods"][replace2]
    #     vis["final_Q"][replace2]==vis["Q2"][replace2]
    #     vis["final_period"][replace3] = -99 # Not bothering with peaks file right now
    #     vis["final_Q"][replace3]==vis["Q3"][replace3]
    #
    # allvis = join(vis1,vis2,keys=["TIC","provenance_name","flux_cols","sequence_number"])
    #
    # x_file = f"{cluster}_crossmatch_xmatch_TIC.csv"
    # xmatch = at.read(x_file,delimiter=",")
