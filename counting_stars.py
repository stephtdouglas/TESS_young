import os, sys
import glob
import itertools

import numpy as np
from astropy.table import join,vstack,Table
import astropy.io.fits as fits
import astropy.io.ascii as at


# import warnings
# with warnings.catch_warnings():
#     from astropy.table import VerifyWarning
#     warnings.filterwarnings("ignore", category=VerifyWarning)
#     import astropy.io.ascii as at
#     import astropy.io.fits as fits
#     from astropy.table import join,vstack,Table


def check_inputs_outputs(cluster,date,allcat):
    matched_file = f"tables/{cluster}_crossmatch.csv"
    tic_file = f"{cluster}_crossmatch_xmatch_TIC.csv"
    cat_init_file2 = os.path.expanduser(f"~/Dropbox/data/MINESweeper/catalog_{cluster}_v0.fits")

    mcat = at.read(matched_file)
    tcat = at.read(tic_file)
    with fits.open(cat_init_file2) as hdu:
        cat_init2 = Table(hdu[1].data)
        # print(cat_init.dtype)
        cat_init2.rename_column("GaiaEDR3_ID","GAIAEDR3_ID")

    print(len(mcat),"initial crossmatch")
    print(len(np.unique(mcat["GAIAEDR3_ID"])))
    print(len(tcat),"TIC crossmatch")
    print(len(np.unique(tcat["GAIAEDR3_ID"])))
    print(len(cat_init2),"MINESweeper")
    print(len(np.unique(cat_init2["GAIAEDR3_ID"])))

    ncat = len(np.where(allcat["Cluster"]==cluster)[0])
    has_periods = (allcat["Prot1"]>0) & (allcat["Prot1"].mask==False)
    nper = len(np.where((allcat["Cluster"]==cluster) & has_periods)[0])

    print(ncat,"in final catalog")
    print(nper,"in final catalog with periods")
    print(len(np.unique(allcat[allcat["Cluster"]==cluster]["GAIAEDR3_ID"])))

    return(ncat,nper)

if __name__=="__main__":
    clusters = ["IC_2391","Collinder_135","NGC_2451A","NGC_2547","IC_2602"]
    dates = ["2021-06-22","2021-06-18","2021-06-21","2021-06-21","2021-07-02"]

    final_file = "tab_all_stars.csv"
    allcat = at.read(final_file)

    ntot, ptot = 0, 0
    for i in range(5):
        print("\n\n",clusters[i])
        ncat, nper = check_inputs_outputs(clusters[i],dates[i],allcat)
        ntot += ncat
        ptot += nper

    print("\n\nALL")
    print(ntot, "in catalog", ptot, "with periods")

    print(len(allcat),"in catalog",
          len(np.where(allcat["Q1"]<=1)[0]), "with periods")
    print(len(np.where(allcat["to_plot"]==1)[0]),"good membership",
          len(np.where((allcat["Q1"]<=1) & (allcat["to_plot"]==1))[0]), "with periods")

    print(len(np.where((allcat["to_plot"]==1) & (allcat["LitPeriod"]>0))[0]), "with lit periods")

    solar = (allcat["Mass"]>=0.9) & (allcat["Mass"]<=1.1)
    print(len(np.where((allcat["to_plot"]==1) & solar)[0]),"good membership solar",
          len(np.where((allcat["Q1"]<=1) & (allcat["to_plot"]==1) & solar)[0]), "solar with periods")
    print(len(np.where((allcat["to_plot"]==1) & (allcat["LitPeriod"]>0) & solar)[0]), "solar with lit periods")
