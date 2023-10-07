import os, sys, pathlib
import glob
import itertools

import numpy as np
from astropy.table import join,vstack,Table
import astropy.io.fits as fits
import astropy.io.ascii as at


from tess_young.get_const import *
import tess_young
_DIR = pathlib.Path(tess_young.__file__).resolve().parent.parent


def check_inputs_outputs(cluster,date,allcat):
    matched_file = os.path.join(_DIR,f"tables/{cluster}_crossmatch.csv")
    tic_file = os.path.join(_DIR,f"{cluster}_crossmatch_xmatch_TIC.csv")
    cat_init_file2 = os.path.expanduser(f"~/Dropbox/data/MINESweeper/catalog_{cluster}_v0.fits")

    mcat = at.read(matched_file)
    tcat = at.read(tic_file)
    with fits.open(cat_init_file2) as hdu:
        cat_init2 = Table(hdu[1].data)
        # print(cat_init.dtype)
        cat_init2.rename_column("GaiaEDR3_ID","GAIAEDR3_ID")

    print(len(mcat),"initial crossmatch")
    print(len(np.unique(mcat["GAIAEDR3_ID"])),"unique Gaia IDs")
    print(len(tcat),"TIC crossmatch")
    print(len(np.unique(tcat["GAIAEDR3_ID"])),"unique Gaia IDs")
    print(len(cat_init2),"MINESweeper")
    print(len(np.unique(cat_init2["GAIAEDR3_ID"])),"unique Gaia IDs")

    ncat = len(np.where(allcat["Cluster"]==cluster)[0])
    has_periods = (allcat["Prot1"]>0) & (allcat["Prot1"].mask==False)
    nper = len(np.where((allcat["Cluster"]==cluster) & has_periods)[0])
    has_lit = (allcat["LitPeriod"].mask==False)
    nlit = len(np.where((allcat["Cluster"]==cluster) & has_lit)[0])
    new_periods = has_periods & (has_lit==False)
    nnew = len(np.where((allcat["Cluster"]==cluster) & new_periods)[0])

    print(ncat,"in final catalog")
    print(nper,"in final catalog with periods")
    print(len(np.unique(allcat[allcat["Cluster"]==cluster]["GAIAEDR3_ID"])),"unique Gaia IDs")
    print(nlit,"literature periods")
    print(nnew,"new periods from TESS")

    return(ncat,nper)

if __name__=="__main__":
    clusters = ["IC_2391","Collinder_135","NGC_2451A","NGC_2547","IC_2602"]
    dates = ["2021-06-22","2021-06-18","2021-06-21","2021-06-21","2021-07-02"]

    final_file = os.path.join(_DIR,"tab_all_stars_orig.csv")
    allcat = at.read(final_file)

    ntot, ptot = 0, 0
    for i in range(5):
        print("\n----------------------------------------\n",clusters[i])
        ncat, nper = check_inputs_outputs(clusters[i],dates[i],allcat)
        ntot += ncat
        ptot += nper

    print("\n----------------------------------------\nALL")
    print(ntot, "in catalog", ptot, "with periods")

    print(len(allcat),"in catalog",
          len(np.where(allcat["Q1"]<=1)[0]), "with periods")
    print(len(np.where(allcat["Q1"]<7)[0]), "with TESS data")
    print(len(np.unique(allcat["GAIAEDR3_ID"])),"unique Gaia IDs")
    print(len(np.where((allcat["Q1"]<7)&(allcat["to_plot"]==1))[0]), "with TESS data and good membership")
    print(len(np.where(allcat["to_plot"]==1)[0]),"good membership",
          len(np.where((allcat["Q1"]<=1) & (allcat["to_plot"]==1))[0]), "with periods")

    print(len(np.where((allcat["to_plot"]==0) & (allcat["LitPeriod"].mask==False))[0]), "with lit periods, only one catalog")
    print(len(np.where((allcat["to_plot"]==1) & (allcat["LitPeriod"].mask==False))[0]), "with lit periods, good for plotting")
    print("\n")

    solar = (allcat["Mass"]>=0.9) & (allcat["Mass"]<=1.1)
    above09 = (allcat["Mass"]>=0.9) 
    print(len(np.where((allcat["to_plot"]==1) & solar)[0]),"good membership solar\n",
          len(np.where((allcat["Q1"]<=1) & (allcat["to_plot"]==1) & solar)[0]), "solar with periods")
    print(len(np.where((allcat["to_plot"]==1) & (allcat["LitPeriod"]>0) & solar)[0]), "solar with lit periods")
    print(len(np.where((allcat["to_plot"]==1) & (allcat["Q1"]<=1) & (allcat["LitPeriod"]<=0) & solar)[0]), "solar new periods")

    print("\n")
    print(len(np.where((allcat["to_plot"]==1) & above09)[0]),"good membership 'high' mass\n",
          len(np.where((allcat["Q1"]<=1) & (allcat["to_plot"]==1) & above09)[0]), "'high' mass with periods")
    print(len(np.where((allcat["to_plot"]==1) & (allcat["LitPeriod"]>0) & above09)[0]), "'high' mass with lit periods")
    print(len(np.where((allcat["to_plot"]==1) & (allcat["Q1"]<=1) & (allcat["LitPeriod"]<=0) & above09)[0]), "'high' mass new periods")


    print("\n")
    print(len(np.where(allcat["Q1"]<=1)[0]), "total new periods, any membership")
    print(len(np.where((allcat["Q1"]<=1) & (allcat["LitPeriod"]<=0))[0]), "total new periods, any membership")
    print(len(np.where((allcat["to_plot"]==1) & (allcat["Q1"]<=1) & (allcat["LitPeriod"]<=0))[0]), "total new periods, good membership")
    print(len(np.where((allcat["to_plot"]==1) & (allcat["Q1"]<=1) & (allcat["LitPeriod"]>0))[0]), "lit overlap periods, good membership")
    print(len(np.where((allcat["Q1"]<=1) & (allcat["LitPeriod"]<=0))[0]), "total new periods")
    print(len(np.where((allcat["Q1"]<=1) & (allcat["LitPeriod"]>0))[0]), "lit overlap periods")


    print("\n")
    print(len(np.where((allcat["Q1"]<=1) & (allcat["Bl?"]=="y"))[0]), "TESS periods, yes blended")
    print(len(np.where((allcat["Q1"]<=1) & (allcat["Bl?"]=="m"))[0]), "TESS periods, maybe blended")
    print(len(np.where((allcat["Q1"]<=1) & (allcat["Bl?"]=="n"))[0]), "TESS periods, not blended")

    print(len(np.where((allcat["Q1"]<=6) & (allcat["Bl?"]=="y"))[0]), "TESS targets, yes blended")
    print(len(np.where((allcat["Q1"]<=6) & (allcat["Bl?"]=="m"))[0]), "TESS targets, maybe blended")
    print(len(np.where((allcat["Q1"]<=6) & (allcat["Bl?"]=="n"))[0]), "TESS targets, not blended")
