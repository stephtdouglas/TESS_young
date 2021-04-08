"""
Script to crossmatch open cluster catalogs
"""

import os, sys
from datetime import date

import matplotlib.pyplot as plt
import numpy as np
import astropy.io.ascii as at
import astropy.io.fits as fits
from astropy.table import join, Table
from astroquery.xmatch import XMatch
from astropy import units as u

def xmatch_ic2391(to_plot=False):
    # Phill's catalog, which also includes all the other Gaia information
    hdbscanfile = os.path.expanduser("~/Dropbox/EDR3/scats/IC_2391.fits")
    with fits.open(hdbscanfile) as hdu:
        hdbscan = Table(hdu[1].data)
        # hdbscan_memb = hdbscan[hdbscan["MemBool"]==1]

    hdbscan.rename_column("GAIADER3_ID","GAIAEDR3_ID")
    # print(hdbscan.dtype)

    # Gaia-ESO Survey catalog (Jackson+2020)
    ges0 = at.read("catalogs/Jackson2020_table2_GaiaEDR3_xmatch.csv")
    ges1 = at.read("catalogs/Jackson2020_table3_GaiaDR2_data.txt",
                  header_start=1,data_start=4)
    ges = join(ges0,ges1,keys=["target","filter","cluster"])
    ges = ges[ges["cluster"]=="IC2391"]
    # print(ges.dtype.names)
    ges.rename_column("source_id","GAIAEDR3_ID")
    # print(ges["GAIADR3_ID"])

    # Cantat-Gaudin et al. 2020
    cantat = at.read("catalogs/cantat-gaudin2020_ic2391_10deg_GaiaEDR3_xmatch.csv",
                     data_start=3)
    cantat.rename_column("source_id","GAIAEDR3_ID")
    cantat = cantat[cantat["Cluster"]=="IC_2391"]
    # print(cantat.dtype.names)
    # print(cantat["GAIAEDR3_ID"])


    # First join the GES and Cantat-Gaudin catalogs
    cg = join(ges,cantat,keys=["GAIAEDR3_ID"],join_type="outer",
              table_names=["GES","Cantat-Gaudin"])
    # print(cg.dtype)
    # print(cg["GAIAEDR3_ID"])

    # Then join those to the HDBScan catalog
    allcat = join(hdbscan,cg,keys=["GAIAEDR3_ID"],join_type="outer")
    # print(allcat.dtype)
    # print(allcat["GAIAEDR3_ID"])


    bp_rp = allcat["GAIAEDR3_BP"] - allcat["GAIAEDR3_RP"]
    hdb_memb = (allcat["MemBool"]==1)  & (allcat["MemBool"].mask==False)
    ges_memb = (allcat["prob_p"]>0) & (allcat["prob_p"].mask==False) #TODO: make sure this is the right cut-off
    can_memb = (allcat["proba"]>0) & (allcat["proba"].mask==False)

    print(allcat["GAIAEDR3_PMRA","GAIAEDR3_PMDEC","GAIAEDR3_PARALLAX"][ges_memb])
    print(allcat["GAIAEDR3_G_CORRECTED","GAIAEDR3_BP","GAIAEDR3_RP"][ges_memb])

    print(allcat["GAIAEDR3_PMRA","GAIAEDR3_PMDEC","GAIAEDR3_PARALLAX"][can_memb])
    print(allcat["GAIAEDR3_G_CORRECTED","GAIAEDR3_BP","GAIAEDR3_RP"][can_memb])

    all_memb = hdb_memb & ges_memb & can_memb
    any_memb = hdb_memb | ges_memb | can_memb
    print(len(np.where(all_memb)[0]),len(np.where(any_memb)[0]))

    hdbscan_lowq = (allcat["HDBscan_Cluster"]==-2) & (allcat["HDBscan_Cluster"].mask==False)

    # Check for "low quality" members in the same region of the diagrams
    # as the confirmed members
    pmra_lims = [np.nanmin(allcat["GAIAEDR3_PMRA"][hdb_memb]),
                 np.nanmax(allcat["GAIAEDR3_PMRA"][hdb_memb])]
    pmdec_lims = [np.nanmin(allcat["GAIAEDR3_PMDEC"][hdb_memb]),
                 np.nanmax(allcat["GAIAEDR3_PMDEC"][hdb_memb])]
    plx_lims = [np.nanmin(allcat["GAIAEDR3_PARALLAX"][hdb_memb]),
                 np.nanmax(allcat["GAIAEDR3_PARALLAX"][hdb_memb])]
    print(pmra_lims,pmdec_lims,plx_lims)
    lowq_close = (hdbscan_lowq & (allcat["GAIAEDR3_PMRA"]>pmra_lims[0])
                & (allcat["GAIAEDR3_PMRA"]<pmra_lims[1])
                & (allcat["GAIAEDR3_PMDEC"]>pmdec_lims[0])
                & (allcat["GAIAEDR3_PMDEC"]<pmdec_lims[1])
                & (allcat["GAIAEDR3_PARALLAX"]>plx_lims[0])
                & (allcat["GAIAEDR3_PARALLAX"]<plx_lims[1])
                & (allcat["GAIAEDR3_PMRA"].mask==False)
                & (allcat["GAIAEDR3_PMDEC"].mask==False)
                & (allcat["GAIAEDR3_PARALLAX"].mask==False)
                )
    # lowq_close = lowq_close & ~any_memb
    print(np.where(lowq_close)[0])
    print(allcat["GAIAEDR3_PMRA","GAIAEDR3_PMDEC","GAIAEDR3_PARALLAX"][lowq_close])
    print(allcat["GAIAEDR3_G_CORRECTED","GAIAEDR3_BP","GAIAEDR3_RP"][lowq_close])

    # Print overlap information
    lists = [hdb_memb,ges_memb,can_memb,hdbscan_lowq]
    names = ["HDBScan","GES","Cantat-Gaudin","HDBscan Low Quality"]

    for i in range(3):
        uniq_list = np.where(lists[i])[0]
        print(f"{len(uniq_list)} members in {names[i]}")
        for j in range(4):
            if i!=j:
                joint_list = lists[i] & lists[j]
                print(f"\t{len(np.where(joint_list)[0])} members also in {names[j]}")
                if j<4:
                    uniq_list = np.setdiff1d(uniq_list,np.where(lists[j])[0])
        print(f"\t{len(uniq_list)} unique members")


    # fix columns in the combined catalog
    col_match = at.read("hdbscan_edr3_cols.csv",delimiter=",")
    print(col_match)
    for i in [1,2]:
        name = names[i]
        memb = np.where(lists[i])[0]
        for k in range(len(col_match)):
            colname = f"{col_match['edr3_col'][k]}_{name}"
            for j in memb:
                hdb_colname = col_match["hdbscan_col"][k]
                allcat[hdb_colname][j] = allcat[colname][j]
                allcat[colname].mask[j] = False
            allcat.remove_column(colname)
    print(allcat.dtype.names)

    # For now, will just have to output the file and run it through xmatch on my own. Sigh.
    # print("Starting Xmatch")
    # # Now crossmatch to the TIC to get TIC IDs
    # allpos = allcat["GAIAEDR3_RA","GAIAEDR3_DEC","GAIAEDR3_ID"]
    # allpos.rename_column("GAIAEDR3_RA","ra")
    # allpos.rename_column("GAIAEDR3_DEC","dec")
    # print(allpos)
    # table = XMatch.query(cat1=allpos,
    #                      cat2='vizier:IV/38/tic',
    #                      max_distance=1 * u.arcsec,
    #                      colRA1='ra',colDec1='dec',cache=True)
    #                      # colRA1='GAIAEDR3_RA',colDec1='GAIAEDR3_DEC')
    # print("Finished Xmatch")
    # print(table.dtype.names)
    # table.rename_column("ra","GAIAEDR3_RA")
    # table.rename_column("dec","GAIAEDR3_DEC")
    # for colname in table.dtype.names:
    #     if ("GAIAEDR3" in colname) or ("TIC" in colname):
    #         continue
    #     else:
    #         table.rename_column(colname,f"{colname}_TIC")
    #
    # finalcat = join(allcat,table,keys=["GAIAEDR3_ID","GAIAEDR3_RA","GAIAEDR3_DEC"])
    # print(len(allcat),len(table),len(finalcat))

    # For now, will just have to output the file and run it through xmatch on my own. Sigh.
    finalcat = allcat

    # sys.exit(0)



    # Write out the catalog to two files: a full catalog with all columns,
    # and a sub-catalog to use for further analysis
    subcat = finalcat['TMASS_ID', 'UKIDSS_ID', 'GAIAEDR3_ID', 'GAIADR2_G',
                    'GAIADR2_BP', 'GAIADR2_RP', 'GAIAEDR3_G', 'GAIAEDR3_BP',
                    'GAIAEDR3_RP', 'GAIADR2_PARALLAX', 'GAIADR2_PARALLAX_ERROR',
                    'GAIADR2_RUWE', 'GAIAEDR3_RA', 'GAIAEDR3_DEC',
                    'GAIAEDR3_PARALLAX', 'GAIAEDR3_PARALLAX_ERROR',
                    'GAIAEDR3_PMRA', 'GAIAEDR3_PMDEC', 'GAIAEDR3_PMRA_ERROR',
                    'GAIAEDR3_PMDEC_ERROR', 'GAIAEDR3_RUWE',
                    'GAIAEDR3_G_CORRECTED', 'GAIAEDR3_PARALLAX_CORRECTED',
                    'HDBscan_MemProb', 'HDBscan_Cluster', 'HDBscan_Stability',
                    'MemBool', 'angDist_GES', 'target', 'filter', 'cluster',
                    'S/N', 'Teff', 'logg', 'gamma', 'Ks', 'RV_1', 'e_RV',
                    'logL', 'ra_epoch2000_GES', 'dec_epoch2000_GES', 'PLX',
                    'e_PLX', 'VRA', 'e_VRA', 'VDec', 'e_VDec', 'RV_2', 'SRV',
                    'Gflag', 'prob_p', 'p_filter', 'angDist_Cantat-Gaudin',
                    'RA_ICRS', 'DE_ICRS', 'GaiaDR2', 'Plx', 'pmRA*', 'pmDE',
                    'RV', 'o_Gmag', 'Gmag', 'BP-RP', 'proba', 'Cluster',
                    'Teff50',
                    'ra_epoch2000_Cantat-Gaudin', 'dec_epoch2000_Cantat-Gaudin']
                    # ,'TIC','angDist_TIC']

    subcat[any_memb].write("IC2391_crossmatch.fits",overwrite=True)
    subcat[any_memb].write("IC2391_crossmatch.csv",overwrite=True)
    # finalcat[any_memb].write("IC2391_crossmatch_allcolumns.fits",overwrite=True)


    if to_plot:
        # Plot astrometry for confirmed members and for "low quality" stars that
        # were left out of the HDBScan run
        plt.figure(figsize=(10,10))
        ax1 = plt.subplot(221)
        for i,memb in enumerate(lists[:3]):
            plt.plot(finalcat["GAIAEDR3_PMDEC"][memb],finalcat["GAIAEDR3_PARALLAX"][memb],'.',
                    label=names[i],alpha=0.75)
        plt.plot(finalcat["GAIAEDR3_PMDEC"][hdbscan_lowq],
                 finalcat["GAIAEDR3_PARALLAX"][hdbscan_lowq],'k.',label="Quality cut",
                 alpha=0.5)
        plt.ylabel("Parallax (mas)")
        ax1.set_ylim(-20,15)
        plt.legend(loc=3)

        ax2 = plt.subplot(223,sharex=ax1)
        for i,memb in enumerate(lists[:3]):
            plt.plot(finalcat["GAIAEDR3_PMDEC"][memb],finalcat["GAIAEDR3_PMRA"][memb],'.',
                    label=names[i],alpha=0.75)
        plt.plot(finalcat["GAIAEDR3_PMDEC"][hdbscan_lowq],
                 finalcat["GAIAEDR3_PMRA"][hdbscan_lowq],'k.',label="Quality cut",
                 alpha=0.25)
        plt.ylabel("PMRA (mas/yr)")
        plt.xlabel("PMDEC (mas/yr)")
        ax1.set_xlim(-50,50)
        ax2.set_ylim(-40,25)

        ax3 = plt.subplot(224,sharey=ax2)
        for i,memb in enumerate(lists[:3]):
            plt.plot(finalcat["GAIAEDR3_PARALLAX"][memb],
                    finalcat["GAIAEDR3_PMRA"][memb],'.',
                    label=names[i],alpha=0.75)
        plt.plot(finalcat["GAIAEDR3_PARALLAX"][hdbscan_lowq],
                 finalcat["GAIAEDR3_PMRA"][hdbscan_lowq],'k.',label="Quality cut",
                 alpha=0.95)
        ax3.set_xlim(-20,15)
        plt.xlabel("Parallax (mas)")
        plt.savefig("ic2391_test_match_astrometry.png")



        # Plot a CMD of all cluster members
        plt.figure()
        # plt.plot(bp_rp,finalcat["GAIAEDR3_G_CORRECTED"],'k.',
        #          label="All stars",alpha=0.05)
        # plt.plot(bp_rp[allmemb],finalcat["GAIAEDR3_G_CORRECTED"][allmemb],'k.',
        #          label="Members",alpha=0.5)
        plt.plot(bp_rp[hdb_memb],finalcat["GAIAEDR3_G_CORRECTED"][hdb_memb],'o',
                 label="HDBScan",alpha=0.75)
        plt.plot(bp_rp[can_memb],finalcat["GAIAEDR3_G_CORRECTED"][can_memb],'o',
                 label="Cantat-Gaudin",alpha=0.75)
        plt.plot(bp_rp[ges_memb],finalcat["GAIAEDR3_G_CORRECTED"][ges_memb],'o',
                 label="GES",alpha=0.75)
        plt.plot(bp_rp[lowq_close],finalcat["GAIAEDR3_G_CORRECTED"][lowq_close],'ko',
                 label="Low Q",alpha=0.75)
        plt.xlim(-0.5,5)
        plt.ylim(25,0)
        plt.legend()
        plt.savefig("ic2391_test_match_cmd.png")
        plt.close()


if __name__=="__main__":

    # test_cantat_gaudin()

    xmatch_ic2391(to_plot=True)
