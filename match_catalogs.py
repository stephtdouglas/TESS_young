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

def xmatch(cluster,hdbscanfile,cgfile,to_plot=False):
    """
    inputs
    ------
    hdbscanfile: filename for the .fits HDBScan catalog from Phill Cargile

    to_plot: bool, whether to create plots to examine crossmatching

    """
    print("\n",cluster,"\n-----------")

    # Phill's catalog, which also includes all the other Gaia information
    with fits.open(hdbscanfile) as hdu:
        hdbscan = Table(hdu[1].data)
        # hdbscan_memb = hdbscan[hdbscan["MemBool"]==1]

    hdbscan.rename_column("GAIADER3_ID","GAIAEDR3_ID")
    # print(hdbscan.dtype)

    # Gaia-ESO Survey catalog (Jackson+2020)
    # This catalog includes all the clusters
    ges0 = at.read("catalogs/Jackson2020_table2_GaiaEDR3_xmatch.csv")
    ges1 = at.read("catalogs/Jackson2020_table3_GaiaDR2_data.txt",
                  header_start=1,data_start=4)
    ges = join(ges0,ges1,keys=["target","filter","cluster"])

    if cluster=="IC_2391":
        gcluster = "IC2391"
    elif cluster=="NGC_2451A":
        gcluster = "NGC2451a"
    else:
        gcluster = cluster.replace("_","")
    ges = ges[ges["cluster"]==gcluster]
    # print(ges.dtype.names)
    ges.rename_column("source_id","GAIAEDR3_ID")
    print(ges["GAIAEDR3_ID"])

    # Cantat-Gaudin et al. 2020
    # With the new Xmatch protocol, this doesn't need to start on 3
    cantat = at.read(cgfile)
    cantat.rename_column("source_id","GAIAEDR3_ID")
    cantat = cantat[cantat["Cluster"]==cluster]
    # print(cantat.dtype.names)
    print(cantat["GAIAEDR3_ID"])


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


    # NOTE: this only looks at "low quality" members, not potential members
    # of other clumps or unassigned members.
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


    if to_plot:

        # Compare G magnitudes
        plt.figure()
        plt.plot(allcat["GAIAEDR3_G_CORRECTED"][hdb_memb],
        allcat["GAIAEDR3_G_CORRECTED"][hdb_memb],'o',
                 label="HDBScan",alpha=0.75)
        plt.plot(allcat["GAIAEDR3_G_CORRECTED"][can_memb],
                allcat["phot_g_mean_mag_Cantat-Gaudin"][can_memb],'o',
                 label="Cantat-Gaudin",alpha=0.75)
        plt.plot(allcat["GAIAEDR3_G_CORRECTED"][ges_memb],
                 allcat["phot_g_mean_mag_GES"][ges_memb],'o',
                 label="GES",alpha=0.75)
        plt.plot(allcat["GAIAEDR3_G_CORRECTED"][can_memb],
                allcat["Gmag"][can_memb],'o',
                 label="Cantat-Gaudin 2",alpha=0.75)
        plt.xlim(25,0)
        plt.ylim(25,0)
        x = np.linspace(0,25,10)
        plt.plot(x,x,'k-')
        plt.xlabel("Gmag final")
        plt.ylabel("Gmag input")
        plt.legend()
        plt.savefig(f"plots/{cluster}_test_match_gmags.png")
        plt.close()

        # Compare astrometry
        plt.figure()
        plt.plot(allcat["GAIAEDR3_PARALLAX"][hdb_memb],
                 allcat["GAIADR2_PARALLAX"][hdb_memb],'o',
                 label="HDBScan",alpha=0.75)
        plt.plot(allcat["GAIAEDR3_PARALLAX"][can_memb],
                allcat["Plx"][can_memb],'o',
                 label="Cantat-Gaudin",alpha=0.75)
        plt.plot(allcat["GAIAEDR3_PARALLAX"][ges_memb],
                 allcat["PLX"][ges_memb],'o',
                 label="GES",alpha=0.75)
        # plt.xlim(25,0)
        # plt.ylim(25,0)
        x = np.linspace(0,9,10)
        plt.plot(x,x,'k-')
        plt.legend()
        plt.xlabel("EDR3 Parallax")
        plt.ylabel("Lit parallax / DR2")
        plt.savefig(f"plots/{cluster}_test_match_parallax.png")
        plt.close()

        # Check angular distance
        plt.figure()
        # plt.plot(allcat["GAIAEDR3_PARALLAX"][hdb_memb],
        #          allcat["GAIADR2_PARALLAX"][hdb_memb],'o',
        #          label="HDBScan",alpha=0.75)
        plt.plot(allcat["GAIAEDR3_PARALLAX"][can_memb],
                 allcat["angDist_Cantat-Gaudin"][can_memb],'o',
                 label="Cantat-Gaudin",alpha=0.75,color="C1")
        plt.plot(allcat["GAIAEDR3_PARALLAX"][ges_memb],
                 allcat["angDist_GES"][ges_memb],'o',
                 label="GES",alpha=0.75,color="C2")
        plt.xlim(-3,10)
        # plt.ylim(25,0)
        plt.legend(loc=3)
        plt.xlabel("EDR3 Parallax")
        plt.ylabel("Xmatch angDist")
        plt.savefig(f"plots/{cluster}_test_match_angdist.png")
        plt.close()

        plt.figure()
        # plt.plot(allcat["GAIAEDR3_PARALLAX"][hdb_memb],
        #          allcat["GAIADR2_PARALLAX"][hdb_memb],'o',
        #          label="HDBScan",alpha=0.75)
        plt.plot(allcat["Gmag"][can_memb],
                 allcat["angDist_Cantat-Gaudin"][can_memb],'o',
                 label="Cantat-Gaudin",alpha=0.75,color="C1")
        plt.plot(allcat["phot_g_mean_mag_GES"][ges_memb],
                 allcat["angDist_GES"][ges_memb],'o',
                 label="GES",alpha=0.75,color="C2")
        plt.xlim(25,0)
        # plt.ylim(25,0)
        plt.legend(loc=3)
        plt.xlabel("Gmag")
        plt.ylabel("Xmatch angDist")
        plt.savefig(f"plots/{cluster}_test_match_angdistG.png")
        plt.close()

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

    subcat[any_memb].write(f"tables/{cluster}_crossmatch.fits",overwrite=True)
    subcat[any_memb].write(f"tables/{cluster}_crossmatch.csv",overwrite=True)
    finalcat[any_memb].write(f"tables/{cluster}_crossmatch_allcolumns.fits",overwrite=True)


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
        plt.savefig(f"plots/{cluster}_test_match_astrometry_lowq.png")

        # Plot astrometry for confirmed members and for "low quality" stars that
        # were left out of the HDBScan run
        plt.figure(figsize=(10,10))
        ax1 = plt.subplot(221)
        for i,memb in enumerate(lists[:3]):
            plt.plot(finalcat["GAIAEDR3_PMDEC"][memb],finalcat["GAIAEDR3_PARALLAX"][memb],'.',
                    label=names[i],alpha=0.75)
        plt.plot(finalcat["GAIAEDR3_PMDEC"],
                 finalcat["GAIAEDR3_PARALLAX"],'k.',label="All",
                 alpha=0.1,zorder=-5)
        plt.ylabel("Parallax (mas)")
        ax1.set_ylim(-20,15)
        plt.legend(loc=3)

        ax2 = plt.subplot(223,sharex=ax1)
        for i,memb in enumerate(lists[:3]):
            plt.plot(finalcat["GAIAEDR3_PMDEC"][memb],finalcat["GAIAEDR3_PMRA"][memb],'.',
                    label=names[i],alpha=0.75)
        plt.plot(finalcat["GAIAEDR3_PMDEC"],
                 finalcat["GAIAEDR3_PMRA"],'k.',label="All",
                 alpha=0.1,zorder=-5)
        plt.ylabel("PMRA (mas/yr)")
        plt.xlabel("PMDEC (mas/yr)")
        ax1.set_xlim(-50,50)
        ax2.set_ylim(-40,25)

        ax3 = plt.subplot(224,sharey=ax2)
        for i,memb in enumerate(lists[:3]):
            plt.plot(finalcat["GAIAEDR3_PARALLAX"][memb],
                    finalcat["GAIAEDR3_PMRA"][memb],'.',
                    label=names[i],alpha=0.75)
        plt.plot(finalcat["GAIAEDR3_PARALLAX"],
                 finalcat["GAIAEDR3_PMRA"],'k.',label="All",
                 alpha=0.1,zorder=-5)
        ax3.set_xlim(-20,15)
        plt.xlabel("Parallax (mas)")
        plt.savefig(f"plots/{cluster}_test_match_astrometry.png")



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
        plt.savefig(f"plots/{cluster}_test_match_cmd.png")
        plt.close()


if __name__=="__main__":

    hdbscanfile = os.path.expanduser("~/Dropbox/EDR3/scats/IC_2391.fits")
    #TODO: switch to the new Xmatch files
    cgfile = "catalogs/cantat-gaudin2020_ic2391_10deg_GaiaEDR3_xmatch_new.csv"
    xmatch("IC_2391",hdbscanfile,cgfile,to_plot=True)

    hdbscanfile = os.path.expanduser("~/Dropbox/EDR3/scats/NGC_2451A.fits")
    #TODO: switch to the new Xmatch files
    cgfile = "catalogs/cantat-gaudin2020_ngc2451A_10deg_GaiaEDR3_xmatch.csv"
    xmatch("NGC_2451A",hdbscanfile,cgfile,to_plot=True)

    hdbscanfile = os.path.expanduser("~/Dropbox/EDR3/scats/NGC_2547.fits")
    #TODO: switch to the new Xmatch files
    cgfile = "catalogs/cantat-gaudin2020_ngc2547_10deg_GaiaEDR3_xmatch.csv"
    xmatch("NGC_2547",hdbscanfile,cgfile,to_plot=True)

    # I'm missing the IC 2602 catalog from Phill
    # and Collinder 135 is missing from Jackson+2020
