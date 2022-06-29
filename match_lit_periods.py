"""
Script to crossmatch literature periods with Gaia/TESS
"""

import os, sys
from datetime import date

import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.cm as cm
import numpy as np
import astropy.io.ascii as at
import astropy.io.fits as fits
from astropy.coordinates import SkyCoord
from astroquery.simbad import Simbad
from astroquery.vizier import Vizier
# from astroquery.gaia import Gaia
from astropy.table import join, Table, vstack
# from astroquery.xmatch import XMatch
from astropy import units as u

from analyze_cluster_output import process_cluster, read_cluster_visual


import get_colors
norm, mapper, cmap2, colors, shapes = get_colors.get_colors()

def vizier_tic(simbad_name,gaia_dr2):
    print("Find TIC for ",simbad_name,gaia_dr2)
    result = Vizier.query_region(simbad_name,radius=2*u.arcsec,
                                 catalog="IV/38/tic")

    tic = "0"

    if len(result)==0:
        print(simbad_name,gaia_dr2,"not found in TIC")
    else:
        for row in result[0]:
            print(row["GAIA","TIC"])
            print(gaia_dr2==str(row["GAIA"]))
            if gaia_dr2==str(row["GAIA"]):
                tic = str(row["TIC"])

    return tic


def vizier_tic_skycoord(simbad_name,gaia_dr2,ra,dec):
    print("Find TIC for ",simbad_name,gaia_dr2)
    result = Vizier.query_region(SkyCoord(ra,dec,unit=u.degree),
                                 radius=2*u.arcsec,
                                 catalog="IV/38/tic")

    tic = "0"


    if len(result)==0:
        print(simbad_name,gaia_dr2,"not found in TIC")
    else:
        for row in result[0]:
            print(row["GAIA","TIC"])
            print(gaia_dr2==str(row["GAIA"]))
            if gaia_dr2==str(row["GAIA"]):
                tic = str(row["TIC"])
            else:
                print(row["Mass"])

    return tic


def xmatch_ic2391_ic2602():

    # Crossmatch Barnes+1999 periods for IC2602
    print("IC2602 Barnes (1999)")
    barnesfile = os.path.expanduser("~/Dropbox/data/catalogs/IC2602_rotation_barnes1999.csv")
    barnes = at.read(barnesfile)

    barnes["TIC"] = np.zeros(len(barnes),"U20")
    barnes["GaiaDR2"] = np.zeros(len(barnes),"U50")
    barnes["SimbadName"] = np.zeros(len(barnes),"U20")

    for i,name in enumerate(barnes["Name"]):

        if "R" in name:
            if "88A" in name:
                simbad_name = "[RSP95] 88"
            elif "24A" in name:
                simbad_name = "[RSP95] 24"
            else:
                simbad_name = name.replace("R","[RSP95] ")
        elif "W" in name:
            simbad_name = "Cl* IC 2602 "+name
        else:
            simbad_name = "IC 2602 "+name[1:]

        print(name,simbad_name)
        # result_table = Simbad.query_object(simbad_name)
        result_table = Simbad.query_objectids(simbad_name)
        # print(result_table)

        barnes["SimbadName"][i] = simbad_name

        if result_table is None:
            print(name,simbad_name,"Not Found")
            continue

        for id in result_table["ID"]:
            if "TIC" in id:
                barnes["TIC"][i] = id[4:]
            elif "DR2" in id:
                barnes["GaiaDR2"][i] = id[9:]

        if (barnes["TIC"][i] == "0") or (barnes["TIC"][i] == ""):
            tic = vizier_tic(simbad_name,barnes["GaiaDR2"][i])
            print(simbad_name,tic)
            barnes["TIC"][i] = tic

    at.write(barnes,"IC2602_rotation_barnes1999_simbad.csv",delimiter=",",
            overwrite=True)


    # Crossmatch Tschape & Rudiger 2001 periods for IC2602
    print("\n\nIC2602 Tschape & Rudiger (2001)")
    trfile = os.path.expanduser("~/Dropbox/data/catalogs/IC2602_rotation_tschape2001.csv")
    tr = at.read(trfile)

    tr["TIC"] = np.zeros(len(tr),"U20")
    tr["GaiaDR2"] = np.zeros(len(tr),"U50")
    tr["SimbadName"] = np.zeros(len(tr),"U20")

    for i,name in enumerate(tr["Name"]):

        if "R" in name:
            simbad_name = name.replace("R","[RSP95] ")
        elif "W" in name:
            simbad_name = "Cl* IC 2602 "+name
        else:
            simbad_name = "IC 2602 "+name[1:]

        print(name,simbad_name)
        # result_table = Simbad.query_object(simbad_name)
        result_table = Simbad.query_objectids(simbad_name)
        # print(result_table)
        tr["SimbadName"][i] = simbad_name

        if result_table is None:
            print(name,simbad_name,"Not Found")
            continue

        for id in result_table["ID"]:
            if "TIC" in id:
                tr["TIC"][i] = id[4:]
            elif "DR2" in id:
                tr["GaiaDR2"][i] = id[9:]

        if (tr["TIC"][i] == "0") or (tr["TIC"][i] == ""):
            tic = vizier_tic(simbad_name,tr["GaiaDR2"][i])
            print(simbad_name,tic)
            tr["TIC"][i] = tic

    at.write(tr,"IC2602_rotation_tschape2001_simbad.csv",delimiter=",",
            overwrite=True)


    # Crossmatch Prosser-Stauffer archival periods for IC2602
    print("\n\nIC2602 Prosser-Stauffer archive (1998)")
    trfile = os.path.expanduser("~/Dropbox/data/prosser_stauffer_archive/2602.spots.csv")
    tr = at.read(trfile)

    tr["TIC"] = np.zeros(len(tr),"U20")
    tr["GaiaDR2"] = np.zeros(len(tr),"U50")
    tr["SimbadName"] = np.zeros(len(tr),"U20")

    for i,name in enumerate(tr["STAR"]):

        if tr["STAR2"].mask[i]==False:
            alt_name = tr["STAR2"][i].replace("(","").replace(")","")
        else:
            alt_name = ""

        if "R" in name:
            simbad_name = name.replace("R","[RSP95] ")
        elif "R" in alt_name:
            simbad_name = alt_name.replace("R","[RSP95] ")
        elif "W" in name:
            simbad_name = "Cl* IC 2602 "+name
        else:
            # This accounts for "B" stars
            simbad_name = "IC 2602 "+name[1:]

        if simbad_name.endswith("A"):
            simbad_name = simbad_name[:-1]

        print(name,simbad_name)
        # result_table = Simbad.query_object(simbad_name)
        result_table = Simbad.query_objectids(simbad_name)
        # print(result_table)
        tr["SimbadName"][i] = simbad_name

        if result_table is None:
            print(name,simbad_name,"Not Found")
            continue

        for id in result_table["ID"]:
            if "TIC" in id:
                tr["TIC"][i] = id[4:]
            elif "DR2" in id:
                tr["GaiaDR2"][i] = id[9:]

        if (tr["TIC"][i] == "0") or (tr["TIC"][i] == ""):
            tic = vizier_tic(simbad_name,tr["GaiaDR2"][i])
            print(simbad_name,tic)
            tr["TIC"][i] = tic

    at.write(tr,"IC2602_rotation_prosserstauffer_simbad.csv",delimiter=",",
            overwrite=True)



    # Crossmatch Patten & Simon 1996 periods for IC2391
    print("\n\nIC2391 Patten & Simon (1996)")
    psfile = os.path.expanduser("~/Dropbox/data/catalogs/IC2391_rotation_patten1996.csv")
    ps = at.read(psfile)
    print(ps.dtype)

    ps["TIC"] = np.zeros(len(ps),"U20")
    ps["GaiaDR2"] = np.zeros(len(ps),"U50")
    ps["SimbadName"] = np.zeros(len(ps),"U20")

    for i,name in enumerate(ps["Name"]):

        if "60a" in name:
            simbad_name = "Cl* IC 2391 SHJM 4"
            ps["TIC"][i] = 812594503
        elif "60b" in name:
            simbad_name = "Cl* IC 2391 SHJM 5"
            ps["TIC"][i] = 94184691
        else:
            simbad_name = name.replace(" "," PSPC ")

        print(name,simbad_name)
        # result_table = Simbad.query_object(simbad_name)
        result_table = Simbad.query_objectids(simbad_name)
        # print(result_table)

        ps["SimbadName"][i] = simbad_name

        if result_table is None:
            print(name,simbad_name,"Not Found")
            continue

        for id in result_table["ID"]:
            if "TIC" in id:
                ps["TIC"][i] = id[4:]
            elif "DR2" in id:
                ps["GaiaDR2"][i] = id[9:]

        if (ps["TIC"][i] == "0") or (ps["TIC"][i] == ""):
            tic = vizier_tic(simbad_name,ps["GaiaDR2"][i])
            print(simbad_name,tic)
            ps["TIC"][i] = tic

    at.write(ps,"IC2391_rotation_patten1996_simbad.csv",delimiter=",",
            overwrite=True)

def xmatch_ngc2547():

    irwinfile = os.path.expanduser("~/Dropbox/data/catalogs/ngc2547_rotation_irwin2008b.tsv")
    ir = at.read(irwinfile,delimiter="|",data_start=3)

    ir["TIC"] = np.zeros(len(ir),"U20")
    ir["GaiaDR2"] = np.zeros(len(ir),"U50")
    ir["SimbadName"] = np.zeros(len(ir),"U30")

    for i,name in enumerate(ir["[IHA2008]"]):

        simbad_name = "[IHA2008] "+ name

        # print(name,simbad_name)
        # result_table = Simbad.query_object(simbad_name)
        result_table = Simbad.query_objectids(simbad_name)
        # print(result_table)
        ir["SimbadName"][i] = simbad_name

        if result_table is None:
            print(name,simbad_name,"Not Found")
            continue

        for id in result_table["ID"]:
            if "TIC" in id:
                ir["TIC"][i] = id[4:]
            elif "DR2" in id:
                ir["GaiaDR2"][i] = id[9:]

        if (ir["TIC"][i] == "0") or (ir["TIC"][i] == ""):
            tic = vizier_tic(simbad_name,ir["GaiaDR2"][i])
            print(simbad_name,tic)
            if tic=="0":
                tic = vizier_tic_skycoord(simbad_name,ir["GaiaDR2"][i],
                                          ir["_RAJ2000"][i],ir["_DEJ2000"][i])
                print(ir["_RAJ2000"][i],ir["_DEJ2000"][i],tic)
                print(ir["Vmag","Icmag","Mass"][i])
            ir["TIC"][i] = tic

    print(len(np.where(ir["TIC"]=="0")[0])," without TIC matches")

    ir.meta = {}
    at.write(ir,"NGC2547_rotation_irwin2008_simbad.csv",delimiter=",",
            overwrite=True)

def catalog_numbers():

    if os.path.exists("tab_all_stars.csv"):
        alldat = at.read("tab_all_stars.csv")
        check_tess = True
    else:
        alldata = None
        check_tess = False

    #### IC 2391 ############
    #########################
    cluster = "IC_2391"
    print(cluster,"\n-------")

    #### Patten+1996
    simbadfile = "IC2391_rotation_patten1996_simbad.csv"
    simbad = at.read(simbadfile, delimiter=",")

    catfile = f"{cluster}_crossmatch_xmatch_TIC.csv"
    cat = at.read(catfile,delimiter=",")

    match = join(simbad,cat,join_type="left",keys=["TIC"],
                 table_names=["lit","new"])
    # print(match.dtype)
    unames = np.unique(match["SimbadName"][match["GAIAEDR3_RA"].mask==False])
    uid = np.unique(match["GAIAEDR3_ID"][match["GAIAEDR3_RA"].mask==False])

    print(len(unames), len(uid),"out of",
          len(simbad),"from Patten & Simon in updated catalog")
    # print(match["Name"],"\n")

    # Check on the TESS data now
    if check_tess:
        match = join(simbad,alldat,join_type="left",keys=["TIC"],
                     table_names=["lit","new"])
        unames = np.unique(match["SimbadName"][(match["Q1"]<9) & (match["Q1"].mask==False)])
        print(len(unames),"out of",
              len(simbad),"from Patten & Simon have TESS data")


    #### IC 2602 ############
    #########################
    cluster = "IC_2602"
    print("\n",cluster,"\n-------")

    #### Barnes+1999
    simbadfile = "IC2602_rotation_barnes1999_simbad.csv"
    simbad = at.read(simbadfile, delimiter=",")

    catfile = f"{cluster}_crossmatch_xmatch_TIC.csv"
    cat = at.read(catfile,delimiter=",")

    match = join(simbad,cat,join_type="left",keys=["TIC"],
                 table_names=["lit","new"])
    # print(match.dtype)

    unames = np.unique(match["SimbadName"][match["GAIAEDR3_RA"].mask==False])
    uid = np.unique(match["GAIAEDR3_ID"][match["GAIAEDR3_RA"].mask==False])

    print(len(unames), len(uid),"out of",
          len(simbad),"from Barnes in updated catalog")
    # print(match["Name"],"\n")
    # uniq, ct = np.unique(match["Name"],return_counts=True)
    # unames = uniq[ct>1]
    # for name in unames:
    #     i = match["Name"]==name
    #     print(match["Name","GaiaDR2_lit","GaiaDR2_new","TIC","angDist",
    #                 "GAIAEDR3_ID","target","filter","angDist_Cantat-Gaudin","proba"][i])

    # Check on the TESS data now
    if check_tess:
        match = join(simbad,alldat,join_type="left",keys=["TIC"],
                     table_names=["lit","new"])
        unames = np.unique(match["SimbadName"][(match["Q1"]<9) & (match["Q1"].mask==False)])
        print(len(unames),"out of",
              len(simbad),"from Barnes have TESS data")

    #### Tschape & Rudiger 2001
    print("\n")
    simbadfile = "IC2602_rotation_tschape2001_simbad.csv"
    simbad = at.read(simbadfile, delimiter=",")

    match = join(simbad,cat,join_type="left",keys=["TIC"],
                 table_names=["lit","new"])
    # print(match.dtype)

    unames = np.unique(match["SimbadName"][match["GAIAEDR3_RA"].mask==False])
    uid = np.unique(match["GAIAEDR3_ID"][match["GAIAEDR3_RA"].mask==False])

    print(len(unames), len(uid), "out of",len(simbad),
          "from Tschape & Rudiger in updated catalog")
    # print(match["Name","TIC","GAIAEDR3_ID","GAIAEDR3_G"])

    # Check on the TESS data now
    if check_tess:
        match = join(simbad,alldat,join_type="left",keys=["TIC"],
                     table_names=["lit","new"])
        unames = np.unique(match["SimbadName"][(match["Q1"]<9) & (match["Q1"].mask==False)])
        print(len(unames),"out of",
              len(simbad),"from Tschape & Rudiger have TESS data")

    #### Prosser & Stauffer Archive
    print("\n")
    simbadfile = "IC2602_rotation_prosserstauffer_simbad.csv"
    simbad = at.read(simbadfile, delimiter=",")

    match = join(simbad,cat,join_type="left",keys=["TIC"],
                 table_names=["lit","new"])
    # print(match.dtype)

    unames = np.unique(match["SimbadName"][match["GAIAEDR3_RA"].mask==False])
    uid = np.unique(match["GAIAEDR3_ID"][match["GAIAEDR3_RA"].mask==False])

    print(len(unames), len(uid),"out of",len(simbad),
          "from the Prosser-Stauffer archive in updated catalog")
    # print(match["Name","TIC","GAIAEDR3_ID","GAIAEDR3_G"])

    patten96 = np.zeros(len(match),bool)
    patten96[(match["REF"]==1)] = True
    for i, row in enumerate(match):
        if "ref 1" in row["NOTES"]:
            patten96[i] = True

    # Now check for just the new periods that aren't in Barnes
    unames = np.unique(match["SimbadName"][((match["GAIAEDR3_RA"].mask==False) &
                                           (patten96))])
    print(len(unames), "out of", 9,
          "from Patten+96 in updated catalog\n")

    print(match["SimbadName","TIC","STAR"][patten96])

    # Check on the TESS data now
    if check_tess:
        match = join(simbad,alldat,join_type="left",keys=["TIC"],
                     table_names=["lit","new"])
        unames = np.unique(match["SimbadName"][match["Q1"]<9])
        print(len(unames),"out of",
              len(simbad),"from the Prosser-Stauffer archive have TESS data")

        patten96 = np.zeros(len(match),bool)
        patten96[(match["REF"]==1)] = True
        for i, row in enumerate(match):
            if "ref 1" in row["NOTES"]:
                patten96[i] = True

        unames = np.unique(match["SimbadName"][((match["Q1"]<9) &
                                                (match["Q1"].mask==False) &
                                                patten96)])
        print(len(unames), "out of",9,
              "from Patten+96 have TESS data")
        print(match["SimbadName","TIC","STAR"][patten96])

    # NGC 2547
    cluster = "NGC_2547"
    print("\n",cluster,"\n-------")
    simbadfile = "NGC2547_rotation_irwin2008_simbad.csv"
    simbad = at.read(simbadfile, delimiter=",")

    catfile = f"{cluster}_crossmatch_xmatch_TIC.csv"
    cat = at.read(catfile,delimiter=",")

    spos = SkyCoord(simbad["_RAJ2000"],simbad["_DEJ2000"],unit=u.degree)
    cpos = SkyCoord(cat["GAIAEDR3_RA"],cat["GAIAEDR3_DEC"],unit=u.degree)
    print(len(spos),len(cpos),len(cat))

    # idx,sep,_ = spos.match_to_catalog_sky(cpos)
    # good_match = np.where(sep<(5*u.arcsec))[0]
    # good_idx = idx[good_match]
    #
    # print(len(idx),len(good_match),len(good_idx))
    # # print(good_match,good_idx)
    #
    # for i in good_match:
    #     if cat["TIC"][idx[i]]!=simbad["TIC"][i]:
    #         print(cat["TIC"][idx[i]],simbad["TIC"][i])

    match = join(simbad[simbad["TIC"]!=0],cat,join_type="left",keys=["TIC"],
                 table_names=["lit","new"])
    # print(match.dtype)
    # print(np.unique(match["SimbadName"]))

    unames = np.unique(match["SimbadName"][match["GAIAEDR3_RA"].mask==False])
    uid = np.unique(match["GAIAEDR3_ID"][match["GAIAEDR3_RA"].mask==False])

    print(len(unames), len(uid),"out of",
          len(simbad),"from Irwin in updated catalog")
    # print(match["Name"],"\n")
    #
    # for row in match:
    #     print(row["GAIAEDR3_RA","_RAJ2000","TIC"])

    # Check on the TESS data now
    if check_tess:
        match = join(simbad,alldat,join_type="left",keys=["TIC"],
                     table_names=["lit","new"])
        unames = np.unique(match["SimbadName"][(match["Q1"]<9) & (match["Q1"].mask==False) & (match["GAIAEDR3_RA"].mask==False)])
        print(len(unames),"out of",
              len(simbad),"from Irwin have TESS data")


def compare_literature():

    if os.path.exists("tab_all_stars.csv"):
        alldat = at.read("tab_all_stars.csv")
        check_tess = True
    else:
        print("WARNING: Result file not found; cannot compare to literature")
        return

    plt.figure(figsize=(11,3.5))

    ############################################################################
    ############################################################################
    # IC 2391
    cluster = "IC_2391"
    print(cluster,"\n-------")
    simbadfile = "IC2391_rotation_patten1996_simbad.csv"
    simbad = at.read(simbadfile, delimiter=",")

    clean = ((alldat["Q1"]==0) | (alldat["Q1"]==1)) & (alldat["Cluster"]==cluster)
    period_col = "Period"

    match = join(simbad,alldat[clean],join_type="left",keys=["TIC"],
                 table_names=["lit","new"])

    ax1 = plt.subplot(131)
    q0 = (match["Q1"]==0) & (match["Q1"].mask==False)
    q1 = (match["Q1"]==1) & (match["Q1"].mask==False)
    ax1.plot(match[period_col][q0],match["Prot1"][q0],shapes[cluster],color=colors[cluster],
             label="Patten & Simon (1996)")
    ax1.plot(match[period_col][q1],match["Prot1"][q1],shapes[cluster],color=colors[cluster],
             mfc="none",label="TESS low-quality")
    ax1.legend(loc=2)

    # Now check for multiperiodic stars
    seconds = (match["Q2"]<=1) & (match["Q2"].mask==False)
    if np.any(seconds):
        q0 = (match["Q2"]==0) & (match["Q2"].mask==False)
        q1 = (match["Q2"]==1) & (match["Q2"].mask==False)
        ax1.plot(match[period_col][q0],match["Prot2"][q0],shapes[cluster],color=colors[cluster])
        ax1.legend(loc=2)
        ax1.plot(match[period_col][q1],match["Prot2"][q1],shapes[cluster],color=colors[cluster],
                 mfc="none")
        for i in np.where(seconds)[0]:
            ax1.plot([match[period_col][i],match[period_col][i]],
                     [match["Prot1"][i],match["Prot2"][i]],'k-')

    # Identify significant discrepancies
    dbl = match["Prot1"]*2.1
    half = match["Prot1"]/2.1
    sig_diff_i = (((match["Q1"]<=1) & (match["Q1"].mask==False))&
                  ((match[period_col]<half) | (match[period_col]>dbl)))
    print(match[sig_diff_i])
    at.write(match[sig_diff_i],f"{cluster}_significant_differences.csv",delimiter=",")


    x = np.linspace(0.1,50,20)
    ax1.plot(x,x,"-",zorder=-5,color="grey")
    ax1.plot(x,x/2,"--",zorder=-5,color="grey")
    ax1.plot(x,x*2,"--",zorder=-5,color="grey")

    ax1.set_xlim(0.1,50)
    ax1.set_ylim(0.1,50)
    ax1.set_xscale("log")
    ax1.set_yscale("log")

    ax1.set_ylabel("TESS Period (d)")
    ax1.set_xlabel("Literature Period (d)")
    ax1.set_title(cluster.replace("_"," "))

    ############################################################################
    ############################################################################
    # IC 2602
    cluster = "IC_2602"
    print("\n",cluster,"\n-------")

    # Barnes periods
    simbadfile = "IC2602_rotation_barnes1999_simbad.csv"
    simbad = at.read(simbadfile, delimiter=",")
    match = join(simbad,alldat,join_type="left",keys=["TIC"],
                 table_names=["lit","new"])

    # Patten 1996 unpublished periods.
    simbadfile2 = "IC2602_rotation_prosserstauffer_simbad.csv"
    simbad2 = at.read(simbadfile2, delimiter=",")
    match2 = join(simbad2,alldat,join_type="left",keys=["TIC"],
                 table_names=["lit","new"])
    patten96 = np.zeros(len(match2),bool)
    patten96[(match2["REF"]==1)] = True
    for i, row in enumerate(match2):
        if "ref 1" in row["NOTES"]:
            patten96[i] = True

    # TODO: separate Q=0 and Q=1 stars, at least for Barnes
    ax2 = plt.subplot(132)
    ax2.plot(match2["P(days)"][patten96],match2["Prot1"][patten96],'D',
             mec="#174f34",mfc="none",
             label="Patten+ (1996)",mew=1.5,zorder=5)
    # ax2.plot(match["Prot"],match["Prot1"],shapes[cluster],color=colors[cluster],
    #          label="Barnes+ (1999)",zorder=4)
    period_col = "Prot"
    q0 = (match["Q1"]==0) & (match["Q1"].mask==False)
    q1 = (match["Q1"]==1) & (match["Q1"].mask==False)
    ax2.plot(match[period_col][q0],match["Prot1"][q0],shapes[cluster],color=colors[cluster],
             label="Barnes+ (1999)")
    ax2.plot(match[period_col][q1],match["Prot1"][q1],shapes[cluster],color=colors[cluster],
             mfc="none",label="TESS low-quality")
    ax2.legend(loc=2)

    x = np.linspace(0.1,50,20)
    ax2.plot(x,x,"-",zorder=-5,color="grey")
    ax2.plot(x,x/2,"--",zorder=-5,color="grey")
    ax2.plot(x,x*2,"--",zorder=-5,color="grey")

    ax2.set_xlim(0.1,50)
    ax2.set_ylim(0.1,50)
    ax2.set_xscale("log")
    ax2.set_yscale("log")

    # x = np.logspace(-1,2,200)
    # for p in [1,2,3,4,5,6,7]:
    #     beat_period = 1/(1/x+1/p)
    #     ax2.plot(x,beat_period,":",zorder=-5,color="lightgrey")
    #     ax2.plot(beat_period,x,":",zorder=-5,color="lightgrey")

    # ax2.set_ylabel("TESS Period (d)")
    ax2.set_xlabel("Literature Period (d)")
    ax2.set_title(cluster.replace("_"," "))


    ############################################################################
    ############################################################################
    # NGC 2547
    cluster = "NGC_2547"
    print("\n",cluster,"\n-------")
    simbadfile = "NGC2547_rotation_irwin2008_simbad.csv"
    simbad = at.read(simbadfile, delimiter=",")

    clean = ((alldat["Q1"]==0) | (alldat["Q1"]==1)) & (alldat["Cluster"]==cluster)
    period_col = "Per"

    match = join(simbad[simbad["TIC"]!=0],alldat[clean],join_type="left",keys=["TIC"],
                 table_names=["lit","new"])

    ax3 = plt.subplot(133)
    q0 = (match["Q1"]==0) & (match["Q1"].mask==False)
    q1 = (match["Q1"]==1) & (match["Q1"].mask==False)
    ax3.plot(match[period_col][q0],match["Prot1"][q0],'v',color=colors[cluster],label="Irwin+ (2008)")
    ax3.plot(match[period_col][q1],match["Prot1"][q1],'v',color=colors[cluster],
             mfc="none",label="TESS low-quality")
    ax3.legend(loc=2)

    # Now check for multiperiodic stars
    seconds = (match["Q2"]<=1) & (match["Q2"].mask==False)
    if np.any(seconds):
        q0 = (match["Q2"]==0) & (match["Q2"].mask==False)
        q1 = (match["Q2"]==1) & (match["Q2"].mask==False)
        ax3.plot(match[period_col][q0],match["Prot2"][q0],shapes[cluster],color=colors[cluster])
        ax3.legend(loc=2)
        ax3.plot(match[period_col][q1],match["Prot2"][q1],shapes[cluster],color=colors[cluster],
                 mfc="none")
        for i in np.where(seconds)[0]:
            ax3.plot([match[period_col][i],match[period_col][i]],
                     [match["Prot1"][i],match["Prot2"][i]],'k-')
    # multi = (match["MP?"]=="y") | (match["MP?"]=="m")
    # ax3.plot(match[period_col][multi],match["Prot1"][multi],'kv',ms=10,mfc="none")
    # blend = (match["Bl?"]=="y") #| (match["Bl?"]=="m")
    # ax3.plot(match[period_col][blend],match["Prot1"][blend],shapes[cluster],
    #          color='k',ms=10,mfc="none",)

    # Identify significant discrepancies
    dbl = match["Prot1"]*2.1
    half = match["Prot1"]/2.1
    sig_diff_i = (((match["Q1"]<=1) & (match["Q1"].mask==False))&
                  ((match[period_col]<half) | (match[period_col]>dbl)))
    print(match[sig_diff_i])
    at.write(match[sig_diff_i],f"{cluster}_significant_differences.csv",delimiter=",")

    x = np.linspace(0.1,50,20)
    ax3.plot(x,x,"-",zorder=-5,color="grey")
    ax3.plot(x,x/2,"--",zorder=-5,color="grey")
    ax3.plot(x,x*2,"--",zorder=-5,color="grey")

    x = np.logspace(-1,2,200)
    for p in [1,2,3,4,5,6,7]:
        beat_period = 1/(1/x+1/p)
        ax3.plot(x,beat_period,":",zorder=-5,color="lightgrey")
        ax3.plot(beat_period,x,":",zorder=-5,color="lightgrey")

    ax3.set_xlim(0.1,50)
    ax3.set_ylim(0.1,50)
    ax3.set_xscale("log")
    ax3.set_yscale("log")

    # ax3.set_ylabel("TESS Period (d)")
    ax3.set_xlabel("Literature Period (d)")
    ax3.set_title(cluster.replace("_"," "))

    plt.savefig(f"plots/literature_comparison.png",
                bbox_inches="tight")
    # plt.show()

if __name__=="__main__":

    xmatch_ic2391_ic2602()
    xmatch_ngc2547()

    catalog_numbers()
    compare_literature()
    # plot_periodcolor_with_literature()
    plt.close("all")
