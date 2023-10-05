"""
Script to crossmatch literature periods with Gaia/TESS
"""

import os, sys
import numpy as np
import astropy.io.ascii as at
from astropy.table import join, Table, vstack

if __name__=="__main__":

    filenames = ["IC2391_rotation_patten1996_simbad.csv",
                 "IC2602_rotation_tschape2001_simbad.csv",
                 "IC2602_rotation_barnes1999_simbad.csv",
                 "NGC2547_rotation_irwin2008_simbad.csv",
                 "IC2391_rotation_messina2011_simbad.csv",
                 "IC2602_rotation_prosserstauffer_simbad.csv"]

    tables = []

    # Read in a given file
    # save the original name, period, and TIC, GaiaDR2, SimbadName
    # change all columns to matching colnames, and add cluster and source columns
    # Columns: Name, Lit Period,
    cluster = "IC 2391"
    filename = filenames[0]
    source = "patten1996"

    tab = at.read(filename,delimiter=",")
    tab.rename_column("Period","LitPeriod")
    tab["Cluster"] = np.empty(len(tab),"U8")
    tab["Cluster"][:] = cluster
    tab["Source"] = np.empty(len(tab),"U22")
    tab["Source"][:] = source
    tab = tab["Name","Cluster","Source","LitPeriod","TIC","GaiaDR2","SimbadName"]
    tables.append(tab)

    ###
    cluster = "IC 2391"
    filename = filenames[4]
    source = "messina2011"

    tab = at.read(filename,delimiter=",")
    tab.rename_column("Per","LitPeriod")
    tab["Cluster"] = np.empty(len(tab),"U8")
    tab["Cluster"][:] = cluster
    tab["Source"] = np.empty(len(tab),"U22")
    tab["Source"][:] = source
    tab = tab["Name","Cluster","Source","LitPeriod","TIC","GaiaDR2","SimbadName"]
    tables.append(tab)

    # ###
    # cluster = "IC 2602"
    # filename = filenames[1]
    # source = "tschape2001"

    # tab = at.read(filename,delimiter=",")
    # tab.rename_column("Prot","LitPeriod")
    # tab["Cluster"] = np.empty(len(tab),"U8")
    # tab["Cluster"][:] = cluster
    # tab["Source"] = np.empty(len(tab),"U22")
    # tab["Source"][:] = source
    # tab = tab["Name","Cluster","Source","LitPeriod","TIC","GaiaDR2","SimbadName"]
    # tables.append(tab)

    ###
    cluster = "IC 2602"
    filename = filenames[5]
    source = "patten1996-1"

    tab = at.read(filename,delimiter=",")
    tab.rename_column("P(days)","LitPeriod")
    tab.rename_column("STAR","Name")
    tab["Cluster"] = np.empty(len(tab),"U8")
    tab["Cluster"][:] = cluster
    tab["Source"] = np.empty(len(tab),"U22")
    tab["Source"][:] = source
    tab = tab["Name","Cluster","Source","LitPeriod","TIC","GaiaDR2","SimbadName"]
    tables.append(tab)


    ###
    cluster = "IC 2602"
    filename = filenames[2]
    source = "barnes1996"

    tab = at.read(filename,delimiter=",")
    tab.rename_column("Prot","LitPeriod")
    tab["Cluster"] = np.empty(len(tab),"U8")
    tab["Cluster"][:] = cluster
    tab["Source"] = np.empty(len(tab),"U22")
    tab["Source"][:] = source
    tab = tab["Name","Cluster","Source","LitPeriod","TIC","GaiaDR2","SimbadName"]
    tables.append(tab)


    ###
    cluster = "NGC 2547"
    filename = filenames[3]
    source = "irwin2008"

    tab = at.read(filename,delimiter=",")
    tab.rename_column("[IHA2008]","Name")
    tab.rename_column("Per","LitPeriod")
    tab["Cluster"] = np.empty(len(tab),"U8")
    tab["Cluster"][:] = cluster
    tab["Source"] = np.empty(len(tab),"U22")
    tab["Source"][:] = source
    tab = tab["Name","Cluster","Source","LitPeriod","TIC","GaiaDR2","SimbadName"]
    tables.append(tab)


    # stack all the new tables together
    all_lit = vstack(tables)
    all_lit.sort("TIC")

    # print out in latex format
    at.write(all_lit,"tab_all_lit_periods.csv",delimiter=",",overwrite=True)

    for i in range(len(all_lit)):
        cite_src = "\\citet{"+all_lit["Source"][i]+"}"
        all_lit["Source"][i] = cite_src

    orig_cols = [column for column in all_lit.columns]
    print(orig_cols)

    for colname in orig_cols:
        new_name = "\\colhead{"+colname+"}"
        all_lit.rename_column(colname,new_name)

    at.write(all_lit,"tab_all_lit_periods.tex",Writer=at.Latex,
             latexdict={'tabletype':'deluxetable','tablealign':'ht',
                        'caption':'Literature periods for cluster stars, crossmatched to current identifiers \\label{tab:lit_periods}',
                        'col_align':'lllDlll'},overwrite=True)
