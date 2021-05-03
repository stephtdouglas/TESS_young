"""
Script to crossmatch literature periods with Gaia/TESS
"""

import os, sys
from datetime import date

import matplotlib.pyplot as plt
import numpy as np
import astropy.io.ascii as at
import astropy.io.fits as fits
from astroquery.simbad import Simbad
from astropy.table import join, Table
# from astroquery.xmatch import XMatch
from astropy import units as u


if __name__=="__main__":

    # Crossmatch Barnes+1999 periods for IC2602
    print("IC2602 Barnes (1999)")
    barnesfile = os.path.expanduser("~/Dropbox/data/catalogs/IC2602_rotation_barnes1999.csv")
    barnes = at.read(barnesfile)

    barnes["TICID"] = np.zeros(len(barnes),"U20")
    barnes["GaiaDR2"] = np.zeros(len(barnes),"U50")

    for i,name in enumerate(barnes["Name"]):

        if "R" in name:
            if "88A" in name:
                simbad_name = "[RSP95] 88"
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

        if result_table is None:
            print(name,simbad_name,"Not Found")
            continue

        for id in result_table["ID"]:
            if "TIC" in id:
                barnes["TICID"][i] = id
            elif "DR2" in id:
                barnes["GaiaDR2"][i] = id

    at.write(barnes,"IC2602_rotation_barnes1999_simbad.csv",delimiter=",")


    # Crossmatch Tschape & Rudiger 2001 periods for IC2602
    print("\n\nIC2602 Tschape & Rudiger (2001)")
    trfile = os.path.expanduser("~/Dropbox/data/catalogs/IC2602_rotation_tschape2001.csv")
    tr = at.read(trfile)

    tr["TICID"] = np.zeros(len(tr),"U20")
    tr["GaiaDR2"] = np.zeros(len(tr),"U50")

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

        if result_table is None:
            print(name,simbad_name,"Not Found")
            continue

        for id in result_table["ID"]:
            if "TIC" in id:
                tr["TICID"][i] = id
            elif "DR2" in id:
                tr["GaiaDR2"][i] = id

    at.write(tr,"IC2602_rotation_tschape2001_simbad.csv",delimiter=",")


    # Crossmatch Patten & Simon 1996 periods for IC2391
    print("\n\nIC2391 Patten & Simon (1996)")
    psfile = os.path.expanduser("~/Dropbox/data/catalogs/IC2391_rotation_patten1996.csv")
    ps = at.read(psfile)
    print(ps.dtype)

    ps["TICID"] = np.zeros(len(ps),"U20")
    ps["GaiaDR2"] = np.zeros(len(ps),"U50")

    for i,name in enumerate(ps["Name"]):

        if "60a" in name:
            simbad_name = "Cl* IC 2391 SHJM 4"
        elif "60b" in name:
            simbad_name = "Cl* IC 2391 SHJM 5"
        else:
            simbad_name = name.replace(" "," PSPC ")

        print(name,simbad_name)
        # result_table = Simbad.query_object(simbad_name)
        result_table = Simbad.query_objectids(simbad_name)
        # print(result_table)

        if result_table is None:
            print(name,simbad_name,"Not Found")
            continue

        for id in result_table["ID"]:
            if "TIC" in id:
                ps["TICID"][i] = id
            elif "DR2" in id:
                ps["GaiaDR2"][i] = id

    at.write(ps,"IC2391_rotation_patten1996_simbad.csv",delimiter=",")
