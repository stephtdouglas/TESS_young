"""
Script to retrieve TESS lightcurves using LightKurve
"""

import os
from datetime import date

import matplotlib.pyplot as plt
import numpy as np
import astropy.io.ascii as at
import astropy.io.fits as fits
from astropy.table import vstack

from lightkurve import search_lightcurve


def download_one_set(ticname,pipeline,sectors=[8,9]):
    """
    Retrieve HLSP light curves from MAST for a particular TIC ID

    Inputs:
        ticname: string, TIC id in the format "TIC ########"
        pipeline: string, HLSP source (e.g., "CDIPS","PATHOS","QLP")

    Returns:
        found: bool, whether a valid light curve was found
        multisector: bool, whether multiple light curves were found

    """

    search = search_lightcurve(ticname, author=pipeline, sector=sectors)
    if len(search)>0:
        lc = search.download_all(download_dir="/data/douglaslab/.lightkurve-cache/")
        return search.table
    else:
        return None

def download_list(ticnames,outfilename,pipelines=["CDIPS","QLP","PATHOS"],
                  sectors=[8,9,10]):
    """
    Retrieve HLSP light curves from MAST for a list of ticnames.

    """

    search_table_list = []

    for ticname in ticnames:
        print(ticname)
        for pipeline in pipelines:
            search_table = download_one_set(f"TIC {ticname}",pipeline,sectors=sectors)
            if search_table is not None:
                search_table_list.append(search_table)
                # print(search_table.dtype.names)
            else:
                print("missing:",ticname,pipeline,sectors)

    out_table = vstack(search_table_list,join_type="outer")
    at.write(out_table,outfilename,delimiter=",")

if __name__=="__main__":

    today = date.today()

    cat_files = ["Collinder_135_crossmatch_xmatch_TIC.csv",
                 "IC_2391_crossmatch_xmatch_TIC.csv",
                 "NGC_2451A_crossmatch_xmatch_TIC.csv",
                 "NGC_2547_crossmatch_xmatch_TIC.csv"
                 ]

    sector_list = [[6,7,8],[8,9,10],[7,8],[7,8,9]]

    for i,filename in enumerate(cat_files):
        print(filename)
        cat = at.read(filename,delimiter=",")
        dl_filename = filename.replace("crossmatch_xmatch_TIC",
                                       f"downloads_{str(today)}")
        print(dl_filename)
        download_list(cat["TIC"],dl_filename,sectors=sector_list[i])
