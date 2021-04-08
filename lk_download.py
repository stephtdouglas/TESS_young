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
        lc = search.download_all()
        return search.table
    else:
        return None

def download_list(ticnames,outfilename,pipelines=["QLP","PATHOS"],sectors=[8,9,10]):
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

    cat = at.read("IC2391_crossmatch_xmatch_TIC.csv",delimiter=",")
    # print(cat["TIC"])
    download_list(cat["TIC"],f"IC2391_downloads_{str(today)}.csv")

    # download_list(["45404408","93270923"],outfilename="test.csv")

    # search_table = download_one_set("TIC 93270923", pipeline="PATHOS")
