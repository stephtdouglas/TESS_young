"""
Script to retrieve TESS lightcurves using LightKurve
"""

import os
from datetime import date
from requests import HTTPError

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

    try:
        search = search_lightcurve(ticname, author=pipeline, sector=sectors)
    except HTTPError:
        print(ticname, pipeline, sectors, "MAST server error")
        return None
        
    if len(search)>0:
        lc = search.download_all(download_dir="/data/douglaslab/.lightkurve-cache/")
        return search.table
    else:
        return None

def download_list(ticnames,outfilename,pipelines=["CDIPS","QLP","PATHOS"],
                  sectors=[8,9,10],start_i=0,start_ct=0):
    """
    Retrieve HLSP light curves from MAST for a list of ticnames.

    """

    search_table_list = []
    ct = start_ct

    for i,ticname in enumerate(ticnames):
        if i < start_i:
            continue
        print(i,ticname,ct)

#        existing_out = glob.glob(f"tables/{outfilename}*")
#        last_out = 
        
#        if (i<(ct+1)*50):
#            last_out = outfilename+f".{ct+1}"
#            if os.path.exists(last_out):
#                print("already done",ct)
#                print(last_out)
#                continue
        
        for pipeline in pipelines:
            search_table = download_one_set(f"TIC {ticname}",pipeline,sectors=sectors)
            if search_table is not None:
                search_table_list.append(search_table)
                # print(search_table.dtype.names)
            #else:
            #    print("missing:",ticname,pipeline,sectors)

        if ((i % 50)==0) and (len(search_table_list)>1):
            ct += 1
            out_table = vstack(search_table_list,join_type="outer")
            print(outfilename+f".{ct}")
            at.write(out_table,outfilename+f".{ct}",delimiter=",")

    out_table = vstack(search_table_list,join_type="outer")
    at.write(out_table,outfilename,delimiter=",")

if __name__=="__main__":

    today = date.today()

    cat_files = ["Collinder_135_crossmatch_xmatch_TIC.csv",
                 "IC_2391_crossmatch_xmatch_TIC.csv",
                 "NGC_2451A_crossmatch_xmatch_TIC.csv",
                 "NGC_2547_crossmatch_xmatch_TIC.csv",
                 "IC_2602_crossmatch_xmatch_TIC.csv"
                ]

    sector_list = [[6,7,8],[8,9,10],[7,8],[7,8,9],[9,10,11]]

    for i,filename in enumerate(cat_files):
        print(i,filename)
        if i<=3:
            continue
        print(filename)
        cat = at.read(filename,delimiter=",")
        dl_filename0 = filename.replace("crossmatch_xmatch_TIC",
                                       f"downloads_{str(today)}")
        dl_filename = os.path.join("tables/",dl_filename0)
        print(dl_filename)
        download_list(cat["TIC"],dl_filename,sectors=sector_list[i])
        
