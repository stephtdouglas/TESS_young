"""
Script to retrieve TESS lightcurves using LightKurve
"""

import os

import matplotlib.pyplot as plt
import numpy as np
import astropy.io.ascii as at
import astropy.io.fits as fits

from lightkurve import search_lightcurve


def download_lightcurves(ticname,pipeline,sectors=[8,9]):
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

if __name__=="__main__":
    download_lightcurves("TIC 93270923", pipeline="PATHOS")
