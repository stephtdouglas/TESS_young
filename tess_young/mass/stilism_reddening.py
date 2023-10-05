"""
Compute reddenings towards stars using Gaia positions and parallaxes, 
and the STILISM extinction maps.

Created 2021
Modified 2023

@author: rmk2432, stephtdouglas
"""
import pathlib, os
from scipy.interpolate import RegularGridInterpolator
import numpy as np
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.table import Table
import astropy.io.ascii as at
import pandas as pd

import tess_young
from tess_young.get_const import *
_DIR = pathlib.Path(tess_young.__file__).resolve().parent.parent


if __name__=="__main__":

    # Collect Ronan's reformated STILISM maps
    EEH = np.load(os.path.join(MODEL_DIR,'STILISM/EBV_SPHERICAL_ERRU.npy'))
    EEL = np.load(os.path.join(MODEL_DIR,'STILISM/EBV_SPHERICAL_ERRL.npy'))
    E   = np.load(os.path.join(MODEL_DIR,'STILISM/EBV_SPHERICAL.npy'))
    larr = np.load(os.path.join(MODEL_DIR,'STILISM/EBV_L.npy'))
    barr = np.load(os.path.join(MODEL_DIR,'STILISM/EBV_B.npy'))
    darr = np.load(os.path.join(MODEL_DIR,'STILISM/EBV_D.npy'))

    HBINS = np.arange(0,5.00001,0.005)
    HBC = (HBINS[1:]+HBINS[:-1])/2.


    interp_red    = RegularGridInterpolator((barr,larr,darr), E, bounds_error=False)
    interp_rederh = RegularGridInterpolator((barr,larr,darr), EEH, bounds_error=False)
    interp_rederl = RegularGridInterpolator((barr,larr,darr), EEL, bounds_error=False)


    ###load files for ra, dec, plx#########
    period_cats = []
    out_cats = []
    for i in range(5):
        print("\n\n",clusters[i])
        cluster = clusters[i]
        tic_file = os.path.join(_DIR,f"{cluster}_crossmatch_xmatch_TIC.csv")

        STARS  = pd.read_csv(tic_file)
        GRA, GDC = STARS['GAIAEDR3_RA'].to_numpy(),STARS['GAIAEDR3_DEC'].to_numpy()
        GPX = STARS['GAIAEDR3_PARALLAX_CORRECTED'].to_numpy()
        GD = 1000./GPX
        # Calculate sky positions
        d = SkyCoord(ra=GRA*u.degree, dec=GDC*u.degree, frame='icrs')
        gd = d.transform_to('galactic')
        Gl, Gb = gd.l, gd.b

        # Calculate reddenings
        reds  = interp_red((Gb, Gl, GD))
        redeh = interp_rederh((Gb, Gl, GD))
        redel = interp_rederl((Gb, Gl, GD))


        # Write out results
        dat = Table.from_pandas(STARS)

        dat["E_BV"] = reds
        dat["eup_E_BV"] = redeh
        dat["edn_E_BV"] = redel
            
        red_table = dat["GAIAEDR3_ID","E_BV","eup_E_BV","edn_E_BV"]

        outfile = os.path.join(_DIR,f"{cluster}_REDS.csv")
        at.write(red_table,outfile,overwrite=True)
