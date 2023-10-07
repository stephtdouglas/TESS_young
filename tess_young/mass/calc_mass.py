"""
De-redden Gaia colors and calculate masses from an isochrone model

Stephanie T. Douglas
"""

import pathlib, os
import numpy as np
from astropy.table import Table, vstack, join
import astropy.io.ascii as at
from scipy.interpolate import interp1d

import tess_young
from tess_young.get_const import *
from tess_young.mass.gaia_reddening import calc_BP_RP0
_DIR = pathlib.Path(tess_young.__file__).resolve().parent.parent

def get_parsec():
    parsec_file = os.path.join(MODEL_DIR,"PARSEC_ZAMS.dat")
    parsec = at.read(parsec_file,header_start=13)

    log_ages = np.unique(parsec["logAge"])

    # Create an interpolator for each model between 25-55 Myr
    # (The file contains 20 and 60 Myr as well)
    isochrones = []
    iso_labels = []
    for age in log_ages[1:-1]:
        loc = (parsec["logAge"]==age) & (parsec["Mass"]<2)
        mod_bprp = parsec["G_BPmag"][loc]-parsec["G_RPmag"][loc]
        mass = parsec["Mass"][loc]
        mod_label=f"{10**(age-6):.0f} Myr"

        isochrones.append(interp1d(mod_bprp,mass,bounds_error=False))
        iso_labels.append(mod_label)

    return isochrones, iso_labels

def calc_mass(bp_rp0, isochrones=None, iso_labels=None):
    if isochrones is None:
        isochrones, iso_labels = get_parsec()

    # Calculate one mass per input isochrone
    mass = np.zeros((len(iso_labels),len(bp_rp0)))

    niso = len(iso_labels)

    for i in range(niso):
        mass[i] = isochrones[i](bp_rp0)

    mass_40 = mass[3]
    mass_err = (np.max(mass,axis=0) - np.min(mass,axis=0))/2

    return mass_40, mass_err


if __name__=="__main__":

    dat0 = at.read("tab_all_stars.csv")

    red_files = []
    for i, cluster in enumerate(clusters):
        red_fname = os.path.join(_DIR,f"{cluster}_REDS.csv")
        red_files.append(at.read(red_fname))
    red = vstack(red_files)

    dat = join(dat0,red,keys=["GAIAEDR3_ID"])
    dat.sort("E_BV")

    bp_rp = dat["GAIAEDR3_BP"]-dat["GAIAEDR3_RP"]
    bp_rp0 = calc_BP_RP0(bp_rp,dat["E_BV"])

    calc_mass(bp_rp0)