import sys, os, pathlib

import numpy as np
from numpy.random import default_rng
import matplotlib.pyplot as plt
import astropy.io.ascii as at
from astropy.table import Table

import tess_young
from tess_young.get_const import *
_DIR = pathlib.Path(tess_young.__file__).resolve().parent.parent
plt.style.use(os.path.join(_DIR,'paper.mplstyle'))

from periodmass import PeriodMassDistribution

class UpSco(PeriodMassDistribution):
    # Replacing: __init__
    # Inheriting: select_obs, plot_obs, plot_period_perc, calc_mass_percentiles

    def __init__(self,mass_limits):
        """
        Read in Angie's UpSco period mass distribution to compare with our models
        """

        # Read in the catalog
        usco_path = os.path.join(_DIR,"catalogs/AngieUpperScoMasses.txt")
        per = at.read(usco_path)
        # per.dtype
        self.cat = Table(per, masked=True, copy=False)
        if mass_limits is not None:
            mass_select = ((self.cat["AngieMass"]>=mass_limits[0]) & 
                           (self.cat["AngieMass"]<=mass_limits[1]))
            self.cat = self.cat[mass_select]

        # Assign the catalog values to attributes
        self.prot_raw = self.cat["P1"]
        self.mass_raw = self.cat["AngieMass"]
        self.mass_err_raw = np.ones_like(self.cat["AngieMass"])
        self.prot_mask = self.cat["P1"].mask
        self.mass_mask = self.cat["AngieMass"].mask

        # Assign masks for consistency, though I'm not sure they matter here
        pmask = ((self.prot_mask==False) & (self.prot_raw>0))
        pmask = pmask.filled(fill_value=False)
        mmask = (self.mass_mask==False)
        self.qmask = pmask & mmask

        # I have to modify the arrays when I run tau-squared, so need to separate the raw data from down-selected
        self.prot = self.prot_raw[self.qmask]
#         # This definitely caused the mask to change, and it added masked values...
#         self.prot.mask = self.qmask
        self.mass = self.mass_raw[self.qmask]
        self.mass_err = self.mass_err_raw[self.qmask]

        self.figsize=(9,9)

        self.param_string = f"UpSco_"

