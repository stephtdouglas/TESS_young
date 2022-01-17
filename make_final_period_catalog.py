import os, sys
import glob
import itertools

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.cm as cm
import astropy.io.ascii as at
from astropy.table import join,vstack,Table
from scipy import stats
from scipy.interpolate import interp1d


def read_validation_results(cluster, date, which=None):
    # Read in my visual inspection results
    if which is None:
        vis_file = f"tables/{cluster}_{date}_results_comments.csv"
    else:
        vis_file = f"tables/{cluster}_{date}_results_comments{which}.csv"
    vis = at.read(vis_file,delimiter=",")
    good = np.where(vis["Select"].mask==False)[0]
    # print(len(good))

    # Limit the table to only the light curves I analyzed
    vis = vis[good]
    vis.rename_column("\ufefftarget_name","TIC")

    # Select the final period based on quality flags
    vis["final_period"] = np.copy(vis["sig_periods"])
    vis["final_Q"] = np.copy(vis["Q"])
    # If I flaggged the highest peak as bad, but selected another peak,
    # Select that one instead
    replace2 = (vis["Q"]==2) & ((vis["Q2"]==1) | (vis["Q2"]==0))
    replace3 = (vis["Q"]==2) & ((vis["Q3"]==1) | (vis["Q3"]==0))
    vis["final_period"][replace2] = vis["sec_periods"][replace2]
    vis["final_Q"][replace2]==vis["Q2"][replace2]
    # TODO: peak files are on other laptop, need to look there for third-
    # highest peaks
    vis["final_period"][replace3] = -99
    vis["final_Q"][replace3]==vis["Q3"][replace3]

    # TODO: identify second periods, for multiperiodic targets

    # Return the output table
    return vis

def make_final_period_catalog(cluster, date, to_plot=False):
    pass

    # Retrieve the two validation catalogs

    # Crossmatch the validation catalogs on TIC IDs

    # Create new columns to hold the final-final values for
    # period, Q, light curve info, peak power, threshold, multi/spot flags

    ##### Identify cases where my results differed between steps
        # If one was flagged as Q=2, flag both as Q=2

        # If the period is the same (and Q=1 or Q=0):

            # but it's just a different light curve, pick the one with the
            # higher periodogram peak and assign the lower Q value

            # If I assigned different Q values, assign the lower Q value

        # If the periods are not the same, print them for now
        # (So I can see what I'm working with)

    # If my results agreed, just assign the first parameters to the final cols

    ##### Crossmatch to the Gaia catalogs, save relevant columns

    # Use the _crossmatch_xmatch_TIC catalogs for their TIC IDs

    # TODO: what to do when two Gaia targets match the same TIC ID?

    # TODO: what to when two TIC IDs match the same Gaia target?

    ##### Add Phill's MINESweeper results


    ##### Check for proximity

    # Check for any nearby stars within ~30 arcsec (pixelish) -
    # likely source confusion.

    # If there is a blended star
        # If the other star has Q=2 or no measured period, do nothing

        # If the other star is >=3 mags fainter, do nothing

        # If I measured two different periods for the stars, do nothing

        # If I measured the same period for both stars
            # If one star has significantly higher periodogram power
            # (threshold?), choose that star

            # Otherwise, assign the period to the brighter star

    # Check for brighter stars within 1 arcmin - possible contamination
        # If I measured the same period as the brighter star, flag the fainter star as having a bad period

        # Otherwise, just flag for possible contamination

    ##### Add literature periods

    # TODO


    ##### Clean up and output the catalog

    # Make sure every star has a value in the "final" columns

    # Select/create columns for output
    # TIC ID
    # Gaia data: EDR3 ID, DR2 ID?, RA, Dec, photometry, RUWE
    # 2MASS ID, J, H, Ks
    # HDBScan membership info
    # Jackson/GES membership info
    # Cantat-Gaudin membership info
    # My membership flag for plotting in this paper
    # Prot1, Pw1, Q1, Sig, Prot2, Pw2, Q2, MP?, SE?
    # Bl? (this will be an automated flag, unlike my previous papers)
    # Literature periods
    # Derived properties (masses, etc, from MINESweeper)
