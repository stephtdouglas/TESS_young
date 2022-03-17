import os, sys
import glob
import itertools

import numpy as np
from numpy.random import default_rng
# import matplotlib.pyplot as plt
# import matplotlib as mpl
# import matplotlib.cm as cm
import astropy.io.ascii as at
import astropy.units as u
from astropy.table import join,vstack,Table
from scipy import stats
from scipy.interpolate import interp1d
from analyze_cluster_output import read_cluster_visual

from analyze_cluster_output import read_cluster_visual

def calc_percentiles(cdat,color_col,period_col,color_name="V-K",e_color_col=None):

    # For 5-30 Myr stars, Pecaut & Mamajek (2013) place G0-G9 stars
    # between 1.35 < V-Ks < 2.10
    if color_name=="V-K":
        color_min = 1.35
        color_max = 2.10
    elif color_name=="Mass":
        color_min = 0.9
        color_max = 1.1
        if e_color_col is None:
            e_color = 0.05 #making this up
    elif color_name=="BP-RP":
        color_min = 0.78
        color_max = 0.90
        if e_color_col is None:
            e_color = 0.01 #making this up

    ##############################################################
    # OK. So let's Monte Carlo this.
    ##############################################################

    # Select all stars with color/mass and errors, and periods
    if e_color_col is not None:
        benchmarks = ((cdat[color_col].mask==False) &
                      (cdat[e_color_col].mask==False) &
                      (cdat[period_col].mask==False))
    else:
        benchmarks = ((cdat[color_col].mask==False) &
                      (cdat[period_col].mask==False))
    raw_solar = benchmarks & (cdat[color_col]<=color_max) & (cdat[color_col]>=color_min)
    nb = len(np.where(benchmarks)[0])


    # # perc = np.percentile(cdat[period_col][solar],[5,25,50,75,95])
    # perc = np.percentile(cdat[period_col][solar],[25,50,90])


    # Number of tests
    ntests = 1000

    # Randomly re-draw colors from within each star’s uncertainty distribution
    rng = default_rng(42)
    if e_color_col is not None:
        new_colors = rng.normal(loc=cdat[color_col][benchmarks],
                                scale=cdat[e_color_col][benchmarks],
                                size=(ntests,nb))
    else:
        new_colors = rng.normal(loc=cdat[color_col][benchmarks],
                                scale=e_color,size=(ntests,nb))

    # print(ntests,nb)
    # print(len(new_colors[0]))

    # Select all stars with re-drawn colors consistent with solar-mass (or whatever) colors
    solar = (new_colors>=color_min) & (new_colors<=color_max)
    # print(np.shape(solar))

    # Recompute 25th, 50th, and 90th percentiles
    p10,p50,p75 = np.zeros(ntests),np.zeros(ntests),np.zeros(ntests)
    for i in range(ntests):
        if len(np.where(solar[i])[0])>0:
            p10[i],p50[i],p75[i] = np.percentile(cdat[period_col][benchmarks][solar[i]],
                                                [10,50,75])


    # Compute the median value for each percentile
    p10_med = np.median(p10)
    p50_med = np.median(p50)
    p75_med = np.median(p75)

    # Compute the error on each percentile using the absolute deviation of the individual estimates around the median
    p10_mad = stats.median_abs_deviation(p10)
    p50_mad = stats.median_abs_deviation(p50)
    p75_mad = stats.median_abs_deviation(p75)

    p10_std = np.std(p10)
    p50_std = np.std(p50)
    p75_std = np.std(p75)

    print("Prot")
    print(f"10% {p10_med:.3f} +/- {p10_mad:.3f} or {p10_std:.3f}")
    print(f"50% {p50_med:.3f} +/- {p50_mad:.3f} or {p50_std:.3f}")
    print(f"75% {p75_med:.3f} +/- {p75_mad:.3f} or {p75_std:.3f}")

    return [p10_med,p50_med,p75_med], cdat[period_col][raw_solar]

def calc_percentiles_omega(cdat,color_col,period_col,color_name="V-K",e_color_col=None):
    """ Gallet & Bouvier (2013, 2015) use angular velocity Omega = 2pi/Prot
    """

    # For 5-30 Myr stars, Pecaut & Mamajek (2013) place G0-G9 stars
    # between 1.35 < V-Ks < 2.10
    if color_name=="V-K":
        color_min = 1.35
        color_max = 2.10
    elif color_name=="Mass":
        color_min = 0.9
        color_max = 1.1
        if e_color_col is None:
            e_color = 0.05 #making this up

    ##############################################################
    # OK. So let's Monte Carlo this.
    ##############################################################

    # Select all stars with color/mass and errors, and periods
    if e_color_col is not None:
        benchmarks = ((cdat[color_col].mask==False) &
                      (cdat[e_color_col].mask==False) &
                      (cdat[period_col].mask==False))
    else:
        benchmarks = ((cdat[color_col].mask==False) &
                      (cdat[period_col].mask==False))
    nb = len(np.where(benchmarks)[0])


    # # perc = np.percentile(cdat[period_col][solar],[5,25,50,75,95])
    # perc = np.percentile(cdat[period_col][solar],[25,50,90])


    # Number of tests
    ntests = 1000

    # Randomly re-draw colors from within each star’s uncertainty distribution
    rng = default_rng(42)
    if e_color_col is not None:
        new_colors = rng.normal(loc=cdat[color_col][benchmarks],
                                scale=cdat[e_color_col][benchmarks],
                                size=(ntests,nb))
    else:
        new_colors = rng.normal(loc=cdat[color_col][benchmarks],
                                scale=e_color,size=(ntests,nb))

    print(ntests,nb)
    print(len(new_colors[0]))

    # Select all stars with re-drawn colors consistent with solar-mass (or whatever) colors
    solar = (new_colors>=color_min) & (new_colors<=color_max)
    print(np.shape(solar))

    # Recompute 25th, 50th, and 90th percentiles
    p25,p50,p90 = np.zeros(ntests),np.zeros(ntests),np.zeros(ntests)
    omega_sun = 2.87e-6 / u.s
    omega = np.pi * 2 / (cdat[period_col]*u.d)
    omega = (omega / omega_sun).decompose().value
    for i in range(ntests):
        if len(np.where(solar[i])[0])>0:
            p25[i],p50[i],p90[i] = np.percentile(omega[benchmarks][solar[i]],
                                                [25,50,90])

    # Compute the median value for each percentile
    p25_med = np.median(p25)
    p50_med = np.median(p50)
    p90_med = np.median(p90)

    # Compute the error on each percentile using the absolute deviation of the individual estimates around the median
    p25_mad = stats.median_abs_deviation(p25)
    p50_mad = stats.median_abs_deviation(p50)
    p90_mad = stats.median_abs_deviation(p90)

    p25_std = np.std(p25)
    p50_std = np.std(p50)
    p90_std = np.std(p90)


    print("Omega")
    print(f"25% {p25_med:.3f} +/- {p25_mad:.3f} or {p25_std:.3f}")
    print(f"50% {p50_med:.3f} +/- {p50_mad:.3f} or {p50_std:.3f}")
    print(f"90% {p90_med:.3f} +/- {p90_mad:.3f} or {p90_std:.3f}")



    print("Converted back to Prot")
    # omega = np.pi * 2 / (cdat[period_col]*u.d)
    # omega = (omega / omega_sun).decompose().value
    print(f"25% {(np.pi*2)/(p25_med*omega_sun.to(1/u.d)):.3f}")
    print(f"50% {(np.pi*2)/(p50_med*omega_sun.to(1/u.d)):.3f}")
    print(f"90% {(np.pi*2)/(p90_med*omega_sun.to(1/u.d)):.3f}")

    return [p25_med,p50_med,p90_med]

def usco_init():
    ## Upper Sco ##############################################
    print("\nUSco")

    usco_file = os.path.expanduser("~/Dropbox/data/catalogs/usco_rhooph_rotation_rebull2018.csv")
    usco = at.read(usco_file,delimiter="|",data_start=3)
    # print(usco.dtype)
    # print(usco[0])

    perc, raw_solar_periods = calc_percentiles(usco,"(V-Ks)0","Per1",color_name="V-K",e_color_col="E(V-Ks)")
    return perc, raw_solar_periods

def young_stars_init():


    ## ONC ##############################################
    print("\nONC")

    cat_file = os.path.expanduser("~/Dropbox/data/catalogs/ONC_rotation_herbst2002_compilation.tsv")
    cat = at.read(cat_file,delimiter="|",data_start=3)
    cat = Table(cat, masked=True, copy=False)

    perc1,raw_solar_periods1 = calc_percentiles(cat,"Mass","Per",color_name="Mass")
    # perc = calc_percentiles_omega(cat,"Mass","Per",color_name="Mass")


    ## Sco Cen ##############################################
    print("\nSco-Cen")

    cat_file = os.path.expanduser("~/Dropbox/data/catalogs/ScoCen_rotation_Mellon2017.tsv")
    cat = at.read(cat_file,delimiter="|",data_start=3)
    cat = Table(cat, masked=True, copy=False)

    perc11,raw_solar_periods11 = calc_percentiles(cat,"Mass","Per",color_name="Mass")
    # perc = calc_percentiles_omega(cat,"Mass","Per",color_name="Mass")



    ## Upper Sco ##############################################
    print("\nUSco")

    usco_file = os.path.expanduser("~/Dropbox/data/catalogs/usco_rhooph_rotation_rebull2018.csv")
    usco = at.read(usco_file,delimiter="|",data_start=3)
    # print(usco.dtype)
    # print(usco[0])


    perc8,raw_solar_periods8 = calc_percentiles(usco,"(V-Ks)0","Per1",color_name="V-K",e_color_col="E(V-Ks)")
    # perc = calc_percentiles_omega(usco,"(V-Ks)0","Per1",color_name="V-K",e_color_col="E(V-Ks)")

    ## h Per ##############################################
    print("\nhPer")

    cat_file = os.path.expanduser("~/Dropbox/data/catalogs/hper_rotation_moraux2013.tsv")
    cat = at.read(cat_file,delimiter="|",data_start=3)
    cat = Table(cat, masked=True, copy=False)
    # print(usco.dtype)
    # print(usco[0])

    perc13,raw_solar_periods13 = calc_percentiles(cat,"Mass","Per",color_name="Mass")
    # perc = calc_percentiles_omega(cat,"Mass","Per",color_name="Mass")

    return [1,8,11,13], [perc1,perc8,perc11,perc13], [raw_solar_periods1,raw_solar_periods8,raw_solar_periods11,raw_solar_periods13]

def zams_percentiles():


    # clusters = ["IC_2391","Collinder_135","NGC_2451A","NGC_2547","IC_2602"]
    # dates = ["2021-06-22","2021-06-18","2021-06-21","2021-06-21","2021-07-02"]

    clusters = ["IC_2391","Collinder_135","NGC_2451A","NGC_2547","IC_2602"]
    all_cat0 = []
    indiv_perc = []
    indiv_solar = []

    dat = at.read("tab_all_stars.csv")
    dat = Table(dat, masked=True, copy=False)
    dat = dat[dat["Q1"]==0]
    perc, solar = calc_percentiles(dat, "Mass","Prot1",color_name="Mass")
    for i in range(5):
        print("\n",clusters[i])
        perc,raw_solar = calc_percentiles(dat[dat["Cluster"]==clusters[i]],
                                          "Mass","Prot1",color_name="Mass")
        indiv_perc.append(perc)
        indiv_solar.append(raw_solar)

    return perc, solar, indiv_perc, indiv_solar

if __name__=="__main__":

    # young_stars_init()
    zams_percentiles()
