"""
Script to measure periods from all downloaded lightkurve files
"""
import sys, os
from datetime import date
import logging
import pickle

logging.basicConfig(level=logging.WARNING)

import matplotlib
matplotlib.use("agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import astropy.io.ascii as at
import astropy.io.fits as fits
from astropy.table import join, Table
from astropy import units as u
from scipy.signal import argrelextrema
import lightkurve as lk

import k2spin
from k2spin import prot
from k2spin import detrend


# def run_one():
#     """
#     Measure periods and associated analysis for a single lightcurve
#     """
#     # Read in the light curve, or correctly assign various arrays
#
#     # If there are multiple apertures to choose from, do aperture selection
#
#     # Assign the appropriate arrays to a k2spin LightCurve object,
#     # or stick with a LightKurve object?
#
#     pass

def run_one(t,f,tic=None,secondary_file=None):
    """Run a lomb-scargle analysis on one light curve.

    Inputs:
    -------
    t: array of epochs/time points

    f: array of fluxes corresponding to time points in t

    tic: object identifier

    Returns:
    --------
    fund_period, fund_power: floats
        period and power corresponding to the highest periodogram peak

    sig_period, sig_power: floats
        period and power corresponding to the highest periodogram peak
        that is a) higher than the bootstrap significance theshold, and
        b) higher than N(=100) nearby points as selected with argrelextrema

    sigma: float
        bootstrap significance threshold
    """
    logging.info(tic)

    ylims = np.percentile(f,[0.5,99.5])

    fig = plt.figure(figsize=(8,10))
    base_grid = gridspec.GridSpec(2,1,height_ratios=[2,3])
    if tic is not None:
        plt.suptitle("TIC {0}".format(tic),fontsize="x-large")

    top_grid = gridspec.GridSpecFromSubplotSpec(2,1,subplot_spec=base_grid[0])
    # Just plot the light curve
    ax = plt.subplot(top_grid[0])
    ax.plot(t,f,'k.')
    print(ax.get_xlim())
    ax.set_ylim(ylims)

    # Run the lomb-scargle periodogram on the light curve
    ls_out = prot.run_ls(t,f,np.ones_like(f),0.1,prot_lims=[0.1,70],
                         run_bootstrap=True)
    # unpack lomb-scargle results
    fund_period, fund_power, periods_to_test, periodogram, aliases, sigmas = ls_out
    logging.info("Prot={0:.3f} Power={1:.3f}".format(fund_period,fund_power))


    # Find all peaks in the periodogram
    peak_locs = argrelextrema(periodogram,np.greater,order=100)
    print(len(peak_locs[0]),periods_to_test[np.argmax(peak_locs[0])])

    # Only keep significant peaks (use bootstrap significance levels)
    sig_locs = peak_locs[0][periodogram[peak_locs[0]]>sigmas[0]]
    sig_periods = periods_to_test[sig_locs]
    sig_powers = periodogram[sig_locs]

    # Plot the periodogram
    ax = plt.subplot(top_grid[1])
    ax.plot(periods_to_test,periodogram,'k-')
    ax.axvline(fund_period,color="r",linestyle=":",linewidth=2)
    ax.axhline(sigmas[0],color="grey",linestyle="-.",linewidth=2)
    ax.set_xscale("log")
    ax.set_xlim(0.1,70)
    # Mark significant peaks (if any) on the periodogram
    num_sig = len(sig_locs)

    if num_sig>0:
        plt.plot(sig_periods,sig_powers*1.1,'kv')
        # What's the most powerful of the significant peaks?
        most_significant = np.argmax(sig_powers)
        most_sig_period = sig_periods[most_significant]
        most_sig_power = sig_powers[most_significant]

        # If there are two or more, also record the second one
        if num_sig>1:
            trim_periods = np.delete(sig_periods, most_significant)
            trim_powers = np.delete(sig_powers, most_significant)
            second_significant = np.argmax(trim_powers)
            sec_period = trim_periods[second_significant]
            sec_power = trim_powers[second_significant]
        else:
            sec_period, sec_power = -9999, -9999

        # Record all secondary peaks, just in case
        if secondary_file is not None:
            for i in range(num_sig):
                secondary_file.write("{0},{1:.4f},{2:.4f},{3:.4f}\n".format(tic,
                                     sig_periods[i],sig_powers[i],sigmas[0]))

    else:
        most_sig_period, most_sig_power = -9999,-9999
        sec_period, sec_power = -9999, -9999

    # Count the number of significant periods besides the first
    # A likely harmonic doesn't count against this
    extra_sig = 0

    # Record type of potential harmonic, if applicable
    harm_type = "-"

    if num_sig>1:
        # Compare the two most significant periods
        period_ratio = sec_period / most_sig_period
        power_ratio = sec_power / most_sig_power

        # Is the secondary peak a 1/2 or 2x harmonic?
        if abs(period_ratio-0.5)<=0.05:
            harm_type = "half"
            extra_sig = num_sig - 2
        elif abs(period_ratio-2.0)<=0.05:
            harm_type = "dbl"
            extra_sig = num_sig - 2
        else:
            extra_sig = num_sig-1

        if (harm_type!="-") and (power_ratio>0.5):
            harm_type = harm_type+"-maybe"

    # plot phase-folded periods
    num_cols = np.int(np.ceil((len(sig_periods)+1) / 2))
    bottom_grid = gridspec.GridSpecFromSubplotSpec(2,num_cols,
                                                   subplot_spec=base_grid[1])

    # Plot the phase-folded light curve corresponding to the max
    # peak in the periodogram
    ax = plt.subplot(bottom_grid[0,0])
    phased_t = t % fund_period / fund_period
    ax.plot(phased_t,f,'r.')
    ax.set_ylim(ylims)
    ax.set_xlim(0,1)
    ax.set_title(r"P$_0$={0:.2f}".format(fund_period))

    # Now plot the phase-folded light curves for all other significant peaks
    row = 0
    for i,per in enumerate(sig_periods[np.argsort(sig_powers)]):
        if (i+1)==num_cols:
            row = 1
        ax = plt.subplot(bottom_grid[row,i+1-num_cols])
        phased_t = t % per / per
        ax.plot(phased_t,f,'k.')

        ax.set_ylim(ylims)
        ax.set_xlim(0,1)
        ax.set_title("P={0:.2f}".format(per))
        ax.tick_params(labelleft=False)

    plt.subplots_adjust(hspace=0.25)

    return (fund_period, fund_power, most_sig_period, most_sig_power,
            sec_period, sec_power, sigmas[0], extra_sig, harm_type)



# def run_list():
#     """
#     Measure periods and associated analysis for a a list of lightcurves,
#     and save results.
#     """
#     pass

def run_list(list_filenames,lc_types,output_filename,data_dir,plot_dir):
    """ Run a list of TESS files through run_one(), and save results.

    Inputs:
    -------
    list_filenames: list or array of filename strings

    output_filenames: string, giving the output filename for a table of results

    data_dir: directory that contains the data files

    plot_dir: directory to save plots in

    """

    n_files = len(list_filenames)
    fund_periods = np.zeros(n_files)
    fund_powers = np.zeros(n_files)
    sig_periods = np.zeros(n_files)
    sig_powers = np.zeros(n_files)
    sec_periods = np.zeros(n_files)
    sec_powers = np.zeros(n_files)
    thresholds = np.zeros(n_files)
    tics = np.zeros(n_files,np.int64)
    num_sig_peaks = np.zeros(n_files,int)
    harm_types = np.empty(n_files,"S10")
    harm_types[:] = "-"

    sec_filename = output_filename.replace(".csv","_allpeaks.csv")
    sec_file = open(sec_filename,"w")
    sec_file.write("TIC,period,power,threshold\n")

    for i,filename in enumerate(list_filenames):
        # The hlsp files don't have consistent filenames from pipeline
        # to pipeline, so I'll have to retrieve them from the header
        # using lightkurve

        full_file = os.path.join(filename,f"{filename}.fits")
        full_path = os.path.join(data_dir,full_file)
        lc_type = lc_types[i].lower()

        if "pathos"==lc_type:
            flux_col = "ap2_flux_cor"
            lc = lk.read(full_path,quality_bitmask="default",flux_column=flux_col)
            lc = lc.remove_outliers()
            tic = lc.meta["TICID"]
            time,flux = lc.time.value,lc.flux.value
            sector = lc.meta["SECTOR"]
            # print(tic,lc_type)
            # print(time)
            # print(flux)
            if sector==8:
                good = ((time>1519) & (time<1530)) | (time>1536.5)
            elif sector==9:
                # These are estimated
                good = ((time>1545) & (time<1556)) | (time>1558)
            elif sector==10:
                good = ((time>1572) & (time<1582)) | (time>1586)
            time, flux = time[good], flux[good]
        elif "cdips" in lc_type:
            flux_col = "TFA2"
            lc = lk.read(full_path,quality_bitmask="default",flux_column=flux_col)
            lc = lc.remove_outliers()
            tic = lc.meta["TICID"]
            time,flux = lc.time.value,lc.flux.value
            sector = lc.meta["SECTOR"]

            # print(tic,lc_type)
            # print(time)
            # print(flux)
            # if sector==8:
            #     good = ((time < 2458530) | (time>2458536.5))
            # elif sector==9:
            #     good = ((time > 2458545) & (time<2458556)) | (time > 2458558)
            # elif sector==10:
            #     good = ((time > 2458573) & (time<2458582.5)) | (time > 2458588)
            # time, flux = time[good], flux[good]
        elif "qlp" in lc_type:
            flux_col = "sap_flux"
            lc = lk.read(full_path,quality_bitmask="default",
                         flux_column=flux_col)
            tic = lc.meta["TICID"]
            sector = lc.meta["SECTOR"]
            time,flux = lc.time.value,lc.flux.value
        else:
            flux_col = "default"
            lc = lk.read(full_path,quality_bitmask="default")
            lc_type = lc.meta["ORIGIN"]
            sector = lc.meta["SECTOR"]
            time,flux = lc.time.value,lc.flux.value
            if "/" in lc_type:
                lc_type = lc_type.replace("/","_")





        one_out = run_one(time,flux,tic,sec_file)

        # Unpack analysis results
        fund_periods[i],fund_powers[i],sig_periods[i] = one_out[:3]
        sig_powers[i],sec_periods[i],sec_powers[i],thresholds[i] = one_out[3:7]
        num_sig_peaks[i],harm_types[i] = one_out[7:]
        tics[i] = tic

        # Save and close the plot files
        print(lc_type,sector)
        plt.savefig("{0}TIC{1}_{2}_{3}_{4}.png".format(plot_dir,tic,lc_type,
                                                       flux_col,sector),
                    bbox_inches="tight")
        plt.close()

        # if i>=10:
        #     break

    sec_file.close()

    data = {"TIC": tics,
            "fund_period": fund_periods,
            "fund_power": fund_powers,
            "sig_period": sig_periods,
            "sig_power": sig_powers,
            "sec_period": sec_periods,
            "sec_power": sec_powers,
            "threshold": thresholds,
            "num_sig": num_sig_peaks,
            "harm_type": harm_types}
    formats = {
            "fund_period": "%0.4f",
            "fund_power": "%0.4f",
            "sig_period": "%0.4f",
            "sig_power": "%0.4f",
            "sec_period": "%0.4f",
            "sec_power": "%0.4f",
            "threshold": "%0.6f"}

    names = ["TIC","fund_period","fund_power",
            "sig_period","sig_power","sec_period","sec_power",
            "num_sig","harm_type","threshold"]

    pickle_file = open(output_filename.replace(".csv",".pkl"),"wb")
    pickle.dump(data,pickle_file)
    pickle_file.close()

    at.write(data,output_filename,names=names,
             formats=formats,delimiter=",")




if __name__=="__main__":

    today = date.isoformat(date.today())
#    logging.basicConfig(level=logging.INFO)

    # Arguments: scriptname, file list, cluster name, output file (optional)

    # Retrieve the input list
    if len(sys.argv)<3:
        print("Please provide a list of light curve files and cluster name")
    else:
        listfile = at.read(sys.argv[1])
        file_list = listfile["obs_id"]
        lc_types = listfile["author"]
        cluster = sys.argv[2]

    if len(sys.argv)>3:
        outfile0 = sys.argv[3]
        outfile = "{0}_{1}.csv".format(outfile0[:-4],today)
    else:
        outfile = "tables/{0}_output_{1}.csv".format(cluster,today)

    # # Break down the data set into subsets for parallel processing
    # arrayid = int(os.getenv("PBS_ARRAYID",0))
    #
    # mini = (arrayid - 1) * 35
    # maxi = min(mini + 35, len(file_list))
    # if arrayid==0:
    #     mini = 0
    #     maxi = len(file_list)
    # else:
    #     outfile = outfile.replace(".csv","_{0}.csv".format(arrayid))
    #
    # print("Array, min(i), max(i)")
    # print(arrayid, mini, maxi)
    mini, maxi = 0,50

    sub_list = file_list[mini:maxi]
    sub_types = lc_types[mini:maxi]

    base_path = "./"# "/vega/astro/users/sd2706/k2/"
    data_path = os.path.expanduser("~/.lightkurve-cache/mastDownload/HLSP/")
    plot_path = base_path+"plots/"

    print(sub_list)

    run_list(sub_list,sub_types,base_path+outfile,data_path,plot_path)
