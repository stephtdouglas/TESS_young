"""
Script to measure periods from all downloaded lightkurve files
"""
import sys, os
from datetime import date
import logging
import pickle

arrayid = int(os.getenv("SLURM_ARRAY_TASK_ID",9999))
jobid = int(os.getenv("SLURM_JOB_ID",9999))
logger = logging.getLogger("measure_periods")
logging.basicConfig(level=logging.DEBUG,
                    filename='/data/douglaslab/script_logs/measure_periods_{0}_{1}.log'.format(jobid,arrayid),
                    format='%(asctime)s %(message)s')
logging.getLogger("matplotlib").setLevel(logging.WARNING)

import matplotlib
matplotlib.use("agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import astropy.io.ascii as at
import astropy.io.fits as fits
from astropy.table import join, Table, unique
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

def run_one(t,f,tic=None,lc_type=None,sector=None,flux_col=None,
            secondary_file=None):
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
    logging.debug(ax.get_xlim())
    ax.set_ylim(ylims)

    # Run the lomb-scargle periodogram on the light curve
    ls_out = prot.run_ls(t,f,np.ones_like(f),0.1,prot_lims=[0.1,70],
                         run_bootstrap=True)
    # unpack lomb-scargle results
    fund_period, fund_power, periods_to_test, periodogram, aliases, sigmas = ls_out
    logging.info("Prot={0:.3f} Power={1:.3f}".format(fund_period,fund_power))


    # Find all peaks in the periodogram
    peak_locs = argrelextrema(periodogram,np.greater,order=100)
    logging.debug(len(peak_locs[0]))
    logging.debug(periods_to_test[np.argmax(peak_locs[0])])

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
                # lc_type,sector,flux_col,
                secondary_file.write("{0},{1},{2},{3},".format(tic,
                                     lc_type,sector,flux_col))
                secondary_file.write("{0:.4f},{1:.4f},{2:.4f}\n".format(
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


def run_list(data_list,output_filename,data_dir,plot_dir,
             flux_cols={"PATHOS":["psf_flux_cor","ap1_flux_cor","ap2_flux_cor","ap3_flux_cor"],
                        "CDIPS":["TFA1","TFA2","TFA3","PCA1","PCA2","PCA3"],
                        "QLP":["sap_flux"]}):
    """ Run a list of TESS files through run_one(), and save results.

    Inputs:
    -------
    list_filenames: list or array of filename strings

    output_filenames: string, giving the output filename for a table of results

    data_dir: string, directory that contains the data files

    plot_dir: string, directory to save plots in

    flux_cols: dictionary containing a list for each lightcurve type in data_list

    """

    n_files = len(data_list)
    raw_data_list = data_list.copy()
    
    data_list["fund_periods"] = np.zeros(n_files)
    data_list["fund_powers"] = np.zeros(n_files)
    data_list["sig_periods"] = np.zeros(n_files)
    data_list["sig_powers"] = np.zeros(n_files)
    data_list["sec_periods"] = np.zeros(n_files)
    data_list["sec_powers"] = np.zeros(n_files)
    data_list["thresholds"] = np.zeros(n_files)
    data_list["tics"] = np.zeros(n_files,np.int64)
    data_list["num_sig_peaks"] = np.zeros(n_files,int)
    data_list["harm_types"] = np.empty(n_files,"S10")
    data_list["harm_types"][:] = "-"
    data_list["flux_cols"] = np.empty(n_files,"S15")
    data_list["flux_cols"][:] = "default"

    # fund_periods = np.zeros(n_files)
    # fund_powers = np.zeros(n_files)
    # sig_periods = np.zeros(n_files)
    # sig_powers = np.zeros(n_files)
    # sec_periods = np.zeros(n_files)
    # sec_powers = np.zeros(n_files)
    # thresholds = np.zeros(n_files)
    # tics = np.zeros(n_files,np.int64)
    # num_sig_peaks = np.zeros(n_files,int)
    # harm_types = np.empty(n_files,"S10")
    # harm_types[:] = "-"

    # TODO: also save sector, lc_type, etc
    sec_filename = output_filename.replace(".csv","_allpeaks.csv")
    sec_file = open(sec_filename,"w")
    sec_file.write("TIC,lc_type,sector,flux_col,period,power,threshold\n")

    for i,row in enumerate(raw_data_list):
        # The hlsp files don't have consistent filenames from pipeline
        # to pipeline, so I'll have to retrieve them from the header
        # using lightkurve

        first = True

        full_file = os.path.join(row["obs_id"],row["productFilename"])
        full_path = os.path.join(data_dir,full_file)
        lc_type = row["provenance_name"]
        sector = row["sequence_number"]
        logging.debug(i)
        logging.debug(lc_type)

        these_flux_cols = flux_cols[lc_type]
        if len(these_flux_cols)==0:
            print(lc_type,"not found in flux_cols")
            continue

        for flux_col in these_flux_cols:
            logging.debug(flux_col)
            lc = lk.read(full_path,quality_bitmask="default",flux_column=flux_col)
            lc = lc.remove_outliers()
            tic = lc.meta["TICID"]
            time,flux = lc.time.value,lc.flux.value


            if "PATHOS"==lc_type:
                if sector==8:
                    good = ((time>1519) & (time<1530)) | (time>1536.5)
                elif sector==9:
                    # These are estimated
                    good = ((time>1545) & (time<1556)) | (time>1558)
                elif sector==10:
                    good = ((time>1572) & (time<1582)) | (time>1586)
                time, flux = time[good], flux[good]


            one_out = run_one(time,flux,tic,lc_type,sector,flux_col,sec_file)

            if first:
                k = i
            else:
                data_list.add_row(data_list[i])
                k = -1

            logging.debug(k)
            
            # Unpack analysis results
            data_list["fund_periods"][k],data_list["fund_powers"][k],data_list["sig_periods"][k] = one_out[:3]
            data_list["sig_powers"][k],data_list["sec_periods"][k],data_list["sec_powers"][k],data_list["thresholds"][k] = one_out[3:7]
            data_list["num_sig_peaks"][k],data_list["harm_types"][k] = one_out[7:]
            data_list["tics"][k] = tic
            data_list["flux_cols"][k] = flux_col

            # Save and close the plot files
            logging.debug(lc_type)
            logging.debug(sector)
            output_plotfile = "{0}TIC{1}_{2}_{3}_{4}.png".format(plot_dir,tic,
                                                           lc_type.upper(),
                                                           flux_col,sector)
            plt.savefig(output_plotfile.replace(" ",""),bbox_inches="tight")
            plt.close()

            logging.debug(output_plotfile.replace(" ","").split("/")[-1])
            
            # If there are multiple aperture types for the same lc file
            # Need to add rows to the output table
            first = False


    sec_file.close()

    # data = {"TIC": tics,
    #         "fund_period": fund_periods,
    #         "fund_power": fund_powers,
    #         "sig_period": sig_periods,
    #         "sig_power": sig_powers,
    #         "sec_period": sec_periods,
    #         "sec_power": sec_powers,
    #         "threshold": thresholds,
    #         "num_sig": num_sig_peaks,
    #         "harm_type": harm_types}
    formats = {
            "fund_periods": "%0.4f",
            "fund_powers": "%0.4f",
            "sig_periods": "%0.4f",
            "sig_powers": "%0.4f",
            "sec_periods": "%0.4f",
            "sec_powers": "%0.4f",
            "thresholds": "%0.6f"}

    # names = ["tics","fund_periods","fund_powers",
    #         "sig_periods","sig_powers","sec_periods","sec_powers",
    #         "num_sig_peaks","harm_types","thresholds"]

    pickle_file = open(output_filename.replace(".csv",".pkl"),"wb")
    pickle.dump(data_list,pickle_file)
    pickle_file.close()

    at.write(data_list,output_filename,#names=names,
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
        sub_info0 = listfile["target_name","provenance_name","sequence_number",
                            "obs_id","productFilename","author"]
        sub_info = unique(sub_info0)
        cluster = sys.argv[2]

    if len(sys.argv)>3:
        outfile0 = sys.argv[3]
        outfile = "{0}_{1}.csv".format(outfile0[:-4],today)
    else:
        outfile = "tables/{0}_output_{1}.csv".format(cluster,today)

    # Break down the data set into subsets for parallel processing
    # TODO: put these back at 50 or 100, longer jobs are easier to run on this cluster
#     arrayid = int(os.getenv("SLURM_ARRAY_TASK_ID",0))

    array_step = 10
    mini = arrayid * array_step
    maxi = min(mini + array_step, len(sub_info))
    if arrayid==9999:
        mini = 0
        maxi = array_step #len(sub_info)
    else:
        outfile = outfile.replace(".csv","_{0}.csv".format(arrayid))

#    logging.basicConfig(filename='/data/douglaslab/script_logs/measure_periods_{0}.log'.format(arrayid))
        
    logging.info("Array, min(i), max(i)")
    logging.info(arrayid)
    logging.info(mini)
    logging.info(maxi)
    
    sub_list = sub_info[mini:maxi]

#    base_path = "./"# "/vega/astro/users/sd2706/k2/"
#    data_path = os.path.expanduser("~/.lightkurve-cache/mastDownload/HLSP/")
#    plot_path = base_path+"plots/"

    base_path = "/data/douglaslab/tess/ic2391/"
    data_path = "/data/douglaslab/.lightkurve-cache/mastDownload/HLSP/"
    plot_path = os.path.join(base_path,"plots/")
    
    logging.info(sub_list)

    run_list(sub_list,base_path+outfile,data_path,plot_path)
