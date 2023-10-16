import os, glob, pathlib

import numpy as np
import astropy.io.ascii as at
from astropy.table import Table, vstack, join
import astropy.units as u
import matplotlib.pyplot as plt

import tess_young
from tess_young.get_const import *
_DIR = pathlib.Path(tess_young.__file__).resolve().parent.parent
plt.style.use(os.path.join(_DIR,'paper.mplstyle'))


def analyze_all_injections(nonvar,pipeline="CDIPS",which="faint"):
    """
    Analyze the results of injection tests on all stars in the nonvar table,
    for the indicated pipeline. 
    """

    table_dir = os.path.join(_DIR,"tables/")
    param_string_wild = f"injection_results_{pipeline}_{which}*.csv"
    search_string = os.path.join(table_dir,param_string_wild)
    # print(search_string)
    output_files = sorted(glob.glob(search_string))

    dat = at.read(os.path.join(_DIR,"tab_all_stars.csv"))
    dat_sub = dat["TIC","GAIAEDR3_G","GAIAEDR3_BP","GAIAEDR3_RP","TIC_Tmag"]

    results_raw = []

    for i, filename in enumerate(output_files):
        tic = filename.split("/")[-1].split("_")[-1].split(".")[0]
        # print(tic)
        # print(int(tic))

        res = at.read(filename)
        res["TIC"] = np.full(len(res),fill_value=tic,dtype="int64")

        results_raw.append(res)

    inj_res00 = vstack(results_raw)


    # Add in brightness
    if "GAIAEDR3_G" in nonvar.dtype.names:
        inj_res0 = join(inj_res00,dat_sub,keys=["TIC"])
        inj_res = join(inj_res0,nonvar,keys=["TIC","GAIAEDR3_G"],table_names=["test","in"])

        colnames = ["Amp","deltaM","TIC_Tmag"]
    else:
        inj_res = join(inj_res00,nonvar,keys="TIC",table_names=["test","in"])
        inj_res.rename_column("Sig","Sig_test")
        colnames = ["Amp","deltaM","Vmag"]

    min_amp,max_amp = 1e-2,2e-1

    bins = {"Amp":np.linspace(min_amp,max_amp,15),
            "deltaM":np.linspace(1,4,13),
            "GAIAEDR3_G":np.linspace(12,18,11),
            "TIC_Tmag":np.linspace(11,19,13),
            "Vmag":np.linspace(8,14,11)
            }

    det = inj_res["Sig_test"]==1
    det_corr = det & (inj_res["Corr"]==1)

    for colname in colnames:

        all_hist, all_edges = np.histogram(inj_res[colname],bins=bins[colname])
        det_hist, det_edges = np.histogram(inj_res[colname][det_corr],bins=bins[colname])

        # plt.figure()
        # ax = plt.subplot(111)
        # ax.hist(inj_res[colname][det],bins=bins[colname])
        # ax.set_title(f"{colname} Raw Detections")
        # ax.set_xlabel(f"{colname}")
        # plt.savefig(os.path.join(_DIR,f"plots/injection_results_{pipeline}_{which}_{colname}_raw.png"))


        plt.figure()
        ax = plt.subplot(111)
        det_frac = det_hist / all_hist
        det_frac[np.isnan(det_frac)] = 0
        ax.step(det_edges[:-1],det_frac,where="post")
        ax.set_title(f"{colname} Correct Detection Fraction")
        ax.set_xlabel(f"{colname}")
        plt.savefig(os.path.join(_DIR,f"plots/injection_results_{pipeline}_{which}_{colname}_fraction.png"))

    col1, col2 = "Amp", "deltaM"
    all_h2d, xedges, yedges = np.histogram2d(inj_res[col1],inj_res[col2],bins=[bins[col1],bins[col2]])
    det_h2d, xedges, yedges = np.histogram2d(inj_res[col1][det],inj_res[col2][det],bins=[bins[col1],bins[col2]])

    frac2d = det_h2d / all_h2d
    frac2d[np.isnan(frac2d)] = 0
    frac2d[frac2d==0] = 0


    left = xedges[0]
    right = xedges[-1]
    bottom = yedges[0]
    top = yedges[-1]
    extents = right, left, bottom, top

    plt.figure()
    ax = plt.subplot(111)
    ax.set_box_aspect(0.9)
    im = ax.imshow(frac2d,extent=extents,origin="lower",aspect="auto")
    # plt.subplots_adjust(left=0.1, right=0.8, top=0.9)
    plt.xlabel(col1)
    plt.ylabel(col2)
    plt.savefig(os.path.join(_DIR,f"plots/injection_results_{pipeline}_{which}_{col1}_{col2}_anysignal.png"))


    # plt.figure()
    # ax = plt.subplot(111)
    # ax.set_box_aspect(0.9)
    # im = ax.imshow(det_h2d,extent=extents,origin="lower",aspect="auto")
    # # plt.subplots_adjust(left=0.1, right=0.8, top=0.9)
    # plt.colorbar(im)
    # plt.savefig(os.path.join(_DIR,f"plots/injection_results_{pipeline}_{which}_{col1}_{col2}_anysignal_raw.png"))


    # Plot where the injected signal is *correctly* detected
    dc_h2d, xedges, yedges = np.histogram2d(inj_res[col1][det_corr],inj_res[col2][det_corr],bins=[bins[col1],bins[col2]])

    fracc2d = dc_h2d / all_h2d
    fracc2d[np.isnan(fracc2d)] = 0
    fracc2d[frac2d==0] = 0

    plt.figure()
    ax = plt.subplot(111)
    ax.set_box_aspect(0.9)
    im = ax.imshow(fracc2d,extent=extents,aspect="auto")
    # plt.subplots_adjust(left=0.1, right=0.8, top=0.9)
    plt.colorbar(im,label="Fraction recovered correctly")
    plt.xlabel(col1)
    plt.ylabel(col2)
    plt.savefig(os.path.join(_DIR,f"plots/injection_results_{pipeline}_{which}_{col1}_{col2}_correct.png"))


    # plt.figure()
    # ax = plt.subplot(111)
    # ax.set_box_aspect(0.9)
    # im = ax.imshow(dc_h2d,extent=extents,aspect="auto")
    # # plt.subplots_adjust(left=0.1, right=0.8, top=0.9)
    # plt.colorbar(im)
    # plt.savefig(os.path.join(_DIR,f"plots/injection_results_{pipeline}_{which}_{col1}_{col2}_correct_raw.png"))


    # plt.figure()
    # plt.plot(inj_res["Pin"],inj_res["Pout"],'.',color="grey",alpha=0.5,zorder=-10)
    # sc = plt.scatter(inj_res["Pin"][good],inj_res["Pout"][good],c=inj_res["Amp"][good],
    #                  zorder=10,vmin=min_amp,vmax=max_amp)
    # plt.colorbar(sc)
    # plt.xlabel("Injected period (d)")
    # plt.ylabel("Detected period (d)")
    # plt.savefig(os.path.join(_DIR,f"plots/injection_results_{pipeline}1all.png"))
    # plt.xlim(0,20)
    # plt.ylim(0,20)
    # plt.savefig(os.path.join(_DIR,f"plots/injection_results_{pipeline}1.png"))

    # plt.figure()
    # plt.plot(inj_res["Pin"],inj_res["Pout"],'.',color="grey",alpha=0.5,zorder=-10)
    # sc = plt.scatter(inj_res["Pin"],inj_res["Pout"],c=inj_res["Amp"],
    #                  zorder=10,vmin=min_amp,vmax=max_amp)
    # plt.plot(inj_res["Pin"][~good],inj_res["Pout"][~good],'kx',zorder=20,alpha=0.5)
    # plt.xlim(0,20)
    # plt.ylim(0,20)
    # plt.xlabel("Injected period (d)")
    # plt.ylabel("Detected period (d)")
    # plt.savefig(os.path.join(_DIR,f"plots/injection_results_{pipeline}2.png"))



if __name__=="__main__":

    # infile = os.path.join(_DIR,"catalogs/nonvar_bright_zhou_vach.csv")
    infile = os.path.join(_DIR,"catalogs/nonvar_faint_douglas.csv")
    nonvar = at.read(infile,delimiter=",")
    # print(nonvar.dtype)

    # analyze_all_injections(nonvar,pipeline="QLP",which="bright")
    analyze_all_injections(nonvar,pipeline="CDIPS",which="faint")
