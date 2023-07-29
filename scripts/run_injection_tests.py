import os, glob, pathlib
import multiprocessing as mp

import numpy as np
import astropy.io.ascii as at
from astropy.table import Table
from astropy.time import Time
import astropy.units as u
from lightkurve import search_lightcurve
import matplotlib.pyplot as plt

import k2spin
from k2spin import prot

import tess_young
from tess_young.get_const import *
_DIR = pathlib.Path(tess_young.__file__).resolve().parent.parent
plt.style.use(os.path.join(_DIR,'paper.mplstyle'))

arrayid = int(os.getenv("SLURM_ARRAY_TASK_ID",9999))
jobid = int(os.getenv("SLURM_JOB_ID",9999))

def test_one_tic(tic,pipeline="CDIPS",check_input=False):

    ticname = f"TIC {tic}"
    wname = f"TIC_{tic}"
    wname2 = str(tic)

    try:
        search = search_lightcurve(ticname, author=pipeline)
    except HTTPError:
        print(ticname, pipeline, sectors, "MAST server error")
        return None
        
    if len(search)>0:
        lc = search.download()#download_dir="/data2/douglaslab/.lightkurve-cache/")
    else:
        print("Search failed")
        return None

    # print(type(lc))
    # print(lc)
    # print(lc.dtype)

    t_raw = lc["time"]
    flat_lc = lc["flux"]
    t = t_raw.btjd

    # CDIPS has everything in magnitudes while QLP is normalized fluxes
    if pipeline=="CDIPS":
        flat_lc_term = 10**(-0.4*flat_lc.to(u.mag).value)

    npts = len(t)

    # print(t)
    # print(flat_lc)

    if check_input:
        # Check the input light curve for signals
        # Run the lomb-scargle periodogram on the light curve
        ls_out = prot.run_ls(t,flat_lc,np.ones_like(flat_lc),0.1,prot_lims=[0.1,70],
                             run_bootstrap=True)
        # unpack lomb-scargle results
        fund_period, fund_power, periods_to_test, periodogram, aliases, sigmas = ls_out
        with open(os.path.join(_DIR,f"tables/injection_input_{pipeline}_{wname}.csv"),"w") as f:
            f.write("TIC,Prot,Pow,threshold,Tmed\n")
            f.write(f"{wname2},{fund_period:.4f},{fund_power:.4f},")
            f.write(f"{sigmas[0]:.4f}")
            tmed = np.nanmedian(flat_lc)
            f.write(f"{tmed:.2f}\n")

    rng = np.random.default_rng(seed=3738449329237479+arrayid)

    min_amp, max_amp = 1e-2, 2e-1

    ntests = 2#100
    inj_res = Table({"Pin":rng.uniform(0.1,20,ntests),
                     "Pout":np.zeros(ntests)*np.nan,
                     "deltaM":rng.uniform(1,4,ntests),
                     "Amp":rng.uniform(min_amp,max_amp,ntests),
                     "Sig":np.zeros(ntests)
                     })

    for i in range(ntests):
        print(i)

        per = inj_res[i]["Pin"]*u.day
        freq = 2*np.pi/per
        # print(per,freq)
        sin_term = t * freq.to(1/u.day).value

        if pipeline=="CDIPS":
            mid = flat_lc + inj_res[i]["deltaM"]*u.mag
            amp = inj_res[i]["Amp"]*u.mag

            faint_lc = mid + amp*np.sin(sin_term)#+rng.normal(size=npts,scale=1e-2)*u.mag

            # # https://www.astro.keele.ac.uk/jkt/pubs/JKTeq-fluxsum.pdf
            add_lc = -2.5*np.log10(flat_lc_term + 
                                   10**(-0.4*faint_lc.to(u.mag).value))
            test_lc = add_lc * u.mag
            # print(test_lc)
        else:
            flux_ratio = 10**(-0.4*inj_res[i]["deltaM"])

            mid = flat_lc * flux_ratio
            amp = 10**(-0.4*inj_res[i]["Amp"])

            faint_lc = mid + amp*np.sin(sin_term)#+rng.normal(size=npts,scale=1e-2)*u.mag

            add_lc = flat_lc + faint_lc
            test_lc = add_lc
            # print(test_lc)

        # Run the lomb-scargle periodogram on the light curve
        ls_out = prot.run_ls(t,test_lc,np.ones_like(test_lc),0.1,prot_lims=[0.1,70],
                             run_bootstrap=True)
        # unpack lomb-scargle results
        fund_period, fund_power, periods_to_test, periodogram, aliases, sigmas = ls_out

        inj_res["Pout"][i] = fund_period

        per_diff = abs(inj_res["Pout"][i]-inj_res["Pin"][i])/inj_res["Pin"][i]
        # A good detection is significant AND matches the input within 5%
        if (fund_power > sigmas[0]) & (per_diff<0.05):
            inj_res["Sig"][i] = 1

    # good = inj_res["Sig"]==1

    # plt.figure()
    # sc = plt.scatter(inj_res["Pin"],inj_res["Pout"],c=inj_res["Amp"],
    #                  zorder=10,vmin=min_amp,vmax=max_amp)
    # plt.plot(inj_res["Pin"][~good],inj_res["Pout"][~good],'kx',zorder=20)
    # plt.xlabel("Injected period (d)")
    # plt.ylabel("Detected period (d)")
    # plt.savefig(os.path.join(_DIR,"plots/test_injection.png"))

    at.write(inj_res,os.path.join(_DIR,f"tables/injection_results_{pipeline}_{wname}.csv"),
             delimiter=",",overwrite=True)

if __name__=="__main__":

    infile = os.path.join(_DIR,"catalogs/nonvar_bright_zhou_vach.csv")
    # infile = os.path.join(_DIR,"catalogs/nonvar_faint_douglas.csv")
    nonvar = at.read(infile,delimiter=",")

    if arrayid==9999:
        i = 0
    else:
        i = arrayid

    tic = nonvar["TIC"][i]

    test_one_tic(tic,pipeline="QLP",check_input=True)