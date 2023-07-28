import os, glob, pathlib
import multiprocessing as mp

import numpy as np
import astropy.io.ascii as at
from astropy.table import Table
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
        lc = search.download()#download_dir="/data/douglaslab/.lightkurve-cache/")
    else:
        print("Search failed")
        return None

    t = lc["time"]
    flat_lc = lc["flux"]
    npts = len(t)

    if check_input:
        # Check the input light curve for signals
        # Run the lomb-scargle periodogram on the light curve
        ls_out = prot.run_ls(t.btjd,flat_lc,np.ones_like(flat_lc),0.1,prot_lims=[0.1,70],
                             run_bootstrap=True)
        # unpack lomb-scargle results
        fund_period, fund_power, periods_to_test, periodogram, aliases, sigmas = ls_out
        with open(f"injection_input_{pipeline}_{wname}.csv","w") as f:
            f.write("TIC,Prot,Pow,threshold,Tmed\n")
            f.write(f"{wname2},{fund_period:.4f},{fund_power:.4f},")
            f.write(f"{sigmas[0]:.4f}")
            tmed = np.nanmedian(flat_lc)
            f.write(f"{tmed:.2f}\n")


    mid = flat_lc + 3*u.mag

    flat_lc_term = 10**(-0.4*flat_lc.to(u.mag).value)

    rng = np.random.default_rng(seed=3738449329237479)

    min_amp, max_amp = 1e-2, 2e-1

    ntests = 100
    inj_res = Table({"Pin":rng.uniform(0.1,20,ntests),
                     "Pout":np.zeros(ntests)*np.nan,
                     "Amp":rng.uniform(min_amp,max_amp,ntests),
                     "Sig":np.zeros(ntests)
                     })

    for i in range(ntests):
        print(i)

        per = inj_res[i]["Pin"]*u.day
        amp = inj_res[i]["Amp"]*u.mag

        freq = 2*np.pi/per
        # print(per,freq)
        sin_term = t.btjd * freq.to(1/u.day).value
        faint_lc = mid + amp*np.sin(sin_term)#+rng.normal(size=npts,scale=1e-2)*u.mag

        # # https://www.astro.keele.ac.uk/jkt/pubs/JKTeq-fluxsum.pdf
        add_lc = -2.5*np.log10(flat_lc_term + 
                               10**(-0.4*faint_lc.to(u.mag).value))
        test_lc = add_lc * u.mag

        # Run the lomb-scargle periodogram on the light curve
        ls_out = prot.run_ls(t.btjd,test_lc,np.ones_like(test_lc),0.1,prot_lims=[0.1,70],
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

    at.write(inj_res,f"injection_{pipeline}_{wname}.csv",
             delimiter=",",overwrite=True)

if __name__=="__main__":

    infile = os.path.join(_DIR,"catalogs/nonvar_faint_douglas.csv")
    nonvar = at.read(infile,delimiter=",")

    array_step = 1
    mini = arrayid * array_step
    maxi = min(mini + array_step, len(nonvar))
    if arrayid==9999:
        mini = 0
        maxi = array_step

    steps = np.arange(mini,maxi)

    for i in steps:
        tic = nonvar["TIC"][i]

        test_one_tic(tic,check_input=True)
