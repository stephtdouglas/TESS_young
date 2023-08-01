import os, glob, pathlib
import multiprocessing as mp

import numpy as np
import astropy.io.ascii as at
from astropy.table import Table
from astropy.time import Time
import astropy.units as u
from lightkurve import LightCurve
import matplotlib.pyplot as plt

import k2spin
from k2spin import prot

import tess_young
from tess_young.get_const import *
_DIR = pathlib.Path(tess_young.__file__).resolve().parent.parent
plt.style.use(os.path.join(_DIR,'paper.mplstyle'))

arrayid = int(os.getenv("SLURM_ARRAY_TASK_ID",9999))
jobid = int(os.getenv("SLURM_JOB_ID",9999))
print(arrayid,jobid)

def test_one_tic(tic,pipeline="CDIPS",which="faint"):

    ticname = f"TIC {tic}"
    wname = f"TIC_{tic}"
    wname2 = str(tic)

    print(ticname)
    

    basefile = os.path.join(_DIR,f"tables/injection_input_{pipeline}_{which}_{wname}.csv")

    if os.path.exists(basefile)==False:
        print("Baseline file not found. Run setup_injection_tests.py first.")
        sys.exit(42)

    with open(basefile,"r") as f:
        fline = f.readline()
        lc_file = fline[2:]

    if os.path.exists(lc_file)==False:
        print("Light curve file not found!")
        print(fline)
        print(lc_file)
        sys.exit(42)

    lc = LightCurve(lc_file[0])
    
    print(type(lc))
    print(lc)
    print(lc.dtype)

    t_raw = lc["time"]
    flat_lc = lc["flux"]
    t = t_raw.btjd

    # CDIPS has everything in magnitudes while QLP is normalized fluxes
    if pipeline=="CDIPS":
        flat_lc_term = 10**(-0.4*flat_lc.to(u.mag).value)

    npts = len(t)

    seed_add = int(tic*1e9)
    print(seed_add)
    rng = np.random.default_rng(seed=3738449329237479+seed_add)

    min_amp, max_amp = 1e-2, 2e-1

    ntests = 5
    inj_res = Table({"Pin":rng.uniform(0.1,20,ntests),
                     "Pout":np.zeros(ntests)*np.nan,
                     "deltaM":rng.uniform(1,4,ntests),
                     "Amp":rng.uniform(min_amp,max_amp,ntests),
                     "Sig":np.zeros(ntests),
                     "Corr":np.zeros(ntests)
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
        if (fund_power > sigmas[0]):
            inj_res["Sig"][i] = 1
        if (per_diff<0.05):
            inj_res["Corr"][i] = 1

    at.write(inj_res,os.path.join(_DIR,f"tables/injection_results_{pipeline}_{which}_{wname}.csv"),
             delimiter=",",overwrite=True)

if __name__=="__main__":

    from argparse import ArgumentParser

    parser = ArgumentParser(description="")
    parser.add_argument("-w", "--which", dest="which", required=True,
                        type=str, help="bright or faint")
    parser.add_argument("-p", "--pipeline", dest="pipeline", required=False,
                        type=str, help="QLP or CDIPS", default="CDIPS")

    args = parser.parse_args()
    print(args)

    which = args.which
    pipeline = args.pipeline

    print(pipeline, which)
    
    if which=="bright":
        infile = os.path.join(_DIR,"catalogs/nonvar_bright_zhou_vach.csv")
    elif which=="faint":
        infile = os.path.join(_DIR,"catalogs/nonvar_faint_douglas.csv")
    else:
        print("selection not found", which)
        sys.exit(42)
    nonvar = at.read(infile,delimiter=",")

    if arrayid==9999:
        i = 0
    else:
        i = arrayid

    tic = nonvar["TIC"][i]
    print(i,tic)
    
    test_one_tic(tic,pipeline=pipeline,which=which)
