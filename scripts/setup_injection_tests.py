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
print(arrayid,jobid)

def setup_one_tic(tic,pipeline="CDIPS",which="faint",check_input=False):

    ticname = f"TIC {tic}"
    wname = f"TIC_{tic}"
    wname2 = str(tic)

    print(ticname)
    
    try:
        search = search_lightcurve(ticname, author=pipeline)
    except HTTPError:
        print(ticname, pipeline, sectors, "MAST server error")
        return None

    print(search)
    
    if len(search)>0:
        lc = search.download(download_dir="/data2/douglaslab/.lightkurve-cache/")
    else:
        print("Search failed")
        return None

    print(lc.meta["FILENAME"])

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
        print("Check the input light curve for signals")
        # Run the lomb-scargle periodogram on the light curve
        ls_out = prot.run_ls(t,flat_lc,np.ones_like(flat_lc),0.1,prot_lims=[0.1,70],
                             run_bootstrap=True)
        # unpack lomb-scargle results
        fund_period, fund_power, periods_to_test, periodogram, aliases, sigmas = ls_out
        with open(os.path.join(_DIR,f"tables/injection_input_{pipeline}_{which}_{wname}.csv"),"w") as f:
            f.write("# "+lc.meta["FILENAME"])
            f.write("\nTIC,Prot,Pow,threshold,Tmed\n")
            f.write(f"{wname2},{fund_period:.4f},{fund_power:.4f},")
            f.write(f"{sigmas[0]:.4f},")
            tmed = np.nanmedian(flat_lc)
            f.write(f"{tmed:.2f}\n")

    return None


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
        for i, tic in enumerate(nonvar["TIC"]):
            print(i,tic)
            setup_one_tic(tic,pipeline=pipeline,which=which,check_input=True)
    else:
        i = arrayid

        tic = nonvar["TIC"][i]
        print(i,tic)
    
        setup_one_tic(tic,pipeline=pipeline,which=which,check_input=True)
