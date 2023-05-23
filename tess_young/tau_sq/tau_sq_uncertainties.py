import os, glob, pathlib

import numpy as np
import matplotlib.pyplot as plt
import astropy.io.ascii as at

from tau_sq_plot import plot_tausq_tracks

import tess_young
from tess_young.get_const import *
_DIR = pathlib.Path(tess_young.__file__).resolve().parent.parent
plt.style.use(os.path.join(_DIR,'paper.mplstyle'))

def analyze_synthetic_obs(outfilebase,compare_model,compare_init):
    compare_colname = f"{compare_model}_{compare_init}"


    _dir = os.path.join(_DIR,"tables/")
    files = glob.glob(f"{_dir}{outfilebase}*.csv")
    files_sorted = np.sort(files)

    # Make a figure and plot all of the individual runs
    fig = plt.figure()
    fig.patch.set_facecolor('w')
    fig.patch.set_alpha(1.0)
    ax = plt.subplot(111)

    n_sets = len(files_sorted)

    best_age = np.zeros(n_sets)*np.nan

    for i,filename in enumerate(files_sorted):
        if filename.endswith("baseline.csv"):
            # This is the reference run
            syn1_file = filename
            del_i = i
        else:
            # These are the synthetic runs
            syn2 = at.read(filename)

            best_loc = np.argmin(syn2[compare_colname])
            best_age[i] = syn2[f"Age_{compare_colname}"][best_loc]
            plot_tausq_tracks(syn2,models_to_plot=[compare_colname],
                              ax=ax)

    best_age = np.delete(best_age,del_i)

    plt.title(outfilebase)
    plt.xlabel("Age [Myr]")
    plt.ylabel("Tau squared")
    plt.savefig(os.path.join(_DIR,f"plots/{outfilebase}_{compare_colname}.png"),dpi=600,bbox_inches="tight")

    # # Histogram the best fit ages (not statistically useful, but interesting)
    # plt.figure()
    # plt.hist(best_age)
    # plt.show()

    # Now do the thing from Naylor & Jeffries where we map the best-fit ages
    # onto their corresponding tau-sq values from the baseline fake dataset
    fig = plt.figure()
    ax = plt.subplot(111)

    if os.path.exists(syn1_file) is False:
        print("Uh oh",syn1_file)
        return
    syn1 = at.read(syn1_file)

    plot_tausq_tracks(syn1,models_to_plot=[compare_colname],ax=ax)

    ax.set_title(f"Baseline {compare_colname}")

    best_tausq = np.zeros(n_sets-1)*np.nan

    # j = np.where(model_names==compare_model)[0][0]
    for i in range(n_sets-1):
        syn_loc = np.where(syn1[f"Age_{compare_colname}"]==best_age[i])[0]
        best_tausq[i] = syn1[compare_colname][syn_loc]
    point_color = "C0" #mapper.to_rgba((j % 3)+1)
    ax.plot(best_age,best_tausq,'.',color=point_color,alpha=0.25,
            label=f"Synthetic results")

    # Now compute the tau-sq value below which 67 (%) of the tests lie
    ts67 = np.percentile(best_tausq,67)

    ax.legend(loc=2)

    ax.axhline(ts67,color="k",ls="--")
    ax.set_xlabel("Model age (Myr)",fontsize=16)
    ax.set_ylabel("tau squared",fontsize=16)

    ax.tick_params(labelsize=12)
    ax.set_xticks(np.arange(0,300,25),minor=True)

    xmin, xmax = np.nanmin(best_age)*0.9, np.nanmax(best_age)*1.1
    good_ages = ((syn1[f"Age_{compare_colname}"]>=xmin) &
                 (syn1[f"Age_{compare_colname}"]<=xmax))
    ymin = np.nanmin(syn1[compare_colname][good_ages])-20
    ymax = np.nanmax(syn1[compare_colname][good_ages])+20
    ax.set_xlim(xmin,xmax)
    ax.set_ylim(ymin,ymax)

    plt.savefig(os.path.join(_DIR,f"plots/{outfilebase}_{compare_colname}_uncertainty.png"),dpi=600,bbox_inches="tight")


    good = syn1[compare_colname]<=ts67
    print(syn1[f"Age_{compare_colname}"][good])
    print(syn1[compare_colname][good])
    print(len(np.where(best_tausq<=ts67)[0]))



if __name__=="__main__":
    from argparse import ArgumentParser
    import logging

    # Define parser object
    parser = ArgumentParser(description="")

    parser.add_argument("outfilebase",help="base string for searching for results")
    parser.add_argument("compare_model",help="base string for searching for results")
    parser.add_argument("compare_init",help="base string for searching for results")

    args = parser.parse_args()

    analyze_synthetic_obs(args.outfilebase,args.compare_model,args.compare_init)