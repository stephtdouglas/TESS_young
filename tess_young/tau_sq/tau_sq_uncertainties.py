import os, glob

import numpy as np
import matplotlib.pyplot as plt
import astropy.io.ascii as at

from tau_sq_plot import plot_tausq_tracks

import ..get_colors
norm, mapper, cmap2, colors, shapes = get_colors.get_colors()
plt.style.use('./paper.mplstyle')
PAPER_DIR = os.path.expanduser("~/my_papers/TESS_young/")


model_names = np.asarray(["UpSco_Mattea2015","UpSco_Mattea2022","UpSco_ZeroTorque",
               "WideHat8Myr_Mattea2015","WideHat8Myr_Mattea2022","WideHat8Myr_ZeroTorque"])
display_names = {"UpSco_Mattea2015":"Matt+15; UpSco initialization",
                 "UpSco_Mattea2022":"Matt+in prep; UpSco initialization",
                 "UpSco_ZeroTorque":"Zero Torque; UpSco initialization",
                 "WideHat8Myr_Mattea2015":"Matt+15; uniform initialization",
                 "WideHat8Myr_Mattea2022":"Matt+in prep; uniform initialization",
                 "WideHat8Myr_ZeroTorque":"Zero Torque; uniform initialization"}

def analyze_synthetic_obs(ref_model,compare_model,n_sets=100,
                          output_filebase="tausq_tracks",best_age_init=80):

    _dir = "./tables/"
    outfilename = f"{output_filebase}_{ref_model}"

    # Make a figure and plot all of the individual runs
    fig = plt.figure()
    fig.patch.set_facecolor('w')
    fig.patch.set_alpha(1.0)
    ax = plt.subplot(111)

    best_age = np.zeros(n_sets)*np.nan

    for i in range(n_sets):
        # syn2_file = os.path.join(_dir,f"tausq_syn_{i:04d}_SYN_{ref_model}_80Myr.csv")
        syn2_file = os.path.join(_dir,f"SYN_binselect_{compare_model}{i:04d}_SYN_{ref_model}_{best_age_init}Myr.csv")
        syn2 = at.read(syn2_file)

        best_loc = np.argmin(syn2[compare_model])
        best_age[i] = syn2[f"Age_{compare_model}"][best_loc]
        plot_tausq_tracks(syn2,models_to_plot=[compare_model],ax=ax)

    plt.title("Synthetic Obs Set 2")
    plt.xlabel("Age [Myr]")
    plt.ylabel("Tau squared")
    plt.savefig(f"plots/{outfilename}_set2.png",dpi=600,bbox_inches="tight")

    # # Histogram the best fit ages (not statistically useful, but interesting)
    # plt.figure()
    # plt.hist(best_age)
    # plt.show()

    # Now do the thing from Naylor & Jeffries where we map the best-fit ages
    # onto their corresponding tau-sq values from the baseline fake dataset
    fig = plt.figure()
    ax = plt.subplot(111)

    syn1_finder = glob.glob(os.path.join(_dir,f"*SYN2*{ref_model}.csv"))
    if len(syn1_finder)==1:
        syn1_file=syn1_finder[0]
    else:
        print("uh oh",syn1_finder)
        return
    # syn2_file = os.path.join(_dir,f"tausq_SYN_bin_{compare_model}{i:04d}_SYN_{ref_model}_{best_age}Myr.csv")
    # syn1_file = os.path.join(_dir,f"tausq_compare_SYN_{ref_model}_80Myr.csv")
    syn1 = at.read(syn1_file)

    plot_tausq_tracks(syn1,models_to_plot=[ref_model],ax=ax)

    ax.set_title(f"Baseline {display_names[ref_model]}")

    best_tausq = np.zeros(n_sets)*np.nan

    j = np.where(model_names==compare_model)[0][0]
    for i in range(n_sets):
        syn_loc = np.where(syn1[f"Age_{compare_model}"]==best_age[i])[0]
        best_tausq[i] = syn1[compare_model][syn_loc]
    point_color = mapper.to_rgba((j % 3)+1)
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
    good_ages = ((syn1[f"Age_{ref_model}"]>=xmin) &
                 (syn1[f"Age_{ref_model}"]<=xmax))
    ymin = np.nanmin(syn1[ref_model][good_ages])-20
    ymax = np.nanmax(syn1[ref_model][good_ages])+20
    ax.set_xlim(xmin,xmax)
    ax.set_ylim(ymin,ymax)

    plt.savefig(f"plots/{outfilename}.png",dpi=600,bbox_inches="tight")


    good = syn1[ref_model]<=ts67
    print(syn1[f"Age_{ref_model}"][good])
    print(syn1[ref_model][good])
    print(len(np.where(best_tausq<=ts67)[0]))



if __name__=="__main__":

    # analyze_synthetic_obs("WideHat8Myr_Mattea2022","WideHat8Myr_Mattea2022")
    analyze_synthetic_obs("WideHat8Myr_Mattea2015","WideHat8Myr_Mattea2015",
                          best_age_init=126,output_filebase="tausq_tessonly")
