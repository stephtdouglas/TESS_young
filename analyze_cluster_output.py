import os, sys
import glob
import itertools

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.cm as cm
import astropy.io.ascii as at
from astropy.table import join,vstack,Table
import astropy.table as table
from scipy import stats
from scipy.interpolate import interp1d

norm = mpl.colors.LogNorm(vmin=0.1, vmax=30)
mapper = cm.ScalarMappable(norm=norm, cmap=cm.viridis)

# Read in and merge the outputs from k2spin

# colors = {"IC_2391": "C0",
#          "IC_2602": "C4",
#          "NGC_2547": "C3",
#          "NGC_2451A": "C2",
#          "Collinder_135": "C1"}

cmap2 = cm.get_cmap("viridis",7)
colors = {"IC_2391": cmap2(0),
         "IC_2602": cmap2(4),
         "NGC_2547": cmap2(3),
         "NGC_2451A": cmap2(2),
         "Collinder_135": cmap2(1)}


shapes= {"IC_2391": "o",
         "IC_2602": "d",
         "NGC_2547": "v",
         "NGC_2451A": "^",
         "Collinder_135": "s"}


mamajek_file = os.path.expanduser("~/Dropbox/Models/mamajek_colors.dat")
mamajek = at.read(mamajek_file,fill_values=[("...","NaN")])
mamajek_mass_to_bp_rp = interp1d(mamajek["Msun"],mamajek["Bp-Rp"],bounds_error=False)

def mass_to_bp_rp(mass):
    return mamajek_mass_to_bp_rp(mass)


def process_cluster(cluster, date, clean_limit=10,
                    return_periodcolor=True, date2=None, to_plot=False):
    """
    Read in results from period measurement, select periods for each star.

    clean_limit: int, 0, 10, 30, 60 (default 10)
        the maximum height of other significant peaks in the periodogram, as a percentage
        of the primary peak
    """
    base_dir = os.path.expanduser(f"~/data/tess/{cluster.lower()}/tables/")
    print(base_dir)
    print(os.path.join(base_dir,f"*{date}*csv"))
    filenames = glob.iglob(os.path.join(base_dir,f"*{date}*csv"))
    if date2 is not None:
        filenames2 = glob.iglob(os.path.join(base_dir,f"*{date2}*csv"))
        # filenames = np.append(filenames,filenames2)
        filenames = itertools.chain(filenames,filenames2)


    all_res = []
    all_peaks0 = []
    for filename in filenames:
        # print(filename)
        if "allpeaks" in filename:
            all_peaks0.append(at.read(filename))
        else:
            all_res.append(at.read(filename))

    all_peaks = vstack(all_peaks0)
    results = vstack(all_res)

    all_peaks.rename_column("lc_type","provenance_name")
    all_peaks.rename_column("sector","sequence_number")
    all_peaks.rename_column("TIC","target_name")

    # print(results.dtype)
    tic_ids = results["target_name"]
    # print(tic_ids)

    # Delete all PATHOS results, since they are dominated by long trends
    pathos = np.where(results["provenance_name"]=="PATHOS")[0]
    results.remove_rows(pathos)


    # Remove that secondary clump in Collinder 135
    # Should probably do this in a smarter fashion in future
    # if cluster=="Collinder_135":
    #     secondary = np.where(results["RA"]<108)

    u_tic = np.unique(results["target_name"])

    N = len(u_tic)
    dtype = [('TIC', 'int64'), ('NProt', 'i4'), ('Nclean', 'i4'), ('Prot', 'f8'), ('Prot_avg', 'f8'), ('Prot_std', 'f8')]
    summary = Table(data=np.zeros(N, dtype=dtype))

    clean_threshold = 0

    for i, tic in enumerate(u_tic):
        loc = results["target_name"]==tic
        summary["TIC"][i] = tic
        summary["NProt"][i] = len(np.where(loc)[0])
        clean = np.where(loc & (results["num_sig_peaks"]==clean_threshold) & (results["sig_periods"]>0))[0]
        summary["Nclean"][i] = len(clean)
        if len(clean)>0:
            clean_short = np.where(loc & (results["num_sig_peaks"]==clean_threshold) & (results["sig_periods"]<=12))[0]

            if len(clean_short)>0:
                summary["Prot"][i] = results["sig_periods"][clean_short][0]
                summary["Prot_avg"][i] = np.mean(results["sig_periods"][clean_short])
                summary["Prot_std"][i] = np.std(results["sig_periods"][clean_short])
            else:
                summary["Prot"][i] = results["sig_periods"][clean][0]
                summary["Prot_avg"][i] = np.mean(results["sig_periods"][clean])
                summary["Prot_std"][i] = np.std(results["sig_periods"][clean])

        else:
            summary["Prot"][i] = -99

    gaia = at.read(f"{cluster}_crossmatch_xmatch_TIC.csv")
    summary_gaia = join(summary,gaia,keys=["TIC"])

    results2 = results.copy()
    results2["clean60"] = np.zeros(len(results2),bool)
    results2["clean30"] = np.zeros(len(results2),bool)
    results2["clean10"] = np.zeros(len(results2),bool)
    results2["clean"] = np.zeros(len(results2),bool)
    results2["half"] = np.zeros(len(results2),int)
    results2["dbl"] = np.zeros(len(results2),int)
    # results2["period_long"] = np.ones(len(results2),float)*-99
    # results2["power_long"] = np.ones(len(results2),float)*-99
    results2["period_short"] = np.ones(len(results2),float)*-99
    results2["power_short"] = np.ones(len(results2),float)*-99

    for i,row in enumerate(results2):

        # print(row)

        # Find all the all_peaks locations, by target_name, provenance_name, and sequence_number
        loc = ((all_peaks["target_name"]==row["target_name"]) &
               (all_peaks["provenance_name"]==row["provenance_name"]) &
               (all_peaks["flux_col"]==row["flux_cols"]) &
               (all_peaks["sequence_number"]==row["sequence_number"]))
        loc = np.where(loc)[0]

    #     print(row["target_name","provenance_name","sequence_number","flux_cols"])

        if len(loc)==0:
    #         print("no peaks found")
    #         print(row["target_name","provenance_name","sequence_number","flux_cols"])
            continue
        elif len(loc)==1:
            results2["clean"][i] = True
            results2["clean10"][i] = True
            results2["clean30"][i] = True
            results2["clean60"][i] = True
            continue

        sub_tab = all_peaks[loc]
    #     print(sub_tab)


        # Are all the other peaks <60% the height of the main peak?
        main = np.argmax(sub_tab["power"])
        max_power = sub_tab["power"][main]
        max_per = sub_tab["period"][main]

        contaminants = np.where(sub_tab["power"]>=(0.6*max_power))[0]
        # one match is just going to be the primary peak itself
        ncontam = len(contaminants) -1

        if ncontam==0:
            results2["clean60"][i] = True
            # print("clean",ncontam,results2["clean"][i])

        else:
            results2["clean60"][i] = False
            # print("not clean",ncontam,results2["clean"][i])

            # Are there potentially half or double harmonics?
            tolerance = 0.05
            contam_diff = abs(max_per-sub_tab["period"][contaminants])/max_per

            if np.any(abs(contam_diff-0.5)<tolerance):
                results2["half"][i] = 1
    #             print("half")
    #             print(contam_diff)

            if np.any(abs(contam_diff-2)<tolerance):
                results2["dbl"][i] = 1
    #             print("double")
    #             print(contam_diff)


        contaminants = np.where(sub_tab["power"]>=(0.3*max_power))[0]
        # one match is just going to be the primary peak itself
        ncontam = len(contaminants) -1

        if ncontam==0:
            results2["clean30"][i] = True
            # print("clean",ncontam,results2["clean"][i])

        else:
            results2["clean30"][i] = False
            # print("not clean",ncontam,results2["clean"][i])

        contaminants = np.where(sub_tab["power"]>=(0.1*max_power))[0]
        # one match is just going to be the primary peak itself
        ncontam = len(contaminants) -1

        if ncontam==0:
            results2["clean10"][i] = True
            # print("clean",ncontam,results2["clean"][i])

        else:
            results2["clean10"][i] = False
            # print("not clean",ncontam,results2["clean"][i])

        # If the highest peak is >13 days, are there other significant periods?
        if max_per>13:
            short = sub_tab["period"]<=13
            if np.any(short):
    #             print("found shorter periods")
                sub_tab2 = sub_tab[short]
                short_max = np.argmax(sub_tab2["power"])
                row["period_short"] = sub_tab2["period"][short_max]
                row["power_short"] = sub_tab2["power"][short_max]
            else:
                continue

    N = len(u_tic)
    dtype = [('TIC', 'int64'), ('NProt', 'i4'), ('Nclean', 'i4'), ('Prot', 'f8'), ('Prot_avg', 'f8'), ('Prot_std', 'f8'),
             ('provenance_name','a30'),('sequence_number','i4'),('flux_cols','a12')]
    summary2 = Table(data=np.zeros(N, dtype=dtype))

    for i, tic in enumerate(u_tic):
        loc = results2["target_name"]==tic
        summary2["TIC"][i] = tic
        summary2["NProt"][i] = len(np.where(loc)[0])
        clean = np.where(loc & (results2[f"clean{clean_limit}"]==True) & (results2["sig_periods"]>0))[0]
        summary2["Nclean"][i] = len(clean)

        clean_cdips = np.intersect1d(clean,np.where(results2["provenance_name"]=="CDIPS")[0])

        if len(clean)>0:
            clean_all = clean
            clean_short = np.where(loc & (results2[f"clean{clean_limit}"]==True) &
                                   (results2["sig_periods"]>0) & (results2["sig_periods"]<=13))[0]

        else:
            summary2["Prot"][i] = -99
            continue

    #     if len(clean_short)==1:
    #         summary2["Prot"][i] = results2["sig_periods"][clean_short][0]
    #         summary2["Prot_avg"][i] = results2["sig_periods"][clean_short][0]
    #         summary2["Prot_std"][i] = 0.0
    #         summary2["provenance_name"][i] = results2["provenance_name"][clean_short][0]
    #         summary2["sequence_number"][i] = results2["sequence_number"][clean_short][0]
    #         summary2["flux_cols"][i] = results2["flux_cols"][clean_short][0]
    #     elif len(clean_short)>1:
    #         sub_tab = results2[clean_short]
    #         best = np.argmax(sub_tab["sig_powers"])
    #         summary2["Prot"][i] = sub_tab["sig_periods"][best]
    #         summary2["Prot_avg"][i] = np.mean(sub_tab["sig_periods"])
    #         summary2["Prot_std"][i] = np.std(sub_tab["sig_periods"])
    # #         if len(np.unique(sub_tab["provenance_name"]))==1:
    # #             summary2["provenance_name"][i] = sub_tab["provenance_name"][0]
    #         summary2["provenance_name"][i] = str(np.asarray(np.unique(sub_tab["provenance_name"])))
        if len(clean)==1:
            summary2["Prot"][i] = results2["sig_periods"][clean][0]
            summary2["Prot_avg"][i] = results2["sig_periods"][clean][0]
            summary2["Prot_std"][i] = 0.0
            summary2["provenance_name"][i] = results2["provenance_name"][clean][0]
            summary2["sequence_number"][i] = results2["sequence_number"][clean][0]
            summary2["flux_cols"][i] = results2["flux_cols"][clean][0]
        elif len(clean)>1:
            sub_tab = results2[clean]
            best = np.argmax(sub_tab["sig_powers"])
            summary2["Prot"][i] = sub_tab["sig_periods"][best]
            summary2["Prot_avg"][i] = np.mean(sub_tab["sig_periods"])
            summary2["Prot_std"][i] = np.std(sub_tab["sig_periods"])
    #         if len(np.unique(sub_tab["provenance_name"]))==1:
    #             summary2["provenance_name"][i] = sub_tab["provenance_name"][0]
            summary2["provenance_name"][i] = str(np.asarray(np.unique(sub_tab["provenance_name"])))
        else:
    #         summary2["Prot"][i] = -99
            summary2["Prot"][i] = results2["sig_periods"][clean][0]
            summary2["Prot_avg"][i] = np.mean(results2["sig_periods"][clean])
            summary2["Prot_std"][i] = np.std(results2["sig_periods"][clean])

    summary2_gaia = join(summary2,gaia,keys=["TIC"])

    # OH DUH. The summary table is already cleaned.
    # This is just telling me the number of clean periods included in the average
    std_avg = summary_gaia["Prot_std"]/summary_gaia["Prot_avg"]
    std_avg2 = summary2_gaia["Prot_std"]/summary2_gaia["Prot_avg"]
    clean = (summary_gaia["Nclean"]>0) & (std_avg<=0.1)
    clean2 = (summary2_gaia["Nclean"]>0) & (std_avg2<=0.1)

    bp_rp = summary_gaia["GAIAEDR3_BP"] - summary_gaia["GAIAEDR3_RP"]
    bp_rp2 = summary2_gaia["GAIAEDR3_BP"] - summary2_gaia["GAIAEDR3_RP"]

    if to_plot:
        plt.figure()
        plt.plot(bp_rp,summary_gaia["Prot"],'.',color="grey",alpha=0.5)
        plt.plot(bp_rp[clean],summary_gaia["Prot"][clean],'d',ms=5,zorder=5)
        plt.plot(bp_rp2[clean2],summary2_gaia["Prot"][clean2],'o')
        print(len(np.where(clean2)[0]))

        plt.ylim(0.1,50)
        plt.xlim(0.5,3.5)
        plt.yscale("log")

        plt.xlabel(r"G$_{BP}$ - G$_{RP}$")
        plt.ylabel("Period (d)")

        ax = plt.gca()
        ax.axhline(12,linestyle="--",color="tab:grey")

        plt.title(cluster)
        plt.savefig(f"plots/periodmass_{cluster}_clean{clean_limit}.png")

    if return_periodcolor:

        # return bp_rp[clean],summary_gaia["Prot"][clean]
        return bp_rp2[clean2],summary2_gaia["Prot"][clean2]

    else:
        return summary2_gaia, clean2, results2

def read_cluster_visual(cluster, date, clean_limit=None,
                        return_periodcolor=True, date2=None, to_plot=False,
                        which=None):
    # Read in my visual inspection results
    if which is None:
        vis_file = f"tables/{cluster}_{date}_results_comments.csv"
    else:
        vis_file = f"tables/{cluster}_{date}_results_comments{which}.csv"
    vis = at.read(vis_file,delimiter=",")
    good = np.where(vis["Select"].mask==False)[0]
    # print(len(good))

    # Limit the table to only the light curves I analyzed
    vis = vis[good]
    vis.rename_column("\ufefftarget_name","TIC")

    vis["final_period"] = np.copy(vis["sig_periods"])
    vis["final_Q"] = np.copy(vis["Q"])
    replace2 = (vis["Q"]==2) & ((vis["Q2"]==1) | (vis["Q2"]==0))
    replace3 = (vis["Q"]==2) & ((vis["Q3"]==1) | (vis["Q3"]==0))
    vis["final_period"][replace2] = vis["sec_periods"][replace2]
    vis["final_Q"][replace2]==vis["Q2"][replace2]
    vis["final_period"][replace3] = -99 # Not bothering with peaks file right now
    vis["final_Q"][replace3]==vis["Q3"][replace3]

    x_file = f"{cluster}_crossmatch_xmatch_TIC.csv"
    xmatch = at.read(x_file,delimiter=",")
    # print(xmatch.dtype)
    if cluster!="Collinder_135":
        xmatch = xmatch["angDist","GAIAEDR3_ID","GAIAEDR3_RA","GAIAEDR3_DEC",
                        "GAIAEDR3_PMRA","GAIAEDR3_PMDEC","GAIAEDR3_G","GAIAEDR3_BP",
                        "GAIAEDR3_RP","GAIAEDR3_RUWE","GAIAEDR3_G_CORRECTED",
                        "MemBool","angDist_GES","prob_p","angDist_Cantat-Gaudin",
                        "proba","TIC"]
    else:
        xmatch = xmatch["angDist","GAIAEDR3_ID","GAIAEDR3_RA","GAIAEDR3_DEC",
                        "GAIAEDR3_PMRA","GAIAEDR3_PMDEC","GAIAEDR3_G","GAIAEDR3_BP",
                        "GAIAEDR3_RP","GAIAEDR3_RUWE","GAIAEDR3_G_CORRECTED",
                        "MemBool","angDist_Cantat-Gaudin",
                        "proba","TIC"]

    match = join(vis,xmatch,join_type="left",keys=["TIC"],table_names=["vis","xmatch"])

    bp_rp = match["GAIAEDR3_BP"] - match["GAIAEDR3_RP"]
    clean = match["final_Q"]==0
    if clean_limit is not None:
        if clean_limit==0:
            all_possible = (match["final_Q"]==0)
        else:
            all_possible = (match["final_Q"]==0) | (match["final_Q"]==1)
    else:
        all_possible = (match["final_Q"]==0) | (match["final_Q"]==1)

    if to_plot:
        plt.figure()
        # plt.plot(bp_rp,match["final_period"],'.',color="grey",alpha=0.5,
        #          label="all detections")
        plt.plot(bp_rp[match["final_Q"]==1],match["final_period"][match["final_Q"]==1],
                 shapes[cluster],color=colors[cluster],ms=6,zorder=5,mfc="none",
                 label="Possible TESS")
        plt.plot(bp_rp[match["final_Q"]==0],match["final_period"][match["final_Q"]==0],
                 shapes[cluster],color=colors[cluster],ms=6,zorder=6,
                 label="Confident TESS")

        plt.legend(loc=2)

        plt.ylim(0.1,50)
        plt.xlim(0.5,3.5)
        plt.yscale("log")

        plt.xlabel(r"G$_{BP}$ - G$_{RP}$")
        plt.ylabel("Period (d)")

        ax = plt.gca()
        ax.axhline(12,linestyle="--",color="tab:grey")

        plt.title(cluster)
        plt.savefig(f"plots/periodmass_{cluster}_visual.png")

    if return_periodcolor:

        return bp_rp[all_possible],match["final_period"][all_possible]

    else:
        return match

def id_solar(bp_rp):
    # return (bp_rp<=1.2) & (bp_rp>=0.8)
    return (bp_rp<=0.9) & (bp_rp>=0.78)


def write_results():

    clusters = ["IC_2391","Collinder_135","NGC_2451A","NGC_2547","IC_2602"]
    dates = ["2021-06-22","2021-06-18","2021-06-21","2021-06-21","2021-07-02"]
    dates2 = [None,None,None,None,"2021-07-03"]

    for i in range(5):
        print(clusters[i],dates[i],dates2[i])
        summary, clean, results = process_cluster(clusters[i],dates[i],
                                                  return_periodcolor=False,
                                                  date2=dates2[i])

        at.write(results,f"tables/{clusters[i]}_{dates[i]}_results_raw.csv",
                 delimiter=",")

def write_peak_results():

    clusters = ["IC_2391","Collinder_135","NGC_2451A","NGC_2547","IC_2602"]
    dates = ["2021-06-22","2021-06-18","2021-06-21","2021-06-21","2021-07-02"]
    dates2 = [None,None,None,None,"2021-07-03"]


    for i in range(5):
        print(clusters[i],dates[i],dates2[i])
        cluster, date, date2 = clusters[i],dates[i],dates2[i]

        base_dir = os.path.expanduser(f"~/data/tess/{cluster.lower()}/tables/")
        print(base_dir)
        print(os.path.join(base_dir,f"*{date}*csv"))
        filenames = glob.iglob(os.path.join(base_dir,f"*{date}*csv"))
        if date2 is not None:
            filenames2 = glob.iglob(os.path.join(base_dir,f"*{date2}*csv"))
            # filenames = np.append(filenames,filenames2)
            filenames = itertools.chain(filenames,filenames2)

        all_peaks0 = []
        for filename in filenames:
            # print(filename)
            if "allpeaks" in filename:
                all_peaks0.append(at.read(filename))
            else:
                continue

        all_peaks = vstack(all_peaks0)

        all_peaks.rename_column("lc_type","provenance_name")
        all_peaks.rename_column("sector","sequence_number")
        all_peaks.rename_column("flux_col","flux_cols")

        all_peaks_final = table.unique(all_peaks)

        # Sort the table. NOTE that this will make the peak with the highest
        # power the last peak in a given section
        all_peaks_final.sort(["TIC","provenance_name","sequence_number",
                             "flux_cols","power"])
        at.write(all_peaks_final,f"tables/{clusters[i]}_{dates[i]}_allpeaks.csv",
                 delimiter=",")

def compare_visual_results(cluster, date):

    plt.figure(figsize=(9,9))

    print(cluster,"\n-------")

    # Read in my visual inspection results
    vis_file = f"tables/{cluster}_{date}_results_comments.csv"
    vis = at.read(vis_file,delimiter=",")
    good = np.where(vis["Select"].mask==False)[0]
    print(len(good))

    # Limit the table to only the light curves I analyzed
    vis = vis[good]
    vis.rename_column("\ufefftarget_name","TIC")

    vis["final_period"] = np.copy(vis["sig_periods"])
    vis["final_Q"] = np.copy(vis["Q"])
    replace2 = (vis["Q"]==2) & ((vis["Q2"]==1) | (vis["Q2"]==0))
    replace3 = (vis["Q"]==2) & ((vis["Q3"]==1) | (vis["Q3"]==0))
    vis["final_period"][replace2] = vis["sec_periods"][replace2]
    vis["final_Q"][replace2]==vis["Q2"][replace2]
    vis["final_period"][replace3] = -99 # Not bothering with peaks file right now
    vis["final_Q"][replace3]==vis["Q3"][replace3]

    # Retrieve the automated results
    for i,clean_limit in enumerate([60,30,10,""]):
        ax1 = plt.subplot(2,2,i+1)
        summary, clean, results = process_cluster(cluster,date,clean_limit=clean_limit,
                                         return_periodcolor=False)

        match = join(vis,summary,join_type="left",keys=["TIC"],
                     table_names=["vis","auto"])

        ax1.plot(match["final_period"][match["final_Q"]==0],
                 match["Prot"][match["final_Q"]==0],'o',color=f"C{i}",
                 label=f"Q=0",ms=5)
        ax1.plot(match["final_period"][match["final_Q"]==1],
                 match["Prot"][match["final_Q"]==1],'o',color=f"C{i}",mfc="none",
                 label=f"Q=1",ms=5)
        ax1.plot(match["final_period"][match["final_Q"]==2],
                 match["Prot"][match["final_Q"]==2],'k+',
                 label=f"Q=2")
        # ax1.plot(match["final_period"][match["final_Q"]==3],
        #          match["Prot"][match["final_Q"]==3],'kx',
        #          label=f"Q=3")
        ax1.legend(loc=4)

        ax1.set_title(f"clean={clean_limit}")

        x = np.linspace(0.1,50,20)
        ax1.plot(x,x,"-",zorder=-5,color="grey")
        ax1.plot(x,x/2,"--",zorder=-5,color="grey")
        ax1.plot(x,x*2,"--",zorder=-5,color="grey")

        ax1.set_xlim(0.1,50)
        ax1.set_ylim(0.1,50)
        ax1.set_xscale("log")
        ax1.set_yscale("log")

        if (i==0) or (i==2):
            ax1.set_ylabel("Automated Period (d)")
        if (i==2) or (i==3):
            ax1.set_xlabel("Visual Inspection Period (d)")

    plt.suptitle(cluster.replace("_"," "),y=0.93)
    plt.savefig(f"plots/visual_compare_{cluster}.png",bbox_inches="tight")
    # plt.show()

def count_targets():
    clusters = ["IC_2391","Collinder_135","NGC_2451A","NGC_2547","IC_2602"]
    dates = ["2021-06-22","2021-06-18","2021-06-21","2021-06-21","2021-07-02"]
    dates2 = [None,None,None,None,"2021-07-03"]

    all_count = 0
    good_count = 0

    for i, cluster in enumerate(clusters):
        date = dates[i]
        # Read in my visual inspection results
        vis_file = f"tables/{cluster}_{date}_results_comments.csv"
        vis = at.read(vis_file,delimiter=",")
        vis.rename_column("\ufefftarget_name","TIC")
        good = np.where((vis["Select"].mask==False) & (vis["Select"]=="y"))[0]

        all = len(np.unique(vis["TIC"]))

        print("\n",cluster)
        print(all,"with TESS data")
        print(len(good),"with periods")

        all_count += all
        good_count += len(good)

    print("\n",all_count,"with TESS data")
    print(good_count,"with periods")




if __name__=="__main__":

    # # write_results()
    # plot_results()
    # # compare_visual_results(cluster="NGC_2451A",date = "2021-06-21")
    # # compare_visual_results(cluster="NGC_2547",date = "2021-06-21")
    # # compare_visual_results(cluster="IC_2391",date = "2021-06-22")
    #
    # plot_all(clean_limit=0)
    # plot_periodcolor_models(clean_limit=0)
    # ax = plot_periodcolor_histogram(clean_limit=0,to_plot_indiv=False)
    # plt.savefig("plots/periodmass_histogram.png",bbox_inches="tight")

    # count_targets()

    write_peak_results()
