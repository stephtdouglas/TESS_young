import os, sys
import glob
import itertools

import numpy as np
from numpy.random import default_rng
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.cm as cm
import astropy.io.ascii as at
import astropy.units as u
from astropy.table import join,vstack,Table
from scipy import stats
from scipy.interpolate import interp1d

norm = mpl.colors.LogNorm(vmin=0.1, vmax=30)
mapper = cm.ScalarMappable(norm=norm, cmap=cm.viridis)

from analyze_cluster_output import colors, shapes
from analyze_cluster_output import mass_to_bp_rp, id_solar
from analyze_cluster_output import read_cluster_visual
from plot_periods import plot_periodcolor_histogram
from model_data_percentiles import young_stars_init, zams_percentiles, zams_percentiles_subset

def plot_periodcolor_models(clean_limit=10):

    ### Matt models, USco init
    ax = plot_periodcolor_histogram(clean_limit)
    # ax = plot_periodcolor(clean_limit)
    model = "Matt et al. (2020), USco initialization"
    model_f = "Matt2020_USco"
    mfile = os.path.expanduser("~/Dropbox/Models/Mattea2020_00040Myr_USco.txt")
    matt = at.read(mfile,names=["mass","prot"])
    matt_bp_rp = mass_to_bp_rp(matt["mass"])
    ax.plot(matt_bp_rp,matt["prot"],'*',color="teal",alpha=0.75)
    ax.axvline(0.9,color="k",linestyle=":")
    ax.axvline(0.78,color="k",linestyle=":")
    ax.set_title(model)
    plt.savefig(f"plots/periodmass_all_clean{clean_limit}_{model_f}.png")
    plt.close("all")

    ### Matt models
    ax = plot_periodcolor_histogram(clean_limit)
    # ax = plot_periodcolor(clean_limit)
    model = "Matt et al. (2020)"
    model_f = "Matt2020"
    mfile = os.path.expanduser("~/Dropbox/Models/Mattea2020_00040Myr.txt")
    matt = at.read(mfile,names=["mass","prot"])
    matt_bp_rp = mass_to_bp_rp(matt["mass"])
    ax.plot(matt_bp_rp,matt["prot"],'*',color="teal",alpha=0.75)
    ax.axvline(0.9,color="k",linestyle=":")
    ax.axvline(0.78,color="k",linestyle=":")
    ax.set_title(model)
    plt.savefig(f"plots/periodmass_all_clean{clean_limit}_{model_f}.png")
    plt.close("all")

    ### Gossage models, M15 wind
    ax = plot_periodcolor_histogram(clean_limit)
    # ax = plot_periodcolor(clean_limit)
    model = "Gossage et al. (2021), M15 Wind"
    model_f = "Gossage21_M15"

    goss_mass = np.zeros(15*11)
    goss_prot = np.zeros(15*11)

    mass = ["010","020","030","040","050","060","070","080","090","100","110","120","130"]
    init_periods = ["0.3d","0.4d","0.6d","0.8d","1.5d","12d","3d","4.5d","6d","8d"]
    desired_age = 40

    ct = 0
    for m in mass:
        for p in init_periods:
            goss_mass[ct] = float(m)/100
            hfile = f"models/MIST_M15_mbraking_feh_p0.00_histfiles/00{m}M_dir/LOGS/{goss_mass[ct]:.1f}M_history_prot{p}0_M15_d02d4_0.5hp_1.4d30K_0.22m_2.6p_3myr_alt3prot_henyey1.82mlta.data"
            # read in the history file (for a particular mass, initial period, and braking formalism)
            data = np.genfromtxt(hfile, skip_header=5, names=True)

            # extract rotation periods at all ages
            periods = 2*np.pi/(86400*data['surf_avg_omega'])

            # mark the index of a desired age [Myr]
            ages_myr = data['star_age']/1e6
            age_diffs = abs(ages_myr - desired_age)
            age_index = np.where(age_diffs == min(age_diffs))[0][0]

            # then extract the rotation period at the desired age
            goss_prot[ct] = periods[age_index]
            ct += 1

    goss_bp_rp = mass_to_bp_rp(goss_mass)
    ax.plot(goss_bp_rp,goss_prot,'*',color="teal",ms=12)
    ax.axvline(0.9,color="k",linestyle=":")
    ax.axvline(0.78,color="k",linestyle=":")
    ax.set_title(model)
    plt.savefig(f"plots/periodmass_all_clean{clean_limit}_{model_f}.png")
    plt.close("all")

    ### Gossage models, G18 wind
    ax = plot_periodcolor_histogram(clean_limit)
    # ax = plot_periodcolor(clean_limit)
    model = "Gossage et al. (2021), G18 Wind"
    model_f = "Gossage21_G18"

    goss_mass = np.zeros(15*11)
    goss_prot = np.zeros(15*11)

    mass = ["010","020","030","040","050","060","070","080","090","100","110","120","130"]
    init_periods = ["0.3d","0.4d","0.6d","0.8d","1.5d","12d","3d","4.5d","6d","8d"]
    desired_age = 40

    ct = 0
    for m in mass:
        for p in init_periods:
            goss_mass[ct] = float(m)/100
            hfile = f"models/MIST_G18_mbraking_feh_p0.00_histfiles/00{m}M_dir/LOGS/{goss_mass[ct]:.1f}M_history_prot{p}0_G18_d02d4_0.5hp_1d_3.2d41c_0.5b_0.03a_3myr_alt3prot_henyey1.82mlta.data"
            # read in the history file (for a particular mass, initial period, and braking formalism)
            data = np.genfromtxt(hfile, skip_header=5, names=True)

            # extract rotation periods at all ages
            periods = 2*np.pi/(86400*data['surf_avg_omega'])

            # mark the index of a desired age [Myr]
            ages_myr = data['star_age']/1e6
            age_diffs = abs(ages_myr - desired_age)
            age_index = np.where(age_diffs == min(age_diffs))[0][0]

            # then extract the rotation period at the desired age
            goss_prot[ct] = periods[age_index]
            ct += 1

    goss_bp_rp = mass_to_bp_rp(goss_mass)
    ax.plot(goss_bp_rp,goss_prot,'*',color="teal",ms=12)
    ax.axvline(0.9,color="k",linestyle=":")
    ax.axvline(0.78,color="k",linestyle=":")
    ax.set_title(model)
    plt.savefig(f"plots/periodmass_all_clean{clean_limit}_{model_f}.png")
    plt.close("all")

    ### Garraffo models
    ax = plot_periodcolor_histogram(clean_limit)
    # ax = plot_periodcolor(clean_limit)
    model = "Garraffo et al. (2018)"
    model_f = "Garr18"

    garr_mass = np.zeros(15*11)
    garr_prot = np.zeros(15*11)

    mass = ["06","07","08","09","10","11","12","13"]
    init_periods = [" 3"," 4"," 6"," 8","15","30","45","60","80","120"]
    desired_age = 40
    gar_dir = os.path.expanduser("~/Dropbox/Models/spintracks_Garraffo_18/")

    ct = 0
    for m in mass:
        for p in init_periods:
            garr_mass[ct] = float(m)/10
            gfile = os.path.join(gar_dir,f"age_p_tau_Mass{m}_P0={p}.txt")

            data = at.read(gfile,names=["age(Myr)","period(d)","unknown"])

            # mark the index of a desired age [Myr]
            ages_myr = data["age(Myr)"]
            age_diffs = abs(ages_myr - desired_age)
            age_index = np.where(age_diffs == min(age_diffs))[0][0]

            # then extract the rotation period at the desired age
            garr_prot[ct] = data["period(d)"][age_index]
            ct += 1

    garr_bp_rp = mass_to_bp_rp(garr_mass)
    ax.plot(garr_bp_rp,garr_prot,'*',color="teal",ms=12)
    ax.axvline(0.9,color="k",linestyle=":")
    ax.axvline(0.78,color="k",linestyle=":")
    ax.set_title(model)
    plt.savefig(f"plots/periodmass_all_clean{clean_limit}_{model_f}.png")
    plt.close("all")

def plot_data_boxes(fig, axes, ages,plot_name="",plot_title="",clean_limit=10,
                      which_plot="individual clusters",subset=False):

    clusters = ["IC_2391","Collinder_135","NGC_2451A","NGC_2547","IC_2602"]
    dates = ["2021-06-22","2021-06-18","2021-06-21","2021-06-21","2021-07-02"]

    dat = at.read("tab_all_stars.csv")
    bp_rp = dat["GAIAEDR3_BP"]-dat["GAIAEDR3_RP"]
    dat = Table(dat, masked=True, copy=False)
    # bp_rp_IC_2391, prot_IC_2391 = read_cluster_visual("IC_2391","2021-06-22",clean_limit,to_plot=False)
    # bp_rp_Collinder_135, prot_Collinder_135 = read_cluster_visual("Collinder_135","2021-06-18",clean_limit,to_plot=False)
    # bp_rp_NGC_2451A, prot_NGC_2451A = read_cluster_visual("NGC_2451A","2021-06-21",clean_limit,to_plot=False)
    # bp_rp_NGC_2547, prot_NGC_2547 = read_cluster_visual("NGC_2547","2021-06-21",clean_limit,to_plot=False)
    # bp_rp_IC_2602, prot_IC_2602 = read_cluster_visual("IC_2602","2021-07-02",clean_limit,to_plot=False)

    if subset:
        prot_raw = dat["Prot1"]
        prot_mask = dat["Prot1"].mask

        max_q=clean_limit
        include_blends=False
        include_lit=True

        qmask = (dat["Q1"]<=max_q)
        pmask = ((prot_mask==False) & (prot_raw>0))
        pmask = pmask.filled(fill_value=False)
        mmask = (dat["Mass"].mask==False)
        if include_blends==False:
            blmask = (dat["Bl?"]=="n") | (dat["Bl?"]=="m")
        else:
            blmask = np.ones(len(dat),bool)
        # If including literature values, have to replace them in the period mask
        if include_lit:
            litmask = (dat["LitPeriod"].mask==False) & (dat["LitPeriod"]>0)

            # I only want to use literature periods when I don't have a valid TESS period
            qmask_init = blmask & qmask & pmask
            use_lit = litmask & (qmask_init==False)
            prot_mask[use_lit] = False
            prot_raw.mask[use_lit] = False
            prot_raw[use_lit] = dat["LitPeriod"][use_lit]

            # Then the final selection should include literature or TESS periods, as appropriate
            lit_or_tess = qmask_init | litmask
            full_qmask = mmask & lit_or_tess
    #             print(np.where(full_qmask)[0])
        else:
            full_qmask = pmask & mmask & qmask & blmask

        dat["Prot_subset"] = prot_raw

        dat = dat[full_qmask]
    else:
        dat = dat[dat["Q1"]<=clean_limit]



    # y_age,y_perc,y_prot = young_stars_init()
    # eightmyr = np.ones_like(usco_perc)*8
    if subset:
        z_perc, z_solar, z_perc_indiv, z_solar_indiv = zams_percentiles_subset()
    else:
        z_perc, z_solar, z_perc_indiv, z_solar_indiv = zams_percentiles()

    clusters = ["IC_2391","Collinder_135","NGC_2451A","NGC_2547","IC_2602"]
    # prot = [prot_IC_2391, prot_Collinder_135, prot_NGC_2451A, prot_NGC_2547, prot_IC_2602]
    # solar = [solar_IC_2391, solar_Collinder_135, solar_NGC_2451A, solar_NGC_2547, solar_IC_2602]

    for i in range(3):
        for j in range(2):
            # if (i==2) and (j==1):
            #     break
            ax = axes[i,j]
            # ax.plot([45,45,45],z_perc,'o',mec="grey",mfc="none",
            #         mew=2,zorder=100)

            if which_plot=="individual clusters":
                # Plot the ZAMS clusters
                boxes = []
                list_ages = []
                for k in range(5):
                    boxes.append({"whislo": np.min(prot[k][solar[k]]),
                                  "q1": z_perc_indiv[k][0],
                                  "med": z_perc_indiv[k][1],
                                  "q3": z_perc_indiv[k][2],
                                  "whishi": np.max(prot[k][solar[k]]),
                                  "fliers": []
                                 })
                    list_ages.append(ages[clusters[k]])
                print(list_ages)
                ax.bxp(bxpstats=boxes,medianprops={"color":"k"},
                       positions=list_ages,widths=np.asarray(list_ages)*0.25,
                       manage_ticks=False,zorder=20)

                # # Plot the younger clusters
                # boxes = []
                # for k in range(4):
                #     boxes.append({"whislo": np.min(y_prot[k]),
                #                   "q1": y_perc[k][0],
                #                   "med": y_perc[k][1],
                #                   "q3": y_perc[k][2],
                #                   "whishi": np.max(y_prot[k]),
                #                   "fliers": []
                #                  })
                # ax.bxp(bxpstats=boxes,medianprops={"color":"k"},
                #        positions=y_age,widths=np.asarray(y_age)*0.25,
                #        manage_ticks=False,zorder=20)
            elif which_plot=="single age":
                z_prot = dat["Prot1"][(dat["Mass"]<=1.1) & (dat["Mass"]>=0.9)]
                # ax.boxplot(prot,sym="*",medianprops={"color":"k"},
                #            positions=[45],widths=[30],
                #            flierprops={"markersize":4},manage_ticks=False,zorder=20,
                #            whis=(10,75))
                #            # whis=(5,95))

                single_box = {"whislo": np.min(z_prot),
                              "q1": z_perc[0],
                              "med": z_perc[1],
                              "q3": z_perc[2],
                              "whishi": np.max(z_prot),
                              "fliers": []
                             }
                print(single_box)

                max_age = 55
                min_age = 30
                avg_age = 45#10**((np.log10(55)+np.log10(30))/2)
                width = (max_age-min_age)
                ax.bxp(bxpstats=[single_box],medianprops={"color":"k"},
                       positions=[avg_age],widths=[width],
                       manage_ticks=False,zorder=20)
                # # Plot the younger clusters
                # boxes = []
                # for k in range(4):
                #     boxes.append({"whislo": np.min(y_prot[k]),
                #                   "q1": y_perc[k][0],
                #                   "med": y_perc[k][1],
                #                   "q3": y_perc[k][2],
                #                   "whishi": np.max(y_prot[k]),
                #                   "fliers": []
                #                  })
                # ax.bxp(bxpstats=boxes,medianprops={"color":"k"},
                #        positions=y_age,widths=np.asarray(y_age)*0.25,
                #        manage_ticks=False,zorder=20)



def plot_model_tracks(ages,plot_name="",plot_title="",clean_limit=10,
                      which_plot="individual clusters",subset=False):

    fig, axes = plt.subplots(nrows=3,ncols=2,sharey=True,figsize=(8,10))
    # plt.suptitle(f"Solar mass, C{clean_limit}{plot_title}",y=0.93)
    plt.suptitle(f"Solar mass{plot_title}",y=0.93)

    plot_data_boxes(fig, axes, ages, plot_name,plot_title,clean_limit,
                    which_plot,subset)

    ########################################################################
    # Sean's models
    ax = axes[0,0]

    linestyle=["--","-",":"]
    for j,model in enumerate(["UpSco_Mattea2015","UpSco_Mattea2022","UpSco_ZeroTorque"]):
        mfiles = glob.glob(f"models/{model}*csv")
        for mfilename in mfiles:
            mfile = at.read(mfilename)
            pinit = float(mfilename.split("/")[-1].split("_")[3][1:5])

            ax.plot(mfile["age(Myr)"],mfile["per"],
                    color=mapper.to_rgba(pinit),linestyle=linestyle[j])
            if j==1:
                ax.plot([0,mfile["age(Myr)"][0]],[mfile["per"][0],mfile["per"][0]],
                        color=mapper.to_rgba(pinit),linestyle=linestyle[j])


    ax.set_xlim(1,1e4)
    ax.tick_params(labelbottom=False)
    ax.set_ylim(0.1,40)
    ax.set_xscale("log")
    ax.set_yscale("log")
    # ax.set_xlabel("Age (Myr)")
    ax.set_ylabel("Period (d)")

    ax.set_title(f"Matt et al. (2015)")#", Solar mass, C{clean_limit}{plot_title}")
    # plt.savefig(f"plots/periodtracks{plot_name}_Matt2015_clean{clean_limit}.png")


    ########################################################################
    # Cecilia's models

    ax = axes[0,1]

    gar_dir = "models/spintracks_Garraffo_18/"
    periods = [" 3"," 4"," 6"," 8","15","30","45","60","80","120"]
    period_label = ["0.3d","0.4d","0.6d","0.8d","1.5d","3.0d","4.5d","6.0d","8.0d","12.0d"]
    for i,p in enumerate(periods):
        gfile = os.path.join(gar_dir,f"age_p_tau_Mass10_P0={p}.txt")

        data = at.read(gfile,names=["age(Myr)","period(d)","unknown"])
        good = (data["age(Myr)"]>=1) & (data["age(Myr)"]<2800)

        color = mapper.to_rgba(float(period_label[i][:-1]))
        ax.plot(data["age(Myr)"][good],data["period(d)"][good],
                label=period_label[i],color=color)

    ax.set_xlim(1,1e4)
    ax.tick_params(labelbottom=False)
    ax.set_ylim(0.1,40)
    ax.set_xscale("log")
    ax.set_yscale("log")
    # ax.set_xlabel("Age (Myr)")
    # ax.set_ylabel("Period (d)")

    ax.set_title(f"Garraffo et al. (2018)")


    ########################################################################
    # Seth's models, Matt wind model
    ax = axes[1,0]

    periods = ["0.3d","0.4d","0.6d","0.8d","1.5d","12d","3d","4.5d","6d","8d"]
    for p in periods:
        hfile = f"models/MIST_M15_mbraking_feh_p0.00_histfiles/00100M_dir/LOGS/1.0M_history_prot{p}0_M15_d02d4_0.5hp_1.4d30K_0.22m_2.6p_3myr_alt3prot_henyey1.82mlta.data"
        # read in the history file (for a particular mass, initial period, and braking formalism)
        data = np.genfromtxt(hfile, skip_header=5, names=True)

        # extract rotation periods at all ages
        periods = 2*np.pi/(86400*data['surf_avg_omega'])

        color = mapper.to_rgba(float(p[:-1]))
        ax.plot(data["star_age"]/1e6,periods,label=p,color=color)

    ax.set_xlim(1,1e4)
    ax.tick_params(labelbottom=False)
    ax.set_xscale("log")
    # plt.ylim(0.1,40)
    # plt.yscale("log")
    # plt.xlabel("Age (Myr)")
    ax.set_ylabel("Period (d)")

    ax.set_title(f"Gossage et al. (2021), M15 wind")#", Solar Mass, C{clean_limit}{plot_title}")
    # plt.savefig(f"plots/periodtracks{plot_name}_Gossage2021_M15_clean{clean_limit}.png")

    ########################################################################
    # Seth's models, Garraffo wind model
    ax = axes[1,1]

    periods = ["0.3d","0.4d","0.6d","0.8d","1.5d","3d","4.5d","6d","8d","12d"]
    for p in periods:
        hfile = f"models/MIST_G18_mbraking_feh_p0.00_histfiles/00100M_dir/LOGS/1.0M_history_prot{p}0_G18_d02d4_0.5hp_1d_3.2d41c_0.5b_0.03a_3myr_alt3prot_henyey1.82mlta.data"
        # read in the history file (for a particular mass, initial period, and braking formalism)
        data = np.genfromtxt(hfile, skip_header=5, names=True)

        # extract rotation periods at all ages
        periods = 2*np.pi/(86400*data['surf_avg_omega'])

        color = mapper.to_rgba(float(p[:-1]))
        ax.plot(data["star_age"]/1e6,periods,label=p,color=color)

    ax.set_xlim(1,1e4)
    ax.set_xscale("log")
    ax.tick_params(labelbottom=False)
    # plt.ylim(0.1,40)
    # plt.xscale("log")
    # plt.yscale("log")
    # ax.set_xlabel("Age (Myr)")
    # plt.ylabel("Period (d)")

    ax.set_title(f"Gossage et al. (2021), G18 wind")#, Solar Mass, C{clean_limit}{plot_title}")
    # plt.savefig(f"plots/periodtracks{plot_name}_Gossage2021_G18_clean{clean_limit}.png")

    ########################################################################
    # Louis's models, 2019 version
    ax = axes[2,0]

    mass = "001p0"
    met = "13"
    periods = ["60","40","20"]
    for p in periods:
        model_dir = os.path.expanduser("~/Dropbox/Models/Amard2019/tracks/")
        tfile = f"M{mass}0Z{met}V{p}.dat"

        data = at.read(os.path.join(model_dir,tfile),data_start=3)


        color = mapper.to_rgba(data["Prot"][0])
        ax.plot(data["time"]/1e6,data["Prot"],label=f"V{p}",color=color)


    ax.set_xlim(1,1e4)
    ax.set_xscale("log")
    # plt.ylim(0.1,40)
    # plt.xscale("log")
    # plt.yscale("log")
    ax.set_xlabel("Age (Myr)")
    ax.set_ylabel("Period (d)")

    ax.set_title(f"Amard et al. (2019)")#", Solar Mass, C{clean_limit}{plot_title}")
    # plt.savefig(f"plots/periodtracks{plot_name}_Amard2019_clean{clean_limit}.png")

    ########################################################################
    # Spada & Lanzafame 2020
    ax = axes[2,1]

    sfile = at.read("models/spadalanzafame2020_tableA1.dat",delimiter="\t")
    i = np.where(sfile["Mass"]==1.00)
    ages_gyr = np.array([0.10,0.12,0.15,0.20,0.22,0.25,0.30,0.40,0.50,
                         0.60,0.70,1.00,1.50,2.00,2.50,4.00,4.57])
    solar_slow = np.array([4.33,4.45,4.64,4.98,5.12,5.33,5.70,6.45,7.21,7.95,
                           8.69,10.79,13.91,16.63,19.06,25.23,27.32])
    color = mapper.to_rgba(float(sfile["0.10Gyr"][i]))
    ax.plot(ages_gyr*1e3,solar_slow,'-',color=color)

    ax.set_xlim(1,1e4)
    ax.tick_params(labelleft=False)
    ax.set_ylim(0.1,40)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Age (Myr)")
    # ax.set_ylabel("Period (d)")

    ax.set_title(f"Spada & Lanzafame (2020)")#", Solar mass, C{clean_limit}{plot_title}")
    # plt.savefig(f"plots/periodtracks{plot_name}_Matt2015_clean{clean_limit}.png")

    plt.savefig(f"plots/periodtracks{plot_name}_page_clean{clean_limit}.png",
                bbox_inches="tight",dpi=600)
    plt.close()


def plot_results():
    # Cantat-Gaudin+2020
    cg20_ages = {"IC_2391": 28.8,
                 "IC_2602": 36.3,
                 "NGC_2547": 32.4,
                 "NGC_2451A": 35.5,
                 "Collinder_135": 26.3}
    # Kharchenko+2013
    khar_ages = {"IC_2391": 112.2,
                 "IC_2602": 221.3,
                 "NGC_2547": 77.6,
                 "NGC_2451A": 57.5,
                 "Collinder_135": 39.8}
    # Ghoza 2012
    g12_ages = {"IC_2391": 46,
                 "IC_2602": 32,
                 "NGC_2547": 38,
                 "NGC_2451A": 60,
                 "Collinder_135": 26}
    for clean_limit in [0,1]:
        plot_model_tracks(cg20_ages,plot_name="_medages",plot_title=", Combined Clusters",
                          clean_limit=clean_limit,which_plot="single age")
        plot_model_tracks(cg20_ages,plot_name="_medages_subset",plot_title=", Combined Clusters, Subset",
                          clean_limit=clean_limit,which_plot="single age",subset=True)
        # plot_model_tracks(cg20_ages,plot_name="_CG20ages",plot_title=", Cantat-Gaudin et al. (2020) Ages",
        # clean_limit=clean_limit)
        # plot_model_tracks(khar_ages,plot_name="_Kharages",plot_title=", Kharchenko et al. Ages",
        # clean_limit=clean_limit)
        # plot_model_tracks(g12_ages,plot_name="_G12ages",plot_title=", Ghoza et al. (2012) Ages",
        # clean_limit=clean_limit)

if __name__=="__main__":

    plot_results()

    # plot_periodcolor_models(clean_limit=0)
    # ax = plot_periodcolor_histogram(clean_limit=0,to_plot_indiv=False)
    # plt.savefig("plots/periodmass_histogram.png",bbox_inches="tight",dpi=600)