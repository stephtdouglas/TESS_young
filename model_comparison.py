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
from model_data_percentiles import young_stars_init, zams_percentiles

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
                      which_plot="individual clusters"):

    bp_rp_IC_2391, prot_IC_2391 = read_cluster_visual("IC_2391","2021-06-22",clean_limit,to_plot=False)
    bp_rp_Collinder_135, prot_Collinder_135 = read_cluster_visual("Collinder_135","2021-06-18",clean_limit,to_plot=False)
    bp_rp_NGC_2451A, prot_NGC_2451A = read_cluster_visual("NGC_2451A","2021-06-21",clean_limit,to_plot=False)
    bp_rp_NGC_2547, prot_NGC_2547 = read_cluster_visual("NGC_2547","2021-06-21",clean_limit,to_plot=False)
    bp_rp_IC_2602, prot_IC_2602 = read_cluster_visual("IC_2602","2021-07-02",clean_limit,to_plot=False)

    solar_IC_2391 = id_solar(bp_rp_IC_2391)
    solar_Collinder_135 = id_solar(bp_rp_Collinder_135)
    solar_NGC_2451A = id_solar(bp_rp_NGC_2451A)
    solar_NGC_2547 = id_solar(bp_rp_NGC_2547)
    solar_IC_2602 = id_solar(bp_rp_IC_2602)

    y_age,y_perc,y_prot = young_stars_init()
    # eightmyr = np.ones_like(usco_perc)*8
    z_perc, z_solar = zams_percentiles()

    for i in range(3):
        for j in range(2):
            # if (i==2) and (j==1):
            #     break
            ax = axes[i,j]
            ax.plot([45,45,45],z_perc,'o',mec="grey",mfc="none",
                    mew=2,zorder=100)

            if which_plot=="individual stars":
                ax.plot(np.ones_like(prot_IC_2391[solar_IC_2391])*ages["IC_2391"],
                         prot_IC_2391[solar_IC_2391],"o",label="IC_2391",
                         color="grey",zorder=20)
                ax.plot(np.ones_like(prot_Collinder_135[solar_Collinder_135])*ages["Collinder_135"],
                         prot_Collinder_135[solar_Collinder_135],"s",
                         label="Collinder_135",color="grey",zorder=20)
                ax.plot(np.ones_like(prot_NGC_2451A[solar_NGC_2451A])*ages["NGC_2451A"],
                         prot_NGC_2451A[solar_NGC_2451A],"^",
                         label="NGC_2451A",color="grey",zorder=20)
                ax.plot(np.ones_like(prot_NGC_2547[solar_NGC_2547])*ages["NGC_2547"],
                         prot_NGC_2547[solar_NGC_2547],"v",label="NGC_2547",
                         color="grey",zorder=20)
                ax.plot(np.ones_like(prot_IC_2602[solar_IC_2602])*ages["IC_2602"],
                         prot_IC_2602[solar_IC_2602],"d",label="IC_2602",
                         color="grey",zorder=20)
                for k in range(4):
                    ax.plot(np.ones_like(y_prot[k]),y_prot[k],'*',color="grey")

            elif which_plot=="individual clusters":
                ax.boxplot(prot_IC_2391[solar_IC_2391],sym="o",medianprops={"color":"grey"},
                           positions=[ages["IC_2391"]],widths=[ages["IC_2391"]*0.25],
                           flierprops={"markersize":4},manage_ticks=False,zorder=20,
                           whis=(10,75))
                           # whis=(5,95))
                ax.boxplot(prot_Collinder_135[solar_Collinder_135],sym="s",medianprops={"color":"grey"},
                           positions=[ages["Collinder_135"]],widths=[ages["Collinder_135"]*0.25],
                           flierprops={"markersize":4},manage_ticks=False,zorder=20,
                           whis=(10,75))
                           # whis=(5,95))
                ax.boxplot(prot_NGC_2451A[solar_NGC_2451A],sym="^",medianprops={"color":"grey"},
                           positions=[ages["NGC_2451A"]],widths=[ages["NGC_2451A"]*0.25],
                           flierprops={"markersize":4},manage_ticks=False,zorder=20,
                           whis=(10,75))
                           # whis=(5,95))
                ax.boxplot(prot_NGC_2547[solar_NGC_2547],sym="v",medianprops={"color":"grey"},
                           positions=[ages["NGC_2547"]],widths=[ages["NGC_2547"]*0.25],
                           flierprops={"markersize":4},manage_ticks=False,zorder=20,
                           whis=(10,75))
                           # whis=(5,95))
                ax.boxplot(prot_IC_2602[solar_IC_2602],sym="d",medianprops={"color":"grey"},
                           positions=[ages["IC_2602"]],widths=[ages["IC_2602"]*0.25],
                           flierprops={"markersize":4},manage_ticks=False,zorder=20,
                           whis=(10,75))
                           # whis=(5,95))
                for k in range(4):
                    ax.plot(np.ones(3)*y_age[k],y_perc[k],'k*')
                    # ax.boxplot(y_prot[k],sym="*",medianprops={"color":"grey"},
                    #            positions=[y_age[k]],widths=[y_age[k]*0.25],
                    #            flierprops={"markersize":4},manage_ticks=False,zorder=20,
                    #            whis=(10,75))
            elif which_plot=="single age":
                prot = np.concatenate((prot_IC_2391[solar_IC_2391],
                                      prot_Collinder_135[solar_Collinder_135],
                                      prot_NGC_2451A[solar_NGC_2451A],
                                      prot_NGC_2547[solar_NGC_2547],
                                      prot_IC_2602[solar_IC_2602]))
                ax.boxplot(prot,sym="*",medianprops={"color":"k"},
                           positions=[45],widths=[30],
                           flierprops={"markersize":4},manage_ticks=False,zorder=20,
                           whis=(10,75))
                           # whis=(5,95))
                for k in range(4):
                    ax.plot(np.ones(3)*y_age[k],y_perc[k],'k*')
                    # ax.boxplot(y_prot[k],sym="*",medianprops={"color":"grey"},
                    #            positions=[y_age[k]],widths=[y_age[k]*0.25],
                    #            flierprops={"markersize":4},manage_ticks=False,zorder=20,
                    #            whis=(10,75))



def plot_model_tracks(ages,plot_name="",plot_title="",clean_limit=10,
                      which_plot="individual clusters"):

    fig, axes = plt.subplots(nrows=3,ncols=2,sharey=True,figsize=(8,10))
    # plt.suptitle(f"Solar mass, C{clean_limit}{plot_title}",y=0.93)
    plt.suptitle(f"Solar mass{plot_title}",y=0.93)

    plot_data_boxes(fig, axes, ages, plot_name,plot_title,clean_limit,
                          which_plot)

    ########################################################################
    # Sean's models
    ax = axes[0,0]

    mfile = at.read("models/spintracks_Matt_ea_15/spintrack10_Mea15.txt",
                    names=["age(yr)","P0_0.7","P0_5","P0_18"])

    ax.plot(mfile["age(yr)"]/1e6,mfile["P0_0.7"],color=mapper.to_rgba(0.7))
    ax.plot(mfile["age(yr)"]/1e6,mfile["P0_5"],color=mapper.to_rgba(5))
    ax.plot(mfile["age(yr)"]/1e6,mfile["P0_18"],color=mapper.to_rgba(18))

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
    # Sean's new models
    ax = axes[2,1]

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
    ax.tick_params(labelleft=False)
    ax.set_ylim(0.1,40)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Age (Myr)")
    # ax.set_ylabel("Period (d)")

    ax.set_title(f"Matt et al. (New)")#", Solar mass, C{clean_limit}{plot_title}")
    # plt.savefig(f"plots/periodtracks{plot_name}_Matt2015_clean{clean_limit}.png")

    plt.savefig(f"plots/periodtracks{plot_name}_page_clean{clean_limit}.png",
                bbox_inches="tight")
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
        # plot_all(clean_limit)
        plot_model_tracks(cg20_ages,plot_name="_CG20ages",plot_title=", Cantat-Gaudin et al. (2020) Ages",
        clean_limit=clean_limit)
        plot_model_tracks(khar_ages,plot_name="_Kharages",plot_title=", Kharchenko et al. Ages",
        clean_limit=clean_limit)
        plot_model_tracks(g12_ages,plot_name="_G12ages",plot_title=", Ghoza et al. (2012) Ages",
        clean_limit=clean_limit)
        break

if __name__=="__main__":

    # # write_results()
    plot_results()
    # # compare_visual_results(cluster="NGC_2451A",date = "2021-06-21")
    # # compare_visual_results(cluster="NGC_2547",date = "2021-06-21")
    # # compare_visual_results(cluster="IC_2391",date = "2021-06-22")
    #
    # plot_all(clean_limit=0)
    # plot_periodcolor_models(clean_limit=0)
    # ax = plot_periodcolor_histogram(clean_limit=0,to_plot_indiv=False)
    # plt.savefig("plots/periodmass_histogram.png",bbox_inches="tight")
