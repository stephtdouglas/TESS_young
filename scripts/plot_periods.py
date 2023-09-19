import os, sys
import glob
import itertools

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.cm as cm
import astropy.io.ascii as at
from astropy.table import join,vstack,Table
from scipy import stats
from scipy.interpolate import interp1d

from analyze_cluster_output import read_cluster_visual
from model_data_percentiles import zams_percentiles
import get_colors
norm, mapper, cmap2, colors, shapes = get_colors.get_colors()

plt.style.use('./paper.mplstyle')
PAPER_DIR = os.path.expanduser("~/my_papers/TESS_young/")

def plot_periodmass(dat, clusters):

    mass = dat["Mass"]
    mass_err = dat["Mass_err"]

    fig2 = plt.figure()
    ax2 = plt.subplot(111)
    ax2.set_ylim(0.1,50)
    ax2.set_xlim(1.3,0.1)
    ax2.set_yscale("log")

    ax2.set_xlabel(r"Mass (M$_{\odot}$)",fontsize=16)
    ax2.set_ylabel("Period (d)",fontsize=16)
    for i, cluster in enumerate(clusters):
        fig = plt.figure()
        ax = plt.subplot(111)
        cloc = (dat["Cluster"]==cluster) & (dat["Q1"]==0)

        ax.errorbar(mass[cloc],dat["Prot1"][cloc],xerr=mass_err[cloc],
                    marker=None,linewidth=0,elinewidth=1,color=colors[cluster])
        ax.plot(mass[cloc],dat["Prot1"][cloc],shapes[cluster],label=cluster,
                 ms=5,color=colors[cluster])
        ax2.plot(mass[cloc],dat["Prot1"][cloc],shapes[cluster],label=cluster,
                 ms=5,color=colors[cluster])
        ax.legend(loc=1)

        ax.set_ylim(0.1,50)
        ax.set_xlim(1.3,0.1)
        ax.set_yscale("log")

        ax.set_xlabel(r"Mass (M$_{\odot}$)",fontsize=16)
        ax.set_ylabel("Period (d)",fontsize=16)

        ax.tick_params(labelsize=14)

        fig.savefig(f"plots/periodmass_{cluster}.png",
                    bbox_inches="tight",dpi=600)


    ax2.legend(loc=1)

    ax2.set_ylim(0.1,50)
    ax2.set_xlim(1.3,0.1)
    ax2.set_yscale("log")

    ax2.set_xlabel(r"Mass (M$_{\odot}$)",fontsize=16)
    ax2.set_ylabel("Period (d)",fontsize=16)

    ax2.tick_params(labelsize=14)
    fig2.savefig(f"plots/periodmass_all_color.png",
                 bbox_inches="tight",dpi=600)



    plt.close("all")
    # ax.axhline(12,linestyle="--",color="tab:grey")

    good = dat["Q1"]<=1
    blend = dat["Bl?"]=="y"
    plt.plot(mass[good & (blend==False)],dat["Prot1"][good & (blend==False)],
             'k.',mfc="k",label="Good periods, no faint blends",zorder=1)#,alpha=0.5)
    plt.plot(mass[good],dat["Prot1"][good],'.',mfc="DarkGrey",mec="DarkGrey",
             label="Good periods, neighbor within 30 arcsec",zorder=0)#,alpha=0.5)
    plt.legend()
    plt.errorbar(mass[good],dat["Prot1"][good],xerr=mass_err[good],
                marker=None,linewidth=0,elinewidth=1,color="DarkGrey",zorder=-10)
    plt.ylim(0.1,50)
    plt.xlim(1.3,0.1)
    plt.yscale("log")

    plt.xlabel(r"Mass (M$_{\odot}$)",fontsize=16)
    plt.ylabel("Period (d)",fontsize=16)

    ax = plt.gca()
    ax.tick_params(labelsize=14)

    plt.savefig("plots/periodmass_all.png",
                bbox_inches="tight",dpi=600)

    lit_loc = (dat["LitPeriod"].mask==False) & (dat["LitPeriod"]>0)
    ax.plot(mass[lit_loc],dat["LitPeriod"][lit_loc],'.',color=colors["NGC_2547"],
            label="Literature")

    # for i, cluster in enumerate(clusters):
    #     cloc = (dat["Cluster"]==cluster) & (dat["Q1"]==0)
    #
    #     ax.errorbar(mass[cloc],dat["Prot1"][cloc],xerr=mass_err[cloc],
    #                 marker=None,linewidth=0,elinewidth=1,color=colors[cluster])
    #     ax.plot(mass[cloc],dat["Prot1"][cloc],shapes[cluster],label=cluster,
    #              ms=5,color=colors[cluster])
    #     ax2.plot(mass[cloc],dat["Prot1"][cloc],shapes[cluster],label=cluster,
    #              ms=5,color=colors[cluster])
    ax.legend(loc=1)
    plt.savefig("plots/periodmass_all_lit.png",
                bbox_inches="tight",dpi=600)


    plt.close()

    plt.errorbar(mass[good],dat["Prot1"][good],xerr=mass_err[good],
                marker=None,linewidth=0,elinewidth=1,color="Grey")
    plt.plot(mass[good],dat["Prot1"][good],'.',alpha=0.5,color="DarkGrey")
    plt.ylim(0.1,50)
    plt.xlim(1.3,0.1)
    plt.yscale("log")

    plt.xlabel(r"Mass (M$_{\odot}$)",fontsize=16)
    plt.ylabel("Period (d)",fontsize=16)

    ax = plt.gca()
    ax.tick_params(labelsize=14)


    z_perc, z_solar, z_perc_indiv, z_solar_indiv, = zams_percentiles()
    z_prot = dat["Prot1"][(dat["Mass"]<=1.1) & (dat["Mass"]>=0.9) & (dat["Q1"]==0)]
    print(np.nanmedian(z_prot),"median")
    single_box = {"whislo": np.min(z_prot),
                  "q1": z_perc[0],
                  "med": z_perc[1],
                  "q3": z_perc[2],
                  "whishi": np.max(z_prot),
                  "fliers": [],
                  "color":mapper.to_rgba(20)
                 }

    colorprop = {"color":mapper.to_rgba(1),"linewidth":3}
    ax.bxp(bxpstats=[single_box],medianprops=colorprop,
           positions=[1.0],widths=[0.2],boxprops=colorprop,
           manage_ticks=False,zorder=20,whiskerprops=colorprop,capprops=colorprop)
    plt.savefig("plots/periodmass_all_bxp.png",bbox_inches="tight")


    plt.close()


def plot_periodcolor(dat, clusters, single_figure=False):

    bp_rp = dat["GAIAEDR3_BP"]-dat["GAIAEDR3_RP"]

    if single_figure:
        fig, axes = plt.subplots(3,2,figsize=(12,12),sharex=True,sharey=True)
        fl_axes = np.asarray(axes).flatten()
        ax_all = fl_axes[-1]

    for i, cluster in enumerate(clusters):

        if single_figure is False:
            plt.figure()
            ax = plt.subplot()
        else:
            ax = fl_axes[i]

        if (single_figure is False) or (i==0):
            ax.set_ylim(0.1,50)
            ax.set_xlim(0.5,3.5)
            ax.set_yscale("log")

        if (single_figure is False) or (i % 2 == 0):
            ax.set_ylabel("Period (d)")#,fontsize=16)
        if (single_figure is False) or (i>=4):
            ax.set_xlabel(r"G$_{BP}$ - G$_{RP}$ (EDR3)")#,fontsize=16)
            # ax.tick_params(labelsize=14)

        # ax.axhline(12,linestyle="--",color="tab:grey")
        ax.set_title(cluster.replace("_"," "))

        cloc = (dat["Cluster"]==cluster) & (dat["Q1"]==0) & (dat["to_plot"]==1)
        cloc1 = (dat["Cluster"]==cluster) & (dat["Q1"]==1) & (dat["to_plot"]==1)

        ax.plot(bp_rp[cloc],dat["Prot1"][cloc],shapes[cluster],
                 label="High Q TESS",
                 ms=5,color=colors[cluster],zorder=6)
        ax.plot(bp_rp[cloc1],dat["Prot1"][cloc1],shapes[cluster],
                 ms=5,color=colors[cluster],mfc="none",zorder=0,
                 label="Low Q TESS")
        ax.legend(loc=1,fontsize=12)

        if single_figure is False:
            plt.savefig(f"plots/periodcolor_{cluster}.png",
                        bbox_inches="tight",dpi=600)

        cl_lit = ((dat["Cluster"]==cluster) & (dat["LitPeriod"].mask==False)
                & (dat["LitPeriod"]>0) & (dat["to_plot"]==1))
        if np.any(cl_lit):
            ax.plot(bp_rp[cl_lit],dat["LitPeriod"][cl_lit],shapes[cluster],
                    mec="DarkGrey",mfc="none",ms=9,label="Literature")
            ax.legend(loc=1,fontsize=12,ncol=2)
            if single_figure is False:
                plt.savefig(f"plots/periodcolor_{cluster}_lit.png",
                            bbox_inches="tight",dpi=600)
        if single_figure is False:
            plt.close()


    if single_figure is False:
        plt.figure()
        ax = plt.subplot()
        ax.set_ylim(0.1,50)
        ax.set_xlim(0.5,3.5)
        ax.set_yscale("log")
        ax.set_xlabel(r"G$_{BP}$ - G$_{RP}$ (EDR3)",fontsize=16)
        ax.set_ylabel("Period (d)",fontsize=16)
    else:
        ax = fl_axes[-1]
        ax.set_xlabel(r"G$_{BP}$ - G$_{RP}$ (EDR3)")

    no_blend = dat["Bl?"]=="n"
    maybe_blend = dat["Bl?"]=="m"
    good = (dat["Q1"]==0)
    ax.plot(bp_rp[good],dat["Prot1"][good],'k.',alpha=0.1,label="High Q TESS")
    ax.plot(bp_rp[good & maybe_blend],dat["Prot1"][good & maybe_blend],'k.',
            mfc="none",label="Possible blend")
    ax.plot(bp_rp[good & no_blend],dat["Prot1"][good & no_blend],'k.',
            label="Not blended")

    lit = dat["LitPeriod"].mask==False
    bad_tess = (dat["Q1"]>=1) | (dat["Bl?"]=="y")
    lit_replace = lit &  bad_tess & (dat["to_plot"]==1)
    ax.plot(bp_rp[lit_replace],dat["LitPeriod"][lit_replace],'*',ms=7,
            label="Lit replaces TESS",color=cmap2(5),mec=cmap2(2))
    ax.legend(fontsize=12,ncol=2)
    ax.set_title("Combined sample")

    if single_figure is False:
        plt.savefig("plots/periodcolor_all.png",bbox_inches="tight")
    else:
        plt.savefig("plots/periocdolor_panel.png",bbox_inches="tight")
        plt.savefig(os.path.join(PAPER_DIR,"fig_periodcolor_panel.pdf"),
                    bbox_inches="tight")
    plt.close()


def plot_periodcolor_membership(dat):

    bp_rp = dat["GAIAEDR3_BP"]-dat["GAIAEDR3_RP"]

    fig, axes = plt.subplots(2,2,figsize=(12,8),sharex=True,sharey=True)
    fl_axes = np.asarray(axes).flatten()
    ax_all = fl_axes[-1]

    # Membership filter for plotting
    hdb_memb = (dat["MemBool"]==1)  #& (dat["MemBool"].mask==False)

    # Jackson+2020 table 4/Section 4 indictes that 0.9 is the membership cutoff
    ges_memb = ((dat["GES_MemProb"]>=0.9) & (dat["GES_MemProb"]<=1))
                # & (dat["GES_MemProb"].mask==False))

    can_memb = ((dat["CG_MemProb"]>=0.7) & (dat["CG_MemProb"]<=1))
                # & (dat["CG_MemProb"].mask==False))

    catalogs = ["HDBScan","Gaia-ESO Survey","Cantat-Gaudin","Combined"]
    memb = [hdb_memb, ges_memb, can_memb, dat["to_plot"]==1]

    for i, catalog in enumerate(catalogs):

        ax = fl_axes[i]

        if (i==0):
            ax.set_ylim(0.1,50)
            ax.set_xlim(0.5,3.5)
            ax.set_yscale("log")

        if (i % 2 == 0):
            ax.set_ylabel("Period (d)")#,fontsize=16)
        if (i>=2):
            ax.set_xlabel(r"G$_{BP}$ - G$_{RP}$ (EDR3)")#,fontsize=16)
            # ax.tick_params(labelsize=14)

        # ax.axhline(12,linestyle="--",color="tab:grey")
        ax.set_title(catalog)

        cloc = memb[i] & (dat["Q1"]==0)
        cloc1 = memb[i] & (dat["Q1"]==1)

        ax.plot(bp_rp[cloc],dat["Prot1"][cloc],'o',
                 label="High Q TESS",ms=5,color="k",zorder=6)
        ax.plot(bp_rp[cloc1],dat["Prot1"][cloc1],'o',
                 label="Low Q TESS",ms=5,color="k",zorder=6,mfc="none")
        ax.legend(loc=1,fontsize=12)


    plt.savefig("plots/periocdolor_membership.png",bbox_inches="tight")
    # plt.savefig(os.path.join(PAPER_DIR,"fig_periodcolor_panel.pdf"),
    #             bbox_inches="tight")
    plt.close()



def plot_periodcolor_ruwe(clean_limit=0,to_plot_indiv=False):
    cat_IC_2391 = read_cluster_visual("IC_2391","2021-06-22",clean_limit,to_plot=to_plot_indiv,return_periodcolor=False)
    cat_Collinder_135 = read_cluster_visual("Collinder_135","2021-06-18",clean_limit,to_plot=to_plot_indiv,return_periodcolor=False)
    cat_NGC_2451A = read_cluster_visual("NGC_2451A","2021-06-21",clean_limit,to_plot=to_plot_indiv,return_periodcolor=False)
    cat_NGC_2547 = read_cluster_visual("NGC_2547","2021-06-21",clean_limit,to_plot=to_plot_indiv,return_periodcolor=False)
    cat_IC_2602 = read_cluster_visual("IC_2602","2021-07-02",clean_limit,to_plot=to_plot_indiv,return_periodcolor=False)
    clusters = ["IC_2391","Collinder_135","NGC_2451A","NGC_2547","IC_2602"]
    cats = [cat_IC_2391, cat_Collinder_135, cat_NGC_2451A, cat_NGC_2547, cat_IC_2602]

    fig1 = plt.figure()
    ax = plt.subplot(111)


    for i, match in enumerate(cats):
        # print(match.dtype)
        bp_rp = match["GAIAEDR3_BP"] - match["GAIAEDR3_RP"]

        fig2 = plt.figure()
        ax2 = plt.subplot(111)


        if clean_limit is not None:
            if clean_limit==0:
                all_possible = (match["final_Q"]==0)
            else:
                all_possible = (match["final_Q"]==0) | (match["final_Q"]==1)
        else:
            all_possible = (match["final_Q"]==0) | (match["final_Q"]==1)

        sc = ax.scatter(bp_rp[all_possible],match["final_period"][all_possible],
                    c=match["GAIAEDR3_RUWE"][all_possible],vmin=0.6,vmax=2,alpha=0.75)
        sc = ax2.scatter(bp_rp[all_possible],match["final_period"][all_possible],
                    c=match["GAIAEDR3_RUWE"][all_possible],vmin=0.6,vmax=2,alpha=0.75)
        plt.colorbar(sc,label="Gaia RUWE - EDR3",ax=ax2)

        ax2.set_ylim(0.1,50)
        ax2.set_xlim(0.5,3.5)
        ax2.set_yscale("log")

        ax2.set_xlabel(r"G$_{BP}$ - G$_{RP}$")
        ax2.set_ylabel("Period (d)")
        ax2.axhline(12,linestyle="--",color="tab:grey")

        ax2.set_title(clusters[i])
        fig2.savefig(f"plots/periodcolor_{clusters[i]}_ruwe.png",bbox_inches="tight",dpi=300)
        plt.close(fig2)


    # TODO: add gemini observations and detections

    plt.colorbar(sc,label="Gaia RUWE - EDR3",ax=ax)

    ax.set_ylim(0.1,50)
    ax.set_xlim(0.5,3.5)
    ax.set_yscale("log")

    ax.set_xlabel(r"G$_{BP}$ - G$_{RP}$")
    ax.set_ylabel("Period (d)")

    ax.axhline(12,linestyle="--",color="tab:grey")

    return ax

def plot_periodcolor_histogram(clean_limit=10,to_plot_indiv=False,ax=None):
    bp_rp_IC_2391, prot_IC_2391 = read_cluster_visual("IC_2391","2021-06-22",clean_limit,to_plot=to_plot_indiv)
    bp_rp_Collinder_135, prot_Collinder_135 = read_cluster_visual("Collinder_135","2021-06-18",clean_limit,to_plot=to_plot_indiv)
    bp_rp_NGC_2451A, prot_NGC_2451A = read_cluster_visual("NGC_2451A","2021-06-21",clean_limit,to_plot=to_plot_indiv)
    bp_rp_NGC_2547, prot_NGC_2547 = read_cluster_visual("NGC_2547","2021-06-21",clean_limit,to_plot=to_plot_indiv)
    bp_rp_IC_2602, prot_IC_2602 = read_cluster_visual("IC_2602","2021-07-02",clean_limit,to_plot=to_plot_indiv)

    bp_rp = np.concatenate([bp_rp_IC_2391,bp_rp_Collinder_135,bp_rp_NGC_2451A,bp_rp_NGC_2547,bp_rp_IC_2602])
    prot = np.concatenate([prot_IC_2391,prot_Collinder_135,prot_NGC_2451A,prot_NGC_2547,prot_IC_2602])

    if ax is None:
        plt.figure()
        ax = plt.subplot(111)
    counts, xedges, yedges, image = ax.hist2d(bp_rp,prot,bins=[np.linspace(0.5,3.5,22),
                                          np.logspace(-1,np.log10(50),20)],cmap="Greys")
    xx = xedges[:-1] + (xedges[1:] - xedges[:-1])/2
    yy = yedges[:-1] + (yedges[1:] - yedges[:-1])/2
    cs = ax.contour(xx,yy,counts.transpose(),linewidths=2,levels=[1,3,10,13,20],colors="k",antialiased=True)
    # plt.clabel(cs,fmt="%1.0f",inline=False)

    ax.set_ylim(0.1,50)
    ax.set_xlim(0.5,3.5)
    ax.set_yscale("log")

    ax.set_xlabel(r"G$_{BP}$ - G$_{RP}$")
    ax.set_ylabel("Period (d)")

    ax.axhline(12,linestyle="--",color="tab:grey")

    return ax



if __name__=="__main__":

    clusters = ["Collinder_135","NGC_2451A","NGC_2547","IC_2391","IC_2602"]
    dat = at.read("tab_all_stars.csv")

    plot_periodcolor(dat, clusters, single_figure=True)
    plot_periodcolor_membership(dat)
    # plot_periodmass(dat, clusters)