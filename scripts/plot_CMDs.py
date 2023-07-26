import os, pathlib
import numpy as np
from numpy.polynomial.polynomial import Polynomial
import astropy.io.ascii as at
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from astropy.stats import sigma_clipped_stats


from tess_young.get_const import *
import tess_young
_DIR = pathlib.Path(tess_young.__file__).resolve().parent.parent
plt.style.use(os.path.join(_DIR,'paper.mplstyle'))

def calc_gaia_extinction(bp_rp0,A0,band="G"):
    """
    Calculate A_band, given an input bp-rp color and A_V
    """
    if band=="G":
        coeffs = [0.9761, -0.1704, 0.0086, 0.0011, -0.0438, 0.0013, 0.0099]
        # poly_color = Polynomial(coeffs[:4])
        # A_G = poly_color(bp_rp0)+coeffs[4]*A0+coeffs[5]*A0**2+coeffs[6]*bp_rp0*A0
        # return A_G
    elif band=="BP":
        coeffs = [1.1517, -0.0871, -0.0333, 0.0173, -0.0230, 0.0006, 0.0043]
        # poly_color = Polynomial(coeffs[:4])
        # A_BP = poly_color(bp_rp0)+coeffs[4]*A0+coeffs[5]*A0**2+coeffs[6]*bp_rp0*A0
        # return A_BP
    elif band=="RP":
        coeffs = [0.6104, -0.0170, -0.0026, -0.0017, -0.0078, 0.00005, 0.0006]

    poly_color = Polynomial(coeffs[:4])
    k_band = poly_color(bp_rp0)+coeffs[4]*A0+coeffs[5]*A0**2+coeffs[6]*bp_rp0*A0
    A_band = k_band * A0
    return k_band, A_band

def calc_gaia_reddening(bp_rp0,A0):
    """
    Calculate E(BP-RP), given an input bp-rp color and A_V
    """
    k_BP, A_BP = calc_gaia_extinction(bp_rp0,A0,"BP")
    k_RP, A_RP = calc_gaia_extinction(bp_rp0,A0,"RP")

    E_BP_RP = A_BP - A_RP
    return E_BP_RP

def plot_gaia_cmd(ax,cluster,gdr="EDR3",**kwargs):
    cat = at.read(os.path.join(_DIR,"tab_all_stars.csv"),delimiter=",")
    cat = cat[cat["Cluster"]==cluster]

    bp_rp = cat[f"GAIA{gdr}_BP"] - cat[f"GAIA{gdr}_RP"]
    dist = 1000/cat["GAIAEDR3_PARALLAX_CORRECTED"]
    abs_G = cat[f"GAIA{gdr}_G"] - 5*np.log10(dist) + 5

    good_memb = cat["to_plot"]==1

    ax.plot(bp_rp[good_memb],abs_G[good_memb],shapes[cluster],color=colors[cluster],label=cluster,**kwargs)


def plot_membership_cmd(ax,cluster,gdr="EDR3",**kwargs):
    cat = at.read(os.path.join(_DIR,"tab_all_stars.csv"),delimiter=",")
    cat = cat[cat["Cluster"]==cluster]

    bp_rp = cat[f"GAIA{gdr}_BP"] - cat[f"GAIA{gdr}_RP"]
    dist = 1000/cat["GAIAEDR3_PARALLAX_CORRECTED"]
    abs_G = cat[f"GAIA{gdr}_G"] - 5*np.log10(dist) + 5

    good_memb = cat["to_plot"]==1

    ax.plot(bp_rp,abs_G,"*",color="grey",label=cluster,ms=2,alpha=0.6)

    ax.plot(bp_rp[good_memb],abs_G[good_memb],shapes[cluster],color=colors[cluster],label=cluster,**kwargs)



def plot_all_cmd():
    plt.figure()
    ax1 = plt.subplot(111)
    ax2 = ax1.inset_axes([0.5, 0.5, 0.47, 0.47])
    for cluster in clusters:
        plot_gaia_cmd(ax1,cluster,alpha=0.75,ms=4)
        plot_gaia_cmd(ax2,cluster,alpha=0.95,ms=4)
    plt.legend(loc=3)
    ax1.set_xlabel(r"G$_{BP}$ - G$_{RP}$ (DR3)")
    ax1.set_ylabel(f"M$_G$ (DR3)")

    ax1.set_ylim(13,-5)
    ax1.set_xlim(-0.6,4)

    ax2.set_xlim(-0.3,0.45)
    ax2.set_ylim(2.5,-3)
    plt.savefig(os.path.join(_DIR,"plots/CMD_all.png"),bbox_inches="tight")
    plt.savefig(os.path.join(PAPER_DIR,"fig_CMD_all.pdf"),bbox_inches="tight")


def plot_panel_cmd():


    # Set up axes: 2 rows 3 cols 
    fig, axes = plt.subplots(nrows=2,ncols=3,sharex=True,sharey=True,
                             figsize=(7.5,7.5))
    fig.subplots_adjust(hspace=0,wspace=0)

    # # Colorbar
    # # left bottom width height
    # fig.subplots_adjust(left=0.1,right=0.9,hspace=0,wspace=0)
    # cbar_ax = fig.add_axes([0.9, 0.11, 0.02, 0.77])

    # Set limits and add labels
    axes[0,0].set_ylim(14,-4)
    axes[0,0].set_xlim(-0.6,4)
    axes[0,0].set_ylabel(f"M$_G$ (DR3)")#, fontsize=16)

    axes[1,0].set_ylabel(f"M$_G$ (DR3)")#, fontsize=16)
    axes[1,0].set_xlabel(r"G$_{BP}$ - G$_{RP}$ (DR3)")#, fontsize=16)
    axes[1,1].set_xlabel(r"G$_{BP}$ - G$_{RP}$ (DR3)")#, fontsize=16)
    axes[1,2].set_xlabel(r"G$_{BP}$ - G$_{RP}$ (DR3)")#, fontsize=16)

    axes[0,0].set_xticks(np.arange(5))
    axes[0,0].set_xticks(np.arange(0,4,0.5),minor=True)

    all_axes = axes.flatten()

    ax_in = all_axes[-1].inset_axes([0.5, 0.53, 0.47, 0.45])


    for i,cluster in enumerate(clusters):
        plot_membership_cmd(all_axes[i],cluster,alpha=0.95,ms=2)
        # all_axes[i].legend(loc=1)
        all_axes[i].text(-0.25,13,cluster,color=colors[cluster])

        plot_gaia_cmd(all_axes[-1],cluster,alpha=0.75,ms=2)
        plot_gaia_cmd(ax_in,cluster,alpha=0.95,ms=2)

    ax_in.set_xlim(-0.3,0.45)
    ax_in.set_ylim(2.5,-3)
    ax_in.set_xticks(np.arange(-0.2,0.6,0.2))
    ax_in.tick_params(labelsize=8)

    plt.savefig(os.path.join(_DIR,"plots/CMD_panel.png"),bbox_inches="tight")
    plt.savefig(os.path.join(PAPER_DIR,"fig_CMD_panel.pdf"),bbox_inches="tight")



if __name__=="__main__":

    # plot_all_cmd()
    plot_panel_cmd()
