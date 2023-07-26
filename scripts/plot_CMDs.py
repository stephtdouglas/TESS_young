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
    cat = at.read(os.path.join(_DIR,f"{cluster}_crossmatch_xmatch_TIC.csv"),delimiter=",")

    cat = at.read(os.path.join(_DIR,"tab_all_stars.csv"),delimiter=",")
    cat = cat[cat["Cluster"]==cluster]

    bp_rp = cat[f"GAIA{gdr}_BP"] - cat[f"GAIA{gdr}_RP"]
    dist = 1000/cat["GAIAEDR3_PARALLAX_CORRECTED"]
    abs_G = cat[f"GAIA{gdr}_G"] - 5*np.log10(dist) + 5
    ax.plot(bp_rp,abs_G,shapes[cluster],color=colors[cluster],label=cluster,**kwargs)


def plot_all_cmd():
    plt.figure()
    ax1 = plt.subplot(111)
    ax2 = ax1.inset_axes([0.5, 0.5, 0.47, 0.47])
    for cluster in clusters:
        plot_gaia_cmd(ax1,cluster,alpha=0.75,ms=4)
        plot_gaia_cmd(ax2,cluster,alpha=0.95,ms=4)
    plt.legend(loc=3)
    ax1.set_xlabel(r"G$_{BP}$ - G$_{RP}$ (EDR3)")
    ax1.set_ylabel(f"M$_G$ (EDR3)")

    ax1.set_ylim(13,-5)
    ax1.set_xlim(-0.6,4)

    ax2.set_xlim(-0.3,0.45)
    ax2.set_ylim(2.5,-3)
    plt.savefig(os.path.join(_DIR,"plots/CMD_all.png"),bbox_inches="tight")
    plt.savefig(os.path.join(PAPER_DIR,"fig_CMD_all.pdf"),bbox_inches="tight")


if __name__=="__main__":

    plot_all_cmd()
