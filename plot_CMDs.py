import os
import numpy as np
from numpy.polynomial.polynomial import Polynomial
import astropy.io.ascii as at
import matplotlib.pyplot as plt

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

def plot_clusters_mist(ages,avs,plot_name="",plot_title="",model_colors=True):
    clusters = ["IC_2391","Collinder_135","NGC_2451A","NGC_2547","IC_2602"]
    symbols = ["o","s","^","v","d"]

    gdr = "EDR3"

    plt.figure(figsize=(8,8))
    ax1 = plt.subplot(111)
    ax2 = ax1.inset_axes([0.5, 0.5, 0.47, 0.47])
    for i,cluster in enumerate(clusters):
        cat = at.read(f"tables/{cluster}_crossmatch.csv",delimiter=",")
        bp_rp = cat[f"GAIA{gdr}_BP"] - cat[f"GAIA{gdr}_RP"]
        dist = 1000/cat["GAIAEDR3_PARALLAX_CORRECTED"]
        # This is a fudge factor because the MIST isochrones aren't lining
        # up with the gaia data
        cmd_off = 0#0.1
        abs_G = cat[f"GAIA{gdr}_G"] - 5*np.log10(dist) + 5 + cmd_off
        ax1.plot(bp_rp,abs_G,symbols[i],label=cluster,ms=3,alpha=0.5)
        ax2.plot(bp_rp,abs_G,symbols[i],label=cluster,ms=3)

    mist_dir = os.path.expanduser("~/Dropbox/Models/MIST_v1.2_vvcrit0.0_UBVRIplus/")
    mist = at.read(os.path.join(mist_dir,"MIST_v1.2_feh_p0.00_afe_p0.0_vvcrit0.0_UBVRIplus.iso.cmd"),
                   header_start=12)
    mist_dir2 = os.path.expanduser("~/Dropbox/Models/MIST_v1.2_vvcrit0.4_UBVRIplus/")
    mist2 = at.read(os.path.join(mist_dir2,"MIST_v1.2_feh_p0.00_afe_p0.0_vvcrit0.4_UBVRIplus.iso.cmd"),
                   header_start=12)

    for i, cluster in enumerate(clusters):
        model_ind = np.argmin(abs(mist["log10_isochrone_age_yr"]-np.log10(ages[cluster]*1e6)))
        model_age = mist["log10_isochrone_age_yr"][model_ind]
        model_age_Myr = 10**(model_age-6)
        model = mist["log10_isochrone_age_yr"]==model_age
        model_bp_rp0 = mist["Gaia_BP_EDR3"][model]-mist["Gaia_RP_EDR3"][model]
        model_g = mist["Gaia_G_EDR3"][model]

        model_ind2 = np.argmin(abs(mist2["log10_isochrone_age_yr"]-np.log10(ages[cluster]*1e6)))
        model_age2 = mist["log10_isochrone_age_yr"][model_ind2]
        model_age_Myr2 = 10**(model_age2-6)
        model2 = mist2["log10_isochrone_age_yr"]==model_age2
        model_bp_rp02 = mist2["Gaia_BP_EDR3"][model2]-mist2["Gaia_RP_EDR3"][model2]
        model_g2 = mist2["Gaia_G_EDR3"][model2]

        if avs[cluster]>0:
            _, A_G = calc_gaia_extinction(model_bp_rp0,avs[cluster],"G")
            E_BP_RP = calc_gaia_reddening(model_bp_rp0,avs[cluster])
            _, A_G2 = calc_gaia_extinction(model_bp_rp02,avs[cluster],"G")
            E_BP_RP2 = calc_gaia_reddening(model_bp_rp02,avs[cluster])
        else:
            A_G, E_BP_RP = 0,0
            A_G2, E_BP_RP2 = 0,0

        if model_colors:
            ax1.plot(model_bp_rp0 + E_BP_RP, model_g + A_G,"-",color=f"C{i}",
                     label=f"{model_age_Myr:.0f} Myr, A_V={avs[cluster]:.2f}, vvcrit0.0")
            ax2.plot(model_bp_rp0 + E_BP_RP, model_g + A_G,"-",color=f"C{i}")

            # ax1.plot(model_bp_rp02 + E_BP_RP2, model_g2 + A_G2,":",color=f"C{i}",
            #          label=f"{model_age_Myr:.0f} Myr, A_V={avs[cluster]:.2f}, vvcrit0.4")
            ax2.plot(model_bp_rp02 + E_BP_RP2, model_g2 + A_G2,":",color=f"C{i}")
        else:
            ax1.plot(model_bp_rp0 + E_BP_RP, model_g + A_G,"-",color=f"grey",
                     label=f"{model_age_Myr:.0f} Myr, A_V={avs[cluster]:.2f}, vvcrit0.0")
            ax2.plot(model_bp_rp0 + E_BP_RP, model_g + A_G,"-",color=f"grey")

            # ax1.plot(model_bp_rp02 + E_BP_RP2, model_g2 + A_G2,":",color=f"grey",
            #          label=f"{model_age_Myr:.0f} Myr, A_V={avs[cluster]:.2f}, vvcrit0.4")
            ax2.plot(model_bp_rp02 + E_BP_RP2, model_g2 + A_G2,":",color=f"grey")

    ax1.plot([],[],":",color="k",label="vvcrit0.4")
    ax1.legend(loc=3)

    ax1.set_xlabel(f"BP-RP ({gdr})")
    ax1.set_ylabel(f"M_G ({gdr})")

    ax1.set_ylim(13,-5)
    ax1.set_xlim(-0.6,4)

    ax2.set_xlim(-0.3,0.35)
    ax2.set_ylim(2,-3)

    ax1.set_title(plot_title)
    plt.savefig(f"plots/CMD_MIST{plot_name}_{gdr}.png")

    # plt.show()
    plt.close()

if __name__=="__main__":

    # Cantat-Gaudin+2020
    cg20_ages = {"IC_2391": 28.8,
                 "IC_2602": 36.3,
                 "NGC_2547": 32.4,
                 "NGC_2451A": 35.5,
                 "Collinder_135": 26.3}
    cg20_av = {"IC_2391": 0.04,
                 "IC_2602": 0.03,
                 "NGC_2547": 0.14,
                 "NGC_2451A": 0,
                 "Collinder_135": 0.01}
    # Kharchenko+2013
    khar_ages = {"IC_2391": 112.2,
                 "IC_2602": 221.3,
                 "NGC_2547": 77.6,
                 "NGC_2451A": 57.5,
                 "Collinder_135": 39.8}
    # They give E(B-V), so A_V = 3.1E(B-V)
    khar_av = {"IC_2391": 3.1*0.052,
                 "IC_2602": 3.1*0.031,
                 "NGC_2547": 3.1*0.04,
                 "NGC_2451A": 0,
                 "Collinder_135": 3.1*0.042}
    # Ghoza 2012
    g12_ages = {"IC_2391": 46,
                 "IC_2602": 32,
                 "NGC_2547": 38,
                 "NGC_2451A": 60,
                 "Collinder_135": 26}

    # Not real ages, I just want to overplot some curves
    age_span = {"IC_2391": 20,
                 "IC_2602": 40,
                 "NGC_2547": 60,
                 "NGC_2451A": 80,
                 "Collinder_135": 120}
    fake_av = {"IC_2391": 0.25,
                 "IC_2602": 0.25,
                 "NGC_2547": 0.25,
                 "NGC_2451A": 0.25,
                 "Collinder_135": 0.25}
    plot_clusters_mist(cg20_ages,cg20_av,"_CG20","Cantat-Gaudin+2020 Ages")
    plot_clusters_mist(khar_ages,khar_av,"_khar","Kharchenko+2013 Ages")
    # plot_clusters_mist(g12_ages,"_G12","Ghoza+2012 Ages")
    plot_clusters_mist(age_span,fake_av,"_20_120","Ages 20-120 Myr",model_colors=False)
