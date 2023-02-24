import os

import numpy as np
from numpy.random import default_rng
import matplotlib.pyplot as plt
import astropy.io.ascii as at
from astropy.table import Table

import ..get_colors
norm, mapper, cmap2, colors, shapes = get_colors.get_colors()
plt.style.use('./paper.mplstyle')
PAPER_DIR = os.path.expanduser("~/my_papers/TESS_young/")



clusters = ["IC_2391","Collinder_135","NGC_2451A","NGC_2547","IC_2602"]
dates = ["2021-06-22","2021-06-18","2021-06-21","2021-06-21","2021-07-02"]


display_names = {"UpSco_Mattea2015":"Matt+15; UpSco initialization",
                 "UpSco_Mattea2022":"Matt+in prep; UpSco initialization",
                 "UpSco_ZeroTorque":"Zero Torque; UpSco initialization",
                 "WideHat8Myr_Mattea2015":"Matt+15; uniform initialization",
                 "WideHat8Myr_Mattea2022":"Matt+in prep; uniform initialization",
                 "WideHat8Myr_ZeroTorque":"Zero Torque; uniform initialization"}
model_names = ["UpSco_Mattea2015","UpSco_Mattea2022","UpSco_ZeroTorque",
               "WideHat8Myr_Mattea2015","WideHat8Myr_Mattea2022","WideHat8Myr_ZeroTorque"]
nmod = 6
nage = 118

class SpinModel:

    def __init__(self,model,model_age,period_scale,init_type="cluster"):
        """
        Inputs
        ------
        model: (string) torque law used. Options:
                             "UpSco_Mattea2015", "UpSco_Mattea2022",
                             "UpSco_"ZeroTorque", "WideHat8Myr_Mattea2015",
                             "WideHat8Myr_Mattea2022", "WideHat8Myr_ZeroTorque"
        model_age: (integer) model age in Myr. Must correspond to a valid age
                   of the chosen model
        period_scale: (string) "log" or "linear"
        init_type: (string) initialization type. "tophat", "cluster", or "kde"
        """

        self.mass_bins = np.arange(0.05,1.4,0.1)

        self.period_scale = period_scale
        self.age = model_age
        self.model_name = model

        if period_scale=="log":
#             self.period_bins = np.logspace(np.log10(0.08),np.log10(40),200)
            self.period_bins = np.logspace(np.log10(0.08),np.log10(40),30)
        else:
            self.period_bins = np.linspace(0,40,30)
#             self.period_bins = np.linspace(0,40,10)

        if init_type in ["tophat","cluster","kde"]:
            self.init_type = init_type
        else:
            print("WARNING: init_type unknown, no normalization will be performed")

        # Read in the model file and calculate the histogram
        if self.init_type=="kde":
            kde_pfile = os.path.expanduser("~/Dropbox/Models/UpSco_KDE_init/all_KDE_corrected.csv")
            kde_prob = at.read(kde_pfile)

            mod_file = os.path.expanduser(f"~/Dropbox/Models/{model}/{model}_{model_age:05d}Myr.txt")
        #     print(mod_file)
            mod = at.read(mod_file,names=["mass","prot"])
            mod["Prob"] = kde_prob["Prob"]
            self.prot_prob = mod["Prob"]

            img_raw, self.xedges, self.yedges = np.histogram2d(mod["mass"],mod["prot"],
                                                               weights=mod["Prob"],density=True,
                                                               bins=[self.mass_bins,self.period_bins])
        else:
            mod_file = os.path.expanduser(f"~/Dropbox/Models/{model}/{model}_{model_age:05d}Myr.txt")
        #     print(mod_file)
            mod = at.read(mod_file,names=["mass","prot"])
            img_raw, self.xedges, self.yedges = np.histogram2d(mod["mass"],mod["prot"],
                                                               bins=[self.mass_bins,self.period_bins])
            self.prot_prob = np.ones_like(mod["prot"])
        # Transpose the image so we can actually plot it
        # by default this is not normalized
        self.img = img_raw.T

        self.mass_array = mod["mass"]
        self.prot_array = mod["prot"]


    def normalize(self):
        # Normalize the histogram if desired
        # calculate dM and dP
        dmass = np.diff(self.mass_bins)
        dper = np.diff(self.period_bins)

        # Each mass bin should contain (1/nmass)*(1/dmass) of the probability
        nmass = len(dmass)

        if (self.init_type=="tophat") or (self.init_type=="cluster"):
            img_nomask = np.copy(self.img)

            # Multiply by dP, to sum correctly
            img_nomask_dp = img_nomask * dper[:,np.newaxis]
            # Then add up the probability in each mass bin
            sum_over_p = np.sum(img_nomask_dp, axis=0)
            # Then normalize each mass bin
            img_nomask_mnorm1 = img_nomask / sum_over_p
            img_nomask_final = img_nomask_mnorm1 / nmass / dmass

            self.img = img_nomask_final
        elif self.init_type=="kde":
            print("init_type KDE; normalization already completed")
        else:
            print("init_type unknown, no normalization performed")

    def add_mask(self):
        img_nomask = np.copy(self.img)

        # mask the image so it doesn't show cells outside the model
        model_exists = np.ones(np.shape(img_nomask),bool)
        mask = np.zeros(np.shape(img_nomask),bool)
        for i in range(len(self.mass_bins)-1):
            mass_loc = ((self.mass_array>=self.mass_bins[i]) &
                        (self.mass_array<self.mass_bins[i+1]))
            # Calculate whether there are periods in each individual bins
            for j in range(len(self.period_bins)-1):
                per_loc = ((self.prot_array>=self.period_bins[j]) &
                           (self.prot_array<self.period_bins[j+1]))
                in_this_bin = np.where(mass_loc & per_loc)[0]
    #             print(mass_bins[i],period_bins[j],in_this_bin)
                if len(in_this_bin)==0:
                    model_exists[j,i] = False
            # Now, for this mass range, define the mask to only exclude bins
            # beyond the range of the model
    #         print(model_exists[:,i])
            mod_min_j = min(np.where(model_exists[:,i]==True)[0])
            mod_max_j = max(np.where(model_exists[:,i]==True)[0])
            mask[:mod_min_j,i] = True
            mask[mod_max_j+1:,i] = True

        self.img = np.ma.masked_array(img_nomask,mask=mask)
        self.mask = mask

    def plot_hist(self,ax=None,fig=None):
        if ax is None:
            if fig is None:
                fig = plt.figure()
            ax = plt.subplot(111)
            ax.set_ylim(min(self.period_bins),max(self.period_bins))
            ax.set_xlim(1.3,0.1)
            ax.tick_params(labelsize=12)
            ax.set_xlabel(r"Mass (M$_\odot$)", fontsize=16)
            ax.set_ylabel("Period (d)", fontsize=16)
            ax.patch.set_facecolor('w')
            ax.patch.set_alpha(1.0)
            fig.patch.set_facecolor('w')
            fig.patch.set_alpha(1.0)

        X, Y = np.meshgrid(self.xedges, self.yedges)
        if self.init_type=="kde":
            ax.pcolormesh(X, Y, self.img,cmap="viridis_r",vmin=1e-4,vmax=0.5)
        else:
            ax.pcolormesh(X, Y, self.img,cmap="viridis_r")
    #     ax.plot(mod["mass"],mod["prot"],'k.',alpha=0.25)
        ax.set_yscale(self.period_scale)

        return ax

    def calc_tau_sq(self, pmd):
        """ calculate tau-squared for an observed PeriodMass distribution
        """
        # area of the region containing the model
        A_pm = ((np.max(self.mass_bins)-np.min(self.mass_bins)) *
                (np.max(self.period_bins)-np.min(self.period_bins)))

        # model weight? value from section 5.2.3
        fscript = 0.7

        # background term in the tau squared sum
        bkgd_i = (1-fscript)/A_pm

        pmd.select_obs(self)

        nprot = len(pmd.prot)
    #     print(nprot)
        sum_tau_sq = 0
        found_count = 0
        for j in range(len(self.period_bins)-1):
            in_p_bin = (pmd.prot>self.period_bins[j]) & (pmd.prot<=self.period_bins[j+1])
            for i in range(len(self.mass_bins)-1):
                if self.mask[j,i]==True:
                    # No model at this index; move on
                    continue
                else:
                    in_m_bin = (pmd.mass>self.mass_bins[i]) & (pmd.mass<=self.mass_bins[i+1])
                    observed = in_p_bin & in_m_bin
                    n_in_bin = len(np.where(observed)[0])
                    if n_in_bin>0:
                        found_count += n_in_bin
                        this_rho_f = fscript * self.img[j,i]
                        this_rho_prime = this_rho_f + bkgd_i
                        sum_tau_sq += np.log(this_rho_prime) * n_in_bin

        n_missed = nprot - found_count
        sum_tau_sq += np.log(bkgd_i) * n_missed

        self.tau_sq = -2 * sum_tau_sq

    def calc_tau_sq_binned(self, pmd):
        """ calculate tau-squared for each individual mass bin in an observed
            PeriodMass distribution
        """
        # area of the region containing the model
        A_pm = ((np.max(self.mass_bins)-np.min(self.mass_bins)) *
                (np.max(self.period_bins)-np.min(self.period_bins)))

        # model weight? value from section 5.2.3
        fscript = 0.7

        # background term in the tau squared sum
        bkgd_i = (1-fscript)/A_pm

        pmd.select_obs(self)

        nprot = len(pmd.prot)
    #     print(nprot)
        sum_tau_sq = np.zeros(len(self.mass_bins)-1)
        for i in range(len(self.mass_bins)-1):
            in_m_bin = (pmd.mass>self.mass_bins[i]) & (pmd.mass<=self.mass_bins[i+1])
            n_in_m_bin = len(np.where(in_m_bin)[0])
            found_count = 0
            for j in range(len(self.period_bins)-1):
                if self.mask[j,i]==True:
                    # No model at this index; move on
                    continue
                else:
                    in_p_bin = (pmd.prot>self.period_bins[j]) & (pmd.prot<=self.period_bins[j+1])
                    observed = in_p_bin & in_m_bin
                    n_in_bin = len(np.where(observed)[0])
                    if n_in_bin>0:
                        found_count += n_in_bin
                        this_rho_f = fscript * self.img[j,i]
                        this_rho_prime = this_rho_f + bkgd_i
                        sum_tau_sq[i] += np.log(this_rho_prime) * n_in_bin

            n_missed = n_in_m_bin - found_count
            sum_tau_sq[i] += np.log(bkgd_i) * n_missed

        self.tau_sq = -2 * sum_tau_sq
