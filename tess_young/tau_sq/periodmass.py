import os, sys, pathlib

import numpy as np
from numpy.random import default_rng
import matplotlib.pyplot as plt
import astropy.io.ascii as at
from astropy.table import Table

import tess_young
from tess_young.get_const import *
_DIR = pathlib.Path(tess_young.__file__).resolve().parent.parent

class PeriodMassDistribution:

    def __init__(self,max_q=0,include_blends=True,include_lit=False,
                 mass_limits=None,cluster="all"):
        """
        max_q: integer, maximum quality flag to include (should be 0 or 1)
        include_blends: boolean, whether or not to include potentially blended targets
        include_lit: boolean, whether or not to include literature values
        mass_limits: tuple or list, minimum and maximum masses to include
        """
        # My crossmatched catalog
        per = at.read(os.path.join(_DIR,"tab_all_stars.csv"))
        if cluster !="all":
            per = per[per["Cluster"]==cluster]
        if len(per)==0:
            print("Error: no matching stars found:",cluster)
            sys.exit(42)
        # per.dtype
        self.cat = Table(per, masked=True, copy=False)
        if mass_limits is not None:
            mass_select = (self.cat["Mass"]>=mass_limits[0]) & (self.cat["Mass"]<=mass_limits[1])
            self.cat = self.cat[mass_select]

        # Assign the catalog values to attributes
        self.prot_raw = self.cat["Prot1"]
        self.mass_raw = self.cat["Mass"]
        self.mass_err_raw = self.cat["Mass_err"]
        self.prot_mask = self.cat["Prot1"].mask
        self.mass_mask = self.cat["Mass"].mask

        # Assign masks for quality and presence of necessary values
        qmask = (self.cat["Q1"]<=max_q) & (self.cat["to_plot"]==1)
        pmask = ((self.prot_mask==False) & (self.prot_raw>0))
        pmask = pmask.filled(fill_value=False)
        mmask = (self.mass_mask==False)
        if include_blends==False:
            blmask = (self.cat["Bl?"]=="n") | (self.cat["Bl?"]=="m")
        else:
            blmask = np.ones(len(self.cat),bool)
        # If including literature values, have to replace them in the period mask
        if include_lit:
            litmask = (self.cat["LitPeriod"].mask==False) & (self.cat["LitPeriod"]>0)

            # I only want to use literature periods when I don't have a valid TESS period
            qmask_init = blmask & qmask & pmask
            use_lit = litmask & (qmask_init==False)
            self.prot_mask[use_lit] = False
            self.prot_raw.mask[use_lit] = False
            self.prot_raw[use_lit] = self.cat["LitPeriod"][use_lit]

            # Then the final selection should include literature or TESS periods, as appropriate
            lit_or_tess = qmask_init | litmask
            self.qmask = mmask & lit_or_tess
#             print(np.where(self.qmask)[0])
        else:
            self.qmask = pmask & mmask & qmask & blmask

        # I have to modify the arrays when I run tau-squared, so need to separate the raw data from down-selected
        self.prot = self.prot_raw[self.qmask]
#         # This definitely caused the mask to change, and it added masked values...
#         self.prot.mask = self.qmask
        self.mass = self.mass_raw[self.qmask]
        self.mass_err = self.mass_err_raw[self.qmask]

        self.figsize=(9,9)

        if cluster=="all":
            self.param_string = f"Qmax{max_q}_blends{include_blends}_lit{include_lit}"
        else:
            self.param_string = f"Qmax{max_q}_blends{include_blends}_lit{include_lit}_{cluster}"


    def select_obs(self, sm):
        use = ((self.mass_raw<max(sm.mass_bins)) &
               (self.mass_raw>=min(sm.mass_bins)) &
               self.qmask)

        self.prot = self.prot_raw[use]
        self.mass = self.mass_raw[use]
        self.mass_err = self.mass_err_raw[use]

    def plot_obs(self,ax=None,plot_errors=False,fig=None):
        if ax is None:
            if fig is None:
                fig = plt.figure(figsize=self.figsize)
            ax = plt.subplot(111)
            ax.set_ylim(0.07,18)
            ax.set_xlim(1.3,0.1)
            ax.tick_params(labelsize=12)
            ax.set_xlabel(r"Mass (M$_\odot$)", fontsize=16)
            ax.set_ylabel("Period (d)", fontsize=16)
            ax.patch.set_facecolor('w')
            ax.patch.set_alpha(1.0)
            fig.patch.set_facecolor('w')
            fig.patch.set_alpha(1.0)

#         # This line is throwing the nan warning
#         # And it's only the errorbar line - the regular plot line doesn't throw one
#         # Weirdly, it seems to be the prot column that has masked values (it shouldn't)
#         print(type(self.mass),np.where(self.mass.mask==True)[0])
#         print(type(self.prot),np.where(self.prot.mask==True)[0])
#         print(type(self.mass_err),np.where(self.mass_err.mask==True)[0])
        ax.errorbar(self.mass,self.prot,xerr=self.mass_err,
                    marker=None,linewidth=0,elinewidth=1,color="DarkGrey",alpha=0.75,
                    zorder=15)

        ax.plot(self.mass,self.prot,"k.",ms=2,alpha=0.75,zorder=16)

        return ax

    def calc_mass_percentiles(self,mass_bins,percentiles=[0,10,50,75,100],ntests=1000):
        self.percentiles = percentiles

        mass_errs = self.mass_err_raw
        nmass = len(mass_bins)-1

        # Define a set of benchmarks that have both periods and masses
        nb = len(np.where(self.qmask)[0])

        # Randomly generate 1000 mass samples
        rng = default_rng(42)
        new_masses = rng.normal(loc=self.mass_raw[self.qmask],
                                scale=mass_errs[self.qmask],
                                size=(ntests,nb))
#         print(np.min(new_masses), np.max(new_masses))

        # Need to generate a set of percentiles for each mass bin
        nperc = len(self.percentiles)
        period_perc = np.zeros((nmass,nperc))

        for i in range(nmass):
            subset = (new_masses>=mass_bins[i]) & (new_masses<mass_bins[i+1])
            bin_perc = np.zeros((ntests,nperc))
            for j in range(ntests):
                if len(self.prot_raw[self.qmask][subset[j]])>0:
                    bin_perc[j] = np.percentile(self.prot_raw[self.qmask][subset[j]],self.percentiles)
                else:
                    bin_perc[j][:] = np.nan
            for k in range(nperc):
                period_perc[i][k] = np.nanmedian(bin_perc[:,k])

        self.period_perc = period_perc
        self.perc_mass_bins = mass_bins

    def plot_period_perc(self,ax=None,fig=None):
        if self.period_perc is None:
            print("Please compute period percentiles first!")
            return None

        if ax is None:
            if fig is None:
                fig = plt.figure(figsize=self.figsize)
            ax = plt.subplot(111)
            ax.set_ylim(0.07,18)
            ax.set_xlim(1.3,0.1)
            ax.tick_params(labelsize=12)
            ax.set_xlabel(r"Mass (M$_\odot$)", fontsize=16)
            ax.set_ylabel("Period (d)", fontsize=16)
            ax.patch.set_facecolor('w')
            ax.patch.set_alpha(1.0)
            fig.patch.set_facecolor('w')
            fig.patch.set_alpha(1.0)

        nmass = len(self.perc_mass_bins)-1
        bin_widths = np.diff(self.perc_mass_bins)
        bin_centers = self.perc_mass_bins[:-1]+bin_widths/2

        boxes = []
        for k in range(nmass):
            boxes.append({"whislo": self.period_perc[k][0],
                  "q1": self.period_perc[k][1],
                  "med": self.period_perc[k][2],
                  "q3": self.period_perc[k][3],
                  "whishi": self.period_perc[k][4],
                  "fliers": []
                 })

        colorprop = {"color":mapper.to_rgba(1),"linewidth":1.5}
        ax.bxp(bxpstats=boxes,positions=bin_centers,widths=bin_widths,
               medianprops=colorprop,boxprops=colorprop,whiskerprops=colorprop,
               capprops=colorprop,manage_ticks=False,zorder=20)
        return ax

class PeriodMassModel(PeriodMassDistribution):
    # Replacing: __init__, calc_mass_percentiles (don't need monte carlo)
    # Inheriting: select_obs, plot_obs, plot_period_perc

    def __init__(self,sm,mass_limits=None,n_select=500,rng_seed=37,
                 id_str=None):
        """
        Inputs
        ------
        sm: SpinModel object
        mass_limits: tuple or list, minimum and maximum masses to include
        n_select: int or array-like, number of synthetic observations to 
                  generate (default 500) n_select should be evenly divisible by 
                  the number of mass_bins for sm otherwise fewer stars may be 
                  returned.
                  If array-like, len must match the number of mass-bins for sm
        """

        # Generate the synthetic observation set
        self.sm = sm

        self.sm.normalize()

        self._generate_sample(n_select,rng_seed)

        # Apply mass limits if needed
        if mass_limits is not None:
            mass_select = (self.mass_raw>=mass_limits[0]) & (self.mass_raw<=mass_limits[1])
            self.mass_raw = self.mass_raw[mass_select]
            self.prot_raw = self.prot_raw[mass_select]

        # Assign the catalog values to attributes - these are fake here since we don't
        # need to exclude any of the modelled values
        n_actual = len(self.prot_raw)
        self.mass_err_raw = np.zeros(n_actual)
        self.prot_mask = np.zeros(n_actual,bool)
        self.mass_mask = np.zeros(n_actual,bool)

        # Assign masks for quality and presence of necessary values
        self.qmask = np.ones(n_actual,bool)
        self.n_select = n_actual

        self.figsize=(9,9)

        # TODO: include model parameters here
        self.param_string = f"SYN_{self.sm.model_name}_{self.sm.age}Myr"
        if id_str is not None:
            self.param_string = self.param_string+id_str

    def _generate_sample(self,n_select,rng_seed):
        # TODO: this should incorporate a mass function

        # For now, just select an even number of stars in every mass bin

        # When the model is normalized, each mass bin has a period distribution
        # that should work as a probability distribution for np.random.choice

        try:
            len_select = len(n_select)
        except:
            len_select = 1
        # print(len_select)

        rng = np.random.default_rng(rng_seed)

        nbins = len(self.sm.mass_bins)-1
        # print(nbins)
        if len_select==1:
            n_per_bin = n_select // nbins
            fake_periods = np.zeros(nbins*n_per_bin).reshape(nbins,n_per_bin)
            fake_masses = np.zeros(nbins*n_per_bin).reshape(nbins,n_per_bin)
        else:
            if len_select!=nbins:
                print("ERROR! The lengths of n_select and sm.mass_bins must match!")
                return
            fake_periods = [np.zeros(n_in_bin) for n_in_bin in n_select]
            fake_masses = [np.zeros(n_in_bin) for n_in_bin in n_select]

        for i in range(nbins):
            bin_center = (self.sm.mass_bins[i]+self.sm.mass_bins[i+1])/2
            if len_select>1:
                n_per_bin = n_select[i]

            fake_masses[i] = np.full(n_per_bin,bin_center)


            model_loc = ((self.sm.mass_array>=self.sm.mass_bins[i]) &
                         (self.sm.mass_array<self.sm.mass_bins[i+1]))

            fake_periods[i] = rng.choice(a=self.sm.prot_array[model_loc],size=n_per_bin,
                                         replace=True,p=self.sm.prot_prob[model_loc])

        if len_select==1:
            self.prot_raw = fake_periods.flatten()
            self.mass_raw = fake_masses.flatten()
        else:
            self.prot_raw = np.concatenate(fake_periods)
            self.mass_raw = np.concatenate(fake_masses)



    def calc_mass_percentiles(self,mass_bins,percentiles=[0,10,50,75,100],ntests=1000):

        self.percentiles = percentiles

        nmass = len(mass_bins)-1

        # Need to generate a set of percentiles for each mass bin
        nperc = len(self.percentiles)
        period_perc = np.zeros((nmass,nperc))

        for i in range(nmass):
            subset = (self.mass>=mass_bins[i]) & (self.mass<mass_bins[i+1])
            period_perc[i] = np.percentile(self.prot[subset],self.percentiles)

        self.period_perc = period_perc
        self.perc_mass_bins = mass_bins



class PeriodMassBootstrap(PeriodMassDistribution):
    # Replacing: __init__
    # Inheriting: select_obs, plot_obs, plot_period_perc, calc_mass_percentiles 

    def __init__(self,pmd,mass_limits=None,rng_seed=37,
                 id_str=None):
        """
        Inputs
        ------
        pmd: another PeriodMassDistribution object to resample from
        mass_limits: tuple or list, minimum and maximum masses to include
        """

        # Generate the synthetic observation set
        self.prot_raw = pmd.prot_raw#[pmd.qmask]
        self.mass_raw = pmd.mass_raw#[pmd.qmask]
        self.mass_err_raw = pmd.mass_err_raw#[pmd.qmask]
        # Ultimately the analysis should use the original mask, 
        # but we can calculate new masses for everything
        self.qmask = pmd.qmask

        # # Apply mass limits if needed
        # if mass_limits is not None:
        #     mass_select = (self.mass_raw>=mass_limits[0]) & (self.mass_raw<=mass_limits[1])
        #     self.mass_raw = self.mass_raw[mass_select]
        #     self.mass_err_raw = self.mass_err_raw[mass_select]
        #     self.prot_raw = self.prot_raw[mass_select]

        # Assign the catalog values to attributes - these are fake here since we don't
        # need to exclude any of the modelled values
        n_actual = len(self.prot_raw)
        self.prot_mask = np.zeros(n_actual,bool)
        self.mass_mask = np.zeros(n_actual,bool)

        self.figsize=(9,9)

        # TODO: include model parameters here
        self.param_string = f"bootstrap_{pmd.param_string}"
        if id_str is not None:
            self.param_string = self.param_string+id_str

        self._generate_sample(rng_seed)

        self.prot = self.prot_raw[self.qmask]
        self.mass = self.mass_raw[self.qmask]
        self.mass_err = self.mass_err_raw[self.qmask]

    def _generate_sample(self,rng_seed):

        rng = np.random.default_rng(rng_seed)
        new_masses = rng.normal(loc=self.mass_raw,
                                scale=self.mass_err_raw)

        self.mass_raw = new_masses
        

