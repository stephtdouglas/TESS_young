import os, sys, time

import numpy as np
import astropy.io.ascii as at

from .tau_sq_run import run_all_models

from .periodmass import PeriodMassDistribution, PeriodMassModel
from .spinmodel import SpinModel


model_names = ["UpSco_Mattea2015","UpSco_Mattea2022","UpSco_ZeroTorque",
               "WideHat8Myr_Mattea2015","WideHat8Myr_Mattea2022","WideHat8Myr_ZeroTorque"]

def generate_synthetic_obs(model,age,period_scale,init_type,
                           n_sets=100,n_per_set=500,id_str=None,
                           start_i=None,end_i=None):
    """
    NOTE: Currently this does not write out the synthetic observations anywhere
    """

    sm = SpinModel(model,age,period_scale,init_type=init_type)
    if init_type!="kde":
        sm.normalize()

    if id_str is None:
        outf = "tausq_syn_"
    else:
        outf = id_str

    if (end_i is None) and (start_i is not None):
        i_iter = range(start_i,n_sets)
    elif (end_i is not None) and (start_i is not None):
        i_iter = range(start_i,end_i)
    else:
        i_iter = range(n_sets)
        
    for i in i_iter:
        pmd = PeriodMassModel(sm,n_select=n_per_set,rng_seed=i)
        pmd.select_obs(sm)

        print(i)
        t = time.localtime()
        current_time = time.strftime("%H:%M:%S", t)
        print(current_time)
        run_all_models(pmd=pmd,output_filebase=f"{outf}{i:04d}",
                       models_to_plot=model_names,to_plot=False)


def count_bins(pmd,sm):
    nmass = len(sm.mass_bins)-1
    n_select = np.zeros(nmass,"int")
    for i in range(nmass):
        subset = (pmd.mass>=sm.mass_bins[i]) & (pmd.mass<sm.mass_bins[i+1])
        n_select[i] = len(np.where(subset)[0])
    return n_select

def one_model(model,age,period_scale,init_type,
              max_q=0,include_blends=True,include_lit=False,
              output_filebase="tausq_SYN_binselect",id_str="SYN_binselect",
              start_i=None,end_i=None):

    # Set up the model to draw fake observations from
    sm = SpinModel(model,age,period_scale,init_type=init_type)
    sm.normalize()

    # Get the number of stars per bin from the observed data set
    pmd_obs = PeriodMassDistribution(max_q=0)
    pmd_obs.select_obs(sm)
    n_select = count_bins(pmd_obs,sm)

    # Generate a fake dataset from the model
    pmd = PeriodMassModel(sm,n_select=n_select,
                      rng_seed=9302,id_str=id_str)
    pmd.select_obs(sm)


    if (start_i is None) or (start_i==0): 
        # Compare the first synthetic set to all models
        run_all_models(pmd=pmd,output_filebase=output_filebase,
                       models_to_plot=model_names)

    # Generate multiple fake model sets and compare to all models
    generate_synthetic_obs(model,age,period_scale,init_type,n_per_set=n_select,
                           id_str=id_str,start_i=start_i,end_i=end_i)


if __name__=="__main__":

    array_id = int(os.getenv("SLURM_ARRAY_TASK_ID",9999))

    model_id = array_id // 10
    start_i = (array_id % 10) * 10
    end_i = start_i + 10
    
    ### Overal params
    period_scale="linear"

    if model_id<=2:
        ### Original run
        max_q=0

        # First, check each model for the best-fit
        filename="tables/tausq_ZAMS_Compare_Widehat_Qmax0_blendsTrue_litFalse.csv"
        ttab = at.read(filename)

        #for model in model_names[3:]:

        model = model_names[model_id+3]
        best_loc = np.argmin(ttab[model])
        best_age = ttab[f"Age_{model}"][best_loc]

        print(filename)
        print(f"best age: {best_age} Myr")

        if "WideHat" in model:
            init_type="kde"
        else:
            init_type="cluster"

        one_model(model,best_age,period_scale,init_type,
                  max_q=max_q,#include_blends=True,include_lit=False,
                  output_filebase=f"tausq_SYN_bin_{model}",id_str=f"SYN_binselect_{model}",
                  start_i=start_i,end_i=end_i)
        
    elif model_id<=5:
        ### Replace blends with literature
        max_q=0
        include_blends=False
        include_lit=True

        # First, check each model for the best-fit
        filename="tables/tausq_ZAMS_Compare_Widehat_Qmax0_blendsFalse_litTrue.csv"
        ttab = at.read(filename)

        #for model in model_names[3:]:
        model = model_names[model_id]
        best_loc = np.argmin(ttab[model])
        best_age = ttab[f"Age_{model}"][best_loc]

        print(filename)
        print(f"best age: {best_age} Myr")

        if "WideHat" in model:
            init_type="kde"
        else:
            init_type="cluster"

        one_model(model,best_age,period_scale,init_type,max_q=max_q,
                  include_blends=include_blends,include_lit=include_lit,
                  output_filebase=f"tausq_SYN2_bin_{model}",id_str=f"SYN2_binselect_{model}",
                  start_i=start_i)
    else:
        print(f"array id {array_id} invalid!")
        print(f"model {model_id}, start-end {start_i} {end_i}")
        sys.exit(42)





    # # 1) Generate a fake model set
    # model="WideHat8Myr_Mattea2022"
    # age=80
    # period_scale="linear"
    # init_type="kde"
    # sm = SpinModel(model,age,period_scale,init_type=init_type)
    #
    # sm.normalize()
    # pmd = PeriodMassModel(sm,n_select=500,rng_seed=9302)
    # pmd.select_obs(sm)
    #
    # # pmd_obs = PeriodMassDistribution(max_q=0)
    # # pmd_obs.select_obs(sm)
    # # n_select = count_bins(pmd_obs,sm)
    # #
    # # pmd = PeriodMassModel(sm,n_select=n_select,
    # #                       rng_seed=9302,id_str="binselect")
    # # pmd.select_obs(sm)
    #
    #
    # # 2) Compare the first synthetic set to all models
    # run_all_models(pmd=pmd,output_filebase="tausq_compare",
    #                models_to_plot=model_names)
    #
    # # 3) Generate multiple fake model sets and compare to all models
    # generate_synthetic_obs(model,age,period_scale,init_type)
