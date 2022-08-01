import os, sys, glob, time

import numpy as np
import matplotlib.pyplot as plt
import astropy.io.ascii as at

from tau_sq import SpinModel, PeriodMassModel
from tau_sq_run import run_all_models

model_names = ["UpSco_Mattea2015","UpSco_Mattea2022","UpSco_ZeroTorque",
               "WideHat8Myr_Mattea2015","WideHat8Myr_Mattea2022","WideHat8Myr_ZeroTorque"]

def generate_synthetic_obs(model,age,period_scale,init_type,
                           n_sets=100,n_per_set=500):

    sm = SpinModel(model,age,period_scale,init_type=init_type)
    if init_type!="kde":
        sm.normalize()

    for i in range(n_sets):
        pmd = PeriodMassModel(sm,n_select=n_per_set,rng_seed=i)
        pmd.select_obs(sm)

        print(i)
        run_all_models(pmd=pmd,output_filebase=f"tausq_syn_{i:04d}",
                       models_to_plot=model_names,to_plot=False)

if __name__=="__main__":

    # 1) Generate a fake model set
    model="WideHat8Myr_Mattea2022"
    age=80
    period_scale="linear"
    init_type="kde"
    sm = SpinModel(model,age,period_scale,init_type=init_type)

    sm.normalize()
    pmd = PeriodMassModel(sm,n_select=500,rng_seed=9302)
    pmd.select_obs(sm)

    # 2) Generate multiple fake model sets (run on the cluster)
    run_all_models(pmd=pmd,output_filebase="tausq_compare",
                   models_to_plot=model_names)

    generate_synthetic_obs(model,age,period_scale,init_type)
