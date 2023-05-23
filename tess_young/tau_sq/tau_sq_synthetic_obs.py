import os, sys, time, pathlib#, datetime
#import multiprocessing as mp

import numpy as np
import astropy.io.ascii as at

from tau_sq_run import run_all_models

from periodmass import PeriodMassDistribution, PeriodMassModel
from spinmodel import SpinModel

from tess_young.get_const import *
import tess_young
_DIR = pathlib.Path(tess_young.__file__).resolve().parent.parent

def generate_synthetic_obs(model,age,period_scale,init_type,
                           n_sets=100,n_per_set=500,id_str=None,
                           start_i=None,end_i=None,mass_limits=None,
                           cluster="all"):
    """
    Generate sets of synthetic observations based on an input spinmodel, 
    and compare to all the models

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
    n_sets: (integer) the number of synthetic datasets to create (default 100)
    n_per_set: (integer) the number of stars to include in each synthetic 
               observation (default 500)
    id_str: (string) run identifier, will be included in output filenames
    start_i: (integer) index within n_sets to start at. 
             If None or 0, will start at 0
    end_i: (integer) index within n_sets to end at. If None, will end at n_sets.

    NOTE: Currently this does not write out the synthetic observations anywhere
    It only writes out the tau-sq values from the fit to each dataset
    (As part of run_all_models)
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
        

    # all_params = []


    for i in i_iter:
        pmd = PeriodMassModel(sm,n_select=n_per_set,rng_seed=i,mass_limits=mass_limits)
        pmd.select_obs(sm)

        # params = [pmd,period_scale,f"{outf}{i:04d}",model_names[3:],
        #           None,mass_limits,False,False,
        #           [init_type,init_type,init_type]]
        # all_params.append(params)

        print(i)
        t = time.localtime()
        current_time = time.strftime("%H:%M:%S", t)
        print(current_time)
        run_all_models(pmd=pmd,output_filebase=f"{outf}{i:04d}",
                       models_to_plot=model_names[3:],to_plot=False,
                       init_types=[init_type,init_type,init_type],
                       mass_limits=mass_limits)
#        break

    # try:
    #     cpu_ct = len(os.sched_getaffinity(0))-1
    #     print(f"Current time: {datetime.datetime.now()} -- AO cpu_count HPC:", cpu_ct)
    # except AttributeError:
    #     cpu_ct = mp.cpu_count()-1
    #     print(f"Current time: {datetime.datetime.now()} -- AO cpu_count PC:", cpu_ct)

    # chunksize = n_sets // cpu_ct

    # with mp.Pool(cpu_ct) as pool:
    #     _ = pool.starmap(run_all_models, all_params, chunksize=chunksize)
       

def count_bins(pmd,sm):
    """
    Count the number of observed stars in each mass bin of the model. 
    """
    nmass = len(sm.mass_bins)-1
    n_select = np.zeros(nmass,"int")
    for i in range(nmass):
        subset = (pmd.mass>=sm.mass_bins[i]) & (pmd.mass<sm.mass_bins[i+1])
        n_select[i] = len(np.where(subset)[0])
    return n_select

def one_model(model,age,period_scale,init_type,
              max_q=0,include_blends=True,include_lit=False,
              output_filebase="tausq_SYN_binselect",id_str="SYN_binselect",
              start_i=None,end_i=None,mass_limits=None,cluster="all"):
    """
    Configure and run one set of synthetic observations, including the number
    of stars and the baseline comparison with a synthetic dataset. 

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
    max_q: integer, maximum quality flag to include (should be 0 or 1)
    include_blends: boolean, whether or not to include potentially blended targets
    include_lit: boolean, whether or not to include literature values
    mass_limits: tuple or list, minimum and maximum masses to include
    n_sets: (integer) the number of synthetic datasets to create (default 100)
    n_per_set: (integer) the number of stars to include in each synthetic 
               observation (default 500)
    output_filebase: (string) identifier for the baseline comparison
    id_str: (string) run identifier, will be included in output filenames w/ an index
    start_i: (integer) index within n_sets to start at. 
             If None or 0, will start at 0
    end_i: (integer) index within n_sets to end at. If None, will end at n_sets.

    """

    # Set up the model to draw fake observations from
    print(model,age,period_scale,init_type)
    sm = SpinModel(model,age,period_scale,init_type=init_type)
    sm.normalize()

    # Get the number of stars per bin from the observed data set
    pmd_obs = PeriodMassDistribution(max_q=max_q,include_blends=include_blends,
                                     include_lit=include_lit,
                                     mass_limits=mass_limits,cluster=cluster)
    pmd_obs.select_obs(sm)
    n_select = count_bins(pmd_obs,sm)

    # Generate a fake dataset from the model
    pmd = PeriodMassModel(sm,n_select=n_select,mass_limits=mass_limits,
                      rng_seed=9302,id_str=id_str)
    pmd.select_obs(sm)

    if (start_i is None) or (start_i==0): 
        # Compare the first synthetic set to all models
        # This produces a baseline tausq vs. age curve
        run_all_models(pmd=pmd,output_filebase=output_filebase,
                       models_to_plot=model_names[3:],
                       init_types=[init_type,init_type,init_type])

    # Generate multiple fake model sets and compare to all models
    # Another script will analyze these results and select the best-fit from each synthetic dataset
    generate_synthetic_obs(model,age,period_scale,init_type,n_per_set=n_select,
                           id_str=id_str,start_i=start_i,end_i=end_i)

if __name__=="__main__":
    from argparse import ArgumentParser
    import logging

    # Define parser object
    parser = ArgumentParser(description="")

    # Define all the necessary arguments
    parser.add_argument("model",help="model to draw synthetic observations from")

    parser.add_argument("init_type", help="init type (tophat, cluster, or kde)")

    parser.add_argument("-p", dest="period_scale", default="linear", required=False,
                        help="linear or log")

    parser.add_argument("-q", dest="max_q", default=0, required=False, help="maximum Q flag to include")

    parser.add_argument("-b", dest="include_blends", default=True, required=False, type=bool,
                        help="whether or not to include potentially blended stars")

    parser.add_argument("-l", dest="include_lit", default=True, required=False, type=bool,
                        help="whether or not to include literature periods")

    parser.add_argument("-o", dest="output_filebase", default="tausq_SYN", required=False,
                        help="identifier for this run")

    parser.add_argument("--high", dest="mass_high", default=1.4, required=False,
                        help="upper limit mass in solar masses",type=float)

    parser.add_argument("--low", dest="mass_low", default=0.05, required=False,
                        help="lower limit mass in solar masses",type=float)

    parser.add_argument("-s", dest="start_i", required=False)

    parser.add_argument("-e", dest="end_i", required=False)

    parser.add_argument("-g", dest="cluster", default="all", required=False,
                        help="Which group/cluster to fit (default is all)")

    age_group = parser.add_mutually_exclusive_group(required=True)
    age_group.add_argument("-a","--age",dest="age", type=int)
    age_group.add_argument("-f", "--filename", dest="filename")


    # parse the input
    args = parser.parse_args()


    print(args)
    # print(model)
    # print(age)

    id_str = args.output_filebase+"_baseline"

    if args.age is not None:
        best_age = args.age
    else:
        try:
            ttab = at.read(args.filename)
        except:
            print("File not found or not readable?")
            raise

        colname = f"{args.model}_{args.init_type}"

        try:
            best_loc = np.argmin(ttab[colname])
        except:
            print("Model/init_type not found in input reference file?",colname)
            raise
        best_age = ttab[f"Age_{colname}"][best_loc]

    mass_limits = [args.mass_low, args.mass_high]
        
    print(args.init_type)
        
    one_model(args.model, best_age, args.period_scale, args.init_type,
              args.max_q, args.include_blends, args.include_lit,
              args.output_filebase, id_str,
              args.start_i, args.end_i, mass_limits, args.cluster)
