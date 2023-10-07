import os, glob, pathlib
import numpy as np
import astropy.io.ascii as at

import tess_young
from tess_young.get_const import *
_DIR = pathlib.Path(tess_young.__file__).resolve().parent.parent
# plt.style.use(os.path.join(_DIR,'paper.mplstyle'))


if __name__=="__main__":
    from argparse import ArgumentParser
    import yaml
    # import logging

    # Define parser object
    parser = ArgumentParser(description="")
    c_group = parser.add_mutually_exclusive_group(required=True)
    c_group.add_argument("-c", "--config", dest="config_file", #required=True,
                        type=str, help="Path to config file that specifies the "
                                       "parameters for this run.")
    c_group.add_argument("-b", "--binned", dest="binned_config", #required=True,
                        type=str, help="Path to config file that specifies the "
                                       "run parameters (binned by mass)")
    args = parser.parse_args()

    if args.config_file is not None:
        # Run a regular fit
        config_file = os.path.abspath(os.path.expanduser(args.config_file))
        with open(config_file, 'r') as f:
            config = yaml.load(f.read())

        max_q = config["max_q"]
        include_blends = config["include_blends"]
        include_lit = config["include_lit"]
        output_filebase = config["output_filebase"]

        param_string_wild = f"{output_filebase}*Qmax{max_q}_blends{include_blends}_lit{include_lit}*"

        table_dir = os.path.join(_DIR,"tables/MINESweeper_v7/")
        # print(os.path.exists(table_dir))

        search_string = os.path.join(table_dir,param_string_wild)
        # print(search_string)
        output_files = sorted(glob.glob(search_string))

        # print(output_files)

        for filename in output_files:
            print("\n")
            print(filename)
            dat = at.read(filename)
            ncol = len(config["models_to_plot"])
            colnames = dat.dtype.names[:ncol]
            for colname in colnames:
                if "tophat" in colname:
                    continue
                print("\t",colname,dat["Age_"+colname][np.argmin(dat[colname])])


    else:
        # Run a fit divided up by mass bins
        config_file = os.path.abspath(os.path.expanduser(args.binned_config))
        with open(config_file, 'r') as f:
            config = yaml.load(f.read())
            config['config_file'] = config_file

        print(config)

        print("\nbinned results not implemented at this time")

