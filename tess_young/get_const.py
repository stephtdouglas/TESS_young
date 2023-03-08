import os

import numpy as np
import matplotlib as mpl
import matplotlib.cm as cm

__all__ = ['norm', 'mapper','cmap2', 'colors', 'shapes',
           'clusters','dates','display_names','model_names','nmod','nage','PAPER_DIR',
           'MODEL_DIR']

PAPER_DIR = os.path.expanduser("~/my_papers/TESS_young/")

MODEL_DIR = os.path.expanduser("~/Dropbox/Models/")

norm = mpl.colors.LogNorm(vmin=0.1, vmax=30)
mapper = cm.ScalarMappable(norm=norm, cmap=cm.viridis)

cmap2 = cm.get_cmap("viridis",7)
colors = {"IC_2391": cmap2(0),
        "IC_2602": cmap2(4),
        "NGC_2547": cmap2(3),
        "NGC_2451A": cmap2(2),
        "Collinder_135": cmap2(1)}


shapes= {"IC_2391": "o",
        "IC_2602": "d",
        "NGC_2547": "v",
        "NGC_2451A": "^",
        "Collinder_135": "s"}

clusters = ["IC_2391","Collinder_135","NGC_2451A","NGC_2547","IC_2602"]
dates = ["2021-06-22","2021-06-18","2021-06-21","2021-06-21","2021-07-02"]


display_names = {"UpSco_Mattea2015":"Matt+15; UpSco initialization",
                 "UpSco_Mattea2022":"Matt+in prep; UpSco initialization",
                 "UpSco_ZeroTorque":"Zero Torque; UpSco initialization",
                 "WideHat8Myr_Mattea2015":"Matt+15; uniform initialization",
                 "WideHat8Myr_Mattea2022":"Matt+in prep; uniform initialization",
                 "WideHat8Myr_ZeroTorque":"Zero Torque; uniform initialization"}
model_names = np.asarray(["UpSco_Mattea2015","UpSco_Mattea2022",
                          "UpSco_ZeroTorque","WideHat8Myr_Mattea2015",
                          "WideHat8Myr_Mattea2022","WideHat8Myr_ZeroTorque"])
nmod = 6
nage = 118
