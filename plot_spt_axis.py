import os, sys

import matplotlib.pyplot as plt
import astropy.io.ascii as at
import numpy as np

def add_spt_gaia(ax,spts=["F8V","G2V","K0V","K4V","K7V","M0V","M2V","M3V","M4V","M5V"]):

    mamajek_file = os.path.expanduser("~/Dropbox/Models/mamajek_colors.dat")
    mamajek = at.read(mamajek_file,fill_values=[("...","NaN")])

    ymax = ax.get_ylim()[1]
    yscale = ax.get_yscale()
    if yscale=="log":
        ytext = ymax * 1.1
    else:
        ytext = ymax + 5

    for spt in spts:
        loc = np.where(mamajek["SpT"]==spt)[0]
        if len(loc)==1:
            text = spt.replace("V","")
            ax.text(mamajek["Bp-Rp"][loc],ytext,text)
