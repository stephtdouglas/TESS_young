

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.cm as cm

def get_colors():
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

    return norm, mapper, cmap2, colors, shapes

if __name__=="__main__":

    norm, mapper, cmap2, colors, shapes = get_colors()

    for cluster in colors.keys():
        print(cluster,mpl.colors.to_hex(colors[cluster]))
