"""
Plot open clusters on the sky with TESS sectors.

Includes code from John Lewis (GitHub @johnarban)
"""

import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord
from astroquery.mast import Observations
import astropy.io.ascii as at
from astropy.wcs import WCS
import matplotlib.pyplot as plt
from matplotlib import cm
# plt.interactive(False)


## helper functions for setting axis limits
def set_xlim(xmin,xmax,w,ax):
    xmin = w.all_world2pix([xmin],[w.wcs.crval[1]],1)[0]
    xmax = w.all_world2pix([xmax],[w.wcs.crval[1]],1)[0]
    ax.set_xlim(xmin,xmax)

def set_ylim(ymin,ymax,w,ax):
    ymin = w.all_world2pix([w.wcs.crval[0]],[ymin],1)[1]
    ymax = w.all_world2pix([w.wcs.crval[0]],[ymax],1)[1]
    ax.set_ylim(ymin,ymax)


def plot_clusters_sterographic():
    clusters = ["IC_2391","Collinder_135","NGC_2451A","NGC_2547"]
    sector_list = [[6,7,8],[8,9,10],[7,8],[7,8,9]]
    for i,cluster in enumerate(clusters):
        cat = at.read(f"tables/{cluster}_crossmatch.csv",delimiter=",")
        cat_pos = SkyCoord(cat["GAIAEDR3_RA"],cat["GAIAEDR3_DEC"],unit=u.degree)

        cluster_skycoord = cat_pos
        x_clu = cluster_skycoord.ra#geocentrictrueecliptic.lon.degree
        y_clu = cluster_skycoord.dec#geocentrictrueecliptic.lat.degree

        w = WCS(naxis=2)
        w.wcs.crpix = [0,0] # center pix doesn't matter if no image
        w.wcs.cdelt = [0.01,0.01] # also doesn't matter if no image

        # These clusters aren't near RA=0 so ignore wrapping issues
        # approximate data center.
        w.wcs.crval = [np.nanmedian(x_clu.degree),np.nanmedian(y_clu.degree)]
        w.wcs.ctype = ["RA---STG", "DEC--STG"]

        ax = plt.subplot(projection=w)
        ax.figure.set_tight_layout(False) # WCSAxes has a tendency to cutoff ticks/labels

        ax.coords[0].set_axislabel('RA')
        ax.coords[1].set_axislabel('Dec')


        ax.plot(x_clu,y_clu,'k.',
                transform=ax.get_transform('world'))
        ax.grid('on')

        set_xlim(np.nanmin(x_clu.degree)*0.95,np.nanmax(x_clu.degree)*1.05,w,ax)
        set_ylim(np.nanmin(y_clu.degree)-5,np.nanmax(y_clu.degree)+5,w,ax)


        ax.set_title(cluster.replace("_"," "))
        plt.show()

if __name__=="__main__":
