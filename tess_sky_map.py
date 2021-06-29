"""
Plot open clusters on the sky with TESS sectors.

Includes code from John Lewis (GitHub @johnarban) and David Rodriguez
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


def parse_s_region(s_region):
    # print(s_region)
    ra = []
    dec = []
    counter = 0
    for elem in s_region.strip().split():
        try:
            value = float(elem)
        except ValueError:
            continue
        if counter % 2 == 0:
            ra.append(value)
        else:
            dec.append(value)
        counter += 1

    return {'ra': ra, 'dec': dec}


def stereographic_map(cluster=None, cluster_skycoord=None, plot_sectors=[8,9,10],
                      save=True, sector_colors=True, legend=True,
                      hemisphere="South"):
    """
    Plot TESS sectors (cycle 1 North or cycle 2 South) on a stereographic
    projection plot, in RA/Dec

    returns
    -------
    ax: matplotlib axis
    """
    obsTable = Observations.query_criteria(dataproduct_type="image", obs_collection='TESS')
    # at.write(obsTable,"map_out.csv",overwrite=False,delimiter=",")
    obsTable = obsTable[obsTable["calib_level"]==3]

    regions = {}
    for line in obsTable:
        regions[line['obs_id']] = (line['s_region'], line['sequence_number'])

    # Sector list for colors
    # For now, I'm only considering Cycles 1 and 2, not the extended mission
    if sector_colors:
        sectors = np.arange(1,27)
        print(sectors)
        colors = cm.get_cmap('viridis', len(sectors))
        colors = colors(range(len(sectors)))

    # Set up the axes
    # Use cluster coordinates to set the center of the WCS
    x_clu = cluster_skycoord.ra#geocentrictrueecliptic.lon.degree
    y_clu = cluster_skycoord.dec#geocentrictrueecliptic.lat.degree

    w = WCS(naxis=2)
    w.wcs.crpix = [0,0] # center pix doesn't matter if no image
    w.wcs.cdelt = [0.01,0.01] # also doesn't matter if no image

    # These clusters aren't near RA=0 so ignore wrapping issues.
    # Approximate data center.
    w.wcs.crval = [np.nanmedian(x_clu.degree),np.nanmedian(y_clu.degree)]
    w.wcs.ctype = ["RA---STG", "DEC--STG"]

    ax = plt.subplot(projection=w)
    ax.figure.set_tight_layout(False) # WCSAxes has a tendency to cutoff ticks/labels

    ax.coords[0].set_axislabel('RA')
    ax.coords[1].set_axislabel('Dec')

    # print("axes set")
    # Parse the edges of the camera footprints, and plot them on the sky
    for k, v in regions.items():

        s = int(k.split("-")[1][1:])

        if (s in plot_sectors)==False:
            continue
        elif (s>26):
            continue
        elif (hemisphere=="South") and (s>=14):
            continue
        elif (hemisphere=="North") and (s<14):
            continue
        try:
            coords = parse_s_region(v[0])
        except:
            print(v[0])
            continue
        xvals = coords['ra'] + [coords['ra'][0]]
        yvals = coords['dec'] + [coords['dec'][0]]
        c = SkyCoord(ra=np.array(xvals) * u.degree, dec=np.array(yvals) * u.degree, frame='icrs')
        ra = c.ra
        dec = c.dec


        if sector_colors:
            ind, = np.where(sectors == v[1])
            sec = sectors[ind]
            ax.fill(ra,dec, color=colors[sec-1][0], alpha=0.2,
                    transform=ax.get_transform('world'))
            ax.plot(ra,dec, color=colors[sec-1][0], alpha=0.5,
                    transform=ax.get_transform('world'))
        else:
            ax.plot(ra,dec, color='b', alpha=0.5,
                    transform=ax.get_transform('world'))

        # set_xlim(120,180,w,ax)
        # set_ylim(-70,-40,w,ax)

    # Now plot the cluster points

    ax.plot(x_clu,y_clu,'k.',ms=4,
            transform=ax.get_transform('world'))
    ax.grid('on')

    set_xlim(np.nanmin(x_clu.degree)*0.97,np.nanmax(x_clu.degree)*1.03,w,ax)
    set_ylim(np.nanmin(y_clu.degree)-3,np.nanmax(y_clu.degree)+3,w,ax)

    ax.set_title(cluster.replace("_"," "))

    if save:
        plt.savefig('plots/TESS_sectors_stereographic_{0}.png'.format(cluster))


    return ax


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
    clusters = ["IC_2391","Collinder_135","NGC_2451A","NGC_2547"]
    sector_list = [[8,9,10],[6,7,8],[7,8],[7,8,9]]
    for i,cluster in enumerate(clusters):
        cat = at.read(f"tables/{cluster}_crossmatch.csv",delimiter=",")
        cat_pos = SkyCoord(cat["GAIAEDR3_RA"],cat["GAIAEDR3_DEC"],unit=u.degree)

        stereographic_map(cluster,cat_pos,sector_list[i])
