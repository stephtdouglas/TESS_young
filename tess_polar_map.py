# Script to make TESS polar map
# David Rodriguez
# Aug 13, 2019

import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord
from astroquery.mast import Observations
import astropy.io.ascii as at
import matplotlib.pyplot as plt
from matplotlib import cm
plt.interactive(False)

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

def polar_map(save=True, sector_colors=True, legend=True,
              hemisphere="South"):
    """
    Plot TESS sectors (cycle 1 North or cycle 2 South) on a polar projection plot

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
    if sector_colors:
        sectors = np.arange(1,27)
        print(sectors)
        colors = cm.get_cmap('viridis', len(sectors))
        colors = colors(range(len(sectors)))

    # Make theplot
    fig = plt.figure(figsize=(12, 8))
    ax = fig.add_subplot(1, 1, 1, projection="polar")
    # if hemisphere=="North":
    #     ax.set_rlim(bottom=90, top=0)
    #     # ax.set_rmin(94)
    #     # ax.set_rmax(0)

    # print("axes set")
    for k, v in regions.items():
        # if ("0014" in k) or ("0015" in k) or ("0016" in k) or ("0017" in k) or ("0018" in k):
        # print(k,v)
        if (int(k.split("-")[1][1:])>26):
            continue
        elif (hemisphere=="South") and (int(k.split("-")[1][1:])>=14):
            continue
        elif (hemisphere=="North") and (int(k.split("-")[1][1:])<14):
            continue
        try:
            coords = parse_s_region(v[0])
        except:
            print(v[0])
            continue
        xvals = coords['ra'] + [coords['ra'][0]]
        yvals = coords['dec'] + [coords['dec'][0]]
        c = SkyCoord(ra=np.array(xvals) * u.degree, dec=np.array(yvals) * u.degree, frame='icrs')
        x_rad = c.geocentrictrueecliptic.lon.wrap_at(180 * u.deg).radian
        y_rad = c.geocentrictrueecliptic.lat.degree
        if hemisphere=="North":
            y_rad = y_rad*-1
        # print("parsed")

        if sector_colors:
            ind, = np.where(sectors == v[1])
            sec = sectors[ind]
            ax.fill(x_rad, y_rad, color=colors[sec-1][0], alpha=0.2)
            ax.plot(x_rad, y_rad, color=colors[sec-1][0], alpha=0.5)
        else:
            ax.plot(x_rad, y_rad, color='b', alpha=0.5)

        # print("plotted")
    # print("done plotting sectors")
    # print(hemisphere,ax.get_rmin(),ax.get_rmax())

    # if sector_colors:
    #     for i, s in enumerate(sectors):
    #         ax.plot([np.nan], [np.nan], color=colors[i], alpha=0.5, label='Sector {}'.format(s))

    # ax.grid(color='grey', ls='dotted')
    # if legend: plt.legend(loc='upper right', bbox_to_anchor=(1.2, 1.0), fontsize='medium')
    # plt.tight_layout()

    # if hemisphere=="North":
    #     ax.set_rticklabels([80,60,40,20,0])

    if save:
        plt.savefig('plots/TESS_sectors_{0}_blank.png'.format(hemisphere.lower()))

    return ax

def polar_cluster(ax,cluster=None,cluster_skycoord=None,save=True,
                  hemisphere="South"):
    """
    Plot members of an open cluster on a polar projection plot

    returns
    -------
    ax: matplotlib axis
    """


    # x_rad = c.geocentrictrueecliptic.lon.wrap_at(180 * u.deg).radian
    # y_rad = c.geocentrictrueecliptic.lat.degree
    x_clu = cluster_skycoord.geocentrictrueecliptic.lon.wrap_at(180*u.deg).radian
    y_clu = cluster_skycoord.geocentrictrueecliptic.lat.degree
    if hemisphere=="North":
        y_clu = y_clu*-1


    ax.plot(x_clu,y_clu,'k.',alpha=0.5,ms=1)

    # ax.set_title(cluster)

    if cluster is not None:
        x_avg = np.average(x_clu)
        y_avg = np.average(y_clu)
        ax.text(x_avg,y_avg+5,cluster,color="k",horizontalalignment="right")

    if save:
        plt.savefig('plots/TESS_sectors_{0}_{1}.png'.format(hemisphere.lower(),
                                                    cluster))
    return ax

def plot_all_clusters():
    ax = polar_map(hemisphere="South",save=False)

    clusters = ["IC_2391","Collinder_135","NGC_2451A","NGC_2547","IC_2602"]
    for cluster in clusters:
        cat = at.read(f"tables/{cluster}_crossmatch.csv",delimiter=",")
        cat_pos = SkyCoord(cat["GAIAEDR3_RA"],cat["GAIAEDR3_DEC"],unit=u.degree)
        ax = polar_cluster(ax,cluster,cat_pos,save=False)

    plt.savefig("plots/TESS_sectors_South_all.png")

    return ax


if __name__=="__main__":
    _ = plot_all_clusters()

    # ax = polar_map(hemisphere="South",save=False)

    # mclusters = ["IC_2391","IC_2602","NGC_2451A","NGC_2547"]
    # cat = at.read("catalogs/Meingast2021_GaiaEDR3_xmatch.csv",delimiter=",")
    # print(np.unique(cat["Cluster"]))
    # cat_pos = SkyCoord(cat["_RAJ2000"],cat["_DEJ2000"],unit=u.degree)
    # for cluster in mclusters:
    #     print(cluster)
    #     loc = (cat["Cluster"]==cluster.replace("_"," ")) & (cat["fc"]<0.1)
    #     print(len(np.where(loc)[0]))
    #     ax = polar_cluster(ax,cluster,cat_pos[loc],save=False)
    #
    # plt.savefig("plots/TESS_sectors_South_Meingast.png")


    # polar_map(hemisphere="North")

    plt.close("all")
