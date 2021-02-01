# Script to make TESS polar map
# David Rodriguez
# Aug 13, 2019

import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord
from astroquery.mast import Observations
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
              hemisphere="South",cluster=None,cluster_skycoord=None):
    obsTable = Observations.query_criteria(dataproduct_type="image", obs_collection='TESS')
    obsTable = obsTable[obsTable["calib_level"]==3]

    regions = {}
    for line in obsTable:
        regions[line['obs_id']] = (line['s_region'], line['sequence_number'])

    # Sector list for colors
    if sector_colors:
        sectors = np.unique([s[1] for s in regions.values()])
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

    if cluster is not None:

        # x_rad = c.geocentrictrueecliptic.lon.wrap_at(180 * u.deg).radian
        # y_rad = c.geocentrictrueecliptic.lat.degree
        x_clu = cluster_skycoord.geocentrictrueecliptic.lon.wrap_at(180*u.deg).radian
        y_clu = cluster_skycoord.geocentrictrueecliptic.lat.degree
        if hemisphere=="North":
            y_clu = y_clu*-1


        ax.plot(x_clu,y_clu,'k.')

        ax.set_title(cluster)

    #
    # if sector_colors:
    #     for i, s in enumerate(sectors):
    #         ax.plot([np.nan], [np.nan], color=colors[i], alpha=0.5, label='Sector {}'.format(s))

    # ax.grid(color='grey', ls='dotted')
    # if legend: plt.legend(loc='upper right', bbox_to_anchor=(1.2, 1.0), fontsize='medium')
    # plt.tight_layout()

    # if hemisphere=="North":
    #     ax.set_rticklabels([80,60,40,20,0])

    if save:
        plt.savefig('TESS_sectors_{0}_{1}.png'.format(hemisphere.lower(),
                                                      cluster))
        plt.close()

    plt.show()

if __name__=="__main__":
    # polar_map(hemisphere="South")

    polar_map(hemisphere="North")
