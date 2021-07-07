import sys,os
import glob
import itertools
import time
import logging
import urllib as urllib

logger = logging.getLogger("make_inspection_plots_tess")
logging.basicConfig(level=logging.DEBUG,format='%(asctime)s %(message)s')
logging.getLogger("matplotlib").setLevel(logging.WARNING)

import numpy as np

import matplotlib
matplotlib.use("agg")
import matplotlib.pyplot as plt
import matplotlib.ticker
import matplotlib.cm as cm
from matplotlib.font_manager import FontProperties
import matplotlib.gridspec as gridspec
from matplotlib import colors
import lightkurve as lk


# import pywcsgrid2
import astroquery
#from astroquery.sdss import SDSS
from astroquery.skyview import SkyView
from astroquery.mast import Observations
from astroquery.mast import Tesscut
from astroquery.vizier import Vizier
from astropy.wcs import WCS
import astropy.io.ascii as at
import astropy.io.fits as fits
import astropy.table as table
from astropy.table import Table, vstack
from astropy import coordinates as coords
from astropy.coordinates import SkyCoord
import astropy.units as u

# import K2fov.projection as proj
# import K2fov.fov as fov
# from K2fov.K2onSilicon import angSepVincenty,getRaDecRollFromFieldnum

# import k2phot
# from k2phot import centroid
# from k2phot import tpf_io
from k2spin import prot
# from hypra.utils import cat_match, cat_io, k2utils
# from hypra.plot import color_mag
# from hypra.plot.plot_k2 import plot_chips, setup_k2_axes
# import convertmass

cmap = plt.cm.viridis_r


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


def stereographic_map(ax, wcs, cluster=None, cluster_skycoord=None,
                      plot_sectors=[8,9,10], sector_colors=True, legend=True,
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

    ax.plot(x_clu,y_clu,'k.',ms=4,alpha=0.25,
            transform=ax.get_transform('world'))
    ax.grid('on')




def set_up_tess_sky(cluster, gs_sky, sky_ind=3):
    clusters = np.array(["IC_2391","Collinder_135","NGC_2451A","NGC_2547"])
    sector_list = [[8,9,10],[6,7,8],[7,8],[7,8,9]]

    i = np.where(clusters==cluster)[0][0]
    logging.debug(cluster)
    logging.debug(i)

    #cat = at.read(f"tables/{cluster}_crossmatch.csv",delimiter=",")
    cat = at.read(f"{cluster}_crossmatch_xmatch_TIC.csv",delimiter=",")
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

    ax = plt.subplot(gs_sky[sky_ind],projection=w)
    ax.figure.set_tight_layout(False) # WCSAxes has a tendency to cutoff ticks/labels

    ax.coords[0].set_axislabel('RA')
    ax.coords[1].set_axislabel('Dec')

    stereographic_map(ax, w, cluster=cluster, cluster_skycoord=cluster_skycoord,
                      plot_sectors=sector_list[i],sector_colors=True, legend=True,
                      hemisphere="South")

    set_xlim(np.nanmin(x_clu.degree)*0.98,np.nanmax(x_clu.degree)*1.02,w,ax)
    set_ylim(np.nanmin(y_clu.degree)-2,np.nanmax(y_clu.degree)+2,w,ax)

    return(ax)


def setup_figure(wcs,cluster):
    """
    Set up a full page figure with axes for a TESS cut out,
    two sky survey images, light curve and periodogram plots.
    """

    fig = plt.figure(figsize=(14,14))

    gs_base = gridspec.GridSpec(1,2,width_ratios=[1,3])

    gs_sky = gridspec.GridSpecFromSubplotSpec(4,1,subplot_spec=gs_base[0])

    sky_axes = [plt.subplot(gs_sky[i],projection=wcs) for i in range(3)]

    sky_axes.append(set_up_tess_sky(cluster, gs_sky, sky_ind=3))

    lc_n = 5

    gs_lc = gridspec.GridSpecFromSubplotSpec(lc_n,1,subplot_spec=gs_base[1])

    lc_axes = [plt.subplot(gs_lc[i]) for i in range(lc_n)]

    return fig, sky_axes, lc_axes


# TODO: function to read in a light curve file and return time, flux, white noise
# (Probably just time and flux in this case)


def plot_phased(ax, t, f, period, power, color):
    """
    Plot the input light curve, phased on the requested period
    """

    phased_t = t % period / period

    ax.plot(phased_t,f,'.',color=color,label="Prot={0:.2f}, Power={1:.2f}".format(float(period),float(power)))
    ax.set_xlim(0,1)
    ax.set_yticklabels([])
    ax.legend(loc="best",numpoints=1,borderaxespad=0)


# TODO: change this to accepting an input filename, because that'll be easier
# than reconstructing it from TIC, author, sector, etc

def read_tess(row,data_dir):
    """
    Read in a lightcurve. row = a row in the output_#.csv table
    """

    print(row["obs_id"],row["productFilename"])
    full_file = os.path.join(row["obs_id"],row["productFilename"])
    full_path = os.path.join(data_dir,full_file)
    lc_type = row["provenance_name"]
    sector = row["sequence_number"]#[0]
    flux_col = row["flux_cols"]

    logging.debug(flux_col)
    lc = lk.read(full_path,quality_bitmask="default",flux_column=flux_col)
    lc = lc.remove_outliers()
    time,flux = lc.time.value,lc.flux.value


    if "PATHOS"==lc_type:
        if sector==8:
            good = ((time>1519) & (time<1530)) | (time>1536.5)
        elif sector==9:
            # These are estimated
            good = ((time>1545) & (time<1556)) | (time>1558)
        elif sector==10:
            good = ((time>1572) & (time<1582)) | (time>1586)
        time, flux = time[good], flux[good]

    return time, flux


def plot_lcs(axes, TIC, row, data_dir, peaks):
    """
    Plot light curves for a given target
    """

    t, f = read_tess(row, data_dir)

    ylims = np.percentile(f[np.isfinite(f)],[0.05,99.05])

    color1 = cmap(0.35)#plt.cm.inferno(0.5)
    color2 = cmap(0.85)#plt.cm.inferno(0.75)
    color3 = cmap(0.6)
    label_fontsize = "small"

    ### Periodogram
    ls_out = prot.run_ls(t,f,np.ones_like(f),0.1,prot_lims=[0.1,70],run_bootstrap=False)
    periods, pgram = ls_out[2], ls_out[3]

    if row["num_sig_peaks"]>=2:
        peak_locs = np.where((peaks["TIC"]==row["target_name"]) &
                             (peaks["sector"]==row["sequence_number"]) &
                             (peaks["lc_type"]==row["provenance_name"]) &
                             (peaks["flux_col"]==row["flux_cols"]))[0]
        logging.debug(peak_locs)
        logging.debug(peaks[peak_locs])
        if len(peak_locs)>2:
            peak_sub = peaks[peak_locs].copy()
            peak_sub.sort("power")
            # sort_locs = peak_locs[np.argsort(peaks["power"][peak_locs])]
            logging.debug(peak_sub)
            # print(peaks["period"][sort_locs], peaks["power"][sort_locs])
            third_period = peak_sub["period"][-3]#[sort_locs[-3]]
            third_power = peak_sub["power"][-3]#[sort_locs[-3]]
    else:
        third_period, third_power = None, None

    # Plot the periodogram
    axes[0].plot(periods,pgram,'k-')
    axes[0].set_xlim(0.1,70)
    axes[0].set_xscale("log")
    axes[0].plot(row["sig_periods"],row["sig_powers"]*1.1,'v',
                mfc=color1,mec="none",ms=11)
    axes[0].plot(row["sec_periods"],row["sec_powers"]*1.1,'v',
                mfc=color2,mec="none",ms=11)
    if third_period is not None:
        axes[0].plot(third_period, third_power*1.1,'v',
                    mfc=color3,mec="none",ms=11)
    if row["sig_powers"]>0:
        ymax = row["sig_powers"]*1.15
    else:
        ymax = max(pgram)
    axes[0].set_ylim((0,ymax))
    axes[0].axhline(row["thresholds"],color='grey',ls="-.",
                    label="threshold {0:2f}".format(float(row["thresholds"])))
    #axes[0].set_xlabel("Period (d)",y=0.95)
    #axes[0].tick_params(labelbottom=False,labeltop=True)
    axes[0].set_xticklabels(["","0.1","1","10"])
    axes[0].set_ylabel("Power",fontsize=label_fontsize)
    axes[0].set_xlabel("Period (d)",fontsize=label_fontsize)

    ### Full light curve
    axes[1].plot(t,f,'k.')
    axes[1].set_ylim(ylims)
    axes[1].set_xlim(t[0],t[-1])
    axes[1].set_ylabel("Light curve",fontsize=label_fontsize)
    axes[1].set_yticklabels([])
    # axes[1].tick_params(labelbottom=False)
    axes[1].set_xlabel("Time (d)",fontsize=label_fontsize)

    ### Phase folded 1
    if row["sig_periods"]>0:
        plot_phased(axes[2],t,f,row["sig_periods"],row["sig_powers"],color=color1)
        axes[2].set_ylim(ylims)
        axes[2].set_ylabel("P1",fontsize=label_fontsize)

        if row["sig_periods"]>=2:
            repeats = np.arange(t[0],t[-1],row["sig_periods"])
            for r in repeats:
                axes[1].axvline(r,ls="--",color=color1,lw=2)
        # axes[4].tick_params(labelbottom=False)
        axes[2].set_xlabel("Phase",fontsize=label_fontsize)
    else:
        axes[2].set_axis_off()


    ### Phase folded 2
    if row["sec_periods"]>0:
        plot_phased(axes[3],t,f,row["sec_periods"],row["sec_powers"],color=color2)
        axes[3].set_ylim(ylims)
        axes[3].set_ylabel("P2",fontsize=label_fontsize)

        # print(row["sec_periods"])
        if row["sec_periods"]>=2:
            repeats = np.arange(t[0],t[-1],row["sec_periods"])
            # print(repeats)
            for r in repeats:
                # print(r)
                axes[1].axvline(r,ls=":",color=color2,lw=2)
        axes[3].set_xlabel("Phase",fontsize=label_fontsize)
    else:
        # axes[5].set_yticks([])
        # axes[5].set_yticklabels([])
        axes[3].set_axis_off()

    ### Phase folded 3
    if third_period is not None:
        plot_phased(axes[4],t,f,third_period,third_power,color=color3)
        axes[4].set_ylim(ylims)
        axes[4].set_ylabel("P3",fontsize=label_fontsize)
        axes[4].set_xlim(0,1)
        axes[4].set_xlabel("Phase",fontsize=label_fontsize)
        #
        # if third_period]>=2:
        #     repeats = np.arange(t[0],t[-1],third_period)
        #     # print(repeats)
        #     for r in repeats:
        #         # print(r)
        #         axes[1].axvline(r,ls="-.",color=color3,lw=2)
    else:
        # axes[6].set_yticks([])
        # axes[6].set_yticklabels([])
        axes[4].set_axis_off()


    plt.subplots_adjust(hspace=0.6)


# TODO: need to actually do the tesscut cutouts for this

def stamp(img, maskmap, ax=None, cmap="cubehelix"):
    """Plot a single pixel stamp."""

    if ax is None:
        fig = plt.figure(figsize=(8,8))
        ax = plt.subplot(111)

    ax.matshow(img, cmap=cmap, origin="lower", norm=colors.LogNorm())
    ax.set_xlim(-0.5,maskmap.shape[1]-0.5)
    ax.set_ylim(-0.5,maskmap.shape[0]-0.5)
    # note from K2: I remain unconvinced that these are labelled right...
    # but the coordinates are plotting right, and this matches the DSS images
    # (except East and West are flipped)!
    ax.set_ylabel("Y")
    ax.set_xlabel("X")

    return ax


def extract_wcs(dataheader):
    """
    Generate WCS for a TPF, since astropy has decided not to cooperate.

    This seems to work for both K2 and TESS TPFs.

    dataheader is the header of extension #1 (zero-indexing!) of a TPF
    """
    w5 = WCS(naxis=2)

    w5.wcs.crpix = [dataheader["1CRPX5"],dataheader["2CRPX5"]]
    w5.wcs.cdelt = np.array([dataheader["1CDLT5"],dataheader["2CDLT5"]])
    w5.wcs.crval = [dataheader["1CRVL5"],dataheader["2CRVL5"]]
    w5.wcs.ctype = [dataheader["1CTYP5"],dataheader["2CTYP5"]]
    w5.wcs.pc = [[dataheader["11PC5"],dataheader["12PC5"]],
                 [dataheader["21PC5"],dataheader["22PC5"]]]
    return w5


def plot_sky(axes, TIC, wcs, pos, tess_img, gaia_pos):
    """
    inputs
    ------
    axes: list of matplotlib.Axes objects
    TIC: TIC ID for the target star
    wcs: world coordinate system for the sky plots
    pos: astropy.coordinates.SkyCoord object, the position of the target star
    tess_img: coadded TESS FFI cutout
    """

    tess_wcs = wcs

    # # Get 2MASS K image
    # twomass_images, pix_2mass, hdr_2mass, wcs_2mass = None, None, None, None
    # try:
    #     twomass_images = SkyView.get_images(position=pos, survey=['2MASS-K'],
    #                                         pixels=500)
    #     # print(twomass_images)
    #     if len(twomass_images)>0:
    #         pix_2mass = twomass_images[0][0].data
    #         hdr_2mass = twomass_images[0][0].header
    #         wcs_2mass = WCS(hdr_2mass)
    # except (astroquery.exceptions.TimeoutError, urllib.HTTPError):
    #     pix_2mass, hdr_2mass, wcs_2mass = None, None, None
    #
    # # # Get DSS2 R image - not SDSS because this is the Hyades!
    # pix,hdr,im,wcs_sdss = None,None,None,None
    # try:
    #     dss_images = SkyView.get_images(position=pos, survey=['DSS2 Red'],
    #                                      pixels=700)
    #     if len(dss_images)>0:
    #         pix = dss_images[0][0].data
    #         hdr = dss_images[0][0].header
    #         wcs_sdss = WCS(hdr)
    # except (astroquery.exceptions.TimeoutError,urllib.HTTPError):
    #     pix, hdr, wcs_sdss = None, None, None


    # Plot TESS image
    # Plot the pixel stamp as usual, except with the WCS
    ax1 = axes[0]
    ax1.matshow(tess_img, origin="lower", cmap=cmap, norm=colors.LogNorm(),
                zorder=-10)

    # Plot TESS image
    # Plot the pixel stamp as usual, except with the WCS
    ax2 = axes[1]
    ax2.matshow(tess_img, origin="lower", cmap="Greys", norm=colors.LogNorm(),
                zorder=-10)
    # tpos = SkyCoord.from_name(f"TIC {TIC}")
    ax2.set_autoscale_on(False)
    if gaia_pos is not None:
        ax2.scatter(gaia_pos.ra,gaia_pos.dec,
                    transform=ax2.get_transform('fk5'), s=6,
                    edgecolor=cmap(0.9),facecolor=cmap(0.9),zorder=10)

    ax2.plot([7.5,9],[10,10],color=cmap(0.5),transform=ax2.get_transform('pixel'))
    ax2.plot([11,12.5],[10,10],color=cmap(0.5),transform=ax2.get_transform('pixel'))
    ax2.plot([10,10],[7.5,9],color=cmap(0.5),transform=ax2.get_transform('pixel'))
    ax2.plot([10,10],[11,12.5],color=cmap(0.5),transform=ax2.get_transform('pixel'))

    # if pix is not None:
    #     # median = np.median(pix)
    #     # stdev = np.std(pix)
    #     # levels = np.linspace(median + stdev, np.max(pix), 20)
    #     # #grey = plt.cm.Greys(0.3)
    #     # ax1[hdr].contour(pix,cmap=plt.cm.Greys, levels=levels,zorder=10)
    #
    #     # Plot the SDSS image rotated into the same frame as the pixel stamp
    #     ax2 = axes[1]
    #     vmax = np.percentile(pix,99.9)
    #     ax2.imshow(pix, origin="lower", norm=colors.LogNorm(), zorder=-10,
    #               transform=ax2.get_transform(wcs_sdss))
    #     # median2 = np.median(coadd)
    #     # stdev2 = np.std(coadd)
    #     # levels2 = np.linspace(median, np.max(coadd), 5)
    #
    # if pix_2mass is not None:
    #
    #     # Matplotlib was choking on negative values when it saved the figure
    #     pos_min = np.min(pix_2mass[pix_2mass>0])
    #     vmin = np.max([pos_min,np.percentile(pix_2mass,0.01)])
    #     vmax = np.percentile(pix_2mass,99.9)
    #     # TODO: figure out how to do this correctly now
    #     # pix_2mass.flags.writeable=True
    #     # pix_2mass[pix_2mass<=0] = pos_min
    #
    #     # Plot the 2MASS image rotated into the same frame as the pixel stamp
    #     ax4 = axes[2]
    #     ax4.imshow(pix_2mass, origin="lower", norm=colors.LogNorm(), zorder=-10,
    #               transform=ax4.get_transform(wcs_2mass))


    # Plot TESS stereographic plot of the cluster
    ax3 = axes[3]
    ax3.plot(pos.ra.value,pos.dec.value,'*',color=cmap(0.01),ms=24,
             mec=cmap(0.99))
    # ax3.set_yticks(np.arange(8,25,4))
    # ax3.set_yticks(np.arange(16,25,2))
    # ax3.set_xticks(np.arange(122,139,4))
    # ax3.set_xticks(np.arange(122,139,2),minor=True)
    ax3.set_xlabel(r"$\alpha_{2000}$",fontsize="large")
    ax3.set_ylabel(r"$\delta_{2000}$",fontsize="large")

    plt.subplots_adjust(hspace=0.15)

    # del(im,pix,hdr)
    # del(twomass_images, pix_2mass, hdr_2mass)


if __name__=="__main__":

    # TODO: make these an input from the commandline
    date = "2021-06-21"
    cluster = "NGC_2451A"
    # base_dir = os.path.expanduser(f"~/data/tess/")
    base_dir = "/data/douglaslab/tess/"
    base_dir = os.path.join(base_dir,f"{cluster.lower()}/tables/")

    filenames = glob.iglob(os.path.join(base_dir,f"*{date}*csv"))

    all_res = []
    all_peaks0 = []
    for filename in filenames:
        if "allpeaks" in filename:
            all_peaks0.append(at.read(filename))
        else:
            all_res.append(at.read(filename))

    all_peaks = vstack(all_peaks0)
    results = vstack(all_res)

    # all_peaks.rename_column("lc_type","provenance_name")
    # all_peaks.rename_column("sector","sequence_number")
    # all_peaks.rename_column("TIC","target_name")


    ###########################################################################
    ###########################################################################
    #### Now a loop over all of them
    # Set up files that will be used for all stars
    # FFI save directory
    dir_ffi = "/data2/douglaslab/tess/ffi/"

    # Cluster catalog with TIC, for ticmags
    catfile = f"{cluster}_crossmatch_xmatch_TIC.csv"
    cat = at.read(catfile,delimiter=",")

    # HDBScan file, used for gaia positions
    # hdbfile = os.path.expanduser("~/Dropbox/EDR3/scats/NGC_2451A.fits")
    hdbfile = "/data/douglaslab/EDR3/scats/NGC_2451A.fits"
    with fits.open(hdbfile) as hdu:
        hdbscan = hdu[1].data
    gaia_pos = SkyCoord(hdbscan["GAIAEDR3_RA"],
                        hdbscan["GAIAEDR3_DEC"],unit=u.degree)
    # TODO: change this all to be relevant to TESS, not K2
    with open("tables/fig_inspect.tbl","w") as f:

        for i,row in enumerate(results):

            tic = row["target_name"]
            sector, lc_type, flux_col = row["sequence_number"], row["provenance_name"], row["flux_cols"]
            output_filename = f"/data2/douglaslab/tess/{cluster.lower()}/inspect_plots/fig_inspect_{tic}_{lc_type}_{flux_col}_{sector}.png"
            if os.path.exists(output_filename):
                continue

            print(row["target_name","provenance_name","sequence_number",
                      "flux_cols"])

            cutout_coord = SkyCoord.from_name(f"TIC {tic}")
            tess_data = Tesscut.download_cutouts(cutout_coord, size=20, path=dir_ffi)
            tess_file = tess_data[-1][0]

            with fits.open(tess_file) as hdu:
                dataheader = hdu[1].header
                pixels = hdu[1].data["FLUX"]
                coadd = np.sum(pixels,axis=0)
                tess_wcs = extract_wcs(dataheader)

            fig, sky_axes, lc_axes = setup_figure(tess_wcs, cluster="NGC_2451A")

            # Define a faint mag limit for stars being plotted over the FFI stamp
            loc = np.where(cat["TIC"]==tic)[0]
            if len(loc)==1:
                faint_limit = cat["Tmag"][loc]+3
                bright = hdbscan["GAIAEDR3_G"]<=faint_limit
            elif len(loc)>1:
                print("Uh oh! Multiple TIC matches for",tic)
                faint_limit = cat["Tmag"][loc][0]+3
                bright = hdbscan["GAIAEDR3_G"]<=faint_limit
            else:
                print("Uh oh! No TIC match for",tic)
                faint_limit=0
                bright = hdbscan["GAIAEDR3_G"]<=faint_limit

            plot_sky(sky_axes, tic, tess_wcs, cutout_coord, coadd, gaia_pos[bright])
            plot_lcs(lc_axes, tic, row,
                     "/data/douglaslab/.lightkurve-cache/mastDownload/HLSP/",
                     all_peaks)
                     # os.path.expanduser("~/.lightkurve-cache/mastDownload/HLSP/"))
            plt.suptitle(f"TIC {tic}: {lc_type}, {flux_col}, Sector {sector}",
                         y=0.91,fontsize=20)

            plt.savefig(output_filename)
            # plt.savefig(f"/data2/douglaslab/tess/cluster.lower()/inspect_plots/fig_inspect_{i}.eps".format(i+700))
            f.write(f"fig_inspect_{i}.eps & TIC {tic}: {lc_type}, {flux_col}, Sector {sector}\n")
#            plt.savefig("/home/stephanie/my_papers/praeK2/fig4.eps".format(ep))
            plt.close("all")

            # Pause 1 second in an attempt to reduce server time-out
            # IF trying to query for images!
            time.sleep(0.5)

#            if i>10:
#                break
