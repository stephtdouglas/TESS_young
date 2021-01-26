import os

import matplotlib.pyplot as plt
import numpy as np
import astropy.io.ascii as at
import astropy.io.fits as fits
from astropy.table import Table, vstack
from astroquery.mast import Observations

def download_lightcurves(ticname,pipeline,sectors=[8,9]):
    """
    Retrieve HLSP light curves from MAST for a particular TIC ID

    Inputs:
        ticname: string, TIC id in the format "TIC ########"
        pipeline: string, HLSP source (e.g., "CDIPS","PATHOS","QLP")

    Returns:
        found: bool, whether a valid light curve was found
        multisector: bool, whether multiple light curves were found
        manifest: Table, with downloaded filenames and input TIC ids

    """

    found = False
    multisector = False

    obs_table = Observations.query_object(ticname,radius="5 arcsec")
#     print(obs_table["obs_id"])
    obsids = obs_table[obs_table["provenance_name"]==pipeline]['obsid']
    if len(obsids)==0:
        print("no light curves found for",ticname,pipeline)
        return found, multisector, []
    data_products_by_id = Observations.get_product_list(obsids)
    data_products = Observations.get_product_list(obsids)
    manifest = Observations.download_products(data_products)
    # files are downloaded to "./mastDownload/HLSP/"

    sec_str0 = [f"{s:04}" for s in sectors]
    print(sec_str0)

    if pipeline=="CDIPS":
        sec_str = [f"-{st}-" for st in sec_str0]
        # Format is -000#-cam#-ccd#
    elif (pipeline=="PATHOS") or (pipeline=="QLP"):
        sec_str = [f"s{st}" for st in sec_str0]
        # format is "s000#"
    else:
        sec_str = sec_str0
        # Just check for the raw sector string, and hope for the best.
        # TODO: Flag with a warning, though

    found_count = 0
    to_delete = []
    manifest["sector"] = np.zeros(len(manifest),"int")
    for i, sector in enumerate(sectors):
        # Check each light curve to make sure it's
        # a) the right star TODO

        # and b) the right sector
        for j, filename in enumerate(manifest["Local Path"]):
            if (".fits" in filename)==False:
                continue
            elif sec_str[i] in filename:
                found = True
                found_count += 1
                manifest["sector"][j] = sector
                break
            else:
                continue

        if found_count>1:
            multisector = True

    for j, filename in enumerate(manifest["Local Path"]):
        if (".fits" in filename)==False:
            to_delete.append(j)
    manifest.remove_rows(to_delete)

    manifest["TIC"] = np.empty(len(manifest),"U14")
    manifest["TIC"][:] = ticname

    manifest["pipeline"] = np.empty(len(manifest),"U6")
    manifest["pipeline"][:] = pipeline

    return found, multisector, manifest


def download_cluster(catalog_filename,cluster_name,sectors,
                     pipelines=["PATHOS","CDIPS","QLP"]):
    """
    Download all MAST HLSP light curves for a given cluster

    Inputs:
    catalog_filename: string, full or relative path to a .csv file
                      must include a TICID column
    cluster_name: string, name of the cluster
    sectors: listlike, sectors to check for light curves for this cluster
    pipelines: listlike, names of HLSP sources (e.g., "CDIPS","PATHOS","QLP")

    Outputs:
    saves a .csv file containing the local filepaths, TIC IDs, and HLSP sources
    for the light curves. Format: {cluster_name}_mast_lightcurves.csv

    """

    memb = at.read(catalog_filename)

    manifest_list = []

    for i,tic in enumerate(memb["TIC"]):
        for source in pipelines:
            found, multisector, manifest = download_lightcurves(f"TIC {tic}",
                                                                source,
                                                                sectors=sectors
                                                                )

            if found is False:
                print("No light curves found for",tic,source)
            else:
                manifest_list.append(manifest)

        # TODO: when running this for real, pause to try to prevent server
        # errors

    all_manifest = vstack(manifest_list)

    at.write(all_manifest,f"{cluster_name}_mast_lightcurves.csv",delimiter=",")

if __name__=="__main__":

    membfile = "IC2391_zorro_DR2memb_EDR3data.csv"

    download_cluster(catalog_filename=membfile,cluster_name="IC2391",
                    sectors=[8,9,10])

    # manifest_list = []
    #
    # found, multisector, manifest = download_lightcurves(ticname="TIC 93270923",
    #                                                     pipeline="PATHOS",
    #                                                     sectors=[8,9])
    # print("PATHOS",found,multisector)
    # print(manifest)
    # manifest_list.append(manifest)
    #
    # found, multisector, manifest = download_lightcurves(ticname="TIC 93270923",
    #                                                     pipeline="CDIPS",
    #                                                     sectors=[8,9])
    # print("CDIPS",found,multisector)
    # print(manifest)
    # manifest_list.append(manifest)
    #
    # all_manifest = vstack(manifest_list)
    # print(all_manifest["TIC","pipeline"])
    # print(all_manifest.dtype)
