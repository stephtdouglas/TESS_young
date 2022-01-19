import os, sys
import glob
import itertools

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.cm as cm
import astropy.io.fits as fits
import astropy.io.ascii as at
from astropy.table import join,vstack,Table
import astropy.table as table
from scipy import stats
from scipy.interpolate import interp1d


def read_validation_results(cluster, date, which=None):
    # Read in my visual inspection results
    if which is None:
        vis_file = f"tables/{cluster}_{date}_results_comments.csv"
    else:
        vis_file = f"tables/{cluster}_{date}_results_comments{which}.csv"
    vis = at.read(vis_file,delimiter=",")
    good = np.where(vis["Select"].mask==False)[0]
    # print(len(good))

    # Read in all the peaks, for ID'ing third periods
    peaks = at.read(f"tables/{cluster}_{date}_allpeaks.csv",delimiter=",")

    # Limit the table to only the light curves I analyzed
    vis = vis[good]
    vis.rename_column("\ufefftarget_name","TIC")

    # Select the final period based on quality flags
    vis["final_period"] = np.copy(vis["sig_periods"])
    vis["final_power"] = np.copy(vis["sig_powers"])
    vis["final_Q"] = np.copy(vis["Q"])
    vis["second_period"] = np.ones_like(vis["final_period"])*-9999
    vis["second_power"] = np.ones_like(vis["final_power"])*-9999
    vis["second_Q"] = np.ones_like(vis["final_Q"])*-9999

    # If I flaggged the highest peak as bad, but selected another peak,
    # Select that one instead
    replace2 = (vis["Q"]==2) & ((vis["Q2"]==1) | (vis["Q2"]==0))
    replace3 = (vis["Q"]==2) & ((vis["Q3"]==1) | (vis["Q3"]==0))
    vis["final_period"][replace2] = vis["sec_periods"][replace2]
    vis["final_Q"][replace2]==vis["Q2"][replace2]
    vis["final_power"][replace2] = vis["sec_powers"][replace2]

    ploc0 = np.ones(len(peaks),bool)
    for i in np.where(replace3)[0]:
        ploc0[:] = True
        for col in ["TIC","provenance_name","sequence_number","flux_cols"]:
            ploc0 = ploc0 & (vis[col][i]==peaks[col])

        ploc = np.where(ploc0)[0]
        if len(ploc)>2:
            # The peaks are sorted in ascending order when I created the peaks
            # file, so I just need to go backwards within the matching set
            vis["final_period"][i] = peaks["period"][-3]
            vis["final_Q"][i]==vis["Q3"][i]
            vis["final_power"][i] = peaks["power"][-3]
        else:
            print(f"Insufficient peaks for {vis['TIC'][i]}!")
            print(vis["TIC","provenance_name","sequence_number","flux_cols",
                      "Q","Q2","Q3","Notes"][i])
            print(peaks[ploc])


    # The MultiProt column didn't necessarily mean that two good periods
    # were detected, so use the Q values to find multiperiodic stars
    good1 = (vis["Q"]==0) | (vis["Q"]==1)
    good2 = (vis["Q2"]==0) | (vis["Q2"]==1)
    good3 = (vis["Q3"]==0) | (vis["Q3"]==1)
    multi_q = (good1 & good2) | (good1 & good3) | (good2 & good3)

    q2 = ["Q2","Q3","Q3"]

    for gi in np.where((good1 & good2))[0]:
        vis["second_period"][gi] = vis["sec_periods"][gi]
        vis["second_power"][gi] = vis["sec_powers"][gi]
        vis["second_Q"][gi] = vis["Q2"][gi]

    for i in np.where((good1 & good3) | (good2 & good3))[0]:
        ploc0[:] = True
        for col in ["TIC","provenance_name","sequence_number","flux_cols"]:
            ploc0 = ploc0 & (vis[col][i]==peaks[col])

        ploc = np.where(ploc0)[0]
        if len(ploc)>2:
            # The peaks are sorted in ascending order when I created the peaks
            # file, so I just need to go backwards within the matching set
            vis["second_period"][i] = peaks["period"][-3]
            vis["second_Q"][i]==vis["Q3"][i]
            vis["second_power"][i] = peaks["power"][-3]
        else:
            print(f"Insufficient peaks for {vis['TIC'][i]}!")
            print(vis["TIC","provenance_name","sequence_number","flux_cols",
                      "Q","Q2","Q3","Notes"][i])
            print(peaks[ploc])

    # Return the output table
    return vis

def make_final_period_catalog(cluster, date, to_plot=False):

    # Retrieve the two validation catalogs
    cat1 = read_validation_results(cluster,date)
    # Temporary fix because I haven't finished re-validating IC 2602
    if cluster!="IC_2602":
        cat2 = read_validation_results(cluster,date,which=2)
    else:
        cat2 = cat1.copy()

    # Crossmatch the validation catalogs on TIC IDs
    allcat = join(cat1,cat2,keys=["TIC"])

    # Create new columns to hold the final-final values for
    # period, Q, light curve info, peak power, threshold, multi/spot flags
    # Colnames from D19:
    # Prot1, Pw1, Q1, Sig, Prot2, Pw2, Q2, MP?, SE?
    # Bl? (this will be an automated flag, unlike my previous papers)
    new_cols_float = [np.zeros(len(allcat))*np.nan for i in range(7)]
    new_cols_char = [np.zeros(len(allcat),"U1") for i in range(3)]
    allcat.add_columns(new_cols_float,names=["Prot1", "Pw1", "Q1", "Sig",
                                             "Prot2", "Pw2", "Q2"])
    allcat.add_columns(new_cols_char,names=["MP?","SE?","Bl?"])
    # Apparently I did not include my notes in the output files for D19
    # So I'm leaving them out for simplicity now
    allcat.add_column(np.zeros(len(allcat),"U5"),name="provenance_name")
    allcat.add_column(np.zeros(len(allcat),"U10"),name="flux_cols")
    allcat.add_column(np.zeros(len(allcat),int),name="sequence_number")
    allcat.add_column(np.zeros(len(allcat),"U200"),name="obs_id")

    ##### Identify cases where my results differed between steps
    # For both primary and possible secondary periods
    check_columns = np.array(["provenance_name","flux_cols","sequence_number",
                              "final_period","final_Q"])
    check_columns2 = np.array(["provenance_name","flux_cols","sequence_number",
                              "second_period","second_Q"])
    diff_results = np.zeros(len(allcat),bool)
    diff_results2 = np.zeros(len(allcat),bool)
    for col in check_columns:
        diff_results = diff_results | (allcat[f"{col}_1"]!=allcat[f"{col}_2"])
    for col in check_columns2:
        diff_results2 = diff_results2 | (allcat[f"{col}_1"]!=allcat[f"{col}_2"])

    ## If my results agreed, just assign the first parameters to the final cols
    good = diff_results==False
    good2 = diff_results2==False
    print(f"Primary: {len(np.where(good)[0])} matches,",
           f"{len(np.where(diff_results)[0])} mismatches")
    print(f"Secondary: {len(np.where(good2)[0])} matches,",
           f"{len(np.where(diff_results2)[0])} mismatches")
    print(f"Combined: {len(np.where(good & good2)[0])} matches,",
           f"{len(np.where(diff_results | diff_results2)[0])} mismatches")
    print(f"Good primary, mismatched secondary:",
           f"{len(np.where(good & diff_results2)[0])}")
    print(f"Mismatched primary, matched secondary:",
           f"{len(np.where(good2 & diff_results)[0])}")
    init_cols = ["provenance_name_1","flux_cols_1","sequence_number_1",
                 "obs_id_1","final_period_1","final_Q_1","final_power_1",
                 "thresholds_1"]
    output_cols = ["provenance_name","flux_cols","sequence_number",
                   "obs_id","Prot1", "Q1", "Pw1", "Sig"]
    init_cols2 = ["second_period_1","second_Q_1","second_power_1"]
    output_cols2 = ["Prot2", "Q2", "Pw2"]
    ncols = len(init_cols)
    ncols2 = len(init_cols2)
    if (len(output_cols)!=ncols) | (len(output_cols2)!=ncols2):
        print("Fix columns!")
        sys.exit(0)

    for i in range(ncols):
        allcat[output_cols[i]][good] = allcat[init_cols[i]][good]
    for i in range(ncols2):
        allcat[output_cols[i]][good2] = allcat[init_cols[i]][good2]

    ## Cycle through all the stars with differing results
    diff_idx = np.where(diff_results)[0]
    print(len(diff_idx))

    primary_still_bad = []
    for i in diff_idx:
        diff_frac = (abs(allcat["final_period_1"][i]-allcat["final_period_2"][i])
                     / min(allcat["final_period_1","final_period_2"][i]))

        half_dbl = allcat["final_period_1"][i]/allcat["final_period_2"][i]

        # If one was flagged as Q=2, flag both as Q=2
        if (allcat["final_Q_1"][i]==2) or (allcat["final_Q_2"][i]==2):
            allcat["Q1"][i] = 2

        # If the period is (almost) the same (and Q=1 or Q=0):
        # elif allcat["final_period_1"][i]==allcat["final_period_2"][i]:
        elif diff_frac<0.05:

            # TODO: but it's just a different light curve, pick the one with the
            # higher periodogram peak and assign the lower Q value
            # higher_peak = np.argmax(allcat["final_"])

            # If I assigned different Q values, assign the lower Q value
            allcat["Q1"][i] = max(allcat["final_Q_1","final_Q_2"][i])
            allcat["Prot1"][i] = allcat["final_period_1"][i]
            allcat["Pw1"][i] = allcat["final_power_1"][i]
            allcat["provenance_name"][i] = allcat["provenance_name_1"][i]
            allcat["flux_cols"][i] = allcat["flux_cols_1"][i]
            allcat["sequence_number"][i] = allcat["sequence_number_1"][i]
            allcat["obs_id"][i] = allcat["obs_id_1"][i]
            allcat["Sig"][i] = allcat["thresholds_1"][i]

        # If they're harmonics and the second period is longer, choose the 2nd
        elif (abs(half_dbl-0.5)<0.02):
            allcat["Q1"][i] = allcat["final_Q_2"][i]
            allcat["Prot1"][i] = allcat["final_period_2"][i]
            allcat["Pw1"][i] = allcat["final_power_2"][i]
            allcat["provenance_name"][i] = allcat["provenance_name_2"][i]
            allcat["flux_cols"][i] = allcat["flux_cols_2"][i]
            allcat["sequence_number"][i] = allcat["sequence_number_2"][i]
            allcat["obs_id"][i] = allcat["obs_id_2"][i]
            allcat["Sig"][i] = allcat["thresholds_2"][i]

        # If they're harmonics and the first period is longer, choose the 1st
        elif (abs(half_dbl-2)<0.08):
            allcat["Q1"][i] = allcat["final_Q_1"][i]
            allcat["Prot1"][i] = allcat["final_period_1"][i]
            allcat["Pw1"][i] = allcat["final_power_1"][i]
            allcat["provenance_name"][i] = allcat["provenance_name_1"][i]
            allcat["flux_cols"][i] = allcat["flux_cols_1"][i]
            allcat["sequence_number"][i] = allcat["sequence_number_1"][i]
            allcat["obs_id"][i] = allcat["obs_id_1"][i]
            allcat["Sig"][i] = allcat["thresholds_1"][i]

        # If the periods are not the same...
        # TODO: so far most of these look impossibly confused, though a few
        # could be salvageable. I'll probably have to have a list of stars
        # that I've resolved by hand at this point.
        else:
            print("\n",allcat["TIC"][i],"primary mismatch")
            print(allcat["TIC","provenance_name_1","flux_cols_1",
                         "sequence_number_1","final_period_1","final_Q_1",
                         "Notes_1"][i])
            print(allcat["TIC","provenance_name_2","flux_cols_2",
                         "sequence_number_2","final_period_2","final_Q_2",
                         "Notes_2"][i])
            primary_still_bad.append(i)


    ## Cycle through all the stars with differing secondary results
    diff_idx2 = np.where(diff_results2)[0]
    print(len(diff_idx2))

    secondary_still_bad = []
    for i in diff_idx2:
        diff_frac = (abs(allcat["second_period_1"][i]-allcat["second_period_2"][i])
                     / min(allcat["second_period_1","second_period_2"][i]))

        half_dbl = allcat["second_period_1"][i]/allcat["second_period_2"][i]

        # If one was flagged as Q=2, flag both as Q=2
        if (allcat["second_Q_1"][i]==2) or (allcat["second_Q_2"][i]==2):
            allcat["Q2"][i] = 2

        # If the period is (almost) the same (and Q=1 or Q=0):
        # elif allcat["final_period_1"][i]==allcat["final_period_2"][i]:
        elif diff_frac<0.05:

            # TODO: but it's just a different light curve, pick the one with the
            # higher periodogram peak and assign the lower Q value
            # higher_peak = np.argmax(allcat["final_"])

            # If I assigned different Q values, assign the lower Q value
            allcat["Q2"][i] = max(allcat["second_Q_1","second_Q_2"][i])
            allcat["Prot2"][i] = allcat["second_period_1"][i]
            allcat["Pw2"][i] = allcat["second_power_1"][i]

        # If they're harmonics and the second period is longer, choose the 2nd
        elif (abs(half_dbl-0.5)<0.02):
            allcat["Q2"][i] = allcat["second_Q_2"][i]
            allcat["Prot2"][i] = allcat["second_period_2"][i]
            allcat["Pw2"][i] = allcat["second_power_2"][i]

        # If they're harmonics and the first period is longer, choose the 1st
        elif (abs(half_dbl-2)<0.08):
            allcat["Q2"][i] = allcat["second_Q_1"][i]
            allcat["Prot2"][i] = allcat["second_period_1"][i]
            allcat["Pw2"][i] = allcat["second_power_1"][i]

        # If the periods are not the same...
        # TODO: so far most of these look impossibly confused, though a few
        # could be salvageable. I'll probably have to have a list of stars
        # that I've resolved by hand at this point.
        else:
            print("\n",allcat["TIC"][i],"secondary mismatch")
            print(allcat["TIC","provenance_name_1","flux_cols_1",
                         "sequence_number_1","final_period_1","final_Q_1",
                         "second_period_1","second_Q_1"][i])
            print(allcat["TIC","provenance_name_2","flux_cols_2",
                         "sequence_number_2","final_period_2","final_Q_2",
                         "second_period_2","second_Q_2"][i])
            secondary_still_bad.append(i)

    primary_still_bad = np.array(primary_still_bad)
    secondary_still_bad = np.array(secondary_still_bad)
    still_bad = np.asarray(np.union1d(primary_still_bad,secondary_still_bad),int)

    if len(still_bad)>0:
        resolved = at.read("resolved_discrepant_validations.dat")
        for i in still_bad:
            loc = np.where(allcat["TIC"][i]==resolved["TIC"])[0]
            if len(loc)!=1:
                print("uh oh!",i,allcat["TIC"][i])
            which = resolved["Which"][loc][0]
            allcat["Q1"][i] = resolved["final_Q"][loc]
            allcat["Prot1"][i] = resolved["final_period"][loc]
            allcat["Pw1"][i] = allcat[f"final_power{which}"][i]
            allcat["provenance_name"][i] = resolved["provenance_name"][loc][0]
            allcat["flux_cols"][i] = resolved["flux_cols"][loc][0]
            allcat["sequence_number"][i] = resolved["sequence_number"][loc][0]
            if resolved["second_Q"][loc]<=2:
                allcat["Q2"][i] = resolved["second_Q"][loc]
            if np.isfinite(resolved["second_period"][loc]):
                allcat["Prot2"][i] = resolved["second_period"][loc]
                allcat["Pw2"][i] = allcat[f"second_power{which}"][i]

            allcat["obs_id"][i] = allcat[f"obs_id{which}"][i]
            allcat["Sig"][i] = allcat[f"thresholds{which}"][i]

    #### Finally, transfer the multiperiodic and spot evolution flags
    allcat["MP?"][:] = "n"
    allcat["MP?"][(allcat["MultiProt_1"]=="y") &
                  (allcat["MultiProt_2"]=="y")] = "y"
    allcat["MP?"][(allcat["MultiProt_1"]=="m") |
                  (allcat["MultiProt_2"]=="m")] = "m"

    allcat["SE?"][:] = "n"
    allcat["SE?"][(allcat["SpotEvol_1"]=="y") &
                  (allcat["SpotEvol_2"]=="y")] = "y"
    allcat["SE?"][(allcat["SpotEvol_1"]=="m") |
                  (allcat["SpotEvol_2"]=="m")] = "m"

    ##### Output a pure period table
    periods = allcat["TIC","provenance_name","flux_cols","sequence_number",
                   "obs_id","Prot1", "Q1", "Pw1","Prot2", "Q2", "Pw2", "Sig",
                   "MP?","SE?"]

    formats = {"Prot1":'%.2f',
               "Pw1":'%.3f',
               "Prot2":'%.2f',
               "Pw2":'%.3f',
               "Sig":'%.3f'}

    # Write out the table. I'm not going to write a tex table - I'll just put
    # the colnames in the manuscript proper
    if os.path.exists("tab_all_tess_periods.csv"):
        print("WARNING: Period table already exists; not re-recreating it")
    else:
        at.write(periods,"tab_all_tess_periods.csv",delimiter=",",formats=formats)

    ##### Crossmatch to the Gaia catalogs, save relevant columns

    # Use the _crossmatch_xmatch_TIC catalogs for their TIC IDs
    tic_file = f"{cluster}_crossmatch_xmatch_TIC.csv"
    tic_match0 = at.read(tic_file,delimiter=",")
    tic_match = tic_match0
    # tic_match = table.unique(tic_match0,silent=True)
    # print(len(tic_match0),len(tic_match))
    if tic_match.masked == False:
        tic_match = Table(tic_match, masked=True, copy=False)

    tic_match = tic_match["TIC","GAIAEDR3_ID","angDist","GAIA","Gmag","Tmag",
                          "Flag"]
    tic_match.rename_column("angDist","TIC_angDist")
    tic_match.rename_column("GAIA","TIC_GAIADR2_ID")
    tic_match.rename_column("Gmag","TIC_Gmag")
    tic_match.rename_column("Tmag","TIC_Tmag")
    tic_match.rename_column("Flag","TIC_Flag")

    ##########################################################################
    ##########################################################################
    ##########################################################################

    # Then match back to my original catalogs (and when I get Phill's
    # MINESweeper results, those should have the same columns, plus extras)
    cat_init_file = f"tables/{cluster}_crossmatch_allcolumns.fits"
    with fits.open(cat_init_file) as hdu:
        cat_init = Table(hdu[1].data)

    if cat_init.masked == False:
        cat_init = Table(cat_init, masked=True, copy=False)

    # Check for doubles in the initial catalog, just in case
    unique_gaia, gct = np.unique(cat_init["GAIAEDR3_ID"], return_counts=True)
    dbl_gaia = unique_gaia[gct>1]
    dbl_gaia_idx = []
    for gid in dbl_gaia:
        loc = np.where(cat_init["GAIAEDR3_ID"]==gid)[0]
        dbl_gaia_idx.append(loc)
    if len(dbl_gaia_idx)>0:
        dbl_gaia_idx = np.sort(np.concatenate(dbl_gaia_idx))


    xmatch = join(cat_init,tic_match,keys="GAIAEDR3_ID",join_type="left")

    # What to do when two Gaia targets match the same TIC ID?
    unique_gaia, gct = np.unique(xmatch["GAIAEDR3_ID"], return_counts=True)
    dbl_gaia = unique_gaia[gct>1]

    # Two TIC IDs matching the same Gaia ID is a bigger problem, because
    # the Gaia targets should all be unique. Let's start there
    print("Removing duplicate Gaia matches")
    to_delete = []
    for gid in dbl_gaia:
        loc = np.where(xmatch["GAIAEDR3_ID"]==gid)[0]
        # Just choosing the closest match would solve a lot of the issues
        # Especially since, in many cases, the TIC also seems to have this
        # issue crossmatching to multiple Gaia targets

        # If all the Gaia duplicates also have the same TIC ID, that's a
        # different problem. After some experimentation, I found two issues
        # 1) GES data was taken through two filters, so sometimes there are
        # duplicate lines of data. Info relevant for us (e.g., probability)
        # are the same across both rows.
        # 2) However, each filter row is also duplicated. I still don't know
        # why these lines were duplicated, but they weren't being removed in
        # table.unique() because GAIADR2_ID (and likely other columns) have
        # nans.
        # For now, I am just saving the first GES line in the table, and
        # deleting the rest.
        if len(np.unique(xmatch["TIC"][loc]))==1:
            # # print(xmatch["TIC","GAIAEDR3_ID","GAIAEDR3_G_CORRECTED","TIC_angDist",
            # #              "TIC_GAIADR2_ID","TIC_Gmag","TIC_Tmag"][loc])
            # print("multi-match for another reason")
            # # print(table.unique(xmatch[loc]))
            # # ucols = ["angDist_GES","target","filter","prob_p",
            # #          "angDist_Cantat-Gaudin","proba","GaiaDR2"]
            # ucols = ["angDist_GES","target","filter","prob_p","cluster",
            #          "S/N","RV_1","panstarrs1_GES","sdssdr13_GES","skymapper2_GES",
            #          "RV_2","p_filter","GAIADR2_ID","TIC"]
            #          # "",""]
            # tu = table.unique(xmatch[loc],keys=ucols)
            # print(tu[ucols])
            # ucols = ["angDist_GES","target","filter","prob_p","cluster",
            #          "S/N","RV_1","panstarrs1_GES","sdssdr13_GES","skymapper2_GES",
            #          "RV_2","p_filter"]
            #          # "",""]
            # tu = table.unique(xmatch[loc],keys=ucols)
            # print(tu[ucols])

            to_delete.append(loc[1:])


        else:
            sort_by_dist = np.argsort(xmatch["TIC_angDist"][loc])
            to_delete.append(loc[sort_by_dist[1:]])

    to_delete = np.concatenate(to_delete)
    xmatch.remove_rows(to_delete)
    print("deleted",len(to_delete))

    # Now re-do the repetition checks, and see if we fixed them
    unique_gaia, gct = np.unique(xmatch["GAIAEDR3_ID"], return_counts=True)
    dbl_gaia = unique_gaia[gct>1]
    dbl_gaia_idx = []
    for gid in dbl_gaia:
        loc = np.where(xmatch["GAIAEDR3_ID"]==gid)[0]
        dbl_gaia_idx.append(loc)
    if len(dbl_gaia_idx)>0:
        print("There are still duplicate EDR3 entries")
        dbl_gaia_idx = np.sort(np.concatenate(dbl_gaia_idx))
        print(dbl_gaia_idx)


    ## What to when two TIC IDs match the same Gaia target?
    print("Removing duplicate TIC matches")
    unique_tic, tct = np.unique(xmatch["TIC"], return_counts=True)
    dbl_tic = unique_tic[(tct>1)]
    dbl_tic_idx = []
    for tid in dbl_tic:
        loc = np.where(xmatch["TIC"]==tid)[0]
        dbl_tic_idx.append(loc)
    if len(dbl_tic_idx)>0:
        dbl_tic_idx = np.sort(np.concatenate(dbl_tic_idx))
        # print(dbl_tic_idx)


    dbl_tic = unique_tic[tct>1]
    to_delete = []
    for tic in dbl_tic:
        loc = np.where((xmatch["TIC"]==tic) & (xmatch["TIC"].mask==False))[0]
        if len(loc)==0:
            print("'duplicates' due to missing TIC matches")
            continue

        # If the TIC has a G mag, keep the Gaia entry with the closer Gmag
        if (np.any(xmatch["TIC_Gmag"].mask==False) and
            np.any(xmatch["GAIAEDR3_G_CORRECTED"].mask==False)):
            gdiff = abs(xmatch["GAIAEDR3_G_CORRECTED"][loc]-xmatch["TIC_Gmag"][loc])
            gsort_idx = np.argsort(gdiff)
            to_delete.append(loc[gsort_idx[1:]])

        # If the TIC doesn't have a G mag (only NGC 2547), just take the closer
        # match
        else:
            sort_by_dist = np.argsort(xmatch["TIC_angDist"][loc])
            to_delete.append(loc[sort_by_dist[1:]])

    if len(to_delete)>0:
        to_delete = np.concatenate(to_delete)
        xmatch.remove_rows(to_delete)
        print("deleted",len(to_delete))

        # Re-do TIC repetition checks
        unique_tic, tct = np.unique(xmatch["TIC"], return_counts=True)
        dbl_tic = unique_tic[(tct>1)]
        dbl_tic_idx = []
        for tid in dbl_tic:
            loc = np.where(xmatch["TIC"]==tid)[0]
            dbl_tic_idx.append(loc)
        if len(dbl_tic_idx)>0:
            dbl_tic_idx = np.sort(np.concatenate(dbl_tic_idx))
            dbl_tic_idx = np.intersect1d(dbl_tic_idx,
                                         np.where(xmatch["TIC"].mask==False)[0])
            if len(dbl_tic_idx)>0:
                print("There are still duplicate TIC IDs")
                print(dbl_tic_idx)

    ##### TODO: Add Phill's MINESweeper results


    ##### Check for proximity

    # Check for any nearby stars within ~30 arcsec (pixelish) -
    # likely source confusion.

    # If there is a blended star
        # If the other star has Q=2 or no measured period, do nothing

        # If the other star is >=3 mags fainter, do nothing

        # If I measured two different periods for the stars, do nothing

        # If I measured the same period for both stars
            # If one star has significantly higher periodogram power
            # (threshold?), choose that star

            # Otherwise, assign the period to the brighter star

    # Check for brighter stars within 1 arcmin - possible contamination
        # If I measured the same period as the brighter star, flag the fainter
        # star as having a bad period

        # Otherwise, just flag for possible contamination

    ##### Add literature periods

    # TODO


    ##### Clean up and output the catalog

    # Make sure every star has a value in the "final" columns

    # Select/create columns for output
    # TIC ID
    # Gaia data: EDR3 ID, DR2 ID?, RA, Dec, photometry, RUWE
    # 2MASS ID, J, H, Ks
    # HDBScan membership info
    # Jackson/GES membership info
    # Cantat-Gaudin membership info
    # My membership flag for plotting in this paper
    # Prot1, Pw1, Q1, Sig, Prot2, Pw2, Q2, MP?, SE?
    # Bl? (this will be an automated flag, unlike my previous papers)
    # Literature periods
    # Derived properties (masses, etc, from MINESweeper)

if __name__=="__main__":
    clusters = ["IC_2391","Collinder_135","NGC_2451A","NGC_2547","IC_2602"]
    dates = ["2021-06-22","2021-06-18","2021-06-21","2021-06-21","2021-07-02"]

    for i in range(5):
        print("\n\n",clusters[i])
        make_final_period_catalog(clusters[i],dates[i])
        # break
