import os, sys

import numpy as np
import astropy.io.fits as fits
import astropy.io.ascii as at
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.table import join,vstack,Table
from astroquery.mast import Catalogs


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
    # These have to be -9999 because comparison with NaNs gets messed up
    vis["second_period"] = np.ones_like(vis["final_period"])*-9999
    vis["second_power"] = np.ones_like(vis["final_power"])*-9999
    vis["second_Q"] = np.ones_like(vis["final_Q"])*9

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

    # q2 = ["Q2","Q3","Q3"] # What was this for?

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
    cat2 = read_validation_results(cluster,date,which=2)
    print(len(np.unique(cat1["TIC"])),"Unique TIC IDs, validation 1")
    print(len(np.unique(cat2["TIC"])),"Unique TIC IDs, validation 2")

    # Crossmatch the validation catalogs on TIC IDs
    allcat = join(cat1,cat2,keys=["TIC"])
    print(len(allcat),"TESS results")
    print(len(np.unique(allcat["TIC"])),"Unique TIC IDs")

    # Create new columns to hold the final-final values for
    # period, Q, light curve info, peak power, threshold, multi/spot flags
    # Colnames from D19:
    # Prot1, Pw1, Q1, Sig, Prot2, Pw2, Q2, MP?, SE?
    # Bl? (this will be an automated flag, unlike my previous papers)
    new_cols_float = [np.ones(len(allcat))*-9999 for i in range(5)]
    new_cols_int = [np.ones(len(allcat))*8 for i in range(2)]
    new_cols_char = [np.zeros(len(allcat),"U1") for i in range(3)]
    allcat.add_columns(new_cols_float,names=["Prot1", "Pw1", "Sig",
                                             "Prot2", "Pw2"])
    allcat.add_columns(new_cols_int,names=["Q1", "Q2"])
    allcat.add_columns(new_cols_char,names=["MP?","SE?","Bl?"])
    # Apparently I did not include my notes in the output files for D19
    # So I'm leaving them out for simplicity now
    allcat.add_column(np.zeros(len(allcat),"U5"),name="provenance_name")
    allcat.add_column(np.zeros(len(allcat),"U10"),name="flux_cols")
    # allcat.add_column(np.zeros(len(allcat),int),name="sequence_number")
    allcat.add_column(np.ones(len(allcat),"int")*-9999,name="sequence_number")
    allcat.add_column(np.zeros(len(allcat),"U200"),name="obs_id")
    allcat["obs_id"][:] = ""

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
        allcat[output_cols2[i]][good2] = allcat[init_cols2[i]][good2]

    ## Cycle through all the stars with differing results
    diff_idx = np.where(diff_results)[0]
    print(len(diff_idx))

    primary_still_bad = []
    for i in diff_idx:
        diff_frac = (abs(allcat["final_period_1"][i]-allcat["final_period_2"][i])
                     / min(allcat["final_period_1","final_period_2"][i]))

        half_dbl = allcat["final_period_1"][i]/allcat["final_period_2"][i]

        # If one was flagged as Q=2, flag both as Q=2
        # AND remove any secondary period
        if (allcat["final_Q_1"][i]==2) or (allcat["final_Q_2"][i]==2):
            allcat["Q1"][i] = 2
            allcat["Prot1"][i] = -9999
            allcat["Pw1"][i] = -9999
            allcat["provenance_name"][i] = allcat["provenance_name_1"][i]
            allcat["flux_cols"][i] = allcat["flux_cols_1"][i]
            allcat["sequence_number"][i] = allcat["sequence_number_1"][i]
            allcat["obs_id"][i] = allcat["obs_id_1"][i]
            allcat["Sig"][i] = allcat["thresholds_1"][i]
            allcat["Q2"][i] = 2
            allcat["Prot2"][i] = -9999
            allcat["Pw2"][i] = -9999

        # If both are -9999, then no need to do anything
        elif ((allcat["final_period_1"][i]<0) and (allcat["final_period_2"][i]<0)):
            continue

        # If the period is (almost) the same (and Q=1 or Q=0):
        # elif allcat["final_period_1"][i]==allcat["final_period_2"][i]:
        elif diff_frac<0.05:

            # TODO: but it's just a different light curve, pick the one with the
            # higher periodogram peak and assign the lower Q value
            # higher_peak = np.argmax(allcat["final_"])

            # If I assigned different Q values, assign the lower Q value
            allcat["Q1"][i] = min(allcat["final_Q_1","final_Q_2"][i])+20
            allcat["Prot1"][i] = allcat["final_period_1"][i]
            allcat["Pw1"][i] = allcat["final_power_1"][i]
            allcat["provenance_name"][i] = allcat["provenance_name_1"][i]
            allcat["flux_cols"][i] = allcat["flux_cols_1"][i]
            allcat["sequence_number"][i] = allcat["sequence_number_1"][i]
            allcat["obs_id"][i] = allcat["obs_id_1"][i]
            allcat["Sig"][i] = allcat["thresholds_1"][i]

        # If they're harmonics and the second period is longer, choose the 2nd
        elif (abs(half_dbl-0.5)<0.02):
            allcat["Q1"][i] = allcat["final_Q_2"][i]+30
            allcat["Prot1"][i] = allcat["final_period_2"][i]
            allcat["Pw1"][i] = allcat["final_power_2"][i]
            allcat["provenance_name"][i] = allcat["provenance_name_2"][i]
            allcat["flux_cols"][i] = allcat["flux_cols_2"][i]
            allcat["sequence_number"][i] = allcat["sequence_number_2"][i]
            allcat["obs_id"][i] = allcat["obs_id_2"][i]
            allcat["Sig"][i] = allcat["thresholds_2"][i]

        # If they're harmonics and the first period is longer, choose the 1st
        elif (abs(half_dbl-2)<0.08):
            allcat["Q1"][i] = allcat["final_Q_1"][i]+40
            allcat["Prot1"][i] = allcat["final_period_1"][i]
            allcat["Pw1"][i] = allcat["final_power_1"][i]
            allcat["provenance_name"][i] = allcat["provenance_name_1"][i]
            allcat["flux_cols"][i] = allcat["flux_cols_1"][i]
            allcat["sequence_number"][i] = allcat["sequence_number_1"][i]
            allcat["obs_id"][i] = allcat["obs_id_1"][i]
            allcat["Sig"][i] = allcat["thresholds_1"][i]

        # If the periods are not the same...
        # so far most of these look impossibly confused, though a few
        # could be salvageable. I have a list of stars
        # that I've resolved by hand
        else:
            print("\n",allcat["TIC"][i],"primary mismatch")
            # print(allcat["TIC","provenance_name_1","flux_cols_1",
            #              "sequence_number_1","final_period_1","final_Q_1",
            #              "Notes_1"][i])
            # print(allcat["TIC","provenance_name_2","flux_cols_2",
            #              "sequence_number_2","final_period_2","final_Q_2",
                         # "Notes_2"][i])
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
            allcat["Prot2"][i] = -9999
            allcat["Pw2"][i] = -9999

        # If both are -9999, then no need to do anything
        elif ((allcat["second_period_1"][i]<0) and (allcat["second_period_2"][i]<0)):
            continue

        # If the period is (almost) the same (and Q=1 or Q=0):
        elif diff_frac<0.05:

            # TODO: but it's just a different light curve, pick the one with the
            # higher periodogram peak and assign the lower Q value
            # higher_peak = np.argmax(allcat["final_"])

            # TODO: Could this be where the random 0/1 values are being included? If the issue is two -9999 values?

            # If I assigned different Q values, assign the lower Q value
            allcat["Q2"][i] = min(allcat["second_Q_1","second_Q_2"][i])+20
            allcat["Prot2"][i] = allcat["second_period_1"][i]
            allcat["Pw2"][i] = allcat["second_power_1"][i]

        # If they're harmonics and the second period is longer, choose the 2nd
        elif (abs(half_dbl-0.5)<0.02):
            allcat["Q2"][i] = allcat["second_Q_2"][i]+30
            allcat["Prot2"][i] = allcat["second_period_2"][i]
            allcat["Pw2"][i] = allcat["second_power_2"][i]

        # If they're harmonics and the first period is longer, choose the 1st
        elif (abs(half_dbl-2)<0.08):
            allcat["Q2"][i] = allcat["second_Q_1"][i]+40
            allcat["Prot2"][i] = allcat["second_period_1"][i]
            allcat["Pw2"][i] = allcat["second_power_1"][i]

        # If the periods are not the same...
        # resolve by hand
        else:
            print("\n",allcat["TIC"][i],"secondary mismatch")
            # print(allcat["TIC","provenance_name_1","flux_cols_1",
            #              "sequence_number_1","final_period_1","final_Q_1",
            #              "second_period_1","second_Q_1"][i])
            # print(allcat["TIC","provenance_name_2","flux_cols_2",
            #              "sequence_number_2","final_period_2","final_Q_2",
            #              "second_period_2","second_Q_2"][i])
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
            if resolved["second_Q"][loc]<=4:
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

    for colname in ["Sig","Pw1","Pw2"]:
        periods[colname].info.format = ".3f"
    for colname in ["Prot1",  "Prot2"]:
        periods[colname].info.format = ".2f"

    ### TODO: Removed to do error tracing - remember to re-add
    # bad_prot = periods["Q1"]>=2
    # periods["Prot1"][bad_prot] = -9999
    # periods["Pw1"][bad_prot] = -9999
    # periods["Prot2"][bad_prot] = -9999
    # periods["Pw2"][bad_prot] = -9999
    # replace_Q2 = bad_prot & (periods["Q2"]<2)
    # periods["Q2"][replace_Q2] = 2

    # bad_prot2 = periods["Q2"]>=2
    # periods["Prot2"][bad_prot2] = -9999
    # periods["Pw2"][bad_prot2] = -9999

    # # Write out the table. I'm not going to write a tex table - I'll just put
    # # the colnames in the manuscript proper
    # if os.path.exists(f"tables/tab_{cluster}_tess_periods.csv"):
    #     print("WARNING: Period table already exists; not re-recreating it")
    # else:
    #     at.write(periods,f"tables/tab_{cluster}_tess_periods.csv",delimiter=",",
    #              formats=formats)

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
    print(len(tic_match),"xmatch_TIC")

    ##########################################################################
    ##########################################################################
    ##########################################################################

    # Then match back to my original catalogs
    # Phill's MINESweeper catalogs only included the GAIA IDs
    cat_init_file1 = f"tables/{cluster}_crossmatch_allcolumns.fits"
    with fits.open(cat_init_file1) as hdu:
        cat_init1 = Table(hdu[1].data)
        # print(cat_init1.dtype)

    cat_init_file2 = os.path.expanduser(f"~/Dropbox/data/MINESweeper/catalog_{cluster}_v0.fits")
    with fits.open(cat_init_file2) as hdu:
        cat_init2 = Table(hdu[1].data)
        # print(cat_init.dtype)
        cat_init2.rename_column("GaiaEDR3_ID","GAIAEDR3_ID")

    cat_init = join(cat_init1,cat_init2,keys="GAIAEDR3_ID")
    print(len(cat_init1),len(cat_init2),"minesweeper inputs")

    if cat_init.masked == False:
        cat_init = Table(cat_init, masked=True, copy=False)

    print(len(cat_init),"input crossmatch/MINESweeper")

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
    print(len(xmatch),"input crossmatch with TIC")

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

    # If there's a star without a TIC ID, introduce a placeholder
    no_tic = xmatch["TIC"].mask==True
    xmatch["TIC"][no_tic] = -9999
    xmatch["TIC"][no_tic].mask = False

    # Join the final period and membership catalogs
    pmcat = join(xmatch,periods,keys=["TIC"],join_type='left')
    print(len(xmatch),len(periods))
    print(len(pmcat),"periods and membership")
    has_periods = np.where((pmcat["Prot1"]>0) & (pmcat["Prot1"].mask==False))[0]
    print(len(has_periods),"periods and membership: actual periods")
    print(np.unique(pmcat["Q1"]))
    unmask = pmcat["Q1"].mask==True
    pmcat["Q1"][unmask] = 7
    pmcat["Q1"].mask[unmask] = False

    ##########################################################################
    ##########################################################################
    ##########################################################################

    ##### Check for proximity with other *cluster* stars

    managed = np.zeros(len(pmcat),bool)
    ppos = SkyCoord(pmcat["GAIAEDR3_RA"],pmcat["GAIAEDR3_DEC"],unit=u.degree)

    neighbor_ct = 0
    removed_ct = 0
    for i,row in enumerate(pmcat):

        if (managed[i]==True) | (row["Q1"]>=2):
            continue

        # Check for any nearby *cluster* stars within ~30 arcsec (pixelish) -
        # likely source confusion.
        sep = ppos[i].separation(ppos)
        pixelish = 30*u.arcsec
        close0 = np.where(sep<pixelish)[0]

        # bright_contam = 1*u.arcmin
        # close0 = np.where(sep<bright_contam)[0]

        if len(close0)==1:
            managed[i] = True
            continue

        # If there is a blended star
        else:
            neighbor_ct += 1
            primary_idx = np.where(close0==i)[0]
            close = np.delete(close0,primary_idx)
            # If the other star has Q=2 or no measured period, do nothing
            neighbors_no_prot = ((pmcat["Q1"][close]>=2))
            # If the other star is >=3 mags fainter, do nothing
            neighbors_faint = ((pmcat["TIC_Tmag"][i]-pmcat["TIC_Tmag"][close])>=-3)

            good_prot = (np.isfinite(pmcat["Prot1"][close]) &
                         (pmcat["Prot1"][close]>0) &
                         (pmcat["Q1"][close]<=1))
            diff_frac = (abs(pmcat["Prot1"][i]-pmcat["Prot1"][close])
                         / pmcat["Prot1"][i])
            # If I measured two different periods for the stars, do nothing
            # neighbors_diff_prot = ((pmcat["Prot1"][close] -
            #                        pmcat["Prot1"][i])!=0)
            neighbors_diff_prot = (diff_frac>=0.05) & good_prot
            same_period = (diff_frac<0.05) & good_prot

            if np.all(neighbors_no_prot | neighbors_diff_prot |
                      neighbors_faint | (managed[close]==True)):
                  continue

            # If I measured the same period for both stars
            elif np.any(same_period):
                # If one star has significantly higher periodogram power
                # (threshold?), choose that star
                # Otherwise, assign the period to the brighter star

                # same period, lower power = remove those periods
                power_diff = pmcat["Pw1"][i]-pmcat["Pw1"][close]
                # Only consider stars with a significantly different power
                higher_power = (power_diff<-0.1) & np.isfinite(pmcat["Pw1"][close])
                # Only consider brighter stars that are at least 1 magnitude brighter
                brighter = (pmcat["TIC_Tmag"][i]-pmcat["TIC_Tmag"][close])>1

                lower_power = (power_diff>0.1) & np.isfinite(pmcat["Pw1"][close])
                same_prot_lower_power = same_period & lower_power
                # print(same_period)
                # print(power_diff>0)
                # print(close)
                if np.any(same_period & (higher_power | brighter)):
                    print("\nSame period, neighbor is brighter/higher power")
                    print(pmcat["TIC","Prot1","Pw1","Q1","TIC_Tmag"][i])
                    print(pmcat["TIC","Prot1","Pw1","Q1","TIC_Tmag"][close[same_period & (higher_power | brighter)]])
                    print(sep.arcminute[close[same_period & (higher_power | brighter)]])
                    pmcat["Prot1"][i] = -9999 # np.nan
                    pmcat["Pw1"][i] = -9999 # np.nan
                    pmcat["Q1"][i] = 5
                    pmcat["Prot2"][i] = -9999 # np.nan
                    pmcat["Pw2"][i] = -9999 # np.nan
                    pmcat["Q2"][i] = 5
                    managed[i] = True
                    removed_ct += 1
                else:
                    print("\nSame period, neighbor is fainter")
                    print(pmcat["TIC","Prot1","Pw1","Q1","TIC_Tmag"][i])
                    print(pmcat["TIC","Prot1","Pw1","Q1","TIC_Tmag"][close[same_prot_lower_power]])
                    print(sep.arcminute[close[same_prot_lower_power]])
                    pmcat["Prot1"][close[same_prot_lower_power]] = -9999 # np.nan
                    pmcat["Pw1"][close[same_prot_lower_power]] = -9999 # np.nan
                    pmcat["Q1"][close[same_prot_lower_power]] = 5
                    pmcat["Prot2"][close[same_prot_lower_power]] = -9999 # np.nan
                    pmcat["Pw2"][close[same_prot_lower_power]] = -9999 # np.nan
                    pmcat["Q2"][close[same_prot_lower_power]] = 5
                    managed[close[same_prot_lower_power]] = True
                    removed_ct += len(np.where(same_prot_lower_power)[0])

    print(neighbor_ct," stars with neighbors")
    print(removed_ct," blended stars removed")

    ##########################################################################
    ##########################################################################
    ##########################################################################

    ##### Check for proximity with any other stars in the Gaia catalog

    # Read in the full regional catalogs from Phill
    cat_file0 = os.path.expanduser(f"~/Dropbox/EDR3/scats/{cluster}.fits")
    with fits.open(cat_file0) as hdu:
        cat0 = Table(hdu[1].data)
        cat0.rename_column("GAIADER3_ID","GAIAEDR3_ID")
        # print(cat0.dtype)

    if cat0.masked == False:
        cat0 = Table(cat0, masked=True, copy=False)

    gpos = SkyCoord(cat0["GAIAEDR3_RA"],cat0["GAIAEDR3_DEC"],unit=u.degree)

    pmcat["Bl?"] = np.empty(len(pmcat),"U1")
    pmcat["Bl?"][:] = "n"

    pmcat["ClosestNeighbor"] = np.zeros(len(pmcat))
    for i,row in enumerate(pmcat):

        sep = ppos[i].separation(gpos)
        pixelish = 30*u.arcsec
        closep = np.where(sep<pixelish)[0]

        bright_contam = 1*u.arcmin
        closeb = np.where(sep<bright_contam)[0]

        if (len(closep)==1) and (len(closeb)==1):
            self_i = np.where(cat0["GAIAEDR3_ID"]==pmcat["GAIAEDR3_ID"][i])[0]
            sep = np.delete(sep,self_i)
            pmcat["ClosestNeighbor"][i] = np.min(sep).to(u.arcsec).value
            continue

        else:

            tic_data = Catalogs.query_region(ppos[i],catalog="TIC",
                                             radius=bright_contam)
            # Check for neighbors that are w/in 3 mag of the target, or brighter
            # If the target has T=10,
            # and the neighbor has T=12, then the difference is -2
            #
            neighbors_faint = ((tic_data["Tmag"][0]-tic_data["Tmag"])>-3)
            # Then check for much brighter neighbors
            neighbors_bright = tic_data["Tmag"]<=12
            # The target itself is not a neighbor
            # Is the target itself always the first line? It seems to be.
            self_i = np.argmin(tic_data["dstArcSec"])
            tic2 = np.delete(tic_data,self_i)
            pmcat["ClosestNeighbor"][i] = np.min(tic2["dstArcSec"])
            neighbors_faint[self_i] = False
            neighbors_bright[self_i] = False

            closept = tic_data["dstArcSec"]<=pixelish.to(u.arcsec).value
            closebt = tic_data["dstArcSec"]<=bright_contam.to(u.arcsec).value

            ncpt = len(np.where(neighbors_faint & closept)[0])
            ncbt = len(np.where(neighbors_bright & closebt)[0])

            if ncpt>0:
                pmcat["Bl?"][i] = "y"
            elif ncbt>0:
                pmcat["Bl?"][i] = "m"


    ##########################################################################
    ##### Add literature periods

    lit = at.read("tab_lit_periods_consolidated.csv",delimiter=",")
    lit.rename_column("Source","LitSource")

    lit2 = lit[(lit["TIC"]>0) & (lit["LitPeriod"].mask==False)]

    pmcat2 = join(pmcat,lit2["TIC","LitPeriod","LitSource"],keys=["TIC"],
                  join_type="left")
    print(len(np.where(pmcat2["LitPeriod"].mask==False)[0]),"Lit periods added")
    print(len(pmcat2),"after xmatch with lit")

    ##### Clean up and output the catalog

    # Make sure every star has a value in the "final" columns

    # Select/create columns for output
    if pmcat2.masked == False:
        pmcat2 = Table(pmcat2, masked=True, copy=False)

    if cluster!="Collinder_135":
        out_cat = pmcat2["TIC","GAIAEDR3_ID","GAIAEDR3_RA","GAIAEDR3_DEC",
                         "GAIAEDR3_PMRA","GAIAEDR3_PMDEC",
                         "GAIAEDR3_PARALLAX","GAIAEDR3_PARALLAX_CORRECTED",
                         "GAIAEDR3_RUWE","GAIAEDR3_G","GAIAEDR3_G_ERR",
                         "GAIAEDR3_G_CORRECTED","GAIAEDR3_BP","GAIAEDR3_BP_ERR",
                         "GAIAEDR3_RP","GAIAEDR3_RP_ERR",
                         "TMASS_ID","TMASS_J","TMASS_J_ERR",
                         "TMASS_H","TMASS_H_ERR","TMASS_K","TMASS_K_ERR",
                         "HDBscan_MemProb","HDBscan_Cluster","HDBscan_Stability",
                         "MemBool",
                         "angDist_GES","target","cluster","prob_p",
                         "angDist_Cantat-Gaudin","proba","Cluster",
                         "av","av_err","dist","dist_err","log(Age)",
                         "log(Age)_err","Mass","Mass_err","log(Teff)",
                         "log(Teff)_err",
                         "Prot1", "Pw1", "Q1", "Sig", "Prot2", "Pw2",
                         "Q2", "MP?", "SE?","Bl?","ClosestNeighbor",
                         # blend??
                         "LitPeriod","LitSource"
                         ]
        out_cat.rename_column("target","GES_Target")
        out_cat.rename_column("cluster","GES_Cluster")
        out_cat.rename_column("prob_p","GES_MemProb")
    else:
        pmcat2.add_column(np.ones(len(pmcat2))*-9999,name="angDist_GES")
        pmcat2.add_column(np.ones(len(pmcat2))*-9999,name="GES_MemProb")
        pmcat2.add_column(np.ones(len(pmcat2),"U16"),name="GES_Target")
        pmcat2.add_column(np.ones(len(pmcat2),"U10"),name="GES_Cluster")
        pmcat2["GES_Target"][:] = ""
        pmcat2["GES_Cluster"][:] = ""
        pmcat2["GES_Target"][:].mask = True
        pmcat2["GES_Cluster"][:].mask = True
        out_cat = pmcat2["TIC","GAIAEDR3_ID","GAIAEDR3_RA","GAIAEDR3_DEC",
                         "GAIAEDR3_PMRA","GAIAEDR3_PMDEC",
                         "GAIAEDR3_PARALLAX","GAIAEDR3_PARALLAX_CORRECTED",
                         "GAIAEDR3_RUWE","GAIAEDR3_G","GAIAEDR3_G_ERR",
                         "GAIAEDR3_G_CORRECTED","GAIAEDR3_BP","GAIAEDR3_BP_ERR",
                         "GAIAEDR3_RP","GAIAEDR3_RP_ERR",
                         "TMASS_ID","TMASS_J","TMASS_J_ERR",
                         "TMASS_H","TMASS_H_ERR","TMASS_K","TMASS_K_ERR",
                         "HDBscan_MemProb","HDBscan_Cluster","HDBscan_Stability",
                         "MemBool",
                         "angDist_GES","GES_MemProb","GES_Target","GES_Cluster",
                         "angDist_Cantat-Gaudin","proba","Cluster",
                         "av","av_err","dist","dist_err","log(Age)",
                         "log(Age)_err","Mass","Mass_err","log(Teff)",
                         "log(Teff)_err",
                         "Prot1", "Pw1", "Q1", "Sig", "Prot2", "Pw2",
                         "Q2", "MP?", "SE?","Bl?","ClosestNeighbor",
                         # blend??
                         "LitPeriod","LitSource"
                         ]

    out_cat.rename_column("proba","CG_MemProb")
    out_cat.rename_column("Cluster","CG_Cluster")

    out_cat.add_column(np.empty(len(out_cat),"U18"),name="Cluster")
    out_cat["Cluster"][:] = cluster
    periods.add_column(np.empty(len(periods),"U18"),name="Cluster")
    periods["Cluster"][:] = cluster


    # Membership filter for plotting
    out_cat.add_column(np.zeros(len(out_cat),int),name="to_plot")
    hdb_memb = (out_cat["MemBool"]==1)  & (out_cat["MemBool"].mask==False)

    # Jackson+2020 Table 4/Section 4 indictes that 0.9 is the membership cutoff
    ges_memb = (((out_cat["GES_MemProb"]>=0.9) & (out_cat["GES_MemProb"]<=1))
                & (out_cat["GES_MemProb"].mask==False))

    can_memb = (((out_cat["CG_MemProb"]>=0.7) & (out_cat["CG_MemProb"]<=1))
                & (out_cat["CG_MemProb"].mask==False))

    hdb_memb1 = np.asarray(hdb_memb,"int")
    ges_memb1 = np.asarray(ges_memb,"int")
    can_memb1 = np.asarray(can_memb,"int")
    sum_memb = hdb_memb1 + ges_memb1 + can_memb1

    for i in range(4):
        print(len(np.where(sum_memb==i)[0]))

    if cluster=="Collinder_135":
        out_cat["to_plot"][hdb_memb | can_memb] = 1
    else:
        out_cat["to_plot"][sum_memb>=2] = 1
    # TODO: I think this is missing some stars! Need to check later
    print("\n",len(np.where(out_cat["to_plot"]==1)[0]),"stars to plot")
    print(len(np.where((out_cat["to_plot"]==1) & (out_cat["Prot1"]>0) &
                       (out_cat["Prot1"].mask==False)
                       )[0]),"with periods")

    # TIC ID
    # Gaia data: EDR3 ID, DR2 ID?, RA, Dec, photometry, RUWE
    out_cat["GAIAEDR3_RA"].info.format = ".6f"
    out_cat["GAIAEDR3_DEC"].info.format = ".6f"
    # 2MASS ID, J, H, Ks
    # HDBScan membership info
    # Jackson/GES membership info
    # Cantat-Gaudin membership info
    # My membership flag for plotting in this paper
    # Prot1, Pw1, Q1, Sig, Prot2, Pw2, Q2, MP?, SE?
    # Literature periods
    for colname in ["GAIAEDR3_RUWE","GAIAEDR3_RA","GAIAEDR3_DEC",
                     "GAIAEDR3_RUWE","GAIAEDR3_G","GAIAEDR3_G_ERR",
                     "GAIAEDR3_G_CORRECTED","GAIAEDR3_BP","GAIAEDR3_BP_ERR",
                     "GAIAEDR3_RP","GAIAEDR3_RP_ERR",
                     "TMASS_ID","TMASS_J","TMASS_J_ERR",
                     "TMASS_H","TMASS_H_ERR","TMASS_K","TMASS_K_ERR",
                     "Sig","Pw1","Pw2"]:
        out_cat[colname].info.format = ".3f"

    for colname in ["Prot1",  "Prot2", "LitPeriod","ClosestNeighbor"]:
        out_cat[colname].info.format = ".2f"

    return periods, out_cat


if __name__=="__main__":
    clusters = ["IC_2391","Collinder_135","NGC_2451A","NGC_2547","IC_2602"]
    dates = ["2021-06-22","2021-06-18","2021-06-21","2021-06-21","2021-07-02"]

    period_cats = []
    out_cats = []
    for i in range(5):
        print("\n\n",clusters[i])
        periods, out_cat = make_final_period_catalog(clusters[i],dates[i])
        period_cats.append(periods)
        out_cats.append(out_cat)
        # break

    all_periods = vstack(period_cats)
    all_outputs = vstack(out_cats)

    at.write(all_outputs,"tab_all_stars.csv", delimiter=",",overwrite=True)
    at.write(all_periods,"tab_all_tess_periods.csv", delimiter=",",overwrite=True)

    with open("tab_cols_all_stars_blank.dat","w") as f:
        for i,colname in enumerate(all_outputs.dtype.names):
            f.write(f"{i+1} & {colname} & \\\\\n")
