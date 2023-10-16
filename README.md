# TESS_young

Code used to produce "Constraining Stellar Rotation at the Zero-Age Main Sequence with TESS" (Douglas et al., submitted to AAS Journals)

The code is licensed under the MIT License. Various catalogs are provided for reference; many of them are from the literature and I make no claim to their copyright. 

The organization of files evolved over the course of the project, and not everything was moved at the end. 
Carrying out the analysis in the paper involved roughly the following course

## Misc Utilities

```get_const.py``` produces various lists and plotting constants used in other scripts

Clone [k2spin][https://github.com/stephtdouglas/k2spin], then add to .bashrc or .bash_profile:
```export PYTHONPATH="$PATH:\<path to k2spin/\>"```


## Catalog matching

Merge catalogs created via HDBScan and from the Gaia-ESO Survey (Jackson+2020) and Cantat-Gaudin+ (2020) 

Jackson+ and Cantat-Gaudin+ catalogs were matched to Gaia IDs using the CDS Xmatch service. 

```tess_young/match_catalogs.py``` merge catalogs

```tess_young/tess_sky_map.py``` plot cluster members in a stereographic projection

## Measuring rotation periods

All scripts were run on a high-performance computing cluster

```tess_young/lk_download.py``` downloads files from MAST. The search was done as a radius search and not as a TIC ID search, so some downloaded files don't correspond to any cluster member. 

```tess_young/measure_periods.py``` takes a list of downloaded light curve files and a cluster name, and runs a periodogram period search on all the input files

```tess_young/make_inspection_plots_tess.py``` Makes a plot showing the periodogram, light curve, and phase-folded light curve. It also shows the TESS pixel stamp and nominally Gaia sources, but the Gaia sources usually failed for some reason. 

At this, point, I went through every light curve and selected the best pipeline/sector for each target. I then did this again in a new file. These results are found under tables/\<cluster\>_\<date\>_results_comments.csv and tables/\<cluster\>_\<date\>_results_comments2.csv

```resolved_discrepant_validations.dat``` Most disagreements between validation steps could be resolved automatically in the catalog script below. This table contains a handful of stars that required manual resolution. 

## Literature rotation periods

```tess_young/match_lit_periods.py``` translate given identifiers for simbad compatibility, then find Gaia and TIC identifiers for each literature target

## Final catalog

```tess_young/make_final_period_catalog.py```


## Tau-squared model fitting

In the ```tess_young/tau_sq/``` folder


