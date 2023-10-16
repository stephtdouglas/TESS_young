# TESS_young

Code used to produce "Constraining Stellar Rotation at the Zero-Age Main Sequence with TESS" (Douglas et al., submitted to AAS Journals)

The code is licensed under the MIT License. Various catalogs are provided for reference; many of them are from the literature and I make no claim to their copyright. 

The organization of files evolved over the course of the project, and not everything was moved at the end. 
Carrying out the analysis in the paper involved roughly the following course

## Misc Utilities

Clone [k2spin][https://github.com/stephtdouglas/k2spin], then add to .bashrc or .bash_profile:
```export PYTHONPATH="$PATH:\<path to k2spin/\>"```

```tess_young/get_const.py``` produces various lists and plotting constants used in other scripts

```tess_young/plot_lit_ages.py``` plots ages found in the literature


## Catalog matching

Merge catalogs created via HDBScan and from the Gaia-ESO Survey (Jackson+2020) and Cantat-Gaudin+ (2020) 

Jackson+ and Cantat-Gaudin+ catalogs were matched to Gaia IDs using the CDS Xmatch service. 

```tess_young/match_catalogs.py``` merge catalogs

```tess_young/tess_sky_map.py``` plot cluster members in a stereographic projection

```tess_young/tess_sky_map.py``` plot cluster members in a stereographic projection

```scripts/plot_CMDs.py``` plot Gaia CMDs for cluster members


## Measuring rotation periods

All scripts were run on a high-performance computing cluster

```tess_young/lk_download.py``` downloads files from MAST. The search was done as a radius search and not as a TIC ID search, so some downloaded files don't correspond to any cluster member. 

```tess_young/measure_periods.py``` takes a list of downloaded light curve files and a cluster name, and runs a periodogram period search on all the input files

```tess_young/make_inspection_plots_tess.py``` Makes a plot showing the periodogram, light curve, and phase-folded light curve. It also shows the TESS pixel stamp and nominally Gaia sources, but the Gaia sources usually failed for some reason. 

At this, point, I went through every light curve and selected the best pipeline/sector for each target. I then did this again in a new file. These results are found under tables/\<cluster\>\_\<date\>\_results\_comments.csv and tables/\<cluster\>\_\<date\>\_results\_comments2.csv

```resolved_discrepant_validations.dat``` Most disagreements between validation steps could be resolved automatically in the catalog script below. This table contains a handful of stars that required manual resolution. 

```tess_young/make_inspection_examples.py``` Plots common light curve features

### Injection tests 

```scripts/setup_injection_tests.py```

```scripts/run_injection_tests.py```

```scripts/analyze_injection_tests.py```


## Literature rotation periods

```tess_young/match_lit_periods.py``` translates given identifiers for simbad compatibility, then find Gaia and TIC identifiers for each literature target

```table_lit_periods.py``` Merges tables of literature periods for different clusters (some manual editing was done afterwards to combine duplicate measurements)

```tess_young/plot_lit_periodmass.py``` Plot the literature rotation periods for cluster members

```tess_young/plot_lit_comparison.py``` Compare TESS periods to literature periods


## Final catalog

```tess_young/make_final_period_catalog.py``` Combines visual validation results with input catalogs, literature rotation periods, and manually resolved validation conflicts, then produces the final catalog tab\_all\_stars.csv. Also searches the Gaia catalog for potential blends, and removes duplicated sources. 

```tess_young/counting_stars.py``` Counts stars with various properties from the final catalog

## Tau-squared model fitting

Scripts are in the ```tess_young/tau_sq/``` folder, yaml input files are in the ```config/``` folder

```tess_young/tau_sq/tau_sq_run.py``` Runs the tau squared fits on the data and model specified in the yaml input file

```tess_young/tau_sq/tau_sq_synthetic_obs.py``` Generates synthetic data sets from the models, and runs tau_sq fits on those synthetic data

```tess_young/tau_sq/tau_sq_bootstrap.py``` Generates new data sets by bootstrap resampling within the inferred stellar mass uncertainties. 

```tess_young/tau_sq/tau_sq_uncertainties.py``` Takes the input from the synthetic or bootstrap runs, and produces uncertainties on the derived model ages

```scripts/plot_tausq_tracks_paper.py``` and ```plot_results_panel.py``` Create paper plots from the tausquared fit results

