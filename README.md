# Kernel-based confidence intervals for median absolute deviations

R code, input data and saved computational results for the manuscript by Hien Duy Nguyen and Luke A. Prendergast.

The estimators, interval interface, simulation drivers and numerical asset generation use native R. The core statistical calculations use base R; JSON reporting uses `jsonlite`. The recorded study used R 4.4.0 and jsonlite 1.8.8. The commands below run from the repository root.

## Use the interval routine

```sh
Rscript computation/R/mad_interval.R --data computation/data/galaxies.csv --column velocity_km_s
```

The routine returns JSON containing the sample MAD, local estimates, variance estimate, confidence intervals and a success or failure result. This command reports in the input's km/s units. The manuscript's galaxy illustration divides velocities by 1000.

## Check the code and saved results

```sh
Rscript scripts/verify_workspace.R
Rscript scripts/check.R NEW_VERIFICATION_DIRECTORY
```

The checks run in a temporary copy, preserving the archived reports. They check deterministic identities, numerical accuracy, interval and command-line contracts, all 30 principal design files, all 270 method summaries, 77 TM boundary refits and 4000 shoulder added-start records. They do not repeat the simulation study. Current results are in `verification/REPOSITORY_CHECKS.md`.

## Rebuild the numerical assets

```sh
Rscript scripts/reproduce_assets.R NEW_OUTPUT_DIRECTORY
```

This reads the completed records and creates 13 table bodies, three figures and their numerical coordinates in a fresh directory. It performs no new simulation, bootstrap resampling or distribution fitting. The original approximate Monte Carlo uncertainty bars are retained. Table filenames retain their study-release numbers; `OUTPUT_MAP.md` maps them to the November 15 manuscript.

For an independent study rerun, use a fresh directory and the settings in `computation/R/RUNNING.md`. Full runs are substantially more expensive than rebuilding the saved assets. The optional comparison with the native MADtests and gld packages remains unexecuted.

## Files

- `computation/R/`: estimation, inference, optimization, simulation and asset generation.
- `computation/data/`: 82 galaxy velocities and eight fixed benchmark vectors, with source notes.
- `computation/generated_r/full/`: compressed replication-level records, summaries, settings, seeds and diagnostics.
- `computation/generated_r/full/inference/`: simulation-median Monte Carlo standard errors.
- `computation/generated_r/full/finalize/assets/`: saved table bodies, figures, hybrid results and selected-fit summaries.
- `RELEASE.json` and `VALIDATION.md`: original study scope and qualifications.
- `MANIFEST_SHA256.json`: repository file hashes.

The principal study uses 1000 datasets per population/sample-size design, TM on the first 200 common datasets and 399 bootstrap resamples per dataset. The code retains failed fits and the declared coverage denominators. Population and optimizer checks describe the specified numerical procedures; they provide no global-optimum certificate or proof of empirical optimizer convergence. The galaxy dataset is a single illustration under an iid working model.
