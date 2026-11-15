# Running the R computation

Use R 4.4.0 with jsonlite 1.8.8. Core methods use base R. `MAD.Rproj` opens the project in RStudio. Run commands from the project root.

## Entry points

```sh
Rscript computation/R/mad_interval.R --data computation/data/galaxies.csv --column velocity_km_s
Rscript computation/R/test_core.R
Rscript computation/tests/test_direct_interface.R
Rscript scripts/verify_workspace.R
```

The interval CLI accepts a CSV column, bandwidth rule, confidence level and reporting scale. It returns strict JSON, including an explicit success flag and failure reason. The study illustration divides galaxy velocities by 1000 internally; the CLI above reports in the input's km/s units.

## Study reproduction

The completed results are already under `computation/generated_r/full/`. Do not repeat them for manuscript or packaging work. To reproduce the entire study independently, use a fresh output directory:

```sh
Rscript computation/R/run.R all --mode full --out computation/reproduced_R
```

Full tasks are principal, population, galaxy, clt, large, resolution, restarts, tm, audit, summary-checks, diagnostic and finalize. Full principal settings are 1000 datasets per design, TM on the first 200, and 399 bootstrap resamples. The CLT task uses 3000 datasets per design; larger-sample, resolution and restart tasks use 1000. A smoke mode uses small counts and cannot validate full coverage.

`run.R finalize` derives all table rows and figures from the completed R records. It checks required counts and shared shoulder datasets and refuses missing or partial inputs. The 13 table bodies are in `full/finalize/assets/tables/`; the three figures and their numeric coordinates are in the corresponding assets directory. The manuscript uses these completed outputs. This repository distributes the numerical assets without a manuscript build. Use `Rscript scripts/reproduce_assets.R NEW_DIRECTORY` to rebuild them from saved records.

## RNG and checkpoints

The master seed is 20260917. All experiments use L'Ecuyer-CMRG, Inversion and Rejection. The seed begins at the master value and folds, in order, experiment, zero-based population index, sample size and zero-based replicate using `h <- (48271*h + key + 1) %% 2147483647`; zero is replaced by one. Bootstrap draws use their own experiment key. Population global search uses the documented scalar seed 783 plus population index. Simulation-median MCSEs used diagnostic seed 20261025 and 400 resamples.

Each local checkpoint is an atomic RDS value with its exact configuration. Source MD5s, numerical settings, RNG and R/platform information are recorded with completed runs. An existing checkpoint is reused only when its configuration is identical. The completed records retain their original execution hashes. Reproduction with the current source should use a fresh directory rather than bypassing configuration guards.

`launch_principal_full.R` supports disjoint sharded execution and merging. A coordinator owns the shard root; do not start a serial principal run in that directory. The archived release omits local checkpoints and process markers. Use the saved-result counts and repository verification report for the completed study; a fresh run writes its own local checkpoints.

## Diagnostics and assets

- `optimizer_checks.R` contains the added-start, wider-bound TM, population and multistart checks. Failed attempts remain in coverage denominators.
- `check_failed_fits.R` reproduces the two retained convergence-limit examples when explicitly needed; the completed diagnostics are already saved.
- `normal_bandwidth_diagnostic.R` computes the first-order Gaussian linear diagnostic using R integration; it is not an exact finite-sample variance or a correction to coverage.
- `benchmark_reference.R` times the eight fixed vectors under `computation/data/benchmark/`. To rerun timings, supply a new output directory as its second argument; hardware and concurrent work affect results.
- `compare_native_packages.R` is optional and unexecuted. It requires MADtests and gld in addition to the ordinary project dependencies.

Current scientific limitations and actual executed validation are in `VALIDATION.md`. No global-optimization or native-package-equivalence claim is made.
