# Archived study validation

This file records validation of the original October 2026 R study. The manuscript paths and manuscript-build claims below are historical; no manuscript is distributed in this repository. Current repository checks are described in README.md and verification/REPOSITORY_CHECKS.md. Scientific calculations and saved study results are retained.

# R computation and validation

The manuscript's 13 tables and three figures use the completed R studies. R 4.4.0 and jsonlite 1.8.8 were used. The statistical implementation uses base R.

## Completed study scope

| Study | Executed scope |
|---|---|
| Principal | Ten populations × three sample sizes × 1000 datasets = 30,000 |
| Bootstrap | 399 resamples per principal dataset = 11,970,000 resamples |
| TM | First 200 datasets per principal design = 6000 attempts |
| Two samples | Four scenarios and seven variance procedures, derived from common principal records |
| Population calibration | Ten populations × J=25,50, ten starts per fit |
| Variance-estimator limit | Two populations × three sizes × 3000 = 18,000 datasets |
| Larger samples | Five populations × n=5000 × 1000 = 5000 datasets |
| Mixture resolution | Shoulder and bimodal × n=50,000 × 1000 = 2000 datasets |
| Added shoulder starts | Four sizes × 1000 common datasets = 4000 |
| TM sensitivity | All 77 flagged exponential/Pareto fits |
| Optimizer audit | Three population searches and 40 datasets with ten starts |
| Galaxy illustration | All six methods on the unchanged 82 observations |
| Gaussian diagnostic | Deterministic integration at n=200,1000,5000 |

The principal study completed at 05:08:41 UTC on 25 October 2026 and the TM refits at 05:12:19 UTC.  Full R asset assembly verified all design counts, method sets and the 4000 shared shoulder datasets, with zero MAD, variance-ratio or coverage disagreements in that dataset-reuse check.

## R checks executed

The following deterministic checks were executed:

```sh
Rscript computation/R/test_core.R
Rscript computation/tests/test_direct_interface.R
```

The 70 interval and CLI checks pass. They exercise finite endpoints, reporting scale, positive affine equivariance, reflection away from endpoint ties, invalid samples, zero MAD, overflow and underflow, strict JSON failure output, and the actual R command-line process. Deterministic core checks pass with maximum FKML quantile/CDF discrepancy 6.28e-12, scaled variance-gradient discrepancy 7.15e-10, population defining-equation residual 3.67e-15, and optimizer quadratic-target parameter error 2.85e-10. Re-evaluating the four kernel rules on the retained galaxy data reproduces the saved R variances within 3.56e-15. These are deterministic and contract checks, not new simulation evidence.

Both sample medians use midpoint conventions, MAD is unnormalised, sample variance has denominator n-1, and quantiles use R type 7. Kernel densities use direct Gaussian sums. The interval API keeps bandwidth and cutpoints on the raw measurement scale. Its ECDF uses closed comparisons at both endpoints; endpoint ties can therefore break exact reflection equivariance. Failure contracts cover zero MAD, nonfinite variance and nonrepresentable interval endpoints.

## Uncertainty and numerical limitations

The master seed is 20260917, with L'Ecuyer-CMRG/Inversion/Rejection and recorded experiment/population/size/replicate keys. Coverage MCSEs use the binomial formula. Within-study method differences use paired indicators on common datasets. TM comparisons use the common first 200 datasets. The retained R-only inference files contain paired principal and two-sample contrasts and simulation-median uncertainty. The latter used 400 diagnostic bootstrap resamples and seed 20261025; this differs from the study's 399-resample interval construction.

One pdQ attempt failed for Cauchy n=50, replicate 648; its hybrid inherits that failure. One TM attempt failed for Pareto n=1000, replicate 110. They reached the specified 700/1600 iteration limits and count as noncoverage. Their input samples and R diagnostics remain under `computation/validation/failed_fits/`.

The wider-bound sensitivity check refitted all 77 flagged records. The maximum absolute relative variance change was 4.5572%; 25 fits reached the wider bounds. One previously nonconverged Pareto fit converged, raising coverage from 0.945 to 0.950. Its original finite interval already covered, so this is failed-fit recovery rather than an endpoint crossing. The check changes bounds, starts and iteration budget; it cannot isolate the effect of bounds alone.

The specified GLD fits retain mixture variance distortion. Kernel consistency does not establish nominal finite-sample calibration: bimodal KDE-R coverage at n=50,000 is 0.969, while its oracle covers 0.957. The lognormal/exponential two-sample ratio coverages are 0.924 for KDE-R, 0.928 for KDE-half, 0.946 for pdQ and 0.945 for the oracle; the paired KDE-R-minus-oracle difference is -0.021 with MCSE 0.00535. For shoulder n=1000, bootstrap coverage 0.960 and KDE-half coverage 0.940 are equally far from 0.95; the bootstrap is longer (median relative lengths 0.203 versus 0.197). These qualifications remain in the manuscript.

The population searches and added-start calculations establish sensitivity, not global optimality. Descriptive closest-method rankings do not establish superiority. Timing batches describe their recorded concurrent workload, not controlled package benchmarks. Neither MADtests nor gld was executed, and no certification of package grid defaults is claimed. The optional native-package driver is clearly marked unexecuted.

## Artifacts and integrity

The modular and expanded manuscript sources are synchronized. The final PDF is compiled with TeX shell escape disabled and inspected. Table values and rankings are checked against the completed R table rows; numerical results and theoretical source are retained. The private Overleaf source export is compared against local publication files. Exact artifact hashes and project verification are recorded under `provenance/`.

`MANIFEST_SHA256.json` and `Rscript scripts/verify_workspace.R` verify the replication files. The R-only source, inputs, result records and build tools are included; local resumable checkpoints and transient logs are excluded from the portable package.
