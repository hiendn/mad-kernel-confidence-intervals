# R runtime measurements

Executed on 25 October 2026 with R 4.4.0 on arm64 macOS. The principal study was configured for four workers; some completed while the study continued. The active-worker count during each batch was not measured. This benchmark used one additional sequential R process. These are descriptive timings under that workload, not a controlled native-package comparison.

The eight inputs are the archived Normal, Lognormal, Shoulder and Bimodal fixtures at sample sizes 200 and 1000. Each method was warmed on each input. The reported value for one dataset is the median of three elapsed-time batch averages. The ranges below are across the four dataset medians at each sample size.

| Method | n = 200 | n = 1000 |
| --- | ---: | ---: |
| KDE-R | 0.074–0.103 ms | 0.168–0.208 ms |
| GLD-pdQ | 2.85–3.28 ms | 2.95–3.75 ms |
| GLD-TM | 0.268–0.340 s | 0.689–1.055 s |
| 399-resample bootstrap | 16.85–19.60 ms | 37.55–44.05 ms |

KDE timing includes the direct local density calculation, algebraic variance and empirical influence-function diagnostic. GLD timings include fitting and evaluating the variance at the sample cutpoints. TM uses the previously computed pdQ candidate in its starting pool, so its time excludes that separate pdQ fit. Bootstrap timing includes new resample-index generation, sample MADs, percentile endpoints and the bootstrap variance.

The timer's smallest observed increment was 1 ms. All 96 retained measurement batches lasted at least 0.1 seconds. Each batch used at least 1000 KDE calls, 20 pdQ calls, three TM calls or 20 bootstrap calls, with larger batches where needed. Short preliminary batches remain in `batch_timings.csv` as calibration or `measurement_short` rows and are excluded from the summaries. Source and fixture hashes were checked before reusing a batch, and the source files were unchanged during the final invocation.

Run from the workspace root with:

```sh
Rscript computation/R/benchmark_reference.R .
```

An existing matching output directory reuses sufficiently long measured batches. Pass a second argument naming a new output directory to collect a fresh set of runtimes. `runtime_summary.csv` gives the per-dataset results, `runtime_ranges.csv` gives the displayed ranges, and `environment.json` plus `sessionInfo.txt` record the execution environment and protocol.
