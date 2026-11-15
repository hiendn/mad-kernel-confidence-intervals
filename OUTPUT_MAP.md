# Numerical outputs for the 15 November 2026 manuscript

Saved display tables are under `computation/generated_r/full/finalize/assets/tables/`. The filenames retain the original study numbering. Match by subject using this map.

| Manuscript table | Subject | Saved display file |
|---|---|---|
|1|Complete additive coverage|`table_01_display.csv`|
|2|Larger-sample mixture comparison|`table_02_display.csv`|
|3|Log-interval coverage|`table_10_display.csv`|
|4|Selected median additive lengths|`table_11_display.csv`|
|5|Population pdQ solutions|`table_03_display.csv`|
|6|Shoulder added starts|`table_04_display.csv`|
|7|Selected shoulder fits|`table_05_display.csv`|
|8|Variance-estimator limit diagnostic|`table_06_display.csv`|
|9|Two-sample coverage|`table_07_display.csv`|
|10|Galaxy MAD intervals|`table_08_display.csv`|
|11|Galaxy local density and CDF estimates|`table_09_display.csv`|
|12|Population truth|`table_12_display.csv`|
|13|TM boundary sensitivity|`table_13_display.csv`|

The corresponding `table_XX_rows.tex` files contain table bodies. Figure inputs and source coordinates are in the adjacent `figures/`, `mixture_figure_data.csv` and `galaxies_density_grid.csv` files.

Additional results cited in the manuscript include:

- `computation/generated_r/full/inference/principal_median_mcse.csv` for simulation-median MCSEs.
- `computation/generated_r/full/finalize/assets/shoulder_selected_fit_distribution.csv` for selected-fit summaries at all four sample sizes.
- `computation/generated_r/full/finalize/assets/tables/hybrid_diagnostic.csv` for the hybrid pdQ-ECDF results.
- `computation/generated_r/full/diagnostic/finite_h_quadrature.csv` for the Gaussian finite-bandwidth diagnostic.
- `computation/generated_r/full/principal/` for replication-level lengths, timings, nuisance estimates, completion and coverage records.
- `computation/generated_r/full/tm/` and `computation/generated_r/full/audit/` for saved boundary and optimization checks.
