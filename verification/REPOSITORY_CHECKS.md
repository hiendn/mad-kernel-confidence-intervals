# Repository verification, 15 November 2026

The original 239-entry native R release was hash-checked before staging. The publication copy retains its statistical methods, saved observations, random streams and numerical results. The original local archive remains unchanged.

Executed checks with R 4.4.0 and jsonlite 1.8.8:

- Deterministic FKML inversion, population-root, variance-gradient, scaling, exceptional-input and saved galaxy-input checks passed.
- All 70 interval and command-line contract checks passed.
- All 30 principal files contain the declared method counts, totalling 246000 procedure rows and 270 method summaries. Reaggregated statistical summaries differ from the saved CSV values by at most 5.33e-15. The largest timing round-trip difference is 5.01e-12 milliseconds.
- The 77 TM boundary refits and 4000 shoulder added-start records were verified.
- Rebuilding from the archived records reproduced all 13 display tables and all 48 mixture-plot points exactly. The largest galaxy-density-grid difference was 3.03e-15.
- All three rebuilt figures match the original saved figures pixel-for-pixel under Poppler at 100 dpi.

The check command was exercised using a new relative output directory and the reproduction command using an absolute script path containing spaces. Tests run in isolated R processes and a temporary copy; their saved reports are included here. No simulation, bootstrap resampling, GLD fit or study optimization was rerun. The deterministic optimizer unit test uses an analytic quadratic target.

Publication adaptations are confined to optional historical-manuscript checking, the saved galaxy-parameter text fallback, two configuration path fields and repository navigation/check scripts. Saved numeric CSVs, galaxy observations, benchmark vectors and original figure PDFs are unchanged. Historical manuscript-specific files, stale progress notes and local path inventories are excluded. Table filenames retain their original numbering; OUTPUT_MAP.md maps them to the 15 November manuscript.

The local compiled manuscript and annotated reviews are maintained separately. This code-and-results repository resolves the missing public locator in AM15-CHK-REPO; it does not apply the new optional readability review.
