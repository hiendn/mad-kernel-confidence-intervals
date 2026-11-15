#!/usr/bin/env Rscript
source('computation/R/nelder_mead.R')
source('computation/R/mad_methods.R')
source('computation/R/studies.R')
base <- 'computation/generated_r/full'
paths <- list.files(file.path(base, 'principal'), '^raw_[0-9]+_[0-9]+[.]csv[.]gz$', recursive = TRUE, full.names = TRUE)
stopifnot(length(paths) == 30L, !anyDuplicated(basename(paths)))
frames <- lapply(paths, read.csv)
stopifnot(sum(vapply(frames, nrow, integer(1))) == 246000L)
for (d in frames) {
  count <- table(d$method)
  stopifnot(count[['GLD-TM']] == 200L, all(count[names(count) != 'GLD-TM'] == 1000L))
}
out <- tempfile('mad-saved-summary-');dir.create(out)
actual <- mad_summarize(mad_bind(frames), out)
expected <- read.csv(file.path(base, 'principal', 'one_sample_summary.csv'))
keys <- c('population', 'n', 'method')
order_rows <- function(d) d[do.call(order, unname(d[keys])), , drop = FALSE]
a <- order_rows(actual); b <- order_rows(expected)
rownames(a) <- rownames(b) <- NULL
stopifnot(nrow(a) == 270L, identical(a[keys], b[keys]))
numeric_columns <- names(b)[vapply(b, is.numeric, logical(1))]
errors <- vapply(numeric_columns, function(n) max(abs(a[[n]] - b[[n]]), na.rm = TRUE), numeric(1))
# CSV decimal round trips can perturb elapsed milliseconds at about1e-12.
tolerance <- ifelse(grepl('_ms$',names(errors)),1e-9,5e-12)
stopifnot(all(errors < tolerance))
tm <- read.csv(file.path(base, 'tm', 'tm_boundary_refit_raw.csv'))
stopifnot(nrow(tm) == 77L)
restart <- list.files(file.path(base, 'restarts'), '^shoulder_restart_raw_[0-9]+[.]csv[.]gz$', recursive = TRUE, full.names = TRUE)
stopifnot(length(restart) == 4L, all(vapply(restart, function(p) nrow(read.csv(p)), integer(1)) == 1000L))
report <- list(scope='Reaggregation of saved records; no sampling, bootstrap resampling, GLD fits or optimizer runs',
 principal_designs=length(paths), principal_procedure_rows=246000L, method_summaries=nrow(a),
 maximum_absolute_statistical_summary_error=max(errors[!grepl('_ms$',names(errors))]),
 maximum_absolute_timing_error_ms=max(errors[grepl('_ms$',names(errors))]), TM_boundary_refits=nrow(tm), shoulder_restart_records=4000L)
print(report)
args <- commandArgs(trailingOnly=TRUE)
if(length(args)) jsonlite::write_json(report,args[1L],pretty=TRUE,auto_unbox=TRUE,digits=17)
invisible(report)
