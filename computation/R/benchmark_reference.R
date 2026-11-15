#!/usr/bin/env Rscript
# Descriptive R runtimes on eight frozen inputs. Timings describe the recorded concurrent workload.
# Run from any directory: Rscript computation/R/benchmark_reference.R [workspace] [output]
args <- commandArgs(trailingOnly = TRUE)
file_arg <- grep("^--file=", commandArgs(), value = TRUE)
script <- normalizePath(sub("^--file=", "", file_arg[1L]), mustWork = TRUE)
root <- if (length(args)) normalizePath(args[1L], mustWork = TRUE) else
  normalizePath(file.path(dirname(script), "../.."), mustWork = TRUE)
out <- if (length(args) >= 2L) args[2L] else file.path(root, "computation/validation/benchmark_reference")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
source(file.path(root, "computation/R/nelder_mead.R"))
source(file.path(root, "computation/R/mad_methods.R"))
source(file.path(root, "computation/R/studies.R"))
stopifnot(requireNamespace("jsonlite", quietly = TRUE))

population_names <- c("Normal", "Lognormal", "Shoulder", "Bimodal")
sizes <- c(200L, 1000L)
minimum_repeats <- c("KDE-R" = 1000L, "GLD-pdQ" = 20L,
                     "GLD-TM" = 3L, "Bootstrap399" = 20L)
minimum_batch_seconds <- 0.1
measured_batches <- 3L
started <- format(Sys.time(), tz = "UTC", usetz = TRUE)
previous <- if (file.exists(file.path(out, "environment.json")))
  jsonlite::fromJSON(file.path(out, "environment.json"), simplifyVector = FALSE) else NULL
source_files <- file.path(root, c("computation/R/nelder_mead.R", "computation/R/mad_methods.R",
                                "computation/R/studies.R", "computation/R/benchmark_reference.R"))
source_hashes <- tools::md5sum(source_files)
fixture_files <- unlist(lapply(population_names, function(pop)
  file.path(root, "computation/data/benchmark", paste0(pop, "_", sizes, ".csv"))))
fixture_hashes <- tools::md5sum(fixture_files)
relative_source_hashes <- as.list(setNames(unname(source_hashes), sub(paste0(root, "/"), "", names(source_hashes), fixed = TRUE)))
relative_fixture_hashes <- as.list(setNames(unname(fixture_hashes), basename(names(fixture_hashes))))
core_names <- c("computation/R/nelder_mead.R", "computation/R/mad_methods.R", "computation/R/studies.R")
if (!is.null(previous)) {
  stopifnot(identical(previous$R, R.version.string),
    identical(unlist(previous$source_md5[core_names]), unlist(relative_source_hashes[core_names])),
    identical(unlist(previous$fixture_md5), unlist(relative_fixture_hashes)))
}
ticks <- numeric(10000L)
for (i in seq_along(ticks)) ticks[i] <- unname(proc.time()[["elapsed"]])
differences <- diff(ticks)
observed_tick <- if (any(differences > 0)) min(differences[differences > 0]) else NA_real_
environment <- list(started_UTC = if (is.null(previous)) started else previous$started_UTC,
  latest_invocation_started_UTC = started, R = R.version.string,
  platform = R.version$platform, machine = unname(Sys.info()[["machine"]]),
  operating_system = paste(Sys.info()[c("sysname", "release")], collapse = " "),
  BLAS = unname(extSoftVersion()[["BLAS"]]),
  RNG = "L'Ecuyer-CMRG / Inversion / Rejection", seed = MAD_SEED,
  timer = "proc.time elapsed/user.self/sys.self; system.time wraps an explicit repeated-call loop",
  smallest_observed_elapsed_increment_seconds = observed_tick,
  concurrent_workload = paste("Coordinator reported four principal R simulation workers active.",
    "This benchmark adds one sequential R process; no controlled idle-machine comparison was made."),
  warmup = "One call per method and dataset before timing; cached weights/grid and R JIT are warm.",
  batch_design = list(minimum_repeats = as.list(minimum_repeats),
    minimum_elapsed_seconds = minimum_batch_seconds, measured_batches = measured_batches),
  method_scope = list(
    "KDE-R" = "kernel_estimate(x), including the local densities, algebraic variance and empirical influence-function diagnostic",
    "GLD-pdQ" = "fit_pdq(x) plus gld_estimate(x, fitted parameters)",
    "GLD-TM" = "fit_tm(x, precomputed pdQ parameters) plus gld_estimate; excludes the precomputed pdQ starting fit",
    "Bootstrap399" = "399 new n-out-of-n resample-index rows, unnormalised midpoint MADs, percentile endpoints and resample variance"),
  resume_policy = "Reuse only measured batches at least 0.1 seconds with the same core-source and fixture hashes; extend shorter batches.",
  previous_benchmark_script_md5 = if (is.null(previous)) NULL else previous$source_md5[["computation/R/benchmark_reference.R"]],
  source_md5 = relative_source_hashes, fixture_md5 = relative_fixture_hashes)
jsonlite::write_json(environment, file.path(out, "environment.json"), pretty = TRUE, auto_unbox = TRUE, digits = 17)
capture.output(sessionInfo(), file = file.path(out, "sessionInfo.txt"))

rows <- if (file.exists(file.path(out, "batch_timings.csv"))) {
  prior <- read.csv(file.path(out, "batch_timings.csv"))
  prior$phase[prior$phase == "measurement" & prior$elapsed_seconds < minimum_batch_seconds] <- "measurement_short"
  lapply(seq_len(nrow(prior)), function(i) prior[i, ])
} else list()
append_timing <- function(row) {
  rows[[length(rows) + 1L]] <<- row
  write.csv(do.call(rbind, rows), file.path(out, "batch_timings.csv"), row.names = FALSE)
}
measure <- function(fn, repeats, population, n, method, phase, batch, seed_index) {
  mad_rng(80L, match(population, population_names) - 1L, n, seed_index)
  gc(verbose = FALSE)
  timing <- system.time(for (i in seq_len(repeats)) value <- fn(), gcFirst = FALSE)
  # Check result after timing; reference fixture fits are required to be usable.
  if (!isTRUE(value$success) || !is.finite(value$value)) stop("Benchmark calculation failed: ", population, "/", n, "/", method)
  row <- data.frame(population = population, n = n, method = method,
    phase = phase, batch = batch, repeats = repeats,
    elapsed_seconds = unname(timing[["elapsed"]]),
    user_seconds = unname(timing[["user.self"]]), system_seconds = unname(timing[["sys.self"]]),
    elapsed_seconds_per_call = unname(timing[["elapsed"]]) / repeats,
    CPU_seconds_per_call = unname(timing[["user.self"]] + timing[["sys.self"]]) / repeats)
  append_timing(row)
  row$elapsed_seconds
}

for (population in population_names) for (n in sizes) {
  path <- file.path(root, "computation/data/benchmark", paste0(population, "_", n, ".csv"))
  x <- read.csv(path)$x
  pdq_start <- fit_pdq(x)$par
  functions <- list(
    "KDE-R" = function() { z <- kernel_estimate(x); list(value = z$value, success = z$info$success) },
    "GLD-pdQ" = function() {
      fit <- fit_pdq(x); z <- gld_estimate(x, fit$par)
      list(value = z$value, success = fit$info$success && z$info$success)
    },
    "GLD-TM" = function() {
      fit <- fit_tm(x, pdq_start); z <- gld_estimate(x, fit$par)
      list(value = z$value, success = fit$info$success && z$info$success)
    },
    "Bootstrap399" = function() {
      indices <- matrix(sample.int(n, 399L * n, replace = TRUE), nrow = 399L, byrow = TRUE)
      z <- bootstrap_from_indices(x, indices)
      list(value = z$root_n_variance, success = all(is.finite(c(z$endpoints, z$root_n_variance))))
    })
  for (method in names(functions)) {
    fn <- functions[[method]]
    prior <- if (length(rows)) do.call(rbind, rows) else NULL
    completed <- if (is.null(prior)) integer() else prior$batch[
      prior$population == population & prior$n == n & prior$method == method & prior$phase == "measurement"]
    if (all(seq_len(measured_batches) %in% completed)) next
    invisible(fn())
    repeats <- minimum_repeats[[method]]
    attempt <- 0L
    repeat {
      attempt <- attempt + 1L
      elapsed <- measure(fn, repeats, population, n, method, "calibration", attempt, attempt)
      if (elapsed >= minimum_batch_seconds) break
      repeats <- repeats * 2L
      if (attempt >= 10L) stop("Could not resolve timer for ", method)
    }
    for (batch in setdiff(seq_len(measured_batches), completed)) {
      repeat {
        elapsed <- measure(fn, repeats, population, n, method, "measurement", batch, 100L + batch)
        if (elapsed >= minimum_batch_seconds) break
        rows[[length(rows)]]$phase <- "measurement_short"
        repeats <- repeats * 2L
      }
    }
    cat(population, n, method, repeats, "calls/batch\n")
    flush.console()
  }
}
data <- do.call(rbind, rows)
measured <- data[data$phase == "measurement", ]
summary <- do.call(rbind, lapply(split(measured,
  interaction(measured$population, measured$n, measured$method, drop = TRUE)), function(g) {
    data.frame(population = g$population[1L], n = g$n[1L], method = g$method[1L],
      min_repeats_per_batch = min(g$repeats), max_repeats_per_batch = max(g$repeats), batches = nrow(g),
      median_elapsed_ms_per_call = 1000 * median(g$elapsed_seconds_per_call),
      min_elapsed_ms_per_call = 1000 * min(g$elapsed_seconds_per_call),
      max_elapsed_ms_per_call = 1000 * max(g$elapsed_seconds_per_call),
      median_CPU_ms_per_call = 1000 * median(g$CPU_seconds_per_call))
  }))
rownames(summary) <- NULL
write.csv(summary, file.path(out, "runtime_summary.csv"), row.names = FALSE)
ranges <- do.call(rbind, lapply(split(summary, interaction(summary$n, summary$method, drop = TRUE)), function(g)
  data.frame(n = g$n[1L], method = g$method[1L], datasets = nrow(g),
    minimum_dataset_median_ms = min(g$median_elapsed_ms_per_call),
    maximum_dataset_median_ms = max(g$median_elapsed_ms_per_call))))
rownames(ranges) <- NULL
write.csv(ranges, file.path(out, "runtime_ranges.csv"), row.names = FALSE)
environment$completed_UTC <- format(Sys.time(), tz = "UTC", usetz = TRUE)
environment$measurement_batches <- nrow(measured)
environment$all_measurement_batches_at_least_0_1_seconds <- all(measured$elapsed_seconds >= minimum_batch_seconds)
environment$source_hashes_unchanged <- identical(source_hashes, tools::md5sum(source_files))
jsonlite::write_json(environment, file.path(out, "environment.json"), pretty = TRUE, auto_unbox = TRUE, digits = 17)
stopifnot(environment$source_hashes_unchanged)
print(ranges)
