#!/usr/bin/env Rscript
# Run from any directory. Outputs default to this durable workspace.
args <- commandArgs(trailingOnly = TRUE)
file_arg <- grep("^--file=", commandArgs(), value = TRUE)
if (!length(file_arg)) stop("Run this driver with Rscript")
script <- normalizePath(sub("^--file=", "", file_arg[1]), mustWork = TRUE)
root <- normalizePath(file.path(dirname(script), "../.."), mustWork = TRUE)
setwd(root)
if (!length(args) || args[1] %in% c("--help", "-h")) {
  cat("Rscript computation/R/run.R TASK [--mode smoke|full] [--out PATH]\n",
      "Tasks: principal, population, galaxy, clt, large, resolution, restarts, tm, audit, summary-checks, finalize, diagnostic, all\n",
      "--mode defaults to smoke. Full uses the manuscript replication counts.\n",
      "--out defaults to computation/generated_r/{MODE}/{TASK}. Checkpoints resume automatically.\n",
      "finalize requires completed R studies.\n", sep = "")
  quit(status = 0L)
}
task <- args[1]; rest <- args[-1]
if (length(rest) %% 2L) stop("Options must have values")
options <- if (length(rest)) setNames(as.list(rest[seq(2, length(rest), 2)]), rest[seq(1, length(rest), 2)]) else list()
if (length(setdiff(names(options), c("--mode", "--out")))) stop("Unknown option")
mode <- if (is.null(options[["--mode"]])) "smoke" else options[["--mode"]]
if (!mode %in% c("smoke", "full")) stop("--mode must be smoke or full")
out <- if (is.null(options[["--out"]])) file.path("computation/generated_r", mode, task) else options[["--out"]]
source("computation/R/nelder_mead.R")
source("computation/R/mad_methods.R")
source("computation/R/studies.R")
full <- mode == "full"
run_task <- function(which, destination) {
  dir.create(destination, recursive = TRUE, showWarnings = FALSE)
  if (which == "principal" && file.exists(file.path(destination,"orchestrator.pid")))
    stop("This directory belongs to the sharded driver. Resume launch_principal_full.R instead.")
  cat("Started", which, "mode", mode, "at", format(Sys.time(), tz = "UTC"), "UTC\n"); flush.console()
  switch(which,
    principal = run_principal(destination, replications = if (full) 1000L else 3L,
      tm_replications = if (full) 200L else 1L, bootstrap = if (full) 399L else 19L,
      pop_indices = 0:9, ns = if (full) c(50L,200L,1000L) else 50L),
    population = run_population(destination, pop_indices = if (full) 0:9 else c(0L,7L,8L),
      Js = if (full) c(25L,50L) else 25L),
    galaxy = run_galaxy(destination),
    clt = run_variance_clt(destination, replications = if (full) 3000L else 20L,
      ns = if (full) c(200L,1000L,5000L) else 200L),
    large = run_long_samples(destination, replications = if (full) 1000L else 1L,
      pop_indices = if (full) NULL else c(0L,7L,8L)),
    resolution = run_long_samples(destination, resolution = TRUE,
      replications = if (full) 1000L else 1L),
    restarts = {
      source("computation/R/optimizer_checks.R")
      run_pdq_restarts(destination, replications = if (full) 1000L else 2L,
        ns = if (full) c(200L,1000L,5000L,50000L) else c(200L,1000L),
        results_dir = dirname(destination))
    },
    `summary-checks` = { source("computation/R/optimizer_checks.R"); run_saved_checks(destination, results_dir = dirname(destination)) },
    tm = {
      source("computation/R/optimizer_checks.R")
      run_tm_boundary(destination, results_dir = file.path(dirname(destination), "principal"),
        ns = if (full) c(50L,200L,1000L) else 50L, max_refits = if (full) Inf else 2L)
    },
    audit = {
      source("computation/R/optimizer_checks.R")
      run_global_search(file.path(destination, "global"), maxiter = if (full) 800L else 60L,
        popsize = if (full) 30L else 6L)
      run_multistart_audit(file.path(destination, "multistart"), reps = if (full) 10L else 1L,
        ns = if (full) c(200L,1000L) else 200L, starts = if (full) 10L else 2L)
    },
    finalize = {
      source("computation/R/optimizer_checks.R")
      source("computation/R/assets.R")
      finalize_r_assets(file.path(dirname(destination)), destination)
    },
    diagnostic = { source("computation/R/normal_bandwidth_diagnostic.R"); normal_bandwidth_diagnostic(destination) },
    stop("Unknown task: ", which))
  cat("Completed", which, "at", format(Sys.time(), tz = "UTC"), "UTC\n"); flush.console()
}
if (task == "all") {
  for (which in c("principal", "population", "galaxy", "clt", "large", "resolution", "restarts", "tm", "audit", "summary-checks", "diagnostic"))
    run_task(which, file.path(out, which))
  if (full) run_task("finalize", file.path(out,"finalize"))
} else run_task(task, out)
