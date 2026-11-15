#!/usr/bin/env Rscript
# Four disjoint workers with one final summary writer. Re-running resumes each
# atomic replicate checkpoint and refuses a changed computation configuration.
args <- commandArgs(trailingOnly = TRUE)
file_arg <- grep("^--file=", commandArgs(), value = TRUE)
script <- normalizePath(sub("^--file=", "", file_arg[1]), mustWork = TRUE)
setwd(normalizePath(file.path(dirname(script), "../.."), mustWork = TRUE))
source("computation/R/nelder_mead.R")
source("computation/R/mad_methods.R")
source("computation/R/studies.R")
out <- if (length(args)) args[1] else "computation/generated_r/full/principal"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
writeLines(as.character(Sys.getpid()), file.path(out, "orchestrator.pid"))
files <- c("computation/R/nelder_mead.R", "computation/R/mad_methods.R", "computation/R/studies.R",
           "computation/R/launch_principal_full.R")
hashes <- tools::md5sum(files)
snapshot <- file.path(out, "source_snapshot")
dir.create(snapshot, showWarnings = FALSE)
for (file in files) {
  dest <- file.path(snapshot, basename(file))
  if (file.exists(dest) && tools::md5sum(dest) != tools::md5sum(file)) stop("Principal source changed; choose a new output directory")
  if (!file.exists(dest)) stopifnot(file.copy(file, dest))
}
write.csv(data.frame(path = files, md5 = unname(hashes)), file.path(out, "source_md5.csv"), row.names = FALSE)
groups <- list(0:2, 3:5, 6:7, 8:9)
jobs <- parallel::mclapply(seq_along(groups)-1L, function(shard) {
  child_out <- file.path(out, paste0("shard_", shard))
  dir.create(child_out, recursive = TRUE, showWarnings = FALSE)
  writeLines(as.character(Sys.getpid()), file.path(child_out, "worker.pid"))
  con <- file(file.path(child_out, "run.log"), "at")
  sink(con); sink(con, type = "message")
  on.exit({sink(type = "message"); sink(); close(con)})
  cat("Started worker", Sys.getpid(), "at", format(Sys.time(), tz = "UTC"), "UTC\n"); flush.console()
  run_principal(child_out, pop_indices = groups[[shard+1L]])
  writeLines(format(Sys.time(), tz = "UTC"), file.path(child_out, "COMPLETE"))
  TRUE
}, mc.cores = 4L, mc.preschedule = FALSE, mc.set.seed = FALSE)
if (any(vapply(jobs, inherits, logical(1), "try-error")) || !all(vapply(jobs, isTRUE, logical(1))))
  stop("At least one shard failed; see shard logs. Existing successful checkpoints remain available.")
paths <- unlist(lapply(seq_along(groups)-1L, function(shard)
  list.files(file.path(out, paste0("shard_", shard)), "^raw_[0-9]+_[0-9]+\\.csv\\.gz$", full.names = TRUE)))
stopifnot(length(paths) == 30L, !anyDuplicated(basename(paths)))
records <- mad_bind(lapply(paths, read.csv))
mad_summarize(records, out); mad_two_sample(records, populations(), out)
mad_csv(mad_bind(lapply(populations(), function(p) c(list(population=p$name, description=p$description),p$true))),file.path(out,"population_truth.csv"))
writeLines(format(Sys.time(), tz = "UTC"), file.path(out, "COMPLETE"))
cat("All 30 principal designs complete and merged at", format(Sys.time(), tz = "UTC"), "UTC\n")
