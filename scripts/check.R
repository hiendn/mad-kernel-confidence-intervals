#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)
script <- sub('^--file=', '', grep('^--file=', commandArgs(), value=TRUE)[1L])
script <- gsub('~+~', ' ', script, fixed=TRUE)
root <- normalizePath(file.path(dirname(script), '..'), mustWork=TRUE)
out <- if(length(args)) normalizePath(args[1L], mustWork=FALSE) else tempfile('mad-verification-')
if(file.exists(out)) stop('Choose a new verification output directory')
dir.create(out, recursive=TRUE)
out <- normalizePath(out, mustWork=TRUE)
work <- tempfile('mad-check-copy-');dir.create(work)
paths <- list.files(root, recursive=TRUE, all.files=FALSE, full.names=FALSE)
paths <- paths[!grepl('^(_reproduced|_verification)/', paths)]
for(p in paths) {
 dir.create(file.path(work,dirname(p)),recursive=TRUE,showWarnings=FALSE)
 if(!file.copy(file.path(root,p),file.path(work,p)))stop('Copy failed: ',p)
}
setwd(work)
for(test in c('computation/R/test_core.R','computation/tests/test_direct_interface.R',
              'scripts/verify_saved_results.R')) {
 log <- file.path(out,paste0(basename(test),'.log'))
 test_args <- test
 if(basename(test)=='verify_saved_results.R')
  test_args <- c(test_args,shQuote(file.path(out,'saved_results_checks.json')))
 status <- system2(file.path(R.home('bin'),'Rscript'),test_args,stdout=log,stderr=log)
 cat(readLines(log,warn=FALSE),sep='\n')
 if(status!=0L)stop('Check failed: ',test,'; see ',log)
}
for(p in c('computation/validation/core_checks.json','computation/validation/core_sessionInfo.txt',
 'computation/tests/direct_interface_checks.json')) file.copy(p,file.path(out,basename(p)))
cat('All checks passed. Reports: ',out,'\n',sep='')
