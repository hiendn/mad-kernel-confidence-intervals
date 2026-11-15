# UNEXECUTED VALIDATION DRIVER. No values from this script appear in the manuscript.
# Run from the project root after installing MADtests and gld.
# Rscript computation/R/compare_native_packages.R [project-root]
args <- commandArgs(trailingOnly=TRUE)
root <- if(length(args)) normalizePath(args[1]) else normalizePath('.')
if(!requireNamespace('MADtests', quietly=TRUE)) stop('Install MADtests and its dependencies first.')
if(!requireNamespace('gld', quietly=TRUE)) stop('Install gld first.')
dir.create(file.path(root,'computation','validation','native_packages'),recursive=TRUE,showWarnings=FALSE)
files <- list.files(file.path(root,'computation','data','benchmark'),pattern='[.]csv$',full.names=TRUE)
if(!length(files)) stop('No fixture CSVs found.')
vm <- function(a,b,c,q) (1+((a-b)/c)^2+4*((a-b)/c)*q)/(4*(a+b)^2)
rows <- list(); k <- 0L
for(path in files) {
 x <- read.csv(path)$x; n <- length(x); m <- median(x); s <- mad(x,constant=1)
 h <- 1.482602218505602*s*n^(-1/5); cuts <- c(m-s,m+s,m)
 ds <- vapply(cuts,function(t) mean(dnorm((t-x)/h))/h,numeric(1))
 q <- 1-mean(x<=cuts[1])-mean(x<=cuts[2])
 funcs <- list(
   direct_gaussian_R=function() vm(ds[1],ds[2],ds[3],q),
   package_kernel_default=function() n*MADtests::varmad(x,constant=1,use.gld=FALSE),
   package_kernel_matched_h=function() n*MADtests::varmad(x,constant=1,use.gld=FALSE,bw=h,n=16384),
   package_GLD_TM=function() n*MADtests::varmad(x,constant=1,use.gld=TRUE,gld.est='TM'))
 for(nm in names(funcs)) {
   msg <- ''; tm <- system.time(val <- tryCatch(funcs[[nm]](),error=function(e){msg<<-conditionMessage(e);NA_real_}))
   k <- k+1L; rows[[k]] <- data.frame(fixture=basename(path),method=nm,n=n,Mhat=m,Shat=s,asymptotic_variance=as.numeric(val),elapsed=unname(tm['elapsed']),error=msg)
 }
}
write.csv(do.call(rbind,rows),file.path(root,'computation','validation','native_packages','comparison.csv'),row.names=FALSE)
capture.output(sessionInfo(),file=file.path(root,'computation','validation','native_packages','sessionInfo.txt'))
# Direct and package calculations use the identical input vectors.
# Package default differs in bandwidth/grid; matched_h still involves package interpolation.
# GLD-TM differences also include optimizer choices. Do not treat any discrepancy automatically as a bug.
