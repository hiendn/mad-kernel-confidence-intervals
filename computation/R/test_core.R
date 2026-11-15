#!/usr/bin/env Rscript
# Deterministic identities, numerical accuracy and failure contracts.
source('computation/R/nelder_mead.R')
source('computation/R/mad_methods.R')
source('computation/R/studies.R')
out <- 'computation/validation'; dir.create(out,recursive=TRUE,showWarnings=FALSE)
pops <- populations(); roundtrip <- gradient_error <- root_error <- 0
for(ab in list(c(.2,.5),c(0,0),c(-.8,-.8),c(.8,-.5),c(4,.1))) {
 par <- c(1,2,ab); p <- c(seq(.001,.999,length.out=1000),1e-6,1-1e-6)
 value <- cdf_pdf(qvec(p,par),par); stopifnot(all(value$pdf>0))
 roundtrip <- max(roundtrip,max(abs(value$cdf-p)))
}
stopifnot(roundtrip<1e-8,identical(qvec(c(0,1),c(1,2,0,0)),c(-Inf,Inf)))
for(pop in pops) {
 u <- unlist(pop$true[c('a','b','c','q')]); g <- do.call(variance_gradient,as.list(u))
 numerical <- vapply(1:4,function(j) {
  h <- 1e-5*max(abs(u[j]),.1); hi <- lo <- u;hi[j]<-hi[j]+h;lo[j]<-lo[j]-h
  (do.call(variance_map,as.list(hi))-do.call(variance_map,as.list(lo)))/(2*h)
 },numeric(1))
 gradient_error <- max(gradient_error,max(abs(numerical-g)/(1+abs(g))))
 root_error <- max(root_error,abs(pop$cdf(pop$true$M)-.5),
  abs(pop$cdf(pop$true$M+pop$true$S)-pop$cdf(pop$true$M-pop$true$S)-.5))
}
stopifnot(gradient_error<1e-6,root_error<1e-9)
x <- c(-5,-2,-1,0,1,4,9,12);base <- kernel_estimate(x)$value
stopifnot(abs(kernel_estimate(7+3*x)$value/(9*base)-1)<1e-10)
stopifnot(!kernel_estimate(rep(1,20))$info$success,!fit_tm(rep(1,20))$info$success,
 !fit_pdq(rep(1,20))$info$success,is.nan(variance_map(1,0,1e-320,0)),
 is.nan(variance_map(1e-320,1e-320,1,0)))
tied <- mad_nuisance(c(0,0,1,1));stopifnot(tied$m==.5,tied$s==.5,
 kernel_estimate(c(0,0,1,1))$info$q==-.5)
# Analytic optimizer targets and an imposed evaluation limit.
f <- nelder_mead(function(z) sum((z-c(2,-3))^2),c(0,0),maxiter=1000,xatol=1e-9,fatol=1e-12)
stopifnot(f$success,max(abs(f$x-c(2,-3)))<1e-7)
# Check the archived R galaxy calculation against its retained input, without refitting.
x <- read.csv('computation/data/galaxies.csv')$velocity_km_s/1000
saved <- read.csv('computation/generated_r/full/galaxy/galaxies_analysis.csv')
errors <- vapply(c('KDE-R','KDE-NR','KDE-half','KDE-double'),function(method) {
 bw <- if(method=='KDE-NR')'nrd0' else 'robust';mult <- switch(method,'KDE-half'=.5,'KDE-double'=2,1)
 abs(kernel_estimate(x,bw,mult)$value-saved$vhat[saved$method==method])
},numeric(1)); stopifnot(max(errors)<1e-10)
report <- list(scope='Deterministic identities and saved R galaxy-input verification; no simulations rerun',
 R=R.version.string,FKML_CDF_quantile_max_error=roundtrip,
 variance_gradient_scaled_max_error=gradient_error,population_equation_max_error=root_error,
 galaxy_max_absolute_variance_errors=as.list(errors),optimizer_quadratic_max_parameter_error=max(abs(f$x-c(2,-3))),
 edge_tests='zero MAD, nonfinite intermediate variance, midpoint medians and closed endpoint ECDF passed')
jsonlite::write_json(report,file.path(out,'core_checks.json'),pretty=TRUE,auto_unbox=TRUE,digits=17)
capture.output(sessionInfo(),file=file.path(out,'core_sessionInfo.txt')); print(report)
