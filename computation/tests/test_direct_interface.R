# R interval contracts and actual command-line execution; no simulations or GLD fits.
source('computation/R/mad_interval.R')
checks <- list();add <- function(name,ok) {
 checks[[length(checks)+1L]] <<- list(case=name,passed=isTRUE(ok));if(!isTRUE(ok))stop(name)
}
close <- function(a,b)length(a)==length(b)&&all(abs(a-b)<=2e-12*(1+abs(b)))
samples <- list(even=c(-4,-2,-1,1,2,4),odd=c(-5,-1,0,2,3,6,9),grid=seq(-3,3,length.out=100),
 galaxies=read.csv('computation/data/galaxies.csv')$velocity_km_s/1000)
for(name in names(samples))for(rule in c('KDE-R','KDE-half','KDE-double','KDE-NR')) {
 x<-samples[[name]];v<-mad_interval(x,bandwidth=rule);s<-mad_interval(x,bandwidth=rule,report_scale=1.4826)
 a<-mad_interval(3*x+7,bandwidth=rule)
 add(paste(name,rule,'finite'),v$success&&all(is.finite(c(v$additive_interval,v$log_interval)))&&all(v$log_interval>0))
 add(paste(name,rule,'report_scale'),identical(v$h,s$h)&&identical(v$cutpoints,s$cutpoints)&&
 close(s$additive_interval,1.4826*v$additive_interval)&&close(s$log_interval,1.4826*v$log_interval))
 add(paste(name,rule,'positive_affine'),close(a$additive_interval,3*v$additive_interval)&&close(a$log_interval,3*v$log_interval))
 if(!any(x %in% v$cutpoints[1:2])) {
  a<-mad_interval(-3*x+7,bandwidth=rule)
  add(paste(name,rule,'reflection_without_endpoint_ties'),close(a$additive_interval,3*v$additive_interval))
 }
}
a<-mad_interval(samples$odd);b<-mad_interval(-3*samples$odd+7)
add('closed_ecdf_endpoint_tie',close(a$local_estimates$q,0)&&close(b$local_estimates$q,-1/7)&&
 !close(b$root_n_variance_raw,9*a$root_n_variance_raw))
failures<-list(mad_interval(rep(2,4)),mad_interval(samples$even,confidence=1-.Machine$double.eps/2),
 mad_interval(samples$even,report_scale=1e308),mad_interval(samples$even,report_scale=1e160),
 mad_interval(samples$even*.1,report_scale=.Machine$double.xmin*.Machine$double.eps))
for(i in seq_along(failures)) {
 v<-failures[[i]];add(paste('finite_failure_contract',i),!v$success&&is.null(v$additive_interval)&&
 is.null(v$log_interval)&&!grepl('NaN|Infinity|Inf',mad_interval_json(v)))
}
for(x in list(0,c(0,Inf),matrix(1:4,2),c(0,NA),numeric()))
 add(paste('invalid_input',length(checks)),inherits(try(mad_interval(x),silent=TRUE),'try-error'))
csv<-tempfile(fileext='.csv');write.csv(data.frame(value=samples$even),csv,row.names=FALSE)
for(name in c('ordinary','confidence_boundary','reporting_overflow')) {
 extra<-switch(name,ordinary=character(),confidence_boundary=c('--confidence',format(1-.Machine$double.eps/2,digits=17)),reporting_overflow=c('--report-scale','1e308'))
 output<-system2(file.path(R.home('bin'),'Rscript'),c('computation/R/mad_interval.R','--data',shQuote(csv),'--column','value',extra),stdout=TRUE)
 parsed<-jsonlite::fromJSON(paste(output,collapse='\n'));add(paste('CLI',name),identical(parsed$success,name=='ordinary'))
}
jsonlite::write_json(list(scope='R contracts and CLI; no study rerun',R=R.version.string,checks=checks),
 'computation/tests/direct_interface_checks.json',pretty=TRUE,auto_unbox=TRUE,digits=17)
cat(length(checks),'R interval checks passed.\n')
