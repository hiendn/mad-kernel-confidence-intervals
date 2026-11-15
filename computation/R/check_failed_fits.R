# Diagnose the two optimizer failures in the complete R principal run.
# Save the exact sample values and convergence diagnostics.
source("computation/R/nelder_mead.R")
source("computation/R/mad_methods.R")
source("computation/R/studies.R")
out <- "computation/validation/failed_fits"
dir.create(out, recursive=TRUE, showWarnings=FALSE)
specs <- data.frame(population=c("Cauchy","Pareto7"), index=c(5L,4L),
                    n=c(50L,1000L), rep=c(648L,110L), method=c("GLD-pdQ","GLD-TM"))
rows <- list()
for(i in seq_len(nrow(specs))) {
  s <- specs[i,]; x <- mad_sample(populations()[[s$index+1L]],1L,s$index,s$n,s$rep)
  id <- paste(s$population,s$n,s$rep,sep="_")
  writeLines(c("x",sprintf("%.17g",x)),file.path(out,paste0(id,".csv")))
  fit <- fit_pdq(x)
  if(s$method=="GLD-TM")fit<-fit_tm(x,fit$par)
  value <- gld_estimate(x,fit$par)$value
  rows[[i]]<-data.frame(id=id,method=s$method,success=fit$info$success,nit=fit$info$nit,
    objective=fit$info$objective,boundary=fit$info$boundary,variance=value,
    loc=fit$par[1],invscale=fit$par[2],shape3=fit$par[3],shape4=fit$par[4])
}
write.csv(do.call(rbind,rows),file.path(out,"R_fits.csv"),row.names=FALSE)
print(do.call(rbind,rows))
