# Optimizer sensitivity studies. Source mad_methods.R, nelder_mead.R and studies.R first.
# These diagnostics leave the frozen principal estimator unchanged.
optimizer_source_md5 <- function()unname(tools::md5sum("computation/R/optimizer_checks.R"))
mad_shape_data <- function(x) {
  xs <- sort(x); pw <- empirical_qd_weights(length(x),if(length(x)<=200)25L else 50L)
  qd <- as.vector(pw$W %*% xs)
  if(any(!is.finite(qd)) || any(qd<=0)) stop("Invalid empirical quantile densities")
  raw <- 1/qd
  list(xs=xs,ps=pw$ps,y=raw/mean(raw))
}

# DE/best/1/bin with immediate updates, Latin-hypercube initial points,
# generation-wise mutation dithering U(.5,1), recombination .7 and bound repair
# by independent uniform draws, followed by Nelder-Mead polishing.
# No global-optimality certificate is implied.
mad_differential_evolution <- function(fn,lower,upper,popsize=30L,maxiter=800L,tol=1e-10) {
  d <- length(lower); N <- popsize*d
  stopifnot(d>=1L,length(upper)==d,N>=4L,all(upper>lower))
  u <- vapply(seq_len(d),function(j) (sample.int(N,N)-1+runif(N))/N,numeric(N))
  pool <- sweep(sweep(u,2L,upper-lower,"*"),2L,lower,"+")
  evaluations <- 0L
  evaluate <- function(x) { evaluations <<- evaluations+1L; fn(x) }
  values <- apply(pool,1L,evaluate); success <- FALSE
  for(iter in seq_len(maxiter)) {
    mutation <- runif(1L,.5,1)
    for(i in seq_len(N)) {
      best <- which.min(values); candidates <- sample(setdiff(seq_len(N),i),2L)
      proposal <- pool[best,]+mutation*(pool[candidates[1],]-pool[candidates[2],])
      bad <- proposal<lower | proposal>upper | !is.finite(proposal)
      if(any(bad)) proposal[bad] <- runif(sum(bad),lower[bad],upper[bad])
      cross <- runif(d)<.7; cross[sample.int(d,1)] <- TRUE
      trial <- pool[i,]; trial[cross] <- proposal[cross]
      v <- evaluate(trial)
      if(is.finite(v) && v<=values[i]) { pool[i,] <- trial; values[i] <- v }
    }
    if(sd(values)*sqrt((N-1)/N)<=tol*abs(mean(values))) { success <- TRUE; break }
  }
  best <- which.min(values)
  polish <- nelder_mead(function(x) if(any(x<lower|x>upper))1e12 else evaluate(x),
                        pool[best,],maxiter=2000L,xatol=1e-10,fatol=1e-12)
  if(polish$fun<values[best]) { pool[best,] <- polish$x; values[best] <- polish$fun }
  list(x=as.numeric(pool[best,]),fun=values[best],success=success,nit=iter,nfev=evaluations)
}

run_global_search <- function(out,pop_indices=c(6L,7L,8L),maxiter=800L,popsize=30L,
                              seed_mode=c("source_scalar","keyed")) {
  seed_mode <- match.arg(seed_mode)
  cfg <- mad_config(out,"global_search",list(implementation_md5=optimizer_source_md5(),pop_indices=pop_indices,maxiter=maxiter,popsize=popsize,
    algorithm="base R DE/best/1/bin; Latin hypercube; immediate updates; NM polish",
    seed_mode=seed_mode,scalar_seed_formula="783 + zero-based population index"))
  pops <- populations(); rows <- list()
  for(pi in pop_indices) {
    row <- mad_checkpoint(file.path(out,"checkpoints",sprintf("de_%d.rds",pi)),cfg,function() {
      pop <- pops[[pi+1L]]; ps <- (seq_len(50)-.5)/50
      y <- pop$pdf(pop$ppf(ps)); y <- y/mean(y)
      if(seed_mode=="source_scalar") {
        RNGkind("L'Ecuyer-CMRG",normal.kind="Inversion",sample.kind="Rejection")
        set.seed(783L+pi)
      } else mad_rng(99L,pi,0L,0L)
      fit <- mad_differential_evolution(function(sh)pdq_obj(sh,y,ps),c(-2,-2),c(5,5),
                                        popsize=popsize,maxiter=maxiter)
      qq <- pop$ppf(c(.25,.5,.75)); gq <- qvec(c(.25,.5,.75),c(0,1,fit$x))
      inv <- diff(gq[c(1,3)])/diff(qq[c(1,3)]); par <- c(qq[2]-gq[2]/inv,inv,fit$x)
      t <- pop$true; fd <- cdf_pdf(c(t$M-t$S,t$M+t$S,t$M),par)
      ratio <- variance_map(fd$pdf[1],fd$pdf[2],fd$pdf[3],1-sum(fd$cdf[1:2]))/t$v
      list(population=pop$name,objective=fit$fun,shape1=fit$x[1],shape2=fit$x[2],ratio=ratio,
           predicted_coverage=2*pnorm(MAD_Z*sqrt(ratio))-1,success=fit$success,nfev=fit$nfev)
    })
    rows[[length(rows)+1L]] <- row
  }
  ans <- mad_bind(rows); mad_csv(ans,file.path(out,"population_global_search.csv")); invisible(ans)
}

run_multistart_audit <- function(out,reps=10L,ns=c(200L,1000L),pop_indices=c(7L,8L),starts=10L) {
  cfg <- mad_config(out,"multistart",list(implementation_md5=optimizer_source_md5(),reps=reps,ns=ns,pop_indices=pop_indices,starts=starts))
  rows <- list(); pops <- populations()
  for(pi in pop_indices) for(n in ns) for(rep in seq_len(reps)-1L) {
    row <- mad_checkpoint(file.path(out,"checkpoints",sprintf("multistart_%d_%d_%d.rds",pi,n,rep)),cfg,function() {
      pop <- pops[[pi+1L]]; x <- mad_sample(pop,1L,pi,n,rep); base <- fit_pdq(x)
      y <- mad_shape_data(x); loss <- rowSums(sweep(grid_pdq(length(y$ps)),2L,y$y,"-")^2)
      best <- base$info$objective; bestpar <- base$par
      for(ix in head(order(loss),starts)) {
        fit <- nelder_mead(pdq_obj,SHAPE_GRID[ix,],args=list(y$y,y$ps),maxiter=2000L,xatol=1e-9,fatol=1e-11)
        if(fit$fun<best) { best <- fit$fun; bestpar <- quartile_params(x,fit$x) }
      }
      v0 <- gld_estimate(x,base$par)$value; v1 <- gld_estimate(x,bestpar)$value
      list(population=pop$name,n=n,rep=rep,base_objective=base$info$objective,
           refined_objective=best,variance_ratio_change=v1/v0-1)
    })
    rows[[length(rows)+1L]] <- row
  }
  ans <- mad_bind(rows); mad_csv(ans,file.path(out,"pdq_multistart_audit.csv")); invisible(ans)
}

run_pdq_restarts <- function(out,ns=c(200L,1000L,5000L,50000L),replications=1000L,results_dir=out) {
  cfg <- mad_config(out,"pdq_restarts",list(implementation_md5=optimizer_source_md5(),ns=ns,replications=replications,results_dir=normalizePath(results_dir,mustWork=FALSE)))
  pop <- populations()[[8L]]; truth <- pop$true; allrows <- list()
  for(n in ns) {
    key <- if(n %in% c(50,200,1000))1L else if(n==5000)4L else 5L
    archive_name <- if(key==1)sprintf("raw_7_%d.csv.gz",n) else if(key==4)"large_n_raw.csv.gz" else "resolution_raw.csv.gz"
    archive_path <- mad_record_path(results_dir,archive_name,required=FALSE)
    archive <- if(!is.null(archive_path)&&file.exists(archive_path))read.csv(archive_path) else data.frame()
    if(nrow(archive)) archive <- archive[archive$population=="Shoulder" & archive$method=="GLD-pdQ" & archive$n==n,,drop=FALSE]
    rows <- vector("list",replications)
    for(rep in seq_len(replications)-1L) rows[[rep+1L]] <- mad_checkpoint(
      file.path(out,"checkpoints",sprintf("restart_%d_%d.rds",n,rep)),cfg,function() {
        x <- mad_sample(pop,key,7L,n,rep); base <- fit_pdq(x); vo <- gld_estimate(x,base$par)$value
        ms <- mad_nuisance(x); y <- mad_shape_data(x)
        error <- NA_real_; rec <- if(nrow(archive))archive[archive$rep==rep,,drop=FALSE] else data.frame()
        if(nrow(rec)) {
          stopifnot(nrow(rec)==1L)
          error <- abs(vo/truth$v-rec$variance_ratio)
          if(error>2e-6 || ("shat" %in% names(rec) && abs(ms$s-rec$shat)>1e-10))
            stop("R regenerated reference mismatch n=",n," rep=",rep)
        }
        bestfun <- base$info$objective; bestshape <- base$par[3:4]; chosen <- "original"
        for(j in 1:2) {
          initial <- if(j==1)c(2.1,4.99) else c(4.99,2.1)
          fit <- nelder_mead(pdq_obj,initial,args=list(y$y,y$ps),maxiter=2000L,xatol=1e-9,fatol=1e-10)
          if(isTRUE(fit$success) && fit$fun<bestfun) {
            bestfun <- fit$fun; bestshape <- fit$x; chosen <- paste0("added_",j)
          }
        }
        par <- quartile_params(x,bestshape); vn <- gld_estimate(x,par)$value
        cover <- function(v) is.finite(v) && v>0 && abs(ms$s-truth$S)<=MAD_Z*sqrt(v/n)
        list(population="Shoulder",n=n,rep=rep,mhat=ms$m,shat=ms$s,
          original_objective=base$info$objective,augmented_objective=bestfun,objective_gain=base$info$objective-bestfun,
          original_v=vo,augmented_v=vn,original_vratio=vo/truth$v,augmented_vratio=vn/truth$v,
          original_cover=isTRUE(base$info$success)&&cover(vo),augmented_cover=cover(vn),
          original_l3=base$par[3],original_l4=base$par[4],augmented_l3=par[3],augmented_l4=par[4],
          selected=chosen,archived_vratio_difference=error)
      })
    d <- mad_bind(rows); mad_csv(d,file.path(out,sprintf("shoulder_restart_raw_%d.csv.gz",n))); allrows[[as.character(n)]] <- d
  }
  summaries <- lapply(allrows,function(d) {
    delta <- as.numeric(d$augmented_cover)-as.numeric(d$original_cover)
    list(n=d$n[1],R=nrow(d),original_coverage=mean(d$original_cover),augmented_coverage=mean(d$augmented_cover),
         paired_difference=mean(delta),paired_mcse=sd(delta)/sqrt(nrow(d)),
         original_median_vratio=median(d$original_vratio),augmented_median_vratio=median(d$augmented_vratio),
         improvement_count=sum(d$objective_gain>1e-6),upper_bound_count=sum(pmax(d$augmented_l3,d$augmented_l4)>4.99),
         max_reproduction_error=if(all(is.na(d$archived_vratio_difference)))NA_real_ else max(d$archived_vratio_difference,na.rm=TRUE))
  })
  summary <- mad_bind(summaries); mad_csv(summary,file.path(out,"shoulder_restart_summary.csv")); invisible(summary)
}

tm_obj_wide <- function(par,z) {
  if(any(!is.finite(par)) || par[2]<=1e-7 || par[2]>1e7 || any(par[3:4]< -6|par[3:4]>12) || abs(par[1])>1e5)return(1e12)
  spacing <- diff(c(0,cdf_pdf(z,par)$cdf,1))
  if(any(!is.finite(spacing)) || any(spacing<=1e-15))return(1e12)
  -sum(log(spacing))
}

mad_tm_wide_refit <- function(x,base=NULL) {
  if(is.null(base)) {pdq<-fit_pdq(x);base<-fit_tm(x,pdq$par)}
  ms <- mad_nuisance(x); iqr <- diff(quantile(x,c(.25,.75),type=7,names=FALSE))
  if(!is.finite(iqr)||iqr<=0||is.null(base$par))stop("Invalid sample or original fit for TM refit")
  xs <- sort((x-ms$m)/iqr); z <- (head(xs,-1)+tail(xs,-1))/2
  initial <- c((base$par[1]-ms$m)/iqr,base$par[2]*iqr,base$par[3:4])
  pushed <- initial; pushed[3:4] <- ifelse(pushed[3:4]>5.99,6.5,ifelse(pushed[3:4]< -2.99,-3.5,pushed[3:4]))
  bestpar <- initial; bestfun <- tm_obj_wide(initial,z); success <- FALSE
  for(st in list(initial,pushed)) {
    fit <- nelder_mead(tm_obj_wide,st,args=list(z),maxiter=10000L,xatol=2e-7,fatol=1e-7)
    if(isTRUE(fit$success)&&fit$fun<=bestfun) {bestfun<-fit$fun;bestpar<-fit$x;success<-TRUE}
  }
  par <- c(ms$m+iqr*bestpar[1],bestpar[2]/iqr,bestpar[3:4])
  list(par=par,info=list(objective=bestfun,improved_converged=success,
                         expanded_bound_hit=any(par[3:4]< -5.99|par[3:4]>11.99)))
}

run_tm_boundary <- function(out,results_dir=out,pop_indices=c(1L,4L),ns=c(50L,200L,1000L),max_refits=Inf) {
  cfg <- mad_config(out,"tm_boundary",list(implementation_md5=optimizer_source_md5(),results_dir=normalizePath(results_dir,mustWork=TRUE),
    pop_indices=pop_indices,ns=ns,max_refits=max_refits))
  pops <- populations(); records <- list(); summaries <- list(); completed <- 0L
  for(pi in pop_indices) for(n in ns) {
    pop <- pops[[pi+1L]]; truth <- pop$true
    path <- mad_record_path(results_dir,sprintf("raw_%d_%d.csv.gz",pi,n),required=TRUE)
    raw <- read.csv(path); raw <- raw[raw$method=="GLD-TM",,drop=FALSE]
    if(!nrow(raw))next
    boundary <- mad_bool(raw$boundary); corrected <- as.numeric(mad_bool(raw$cover)); local <- list()
    for(i in which(boundary)) {
      if(completed>=max_refits)break
      rep <- raw$rep[i]
      row <- mad_checkpoint(file.path(out,"checkpoints",sprintf("tmwide_%d_%d_%d.rds",pi,n,rep)),cfg,function() {
        x <- mad_sample(pop,1L,pi,n,rep); pdq <- fit_pdq(x); base <- fit_tm(x,pdq$par)
        vo <- gld_estimate(x,base$par)$value; error <- abs(vo-raw$vhat[i])/max(abs(raw$vhat[i]),1e-10)
        if(!is.finite(error) || error>2e-5)stop("R TM reproduction mismatch ",pi,",",n,",",rep)
        ms <- mad_nuisance(x); widened <- mad_tm_wide_refit(x,base)
        par <- widened$par;bestfun<-widened$info$objective;success<-widened$info$improved_converged
        vn <- gld_estimate(x,par)$value
        list(population=pop$name,n=n,rep=rep,original_objective=base$info$objective,
          expanded_objective=bestfun,gain=base$info$objective-bestfun,original_v=vo,expanded_v=vn,
          relative_variance_change=(vn-vo)/vo,original_cover=mad_bool(raw$cover[i]),
          expanded_cover=is.finite(vn)&&vn>0&&abs(ms$s-truth$S)<=MAD_Z*sqrt(vn/n),
          original_l3=base$par[3],original_l4=base$par[4],expanded_l3=par[3],expanded_l4=par[4],
          expanded_bound_hit=any(par[3:4]< -5.99|par[3:4]>11.99),improved_converged=success,reproduction_error=error)
      })
      local[[length(local)+1L]] <- row; records[[length(records)+1L]] <- row
      corrected[i] <- as.numeric(row$expanded_cover); completed <- completed+1L
    }
    g <- mad_bind(local)
    summaries[[length(summaries)+1L]] <- list(population=pop$name,n=n,R=nrow(raw),boundary_count=sum(boundary),
      refit_count=nrow(g),all_flagged_refit=nrow(g)==sum(boundary),original_coverage=mean(mad_bool(raw$cover)),
      boundary_refit_coverage=mean(corrected),changed_coverage_count=sum(corrected!=as.numeric(mad_bool(raw$cover))),
      median_relative_variance_change=if(nrow(g))median(g$relative_variance_change) else NA_real_,
      max_absolute_relative_variance_change=if(nrow(g))max(abs(g$relative_variance_change)) else NA_real_,
      max_gain=if(nrow(g))max(g$gain) else NA_real_,expanded_bound_hits=if(nrow(g))sum(g$expanded_bound_hit) else 0L)
    mad_csv(mad_bind(records),file.path(out,"tm_boundary_refit_raw.csv"))
  }
  ans <- mad_bind(summaries);mad_csv(ans,file.path(out,"tm_boundary_refit_summary.csv"));invisible(ans)
}

run_saved_checks <- function(out,results_dir=out) {
  paths <- list.files(results_dir,"^raw_[0-9]+_[0-9]+\\.csv\\.gz$",full.names=TRUE,recursive=TRUE)
  if(anyDuplicated(basename(paths)))stop("Ambiguous duplicate principal records")
  if(!length(paths))stop("No principal-study raw files in ",results_dir)
  summary <- hybrid <- list()
  for(path in paths) {
    d <- read.csv(path); tm <- d[d$method=="GLD-TM",,drop=FALSE]
    if(nrow(tm) && tm$population[1] %in% c("Exponential","Pareto7")) {
      boundary <- mad_bool(tm$boundary);cover<-mad_bool(tm$cover)
      summary[[length(summary)+1L]] <- list(population=tm$population[1],n=tm$n[1],R=nrow(tm),
        boundary_count=sum(boundary),boundary_fraction=mean(boundary),coverage=mean(cover),
        boundary_coverage=mad_mean(cover[boundary]),interior_coverage=mad_mean(cover[!boundary]),
        boundary_median_vratio=mad_median(tm$variance_ratio[boundary]),boundary_replicates=paste(tm$rep[boundary],collapse=";"))
    }
    for(method in c("GLD-pdQ","GLD-pdQ-ECDF")) {
      g<-d[d$method==method,,drop=FALSE]
      if(nrow(g))hybrid[[length(hybrid)+1L]]<-list(population=g$population[1],n=g$n[1],method=method,R=nrow(g),
        coverage=mean(mad_bool(g$cover)),median_vratio=mad_median(g$variance_ratio),
        median_a=mad_median(g$a),median_b=mad_median(g$b),median_c=mad_median(g$c),median_q=mad_median(g$q))
    }
  }
  mad_csv(mad_bind(summary),file.path(out,"tm_boundary_counts.csv"));mad_csv(mad_bind(hybrid),file.path(out,"hybrid_diagnostic.csv"))
  localpath<-mad_record_path(results_dir,"population_pdq_calibration.csv",required=FALSE)
  altpath<-mad_record_path(results_dir,"population_global_search.csv",required=FALSE)
  if(!is.null(localpath)&&!is.null(altpath)&&file.exists(localpath)&&file.exists(altpath)) {
    pop<-populations()[[8L]];v<-pop$true$v;local<-read.csv(localpath);local<-subset(local,population=="Shoulder"&J==50)
    alt<-read.csv(altpath);alt<-subset(alt,population=="Shoulder")
    if(nrow(local)==1L&&nrow(alt)==1L) {
      qq<-pop$ppf(c(.25,.5,.75));gq<-qvec(c(.25,.5,.75),c(0,1,alt$shape1,alt$shape2))
      iv<-diff(gq[c(1,3)])/diff(qq[c(1,3)]);loc<-qq[2]-gq[2]/iv
      rows<-data.frame(solution=c("Retained ten-start local fit","Differential-evolution alternative"),
        criterion="50-point population pdQ squared discrepancy",objective=c(local$objective,alt$objective),
        lambda1=c(local$lambda1,loc),lambda2=c(local$lambda2,iv),lambda3=c(local$lambda3,alt$shape1),
        lambda4=c(local$lambda4,alt$shape2),variance_ratio=c(local$variance_ratio,alt$ratio),
        fitted_variance=c(local$variance_ratio,alt$ratio)*v,predicted_coverage=c(local$limiting_coverage,alt$predicted_coverage))
      mad_csv(rows,file.path(out,"shoulder_population_solutions.csv"))
    }
  }
  invisible(list(boundary=mad_bind(summary),hybrid=mad_bind(hybrid)))
}
