# Paired Monte Carlo uncertainty from saved R study records; no model refitting.
source('computation/R/studies.R')
paired <- function(d,groups,methods,quantities) {
 result<-list()
 keys<-unique(d[groups])
 for(i in seq_len(nrow(keys))) {
  use<-rep(TRUE,nrow(d));for(k in groups)use<-use & d[[k]]==keys[[k]][i]
  sub<-d[use,,drop=FALSE]
  for(pair in combn(methods,2,simplify=FALSE)) {
   a<-sub[sub$method==pair[1],,drop=FALSE];b<-sub[sub$method==pair[2],,drop=FALSE]
   common<-merge(a,b,by='rep',suffixes=c('_1','_2'));if(!nrow(common))next
   for(q in quantities) {
    x<-as.numeric(mad_bool(common[[paste0(q,'_1')]]))-as.numeric(mad_bool(common[[paste0(q,'_2')]]))
    result[[length(result)+1L]]<-c(as.list(keys[i,,drop=FALSE]),list(quantity=q,left=pair[1],right=pair[2],
      R=length(x),difference=mean(x),mcse=sd(x)/sqrt(length(x))))
   }
  }
 }
 mad_bind(result)
}
run_inference <- function(full='computation/generated_r/full',out=file.path(full,'inference'),medians=FALSE) {
 paths<-list.files(file.path(full,'principal'),'^raw_[0-9]+_[0-9]+[.]csv[.]gz$',recursive=TRUE,full.names=TRUE)
 stopifnot(length(paths)==30L);d<-mad_bind(lapply(sort(paths),read.csv));d$success<-mad_bool(d$success)
 methods<-c('KDE-R','KDE-NR','KDE-half','GLD-pdQ','GLD-TM','Oracle','Bootstrap')
 mad_csv(paired(d,c('population','n'),methods,c('cover','logcover')),file.path(out,'paired_principal_contrasts.csv'))
 truth<-read.csv(file.path(full,'principal/population_truth.csv'));truth<-setNames(truth$S,truth$population)
 specs<-list(list('Normal scale contrast','Normal','Normal',1,.5,50L,200L),
  list('Lognormal equal MAD','Lognormal','Exponential',1,1,200L,200L),
  list('Shoulder equal MAD','Shoulder','Normal',1,1,200L,1000L),
  list('Bimodal equal MAD','Bimodal','Normal',1,1,1000L,1000L))
 methods2<-c('KDE-R','KDE-NR','KDE-half','KDE-EIF','GLD-pdQ','GLD-TM','Oracle');indicators<-list()
 for(s in specs)for(method in methods2) {
  a<-d[d$population==s[[2]]&d$n==s[[6]]&d$method==method,,drop=FALSE]
  b<-d[d$population==s[[3]]&d$n==s[[7]]&d$method==method,,drop=FALSE]
  z<-merge(a,b,by='rep',suffixes=c('_1','_2'));k1<-s[[4]]/truth[[s[[2]]]];k2<-s[[5]]/truth[[s[[3]]]]
  e1<-z$shat_1*k1;e2<-z$shat_2*k2;v1<-z$vhat_1*k1^2;v2<-z$vhat_2*k2^2;ok<-z$success_1&z$success_2
  cd<-ok & abs(e1-e2-(s[[4]]-s[[5]]))<=MAD_Z*sqrt(v1/s[[6]]+v2/s[[7]])
  cr<-ok & abs(log(e1/e2)-log(s[[4]]/s[[5]]))<=MAD_Z*sqrt(v1/(s[[6]]*e1^2)+v2/(s[[7]]*e2^2))
  cd[is.na(cd)]<-FALSE;cr[is.na(cr)]<-FALSE
  indicators[[length(indicators)+1L]]<-data.frame(scenario=s[[1]],method=method,rep=z$rep,difference_coverage=cd,ratio_coverage=cr)
 }
 two<-mad_bind(indicators)
 mad_csv(paired(two,'scenario',methods2,c('difference_coverage','ratio_coverage')),file.path(out,'paired_two_sample_contrasts.csv'))
 if(medians) {
  RNGkind("L'Ecuyer-CMRG",normal.kind='Inversion',sample.kind='Rejection');set.seed(20261025L)
  ans<-list()
  for(g in split(d,interaction(d$population,d$n,d$method,drop=TRUE)))for(q in c('length','variance_ratio')) {
   x<-g[[q]][g$success];x<-x[is.finite(x)]
   se<-if(!length(x))NA_real_ else if(all(x==x[1]))0 else sd(replicate(400,median(x[sample.int(length(x),length(x),replace=TRUE)])))
   ans[[length(ans)+1L]]<-list(population=g$population[1],n=g$n[1],method=g$method[1],quantity=q,
    successful=length(x),median=mad_median(x),mcse=se,bootstrap_replications=400L)
  }
  mad_csv(mad_bind(ans),file.path(out,'principal_median_mcse.csv'))
 }
 invisible(out)
}
if(sys.nframe()==0L) {
 args<-commandArgs(trailingOnly=TRUE);run_inference(if(length(args))args[1] else 'computation/generated_r/full',
  if(length(args)>=2)args[2] else 'computation/generated_r/full/inference',medians='--medians'%in%args)
}
