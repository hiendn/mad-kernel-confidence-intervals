# Rebuild derived numeric summaries and figures in R from saved study records.
# Input records retain their provenance. Running this script does not rerun a study.
asset_bool <- function(x, name = "binary column") {
  if (is.logical(x)) return(x)
  y <- tolower(as.character(x))
  if (anyNA(y) || any(!y %in% c("true", "false", "1", "0")))
    stop("Nonbinary values in ", name)
  y %in% c("true", "1")
}

asset_rbind <- function(frames) {
  columns <- unique(unlist(lapply(frames, names)))
  do.call(rbind, lapply(frames, function(d) {
    for (nm in setdiff(columns, names(d))) d[[nm]] <- NA
    d[columns]
  }))
}

asset_tex <- function(d, path) {
  # These are numeric/data rows only, with no generated scientific narrative.
  escape <- function(x) gsub("_", "\\_", as.character(x), fixed = TRUE)
  fmt <- function(x) if (is.numeric(x)) {
    ifelse(is.na(x), "", formatC(x, digits = 6, format = "fg", flag = "#"))
  } else escape(x)
  rows <- apply(as.data.frame(lapply(d, fmt), stringsAsFactors = FALSE), 1,
                function(row) paste0(paste(row, collapse = " & "), " \\\\"))
  writeLines(c(paste0("\\begin{tabular}{", paste(rep("l", ncol(d)), collapse = ""), "}"),
               "\\toprule", paste0(paste(escape(names(d)), collapse = " & "), " \\\\") ,
               "\\midrule", rows, "\\bottomrule", "\\end{tabular}"), path)
}

asset_rank_table <- function(number,d) {
  highlight<-function(rows,columns,target) {
    for(row in rows) {
      values<-round(as.numeric(d[row,columns])*1000)
      if(all(is.na(values)))next
      closest<-which(abs(values-target)==min(abs(values-target),na.rm=TRUE))
      for(column in columns[closest])d[row,column]<<-paste0("\\textbf{",d[row,column],"}")
    }
  }
  if(number %in% c(1L,10L))highlight(seq_len(nrow(d)),c(3L,4L,5L,6L,7L,9L),950)
  if(number==2L) {
    keys<-unique(d[c(1,2)])
    for(i in seq_len(nrow(keys))) {
      rows<-which(d[[1]]==keys[[1]][i]&d[[2]]==keys[[2]][i]&d[[3]]!="Oracle")
      for(column in c(4L,6L)) {
        target<-if(column==4L)950 else 1000
        values<-round(as.numeric(d[rows,column])*1000)
        closest<-rows[which(abs(values-target)==min(abs(values-target),na.rm=TRUE))]
        d[closest,column]<-paste0("\\textbf{",d[closest,column],"}")
      }
    }
  }
  if(number==7L)for(scenario in unique(d[[1]]))for(column in c(4L,5L)) {
    rows<-which(d[[1]]==scenario&d[[2]]!="Oracle")
    values<-round(as.numeric(d[rows,column])*1000)
    closest<-rows[which(abs(values-950)==min(abs(values-950),na.rm=TRUE))]
    d[closest,column]<-paste0("\\textbf{",d[closest,column],"}")
  }
  d
}

asset_check_table_maps <- function(frames,output_dir,verify_values=FALSE,
                                   manuscript_root="manuscript") {
  source_files<-c(rep("tables/results_text.tex",4),"tables/shoulder_selection.tex",
    rep("tables/results_text.tex",2),rep("tables/application_text.tex",2),rep("tables/appendix_tables.tex",4))
  labels<-c("tab:maincover","tab:large","tab:popsolutions","tab:shoulderrestart","tab:shoulderselection",
    "tab:clt","tab:two","tab:galaxy","tab:galaxydensities","tab:logcoverage","tab:length","tab:truth","tab:tmbounds")
  expected_rows<-c(30L,14L,2L,4L,3L,6L,20L,6L,6L,30L,20L,10L,6L)
  expected_columns<-c(9L,6L,7L,7L,6L,6L,5L,5L,5L,9L,8L,8L,7L)
  scope_columns<-list(c(1L,2L),c(1L,2L,3L),1L,1L,1L,c(1L,2L),c(1L,2L,3L),1L,1L,
                      c(1L,2L),c(1L,2L),1L,c(1L,2L))
  clean<-function(x) {
    x<-gsub("\\textbf{","",x,fixed=TRUE);x<-gsub("}","",x,fixed=TRUE)
    x<-gsub("$","",x,fixed=TRUE);x<-gsub("\\chi^2_5","Chi-square5",x,fixed=TRUE)
    x<-gsub("Pareto(7)","Pareto7",x,fixed=TRUE);gsub("[[:space:]]+","",x)
  }
  checks<-lapply(seq_along(labels),function(i) {
    path<-file.path(manuscript_root,source_files[i]);lines<-readLines(path,warn=FALSE)
    start<-which(grepl(paste0("\\label{",labels[i],"}"),lines,fixed=TRUE))[1]
    if(is.na(start))stop("Current manuscript table label missing: ",labels[i])
    endings<-which(grepl("\\end{tabular}",lines,fixed=TRUE)|grepl("\\end{longtable}",lines,fixed=TRUE))
    ending<-endings[endings>=start][1]
    rows<-lapply(lines[start:ending],function(line) {
      if(!grepl(" & ",line,fixed=TRUE))return(NULL)
      line<-strsplit(line,"\\\\",fixed=TRUE)[[1]][1]
      cells<-trimws(strsplit(line,"&",fixed=TRUE)[[1]])
      last<-clean(tail(cells,1))
      if(!grepl("^-?[0-9]+(\\.[0-9]+)?$|^\\[[-0-9.,]+\\]$",last))return(NULL)
      cells
    })
    rows<-Filter(Negate(is.null),rows)
    actual<-frames[[i]]
    if(length(rows)!=expected_rows[i]||nrow(actual)!=expected_rows[i]||
       any(lengths(rows)!=expected_columns[i])||ncol(actual)!=expected_columns[i])
      stop("Current table scope mismatch for ",labels[i],": source ",length(rows)," generated ",nrow(actual))
    old<-do.call(rbind,rows);generated<-as.matrix(actual)
    # Ordering and first-column scientific objects must match even for new R runs.
    if(!identical(unname(clean(old[,scope_columns[[i]],drop=FALSE])),
                  unname(clean(generated[,scope_columns[[i]],drop=FALSE]))))
      stop("Current table row objects/order differ for ",labels[i])
    matches<-all(clean(old)==clean(generated))
    source_bold<-matrix(grepl("\\textbf{",old,fixed=TRUE),nrow=nrow(old))
    ranked<-as.matrix(asset_rank_table(i,actual))
    generated_bold<-matrix(grepl("\\textbf{",ranked,fixed=TRUE),nrow=nrow(ranked))
    bold_matches<-identical(source_bold,generated_bold)
    if(verify_values&&!matches) {
      where<-which(matrix(clean(old)!=clean(generated),nrow=nrow(old)),arr.ind=TRUE)
      stop("Current-table values differ for ",labels[i]," at row ",where[1,1]," column ",where[1,2],
        ": source ",old[where[1,1],where[1,2]],"; R ",generated[where[1,1],where[1,2]])
    }
    if(verify_values&&!bold_matches)stop("Current-table ranking differs for ",labels[i])
    data.frame(table=i,label=labels[i],source=path,source_rows=length(rows),generated_rows=nrow(actual),
      source_columns=expected_columns[i],generated_columns=ncol(actual),
      scope_and_row_order_match=TRUE,display_values_match=if(verify_values)matches else NA,
      highlight_positions_match=if(verify_values)bold_matches else NA,
      validation=if(verify_values)"All displayed scientific values and rankings" else "Source scope and row order")
  })
  write.csv(do.call(rbind,checks),file.path(output_dir,"table_map_checks.csv"),row.names=FALSE)
  invisible(checks)
}

rebuild_assets <- function(root = "computation/generated_r/full/finalize/canonical",
                           output_dir = "computation/generated_r/assets",
                           methods_file = "computation/R/mad_methods.R",
                           source_label = "R study records",
                           check_manuscript = FALSE) {
  root <- normalizePath(root, mustWork = TRUE)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  figures <- file.path(output_dir, "figures")
  tables <- file.path(output_dir, "tables")
  dir.create(figures, showWarnings = FALSE)
  dir.create(tables, showWarnings = FALSE)
  inputs <- character()
  load <- function(rel) {
    path <- file.path(root, rel)
    inputs <<- unique(c(inputs, path))
    read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  }
  selection <- do.call(rbind, lapply(c(200, 1000, 5000, 50000), function(n) {
    d <- load(sprintf("revision_results/shoulder_restart_raw_%s.csv.gz", n))
    stopifnot(nrow(d)>0L, !anyDuplicated(d$rep))
    stopifnot(nrow(d)==1000L)
    high <- pmax(d$augmented_l3, d$augmented_l4) > 4.99
    groups <- list(All = rep(TRUE, nrow(d)), `Upper shape region` = high,
                   `Other selected fits` = !high)
    do.call(rbind, lapply(names(groups), function(group) {
      sub <- d[groups[[group]], , drop = FALSE]
      qs <- if(nrow(sub))quantile(sub$augmented_vratio, c(.1,.25,.5,.75,.9), type = 7, names = FALSE) else rep(NA_real_,5)
      covered <- asset_bool(sub$augmented_cover, "augmented_cover")
      data.frame(n = n, group = group, count = nrow(sub), fraction = nrow(sub)/nrow(d),
                 q10 = qs[1], q25 = qs[2], median = qs[3], q75 = qs[4], q90 = qs[5],
                 descriptive_covered = sum(covered), descriptive_coverage = mean(covered))
    }))
  }))
  rownames(selection) <- NULL
  write.csv(selection, file.path(output_dir, "shoulder_selected_fit_distribution.csv"), row.names = FALSE)
  selected_50000 <- subset(selection, n == 50000)
  asset_tex(selected_50000[c("group","count","fraction","q10","median","q90")],
            file.path(tables, "shoulder_selection_rows.tex"))

  specs <- data.frame(population = c("Shoulder","Shoulder","Shoulder","Bimodal"),
                      index = c(7,7,7,8), n = c(200,1000,1000,50),
                      other = c("KDE-R","KDE-R","KDE-half","Oracle"))
  paired <- do.call(rbind, lapply(seq_len(nrow(specs)), function(i) {
    sp <- specs[i, ]
    d <- load(sprintf("results/raw_%s_%s.csv.gz", sp$index, sp$n))
    a <- d[d$method == "Bootstrap", c("rep","cover")]
    b <- d[d$method == sp$other, c("rep","cover")]
    stopifnot(!anyDuplicated(a$rep), !anyDuplicated(b$rep))
    common <- merge(a,b,by="rep",suffixes=c("_boot","_other"))
    stopifnot(nrow(common)>0L)
    stopifnot(nrow(common)==1000L)
    delta <- as.numeric(asset_bool(common$cover_boot)) - as.numeric(asset_bool(common$cover_other))
    data.frame(population = sp$population, n = sp$n,
                contrast = paste("Bootstrap minus", sp$other), R = length(delta),
                difference = mean(delta), paired_mcse = sd(delta)/sqrt(length(delta)))
  }))
  write.csv(paired, file.path(output_dir, "paired_bootstrap_comparisons.csv"), row.names = FALSE)

  primary <- load("results/one_sample_summary.csv")
  large <- load("results/large_n_summary.csv")
  resolution <- load("results/resolution_summary.csv")
  restart <- load("revision_results/shoulder_restart_summary.csv")
  combined <- asset_rbind(list(primary, large, resolution))
  methods <- c("KDE-R","KDE-half","GLD-pdQ","Oracle","Bootstrap")
  plotdata <- combined[combined$population %in% c("Shoulder","Bimodal") &
                        combined$method %in% methods,
                      c("population","n","method","R","coverage","mcse")]
  stopifnot(!anyDuplicated(plotdata[c("population","n","method")]))
  extra <- data.frame(population = "Shoulder", n = restart$n,
                      method = "GLD-pdQ + starts", R = restart$R,
                      coverage = restart$augmented_coverage)
  extra$mcse <- sqrt(extra$coverage*(1-extra$coverage)/extra$R)
  plotdata <- rbind(plotdata, extra)
  plotdata <- plotdata[order(plotdata$population, plotdata$method, plotdata$n), ]
  rownames(plotdata) <- NULL
  write.csv(plotdata, file.path(output_dir, "mixture_figure_data.csv"), row.names = FALSE)
  cols <- c("#1f77b4","#ff7f0e","#2ca02c","#d62728","#9467bd","#8c564b")
  points <- c(16,15,17,18,4,8)
  labels <- c("KDE-R", "KDE-half", "pdQ (original)", "Oracle", "Bootstrap (399)", "pdQ + starts")
  for (pop in c("Shoulder","Bimodal")) {
    pdf(file.path(figures, paste0("coverage_", tolower(pop), ".pdf")), width=7.2, height=4.35,
        family="Times", title=paste(pop, "coverage from",source_label))
    par(mar=c(4,4,1,1), cex=.88)
    plot(NA, xlim=c(38,65000), ylim=c(.745,1.008), log="x", xaxt="n", yaxt="n",
         xlab="Sample size (log scale)", ylab="Estimated coverage")
    axis(1, at=c(50,200,1000,5000,50000), labels=c("50","200","1000","5000","50000"))
    axis(2, at=seq(.75,1,.05), las=1)
    abline(h=.95, lty=2, col="grey40")
    chosen <- c(methods, if(pop == "Shoulder") "GLD-pdQ + starts")
    for (j in seq_along(chosen)) {
      d <- plotdata[plotdata$population==pop & plotdata$method==chosen[j], ]
      d <- d[order(d$n), ]
      radius <- qnorm(.975)*d$mcse
      low <- pmax(0,d$coverage-radius); high <- pmin(1,d$coverage+radius)
      if (!chosen[j] %in% c("Bootstrap","GLD-pdQ + starts","KDE-half"))
        lines(d$n, d$coverage, col=cols[j], lwd=1.1)
      visible_bar <- is.finite(low) & is.finite(high) & high > low
      if(any(visible_bar))arrows(d$n[visible_bar],low[visible_bar],d$n[visible_bar],high[visible_bar],
                                angle=90,code=3,length=.035,col=cols[j],lwd=.7)
      points(d$n,d$coverage,col=cols[j],pch=points[j],cex=.8)
    }
    legend("bottom", legend=labels[seq_along(chosen)], col=cols[seq_along(chosen)],
           pch=points[seq_along(chosen)], lty=c(1,NA,1,1,NA,NA)[seq_along(chosen)],
           ncol=3, cex=.77, bg="white", inset=.01)
    dev.off()
  }

  # Numeric table source exports preserve current scientific values without
  # recreating the older generator's obsolete narrative and removed sections.
  sources <- c("results/one_sample_summary.csv", "results/large_n_summary.csv",
    "results/resolution_summary.csv", "revision_results/shoulder_population_solutions.csv",
    "revision_results/shoulder_restart_summary.csv", "results/variance_clt_check.csv",
    "results/two_sample_summary.csv", "results/galaxies_analysis.csv",
    "results/population_truth.csv", "revision_results/tm_boundary_counts.csv",
    "revision_results/tm_boundary_refit_summary.csv", "revision_results/hybrid_diagnostic.csv")
  for (rel in sources) {
    d <- load(rel)
    stem <- sub("\\.csv$", "", basename(rel))
    write.csv(d, file.path(tables, paste0(stem,".csv")), row.names=FALSE)
  }
  # Table row files follow the thirteen current-paper computations. Captions and
  # surrounding interpretation remain manuscript-owned.
  publication_tables<-list()
  f<-function(x,digits=3)ifelse(is.na(x),"---",formatC(x,digits=digits,format="f"))
  integer_text<-function(x)formatC(x,digits=0,format="f")
  pop_order<-c("Normal","Exponential","Lognormal","Chi-square5","Pareto7","Cauchy","GLD","Shoulder","Bimodal","Contaminated")
  table_rows <- function(number,d) {
    publication_tables[[number]]<<-d
    asset_tex(asset_rank_table(number,d),file.path(tables,sprintf("table_%02d_rows.tex",number)))
    write.csv(d,file.path(tables,sprintf("table_%02d_display.csv",number)),row.names=FALSE)
  }
  cover_rows <- function(value,selected=primary,methods=c("KDE-R","KDE-NR","KDE-half","GLD-pdQ","GLD-TM","Oracle","Bootstrap")) {
    keys<-unique(selected[c("population","n")]);keys<-keys[order(match(keys$population,pop_order),keys$n),]
    do.call(rbind,lapply(seq_len(nrow(keys)),function(i) {
      k<-keys[i,];s<-selected[selected$population==k$population&selected$n==k$n,]
      values<-vapply(methods,
        function(m){z<-s[s$method==m,value];if(length(z))z[1] else NA_real_},numeric(1))
      data.frame(Population=k$population,n=integer_text(k$n),as.list(f(values)),check.names=FALSE)
    }))
  }
  table_rows(1,cover_rows("coverage"))
  d<-combined[combined$population %in% c("Shoulder","Bimodal")&combined$n %in% c(5000,50000)&
    combined$method %in% c("KDE-R","KDE-half","GLD-pdQ","Oracle"),]
  d<-d[order(match(d$population,c("Shoulder","Bimodal")),d$n,match(d$method,c("KDE-R","KDE-half","GLD-pdQ","Oracle"))),]
  table_rows(2,data.frame(Population=d$population,n=integer_text(d$n),Method=d$method,
    Coverage=f(d$coverage),MCSE=f(d$mcse),`Median variance ratio`=f(d$median_vratio),check.names=FALSE))
  d<-load("revision_results/shoulder_population_solutions.csv")
  table_rows(3,data.frame(Solution=c("Retained local","Alternative"),Objective=f(d$objective,4),
    Shape3=f(d$lambda3,4),Shape4=f(d$lambda4,4),`Fitted variance`=f(d$fitted_variance,4),
    `Variance ratio`=f(d$variance_ratio,4),Transform=f(d$predicted_coverage,4),check.names=FALSE))
  restart<-restart[order(restart$n),]
  table_rows(4,data.frame(n=integer_text(restart$n),Original=f(restart$original_coverage),
    `Added starts`=f(restart$augmented_coverage),Difference=f(restart$paired_difference),MCSE=f(restart$paired_mcse),
    `Original median ratio`=f(restart$original_median_vratio),`Added-start median ratio`=f(restart$augmented_median_vratio),check.names=FALSE))
  table_rows(5,data.frame(`Selected fits`=selected_50000$group,Count=integer_text(selected_50000$count),
    Fraction=f(selected_50000$fraction),`10th percentile`=f(selected_50000$q10),Median=f(selected_50000$median),
    `90th percentile`=f(selected_50000$q90),check.names=FALSE))
  d<-load("results/variance_clt_check.csv")
  d<-d[order(match(d$population,c("Normal","Lognormal")),d$n),]
  table_rows(6,data.frame(Population=d$population,n=integer_text(d$n),`Mean Z`=f(d$z_mean),`SD Z`=f(d$z_sd),
    `Central 95 percent`=f(d$central95),MSE=f(d$mse,5),check.names=FALSE))
  d<-load("results/two_sample_summary.csv")
  scenario_names<-c("Normal scale contrast","Lognormal equal MAD","Shoulder equal MAD","Bimodal equal MAD")
  scenario_labels<-c("Normal ratio 2","Lognormal/Exp.","Shoulder/Normal","Bimodal/Normal")
  methods_two<-c("KDE-R","KDE-half","GLD-pdQ","GLD-TM","Oracle")
  d<-d[d$method %in% methods_two,];d<-d[order(match(d$scenario,scenario_names),match(d$method,methods_two)),]
  table_rows(7,data.frame(Scenario=scenario_labels[match(d$scenario,scenario_names)],Method=d$method,R=integer_text(d$R),
    Difference=f(d$difference_coverage),Ratio=f(d$ratio_coverage),check.names=FALSE))
  d<-load("results/galaxies_analysis.csv")
  d<-d[order(match(d$method,c("KDE-R","KDE-NR","KDE-half","KDE-double","GLD-pdQ","GLD-TM"))),]
  table_rows(8,data.frame(Method=d$method,Variance=f(d$vhat),SE=f(d$se),
    `Additive interval`=paste0("[",f(d$lower),", ",f(d$upper),"]"),
    `Log interval`=paste0("[",f(d$loglower),", ",f(d$logupper),"]"),check.names=FALSE))
  table_rows(9,data.frame(Method=d$method,a=f(d$a,4),b=f(d$b,4),c=f(d$c,4),q=f(d$q,4),check.names=FALSE))
  table_rows(10,cover_rows("log_coverage"))
  table_rows(11,cover_rows("median_length",primary[primary$n %in% c(200,1000),],
    methods=c("KDE-R","KDE-half","GLD-pdQ","GLD-TM","Oracle","Bootstrap")))
  d<-load("results/population_truth.csv")
  d<-d[order(match(d$population,pop_order)),]
  table_rows(12,data.frame(Population=d$population,lapply(d[c("M","S","a","b","c","q","v")],f,digits=4),check.names=FALSE))
  d<-load("revision_results/tm_boundary_refit_summary.csv")
  d<-d[order(match(d$population,c("Exponential","Pareto7")),d$n),]
  table_rows(13,data.frame(Population=d$population,n=integer_text(d$n),Flagged=paste0(d$boundary_count,"/",d$R),
    `Median change (percent)`=f(100*d$median_relative_variance_change,2),
    `Max absolute change (percent)`=f(100*d$max_absolute_relative_variance_change,2),
    `New hits`=integer_text(d$expanded_bound_hits),Coverage=f(d$boundary_refit_coverage),check.names=FALSE))
  if (check_manuscript) asset_check_table_maps(publication_tables,tables)
  app <- load("results/galaxies_analysis.csv")
  gx <- load("sources/galaxies.csv")[[1]]/1000
  parpath <- file.path(root,"results/galaxies_gld_parameters.json")
  inputs <- unique(c(inputs,parpath))
  if (!requireNamespace("jsonlite",quietly=TRUE)) stop("jsonlite is needed to read saved GLD parameters")
  pars <- jsonlite::fromJSON(parpath)
  if (!exists("cdf_pdf", mode="function")) source(methods_file, local=environment())
  grid <- seq(8,36,length.out=1000)
  gaussian_at <- function(h) vapply(grid, function(t) mean(dnorm((t-gx)/h))/h, numeric(1))
  density_grid <- data.frame(velocity=grid)
  density_grid[["KDE-R"]] <- gaussian_at(app$h[app$method=="KDE-R"])
  density_grid[["KDE-half"]] <- gaussian_at(app$h[app$method=="KDE-half"])
  for (method in c("GLD-pdQ","GLD-TM")) {
    fit <- cdf_pdf(grid, as.numeric(pars[[method]]))
    density_grid[[method]] <- if(is.list(fit) && !is.null(fit$pdf)) fit$pdf else fit[[2]]
  }
  write.csv(density_grid,file.path(output_dir,"galaxies_density_grid.csv"),row.names=FALSE)
  pdf(file.path(figures,"galaxies_density.pdf"),width=6.6,height=3.6,family="Times",
      title=paste("Galaxy density evaluation from",source_label))
  par(mar=c(4,4,1,1),cex=.88)
  histogram <- hist(gx, breaks=seq(min(gx),max(gx),length.out=17), right=FALSE, plot=FALSE)
  ymax <- 1.04*max(histogram$density,unlist(density_grid[-1]))
  hist(gx, breaks=seq(min(gx),max(gx),length.out=17), freq=FALSE, right=FALSE,
       col=adjustcolor("grey60",.3),border="white",xlim=c(8,36),ylim=c(0,ymax),
       xlab="Velocity (1000 km/s)",ylab="Density",main="")
  for (j in seq_len(4)) lines(grid,density_grid[[j+1]],col=cols[j],lty=j,lwd=1.1)
  m <- median(gx); s <- median(abs(gx-m)); abline(v=c(m-s,m,m+s),lty=3,lwd=.7)
  legend("topright",legend=c("Observed histogram","KDE-R","KDE-half","GLD-pdQ","GLD-TM"),
         col=c("grey60",cols[1:4]),lty=c(NA,1:4),pch=c(15,NA,NA,NA,NA),cex=.77,ncol=2,bty="n")
  dev.off()
  manifest <- data.frame(path=inputs,md5=unname(tools::md5sum(inputs)))
  write.csv(manifest,file.path(output_dir,"input_md5.csv"),row.names=FALSE)
  writeLines(c(paste("Inputs:",source_label,"and fitted parameters."),
    "Computed in R: summaries, paired comparisons, figure coordinates and density evaluations.",
    "No simulation, bootstrap resampling or GLD optimization is run by this script.",
    "Quantiles: R type=7; group membership follows maximum augmented shape >4.99.",
    "Requires: base R, jsonlite, cdf_pdf from mad_methods.R.", paste("R:",R.version.string)),
    file.path(output_dir,"PROVENANCE.txt"))
  invisible(publication_tables)
}

asset_restart_base_checks <- function(canonical,output_dir) {
  checks<-lapply(c(200L,1000L,5000L,50000L),function(n) {
    restart<-read.csv(file.path(canonical,"revision_results",sprintf("shoulder_restart_raw_%d.csv.gz",n)))
    basename<-if(n %in% c(200,1000))sprintf("raw_7_%d.csv.gz",n) else
      if(n==5000)"large_n_raw.csv.gz" else "resolution_raw.csv.gz"
    base<-read.csv(file.path(canonical,"results",basename))
    base<-base[base$population=="Shoulder"&base$method=="GLD-pdQ"&base$n==n,]
    if(nrow(base)!=1000L||nrow(restart)!=1000L||anyDuplicated(base$rep)||anyDuplicated(restart$rep)||
       !setequal(base$rep,restart$rep))stop("Incomplete shared R shoulder datasets at n=",n)
    common<-merge(restart,base,by="rep",suffixes=c("_restart","_base"),sort=TRUE)
    valid<-is.finite(common$original_vratio)&is.finite(common$variance_ratio)
    max_ratio<-if(any(valid))max(abs(common$original_vratio[valid]-common$variance_ratio[valid])) else NA_real_
    coverage_disagreements<-sum(asset_bool(common$original_cover)!=asset_bool(common$cover))
    mad_difference<-if(all(c("shat_restart","shat_base") %in% names(common)))
      max(abs(common$shat_restart-common$shat_base)) else NA_real_
    if((is.finite(max_ratio)&&max_ratio>2e-6)||coverage_disagreements>0L||
       (is.finite(mad_difference)&&mad_difference>1e-10))stop("Shoulder restart/base mismatch at n=",n)
    data.frame(n=n,shared_datasets=nrow(common),finite_ratio_pairs=sum(valid),
      max_absolute_variance_ratio_difference=max_ratio,coverage_disagreements=coverage_disagreements,
      max_absolute_MAD_difference=mad_difference)
  })
  dir.create(output_dir,recursive=TRUE,showWarnings=FALSE)
  result<-do.call(rbind,checks)
  write.csv(result,file.path(output_dir,"restart_base_checks.csv"),row.names=FALSE)
  result
}

finalize_r_assets <- function(full_base,output_dir="computation/generated_r/finalized_R") {
  if(!exists("mad_record_path",mode="function"))stop("Source studies.R before finalizing")
  full_base<-normalizePath(full_base,mustWork=TRUE)
  dir.create(output_dir,recursive=TRUE,showWarnings=FALSE)
  output_dir<-normalizePath(output_dir,mustWork=TRUE)
  for(task in c("principal","population","galaxy","clt","large","resolution","restarts","tm","audit")) {
    taskroot<-file.path(full_base,task)
    if(identical(output_dir,taskroot)||startsWith(output_dir,paste0(taskroot,"/")))
      stop("Finalize output must be outside native task directories: ",taskroot)
  }
  canonical<-file.path(output_dir,"canonical")
  for(d in c("results","revision_results","sources"))dir.create(file.path(canonical,d),recursive=TRUE,showWarnings=FALSE)
  links<-list()
  link_record<-function(source,destination) {
    source<-normalizePath(source,mustWork=TRUE)
    if(file.exists(destination)) {
      existing<-Sys.readlink(destination)
      if(is.na(existing)||!nzchar(existing)||normalizePath(destination)!=source)
        stop("Existing canonical asset target differs: ",destination)
    } else if(!file.symlink(source,destination))stop("Could not link ",source)
    links[[length(links)+1L]]<<-data.frame(source=source,target=destination)
  }
  principal<-file.path(full_base,"principal")
  raw<-list.files(principal,"^raw_[0-9]+_[0-9]+\\.csv\\.gz$",recursive=TRUE,full.names=TRUE)
  expected<-as.vector(outer(0:9,c(50,200,1000),function(i,n)sprintf("raw_%d_%d.csv.gz",i,n)))
  if(anyDuplicated(basename(raw))||!setequal(basename(raw),expected))stop("Full principal run requires all30 distinct design files")
  frames<-lapply(raw,function(p) {
    d<-read.csv(p);counts<-table(d$method)
    if(!all(c("GLD-TM","Bootstrap","KDE-R","GLD-pdQ","Oracle") %in% names(counts)) ||
       any(counts[names(counts)!="GLD-TM"]!=1000L)||counts[["GLD-TM"]]!=200L)
      stop("Incomplete principal records: ",p)
    link_record(p,file.path(canonical,"results",basename(p)));d
  })
  df<-mad_bind(frames);mad_summarize(df,file.path(canonical,"results"))
  mad_two_sample(df,populations(),file.path(canonical,"results"))
  maps<-list(population=c("population_pdq_calibration.csv","population_truth.csv"),
    galaxy=c("galaxies_analysis.csv","galaxies_gld_parameters.json"),clt="variance_clt_check.csv",
    large=c("large_n_summary.csv","large_n_raw.csv.gz"),resolution=c("resolution_summary.csv","resolution_raw.csv.gz"))
  for(task in names(maps))for(nm in maps[[task]]) {
    source<-mad_record_path(file.path(full_base,task),nm,required=FALSE)
    if(!file.exists(source)&&nm=="population_truth.csv") {
      truth<-mad_bind(lapply(populations(),function(p)c(list(population=p$name,description=p$description),p$true)))
      mad_csv(truth,file.path(canonical,"results",nm))
    } else if(!file.exists(source)&&task=="galaxy"&&nm=="galaxies_gld_parameters.json") {
      # Portable releases retain the saved fit vectors as dput text.
      # Read those values without refitting, then check the JSON round trip.
      if(!requireNamespace("jsonlite",quietly=TRUE))stop("jsonlite is needed to serialize native R galaxy parameters")
      parameter_path<-mad_record_path(file.path(full_base,task),"galaxies_gld_parameters.rds",required=FALSE)
      if(file.exists(parameter_path)) {
        parameters<-readRDS(parameter_path)
      } else {
        parameter_path<-mad_record_path(file.path(full_base,task),"galaxies_gld_parameters.R",required=TRUE)
        parameters<-dget(parameter_path)
      }
      destination<-file.path(canonical,"results",nm)
      if(file.exists(destination)) {
        if(!isTRUE(all.equal(jsonlite::fromJSON(destination),parameters,tolerance=0,check.attributes=FALSE)))
          stop("Existing canonical galaxy parameters differ from the saved native R values: ",destination)
      } else jsonlite::write_json(parameters,destination,pretty=TRUE,auto_unbox=TRUE,digits=17)
      if(!isTRUE(all.equal(jsonlite::fromJSON(destination),parameters,tolerance=0,check.attributes=FALSE)))
        stop("Native galaxy parameter JSON failed exact round-trip verification")
      links[[length(links)+1L]]<-data.frame(source=normalizePath(parameter_path),target=destination)
    } else link_record(source,file.path(canonical,"results",nm))
  }
  for(spec in list(list(file="large_n_summary.csv",populations=5L,methods=3L),
                   list(file="resolution_summary.csv",populations=2L,methods=4L))) {
    d<-read.csv(file.path(canonical,"results",spec$file))
    if(nrow(d)!=spec$populations*spec$methods || length(unique(d$population))!=spec$populations || any(d$R!=1000L))
      stop("Incomplete native R long-sample study: ",spec$file)
  }
  clt<-read.csv(file.path(canonical,"results","variance_clt_check.csv"))
  if(nrow(clt)!=6L||any(clt$R!=3000L)||!setequal(clt$n,c(200,1000,5000)))
    stop("Incomplete native R variance-limit study")
  for(nm in c("shoulder_restart_summary.csv",sprintf("shoulder_restart_raw_%d.csv.gz",c(200,1000,5000,50000)))) {
    source<-mad_record_path(file.path(full_base,"restarts"),nm,required=TRUE)
    d<-read.csv(source)
    if(grepl("_raw_",nm)&&nrow(d)!=1000L)stop("Incomplete shoulder restart records: ",source)
    link_record(source,file.path(canonical,"revision_results",nm))
  }
  asset_restart_base_checks(canonical,output_dir)
  for(nm in c("tm_boundary_refit_summary.csv","tm_boundary_refit_raw.csv"))
    link_record(mad_record_path(file.path(full_base,"tm"),nm,required=TRUE),file.path(canonical,"revision_results",nm))
  tm<-read.csv(file.path(canonical,"revision_results","tm_boundary_refit_summary.csv"))
  if(nrow(tm)!=6L||any(tm$R!=200L)||("all_flagged_refit" %in% names(tm)&&!all(tm$all_flagged_refit)))
    stop("Incomplete native R TM boundary follow-up")
  # Derive optimizer and hybrid summaries from the saved R records.
  link_record(mad_record_path(file.path(full_base,"audit"),"population_global_search.csv",required=TRUE),
    file.path(canonical,"results","population_global_search.csv"))
  run_saved_checks(file.path(canonical,"revision_results"),file.path(canonical,"results"))
  link_record("computation/data/galaxies.csv",file.path(canonical,"sources","galaxies.csv"))
  mad_csv(mad_bind(links),file.path(output_dir,"source_links.csv"))
  rebuild_assets(canonical,file.path(output_dir,"assets"),
                 source_label="R simulation and fitting outputs; galaxy observations")
}

if (sys.nframe()==0L) {
  args <- commandArgs(trailingOnly=TRUE)
  root <- if(length(args)>=1) args[1] else "computation/generated_r/full/finalize/canonical"
  out <- if(length(args)>=2) args[2] else "computation/generated_r/assets"
  print(rebuild_assets(root,out))
}
