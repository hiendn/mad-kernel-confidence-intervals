# Reproducible study drivers for the MAD manuscript.
# See RUNNING.md for seeds, checkpoint guards and completed study counts.

MAD_SEED <- 20260917L
MAD_Z <- qnorm(.975)
`%or%` <- function(x, y) if (is.null(x)) y else x
mad_elapsed <- function() unname(proc.time()[["elapsed"]])
mad_bool <- function(x) {
  if (is.logical(x)) return(x)
  y <- tolower(trimws(as.character(x))); absent <- is.na(y) | y %in% c("", "na", "nan", "none")
  if (any(!absent & !y %in% c("true", "false", "1", "0"))) stop("Invalid boolean value")
  result <- y %in% c("true", "1"); result[absent] <- NA; result
}
mad_mean <- function(x) if (length(x) && any(!is.na(x))) mean(x, na.rm = TRUE) else NA_real_
mad_median <- function(x) if (length(x) && any(!is.na(x))) median(x, na.rm = TRUE) else NA_real_
mad_bind <- function(rows) {
  if (!length(rows)) return(data.frame())
  rows <- lapply(rows, function(x) if (is.data.frame(x)) x else as.data.frame(x, stringsAsFactors = FALSE))
  fields <- unique(unlist(lapply(rows, names)))
  do.call(rbind, lapply(rows, function(x) { for (n in setdiff(fields, names(x))) x[[n]] <- NA; x[fields] }))
}
mad_csv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  con <- if (grepl("\\.gz$", path)) gzfile(path, "wt") else file(path, "wt")
  on.exit(close(con)); write.csv(x, con, row.names = FALSE, na = "")
  invisible(path)
}
mad_record_path <- function(root, basename, required = FALSE) {
  direct <- file.path(root, basename)
  if (file.exists(direct)) return(direct)
  candidates <- list.files(root, recursive = TRUE, full.names = TRUE)
  candidates <- candidates[base::basename(candidates) == basename]
  if (length(candidates) > 1L) stop("Ambiguous saved records for ", basename, " under ", root)
  if (!length(candidates)) {
    if (required) stop("Missing saved records ", basename, " under ", root)
    return(direct)
  }
  candidates
}

# Integer keys are folded in this documented order.
# The L'Ecuyer-CMRG algorithm is fixed explicitly.
mad_seed <- function(experiment, population, n, rep, master = MAD_SEED) {
  h <- as.double(master)
  for (key in c(experiment, population, n, rep)) h <- (48271 * h + key + 1) %% 2147483647
  as.integer(if (h == 0) 1 else h)
}
mad_rng <- function(experiment, population, n, rep, master = MAD_SEED) {
  RNGkind("L'Ecuyer-CMRG", normal.kind = "Inversion", sample.kind = "Rejection")
  set.seed(mad_seed(experiment, population, n, rep, master))
}
mad_sample <- function(pop, experiment, population, n, rep) {
  mad_rng(experiment, population, n, rep); pop$rvs(n)
}

# Each row is a resample of one-based indices; supplying indices makes
# bootstrap arithmetic deterministic.
bootstrap_from_indices <- function(x, indices) {
  indices <- as.matrix(indices); n <- length(x)
  if (ncol(indices) != n || nrow(indices) < 2L || anyNA(indices) ||
      any(indices != as.integer(indices)) || any(indices < 1L | indices > n))
    stop("indices must have >=2 rows and n valid one-based integer indices per row")
  sb <- apply(indices, 1L, function(index) { y <- x[index]; median(abs(y-median(y))) })
  list(mads = sb, endpoints = quantile(sb, c(.025, .975), type = 7, names = FALSE), root_n_variance = var(sb)*n)
}

# Checkpoints are one atomic RDS per replicate. Existing checkpoints are used
# only if their exact computation configuration agrees. Failures remain rows.
mad_checkpoint <- function(path, config, compute) {
  if (file.exists(path)) {
    old <- readRDS(path)
    if (!identical(old$config, config)) stop("Checkpoint configuration mismatch: ", path)
    return(old$value)
  }
  value <- compute(); dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- paste0(path, ".tmp-", Sys.getpid())
  saveRDS(list(config = config, value = value), tmp)
  if (!file.rename(tmp, path)) stop("Could not persist checkpoint ", path)
  value
}
mad_config <- function(out, task, args) {
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  cfg <- list(task = task, seed = MAD_SEED, rng = "L'Ecuyer-CMRG/Inversion/Rejection",
              seed_key_order = "master, experiment, zero-based population, n, zero-based replication",
              seed_fold = "h <- (48271*h + key + 1) mod 2147483647",
              language = R.version.string, platform = R.version$platform,
              source_md5 = tools::md5sum(c("computation/R/nelder_mead.R", "computation/R/mad_methods.R", "computation/R/studies.R")),
              args = args)
  path <- file.path(out, paste0(task, "_config.RDS"))
  if (file.exists(path) && !identical(readRDS(path), cfg)) stop("Output configuration differs: ", path)
  saveRDS(cfg, path); dput(cfg, file = file.path(out, paste0(task, "_config.R")))
  capture.output(sessionInfo(), file = file.path(out, paste0(task, "_session.txt")))
  cfg
}

mad_record <- function(method, v, x, truth, rep, seconds = 0, info = list()) {
  ms <- mad_nuisance(x); n <- length(x); s <- ms$s
  ok <- is.finite(v) && v > 0 && isTRUE(info$success %or% TRUE) && is.finite(s) && s > 0
  row <- list(rep = rep, method = method, n = n, shat = s, mhat = ms$m,
              true_s = truth$S, vhat = v, true_v = truth$v, success = ok, seconds = seconds)
  for (key in c("h", "a", "b", "c", "q", "nit", "boundary", "floor")) row[[key]] <- info[[key]] %or% NA
  if (ok) {
    se <- sqrt(v/n); lo <- s - MAD_Z*se; hi <- s + MAD_Z*se; rad <- MAD_Z*se/s
    loglo <- s*exp(-rad); loghi <- if (rad < 700) s*exp(rad) else Inf
    row <- c(row, list(cover = lo <= truth$S && truth$S <= hi,
      logcover = loglo <= truth$S && truth$S <= loghi,
      length = (hi-lo)/truth$S, loglength = (loghi-loglo)/truth$S,
      variance_ratio = v/truth$v, lower_negative = lo < 0))
  } else row <- c(row, list(cover = FALSE, logcover = FALSE, length = NA_real_,
    loglength = NA_real_, variance_ratio = NA_real_, lower_negative = FALSE))
  row
}

mad_principal_rep <- function(pop, pi, n, rep, rtm, B, include_gld = TRUE) {
  x <- mad_sample(pop, 1, pi, n, rep); rows <- list(); add <- function(row) rows[[length(rows)+1L]] <<- row
  for (spec in list(c("KDE-R", "robust", 1), c("KDE-NR", "nrd0", 1), c("KDE-half", "robust", .5))) {
    t <- mad_elapsed(); v <- kernel_estimate(x, spec[2], as.numeric(spec[3])); elapsed <- mad_elapsed()-t
    add(mad_record(spec[1], v$value, x, pop$true, rep, elapsed, v$info))
    if (spec[1] == "KDE-R") add(mad_record("KDE-EIF", v$info$eif, x, pop$true, rep, elapsed, v$info))
  }
  add(mad_record("Oracle", pop$true$v, x, pop$true, rep))
  if (include_gld) {
    t <- mad_elapsed(); fit <- fit_pdq(x); v <- gld_estimate(x, fit$par)
    meta <- modifyList(v$info, fit$info); elapsed <- mad_elapsed()-t
    add(mad_record("GLD-pdQ", v$value, x, pop$true, rep, elapsed, meta))
    vh <- gld_estimate(x, fit$par, TRUE)
    add(mad_record("GLD-pdQ-ECDF", vh$value, x, pop$true, rep, elapsed, modifyList(vh$info, fit$info)))
    if (rep < rtm) {
      t <- mad_elapsed(); ft <- fit_tm(x, fit$par); vt <- gld_estimate(x, ft$par)
      add(mad_record("GLD-TM", vt$value, x, pop$true, rep, mad_elapsed()-t, modifyList(vt$info, ft$info)))
    }
  }
  if (B > 0) {
    mad_rng(2, pi, n, rep); t <- mad_elapsed()
    index <- matrix(sample.int(n, B*n, replace = TRUE), nrow = B, byrow = TRUE)
    boot <- bootstrap_from_indices(x, index); sb <- boot$mads; ends <- boot$endpoints; ms <- mad_nuisance(x)
    add(list(rep = rep, method = "Bootstrap", n = n, shat = ms$s, mhat = ms$m,
      true_s = pop$true$S, vhat = var(sb)*n, true_v = pop$true$v, success = TRUE,
      seconds = mad_elapsed()-t, cover = ends[1] <= pop$true$S && pop$true$S <= ends[2],
      logcover = ends[1] <= pop$true$S && pop$true$S <= ends[2], length = diff(ends)/pop$true$S,
      loglength = diff(ends)/pop$true$S, variance_ratio = var(sb)*n/pop$true$v, lower_negative = FALSE))
  }
  result <- mad_bind(rows); result$population <- pop$name
  result[c("population", setdiff(names(result), "population"))]
}

mad_summarize <- function(df, out) {
  if (!nrow(df)) return(invisible(data.frame()))
  for (k in intersect(c("success", "cover", "logcover", "boundary"), names(df))) df[[k]] <- mad_bool(df[[k]])
  groups <- split(df, interaction(df$population, df$n, df$method, drop = TRUE))
  summary <- mad_bind(lapply(groups, function(g) {
    ok <- g$success; R <- nrow(g); cp <- mean(g$cover); lc <- mean(g$logcover)
    list(population = g$population[1], n = g$n[1], method = g$method[1], R = R,
      completed = sum(ok), coverage = cp, mcse = sqrt(cp*(1-cp)/R),
      log_coverage = lc, log_mcse = sqrt(lc*(1-lc)/R), conditional_coverage = mad_mean(g$cover[ok]),
      failed = 1-mean(ok), failsafe_coverage = cp + 1-mean(ok),
      mean_length = mad_mean(g$length[ok]), median_length = mad_median(g$length[ok]),
      median_loglength = mad_median(g$loglength[ok]), median_vratio = mad_median(g$variance_ratio[ok]),
      mean_vratio = mad_mean(g$variance_ratio[ok]), median_ms = 1000*median(g$seconds),
      mean_ms = 1000*mean(g$seconds), boundary = mad_mean(g$boundary))
  }))
  mad_csv(summary, file.path(out, "one_sample_summary.csv")); pairs <- list()
  for (g in split(df, interaction(df$population, df$n, drop = TRUE))) {
    a <- g[g$method == "KDE-R", c("rep", "cover")]
    for (other in c("GLD-pdQ", "GLD-TM", "Oracle")) {
      b <- g[g$method == other, c("rep", "cover")]; z <- merge(a, b, by = "rep")
      if (!nrow(z)) next
      d <- as.numeric(z$cover.x)-as.numeric(z$cover.y)
      pairs[[length(pairs)+1]] <- list(population = g$population[1], n = g$n[1],
        comparison = paste("KDE-R minus", other), R = length(d), difference = mean(d), mcse = sd(d)/sqrt(length(d)))
    }
  }
  mad_csv(mad_bind(pairs), file.path(out, "paired_coverage_differences.csv")); invisible(summary)
}

mad_two_sample <- function(df, pops, out) {
  truth <- setNames(lapply(pops, `[[`, "true"), vapply(pops, `[[`, "", "name"))
  scenarios <- list(list("Normal scale contrast", "Normal", "Normal", 1, .5, 50, 200),
    list("Lognormal equal MAD", "Lognormal", "Exponential", 1, 1, 200, 200),
    list("Shoulder equal MAD", "Shoulder", "Normal", 1, 1, 200, 1000),
    list("Bimodal equal MAD", "Bimodal", "Normal", 1, 1, 1000, 1000))
  rows <- list()
  for (s in scenarios) for (method in c("KDE-R", "KDE-NR", "KDE-half", "KDE-EIF", "GLD-pdQ", "GLD-TM", "Oracle")) {
    a <- df[df$population == s[[2]] & df$n == s[[6]] & df$method == method, ]
    b <- df[df$population == s[[3]] & df$n == s[[7]] & df$method == method, ]
    a <- a[order(a$rep), ]; b <- b[order(b$rep), ]; R <- min(nrow(a), nrow(b)); if (!R) next
    a <- a[seq_len(R), ]; b <- b[seq_len(R), ]
    if (!identical(a$rep, b$rep)) stop("Two-sample common-prefix replication IDs differ")
    k1 <- s[[4]]/truth[[s[[2]]]]$S; k2 <- s[[5]]/truth[[s[[3]]]]$S
    e1 <- a$shat*k1; e2 <- b$shat*k2; v1 <- a$vhat*k1^2; v2 <- b$vhat*k2^2
    ok <- mad_bool(a$success) & mad_bool(b$success); se <- sqrt(v1/s[[6]]+v2/s[[7]])
    logse <- sqrt(v1/(s[[6]]*e1^2)+v2/(s[[7]]*e2^2))
    cd <- ok & abs(e1-e2-(s[[4]]-s[[5]])) <= MAD_Z*se
    cr <- ok & abs(log(e1/e2)-log(s[[4]]/s[[5]])) <= MAD_Z*logse
    cd[is.na(cd)] <- FALSE; cr[is.na(cr)] <- FALSE
    rows[[length(rows)+1]] <- list(scenario = s[[1]], pop1 = s[[2]], pop2 = s[[3]], n1 = s[[6]], n2 = s[[7]],
      method = method, R = R, true_ratio = s[[4]]/s[[5]], difference_coverage = mean(cd), ratio_coverage = mean(cr),
      ratio_mcse = sqrt(mean(cr)*(1-mean(cr))/R), completed = sum(ok),
      median_difference_width = mad_median(2*MAD_Z*se[ok]), median_logratio_width = mad_median(2*MAD_Z*logse[ok]))
  }
  ans <- mad_bind(rows); mad_csv(ans, file.path(out, "two_sample_summary.csv")); invisible(ans)
}

run_principal <- function(out, replications = 1000L, tm_replications = 200L, bootstrap = 399L,
                          pop_indices = 0:9, ns = c(50L, 200L, 1000L), include_gld = TRUE) {
  stopifnot(replications > 0, tm_replications >= 0, bootstrap == 0 || bootstrap >= 2,
    all(pop_indices %in% 0:9), all(ns > 1))
  cfg <- mad_config(out, "principal", as.list(environment())[setdiff(ls(), c("out"))])
  pops <- populations(); truth <- mad_bind(lapply(pops, function(p) c(list(population = p$name, description = p$description), p$true)))
  mad_csv(truth, file.path(out, "population_truth.csv")); paths <- character()
  for (pi in pop_indices) for (n in ns) {
    t <- mad_elapsed(); pop <- pops[[pi+1L]]; all <- vector("list", replications)
    for (rep in seq_len(replications)-1L) {
      path <- file.path(out, "checkpoints", sprintf("principal_%02d_%d_%04d.rds", pi, n, rep))
      all[[rep+1L]] <- mad_checkpoint(path, cfg, function() mad_principal_rep(pop, pi, n, rep, tm_replications, bootstrap, include_gld))
    }
    path <- file.path(out, sprintf("raw_%d_%d.csv.gz", pi, n)); mad_csv(mad_bind(all), path); paths <- c(paths, path)
    cat(sprintf("Principal %s n=%d R=%d completed in %.2fs\n", pop$name, n, replications, mad_elapsed()-t)); flush.console()
  }
  df <- mad_bind(lapply(paths, read.csv)); mad_summarize(df, out); mad_two_sample(df, pops, out); invisible(df)
}

population_pdq <- function(pop, J = 50L) {
  ps <- (seq_len(J)-.5)/J; y <- pop$pdf(pop$ppf(ps)); y <- y/mean(y)
  losses <- rowSums((grid_pdq(J)-matrix(y, nrow(SHAPE_GRID), J, byrow = TRUE))^2); best <- NULL
  for (ix in order(losses)[seq_len(10L)]) {
    fit <- nelder_mead(function(sh) pdq_obj(sh, y, ps, QNODE, QWEIGHT), SHAPE_GRID[ix, ],
      maxiter = 2000L, xatol = 1e-10, fatol = 1e-12)
    if (is.null(best) || fit$fun < best$fun) best <- fit
  }
  qq <- pop$ppf(c(.25, .5, .75)); gq <- qvec(c(.25, .5, .75), c(0, 1, best$x))
  inv <- diff(gq[c(1,3)])/diff(qq[c(1,3)]); par <- c(qq[2]-gq[2]/inv, inv, best$x)
  t <- pop$true; fd <- cdf_pdf(c(t$M-t$S, t$M+t$S, t$M), par)
  v <- variance_map(fd$pdf[1], fd$pdf[2], fd$pdf[3], 1-sum(fd$cdf[1:2])); vh <- variance_map(fd$pdf[1], fd$pdf[2], fd$pdf[3], t$q)
  list(par = par, info = list(population = pop$name, J = J, lambda1 = par[1], lambda2 = par[2],
    lambda3 = par[3], lambda4 = par[4], objective = best$fun, success = best$success,
    variance_ratio = v/t$v, hybrid_ratio = vh/t$v, limiting_coverage = 2*pnorm(MAD_Z*sqrt(v/t$v))-1,
    a_fit = fd$pdf[1], b_fit = fd$pdf[2], c_fit = fd$pdf[3], q_fit = 1-sum(fd$cdf[1:2])))
}
run_population <- function(out, pop_indices = 0:9, Js = c(25L, 50L)) {
  cfg <- mad_config(out, "population", list(pop_indices = pop_indices, Js = Js)); pops <- populations(); rows <- list()
  for (pi in pop_indices) for (J in Js) {
    value <- mad_checkpoint(file.path(out, "checkpoints", sprintf("population_%d_%d.rds", pi, J)), cfg,
      function() population_pdq(pops[[pi+1L]], J))
    rows[[length(rows)+1L]] <- value$info; cat("Population", pops[[pi+1L]]$name, "J", J, "completed\n")
  }
  ans <- mad_bind(rows); mad_csv(ans, file.path(out, "population_pdq_calibration.csv")); invisible(ans)
}

run_galaxy <- function(out, input = "computation/data/galaxies.csv") {
  xx <- read.csv(input)$velocity_km_s/1000; ms <- mad_nuisance(xx); n <- length(xx); rows <- list(); pars <- list()
  add <- function(name, v, info) {
    se <- sqrt(v/n); rad <- MAD_Z*se/ms$s
    rows[[length(rows)+1L]] <<- c(list(method = name, mhat = ms$m, shat = ms$s, vhat = v, se = se,
      lower = ms$s-MAD_Z*se, upper = ms$s+MAD_Z*se, loglower = ms$s*exp(-rad), logupper = ms$s*exp(rad)), info)
  }
  for (spec in list(c("KDE-R", "robust", 1), c("KDE-NR", "nrd0", 1), c("KDE-half", "robust", .5), c("KDE-double", "robust", 2))) {
    v <- kernel_estimate(xx, spec[2], as.numeric(spec[3])); add(spec[1], v$value, v$info)
  }
  pd <- fit_pdq(xx); tm <- fit_tm(xx, pd$par)
  for (name in c("GLD-pdQ", "GLD-TM")) {
    fit <- if (name == "GLD-pdQ") pd else tm; v <- gld_estimate(xx, fit$par)
    add(name, v$value, c(v$info, list(fit_converged = fit$info$success))); pars[[name]] <- fit$par
  }
  ans <- mad_bind(rows); mad_csv(ans, file.path(out, "galaxies_analysis.csv"))
  saveRDS(pars, file.path(out, "galaxies_gld_parameters.rds")); dput(pars, file = file.path(out, "galaxies_gld_parameters.R"))
  invisible(list(analysis = ans, parameters = pars))
}

run_variance_clt <- function(out, replications = 3000L, ns = c(200L, 1000L, 5000L), pop_indices = c(0L, 2L)) {
  cfg <- mad_config(out, "variance_clt", list(replications = replications, ns = ns, pop_indices = pop_indices)); pops <- populations(); rows <- list()
  for (pi in pop_indices) for (n in ns) {
    pop <- pops[[pi+1L]]; t <- pop$true; points <- c(t$M-t$S, t$M+t$S, t$M); density <- c(t$a, t$b, t$c)
    g <- variance_gradient(t$a, t$b, t$c, t$q)[1:3]
    f2 <- if (pi == 0L) (points^2-1)*density else { ll <- log(points); density*((1+ll)^2+ll)/points^2 }
    bias <- .5*sum(g*f2); h <- 1.482602218505602*t$S*n^(-.2)
    # pi is the population key, so refer to the mathematical constant explicitly.
    avar <- sum(g*g*density)/(2*sqrt(base::pi))
    vals <- vector("list", replications)
    for (rep in seq_len(replications)-1L) vals[[rep+1L]] <- mad_checkpoint(
      file.path(out, "checkpoints", sprintf("clt_%d_%d_%04d.rds", pi, n, rep)), cfg, function() {
        x <- mad_sample(pop, 3, pi, n, rep); ms <- mad_nuisance(x); ds <- gaussian_kde_at(x, ms$points, h)
        q <- 1-mean(x <= ms$points[1])-mean(x <= ms$points[2]); v <- variance_map(ds[1], ds[2], ds[3], q)
        list(rep = rep, population = pop$name, n = n, h = h, vhat = v,
          z = sqrt(n*h)*(v-t$v-h*h*bias)/sqrt(avar))
      })
    raw <- mad_bind(vals); mad_csv(raw, file.path(out, sprintf("variance_clt_raw_%d_%d.csv.gz", pi, n)))
    zs <- raw$z; qs <- quantile(zs, c(.025, .975), names = FALSE, type = 7)
    rows[[length(rows)+1]] <- list(population = pop$name, n = n, R = replications, h = h,
      bias_coefficient = bias, variance_coefficient = avar, z_mean = mean(zs), z_sd = sd(zs), z_q025 = qs[1], z_q975 = qs[2],
      central95 = mean(abs(zs) <= MAD_Z), mse = mean((raw$vhat-t$v)^2), amse = h^4*bias^2+avar/(n*h))
    mad_csv(mad_bind(rows), file.path(out, "variance_clt_check.csv")); cat("CLT", pop$name, n, "completed\n")
  }
  invisible(mad_bind(rows))
}

mad_long_summary <- function(df, out, stem) {
  ss <- mad_bind(lapply(split(df, interaction(df$population, df$n, df$method, drop = TRUE)), function(g) {
    cp <- mean(g$cover); list(population = g$population[1], n = g$n[1], method = g$method[1], R = nrow(g),
      coverage = cp, completed = sum(g$success), median_vratio = mad_median(g$variance_ratio),
      mean_vratio = mad_mean(g$variance_ratio), mcse = sqrt(cp*(1-cp)/nrow(g)))
  }))
  mad_csv(df, file.path(out, paste0(stem, "_raw.csv.gz"))); mad_csv(ss, file.path(out, paste0(stem, "_summary.csv"))); invisible(ss)
}
run_long_samples <- function(out, resolution = FALSE, replications = 1000L, pop_indices = NULL) {
  if (is.null(pop_indices)) pop_indices <- if (resolution) c(7L, 8L) else c(0L, 2L, 6L, 7L, 8L)
  n <- if (resolution) 50000L else 5000L; key <- if (resolution) 5L else 4L; stem <- if (resolution) "resolution" else "large_n"
  cfg <- mad_config(out, stem, list(replications = replications, pop_indices = pop_indices, n = n)); pops <- populations(); rows <- list()
  for (pi in pop_indices) {
    pop <- pops[[pi+1L]]
    for (rep in seq_len(replications)-1L) rows[[length(rows)+1L]] <- mad_checkpoint(
      file.path(out, "checkpoints", sprintf("%s_%d_%04d.rds", stem, pi, rep)), cfg, function() {
        x <- mad_sample(pop, key, pi, n, rep); fit <- fit_pdq(x); vg <- gld_estimate(x, fit$par)
        a <- kernel_estimate(x); ans <- list(mad_record("KDE-R", a$value, x, pop$true, rep, info = a$info),
          mad_record("GLD-pdQ", vg$value, x, pop$true, rep, info = modifyList(vg$info, fit$info)),
          mad_record("Oracle", pop$true$v, x, pop$true, rep))
        if (resolution) { b <- kernel_estimate(x, multiplier = .5); ans[[length(ans)+1L]] <- mad_record("KDE-half", b$value, x, pop$true, rep, info = b$info) }
        ans <- mad_bind(ans)
        # The original long-sample drivers retain a finite fitted variance ratio
        # even when optimization reports nonconvergence. Coverage still counts
        # that fit as failed. Preserve that diagnostic contract here.
        ans$variance_ratio <- ans$vhat/pop$true$v
        ans$population <- pop$name; ans
      })
    mad_long_summary(mad_bind(rows), out, stem); cat(stem, pop$name, "completed\n")
  }
  invisible(mad_bind(rows))
}
