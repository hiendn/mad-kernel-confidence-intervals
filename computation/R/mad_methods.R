# Direct MAD variance and FKML reference procedures for MAD-Kernel-2026-10-25.
# All sample quantiles use type 7; all MADs are unnormalised (constant = 1).
# Numerical fitting follows the frozen reference protocol, not gld package defaults.

SQRT2PI <- sqrt(2 * pi)
GRID <- c(-0.9, -0.5, -0.1, 0, 0.1, 0.2, 0.4, 0.8, 1, 1.5)
SHAPE_GRID <- cbind(rep(GRID, each = length(GRID)), rep(GRID, length(GRID)))

gauss_legendre <- function(n = 96L) {
  if (n < 1L || n != as.integer(n)) stop("n must be a positive integer")
  if (n == 1L) return(list(nodes = 0.5, weights = 1))
  j <- seq_len(n - 1L)
  b <- j / sqrt(4 * j^2 - 1)
  jacobi <- matrix(0, n, n)
  jacobi[cbind(j, j + 1L)] <- b
  jacobi[cbind(j + 1L, j)] <- b
  eig <- eigen(jacobi, symmetric = TRUE)
  ix <- order(eig$values)
  list(nodes = (eig$values[ix] + 1) / 2,
       weights = eig$vectors[1L, ix]^2)
}
.quadrature96 <- gauss_legendre(96L)
QNODE <- .quadrature96$nodes
QWEIGHT <- .quadrature96$weights

.powerpart <- function(u, a) {
  if (abs(a) < 1e-8) {
    l <- log(u)
    l + 0.5 * a * l^2 + (a^2 / 6) * l^3
  } else expm1(a * log(u)) / a
}
q_scalar <- function(p, loc, invscale, a, b) {
  loc + (.powerpart(p, a) - .powerpart(1 - p, b)) / invscale
}
qprime_scalar <- function(p, invscale, a, b) {
  (p^(a - 1) + (1 - p)^(b - 1)) / invscale
}
qvec <- function(p, par) {
  if (length(par) != 4L || any(!is.finite(par)) || par[2L] <= 0)
    stop("par must contain finite location, positive inverse scale and two shapes")
  if (any(!is.finite(p)) || any(p < 0 | p > 1)) stop("p must lie in [0, 1]")
  # Evaluate endpoints by their analytic limits, avoiding 0 * Inf in the series.
  ans <- rep(NA_real_, length(p))
  inside <- p > 0 & p < 1
  ans[inside] <- q_scalar(p[inside], par[1L], par[2L], par[3L], par[4L])
  ans[p == 0] <- if (par[3L] > 0) par[1L] - 1 / (par[2L] * par[3L]) else -Inf
  ans[p == 1] <- if (par[4L] > 0) par[1L] + 1 / (par[2L] * par[4L]) else Inf
  ans
}

cdf_pdf <- function(x, par) {
  if (length(par) != 4L || any(!is.finite(par)) || par[2L] <= 0)
    stop("invalid FKML parameters")
  loc <- par[1L]; invscale <- par[2L]; a <- par[3L]; b <- par[4L]
  lower <- if (a > 0) loc - 1 / (invscale * a) else -Inf
  upper <- if (b > 0) loc + 1 / (invscale * b) else Inf
  cdf <- pdf <- rep(NA_real_, length(x))
  below <- !is.na(x) & x <= lower
  above <- !is.na(x) & x >= upper
  cdf[below] <- 0; pdf[below] <- 0
  cdf[above] <- 1; pdf[above] <- 0
  # Finite support endpoints have density zero by convention.
  active <- which(!is.na(x) & !below & !above)
  if (length(active)) {
    xx <- x[active]
    lo <- rep(1e-14, length(xx)); hi <- rep(1 - 1e-14, length(xx))
    p <- rep(0.5, length(xx)); remaining <- seq_along(xx)
    for (k in seq_len(70L)) {
      if (!length(remaining)) break
      ii <- remaining
      err <- q_scalar(p[ii], loc, invscale, a, b) - xx[ii]
      done <- abs(err) < 2e-12 * (1 + abs(xx[ii]))
      ii <- ii[!done]; err <- err[!done]
      if (!length(ii)) break
      pos <- err > 0
      hi[ii[pos]] <- p[ii[pos]]; lo[ii[!pos]] <- p[ii[!pos]]
      newp <- p[ii] - err / qprime_scalar(p[ii], invscale, a, b)
      bad <- !is.finite(newp) | newp <= lo[ii] | newp >= hi[ii]
      newp[bad] <- (lo[ii[bad]] + hi[ii[bad]]) / 2
      done <- abs(newp - p[ii]) < 2e-14
      p[ii] <- newp
      remaining <- ii[!done]
    }
    cdf[active] <- p
    pdf[active] <- 1 / qprime_scalar(p, invscale, a, b)
  }
  list(cdf = cdf, pdf = pdf)
}

pdq_values <- function(a, b, ps, nodes = QNODE, weights = QWEIGHT) {
  integ <- sum(weights / (nodes^(a - 1) + (1 - nodes)^(b - 1)))
  1 / (ps^(a - 1) + (1 - ps)^(b - 1)) / integ
}
pdq_obj <- function(par, y, ps, nodes = QNODE, weights = QWEIGHT) {
  if (any(!is.finite(par)) || any(par < -2 | par > 5)) return(1e12)
  sum((pdq_values(par[1L], par[2L], ps, nodes, weights) - y)^2)
}
.qd_cache <- new.env(parent = emptyenv())
empirical_qd_weights <- function(n, J) {
  key <- paste(n, J, sep = ":")
  if (exists(key, .qd_cache, inherits = FALSE)) return(.qd_cache[[key]])
  ps <- (seq_len(J) - 0.5) / J
  z <- qnorm(ps); dens <- dnorm(z)
  qor <- dens^2 / (1 + 2 * z^2)
  bw <- pmin((15 / n)^0.2 * qor^0.4, 0.95 * pmin(ps, 1 - ps))
  lower <- (seq_len(n) - 1) / n; upper <- seq_len(n) / n
  u <- sweep(outer(ps, lower, "-"), 1L, bw, "/")
  v <- sweep(outer(ps, upper, "-"), 1L, bw, "/")
  kk <- function(z) 0.75 * pmax(1 - z^2, 0)
  W <- sweep(kk(u) - kk(v), 1L, bw, "/")
  ans <- list(ps = ps, W = W)
  # Bound retained matrices: long studies revisit the same few (n, J) designs.
  if (length(ls(.qd_cache)) >= 32L) rm(list = ls(.qd_cache)[1L], envir = .qd_cache)
  .qd_cache[[key]] <- ans
  ans
}
.grid_cache <- new.env(parent = emptyenv())
grid_pdq <- function(J) {
  key <- as.character(J)
  if (!exists(key, .grid_cache, inherits = FALSE)) {
    ps <- (seq_len(J) - 0.5) / J
    .grid_cache[[key]] <- t(vapply(seq_len(nrow(SHAPE_GRID)), function(i)
      pdq_values(SHAPE_GRID[i, 1L], SHAPE_GRID[i, 2L], ps), numeric(J)))
  }
  .grid_cache[[key]]
}
quartile_params <- function(x, shapes) {
  qq <- as.numeric(quantile(x, c(0.25, 0.5, 0.75), type = 7, names = FALSE))
  gq <- qvec(c(0.25, 0.5, 0.75), c(0, 1, shapes))
  inv <- (gq[3L] - gq[1L]) / (qq[3L] - qq[1L])
  c(qq[2L] - gq[2L] / inv, inv, shapes)
}
.check_sample <- function(x) {
  if (!is.numeric(x) || length(x) < 2L || any(!is.finite(x)))
    stop("x must contain at least two finite numeric observations")
  as.numeric(x)
}
.require_optimizer <- function() {
  if (!exists("nelder_mead", mode = "function"))
    stop("Source computation/R/nelder_mead.R before fitting GLDs")
}
fit_pdq <- function(x) {
  x <- .check_sample(x); xs <- sort(x); n <- length(xs)
  J <- if (n <= 200L) 25L else 50L
  pw <- empirical_qd_weights(n, J); qd <- as.vector(pw$W %*% xs)
  if (any(qd <= 0) || any(!is.finite(qd)))
    return(list(par = NULL, info = list(success = FALSE,
      reason = "nonpositive quantile-density estimate")))
  raw <- 1 / qd; y <- raw / mean(raw)
  losses <- rowSums(sweep(grid_pdq(J), 2L, y, "-")^2)
  start <- SHAPE_GRID[which.min(losses), ]
  .require_optimizer()
  res <- nelder_mead(pdq_obj, start, args = list(y, pw$ps, QNODE, QWEIGHT),
    maxiter = 700L, xatol = 1e-7, fatol = 1e-8)
  par <- quartile_params(xs, res$x)
  list(par = par, info = list(success = res$success, objective = res$fun,
    nit = res$nit, boundary = any(res$x < -1.99 | res$x > 4.99), J = J))
}
tm_obj <- function(par, z) {
  loc <- par[1L]; invscale <- par[2L]; a <- par[3L]; b <- par[4L]
  if (any(!is.finite(par)) || invscale <= 1e-7 || invscale > 1e7 ||
      a < -3 || a > 6 || b < -3 || b > 6 || abs(loc) > 1e5) return(1e12)
  p <- cdf_pdf(z, par)$cdf
  spacing <- diff(c(0, p, 1))
  if (any(!is.finite(spacing)) || any(spacing <= 1e-15)) return(1e12)
  -sum(log(spacing))
}
fit_tm <- function(x, pdq_par = NULL) {
  x <- .check_sample(x); med <- median(x)
  iqr <- diff(quantile(x, c(0.25, 0.75), type = 7, names = FALSE))
  if (iqr <= 0) return(list(par = NULL,
    info = list(success = FALSE, reason = "zero IQR")))
  xs <- sort((x - med) / iqr)
  z <- (xs[-length(xs)] + xs[-1L]) / 2
  candidates <- t(vapply(seq_len(nrow(SHAPE_GRID)), function(i)
    quartile_params(xs, SHAPE_GRID[i, ]), numeric(4L)))
  if (!is.null(pdq_par)) candidates <- rbind(candidates,
    c((pdq_par[1L] - med) / iqr, pdq_par[2L] * iqr, pdq_par[3:4]))
  vals <- apply(candidates, 1L, tm_obj, z = z)
  # Stable candidate ordering is deliberate; equal penalties retain grid order.
  ix <- order(vals, seq_along(vals))[1:2]
  .require_optimizer(); best <- NULL
  for (i in ix) {
    res <- nelder_mead(tm_obj, candidates[i, ], args = list(z),
      maxiter = 1600L, xatol = 2e-6, fatol = 1e-6)
    if (is.null(best) || res$fun < best$fun) best <- res
  }
  par <- c(med + iqr * best$x[1L], best$x[2L] / iqr, best$x[3:4])
  list(par = par, info = list(success = best$success && best$fun < 1e11,
    objective = best$fun, nit = best$nit,
    boundary = any(best$x[3:4] < -2.99 | best$x[3:4] > 5.99)))
}

variance_map <- function(a, b, c, q) {
  if (any(!is.finite(c(a, b, c, q))) || c <= 0 || a + b <= 0) return(NaN)
  lambda <- (a - b) / c
  H <- 1 + lambda^2 + 4 * lambda * q
  if (!is.finite(H) || H <= 0) return(NaN)
  value <- H / (4 * (a + b)^2)
  if (is.finite(value)) value else NaN
}
variance_gradient <- function(a, b, c, q) {
  A <- a + b; lambda <- (a - b) / c
  H <- 1 + lambda^2 + 4 * lambda * q; T <- 2 * lambda + 4 * q
  c(T / (4 * A^2 * c) - H / (2 * A^3),
    -T / (4 * A^2 * c) - H / (2 * A^3),
    -T * lambda / (4 * A^2 * c), lambda / A^2)
}
mad_nuisance <- function(x) {
  x <- .check_sample(x); m <- median(x); s <- median(abs(x - m))
  list(m = m, s = s, points = c(m - s, m + s, m))
}
gaussian_kde_at <- function(x, points, h) {
  if (length(h) != 1L || !is.finite(h) || h <= 0) stop("h must be finite and positive")
  vapply(points, function(t) mean(dnorm((t - x) / h)) / h, numeric(1L))
}
kernel_estimate <- function(x, bandwidth = "robust", multiplier = 1) {
  x <- .check_sample(x); ns <- mad_nuisance(x); n <- length(x)
  if (ns$s <= 0) return(list(value = NaN,
    info = list(success = FALSE, reason = "zero MAD")))
  if (length(multiplier) != 1L || !is.finite(multiplier) || multiplier <= 0)
    stop("multiplier must be finite and positive")
  if (bandwidth == "robust") {
    h <- 1.482602218505602 * ns$s * n^(-0.2) * multiplier
  } else if (bandwidth == "nrd0") {
    scale <- min(sd(x), diff(quantile(x, c(0.25, 0.75), type = 7, names = FALSE)) / 1.34)
    h <- 0.9 * scale * n^(-0.2) * multiplier
  } else stop("Unknown bandwidth")
  dens <- gaussian_kde_at(x, ns$points, h)
  a <- dens[1L]; b <- dens[2L]; c <- dens[3L]
  q <- 1 - mean(x <= ns$points[1L]) - mean(x <= ns$points[2L])
  v <- variance_map(a, b, c, q)
  floor <- 1e-10 / (n * ns$s); ct <- max(c, floor); At <- max(a + b, floor)
  phi <- (0.5 - (abs(x - ns$m) <= ns$s)) / At +
    (a - b) / (At * ct) * (0.5 - (x <= ns$m))
  eif <- mean((phi - mean(phi))^2)
  list(value = v, info = list(success = is.finite(v) && v > 0,
    h = h, a = a, b = b, c = c, q = q, eif = eif,
    floor = c < floor || a + b < floor))
}
gld_estimate <- function(x, par, empirical_cdf = FALSE) {
  ns <- mad_nuisance(x)
  if (is.null(par)) return(list(value = NaN, info = list()))
  fitted <- cdf_pdf(ns$points, par)
  q <- 1 - fitted$cdf[1L] - fitted$cdf[2L]
  if (empirical_cdf) q <- 1 - mean(x <= ns$points[1L]) - mean(x <= ns$points[2L])
  a <- fitted$pdf[1L]; b <- fitted$pdf[2L]; c <- fitted$pdf[3L]
  v <- variance_map(a, b, c, q)
  list(value = v, info = list(a = a, b = b, c = c, q = q,
    success = is.finite(v) && v > 0))
}

compute_truth <- function(pop) {
  m <- pop$ppf(0.5)
  hi <- max(pop$ppf(0.75) - m, m - pop$ppf(0.25)) * 1.01
  while (pop$cdf(m + hi) - pop$cdf(m - hi) < 0.5) hi <- hi * 2
  s <- uniroot(function(s) pop$cdf(m + s) - pop$cdf(m - s) - 0.5,
               c(0, hi), tol = 1e-13)$root
  d <- pop$pdf(c(m - s, m + s, m))
  q <- 1 - pop$cdf(m - s) - pop$cdf(m + s)
  list(M = m, S = s, a = d[1L], b = d[2L], c = d[3L], q = q,
       v = variance_map(d[1L], d[2L], d[3L], q))
}
new_population <- function(name, description, cdf, pdf, ppf, rvs) {
  pop <- list(name = name, description = description, cdf = cdf,
              pdf = pdf, ppf = ppf, rvs = rvs)
  pop$true <- compute_truth(pop)
  pop
}
mixture_population <- function(name, weights, means, sds) {
  stopifnot(length(weights) == length(means), length(means) == length(sds),
    all(weights >= 0), abs(sum(weights) - 1) < 1e-12, all(sds > 0))
  cdf <- function(x) vapply(x, function(t) sum(weights * pnorm((t - means) / sds)), numeric(1))
  pdf <- function(x) vapply(x, function(t) sum(weights * dnorm((t - means) / sds) / sds), numeric(1))
  ppf <- function(p) vapply(p, function(u) {
    if (u == 0) return(-Inf)
    if (u == 1) return(Inf)
    uniroot(function(t) cdf(t) - u,
      c(min(means - 15 * sds), max(means + 15 * sds)), tol = 1e-13)$root
  }, numeric(1))
  rvs <- function(n) {
    ix <- sample.int(length(weights), n, replace = TRUE, prob = weights)
    means[ix] + sds[ix] * rnorm(n)
  }
  new_population(name, paste0("normal mixture weights=", paste(weights, collapse = ","),
    "; means=", paste(means, collapse = ","), "; sds=", paste(sds, collapse = ",")),
    cdf, pdf, ppf, rvs)
}
gld_population <- function(name, par) {
  cdf <- function(x) cdf_pdf(x, par)$cdf
  pdf <- function(x) cdf_pdf(x, par)$pdf
  ppf <- function(p) qvec(p, par)
  new_population(name, paste("FKML GLD", paste(par, collapse = ",")),
    cdf, pdf, ppf, function(n) qvec(runif(n), par))
}
populations <- function() {
  list(
    new_population("Normal", "N(0,1)", pnorm, dnorm, qnorm, rnorm),
    new_population("Exponential", "Exp(rate=1)", pexp, dexp, qexp, rexp),
    new_population("Lognormal", "lognormal(0,1)", plnorm, dlnorm, qlnorm, rlnorm),
    new_population("Chi-square5", "chi-square(5)", function(x) pchisq(x, 5),
      function(x) dchisq(x, 5), function(p) qchisq(p, 5), function(n) rchisq(n, 5)),
    new_population("Pareto7", "Pareto(x_min=1, shape=7)",
      function(x) ifelse(x < 1, 0, -expm1(-7 * log(pmax(x, 1)))),
      function(x) ifelse(x < 1, 0, 7 * pmax(x, 1)^(-8)),
      function(p) (1 - p)^(-1/7), function(n) (1 - runif(n))^(-1/7)),
    new_population("Cauchy", "standard Cauchy", pcauchy, dcauchy, qcauchy, rcauchy),
    gld_population("GLD", c(0, 1, 0.2, 0.5)),
    mixture_population("Shoulder", c(0.75, 0.25), c(0, 2.5), c(1, 0.35)),
    mixture_population("Bimodal", c(0.5, 0.5), c(-1.25, 1.25), c(0.65, 0.65)),
    mixture_population("Contaminated", c(0.9, 0.1), c(0, 0), c(1, 5)))
}
