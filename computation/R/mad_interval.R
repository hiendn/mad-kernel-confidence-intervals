# Direct Gaussian-kernel intervals for the unnormalised sample MAD.
# The interval interface uses the raw
# algebraic variance, with no denominator floor or theoretical fallback.

mad_interval <- function(x, bandwidth = "KDE-R", confidence = 0.95,
                         report_scale = 1) {
  if (!is.atomic(x) || !is.null(dim(x)) ||
      !(is.numeric(x) || is.logical(x)) || length(x) < 2L) {
    stop("x must be a one-dimensional numeric sample with at least two observations",
         call. = FALSE)
  }
  x <- as.double(x)
  if (any(!is.finite(x))) stop("x must contain only finite observations", call. = FALSE)
  if (length(confidence) != 1L || !is.numeric(confidence) ||
      !is.finite(confidence) || confidence <= 0 || confidence >= 1) {
    stop("confidence must lie strictly between zero and one", call. = FALSE)
  }
  if (length(report_scale) != 1L || !is.numeric(report_scale) ||
      !is.finite(report_scale) || report_scale <= 0) {
    stop("report_scale must be finite and positive", call. = FALSE)
  }
  if (length(bandwidth) != 1L || is.na(bandwidth) ||
      !bandwidth %in% c("KDE-R", "KDE-half", "KDE-double", "KDE-NR")) {
    stop("Unrecognised reference bandwidth", call. = FALSE)
  }

  n <- length(x)
  m <- stats::median(x)
  s <- stats::median(abs(x - m))
  cuts <- c(m - s, m + s, m)
  result <- list(n = n, median = m, mad_raw = s,
                 mad_reported = report_scale * s, report_scale = report_scale,
                 confidence = confidence, bandwidth_rule = bandwidth,
                 success = FALSE, reason = NULL, h = NULL, cutpoints = cuts,
                 root_n_variance_raw = NULL, root_n_variance_reported = NULL,
                 sampling_variance = NULL, standard_error = NULL,
                 additive_interval = NULL, log_interval = NULL)

  finite_or_null <- function(value) {
    if (is.list(value)) return(lapply(value, finite_or_null))
    if (is.numeric(value) && any(!is.finite(value))) {
      if (length(value) == 1L) return(NULL)
      return(lapply(value, finite_or_null))
    }
    value
  }
  failure <- function(reason) {
    result["reason"] <- list(reason)
    finite_or_null(result)
  }

  if (!is.finite(m)) return(failure("nonfinite sample median"))
  if (!is.finite(s) || s <= 0) return(failure("zero or nonfinite raw sample MAD"))
  if (bandwidth == "KDE-NR") {
    quartiles <- stats::quantile(x, c(0.25, 0.75), type = 7, names = FALSE)
    h <- 0.9 * min(stats::sd(x), diff(quartiles) / 1.34) * n^(-0.2)
  } else {
    multiplier <- switch(bandwidth, "KDE-R" = 1, "KDE-half" = 0.5,
                         "KDE-double" = 2)
    h <- s * n^(-0.2) * multiplier / stats::qnorm(0.75)
  }
  if (!is.finite(h) || h <= 0) return(failure("nonpositive or nonfinite bandwidth"))
  result["h"] <- list(h)
  if (any(!is.finite(cuts))) return(failure("nonfinite empirical cutpoint"))

  # Direct sums preserve the estimator used in the paper; density() would
  # introduce a grid and interpolation.
  dens <- vapply(cuts, function(t) mean(stats::dnorm((t - x) / h)) / h,
                 numeric(1))
  a <- dens[1L]
  b <- dens[2L]
  c <- dens[3L]
  q <- 1 - mean(x <= cuts[1L]) - mean(x <= cuts[2L])
  result["local_estimates"] <- list(list(a = a, b = b, c = c, q = q))
  if (any(!is.finite(dens)) || c <= 0 || a + b <= 0) {
    return(failure("invalid local density estimate"))
  }
  lambda <- (a - b) / c
  density_sum_square <- (a + b)^2
  if (!is.finite(density_sum_square) || density_sum_square == 0) {
    return(failure("raw algebraic variance is not representable"))
  }
  v <- (1 + lambda^2 + 4 * lambda * q) / (4 * density_sum_square)
  if (!is.finite(v) || v <= 0) {
    return(failure("nonpositive or nonfinite raw algebraic variance"))
  }
  z <- stats::qnorm((1 + confidence) / 2)
  se <- sqrt(v / n)
  additive <- report_scale * c(s - z * se, s + z * se)
  logarithmic <- report_scale * s * exp(c(-z, z) * se / s)
  if (any(!is.finite(additive)) || any(!is.finite(logarithmic))) {
    return(failure("nonfinite computed interval endpoint"))
  }
  if (any(logarithmic <= 0)) return(failure("log-interval endpoint underflowed to zero"))
  scale_square <- report_scale^2
  if (!is.finite(scale_square)) {
    return(failure("reported estimate or variance is not representable"))
  }
  reported_variance <- scale_square * v
  sampling_variance <- scale_square * v / n
  reported_se <- report_scale * se
  if (any(!is.finite(c(result$mad_reported, reported_variance,
                       sampling_variance, reported_se)))) {
    return(failure("nonfinite reported estimate or variance"))
  }
  result$success <- TRUE
  result$root_n_variance_raw <- v
  result$root_n_variance_reported <- reported_variance
  result$sampling_variance <- sampling_variance
  result$standard_error <- reported_se
  result$additive_interval <- additive
  result$log_interval <- logarithmic
  result
}

mad_interval_json <- function(result) {
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop("Package 'jsonlite' is required only for JSON output", call. = FALSE)
  }
  paste0(jsonlite::toJSON(result, auto_unbox = TRUE, null = "null",
                         digits = 17, pretty = TRUE), "\n")
}

mad_interval_main <- function(args = commandArgs(trailingOnly = TRUE)) {
  usage <- paste(
    "Usage: Rscript computation/R/mad_interval.R --data FILE --column NAME",
    "       [--divide-by 1] [--bandwidth KDE-R] [--confidence 0.95]",
    "       [--report-scale 1] [--output FILE]", sep = "\n")
  if (any(args %in% c("--help", "-h"))) {
    cat(usage, "\n")
    return(invisible(NULL))
  }
  known <- c("--data", "--column", "--divide-by", "--bandwidth",
             "--confidence", "--report-scale", "--output")
  if (!length(args) || length(args) %% 2L ||
      any(!args[seq.int(1L, length(args), by = 2L)] %in% known)) {
    stop(usage, call. = FALSE)
  }
  opts <- list("--divide-by" = "1", "--bandwidth" = "KDE-R",
               "--confidence" = "0.95", "--report-scale" = "1")
  if (length(args)) {
    for (i in seq.int(1L, length(args), by = 2L)) opts[args[i]] <- list(args[i + 1L])
  }
  if (is.null(opts[["--data"]]) || is.null(opts[["--column"]])) {
    stop(usage, call. = FALSE)
  }
  # R 4.4's character-to-double conversion can round the decimal
  # 0.9999999999999999 up to 1. jsonlite preserves this representable
  # probability, which matters for the explicit endpoint-overflow check.
  parse_number <- function(value) {
    if (!requireNamespace("jsonlite", quietly = TRUE)) {
      stop("Package 'jsonlite' is required for command-line JSON input and output", call. = FALSE)
    }
    number <- jsonlite::fromJSON(value)
    if (!is.numeric(number) || length(number) != 1L) {
      stop("Numeric options must contain one number", call. = FALSE)
    }
    as.double(number)
  }
  divide_by <- parse_number(opts[["--divide-by"]])
  if (!is.finite(divide_by) || divide_by <= 0) {
    stop("--divide-by must be finite and positive", call. = FALSE)
  }
  dat <- utils::read.csv(opts[["--data"]], check.names = FALSE,
                         stringsAsFactors = FALSE)
  column <- opts[["--column"]]
  if (!column %in% names(dat)) stop(sprintf("Column '%s' not found", column), call. = FALSE)
  values <- suppressWarnings(as.double(dat[[column]])) / divide_by
  result <- mad_interval(values, bandwidth = opts[["--bandwidth"]],
                         confidence = parse_number(opts[["--confidence"]]),
                         report_scale = parse_number(opts[["--report-scale"]]))
  output <- mad_interval_json(result)
  if (is.null(opts[["--output"]])) cat(output) else {
    dir.create(dirname(opts[["--output"]]), recursive = TRUE, showWarnings = FALSE)
    writeLines(output, opts[["--output"]], useBytes = TRUE)
  }
  invisible(result)
}

if (sys.nframe() == 0L) {
  tryCatch(mad_interval_main(), error = function(e) {
    message(conditionMessage(e))
    quit(status = 1L)
  })
}
