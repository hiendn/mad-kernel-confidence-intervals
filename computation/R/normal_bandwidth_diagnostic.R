# Deterministic first-order Gaussian finite-bandwidth diagnostic.
# This is neither a finite-sample variance nor a coverage correction.
normal_bandwidth_diagnostic <- function(output_dir = "computation/generated_r/diagnostics") {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  s <- qnorm(.75)
  f <- dnorm(s)
  g <- -1 / (2 * (2 * f)^3)
  k <- -2 * s * f * g
  af <- 2 * g^2 * f / (2 * sqrt(pi))
  segments <- list(c(-Inf, -s), c(-s, s), c(s, Inf))
  rows <- lapply(c(200L, 1000L, 5000L), function(n) {
    h <- n^(-.2)
    mu <- dnorm(s, sd = sqrt(1 + h^2))
    w <- function(x) g * (dnorm((s-x)/h)/h + dnorm((-s-x)/h)/h - 2*mu) +
      k * (.5 - as.numeric(abs(x) <= s)) / (2*f)
    moment <- function(power) sum(vapply(segments, function(ab) {
      integrate(function(x) w(x)^power * dnorm(x), ab[1], ab[2],
                abs.tol = 1e-10, rel.tol = 1e-10, subdivisions = 500L)$value
    }, numeric(1)))
    data.frame(n = n, E_W = moment(1), quadrature_SD = sqrt(h*moment(2)/af))
  })
  result <- do.call(rbind, rows)
  write.csv(result, file.path(output_dir, "finite_h_quadrature.csv"), row.names = FALSE)
  if (requireNamespace("jsonlite", quietly = TRUE)) {
    jsonlite::write_json(result, file.path(output_dir, "finite_h_quadrature.json"),
                         dataframe = "rows", pretty = TRUE, digits = 17)

  }
  writeLines(c("Computation: native R integrate(), split at both MAD endpoints.",
               "RNG: none; deterministic Gaussian integration.",
               "Scope: first-order linear diagnostic, not full estimator variance.",
               paste("R:", R.version.string)), file.path(output_dir, "PROVENANCE.txt"))
  result
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  output <- if (length(args)) args[1] else "computation/generated_r/diagnostics"
  print(normal_bandwidth_diagnostic(output))
}
