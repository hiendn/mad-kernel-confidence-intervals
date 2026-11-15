# Nonadaptive Nelder-Mead in base R, with separate absolute parameter and
# objective stopping tolerances. Objective ties retain the vertex order.
# Objectives enforce their own domains; vertices are not clipped to bounds.
nelder_mead <- function(fn, x0, args = list(), maxiter = 700,
                        xatol = 1e-7, fatol = 1e-8, maxfev = Inf,
                        initial_simplex = NULL) {
  x0 <- as.numeric(x0)
  n <- length(x0)
  if (!n || any(!is.finite(x0))) stop("x0 must be a finite nonempty vector")
  if (!is.function(fn) || !is.list(args)) stop("fn must be a function and args a list")
  for (setting in list(maxiter, maxfev, xatol, fatol)) {
    if (length(setting) != 1L || is.na(setting) || setting < 0)
      stop("limits and tolerances must be nonnegative scalars")
  }
  if (is.null(initial_simplex)) {
    vertices <- matrix(x0, nrow = n + 1L, ncol = n, byrow = TRUE)
    for (j in seq_len(n))
      vertices[j + 1L, j] <- if (x0[j] == 0) 0.00025 else 1.05 * x0[j]
  } else {
    vertices <- as.matrix(initial_simplex)
    storage.mode(vertices) <- "double"
    if (!identical(dim(vertices), c(n + 1L, n)) || any(!is.finite(vertices)))
      stop("initial_simplex must be a finite (length(x0)+1) by length(x0) matrix")
  }
  evaluations <- 0L
  limited <- structure(list(message = "function evaluation limit"),
                       class = c("nm_evaluation_limit", "error", "condition"))
  evaluate <- function(x) {
    if (evaluations >= maxfev) stop(limited)
    evaluations <<- evaluations + 1L
    value <- do.call(fn, c(list(as.numeric(x)), args))
    if (length(value) != 1L || !is.numeric(value))
      stop("objective must return one numeric value")
    as.numeric(value)
  }
  values <- rep(Inf, n + 1L)
  tryCatch({
    for (j in seq_len(n + 1L)) values[j] <- evaluate(vertices[j, ])
  }, nm_evaluation_limit = function(e) NULL)
  rank_vertices <- function() {
    # Radix ordering is stable, including the penalties used by GLD fits.
    index <- order(values, na.last = TRUE, method = "radix")
    vertices <<- vertices[index, , drop = FALSE]
    values <<- values[index]
  }
  rank_vertices()
  iterations <- 1L
  while (evaluations < maxfev && iterations < maxiter) {
    spread_x <- max(abs(sweep(vertices[-1L, , drop = FALSE], 2L,
                              vertices[1L, ], "-")))
    spread_f <- max(abs(values[-1L] - values[1L]))
    if (isTRUE(spread_x <= xatol && spread_f <= fatol)) break
    tryCatch({
      centroid <- colSums(vertices[seq_len(n), , drop = FALSE]) / n
      worst <- vertices[n + 1L, ]
      reflected <- 2 * centroid - worst
      reflected_value <- evaluate(reflected)
      shrink <- FALSE
      if (isTRUE(reflected_value < values[1L])) {
        expanded <- 3 * centroid - 2 * worst
        expanded_value <- evaluate(expanded)
        if (isTRUE(expanded_value < reflected_value)) {
          vertices[n + 1L, ] <- expanded
          values[n + 1L] <- expanded_value
        } else {
          vertices[n + 1L, ] <- reflected
          values[n + 1L] <- reflected_value
        }
      } else if (isTRUE(reflected_value < values[n])) {
        vertices[n + 1L, ] <- reflected
        values[n + 1L] <- reflected_value
      } else if (isTRUE(reflected_value < values[n + 1L])) {
        contracted <- 1.5 * centroid - 0.5 * worst
        contracted_value <- evaluate(contracted)
        if (isTRUE(contracted_value <= reflected_value)) {
          vertices[n + 1L, ] <- contracted
          values[n + 1L] <- contracted_value
        } else shrink <- TRUE
      } else {
        contracted <- 0.5 * centroid + 0.5 * worst
        contracted_value <- evaluate(contracted)
        if (isTRUE(contracted_value < values[n + 1L])) {
          vertices[n + 1L, ] <- contracted
          values[n + 1L] <- contracted_value
        } else shrink <- TRUE
      }
      if (shrink) {
        for (j in 2L:(n + 1L)) {
          vertices[j, ] <- vertices[1L, ] + 0.5 * (vertices[j, ] - vertices[1L, ])
          values[j] <- evaluate(vertices[j, ])
        }
      }
      iterations <- iterations + 1L
    }, nm_evaluation_limit = function(e) NULL)
    rank_vertices()
  }
  status <- if (evaluations >= maxfev) 1L else if (iterations >= maxiter) 2L else 0L
  message <- switch(as.character(status),
                    "0" = "Optimization terminated successfully.",
                    "1" = "Maximum number of function evaluations reached.",
                    "2" = "Maximum number of iterations reached.")
  list(x = as.numeric(vertices[1L, ]), fun = values[1L],
       nit = iterations, nfev = evaluations, success = status == 0L,
       status = status, message = message,
       final_simplex = list(x = vertices, fun = values))
}
