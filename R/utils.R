.as_data_matrix <- function(x, name = "x", min_rows = 1L) {
  if (is.data.frame(x)) {
    x <- data.matrix(x)
  }
  if (is.vector(x) && is.numeric(x)) {
    x <- matrix(x, nrow = 1L)
  }
  if (!is.matrix(x) || !is.numeric(x)) {
    stop(sprintf("`%s` must be a numeric matrix or data frame.", name),
         call. = FALSE)
  }
  storage.mode(x) <- "double"
  if (nrow(x) < min_rows || ncol(x) < 1L) {
    stop(sprintf("`%s` must have at least %d row(s) and one column.",
                 name, min_rows), call. = FALSE)
  }
  if (anyNA(x) || any(!is.finite(x))) {
    stop(sprintf("`%s` must contain only finite values.", name),
         call. = FALSE)
  }
  x
}

.as_location <- function(center, p, name = "center") {
  center_names <- names(center)
  center <- as.numeric(center)
  if (length(center) != p || anyNA(center) || any(!is.finite(center))) {
    stop(sprintf("`%s` must be a finite numeric vector of length %d.",
                 name, p), call. = FALSE)
  }
  names(center) <- center_names
  center
}

.weighted_median <- function(x, weights) {
  ordering <- order(x)
  x <- x[ordering]
  weights <- weights[ordering]
  x[which(cumsum(weights) >= sum(weights) / 2)[1L]]
}

.validate_iteration_controls <- function(tol, max_iter, zero_tol) {
  tol <- as.numeric(tol)
  if (length(tol) != 1L || is.na(tol) || !is.finite(tol) || tol <= 0) {
    stop("`tol` must be a finite positive number.", call. = FALSE)
  }

  if (!is.numeric(max_iter) || length(max_iter) != 1L ||
      is.na(max_iter) || !is.finite(max_iter) || max_iter < 1 ||
      max_iter != floor(max_iter) || max_iter > .Machine$integer.max) {
    stop("`max_iter` must be a positive integer.", call. = FALSE)
  }

  zero_tol <- as.numeric(zero_tol)
  if (length(zero_tol) != 1L || is.na(zero_tol) ||
      !is.finite(zero_tol) || zero_tol < 0) {
    stop("`zero_tol` must be a finite non-negative number.", call. = FALSE)
  }

  list(tol = tol, max_iter = as.integer(max_iter), zero_tol = zero_tol)
}

.validate_zero_tol <- function(zero_tol) {
  zero_tol <- as.numeric(zero_tol)
  if (length(zero_tol) != 1L || is.na(zero_tol) ||
      !is.finite(zero_tol) || zero_tol < 0) {
    stop("`zero_tol` must be a finite non-negative number.", call. = FALSE)
  }
  zero_tol
}

.center_and_scale <- function(x, center, zero_tol = 0) {
  zero_tol <- .validate_zero_tol(zero_tol)

  centered <- suppressWarnings(
    sweep(x, 2L, center, check.margin = FALSE)
  )
  overflow_fallback <- any(!is.finite(centered))
  if (overflow_fallback) {
    scale <- max(abs(x), abs(center))
    if (!is.finite(scale) || scale <= 0) {
      scale <- 1
    }
    centered <- sweep(
      x / scale, 2L, center / scale, check.margin = FALSE
    )
  } else {
    scale <- max(abs(centered))
    if (!is.finite(scale) || scale <= 0) {
      scale <- 1
    } else {
      centered <- centered / scale
    }
  }

  scaled_zero_tol <- zero_tol / scale
  if (!is.finite(scaled_zero_tol)) {
    scaled_zero_tol <- .Machine$double.xmax
  }
  list(
    x = centered,
    scale = scale,
    zero_tol = scaled_zero_tol,
    overflow_fallback = overflow_fallback
  )
}

.has_zero_residual <- function(x, center, zero_tol = 0) {
  scaled <- .center_and_scale(x, center, zero_tol)
  cpp_spatial_sign(scaled$x, scaled$zero_tol)$n_zero > 0
}

.hr_default_initial <- function(x, tol, max_iter, zero_tol, warn) {
  spatial <- as.numeric(spatial_median(
    x, tol = tol, max_iter = max_iter, zero_tol = zero_tol, warn = warn
  ))
  if (!.has_zero_residual(x, spatial, zero_tol)) {
    return(spatial)
  }

  sample_mean <- colMeans(x)
  if (!.has_zero_residual(x, sample_mean, zero_tol)) {
    return(sample_mean)
  }

  if (nrow(x) > 2L) {
    for (i in seq_len(nrow(x))) {
      candidate <- colMeans(x[-i, , drop = FALSE])
      if (!.has_zero_residual(x, candidate, zero_tol)) {
        return(candidate)
      }
    }
  }

  for (fraction in 1 / seq.int(2, 33)) {
    for (i in seq_len(nrow(x))) {
      candidate <- (1 - fraction) * sample_mean + fraction * x[i, ]
      if (!.has_zero_residual(x, candidate, zero_tol)) {
        return(candidate)
      }
    }
  }
  stop(
    "The HR equations need a noncoincident starting location, but none " %+%
      "could be constructed from these data.",
    call. = FALSE
  )
}

.resolve_center <- function(x, center, tol, max_iter, zero_tol) {
  p <- ncol(x)
  if (is.numeric(center)) {
    return(.as_location(center, p))
  }
  center <- match.arg(center, c("spatial", "mean", "none"))
  switch(
    center,
    spatial = as.numeric(spatial_median(
      x, tol = tol, max_iter = max_iter, zero_tol = zero_tol
    )),
    mean = colMeans(x),
    none = rep.int(0, p)
  )
}

.shape_dimnames <- function(estimate, x) {
  variable_names <- colnames(x)
  if (!is.null(variable_names)) {
    dimnames(estimate) <- list(variable_names, variable_names)
  }
  estimate
}
