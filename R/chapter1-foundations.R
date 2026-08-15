#' Spatial signs
#'
#' Computes the spatial sign map \eqn{U(x)=x/\lVert x\rVert_2}, with the zero
#' vector mapped to zero. For a matrix, rows are observations and columns are
#' variables.
#'
#' @param x A numeric vector, matrix, or data frame.
#' @param center Optional center to subtract from every observation.
#' @param zero_tol Non-negative tolerance below which a norm is treated as zero.
#'   The default detects exact zeros only and therefore preserves scale
#'   equivariance.
#'
#' @return A vector when `x` is a vector, otherwise a matrix. The result carries
#'   `norms` and `n_zero` attributes.
#' @references
#' Oja, H. (2010). *Multivariate Nonparametric Methods with R: An Approach
#' Based on Spatial Signs and Ranks*. Springer.
#' \doi{10.1007/978-1-4419-0468-3}.
#' @export
#'
#' @examples
#' spatial_sign(c(3, 4))
#' spatial_sign(rbind(c(3, 4), c(0, 0)))
spatial_sign <- function(x, center = NULL, zero_tol = 0) {
  vector_input <- is.null(dim(x))
  if (vector_input) {
    if (!is.numeric(x)) {
      stop("`x` must be numeric.", call. = FALSE)
    }
    x <- matrix(x, nrow = 1L)
  }
  x <- .as_data_matrix(x)
  norm_scale <- 1
  if (!is.null(center)) {
    center <- .as_location(center, ncol(x))
    scaled <- .center_and_scale(x, center, zero_tol)
    x <- scaled$x
    zero_tol <- scaled$zero_tol
    norm_scale <- scaled$scale
  }
  result <- cpp_spatial_sign(x, zero_tol)
  result$norms <- result$norms * norm_scale
  signs <- result$signs
  attr(signs, "norms") <- as.numeric(result$norms)
  attr(signs, "n_zero") <- as.numeric(result$n_zero)
  if (vector_input) {
    signs <- as.numeric(signs)
    attr(signs, "norms") <- as.numeric(result$norms)
    attr(signs, "norm") <- as.numeric(result$norms)
    attr(signs, "n_zero") <- as.numeric(result$n_zero)
  }
  signs
}


#' Spatial median
#'
#' Minimizes the weighted average Euclidean distance to the observations. The
#' implementation uses a modified Weiszfeld iteration that remains valid when
#' an iterate coincides with an observation.
#'
#' @param x A numeric matrix or data frame with observations in rows.
#' @param weights Optional non-negative observation weights.
#' @param initial Optional starting location. Coordinatewise weighted medians
#'   are used by default.
#' @param tol Tolerance for the normalized subgradient-equation residual.
#' @param max_iter Maximum number of modified Weiszfeld iterations.
#' @param zero_tol Non-negative tolerance for coincident observations.
#' @param warn If `TRUE`, warn when the iteration does not converge.
#'
#' @return A numeric location vector. Convergence diagnostics are stored in the
#'   `objective`, `iterations`, `converged`, `relative_change`, and
#'   `equation_residual` attributes.
#' @export
#'
#' @references
#' Oja, H. (2010). *Multivariate Nonparametric Methods with R*. Springer.
#'
#' @examples
#' x <- rbind(c(0, 0), c(1, 0), c(0, 1), c(20, 20))
#' spatial_median(x)
spatial_median <- function(x, weights = NULL, initial = NULL,
                           tol = 1e-8, max_iter = 500L,
                           zero_tol = 0, warn = TRUE) {
  x <- .as_data_matrix(x)
  n <- nrow(x)
  p <- ncol(x)
  if (is.null(weights)) {
    weights <- rep.int(1, n)
  }
  weights <- as.numeric(weights)
  weight_sum <- sum(weights)
  if (length(weights) != n || anyNA(weights) || any(!is.finite(weights)) ||
      any(weights < 0) || !is.finite(weight_sum) || weight_sum <= 0) {
    stop("`weights` must contain one finite non-negative value per row and " %+%
           "have a positive sum.", call. = FALSE)
  }
  if (is.null(initial)) {
    initial <- vapply(
      seq_len(p), function(j) .weighted_median(x[, j], weights), numeric(1)
    )
  } else {
    initial <- .as_location(initial, p, "initial")
  }
  origin <- initial
  scaled <- .center_and_scale(x, origin, zero_tol)
  result <- cpp_spatial_median(
    scaled$x, weights, numeric(p), as.numeric(tol), as.integer(max_iter),
    scaled$zero_tol
  )
  if (isTRUE(warn) && !isTRUE(result$converged)) {
    warning("The spatial median iteration did not converge in `max_iter` steps.",
            call. = FALSE)
  }
  displacement <- as.numeric(result$location)
  if (scaled$overflow_fallback) {
    estimate <- (origin / scaled$scale + displacement) * scaled$scale
  } else {
    estimate <- origin + displacement * scaled$scale
  }
  result$objective <- result$objective * scaled$scale
  names(estimate) <- colnames(x)
  for (field in c("objective", "iterations", "converged",
                  "relative_change", "equation_residual")) {
    attr(estimate, field) <- result[[field]]
  }
  estimate
}


#' Empirical spatial ranks
#'
#' For every row \eqn{x_i}, computes
#' \deqn{m^{-1}\sum_{j=1}^m U(x_i-y_j)}
#' with respect to the rows \eqn{y_j} of a reference sample. When `reference` is
#' omitted, the reference is `x`, including the zero self-difference exactly as
#' in the book definition.
#'
#' @param x Evaluation observations in rows.
#' @param reference Optional reference observations in rows.
#' @param zero_tol Non-negative tolerance below which differences are zero.
#'
#' @return A matrix of empirical spatial ranks.
#' @references
#' Oja, H. (2010). *Multivariate Nonparametric Methods with R: An Approach
#' Based on Spatial Signs and Ranks*. Springer.
#' \doi{10.1007/978-1-4419-0468-3}.
#' @export
#'
#' @examples
#' spatial_rank(rbind(c(0, 0), c(1, 0), c(0, 1)))
spatial_rank <- function(x, reference = NULL, zero_tol = 0) {
  x <- .as_data_matrix(x)
  if (is.null(reference)) {
    reference <- x
  } else {
    reference <- .as_data_matrix(reference, "reference")
  }
  zero_tol <- .validate_zero_tol(zero_tol)
  ranks <- cpp_spatial_rank(x, reference, zero_tol)
  dimnames(ranks) <- list(rownames(x), colnames(x))
  ranks
}


#' Spatial sign covariance matrix
#'
#' Computes the average outer product of centered spatial signs. By default the
#' sample is centered at its spatial median, matching Chapter 1 of the book.
#'
#' @param x Observations in rows.
#' @param center Either a numeric center or one of `"spatial"`, `"mean"`, and
#'   `"none"`.
#' @param tol,max_iter Convergence controls used when estimating a spatial
#'   center.
#' @param zero_tol Tolerance for zero residuals.
#'
#' @return A positive semidefinite matrix with `center` and `n_zero` attributes.
#'   Its trace is one when there are no zero residuals.
#' @export
#'
#' @references
#' Visuri, S., Koivunen, V., and Oja, H. (2000). Sign and rank covariance
#' matrices. *Journal of Statistical Planning and Inference*, 91, 557-575.
#'
#' @examples
#' set.seed(1)
#' sscm(matrix(rnorm(60), 20, 3))
sscm <- function(x, center = "spatial", tol = 1e-8, max_iter = 500L,
                 zero_tol = 0) {
  x <- .as_data_matrix(x)
  center_value <- .resolve_center(x, center, tol, max_iter, zero_tol)
  scaled <- .center_and_scale(x, center_value, zero_tol)
  sign_result <- cpp_spatial_sign(scaled$x, scaled$zero_tol)
  estimate <- crossprod(sign_result$signs) / nrow(x)
  estimate <- .shape_dimnames(estimate, x)
  attr(estimate, "center") <- center_value
  attr(estimate, "n_zero") <- as.numeric(sign_result$n_zero)
  estimate
}


#' Multivariate spatial Kendall matrix
#'
#' Computes the exact U-statistic average of outer products of spatial signs of
#' all pairwise differences. Pairwise differencing makes the estimator
#' translation invariant and removes the need to estimate location.
#'
#' @param x Observations in rows.
#' @param zero_tol Tolerance for tied pairwise differences.
#'
#' @return A positive semidefinite matrix with `n_pairs` and `n_zero_pairs`
#'   attributes. Its trace is `1 - n_zero_pairs / n_pairs`; in particular, it
#'   has unit trace when there are no tied pairs.
#' @export
#'
#' @references
#' Han, F. and Liu, H. (2018). ECA: High-dimensional elliptical component
#' analysis in non-Gaussian distributions. *Journal of the American Statistical
#' Association*, 113, 252-268.
#'
#' @examples
#' set.seed(2)
#' spatial_kendall(matrix(rnorm(40), 10, 4))
spatial_kendall <- function(x, zero_tol = 0) {
  x <- .as_data_matrix(x, min_rows = 2L)
  zero_tol <- .validate_zero_tol(zero_tol)
  result <- cpp_spatial_kendall(x, zero_tol)
  estimate <- .shape_dimnames(result$matrix, x)
  attr(estimate, "n_pairs") <- as.numeric(result$n_pairs)
  attr(estimate, "n_zero_pairs") <- as.numeric(result$n_zero_pairs)
  estimate
}


#' Spatial rank covariance matrix
#'
#' Computes the average outer product of the empirical spatial ranks defined
#' with denominator `n`, including zero self-differences.
#'
#' @inheritParams spatial_rank
#'
#' @return A positive semidefinite matrix.
#' @references
#' Oja, H. (2010). *Multivariate Nonparametric Methods with R: An Approach
#' Based on Spatial Signs and Ranks*. Springer.
#' \doi{10.1007/978-1-4419-0468-3}.
#' @export
#'
#' @examples
#' set.seed(3)
#' spatial_rank_covariance(matrix(rnorm(30), 10, 3))
spatial_rank_covariance <- function(x, zero_tol = 0) {
  x <- .as_data_matrix(x)
  ranks <- spatial_rank(x, zero_tol = zero_tol)
  estimate <- crossprod(ranks) / nrow(x)
  .shape_dimnames(estimate, x)
}


#' Normalize a shape matrix
#'
#' Removes the unidentified scalar from a positive definite shape matrix.
#'
#' @param shape A finite square numeric matrix.
#' @param method `"trace"` gives trace equal to the dimension; `"determinant"`
#'   gives determinant one.
#'
#' @return The symmetrized, normalized matrix.
#' @references
#' Fang, K.-T. and Anderson, T. W. (1990). *Statistical Inference in
#' Elliptically Contoured and Related Distributions*. Allerton Press.
#' @export
#'
#' @examples
#' shape <- matrix(c(4, 1, 1, 1), 2, 2)
#' normalize_shape(shape, method = "determinant")
normalize_shape <- function(shape, method = c("trace", "determinant")) {
  method <- match.arg(method)
  shape <- as.matrix(shape)
  storage.mode(shape) <- "double"
  if (!is.numeric(shape) || nrow(shape) != ncol(shape) || nrow(shape) < 1L ||
      anyNA(shape) || any(!is.finite(shape))) {
    stop("`shape` must be a finite square numeric matrix.", call. = FALSE)
  }
  matrix_scale <- max(abs(shape))
  if (matrix_scale <= 0) {
    stop("`shape` must be positive definite.", call. = FALSE)
  }
  scaled_shape <- shape / matrix_scale
  shape <- (scaled_shape + t(scaled_shape)) / 2
  eigenvalues <- eigen(shape, symmetric = TRUE, only.values = TRUE)$values
  if (min(eigenvalues) <= 0) {
    stop("`shape` must be positive definite.", call. = FALSE)
  }
  if (method == "trace") {
    shape * nrow(shape) / sum(diag(shape))
  } else {
    log_determinant <- as.numeric(determinant(shape, logarithm = TRUE)$modulus)
    shape / exp(log_determinant / nrow(shape))
  }
}


#' Tyler's shape estimator
#'
#' Solves Tyler's fixed-point equation and applies the book's trace
#' normalization, \eqn{\mathrm{tr}(V)=p}. The exact unregularized estimator
#' requires more observations than variables and data in general position.
#'
#' @param x Observations in rows.
#' @param center Either a numeric center or one of `"mean"`, `"spatial"`, and
#'   `"none"`. The default sample mean preserves affine equivariance of the
#'   resulting Tyler shape. A spatial-median center is more resistant to
#'   magnitude outliers but is only orthogonally equivariant.
#' @param initial Optional positive definite starting shape.
#' @param tol Relative Frobenius convergence tolerance.
#' @param max_iter Maximum number of fixed-point iterations.
#' @param zero_tol Tolerance used to detect undefined zero residuals. The
#'   default detects exact zeros only.
#' @param warn If `TRUE`, warn when the iteration does not converge.
#'
#' @return A trace-normalized shape matrix. Convergence diagnostics and the
#'   center are stored as attributes.
#' @export
#'
#' @references
#' Tyler, D. E. (1987). A distribution-free M-estimator of multivariate
#' scatter. *The Annals of Statistics*, 15, 234-251.
#'
#' @examples
#' set.seed(4)
#' x <- matrix(rt(200, df = 3), 50, 4)
#' tyler_shape(x)
tyler_shape <- function(x, center = "mean", initial = NULL,
                        tol = 1e-8, max_iter = 500L,
                        zero_tol = 0, warn = TRUE) {
  x <- .as_data_matrix(x)
  n <- nrow(x)
  p <- ncol(x)
  if (n <= p) {
    stop("The unregularized Tyler estimator requires `n > p`.",
         call. = FALSE)
  }
  center_value <- .resolve_center(x, center, tol, max_iter, zero_tol)
  scaled <- .center_and_scale(x, center_value, zero_tol)
  if (is.null(initial)) {
    initial <- diag(p)
  } else {
    initial <- normalize_shape(initial, "trace")
    if (!identical(dim(initial), c(p, p))) {
      stop("`initial` must match the data dimension.", call. = FALSE)
    }
  }
  result <- cpp_tyler_shape(
    scaled$x, initial, as.numeric(tol), as.integer(max_iter),
    scaled$zero_tol
  )
  if (isTRUE(warn) && !isTRUE(result$converged)) {
    warning("Tyler's fixed-point iteration did not converge in `max_iter` steps.",
            call. = FALSE)
  }
  estimate <- .shape_dimnames(result$shape, x)
  attr(estimate, "center") <- center_value
  for (field in c("iterations", "converged", "relative_change",
                  "equation_residual")) {
    attr(estimate, field) <- result[[field]]
  }
  estimate
}


#' Hettmansperger-Randles affine-equivariant location and shape
#'
#' Alternates the location and Tyler-type shape updates stated in Chapter 1.
#' The returned shape has trace equal to the data dimension.
#'
#' @param x Observations in rows.
#' @param initial_location Optional starting location. The default begins with
#'   the spatial median and deterministically falls back to a noncoincident
#'   sample mean or leave-one-out mean when the median equals an observation.
#' @param initial_shape Optional positive definite starting shape.
#' @param tol Tolerance for both normalized estimating-equation residuals.
#' @param max_iter Maximum number of alternating iterations.
#' @param zero_tol Tolerance used to detect zero whitened residuals. The
#'   default detects exact zeros only.
#' @param warn If `TRUE`, warn when the iteration does not converge.
#'
#' @return An object of class `hd_hr` containing `location`, `shape`, and
#'   convergence diagnostics. In one dimension the method is defined by the
#'   conventional sample median and unit shape. This univariate extension does
#'   not report the multivariate joint shape-equation residual, which is
#'   returned as `NA`.
#' @export
#'
#' @references
#' Hettmansperger, T. P. and Randles, R. H. (2002). A practical affine
#' equivariant multivariate median. *Biometrika*, 89, 851-860.
#'
#' @examples
#' set.seed(5)
#' x <- matrix(rt(240, df = 4), 60, 4)
#' hr_estimator(x)
hr_estimator <- function(x, initial_location = NULL, initial_shape = NULL,
                         tol = 1e-8, max_iter = 500L,
                         zero_tol = 0, warn = TRUE) {
  x <- .as_data_matrix(x)
  controls <- .validate_iteration_controls(tol, max_iter, zero_tol)
  tol <- controls$tol
  max_iter <- controls$max_iter
  zero_tol <- controls$zero_tol
  n <- nrow(x)
  p <- ncol(x)
  if (n <= p) {
    stop("The unregularized HR estimator requires `n > p`.", call. = FALSE)
  }
  if (p == 1L) {
    if (!is.null(initial_location)) {
      .as_location(initial_location, 1L, "initial_location")
    }
    if (!is.null(initial_shape)) {
      checked_shape <- normalize_shape(initial_shape, "trace")
      if (!identical(dim(checked_shape), c(1L, 1L))) {
        stop("`initial_shape` must match the data dimension.", call. = FALSE)
      }
    }
    location <- stats::median(x[, 1L])
    signs <- sign(x[, 1L] - location)
    n_zero <- sum(signs == 0)
    location_residual <- max(0, abs(sum(signs)) - n_zero) / n
    result <- list(
      location = stats::setNames(location, colnames(x)),
      shape = .shape_dimnames(matrix(1, 1L, 1L), x),
      iterations = 0L,
      best_iteration = 0L,
      converged = TRUE,
      location_change = 0,
      shape_change = 0,
      location_equation_residual = location_residual,
      shape_equation_residual = NA_real_,
      equation_residual = location_residual,
      n_zero_residuals = n_zero,
      equation_convention = "univariate sample median; shape equation not applicable"
    )
    class(result) <- "hd_hr"
    return(result)
  }
  if (is.null(initial_location)) {
    initial_location <- .hr_default_initial(
      x, tol = tol, max_iter = max_iter, zero_tol = zero_tol, warn = warn
    )
  } else {
    initial_location <- .as_location(
      initial_location, p, "initial_location"
    )
  }
  if (is.null(initial_shape)) {
    initial_shape <- diag(p)
  } else {
    initial_shape <- normalize_shape(initial_shape, "trace")
    if (!identical(dim(initial_shape), c(p, p))) {
      stop("`initial_shape` must match the data dimension.", call. = FALSE)
    }
  }
  origin <- initial_location
  scaled <- .center_and_scale(x, origin, zero_tol)
  result <- cpp_hr_estimator(
    scaled$x, numeric(p), initial_shape, as.numeric(tol),
    as.integer(max_iter), scaled$zero_tol
  )
  if (isTRUE(warn) && !isTRUE(result$converged)) {
    warning("The HR iteration did not converge in `max_iter` steps; " %+%
              "the returned iterate had the smallest observed equation residual.",
            call. = FALSE)
  }
  displacement <- as.numeric(result$location)
  if (scaled$overflow_fallback) {
    result$location <- (origin / scaled$scale + displacement) * scaled$scale
  } else {
    result$location <- origin + displacement * scaled$scale
  }
  names(result$location) <- colnames(x)
  result$shape <- .shape_dimnames(result$shape, x)
  class(result) <- "hd_hr"
  result
}


#' Angular central Gaussian log-likelihood
#'
#' Evaluates the angular central Gaussian log-likelihood for nonzero row
#' vectors. Row lengths need not equal one because they cancel from the
#' quadratic contribution.
#'
#' @param x Nonzero directions in rows.
#' @param shape A positive definite shape matrix.
#' @param include_constant Include the density normalizing constant.
#'
#' @return A scalar log-likelihood.
#' @references
#' Tyler, D. E. (1987). A distribution-free M-estimator of multivariate
#' scatter. *Annals of Statistics*, 15, 234--251.
#' \doi{10.1214/aos/1176350263}.
#' @export
#'
#' @examples
#' directions <- rbind(c(1, 0), c(0, 1), c(1, 1))
#' acg_loglik(directions, diag(2))
acg_loglik <- function(x, shape, include_constant = FALSE) {
  x <- .as_data_matrix(x)
  p <- ncol(x)
  shape <- normalize_shape(shape, "trace")
  if (!identical(dim(shape), c(p, p))) {
    stop("`shape` must match the dimension of `x`.", call. = FALSE)
  }
  value <- cpp_acg_loglik(x, shape)
  if (isTRUE(include_constant)) {
    value <- value + nrow(x) * (
      lgamma(p / 2) - log(2) - p * log(pi) / 2
    )
  }
  value
}


#' Simulate an elliptically symmetric sample
#'
#' Uses the stochastic representation \eqn{X=\mu+\xi A U}, where `shape` is
#' \eqn{AA^T}, \eqn{U} is uniform on the unit sphere, and `radial` supplies the
#' non-negative radii. The default chi radius yields a multivariate Gaussian
#' sample with covariance `shape`.
#'
#' @param n Number of observations.
#' @param location Location vector. A scalar zero is expanded when `shape` has
#'   dimension greater than one.
#' @param shape Positive semidefinite scatter matrix.
#' @param radial Either `NULL`, a function accepting `n`, or a non-negative
#'   numeric vector of length `n`.
#'
#' @return An `n` by `p` numeric matrix.
#' @references
#' Fang, K.-T. and Anderson, T. W. (1990). *Statistical Inference in
#' Elliptically Contoured and Related Distributions*. Allerton Press.
#' @export
#'
#' @examples
#' set.seed(6)
#' relliptical(5, location = c(1, -1), shape = diag(c(1, 4)))
relliptical <- function(n, location = 0,
                        shape = diag(length(location)), radial = NULL) {
  if (!is.numeric(n) || length(n) != 1L || is.na(n) || !is.finite(n) ||
      n < 1 || n != floor(n) || n > .Machine$integer.max) {
    stop("`n` must be a positive integer.", call. = FALSE)
  }
  n <- as.integer(n)
  shape <- as.matrix(shape)
  storage.mode(shape) <- "double"
  if (nrow(shape) != ncol(shape) || nrow(shape) < 1L ||
      anyNA(shape) || any(!is.finite(shape))) {
    stop("`shape` must be a finite square numeric matrix.", call. = FALSE)
  }
  p <- nrow(shape)
  if (length(location) == 1L && isTRUE(location == 0) && p > 1L) {
    location <- rep.int(0, p)
  }
  location <- .as_location(location, p, "location")
  matrix_scale <- max(abs(shape))
  if (matrix_scale == 0) {
    square_root <- matrix(0, p, p)
  } else {
    scaled_shape <- shape / matrix_scale
    shape <- (scaled_shape + t(scaled_shape)) / 2
    decomposition <- eigen(shape, symmetric = TRUE)
    scale_tol <- 100 * .Machine$double.eps *
      max(1, max(abs(decomposition$values)))
    if (min(decomposition$values) < -scale_tol) {
      stop("`shape` must be positive semidefinite.", call. = FALSE)
    }
    square_root <- sqrt(matrix_scale) * decomposition$vectors %*%
      (sqrt(pmax(decomposition$values, 0)) * t(decomposition$vectors))
  }

  directions <- matrix(stats::rnorm(n * p), n, p)
  lengths <- sqrt(rowSums(directions^2))
  while (any(lengths == 0)) {
    zero <- lengths == 0
    directions[zero, ] <- matrix(
      stats::rnorm(sum(zero) * p), sum(zero), p
    )
    lengths[zero] <- sqrt(rowSums(directions[zero, , drop = FALSE]^2))
  }
  directions <- directions / lengths

  if (is.null(radial)) {
    radii <- sqrt(stats::rchisq(n, df = p))
  } else if (is.function(radial)) {
    radii <- radial(n)
  } else {
    radii <- radial
  }
  radii <- as.numeric(radii)
  if (length(radii) != n || anyNA(radii) || any(!is.finite(radii)) ||
      any(radii < 0)) {
    stop("`radial` must produce `n` finite non-negative radii.",
         call. = FALSE)
  }
  sample <- (directions * radii) %*% t(square_root)
  sample <- sweep(sample, 2L, location, "+", check.margin = FALSE)
  colnames(sample) <- names(location)
  sample
}


#' Simulate a spherically symmetric sample
#'
#' @param n Number of observations.
#' @param p Dimension.
#' @param location Location vector.
#' @param radial Radial specification passed to [relliptical()].
#'
#' @return An `n` by `p` numeric matrix.
#' @references
#' Fang, K.-T. and Anderson, T. W. (1990). *Statistical Inference in
#' Elliptically Contoured and Related Distributions*. Allerton Press.
#' @export
#'
#' @examples
#' set.seed(7)
#' rspherical(4, p = 3, radial = rep(1, 4))
rspherical <- function(n, p, location = rep.int(0, p), radial = NULL) {
  if (!is.numeric(p) || length(p) != 1L || is.na(p) || !is.finite(p) ||
      p < 1 || p != floor(p) || p > .Machine$integer.max) {
    stop("`p` must be a positive integer.", call. = FALSE)
  }
  p <- as.integer(p)
  relliptical(n, location = location, shape = diag(p), radial = radial)
}
