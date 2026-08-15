.clx_extreme_value_cdf <- function(q) {
  q <- as.numeric(q)
  log_intensity <- -0.5 * log(pi) - q / 2
  answer <- rep.int(NA_real_, length(q))
  finite_or_infinite <- !is.na(log_intensity)
  overflow <- finite_or_infinite &
    log_intensity > log(.Machine$double.xmax)
  answer[overflow] <- 0
  ordinary <- finite_or_infinite & !overflow
  answer[ordinary] <- exp(-exp(log_intensity[ordinary]))
  answer
}

.clx_extreme_value_survival <- function(q) {
  q <- as.numeric(q)
  log_intensity <- -0.5 * log(pi) - q / 2
  answer <- rep.int(NA_real_, length(q))
  finite_or_infinite <- !is.na(log_intensity)
  overflow <- finite_or_infinite &
    log_intensity > log(.Machine$double.xmax)
  answer[overflow] <- 1
  ordinary <- finite_or_infinite & !overflow
  intensity <- exp(log_intensity[ordinary])
  answer[ordinary] <- -expm1(-intensity)
  answer
}

.clx_extreme_value_quantile <- function(probability) {
  probability <- as.numeric(probability)
  if (anyNA(probability) || any(probability < 0 | probability > 1)) {
    stop("`probability` must contain values in [0, 1].", call. = FALSE)
  }
  answer <- rep.int(NA_real_, length(probability))
  answer[probability == 0] <- -Inf
  answer[probability == 1] <- Inf
  interior <- probability > 0 & probability < 1
  answer[interior] <- -log(pi) -
    2 * log(-log(probability[interior]))
  answer
}

.clx_critical_value <- function(alpha, p) {
  alpha <- as.numeric(alpha)
  p <- as.numeric(p)
  if (length(alpha) != 1L || is.na(alpha) || !is.finite(alpha) ||
      alpha <= 0 || alpha >= 1) {
    stop("`alpha` must be one finite number in (0, 1).", call. = FALSE)
  }
  if (length(p) != 1L || is.na(p) || !is.finite(p) ||
      p < 2 || p != floor(p)) {
    stop("`p` must be one integer at least 2.", call. = FALSE)
  }
  2 * log(p) - log(log(p)) +
    .clx_extreme_value_quantile(1 - alpha)
}

.clx_positive_scalar <- function(value, name, upper = Inf,
                                 allow_zero = FALSE) {
  value <- as.numeric(value)
  lower_ok <- if (allow_zero) value >= 0 else value > 0
  if (length(value) != 1L || is.na(value) || !is.finite(value) ||
      !lower_ok || value >= upper) {
    interval <- if (allow_zero) "[0, 1)" else "positive"
    if (is.finite(upper) && !allow_zero) {
      interval <- sprintf("(0, %s)", format(upper))
    }
    stop(sprintf("`%s` must be one finite %s number.", name, interval),
         call. = FALSE)
  }
  value
}

.clx_validate_precision <- function(precision, p, oracle) {
  if (!is.matrix(precision) || !is.numeric(precision) ||
      !identical(dim(precision), c(p, p))) {
    stop(sprintf("`precision` must be a numeric %d by %d matrix.", p, p),
         call. = FALSE)
  }
  storage.mode(precision) <- "double"
  if (anyNA(precision) || any(!is.finite(precision))) {
    stop("`precision` must contain only finite values.", call. = FALSE)
  }

  matrix_scale <- max(abs(precision))
  safe_scale <- max(matrix_scale, .Machine$double.xmin)
  relative_asymmetry <- max(abs(
    precision / safe_scale - t(precision) / safe_scale
  ))
  maximum_asymmetry <- relative_asymmetry * safe_scale
  symmetry_tolerance <- sqrt(.Machine$double.eps) *
    safe_scale
  if (relative_asymmetry > sqrt(.Machine$double.eps)) {
    stop("`precision` must be symmetric within numerical tolerance.",
         call. = FALSE)
  }
  symmetrized <- relative_asymmetry > 0
  precision <- precision / 2 + t(precision) / 2
  eigenvalues <- tryCatch(
    eigen(precision, symmetric = TRUE, only.values = TRUE)$values,
    error = function(error) {
      stop("The eigendecomposition of `precision` failed: ",
           conditionMessage(error), call. = FALSE)
    }
  )
  if (any(!is.finite(eigenvalues))) {
    stop("The eigenvalues of `precision` must be finite.", call. = FALSE)
  }
  positive_definite <- min(eigenvalues) > 0
  if (oracle && !positive_definite) {
    stop("Oracle `precision` must be positive definite.", call. = FALSE)
  }

  list(
    precision = precision,
    diagnostics = list(
      maximum.asymmetry = maximum_asymmetry,
      symmetry.tolerance = symmetry_tolerance,
      symmetrized = symmetrized,
      positive.definite = positive_definite,
      minimum.eigenvalue = min(eigenvalues),
      maximum.eigenvalue = max(eigenvalues),
      eigenvalue.adjustment = "none"
    )
  )
}

#' Cai--Liu--Xia precision-adjusted two-sample maximum test
#'
#' Tests equality of two high-dimensional mean vectors using the
#' precision-adjusted maximum statistic of Cai, Liu, and Xia (2014).  The
#' calibration is the type-I extreme-value limit
#' \deqn{F(g)=\exp\{-\pi^{-1/2}\exp(-g/2)\},}
#' for
#' \deqn{G=M-2\log(p)+\log\{\log(p)\}.}
#'
#' `precision_source` deliberately separates three statistically different
#' paths.  With `"oracle"`, `precision` is a known population precision matrix
#' and the denominator is its diagonal.  With `"feasible"`, `precision` is a
#' user-supplied estimate and the denominator is formed from empirical
#' within-group variances of the transformed observations, each with divisor
#' \eqn{n_k}, as in equations (6)--(7) of the paper.  With `"adaptive"`, the
#' function estimates a precision matrix by hard-thresholding the pooled
#' covariance and then uses that same feasible denominator.  It never uses
#' `diag(precision)` to standardise an estimated-precision statistic.
#'
#' For the adaptive backend, the pooled covariance and
#' \eqn{\widehat\theta_{ij}} both use divisor \eqn{N=n_1+n_2}, and
#' \deqn{\lambda_{ij}=\delta
#'   \sqrt{\widehat\theta_{ij}\log(p)/N}.}
#' If the thresholded covariance is not numerically positive definite, its
#' eigenvalues are raised to `eigen_floor` times a data-scale reference before
#' inversion.  The floor, the unadjusted and adjusted extreme eigenvalues, and
#' the number of modified eigenvalues are all returned in `diagnostics`.
#' Both samples are centered relative to observed anchors before precision
#' transformation, avoiding cancellation under a large common translation.
#' The adaptive covariance and fourth-order threshold moments are also formed
#' on an internal common scale and then mapped back to the original units.
#'
#' @param x,y Numeric matrices or data frames with observations in rows and
#'   the same variables in columns.  Each group needs at least two rows and
#'   the dimension must satisfy \eqn{p\geq 2}.
#' @param precision A finite symmetric \eqn{p\times p} matrix.  It is required
#'   for `precision_source = "oracle"` or `"feasible"` and must be `NULL` for
#'   the adaptive path.  An oracle matrix must be positive definite.  A
#'   feasible estimate may be indefinite, as allowed by the original paper,
#'   provided all transformed empirical variances are strictly positive.
#' @param precision_source One of `"adaptive"`, `"oracle"`, or `"feasible"`.
#'   This argument must identify whether a supplied matrix is known or
#'   estimated because the two cases have different denominators.
#' @param threshold_delta Positive multiplier \eqn{\delta} for adaptive
#'   thresholding.  The paper recommends 2 as a fixed choice.
#' @param eigen_floor A non-negative relative eigenvalue floor below 1 for the
#'   adaptive thresholded covariance.  Set it to zero to prohibit adjustment;
#'   a non-positive thresholded eigenvalue then produces an error.
#'
#' @return An object of class `c("hd_location_test", "htest")`.  The reported
#'   statistic `G` has the extreme-value calibration, while `raw.statistic`
#'   contains `M`.  Coordinatewise transformed scores, their denominators,
#'   the maximizing coordinate index and name, and the precision matrix used
#'   by the test are retained in `components`.  `diagnostics` identifies the
#'   precision path, denominator construction, and every adaptive-threshold
#'   eigenvalue adjustment.
#'
#' @references
#' Cai, T. T., Liu, W., and Xia, Y. (2014). Two-sample test of high
#' dimensional means under dependence. *Journal of the Royal Statistical
#' Society: Series B*, **76**, 349--372.
#'
#' Cai, T. T. and Liu, W. (2011). Adaptive thresholding for sparse covariance
#' matrix estimation. *Journal of the American Statistical Association*,
#' **106**, 672--684.
#'
#' @examples
#' set.seed(29)
#' x <- matrix(rnorm(80), 20, 4)
#' y <- matrix(rnorm(96, 0.2), 24, 4)
#' cai_liu_xia_two_sample_test(x, y)
#'
#' omega <- diag(4)
#' cai_liu_xia_two_sample_test(
#'   x, y, precision = omega, precision_source = "oracle"
#' )
#'
#' @export
cai_liu_xia_two_sample_test <- function(
    x, y, precision = NULL,
    precision_source = c("adaptive", "oracle", "feasible"),
    threshold_delta = 2,
    eigen_floor = sqrt(.Machine$double.eps)) {
  call <- match.call()
  x_name <- deparse1(substitute(x))
  y_name <- deparse1(substitute(y))
  x <- .as_data_matrix(x, "x", min_rows = 2L)
  y <- .as_data_matrix(y, "y", min_rows = 2L)
  .check_two_sample_variables(x, y)
  n1 <- nrow(x)
  n2 <- nrow(y)
  p <- ncol(x)
  if (p < 2L) {
    stop("The CLX extreme-value calibration requires `p >= 2`.",
         call. = FALSE)
  }

  precision_source <- match.arg(precision_source)
  threshold_delta <- .clx_positive_scalar(
    threshold_delta, "threshold_delta"
  )
  eigen_floor <- .clx_positive_scalar(
    eigen_floor, "eigen_floor", upper = 1, allow_zero = TRUE
  )
  variable_names <- .hotelling_variable_names(x)

  adaptive <- identical(precision_source, "adaptive")
  oracle <- identical(precision_source, "oracle")
  adaptive_fit <- NULL
  if (adaptive) {
    if (!is.null(precision)) {
      stop(
        "`precision` must be `NULL` when `precision_source = \"adaptive\"`.",
        call. = FALSE
      )
    }
    adaptive_fit <- cpp_clx_adaptive_precision(
      x, y, threshold_delta, eigen_floor
    )
    precision <- as.matrix(adaptive_fit$precision)
    precision_diagnostics <- list(
      maximum.asymmetry = max(abs(precision - t(precision))),
      symmetry.tolerance = sqrt(.Machine$double.eps) *
        max(abs(precision), .Machine$double.xmin),
      symmetrized = FALSE,
      positive.definite = TRUE,
      minimum.eigenvalue =
        1 / as.numeric(adaptive_fit$maximum_eigenvalue_after),
      maximum.eigenvalue =
        1 / as.numeric(adaptive_fit$minimum_eigenvalue_after),
      eigenvalue.adjustment = if (isTRUE(adaptive_fit$adjustment_applied)) {
        "thresholded covariance eigenvalue floor"
      } else {
        "none"
      }
    )
  } else {
    if (is.null(precision)) {
      stop(sprintf(
        "`precision` is required when `precision_source = \"%s\"`.",
        precision_source
      ), call. = FALSE)
    }
    checked_precision <- .clx_validate_precision(precision, p, oracle)
    precision <- checked_precision$precision
    precision_diagnostics <- checked_precision$diagnostics
  }
  dimnames(precision) <- list(variable_names, variable_names)

  fit <- cpp_clx_two_sample(x, y, precision, oracle)
  maximum <- as.numeric(fit$M)
  centered_maximum <- as.numeric(fit$G)
  p_value <- .clx_extreme_value_survival(centered_maximum)
  argmax_index <- as.integer(fit$argmax)
  argmax_name <- variable_names[[argmax_index]]

  mean_x <- stats::setNames(as.numeric(fit$mean_x), variable_names)
  mean_y <- stats::setNames(as.numeric(fit$mean_y), variable_names)
  difference <- stats::setNames(
    as.numeric(fit$difference), variable_names
  )
  transformed_difference <- stats::setNames(
    as.numeric(fit$transformed_difference), variable_names
  )
  denominator <- stats::setNames(
    as.numeric(fit$denominator), variable_names
  )
  scores <- stats::setNames(as.numeric(fit$scores), variable_names)
  squared_scores <- stats::setNames(
    as.numeric(fit$squared_scores), variable_names
  )

  components <- list(
    mean.x = mean_x,
    mean.y = mean_y,
    difference = difference,
    transformed.difference = transformed_difference,
    denominator = denominator,
    standardized.scores = scores,
    squared.scores = squared_scores,
    maximum.coordinate = list(index = argmax_index, name = argmax_name),
    precision = precision,
    effective.sample.size = as.numeric(fit$effective_n),
    M = maximum,
    G = centered_maximum
  )
  if (!oracle) {
    components$transformed.variance.x <- stats::setNames(
      as.numeric(fit$variance_x), variable_names
    )
    components$transformed.variance.y <- stats::setNames(
      as.numeric(fit$variance_y), variable_names
    )
  }

  adaptive_diagnostics <- NULL
  if (adaptive) {
    adaptive_diagnostics <- list(
      pooled.covariance.denominator = n1 + n2,
      theta.denominator = n1 + n2,
      threshold.sample.size = n1 + n2,
      threshold.formula =
        "delta * sqrt(theta_hat_ij * log(p) / (n1 + n2))",
      threshold.delta = as.numeric(adaptive_fit$threshold_delta),
      internal.data.scale = as.numeric(adaptive_fit$internal_data_scale),
      minimum.threshold = as.numeric(adaptive_fit$minimum_threshold),
      maximum.threshold = as.numeric(adaptive_fit$maximum_threshold),
      retained.nonzero.entries =
        as.integer(adaptive_fit$retained_nonzero_entries),
      retained.off.diagonal.entries =
        as.integer(adaptive_fit$retained_off_diagonal_entries),
      eigen.floor.relative =
        as.numeric(adaptive_fit$relative_eigen_floor),
      eigen.floor.absolute =
        as.numeric(adaptive_fit$absolute_eigen_floor),
      eigen.floor.applied = isTRUE(adaptive_fit$adjustment_applied),
      adjusted.eigenvalues =
        as.integer(adaptive_fit$adjusted_eigenvalues),
      minimum.eigenvalue.before =
        as.numeric(adaptive_fit$minimum_eigenvalue_before),
      maximum.eigenvalue.before =
        as.numeric(adaptive_fit$maximum_eigenvalue_before),
      minimum.eigenvalue.after =
        as.numeric(adaptive_fit$minimum_eigenvalue_after),
      maximum.eigenvalue.after =
        as.numeric(adaptive_fit$maximum_eigenvalue_after)
    )
  }

  .new_hd_location_test(
    statistic = c(G = centered_maximum),
    parameter = c(p = p),
    p.value = p_value,
    method = paste(
      "Cai-Liu-Xia two-sample precision-adjusted maximum test",
      sprintf("(%s precision path)", precision_source)
    ),
    data.name = paste(x_name, "and", y_name),
    alternative = "two.sided",
    raw.statistic = c(M = maximum),
    estimate = difference,
    null.value = stats::setNames(numeric(p), variable_names),
    null.distribution = list(
      family = "CLX type-I extreme-value",
      parameters = c(intensity = 1 / sqrt(pi), scale = 2),
      cdf = "exp(-pi^(-1/2) * exp(-g / 2))",
      centering = "M - 2 * log(p) + log(log(p))",
      exact = FALSE,
      tail = "upper",
      assumptions = paste(
        "Cai-Liu-Xia high-dimensional Gaussian common-covariance",
        "and weak transformed-dependence conditions"
      )
    ),
    variance = denominator,
    components = components,
    diagnostics = list(
      precision.source = precision_source,
      statistic.path = if (oracle) "oracle" else "feasible",
      denominator.source = if (oracle) {
        "diagonal of known population precision"
      } else {
        paste(
          "transformed empirical within-group variances",
          "with divisors n1 and n2"
        )
      },
      transformed.variance.divisors = if (oracle) {
        NULL
      } else {
        c(x = n1, y = n2)
      },
      argmax.index = argmax_index,
      argmax.name = argmax_name,
      precision = precision_diagnostics,
      adaptive.thresholding = adaptive_diagnostics,
      regularization = if (adaptive) {
        "adaptive hard thresholding, with diagnosed eigen-floor if needed"
      } else {
        "none"
      },
      variance.repair = "none"
    ),
    n = c(x = n1, y = n2),
    p = p,
    call = call
  )
}
