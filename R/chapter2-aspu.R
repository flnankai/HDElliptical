.aspu_validate_powers <- function(powers) {
  if (!is.numeric(powers) || length(powers) < 1L || anyNA(powers)) {
    stop("`powers` must be a non-empty numeric vector without missing values.",
         call. = FALSE)
  }
  finite <- is.finite(powers)
  if (any(!finite & powers != Inf) ||
      any(powers[finite] <= 0) ||
      any(powers[finite] != floor(powers[finite])) ||
      any(powers[finite] > .Machine$integer.max)) {
    stop(
      "Finite `powers` must be positive integers; the only allowed " %+%
        "infinite value is `Inf`.",
      call. = FALSE
    )
  }
  if (anyDuplicated(powers)) {
    stop("`powers` must not contain duplicates.", call. = FALSE)
  }
  finite_powers <- sort(as.integer(powers[finite]))
  c(as.numeric(finite_powers), if (any(!finite)) Inf else numeric())
}

.aspu_validate_integer <- function(value, name, minimum = 1L,
                                   maximum = .Machine$integer.max) {
  if (!is.numeric(value) || length(value) != 1L || is.na(value) ||
      !is.finite(value) || value != floor(value) ||
      value < minimum || value > maximum) {
    stop(sprintf(
      "`%s` must be one integer in %d, ..., %d.",
      name, minimum, maximum
    ), call. = FALSE)
  }
  as.integer(value)
}

.aspu_validate_bandwidth <- function(bandwidth, p, source) {
  if (is.null(bandwidth)) {
    return(c(first = -1L, second = -1L))
  }
  if (identical(source, "supplied")) {
    stop(
      "`bandwidth` must be `NULL` when `correlation_source = \"supplied\"`; " %+%
        "band the covariance estimate before supplying its correlation.",
      call. = FALSE
    )
  }
  expected <- if (identical(source, "unequal")) c(1L, 2L) else 1L
  if (!is.numeric(bandwidth) || anyNA(bandwidth) ||
      any(!is.finite(bandwidth)) ||
      !(length(bandwidth) %in% expected) ||
      any(bandwidth != floor(bandwidth)) ||
      any(bandwidth < 0) || any(bandwidth >= p)) {
    length_text <- if (identical(source, "unequal")) {
      "one value (recycled) or two values"
    } else {
      "one value"
    }
    stop(sprintf(
      "`bandwidth` must contain %s in 0, ..., p - 1.", length_text
    ), call. = FALSE)
  }
  if (length(bandwidth) == 1L) {
    bandwidth <- rep.int(bandwidth, 2L)
  }
  stats::setNames(as.integer(bandwidth[seq_len(2L)]),
                  c("first", "second"))
}

.aspu_validate_supplied_correlation <- function(correlation, p,
                                                variable_names) {
  if (!is.matrix(correlation) || !is.numeric(correlation) ||
      !identical(dim(correlation), c(p, p))) {
    stop(sprintf(
      "`supplied_correlation` must be a numeric %d by %d matrix.", p, p
    ), call. = FALSE)
  }
  storage.mode(correlation) <- "double"
  if (anyNA(correlation) || any(!is.finite(correlation))) {
    stop("`supplied_correlation` must contain only finite values.",
         call. = FALSE)
  }
  supplied_names <- dimnames(correlation)
  for (margin_names in supplied_names) {
    if (!is.null(margin_names) &&
        !identical(margin_names, variable_names)) {
      stop(
        "Dimnames of `supplied_correlation` must match the sample " %+%
          "variable names in order.",
        call. = FALSE
      )
    }
  }
  dimnames(correlation) <- list(variable_names, variable_names)
  correlation
}

.aspu_validate_standard_errors <- function(standard_errors, p,
                                           variable_names) {
  if (is.null(standard_errors)) {
    return(numeric())
  }
  input_names <- names(standard_errors)
  standard_errors <- as.numeric(standard_errors)
  if (length(standard_errors) != p || anyNA(standard_errors) ||
      any(!is.finite(standard_errors)) || any(standard_errors <= 0)) {
    stop(sprintf(
      "`standard_errors` must contain %d finite, strictly positive values.",
      p
    ), call. = FALSE)
  }
  if (!is.null(input_names) && !identical(input_names, variable_names)) {
    stop(
      "Named `standard_errors` must use the sample variable names in order.",
      call. = FALSE
    )
  }
  standard_errors
}

.aspu_psd_correlation <- function(correlation, psd_adjust, psd_tol,
                                  label) {
  correlation_dimnames <- dimnames(correlation)
  dimension <- nrow(correlation)
  matrix_scale <- max(abs(correlation), .Machine$double.xmin)
  relative_asymmetry <- max(abs(
    correlation / matrix_scale - t(correlation) / matrix_scale
  ))
  symmetry_tolerance <- sqrt(.Machine$double.eps)
  if (relative_asymmetry > symmetry_tolerance) {
    stop(sprintf("The %s correlation matrix must be symmetric.", label),
         call. = FALSE)
  }
  maximum_diagonal_error <- max(abs(diag(correlation) - 1))
  if (maximum_diagonal_error > symmetry_tolerance) {
    stop(sprintf("The %s correlation matrix must have a unit diagonal.",
                 label), call. = FALSE)
  }
  correlation <- correlation / 2 + t(correlation) / 2
  diag(correlation) <- 1
  if (any(abs(correlation) > 1 + symmetry_tolerance)) {
    stop(sprintf(
      "Every entry of the %s correlation matrix must lie in [-1, 1].",
      label
    ), call. = FALSE)
  }

  decomposition <- tryCatch(
    eigen(correlation, symmetric = TRUE),
    error = function(error) {
      stop(sprintf(
        "The eigendecomposition of the %s correlation matrix failed: %s",
        label, conditionMessage(error)
      ), call. = FALSE)
    }
  )
  eigenvalues_before <- decomposition$values
  minimum_before <- min(eigenvalues_before)
  adjustment_applied <- FALSE
  frobenius_change <- 0
  if (identical(psd_adjust, "error")) {
    if (minimum_before < -psd_tol) {
      stop(sprintf(
        paste0(
          "The %s correlation matrix is not positive semidefinite ",
          "(minimum eigenvalue %.6g). Set `psd_adjust = \"eigen_clip\"` ",
          "to request an explicit repair."
        ),
        label, minimum_before
      ), call. = FALSE)
    }
  } else if (minimum_before < psd_tol) {
    clipped <- pmax(eigenvalues_before, psd_tol)
    adjusted <- decomposition$vectors %*%
      (clipped * t(decomposition$vectors))
    marginal_scale <- sqrt(diag(adjusted))
    adjusted <- adjusted / tcrossprod(marginal_scale)
    adjusted <- adjusted / 2 + t(adjusted) / 2
    diag(adjusted) <- 1
    frobenius_change <- sqrt(sum((adjusted - correlation)^2))
    correlation <- adjusted
    adjustment_applied <- TRUE
  }
  eigenvalues_after <- eigen(
    correlation, symmetric = TRUE, only.values = TRUE
  )$values
  dimnames(correlation) <- correlation_dimnames
  list(
    correlation = correlation,
    diagnostics = list(
      label = label,
      action = psd_adjust,
      tolerance = psd_tol,
      adjusted = adjustment_applied,
      minimum.eigenvalue.before = minimum_before,
      maximum.eigenvalue.before = max(eigenvalues_before),
      minimum.eigenvalue.after = min(eigenvalues_after),
      maximum.eigenvalue.after = max(eigenvalues_after),
      frobenius.adjustment = frobenius_change,
      maximum.relative.asymmetry = relative_asymmetry,
      maximum.diagonal.error = maximum_diagonal_error
    )
  )
}

.aspu_miwa_group <- function(z, correlation, tail, miwa_steps,
                             psd_adjust, psd_tol, power_names) {
  count <- length(z)
  if (count == 1L) {
    p_value <- if (identical(tail, "two.sided")) {
      2 * stats::pnorm(abs(z), lower.tail = FALSE)
    } else {
      stats::pnorm(z, lower.tail = FALSE)
    }
    return(list(
      statistic = if (identical(tail, "two.sided")) abs(z) else z,
      p.value = p_value,
      correlation = matrix(1, 1L, 1L,
                           dimnames = list(power_names, power_names)),
      integration = list(
        algorithm = "univariate normal", deterministic = TRUE,
        absolute.error = 0, message = "exact univariate normal tail"
      ),
      psd = list(
        label = paste(tail, "power family"), action = "not needed",
        tolerance = psd_tol, adjusted = FALSE,
        minimum.eigenvalue.before = 1,
        maximum.eigenvalue.before = 1,
        minimum.eigenvalue.after = 1,
        maximum.eigenvalue.after = 1,
        frobenius.adjustment = 0,
        maximum.relative.asymmetry = 0,
        maximum.diagonal.error = 0
      )
    ))
  }
  if (count > 20L) {
    stop(
      "Deterministic `mvtnorm::Miwa` calibration supports at most 20 " %+%
        "powers in each odd or even family.",
      call. = FALSE
    )
  }
  dimnames(correlation) <- list(power_names, power_names)
  checked <- .aspu_psd_correlation(
    correlation, psd_adjust, psd_tol,
    paste(tail, "finite-power")
  )
  correlation <- checked$correlation
  threshold <- if (identical(tail, "two.sided")) max(abs(z)) else max(z)
  lower <- if (identical(tail, "two.sided")) {
    rep.int(-threshold, count)
  } else {
    rep.int(-Inf, count)
  }
  upper <- rep.int(threshold, count)
  if (!requireNamespace("mvtnorm", quietly = TRUE)) {
    stop(
      "Package `mvtnorm` is required for multivariate analytical aSPU " %+%
        "calibration.",
      call. = FALSE
    )
  }
  probability <- tryCatch(
    mvtnorm::pmvnorm(
      lower = lower, upper = upper, mean = numeric(count),
      corr = correlation,
      algorithm = mvtnorm::Miwa(steps = miwa_steps, checkCorr = TRUE),
      keepAttr = TRUE
    ),
    error = function(error) {
      stop(sprintf(
        paste0(
          "Deterministic Miwa integration failed for the %s family: %s. ",
          "A singular numerical correlation can be handled explicitly with ",
          "`psd_adjust = \"eigen_clip\"`."
        ),
        tail, conditionMessage(error)
      ), call. = FALSE)
    }
  )
  probability_value <- as.numeric(probability)
  if (length(probability_value) != 1L || !is.finite(probability_value) ||
      probability_value < -1e-10 || probability_value > 1 + 1e-10) {
    stop("Deterministic Miwa integration returned an invalid probability.",
         call. = FALSE)
  }
  probability_value <- min(1, max(0, probability_value))
  p_value <- min(1, max(0, 1 - probability_value))
  list(
    statistic = threshold,
    p.value = p_value,
    correlation = correlation,
    integration = list(
      algorithm = "mvtnorm::Miwa", deterministic = TRUE,
      steps = miwa_steps,
      absolute.error = as.numeric(attr(probability, "error")),
      message = as.character(attr(probability, "msg"))
    ),
    psd = checked$diagnostics
  )
}

.aspu_extreme_survival <- function(g) {
  log_intensity <- -0.5 * log(pi) - g / 2
  if (log_intensity > log(.Machine$double.xmax)) {
    return(1)
  }
  -expm1(-exp(log_intensity))
}

.aspu_rescale_power_values <- function(values, log_multiplier) {
  if (!identical(length(values), length(log_multiplier))) {
    stop("Internal aSPU power-rescaling dimensions do not agree.",
         call. = FALSE)
  }
  result <- values
  nonzero <- values != 0
  if (!any(nonzero)) {
    return(result)
  }
  log_magnitude <- log(abs(values[nonzero])) + log_multiplier[nonzero]
  sign_value <- sign(values[nonzero])
  maximum_log <- log(.Machine$double.xmax)
  minimum_log <- log(.Machine$double.xmin) + log(.Machine$double.eps)
  rescaled <- numeric(length(log_magnitude))
  rescaled[log_magnitude > maximum_log] <- Inf
  ordinary <- log_magnitude <= maximum_log & log_magnitude >= minimum_log
  rescaled[ordinary] <- exp(log_magnitude[ordinary])
  result[nonzero] <- sign_value * rescaled
  result
}

#' Xu--Lin--Wei--Pan analytical adaptive sum-of-powers test
#'
#' Tests equality of two high-dimensional mean vectors with an analytical
#' version of the adaptive sum-of-powers (aSPU) test of Xu, Lin, Wei, and Pan
#' (2016). Observations are rows and variables are columns. Write the
#' inverse-variance-standardised differences as
#' \deqn{W_j={\bar X_{1j}-\bar X_{2j}\over
#' \{\widehat\sigma_{1,jj}/n_1+\widehat\sigma_{2,jj}/n_2\}^{1/2}},}
#' with the common-covariance version replacing the denominator by the pooled
#' marginal variance times \eqn{1/n_1+1/n_2}. With the default
#' `score_scale = "paper_raw"`, the primary paper's finite-power statistic is
#' \eqn{L_\gamma=\sum_j(\bar X_{1j}-\bar X_{2j})^\gamma}; its Gaussian
#' moments use the estimated covariance matrix of the mean difference. With
#' `score_scale = "book_studentized"`, the scale-invariant book variant is
#' \eqn{T_\gamma=\sum_j W_j^\gamma}, whose moments use the correlation matrix
#' of \eqn{W}. Both paths use \eqn{M=\max_j|W_j|} for the infinite power.
#'
#' The two finite-power paths are deliberately explicit because they are not
#' generally numerically equivalent: the primary-paper path is sensitive to
#' coordinate units, whereas the book path is invariant to positive diagonal
#' rescaling. For numerical stability, raw scores and their covariance are
#' divided by one common reported scale before the C++ moment calculation;
#' reported finite-power moments are mapped back to the original units, while
#' standardized statistics and p-values are unaffected by that normalization.
#'
#' Every finite positive integer \eqn{\gamma} is allowed. In particular, the
#' original method's default \code{1:6} includes odd powers; a restriction to
#' even powers in the accompanying book draft is a transcription error. Under
#' the Gaussian working limit, the standardized odd-power statistics use a
#' joint two-sided normal tail and the standardized even-power statistics use
#' a joint upper normal tail. The maximum statistic uses
#' \deqn{G=M^2-2\log(p)+\log\{\log(p)\},\qquad
#' P(G\leq g)=\exp\{-\pi^{-1/2}\exp(-g/2)\}.}
#' The odd, even, and maximum groups are asymptotically independent. If
#' \eqn{m} non-empty groups were requested, the final p-value is
#' \eqn{1-\{1-\min(P_O,P_E,P_\infty)\}^m}. Thus a subset of power families is
#' combined with its actual group count rather than always using exponent 3.
#'
#' Finite-power Gaussian moments are evaluated exactly from Isserlis pairings
#' in C++, and multivariate normal rectangles are evaluated deterministically
#' with `mvtnorm::Miwa`. This function deliberately implements no permutation
#' or parametric-bootstrap calibration. `bandwidth` applies the paper's hard
#' covariance band, setting entries with \eqn{|j-k|} larger than the supplied
#' width to zero. A banded estimate need not be positive semidefinite:
#' `psd_adjust = "error"` rejects it, while `"eigen_clip"` performs and reports
#' an explicit eigenvalue clipping repair. No silent ridge, absolute-value
#' variance repair, or variance floor is used.
#'
#' @param x,y Numeric matrices or data frames with observations in rows and
#'   the same variables in columns. Each group needs at least two rows.
#' @param powers Distinct positive integer powers, optionally including
#'   `Inf`. The default is `c(1:6, Inf)`.
#' @param score_scale Finite-power coordinate definition. `"paper_raw"`
#'   (the default) uses the primary paper's unstandardised sample-mean
#'   differences. `"book_studentized"` uses the scale-invariant
#'   inverse-variance-standardised coordinates in the accompanying book.
#' @param correlation_source How to estimate the null correlation of the
#'   coordinate contrasts: pooled common covariance, unequal group
#'   covariances, or a supplied correlation matrix.
#' @param supplied_correlation A finite \eqn{p\times p} correlation matrix,
#'   required only for `correlation_source = "supplied"`.
#' @param standard_errors Optional finite positive coordinate standard errors
#'   for the supplied-correlation path. If `NULL`, unequal-covariance sample
#'   standard errors are used.
#' @param bandwidth Optional fixed hard-band width. For a common covariance it
#'   is one integer in \eqn{0,\ldots,p-1}; for unequal covariances it may be one
#'   recycled width or two group-specific widths. `NULL` leaves estimates
#'   unbanded. It is unavailable for a supplied correlation.
#' @param psd_adjust Either `"error"` or `"eigen_clip"`. The latter explicitly
#'   clips eigenvalues below `psd_tol` and reports the adjustment.
#' @param psd_tol Positive eigenvalue tolerance used by `psd_adjust`.
#' @param miwa_steps Positive integer grid size for deterministic Miwa
#'   integration. Miwa supports at most 20 selected powers per parity family.
#'
#' @return An object of class `c("hd_location_test", "htest")`. `components`
#'   retains coordinate contrasts, all SPU statistics and individual p-values,
#'   Gaussian means/covariances/correlations, group statistics and p-values,
#'   and the maximum-statistic calibration. `diagnostics` records covariance,
#'   banding, PSD, numerical-integration, and group-combination choices.
#'
#' @references
#' Xu, G., Lin, L., Wei, P., and Pan, W. (2016). An adaptive two-sample test
#' for high-dimensional means. *Biometrika*, **103**, 609--624.
#' \doi{10.1093/biomet/asw029}
#'
#' @examples
#' set.seed(83)
#' x <- matrix(rnorm(80), 20, 4)
#' y <- matrix(rnorm(88, 0.15), 22, 4)
#' xu_lin_wei_pan_aspu_test(x, y)
#'
#' @export
xu_lin_wei_pan_aspu_test <- function(
    x, y, powers = c(1:6, Inf),
    score_scale = c("paper_raw", "book_studentized"),
    correlation_source = c("common", "unequal", "supplied"),
    supplied_correlation = NULL, standard_errors = NULL,
    bandwidth = NULL,
    psd_adjust = c("error", "eigen_clip"),
    psd_tol = sqrt(.Machine$double.eps), miwa_steps = 128L) {
  call <- match.call()
  x_name <- deparse1(substitute(x))
  y_name <- deparse1(substitute(y))
  x <- .as_data_matrix(x, "x", min_rows = 2L)
  y <- .as_data_matrix(y, "y", min_rows = 2L)
  .check_two_sample_variables(x, y)
  n1 <- nrow(x)
  n2 <- nrow(y)
  p <- ncol(x)
  variable_names <- .hotelling_variable_names(x)

  powers <- .aspu_validate_powers(powers)
  score_scale <- match.arg(score_scale)
  finite_powers <- powers[is.finite(powers)]
  has_infinite_power <- any(is.infinite(powers))
  if (has_infinite_power && p < 2L) {
    stop("The infinite-power extreme-value calibration requires `p >= 2`.",
         call. = FALSE)
  }
  correlation_source <- match.arg(correlation_source)
  psd_adjust <- match.arg(psd_adjust)
  psd_tol <- as.numeric(psd_tol)
  if (length(psd_tol) != 1L || is.na(psd_tol) || !is.finite(psd_tol) ||
      psd_tol <= 0 || psd_tol >= 1) {
    stop("`psd_tol` must be one finite number strictly between zero and one.",
         call. = FALSE)
  }
  miwa_steps <- .aspu_validate_integer(miwa_steps, "miwa_steps")
  widths <- .aspu_validate_bandwidth(bandwidth, p, correlation_source)

  if (identical(correlation_source, "supplied")) {
    if (is.null(supplied_correlation)) {
      stop(
        "`supplied_correlation` is required when " %+%
          "`correlation_source = \"supplied\"`.",
        call. = FALSE
      )
    }
    supplied_correlation <- .aspu_validate_supplied_correlation(
      supplied_correlation, p, variable_names
    )
    standard_errors <- .aspu_validate_standard_errors(
      standard_errors, p, variable_names
    )
  } else {
    if (!is.null(supplied_correlation) || !is.null(standard_errors)) {
      stop(
        "`supplied_correlation` and `standard_errors` are only available " %+%
          "when `correlation_source = \"supplied\"`.",
        call. = FALSE
      )
    }
    supplied_correlation <- matrix(numeric(), 0L, 0L)
    standard_errors <- numeric()
  }

  standardized <- cpp_aspu_standardize(
    x, y, correlation_source,
    unname(widths[["first"]]), unname(widths[["second"]]),
    supplied_correlation, standard_errors
  )
  coordinate_correlation <- as.matrix(standardized$correlation)
  dimnames(coordinate_correlation) <- list(variable_names, variable_names)
  checked_coordinate <- .aspu_psd_correlation(
    coordinate_correlation, psd_adjust, psd_tol, "coordinate"
  )
  coordinate_correlation <- checked_coordinate$correlation

  coordinate_names <- variable_names
  mean_x <- stats::setNames(as.numeric(standardized$mean_x), coordinate_names)
  mean_y <- stats::setNames(as.numeric(standardized$mean_y), coordinate_names)
  difference <- stats::setNames(
    as.numeric(standardized$difference), coordinate_names
  )
  standard_error <- stats::setNames(
    as.numeric(standardized$standard_error), coordinate_names
  )
  w <- stats::setNames(as.numeric(standardized$W), coordinate_names)

  finite_fit <- NULL
  finite_components <- list()
  individual_p_values <- numeric()
  group_results <- list()
  family_psd <- list()
  family_integration <- list()
  if (length(finite_powers) > 0L) {
    if (identical(score_scale, "paper_raw")) {
      if (any(!is.finite(standard_error)) || any(standard_error <= 0)) {
        stop(
          paste(
            "The primary-paper raw finite-power path requires marginal",
            "standard errors representable in double precision; use",
            "`score_scale = \"book_studentized\"` for more extreme mixed units."
          ),
          call. = FALSE
        )
      }
      finite_normalizer <- max(standard_error)
      relative_standard_error <- standard_error / finite_normalizer
      if (any(!is.finite(relative_standard_error)) ||
          any(relative_standard_error <= 0)) {
        stop(
          paste(
            "The primary-paper raw finite-power covariance has a marginal",
            "scale outside double precision after common normalization; use",
            "`score_scale = \"book_studentized\"` for these mixed units."
          ),
          call. = FALSE
        )
      }
      finite_score <- difference
      finite_score_normalized <- as.numeric(w) *
        as.numeric(relative_standard_error)
      finite_null_coordinate_covariance_normalized <-
        coordinate_correlation * tcrossprod(relative_standard_error)
      finite_score_definition <- "unstandardised sample-mean difference"
    } else {
      finite_normalizer <- 1
      relative_standard_error <- stats::setNames(
        rep.int(1, p), coordinate_names
      )
      finite_score <- w
      finite_score_normalized <- as.numeric(w)
      finite_null_coordinate_covariance_normalized <-
        coordinate_correlation
      finite_score_definition <-
        "inverse-variance standardised mean difference W"
    }
    finite_fit <- cpp_aspu_power_moments(
      finite_score_normalized,
      finite_null_coordinate_covariance_normalized,
      as.integer(finite_powers)
    )
    finite_labels <- as.character(finite_powers)
    log_normalizer <- log(finite_normalizer)
    observed_normalized <- as.numeric(finite_fit$observed)
    mean_normalized <- as.numeric(finite_fit$null_mean)
    variance_normalized <- as.numeric(finite_fit$null_variance)
    finite_observed <- stats::setNames(
      .aspu_rescale_power_values(
        observed_normalized, finite_powers * log_normalizer
      ),
      finite_labels
    )
    finite_mean <- stats::setNames(
      .aspu_rescale_power_values(
        mean_normalized, finite_powers * log_normalizer
      ),
      finite_labels
    )
    finite_variance <- stats::setNames(
      .aspu_rescale_power_values(
        variance_normalized, 2 * finite_powers * log_normalizer
      ),
      finite_labels
    )
    finite_z <- stats::setNames(
      as.numeric(finite_fit$standardized), finite_labels
    )
    finite_covariance_normalized <- as.matrix(finite_fit$null_covariance)
    covariance_exponents <- outer(finite_powers, finite_powers, "+")
    finite_covariance <- .aspu_rescale_power_values(
      finite_covariance_normalized,
      covariance_exponents * log_normalizer
    )
    finite_correlation <- as.matrix(finite_fit$null_correlation)
    dimnames(finite_covariance) <- list(finite_labels, finite_labels)
    dimnames(finite_correlation) <- list(finite_labels, finite_labels)
    odd <- finite_powers %% 2 == 1
    individual_p_values <- stats::setNames(
      ifelse(
        odd,
        2 * stats::pnorm(abs(finite_z), lower.tail = FALSE),
        stats::pnorm(finite_z, lower.tail = FALSE)
      ),
      paste0("SPU_", finite_labels)
    )
    if (any(odd)) {
      odd_names <- finite_labels[odd]
      odd_result <- .aspu_miwa_group(
        finite_z[odd], finite_correlation[odd, odd, drop = FALSE],
        "two.sided", miwa_steps, psd_adjust, psd_tol, odd_names
      )
      group_results$odd <- list(
        powers = finite_powers[odd], tail = "joint two-sided normal",
        statistic = odd_result$statistic,
        p.value = odd_result$p.value,
        correlation = odd_result$correlation
      )
      family_psd$odd <- odd_result$psd
      family_integration$odd <- odd_result$integration
    }
    if (any(!odd)) {
      even_names <- finite_labels[!odd]
      even_result <- .aspu_miwa_group(
        finite_z[!odd], finite_correlation[!odd, !odd, drop = FALSE],
        "upper", miwa_steps, psd_adjust, psd_tol, even_names
      )
      group_results$even <- list(
        powers = finite_powers[!odd], tail = "joint upper normal",
        statistic = even_result$statistic,
        p.value = even_result$p.value,
        correlation = even_result$correlation
      )
      family_psd$even <- even_result$psd
      family_integration$even <- even_result$integration
    }
    finite_components <- list(
      powers = finite_powers,
      score.scale = score_scale,
      score.definition = finite_score_definition,
      score = finite_score,
      common.normalization.scale = finite_normalizer,
      normalized.score = stats::setNames(
        finite_score_normalized, coordinate_names
      ),
      null.coordinate.covariance =
        .aspu_rescale_power_values(
          finite_null_coordinate_covariance_normalized,
          rep.int(2 * log_normalizer, p * p)
        ),
      normalized.null.coordinate.covariance =
        finite_null_coordinate_covariance_normalized,
      statistic = finite_observed,
      normalized.statistic = stats::setNames(
        observed_normalized, finite_labels
      ),
      null.mean = finite_mean,
      normalized.null.mean = stats::setNames(
        mean_normalized, finite_labels
      ),
      null.variance = finite_variance,
      normalized.null.variance = stats::setNames(
        variance_normalized, finite_labels
      ),
      standardized.statistic = finite_z,
      null.covariance = finite_covariance,
      normalized.null.covariance = finite_covariance_normalized,
      null.correlation = finite_correlation,
      p.value = individual_p_values
    )
  }

  maximum_components <- NULL
  if (has_infinite_power) {
    maximum <- max(abs(as.numeric(standardized$W)))
    g <- maximum^2 - 2 * log(p) + log(log(p))
    maximum_p_value <- .aspu_extreme_survival(g)
    individual_p_values <- c(
      individual_p_values,
      stats::setNames(maximum_p_value, "SPU_Inf")
    )
    group_results$infinite <- list(
      powers = Inf, tail = "upper extreme-value",
      statistic = g, p.value = maximum_p_value
    )
    maximum_components <- list(
      M = maximum,
      G = g,
      p.value = maximum_p_value,
      cdf = "exp(-pi^(-1/2) * exp(-g / 2))"
    )
  }

  group_p_values <- vapply(
    group_results, function(group) as.numeric(group$p.value), numeric(1L)
  )
  group_count <- length(group_p_values)
  minimum_group_p <- min(group_p_values)
  combined_p <- if (minimum_group_p >= 1) {
    1
  } else {
    -expm1(group_count * log1p(-minimum_group_p))
  }
  combined_p <- min(1, max(0, combined_p))

  components <- list(
    mean.x = mean_x,
    mean.y = mean_y,
    difference = difference,
    standard.error = standard_error,
    W = w,
    coordinate.correlation = coordinate_correlation,
    finite = finite_components,
    maximum = maximum_components,
    individual.p.value = individual_p_values,
    groups = group_results,
    group.p.value = group_p_values,
    minimum.group.p.value = minimum_group_p,
    combination.exponent = group_count
  )

  .new_hd_location_test(
    statistic = c(`minimum group p` = minimum_group_p),
    parameter = c(p = p, groups = group_count),
    p.value = combined_p,
    method = paste(
      "Xu-Lin-Wei-Pan analytical adaptive sum-of-powers two-sample test",
      if (identical(score_scale, "paper_raw")) {
        "(primary-paper raw finite-power coordinates)"
      } else {
        "(studentised-coordinate book variant)"
      },
      sprintf("(%s covariance path)", correlation_source)
    ),
    data.name = paste(x_name, "and", y_name),
    alternative = "two.sided",
    raw.statistic = c(
      stats::setNames(
        if (length(finite_powers)) as.numeric(finite_observed) else
          numeric(),
        if (length(finite_powers)) paste0("SPU_", finite_powers) else
          character()
      ),
      if (has_infinite_power) c(SPU_Inf = maximum_components$M) else numeric()
    ),
    estimate = difference,
    null.value = stats::setNames(numeric(p), coordinate_names),
    null.distribution = list(
      family = "analytical aSPU Gaussian and extreme-value limits",
      exact = FALSE,
      finite.powers = paste(
        "Isserlis Gaussian moments; odd joint two-sided and",
        "even joint upper normal tails"
      ),
      infinite.power = if (has_infinite_power) {
        "exp(-pi^(-1/2) * exp(-g / 2))"
      } else {
        NULL
      },
      combination = sprintf(
        "1 - (1 - minimum group p)^%d", group_count
      ),
      assumptions = paste(
        "Xu-Lin-Wei-Pan high-dimensional Gaussian power-sum limits,",
        "weak coordinate dependence, and consistent covariance estimation"
      )
    ),
    variance = standard_error^2,
    components = components,
    diagnostics = list(
      calibration = "analytical; no permutation or bootstrap",
      score.scale = score_scale,
      finite.power.coordinate.definition = if (
        identical(score_scale, "paper_raw")
      ) {
        "unstandardised sample-mean differences (primary paper)"
      } else {
        "inverse-variance standardised W (scale-invariant book variant)"
      },
      finite.power.common.normalization.scale = if (
        length(finite_powers)
      ) finite_normalizer else NULL,
      correlation.source = correlation_source,
      standard.error.source = as.character(
        standardized$standard_error_source
      ),
      bandwidth = if (is.null(bandwidth)) NULL else widths,
      banding.rule = if (is.null(bandwidth)) {
        "none"
      } else {
        "hard zero for abs(j - k) greater than bandwidth"
      },
      coordinate.psd = checked_coordinate$diagnostics,
      finite.family.psd = family_psd,
      integration = family_integration,
      psd.adjust = psd_adjust,
      variance.repair = "none",
      coordinate.internal.scale = stats::setNames(
        as.numeric(standardized$column_scale), coordinate_names
      ),
      combination.groups = names(group_results),
      combination.exponent = group_count,
      group.asymptotic.independence = TRUE,
      maximum.requires.p.at.least.two = has_infinite_power
    ),
    n = c(x = n1, y = n2),
    p = p,
    call = call
  )
}
