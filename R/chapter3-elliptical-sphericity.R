.sph_validate_alpha <- function(alpha) {
  alpha <- as.numeric(alpha)
  if (length(alpha) != 1L || is.na(alpha) || !is.finite(alpha) ||
      alpha <= 0 || alpha >= 1) {
    stop("`alpha` must be one finite number strictly between zero and one.",
         call. = FALSE)
  }
  alpha
}

.sph_validate_logical <- function(value, name) {
  if (!is.logical(value) || length(value) != 1L || is.na(value)) {
    stop(sprintf("`%s` must be TRUE or FALSE.", name), call. = FALSE)
  }
  value
}

.sph_variable_names <- function(x) {
  if (is.null(colnames(x))) paste0("variable", seq_len(ncol(x))) else colnames(x)
}

.sph_observation_names <- function(x) {
  if (is.null(rownames(x))) paste0("observation", seq_len(nrow(x))) else rownames(x)
}

.sph_new_test <- function(statistic, p.value, method, data.name,
                          raw.statistic, estimate = NULL,
                          components = list(), diagnostics = list(),
                          call = NULL) {
  result <- list(
    statistic = statistic,
    p.value = as.numeric(p.value),
    alternative = "greater",
    method = method,
    data.name = data.name,
    raw.statistic = raw.statistic,
    null.value = c(sphericity = 1),
    components = components,
    diagnostics = diagnostics,
    call = call
  )
  if (!is.null(estimate)) result$estimate <- estimate
  class(result) <- c("hd_sphericity_test", "htest")
  result
}

.sph_fit_spatial_center <- function(x, tol, max_iter, strict) {
  controls <- .validate_iteration_controls(tol, max_iter, zero_tol = 0)
  strict <- .sph_validate_logical(strict, "strict")
  location <- spatial_median(
    x, tol = controls$tol, max_iter = controls$max_iter,
    zero_tol = 0, warn = FALSE
  )
  diagnostics <- list(
    converged = isTRUE(attr(location, "converged")),
    iterations = as.integer(attr(location, "iterations")),
    tolerance = controls$tol,
    max.iterations.allowed = controls$max_iter,
    objective = as.numeric(attr(location, "objective")),
    relative.change = as.numeric(attr(location, "relative_change")),
    equation.residual = as.numeric(attr(location, "equation_residual")),
    strict = strict,
    algorithm = "Chapter 1 modified Weiszfeld spatial_median()"
  )
  if (!diagnostics$converged) {
    text <- sprintf(
      paste0(
        "The spatial-median iteration did not converge within %d steps ",
        "(equation residual %.6g, relative change %.6g)."
      ),
      controls$max_iter, diagnostics$equation.residual,
      diagnostics$relative.change
    )
    if (strict) stop(text, call. = FALSE)
    warning(paste0(text, " Returning the last finite iterate."),
            call. = FALSE)
  }
  location <- as.numeric(location)
  if (anyNA(location) || any(!is.finite(location))) {
    stop("The fitted spatial median is non-finite; no repair is applied.",
         call. = FALSE)
  }
  list(location = location, diagnostics = diagnostics)
}

.sph_centered_sign_core <- function(x, center, compute_second_order,
                                    keep_sscm) {
  scaled <- .center_and_scale(x, center, zero_tol = 0)
  if (compute_second_order) {
    center.scaled <- center / scaled$scale
    if (anyNA(center.scaled) || any(!is.finite(center.scaled))) {
      stop(
        paste0(
          "The origin-dependent second-order radius correction is ",
          "non-finite at the data scale; no truncation is applied."
        ),
        call. = FALSE
      )
    }
  } else {
    center.scaled <- numeric(ncol(x))
  }
  core <- cpp_elliptical_sphericity_sign_core(
    scaled$x, center.scaled, compute_second_order, keep_sscm
  )
  variable.names <- .sph_variable_names(x)
  observation.names <- .sph_observation_names(x)
  dimnames(core$signs) <- list(observation.names, variable.names)
  names(core$radii_scaled) <- observation.names
  names(core$log_radii_scaled) <- observation.names
  names(core$sscm_diagonal) <- variable.names
  if (!is.null(core$corrected_radii_scaled)) {
    names(core$corrected_radii_scaled) <- observation.names
  }
  if (!is.null(core$sscm)) {
    dimnames(core$sscm) <- list(variable.names, variable.names)
  }
  core$centering_scale <- scaled$scale
  core$centering_overflow_fallback <- scaled$overflow_fallback
  core
}

.sph_log_mean_exp <- function(value) {
  maximum <- max(value)
  maximum + log(mean(exp(value - maximum)))
}

.sph_inverse_moment_ratios <- function(log_radii) {
  if (length(log_radii) < 1L || anyNA(log_radii) ||
      any(!is.finite(log_radii))) {
    stop(
      paste0(
        "Inverse radial moments are undefined because at least one fitted ",
        "residual radius is zero or non-finite; no deletion or floor is applied."
      ),
      call. = FALSE
    )
  }
  log.mean.inverse <- .sph_log_mean_exp(-log_radii)
  log.ratio2 <- .sph_log_mean_exp(-2 * log_radii) -
    2 * log.mean.inverse
  log.ratio3 <- .sph_log_mean_exp(-3 * log_radii) -
    3 * log.mean.inverse
  ratios <- c(
    inverse.second.ratio = exp(log.ratio2),
    inverse.third.ratio = exp(log.ratio3)
  )
  if (anyNA(ratios) || any(!is.finite(ratios))) {
    stop("The inverse radial-moment ratios are non-finite; no repair is applied.",
         call. = FALSE)
  }
  list(
    ratio2 = unname(ratios[[1L]]),
    ratio3 = unname(ratios[[2L]]),
    log.mean.inverse.radius = log.mean.inverse,
    log.ratios = c(second = log.ratio2, third = log.ratio3)
  )
}

.sph_delta_from_ratios <- function(n, ratio2, ratio3) {
  n <- as.numeric(n)
  delta <- (2 - 2 * ratio2 + ratio2^2) / n^2 +
    (8 * ratio2 - 6 * ratio2^2 + 2 * ratio2 * ratio3 -
       2 * ratio3) / n^3
  if (length(delta) != 1L || is.na(delta) || !is.finite(delta)) {
    stop("The feasible bias estimate is non-finite; no repair is applied.",
         call. = FALSE)
  }
  delta
}

.sph_sign_sum_components <- function(core, n, p, bias_estimator) {
  if (bias_estimator == "normal_limit") {
    ratios <- list(
      ratio2 = NA_real_, ratio3 = NA_real_,
      log.mean.inverse.radius = NA_real_,
      log.ratios = c(second = NA_real_, third = NA_real_)
    )
    delta <- n^(-2) + 2 * n^(-3)
    radii.used <- NULL
    radius.definition <- "large-p normal-limit shortcut"
  } else {
    if (bias_estimator == "second_order") {
      corrected <- as.numeric(core$corrected_radii_scaled)
      if (length(corrected) != n || anyNA(corrected) ||
          any(!is.finite(corrected)) || any(corrected <= 0)) {
        stop(
          paste0(
            "The 2014 second-order corrected radii contain a non-positive ",
            "or non-finite value; no absolute value, deletion, or floor is applied."
          ),
          call. = FALSE
        )
      }
      log.radii <- log(corrected)
      radii.used <- corrected
      radius.definition <- paste(
        "R_i* = Rhat_i + thetahat' Uhat_i -",
        "||thetahat||^2 / (2 Rhat_i), on a common scale"
      )
    } else {
      log.radii <- as.numeric(core$log_radii_scaled)
      radii.used <- as.numeric(core$radii_scaled)
      radius.definition <- "Rhat_i = ||X_i - thetahat||, on a common scale"
    }
    ratios <- .sph_inverse_moment_ratios(log.radii)
    delta <- .sph_delta_from_ratios(n, ratios$ratio2, ratios$ratio3)
  }

  variance <- 4 * (p - 1) / (n * (n - 1) * (p + 2))
  if (!is.finite(variance) || variance <= 0) {
    stop("The null variance is not finite and positive.", call. = FALSE)
  }
  standard.error <- sqrt(variance)
  q.tilde <- as.numeric(core$q_tilde)
  z <- (q.tilde - p * delta) / standard.error
  if (!is.finite(z)) {
    stop("The standardized spatial-sign statistic is non-finite.",
         call. = FALSE)
  }
  list(
    q.tilde = q.tilde,
    ordered.sign.sum = as.numeric(core$ordered_sign_sum),
    delta = delta,
    bias = p * delta,
    variance = variance,
    standard.error = standard.error,
    z = z,
    p.value = stats::pnorm(z, lower.tail = FALSE),
    ratios = ratios,
    radii.used.scaled = radii.used,
    radius.definition = radius.definition
  )
}

.sph_sign_max_components <- function(core, n, p, alpha) {
  comparison.count <- p * (p + 1) / 2
  maximum <- as.numeric(core$max_standardized_score)
  statistic <- maximum - 2 * log(comparison.count) +
    log(log(comparison.count))
  exponential.argument <- exp(-statistic / 2) / sqrt(pi)
  p.value <- -expm1(-exponential.argument)
  critical <- -log(pi) - 2 * log(-log1p(-alpha))
  list(
    statistic = statistic,
    p.value = p.value,
    critical.value = critical,
    reject = statistic >= critical,
    comparison.count = comparison.count,
    maximum.standardized.square = maximum,
    maximum.type = if (isTRUE(core$maximum_is_diagonal)) {
      "diagonal"
    } else {
      "off-diagonal"
    },
    maximum.diagonal.score = as.numeric(core$max_diagonal_score),
    maximum.diagonal.index = as.integer(core$max_diagonal_index),
    maximum.off.diagonal.score = as.numeric(core$max_off_diagonal_score),
    maximum.off.diagonal.value = as.numeric(core$max_off_diagonal_value),
    maximum.off.diagonal.indices =
      as.integer(core$max_off_diagonal_indices),
    gumbel.cdf = exp(-exponential.argument)
  )
}

.sph_truncated_cauchy <- function(p.values) {
  p.values <- as.numeric(p.values)
  if (length(p.values) < 1L || anyNA(p.values) ||
      any(!is.finite(p.values)) || any(p.values < 0) ||
      any(p.values > 1)) {
    stop("Component p-values must be finite numbers in [0, 1].",
         call. = FALSE)
  }
  active <- p.values < 0.5
  terms <- numeric(length(p.values))
  terms[active] <- 0.5 / tan(pi * p.values[active])
  score <- sum(terms)
  if (is.nan(score) || score < 0) {
    stop("The truncated Cauchy score is invalid.", call. = FALSE)
  }
  combined <- atan2(1, score) / pi
  list(
    score = score,
    p.value = combined,
    active = active,
    terms = terms,
    empty.active.set = !any(active)
  )
}


#' Zou--Peng--Feng--Wang spatial-sign sphericity test
#'
#' Tests sphericity of an elliptical distribution with the bias-corrected
#' spatial-sign sum statistic of Zou et al. (2014). Rows of `x` are
#' observations. The raw statistic is
#' \deqn{\widetilde Q = \frac{p}{n(n-1)}\sum_{i\ne j}
#' (\widehat U_i^\mathsf{T}\widehat U_j)^2-1,}
#' where \eqn{\widehat U_i=U(X_i-\widehat\theta)} and
#' \eqn{\widehat\theta} is the spatial median. The standardized statistic is
#' \eqn{(\widetilde Q-p\widehat\delta)/\sigma_0}, with
#' \eqn{\sigma_0^2=4(p-1)/\{n(n-1)(p+2)\}}.
#'
#' `bias_estimator = "residual"` uses the translation-invariant inverse
#' residual moments adopted for the feasible sign-sum component in Zhao et al.
#' (2026), and is the practical default. `"second_order"` reproduces the
#' origin-dependent corrected-radius prescription in Zou et al. (2014), after
#' that paper sets the unknown population centre to zero. It is retained for
#' formula auditing and is not silently presented as translation invariant.
#' `"normal_limit"` uses their large-dimensional shortcut
#' \eqn{\widehat\delta=n^{-2}+2n^{-3}}. Exact zero residuals make inverse
#' moments undefined and are rejected unless the normal-limit shortcut is
#' requested. No ridge, deletion, absolute value, or numerical floor is used.
#'
#' @param x A finite numeric matrix or data frame with observations in rows;
#'   at least two rows and two columns are required.
#' @param alpha A finite test level strictly between zero and one.
#' @param bias_estimator Feasible bias calibration: `"residual"`, the literal
#'   2014 `"second_order"` prescription, or `"normal_limit"`.
#' @param tol,max_iter Convergence controls passed to [spatial_median()].
#' @param strict If `TRUE`, spatial-median nonconvergence is an error. If
#'   `FALSE`, the last finite iterate is used with a warning.
#'
#' @return An object of class `c("hd_sphericity_test", "htest")`. Raw
#'   components include the fitted centre and signs, radii, inverse-moment
#'   ratios, feasible bias, null variance, and convergence diagnostics.
#' @export
#'
#' @references
#' Zou, C., Peng, L., Feng, L., and Wang, Z. (2014). Multivariate-sign-based
#' high-dimensional tests for sphericity. *Biometrika*, 101, 229--236.
#' \doi{10.1093/biomet/ast040}.
#'
#' @examples
#' x <- rbind(c(-2, 0), c(-1, 1), c(0, -2), c(1, 2), c(3, -1), c(2, 1))
#' zou_peng_feng_wang_sphericity_test(x)
zou_peng_feng_wang_sphericity_test <- function(
    x, alpha = 0.05,
    bias_estimator = c("residual", "second_order", "normal_limit"),
    tol = 1e-8, max_iter = 1000L, strict = TRUE) {
  call <- match.call()
  data.name <- deparse1(substitute(x))
  x <- .as_data_matrix(x, "x", min_rows = 2L)
  if (ncol(x) < 2L) {
    stop("`x` must have at least two columns for sphericity testing.",
         call. = FALSE)
  }
  alpha <- .sph_validate_alpha(alpha)
  bias_estimator <- match.arg(bias_estimator)
  fitted <- .sph_fit_spatial_center(x, tol, max_iter, strict)
  names(fitted$location) <- .sph_variable_names(x)
  core <- .sph_centered_sign_core(
    x, fitted$location,
    compute_second_order = identical(bias_estimator, "second_order"),
    keep_sscm = FALSE
  )
  sum.component <- .sph_sign_sum_components(
    core, nrow(x), ncol(x), bias_estimator
  )
  critical <- stats::qnorm(1 - alpha)

  .sph_new_test(
    statistic = c(Z = sum.component$z),
    p.value = sum.component$p.value,
    method = paste(
      "Zou-Peng-Feng-Wang spatial-sign sphericity test",
      sprintf("(%s bias calibration)", bias_estimator)
    ),
    data.name = data.name,
    raw.statistic = c(Q.tilde = sum.component$q.tilde),
    estimate = fitted$location,
    components = list(
      Q.tilde = sum.component$q.tilde,
      ordered.sign.sum = sum.component$ordered.sign.sum,
      delta.hat = sum.component$delta,
      estimated.null.bias = sum.component$bias,
      sigma0.squared = sum.component$variance,
      sigma0 = sum.component$standard.error,
      inverse.moment.ratios = sum.component$ratios,
      radius.definition = sum.component$radius.definition,
      radii.used.scaled = sum.component$radii.used.scaled,
      fitted.radii.scaled = core$radii_scaled,
      second.order.corrected.radii.scaled =
        core$corrected_radii_scaled,
      fitted.signs = core$signs,
      fitted.location = fitted$location,
      critical.value = critical,
      reject = sum.component$z > critical
    ),
    diagnostics = list(
      bias.estimator = bias_estimator,
      bias.source = switch(
        bias_estimator,
        residual = paste(
          "translation-invariant residual inverse moments used by the",
          "2026 adaptive sign-sum component"
        ),
        second_order = paste(
          "literal 2014 corrected-radius audit path after setting the",
          "population centre to the coordinate origin"
        ),
        normal_limit = "2014 large-dimensional normal radial shortcut"
      ),
      translation.invariant.calibration =
        !identical(bias_estimator, "second_order"),
      zero.residuals = as.numeric(core$n_zero),
      centering.scale = core$centering_scale,
      centering.overflow.fallback = core$centering_overflow_fallback,
      second.order.extra.scale = core$second_order_scale,
      second.order.nonpositive = core$corrected_nonpositive,
      second.order.nonfinite = core$corrected_nonfinite,
      spatial.median = fitted$diagnostics,
      p.growth.condition = "primary theorem assumes p = O(n^2)",
      tail = "upper normal",
      zero.sign.convention = "U(0) = 0",
      regularization = "none",
      ridge = "none",
      numerical.floor = "none",
      deletion = "none",
      perturbation = "none"
    ),
    call = call
  )
}


#' Feng--Liu spatial-rank sphericity tests
#'
#' Implements the Spearman- and Kendall-type high-dimensional sphericity tests
#' of Feng and Liu (2017). `method = "spearman"` uses
#' \deqn{\widehat{\operatorname{tr}(\Omega^2)}=
#' \{2n(n-1)(n-2)(n-3)\}^{-1}
#' \sum^*(U_{ij}^\mathsf{T}U_{kl})(U_{kj}^\mathsf{T}U_{il})}
#' and \eqn{\widetilde Q=4p\widehat{\operatorname{tr}(\Omega^2)}-1}.
#' `method = "kendall"` removes the factor two, replaces the summand by
#' \eqn{(U_{ij}^\mathsf{T}U_{kl})^2}, and uses
#' \eqn{\widetilde Q=p\widehat{\operatorname{tr}(\Omega^2)}-1}. The star
#' means that all four ordered indices are distinct. Both divide by the same
#' \eqn{\sigma_0} used by the spatial-sign test.
#'
#' Pair directions are computed once in C++. For each unordered four-set, its
#' 24 ordered terms are reduced exactly to the three disjoint pairings. A tied
#' pair has the literal direction `U(0) = 0`, is counted in diagnostics, and is
#' never dropped from a denominator.
#'
#' @inheritParams zou_peng_feng_wang_sphericity_test
#' @param method Either `"spearman"` or `"kendall"`.
#' @param keep_pair_signs Whether to retain the
#'   \eqn{\binom{n}{2}\times p} pair-sign matrix and its endpoint table.
#'
#' @return An object of class `c("hd_sphericity_test", "htest")` containing
#'   the exact ordered sum, denominator, trace estimate, pair-tie diagnostics,
#'   and optional pair signs.
#' @export
#'
#' @references
#' Feng, L. and Liu, B. (2017). High-dimensional rank tests for sphericity.
#' *Journal of Multivariate Analysis*, 155, 217--233.
#' \doi{10.1016/j.jmva.2017.01.005}.
#'
#' @examples
#' x <- rbind(c(-2, 0), c(-1, 1), c(0, -2), c(1, 2), c(3, -1), c(2, 1))
#' feng_liu_rank_sphericity_test(x, method = "kendall")
feng_liu_rank_sphericity_test <- function(
    x, method = c("spearman", "kendall"), alpha = 0.05,
    keep_pair_signs = FALSE) {
  call <- match.call()
  data.name <- deparse1(substitute(x))
  x <- .as_data_matrix(x, "x", min_rows = 4L)
  if (ncol(x) < 2L) {
    stop("`x` must have at least two columns for sphericity testing.",
         call. = FALSE)
  }
  method <- match.arg(method)
  alpha <- .sph_validate_alpha(alpha)
  keep_pair_signs <- .sph_validate_logical(
    keep_pair_signs, "keep_pair_signs"
  )
  core <- cpp_elliptical_sphericity_rank_core(
    x, method, keep_pair_signs
  )
  n <- nrow(x)
  p <- ncol(x)
  variance <- 4 * (p - 1) / (n * (n - 1) * (p + 2))
  sigma0 <- sqrt(variance)
  z <- as.numeric(core$q_tilde) / sigma0
  p.value <- stats::pnorm(z, lower.tail = FALSE)
  critical <- stats::qnorm(1 - alpha)
  variable.names <- .sph_variable_names(x)
  if (!is.null(core$pair_signs)) {
    colnames(core$pair_signs) <- variable.names
    colnames(core$pair_endpoints) <- c("first", "second")
  }

  .sph_new_test(
    statistic = c(Z = z),
    p.value = p.value,
    method = sprintf(
      "Feng-Liu %s-type spatial-rank sphericity test", method
    ),
    data.name = data.name,
    raw.statistic = c(Q.tilde = as.numeric(core$q_tilde)),
    components = list(
      Q.tilde = as.numeric(core$q_tilde),
      trace.estimate = as.numeric(core$trace_estimate),
      ordered.sum = as.numeric(core$ordered_sum),
      ordered.denominator = as.numeric(core$ordered_denominator),
      trace.denominator = if (method == "spearman") {
        2 * as.numeric(core$ordered_denominator)
      } else {
        as.numeric(core$ordered_denominator)
      },
      sigma0.squared = variance,
      sigma0 = sigma0,
      pair.signs = core$pair_signs,
      pair.endpoints = core$pair_endpoints,
      critical.value = critical,
      reject = z > critical
    ),
    diagnostics = list(
      rank.method = method,
      ordered.indices = "all four indices distinct",
      quadruple.reduction = core$quadruple_reduction,
      pair.count = as.numeric(core$n_pairs),
      zero.pair.directions = as.numeric(core$n_zero_pairs),
      zero.sign.convention = "U(0) = 0; tied pairs remain in denominators",
      p.growth.condition = "primary theorem assumes p = O(n^2)",
      tail = "upper normal",
      regularization = "none",
      ridge = "none",
      numerical.floor = "none",
      deletion = "none",
      perturbation = "none"
    ),
    call = call
  )
}


#' Zhao--Yang--Zhang--Feng--Wang spatial-sign max sphericity test
#'
#' Computes the sparse-alternative max statistic from the fitted spatial-sign
#' covariance matrix \eqn{\widehat\Omega=(\widehat\psi_{ij})}. The diagonal
#' and off-diagonal entries are standardized with their distinct null
#' variances, maximized, and centered by
#' \eqn{-2\log\{p(p+1)/2\}+\log\log\{p(p+1)/2\}}. The upper-tail limiting cdf
#' is \eqn{G(t)=\exp\{-\pi^{-1/2}\exp(-t/2)\}}.
#'
#' Exact zero residuals retain `U(0) = 0`, so the empirical SSCM can have trace
#' below one; the count and trace are returned. `keep_sscm = FALSE` avoids
#' allocating a \eqn{p\times p} return matrix while still computing the exact
#' maximum. No covariance regularization or numerical repair is performed.
#'
#' @inheritParams zou_peng_feng_wang_sphericity_test
#' @param keep_sscm Whether to retain the full fitted spatial-sign covariance
#'   matrix. Its diagonal and maximizing entries are always returned.
#'
#' @return An object of class `c("hd_sphericity_test", "htest")` with the
#'   Gumbel statistic, upper-tail p-value, fitted signs, SSCM diagnostics, and
#'   exact maximizing entry.
#' @export
#'
#' @references
#' Zhao, P., Yang, F., Zhang, X., Feng, L., and Wang, Z. (2026). Adaptive tests
#' for high-dimensional sphericity under different distribution types.
#' *Journal of Multivariate Analysis*, 214, 105634.
#' \doi{10.1016/j.jmva.2026.105634}.
#'
#' @examples
#' x <- rbind(c(-2, 0), c(-1, 1), c(0, -2), c(1, 2), c(3, -1), c(2, 1))
#' zhao_yang_zhang_feng_wang_sign_max_test(x)
zhao_yang_zhang_feng_wang_sign_max_test <- function(
    x, alpha = 0.05, tol = 1e-8, max_iter = 1000L,
    strict = TRUE, keep_sscm = FALSE) {
  call <- match.call()
  data.name <- deparse1(substitute(x))
  x <- .as_data_matrix(x, "x", min_rows = 2L)
  if (ncol(x) < 2L) {
    stop("`x` must have at least two columns for sphericity testing.",
         call. = FALSE)
  }
  alpha <- .sph_validate_alpha(alpha)
  keep_sscm <- .sph_validate_logical(keep_sscm, "keep_sscm")
  fitted <- .sph_fit_spatial_center(x, tol, max_iter, strict)
  names(fitted$location) <- .sph_variable_names(x)
  core <- .sph_centered_sign_core(
    x, fitted$location, compute_second_order = FALSE,
    keep_sscm = keep_sscm
  )
  maximum <- .sph_sign_max_components(core, nrow(x), ncol(x), alpha)

  .sph_new_test(
    statistic = c(T.SM = maximum$statistic),
    p.value = maximum$p.value,
    method = paste(
      "Zhao-Yang-Zhang-Feng-Wang spatial-sign max sphericity test",
      "(Gumbel calibration)"
    ),
    data.name = data.name,
    raw.statistic = c(maximum.standardized.square =
                        maximum$maximum.standardized.square),
    estimate = fitted$location,
    components = list(
      T.SM = maximum$statistic,
      p.SM = maximum$p.value,
      maximum.standardized.square =
        maximum$maximum.standardized.square,
      maximum.type = maximum$maximum.type,
      maximum.diagonal.score = maximum$maximum.diagonal.score,
      maximum.diagonal.index = maximum$maximum.diagonal.index,
      maximum.off.diagonal.score = maximum$maximum.off.diagonal.score,
      maximum.off.diagonal.value = maximum$maximum.off.diagonal.value,
      maximum.off.diagonal.indices =
        maximum$maximum.off.diagonal.indices,
      comparison.count = maximum$comparison.count,
      gumbel.cdf = maximum$gumbel.cdf,
      critical.value = maximum$critical.value,
      reject = maximum$reject,
      fitted.location = fitted$location,
      fitted.signs = core$signs,
      sscm.diagonal = core$sscm_diagonal,
      sscm = core$sscm
    ),
    diagnostics = list(
      zero.residuals = as.numeric(core$n_zero),
      sscm.trace = as.numeric(core$sscm_trace),
      spatial.median = fitted$diagnostics,
      centering.scale = core$centering_scale,
      centering.overflow.fallback = core$centering_overflow_fallback,
      tail = "upper Gumbel; primary critical rule uses T.SM >= q_alpha",
      asymptotic.condition = "primary theorem assumes log^5(p^2 n) / n -> 0",
      zero.sign.convention = "U(0) = 0",
      regularization = "none",
      ridge = "none",
      numerical.floor = "none",
      deletion = "none",
      perturbation = "none"
    ),
    call = call
  )
}


#' Adaptive dense--sparse elliptical sphericity test
#'
#' Combines the feasible spatial-sign sum p-value (`p.SS`) with the
#' spatial-sign max p-value (`p.SM`) using the truncated Cauchy rule of Zhao et
#' al. (2026). The spatial median and fitted signs are computed once. The
#' returned `Cauchy.score` is the inner non-negative score,
#' \deqn{\frac12\tan\{\pi(1/2-p_{SS})\}I(p_{SS}<1/2)+
#'       \frac12\tan\{\pi(1/2-p_{SM})\}I(p_{SM}<1/2),}
#' while `p.value` is the distinct final quantity
#' \eqn{1-F_C(\mathrm{Cauchy.score})}. If neither component p-value is below
#' one half, the score is zero and the combined p-value is exactly one half.
#'
#' This function uses the translation-invariant residual inverse-moment
#' estimate of the spatial-sign sum bias specified by the adaptive-test
#' primary source. An exact fitted zero residual makes that feasible bias
#' undefined and causes an explicit error. No floor, deletion, ridge, or
#' perturbation is applied.
#'
#' @inheritParams zhao_yang_zhang_feng_wang_sign_max_test
#'
#' @return An object of class `c("hd_sphericity_test", "htest")`. The scalar
#'   `statistic` is the Cauchy score and `p.value` is the final combined
#'   probability. `components$sum` and `components$max` retain both complete
#'   component calibrations.
#' @export
#'
#' @references
#' Zhao, P., Yang, F., Zhang, X., Feng, L., and Wang, Z. (2026). Adaptive tests
#' for high-dimensional sphericity under different distribution types.
#' *Journal of Multivariate Analysis*, 214, 105634.
#' \doi{10.1016/j.jmva.2026.105634}.
#'
#' @examples
#' x <- rbind(c(-2, 0), c(-1, 1), c(0, -2), c(1, 2), c(3, -1), c(2, 1))
#' zhao_yang_zhang_feng_wang_adaptive_sphericity_test(x)
zhao_yang_zhang_feng_wang_adaptive_sphericity_test <- function(
    x, alpha = 0.05, tol = 1e-8, max_iter = 1000L,
    strict = TRUE, keep_sscm = FALSE) {
  call <- match.call()
  data.name <- deparse1(substitute(x))
  x <- .as_data_matrix(x, "x", min_rows = 2L)
  if (ncol(x) < 2L) {
    stop("`x` must have at least two columns for sphericity testing.",
         call. = FALSE)
  }
  alpha <- .sph_validate_alpha(alpha)
  keep_sscm <- .sph_validate_logical(keep_sscm, "keep_sscm")
  fitted <- .sph_fit_spatial_center(x, tol, max_iter, strict)
  names(fitted$location) <- .sph_variable_names(x)
  core <- .sph_centered_sign_core(
    x, fitted$location, compute_second_order = FALSE,
    keep_sscm = keep_sscm
  )
  sum.component <- .sph_sign_sum_components(
    core, nrow(x), ncol(x), "residual"
  )
  max.component <- .sph_sign_max_components(
    core, nrow(x), ncol(x), alpha
  )
  combination <- .sph_truncated_cauchy(
    c(p.SS = sum.component$p.value, p.SM = max.component$p.value)
  )

  .sph_new_test(
    statistic = c(Cauchy.score = combination$score),
    p.value = combination$p.value,
    method = paste(
      "Zhao-Yang-Zhang-Feng-Wang adaptive elliptical sphericity test",
      "(spatial-sign sum/max truncated-Cauchy combination)"
    ),
    data.name = data.name,
    raw.statistic = c(Cauchy.score = combination$score),
    estimate = fitted$location,
    components = list(
      Cauchy.score = combination$score,
      combined.p.value = combination$p.value,
      active.components = stats::setNames(
        combination$active, c("p.SS", "p.SM")
      ),
      Cauchy.terms = stats::setNames(
        combination$terms, c("p.SS", "p.SM")
      ),
      empty.active.set = combination$empty.active.set,
      reject = combination$p.value <= alpha,
      sum = list(
        T.SS = sum.component$z,
        p.SS = sum.component$p.value,
        Q.tilde = sum.component$q.tilde,
        delta.hat = sum.component$delta,
        estimated.null.bias = sum.component$bias,
        sigma0.squared = sum.component$variance,
        sigma0 = sum.component$standard.error,
        inverse.moment.ratios = sum.component$ratios,
        fitted.radii.scaled = core$radii_scaled
      ),
      max = list(
        T.SM = max.component$statistic,
        p.SM = max.component$p.value,
        maximum.standardized.square =
          max.component$maximum.standardized.square,
        maximum.type = max.component$maximum.type,
        maximum.diagonal.score = max.component$maximum.diagonal.score,
        maximum.diagonal.index = max.component$maximum.diagonal.index,
        maximum.off.diagonal.score =
          max.component$maximum.off.diagonal.score,
        maximum.off.diagonal.value =
          max.component$maximum.off.diagonal.value,
        maximum.off.diagonal.indices =
          max.component$maximum.off.diagonal.indices,
        comparison.count = max.component$comparison.count,
        gumbel.cdf = max.component$gumbel.cdf,
        critical.value = max.component$critical.value,
        reject = max.component$reject
      ),
      fitted.location = fitted$location,
      fitted.signs = core$signs,
      sscm.diagonal = core$sscm_diagonal,
      sscm = core$sscm
    ),
    diagnostics = list(
      sum.bias.estimator = "translation-invariant residual inverse moments",
      cauchy.output.contract = paste(
        "statistic is the inner score; p.value is 1 minus the standard",
        "Cauchy cdf of that score"
      ),
      cauchy.truncation = "a component contributes only when p < 0.5",
      empty.active.set.rule = "score = 0 and combined p-value = 0.5",
      zero.residuals = as.numeric(core$n_zero),
      sscm.trace = as.numeric(core$sscm_trace),
      spatial.median = fitted$diagnostics,
      centering.scale = core$centering_scale,
      centering.overflow.fallback = core$centering_overflow_fallback,
      asymptotic.condition = paste(
        "primary joint null result assumes p = O(n^2) in addition to",
        "the component regularity conditions"
      ),
      zero.sign.convention = "U(0) = 0",
      regularization = "none",
      ridge = "none",
      numerical.floor = "none",
      deletion = "none",
      perturbation = "none"
    ),
    call = call
  )
}
