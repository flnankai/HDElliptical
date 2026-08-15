.erht_alpha <- function(alpha) {
  alpha <- as.numeric(alpha)
  if (length(alpha) != 1L || is.na(alpha) || !is.finite(alpha) ||
      alpha <= 0 || alpha >= 1) {
    stop("`alpha` must be one finite number strictly between zero and one.",
         call. = FALSE)
  }
  alpha
}

.erht_strict <- function(strict) {
  if (!is.logical(strict) || length(strict) != 1L || is.na(strict)) {
    stop("`strict` must be `TRUE` or `FALSE`.", call. = FALSE)
  }
  strict
}

.erht_keep_companion <- function(value) {
  if (!is.logical(value) || length(value) != 1L || is.na(value)) {
    stop("`keep_companion` must be `TRUE` or `FALSE`.", call. = FALSE)
  }
  value
}

.erht_rho <- function(rho, multiple = TRUE) {
  rho <- as.numeric(rho)
  if (!length(rho) || anyNA(rho) || any(!is.finite(rho)) || any(rho <= 0)) {
    stop("`rho` must contain finite, strictly positive ridge values.",
         call. = FALSE)
  }
  if (!multiple && length(rho) != 1L) {
    stop("The fixed-ridge ERHT function requires exactly one `rho` value.",
         call. = FALSE)
  }
  if (multiple && length(rho) < 2L) {
    stop("ERHT--CC requires at least two ridge values.", call. = FALSE)
  }
  if (length(rho) > 1L && any(diff(rho) <= 0)) {
    stop("`rho` must be strictly increasing with no duplicates.",
         call. = FALSE)
  }
  rho
}

.erht_rho_labels <- function(rho) {
  paste0("rho=", sprintf("%.17g", rho))
}

.erht_scale_value <- function(value, scale, power) {
  value <- as.numeric(value)
  if (scale == 1 || power == 0) return(value)
  answer <- numeric(length(value))
  zero <- value == 0
  answer[zero] <- 0
  use <- !zero & !is.na(value)
  if (any(use)) {
    log.absolute <- log(abs(value[use])) + power * log(scale)
    answer[use] <- sign(value[use]) * exp(log.absolute)
  }
  answer[is.na(value)] <- NA_real_
  answer
}

.erht_two_scale_z <- function(raw.null.scaled, center.residual.scaled,
                              variance.residual.scaled,
                              null.scale, residual.scale) {
  raw.null.scaled <- as.numeric(raw.null.scaled)
  center.residual.scaled <- as.numeric(center.residual.scaled)
  variance.residual.scaled <- as.numeric(variance.residual.scaled)
  if (length(raw.null.scaled) != length(center.residual.scaled) ||
      length(raw.null.scaled) != length(variance.residual.scaled) ||
      anyNA(raw.null.scaled) || any(raw.null.scaled < 0) ||
      anyNA(center.residual.scaled) || any(center.residual.scaled < 0) ||
      anyNA(variance.residual.scaled) ||
      any(!is.finite(variance.residual.scaled)) ||
      any(variance.residual.scaled <= 0)) {
    stop("Invalid internally scaled ERHT components.", call. = FALSE)
  }

  log.multiplier <- 2 * (log(null.scale) - log(residual.scale))
  answer <- numeric(length(raw.null.scaled))
  for (i in seq_along(answer)) {
    log.quadratic <- if (raw.null.scaled[i] > 0) {
      log(raw.null.scaled[i]) + log.multiplier
    } else {
      -Inf
    }
    log.center <- if (center.residual.scaled[i] > 0) {
      log(center.residual.scaled[i])
    } else {
      -Inf
    }
    log.denominator <- 0.5 * log(variance.residual.scaled[i])

    if (is.infinite(log.quadratic) && log.quadratic < 0 &&
        is.infinite(log.center) && log.center < 0) {
      answer[i] <- 0
      next
    }
    if (identical(log.quadratic, log.center)) {
      answer[i] <- 0
      next
    }

    if (log.quadratic > log.center) {
      sign.value <- 1
      log.maximum <- log.quadratic
      log.factor <- log(-expm1(log.center - log.quadratic))
    } else {
      sign.value <- -1
      log.maximum <- log.center
      log.factor <- log(-expm1(log.quadratic - log.center))
    }
    log.absolute.z <- log.maximum + log.factor - log.denominator
    answer[i] <- if (log.absolute.z > log(.Machine$double.xmax)) {
      sign.value * Inf
    } else {
      sign.value * exp(log.absolute.z)
    }
  }
  answer
}

.erht_tail <- function(z) {
  list(
    p.value = stats::pnorm(z, lower.tail = FALSE),
    log.p.value = stats::pnorm(z, lower.tail = FALSE, log.p = TRUE),
    log.one.minus.p.value = stats::pnorm(z, lower.tail = TRUE, log.p = TRUE)
  )
}

.erht_cot_component <- function(p.value, log.p.value,
                                log.one.minus.p.value, weight) {
  if (p.value == 0.5) {
    return(list(sign = 0, log.abs = -Inf))
  }
  if (p.value > 1e-8 && p.value < 1 - 1e-8) {
    transformed <- tan(pi * (0.5 - p.value))
    return(list(
      sign = sign(transformed),
      log.abs = log(weight) + log(abs(transformed))
    ))
  }
  if (p.value < 0.5) {
    return(list(
      sign = 1,
      log.abs = log(weight) - log(pi) - log.p.value
    ))
  }
  list(
    sign = -1,
    log.abs = log(weight) - log(pi) - log.one.minus.p.value
  )
}

.erht_signed_log_add <- function(first, second) {
  if (first$sign == 0) return(second)
  if (second$sign == 0) return(first)
  if (first$sign == second$sign) {
    maximum <- max(first$log.abs, second$log.abs)
    minimum <- min(first$log.abs, second$log.abs)
    if (is.infinite(maximum) && maximum > 0) {
      return(list(sign = first$sign, log.abs = Inf))
    }
    return(list(
      sign = first$sign,
      log.abs = maximum + log1p(exp(minimum - maximum))
    ))
  }
  if (is.infinite(first$log.abs) && first$log.abs > 0 &&
      is.infinite(second$log.abs) && second$log.abs > 0) {
    stop(
      paste0(
        "The ERHT Cauchy combination is indeterminate because exact ",
        "opposite infinite marginal transforms are present."
      ),
      call. = FALSE
    )
  }
  if (identical(first$log.abs, second$log.abs)) {
    return(list(sign = 0, log.abs = -Inf))
  }
  larger <- if (first$log.abs > second$log.abs) first else second
  smaller <- if (first$log.abs > second$log.abs) second else first
  ratio <- exp(smaller$log.abs - larger$log.abs)
  if (ratio == 1) return(list(sign = 0, log.abs = -Inf))
  list(
    sign = larger$sign,
    log.abs = larger$log.abs + log1p(-ratio)
  )
}

.erht_cauchy_grid <- function(tails, weights) {
  combined <- list(sign = 0, log.abs = -Inf)
  for (i in seq_along(weights)) {
    combined <- .erht_signed_log_add(
      combined,
      .erht_cot_component(
        tails$p.value[i], tails$log.p.value[i],
        tails$log.one.minus.p.value[i], weights[i]
      )
    )
  }
  if (combined$sign == 0) {
    return(list(
      statistic = 0, statistic.sign = 0,
      log.absolute.statistic = -Inf, angle = 0, p.value = 0.5,
      log.p.value = log(0.5), log.one.minus.p.value = log(0.5)
    ))
  }
  inverse.magnitude <- if (-combined$log.abs >
      log(.Machine$double.xmax)) Inf else exp(-combined$log.abs)
  complement.angle <- atan(inverse.magnitude)
  angle <- combined$sign * (pi / 2 - complement.angle)
  if (combined$sign > 0) {
    p.value <- complement.angle / pi
    log.p.value <- if (p.value > 0) {
      log(p.value)
    } else {
      -combined$log.abs - log(pi)
    }
    log.one.minus.p.value <- if (p.value < 1) log1p(-p.value) else -Inf
  } else {
    one.minus.p <- complement.angle / pi
    p.value <- 1 - one.minus.p
    log.one.minus.p.value <- if (one.minus.p > 0) {
      log(one.minus.p)
    } else {
      -combined$log.abs - log(pi)
    }
    log.p.value <- if (p.value > 0) log(p.value) else -Inf
  }
  if (!is.finite(p.value) || p.value < 0 || p.value > 1) {
    stop("Stable ERHT Cauchy aggregation produced an invalid probability.",
         call. = FALSE)
  }
  statistic <- if (combined$log.abs > log(.Machine$double.xmax)) {
    combined$sign * Inf
  } else {
    combined$sign * exp(combined$log.abs)
  }
  list(
    statistic = statistic,
    statistic.sign = combined$sign,
    log.absolute.statistic = combined$log.abs,
    angle = angle,
    p.value = p.value,
    log.p.value = log.p.value,
    log.one.minus.p.value = log.one.minus.p.value
  )
}

.erht_fit_grid <- function(x, mu, rho, tol, max_iter, strict,
                           keep_companion) {
  x <- .as_data_matrix(x, "x", min_rows = 3L)
  if (ncol(x) < 2L) {
    stop("ERHT requires at least two variables.", call. = FALSE)
  }
  p <- ncol(x)
  mu <- if (is.null(mu)) numeric(p) else .as_location(mu, p, "mu")
  controls <- .validate_iteration_controls(tol, max_iter, 0)
  strict <- .erht_strict(strict)
  keep_companion <- .erht_keep_companion(keep_companion)

  location <- spatial_median(
    x, tol = controls$tol, max_iter = controls$max_iter,
    zero_tol = 0, warn = FALSE
  )
  median.diagnostics <- list(
    objective = attr(location, "objective"),
    iterations = attr(location, "iterations"),
    converged = isTRUE(attr(location, "converged")),
    relative.change = attr(location, "relative_change"),
    equation.residual = attr(location, "equation_residual")
  )
  if (!median.diagnostics$converged) {
    message <- paste0(
      "The ERHT sample spatial median did not satisfy its convergence ",
      "criterion in `max_iter` iterations."
    )
    if (strict) stop(message, call. = FALSE) else warning(message, call. = FALSE)
  }

  residual.scale <- .center_and_scale(x, as.numeric(location), 0)
  null.scale <- .center_and_scale(
    matrix(as.numeric(location), nrow = 1L), mu, 0
  )
  # The inverse-radius calibration and the null displacement have different
  # homogeneities.  Scaling them together can make fourth powers of the
  # inverse radii overflow when the null is far from an otherwise ordinary
  # sample.  Normalize each part by its own positive scale and combine the
  # two scales only in the final standardized statistic.
  residuals.scaled <- residual.scale$x
  null.difference.scaled <- as.numeric(null.scale$x)
  native <- cpp_erht_grid(
    residuals.scaled, null.difference.scaled, rho, keep_companion
  )

  native$z <- .erht_two_scale_z(
    native$raw_statistic_scaled,
    native$center_scaled,
    native$variance_scaled,
    null.scale$scale,
    residual.scale$scale
  )

  psi.exponents <- c(0, 1, 2, 2, 3, 4)
  psi.raw <- native$psi_scaled
  for (j in seq_along(psi.exponents)) {
    psi.raw[, j] <- .erht_scale_value(
      psi.raw[, j], residual.scale$scale, -psi.exponents[j]
    )
  }
  g.exponents <- c(2, 3, 4)
  g.raw <- native$g_scaled
  for (j in seq_along(g.exponents)) {
    g.raw[, j] <- .erht_scale_value(
      g.raw[, j], residual.scale$scale, g.exponents[j]
    )
  }

  variable.names <- .hotelling_variable_names(x)
  names(location) <- names(mu) <- variable.names
  row.labels <- rownames(x)
  if (is.null(row.labels)) row.labels <- paste0("x", seq_len(nrow(x)))
  rho.labels <- .erht_rho_labels(rho)
  names(native$radii_scaled) <- names(native$weights_scaled) <- row.labels
  colnames(native$psi_scaled) <- colnames(psi.raw) <-
    c("psi00", "psi01", "psi02", "psi11", "psi12", "psi22")
  colnames(native$g_scaled) <- colnames(g.raw) <- c("g1", "g2", "g3")
  rownames(native$psi_scaled) <- rownames(psi.raw) <-
    rownames(native$g_scaled) <- rownames(g.raw) <- rho.labels
  rownames(native$companion_eigenvalues) <- rho.labels
  if (!is.null(native$companion_matrices)) {
    names(native$companion_matrices) <- rho.labels
    native$companion_matrices <- lapply(native$companion_matrices, function(a) {
      dimnames(a) <- list(row.labels, row.labels)
      a
    })
  }

  list(
    x = x,
    mu = mu,
    location = location,
    rho = rho,
    native = native,
    scale = residual.scale$scale,
    null.scale = null.scale$scale,
    raw.statistic = .erht_scale_value(
      native$raw_statistic_scaled, null.scale$scale, 2
    ),
    center = .erht_scale_value(
      native$center_scaled, residual.scale$scale, 2
    ),
    variance = .erht_scale_value(
      native$variance_scaled, residual.scale$scale, 4
    ),
    D = .erht_scale_value(native$D_scaled, residual.scale$scale, -2),
    mu.functional = .erht_scale_value(
      native$mu_scaled, residual.scale$scale, 2
    ),
    sigma.D.squared = .erht_scale_value(
      native$sigma_D_squared_scaled, residual.scale$scale, 4
    ),
    e = .erht_scale_value(native$e_scaled, residual.scale$scale, -1),
    t = .erht_scale_value(native$t_scaled, residual.scale$scale, -2),
    b1 = .erht_scale_value(native$b1_scaled, residual.scale$scale, -1),
    b2 = .erht_scale_value(native$b2_scaled, residual.scale$scale, -2),
    radii = .erht_scale_value(
      native$radii_scaled, residual.scale$scale, 1
    ),
    weights = .erht_scale_value(
      native$weights_scaled, residual.scale$scale, -1
    ),
    psi = psi.raw,
    g = g.raw,
    median.diagnostics = median.diagnostics,
    overflow.fallback = residual.scale$overflow_fallback ||
      null.scale$overflow_fallback
  )
}

.erht_components <- function(fit, tails, alpha) {
  native <- fit$native
  marginal <- data.frame(
    rho = fit$rho,
    T = fit$raw.statistic,
    center = fit$center,
    variance = fit$variance,
    Z = as.numeric(native$z),
    p.value = tails$p.value,
    log.p.value = tails$log.p.value,
    critical.Z = stats::qnorm(1 - alpha),
    reject = as.numeric(native$z) > stats::qnorm(1 - alpha),
    row.names = .erht_rho_labels(fit$rho),
    check.names = FALSE
  )
  list(
    location = fit$location,
    null.location = fit$mu,
    location.difference = fit$location - fit$mu,
    marginal = marginal,
    e = fit$e,
    t = fit$t,
    kappa = stats::setNames(
      as.numeric(native$kappa), .erht_rho_labels(fit$rho)
    ),
    b1 = fit$b1,
    b2 = fit$b2,
    D = fit$D,
    mu.functional = fit$mu.functional,
    sigma.D.squared = fit$sigma.D.squared,
    psi = fit$psi,
    g = fit$g,
    radii = fit$radii,
    inverse.distances = fit$weights,
    scaled = list(
      T = native$raw_statistic_scaled,
      center = native$center_scaled,
      variance = native$variance_scaled,
      e = native$e_scaled,
      t = native$t_scaled,
      b1 = native$b1_scaled,
      b2 = native$b2_scaled,
      D = native$D_scaled,
      mu.functional = native$mu_scaled,
      sigma.D.squared = native$sigma_D_squared_scaled,
      psi = native$psi_scaled,
      g = native$g_scaled,
      radii = native$radii_scaled,
      inverse.distances = native$weights_scaled,
      residual.scale = fit$scale,
      null.displacement.scale = fit$null.scale
    ),
    companion = list(
      gram.eigenvalues = native$gram_eigenvalues,
      eigenvalues = native$companion_eigenvalues,
      matrices = native$companion_matrices
    )
  )
}

#' Elliptical regularized Hotelling fixed-ridge test
#'
#' Tests a high-dimensional location vector under elliptical symmetry and
#' pervasive dependence using the fixed-ridge ERHT statistic of Feng, Zhou,
#' and Wang (2026). Observations are rows and variables are columns.
#'
#' Let \eqn{\widehat\theta} be the sample spatial median,
#' \eqn{\widehat Y_i=\sqrt p\,U(X_i-\widehat\theta)}, and
#' \deqn{\widehat R_n=n^{-1}\sum_i\widehat Y_i\widehat Y_i^{\mathsf T}.}
#' For \eqn{\rho>0}, the raw statistic is
#' \deqn{T_n(\rho)=n(\widehat\theta-\theta_0)^{\mathsf T}
#' (\widehat R_n+\rho I)^{-1}(\widehat\theta-\theta_0).}
#' The function implements the paper's direct feasible companion-matrix
#' centering \eqn{n\widehat\mu_n(\rho)} and variance
#' \eqn{n\widehat\sigma_{D,n}^2(\rho)}, returning the upper-tail normal
#' calibration of
#' \deqn{Z_n(\rho)=\{T_n(\rho)-n\widehat\mu_n(\rho)\}/
#' \{n\widehat\sigma_{D,n}^2(\rho)\}^{1/2}.}
#'
#' The C++ kernel uses one economy SVD of the sign matrix. It evaluates the
#' ridge quadratic in row and null spaces and switches between the companion
#' matrix \eqn{A} and its complement \eqn{I-A} for the feasible calibration.
#' Thus it never forms a persistent \eqn{p\times p} inverse or subtracts
#' nearly equal Woodbury, diagonal-weight, or companion terms. A
#' zero fitted residual makes the required inverse distance undefined and is
#' reported as an error; no observation is omitted or perturbed and no
#' additional ridge, pseudoinverse, absolute-value repair, or variance floor
#' is used.
#'
#' The paper's numerical section mentions a separate Bartlett center
#' correction with \eqn{a=8n/p}, but neither the article nor its public arXiv
#' source specifies the formula that maps \eqn{a} into the feasible center.
#' This function therefore implements the fully stated equations only and
#' deliberately does not guess that simulation-program modification.
#'
#' @param x Numeric matrix or data frame with observations in rows. At least
#'   three observations and two variables are required.
#' @param mu Null location vector. `NULL` means the zero vector.
#' @param rho One finite, strictly positive ridge value.
#' @param alpha Significance level for the upper-tail rejection rule.
#' @param tol,max_iter Convergence controls for the sample spatial median.
#' @param strict If `TRUE`, fail when the spatial median does not converge;
#'   otherwise warn and continue only when all ERHT quantities remain defined.
#' @param keep_companion If `TRUE`, retain the full \eqn{n\times n} companion
#'   matrix. The default retains only its eigenvalues and scalar functionals.
#'
#' @return An object of class `c("hd_location_test", "htest")`. The reported
#'   statistic is \eqn{Z_n(\rho)}; `raw.statistic` is \eqn{T_n(\rho)}.
#'   Components include all feasible functionals, raw and internally scaled
#'   values, inverse distances, spatial-median diagnostics, and companion
#'   eigenvalues. Raw dimensional fields can underflow to zero or overflow to
#'   infinity when their physical units exceed the double-precision range;
#'   the separately normalized fields in `components$scaled` remain the
#'   numerical-audit representation used to compute the reported statistic.
#'
#' @references
#' Feng, L., Zhou, L., and Wang, X. (2026). Elliptical regularized Hotelling
#' testing for high dimensional data. arXiv:2606.25942.
#' \doi{10.48550/arXiv.2606.25942}.
#'
#' @examples
#' set.seed(1)
#' x <- matrix(rnorm(120), 30, 4)
#' elliptical_regularized_hotelling_test(x, rho = 0.5)
#'
#' @export
elliptical_regularized_hotelling_test <- function(
    x, mu = NULL, rho = 0.5, alpha = 0.05,
    tol = 1e-8, max_iter = 1000L, strict = TRUE,
    keep_companion = FALSE) {
  call <- match.call()
  data.name <- deparse1(substitute(x))
  rho <- .erht_rho(rho, multiple = FALSE)
  alpha <- .erht_alpha(alpha)
  fit <- .erht_fit_grid(
    x, mu, rho, tol, max_iter, strict, keep_companion
  )
  z <- as.numeric(fit$native$z[1L])
  tails <- .erht_tail(z)
  components <- .erht_components(fit, tails, alpha)

  .new_hd_location_test(
    statistic = c(Z.ERHT = z),
    parameter = c(rho = rho),
    p.value = tails$p.value,
    method = "Feng-Zhou-Wang elliptical regularized Hotelling test",
    data.name = data.name,
    alternative = "two.sided",
    raw.statistic = c(T.ERHT = fit$raw.statistic[1L]),
    estimate = fit$location,
    null.value = fit$mu,
    null.distribution = list(
      family = "asymptotic standard normal",
      exact = FALSE,
      tail = "upper",
      assumptions = paste(
        "elliptical symmetry, positive radii with controlled inverse",
        "moments, proportional p/n growth, and a fixed pervasive-plus-bulk",
        "shape spectrum"
      )
    ),
    variance = c(estimated = fit$variance[1L]),
    components = components,
    diagnostics = list(
      spatial.median = fit$median.diagnostics,
      companion.route = "economy SVD with stable A / (I - A) representation",
      ridge.quadratic.route = fit$native$quadratic_route,
      ridge.quadratic.spectral.rank =
        fit$native$quadratic_spectral_rank,
      ridge.quadratic.spectral.terms =
        fit$native$quadratic_spectral_terms,
      ridge.quadratic.singular.tolerance =
        fit$native$quadratic_singular_tolerance,
      aspect.ratio = fit$native$aspect_ratio,
      raw.minimum.gram.eigenvalue =
        fit$native$raw_minimum_gram_eigenvalue,
      roundoff.clipped.gram.eigenvalues =
        fit$native$roundoff_clipped_gram_eigenvalues,
      internal.residual.scale = fit$scale,
      internal.null.displacement.scale = fit$null.scale,
      subtraction.overflow.fallback = fit$overflow.fallback,
      quadratic.roundoff.clipped =
        fit$native$quadratic_roundoff_clipped,
      zero.fitted.residuals.allowed = FALSE,
      extra.regularization = FALSE,
      bartlett.center.correction = paste(
        "not applied: the public paper states a=8n/p for its numerical",
        "program but does not specify the correction formula"
      ),
      alpha = alpha,
      critical.Z = stats::qnorm(1 - alpha),
      rejected = z > stats::qnorm(1 - alpha),
      log.p.value = tails$log.p.value
    ),
    n = nrow(fit$x), p = ncol(fit$x), call = call
  )
}

#' Elliptical regularized Hotelling test with Cauchy aggregation
#'
#' Computes the feasible ERHT statistic over a fixed deterministic ridge grid
#' and combines its upper-tail normal p-values with the analytic Cauchy rule
#' of Feng, Zhou, and Wang (2026). Equal weights are used by default. The
#' default grid \eqn{\{0.1,0.2,\ldots,1\}} is the grid used in the paper's
#' numerical implementation; this function does not reproduce any simulation
#' design from that paper.
#'
#' For positive weights \eqn{\varpi_k} summing to one,
#' \deqn{T_{CC}=\sum_k\varpi_k\tan[\pi\{1/2-p_k\}],\qquad
#' p_{CC}=1/2-\pi^{-1}\arctan(T_{CC}).}
#' The computation uses signed-log cotangents and an angle representation to
#' avoid overflow and cancellation in extreme tails. The article proves
#' dependence-robust small-tail validity for this analytic p-value; except in
#' special dependence structures it does not claim exact fixed-level
#' calibration at an ordinary level such as 0.05.
#'
#' @inheritParams elliptical_regularized_hotelling_test
#' @param rho A strictly increasing vector of at least two positive ridge
#'   values.
#' @param weights Optional prespecified positive Cauchy weights, one per ridge.
#'   They must be deterministic and chosen independently of the observed test
#'   statistics, as required by the fixed-grid result. They are normalized to
#'   sum to one. `NULL` gives equal weights.
#'
#' @return An object of class `c("hd_location_test", "htest")`. Its p-value is
#'   the analytic Cauchy combination. The marginal table contains every raw
#'   ERHT statistic, feasible center and variance, standardized statistic,
#'   p-value, log-p-value, and fixed-ridge decision.
#'
#' @references
#' Feng, L., Zhou, L., and Wang, X. (2026). Elliptical regularized Hotelling
#' testing for high dimensional data. arXiv:2606.25942.
#' \doi{10.48550/arXiv.2606.25942}.
#'
#' @examples
#' set.seed(2)
#' x <- matrix(rnorm(160), 40, 4)
#' elliptical_regularized_hotelling_cauchy_test(
#'   x, rho = c(0.2, 0.5, 1)
#' )
#'
#' @export
elliptical_regularized_hotelling_cauchy_test <- function(
    x, mu = NULL, rho = seq(0.1, 1, by = 0.1), weights = NULL,
    alpha = 0.05, tol = 1e-8, max_iter = 1000L, strict = TRUE,
    keep_companion = FALSE) {
  call <- match.call()
  data.name <- deparse1(substitute(x))
  rho <- .erht_rho(rho, multiple = TRUE)
  alpha <- .erht_alpha(alpha)
  if (is.null(weights)) weights <- rep.int(1 / length(rho), length(rho))
  weights <- as.numeric(weights)
  if (length(weights) != length(rho) || anyNA(weights) ||
      any(!is.finite(weights)) || any(weights <= 0)) {
    stop("`weights` must contain one finite, strictly positive value per ridge.",
         call. = FALSE)
  }
  weight.scale <- max(weights)
  weights.scaled <- weights / weight.scale
  if (any(weights.scaled == 0)) {
    stop(
      paste0(
        "The positive `weights` span more than the representable ",
        "double-precision range after normalization."
      ),
      call. = FALSE
    )
  }
  weights <- weights.scaled / sum(weights.scaled)
  if (any(weights == 0)) {
    stop(
      paste0(
        "The positive `weights` span more than the representable ",
        "double-precision range after normalization."
      ),
      call. = FALSE
    )
  }
  fit <- .erht_fit_grid(
    x, mu, rho, tol, max_iter, strict, keep_companion
  )
  z <- as.numeric(fit$native$z)
  tails <- .erht_tail(z)
  combined <- .erht_cauchy_grid(tails, weights)
  components <- .erht_components(fit, tails, alpha)
  components$cauchy <- c(
    statistic = combined$statistic,
    sign = combined$statistic.sign,
    log.absolute.statistic = combined$log.absolute.statistic,
    angle = combined$angle,
    p.value = combined$p.value,
    log.p.value = combined$log.p.value,
    log.one.minus.p.value = combined$log.one.minus.p.value
  )
  components$cauchy.weights <- stats::setNames(
    weights, .erht_rho_labels(rho)
  )

  display.statistic <- if (is.finite(combined$statistic)) {
    combined$statistic
  } else {
    combined$statistic.sign * .Machine$double.xmax
  }
  .new_hd_location_test(
    statistic = c(T.CC = display.statistic),
    parameter = c(K = length(rho)),
    p.value = combined$p.value,
    method = paste(
      "Feng-Zhou-Wang elliptical regularized Hotelling test",
      "with fixed-grid Cauchy aggregation"
    ),
    data.name = data.name,
    alternative = "two.sided",
    raw.statistic = stats::setNames(
      fit$raw.statistic, .erht_rho_labels(rho)
    ),
    estimate = fit$location,
    null.value = fit$mu,
    null.distribution = list(
      family = "analytic standard-Cauchy combination",
      exact = FALSE,
      tail = "upper",
      fixed.level.exact = FALSE,
      validity = "dependence-robust asymptotic small-tail calibration",
      assumptions = paste(
        "the fixed-grid joint ERHT limit under elliptical symmetry,",
        "proportional growth, and pervasive-plus-bulk dependence"
      )
    ),
    variance = fit$variance,
    components = components,
    diagnostics = list(
      spatial.median = fit$median.diagnostics,
      ridge.grid = rho,
      weights = weights,
      grid.fixed = TRUE,
      companion.route = "one economy SVD with stable A / (I - A) representation",
      ridge.quadratic.route = fit$native$quadratic_route,
      ridge.quadratic.spectral.rank =
        fit$native$quadratic_spectral_rank,
      ridge.quadratic.spectral.terms =
        fit$native$quadratic_spectral_terms,
      ridge.quadratic.singular.tolerance =
        fit$native$quadratic_singular_tolerance,
      aspect.ratio = fit$native$aspect_ratio,
      internal.residual.scale = fit$scale,
      internal.null.displacement.scale = fit$null.scale,
      subtraction.overflow.fallback = fit$overflow.fallback,
      quadratic.roundoff.clipped =
        fit$native$quadratic_roundoff_clipped,
      cauchy.evaluation = "signed-log weighted cotangents and stable angle",
      analytic.cauchy.fixed.level.exact = FALSE,
      extra.regularization = FALSE,
      bartlett.center.correction = paste(
        "not applied: the public paper states a=8n/p for its numerical",
        "program but does not specify the correction formula"
      ),
      alpha = alpha,
      rejected = combined$p.value <= alpha
    ),
    n = nrow(fit$x), p = ncol(fit$x), call = call
  )
}
