.os_validate_logical <- function(value, name) {
  if (!is.logical(value) || length(value) != 1L || is.na(value)) {
    stop(sprintf("`%s` must be TRUE or FALSE.", name), call. = FALSE)
  }
  value
}

.os_validate_positive <- function(value, name, allow_zero = FALSE) {
  value <- as.numeric(value)
  valid <- length(value) == 1L && !is.na(value) && is.finite(value)
  valid <- valid && if (allow_zero) value >= 0 else value > 0
  if (!valid) {
    qualifier <- if (allow_zero) "non-negative" else "positive"
    stop(sprintf("`%s` must be one finite %s number.", name, qualifier),
         call. = FALSE)
  }
  value
}

.os_validate_integer <- function(value, name, minimum = 1L,
                                 maximum = .Machine$integer.max) {
  if (!is.numeric(value) || length(value) != 1L || is.na(value) ||
      !is.finite(value) || value != floor(value) || value < minimum ||
      value > maximum) {
    stop(sprintf("`%s` must be an integer from %d through %d.",
                 name, minimum, maximum), call. = FALSE)
  }
  as.integer(value)
}

.os_matrix_names <- function(matrix, variable.names) {
  dimnames(matrix) <- list(variable.names, variable.names)
  matrix
}

.os_covariance_result <- function(estimate, method, call,
                                  sample.covariance = NULL,
                                  components = list(), diagnostics = list()) {
  if (!is.matrix(estimate) || nrow(estimate) != ncol(estimate) ||
      anyNA(estimate) || any(!is.finite(estimate))) {
    stop("The covariance estimate is not a finite square matrix.",
         call. = FALSE)
  }
  result <- list(
    estimate = estimate,
    sample.covariance = sample.covariance,
    method = method,
    components = components,
    diagnostics = diagnostics,
    call = call
  )
  class(result) <- c("hd_covariance_estimator", "list")
  result
}

.os_shape_result <- function(estimate, method, call, center,
                             components = list(), diagnostics = list()) {
  if (!is.matrix(estimate) || nrow(estimate) != ncol(estimate) ||
      anyNA(estimate) || any(!is.finite(estimate))) {
    stop("The shape estimate is not a finite square matrix.", call. = FALSE)
  }
  result <- list(
    estimate = estimate,
    center = center,
    method = method,
    components = components,
    diagnostics = diagnostics,
    call = call
  )
  class(result) <- c("hd_shape_estimator", "list")
  result
}

.os_restore_power <- function(value, scale, power, what) {
  result <- value
  source.nonzero <- value != 0
  for (iteration in seq_len(power)) result <- result * scale
  if (anyNA(result) || any(!is.finite(result))) {
    stop(sprintf("The %s overflows in the original measurement units.", what),
         call. = FALSE)
  }
  if (any(source.nonzero & result == 0)) {
    stop(sprintf("The %s underflows in the original measurement units.", what),
         call. = FALSE)
  }
  result
}

.os_kurtosis_from_g2 <- function(g2, n, p) {
  if (length(g2) != p || anyNA(g2) || any(!is.finite(g2))) {
    stop(
      paste0(
        "Every marginal variable must have positive sample variation for ",
        "the corrected kurtosis; no variable is silently omitted."
      ),
      call. = FALSE
    )
  }
  corrected <- (n - 1) / ((n - 2) * (n - 3)) * ((n + 1) * g2 + 6)
  raw <- mean(corrected) / 3
  lower <- -2 / (p + 2)
  estimate <- max(lower, raw)
  list(
    estimate = estimate,
    raw = raw,
    lower.bound = lower,
    at.lower.bound = estimate == lower,
    below.lower.before.projection = raw < lower,
    marginal.g2 = g2,
    corrected.marginal.kurtosis = corrected,
    official.code.boundary.repair = "not used (no 0.99 multiplier)"
  )
}

.os_sample_components <- function(x, min_rows = 4L) {
  x <- .as_data_matrix(x, "x", min_rows = min_rows)
  n <- nrow(x)
  p <- ncol(x)
  moment <- cpp_ollila_sample_moments(x)
  if (any(as.logical(moment$zero_variance))) {
    stop(
      paste0(
        "Every marginal variable must have positive sample variation; ",
        "zero-variance variables are not dropped."
      ),
      call. = FALSE
    )
  }
  scale2 <- .os_restore_power(1, moment$data_scale, 2,
                              "covariance scale")
  variable.names <- if (is.null(colnames(x))) {
    paste0("variable", seq_len(p))
  } else {
    colnames(x)
  }
  covariance <- .os_matrix_names(.os_restore_power(
    as.matrix(moment$covariance_scaled), as.numeric(moment$data_scale), 2,
    "sample covariance"
  ), variable.names)
  mean <- as.numeric(moment$mean_scaled) * moment$data_scale
  names(mean) <- variable.names
  kurtosis <- .os_kurtosis_from_g2(as.numeric(moment$g2), n, p)
  list(
    x = x, n = n, p = p, variable.names = variable.names,
    data.scale = as.numeric(moment$data_scale),
    scale.squared = scale2,
    mean = mean,
    covariance = covariance,
    covariance.scaled = as.matrix(moment$covariance_scaled),
    trace.scaled = as.numeric(moment$trace_covariance_scaled),
    trace.square.scaled =
      as.numeric(moment$trace_covariance_squared_scaled),
    kurtosis = kurtosis
  )
}

.os_spatial_median_diagnostics <- function(location, tol, max_iter, strict) {
  list(
    converged = isTRUE(attr(location, "converged")),
    iterations = as.integer(attr(location, "iterations")),
    tolerance = tol,
    max.iterations.allowed = max_iter,
    objective = as.numeric(attr(location, "objective")),
    relative.change = as.numeric(attr(location, "relative_change")),
    equation.residual = as.numeric(attr(location, "equation_residual")),
    strict = strict,
    algorithm = "Chapter 1 modified Weiszfeld spatial_median()"
  )
}

.os_resolve_center <- function(x, center, tol, max_iter, strict) {
  controls <- .validate_iteration_controls(tol, max_iter, zero_tol = 0)
  strict <- .os_validate_logical(strict, "strict")
  p <- ncol(x)
  if (is.numeric(center)) {
    value <- .as_location(center, p, "center")
    return(list(
      value = value, label = "supplied numeric center",
      diagnostics = list(
        estimated = FALSE, converged = NA, strict = strict,
        tolerance = controls$tol,
        max.iterations.allowed = controls$max_iter
      )
    ))
  }
  center <- match.arg(center, c("spatial", "mean", "none"))
  if (center == "spatial") {
    value <- spatial_median(
      x, tol = controls$tol, max_iter = controls$max_iter,
      zero_tol = 0, warn = FALSE
    )
    diagnostics <- .os_spatial_median_diagnostics(
      value, controls$tol, controls$max_iter, strict
    )
    if (!diagnostics$converged) {
      message <- sprintf(
        paste0(
          "The spatial-median iteration did not converge within %d steps ",
          "(equation residual %.6g, relative change %.6g)."
        ),
        controls$max_iter, diagnostics$equation.residual,
        diagnostics$relative.change
      )
      if (strict) stop(message, call. = FALSE)
      warning(paste0(message, " Returning the last finite iterate."),
              call. = FALSE)
    }
    value <- as.numeric(value)
    label <- "sample spatial median"
  } else if (center == "mean") {
    moment <- cpp_ollila_sample_moments(x)
    value <- as.numeric(moment$mean_scaled) * as.numeric(moment$data_scale)
    diagnostics <- list(
      estimated = TRUE, converged = TRUE, iterations = 0L,
      strict = strict, algorithm = "overflow-safe arithmetic sample mean"
    )
    label <- "sample arithmetic mean"
  } else {
    value <- numeric(p)
    diagnostics <- list(
      estimated = FALSE, converged = NA, iterations = 0L,
      strict = strict, algorithm = "no centering"
    )
    label <- "no centering (data treated as already centered)"
  }
  if (anyNA(value) || any(!is.finite(value))) {
    stop("The selected center is non-finite; no repair is applied.",
         call. = FALSE)
  }
  names(value) <- if (is.null(colnames(x))) NULL else colnames(x)
  list(value = value, label = label, diagnostics = diagnostics)
}

.os_sign_components <- function(x, center, keep_signs = FALSE,
                                require_positive_radii = TRUE) {
  scaled <- .center_and_scale(x, center, zero_tol = 0)
  core <- cpp_ollila_sign_components(scaled$x, keep_signs)
  if (require_positive_radii && as.numeric(core$n_zero) > 0) {
    stop(
      paste0(
        "A fitted residual radius is exactly zero, so the spatial-sign ",
        "estimator or inverse radial moments are undefined; no observation ",
        "is deleted and no radius floor is applied."
      ),
      call. = FALSE
    )
  }
  variable.names <- if (is.null(colnames(x))) {
    paste0("variable", seq_len(ncol(x)))
  } else {
    colnames(x)
  }
  observation.names <- if (is.null(rownames(x))) {
    paste0("observation", seq_len(nrow(x)))
  } else {
    rownames(x)
  }
  core$sscm <- .os_matrix_names(as.matrix(core$sscm), variable.names)
  core$shape <- .os_matrix_names(as.matrix(core$shape), variable.names)
  dimnames(core$eigenvectors) <- list(variable.names, NULL)
  names(core$sscm_eigenvalues) <- paste0("component", seq_len(ncol(x)))
  names(core$radii_scaled) <- observation.names
  names(core$log_radii_scaled) <- observation.names
  if (!is.null(core$signs)) {
    dimnames(core$signs) <- list(observation.names, variable.names)
  }
  core$residual.scale <- scaled$scale
  core$centering.overflow.fallback <- scaled$overflow_fallback
  core
}

.os_log_mean_exp <- function(value) {
  maximum <- max(value)
  maximum + log(mean(exp(value - maximum)))
}

.os_radial_bias <- function(log_radii) {
  if (length(log_radii) < 1L || anyNA(log_radii) ||
      any(!is.finite(log_radii))) {
    stop(
      "All residual radii must be finite and positive; no floor is applied.",
      call. = FALSE
    )
  }
  n <- length(log_radii)
  log.m1 <- .os_log_mean_exp(-log_radii)
  log.r2 <- .os_log_mean_exp(-2 * log_radii) - 2 * log.m1
  log.r3 <- .os_log_mean_exp(-3 * log_radii) - 3 * log.m1
  r2 <- exp(log.r2)
  r3 <- exp(log.r3)
  delta <- (2 - 2 * r2 + r2^2) / n^2 +
    (8 * r2 - 6 * r2^2 + 2 * r2 * r3 - 2 * r3) / n^3
  if (any(!is.finite(c(r2, r3, delta)))) {
    stop("The radial bias correction is non-finite; no repair is applied.",
         call. = FALSE)
  }
  list(ratio2 = r2, ratio3 = r3, delta = delta,
       log.mean.inverse.radius = log.m1)
}

.os_rsscm_components <- function(core, n, p) {
  a <- as.numeric(core$a)
  if (!is.finite(a)) {
    stop("The RSSCM signal statistic `a` is non-finite.", call. = FALSE)
  }
  if (a < 1) {
    stop(
      paste0(
        "The computed RSSCM signal statistic satisfies a < 1, contrary to ",
        "the unit-trace SSCM identity; no numerical floor is applied."
      ),
      call. = FALSE
    )
  }
  if (a == 1) {
    alpha.raw <- NA_real_
    alpha.data <- 0
    signal.free <- TRUE
  } else {
    alpha.raw <- (n / (n - 1) * (a - p / n) - 1) / (a - 1)
    if (!is.finite(alpha.raw)) {
      stop("The RSSCM data weight is non-finite; no repair is applied.",
           call. = FALSE)
    }
    alpha.data <- min(1, max(0, alpha.raw))
    signal.free <- FALSE
  }
  estimate <- alpha.data * core$shape + (1 - alpha.data) * diag(p)
  list(
    estimate = estimate,
    alpha.data = alpha.data,
    alpha.raw = alpha.raw,
    shrinkage.intensity = 1 - alpha.data,
    a = a,
    signal.free = signal.free,
    projection = "published projection of alpha_data onto [0, 1]"
  )
}

.os_validate_basic_controls <- function(quadrature_order, inversion_tol,
                                        integration_tol, inversion_max_iter,
                                        max_bracket_iter) {
  list(
    quadrature.order = .os_validate_integer(
      quadrature_order, "quadrature_order", 8L, 2048L
    ),
    inversion.tolerance = .os_validate_positive(
      inversion_tol, "inversion_tol"
    ),
    integration.tolerance = .os_validate_positive(
      integration_tol, "integration_tol"
    ),
    max.iterations = .os_validate_integer(
      inversion_max_iter, "inversion_max_iter"
    ),
    max.bracket.iterations = .os_validate_integer(
      max_bracket_iter, "max_bracket_iter"
    )
  )
}

.os_basic_from_eigensystem <- function(eigenvectors, delta, p, controls,
                                       variable.names) {
  inversion <- cpp_ollila_basic_inverse(
    as.numeric(delta), p, controls$quadrature.order,
    controls$inversion.tolerance, controls$integration.tolerance,
    controls$max.iterations, controls$max.bracket.iterations
  )
  lambda <- as.numeric(inversion$lambda_normalized)
  estimate <- eigenvectors %*% (lambda * t(eigenvectors))
  estimate <- 0.5 * (estimate + t(estimate))
  estimate <- .os_matrix_names(estimate, variable.names)
  list(estimate = estimate, inversion = inversion)
}


#' Ollila--Raninen elliptical shrinkage covariance estimator
#'
#' Estimates a covariance matrix by
#' \deqn{\widehat\Sigma=\widehat\beta S+
#' (1-\widehat\beta)\widehat\eta I_p,}
#' where `S` is the unknown-mean unbiased sample covariance (divisor
#' \eqn{n-1}), \eqn{\widehat\eta=\operatorname{tr}(S)/p}, and
#' \deqn{\widehat\beta=\frac{\widehat\gamma-1}
#' {\widehat\gamma-1+\widehat\kappa(2\widehat\gamma+p)/n+
#' (\widehat\gamma+p)/(n-1)}}.
#'
#' The corrected marginal excess kurtoses are
#' \deqn{K_j=\frac{n-1}{(n-2)(n-3)}\{(n+1)g_{2j}+6\},}
#' and \eqn{\widehat\kappa=\max\{-2/(p+2),p^{-1}
#' \sum_jK_j/3\}}. Equality at the theoretical lower bound is retained and
#' diagnosed; the `0.99` boundary modification in the authors' MATLAB code is
#' deliberately not used.
#'
#' `sphericity = "ell1"` uses the spatial-median SSCM estimator
#' \deqn{\widehat\gamma_1^*=\frac{n}{n-1}
#' \{p\operatorname{tr}(S_{sgn}^2)-p/n\}.}
#' `"ell2"` uses
#' \deqn{\widehat\gamma_2^*=b_n\left\{
#' \frac{p\operatorname{tr}(S^2)}{\operatorname{tr}(S)^2}
#' -a_n\frac pn\right\},}
#' where \eqn{a_n=n(n/(n-1)+\widehat\kappa)/(n+\widehat\kappa)} and
#' \eqn{b_n=(n+\widehat\kappa)(n-1)^2/
#' [(n-2)\{3\widehat\kappa(n-1)+n(n+1)\}]}.
#' Both are projected to `[1, p]` as prescribed in the publication.
#' `"ell3"` selects the smaller projected sphericity, with ties assigned to
#' Ell2.
#'
#' @param x A finite real numeric matrix or data frame; observations are rows.
#'   At least four observations are required and every marginal variable must
#'   have positive sample variation.
#' @param sphericity One of `"ell1"`, `"ell2"`, or `"ell3"`.
#' @param tol,max_iter Convergence controls passed to [spatial_median()] when
#'   Ell1 is needed.
#' @param strict Whether spatial-median nonconvergence is an error. With
#'   `FALSE`, the last finite iterate is used with a warning.
#' @param keep_signs Whether to retain fitted spatial signs for Ell1/Ell3.
#'
#' @return An `hd_covariance_estimator` containing the estimate, unbiased SCM,
#'   scale, both raw/projected sphericities, corrected kurtosis, data weight,
#'   shrinkage intensity, and no-repair diagnostics.
#' @export
#'
#' @references
#' Ollila, E. and Raninen, E. (2019). Optimal shrinkage covariance matrix
#' estimation under random sampling from elliptical distributions.
#' *IEEE Transactions on Signal Processing*, 67, 2707--2719.
#' \doi{10.1109/TSP.2019.2908144}.
#'
#' @examples
#' x <- matrix(c(-2, 0, 1, -1, 2, 0, 0, -2, 1, 1, 1, -1,
#'               2, 1, 0, 3, -1, 2, -2, -1, -1, 1, 3, 1),
#'             ncol = 3, byrow = TRUE)
#' ollila_raninen_shrinkage_covariance(x, sphericity = "ell3")
ollila_raninen_shrinkage_covariance <- function(
    x, sphericity = c("ell1", "ell2", "ell3"),
    tol = 1e-8, max_iter = 1000L, strict = TRUE,
    keep_signs = FALSE) {
  call <- match.call()
  sphericity <- match.arg(sphericity)
  strict <- .os_validate_logical(strict, "strict")
  keep_signs <- .os_validate_logical(keep_signs, "keep_signs")
  sample <- .os_sample_components(x, min_rows = 4L)
  n <- sample$n
  p <- sample$p
  trace.s <- sample$trace.scaled
  if (!is.finite(trace.s) || trace.s <= 0) {
    stop("The sample covariance has non-positive trace.", call. = FALSE)
  }

  kappa <- sample$kurtosis$estimate
  a.n <- n / (n + kappa) * (n / (n - 1) + kappa)
  b.denominator <- (n - 2) * (3 * kappa * (n - 1) + n * (n + 1))
  if (!is.finite(b.denominator) || b.denominator == 0) {
    stop("The Ell2 coefficient denominator is zero or non-finite.",
         call. = FALSE)
  }
  b.n <- (kappa + n) * (n - 1)^2 / b.denominator
  gamma2.raw <- b.n * (
    p * sample$trace.square.scaled / trace.s^2 - a.n * p / n
  )
  if (!is.finite(gamma2.raw)) {
    stop("The Ell2 sphericity estimate is non-finite.", call. = FALSE)
  }
  gamma2 <- min(p, max(1, gamma2.raw))

  center.fit <- NULL
  sign <- NULL
  gamma1.raw <- gamma1 <- NA_real_
  if (sphericity != "ell2" && p == 1L) {
    gamma1.raw <- gamma1 <- 1
  } else if (sphericity != "ell2") {
    center.fit <- .os_resolve_center(
      sample$x, "spatial", tol, max_iter, strict
    )
    sign <- .os_sign_components(
      sample$x, center.fit$value, keep_signs = keep_signs,
      require_positive_radii = TRUE
    )
    gamma1.raw <- n / (n - 1) * (p * sum(sign$sscm^2) - p / n)
    if (!is.finite(gamma1.raw)) {
      stop("The Ell1 sphericity estimate is non-finite.", call. = FALSE)
    }
    gamma1 <- min(p, max(1, gamma1.raw))
  }

  if (sphericity == "ell1") {
    gamma <- gamma1
    selected <- "ell1"
  } else if (sphericity == "ell2") {
    gamma <- gamma2
    selected <- "ell2"
  } else if (gamma1 < gamma2) {
    gamma <- gamma1
    selected <- "ell1"
  } else {
    gamma <- gamma2
    selected <- "ell2"
  }
  beta.denominator <- (gamma - 1) + kappa * (2 * gamma + p) / n +
    (gamma + p) / (n - 1)
  if (!is.finite(beta.denominator) || beta.denominator <= 0) {
    stop(
      "The published SCM data-weight denominator is not finite and positive.",
      call. = FALSE
    )
  }
  beta <- (gamma - 1) / beta.denominator
  if (!is.finite(beta) || beta < 0 || beta > 1) {
    stop(
      "The published SCM data weight lies outside [0, 1]; no clamp is applied.",
      call. = FALSE
    )
  }
  eta <- .os_restore_power(
    trace.s / p, sample$data.scale, 2, "average sample variance"
  )
  estimate <- beta * sample$covariance + (1 - beta) * eta * diag(p)
  estimate <- .os_matrix_names(estimate, sample$variable.names)

  .os_covariance_result(
    estimate = estimate,
    method = sprintf(
      "Ollila-Raninen elliptical shrinkage covariance (%s)", sphericity
    ),
    call = call,
    sample.covariance = sample$covariance,
    components = list(
      eta = eta,
      kappa = kappa,
      corrected.marginal.kurtosis =
        sample$kurtosis$corrected.marginal.kurtosis,
      marginal.g2 = sample$kurtosis$marginal.g2,
      sphericity = gamma,
      sphericity.selected = selected,
      ell1 = list(raw = gamma1.raw, projected = gamma1),
      ell2 = list(raw = gamma2.raw, projected = gamma2,
                  a.n = a.n, b.n = b.n),
      scm.data.weight = beta,
      shrinkage.intensity = 1 - beta,
      fitted.location = if (is.null(center.fit)) NULL else center.fit$value,
      sscm = if (is.null(sign)) NULL else sign$sscm,
      fitted.signs = if (is.null(sign)) NULL else sign$signs
    ),
    diagnostics = list(
      n = n, p = p,
      p.one.sphericity.known = p == 1L,
      covariance.divisor = n - 1,
      covariance.mean = "unknown; arithmetic sample mean removed",
      internal.data.scale = sample$data.scale,
      kappa.raw = sample$kurtosis$raw,
      kappa.lower.bound = sample$kurtosis$lower.bound,
      kappa.at.lower.bound = sample$kurtosis$at.lower.bound,
      kappa.below.lower.before.projection =
        sample$kurtosis$below.lower.before.projection,
      official.kappa.boundary.repair =
        sample$kurtosis$official.code.boundary.repair,
      ell3.tie.rule = "ties select Ell2",
      spatial.median = if (is.null(center.fit)) NULL else
        center.fit$diagnostics,
      zero.residuals = if (is.null(sign)) NA_real_ else
        as.numeric(sign$n_zero),
      regularization = "published scalar identity shrinkage only",
      ridge = "none",
      pseudoinverse = "none",
      numerical.floor = "none",
      deletion = "none"
    )
  )
}


#' Regularized spatial-sign covariance (RSSCM)
#'
#' Forms \eqn{V=pS_{sgn}}, lets \eqn{a=\operatorname{tr}(V^2)/p}, and
#' estimates the data weight by
#' \deqn{\widehat\alpha_{raw}=
#' \frac{\{n/(n-1)\}(a-p/n)-1}{a-1}.}
#' The published estimator projects this value to `[0, 1]` and returns
#' \eqn{\widehat\alpha V+(1-\widehat\alpha)I_p}. When `a` is exactly one,
#' the formula is `0/0`; this implementation returns data weight zero and
#' explicitly labels the signal-free case instead of producing `NaN`.
#'
#' Spatial-median centering is the paper-facing default. `center = "mean"`,
#' `"none"`, or a numeric vector are practical alternatives and are labelled
#' as such in diagnostics. Exact zero residuals are errors: the official code's
#' observation deletion is not reproduced.
#'
#' @param x A finite numeric matrix or data frame, observations in rows.
#' @param center One of `"spatial"`, `"mean"`, `"none"`, or a supplied finite
#'   numeric center. `"none"` treats the data as already centered.
#' @param tol,max_iter,strict Spatial-median convergence controls.
#' @param keep_signs Whether to retain the fitted sign matrix.
#'
#' @return An `hd_covariance_estimator` whose estimate has trace `p`, with
#'   SSCM, raw/projected data weight, shrinkage intensity, centre, and
#'   diagnostics.
#' @export
#'
#' @references
#' Raninen, E. and Ollila, E. (2022). Bias adjusted sign covariance matrix.
#' *IEEE Signal Processing Letters*, 29, 339--343.
#' \doi{10.1109/LSP.2021.3134940}.
#'
#' @examples
#' x <- matrix(c(-2, 0, -1, 2, 0, -2, 1, 1, 2, -1, 3, 2),
#'             ncol = 2, byrow = TRUE)
#' regularized_spatial_sign_covariance(x)
regularized_spatial_sign_covariance <- function(
    x, center = c("spatial", "mean", "none"),
    tol = 1e-8, max_iter = 1000L, strict = TRUE,
    keep_signs = FALSE) {
  call <- match.call()
  x <- .as_data_matrix(x, "x", min_rows = 2L)
  keep_signs <- .os_validate_logical(keep_signs, "keep_signs")
  variable.names <- if (is.null(colnames(x))) {
    paste0("variable", seq_len(ncol(x)))
  } else {
    colnames(x)
  }
  if (ncol(x) == 1L) {
    estimate <- matrix(1, 1, 1, dimnames = list(variable.names, variable.names))
    return(.os_covariance_result(
      estimate, "Regularized spatial-sign covariance (p = 1)", call,
      sample.covariance = estimate,
      components = list(
        sscm = estimate, unregularized.shape = estimate,
        a = 1, alpha.raw = NA_real_, alpha.data = 0,
        shrinkage.intensity = 1, fitted.location = NA_real_,
        fitted.signs = NULL
      ),
      diagnostics = list(
        p.one.special.case = TRUE, signal.free = TRUE,
        center = "irrelevant for the one-dimensional shape",
        regularization = "identity", ridge = "none",
        numerical.floor = "none", deletion = "none"
      )
    ))
  }
  center.fit <- .os_resolve_center(x, center, tol, max_iter, strict)
  sign <- .os_sign_components(
    x, center.fit$value, keep_signs = keep_signs,
    require_positive_radii = TRUE
  )
  rsscm <- .os_rsscm_components(sign, nrow(x), ncol(x))
  rsscm$estimate <- .os_matrix_names(rsscm$estimate, variable.names)

  .os_covariance_result(
    rsscm$estimate,
    "Regularized spatial-sign covariance (RSSCM)", call,
    sample.covariance = sign$sscm,
    components = list(
      sscm = sign$sscm,
      unregularized.shape = sign$shape,
      a = rsscm$a,
      alpha.raw = rsscm$alpha.raw,
      alpha.data = rsscm$alpha.data,
      shrinkage.intensity = rsscm$shrinkage.intensity,
      projection = rsscm$projection,
      fitted.location = center.fit$value,
      fitted.radii.scaled = sign$radii_scaled,
      fitted.signs = sign$signs
    ),
    diagnostics = list(
      p.one.special.case = FALSE,
      signal.free = rsscm$signal.free,
      center = center.fit$label,
      center.details = center.fit$diagnostics,
      zero.residuals = as.numeric(sign$n_zero),
      sscm.trace = as.numeric(sign$trace_sscm),
      residual.scale = sign$residual.scale,
      centering.overflow.fallback = sign$centering.overflow.fallback,
      official.zero.deletion = "not used",
      regularization = "published scalar identity shrinkage only",
      ridge = "none", pseudoinverse = "none",
      numerical.floor = "none", deletion = "none"
    )
  )
}


#' BASIC bias-adjusted spatial-sign shape estimator
#'
#' Applies the real-valued BASIC approximate inverse map independently to the
#' eigenvalues of the fitted SSCM. For a shape eigenvalue \eqn{\lambda}, the
#' map used by Raninen and Ollila is
#' \deqn{\widetilde\delta(\lambda)=\frac12\int_0^1
#' \frac{\lambda}{1-t+t\lambda}(1-t)^{p/2-1}\,dt.}
#' The implementation evaluates the equivalent smooth integral after
#' \eqn{1-t=v^2} with paired Gauss--Legendre rules, brackets each inverse root,
#' and uses bisection. It is independent of the official lookup table:
#' no interpolation, spline extrapolation, eigenvalue floor, or out-of-range
#' repair is used. Raw inverted eigenvalues are finally normalized to sum `p`.
#'
#' Boundary SSCM eigenvalues equal to zero map exactly to zero and are reported.
#' For `p = 1`, shape is identically one and no numerical inversion is needed.
#'
#' @inheritParams regularized_spatial_sign_covariance
#' @param quadrature_order Order of the coarse Gauss--Legendre rule; a rule of
#'   twice this order supplies the returned value and error estimate.
#' @param inversion_tol Positive absolute/relative bisection tolerance.
#' @param integration_tol Positive maximum coarse/fine quadrature discrepancy.
#' @param inversion_max_iter Maximum bisection iterations per eigenvalue.
#' @param max_bracket_iter Maximum upper-bracket doublings per eigenvalue.
#'
#' @return An `hd_shape_estimator` with the BASIC shape, SSCM eigensystem,
#'   raw and normalized inverse eigenvalues, mapped values, integration and
#'   inversion errors, brackets, iterations, and boundary diagnostics.
#' @export
#'
#' @references
#' Raninen, E. and Ollila, E. (2022). Bias adjusted sign covariance matrix.
#' *IEEE Signal Processing Letters*, 29, 339--343.
#' \doi{10.1109/LSP.2021.3134940}.
#'
#' @examples
#' x <- matrix(c(-2, 0, -1, 2, 0, -2, 1, 1, 2, -1, 3, 2),
#'             ncol = 2, byrow = TRUE)
#' basic_shape(x, quadrature_order = 32)
basic_shape <- function(
    x, center = c("spatial", "mean", "none"),
    tol = 1e-8, max_iter = 1000L, strict = TRUE,
    keep_signs = FALSE, quadrature_order = 64L,
    inversion_tol = 1e-10, integration_tol = 1e-10,
    inversion_max_iter = 100L,
    max_bracket_iter = 100L) {
  call <- match.call()
  x <- .as_data_matrix(x, "x", min_rows = 2L)
  keep_signs <- .os_validate_logical(keep_signs, "keep_signs")
  variable.names <- if (is.null(colnames(x))) {
    paste0("variable", seq_len(ncol(x)))
  } else {
    colnames(x)
  }
  if (ncol(x) == 1L) {
    estimate <- matrix(1, 1, 1, dimnames = list(variable.names, variable.names))
    return(.os_shape_result(
      estimate, "BASIC shape estimator (p = 1)", call,
      center = NA_real_,
      components = list(
        sscm = estimate, delta = 1, lambda.raw = 1,
        lambda.normalized = 1
      ),
      diagnostics = list(
        p.one.special.case = TRUE, inversion = "not required",
        boundary.zero = FALSE, extrapolation = "none",
        ridge = "none", numerical.floor = "none", deletion = "none"
      )
    ))
  }
  controls <- .os_validate_basic_controls(
    quadrature_order, inversion_tol, integration_tol,
    inversion_max_iter, max_bracket_iter
  )
  center.fit <- .os_resolve_center(x, center, tol, max_iter, strict)
  sign <- .os_sign_components(
    x, center.fit$value, keep_signs = keep_signs,
    require_positive_radii = TRUE
  )
  corrected <- .os_basic_from_eigensystem(
    sign$eigenvectors, sign$sscm_eigenvalues, ncol(x), controls,
    variable.names
  )
  inversion <- corrected$inversion

  .os_shape_result(
    corrected$estimate,
    "BASIC bias-adjusted spatial-sign shape estimator", call,
    center = center.fit$value,
    components = list(
      sscm = sign$sscm,
      sscm.eigenvalues = sign$sscm_eigenvalues,
      eigenvectors = sign$eigenvectors,
      lambda.raw = inversion$lambda_raw,
      lambda.normalized = inversion$lambda_normalized,
      delta.mapped = inversion$delta_mapped,
      fitted.radii.scaled = sign$radii_scaled,
      fitted.signs = sign$signs
    ),
    diagnostics = list(
      p.one.special.case = FALSE,
      center = center.fit$label,
      center.details = center.fit$diagnostics,
      zero.residuals = as.numeric(sign$n_zero),
      boundary.zero = as.logical(inversion$boundary_zero),
      integration.error = inversion$integration_error,
      inversion.error = inversion$inversion_error,
      maximum.integration.error = inversion$max_integration_error,
      maximum.inversion.error = inversion$max_inversion_error,
      brackets = inversion$brackets,
      inversion.iterations = inversion$iterations,
      quadrature.order.coarse = inversion$quadrature_order_coarse,
      quadrature.order.fine = inversion$quadrature_order_fine,
      inversion.converged = isTRUE(inversion$converged),
      extrapolation = inversion$extrapolation,
      map = paste(
        "real BASIC approximate map with one-half integral factor;",
        "independent Gauss-Legendre inversion"
      ),
      ridge = "none", pseudoinverse = "none",
      numerical.floor = "none", deletion = "none"
    )
  )
}


#' BASICS shrinkage plus bias-adjusted shape estimator
#'
#' Computes the RSSCM data weight, forms
#' \eqn{\widehat\Lambda_{RSSCM}=\widehat\alpha pS_{sgn}+
#' (1-\widehat\alpha)I_p}, and applies the same self-contained BASIC inverse
#' map as [basic_shape()] to the eigenvalues of
#' \eqn{\widehat\Lambda_{RSSCM}/p}. This follows the BASICS construction while
#' avoiding the official code's spline extrapolation. The RSSCM and BASIC
#' numerical diagnostics are both retained.
#'
#' @inheritParams basic_shape
#'
#' @return An `hd_shape_estimator` containing the BASICS shape, RSSCM, its
#'   raw/projected shrinkage weight, inverse-map details, and full diagnostics.
#' @export
#'
#' @references
#' Raninen, E. and Ollila, E. (2022). Bias adjusted sign covariance matrix.
#' *IEEE Signal Processing Letters*, 29, 339--343.
#' \doi{10.1109/LSP.2021.3134940}.
#'
#' @examples
#' x <- matrix(c(-2, 0, -1, 2, 0, -2, 1, 1, 2, -1, 3, 2),
#'             ncol = 2, byrow = TRUE)
#' basics_shape(x, quadrature_order = 32)
basics_shape <- function(
    x, center = c("spatial", "mean", "none"),
    tol = 1e-8, max_iter = 1000L, strict = TRUE,
    keep_signs = FALSE, quadrature_order = 64L,
    inversion_tol = 1e-10, integration_tol = 1e-10,
    inversion_max_iter = 100L,
    max_bracket_iter = 100L) {
  call <- match.call()
  x <- .as_data_matrix(x, "x", min_rows = 2L)
  keep_signs <- .os_validate_logical(keep_signs, "keep_signs")
  variable.names <- if (is.null(colnames(x))) {
    paste0("variable", seq_len(ncol(x)))
  } else {
    colnames(x)
  }
  if (ncol(x) == 1L) {
    estimate <- matrix(1, 1, 1, dimnames = list(variable.names, variable.names))
    return(.os_shape_result(
      estimate, "BASICS shape estimator (p = 1)", call,
      center = NA_real_,
      components = list(
        sscm = estimate, rsscm = estimate, a = 1,
        alpha.raw = NA_real_, alpha.data = 0,
        shrinkage.intensity = 1, delta = 1,
        lambda.raw = 1, lambda.normalized = 1
      ),
      diagnostics = list(
        p.one.special.case = TRUE, signal.free = TRUE,
        inversion = "not required", extrapolation = "none",
        ridge = "none", numerical.floor = "none", deletion = "none"
      )
    ))
  }
  controls <- .os_validate_basic_controls(
    quadrature_order, inversion_tol, integration_tol,
    inversion_max_iter, max_bracket_iter
  )
  center.fit <- .os_resolve_center(x, center, tol, max_iter, strict)
  sign <- .os_sign_components(
    x, center.fit$value, keep_signs = keep_signs,
    require_positive_radii = TRUE
  )
  rsscm <- .os_rsscm_components(sign, nrow(x), ncol(x))
  delta <- rsscm$alpha.data * as.numeric(sign$sscm_eigenvalues) +
    (1 - rsscm$alpha.data) / ncol(x)
  corrected <- .os_basic_from_eigensystem(
    sign$eigenvectors, delta, ncol(x), controls, variable.names
  )
  rsscm$estimate <- .os_matrix_names(rsscm$estimate, variable.names)
  inversion <- corrected$inversion

  .os_shape_result(
    corrected$estimate,
    "BASICS regularized bias-adjusted spatial-sign shape estimator", call,
    center = center.fit$value,
    components = list(
      sscm = sign$sscm,
      unregularized.shape = sign$shape,
      rsscm = rsscm$estimate,
      a = rsscm$a,
      alpha.raw = rsscm$alpha.raw,
      alpha.data = rsscm$alpha.data,
      shrinkage.intensity = rsscm$shrinkage.intensity,
      rsscm.eigenvalues.divided.by.p = delta,
      eigenvectors = sign$eigenvectors,
      lambda.raw = inversion$lambda_raw,
      lambda.normalized = inversion$lambda_normalized,
      delta.mapped = inversion$delta_mapped,
      fitted.radii.scaled = sign$radii_scaled,
      fitted.signs = sign$signs
    ),
    diagnostics = list(
      p.one.special.case = FALSE,
      signal.free = rsscm$signal.free,
      alpha.projection = rsscm$projection,
      center = center.fit$label,
      center.details = center.fit$diagnostics,
      zero.residuals = as.numeric(sign$n_zero),
      boundary.zero = as.logical(inversion$boundary_zero),
      integration.error = inversion$integration_error,
      inversion.error = inversion$inversion_error,
      maximum.integration.error = inversion$max_integration_error,
      maximum.inversion.error = inversion$max_inversion_error,
      brackets = inversion$brackets,
      inversion.iterations = inversion$iterations,
      quadrature.order.coarse = inversion$quadrature_order_coarse,
      quadrature.order.fine = inversion$quadrature_order_fine,
      inversion.converged = isTRUE(inversion$converged),
      extrapolation = inversion$extrapolation,
      map = paste(
        "real BASIC approximate map with one-half integral factor;",
        "independent Gauss-Legendre inversion"
      ),
      ridge = "none", pseudoinverse = "none",
      numerical.floor = "none", deletion = "none"
    )
  )
}


.os_validate_pool_data <- function(data) {
  if (is.matrix(data) || is.data.frame(data)) data <- list(data)
  if (!is.list(data) || length(data) < 1L) {
    stop("`data` must be a non-empty list of sample matrices.",
         call. = FALSE)
  }
  if (is.null(names(data))) names(data) <- paste0("group", seq_along(data))
  empty.names <- names(data) == "" | is.na(names(data))
  names(data)[empty.names] <- paste0("group", which(empty.names))
  if (anyDuplicated(names(data))) {
    stop("`data` group names must be unique.", call. = FALSE)
  }
  group.names <- names(data)
  data <- lapply(seq_along(data), function(index) {
    .as_data_matrix(data[[index]], sprintf("data[[%d]]", index), min_rows = 4L)
  })
  names(data) <- group.names
  dimensions <- vapply(data, ncol, integer(1))
  if (length(unique(dimensions)) != 1L) {
    stop("All groups in `data` must have the same number of columns.",
         call. = FALSE)
  }
  column.names <- lapply(data, colnames)
  has.names <- !vapply(column.names, is.null, logical(1))
  if (any(has.names) && !all(has.names)) {
    stop("Either every group or no group must supply column names.",
         call. = FALSE)
  }
  if (all(has.names)) {
    reference <- column.names[[1L]]
    if (any(!vapply(column.names, identical, logical(1), reference))) {
      stop("Column names and their order must agree across all groups.",
           call. = FALSE)
    }
    variable.names <- reference
  } else {
    variable.names <- paste0("variable", seq_len(dimensions[[1L]]))
  }
  names(data) <- group.names
  list(data = data, group.names = group.names,
       variable.names = variable.names, p = dimensions[[1L]])
}

.os_pool_group_components <- function(x.scaled, group.name,
                                      variable.names, tol, max_iter,
                                      strict, keep_signs) {
  n <- nrow(x.scaled)
  p <- ncol(x.scaled)
  moment <- cpp_ollila_sample_moments(x.scaled)
  if (any(as.logical(moment$zero_variance))) {
    stop(sprintf(
      paste0(
        "Group `%s` has a zero-variance marginal variable; no variable is ",
        "silently omitted."
      ), group.name
    ), call. = FALSE)
  }
  covariance.scaled <- .os_restore_power(
    as.matrix(moment$covariance_scaled), as.numeric(moment$data_scale), 2,
    sprintf("internally scaled covariance for group `%s`", group.name)
  )
  covariance.scaled <- .os_matrix_names(covariance.scaled, variable.names)
  trace.s <- sum(diag(covariance.scaled))
  eta <- trace.s / p
  kurtosis <- .os_kurtosis_from_g2(as.numeric(moment$g2), n, p)

  center.fit <- .os_resolve_center(
    x.scaled, "spatial", tol, max_iter, strict
  )
  sign <- .os_sign_components(
    x.scaled, center.fit$value, keep_signs = keep_signs,
    require_positive_radii = TRUE
  )
  gamma0.raw <- p * n / (n - 1) * (sum(sign$sscm^2) - 1 / n)
  radial <- .os_radial_bias(as.numeric(sign$log_radii_scaled))
  gamma.raw <- gamma0.raw - p * radial$delta
  if (any(!is.finite(c(gamma0.raw, gamma.raw, eta)))) {
    stop(sprintf("Group `%s` produced a non-finite pooling parameter.",
                 group.name), call. = FALSE)
  }
  gamma <- min(p, max(1, gamma.raw))
  tr.sigma.squared <- p * eta^2 * gamma
  tau1 <- 1 / (n - 1) + kurtosis$estimate / n
  tau2 <- kurtosis$estimate / n
  delta.diagonal <- (
    tau1 * trace.s^2 + (tau1 + tau2) * tr.sigma.squared
  ) / p
  if (!is.finite(delta.diagonal) || delta.diagonal < 0) {
    stop(sprintf(
      "Group `%s` has a negative or non-finite estimated SCM MSE.",
      group.name
    ), call. = FALSE)
  }
  list(
    n = n,
    covariance.scaled = covariance.scaled,
    sample.mean.scaled = as.numeric(moment$mean_scaled) *
      as.numeric(moment$data_scale),
    spatial.center.scaled = center.fit$value,
    spatial.center.diagnostics = center.fit$diagnostics,
    sscm = sign$sscm,
    fitted.signs = sign$signs,
    radii.scaled = sign$radii_scaled,
    eta.scaled = eta,
    trace.scaled = trace.s,
    kappa = kurtosis$estimate,
    kurtosis = kurtosis,
    gamma0.raw = gamma0.raw,
    gamma.raw = gamma.raw,
    gamma = gamma,
    radial = radial,
    tr.sigma.squared.scaled = tr.sigma.squared,
    tau1 = tau1,
    tau2 = tau2,
    delta.diagonal.scaled = delta.diagonal,
    zero.residuals = as.numeric(sign$n_zero),
    residual.scale = sign$residual.scale
  )
}


#' Linear or convex pooling of covariance matrices
#'
#' Estimates every group covariance as a non-negative linear combination of
#' all unbiased group SCMs and, optionally, the identity. For target group
#' `k`, coefficients solve
#' \deqn{\min_a\;\frac12a^\mathsf{T}(C+\Delta)a-C_{\cdot k}^\mathsf{T}a,}
#' subject to the requested lower bounds. `method = "convex"` additionally
#' imposes \eqn{1^\mathsf{T}a=1}.
#'
#' Each SCM uses divisor \eqn{n_k-1}. For group `k`, the spatial-median SSCM
#' and positive residual radii give
#' \deqn{\widehat\gamma_{0k}=\frac{pn_k}{n_k-1}
#' \{\operatorname{tr}(S_{sgn,k}^2)-1/n_k\},\qquad
#' \widehat\gamma_k=\Pi_{[1,p]}
#' (\widehat\gamma_{0k}-p\widehat\delta_k),}
#' with the Zou et al. inverse-radial bias correction. Corrected marginal
#' kurtosis follows [ollila_raninen_shrinkage_covariance()]. The real-valued
#' diagonal MSE term is
#' \deqn{\Delta_{kk}=p^{-1}\{\tau_{1k}\operatorname{tr}(S_k)^2+
#' (\tau_{1k}+\tau_{2k})p\eta_k^2\gamma_k\},}
#' where \eqn{\tau_{1k}=1/(n_k-1)+\kappa_k/n_k} and
#' \eqn{\tau_{2k}=\kappa_k/n_k}. The diagonal of `C` is
#' \eqn{\eta_k^2\gamma_k}; off-diagonal entries use
#' \eqn{\operatorname{tr}(S_{sgn,k}S_{sgn,l})
#' \operatorname{tr}(S_k)\operatorname{tr}(S_l)/p}.
#'
#' The independent projected-FISTA solver adds no ridge or pseudoinverse. It
#' returns feasibility, projected-gradient, KKT, complementarity-gap,
#' eigenvalue and restart diagnostics. Singular convex objectives are allowed.
#' Nonconvergence is an error by default; `strict = FALSE` returns the last
#' feasible iterate with a warning and labels it uncalibrated.
#'
#' @param data A non-empty list of finite numeric sample matrices. A single
#'   matrix is accepted and wrapped as one group. Rows are observations; every
#'   group needs at least four rows, a common column dimension, matching column
#'   names, positive marginal variation, and positive spatial residual radii.
#' @param method Either non-negative `"linear"` pooling or `"convex"` pooling
#'   with coefficients summing to one.
#' @param identity Whether to add the identity as an additional target.
#' @param identity_lower A non-negative scalar or one value per target group,
#'   giving the identity coefficient lower bound. The default zero does not
#'   reproduce the official code's numerical `1e-8` choice.
#' @param tol,max_iter Spatial-median convergence controls.
#' @param solver_tol Positive projected-gradient/KKT solver tolerance.
#' @param solver_max_iter Positive QP iteration limit.
#' @param strict Whether spatial-median or QP nonconvergence is an error.
#' @param keep_signs Whether to retain each group's fitted sign matrix.
#'
#' @return A `linear_pool_covariance` list with coefficient matrix, pooled
#'   estimates, SCMs, SSCMs, centers, radii, eta/gamma/kappa, `C`, `Delta`,
#'   scaled objective matrices, and per-target solver diagnostics.
#' @export
#'
#' @references
#' Raninen, E., Tyler, D. E., and Ollila, E. (2022). Linear pooling of sample
#' covariance matrices. *IEEE Transactions on Signal Processing*, 70,
#' 659--672. \doi{10.1109/TSP.2021.3139207}.
#'
#' @examples
#' x1 <- matrix(c(-2, 0, -1, 2, 0, -2, 1, 1, 2, -1, 3, 2,
#'                -1, -2, 2, 0), ncol = 2, byrow = TRUE)
#' x2 <- matrix(c(-1, 1, 0, 3, 1, -2, 2, 2, 3, 0, -2, -1,
#'                1, 2, 0, -3), ncol = 2, byrow = TRUE)
#' linear_pool_covariance(list(first = x1, second = x2))
linear_pool_covariance <- function(
    data, method = c("linear", "convex"), identity = TRUE,
    identity_lower = 0, tol = 1e-8, max_iter = 1000L,
    solver_tol = 1e-10, solver_max_iter = 10000L,
    strict = TRUE, keep_signs = FALSE) {
  call <- match.call()
  method <- match.arg(method)
  identity <- .os_validate_logical(identity, "identity")
  strict <- .os_validate_logical(strict, "strict")
  keep_signs <- .os_validate_logical(keep_signs, "keep_signs")
  solver_tol <- .os_validate_positive(solver_tol, "solver_tol")
  solver_max_iter <- .os_validate_integer(
    solver_max_iter, "solver_max_iter"
  )
  groups <- .os_validate_pool_data(data)
  K <- length(groups$data)
  p <- groups$p
  if (!is.numeric(identity_lower) || anyNA(identity_lower) ||
      any(!is.finite(identity_lower)) || any(identity_lower < 0) ||
      !(length(identity_lower) %in% c(1L, K))) {
    stop(
      "`identity_lower` must be a finite non-negative scalar or one value per group.",
      call. = FALSE
    )
  }
  identity.lower <- rep(as.numeric(identity_lower), length.out = K)
  names(identity.lower) <- groups$group.names
  if (!identity && any(identity.lower != 0)) {
    stop("`identity_lower` must be zero when `identity = FALSE`.",
         call. = FALSE)
  }
  if (method == "convex" && any(identity.lower > 1)) {
    stop("Convex pooling requires every identity lower bound to be at most one.",
         call. = FALSE)
  }

  global.scale <- max(vapply(
    groups$data, function(value) max(abs(value)), numeric(1)
  ))
  if (!is.finite(global.scale) || global.scale <= 0) {
    stop("The pooled samples must contain nonzero finite data.", call. = FALSE)
  }
  scale2 <- .os_restore_power(1, global.scale, 2,
                              "pooled covariance scale")
  scale4 <- .os_restore_power(1, global.scale, 4,
                              "pooled fourth-order scale")
  scaled.data <- lapply(groups$data, function(value) value / global.scale)
  group.components <- lapply(seq_len(K), function(index) {
    .os_pool_group_components(
      scaled.data[[index]], groups$group.names[[index]],
      groups$variable.names, tol, max_iter, strict, keep_signs
    )
  })
  names(group.components) <- groups$group.names

  eta.scaled <- vapply(group.components, `[[`, numeric(1), "eta.scaled")
  gamma <- vapply(group.components, `[[`, numeric(1), "gamma")
  gamma.raw <- vapply(group.components, `[[`, numeric(1), "gamma.raw")
  gamma0.raw <- vapply(group.components, `[[`, numeric(1), "gamma0.raw")
  kappa <- vapply(group.components, `[[`, numeric(1), "kappa")
  sample.sizes <- vapply(group.components, `[[`, integer(1), "n")
  names(eta.scaled) <- names(gamma) <- names(gamma.raw) <-
    names(gamma0.raw) <- names(kappa) <- names(sample.sizes) <-
    groups$group.names

  C.scaled <- matrix(0, K, K,
                     dimnames = list(groups$group.names, groups$group.names))
  diag(C.scaled) <- eta.scaled^2 * gamma
  for (left in seq_len(K - 1L)) {
    for (right in seq.int(left + 1L, K)) {
      value <- sum(
        group.components[[left]]$sscm * group.components[[right]]$sscm
      ) * group.components[[left]]$trace.scaled *
        group.components[[right]]$trace.scaled / p
      C.scaled[left, right] <- C.scaled[right, left] <- value
    }
  }
  Delta.scaled <- diag(
    vapply(group.components, `[[`, numeric(1), "delta.diagonal.scaled"),
    K, K
  )
  dimnames(Delta.scaled) <- list(groups$group.names, groups$group.names)
  if (any(!is.finite(C.scaled)) || any(!is.finite(Delta.scaled))) {
    stop("The estimated pooling objective is non-finite.", call. = FALSE)
  }

  if (identity) {
    identity.scaled <- 1 / scale2
    objective.C <- rbind(
      cbind(C.scaled, eta.scaled * identity.scaled),
      c(eta.scaled * identity.scaled, identity.scaled^2)
    )
    objective.Delta <- rbind(
      cbind(Delta.scaled, numeric(K)), numeric(K + 1L)
    )
    coefficient.names <- c(groups$group.names, "identity")
  } else {
    identity.scaled <- NULL
    objective.C <- C.scaled
    objective.Delta <- Delta.scaled
    coefficient.names <- groups$group.names
  }
  dimnames(objective.C) <- dimnames(objective.Delta) <-
    list(coefficient.names, coefficient.names)
  coefficients <- matrix(
    NA_real_, length(coefficient.names), K,
    dimnames = list(coefficient.names, groups$group.names)
  )
  solver <- vector("list", K)
  names(solver) <- groups$group.names
  for (target in seq_len(K)) {
    lower <- numeric(length(coefficient.names))
    if (identity) lower[[length(lower)]] <- identity.lower[[target]]
    fit <- cpp_ollila_pool_qp(
      objective.C + objective.Delta, objective.C[, target], lower,
      identical(method, "convex"), solver_tol, solver_max_iter
    )
    coefficients[, target] <- as.numeric(fit$solution)
    solver[[target]] <- fit$diagnostics
    if (!isTRUE(fit$diagnostics$converged)) {
      message <- sprintf(
        paste0(
          "Pooling QP for group `%s` did not converge within %d iterations ",
          "(scaled projected-gradient residual %.6g)."
        ),
        groups$group.names[[target]], solver_max_iter,
        fit$diagnostics$projected_gradient_residual_scaled
      )
      if (strict) stop(message, call. = FALSE)
      warning(paste0(message, " Returning the last feasible, uncalibrated iterate."),
              call. = FALSE)
    }
  }

  sample.covariances <- lapply(group.components, function(component) {
    .os_matrix_names(.os_restore_power(
      component$covariance.scaled, global.scale, 2,
      "pooled sample covariance"
    ), groups$variable.names)
  })
  spatial.sign.covariances <- lapply(group.components, `[[`, "sscm")
  pooled <- lapply(seq_len(K), function(target) {
    estimate <- matrix(0, p, p)
    for (source in seq_len(K)) {
      estimate <- estimate + coefficients[source, target] *
        sample.covariances[[source]]
    }
    if (identity) {
      estimate <- estimate + coefficients[K + 1L, target] * diag(p)
    }
    if (anyNA(estimate) || any(!is.finite(estimate))) {
      stop(sprintf("The pooled estimate for group `%s` is non-finite.",
                   groups$group.names[[target]]), call. = FALSE)
    }
    .os_matrix_names(0.5 * (estimate + t(estimate)), groups$variable.names)
  })
  names(sample.covariances) <- names(spatial.sign.covariances) <-
    names(pooled) <- groups$group.names

  C <- .os_restore_power(C.scaled, global.scale, 4, "pooling C matrix")
  Delta <- .os_restore_power(
    Delta.scaled, global.scale, 4, "pooling Delta matrix"
  )
  C.extended <- if (identity) {
    rbind(cbind(C, eta.scaled * scale2), c(eta.scaled * scale2, 1))
  } else C
  Delta.extended <- if (identity) {
    rbind(cbind(Delta, numeric(K)), numeric(K + 1L))
  } else Delta
  dimnames(C.extended) <- dimnames(Delta.extended) <-
    list(coefficient.names, coefficient.names)
  centers <- lapply(group.components, function(component) {
    value <- component$spatial.center.scaled * global.scale
    names(value) <- groups$variable.names
    value
  })
  radii <- lapply(group.components, `[[`, "radii.scaled")
  names(centers) <- names(radii) <- groups$group.names
  eta <- .os_restore_power(
    eta.scaled, global.scale, 2, "pooled average variances"
  )
  names(eta) <- groups$group.names
  tr.sigma.squared <- .os_restore_power(
    vapply(group.components, `[[`, numeric(1), "tr.sigma.squared.scaled"),
    global.scale, 4, "pooled squared covariance traces"
  )
  names(tr.sigma.squared) <- groups$group.names

  result <- list(
    estimates = pooled,
    coefficients = coefficients,
    sample.covariances = sample.covariances,
    spatial.sign.covariances = spatial.sign.covariances,
    spatial.centers = centers,
    residual.radii.scaled = radii,
    eta = eta,
    gamma = gamma,
    gamma.raw = gamma.raw,
    gamma0.raw = gamma0.raw,
    kappa = kappa,
    C = C,
    Delta = Delta,
    C.extended = C.extended,
    Delta.extended = Delta.extended,
    solver = solver,
    method = sprintf("%s covariance pooling", method),
    components = list(
      radial.bias = lapply(group.components, `[[`, "radial"),
      corrected.marginal.kurtosis = lapply(
        group.components,
        function(value) value$kurtosis$corrected.marginal.kurtosis
      ),
      kappa.raw = vapply(
        group.components, function(value) value$kurtosis$raw, numeric(1)
      ),
      kappa.lower.bound = vapply(
        group.components,
        function(value) value$kurtosis$lower.bound, numeric(1)
      ),
      kappa.at.lower.bound = vapply(
        group.components,
        function(value) value$kurtosis$at.lower.bound, logical(1)
      ),
      tau1 = vapply(group.components, `[[`, numeric(1), "tau1"),
      tau2 = vapply(group.components, `[[`, numeric(1), "tau2"),
      tr.Sigma.squared = tr.sigma.squared,
      fitted.signs = if (keep_signs) {
        lapply(group.components, `[[`, "fitted.signs")
      } else NULL,
      objective.C.scaled = objective.C,
      objective.Delta.scaled = objective.Delta
    ),
    diagnostics = list(
      groups = groups$group.names,
      sample.sizes = sample.sizes,
      p = p,
      covariance.divisor = "n_k - 1 for every group",
      covariance.mean = "unknown; arithmetic sample mean removed",
      spatial.center = "ordinary sample spatial median",
      spatial.median = lapply(
        group.components, `[[`, "spatial.center.diagnostics"
      ),
      zero.residuals = vapply(
        group.components, `[[`, numeric(1), "zero.residuals"
      ),
      residual.radius.scale.relative.to.global = vapply(
        group.components, `[[`, numeric(1), "residual.scale"
      ),
      global.data.scale = global.scale,
      fourth.order.scale = scale4,
      identity = identity,
      identity.lower = identity.lower,
      identity.default.difference = paste(
        "default is zero; official MATLAB convenience default 1e-8",
        "is not imposed"
      ),
      convex = identical(method, "convex"),
      coefficient.constraint = if (method == "convex") {
        "a >= lower and sum(a) = 1"
      } else {
        "a >= lower"
      },
      strict = strict,
      solver.tolerance = solver_tol,
      solver.max.iterations = solver_max_iter,
      calibrated = all(vapply(solver, function(value) {
        isTRUE(value$converged)
      }, logical(1))),
      official.zero.deletion = "not used",
      official.radius.floor = "not used",
      official.kappa.boundary.repair = "not used",
      ridge = "none", pseudoinverse = "none",
      numerical.floor = "none", perturbation = "none"
    ),
    call = call
  )
  class(result) <- c("linear_pool_covariance", "list")
  result
}
