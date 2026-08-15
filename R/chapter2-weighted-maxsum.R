.yzf_validate_m <- function(m) {
  m <- as.numeric(m)
  if (length(m) != 1L || is.na(m) || !is.finite(m) || m > 1) {
    stop("`m` must be one finite number no greater than one.",
         call. = FALSE)
  }
  m
}


.yzf_validate_strict <- function(strict) {
  if (!is.logical(strict) || length(strict) != 1L || is.na(strict)) {
    stop("`strict` must be `TRUE` or `FALSE`.", call. = FALSE)
  }
  strict
}


.yzf_validate_alpha <- function(alpha) {
  alpha <- as.numeric(alpha)
  if (length(alpha) != 1L || is.na(alpha) || !is.finite(alpha) ||
      alpha <= 0 || alpha >= 1) {
    stop("`alpha` must be one finite number strictly between zero and one.",
         call. = FALSE)
  }
  alpha
}


.yzf_stability_message <- function(full.stable, pair.stable = TRUE,
                                   pair.failures = 0L, max.iter,
                                   context, strict) {
  messages <- character()
  if (!isTRUE(full.stable)) {
    messages <- c(
      messages,
      sprintf("the full-sample weighted fit did not stabilize in %d updates",
              max.iter)
    )
  }
  if (!isTRUE(pair.stable)) {
    messages <- c(
      messages,
      sprintf("%d leave-two-out diagonal fit(s) did not stabilize in %d updates",
              as.integer(pair.failures), max.iter)
    )
  }
  if (!length(messages)) {
    return(invisible(NULL))
  }
  message <- sprintf(
    "%s: %s. Last iterates are available with complete diagnostics; no ridge, floor, or iterate substitution was applied.",
    context, paste(messages, collapse = "; ")
  )
  if (isTRUE(strict)) {
    stop(message, call. = FALSE)
  }
  warning(message, call. = FALSE)
  invisible(NULL)
}


.yzf_gumbel_tail <- function(statistic) {
  log.intensity <- -0.5 * log(pi) - statistic / 2
  if (log.intensity > log(.Machine$double.xmax)) {
    return(list(
      p.value = 1,
      log.p.value = 0,
      cdf = 0,
      log.cdf = -Inf,
      log.intensity = log.intensity
    ))
  }
  intensity <- exp(log.intensity)
  if (intensity == 0) {
    return(list(
      p.value = 0,
      log.p.value = log.intensity,
      cdf = 1,
      log.cdf = 0,
      log.intensity = log.intensity
    ))
  }
  p.value <- -expm1(-intensity)
  list(
    p.value = p.value,
    log.p.value = if (p.value > 0) log(p.value) else log.intensity,
    cdf = exp(-intensity),
    log.cdf = -intensity,
    log.intensity = log.intensity
  )
}


.yzf_cot_component <- function(p.value, log.p.value,
                               log.one.minus.p.value) {
  if (p.value == 0.5) {
    return(list(sign = 0, log.abs = -Inf))
  }
  if (p.value > 1e-8 && p.value < 1 - 1e-8) {
    value <- tan(pi * (0.5 - p.value))
    return(list(sign = sign(value), log.abs = log(abs(value))))
  }
  if (p.value < 0.5) {
    return(list(sign = 1, log.abs = -log(pi) - log.p.value))
  }
  list(sign = -1, log.abs = -log(pi) - log.one.minus.p.value)
}


.yzf_signed_log_average <- function(first, second) {
  log.two <- log(2)
  first$log.abs <- first$log.abs - log.two
  second$log.abs <- second$log.abs - log.two
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
  if (identical(first$log.abs, second$log.abs)) {
    return(list(sign = 0, log.abs = -Inf))
  }
  larger <- if (first$log.abs > second$log.abs) first else second
  smaller <- if (first$log.abs > second$log.abs) second else first
  difference <- exp(smaller$log.abs - larger$log.abs)
  if (difference == 1) {
    return(list(sign = 0, log.abs = -Inf))
  }
  list(
    sign = larger$sign,
    log.abs = larger$log.abs + log1p(-difference)
  )
}


.yzf_cauchy_combine <- function(p.max, log.p.max, log.q.max,
                                p.sum, log.p.sum, log.q.sum) {
  reverse.endpoints <-
    (p.max == 0 && p.sum == 1) || (p.max == 1 && p.sum == 0)
  if (isTRUE(reverse.endpoints)) {
    stop(
      paste0(
        "The Cauchy combination is indeterminate for exact reverse ",
        "component endpoints 0 and 1 (+Inf - Inf)."
      ),
      call. = FALSE
    )
  }
  combined <- .yzf_signed_log_average(
    .yzf_cot_component(p.max, log.p.max, log.q.max),
    .yzf_cot_component(p.sum, log.p.sum, log.q.sum)
  )
  if (combined$sign == 0) {
    return(list(
      statistic = 0,
      statistic.sign = 0,
      log.absolute.statistic = -Inf,
      angle = 0,
      p.value = 0.5
    ))
  }
  inverse.magnitude <- if (
    -combined$log.abs > log(.Machine$double.xmax)
  ) Inf else exp(-combined$log.abs)
  complement.angle <- atan(inverse.magnitude)
  angle <- combined$sign * (pi / 2 - complement.angle)
  p.value <- 0.5 - angle / pi
  if (!is.finite(p.value) || p.value < 0 || p.value > 1) {
    stop("Stable Cauchy combination produced an invalid probability.",
         call. = FALSE)
  }
  statistic <- if (
    combined$log.abs > log(.Machine$double.xmax)
  ) combined$sign * Inf else combined$sign * exp(combined$log.abs)
  list(
    statistic = statistic,
    statistic.sign = combined$sign,
    log.absolute.statistic = combined$log.abs,
    angle = angle,
    p.value = p.value
  )
}


.yzf_full_components <- function(fit, variable.names) {
  direction <- as.matrix(fit$full_direction)
  colnames(direction) <- variable.names
  if (!is.null(rownames(direction))) {
    rownames(direction) <- rownames(direction)
  }
  list(
    weighted.location = stats::setNames(
      as.numeric(fit$location), variable.names
    ),
    weighted.location.standardized = stats::setNames(
      as.numeric(fit$location_standardized), variable.names
    ),
    scale.diagonal = stats::setNames(
      as.numeric(fit$diagonal_input_canonical), variable.names
    ),
    log.scale.diagonal = stats::setNames(
      as.numeric(fit$log_diagonal_input_canonical), variable.names
    ),
    scale.diagonal.standardized = stats::setNames(
      as.numeric(fit$diagonal_standardized), variable.names
    ),
    radius = as.numeric(fit$full_radius),
    radial.weight = as.numeric(fit$full_weight),
    direction = direction,
    fit = fit$fit_diagnostics
  )
}


.yzf_fit_diagnostics <- function(fit, controls, m, strict) {
  list(
    iteration.stable = isTRUE(fit$iteration_stable),
    iterations = as.integer(fit$fit_diagnostics$iterations),
    relative.update = as.numeric(fit$fit_diagnostics$relative.update),
    location.relative.update = as.numeric(
      fit$fit_diagnostics$location.relative.update
    ),
    log.diagonal.relative.update = as.numeric(
      fit$fit_diagnostics$log.diagonal.relative.update
    ),
    score.residual = as.numeric(fit$fit_diagnostics$score.residual),
    location.score.residual = as.numeric(
      fit$fit_diagnostics$location.score.residual
    ),
    diagonal.score.residual = as.numeric(
      fit$fit_diagnostics$diagonal.score.residual
    ),
    minimum.residual.distance = as.numeric(
      fit$fit_diagnostics$minimum.residual.distance
    ),
    zero.residual.count = as.integer(
      fit$fit_diagnostics$zero.residual.count
    ),
    convergence.basis = fit$fit_diagnostics$convergence.basis,
    score.residual.role = paste(
      "weighted estimating-equation diagnostic normalized by total radial",
      "weight; reported separately from iterate stability"
    ),
    m = m,
    tolerance = controls$tol,
    max.iterations.allowed = controls$max_iter,
    zero.tolerance = controls$zero_tol,
    strict = strict,
    estimating.equations = paste(
      "sum r^m U = 0 and p * diag(mean(U U')) = I"
    ),
    corrected.location.update = paste(
      "theta <- theta + D^(1/2) sum(r^m U) / sum(r^(m-1));",
      "the D^(-1/2) premultiplier printed in arXiv:2501.14168 is",
      "dimensionally inconsistent"
    ),
    corrected.scale.update = paste(
      "D <- p D^(1/2) diag(mean(U U')) D^(1/2);",
      "the final D^(-1/2) printed in arXiv:2501.14168 is a typo"
    ),
    scale.identification = paste(
      "the unnormalised internal diagonal retains the common scale of the",
      "marginal-variance initialization; displayed input-coordinate",
      "diagonals are canonicalized to largest entry one"
    ),
    U.at.zero = "U(0) = 0; negative radial powers remain undefined",
    internal.column.log.scale = as.numeric(
      fit$column_log_residual_scale
    ),
    subtraction.overflow.fallback.columns = as.numeric(
      fit$subtraction_overflow_columns
    ),
    regularization = "none",
    radius.perturbation = "none",
    weight.cap = "none",
    variance.repair = "none"
  )
}


#' Weighted scaled spatial median and diagonal HR scale
#'
#' Fits the weighted diagonal Hettmansperger--Randles equations used by Yan,
#' Zhao, and Feng.  With
#' \eqn{e_i=D^{-1/2}(X_i-\theta)}, \eqn{r_i=\lVert e_i\rVert}, and
#' \eqn{U_i=U(e_i)}, the fitted pair solves
#' \deqn{\sum_i r_i^m U_i=0,\qquad
#' p\,\operatorname{diag}\{n^{-1}\sum_iU_iU_i^{\mathsf T}\}=I_p.}
#' The default \eqn{m=-1} is the inverse-norm weighted estimator; the paper's
#' theoretical family permits any finite \eqn{m\leq 1}.
#'
#' The implementation uses the dimensionally coherent iterations
#' \deqn{\theta^+=\theta+D^{1/2}
#' \frac{\sum_i r_i^mU_i}{\sum_i r_i^{m-1}},\qquad
#' D^+=pD^{1/2}\operatorname{diag}\{n^{-1}\sum_iU_iU_i^{\mathsf T}\}D^{1/2}.}
#' These correct two incompatible \eqn{D^{-1/2}} factors printed in the
#' current arXiv source.  The estimating equations, units, and Bahadur
#' representation all require \eqn{D^{1/2}} in those positions.
#'
#' @param x A numeric matrix or data frame with observations in rows and at
#'   least two rows.
#' @param m Finite radial power no greater than one.  The default is `-1`.
#' @param tol Positive relative iteration tolerance.
#' @param max_iter Positive maximum number of iterations.
#' @param zero_tol Non-negative threshold at which a standardized residual is
#'   treated as zero.  The default detects exact zeros only.
#' @param strict If `TRUE`, non-stabilization is an error.  If `FALSE`, the
#'   final iterate is returned with a warning and full diagnostics.
#'
#' @return A list of class `weighted_scaled_spatial_median` containing the
#'   fitted location, canonical diagonal scale, standardized radii and
#'   directions, radial weights, and convergence diagnostics.
#'
#' @references
#' Yan, G., Zhao, P., and Feng, L. (2025). Inverse norm weighted maxsum test
#' for high dimensional location parameters. arXiv:2501.14168.
#' \url{https://arxiv.org/abs/2501.14168}.
#'
#' @examples
#' x <- matrix(c(-2, 1, 0, 3, -1, 2, 1, -3, 2, 0, 4, -2), ncol = 2)
#' weighted_scaled_spatial_median(x, m = 1, tol = 1e-6)
#'
#' @export
weighted_scaled_spatial_median <- function(x, m = -1, tol = 1e-8,
                                           max_iter = 500L,
                                           zero_tol = 0,
                                           strict = TRUE) {
  call <- match.call()
  x <- .as_data_matrix(x, "x", min_rows = 2L)
  m <- .yzf_validate_m(m)
  strict <- .yzf_validate_strict(strict)
  controls <- .validate_iteration_controls(tol, max_iter, zero_tol)
  fit <- cpp_weighted_scaled_spatial_median(
    x, m, controls$tol, controls$max_iter, controls$zero_tol
  )
  .yzf_stability_message(
    fit$iteration_stable, max.iter = controls$max_iter,
    context = "Weighted scaled spatial median", strict = strict
  )
  variable.names <- .hotelling_variable_names(x)
  components <- .yzf_full_components(fit, variable.names)
  diagnostics <- .yzf_fit_diagnostics(fit, controls, m, strict)
  names(diagnostics$internal.column.log.scale) <- variable.names
  result <- c(
    components,
    list(
      m = m,
      diagnostics = diagnostics,
      n = nrow(x),
      p = ncol(x),
      call = call
    )
  )
  class(result) <- c("weighted_scaled_spatial_median", "list")
  result
}


#' Yan--Zhao--Feng weighted max test
#'
#' Tests \eqn{H_0:\theta=\mu_0} with the weighted max statistic of Yan,
#' Zhao, and Feng.  The full-sample weighted estimator is obtained with
#' `weighted_scaled_spatial_median()` estimating equations.  Define
#' \eqn{\widehat\zeta_k=n^{-1}\sum_i\widehat r_i^k}.  The uncentred and
#' centred statistics are
#' \deqn{M_m=n\lVert\widehat D^{-1/2}(\widehat\theta-\mu_0)\rVert_\infty^2
#' p(1-n^{-1/2})\widehat\zeta_{m-1}^2/\widehat\zeta_{2m},}
#' \deqn{T_{MAX}^{(m)}=M_m-2\log p+\log\log p.}
#' Its feasible null cdf is
#' \eqn{F(t)=\exp\{-\pi^{-1/2}\exp(-t/2)\}}.  The returned upper-tail
#' p-value is evaluated with stable exponential tails.
#'
#' @inheritParams weighted_scaled_spatial_median
#' @param mu A finite null-location vector.  The default is zero.
#' @param alpha Significance level for the reported upper-tail rejection rule.
#'
#' @return An object of class `c("hd_location_test", "htest")`.  Raw
#'   components include both max statistics, radial moments and their logs,
#'   the complete fitted estimator, and convergence/no-repair diagnostics.
#'
#' @references
#' Yan, G., Zhao, P., and Feng, L. (2025). Inverse norm weighted maxsum test
#' for high dimensional location parameters. arXiv:2501.14168.
#' \url{https://arxiv.org/abs/2501.14168}.
#'
#' @examples
#' set.seed(2501)
#' x <- matrix(stats::rnorm(80), 10, 8)
#' yan_zhao_feng_weighted_max_test(x, m = 1, tol = 1e-6)
#'
#' @export
yan_zhao_feng_weighted_max_test <- function(
    x, mu = NULL, m = -1, alpha = 0.05, tol = 1e-8,
    max_iter = 500L, zero_tol = 0, strict = TRUE) {
  call <- match.call()
  data.name <- deparse1(substitute(x))
  x <- .as_data_matrix(x, "x", min_rows = 2L)
  n <- nrow(x)
  p <- ncol(x)
  if (p < 2L) {
    stop("The weighted max Gumbel calibration requires at least two variables.",
         call. = FALSE)
  }
  if (is.null(mu)) mu <- numeric(p) else mu <- .as_location(mu, p, "mu")
  m <- .yzf_validate_m(m)
  alpha <- .yzf_validate_alpha(alpha)
  strict <- .yzf_validate_strict(strict)
  controls <- .validate_iteration_controls(tol, max_iter, zero_tol)
  fit <- cpp_yzf_weighted_max(
    x, mu, m, controls$tol, controls$max_iter, controls$zero_tol
  )
  .yzf_stability_message(
    fit$iteration_stable, max.iter = controls$max_iter,
    context = "Yan--Zhao--Feng weighted max test", strict = strict
  )
  tail <- .yzf_gumbel_tail(as.numeric(fit$max_centred))
  variable.names <- .hotelling_variable_names(x)
  full <- .yzf_full_components(fit, variable.names)
  mu.named <- stats::setNames(as.numeric(mu), variable.names)
  fit.diagnostics <- .yzf_fit_diagnostics(fit, controls, m, strict)
  names(fit.diagnostics$internal.column.log.scale) <- variable.names
  critical.value <- -log(pi) - 2 * log(-log1p(-alpha))
  reject <- isTRUE(as.numeric(fit$max_centred) > critical.value)

  .new_hd_location_test(
    statistic = c(T.MAX = as.numeric(fit$max_centred)),
    parameter = c(m = m, p = p),
    p.value = tail$p.value,
    method = "Yan-Zhao-Feng inverse-norm weighted max test",
    data.name = data.name,
    alternative = "two.sided",
    raw.statistic = c(M.MAX = as.numeric(fit$max_raw)),
    estimate = full$weighted.location,
    null.value = mu.named,
    null.distribution = list(
      family = "Gumbel exp(-pi^(-1/2) * exp(-t/2))",
      exact = FALSE,
      tail = "upper",
      finite.sample.correction = "1 - n^(-1/2)",
      assumptions = paste(
        "Yan-Zhao-Feng high-dimensional weak-dependence asymptotics;",
        "m <= 1 and required radial moments exist; stable weighted fit"
      )
    ),
    components = c(
      list(
        M.MAX = as.numeric(fit$max_raw),
        T.MAX = as.numeric(fit$max_centred),
        p.MAX = tail$p.value,
        log.p.MAX = tail$log.p.value,
        log.F.MAX = tail$log.cdf,
        maximum.standardized.location = as.numeric(
          fit$maximum_standardized_location
        ),
        zeta.m.minus.1.hat = as.numeric(fit$zeta_m_minus_1),
        zeta.2m.hat = as.numeric(fit$zeta_2m),
        log.zeta.m.minus.1.hat = as.numeric(
          fit$log_zeta_m_minus_1
        ),
        log.zeta.2m.hat = as.numeric(fit$log_zeta_2m),
        moment.ratio.hat = as.numeric(fit$moment_ratio),
        log.moment.ratio.hat = as.numeric(fit$log_moment_ratio),
        finite.sample.correction = as.numeric(
          fit$finite_sample_correction
        ),
        null.location = mu.named,
        sample.mean = stats::setNames(
          as.numeric(fit$sample_mean), variable.names
        )
      ),
      full
    ),
    diagnostics = c(
      fit.diagnostics,
      list(
        calibration = paste(
          "feasible Gumbel calibration using full-sample radial moments;",
          "no bootstrap or simulation calibration"
        ),
        statistic.scale = paste(
          "T.MAX is centred by -2 log(p) + log(log(p)); M.MAX is",
          "the uncentred statistic"
        ),
        moment.convention = "zeta_k = mean(r^k)",
        applicability = paste(
          "asymptotic high-dimensional p,n growth and the paper's",
          "weak-coordinate-dependence conditions are required"
        ),
        rejection = list(
          alpha = alpha,
          critical.value = critical.value,
          rule = "reject when T.MAX > -log(pi) - 2*log(-log(1-alpha))",
          reject = reject
        )
      )
    ),
    n = n,
    p = p,
    call = call
  )
}


#' Yan--Zhao--Feng weighted max-sum test
#'
#' Combines the weighted max test with the general weighted sum statistic
#' \deqn{T_{SUM}^{(m)}=\frac{2}{n(n-1)}\sum_{i<j}
#' r_{ij,i}^m r_{ij,j}^m U_{ij,i}^{\mathsf T}U_{ij,j}.}
#' Here \eqn{D_{ij}} is the *unweighted* diagonal HR scale fitted after
#' deleting observations \eqn{i,j}; endpoints are centred at the null
#' \eqn{\mu_0}, not at the fitted leave-two-out location.
#'
#' The operational sum calibration is the direct feasible estimator
#' \deqn{\widehat\sigma_m^2=2n^{-4}\sum_{i\ne j}r_{ij,i}^{2m}r_{ij,j}^{2m}
#' \{(U_{ij,i}-\widetilde\mu_{ij})^{\mathsf T}U_{ij,j}\}
#' \{(U_{ij,j}-\widetilde\mu_{ij})^{\mathsf T}U_{ij,i}\},}
#' where
#' \eqn{\widetilde\mu_{ij}=(n-2)^{-1}\sum_{k\ne i,j}U_{ij,k}} is the
#' unweighted, null-centred leave-two sign mean established in Section S.3 of
#' the official Feng--Liu--Ma supplement.  The sum and max upper-tail p-values
#' are combined with equal Cauchy weights.  At \eqn{m=-1}, every sum component
#' reduces to `inst_one_sample_test()`.
#'
#' @inheritParams yan_zhao_feng_weighted_max_test
#'
#' @return An object of class `c("hd_location_test", "htest")`.  Its
#'   `components` field includes max, sum, direct variance, Cauchy, every
#'   leave-two kernel, fitted scales and locations, and radial quantities.
#'   `diagnostics` records all full/leave-out updates, residuals, zero counts,
#'   asymptotic scope, and the no-repair contract.
#'
#' @references
#' Yan, G., Zhao, P., and Feng, L. (2025). Inverse norm weighted maxsum test
#' for high dimensional location parameters. arXiv:2501.14168.
#' \url{https://arxiv.org/abs/2501.14168}.
#'
#' Feng, L., Liu, B., and Ma, Y. (2020). Supplementary material for *An
#' inverse norm sign test of location parameter for high-dimensional data*.
#' \doi{10.6084/m9.figshare.11914095.v2}.
#'
#' @examples
#' set.seed(2502)
#' x <- matrix(stats::rnorm(80), 10, 8)
#' \donttest{
#' yan_zhao_feng_weighted_maxsum_test(x, m = 0, tol = 1e-6)
#' }
#'
#' @export
yan_zhao_feng_weighted_maxsum_test <- function(
    x, mu = NULL, m = -1, alpha = 0.05, tol = 1e-8,
    max_iter = 500L, zero_tol = 0, strict = TRUE) {
  call <- match.call()
  data.name <- deparse1(substitute(x))
  x <- .as_data_matrix(x, "x", min_rows = 4L)
  n <- nrow(x)
  p <- ncol(x)
  if (p < 2L) {
    stop("The weighted max-sum Gumbel calibration requires at least two variables.",
         call. = FALSE)
  }
  if (is.null(mu)) mu <- numeric(p) else mu <- .as_location(mu, p, "mu")
  m <- .yzf_validate_m(m)
  alpha <- .yzf_validate_alpha(alpha)
  strict <- .yzf_validate_strict(strict)
  controls <- .validate_iteration_controls(tol, max_iter, zero_tol)
  fit <- cpp_yzf_weighted_maxsum(
    x, mu, m, controls$tol, controls$max_iter, controls$zero_tol
  )
  sum.fit <- fit$sum
  .yzf_stability_message(
    fit$iteration_stable,
    pair.stable = sum.fit$all_iteration_stable,
    pair.failures = sum.fit$stability_failures,
    max.iter = controls$max_iter,
    context = "Yan--Zhao--Feng weighted max-sum test",
    strict = strict
  )

  max.tail <- .yzf_gumbel_tail(as.numeric(fit$max_centred))
  z.sum <- as.numeric(sum.fit$z)
  p.sum <- stats::pnorm(z.sum, lower.tail = FALSE)
  log.p.sum <- stats::pnorm(z.sum, lower.tail = FALSE, log.p = TRUE)
  log.q.sum <- stats::pnorm(z.sum, lower.tail = TRUE, log.p = TRUE)
  combined <- .yzf_cauchy_combine(
    max.tail$p.value, max.tail$log.p.value, max.tail$log.cdf,
    p.sum, log.p.sum, log.q.sum
  )

  variable.names <- .hotelling_variable_names(x)
  full <- .yzf_full_components(fit, variable.names)
  mu.named <- stats::setNames(as.numeric(mu), variable.names)
  fit.diagnostics <- .yzf_fit_diagnostics(fit, controls, m, strict)
  names(fit.diagnostics$internal.column.log.scale) <- variable.names
  pair.labels <- paste0(sum.fit$pair_i, ",", sum.fit$pair_j)
  pair.index <- data.frame(
    i = as.integer(sum.fit$pair_i),
    j = as.integer(sum.fit$pair_j),
    row.names = pair.labels
  )
  name.pairs <- function(value) {
    stats::setNames(as.numeric(value), pair.labels)
  }
  name.pair.matrix <- function(value, columns = variable.names) {
    value <- as.matrix(value)
    dimnames(value) <- list(pair.labels, columns)
    value
  }
  endpoint.columns <- c("i", "j")
  factor.columns <- c(
    "(U_i-mu_tilde)'U_j", "(U_j-mu_tilde)'U_i"
  )
  pair.diagnostics <- sum.fit$fit_diagnostics
  for (field in c(
    "iterations", "iteration.stable", "relative.update",
    "location.relative.update", "log.diagonal.relative.update",
    "score.residual", "location.score.residual",
    "diagonal.score.residual", "minimum.residual.distance",
    "zero.residual.count"
  )) {
    names(pair.diagnostics[[field]]) <- pair.labels
  }
  critical.max <- -log(pi) - 2 * log(-log1p(-alpha))
  critical.sum <- stats::qnorm(1 - alpha)

  .new_hd_location_test(
    statistic = c(Cauchy.angle = combined$angle),
    parameter = c(m = m, p = p),
    p.value = combined$p.value,
    method = "Yan-Zhao-Feng inverse-norm weighted max-sum test",
    data.name = data.name,
    alternative = "two.sided",
    raw.statistic = c(
      M.MAX = as.numeric(fit$max_raw),
      T.MAX = as.numeric(fit$max_centred),
      T.SUM = as.numeric(sum.fit$T),
      Cauchy = combined$statistic
    ),
    estimate = full$weighted.location,
    null.value = mu.named,
    null.distribution = list(
      family = "equal-weight Cauchy combination of Gumbel and normal limits",
      exact = FALSE,
      tail = "combined upper tails",
      assumptions = paste(
        "Yan-Zhao-Feng asymptotic independence and high-dimensional",
        "weak-dependence conditions; required radial moments and stable fits"
      )
    ),
    variance = c(sum.primary.direct = as.numeric(sum.fit$variance)),
    components = c(
      list(
        M.MAX = as.numeric(fit$max_raw),
        T.MAX = as.numeric(fit$max_centred),
        p.MAX = max.tail$p.value,
        log.p.MAX = max.tail$log.p.value,
        T.SUM = as.numeric(sum.fit$T),
        Z.SUM = z.sum,
        p.SUM = p.sum,
        log.p.SUM = log.p.sum,
        sigma2.direct.hat = as.numeric(sum.fit$variance),
        sigma.direct.hat = as.numeric(sum.fit$standard_error),
        Cauchy.statistic = combined$statistic,
        Cauchy.statistic.sign = combined$statistic.sign,
        log.absolute.Cauchy.statistic = combined$log.absolute.statistic,
        Cauchy.angle = combined$angle,
        p.CC = combined$p.value,
        zeta.m.minus.1.hat = as.numeric(fit$zeta_m_minus_1),
        zeta.2m.hat = as.numeric(fit$zeta_2m),
        log.zeta.m.minus.1.hat = as.numeric(
          fit$log_zeta_m_minus_1
        ),
        log.zeta.2m.hat = as.numeric(fit$log_zeta_2m),
        moment.ratio.hat = as.numeric(fit$moment_ratio),
        log.moment.ratio.hat = as.numeric(fit$log_moment_ratio),
        finite.sample.correction = as.numeric(
          fit$finite_sample_correction
        ),
        maximum.standardized.location = as.numeric(
          fit$maximum_standardized_location
        ),
        ordered.variance.kernel.sum = as.numeric(
          sum.fit$ordered_variance_kernel_sum
        ),
        zeta.2m.crossfit.hat = as.numeric(
          sum.fit$zeta_2m_crossfit
        ),
        trace.R2.hat = as.numeric(sum.fit$trace_R2_hat),
        sigma2.oracle.pair.adjusted.diagnostic = as.numeric(
          sum.fit$oracle_variance_pair_adjusted
        ),
        sigma2.oracle.asymptotic.diagnostic = as.numeric(
          sum.fit$oracle_variance_asymptotic
        ),
        trace.weighted.score.covariance.squared.hat = as.numeric(
          sum.fit$weighted_score_trace_hat
        ),
        pair.index = pair.index,
        pair.test.inner.product = name.pairs(
          sum.fit$pair_test_inner_product
        ),
        pair.statistic.kernel = name.pairs(
          sum.fit$pair_statistic_kernel
        ),
        pair.leaveout.sign.mean.norm = name.pairs(
          sum.fit$pair_sign_mean_norm
        ),
        pair.variance.factor = name.pair.matrix(
          sum.fit$pair_variance_factor, factor.columns
        ),
        pair.variance.kernel = name.pairs(
          sum.fit$pair_variance_kernel
        ),
        pair.trace.inner.product.squared = name.pairs(
          sum.fit$pair_trace_inner_product_squared
        ),
        pair.trace.zero.signs = name.pairs(
          sum.fit$pair_trace_zero_signs
        ),
        pair.null.zero.signs = name.pairs(
          sum.fit$pair_null_zero_signs
        ),
        endpoint.radius = name.pair.matrix(
          sum.fit$endpoint_radius, endpoint.columns
        ),
        endpoint.radial.weight = name.pair.matrix(
          sum.fit$endpoint_weight, endpoint.columns
        ),
        leaveout.location = name.pair.matrix(
          sum.fit$leaveout_location
        ),
        leaveout.location.standardized = name.pair.matrix(
          sum.fit$leaveout_location_standardized
        ),
        leaveout.scale.diagonal = name.pair.matrix(
          sum.fit$leaveout_diagonal_input_canonical
        ),
        leaveout.log.scale.diagonal = name.pair.matrix(
          sum.fit$leaveout_log_diagonal_input_canonical
        ),
        leaveout.scale.diagonal.standardized = name.pair.matrix(
          sum.fit$leaveout_diagonal_standardized
        ),
        pair.count = as.numeric(sum.fit$pair_count),
        ordered.pair.count = as.numeric(sum.fit$ordered_pair_count),
        null.location = mu.named,
        sample.mean = stats::setNames(
          as.numeric(fit$sample_mean), variable.names
        )
      ),
      full
    ),
    diagnostics = c(
      fit.diagnostics,
      list(
        pair.fit = pair.diagnostics,
        pair.iteration.stable = isTRUE(
          sum.fit$all_iteration_stable
        ),
        pair.stability.failures = as.integer(
          sum.fit$stability_failures
        ),
        endpoint.center = paste(
          "null location mu only; the leave-two-out fitted location",
          "participates only in estimating D_ij"
        ),
        variance.sign.center = paste(
          "unweighted mean of null-centred U_ij,k over k != i,j;",
          "not radially weighted and not the fitted leave-two location"
        ),
        primary.sum.variance = paste(
          "Feng-Liu-Ma official supplement S.3 generalized by r^(2m):",
          "2*n^(-4) times an ordered i != j sum"
        ),
        primary.sum.calibration = "direct feasible variance only",
        oracle.factorized.variance.role = "diagnostic only",
        leaveout.fit.weight = "unweighted diagonal HR (m = 0)",
        pair.convention = paste(
          "unordered leave-two fits; direct variance reconstructed as",
          "the published ordered-pair sum"
        ),
        cauchy.evaluation = paste(
          "signed-log cotangent average with an angle representation;",
          "no p-value clipping"
        ),
        calibration = paste(
          "Gumbel max plus direct-feasible normal sum, combined by equal",
          "Cauchy weights; no bootstrap or simulation calibration"
        ),
        applicability = paste(
          "asymptotic high-dimensional p,n growth, required radial moments,",
          "and Yan-Zhao-Feng weak-dependence/independence conditions"
        ),
        rejection = list(
          alpha = alpha,
          rule = "reject when p.CC < alpha",
          reject = isTRUE(combined$p.value < alpha),
          component.rules = c(
            max = sprintf("T.MAX > %.17g", critical.max),
            sum = sprintf("Z.SUM > %.17g", critical.sum)
          )
        )
      )
    ),
    n = n,
    p = p,
    call = call
  )
}
