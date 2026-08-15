.ssmax_validate_strict <- function(strict) {
  if (!is.logical(strict) || length(strict) != 1L || is.na(strict)) {
    stop("`strict` must be TRUE or FALSE.", call. = FALSE)
  }
  strict
}

.ssmax_validate_alpha <- function(alpha) {
  alpha <- as.numeric(alpha)
  if (length(alpha) != 1L || is.na(alpha) || !is.finite(alpha) ||
      alpha <= 0 || alpha >= 1) {
    stop("`alpha` must be one finite number strictly between zero and one.",
         call. = FALSE)
  }
  alpha
}

.ssmax_check_fit <- function(fit, strict, context) {
  if (isTRUE(fit$iteration_stable)) {
    return(invisible(fit))
  }
  message <- sprintf(
    paste0(
      "%s did not stabilize by relative location/log-diagonal update within ",
      "%d updates (score residual %.6g, relative update %.6g)."
    ),
    context, as.integer(fit$max_iterations),
    unname(fit$score_residual), unname(fit$relative_update)
  )
  if (strict) {
    stop(message, call. = FALSE)
  }
  warning(paste0(message, " Returning the last finite iterate."),
          call. = FALSE)
  invisible(fit)
}

.ssmax_gumbel_tail <- function(centered) {
  centered <- as.numeric(centered)
  if (length(centered) != 1L || is.na(centered) || !is.finite(centered)) {
    stop("The centered max statistic must be one finite number.",
         call. = FALSE)
  }

  # F(x) = exp{-a(x)}, a(x) = pi^(-1/2) exp(-x/2).
  log.intensity <- -0.5 * log(pi) - 0.5 * centered
  if (log.intensity > log(.Machine$double.xmax)) {
    return(list(
      p.value = 1,
      log.p.value = 0,
      log.one.minus.p.value = -Inf,
      log.intensity = log.intensity
    ))
  }

  intensity <- exp(log.intensity)
  if (intensity == 0) {
    return(list(
      p.value = 0,
      log.p.value = log.intensity,
      log.one.minus.p.value = 0,
      log.intensity = log.intensity
    ))
  }

  log.p <- if (intensity <= log(2)) {
    log(-expm1(-intensity))
  } else {
    log1p(-exp(-intensity))
  }
  list(
    p.value = exp(log.p),
    log.p.value = log.p,
    log.one.minus.p.value = -intensity,
    log.intensity = log.intensity
  )
}

.ssmax_gumbel_quantile <- function(probability) {
  probability <- as.numeric(probability)
  if (length(probability) != 1L || is.na(probability) ||
      !is.finite(probability) || probability <= 0 || probability >= 1) {
    stop("`probability` must be one finite number in (0, 1).",
         call. = FALSE)
  }
  -log(pi) - 2 * log(-log(probability))
}

.ssmax_signed_log_cot <- function(log.p, log.one.minus.p) {
  values <- c(log.p, log.one.minus.p)
  if (length(log.p) != 1L || length(log.one.minus.p) != 1L ||
      anyNA(values) || any(values > 0) ||
      any(!is.finite(values) & values != -Inf)) {
    stop("Cauchy component log tails must be non-positive or -Inf.",
         call. = FALSE)
  }
  p.is.zero <- is.infinite(log.p) && log.p < 0
  q.is.zero <- is.infinite(log.one.minus.p) && log.one.minus.p < 0
  if (p.is.zero && q.is.zero) {
    stop("A Cauchy component cannot have both probability tails equal to zero.",
         call. = FALSE)
  }

  if (identical(log.p, log.one.minus.p)) {
    return(list(sign = 0, log.absolute = -Inf))
  }
  sign <- if (log.p < log.one.minus.p) 1 else -1
  log.u <- min(log.p, log.one.minus.p)
  if (is.infinite(log.u) && log.u < 0) {
    return(list(sign = sign, log.absolute = Inf))
  }

  # For a very small tail, log(cot(pi*u)) = -log(pi)-log(u)+O(u^2).
  # Avoid materialising u when only its log remains representable.
  log.cot <- if (log.u < log(1e-8)) {
    -log(pi) - log.u
  } else {
    u <- exp(log.u)
    log(cospi(u)) - log(sinpi(u))
  }
  list(sign = sign, log.absolute = log.cot - log(2))
}

.ssmax_cauchy_combine_logtails <- function(log.p1, log.q1,
                                            log.p2, log.q2) {
  is.zero.tail <- function(value) {
    length(value) == 1L && !is.na(value) &&
      is.infinite(value) && value < 0
  }
  reverse.endpoints <-
    (is.zero.tail(log.p1) && is.zero.tail(log.q2)) ||
    (is.zero.tail(log.q1) && is.zero.tail(log.p2))
  if (reverse.endpoints) {
    stop(
      paste0(
        "The Cauchy combination is indeterminate for exact reverse ",
        "component endpoints 0 and 1 (+Inf - Inf)."
      ),
      call. = FALSE
    )
  }

  first <- .ssmax_signed_log_cot(log.p1, log.q1)
  second <- .ssmax_signed_log_cot(log.p2, log.q2)
  if (first$sign == 0) {
    combined <- second
  } else if (second$sign == 0) {
    combined <- first
  } else if (first$sign == second$sign) {
    largest <- max(first$log.absolute, second$log.absolute)
    if (is.infinite(largest)) {
      combined <- list(sign = first$sign, log.absolute = largest)
    } else {
      combined <- list(
        sign = first$sign,
        log.absolute = largest + log1p(exp(
          min(first$log.absolute, second$log.absolute) - largest
        ))
      )
    }
  } else if (identical(first$log.absolute, second$log.absolute)) {
    combined <- list(sign = 0, log.absolute = -Inf)
  } else {
    first.is.larger <- first$log.absolute > second$log.absolute
    larger <- if (first.is.larger) first else second
    smaller <- if (first.is.larger) second else first
    combined <- list(
      sign = larger$sign,
      log.absolute = larger$log.absolute +
        log(-expm1(smaller$log.absolute - larger$log.absolute))
    )
  }

  if (combined$sign == 0) {
    p.value <- 0.5
    log.p.value <- log(0.5)
    log.one.minus.p.value <- log(0.5)
    raw.statistic <- 0
  } else if (is.infinite(combined$log.absolute)) {
    p.value <- if (combined$sign > 0) 0 else 1
    log.p.value <- if (combined$sign > 0) -Inf else 0
    log.one.minus.p.value <- if (combined$sign > 0) 0 else -Inf
    raw.statistic <- combined$sign * Inf
  } else {
    magnitude.scale <- max(0, combined$log.absolute)
    y <- exp(-magnitude.scale)
    x <- combined$sign * exp(combined$log.absolute - magnitude.scale)
    p.value <- atan2(y, x) / pi
    if (combined$log.absolute > 20) {
      small.log.tail <- -combined$log.absolute - log(pi)
      if (combined$sign > 0) {
        log.p.value <- small.log.tail
        log.one.minus.p.value <- log1p(-exp(small.log.tail))
      } else {
        log.one.minus.p.value <- small.log.tail
        log.p.value <- log1p(-exp(small.log.tail))
      }
    } else {
      log.p.value <- log(p.value)
      log.one.minus.p.value <- log1p(-p.value)
    }
    raw.statistic <- if (combined$log.absolute >
                         log(.Machine$double.xmax)) {
      combined$sign * Inf
    } else {
      combined$sign * exp(combined$log.absolute)
    }
  }

  list(
    p.value = p.value,
    log.p.value = log.p.value,
    log.one.minus.p.value = log.one.minus.p.value,
    statistic = raw.statistic,
    angle = pi * (0.5 - p.value),
    sign = combined$sign,
    log.absolute = combined$log.absolute,
    component.signed.log = list(first = first, second = second)
  )
}

.ssmax_name_fit <- function(fit, variable.names, observation.names) {
  name.vector <- function(value) {
    stats::setNames(as.numeric(value), variable.names)
  }
  directions <- as.matrix(fit$directions)
  dimnames(directions) <- list(observation.names, variable.names)
  fit$location <- name.vector(fit$location)
  fit$location_relative_origin <- name.vector(fit$location_relative_origin)
  fit$location_standardized <- name.vector(fit$location_standardized)
  fit$standardized_location_coordinate <- name.vector(
    fit$standardized_location_coordinate
  )
  fit$origin <- name.vector(fit$origin)
  fit$sample_mean <- name.vector(fit$sample_mean)
  fit$diagonal_standardized <- name.vector(fit$diagonal_standardized)
  fit$diagonal_input_canonical <- name.vector(
    fit$diagonal_input_canonical
  )
  fit$log_diagonal_input_canonical <- name.vector(
    fit$log_diagonal_input_canonical
  )
  fit$coordinate_statistic <- name.vector(fit$coordinate_statistic)
  fit$column_log_scale <- name.vector(fit$column_log_scale)
  fit$directions <- directions
  fit$radii <- stats::setNames(as.numeric(fit$radii), observation.names)
  fit$inverse_radii <- stats::setNames(
    as.numeric(fit$inverse_radii), observation.names
  )
  fit
}

.ssmax_fit_diagnostics <- function(fit, strict) {
  list(
    iteration.stable = isTRUE(fit$iteration_stable),
    iterations = as.integer(fit$iterations),
    tolerance = unname(fit$tolerance),
    max.iterations.allowed = as.integer(fit$max_iterations),
    relative.update = unname(fit$relative_update),
    location.relative.update = unname(fit$location_relative_update),
    log.diagonal.update = unname(fit$log_diagonal_update),
    score.residual = unname(fit$score_residual),
    location.equation.residual = unname(fit$location_score_residual),
    diagonal.equation.residual = unname(fit$diagonal_score_residual),
    minimum.residual.distance = unname(fit$minimum_residual_distance),
    convergence.basis =
      "relative location and log-diagonal iterate update no larger than tol",
    score.residual.role = paste(
      "the maximum full-sample mean-sign/diagonal-sign equation residual",
      "is reported separately from iterate stability"
    ),
    strict = strict,
    scale.identification = paste(
      "internal safely standardized D has unit geometric mean; returned",
      "input-coordinate D has maximum diagonal one because the HR system",
      "identifies only relative diagonal scale"
    ),
    radial.scale = paste(
      "radii and inverse-radius moments use the unit-geometric-mean D in",
      "the safely standardized coordinate system"
    ),
    internal.column.log.scale = fit$column_log_scale,
    subtraction.overflow.fallback.columns =
      unname(fit$subtraction_overflow_columns),
    zero.tolerance = unname(fit$zero_tolerance),
    regularization = "none",
    numerical.floor = "none",
    perturbation = "none",
    pseudoinverse = "none"
  )
}

#' Scaled spatial median with diagonal HR standardisation
#'
#' Fits the full-sample diagonal Hettmansperger--Randles system used by
#' Liu, Feng, Zhao and Wang:
#' \deqn{n^{-1}\sum_i U\{D^{-1/2}(X_i-\theta)\}=0,}
#' \deqn{(p/n)\operatorname{diag}\sum_i
#' U\{D^{-1/2}(X_i-\theta)\}U\{D^{-1/2}(X_i-\theta)\}^{\mathsf T}=I_p.}
#' Starting with the sample mean and marginal sample variances, the function
#' applies the paper's simultaneous location and diagonal-scale recursion.
#'
#' The equations identify `D` only up to a common positive multiplier.
#' Internally, `D` has unit geometric mean after safe coordinate
#' standardisation; the returned input-coordinate diagonal is divided by its
#' largest entry. Both identifications leave the location and signs unchanged.
#' Radial quantities are explicitly labelled as belonging to the internal
#' unit-geometric-mean identification. The original paper notes that general
#' existence, uniqueness, and convergence of this recursion are not proved.
#'
#' With `strict = TRUE`, failure of the relative location/log-diagonal update
#' to reach `tol` within `max_iter` updates is an error. With `strict = FALSE`,
#' the last finite iterate is returned with a warning and
#' `iteration.stable = FALSE`. Both equation residuals are reported separately
#' and are not mislabelled as the paper's unspecified stopping rule. Coincident
#' residuals, non-positive marginal variation, and non-finite updates are
#' errors; no ridge, floor, perturbation, or pseudoinverse is used.
#'
#' @param x A numeric matrix or data frame with observations in rows and
#'   variables in columns. At least two observations and positive marginal
#'   sample variation are required.
#' @param tol A finite positive estimating-equation tolerance.
#' @param max_iter A positive integer maximum number of recursion updates.
#' @param zero_tol A non-negative threshold for declaring an internally
#'   standardised residual radius singular. The literal default is zero.
#' @param strict Whether non-convergence is an error (`TRUE`) or warning with
#'   the last finite iterate returned (`FALSE`).
#'
#' @return A list containing `location`, the identified scale diagonal and its
#'   logarithm, spatial directions, radial quantities, and convergence and
#'   equation diagnostics.
#'
#' @references
#' Liu, B., Feng, L., Zhao, P. and Wang, Z. Spatial-sign based maxsum test for
#' high-dimensional location parameters. *Statistica Sinica*, accepted.
#' \doi{10.5705/ss.202024.0051}.
#'
#' @examples
#' x <- matrix(c(0.7, -0.4, 1.1, -1.2, 0.2, 1.5, -0.8,
#'               -1.1, 0.8, 0.3, -0.5, 1.4, -0.2, 0.6), 7, 2)
#' scaled_spatial_median(x, tol = 1e-6)
#'
#' @export
scaled_spatial_median <- function(x, tol = 1e-7, max_iter = 500L,
                                  zero_tol = 0, strict = TRUE) {
  call <- match.call()
  data.name <- deparse1(substitute(x))
  x <- .as_data_matrix(x, "x", min_rows = 2L)
  controls <- .validate_iteration_controls(tol, max_iter, zero_tol)
  strict <- .ssmax_validate_strict(strict)
  variable.names <- .hotelling_variable_names(x)
  observation.names <- if (is.null(rownames(x))) {
    paste0("observation", seq_len(nrow(x)))
  } else {
    rownames(x)
  }

  fit <- cpp_scaled_spatial_median(
    x, numeric(ncol(x)), FALSE, controls$tol, controls$max_iter,
    controls$zero_tol
  )
  .ssmax_check_fit(fit, strict, "The scaled spatial median recursion")
  fit <- .ssmax_name_fit(fit, variable.names, observation.names)

  result <- list(
    location = fit$location,
    scale.diagonal = fit$diagonal_input_canonical,
    log.scale.diagonal = fit$log_diagonal_input_canonical,
    scale.diagonal.standardized = fit$diagonal_standardized,
    location.standardized = fit$location_standardized,
    directions = fit$directions,
    radii.standardized = fit$radii,
    inverse.radii.standardized = fit$inverse_radii,
    inverse.radius.moment.standardized = unname(fit$zeta1_hat),
    sample.mean = fit$sample_mean,
    diagnostics = .ssmax_fit_diagnostics(fit, strict),
    n = nrow(x),
    p = ncol(x),
    data.name = data.name,
    call = call
  )
  class(result) <- c("scaled_spatial_median", "list")
  result
}

#' Spatial-sign max test for a high-dimensional location
#'
#' Tests \eqn{H_0:\theta=\mu} with the feasible max statistic of Liu, Feng,
#' Zhao and Wang. The full-sample scaled spatial median and diagonal HR scale
#' are fitted first. With
#' \eqn{\widehat r_i=\|\widehat D^{-1/2}
#' (X_i-\widehat\theta)\|} and
#' \eqn{\widehat\zeta_1=n^{-1}\sum_i\widehat r_i^{-1}},
#' \deqn{T_{\rm MAX}=n\|\widehat D^{-1/2}
#' (\widehat\theta-\mu)\|_\infty^2\,p\widehat\zeta_1^2
#' (1-n^{-1/2}).}
#' Every factor shown is multiplicative, including the unsquared finite-sample
#' factor. Under the paper's high-dimensional conditions,
#' \deqn{T_{\rm MAX}-2\log p+\log\log p}
#' has limiting cdf
#' \eqn{F(t)=\exp\{-\pi^{-1/2}\exp(-t/2)\}}. The vector alternative is
#' scientifically two-sided, while large max statistics use the upper tail.
#'
#' `strict`, singular residual handling, and diagonal-scale identification are
#' as described in [scaled_spatial_median()]. At least two variables are
#' required for the literal \eqn{\log\log p} centering.
#'
#' @inheritParams scaled_spatial_median
#' @param mu A finite numeric null-location vector. The default is zero.
#' @param alpha A finite test level strictly between zero and one.
#'
#' @return An object of class `c("hd_location_test", "htest")`. It contains
#'   the raw max statistic, centered Gumbel statistic, upper-tail probability
#'   and log tails, critical values, the full-sample fit, radial moment, and
#'   convergence diagnostics.
#'
#' @references
#' Liu, B., Feng, L., Zhao, P. and Wang, Z. Spatial-sign based maxsum test for
#' high-dimensional location parameters. *Statistica Sinica*, accepted.
#' \doi{10.5705/ss.202024.0051}.
#'
#' @examples
#' x <- matrix(c(-2, -1, 0, 1, 2, 3, -1, 2, 1, -2, 3, 0,
#'               1, 0, -1, 2, -2, 1), 6, 3)
#' spatial_sign_max_test(x, tol = 1e-6)
#'
#' @export
spatial_sign_max_test <- function(x, mu = NULL, alpha = 0.05,
                                  tol = 1e-7, max_iter = 500L,
                                  zero_tol = 0, strict = TRUE) {
  call <- match.call()
  data.name <- deparse1(substitute(x))
  x <- .as_data_matrix(x, "x", min_rows = 2L)
  n <- nrow(x)
  p <- ncol(x)
  if (p < 2L) {
    stop("The spatial-sign max test requires at least two variables.",
         call. = FALSE)
  }
  if (is.null(mu)) {
    mu <- numeric(p)
  } else {
    mu <- .as_location(mu, p, "mu")
  }
  alpha <- .ssmax_validate_alpha(alpha)
  controls <- .validate_iteration_controls(tol, max_iter, zero_tol)
  strict <- .ssmax_validate_strict(strict)
  variable.names <- .hotelling_variable_names(x)
  observation.names <- if (is.null(rownames(x))) {
    paste0("observation", seq_len(n))
  } else {
    rownames(x)
  }

  fit <- cpp_scaled_spatial_median(
    x, mu, TRUE, controls$tol, controls$max_iter, controls$zero_tol
  )
  .ssmax_check_fit(fit, strict, "The spatial-sign max fit")
  fit <- .ssmax_name_fit(fit, variable.names, observation.names)

  T.max <- unname(fit$T_max)
  centered <- T.max - 2 * log(p) + log(log(p))
  tail <- .ssmax_gumbel_tail(centered)
  quantile <- .ssmax_gumbel_quantile(1 - alpha)
  raw.critical <- 2 * log(p) - log(log(p)) + quantile
  variable.mu <- stats::setNames(as.numeric(mu), variable.names)
  diagnostics <- .ssmax_fit_diagnostics(fit, strict)
  diagnostics$tail <- "upper"
  diagnostics$alpha <- alpha
  diagnostics$rejected <- tail$p.value < alpha
  diagnostics$finite.sample.factor <- unname(fit$finite_sample_factor)
  diagnostics$gumbel.log.intensity <- tail$log.intensity

  .new_hd_location_test(
    statistic = c(Gumbel.centered = centered),
    p.value = tail$p.value,
    method = paste(
      "Liu-Feng-Zhao-Wang spatial-sign max test",
      "(full-sample diagonal HR fit; feasible Gumbel calibration)"
    ),
    data.name = data.name,
    alternative = "two.sided",
    raw.statistic = c(T.MAX = T.max),
    estimate = fit$location,
    null.value = variable.mu,
    null.distribution = list(
      family = "Gumbel-type extreme-value",
      cdf = "exp{-pi^(-1/2) exp(-x/2)}",
      exact = FALSE,
      tail = "upper",
      centered.quantile = c(`1-alpha` = quantile),
      raw.critical.value = c(`1-alpha` = raw.critical),
      assumptions = paste(
        "Liu-Feng-Zhao-Wang high-dimensional radial-moment, Bahadur",
        "remainder, growth, and weak-correlation conditions"
      )
    ),
    components = list(
      location = fit$location,
      location.minus.null = fit$location_relative_origin,
      null.location = variable.mu,
      sample.mean = fit$sample_mean,
      scale.diagonal = fit$diagonal_input_canonical,
      log.scale.diagonal = fit$log_diagonal_input_canonical,
      scale.diagonal.standardized = fit$diagonal_standardized,
      standardized.location.coordinate =
        fit$standardized_location_coordinate,
      coordinate.statistic = fit$coordinate_statistic,
      directions = fit$directions,
      radii.standardized = fit$radii,
      inverse.radii.standardized = fit$inverse_radii,
      zeta1.hat.standardized = unname(fit$zeta1_hat),
      finite.sample.factor = unname(fit$finite_sample_factor),
      T.MAX = T.max,
      centered.max = centered,
      p.MAX = tail$p.value,
      log.p.MAX = tail$log.p.value,
      log.one.minus.p.MAX = tail$log.one.minus.p.value,
      critical.value.centered = quantile,
      critical.value.raw = raw.critical,
      n = n,
      p = p
    ),
    diagnostics = diagnostics,
    n = n,
    p = p,
    call = call
  )
}

#' Spatial-sign max-sum test for a high-dimensional location
#'
#' Combines the feasible spatial-sign max p-value from
#' [spatial_sign_max_test()] with the feasible Feng--Sun sum p-value returned
#' directly by [feng_sun_one_sample_test()]. For component p-values
#' \eqn{p_{\rm MAX}} and \eqn{p_{\rm SUM}}, the published combination is
#' \deqn{1-G[0.5\tan\{\pi(0.5-p_{\rm MAX})\}+
#' 0.5\tan\{\pi(0.5-p_{\rm SUM})\}],}
#' where \eqn{G} is the standard Cauchy cdf.
#'
#' The implementation retains both component log tails and evaluates the
#' combination in signed-log form, finishing with `atan2`. It therefore does
#' not clip component probabilities. Exact equal endpoints map to the same
#' endpoint; exact reverse endpoints 0 and 1 are reported as mathematically
#' indeterminate rather than silently set to one half. The sum component uses
#' literal leave-two-out fits and consequently requires at least four
#' observations.
#'
#' @inheritParams spatial_sign_max_test
#'
#' @return An object of class `c("hd_location_test", "htest")` containing the
#'   combined p-value, both complete component-test objects, their ordinary
#'   and logarithmic p-value tails, and signed-log Cauchy diagnostics.
#'
#' @references
#' Liu, B., Feng, L., Zhao, P. and Wang, Z. Spatial-sign based maxsum test for
#' high-dimensional location parameters. *Statistica Sinica*, accepted.
#' \doi{10.5705/ss.202024.0051}.
#'
#' @examples
#' x <- matrix(c(-2, -1, 0, 1, 2, 3, 1, -1,
#'               -1, 2, 1, -2, 3, 0, -2, 1,
#'               1, 0, -1, 2, -2, 1, 3, -1), 8, 3)
#' spatial_sign_maxsum_test(x, tol = 1e-6)
#'
#' @export
spatial_sign_maxsum_test <- function(x, mu = NULL, alpha = 0.05,
                                     tol = 1e-7, max_iter = 500L,
                                     zero_tol = 0, strict = TRUE) {
  call <- match.call()
  data.name <- deparse1(substitute(x))
  x <- .as_data_matrix(x, "x", min_rows = 4L)
  p <- ncol(x)
  if (p < 2L) {
    stop("The spatial-sign max-sum test requires at least two variables.",
         call. = FALSE)
  }
  if (is.null(mu)) {
    mu <- numeric(p)
  } else {
    mu <- .as_location(mu, p, "mu")
  }
  alpha <- .ssmax_validate_alpha(alpha)
  controls <- .validate_iteration_controls(tol, max_iter, zero_tol)
  strict <- .ssmax_validate_strict(strict)

  max.component <- spatial_sign_max_test(
    x, mu = mu, alpha = alpha, tol = controls$tol,
    max_iter = controls$max_iter, zero_tol = controls$zero_tol,
    strict = strict
  )
  # This direct call is deliberate: the max-sum procedure uses the feasible
  # Feng--Sun leave-two-out sum statistic, not a full-sample approximation.
  sum.component <- feng_sun_one_sample_test(
    x, mu = mu, tol = controls$tol, max_iter = controls$max_iter
  )

  max.log.p <- max.component$components$log.p.MAX
  max.log.q <- max.component$components$log.one.minus.p.MAX
  sum.z <- unname(sum.component$statistic)
  sum.log.p <- stats::pnorm(sum.z, lower.tail = FALSE, log.p = TRUE)
  sum.log.q <- stats::pnorm(sum.z, lower.tail = TRUE, log.p = TRUE)
  combination <- .ssmax_cauchy_combine_logtails(
    max.log.p, max.log.q, sum.log.p, sum.log.q
  )
  variable.names <- .hotelling_variable_names(x)
  variable.mu <- stats::setNames(as.numeric(mu), variable.names)

  .new_hd_location_test(
    statistic = c(Cauchy.angle = combination$angle),
    p.value = combination$p.value,
    method = paste(
      "Liu-Feng-Zhao-Wang spatial-sign max-sum test",
      "(stable Cauchy combination of feasible max and sum components)"
    ),
    data.name = data.name,
    alternative = "two.sided",
    raw.statistic = c(Cauchy = combination$statistic),
    estimate = max.component$estimate,
    null.value = variable.mu,
    null.distribution = list(
      family = "standard Cauchy combination",
      parameters = c(location = 0, scale = 1),
      exact = FALSE,
      tail = "upper",
      assumptions = paste(
        "asymptotic independence of the Liu-Feng-Zhao-Wang max and",
        "Feng-Sun sum components under the paper's conditions"
      )
    ),
    components = list(
      max.test = max.component,
      sum.test = sum.component,
      p.MAX = max.component$p.value,
      log.p.MAX = max.log.p,
      log.one.minus.p.MAX = max.log.q,
      p.SUM = sum.component$p.value,
      log.p.SUM = sum.log.p,
      log.one.minus.p.SUM = sum.log.q,
      Cauchy.statistic = combination$statistic,
      Cauchy.angle = combination$angle,
      log.p.Cauchy = combination$log.p.value,
      log.one.minus.p.Cauchy = combination$log.one.minus.p.value,
      Cauchy.sign = combination$sign,
      Cauchy.log.absolute = combination$log.absolute,
      alpha = alpha,
      rejected = combination$p.value < alpha
    ),
    diagnostics = list(
      max.iteration.stable =
        max.component$diagnostics$iteration.stable,
      max.score.residual = max.component$diagnostics$score.residual,
      sum.leaveout.fits.converged =
        sum.component$diagnostics$all.leaveout.fits.converged,
      combination = "equal-weight signed-log cotangents with atan2 tail",
      probability.clipping = "none",
      reverse.endpoint.policy =
        "exact 0/1 or 1/0 component endpoints are indeterminate errors",
      regularization = "none"
    ),
    n = nrow(x),
    p = p,
    call = call
  )
}
