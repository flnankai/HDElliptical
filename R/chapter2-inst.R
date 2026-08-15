#' Feng--Liu--Ma one-sample inverse norm sign test
#'
#' Tests \eqn{H_0:\theta=\mu_0} against an unrestricted multivariate
#' location alternative with the inverse norm sign test (INST) of Feng, Liu,
#' and Ma (2021).  For every unordered pair \eqn{i<j}, a location and positive
#' diagonal scale are jointly fitted from the sample with observations
#' \eqn{i,j} removed.  Writing the fitted scale as \eqn{\widehat D_{ij}},
#' \eqn{r_{ij,k}=\lVert\widehat D_{ij}^{-1/2}(X_k-\mu_0)\rVert}, and
#' \eqn{U_{ij,k}=U\{\widehat D_{ij}^{-1/2}(X_k-\mu_0)\}}, the primary raw
#' statistic is
#' \deqn{T_{INST}=\frac{2}{n(n-1)}\sum_{i<j}
#' r_{ij,i}^{-1}r_{ij,j}^{-1}U_{ij,i}^{\mathsf T}U_{ij,j}.}
#' The pair-specific fitted location is used only to estimate
#' \eqn{\widehat D_{ij}}.  It is deliberately not subtracted from either
#' endpoint in this statistic.
#'
#' The operational calibration is the direct feasible variance in Section
#' S.3 of the paper's official supplement.  If
#' \eqn{\widetilde\mu_{ij}=(n-2)^{-1}\sum_{k\ne i,j}U_{ij,k}}, then
#' \deqn{\widehat\sigma_n^2=2n^{-4}\sum_{i\ne j}
#' r_{ij,i}^{-2}r_{ij,j}^{-2}
#' \{(U_{ij,i}-\widetilde\mu_{ij})^{\mathsf T}U_{ij,j}\}
#' \{(U_{ij,j}-\widetilde\mu_{ij})^{\mathsf T}U_{ij,i}\}.}
#' This is an ordered-pair sum with the published \eqn{2n^{-4}} multiplier;
#' it is not replaced by a finite-sample pair-count denominator.  The reported
#' p-value is based only on \eqn{T_{INST}/\widehat\sigma_n}.
#'
#' For interpretation, the result additionally reports cross-fitted estimates
#' of \eqn{E(r^{-1})}, \eqn{\nu_{2,IN}=E(r^{-2})},
#' \eqn{c_{0,IN}=E(r^{-2})}, and \eqn{\operatorname{tr}(R^2)}.  The displayed
#' oracle-factorized variance
#' \eqn{2\widehat\nu_{2,IN}^2\widehat{\operatorname{tr}(R^2)}/
#' \{n(n-1)p^2\}} is diagnostic only and never replaces the primary direct
#' variance.
#'
#' Each diagonal HR recursion starts at the leave-two-out sample mean and
#' marginal sample variances.  The common scale from that initialization is
#' retained because the fixed-point equations identify only relative diagonal
#' scales.  `iteration.stable` means that both the relative location update and
#' log-diagonal update are at most `tol`; the estimating-equation
#' `score.residual` is returned separately.  With `strict = TRUE`, any required
#' fit that is not stable after `max_iter` updates is an error.  With
#' `strict = FALSE`, the last iterates are used with a warning and complete
#' diagnostics.
#'
#' Exact or tolerance-defined zero radii make inverse norm weighting singular.
#' Such a case is an error.  A non-finite or non-positive published variance is
#' also an error.  No radius perturbation, weight cap, ridge, absolute value,
#' numerical floor, or variance substitution is applied.
#'
#' Null-centred residuals are internally divided by a positive scale in each
#' coordinate.  This is algebraically neutral under a common translation and
#' nonsingular diagonal change of units, and protects the calculation from
#' avoidable overflow.  Returned input-coordinate diagonals are canonical
#' displays with largest entry one; their canonical log values retain scale
#' ratio information when a displayed entry underflows.  The unnormalised
#' internal diagonals retain the initialization's common scale.
#'
#' @param x A numeric matrix or data frame with observations in rows and
#'   variables in columns.  At least four observations are required.
#' @param mu A finite numeric vector giving the null location.  The default is
#'   the zero vector.
#' @param alpha Significance level for the returned upper-tail rejection rule.
#' @param tol Positive tolerance for relative location and log-diagonal
#'   updates in every leave-two-out fit.
#' @param max_iter Positive maximum number of updates for each fit.
#' @param zero_tol Non-negative tolerance below which a standardized radius is
#'   treated as singular.  The default detects exact zeros only.
#' @param strict If `TRUE`, error when any leave-two-out iteration is not
#'   stable.  If `FALSE`, warn and return the last iterates.
#'
#' @return An object of class `c("hd_location_test", "htest")`.  Its
#'   `components` field contains the primary statistic and direct variance,
#'   cross-fitted radial moments and trace diagnostics, every pair's kernels,
#'   radii, fitted locations, and diagonal scales.  Its `diagnostics` field
#'   records iteration stability, score residuals, centring roles, rejection
#'   rule, numerical scaling, and the no-repair contract.
#'
#' @references
#' Feng, L., Liu, B., and Ma, Y. (2021). An inverse norm sign test of location
#' parameter for high-dimensional data. *Journal of Business & Economic
#' Statistics*, **39**, 807--815.
#' \doi{10.1080/07350015.2020.1736084}.
#'
#' Feng, L., Liu, B., and Ma, Y. (2020). Supplementary material for *An
#' inverse norm sign test of location parameter for high-dimensional data*.
#' \doi{10.6084/m9.figshare.11914095.v2}.
#'
#' @examples
#' set.seed(2021)
#' x <- matrix(stats::rt(40, df = 5), 10, 4)
#' inst_one_sample_test(x, tol = 1e-6)
#'
#' @export
inst_one_sample_test <- function(x, mu = NULL, alpha = 0.05,
                                 tol = 1e-8, max_iter = 500L,
                                 zero_tol = 0, strict = TRUE) {
  call <- match.call()
  data.name <- deparse1(substitute(x))
  x <- .as_data_matrix(x, "x", min_rows = 4L)
  n <- nrow(x)
  p <- ncol(x)
  if (is.null(mu)) {
    mu <- numeric(p)
  } else {
    mu <- .as_location(mu, p, "mu")
  }
  alpha <- as.numeric(alpha)
  if (length(alpha) != 1L || is.na(alpha) || !is.finite(alpha) ||
      alpha <= 0 || alpha >= 1) {
    stop("`alpha` must be one finite number strictly between zero and one.",
         call. = FALSE)
  }
  controls <- .validate_iteration_controls(tol, max_iter, zero_tol)
  if (!is.logical(strict) || length(strict) != 1L || is.na(strict)) {
    stop("`strict` must be `TRUE` or `FALSE`.", call. = FALSE)
  }

  fit <- cpp_inst_one_sample(
    x, mu, controls$tol, controls$max_iter, controls$zero_tol
  )
  if (!isTRUE(fit$all_iteration_stable)) {
    message <- sprintf(
      paste0(
        "INST did not become iteration-stable for %d leave-two-out fit(s) ",
        "in `max_iter = %d` updates. No ridge or iterate substitution was ",
        "applied."
      ),
      as.integer(fit$stability_failures), controls$max_iter
    )
    if (isTRUE(strict)) {
      stop(message, call. = FALSE)
    }
    warning(message, call. = FALSE)
  }

  variable.names <- .hotelling_variable_names(x)
  mu.named <- stats::setNames(as.numeric(mu), variable.names)
  sample.mean <- stats::setNames(
    as.numeric(fit$sample_mean), variable.names
  )
  z <- as.numeric(fit$z)
  critical.value <- stats::qnorm(1 - alpha)
  reject <- isTRUE(z > critical.value)

  pair.labels <- paste0(fit$pair_i, ",", fit$pair_j)
  pair.index <- data.frame(
    i = as.integer(fit$pair_i),
    j = as.integer(fit$pair_j),
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
  variance.factor.columns <- c(
    "(U_i-mu_tilde)'U_j", "(U_j-mu_tilde)'U_i"
  )
  pair.fit <- fit$fit_diagnostics
  names(pair.fit$iterations) <- pair.labels
  names(pair.fit$iteration.stable) <- pair.labels
  names(pair.fit$relative.update) <- pair.labels
  names(pair.fit$location.relative.update) <- pair.labels
  names(pair.fit$log.diagonal.relative.update) <- pair.labels
  names(pair.fit$score.residual) <- pair.labels
  names(pair.fit$location.score.residual) <- pair.labels
  names(pair.fit$diagonal.score.residual) <- pair.labels
  names(pair.fit$minimum.residual.distance) <- pair.labels

  .new_hd_location_test(
    statistic = c(Z = z),
    p.value = stats::pnorm(z, lower.tail = FALSE),
    method = paste(
      "Feng-Liu-Ma one-sample inverse norm sign test",
      "(INST; primary S.3 direct feasible variance)"
    ),
    data.name = data.name,
    alternative = "two.sided",
    raw.statistic = c(T.INST = as.numeric(fit$T)),
    estimate = sample.mean,
    null.value = mu.named,
    null.distribution = list(
      family = "normal",
      parameters = c(mean = 0, sd = 1),
      exact = FALSE,
      tail = "upper",
      assumptions = paste(
        "Feng-Liu-Ma high-dimensional elliptical asymptotics;",
        "independent observations; finite inverse radial moments;",
        "stable leave-two-out diagonal HR fits"
      )
    ),
    variance = c(primary.direct = as.numeric(fit$variance)),
    components = list(
      sample.mean = sample.mean,
      null.location = mu.named,
      T.INST = as.numeric(fit$T),
      sigma2.direct.hat = as.numeric(fit$variance),
      sigma.direct.hat = as.numeric(fit$standard_error),
      ordered.variance.kernel.sum =
        as.numeric(fit$ordered_variance_kernel_sum),
      inverse.radius.mean.hat = as.numeric(fit$inverse_radius_mean),
      nu2.IN.hat = as.numeric(fit$nu2),
      c0.IN.hat = as.numeric(fit$c0),
      trace.R2.hat = as.numeric(fit$trace_R2_hat),
      sigma2.oracle.factorized.diagnostic =
        as.numeric(fit$oracle_factorized_variance),
      trace.weighted.score.covariance.squared.hat =
        as.numeric(fit$weighted_score_trace_hat),
      pair.index = pair.index,
      pair.test.inner.product = name.pairs(
        fit$pair_test_inner_product
      ),
      pair.statistic.kernel = name.pairs(fit$pair_statistic_kernel),
      pair.leaveout.sign.mean.norm = name.pairs(
        fit$pair_sign_mean_norm
      ),
      pair.variance.factor = name.pair.matrix(
        fit$pair_variance_factor, variance.factor.columns
      ),
      pair.variance.kernel = name.pairs(fit$pair_variance_kernel),
      pair.trace.inner.product.squared = name.pairs(
        fit$pair_trace_inner_product_squared
      ),
      pair.trace.zero.signs = name.pairs(fit$pair_trace_zero_signs),
      endpoint.radius = name.pair.matrix(
        fit$endpoint_radius, endpoint.columns
      ),
      endpoint.inverse.radius = name.pair.matrix(
        fit$endpoint_inverse_radius, endpoint.columns
      ),
      leaveout.location = name.pair.matrix(fit$leaveout_location),
      leaveout.location.standardized = name.pair.matrix(
        fit$leaveout_location_standardized
      ),
      leaveout.scale.diagonal.standardized = name.pair.matrix(
        fit$leaveout_diagonal_standardized
      ),
      leaveout.scale.diagonal = name.pair.matrix(
        fit$leaveout_diagonal_input_canonical
      ),
      leaveout.log.scale.diagonal = name.pair.matrix(
        fit$leaveout_log_diagonal_input_canonical
      ),
      pair.count = as.numeric(fit$pair_count),
      ordered.pair.count = as.numeric(fit$ordered_pair_count),
      n = as.numeric(fit$n),
      p = as.numeric(fit$p)
    ),
    diagnostics = list(
      iteration.stable = isTRUE(fit$all_iteration_stable),
      stability.failures = as.integer(fit$stability_failures),
      relative.update = as.numeric(pair.fit$worst.relative.update),
      score.residual = as.numeric(pair.fit$worst.score.residual),
      minimum.residual.distance =
        as.numeric(pair.fit$smallest.residual.distance),
      convergence.basis = pair.fit$convergence.basis,
      score.residual.role = paste(
        "estimating-equation diagnostic reported separately from relative",
        "iterate stability"
      ),
      pair.fit = pair.fit,
      endpoint.center = paste(
        "null location mu only; the leave-two-out fitted location",
        "participates only in estimating D_ij"
      ),
      variance.sign.center = paste(
        "unweighted mean of null-centred U_ij,k over k != i,j;",
        "not the fitted leave-two-out location"
      ),
      primary.variance.estimator = paste(
        "Feng-Liu-Ma official supplement Section S.3:",
        "2 * n^(-4) times an ordered i != j sum"
      ),
      primary.calibration = "direct feasible variance only",
      oracle.factorized.variance.role = "diagnostic only",
      radial.moment.convention = paste(
        "cross-fitted endpoint average over n(n-1) pair uses;",
        "nu2.IN = c0.IN = mean(r^(-2)) for omega(r)=1/r"
      ),
      trace.estimator = paste(
        "pair-specific leave-two-out fitted-location directions;",
        "p^2/[n(n-1)] times the ordered squared-inner-product sum"
      ),
      pair.convention = "unordered fits; primary variance reconstructed as ordered pairs",
      rejection = list(
        alpha = alpha,
        critical.value = critical.value,
        rule = "reject when Z > qnorm(1 - alpha)",
        reject = reject
      ),
      tolerance = controls$tol,
      max.iterations.allowed = controls$max_iter,
      zero.tolerance = controls$zero_tol,
      inverse.norm.at.zero = "undefined; error without perturbation",
      internal.scaling = paste(
        "null-centred positive per-coordinate scaling; common translation",
        "and nonsingular diagonal changes of units are neutral"
      ),
      internal.column.log.scale = stats::setNames(
        as.numeric(fit$column_log_residual_scale), variable.names
      ),
      subtraction.overflow.fallback.columns =
        as.numeric(fit$subtraction_overflow_columns),
      scale.identification = paste(
        "the internal diagonal retains the common scale of the marginal",
        "sample-variance initialization; displayed input diagonals are",
        "canonicalized only for reporting"
      ),
      regularization = "none",
      radius.perturbation = "none",
      weight.cap = "none",
      variance.repair = "none"
    ),
    n = n,
    p = p,
    call = call
  )
}
