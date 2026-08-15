#' Feng--Sun scalar-invariant one-sample spatial-sign test
#'
#' Tests \eqn{H_0:\theta=\mu_0} against an unrestricted location alternative
#' with the scalar-invariant high-dimensional spatial-sign procedure of Feng
#' and Sun (2016).  For every unordered pair \eqn{i<j}, the function jointly
#' fits a leave-two-out location \eqn{\widehat\theta_{ij}} and positive diagonal
#' scale \eqn{\widehat D_{ij}} from all observations except \eqn{i,j}.  The raw
#' statistic is
#' \deqn{T_{SS}=\frac{2}{n(n-1)}\sum_{i<j}
#' U\{\widehat D_{ij}^{-1/2}(X_i-\mu_0)\}^\mathsf{T}
#' U\{\widehat D_{ij}^{-1/2}(X_j-\mu_0)\}.}
#' In particular, \eqn{\widehat\theta_{ij}} does **not** enter this numerator.
#'
#' The implementation also supplies the feasible null calibration from the
#' original paper, which is needed for an operational test but is omitted from
#' the short presentation in the book chapter:
#' \deqn{\widehat{\operatorname{tr}(R^2)}=
#' \frac{p^2}{n(n-1)}\sum_{i\ne j}
#' \left[U\{\widehat D_{ij}^{-1/2}(X_i-\widehat\theta_{ij})\}^\mathsf{T}
#' U\{\widehat D_{ij}^{-1/2}(X_j-\widehat\theta_{ij})\}\right]^2,}
#' \deqn{\widehat\sigma_n^2=
#' \frac{2\widehat{\operatorname{tr}(R^2)}}{n(n-1)p^2},
#' \qquad Z_{SS}=T_{SS}/\widehat\sigma_n.}
#' The p-value is the upper tail of the asymptotic standard normal law.  The
#' vector alternative is conventionally labelled `two.sided`, while large
#' positive quadratic evidence is the rejection direction.
#'
#' Each leave-out fit follows the paper's diagonal HR recursion.  Its initial
#' location and scale are the leave-out sample mean and marginal sample
#' variances.  Iteration stops only when both the infinity norm of the mean
#' sign equation and the maximum diagonal-scale equation residual are at most
#' `tol`.  Failure to converge, a coincident training observation (for which
#' the inverse-radius location update is undefined), a non-positive marginal
#' scale, or a non-positive feasible variance raises an error.  No ridge,
#' absolute-value repair, or numerical floor is used.  The paper itself notes
#' that general existence, uniqueness, and convergence of this diagonal HR
#' recursion are not established.
#'
#' Before fitting, null-centred data are divided by a safe scale separately in
#' every variable.  This is an algebraically neutral coordinatewise change of
#' units for the Feng--Sun statistic and protects both very large and very
#' small finite inputs.  Returned `leaveout.scale.diagonal` values are mapped
#' back to input-coordinate geometry and canonically normalised so their
#' largest diagonal entry is one; the scale equations identify the diagonal
#' only up to a common multiplier.  The corresponding canonical log diagonals
#' are also returned, retaining information if a display-scale entry
#' underflows.  `leaveout.location` is in the original input units.
#'
#' The literal algorithm requires \eqn{n(n-1)/2} separate iterative fits and
#' is therefore computationally intensive.  Its compiled implementation is
#' intended for faithful, moderate-sample use rather than silently replacing
#' the published statistic by a full-sample approximation.
#'
#' @param x A numeric matrix or data frame with observations in rows and
#'   variables in columns.  At least four observations are required, and every
#'   leave-two-out fit must have positive variation in every variable.
#' @param mu A finite numeric vector giving the null location.  The default is
#'   the zero vector.
#' @param tol A finite positive tolerance for both leave-out estimating
#'   equations.  The default is `1e-7`.
#' @param max_iter A positive integer giving the maximum number of diagonal HR
#'   updates for each pair.  Non-convergence is an error.
#'
#' @return An object of class `c("hd_location_test", "htest")`.  Components
#'   include `T.SS`, the ordered-pair trace estimate, `sigma2.hat`, pairwise
#'   numerator and variance contributions, and every leave-two-out location
#'   and canonical diagonal-scale fit.  Diagnostics report equation residuals,
#'   iteration counts, exact zero signs, internal scaling, and the no-repair
#'   policy.
#'
#' @references
#' Feng, L. and Sun, F. (2016). Spatial-sign based high-dimensional location
#' test. *Electronic Journal of Statistics*, **10**, 2420--2434.
#' \doi{10.1214/16-EJS1176}.
#'
#' @examples
#' set.seed(2608)
#' x <- matrix(stats::rt(48, df = 5), 8, 6)
#' feng_sun_one_sample_test(x, tol = 1e-6)
#'
#' @export
feng_sun_one_sample_test <- function(x, mu = NULL, tol = 1e-7,
                                     max_iter = 500L) {
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
  controls <- .validate_iteration_controls(tol, max_iter, zero_tol = 0)

  fit <- cpp_feng_sun_one_sample(
    x, mu, controls$tol, controls$max_iter
  )
  z <- unname(fit$z)
  p.value <- stats::pnorm(z, lower.tail = FALSE)
  variable.names <- .hotelling_variable_names(x)
  mu.named <- stats::setNames(as.numeric(mu), variable.names)
  sample.mean <- stats::setNames(
    as.numeric(fit$sample_mean), variable.names
  )

  pair.labels <- paste0(fit$pair_i, ",", fit$pair_j)
  name.pairs <- function(value) {
    stats::setNames(as.numeric(value), pair.labels)
  }
  name.pair.matrix <- function(value) {
    value <- as.matrix(value)
    dimnames(value) <- list(pair.labels, variable.names)
    value
  }
  pair.index <- data.frame(
    i = as.integer(fit$pair_i),
    j = as.integer(fit$pair_j),
    row.names = pair.labels
  )

  pair.iterations <- name.pairs(fit$pair_iterations)
  pair.location.residual <- name.pairs(fit$pair_location_residual)
  pair.scale.residual <- name.pairs(fit$pair_scale_residual)
  pair.test.zero.signs <- name.pairs(fit$pair_test_zero_signs)
  pair.variance.zero.signs <- name.pairs(fit$pair_variance_zero_signs)

  .new_hd_location_test(
    statistic = c(Z = z),
    p.value = p.value,
    method = paste(
      "Feng-Sun one-sample scalar-invariant spatial-sign test",
      "(literal leave-two-out fits; feasible asymptotic normal calibration)"
    ),
    data.name = data.name,
    alternative = "two.sided",
    raw.statistic = c(T.SS = unname(fit$T_SS)),
    estimate = sample.mean,
    null.value = mu.named,
    null.distribution = list(
      family = "normal",
      parameters = c(mean = 0, sd = 1),
      exact = FALSE,
      tail = "upper",
      assumptions = paste(
        "Feng-Sun elliptical high-dimensional conditions (C1)-(C3);",
        "independent observations; converged leave-two-out diagonal HR fits"
      )
    ),
    variance = c(estimated = unname(fit$sigma2_hat)),
    components = list(
      sample.mean = sample.mean,
      null.location = mu.named,
      T.SS = unname(fit$T_SS),
      test.inner.product.sum = unname(fit$test_inner_product_sum),
      variance.inner.product.squared.sum =
        unname(fit$variance_inner_product_squared_sum),
      trace.R2.hat = unname(fit$trace_R2_hat),
      sigma2.hat = unname(fit$sigma2_hat),
      sigma.hat = unname(fit$sigma_hat),
      pair.index = pair.index,
      pair.test.inner.product = name.pairs(
        fit$pair_test_inner_product
      ),
      pair.variance.inner.product.squared = name.pairs(
        fit$pair_variance_inner_product_squared
      ),
      leaveout.location = name.pair.matrix(fit$leaveout_theta),
      leaveout.location.standardized = name.pair.matrix(
        fit$leaveout_theta_standardised
      ),
      leaveout.scale.diagonal = name.pair.matrix(
        fit$leaveout_D_input_canonical
      ),
      leaveout.log.scale.diagonal = name.pair.matrix(
        fit$leaveout_log_D_input_canonical
      ),
      leaveout.scale.diagonal.standardized = name.pair.matrix(
        fit$leaveout_D_standardised
      ),
      pair.count = unname(fit$pair_count),
      ordered.pair.count = unname(fit$ordered_pair_count),
      n = unname(fit$n),
      p = unname(fit$p)
    ),
    diagnostics = list(
      all.leaveout.fits.converged = TRUE,
      tolerance = unname(fit$tolerance),
      max.iterations.allowed = unname(fit$max_iterations),
      pair.iterations = pair.iterations,
      iteration.summary = c(
        minimum = min(pair.iterations),
        median = stats::median(pair.iterations),
        mean = mean(pair.iterations),
        maximum = max(pair.iterations)
      ),
      pair.location.equation.residual = pair.location.residual,
      pair.scale.equation.residual = pair.scale.residual,
      maximum.location.equation.residual =
        max(pair.location.residual),
      maximum.scale.equation.residual = max(pair.scale.residual),
      test.location.center =
        "null location mu; leave-two-out theta is not used in T.SS",
      variance.location.center =
        "pair-specific leave-two-out theta",
      pair.test.zero.signs = pair.test.zero.signs,
      pair.variance.zero.signs = pair.variance.zero.signs,
      test.zero.sign.pair.uses = sum(pair.test.zero.signs),
      variance.zero.sign.pair.uses = sum(pair.variance.zero.signs),
      sign.at.zero = "U(0) = 0 outside the inverse-radius fit update",
      scale.normalization = paste(
        "each internal D has max diagonal 1; returned input-coordinate D",
        "also has max diagonal 1 because the equations identify only ratios"
      ),
      internal.scaling =
        "safe null-centred per-column max-absolute scaling",
      internal.column.log.scale = stats::setNames(
        as.numeric(fit$column_log_residual_scale), variable.names
      ),
      subtraction.overflow.fallback.columns =
        unname(fit$subtraction_overflow_columns),
      computation = "literal pair-specific leave-two-out diagonal HR fits",
      regularization = "none",
      variance.repair = "none"
    ),
    n = n,
    p = p,
    call = call
  )
}
