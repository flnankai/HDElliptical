#' Park--Ayyala one-sample high-dimensional mean test
#'
#' Tests \eqn{H_0: \mu=\mu_0} with the diagonal leave-two-out statistic of
#' Park and Ayyala (2013). For every ordered pair of observations, the
#' diagonal covariance used in the quadratic product is estimated after
#' removing that pair. The estimated null variance uses the matching
#' leave-two-out means. At least six observations are required: below this
#' threshold the inverse-chi-square bias correction is not valid (and at
#' \eqn{n=5} the corrected numerator is identically zero). The algebraic API
#' therefore permits \eqn{n=6,7}, but squared inverse leave-out variances have
#' no finite Gaussian expectation at those sizes; their asymptotic normal
#' calibration is especially fragile, so \eqn{n\geq 8} is strongly preferred.
#'
#' Each residual `x - mu` is formed in extended precision and divided by a
#' common scale within its variable. This is algebraically neutral for the
#' Park--Ayyala statistic and protects its leave-out moments from avoidable
#' overflow, underflow, and cancellation under a large common translation. No
#' zero variance is replaced, no ridge is added, and a non-positive estimated
#' null variance is reported as an error.
#'
#' @param x A numeric matrix or data frame with observations in rows and at
#'   least six rows.
#' @param mu A finite null-mean vector. `NULL` uses the zero vector.
#'
#' @return An object of class `c("hd_location_test", "htest")`. The
#'   `components` field contains the two ordered-pair sums, bias-corrected raw
#'   statistic, and estimated null variance. `finite.sample.correction` is
#'   \eqn{(n-5)/(n-3)}, whereas `bias.factor` is the full coefficient
#'   \eqn{(n-5)/\{n(n-1)(n-3)\}} multiplying the ordered-pair sum.
#'
#' @references
#' Park, J. and Ayyala, D. N. (2013). A test for the mean vector in large
#' dimension and small samples. *Journal of Statistical Planning and
#' Inference*, 143, 929--943.
#'
#' @examples
#' set.seed(23)
#' x <- matrix(rnorm(24), nrow = 8, ncol = 3)
#' park_ayyala_one_sample_test(x)
#'
#' @export
park_ayyala_one_sample_test <- function(x, mu = NULL) {
  call <- match.call()
  data.name <- deparse1(substitute(x))
  x <- .as_data_matrix(x, "x", min_rows = 6L)
  n <- nrow(x)
  p <- ncol(x)
  if (is.null(mu)) {
    mu <- numeric(p)
  } else {
    mu <- .as_location(mu, p, "mu")
  }

  out <- cpp_park_ayyala_one_sample(x, mu)
  z <- as.numeric(out$z)
  p.value <- stats::pnorm(z, lower.tail = FALSE)
  variable.names <- .hotelling_variable_names(x)
  sample.mean <- stats::setNames(as.numeric(out$sample_mean), variable.names)
  difference <- stats::setNames(as.numeric(out$difference), variable.names)
  mu <- stats::setNames(as.numeric(mu), variable.names)

  .new_hd_location_test(
    statistic = c(Z = z),
    p.value = p.value,
    method = paste(
      "Park-Ayyala one-sample high-dimensional mean test",
      "(leave-two-out diagonal covariance; asymptotic normal calibration)"
    ),
    data.name = data.name,
    alternative = "two.sided",
    raw.statistic = c(U.PA = as.numeric(out$U)),
    estimate = sample.mean,
    null.value = mu,
    null.distribution = list(
      family = "normal",
      parameters = c(mean = 0, sd = 1),
      exact = FALSE,
      tail = "upper",
      assumptions = paste(
        "Park-Ayyala high-dimensional asymptotics;",
        "strictly positive leave-two-out marginal variances"
      )
    ),
    variance = c(denominator = as.numeric(out$variance)),
    components = list(
      sample.mean = sample.mean,
      difference = difference,
      N = as.numeric(out$N),
      p = as.numeric(out$p),
      ordered.pair.sum = as.numeric(out$T1),
      variance.pair.sum = as.numeric(out$T2),
      finite.sample.correction = as.numeric(out$finite_sample_correction),
      bias.factor = as.numeric(out$bias_factor),
      U.PA = as.numeric(out$U),
      denominator.variance = as.numeric(out$variance)
    ),
    diagnostics = list(
      leaveout.pairs = as.numeric(out$leaveout_pairs),
      minimum.scaled.leaveout.variance =
        as.numeric(out$min_leaveout_variance),
      internal.scaling =
        "extended-precision x - mu, then per-column residual scaling",
      internal.scale.factors = stats::setNames(
        as.numeric(out$column_scale), variable.names
      ),
      regularization = "none",
      variance.repair = "none"
    ),
    n = c(x = n),
    p = p,
    call = call
  )
}


#' Chen--Qin two-sample high-dimensional mean test
#'
#' Tests equality of two mean vectors with the covariance-unrestricted
#' U-statistic of Chen and Qin (2010). Its numerator is
#' \deqn{T_{CQ}=\|\bar X-\bar Y\|^2-
#'   \operatorname{tr}(S_X)/n_1-\operatorname{tr}(S_Y)/n_2.}
#' The variance estimator is the original paper's leave-out estimator of the
#' two within-sample covariance traces and their cross trace. This is distinct
#' from the common-covariance Bai--Saranadasa test.
#'
#' The numerator is exactly invariant to a common translation. The original
#' finite-sample leave-out variance estimator, however, is not translation
#' invariant; this known property is retained rather than silently replacing
#' the published estimator. Both samples are internally divided by one common
#' global scale. The returned `components` include both the original-unit
#' quantities and their finite internally scaled counterparts. An
#' original-unit fourth-order quantity can legitimately overflow to `Inf` or
#' underflow to zero when the input itself spans the limits of double
#' precision; the reported statistic is reconstructed from the scaled
#' quantities and remains well defined. No ridge, absolute-value repair, or
#' variance flooring is applied.
#'
#' @param x,y Numeric matrices or data frames with observations in rows and
#'   the same variables in columns. Each group must contain at least three
#'   observations.
#'
#' @return An object of class `c("hd_location_test", "htest")`. The
#'   `components` field includes the U-statistic numerator, leave-out trace
#'   estimates `A1`, `A2`, and `A12`, the estimated null variance, and their
#'   `.scaled` counterparts used for stable calibration.
#'
#' @references
#' Chen, S. X. and Qin, Y.-L. (2010). A two-sample test for high-dimensional
#' data with applications to gene-set testing. *Annals of Statistics*, 38,
#' 808--835.
#'
#' @examples
#' set.seed(24)
#' x <- matrix(rnorm(21), nrow = 7, ncol = 3)
#' y <- matrix(rnorm(24, mean = 0.2), nrow = 8, ncol = 3)
#' chen_qin_two_sample_test(x, y)
#'
#' @export
chen_qin_two_sample_test <- function(x, y) {
  call <- match.call()
  x.name <- deparse1(substitute(x))
  y.name <- deparse1(substitute(y))
  x <- .as_data_matrix(x, "x", min_rows = 3L)
  y <- .as_data_matrix(y, "y", min_rows = 3L)
  .check_two_sample_variables(x, y)

  n1 <- nrow(x)
  n2 <- nrow(y)
  p <- ncol(x)
  out <- cpp_chen_qin_two_sample(x, y)
  z <- as.numeric(out$z)
  p.value <- stats::pnorm(z, lower.tail = FALSE)
  variable.names <- .hotelling_variable_names(x)
  mean.x <- stats::setNames(as.numeric(out$mean_x), variable.names)
  mean.y <- stats::setNames(as.numeric(out$mean_y), variable.names)
  difference <- stats::setNames(as.numeric(out$difference), variable.names)

  .new_hd_location_test(
    statistic = c(Z = z),
    p.value = p.value,
    method = paste(
      "Chen-Qin two-sample high-dimensional mean test",
      "(unequal covariance; asymptotic normal calibration)"
    ),
    data.name = paste(x.name, "and", y.name),
    alternative = "two.sided",
    raw.statistic = c(T.CQ = as.numeric(out$T)),
    estimate = difference,
    null.value = stats::setNames(numeric(p), variable.names),
    null.distribution = list(
      family = "normal",
      parameters = c(mean = 0, sd = 1),
      exact = FALSE,
      tail = "upper",
      assumptions = paste(
        "Chen-Qin high-dimensional asymptotics; independent samples;",
        "positive original leave-out variance estimate"
      )
    ),
    variance = c(denominator = as.numeric(out$variance)),
    components = list(
      mean.x = mean.x,
      mean.y = mean.y,
      difference = difference,
      n1 = as.numeric(out$n1),
      n2 = as.numeric(out$n2),
      p = as.numeric(out$p),
      trace.S1 = as.numeric(out$trace_S1),
      trace.S2 = as.numeric(out$trace_S2),
      T.CQ = as.numeric(out$T),
      A1 = as.numeric(out$A1),
      A2 = as.numeric(out$A2),
      A12 = as.numeric(out$A12),
      denominator.variance = as.numeric(out$variance),
      T.CQ.scaled = as.numeric(out$T_scaled),
      A1.scaled = as.numeric(out$A1_scaled),
      A2.scaled = as.numeric(out$A2_scaled),
      A12.scaled = as.numeric(out$A12_scaled),
      denominator.variance.scaled = as.numeric(out$variance_scaled)
    ),
    diagnostics = list(
      covariance.model = "unrestricted",
      variance.estimator = "original Chen-Qin leave-out estimator",
      numerator.translation.invariant = TRUE,
      finite.sample.variance.translation.invariant = FALSE,
      internal.scaling = "common global scaling of both groups",
      internal.scale.factor = as.numeric(out$global_scale),
      statistic.reconstruction =
        "T.CQ.scaled / sqrt(denominator.variance.scaled)",
      original.unit.components = paste(
        "may overflow or underflow outside the representable double range;",
        "scaled components remain the computational source of truth"
      ),
      regularization = "none",
      variance.repair = "none"
    ),
    n = c(x = n1, y = n2),
    p = p,
    call = call
  )
}
