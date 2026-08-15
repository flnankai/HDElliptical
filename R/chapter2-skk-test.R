#' Srivastava--Katayama--Kano two-sample mean test
#'
#' Tests equality of two high-dimensional mean vectors using the
#' diagonal-standardised statistic of Srivastava, Katayama, and Kano (SKK).
#' The two samples may have different covariance matrices.  With unbiased
#' sample covariance matrices \eqn{S_1} and \eqn{S_2}, define
#' \deqn{\widehat D = \operatorname{diag}(S_1)/n_1 +
#'                    \operatorname{diag}(S_2)/n_2,}
#' \deqn{Q = (\bar X_1-\bar X_2)^\mathsf{T}\widehat D^{-1}
#'             (\bar X_1-\bar X_2), \qquad
#'       \widehat q = (Q-p)/\sqrt p.}
#'
#' This implementation uses the correction in the published corrigendum.  If
#' \eqn{A_k=\widehat D^{-1/2}S_k\widehat D^{-1/2}}, it computes
#' \deqn{\widehat F_k = p^{-1}\left\{
#'       \operatorname{tr}(A_k^2) -
#'       \operatorname{tr}(A_k)^2/(n_k-1)\right\},}
#' and
#' \deqn{\widehat G=p^{-1}\operatorname{tr}(A_1A_2).}
#' Thus
#' \deqn{\widehat V_q = 2\widehat F_1/n_1^2 +
#'       2\widehat F_2/n_2^2 + 4\widehat G/(n_1n_2).}
#' The feasible finite-sample factor is calculated from the *sample* matrix
#' \deqn{\widehat R=A_1/n_1+A_2/n_2, \qquad
#'       \widehat c=1+\operatorname{tr}(\widehat R^2)/p^{3/2},}
#' and the reported statistic is
#' \deqn{Z_{\mathrm{SKK}}=\widehat q/
#'       \sqrt{\widehat V_q\widehat c}.}
#'
#' Trace functionals are formed through a \eqn{p\times p} primal Gram matrix
#' only when \eqn{p\le n_1+n_2}; otherwise an
#' \eqn{(n_1+n_2)\times(n_1+n_2)} dual Gram matrix is used, so the
#' high-dimensional path does not construct a \eqn{p\times p} object.  Before
#' centring, both samples are divided by the same scale within each variable.
#' This algebraically neutral change of units protects very large and small
#' inputs while preserving componentwise diagonal-scale invariance.  Every
#' entry of \eqn{\widehat D} and the final normalising variance must be strictly
#' positive.  No ridge, absolute-value repair, or variance flooring is used.
#'
#' @param x,y Numeric matrices or data frames with observations in rows and
#'   the same variables in columns.  Each group must have at least two
#'   observations.
#'
#' @return An object of class `c("hd_location_test", "htest")`.  Besides the
#'   upper-tail asymptotic normal p-value, `components` contains `Q.raw`,
#'   `q.scaled`, `F1`, `F2`, `G`, the variance estimate, the sample
#'   `trace.R.hat2` and `c.hat`, sample means, and marginal variance
#'   estimates.  `diagnostics` records the primal/dual path, internal column
#'   scales, covariance denominators, and the absence of regularisation.
#'
#' @references
#' Srivastava, M. S., Katayama, S., and Kano, Y. (2013). A two sample test in
#' high dimensional data. *Journal of Multivariate Analysis*, **114**,
#' 349--358. \doi{10.1016/j.jmva.2012.08.014}.
#'
#' Srivastava, M. S., Katayama, S., and Kano, Y. (2013). Corrigendum to "A two
#' sample test in high dimensional data". *Journal of Multivariate Analysis*,
#' **120**, 251. \doi{10.1016/j.jmva.2013.04.016}.
#'
#' @examples
#' set.seed(2301)
#' x <- matrix(rnorm(60), 12, 5)
#' y <- matrix(rnorm(70, 0.2), 14, 5)
#' srivastava_katayama_kano_two_sample_test(x, y)
#'
#' @export
srivastava_katayama_kano_two_sample_test <- function(x, y) {
  call <- match.call()
  x.name <- deparse1(substitute(x))
  y.name <- deparse1(substitute(y))
  x <- .as_data_matrix(x, "x", min_rows = 2L)
  y <- .as_data_matrix(y, "y", min_rows = 2L)
  .check_two_sample_variables(x, y)

  n1 <- nrow(x)
  n2 <- nrow(y)
  p <- ncol(x)
  fit <- cpp_skk_two_sample(x, y)
  z <- unname(fit$z)
  p.value <- stats::pnorm(z, lower.tail = FALSE)

  variable.names <- .hotelling_variable_names(x)
  mean.x <- stats::setNames(as.numeric(fit$mean_x), variable.names)
  mean.y <- stats::setNames(as.numeric(fit$mean_y), variable.names)
  difference <- stats::setNames(
    as.numeric(fit$difference), variable.names
  )
  variance1.diagonal <- stats::setNames(
    as.numeric(fit$variance1_diagonal), variable.names
  )
  variance2.diagonal <- stats::setNames(
    as.numeric(fit$variance2_diagonal), variable.names
  )
  d.hat.diagonal <- stats::setNames(
    as.numeric(fit$D_hat_diagonal), variable.names
  )

  .new_hd_location_test(
    statistic = c(Z = z),
    p.value = p.value,
    method = paste(
      "Srivastava-Katayama-Kano two-sample high-dimensional mean test",
      "(unequal covariance; asymptotic normal calibration)"
    ),
    data.name = paste(x.name, "and", y.name),
    alternative = "two.sided",
    raw.statistic = c(
      Q = unname(fit$Q_raw),
      q.SKK = unname(fit$q_scaled)
    ),
    estimate = difference,
    null.value = stats::setNames(numeric(p), variable.names),
    null.distribution = list(
      family = "normal",
      parameters = c(mean = 0, sd = 1),
      exact = FALSE,
      tail = "upper",
      assumptions = paste(
        "Srivastava-Katayama-Kano high-dimensional Gaussian asymptotics;",
        "independent samples; strictly positive combined marginal variances"
      )
    ),
    variance = c(denominator = unname(fit$denominator_variance)),
    components = list(
      mean.x = mean.x,
      mean.y = mean.y,
      difference = difference,
      variance1.diagonal = variance1.diagonal,
      variance2.diagonal = variance2.diagonal,
      D.hat.diagonal = d.hat.diagonal,
      n1 = unname(fit$n1),
      n2 = unname(fit$n2),
      p = unname(fit$p),
      Q.raw = unname(fit$Q_raw),
      q.scaled = unname(fit$q_scaled),
      F1 = unname(fit$F1),
      F2 = unname(fit$F2),
      G = unname(fit$G),
      q.variance = unname(fit$variance_q),
      trace.R.hat2 = unname(fit$trace_R_hat2),
      c.hat = unname(fit$c_hat),
      denominator.variance = unname(fit$denominator_variance),
      trace.A1 = unname(fit$trace_A1),
      trace.A2 = unname(fit$trace_A2),
      trace.A1.squared = unname(fit$trace_A1_squared),
      trace.A2.squared = unname(fit$trace_A2_squared),
      trace.A1.A2 = unname(fit$trace_A1_A2)
    ),
    diagnostics = list(
      covariance.model = "unequal",
      sample.covariance.denominators = c(
        group1 = n1 - 1,
        group2 = n2 - 1
      ),
      F.bias.correction.denominators = c(
        group1 = n1 - 1,
        group2 = n2 - 1
      ),
      q.standardisation = "(Q - p) / sqrt(p)",
      finite.sample.correction =
        "1 + trace(sample R-hat squared) / p^(3/2)",
      trace.computation = paste(unname(fit$gram_type), "Gram matrix"),
      gram.type = unname(fit$gram_type),
      gram.dimension = unname(fit$gram_dimension),
      constructs.p.by.p.matrix = identical(unname(fit$gram_type), "primal"),
      internal.scaling = "per-column common scaling of x and y",
      internal.scale.factors = stats::setNames(
        as.numeric(fit$column_scale), variable.names
      ),
      regularization = "none",
      variance.repair = "none"
    ),
    n = c(group1 = n1, group2 = n2),
    p = p,
    call = call
  )
}
