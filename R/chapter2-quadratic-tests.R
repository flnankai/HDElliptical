#' Srivastava--Du one-sample high-dimensional mean test
#'
#' Tests \eqn{H_0: \mu = \mu_0} using the diagonal-standardised
#' Srivastava--Du statistic.  The alternative is unrestricted, while the
#' non-negative quadratic evidence is calibrated in the upper tail of its
#' asymptotic standard normal distribution.
#'
#' The implementation follows the finite-sample centring and scale correction
#' in Srivastava and Du (2008).  It computes \eqn{\operatorname{tr}(R^2)}
#' through the smaller of a primal \eqn{p \times p} matrix and a dual
#' \eqn{N \times N} Gram matrix.  Before centring, each column of `x` and the
#' matching entry of `mu` are divided by one common column scale.  This is an
#' algebraically neutral change of units that protects the statistic at very
#' large or small numerical scales; reported means, differences, and marginal
#' variances are converted back to the input units.  Every sample variance must
#' be strictly positive; no ridge, absolute-value repair, or variance flooring
#' is applied.
#'
#' @param x A numeric matrix or data frame with observations in rows and
#'   variables in columns.  At least four observations are required.
#' @param mu A numeric vector giving the null mean.  The default is the zero
#'   vector.
#'
#' @return An object of class `c("hd_location_test", "htest")`.  In addition
#'   to the standardised statistic and upper-tail p-value, `components`
#'   contains `A`, `nu`, `trace.R2`, `c`, the null centre, the denominator
#'   variance, the sample mean, the mean difference, and the marginal sample
#'   variances.
#'
#' @references
#' Srivastava, M. S. and Du, M. (2008). A test for the mean vector with fewer
#' observations than the dimension. *Journal of Multivariate Analysis*, 99,
#' 386--402.
#'
#' @examples
#' set.seed(21)
#' x <- matrix(rnorm(32), nrow = 8, ncol = 4)
#' srivastava_du_one_sample_test(x)
#'
#' @export
srivastava_du_one_sample_test <- function(x, mu = NULL) {
  call <- match.call()
  data.name <- deparse(substitute(x))
  x <- .as_data_matrix(x, "x", min_rows = 4L)

  n <- nrow(x)
  p <- ncol(x)
  if (is.null(mu)) {
    mu <- numeric(p)
  } else {
    mu <- .as_location(mu, p, "mu")
  }

  out <- cpp_srivastava_du_one_sample(x, mu)
  z <- unname(out$z)
  p.value <- stats::pnorm(z, lower.tail = FALSE)

  variable.names <- .hotelling_variable_names(x)
  sample.mean <- stats::setNames(as.numeric(out$sample_mean), variable.names)
  difference <- stats::setNames(as.numeric(out$difference), variable.names)
  mu <- stats::setNames(as.numeric(mu), variable.names)
  variance.diagonal <- stats::setNames(
    as.numeric(out$variance_diagonal),
    variable.names
  )

  .new_hd_location_test(
    statistic = c(Z = z),
    p.value = p.value,
    method = paste(
      "Srivastava-Du one-sample high-dimensional mean test",
      "(asymptotic normal calibration)"
    ),
    data.name = data.name,
    alternative = "two.sided",
    raw.statistic = c(A = unname(out$A)),
    estimate = sample.mean,
    null.value = mu,
    null.distribution = list(
      family = "normal",
      parameters = c(mean = 0, sd = 1),
      exact = FALSE,
      tail = "upper",
      assumptions = paste(
        "Srivastava-Du high-dimensional asymptotics;",
        "strictly positive marginal sample variances"
      )
    ),
    variance = c(denominator = unname(out$variance)),
    components = list(
      sample.mean = sample.mean,
      difference = difference,
      variance.diagonal = variance.diagonal,
      A = unname(out$A),
      N = unname(out$N),
      nu = unname(out$nu),
      p = unname(out$p),
      trace.R2 = unname(out$trace_R2),
      c = unname(out$c),
      null.center = unname(out$null_center),
      denominator.variance = unname(out$variance)
    ),
    diagnostics = list(
      trace.computation = paste(unname(out$gram_type), "Gram matrix"),
      gram.type = unname(out$gram_type),
      gram.dimension = unname(out$gram_dimension),
      constructs.p.by.p.matrix = identical(unname(out$gram_type), "primal"),
      internal.scaling = "per-column common scaling of x and mu",
      internal.scale.factors = stats::setNames(
        as.numeric(out$column_scale),
        variable.names
      ),
      regularization = "none",
      variance.repair = "none"
    ),
    n = n,
    p = p,
    call = call
  )
}

#' Original Bai--Saranadasa two-sample mean test
#'
#' Tests equality of two multivariate means under a common covariance matrix
#' using the original Bai--Saranadasa statistic.  This is deliberately a
#' separate API from covariance-unrestricted Chen--Qin tests: its numerator is
#' \deqn{\xi = \|\bar X_1 - \bar X_2\|^2
#'              - \frac{N}{n_1 n_2}\operatorname{tr}(S_p),}
#' where \eqn{S_p} is the unbiased pooled covariance matrix and
#' \eqn{N=n_1+n_2}.
#'
#' The traces of \eqn{S_p} and \eqn{S_p^2} are computed using the smaller of a
#' primal \eqn{p \times p} matrix and a dual \eqn{N \times N} Gram matrix.
#' Both groups are divided by the same global numerical scale before any
#' moments are formed.  The standardised statistic is computed on that safe
#' scale, while reported means, traces, numerator, and variance are converted
#' back according to their physical dimensions.  Consequently a reported
#' fourth-order component can be infinite when its value exceeds double range,
#' or zero when it is below the representable range, without making the
#' standardised statistic or p-value non-finite.  The finite internal
#' normalising variance is retained in `diagnostics`.  The estimated null
#' variance must be strictly positive; no ridge, absolute-value repair, or
#' variance flooring is applied.
#'
#' @param x,y Numeric matrices or data frames, with observations in rows and
#'   the same variables in columns.  Each group must have at least two
#'   observations.
#'
#' @return An object of class `c("hd_location_test", "htest")`.  Its
#'   `components` include the original numerator `xi`, pooled covariance
#'   traces, the bias-corrected estimate of \eqn{\operatorname{tr}(\Sigma^2)},
#'   and the estimated null variance.
#'
#' @references
#' Bai, Z. and Saranadasa, H. (1996). Effect of high dimension: by an example
#' of a two sample problem. *Statistica Sinica*, 6, 311--329.
#'
#' @examples
#' set.seed(22)
#' x <- matrix(rnorm(24), nrow = 8, ncol = 3)
#' y <- matrix(rnorm(27, mean = 0.2), nrow = 9, ncol = 3)
#' bai_saranadasa_two_sample_test(x, y)
#'
#' @export
bai_saranadasa_two_sample_test <- function(x, y) {
  call <- match.call()
  x.name <- deparse(substitute(x))
  y.name <- deparse(substitute(y))
  x <- .as_data_matrix(x, "x", min_rows = 2L)
  y <- .as_data_matrix(y, "y", min_rows = 2L)
  .check_two_sample_variables(x, y)

  n1 <- nrow(x)
  n2 <- nrow(y)
  p <- ncol(x)
  out <- cpp_bai_saranadasa_two_sample(x, y)
  z <- unname(out$z)
  p.value <- stats::pnorm(z, lower.tail = FALSE)

  variable.names <- .hotelling_variable_names(x)
  mean.x <- stats::setNames(as.numeric(out$mean_x), variable.names)
  mean.y <- stats::setNames(as.numeric(out$mean_y), variable.names)
  difference <- stats::setNames(as.numeric(out$difference), variable.names)

  .new_hd_location_test(
    statistic = c(Z = z),
    p.value = p.value,
    method = paste(
      "Bai-Saranadasa two-sample high-dimensional mean test",
      "(common covariance; asymptotic normal calibration)"
    ),
    data.name = paste(x.name, "and", y.name),
    alternative = "two.sided",
    raw.statistic = c(xi = unname(out$xi)),
    estimate = difference,
    null.value = stats::setNames(numeric(p), variable.names),
    null.distribution = list(
      family = "normal",
      parameters = c(mean = 0, sd = 1),
      exact = FALSE,
      tail = "upper",
      assumptions = paste(
        "Bai-Saranadasa high-dimensional asymptotics;",
        "independent samples with a common covariance matrix"
      )
    ),
    variance = c(denominator = unname(out$variance)),
    components = list(
      mean.x = mean.x,
      mean.y = mean.y,
      difference = difference,
      n1 = unname(out$n1),
      n2 = unname(out$n2),
      N = unname(out$N),
      nu = unname(out$nu),
      p = unname(out$p),
      a = unname(out$a),
      trace.Sp = unname(out$trace_Sp),
      trace.Sp2 = unname(out$trace_Sp2),
      B = unname(out$B),
      trace.Sigma2.hat = unname(out$trace_sigma2_hat),
      xi = unname(out$xi),
      denominator.variance = unname(out$variance)
    ),
    diagnostics = list(
      covariance.model = "common",
      pooled.covariance.denominator = unname(out$nu),
      trace.computation = paste(unname(out$gram_type), "Gram matrix"),
      gram.type = unname(out$gram_type),
      gram.dimension = unname(out$gram_dimension),
      constructs.p.by.p.matrix = identical(unname(out$gram_type), "primal"),
      internal.scaling = "common global scaling of both groups",
      internal.scale.factor = unname(out$global_scale),
      internal.scaled.xi = unname(out$xi_scaled),
      internal.scaled.denominator.variance = unname(out$variance_scaled),
      regularization = "none",
      variance.repair = "none"
    ),
    n = c(group1 = n1, group2 = n2),
    p = p,
    call = call
  )
}
