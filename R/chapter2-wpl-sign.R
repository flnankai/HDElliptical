#' Wang--Peng--Li one-sample high-dimensional spatial-sign test
#'
#' Tests \eqn{H_0:\mu=\mu_0} against an unrestricted location alternative
#' using the raw spatial-sign statistic of Wang, Peng, and Li (2015).  For
#' \eqn{Z_i=U(X_i-\mu_0)}, where \eqn{U(v)=v/\lVert v\rVert} for nonzero
#' \eqn{v} and \eqn{U(0)=0}, the unstandardised statistic is
#' \deqn{T_n=\sum_{1\leq i<j\leq n} Z_i^\mathsf{T}Z_j.}
#'
#' The unknown \eqn{\operatorname{tr}(B^2)} in
#' \eqn{\operatorname{Var}(T_n)=n(n-1)\operatorname{tr}(B^2)/2} is estimated
#' by the feasible cross-validation estimator in equation (7) of the original
#' paper:
#' \deqn{\widehat{\operatorname{tr}(B^2)}=\frac{1}{n(n-1)}
#' \sum_{j\ne k}\operatorname{tr}\{(Z_j-\bar Z_{(j,k)})Z_j^\mathsf{T}
#' (Z_k-\bar Z_{(j,k)})Z_k^\mathsf{T}\},}
#' where \eqn{\bar Z_{(j,k)}} is the sign mean after deleting observations
#' \eqn{j} and \eqn{k}.  The implementation evaluates this literal
#' leave-two-out formula.  It does not substitute the paper's equation (8),
#' whose shortcut uses \eqn{\lVert Z_i\rVert^2=1}; consequently the documented
#' \eqn{U(0)=0} convention remains coherent even when an observation equals
#' the null location exactly.
#'
#' Residual directions are normalised after max-absolute scaling, and
#' subtraction is prescaled only when direct finite subtraction overflows.
#' Thus common nonzero rescaling of all residuals leaves the calculation
#' unchanged even at extreme finite units.  Exact null residuals contribute a
#' zero sign and are counted in `diagnostics`.  A non-finite or non-positive
#' cross-validation variance estimate is a genuine degenerate calibration and
#' raises an error; no absolute value, ridge, or numerical floor is used.
#'
#' The p-value uses the upper tail of the WPL asymptotic standard normal law.
#' Although the vector alternative is conventionally labelled `two.sided`,
#' large positive quadratic evidence is the rejection direction.  The
#' calibration requires the high-dimensional trace and concentration
#' conditions in Wang, Peng, and Li (2015); it is not an exact finite-sample
#' test.
#'
#' @param x A numeric matrix or data frame with observations in rows and
#'   variables in columns.  At least three observations are required.
#' @param mu A finite numeric vector giving the null location.  The default is
#'   the zero vector.
#'
#' @return An object of class `c("hd_location_test", "htest")`.  Its raw
#'   components include `T.WPL`, the literal cross-validation numerator,
#'   `trace.B2.hat`, its estimated variance and standard error, and the sum and
#'   mean of the spatial signs.  `diagnostics` records exact zero signs,
#'   overflow-safe subtraction fallbacks, and the no-repair variance policy.
#'
#' @references
#' Wang, L., Peng, B., and Li, R. (2015). A high-dimensional nonparametric
#' multivariate test for mean vector. *Journal of the American Statistical
#' Association*, **110**, 1658--1669.
#' \doi{10.1080/01621459.2014.988215}.
#'
#' @examples
#' set.seed(2601)
#' x <- matrix(stats::rt(240, df = 4), 20, 12)
#' wang_peng_li_one_sample_test(x)
#'
#' @export
wang_peng_li_one_sample_test <- function(x, mu = NULL) {
  call <- match.call()
  data.name <- deparse1(substitute(x))
  x <- .as_data_matrix(x, "x", min_rows = 3L)
  n <- nrow(x)
  p <- ncol(x)
  if (is.null(mu)) {
    mu <- numeric(p)
  } else {
    mu <- .as_location(mu, p, "mu")
  }

  fit <- cpp_wang_peng_li_one_sample(x, mu)
  z <- unname(fit$z)
  p.value <- stats::pnorm(z, lower.tail = FALSE)

  variable.names <- .hotelling_variable_names(x)
  mu.named <- stats::setNames(as.numeric(mu), variable.names)
  sample.mean <- stats::setNames(
    as.numeric(fit$sample_mean), variable.names
  )
  direction.sum <- stats::setNames(
    as.numeric(fit$direction_sum), variable.names
  )
  direction.mean <- stats::setNames(
    as.numeric(fit$direction_mean), variable.names
  )

  .new_hd_location_test(
    statistic = c(Z = z),
    p.value = p.value,
    method = paste(
      "Wang-Peng-Li one-sample high-dimensional spatial-sign test",
      "(literal leave-two-out variance; asymptotic normal calibration)"
    ),
    data.name = data.name,
    alternative = "two.sided",
    raw.statistic = c(T.WPL = unname(fit$T_raw)),
    estimate = sample.mean,
    null.value = mu.named,
    null.distribution = list(
      family = "normal",
      parameters = c(mean = 0, sd = 1),
      exact = FALSE,
      tail = "upper",
      assumptions = paste(
        "Wang-Peng-Li high-dimensional trace/concentration conditions;",
        "independent observations; positive leave-two-out variance estimate"
      )
    ),
    variance = c(estimated = unname(fit$variance_hat)),
    components = list(
      sample.mean = sample.mean,
      null.location = mu.named,
      direction.sum = direction.sum,
      direction.mean = direction.mean,
      sum.direction.norm.squared =
        unname(fit$sum_direction_norm_squared),
      T.WPL = unname(fit$T_raw),
      cross.validation.numerator =
        unname(fit$cross_validation_numerator),
      trace.B2.hat = unname(fit$trace_B2_hat),
      variance.hat = unname(fit$variance_hat),
      standard.error = unname(fit$standard_error),
      pair.count = unname(fit$pair_count),
      ordered.pair.count = unname(fit$ordered_pair_count),
      n = unname(fit$n),
      p = unname(fit$p)
    ),
    diagnostics = list(
      sign.at.zero = "U(0) = 0",
      zero.signs = unname(fit$zero_signs),
      nonzero.signs = unname(fit$nonzero_signs),
      zero.sign.proportion = unname(fit$zero_signs) / n,
      variance.estimator = paste(
        "Wang-Peng-Li (2015), equation (7), literal ordered-pair",
        "leave-two-out cross-validation formula"
      ),
      equation.8.shortcut.used = FALSE,
      subtraction.overflow.fallback.rows =
        unname(fit$overflow_fallback_rows),
      direction.normalization = "rowwise max-absolute scaling",
      regularization = "none",
      variance.repair = "none"
    ),
    n = n,
    p = p,
    call = call
  )
}
