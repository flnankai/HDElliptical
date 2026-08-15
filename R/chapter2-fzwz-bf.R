#' Feng--Zou--Wang--Zhu scale-invariant Behrens--Fisher test
#'
#' Tests equality of two high-dimensional mean vectors without requiring equal
#' covariance matrices.  Observations are rows and variables are columns.  Let
#' \eqn{\gamma=n_1/n_2}, let \eqn{\widehat\sigma_{sk}^2} be the unbiased
#' sample variance of variable \eqn{k} in group \eqn{s}, and define
#' \deqn{A_k=(\bar X_{1k}-\bar X_{2k})^2-
#'   \widehat\sigma_{1k}^2/n_1-\widehat\sigma_{2k}^2/n_2.}
#' The initial statistic is
#' \deqn{T_{BF}=\sum_k
#'   \frac{A_k}{\widehat\sigma_{1k}^2+
#'   \gamma\widehat\sigma_{2k}^2}.}
#'
#' The reported statistic subtracts the Feng--Zou--Wang--Zhu plug-in estimate
#' of the asymptotic null centering and divides by their ratio-consistent null
#' standard deviation.  If
#' \eqn{\widehat\kappa_{sk}=n_s^{-1}\sum_i
#' (X_{sik}-\bar X_{sk})^3} and
#' \eqn{D_k=\widehat\sigma_{1k}^2+
#' \gamma\widehat\sigma_{2k}^2}, the centering is
#' \deqn{\widehat\mu_n=\widehat b_1+\widehat b_2,}
#' with
#' \deqn{\widehat b_1=\sum_k\left\{
#' \frac{2\widehat\sigma_{1k}^4}{n_1(n_1-1)D_k^2}+
#' \frac{2\gamma\widehat\sigma_{2k}^4}{n_2(n_2-1)D_k^2}
#' \right\},}
#' and
#' \deqn{\widehat b_2=\sum_k\frac{2}{D_k^3}
#' \left(\frac{\widehat\kappa_{1k}}{n_1}-
#' \frac{\gamma\widehat\kappa_{2k}}{n_2}\right)^2.}
#' In particular, the second term of \eqn{\widehat b_1} contains one power of
#' \eqn{\gamma}.  This follows the published paper; the current book draft at
#' `ch2_location.tex:685` has an extra power of \eqn{\gamma}.  The displayed
#' population expansion in that draft is also an asymptotic centering up to a
#' smaller-order remainder, not an exact finite-sample expectation.
#'
#' The two within-group covariance trace estimates use exactly the paper's
#' leave-four-out marginal variances and denominator \eqn{2P_{n_s}^4}.  The
#' cross trace removes two observations from each group and uses denominator
#' \eqn{4P_{n_1}^2P_{n_2}^2}.  At least six observations are consequently
#' required in each group.  The C++ kernel sums the 24 orderings associated
#' with each unordered within-group quadruple from a local \eqn{4\times4}
#' dual Gram matrix.  This reduces computation without changing the published
#' estimator and never constructs a \eqn{p\times p} matrix.
#'
#' Both samples are internally translated by a common anchor and divided by a
#' common scale within each variable.  This is algebraically neutral and
#' protects the statistic under very large or small measurement units.  Every
#' full-sample and leave-out combined marginal variance, and the final null
#' variance estimate, must have the sign required by the original formulas.
#' No ridge, absolute-value repair, or variance floor is applied.
#'
#' @param x,y Numeric matrices or data frames with observations in rows and
#'   the same variables in columns.  Each group must contain at least six
#'   observations.
#'
#' @return An object of class `c("hd_location_test", "htest")`.  The
#'   upper-tail normal p-value is based on
#'   \eqn{Z=(T_{BF}-\widehat\mu_n)/\widehat\sigma_n}.  `components` contains
#'   `T.BF`, `Q3 = n1 * T.BF`, the two centering terms, all three published
#'   leave-out trace estimates, their sample-size coefficients and permutation
#'   counts, and coordinate-level quantities.  `diagnostics` records the
#'   exact leave-out policy and the absence of numerical repair.
#'
#' @references
#' Feng, L., Zou, C., Wang, Z., and Zhu, L. (2015). Two-sample
#' Behrens--Fisher problem for high-dimensional data. *Statistica Sinica*,
#' **25**, 1297--1312. \doi{10.5705/ss.2014.048}.
#'
#' @examples
#' set.seed(2411)
#' x <- matrix(rnorm(56), 8, 7)
#' y <- matrix(rnorm(63, 0.15), 9, 7)
#' feng_zou_wang_zhu_two_sample_test(x, y)
#'
#' @export
feng_zou_wang_zhu_two_sample_test <- function(x, y) {
  call <- match.call()
  x.name <- deparse1(substitute(x))
  y.name <- deparse1(substitute(y))
  x <- .as_data_matrix(x, "x", min_rows = 6L)
  y <- .as_data_matrix(y, "y", min_rows = 6L)
  .check_two_sample_variables(x, y)

  n1 <- nrow(x)
  n2 <- nrow(y)
  p <- ncol(x)
  fit <- cpp_fzwz_bf_two_sample(x, y)
  z <- as.numeric(fit$z)
  p.value <- stats::pnorm(z, lower.tail = FALSE)
  variable.names <- .hotelling_variable_names(x)

  named.vector <- function(value) {
    stats::setNames(as.numeric(value), variable.names)
  }
  mean.x <- named.vector(fit$mean_x)
  mean.y <- named.vector(fit$mean_y)
  difference <- named.vector(fit$difference)
  variance1.diagonal <- named.vector(fit$variance1_diagonal)
  variance2.diagonal <- named.vector(fit$variance2_diagonal)
  d.hat.diagonal <- named.vector(fit$D_hat_diagonal)
  a.coordinate <- named.vector(fit$A_coordinate)
  a.coordinate.scaled <- named.vector(fit$A_coordinate_scaled)
  coordinate.contribution <- named.vector(fit$coordinate_contribution)

  .new_hd_location_test(
    statistic = c(Z = z),
    p.value = p.value,
    method = paste(
      "Feng-Zou-Wang-Zhu scale-invariant high-dimensional",
      "Behrens-Fisher test (asymptotic normal calibration)"
    ),
    data.name = paste(x.name, "and", y.name),
    alternative = "two.sided",
    raw.statistic = c(T.BF = as.numeric(fit$T_BF)),
    estimate = difference,
    null.value = stats::setNames(numeric(p), variable.names),
    null.distribution = list(
      family = "normal",
      parameters = c(mean = 0, sd = 1),
      exact = FALSE,
      tail = "upper",
      assumptions = paste(
        "Feng-Zou-Wang-Zhu high-dimensional asymptotics; independent",
        "samples; finite eighth moments; strictly positive full and",
        "leave-out combined marginal variances"
      )
    ),
    variance = c(denominator = as.numeric(fit$variance)),
    components = list(
      mean.x = mean.x,
      mean.y = mean.y,
      difference = difference,
      variance1.diagonal = variance1.diagonal,
      variance2.diagonal = variance2.diagonal,
      D.hat.diagonal = d.hat.diagonal,
      A.coordinate = a.coordinate,
      A.coordinate.scaled = a.coordinate.scaled,
      coordinate.contribution = coordinate.contribution,
      gamma = as.numeric(fit$gamma),
      n1 = as.numeric(fit$n1),
      n2 = as.numeric(fit$n2),
      p = as.numeric(fit$p),
      Q3 = as.numeric(fit$Q3),
      T.BF = as.numeric(fit$T_BF),
      null.centering = as.numeric(fit$null_centering),
      b1 = as.numeric(fit$bias1),
      b2 = as.numeric(fit$bias2),
      trace1 = as.numeric(fit$trace1),
      trace2 = as.numeric(fit$trace2),
      trace12 = as.numeric(fit$trace12),
      variance.coefficients = fit$variance_coefficients,
      denominator.variance = as.numeric(fit$variance),
      P4.n1 = as.numeric(fit$P4_n1),
      P4.n2 = as.numeric(fit$P4_n2),
      P2.n1 = as.numeric(fit$P2_n1),
      P2.n2 = as.numeric(fit$P2_n2),
      within1.combinations = as.numeric(fit$within1_combinations),
      within2.combinations = as.numeric(fit$within2_combinations),
      cross.pair.combinations =
        as.numeric(fit$cross_pair_combinations)
    ),
    diagnostics = list(
      covariance.model = "unequal",
      sample.variance.denominators = c(
        group1 = n1 - 1,
        group2 = n2 - 1
      ),
      sample.third.moment.denominators = c(
        group1 = n1,
        group2 = n2
      ),
      b1.group2.gamma.power = 1,
      within.trace.estimator = paste(
        "published leave-four-out estimator; denominator 2 * P_n^4;",
        "24 ordered kernels per unordered quadruple"
      ),
      cross.trace.estimator = paste(
        "published two-plus-two leave-out estimator; denominator",
        "4 * P_n1^2 * P_n2^2"
      ),
      trace.computation = "local 4 by 4 dual Gram symmetry reduction",
      constructs.p.by.p.matrix = FALSE,
      minimum.full.D.scaled =
        as.numeric(fit$minimum_full_D_scaled),
      minimum.leaveout.D.scaled =
        as.numeric(fit$minimum_leaveout_D_scaled),
      internal.scaling = paste(
        "common per-column anchor and scale across both groups;",
        "long-double computational moments"
      ),
      internal.scale.factors = named.vector(fit$column_scale),
      regularization = "none",
      variance.repair = "none"
    ),
    n = c(group1 = n1, group2 = n2),
    p = p,
    call = call
  )
}
