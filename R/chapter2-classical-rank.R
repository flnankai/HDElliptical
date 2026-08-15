.classical_location_chisq <- function(fit, method, data.name, estimate,
                                      null.value, variance, components,
                                      diagnostics, n, p, call,
                                      assumptions) {
  statistic <- as.numeric(fit$statistic)
  df <- c(df = p)
  p_value <- stats::pchisq(statistic, df = p, lower.tail = FALSE)

  .new_hd_location_test(
    statistic = c(`X-squared` = statistic),
    parameter = df,
    p.value = p_value,
    method = method,
    data.name = data.name,
    raw.statistic = c(Q = statistic),
    estimate = estimate,
    null.value = null.value,
    null.distribution = list(
      family = "chi-squared",
      parameters = df,
      exact = FALSE,
      approximation = "fixed-p asymptotic",
      assumptions = assumptions
    ),
    variance = variance,
    components = components,
    diagnostics = c(
      list(
        rank = as.integer(fit$rank),
        rcond = as.numeric(fit$rcond),
        solver = as.character(fit$solver),
        calibration = "fixed-p asymptotic"
      ),
      diagnostics
    ),
    n = n,
    p = p,
    call = call
  )
}

.classical_zero_tolerance <- function(zero_tol) {
  zero_tol <- as.numeric(zero_tol)
  if (length(zero_tol) != 1L || is.na(zero_tol) ||
      !is.finite(zero_tol) || zero_tol < 0) {
    stop("`zero_tol` must be a finite non-negative number.", call. = FALSE)
  }
  zero_tol
}


#' Classical spatial sign and rank location tests
#'
#' Implements the fixed-dimensional spatial-sign test, the one-sample spatial
#' signed-rank test, and the pooled two-sample spatial-rank test. All three use
#' the zero-direction convention \eqn{U(0)=0} and a Cholesky solve; a singular
#' studentizing matrix is reported rather than replaced by a generalized
#' inverse.
#'
#' For centered observations \eqn{Z_i=X_i-\mu_0}, the sign statistic uses
#' \deqn{\bar U=n^{-1}\sum_i U(Z_i),\quad
#' B_U=n^{-1}\sum_i U(Z_i)U(Z_i)^T,\quad
#' Q_{sign}=n\bar U^T B_U^{-1}\bar U.}
#' The signed ranks are
#' \deqn{R_i=n^{-1}\sum_j U(Z_i+Z_j).}
#' Since the first-order projection of their average has covariance
#' \eqn{4B_R}, where \eqn{B_R=n^{-1}\sum_i R_iR_i^T}, its statistic is
#' \deqn{Q_{SR}=\frac{n}{4}\bar R^T B_R^{-1}\bar R.}
#' The factor `1/4` is essential.
#'
#' For two samples, every observation is ranked against the pooled sample:
#' \deqn{R(Y_i)=N^{-1}\sum_j U(Y_i-Y_j).}
#' If \eqn{C=(N-1)^{-1}\sum_i R(Y_i)R(Y_i)^T}, the sample covariance of
#' the pooled ranks (whose average is exactly zero), the statistic is
#' \deqn{Q_{2SR}=\frac{n_1n_2}{N}(\bar R_1-\bar R_2)^T
#' C^{-1}(\bar R_1-\bar R_2).}
#' These calibrations are asymptotic chi-squared laws for fixed dimension, not
#' high-dimensional approximations. Euclidean spatial signs make the tests
#' invariant to translations, common positive rescaling, and orthogonal
#' transformations, but not to arbitrary nonspherical affine transformations.
#'
#' @param x A numeric matrix or data frame with observations in rows.
#' @param mu For a one-sample test, a finite null location with one value per
#'   column of `x`. `NULL` uses the zero vector.
#' @param y For `spatial_rank_test()`, a second numeric matrix or data frame
#'   with the same variables as `x`.
#' @param tol Reciprocal-condition-number tolerance for the Cholesky solve. It
#'   must be strictly between zero and one.
#' @param zero_tol Non-negative tolerance, in the original data units, below
#'   which a residual, pair sum, or pair difference is mapped to zero.
#'
#' @return An object of classes `hd_location_test` and `htest`. `statistic` is
#'   the chi-squared statistic, `parameter` is its asymptotic degrees of
#'   freedom, `variance` is the studentizing second-moment matrix (or `4 B_R`
#'   for the signed-rank test), and `components` contains the directional mean
#'   and unscaled moment matrices. Numerical rank, reciprocal condition number,
#'   solver, zero counts, and calibration type are in `diagnostics`.
#'
#' @references
#' Mottonen, J. and Oja, H. (1995). Multivariate spatial sign and rank methods.
#' *Journal of Nonparametric Statistics*, **5**, 201--213.
#'
#' Oja, H. (2010). *Multivariate Nonparametric Methods with R*. Springer.
#'
#' @examples
#' set.seed(12)
#' x <- matrix(rnorm(80), 20, 4)
#' spatial_sign_test(x)
#' spatial_signed_rank_test(x)
#'
#' y <- matrix(rnorm(96, 0.25), 24, 4)
#' spatial_rank_test(x, y)
#'
#' @name classical_spatial_tests
NULL


#' @rdname classical_spatial_tests
#' @export
spatial_sign_test <- function(x, mu = NULL,
                              tol = sqrt(.Machine$double.eps),
                              zero_tol = 0) {
  call <- match.call()
  data_name <- deparse1(substitute(x))
  x <- .as_data_matrix(x, min_rows = 2L)
  n <- nrow(x)
  p <- ncol(x)
  if (n <= p) {
    stop(sprintf(
      "The classical spatial-sign test requires `n > p`; received n = %d and p = %d.",
      n, p
    ), call. = FALSE)
  }
  tol <- .hotelling_tolerance(tol)
  if (is.null(mu)) {
    mu <- rep.int(0, p)
  } else {
    mu <- .as_location(mu, p, "mu")
  }
  zero_tol <- .classical_zero_tolerance(zero_tol)
  fit <- cpp_spatial_sign_test(x, mu, tol, zero_tol)

  variable_names <- .hotelling_variable_names(x)
  names(mu) <- variable_names
  mean_sign <- as.numeric(fit$mean_sign)
  names(mean_sign) <- variable_names
  second_moment <- as.matrix(fit$second_moment)
  dimnames(second_moment) <- list(variable_names, variable_names)

  .classical_location_chisq(
    fit = fit,
    method = "One-sample classical spatial-sign test",
    data.name = data_name,
    estimate = mean_sign,
    null.value = mu,
    variance = second_moment,
    components = list(
      mean.sign = mean_sign,
      sign.second.moment = second_moment,
      hypothesized.location = mu
    ),
    diagnostics = list(
      tolerance = tol,
      zero.tolerance = zero_tol,
      n.zero.signs = as.numeric(fit$n_zero)
    ),
    n = c(x = n),
    p = p,
    call = call,
    assumptions = c(
      "independent observations",
      "central symmetry about the hypothesized location",
      "fixed dimension with nonsingular sign second moment"
    )
  )
}


#' @rdname classical_spatial_tests
#' @export
spatial_signed_rank_test <- function(x, mu = NULL,
                                     tol = sqrt(.Machine$double.eps),
                                     zero_tol = 0) {
  call <- match.call()
  data_name <- deparse1(substitute(x))
  x <- .as_data_matrix(x, min_rows = 2L)
  n <- nrow(x)
  p <- ncol(x)
  if (n <= p) {
    stop(sprintf(
      paste0(
        "The classical spatial signed-rank test requires `n > p`; ",
        "received n = %d and p = %d."
      ),
      n, p
    ), call. = FALSE)
  }
  tol <- .hotelling_tolerance(tol)
  if (is.null(mu)) {
    mu <- rep.int(0, p)
  } else {
    mu <- .as_location(mu, p, "mu")
  }
  zero_tol <- .classical_zero_tolerance(zero_tol)
  fit <- cpp_spatial_signed_rank_test(x, mu, tol, zero_tol)

  variable_names <- .hotelling_variable_names(x)
  names(mu) <- variable_names
  mean_rank <- as.numeric(fit$mean_rank)
  names(mean_rank) <- variable_names
  rank_second_moment <- as.matrix(fit$rank_second_moment)
  dimnames(rank_second_moment) <- list(variable_names, variable_names)
  asymptotic_covariance <- 4 * rank_second_moment

  .classical_location_chisq(
    fit = fit,
    method = "One-sample classical spatial signed-rank test",
    data.name = data_name,
    estimate = mean_rank,
    null.value = mu,
    variance = asymptotic_covariance,
    components = list(
      mean.signed.rank = mean_rank,
      rank.second.moment = rank_second_moment,
      asymptotic.covariance = asymptotic_covariance,
      quadratic.multiplier = n / 4,
      hypothesized.location = mu
    ),
    diagnostics = list(
      tolerance = tol,
      zero.tolerance = zero_tol,
      n.zero.ordered.pairs = as.numeric(fit$n_zero_ordered_pairs),
      n.ordered.pairs = as.numeric(n)^2,
      hoeffding.factor = 4
    ),
    n = c(x = n),
    p = p,
    call = call,
    assumptions = c(
      "independent observations",
      "central symmetry about the hypothesized location",
      "fixed dimension with nonsingular signed-rank second moment"
    )
  )
}


#' @rdname classical_spatial_tests
#' @export
spatial_rank_test <- function(x, y,
                              tol = sqrt(.Machine$double.eps),
                              zero_tol = 0) {
  call <- match.call()
  x_name <- deparse1(substitute(x))
  y_name <- deparse1(substitute(y))
  x <- .as_data_matrix(x, "x", min_rows = 2L)
  y <- .as_data_matrix(y, "y", min_rows = 2L)
  .check_two_sample_variables(x, y)
  n1 <- nrow(x)
  n2 <- nrow(y)
  total <- n1 + n2
  p <- ncol(x)
  if (total <= p) {
    stop(sprintf(
      paste0(
        "The pooled spatial-rank test requires `n1 + n2 > p`; ",
        "received n1 + n2 = %d and p = %d."
      ),
      total, p
    ), call. = FALSE)
  }
  tol <- .hotelling_tolerance(tol)
  zero_tol <- .classical_zero_tolerance(zero_tol)
  fit <- cpp_spatial_rank_test(x, y, tol, zero_tol)

  variable_names <- .hotelling_variable_names(x)
  mean_rank_x <- as.numeric(fit$mean_rank_x)
  mean_rank_y <- as.numeric(fit$mean_rank_y)
  difference <- as.numeric(fit$difference)
  names(mean_rank_x) <- names(mean_rank_y) <- names(difference) <-
    variable_names
  zero_difference <- stats::setNames(rep.int(0, p), variable_names)
  rank_covariance <- as.matrix(fit$rank_covariance)
  dimnames(rank_covariance) <- list(variable_names, variable_names)

  .classical_location_chisq(
    fit = fit,
    method = "Two-sample pooled spatial-rank test",
    data.name = paste(x_name, "and", y_name),
    estimate = difference,
    null.value = zero_difference,
    variance = rank_covariance,
    components = list(
      mean.rank.x = mean_rank_x,
      mean.rank.y = mean_rank_y,
      difference = difference,
      pooled.rank.covariance = rank_covariance,
      quadratic.multiplier = n1 * n2 / total
    ),
    diagnostics = list(
      tolerance = tol,
      zero.tolerance = zero_tol,
      n.zero.ordered.pairs = as.numeric(fit$n_zero_ordered_pairs),
      n.ordered.pairs = as.numeric(total)^2,
      includes.zero.self.differences = TRUE,
      covariance.denominator = total - 1
    ),
    n = c(x = n1, y = n2),
    p = p,
    call = call,
    assumptions = c(
      "independent samples from a common continuous distribution under the null",
      "fixed dimension with nonsingular pooled-rank covariance"
    )
  )
}
