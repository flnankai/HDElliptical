.new_hd_location_test <- function(statistic, parameter = numeric(), p.value,
                                  method, data.name,
                                  alternative = "two.sided",
                                  raw.statistic = NULL, estimate = NULL,
                                  null.value = NULL,
                                  null.distribution = list(),
                                  variance = NULL, components = list(),
                                  diagnostics = list(), n, p, call) {
  statistic_name <- names(statistic)
  statistic <- as.numeric(statistic)
  if (length(statistic) != 1L || !is.finite(statistic)) {
    stop("`statistic` must be one finite numeric value.", call. = FALSE)
  }
  if (is.null(statistic_name) || !nzchar(statistic_name[[1L]])) {
    statistic_name <- "statistic"
  }
  names(statistic) <- statistic_name[[1L]]

  parameter_names <- names(parameter)
  parameter <- as.numeric(parameter)
  if (anyNA(parameter) || any(!is.finite(parameter))) {
    stop("`parameter` must contain only finite numeric values.", call. = FALSE)
  }
  names(parameter) <- parameter_names

  p.value <- as.numeric(p.value)
  if (length(p.value) != 1L || !is.finite(p.value) ||
      p.value < 0 || p.value > 1) {
    stop("`p.value` must be one finite number in [0, 1].", call. = FALSE)
  }
  if (!is.character(method) || length(method) != 1L || !nzchar(method) ||
      !is.character(data.name) || length(data.name) != 1L ||
      !nzchar(data.name)) {
    stop("`method` and `data.name` must be non-empty strings.", call. = FALSE)
  }
  if (!identical(alternative, "two.sided")) {
    stop("Multivariate location tests currently use `alternative = \"two.sided\"`.",
         call. = FALSE)
  }
  if (!is.list(null.distribution) ||
      !is.character(null.distribution$family) ||
      length(null.distribution$family) != 1L) {
    stop("`null.distribution` must be a list with a scalar `family` field.",
         call. = FALSE)
  }
  if (!is.list(components) || !is.list(diagnostics)) {
    stop("`components` and `diagnostics` must be lists.", call. = FALSE)
  }
  n_names <- names(n)
  n <- as.numeric(n)
  names(n) <- n_names
  p <- as.integer(p)
  if (length(n) < 1L || anyNA(n) || any(!is.finite(n)) || any(n < 1) ||
      length(p) != 1L || is.na(p) || p < 1L) {
    stop("`n` and `p` must describe positive sample and dimension sizes.",
         call. = FALSE)
  }

  result <- list(
    statistic = statistic,
    parameter = parameter,
    p.value = p.value,
    method = method,
    data.name = data.name,
    alternative = alternative,
    raw.statistic = raw.statistic,
    estimate = estimate,
    null.value = null.value,
    null.distribution = null.distribution,
    variance = variance,
    components = components,
    diagnostics = diagnostics,
    n = n,
    p = p,
    call = call
  )
  class(result) <- c("hd_location_test", "htest")
  result
}

.hotelling_tolerance <- function(tol) {
  tol <- as.numeric(tol)
  if (length(tol) != 1L || !is.finite(tol) || tol <= 0 || tol >= 1) {
    stop("`tol` must be a finite number strictly between zero and one.",
         call. = FALSE)
  }
  tol
}

.hotelling_variable_names <- function(x) {
  if (is.null(colnames(x))) {
    paste0("V", seq_len(ncol(x)))
  } else {
    colnames(x)
  }
}

.check_two_sample_variables <- function(x, y) {
  if (ncol(x) != ncol(y)) {
    stop("`x` and `y` must have the same number of columns.", call. = FALSE)
  }
  if (!is.null(colnames(x)) && !is.null(colnames(y)) &&
      !identical(colnames(x), colnames(y))) {
    stop("When both samples have column names, `x` and `y` must use the " %+%
           "same names in the same order.", call. = FALSE)
  }
  invisible(NULL)
}

#' Classical Hotelling location tests
#'
#' Performs the exact Gaussian one- or two-sample Hotelling
#' \eqn{T^2} test. Observations are rows and variables are columns. The
#' two-sample test uses the pooled covariance matrix and therefore assumes a
#' common population covariance matrix.
#'
#' The covariance system is solved by a Cholesky factorization. A generalized
#' inverse is deliberately not used: singular or numerically ill-conditioned
#' covariance matrices invalidate the exact \eqn{F} calibration and produce an
#' error directing the user to a high-dimensional test.
#'
#' @param x A numeric matrix or data frame with observations in rows.
#' @param mu For the one-sample test, a finite null mean vector with one value
#'   per column of `x`. `NULL` uses the zero vector.
#' @param y For the two-sample test, a second numeric matrix or data frame with
#'   observations in rows and the same variables as `x`.
#' @param tol Reciprocal-condition-number tolerance used after Cholesky
#'   factorization. It must be strictly between zero and one.
#'
#' @return An object of classes `hd_location_test` and `htest`. In addition to
#'   the standard `htest` fields, it contains the raw Hotelling `T2`, the
#'   covariance estimate, exact null-distribution metadata, and numerical
#'   diagnostics including the reciprocal condition number and solver.
#'
#' @references
#' Hotelling, H. (1931). The generalization of Student's ratio. *Annals of
#' Mathematical Statistics*, **2**, 360--378.
#'
#' @examples
#' set.seed(11)
#' x <- matrix(rnorm(60), 20, 3)
#' hotelling_one_sample_test(x)
#'
#' y <- matrix(rnorm(75, 0.25), 25, 3)
#' hotelling_two_sample_test(x, y)
#'
#' @name hotelling_tests
NULL

#' @rdname hotelling_tests
#' @export
hotelling_one_sample_test <- function(x, mu = NULL,
                                      tol = sqrt(.Machine$double.eps)) {
  call <- match.call()
  data_name <- deparse1(substitute(x))
  x <- .as_data_matrix(x, min_rows = 2L)
  n <- nrow(x)
  p <- ncol(x)
  if (n <= p) {
    stop(sprintf(
      "One-sample Hotelling's test requires `n > p`; received n = %d and p = %d.",
      n, p
    ), call. = FALSE)
  }
  tol <- .hotelling_tolerance(tol)
  if (is.null(mu)) {
    mu <- rep.int(0, p)
  } else {
    mu <- .as_location(mu, p, "mu")
  }

  fit <- cpp_hotelling_one_sample(x, mu, tol)
  t_squared <- as.numeric(fit$t_squared)
  scale_factor <- (n - p) / (p * (n - 1))
  f_statistic <- scale_factor * t_squared
  df <- c(df1 = p, df2 = n - p)
  p_value <- stats::pf(f_statistic, df[["df1"]], df[["df2"]],
                       lower.tail = FALSE)

  variable_names <- .hotelling_variable_names(x)
  sample_mean <- as.numeric(fit$mean)
  difference <- as.numeric(fit$difference)
  names(sample_mean) <- variable_names
  names(difference) <- variable_names
  names(mu) <- variable_names
  covariance <- as.matrix(fit$covariance)
  dimnames(covariance) <- list(variable_names, variable_names)

  .new_hd_location_test(
    statistic = c(F = f_statistic),
    parameter = df,
    p.value = p_value,
    method = "One-sample Hotelling's T-squared test",
    data.name = data_name,
    raw.statistic = c(T2 = t_squared),
    estimate = sample_mean,
    null.value = mu,
    null.distribution = list(
      family = "F",
      parameters = df,
      exact = TRUE,
      assumptions = c("multivariate normality", "positive definite covariance")
    ),
    variance = covariance,
    components = list(
      sample.mean = sample_mean,
      difference = difference,
      covariance = covariance,
      f.scale = scale_factor
    ),
    diagnostics = list(
      rank = as.integer(fit$rank),
      rcond = as.numeric(fit$rcond),
      solver = as.character(fit$solver),
      tolerance = tol,
      exact.calibration = TRUE
    ),
    n = c(x = n),
    p = p,
    call = call
  )
}

#' @rdname hotelling_tests
#' @export
hotelling_two_sample_test <- function(x, y,
                                      tol = sqrt(.Machine$double.eps)) {
  call <- match.call()
  x_name <- deparse1(substitute(x))
  y_name <- deparse1(substitute(y))
  x <- .as_data_matrix(x, "x", min_rows = 2L)
  y <- .as_data_matrix(y, "y", min_rows = 2L)
  .check_two_sample_variables(x, y)
  n1 <- nrow(x)
  n2 <- nrow(y)
  p <- ncol(x)
  total <- n1 + n2
  if (total <= p + 1L) {
    stop(sprintf(
      paste0(
        "Two-sample Hotelling's test requires `n1 + n2 > p + 1`; ",
        "received n1 + n2 = %d and p = %d."
      ),
      total, p
    ), call. = FALSE)
  }
  tol <- .hotelling_tolerance(tol)

  fit <- cpp_hotelling_two_sample(x, y, tol)
  t_squared <- as.numeric(fit$t_squared)
  scale_factor <- (total - p - 1) / (p * (total - 2))
  f_statistic <- scale_factor * t_squared
  df <- c(df1 = p, df2 = total - p - 1)
  p_value <- stats::pf(f_statistic, df[["df1"]], df[["df2"]],
                       lower.tail = FALSE)

  variable_names <- .hotelling_variable_names(x)
  mean_x <- as.numeric(fit$mean_x)
  mean_y <- as.numeric(fit$mean_y)
  difference <- as.numeric(fit$difference)
  names(mean_x) <- names(mean_y) <- names(difference) <- variable_names
  null_difference <- stats::setNames(rep.int(0, p), variable_names)
  covariance <- as.matrix(fit$covariance)
  dimnames(covariance) <- list(variable_names, variable_names)

  .new_hd_location_test(
    statistic = c(F = f_statistic),
    parameter = df,
    p.value = p_value,
    method = "Two-sample Hotelling's T-squared test (equal covariance)",
    data.name = paste(x_name, "and", y_name),
    raw.statistic = c(T2 = t_squared),
    estimate = difference,
    null.value = null_difference,
    null.distribution = list(
      family = "F",
      parameters = df,
      exact = TRUE,
      assumptions = c(
        "independent multivariate normal samples",
        "common positive definite covariance"
      )
    ),
    variance = covariance,
    components = list(
      mean.x = mean_x,
      mean.y = mean_y,
      difference = difference,
      pooled.covariance = covariance,
      f.scale = scale_factor
    ),
    diagnostics = list(
      rank = as.integer(fit$rank),
      rcond = as.numeric(fit$rcond),
      solver = as.character(fit$solver),
      tolerance = tol,
      common.covariance = TRUE,
      exact.calibration = TRUE
    ),
    n = c(x = n1, y = n2),
    p = p,
    call = call
  )
}
