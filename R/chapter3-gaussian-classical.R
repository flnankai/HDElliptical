.ch3g_new_test <- function(statistic, p.value, method, data.name,
                           parameter = numeric(), alternative = "greater",
                           raw.statistic = NULL, estimate = NULL,
                           null.value = NULL, null.distribution = list(),
                           components = list(), diagnostics = list(), n, p,
                           call) {
  statistic.names <- names(statistic)
  statistic <- as.numeric(statistic)
  if (length(statistic) != 1L || !is.finite(statistic)) {
    stop("The calibrated test statistic is not finite.", call. = FALSE)
  }
  names(statistic) <- if (length(statistic.names) == 1L &&
    nzchar(statistic.names)) statistic.names else "statistic"
  p.value <- as.numeric(p.value)
  if (length(p.value) != 1L || !is.finite(p.value) ||
      p.value < 0 || p.value > 1) {
    stop("The computed p-value is not in [0, 1].", call. = FALSE)
  }
  parameter.names <- names(parameter)
  parameter <- as.numeric(parameter)
  if (length(parameter)) {
    names(parameter) <- parameter.names
  } else {
    parameter <- NULL
  }
  answer <- list(
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
    components = components,
    diagnostics = diagnostics,
    n = n,
    p = p,
    call = call
  )
  class(answer) <- c("hd_covariance_test", "htest")
  answer
}


.ch3g_validate_center <- function(center) {
  if (!is.logical(center) || length(center) != 1L || is.na(center)) {
    stop("`center` must be TRUE or FALSE.", call. = FALSE)
  }
  center
}


.ch3g_centered_scaled <- function(x, center) {
  center <- .ch3g_validate_center(center)
  scale <- max(abs(x))
  if (!is.finite(scale) || scale <= 0) {
    stop("`x` must contain nonconstant, nonzero finite data.", call. = FALSE)
  }
  scaled <- x / scale
  location <- if (center) colMeans(scaled) else numeric(ncol(x))
  residual <- sweep(scaled, 2L, location, "-")
  list(
    residual = residual,
    scale = scale,
    location.scaled = location,
    effective.df = nrow(x) - as.integer(center)
  )
}


.ch3g_spd <- function(x, name, p = NULL) {
  x <- as.matrix(x)
  storage.mode(x) <- "double"
  if (length(dim(x)) != 2L || nrow(x) != ncol(x) ||
      (!is.null(p) && nrow(x) != p) || anyNA(x) || any(!is.finite(x))) {
    stop(sprintf("`%s` must be a finite square matrix of the required dimension.",
                 name), call. = FALSE)
  }
  tolerance <- 100 * .Machine$double.eps * max(1, max(abs(x)))
  if (max(abs(x - t(x))) > tolerance) {
    stop(sprintf("`%s` must be symmetric.", name), call. = FALSE)
  }
  x <- (x + t(x)) / 2
  factor <- tryCatch(chol(x), error = function(e) NULL)
  if (is.null(factor)) {
    stop(sprintf("`%s` must be strictly positive definite.", name),
         call. = FALSE)
  }
  list(matrix = x, chol = factor)
}


#' Classical Gaussian likelihood-ratio test for a covariance matrix
#'
#' Tests \eqn{H_0:\Sigma=\Sigma_0} with the Gaussian likelihood-ratio
#' statistic
#' \deqn{n\{\operatorname{tr}(\Sigma_0^{-1}S_n)
#' -\log|\Sigma_0^{-1}S_n|-p\},}
#' where \eqn{S_n=n^{-1}\sum_i(X_i-\bar X)(X_i-\bar X)'} when
#' `center = TRUE`. Its fixed-dimensional reference distribution is
#' chi-squared with \eqn{p(p+1)/2} degrees of freedom. The sample covariance
#' must be positive definite; no pseudo-determinant or ridge is substituted.
#'
#' @param x Numeric matrix with observations in rows.
#' @param sigma0 Finite symmetric positive-definite null covariance matrix.
#' @param center Whether to estimate and remove the mean. `FALSE` implements
#'   the known-zero-mean model.
#'
#' @return An object of class `c("hd_covariance_test", "htest")`.
#'
#' @references
#' Anderson, T. W. (2003). *An Introduction to Multivariate Statistical
#' Analysis*, 3rd ed. Wiley.
#'
#' @examples
#' x <- matrix(rnorm(80), 20, 4)
#' gaussian_covariance_lrt(x, diag(4))
#'
#' @export
gaussian_covariance_lrt <- function(x, sigma0, center = TRUE) {
  call <- match.call()
  data.name <- deparse1(substitute(x))
  x <- .as_data_matrix(x, "x", min_rows = 2L)
  center <- .ch3g_validate_center(center)
  n <- nrow(x)
  p <- ncol(x)
  if (n - as.integer(center) < p) {
    stop("The sample covariance must have at least `p` residual degrees of freedom.",
         call. = FALSE)
  }
  sigma0 <- .ch3g_spd(sigma0, "sigma0", p)$matrix

  common.scale <- max(max(abs(x)), sqrt(max(abs(sigma0))))
  if (!is.finite(common.scale) || common.scale <= 0) {
    stop("A finite positive common scale could not be constructed.",
         call. = FALSE)
  }
  xs <- x / common.scale
  sigma0s <- sigma0 / common.scale / common.scale
  sigma.fit <- .ch3g_spd(sigma0s, "sigma0", p)
  residual <- if (center) sweep(xs, 2L, colMeans(xs), "-") else xs
  whitened <- t(backsolve(sigma.fit$chol, t(residual), transpose = TRUE))
  relative <- crossprod(whitened) / n
  eigenvalues <- eigen(relative, symmetric = TRUE, only.values = TRUE)$values
  if (any(!is.finite(eigenvalues)) || any(eigenvalues <= 0)) {
    stop("The sample covariance relative to `sigma0` is not positive definite.",
         call. = FALSE)
  }
  trace.term <- sum(eigenvalues)
  logdet.term <- sum(log(eigenvalues))
  statistic <- n * (trace.term - logdet.term - p)
  df <- p * (p + 1) / 2

  .ch3g_new_test(
    statistic = c(X.squared = statistic),
    p.value = stats::pchisq(statistic, df = df, lower.tail = FALSE),
    method = "Classical Gaussian covariance likelihood-ratio test",
    data.name = data.name,
    parameter = c(df = df),
    raw.statistic = c(LRT = statistic),
    null.value = sigma0,
    null.distribution = list(
      family = "chi-squared", parameters = c(df = df), exact = FALSE,
      tail = "upper", assumptions = "Gaussian sampling with fixed p"
    ),
    components = list(
      relative.eigenvalues = eigenvalues,
      trace.term = trace.term,
      log.determinant.term = logdet.term
    ),
    diagnostics = list(
      center.estimated = center,
      covariance.divisor = n,
      internal.common.scale = common.scale,
      regularization = "none",
      determinant = "ordinary positive-definite determinant"
    ),
    n = n, p = p, call = call
  )
}


#' Mauchly--Box likelihood-ratio test of Gaussian sphericity
#'
#' With residual Wishart degrees of freedom \eqn{m}, this function computes
#' \deqn{V=|S|/(\operatorname{tr}(S)/p)^p}
#' and the Box--Bartlett statistic
#' \deqn{-\left\{m-\frac{2p^2+p+2}{6p}\right\}\log V.}
#' The denominator in the correction is \eqn{p}, not \eqn{p+1}. The latter
#' appears in the accompanying book draft and is a transcription error.
#'
#' @inheritParams gaussian_covariance_lrt
#' @return An `hd_covariance_test` object.
#' @references Mauchly, J. W. (1940). *Annals of Mathematical Statistics*,
#'   11, 204--209. Wang, Q. and Yao, J. (2013). *Electronic Journal of
#'   Statistics*, 7, 2164--2192.
#' @examples
#' set.seed(35)
#' x <- matrix(rnorm(60), nrow = 20, ncol = 3)
#' mauchly_sphericity_test(x)
#' @export
mauchly_sphericity_test <- function(x, center = TRUE) {
  call <- match.call()
  data.name <- deparse1(substitute(x))
  x <- .as_data_matrix(x, "x", min_rows = 3L)
  p <- ncol(x)
  if (p < 2L) stop("Mauchly's test requires at least two variables.", call. = FALSE)
  work <- .ch3g_centered_scaled(x, center)
  m <- work$effective.df
  if (m < p) {
    stop("Mauchly's determinant requires residual degrees of freedom at least p.",
         call. = FALSE)
  }
  S <- crossprod(work$residual) / m
  values <- eigen(S, symmetric = TRUE, only.values = TRUE)$values
  if (any(values <= 0) || any(!is.finite(values))) {
    stop("The residual sample covariance is not positive definite.", call. = FALSE)
  }
  log.V <- sum(log(values)) - p * log(mean(values))
  correction.constant <- (2 * p^2 + p + 2) / (6 * p)
  rho <- 1 - correction.constant / m
  if (rho <= 0) {
    stop("The Box--Bartlett correction is non-positive at this sample size.",
         call. = FALSE)
  }
  statistic <- -rho * m * log.V
  df <- (p - 1) * (p + 2) / 2
  .ch3g_new_test(
    statistic = c(X.squared = statistic),
    p.value = stats::pchisq(statistic, df, lower.tail = FALSE),
    method = "Mauchly-Box Gaussian sphericity test",
    data.name = data.name,
    parameter = c(df = df),
    raw.statistic = c(log.V = log.V, V = exp(log.V)),
    null.distribution = list(
      family = "chi-squared", parameters = c(df = df), exact = FALSE,
      tail = "upper", assumptions = "Gaussian sampling with fixed p"
    ),
    components = list(
      log.V = log.V, V = exp(log.V), rho = rho,
      correction.constant = correction.constant
    ),
    diagnostics = list(
      residual.df = m, center.estimated = center,
      correction.denominator = "6 * p * residual.df",
      internal.data.scale = work$scale,
      regularization = "none"
    ),
    n = nrow(x), p = p, call = call
  )
}


#' John's classical trace test of Gaussian sphericity
#'
#' Computes \eqn{U=p\operatorname{tr}(S^2)/\operatorname{tr}^2(S)-1}
#' and uses \eqn{mpU/2}, with residual degrees of freedom \eqn{m}, as a
#' chi-squared statistic with \eqn{(p-1)(p+2)/2} degrees of freedom. The
#' multiplier is \eqn{p}; the \eqn{p+2} multiplier in the book draft belongs
#' to the spatial-sign analogue, not John's covariance statistic.
#'
#' @inheritParams gaussian_covariance_lrt
#' @return An `hd_covariance_test` object.
#' @references
#' John, S. (1971). Some optimal multivariate tests.
#' *Biometrika*, 58, 123--127.
#' @examples
#' set.seed(36)
#' x <- matrix(rnorm(60), nrow = 20, ncol = 3)
#' john_sphericity_test(x)
#' @export
john_sphericity_test <- function(x, center = TRUE) {
  call <- match.call()
  data.name <- deparse1(substitute(x))
  x <- .as_data_matrix(x, "x", min_rows = 3L)
  p <- ncol(x)
  if (p < 2L) stop("John's test requires at least two variables.", call. = FALSE)
  work <- .ch3g_centered_scaled(x, center)
  m <- work$effective.df
  S <- crossprod(work$residual) / m
  trace.S <- sum(diag(S))
  if (!is.finite(trace.S) || trace.S <= 0) {
    stop("The residual covariance has non-positive trace.", call. = FALSE)
  }
  U <- p * sum(S * S) / trace.S^2 - 1
  statistic <- m * p * U / 2
  df <- (p - 1) * (p + 2) / 2
  .ch3g_new_test(
    statistic = c(X.squared = statistic),
    p.value = stats::pchisq(statistic, df, lower.tail = FALSE),
    method = "John's classical Gaussian sphericity test",
    data.name = data.name,
    parameter = c(df = df),
    raw.statistic = c(U = U),
    null.distribution = list(
      family = "chi-squared", parameters = c(df = df), exact = FALSE,
      tail = "upper", assumptions = "Gaussian sampling with fixed p"
    ),
    components = list(U = U, trace.S = trace.S, trace.S2 = sum(S * S)),
    diagnostics = list(
      residual.df = m, center.estimated = center,
      multiplier = "residual.df * p / 2",
      internal.data.scale = work$scale,
      regularization = "none"
    ),
    n = nrow(x), p = p, call = call
  )
}


#' Nagao's classical Gaussian identity-covariance test
#'
#' Tests \eqn{H_0:\Sigma=I_p} with
#' \deqn{\frac{m}{2}\operatorname{tr}(S-I_p)^2,}
#' using the covariance with divisor equal to the residual Wishart degrees of
#' freedom \eqn{m}. Its fixed-dimensional reference has \eqn{p(p+1)/2}
#' degrees of freedom. The factor one-half is missing from the book draft.
#'
#' @inheritParams gaussian_covariance_lrt
#' @return An `hd_covariance_test` object.
#' @references
#' Nagao, H. (1973). On some test criteria for covariance matrix.
#' *Annals of Statistics*, 1, 700--709.
#' @examples
#' set.seed(37)
#' x <- matrix(rnorm(60), nrow = 20, ncol = 3)
#' nagao_identity_test(x)
#' @export
nagao_identity_test <- function(x, center = TRUE) {
  call <- match.call()
  data.name <- deparse1(substitute(x))
  x <- .as_data_matrix(x, "x", min_rows = 3L)
  center <- .ch3g_validate_center(center)
  n <- nrow(x)
  p <- ncol(x)
  m <- n - as.integer(center)
  residual <- if (center) sweep(x, 2L, colMeans(x), "-") else x
  S <- crossprod(residual) / m
  departure <- S - diag(p)
  distance <- sum(departure * departure)
  statistic <- m * distance / 2
  df <- p * (p + 1) / 2
  if (!is.finite(statistic)) {
    stop("Nagao's statistic overflowed; rescale the measurement units and null covariance together.",
         call. = FALSE)
  }
  .ch3g_new_test(
    statistic = c(X.squared = statistic),
    p.value = stats::pchisq(statistic, df, lower.tail = FALSE),
    method = "Nagao's classical Gaussian identity-covariance test",
    data.name = data.name,
    parameter = c(df = df),
    raw.statistic = c(V.N = distance / p),
    null.value = diag(p),
    null.distribution = list(
      family = "chi-squared", parameters = c(df = df), exact = FALSE,
      tail = "upper", assumptions = "Gaussian sampling with fixed p"
    ),
    components = list(V.N = distance / p, frobenius.distance.squared = distance),
    diagnostics = list(
      residual.df = m, center.estimated = center,
      covariance.divisor = m,
      multiplier = "residual.df / 2",
      regularization = "none"
    ),
    n = n, p = p, call = call
  )
}


#' Hallin--Paindaveine signed-rank test for elliptical shape
#'
#' For null shape \eqn{V_0}, the observations are sphericized, their radii
#' are ranked, and their directions are weighted by a score `K`. The statistic
#' is Hallin and Paindaveine's equation (4.3),
#' \deqn{\frac{p(p+2)}{2nE\{K^2(U)\}}
#' \sum_{i,j}K(R_i/(n+1))K(R_j/(n+1))
#' \{(U_i'U_j)^2-1/p\}.}
#' It has an asymptotic chi-squared reference with
#' \eqn{p(p+1)/2-1} degrees of freedom. Scores are `"sign"` (constant),
#' `"wilcoxon"` (power one), `"spearman"` (power two), and `"vdw"`
#' (chi-squared normal scores). Non-sign scores require distinct radii, as in
#' the paper's continuous elliptical model. No random tie breaking is used.
#'
#' @param x Numeric matrix with observations in rows.
#' @param shape0 Null shape matrix. Multiplying it by a positive scalar has no
#'   effect.
#' @param center A supplied finite center vector. `NULL` estimates the spatial
#'   median.
#' @param score Signed-rank score family.
#' @param tol,max_iter,zero_tol Controls passed to \code{spatial_median()} when
#'   the center is estimated.
#'
#' @return An `hd_covariance_test` object. `components` retains scores, ranks,
#'   signs, radii on the log scale, and the weighted sign scatter matrix.
#' @references Hallin, M. and Paindaveine, D. (2006). *Annals of Statistics*,
#'   34, 2707--2756. \doi{10.1214/009053606000000731}.
#' @examples
#' set.seed(38)
#' x <- matrix(rnorm(30), nrow = 15, ncol = 2)
#' hallin_paindaveine_shape_test(
#'   x, diag(2), center = c(0, 0), score = "sign"
#' )
#' @export
hallin_paindaveine_shape_test <- function(
    x, shape0 = NULL, center = NULL,
    score = c("sign", "wilcoxon", "spearman", "vdw"),
    tol = 1e-8, max_iter = 500L, zero_tol = 0) {
  call <- match.call()
  data.name <- deparse1(substitute(x))
  x <- .as_data_matrix(x, "x", min_rows = 3L)
  n <- nrow(x)
  p <- ncol(x)
  if (p < 2L) stop("Shape testing requires at least two variables.", call. = FALSE)
  score <- match.arg(score)
  if (is.null(shape0)) shape0 <- diag(p)
  shape.fit <- .ch3g_spd(shape0, "shape0", p)
  if (is.null(center)) {
    fitted.center <- spatial_median(
      x, tol = tol, max_iter = max_iter, zero_tol = zero_tol
    )
    center.value <- as.numeric(fitted.center)
    center.source <- "spatial median"
    center.diagnostics <- unclass(fitted.center)
  } else {
    center.value <- .as_location(center, p, "center")
    center.source <- "supplied"
    center.diagnostics <- NULL
  }
  residual <- sweep(x, 2L, center.value, "-")
  sphericized <- t(backsolve(shape.fit$chol, t(residual), transpose = TRUE))
  geometry <- cpp_ch3_gaussian_sign_geometry(sphericized)
  log.radii <- as.numeric(geometry$log_radii)
  if (score != "sign" && anyDuplicated(log.radii)) {
    stop("Non-sign Hallin--Paindaveine scores require distinct radii; ties are not broken silently.",
         call. = FALSE)
  }
  ranks <- if (score == "sign") rep.int(NA_integer_, n) else
    rank(log.radii, ties.method = "first")
  u <- if (score == "sign") rep(NA_real_, n) else ranks / (n + 1)
  scores <- switch(
    score,
    sign = rep(1, n),
    wilcoxon = u,
    spearman = u^2,
    vdw = stats::qchisq(u, df = p)
  )
  score.second <- switch(
    score,
    sign = 1,
    wilcoxon = 1 / 3,
    spearman = 1 / 5,
    vdw = p * (p + 2)
  )
  out <- cpp_ch3_gaussian_rank_shape(
    geometry$directions, scores, score.second
  )
  statistic <- as.numeric(out$statistic)
  numerical.tolerance <- 1000 * .Machine$double.eps *
    max(1, abs(as.numeric(out$trace_square)))
  if (statistic < -numerical.tolerance) {
    stop("The signed-rank quadratic contrast is materially negative.", call. = FALSE)
  }
  if (statistic < 0) statistic <- 0
  df <- p * (p + 1) / 2 - 1

  .ch3g_new_test(
    statistic = c(X.squared = statistic),
    p.value = stats::pchisq(statistic, df, lower.tail = FALSE),
    method = paste0("Hallin-Paindaveine ", score,
                    "-score signed-rank shape test"),
    data.name = data.name,
    parameter = c(df = df),
    raw.statistic = c(Q.K = statistic),
    null.value = shape.fit$matrix,
    null.distribution = list(
      family = "chi-squared", parameters = c(df = df), exact = FALSE,
      tail = "upper",
      assumptions = paste(
        "continuous elliptical sampling; fixed p;",
        "root-n consistent center when estimated"
      )
    ),
    components = list(
      score = score,
      score.second.moment = score.second,
      scores = scores,
      ranks = ranks,
      log.radii = log.radii,
      directions = geometry$directions,
      weighted.sign.scatter = out$scatter,
      quadratic.contrast = as.numeric(out$contrast)
    ),
    diagnostics = list(
      center = center.value,
      center.source = center.source,
      center.fit = center.diagnostics,
      shape.normalization = "irrelevant positive scalar",
      tied.radii = anyDuplicated(log.radii) > 0L,
      regularization = "none"
    ),
    n = n, p = p, call = call
  )
}


#' Spatial-sign sphericity test
#'
#' This is the constant-score special case of
#' \code{hallin_paindaveine_shape_test()}. It reports
#' \eqn{Q_S=p\operatorname{tr}(\Omega-I_p/p)^2} and calibrates
#' \eqn{n(p+2)Q_S/2} by chi-squared with
#' \eqn{(p-1)(p+2)/2} degrees of freedom.
#'
#' @inheritParams hallin_paindaveine_shape_test
#' @return An `hd_covariance_test` object.
#' @references
#' Hallin, M. and Paindaveine, D. (2006). *Annals of Statistics*,
#' 34, 2707--2756. \doi{10.1214/009053606000000731}.
#' @examples
#' set.seed(39)
#' x <- matrix(rnorm(30), nrow = 15, ncol = 2)
#' spatial_sign_sphericity_test(x, center = c(0, 0))
#' @export
spatial_sign_sphericity_test <- function(
    x, center = NULL, tol = 1e-8, max_iter = 500L, zero_tol = 0) {
  call <- match.call()
  answer <- hallin_paindaveine_shape_test(
    x = x, shape0 = diag(ncol(as.matrix(x))), center = center,
    score = "sign", tol = tol, max_iter = max_iter, zero_tol = zero_tol
  )
  omega <- answer$components$weighted.sign.scatter
  p <- answer$p
  n <- answer$n
  Q.S <- p * (sum(omega * omega) - 1 / p)
  answer$raw.statistic <- c(Q.S = Q.S)
  answer$components$Q.S <- Q.S
  answer$method <- "Fixed-p spatial-sign sphericity test"
  answer$call <- call
  answer
}
