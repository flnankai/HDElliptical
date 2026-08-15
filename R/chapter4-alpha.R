.ch4_alpha_matrix <- function(x, name, min_rows = 2L, min_cols = 1L) {
  if (is.data.frame(x)) {
    if (!all(vapply(x, is.numeric, logical(1)))) {
      stop(sprintf("%s must contain only numeric columns.", name),
           call. = FALSE)
    }
    x <- as.matrix(x)
  }
  if (!is.matrix(x) || !is.numeric(x) || length(dim(x)) != 2L) {
    stop(sprintf("%s must be a numeric matrix or numeric data frame.", name),
         call. = FALSE)
  }
  storage.mode(x) <- "double"
  if (nrow(x) < min_rows || ncol(x) < min_cols) {
    stop(sprintf(
      "%s must have at least %d rows and %d columns.",
      name, min_rows, min_cols
    ), call. = FALSE)
  }
  if (anyNA(x) || any(!is.finite(x))) {
    stop(sprintf("%s must contain only finite values.", name),
         call. = FALSE)
  }
  x
}

.ch4_alpha_factors <- function(factors, n) {
  if (is.null(factors)) {
    return(matrix(numeric(), nrow = n, ncol = 0L))
  }
  if (is.vector(factors) && is.numeric(factors) && is.null(dim(factors))) {
    factors <- matrix(factors, ncol = 1L)
  }
  factors <- .ch4_alpha_matrix(
    factors, "factors", min_rows = n, min_cols = 1L
  )
  if (nrow(factors) != n) {
    stop("returns and factors must have the same number of rows.",
         call. = FALSE)
  }
  factors
}

.ch4_alpha_column_scale <- function(x, name) {
  if (!ncol(x)) {
    return(list(data = x, scale = numeric()))
  }
  scale <- apply(abs(x), 2L, max)
  if (any(!is.finite(scale)) || any(scale <= 0)) {
    stop(sprintf("Every column of %s must have positive finite magnitude.",
                 name), call. = FALSE)
  }
  list(data = sweep(x, 2L, scale, "/"), scale = as.numeric(scale))
}

.ch4_alpha_has_full_column_rank <- function(x) {
  if (!ncol(x)) {
    return(TRUE)
  }
  rank <- qr(x, tol = sqrt(.Machine$double.eps), LAPACK = FALSE)$rank
  isTRUE(rank == ncol(x))
}

.ch4_alpha_prepare <- function(returns, factors, min_assets = 1L) {
  returns <- .ch4_alpha_matrix(
    returns, "returns", min_rows = 2L, min_cols = min_assets
  )
  factors <- .ch4_alpha_factors(factors, nrow(returns))
  y <- .ch4_alpha_column_scale(returns, "returns")
  f <- .ch4_alpha_column_scale(factors, "factors")
  if (nrow(returns) <= ncol(factors) + 1L) {
    stop(
      "The intercept-plus-factor regression needs positive residual degrees of freedom.",
      call. = FALSE
    )
  }
  design <- cbind(intercept = 1, f$data)
  design.cross <- crossprod(design)
  design.root <- if (.ch4_alpha_has_full_column_rank(design)) {
    tryCatch(chol(design.cross), error = function(e) NULL)
  } else {
    NULL
  }
  if (is.null(design.root)) {
    stop(
      "The intercept-plus-factor design is rank deficient; no generalized inverse is used.",
      call. = FALSE
    )
  }
  list(
    returns = returns,
    factors = factors,
    y = y$data,
    f = f$data,
    y.scale = y$scale,
    f.scale = f$scale,
    T = nrow(returns),
    N = ncol(returns),
    K = ncol(factors),
    asset.names = if (is.null(colnames(returns))) {
      paste0("asset", seq_len(ncol(returns)))
    } else {
      colnames(returns)
    },
    factor.names = if (is.null(colnames(factors))) {
      paste0("factor", seq_len(ncol(factors)))
    } else {
      colnames(factors)
    }
  )
}

.ch4_alpha_ols <- function(prepared) {
  fit <- cpp_ch4_alpha_ols(prepared$y, prepared$f)
  fit$alpha <- as.numeric(fit$alpha)
  fit$t_squared <- as.numeric(fit$t_squared)
  fit$h <- as.numeric(fit$h)
  names(fit$alpha) <- prepared$asset.names
  names(fit$t_squared) <- prepared$asset.names
  colnames(fit$residuals) <- prepared$asset.names
  dimnames(fit$residual_correlation) <- list(
    prepared$asset.names, prepared$asset.names
  )
  alpha.input <- fit$alpha * prepared$y.scale
  names(alpha.input) <- prepared$asset.names
  slopes.scaled <- if (prepared$K) {
    fit$coefficients[-1L, , drop = FALSE]
  } else {
    matrix(numeric(), nrow = 0L, ncol = prepared$N)
  }
  slopes.input <- slopes.scaled
  if (prepared$K) {
    slopes.input <- sweep(slopes.input, 1L, prepared$f.scale, "/")
    slopes.input <- sweep(slopes.input, 2L, prepared$y.scale, "*")
    dimnames(slopes.input) <- list(
      prepared$factor.names, prepared$asset.names
    )
  }
  list(
    native = fit,
    alpha = alpha.input,
    slopes = slopes.input,
    slopes.scaled = slopes.scaled,
    residuals.scaled = fit$residuals,
    residuals = sweep(fit$residuals, 2L, prepared$y.scale, "*")
  )
}

.ch4_alpha_probability <- function(x, name) {
  x <- as.numeric(x)
  if (length(x) != 1L || is.na(x) || !is.finite(x) ||
      x <= 0 || x >= 1) {
    stop(sprintf("%s must be one finite number strictly between zero and one.",
                 name), call. = FALSE)
  }
  x
}

.ch4_alpha_positive <- function(x, name, integer = FALSE) {
  raw <- x
  x <- as.numeric(x)
  if (length(x) != 1L || is.na(x) || !is.finite(x) || x <= 0 ||
      (integer && (x != floor(x) || x > .Machine$integer.max))) {
    stop(sprintf("%s must be one finite positive%s number.",
                 name, if (integer) " integer" else ""), call. = FALSE)
  }
  if (integer) as.integer(x) else as.numeric(raw)
}

.ch4_alpha_nonnegative <- function(x, name) {
  x <- as.numeric(x)
  if (length(x) != 1L || is.na(x) || !is.finite(x) || x < 0) {
    stop(sprintf("%s must be one finite non-negative number.", name),
         call. = FALSE)
  }
  x
}

.ch4_alpha_flag <- function(x, name) {
  if (!is.logical(x) || length(x) != 1L || is.na(x)) {
    stop(sprintf("%s must be TRUE or FALSE.", name), call. = FALSE)
  }
  x
}

.ch4_alpha_gumbel_tail <- function(centered) {
  centered <- as.numeric(centered)
  if (length(centered) != 1L || is.na(centered) ||
      !is.finite(centered)) {
    stop("The centered maximum statistic is not finite.", call. = FALSE)
  }
  log.intensity <- -0.5 * log(pi) - 0.5 * centered
  if (log.intensity > log(.Machine$double.xmax)) {
    return(list(p.value = 1, log.p.value = 0,
                log.cdf = -Inf, log.intensity = log.intensity))
  }
  intensity <- exp(log.intensity)
  if (intensity == 0) {
    return(list(p.value = 0, log.p.value = log.intensity,
                log.cdf = 0, log.intensity = log.intensity))
  }
  log.p <- if (intensity <= log(2)) {
    log(-expm1(-intensity))
  } else {
    log1p(-exp(-intensity))
  }
  list(
    p.value = exp(log.p),
    log.p.value = log.p,
    log.cdf = -intensity,
    log.intensity = log.intensity
  )
}

.ch4_alpha_htest <- function(statistic, p.value, method, data.name,
                              estimate, raw.statistic = NULL,
                              parameter = NULL, components = list(),
                              diagnostics = list(), call = NULL) {
  answer <- list(
    statistic = statistic,
    parameter = parameter,
    p.value = as.numeric(p.value),
    method = method,
    data.name = data.name,
    alternative = "two.sided",
    estimate = estimate,
    null.value = stats::setNames(0, "alpha"),
    raw.statistic = raw.statistic,
    components = components,
    diagnostics = diagnostics,
    call = call
  )
  answer <- answer[!vapply(answer, is.null, logical(1))]
  class(answer) <- c("hd_alpha_test", "htest")
  answer
}


#' Gibbons--Ross--Shanken exact alpha test
#'
#' Tests whether all intercepts are zero in an unconditional linear
#' factor-pricing model. Rows of returns and factors are observations and
#' columns are assets and factors. The joint regression contains both an
#' intercept and all factor columns. With unrestricted OLS residual matrix E,
#' Vhat = E' E / T, and h the residual from regressing the all-ones vector on
#' the factors, the implemented primary statistic is
#' \deqn{[(T-N-K)/N][h'h/T]\hat\alpha'\widehat V^{-1}\hat\alpha,}
#' with the exact F distribution having N and T-N-K degrees of freedom under
#' the Gaussian GRS assumptions.
#'
#' This differs from two finite-sample expressions in the book draft. Its
#' displayed slope estimator omits the intercept although its residual formula
#' includes one. It also divides the residual covariance by T-K-1 while
#' retaining the primary statistic scaling for a divisor-T covariance, missing
#' the compensating factor T/(T-K-1). This function follows GRS exactly and
#' reports both divisors for audit.
#'
#' @param returns Finite observation-by-asset numeric matrix or data frame.
#' @param factors NULL, a finite length-T numeric vector, or a finite
#'   observation-by-factor numeric matrix or data frame.
#' @return An object of class c("hd_alpha_test", "htest") with the exact
#'   statistic, p-value, fitted intercepts, and OLS diagnostics.
#' @references
#' Gibbons, M. R., Ross, S. A., and Shanken, J. (1989). A test of the
#' efficiency of a given portfolio. Econometrica, 57, 1121--1152.
#' \doi{10.2307/1913625}.
#' @examples
#' f <- cbind(seq(-1, 1, length.out = 10))
#' y <- cbind(0.2 + f[, 1] + sin(seq_len(10)),
#'            -0.1 - 0.5 * f[, 1] + cos(seq_len(10)))
#' grs_alpha_test(y, f)
#' @export
grs_alpha_test <- function(returns, factors = NULL) {
  call <- match.call()
  data.name <- deparse1(substitute(returns))
  prepared <- .ch4_alpha_prepare(returns, factors)
  if (prepared$T <= prepared$N + prepared$K) {
    stop("The exact GRS test requires T > N + K.", call. = FALSE)
  }
  ols <- .ch4_alpha_ols(prepared)
  scatter <- crossprod(ols$residuals.scaled) / prepared$T
  root <- if (.ch4_alpha_has_full_column_rank(ols$residuals.scaled)) {
    tryCatch(chol(scatter), error = function(e) NULL)
  } else {
    NULL
  }
  if (is.null(root)) {
    stop(
      "The divisor-T OLS residual covariance must be positive definite; no generalized inverse or ridge is used.",
      call. = FALSE
    )
  }
  whitened <- backsolve(root, ols$native$alpha, transpose = TRUE)
  quadratic <- sum(whitened^2)
  statistic <- (prepared$T - prepared$N - prepared$K) /
    prepared$N * ols$native$h2 / prepared$T * quadratic
  if (!is.finite(statistic) || statistic < 0) {
    stop("The GRS quadratic form is invalid.", call. = FALSE)
  }
  denominator.df <- prepared$T - prepared$N - prepared$K
  p.value <- stats::pf(
    statistic, prepared$N, denominator.df, lower.tail = FALSE
  )

  .ch4_alpha_htest(
    statistic = stats::setNames(statistic, "GRS F"),
    p.value = p.value,
    method = "Gibbons-Ross-Shanken exact alpha test",
    data.name = data.name,
    estimate = ols$alpha,
    raw.statistic = stats::setNames(quadratic, "alpha quadratic"),
    parameter = c(df1 = prepared$N, df2 = denominator.df),
    components = list(
      alpha = ols$alpha,
      slopes = ols$slopes,
      residuals = ols$residuals,
      residual.covariance.scaled.divisor.T = scatter,
      residual.covariance.scaled.divisor.v =
        crossprod(ols$residuals.scaled) / ols$native$residual_df,
      asset.column.scales = prepared$y.scale,
      factor.column.scales = prepared$f.scale,
      h = as.numeric(ols$native$h),
      h2 = unname(ols$native$h2),
      quadratic.scaled = quadratic,
      T = prepared$T,
      N = prepared$N,
      K = prepared$K
    ),
    diagnostics = list(
      calibration = "exact upper-tail F under primary Gaussian assumptions",
      covariance.divisor = "T",
      book.slope.conflict = TRUE,
      book.covariance.divisor.conflict = TRUE,
      generalized.inverse = "none",
      ridge = "none"
    ),
    call = call
  )
}


#' Pesaran--Yamagata large-N sum alpha test
#'
#' Implements the feasible statistic of Pesaran and Yamagata. Let v=T-K-1
#' and t_i be the usual unrestricted OLS intercept t statistic. Residual
#' correlations are retained only when
#' \eqn{|\sqrt v\,\hat\rho_{ij}|>
#' \Phi^{-1}\{1-p_0/(2N^\delta)\}}. If \eqn{\widetilde\rho^2} is twice the
#' sum of squared retained upper-triangular correlations divided by
#' N(N-1), the primary statistic is
#' \deqn{\frac{N^{-1/2}\sum_i\{t_i^2-v/(v-2)\}}
#' {[v/(v-2)]\{2(v-1)(1+(N-1)\widetilde\rho^2)/(v-4)\}^{1/2}}.}
#'
#' The generic dense statistic in the book draft is not this feasible
#' finite-v test: it omits the t-square centering and the multiple-testing
#' threshold estimator. This function follows the primary paper.
#'
#' @inheritParams grs_alpha_test
#' @param p0 Finite thresholding probability in (0,1), default 0.1.
#' @param delta Finite positive threshold exponent, default 1.
#' @return An asymptotic upper-tail normal alpha-test object.
#' @references
#' Pesaran, M. H. and Yamagata, T. (2023). Testing for alpha in linear factor
#' pricing models with a large number of securities. Journal of Financial
#' Econometrics. \doi{10.1093/jjfinec/nbad002}.
#' @examples
#' f <- cbind(seq(-1, 1, length.out = 12))
#' y <- outer(seq_len(12), 1:3, function(i, j) sin(i + j)) +
#'   f %*% matrix(1:3, nrow = 1)
#' pesaran_yamagata_alpha_test(y, f)
#' @export
pesaran_yamagata_alpha_test <- function(returns, factors = NULL,
                                         p0 = 0.1, delta = 1) {
  call <- match.call()
  data.name <- deparse1(substitute(returns))
  p0 <- .ch4_alpha_probability(p0, "p0")
  delta <- .ch4_alpha_positive(delta, "delta")
  prepared <- .ch4_alpha_prepare(returns, factors, min_assets = 2L)
  ols <- .ch4_alpha_ols(prepared)
  v <- unname(ols$native$residual_df)
  if (v <= 4) {
    stop("The Pesaran-Yamagata finite-v standardization requires v > 4.",
         call. = FALSE)
  }
  threshold <- stats::qnorm(
    1 - p0 / (2 * prepared$N^delta)
  )
  upper <- upper.tri(ols$native$residual_correlation)
  rho <- ols$native$residual_correlation[upper]
  keep <- abs(sqrt(v) * rho) > threshold
  rho.tilde.sq <- 2 * sum(rho[keep]^2) /
    (prepared$N * (prepared$N - 1))
  center <- v / (v - 2)
  scale <- center * sqrt(
    2 * (v - 1) / (v - 4) *
      (1 + (prepared$N - 1) * rho.tilde.sq)
  )
  numerator <- sum(ols$native$t_squared - center) / sqrt(prepared$N)
  statistic <- numerator / scale
  if (!is.finite(statistic)) {
    stop("The Pesaran-Yamagata statistic is not finite.", call. = FALSE)
  }
  p.value <- stats::pnorm(statistic, lower.tail = FALSE)

  .ch4_alpha_htest(
    statistic = stats::setNames(statistic, "J.hat"),
    p.value = p.value,
    method = "Pesaran-Yamagata feasible large-N sum alpha test",
    data.name = data.name,
    estimate = ols$alpha,
    raw.statistic = stats::setNames(numerator, "centered t-square sum"),
    parameter = stats::setNames(v, "residual df"),
    components = list(
      alpha = ols$alpha,
      t.squared = ols$native$t_squared,
      residual.correlation = ols$native$residual_correlation,
      correlation.threshold = threshold,
      retained.upper.triangle = keep,
      retained.count = sum(keep),
      rho.tilde.squared = rho.tilde.sq,
      finite.v.center = center,
      finite.v.scale = scale,
      p0 = p0,
      delta = delta,
      T = prepared$T,
      N = prepared$N,
      K = prepared$K
    ),
    diagnostics = list(
      calibration = "asymptotic upper-tail standard normal",
      book.generic.sum.is.primary.PY = FALSE,
      covariance.repair = "none"
    ),
    call = call
  )
}


#' Feng--Lan--Liu--Ma high-dimensional maximum alpha test
#'
#' Computes the maximum squared unrestricted OLS intercept t statistic,
#' \eqn{M=\max_i t_i^2}, and its centered value
#' \eqn{M-2\log N+\log\log N}. The limiting cdf is
#' \eqn{\exp\{-\pi^{-1/2}\exp(-x/2)\}}. The upper tail is evaluated without
#' subtracting a cdf numerically.
#'
#' @inheritParams grs_alpha_test
#' @return An asymptotic Gumbel-calibrated upper-tail alpha-test object.
#' @references
#' Feng, L., Lan, W., Liu, B., and Ma, Y. (2022). High-dimensional test for
#' alpha in linear factor pricing models. Journal of Econometrics.
#' \doi{10.1016/j.jeconom.2021.07.011}.
#' @examples
#' f <- cbind(seq(-1, 1, length.out = 12))
#' y <- outer(seq_len(12), 1:3, function(i, j) cos(i + 2 * j))
#' feng_lan_liu_ma_alpha_max_test(y, f)
#' @export
feng_lan_liu_ma_alpha_max_test <- function(returns, factors = NULL) {
  call <- match.call()
  data.name <- deparse1(substitute(returns))
  prepared <- .ch4_alpha_prepare(returns, factors, min_assets = 2L)
  ols <- .ch4_alpha_ols(prepared)
  raw <- max(ols$native$t_squared)
  centered <- raw - 2 * log(prepared$N) + log(log(prepared$N))
  tail <- .ch4_alpha_gumbel_tail(centered)

  .ch4_alpha_htest(
    statistic = stats::setNames(centered, "Gumbel centered"),
    p.value = tail$p.value,
    method = "Feng-Lan-Liu-Ma high-dimensional maximum alpha test",
    data.name = data.name,
    estimate = ols$alpha,
    raw.statistic = stats::setNames(raw, "maximum t squared"),
    components = list(
      alpha = ols$alpha,
      t.squared = ols$native$t_squared,
      maximum.t.squared = raw,
      centered.maximum = centered,
      log.p.value = tail$log.p.value,
      log.gumbel.cdf = tail$log.cdf,
      T = prepared$T,
      N = prepared$N,
      K = prepared$K
    ),
    diagnostics = list(
      calibration = "upper tail of primary type-I extreme-value limit",
      gumbel.cdf = "exp{-pi^(-1/2) exp(-x/2)}",
      numerical.subtraction.of.cdf = FALSE
    ),
    call = call
  )
}


#' Feng--Lan--Liu--Ma adaptive Gaussian alpha test
#'
#' Combines the Pesaran--Yamagata sum p-value and the Feng--Lan--Liu--Ma max
#' p-value with the primary paper's two-test Bonferroni rule
#' \deqn{p_{\rm COM}=\min\{1,2\min(p_{\rm PY},p_{\rm MAX})\}.}
#' It deliberately does not implement the generic Cauchy benchmark introduced
#' in the book draft, which is not the combination proposed in the cited
#' Feng--Lan--Liu--Ma method.
#'
#' @inheritParams pesaran_yamagata_alpha_test
#' @return An alpha-test object containing both complete component tests.
#' @references
#' Feng, L., Lan, W., Liu, B., and Ma, Y. (2022). High-dimensional test for
#' alpha in linear factor pricing models. Journal of Econometrics.
#' \doi{10.1016/j.jeconom.2021.07.011}.
#' @examples
#' f <- cbind(seq(-1, 1, length.out = 12))
#' y <- outer(seq_len(12), 1:3, function(i, j) sin(i * j))
#' gaussian_alpha_combination_test(y, f)
#' @export
gaussian_alpha_combination_test <- function(returns, factors = NULL,
                                             p0 = 0.1, delta = 1) {
  call <- match.call()
  data.name <- deparse1(substitute(returns))
  sum.component <- pesaran_yamagata_alpha_test(
    returns, factors, p0 = p0, delta = delta
  )
  max.component <- feng_lan_liu_ma_alpha_max_test(returns, factors)
  component.p <- c(
    PY = sum.component$p.value,
    MAX = max.component$p.value
  )
  p.value <- min(1, 2 * min(component.p))

  .ch4_alpha_htest(
    statistic = stats::setNames(min(component.p), "minimum component p"),
    p.value = p.value,
    method = paste(
      "Feng-Lan-Liu-Ma adaptive alpha test",
      "(primary two-component Bonferroni combination)"
    ),
    data.name = data.name,
    estimate = sum.component$estimate,
    raw.statistic = component.p,
    components = list(
      p.values = component.p,
      sum.test = sum.component,
      max.test = max.component,
      multiplicity.factor = 2
    ),
    diagnostics = list(
      calibration = "primary Bonferroni-adjusted minimum p-value",
      book.generic.Cauchy.implemented = FALSE
    ),
    call = call
  )
}


.ch4_alpha_controls <- function(tol, max_iter, zero_tol) {
  list(
    tol = .ch4_alpha_positive(tol, "tol"),
    max_iter = .ch4_alpha_positive(max_iter, "max_iter", integer = TRUE),
    zero_tol = .ch4_alpha_nonnegative(zero_tol, "zero_tol")
  )
}

.ch4_alpha_seed <- function(seed) {
  if (is.null(seed)) {
    return(NULL)
  }
  seed <- as.numeric(seed)
  if (length(seed) != 1L || is.na(seed) || !is.finite(seed) ||
      seed < 0 || seed != floor(seed) ||
      seed > .Machine$integer.max) {
    stop(
      "seed must be NULL or one integer from zero through R's integer limit.",
      call. = FALSE
    )
  }
  as.integer(seed)
}

.ch4_alpha_with_seed <- function(seed, fun) {
  if (is.null(seed)) {
    return(fun())
  }
  had.state <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had.state) {
    old.state <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  }
  on.exit({
    if (had.state) {
      assign(".Random.seed", old.state, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv,
                      inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  set.seed(seed)
  fun()
}

.ch4_alpha_canonical_diagonal <- function(diagonal, column.scale) {
  log.diagonal <- log(as.numeric(diagonal)) + 2 * log(column.scale)
  log.diagonal <- log.diagonal - max(log.diagonal)
  list(value = exp(log.diagonal), log = log.diagonal)
}

.ch4_alpha_check_lfm_fit <- function(fit, context) {
  if (!isTRUE(fit$converged)) {
    stop(sprintf(
      paste0(
        "%s did not satisfy the diagonal spatial-sign equation within ",
        "%d updates (equation residual %.6g, relative update %.6g)."
      ),
      context, as.integer(fit$iterations),
      as.numeric(fit$equation_residual),
      as.numeric(fit$relative_update)
    ), call. = FALSE)
  }
  invisible(fit)
}

.ch4_alpha_lfm_core <- function(prepared, controls, compute_trace = TRUE) {
  fit <- cpp_ch4_lfm_spatial_sign_core(
    prepared$y, prepared$f, controls$tol, controls$max_iter,
    controls$zero_tol, compute_trace
  )
  .ch4_alpha_check_lfm_fit(fit, "The Liu-Feng-Ma diagonal-scale fit")
  fit
}

.ch4_alpha_truncated_cauchy <- function(p.values) {
  p.values <- as.numeric(p.values)
  if (!length(p.values) || anyNA(p.values) ||
      any(!is.finite(p.values)) || any(p.values < 0) ||
      any(p.values > 1)) {
    stop("Truncated-Cauchy component p-values must lie in [0,1].",
         call. = FALSE)
  }
  active <- p.values < 0.5
  terms <- numeric(length(p.values))
  terms[active] <- 0.5 / tan(pi * p.values[active])
  score <- sum(terms)
  if (is.nan(score) || score < 0) {
    stop("The truncated-Cauchy score is invalid.", call. = FALSE)
  }
  list(
    statistic = score,
    p.value = atan2(1, score) / pi,
    terms = terms,
    active = active
  )
}

.ch4_alpha_restricted_residuals <- function(prepared) {
  if (!prepared$K) {
    return(list(
      slopes = matrix(numeric(), nrow = 0L, ncol = prepared$N),
      residuals = prepared$y,
      fitted = matrix(0, nrow = prepared$T, ncol = prepared$N),
      h = rep(1, prepared$T),
      eta = 0,
      omega = 1
    ))
  }
  root <- tryCatch(chol(crossprod(prepared$f)),
                   error = function(e) NULL)
  if (is.null(root)) {
    stop(
      "The restricted factor cross-product is singular; no generalized inverse is used.",
      call. = FALSE
    )
  }
  rhs <- crossprod(prepared$f, prepared$y)
  slopes <- backsolve(root, forwardsolve(t(root), rhs))
  fitted <- prepared$f %*% slopes
  intercept.rhs <- crossprod(prepared$f, rep(1, prepared$T))
  intercept.coef <- backsolve(
    root, forwardsolve(t(root), intercept.rhs)
  )
  projection <- as.numeric(prepared$f %*% intercept.coef)
  h <- 1 - projection
  list(
    slopes = slopes,
    residuals = prepared$y - fitted,
    fitted = fitted,
    h = h,
    eta = sum(projection^2) / prepared$T,
    omega = sum(h^2) / prepared$T
  )
}


#' Liu--Feng--Ma robust spatial-sign sum alpha test
#'
#' Implements the primary heavy-tailed alpha test. Restricted factor residuals
#' are standardized by the full-sample diagonal spatial-sign fixed point. If
#' U_t denotes the resulting sign and h is the residualized intercept, the
#' quadratic form is
#' \deqn{Q=N(h'h)^{-1}\sum_{s\ne t}h_s h_t U_s'U_t.}
#' The variance component is the primary split leave-two-out estimator. After
#' deleting a pair, the remaining rows are split into chronological first and
#' second halves, separate restricted slopes are fitted, and the two deleted
#' residuals use the corresponding slopes. Its denominator is
#' \eqn{(h'h)(h'h-1)}, not \eqn{(h'h)^2}.
#'
#' The recommended feasible statistic is
#' \deqn{(Q-\widehat\delta_Q)/
#' \{2\widehat{\operatorname{tr}(R^2)}\}^{1/2}.}
#' By default, \eqn{\widehat\delta_Q} is the mean of bootstrap Q values formed
#' by multiplying every unrestricted OLS residual entry by an independent
#' asset-time Rademacher sign and adding the restricted fitted values. This is
#' the method's
#' inferential bias calibration, not a reproduction of its simulation study.
#' An explicit seed is locally isolated: the caller's random-number state is
#' restored even after an error. With a NULL seed the current stream advances.
#'
#' The book draft omits the primary wild-bootstrap bias term, uses
#' \eqn{(h'h)^2} in the trace denominator, and describes pair-deleted scale
#' fits rather than the paper's full-sample diagonal scale plus split
#' pair-deleted slope fits. These choices are not silently mixed here.
#' The book's weighted/INST construction splices this incompatible trace into
#' a paper whose complete feasible calibration is not given in the cited
#' material, so weighted/INST is intentionally review-only. Likewise, the
#' dependent and Lq paragraphs do not specify a complete feasible
#' centering, long-run variance, and calibration and are not reconstructed.
#'
#' @inheritParams grs_alpha_test
#' @param bias Bias calibration: primary "wild_bootstrap", literal zero
#'   "none" for formula auditing, or a user "supplied" value.
#' @param delta_q Finite supplied bias when bias is "supplied"; otherwise
#'   must be NULL.
#' @param bootstrap_reps Positive number of Rademacher replicates. The primary
#'   recommendation is 100.
#' @param seed NULL to use and advance the current random stream, or a
#'   non-negative integer for an isolated deterministic bootstrap.
#' @param keep_bootstrap Whether to retain all bootstrap Q values.
#' @param tol Positive diagonal fixed-point tolerance.
#' @param max_iter Positive integer update limit.
#' @param zero_tol Non-negative standardized zero-radius threshold; default
#'   zero uses the literal spatial-sign convention.
#' @return An asymptotic upper-tail normal alpha-test object with every
#'   quadratic, trace, bias, scale, convergence, and bootstrap diagnostic.
#' @references
#' Liu, B., Feng, L., and Ma, Y. (2023). High-dimensional alpha test of linear
#' factor pricing models with heavy-tailed distributions. Statistica Sinica,
#' 33, 1389--1410. \doi{10.5705/ss.202021.0134}.
#' @examples
#' f <- cbind(seq(-1, 1, length.out = 10))
#' y <- outer(seq_len(10), 1:3, function(i, j) sin(i + j / 2))
#' liu_feng_ma_spatial_sign_alpha_test(
#'   y, f, bootstrap_reps = 2, seed = 1
#' )
#' @export
liu_feng_ma_spatial_sign_alpha_test <- function(
    returns, factors = NULL,
    bias = c("wild_bootstrap", "none", "supplied"),
    delta_q = NULL, bootstrap_reps = 100L, seed = NULL,
    keep_bootstrap = FALSE, tol = 1e-7, max_iter = 500L,
    zero_tol = 0) {
  call <- match.call()
  data.name <- deparse1(substitute(returns))
  bias <- match.arg(bias)
  bootstrap_reps <- .ch4_alpha_positive(
    bootstrap_reps, "bootstrap_reps", integer = TRUE
  )
  seed <- .ch4_alpha_seed(seed)
  keep_bootstrap <- .ch4_alpha_flag(keep_bootstrap, "keep_bootstrap")
  controls <- .ch4_alpha_controls(tol, max_iter, zero_tol)
  if (bias == "supplied") {
    delta.q <- as.numeric(delta_q)
    if (length(delta.q) != 1L || is.na(delta.q) ||
        !is.finite(delta.q)) {
      stop("delta_q must be one finite number when bias is supplied.",
           call. = FALSE)
    }
  } else {
    if (!is.null(delta_q)) {
      stop("delta_q must be NULL unless bias is supplied.",
           call. = FALSE)
    }
    delta.q <- if (bias == "none") 0 else NA_real_
  }

  prepared <- .ch4_alpha_prepare(returns, factors, min_assets = 2L)
  if (prepared$T < 2L * prepared$K + 2L) {
    stop(
      "The primary split trace estimator requires T >= 2*K + 2.",
      call. = FALSE
    )
  }
  core <- .ch4_alpha_lfm_core(prepared, controls, compute_trace = TRUE)
  bootstrap.q <- NULL
  if (bias == "wild_bootstrap") {
    unrestricted <- .ch4_alpha_ols(prepared)
    restricted.fitted <- prepared$y - core$restricted_residuals
    bootstrap.q <- .ch4_alpha_with_seed(seed, function() {
      values <- numeric(bootstrap_reps)
      for (b in seq_len(bootstrap_reps)) {
        multipliers <- matrix(
          sample(c(-1, 1), prepared$T * prepared$N, replace = TRUE),
          nrow = prepared$T, ncol = prepared$N
        )
        y.star <- restricted.fitted +
          unrestricted$residuals.scaled * multipliers
        fit.star <- cpp_ch4_lfm_spatial_sign_core(
          y.star, prepared$f, controls$tol, controls$max_iter,
          controls$zero_tol, FALSE
        )
        .ch4_alpha_check_lfm_fit(
          fit.star, sprintf("Wild-bootstrap replicate %d", b)
        )
        values[b] <- fit.star$Q
      }
      values
    })
    delta.q <- mean(bootstrap.q)
  }
  statistic <- (core$Q - delta.q) / sqrt(2 * core$trace_R2)
  if (!is.finite(statistic)) {
    stop("The Liu-Feng-Ma standardized statistic is not finite.",
         call. = FALSE)
  }
  p.value <- stats::pnorm(statistic, lower.tail = FALSE)
  diagonal <- .ch4_alpha_canonical_diagonal(
    core$diagonal, prepared$y.scale
  )
  names(diagonal$value) <- names(diagonal$log) <- prepared$asset.names
  directions <- as.matrix(core$directions)
  dimnames(directions) <- list(rownames(prepared$returns),
                               prepared$asset.names)

  .ch4_alpha_htest(
    statistic = stats::setNames(statistic, "T.SS"),
    p.value = p.value,
    method = paste(
      "Liu-Feng-Ma robust spatial-sign sum alpha test",
      sprintf("(%s bias calibration)", bias)
    ),
    data.name = data.name,
    estimate = stats::setNames(
      .ch4_alpha_ols(prepared)$alpha, prepared$asset.names
    ),
    raw.statistic = c(Q = core$Q, delta.Q = delta.q),
    components = list(
      Q = unname(core$Q),
      delta.Q = delta.q,
      trace.R.squared = unname(core$trace_R2),
      trace.denominator = unname(core$trace_denominator),
      trace.split = core$trace_split,
      h = as.numeric(core$h),
      h2 = unname(core$h2),
      scale.diagonal = diagonal$value,
      log.scale.diagonal = diagonal$log,
      scale.diagonal.scaled.coordinates = as.numeric(core$diagonal),
      directions = directions,
      radii.scaled.coordinates = as.numeric(core$radii),
      restricted.slopes.scaled = core$restricted_slopes,
      restricted.residuals.scaled = core$restricted_residuals,
      bootstrap.Q = if (keep_bootstrap) bootstrap.q else NULL,
      bootstrap.reps = if (bias == "wild_bootstrap") {
        bootstrap_reps
      } else {
        0L
      },
      seed = seed,
      T = prepared$T,
      N = prepared$N,
      K = prepared$K
    ),
    diagnostics = list(
      calibration = "asymptotic upper-tail standard normal",
      bias.method = bias,
      rng.isolated = !is.null(seed) && bias == "wild_bootstrap",
      converged = isTRUE(core$converged),
      iterations = as.integer(core$iterations),
      equation.residual = unname(core$equation_residual),
      relative.update = unname(core$relative_update),
      zero.residuals = as.integer(core$zero_residuals),
      trace.denominator.primary = "(h'h)(h'h-1)",
      book.trace.denominator.conflict = TRUE,
      book.bias.omission = TRUE,
      book.leaveout.scale.conflict = TRUE,
      weighted.INST.status = "review-only: complete feasible calibration unavailable",
      dependent.Lq.status = "review-only: complete feasible calibration unavailable",
      ridge = "none",
      scale.floor = "none",
      generalized.inverse = "none"
    ),
    call = call
  )
}


#' Zhao--Feng--Wang--Wang robust maximum and adaptive alpha test
#'
#' Fits the simultaneous scaled spatial median to the restricted factor
#' residuals and computes the primary robust maximum statistic. Let r denote
#' the fitted standardized radii and
#' \eqn{\eta=T^{-1}1'P_F1}. The implemented nuisance constant is
#' \deqn{\widehat\zeta =
#' \frac{N\{\overline{r^{-1}}\}^2}
#' {1-2\eta\overline{r^{-1}}\bar r+
#' \eta\overline{r^{-2}}\,\overline{r^2}}.}
#' The centered maximum is
#' \deqn{T\widehat\zeta
#' \|\widehat D^{-1/2}\widehat\theta\|_\infty^2
#' -2\log N+\log\log N,}
#' calibrated by the upper tail of
#' \eqn{\exp\{-\pi^{-1/2}\exp(-x/2)\}}.
#'
#' The book draft reverses the squared inverse-radius moment in zeta and
#' replaces the primary final coefficient eta by eta squared. Both changes
#' alter the statistic; this implementation uses the primary formula and
#' returns every empirical radial moment. For component "combined", the
#' robust maximum p-value and [liu_feng_ma_spatial_sign_alpha_test()] p-value
#' are combined by the paper's truncated Cauchy rule: only component p-values
#' below one half contribute. This is distinct from the untruncated generic
#' Cauchy construction introduced elsewhere in the book.
#'
#' @inheritParams liu_feng_ma_spatial_sign_alpha_test
#' @param component "max" for the robust maximum alone or "combined" for
#'   the primary truncated-Cauchy max/sum procedure.
#' @return An upper-tail alpha-test object containing the complete scaled
#'   spatial-median fit, radial moments, corrected zeta, maximum calibration,
#'   and, for the combined procedure, the full LFM component.
#' @references
#' Zhao, P., Feng, L., Wang, Z., and Wang, X. Robust high-dimensional alpha
#' testing for linear factor pricing models. Oxford Bulletin of Economics and
#' Statistics (2026). \doi{10.1111/obes.70080}.
#' Preprint: \url{https://arxiv.org/abs/2408.06612}.
#' @examples
#' f <- cbind(seq(-1, 1, length.out = 12))
#' y <- outer(seq_len(12), 1:3, function(i, j) cos(i + j / 3))
#' zhao_feng_wang_wang_robust_alpha_test(y, f)
#' @export
zhao_feng_wang_wang_robust_alpha_test <- function(
    returns, factors = NULL, component = c("max", "combined"),
    bias = c("wild_bootstrap", "none", "supplied"),
    delta_q = NULL, bootstrap_reps = 100L, seed = NULL,
    keep_bootstrap = FALSE, tol = 1e-7, max_iter = 500L,
    zero_tol = 0) {
  call <- match.call()
  data.name <- deparse1(substitute(returns))
  component <- match.arg(component)
  bias <- match.arg(bias)
  controls <- .ch4_alpha_controls(tol, max_iter, zero_tol)
  prepared <- .ch4_alpha_prepare(returns, factors, min_assets = 2L)
  restricted <- .ch4_alpha_restricted_residuals(prepared)

  fit <- cpp_scaled_spatial_median(
    restricted$residuals, numeric(prepared$N), TRUE,
    controls$tol, controls$max_iter, controls$zero_tol
  )
  if (!isTRUE(fit$iteration_stable)) {
    stop(sprintf(
      paste0(
        "The Zhao-Feng-Wang-Wang scaled spatial median did not stabilize ",
        "within %d updates (score residual %.6g, relative update %.6g)."
      ),
      as.integer(fit$max_iterations), as.numeric(fit$score_residual),
      as.numeric(fit$relative_update)
    ), call. = FALSE)
  }
  radii <- as.numeric(fit$radii)
  if (anyNA(radii) || any(!is.finite(radii)) || any(radii <= 0)) {
    stop(
      "Every fitted radial distance must be strictly positive and finite; no deletion or floor is used.",
      call. = FALSE
    )
  }
  radial <- c(
    mean.r = mean(radii),
    mean.inverse.r = mean(1 / radii),
    mean.r.squared = mean(radii^2),
    mean.inverse.r.squared = mean(1 / radii^2)
  )
  if (any(!is.finite(radial))) {
    stop("A required empirical radial moment is not finite.",
         call. = FALSE)
  }
  zeta.denominator <- unname(
    1 - 2 * restricted$eta * radial["mean.inverse.r"] *
      radial["mean.r"] +
      restricted$eta * radial["mean.inverse.r.squared"] *
      radial["mean.r.squared"]
  )
  if (!is.finite(zeta.denominator) || zeta.denominator <= 0) {
    stop(
      "The primary zeta denominator is not strictly positive and finite; no absolute value or floor is used.",
      call. = FALSE
    )
  }
  zeta <- prepared$N * unname(radial["mean.inverse.r"])^2 /
    zeta.denominator
  coordinates <- as.numeric(fit$standardized_location_coordinate)
  maximum.square <- max(coordinates^2)
  raw.maximum <- prepared$T * zeta * maximum.square
  centered <- raw.maximum - 2 * log(prepared$N) +
    log(log(prepared$N))
  if (!is.finite(zeta) || zeta <= 0 || !is.finite(raw.maximum) ||
      !is.finite(centered)) {
    stop("The robust maximum statistic is not finite and positive.",
         call. = FALSE)
  }
  max.tail <- .ch4_alpha_gumbel_tail(centered)

  log.diagonal <- as.numeric(fit$log_diagonal_input_canonical) +
    2 * log(prepared$y.scale)
  log.diagonal <- log.diagonal - max(log.diagonal)
  diagonal <- exp(log.diagonal)
  names(log.diagonal) <- names(diagonal) <- prepared$asset.names
  theta <- as.numeric(fit$location_relative_origin) * prepared$y.scale
  names(theta) <- prepared$asset.names
  names(coordinates) <- prepared$asset.names
  directions <- as.matrix(fit$directions)
  dimnames(directions) <- list(rownames(prepared$returns),
                               prepared$asset.names)

  sum.component <- NULL
  combination <- NULL
  final.statistic <- stats::setNames(centered, "Gumbel centered")
  final.p <- max.tail$p.value
  final.method <- "Zhao-Feng-Wang-Wang robust maximum alpha test"
  if (component == "combined") {
    sum.component <- liu_feng_ma_spatial_sign_alpha_test(
      returns, factors, bias = bias, delta_q = delta_q,
      bootstrap_reps = bootstrap_reps, seed = seed,
      keep_bootstrap = keep_bootstrap, tol = controls$tol,
      max_iter = controls$max_iter, zero_tol = controls$zero_tol
    )
    combination <- .ch4_alpha_truncated_cauchy(c(
      SS = sum.component$p.value,
      SM = max.tail$p.value
    ))
    names(combination$terms) <- names(combination$active) <- c("SS", "SM")
    final.statistic <- stats::setNames(
      combination$statistic, "truncated Cauchy"
    )
    final.p <- combination$p.value
    final.method <- paste(
      "Zhao-Feng-Wang-Wang robust adaptive alpha test",
      "(primary truncated-Cauchy SS/SM combination)"
    )
  }

  .ch4_alpha_htest(
    statistic = final.statistic,
    p.value = final.p,
    method = final.method,
    data.name = data.name,
    estimate = theta,
    raw.statistic = c(
      robust.maximum = raw.maximum,
      Gumbel.centered = centered
    ),
    components = list(
      theta = theta,
      theta.scaled.coordinates = as.numeric(fit$location_relative_origin),
      standardized.theta = coordinates,
      scale.diagonal = diagonal,
      log.scale.diagonal = log.diagonal,
      scale.diagonal.internal = as.numeric(fit$diagonal_standardized),
      directions = directions,
      radii.internal = radii,
      radial.moments = radial,
      eta = restricted$eta,
      omega = restricted$omega,
      zeta.denominator = zeta.denominator,
      zeta = zeta,
      maximum.standardized.square = maximum.square,
      robust.maximum = raw.maximum,
      centered.maximum = centered,
      p.SM = max.tail$p.value,
      log.p.SM = max.tail$log.p.value,
      restricted.slopes.scaled = restricted$slopes,
      restricted.residuals.scaled = restricted$residuals,
      sum.test = sum.component,
      p.SS = if (is.null(sum.component)) NULL else sum.component$p.value,
      truncated.cauchy = combination,
      T = prepared$T,
      N = prepared$N,
      K = prepared$K
    ),
    diagnostics = list(
      component = component,
      max.calibration = "upper tail of primary type-I extreme-value limit",
      zeta.formula = paste(
        "N*mean(r^-1)^2/[1-2*eta*mean(r^-1)*mean(r)",
        "+eta*mean(r^-2)*mean(r^2)]"
      ),
      book.zeta.inverse.moment.conflict = TRUE,
      book.zeta.eta.squared.conflict = TRUE,
      book.generic.Cauchy.implemented = FALSE,
      combined.calibration = if (component == "combined") {
        "primary truncated-Cauchy combination"
      } else {
        NULL
      },
      iteration.stable = isTRUE(fit$iteration_stable),
      iterations = as.integer(fit$iterations),
      score.residual = as.numeric(fit$score_residual),
      relative.update = as.numeric(fit$relative_update),
      zero.radius.handling = "strict failure",
      radius.floor = "none",
      ridge = "none",
      generalized.inverse = "none"
    ),
    call = call
  )
}
