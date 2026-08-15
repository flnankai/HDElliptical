.ch4ind_validate_matrix <- function(x, name, min_rows = 2L,
                                    min_cols = 1L) {
  if (is.data.frame(x)) x <- as.matrix(x)
  if (is.atomic(x) && is.null(dim(x))) x <- matrix(x, ncol = 1L)
  if (!is.matrix(x) || !is.numeric(x)) {
    stop(sprintf("`%s` must be a numeric matrix.", name), call. = FALSE)
  }
  storage.mode(x) <- "double"
  if (nrow(x) < min_rows || ncol(x) < min_cols) {
    stop(sprintf(
      "`%s` must have at least %d rows and %d column%s.",
      name, min_rows, min_cols, if (min_cols == 1L) "" else "s"
    ), call. = FALSE)
  }
  if (anyNA(x) || any(!is.finite(x))) {
    stop(sprintf("`%s` must contain only finite values.", name),
         call. = FALSE)
  }
  x
}

.ch4ind_validate_logical <- function(value, name) {
  if (!is.logical(value) || length(value) != 1L || is.na(value)) {
    stop(sprintf("`%s` must be TRUE or FALSE.", name), call. = FALSE)
  }
  value
}

.ch4ind_validate_integer <- function(value, name, minimum = 0L) {
  if (!is.numeric(value) || length(value) != 1L || is.na(value) ||
      !is.finite(value) || value != floor(value) || value < minimum ||
      value > .Machine$integer.max) {
    stop(sprintf(
      "`%s` must be one integer between %d and .Machine$integer.max.",
      name, minimum
    ), call. = FALSE)
  }
  as.integer(value)
}

.ch4ind_scale_columns <- function(x, name) {
  scales <- apply(abs(x), 2L, max)
  if (any(!is.finite(scales)) || any(scales <= 0)) {
    stop(sprintf("Every column of `%s` must have a positive finite scale.",
                 name), call. = FALSE)
  }
  list(value = sweep(x, 2L, scales, "/"), scales = scales)
}

.ch4ind_scale_global <- function(x, name) {
  scale <- max(abs(x))
  if (!is.finite(scale) || scale <= 0) {
    stop(sprintf("`%s` must have a positive finite global scale.", name),
         call. = FALSE)
  }
  list(value = x / scale, scale = scale)
}

.ch4ind_names <- function(x, prefix) {
  if (is.null(colnames(x))) paste0(prefix, seq_len(ncol(x))) else colnames(x)
}

.ch4ind_time_names <- function(x) {
  if (is.null(rownames(x))) paste0("time", seq_len(nrow(x))) else rownames(x)
}

.ch4ind_extreme_tail <- function(statistic, log_constant) {
  log.lambda <- log_constant - statistic / 2
  if (log.lambda >= log(.Machine$double.xmax)) {
    return(list(p.value = 1, log.p.value = 0, cdf = 0,
                log.lambda = log.lambda))
  }
  if (log.lambda <= log(.Machine$double.xmin)) {
    return(list(p.value = exp(log.lambda), log.p.value = log.lambda,
                cdf = 1, log.lambda = log.lambda))
  }
  lambda <- exp(log.lambda)
  p.value <- -expm1(-lambda)
  list(
    p.value = p.value,
    log.p.value = log(p.value),
    cdf = exp(-lambda),
    log.lambda = log.lambda
  )
}

.ch4ind_new_test <- function(statistic, p.value, alternative, method,
                             data.name, raw.statistic, components,
                             diagnostics, call, estimate = NULL,
                             parameter = NULL) {
  answer <- list(
    statistic = statistic,
    p.value = as.numeric(p.value),
    alternative = alternative,
    method = method,
    data.name = data.name,
    raw.statistic = raw.statistic,
    null.value = c(dependence = 0),
    components = components,
    diagnostics = diagnostics,
    call = call
  )
  if (!is.null(estimate)) answer$estimate <- estimate
  if (!is.null(parameter)) answer$parameter <- parameter
  class(answer) <- c("hd_independence_test", "htest")
  answer
}

.ch4ind_design_qr <- function(design, time, label) {
  design <- .ch4ind_validate_matrix(design, label, min_rows = time,
                                    min_cols = 1L)
  if (nrow(design) != time) {
    stop(sprintf("`%s` must have exactly T rows.", label), call. = FALSE)
  }
  if (ncol(design) >= time) {
    stop(sprintf("`%s` must have fewer than T columns.", label),
         call. = FALSE)
  }
  scaled <- .ch4ind_scale_columns(design, label)$value
  fit <- qr(scaled, tol = sqrt(.Machine$double.eps), LAPACK = FALSE)
  if (fit$rank != ncol(scaled)) {
    stop(sprintf("`%s` must have full column rank.", label),
         call. = FALSE)
  }
  list(
    qr = fit,
    q = qr.Q(fit, complete = FALSE)[, seq_len(ncol(scaled)), drop = FALSE],
    rank = ncol(scaled)
  )
}

.ch4ind_residualize <- function(panel, regressors, scaling) {
  panel <- .ch4ind_validate_matrix(panel, "panel", min_rows = 2L,
                                   min_cols = 2L)
  time <- nrow(panel)
  units <- ncol(panel)
  scaled <- if (scaling == "column") {
    .ch4ind_scale_columns(panel, "panel")$value
  } else {
    .ch4ind_scale_global(panel, "panel")$value
  }
  if (is.null(regressors)) {
    return(list(
      residuals = scaled, rank = 0L, q = NULL, common = TRUE,
      design.route = "supplied residuals (no regression projection)"
    ))
  }

  if (is.matrix(regressors) || is.data.frame(regressors)) {
    design <- .ch4ind_design_qr(regressors, time, "regressors")
    residuals <- qr.resid(design$qr, scaled)
    return(list(
      residuals = residuals, rank = design$rank,
      q = list(design$q), common = TRUE,
      design.route = "common full-rank supplied design"
    ))
  }

  if (!is.list(regressors) || length(regressors) != units) {
    stop("`regressors` must be NULL, one T by p matrix, or a list of " %+%
           "one T by p matrix per panel unit.", call. = FALSE)
  }
  fits <- vector("list", units)
  residuals <- matrix(0, time, units, dimnames = dimnames(panel))
  ranks <- integer(units)
  for (unit in seq_len(units)) {
    fits[[unit]] <- .ch4ind_design_qr(
      regressors[[unit]], time, sprintf("regressors[[%d]]", unit)
    )
    ranks[unit] <- fits[[unit]]$rank
    residuals[, unit] <- qr.resid(fits[[unit]]$qr, scaled[, unit])
  }
  if (length(unique(ranks)) != 1L) {
    stop("All unit-specific designs must have the same column rank p.",
         call. = FALSE)
  }
  list(
    residuals = residuals, rank = ranks[1L],
    q = lapply(fits, `[[`, "q"), common = FALSE,
    design.route = "unit-specific full-rank supplied designs"
  )
}

.ch4ind_projection_overlap <- function(fit, time, units) {
  pairs <- units * (units - 1) / 2
  rank <- fit$rank
  if (rank == 0L) return(pairs * time)
  if (fit$common) return(pairs * (time - rank))
  answer <- 0
  for (first in seq_len(units - 1L)) {
    for (second in (first + 1L):units) {
      overlap <- crossprod(fit$q[[first]], fit$q[[second]])
      answer <- answer + time - 2 * rank + sum(overlap^2)
    }
  }
  answer
}

.ch4ind_with_local_seed <- function(seed, expression) {
  if (is.null(seed)) return(force(expression))
  had.state <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had.state) old.state <- get(".Random.seed", envir = .GlobalEnv)
  on.exit({
    if (had.state) {
      assign(".Random.seed", old.state, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv,
                      inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  set.seed(seed)
  force(expression)
}


#' Gaussian Wilks block-independence test
#'
#' Computes the classical Gaussian likelihood-ratio test for independence of
#' two vector blocks. If \eqn{\hat\rho_1,\ldots,\hat\rho_m} are the sample
#' canonical correlations, \eqn{m=\min(p,q)}, then
#' \deqn{\Lambda=\prod_{j=1}^m(1-\hat\rho_j^2)}
#' and the Bartlett statistic is
#' \deqn{-\{n-1-(p+q+1)/2\}\log\Lambda,}
#' calibrated against \eqn{\chi^2_{pq}}. This is a fixed-dimensional,
#' Gaussian benchmark; it is not a high-dimensional repair of Wilks' test.
#'
#' @param x Numeric \eqn{n} by \eqn{p} matrix for the first block.
#' @param y Numeric \eqn{n} by \eqn{q} matrix for the second block.
#'
#' @return An object inheriting from `htest`. `components` contains Wilks'
#'   Lambda, its log value, canonical correlations, and the three sample
#'   covariance blocks.
#' @references Anderson, T. W. (2003). *An Introduction to Multivariate
#'   Statistical Analysis*, 3rd edition, Chapter 8. Wiley.
#' @examples
#' x <- matrix(c(-2, 0, -1, 2, 0, -1, 1, 1, 2, -2, 3, 0), 6, 2)
#' y <- matrix(c(1, -2, 0, 2, -1, 3), 6, 1)
#' gaussian_wilks_independence_test(x, y)
#' @export
gaussian_wilks_independence_test <- function(x, y) {
  call <- match.call()
  x.name <- deparse(substitute(x))
  y.name <- deparse(substitute(y))
  x <- .ch4ind_validate_matrix(x, "x", min_rows = 3L, min_cols = 1L)
  y <- .ch4ind_validate_matrix(y, "y", min_rows = 3L, min_cols = 1L)
  if (nrow(x) != nrow(y)) {
    stop("`x` and `y` must have the same number of rows.", call. = FALSE)
  }
  n <- nrow(x)
  p <- ncol(x)
  q <- ncol(y)
  if (n <= p + q) {
    stop("Classical Wilks calibration requires n > p + q.",
         call. = FALSE)
  }
  x.scaled <- .ch4ind_scale_columns(x, "x")
  y.scaled <- .ch4ind_scale_columns(y, "y")
  core <- cpp_ch4ind_wilks_core(x.scaled$value, y.scaled$value)
  multiplier <- n - 1 - (p + q + 1) / 2
  if (multiplier <= 0) {
    stop("The Bartlett multiplier is not positive.", call. = FALSE)
  }
  statistic <- -multiplier * core$log_lambda
  p.value <- stats::pchisq(statistic, df = p * q, lower.tail = FALSE)
  names(core$canonical_correlations) <- paste0(
    "canonical", seq_along(core$canonical_correlations)
  )
  dimnames(core$covariance_x) <- list(.ch4ind_names(x, "x"),
                                      .ch4ind_names(x, "x"))
  dimnames(core$covariance_y) <- list(.ch4ind_names(y, "y"),
                                      .ch4ind_names(y, "y"))
  dimnames(core$cross_covariance) <- list(.ch4ind_names(x, "x"),
                                          .ch4ind_names(y, "y"))
  .ch4ind_new_test(
    statistic = stats::setNames(statistic, "Bartlett Wilks"),
    p.value = p.value,
    alternative = "greater",
    method = "Gaussian Wilks test of block independence",
    data.name = paste(x.name, "and", y.name),
    raw.statistic = c(lambda = core$lambda,
                      log.lambda = core$log_lambda),
    components = c(core, list(bartlett.multiplier = multiplier)),
    diagnostics = list(
      calibration = "fixed-dimensional Gaussian chi-square",
      degrees.of.freedom = p * q,
      numerical.repair = "none; singular cross-products fail explicitly",
      x.column.scales = x.scaled$scales,
      y.column.scales = y.scaled$scales
    ),
    call = call,
    estimate = stats::setNames(core$lambda, "Wilks Lambda"),
    parameter = c(df = p * q)
  )
}


#' Pesaran CD test for cross-sectional independence
#'
#' Given a \eqn{T} by \eqn{N} residual matrix, computes
#' \deqn{CD=\sqrt{2T/\{N(N-1)\}}\sum_{i<j}\hat\rho_{ij}.}
#' The standard-normal calibration is two-sided. Residuals are used exactly as
#' supplied: include an intercept in the preceding unit regressions when
#' centering is required by the model.
#'
#' @param residuals Numeric \eqn{T} by \eqn{N} matrix, with time in rows and
#'   panel units in columns.
#' @param keep_correlations Whether to retain the \eqn{N} by \eqn{N} sample
#'   residual-correlation matrix.
#'
#' @return An object inheriting from `htest`.
#' @references Pesaran, M. H. (2004). General Diagnostic Tests for Cross
#'   Section Dependence in Panels. IZA Discussion Paper 1240.
#'   \url{https://docs.iza.org/dp1240.pdf}
#' @examples
#' e <- matrix(c(-2, 1, 0, 2, -1, 3, 1, -2, 2, 1, -3, 1), 4, 3)
#' pesaran_cd_test(e)
#' @export
pesaran_cd_test <- function(residuals, keep_correlations = FALSE) {
  call <- match.call()
  data.name <- deparse(substitute(residuals))
  residuals <- .ch4ind_validate_matrix(
    residuals, "residuals", min_rows = 2L, min_cols = 2L
  )
  keep_correlations <- .ch4ind_validate_logical(
    keep_correlations, "keep_correlations"
  )
  scaled <- .ch4ind_scale_columns(residuals, "residuals")
  core <- cpp_ch4ind_pairwise_core(scaled$value, keep_correlations)
  time <- nrow(residuals)
  units <- ncol(residuals)
  statistic <- sqrt(2 * time / (units * (units - 1))) * core$sum
  log.half <- stats::pnorm(abs(statistic), lower.tail = FALSE, log.p = TRUE)
  p.value <- min(1, 2 * exp(log.half))
  unit.names <- .ch4ind_names(residuals, "unit")
  names(core$maximum_pair) <- c("first", "second")
  core$maximum_pair.names <- unit.names[core$maximum_pair]
  if (keep_correlations) {
    dimnames(core$correlations) <- list(unit.names, unit.names)
  }
  .ch4ind_new_test(
    statistic = stats::setNames(statistic, "CD"),
    p.value = p.value,
    alternative = "two.sided",
    method = "Pesaran CD test for cross-sectional independence",
    data.name = data.name,
    raw.statistic = c(sum.pairwise.correlations = core$sum),
    components = c(core, list(
      time.points = time,
      units = units,
      standardization = sqrt(2 * time / (units * (units - 1)))
    )),
    diagnostics = list(
      calibration = "standard normal",
      residuals.centered.by.function = FALSE,
      cancellation.possible = TRUE,
      numerical.repair = "none"
    ),
    call = call,
    estimate = stats::setNames(
      core$sum / (units * (units - 1) / 2),
      "mean pairwise residual correlation"
    )
  )
}


#' Feng-Jiang-Liu-Xiong max-sum panel-independence test
#'
#' Implements the three procedures in Feng, Jiang, Liu and Xiong (2022) for
#' serially uncorrelated panel errors. For OLS residual correlations
#' \eqn{\hat\rho_{ij}},
#' \deqn{S_N=\sum_{i<j}T\hat\rho_{ij}^2,\qquad
#' L_N=\max_{i<j}|\hat\rho_{ij}|.}
#' The sum component is \eqn{(S_N-\mu_N)/N}, where
#' \deqn{\mu_N=\frac{T}{(T-p)^2}\sum_{i<j}\mathrm{tr}(P_iP_j),}
#' and the max component is
#' \eqn{TL_N^2-4\log N+\log\log N}. Their primary max-sum statistic is
#' \eqn{C_N=\min(p_L,p_S)} with calibrated p-value
#' \eqn{2C_N-C_N^2}.
#'
#' @param panel Numeric \eqn{T} by \eqn{N} matrix. With `regressors = NULL`,
#'   columns are treated as already-computed residual vectors. Otherwise they
#'   are unit outcomes and OLS residuals are computed internally.
#' @param regressors `NULL`, one common \eqn{T} by \eqn{p} design matrix, or a
#'   list of \eqn{N} unit-specific \eqn{T} by \eqn{p} full-rank designs. The
#'   designs must have the same \eqn{p}; no intercept is added implicitly.
#' @param component One of `"max-sum"`, `"max"`, or `"sum"`.
#' @param keep_correlations Whether to retain the residual-correlation matrix.
#'
#' @return An object inheriting from `htest`, with all three component
#'   statistics and p-values retained in `components`.
#' @references Feng, L., Jiang, T., Liu, B. and Xiong, W. (2022). Max-Sum
#'   Tests for Cross-Sectional Independence of High-Dimensional Panel Data.
#'   *Annals of Statistics*, 50, 1124-1143.
#'   \doi{10.1214/21-AOS2142}
#' @examples
#' panel <- matrix(c(-2, 1, 0, 2, -1, 3, 1, -2,
#'                   2, 1, -3, 1, 3, -1, 2, -2), 4, 4)
#' feng_jiang_liu_xiong_panel_independence_test(panel)
#' @export
feng_jiang_liu_xiong_panel_independence_test <- function(
    panel, regressors = NULL,
    component = c("max-sum", "max", "sum"),
    keep_correlations = FALSE) {
  call <- match.call()
  data.name <- deparse(substitute(panel))
  component <- match.arg(component)
  keep_correlations <- .ch4ind_validate_logical(
    keep_correlations, "keep_correlations"
  )
  fit <- .ch4ind_residualize(panel, regressors, scaling = "column")
  core <- cpp_ch4ind_pairwise_core(fit$residuals, keep_correlations)
  time <- nrow(fit$residuals)
  units <- ncol(fit$residuals)
  residual.df <- time - fit$rank
  if (residual.df <= 0L) {
    stop("The residual degrees of freedom T - p must be positive.",
         call. = FALSE)
  }
  overlap.sum <- .ch4ind_projection_overlap(fit, time, units)
  mu <- time * overlap.sum / residual.df^2
  sum.raw <- time * core$sum_squares
  sum.z <- (sum.raw - mu) / units
  max.gumbel <- time * core$maximum^2 - 4 * log(units) + log(log(units))
  max.tail <- .ch4ind_extreme_tail(
    max.gumbel, -0.5 * log(8 * pi)
  )
  log.p.sum <- stats::pnorm(sum.z, lower.tail = FALSE, log.p = TRUE)
  p.sum <- exp(log.p.sum)
  c.statistic <- min(max.tail$p.value, p.sum)
  p.combined <- c.statistic * (2 - c.statistic)

  unit.names <- .ch4ind_names(fit$residuals, "unit")
  names(core$maximum_pair) <- c("first", "second")
  core$maximum_pair.names <- unit.names[core$maximum_pair]
  if (keep_correlations) {
    dimnames(core$correlations) <- list(unit.names, unit.names)
  }
  components <- c(core, list(
    time.points = time,
    units = units,
    regression.rank = fit$rank,
    residual.degrees.of.freedom = residual.df,
    projection.overlap.sum = overlap.sum,
    mu = mu,
    sum.raw = sum.raw,
    sum.z = sum.z,
    p.sum = p.sum,
    log.p.sum = log.p.sum,
    max.gumbel = max.gumbel,
    p.max = max.tail$p.value,
    log.p.max = max.tail$log.p.value,
    max.cdf = max.tail$cdf,
    C.N = c.statistic,
    p.max.sum = p.combined
  ))
  selected <- switch(
    component,
    `max-sum` = list(statistic = c(C.N = c.statistic),
                     p.value = p.combined, alternative = "less"),
    max = list(statistic = c(`max Gumbel` = max.gumbel),
               p.value = max.tail$p.value, alternative = "greater"),
    sum = list(statistic = c(`sum Z` = sum.z),
               p.value = p.sum, alternative = "greater")
  )
  .ch4ind_new_test(
    statistic = selected$statistic,
    p.value = selected$p.value,
    alternative = selected$alternative,
    method = paste0(
      "Feng-Jiang-Liu-Xiong panel independence test (", component, ")"
    ),
    data.name = data.name,
    raw.statistic = c(S.N = sum.raw, L.N = core$maximum),
    components = components,
    diagnostics = list(
      selected.component = component,
      design.route = fit$design.route,
      serial.correlation.allowed = FALSE,
      sum.calibration = "Gaussian CLT with exact projection centering",
      max.calibration = "type-I extreme value; constant 1/sqrt(8*pi)",
      combination = "minimum p-value with cdf 2*c-c^2",
      numerical.repair = "none"
    ),
    call = call,
    estimate = c(maximum.absolute.correlation = core$maximum)
  )
}


#' Wang-Liu-Feng-Ma serial-panel Fisher independence test
#'
#' Implements the serial-correlation-aware panel procedure of Wang, Liu, Feng
#' and Ma. The dense component is the signed-correlation sum
#' \deqn{S_N=\sqrt{2/\{N(N-1)\}}\sum_{i<j}\hat\rho_{ij}}
#' with the published leave-two-units-out variance estimator. The sparse
#' component uses \eqn{L_N=\max_{i<j}\hat\rho_{ij}^2} and the temporal effective
#' dimension \eqn{\mathrm{tr}^2(\widetilde\Sigma)/
#' \|\widetilde\Sigma\|_F^2}. Fisher's statistic combines the two upper-tail
#' p-values and is calibrated by \eqn{\chi^2_4}.
#'
#' The internal temporal estimator is implemented literally: the cross-unit
#' sample covariance is hard-thresholded using the paper's \eqn{\hat P_N} and
#' \eqn{\nu>\sqrt 2}. No positive-definite projection, ridge, or replacement of
#' a non-positive \eqn{\hat P_N} is applied. A scientifically justified
#' temporal covariance estimate may instead be supplied explicitly.
#'
#' @param panel Numeric \eqn{T} by \eqn{N} residual or outcome matrix.
#' @param regressors As in
#'   [feng_jiang_liu_xiong_panel_independence_test()].
#' @param temporal_covariance Optional supplied positive-definite \eqn{T} by
#'   \eqn{T} temporal covariance. If `NULL`, use the paper's thresholded sample
#'   estimator.
#' @param nu Threshold constant, strictly greater than \eqn{\sqrt 2}; used only
#'   by the internal temporal estimator. The default `1.42` is the value used
#'   in the primary paper's applications.
#' @param component One of `"fisher"`, `"max"`, or `"sum"`.
#' @param keep_matrices Whether to retain correlation and temporal-estimation
#'   matrices.
#'
#' @return An object inheriting from `htest`.
#' @references Wang, H., Liu, B., Feng, L. and Ma, Y. (2026). Fisher's
#'   Combined Probability Test for Cross-Sectional Independence in Panel Data
#'   Models with Serial Correlation. *Statistica Sinica*, 36, 1-21.
#'   \doi{10.5705/ss.202023.0348}
#' @examples
#' panel <- matrix(c(-2, 1, 0, 2, -1, 3, 1, -2,
#'                   2, 1, -3, 1, 3, -1, 2, -2,
#'                   1, 3, -2, 1), 5, 4)
#' wang_liu_feng_ma_serial_panel_test(
#'   panel, temporal_covariance = diag(nrow(panel))
#' )
#' @export
wang_liu_feng_ma_serial_panel_test <- function(
    panel, regressors = NULL, temporal_covariance = NULL, nu = 1.42,
    component = c("fisher", "max", "sum"), keep_matrices = FALSE) {
  call <- match.call()
  data.name <- deparse(substitute(panel))
  component <- match.arg(component)
  keep_matrices <- .ch4ind_validate_logical(keep_matrices, "keep_matrices")
  if (!is.numeric(nu) || length(nu) != 1L || is.na(nu) ||
      !is.finite(nu) || nu <= sqrt(2)) {
    stop("`nu` must be one finite number strictly greater than sqrt(2).",
         call. = FALSE)
  }
  fit <- .ch4ind_residualize(panel, regressors, scaling = "global")
  if (ncol(fit$residuals) < 3L) {
    stop("Serial-panel calibration requires at least three panel units.",
         call. = FALSE)
  }
  estimate.temporal <- is.null(temporal_covariance)
  if (estimate.temporal) {
    temporal <- matrix(numeric(), 0L, 0L)
  } else {
    temporal <- .ch4ind_validate_matrix(
      temporal_covariance, "temporal_covariance",
      min_rows = nrow(fit$residuals), min_cols = nrow(fit$residuals)
    )
    if (!identical(dim(temporal), rep(nrow(fit$residuals), 2L))) {
      stop("`temporal_covariance` must be a T by T matrix.",
           call. = FALSE)
    }
  }
  core <- cpp_ch4ind_serial_panel_core(
    fit$residuals, nu, estimate.temporal, temporal, keep_matrices
  )
  units <- ncol(fit$residuals)
  sum.z <- core$sum_statistic / sqrt(core$sum_variance)
  log.p.sum <- stats::pnorm(sum.z, lower.tail = FALSE, log.p = TRUE)
  p.sum <- exp(log.p.sum)
  max.gumbel <- core$effective_dimension * core$maximum_square -
    4 * log(units) + log(log(units))
  max.tail <- .ch4ind_extreme_tail(
    max.gumbel, -0.5 * log(8 * pi)
  )
  fisher <- -2 * (log.p.sum + max.tail$log.p.value)
  p.fisher <- stats::pchisq(fisher, df = 4, lower.tail = FALSE)

  unit.names <- .ch4ind_names(fit$residuals, "unit")
  time.names <- .ch4ind_time_names(fit$residuals)
  names(core$maximum_pair) <- c("first", "second")
  core$maximum_pair.names <- unit.names[core$maximum_pair]
  if (keep_matrices) {
    dimnames(core$correlations) <- list(unit.names, unit.names)
    if (!is.null(core$sigma_hat)) {
      dimnames(core$sigma_hat) <- list(time.names, time.names)
    }
    dimnames(core$sigma_tilde) <- list(time.names, time.names)
    if (!is.null(core$u_hat)) {
      dimnames(core$u_hat) <- list(unit.names, unit.names)
    }
  }
  components <- c(core, list(
    time.points = nrow(fit$residuals),
    units = units,
    regression.rank = fit$rank,
    sum.z = sum.z,
    p.sum = p.sum,
    log.p.sum = log.p.sum,
    max.gumbel = max.gumbel,
    p.max = max.tail$p.value,
    log.p.max = max.tail$log.p.value,
    max.cdf = max.tail$cdf,
    fisher = fisher,
    p.fisher = p.fisher
  ))
  selected <- switch(
    component,
    fisher = list(statistic = c(Fisher = fisher), p.value = p.fisher),
    max = list(statistic = c(`max Gumbel` = max.gumbel),
               p.value = max.tail$p.value),
    sum = list(statistic = c(`sum Z` = sum.z), p.value = p.sum)
  )
  .ch4ind_new_test(
    statistic = selected$statistic,
    p.value = selected$p.value,
    alternative = "greater",
    method = paste0(
      "Wang-Liu-Feng-Ma serial-panel independence test (", component, ")"
    ),
    data.name = data.name,
    raw.statistic = c(S.N = core$sum_statistic,
                      L.N = core$maximum_square),
    components = components,
    diagnostics = list(
      selected.component = component,
      design.route = fit$design.route,
      temporal.route = if (estimate.temporal) {
        "literal primary hard-thresholded sample estimator"
      } else {
        "supplied positive-definite temporal covariance"
      },
      threshold.nu = if (estimate.temporal) nu else NA_real_,
      serial.correlation.allowed = TRUE,
      max.calibration = "type-I extreme value; constant 1/sqrt(8*pi)",
      sum.calibration = "Baltagi-Kao-Peng signed-sum variance",
      combination = "Fisher chi-square with 4 degrees of freedom",
      numerical.repair = "none; no ridge or positive-definite projection"
    ),
    call = call,
    estimate = c(maximum.squared.correlation = core$maximum_square),
    parameter = if (component == "fisher") c(df = 4) else NULL
  )
}


#' Wang-Liu-Feng rank max-sum test for vector independence
#'
#' Tests independence between two high-dimensional random vectors using the
#' Spearman or Kendall procedures of Wang, Liu and Feng. For every cross-block
#' coordinate pair, the function computes the exact primary rank correlation.
#' It then forms the max statistic, the centered sum of squared correlations,
#' and their Fisher combination. The max null variance is \eqn{1/(n-1)} for
#' Spearman and \eqn{2(2n+5)/\{9n(n-1)\}} for Kendall. The sum variance is
#' estimated by the paper's intrinsic permutation calibration; this is part of
#' the callable method, not a replication simulation.
#'
#' Exact ties are rejected because the primary finite-sample null moments and
#' distribution-free calibration assume continuous margins. The chapter's
#' general mutual-independence construction and the degenerate Hoeffding D,
#' Blum-Kiefer-Rosenblatt R, and Bergsma-Dassios-Yanagimoto tau-star families
#' remain review-only: this function does not invent missing executable
#' variance/eigenspectrum contracts for them.
#'
#' @param x Numeric \eqn{n} by \eqn{p} matrix.
#' @param y Numeric \eqn{n} by \eqn{q} matrix with the same rows as `x`.
#' @param measure Either `"spearman"` or `"kendall"`.
#' @param component One of `"fisher"`, `"max"`, or `"sum"`.
#' @param B Number of intrinsic permutations used to estimate the sum
#'   variance; at least two.
#' @param seed Optional non-negative integer. An explicit seed is localized and
#'   preserves the caller's R random-number state. `NULL` uses the caller's
#'   stream normally.
#' @param keep_correlations Whether to retain the \eqn{p} by \eqn{q} matrix of
#'   pairwise rank correlations.
#' @param keep_permutation Whether to retain the intrinsic permutation
#'   statistics and permutation index matrix.
#'
#' @return An object inheriting from `htest`, with max, sum, and Fisher
#'   components.
#' @references Wang, H., Liu, B. and Feng, L. (2026). Testing Independence
#'   Between High-Dimensional Random Vectors Using Rank-Based Max-Sum Tests.
#'   *Scandinavian Journal of Statistics*, 53, 821-847.
#'   \doi{10.1111/sjos.70063}
#' @examples
#' x <- matrix(c(1, 4, 2, 6, 3, 5, 2, 6, 1, 5, 3, 4), 6, 2)
#' y <- matrix(c(6, 2, 5, 1, 4, 3), 6, 1)
#' wang_liu_feng_vector_independence_test(x, y, B = 19, seed = 7)
#' @export
wang_liu_feng_vector_independence_test <- function(
    x, y, measure = c("spearman", "kendall"),
    component = c("fisher", "max", "sum"), B = 199L, seed = NULL,
    keep_correlations = FALSE, keep_permutation = FALSE) {
  call <- match.call()
  x.name <- deparse(substitute(x))
  y.name <- deparse(substitute(y))
  measure <- match.arg(measure)
  component <- match.arg(component)
  x <- .ch4ind_validate_matrix(x, "x", min_rows = 3L, min_cols = 1L)
  y <- .ch4ind_validate_matrix(y, "y", min_rows = 3L, min_cols = 1L)
  if (nrow(x) != nrow(y)) {
    stop("`x` and `y` must have the same number of rows.", call. = FALSE)
  }
  comparisons <- ncol(x) * ncol(y)
  if (comparisons < 2L) {
    stop("The high-dimensional max calibration requires p * q >= 2.",
         call. = FALSE)
  }
  B <- .ch4ind_validate_integer(B, "B", minimum = 2L)
  if (!is.null(seed)) {
    seed <- .ch4ind_validate_integer(seed, "seed", minimum = 0L)
  }
  keep_correlations <- .ch4ind_validate_logical(
    keep_correlations, "keep_correlations"
  )
  keep_permutation <- .ch4ind_validate_logical(
    keep_permutation, "keep_permutation"
  )
  draw.permutations <- function() {
    answer <- replicate(B, sample.int(nrow(x)), simplify = "matrix")
    if (B == 1L) answer <- matrix(answer, ncol = 1L)
    storage.mode(answer) <- "integer"
    answer
  }
  permutations <- .ch4ind_with_local_seed(seed, draw.permutations())
  core <- cpp_ch4ind_rank_vector_core(
    x, y, match(measure, c("spearman", "kendall")) - 1L,
    permutations, keep_correlations, keep_permutation
  )
  sum.z <- core$sum_statistic / sqrt(core$permutation_variance)
  log.p.sum <- stats::pnorm(sum.z, lower.tail = FALSE, log.p = TRUE)
  p.sum <- exp(log.p.sum)
  max.gumbel <- core$maximum^2 / core$null_second_moment -
    2 * log(comparisons) + log(log(comparisons))
  max.tail <- .ch4ind_extreme_tail(max.gumbel, -0.5 * log(pi))
  fisher <- -2 * (log.p.sum + max.tail$log.p.value)
  p.fisher <- stats::pchisq(fisher, df = 4, lower.tail = FALSE)

  x.names <- .ch4ind_names(x, "x")
  y.names <- .ch4ind_names(y, "y")
  names(core$maximum_index) <- c("x", "y")
  core$maximum.names <- c(
    x.names[core$maximum_index[1L]], y.names[core$maximum_index[2L]]
  )
  if (keep_correlations) {
    dimnames(core$correlations) <- list(x.names, y.names)
  }
  if (keep_permutation) {
    colnames(permutations) <- paste0("permutation", seq_len(B))
    core$permutation.indices <- permutations
  }
  components <- c(core, list(
    observations = nrow(x),
    x.dimension = ncol(x),
    y.dimension = ncol(y),
    comparisons = comparisons,
    sum.z = sum.z,
    p.sum = p.sum,
    log.p.sum = log.p.sum,
    max.gumbel = max.gumbel,
    p.max = max.tail$p.value,
    log.p.max = max.tail$log.p.value,
    max.cdf = max.tail$cdf,
    fisher = fisher,
    p.fisher = p.fisher
  ))
  selected <- switch(
    component,
    fisher = list(statistic = c(Fisher = fisher), p.value = p.fisher),
    max = list(statistic = c(`max Gumbel` = max.gumbel),
               p.value = max.tail$p.value),
    sum = list(statistic = c(`sum Z` = sum.z), p.value = p.sum)
  )
  .ch4ind_new_test(
    statistic = selected$statistic,
    p.value = selected$p.value,
    alternative = "greater",
    method = paste0(
      "Wang-Liu-Feng ", tools::toTitleCase(measure),
      " vector-independence test (", component, ")"
    ),
    data.name = paste(x.name, "and", y.name),
    raw.statistic = c(
      maximum.absolute.rank.correlation = core$maximum,
      centered.sum.squares = core$sum_statistic
    ),
    components = components,
    diagnostics = list(
      selected.component = component,
      measure = measure,
      tie.policy = "error; primary continuous-margin calibration",
      sum.calibration = "intrinsic X-row permutation variance",
      permutations = B,
      seed = seed,
      explicit.seed.preserves.caller.RNG.state = !is.null(seed),
      max.calibration = "type-I extreme value; constant 1/sqrt(pi)",
      combination = "Fisher chi-square with 4 degrees of freedom",
      review.only = c(
        "mutual independence across more than two vector blocks",
        "Hoeffding D / Blum-Kiefer-Rosenblatt R / tau-star degenerate kernels"
      )
    ),
    call = call,
    estimate = c(maximum.absolute.rank.correlation = core$maximum),
    parameter = if (component == "fisher") c(df = 4) else NULL
  )
}
