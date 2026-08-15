.ch4wn_validate_matrix <- function(x, name = "x") {
  if (is.data.frame(x)) x <- as.matrix(x)
  if (is.atomic(x) && is.null(dim(x))) x <- matrix(x, ncol = 1L)
  if (!is.matrix(x) || !is.numeric(x)) {
    stop(sprintf("`%s` must be a numeric matrix (rows are time points).", name),
         call. = FALSE)
  }
  storage.mode(x) <- "double"
  if (nrow(x) < 3L || ncol(x) < 1L) {
    stop(sprintf("`%s` must have at least three rows and one column.", name),
         call. = FALSE)
  }
  if (anyNA(x) || any(!is.finite(x))) {
    stop(sprintf("`%s` must contain only finite values.", name),
         call. = FALSE)
  }
  x
}

.ch4wn_validate_lag <- function(lag, n) {
  if (!is.numeric(lag) || length(lag) != 1L || is.na(lag) ||
      !is.finite(lag) || lag != floor(lag)) {
    stop("`lag` must be one finite integer.", call. = FALSE)
  }
  if (lag < 1 || lag > n - 2L) {
    stop("`lag` must be between 1 and n - 2.", call. = FALSE)
  }
  as.integer(lag)
}

.ch4wn_validate_logical <- function(value, name) {
  if (!is.logical(value) || length(value) != 1L || is.na(value)) {
    stop(sprintf("`%s` must be TRUE or FALSE.", name), call. = FALSE)
  }
  value
}

.ch4wn_variable_names <- function(x) {
  if (is.null(colnames(x))) paste0("variable", seq_len(ncol(x))) else colnames(x)
}

.ch4wn_stable_residuals <- function(x, center) {
  center <- match.arg(center, c("none", "mean"))
  data.scale <- max(abs(x))
  if (!is.finite(data.scale) || data.scale <= 0) {
    stop("The data matrix is identically zero; no calibration is defined.",
         call. = FALSE)
  }
  scaled <- x / data.scale
  if (center == "mean") {
    center.scaled <- colMeans(scaled)
    residuals <- sweep(scaled, 2L, center.scaled, FUN = "-")
  } else {
    center.scaled <- numeric(ncol(x))
    residuals <- scaled
  }
  if (anyNA(residuals) || any(!is.finite(residuals))) {
    stop("Stable common-scale preprocessing produced a non-finite residual.",
         call. = FALSE)
  }
  list(
    residuals = residuals,
    center = center,
    center.scaled = center.scaled,
    center.original = center.scaled * data.scale,
    data.scale = data.scale,
    fourth.power.log.scale = 4 * log(data.scale)
  )
}

.ch4wn_extreme_tail <- function(statistic) {
  log.lambda <- -0.5 * statistic - 0.5 * log(pi)
  log.maximum <- log(.Machine$double.xmax)
  log.minimum <- log(.Machine$double.xmin)
  if (log.lambda >= log.maximum) {
    return(list(p.value = 1, log.p.value = 0, cdf = 0,
                log.lambda = log.lambda))
  }
  if (log.lambda <= log.minimum) {
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

.ch4wn_new_test <- function(statistic, p.value, method, data.name,
                            raw.statistic, components, diagnostics,
                            parameter = NULL, estimate = NULL,
                            call = NULL) {
  answer <- list(
    statistic = statistic,
    p.value = as.numeric(p.value),
    alternative = "greater",
    method = method,
    data.name = data.name,
    raw.statistic = raw.statistic,
    null.value = c(`serial dependence` = 0),
    components = components,
    diagnostics = diagnostics,
    call = call
  )
  if (!is.null(parameter)) answer$parameter <- parameter
  if (!is.null(estimate)) answer$estimate <- estimate
  class(answer) <- c("hd_white_noise_test", "htest")
  answer
}

#' Classical Box--Pierce or Ljung--Box white-noise test
#'
#' Computes the univariate portmanteau statistic written in Chapter 4. This
#' function is intentionally univariate; the high-dimensional procedures are
#' provided by the other functions in this file.
#'
#' @param x Numeric vector, ordered in time.
#' @param lag Positive truncation lag. It must be smaller than `length(x)`.
#' @param type Either `"Ljung-Box"` or `"Box-Pierce"`.
#' @param center Whether to subtract the sample mean before computing sample
#'   autocorrelations.
#'
#' @return An object inheriting from `htest`. The `components` field contains
#'   every lagged sample autocorrelation and its contribution.
#' @references
#' Box, G. E. P. and Pierce, D. A. (1970). Distribution of residual
#' autocorrelations in autoregressive-integrated moving average time series
#' models. *Journal of the American Statistical Association*, 65, 1509--1526.
#' Ljung, G. M. and Box, G. E. P. (1978). On a measure of lack of fit in time
#' series models. *Biometrika*, 65, 297--303.
#' @examples
#' x <- c(0.2, -0.1, 0.3, 0.05, -0.2, 0.1, 0.4, -0.3)
#' white_noise_portmanteau_test(x, lag = 2)
#' @export
white_noise_portmanteau_test <- function(
    x, lag = 1L, type = c("Ljung-Box", "Box-Pierce"), center = TRUE) {
  call <- match.call()
  data.name <- deparse(substitute(x))
  type <- match.arg(type)
  center <- .ch4wn_validate_logical(center, "center")
  if (!is.numeric(x) || !is.null(dim(x)) || length(x) < 2L ||
      anyNA(x) || any(!is.finite(x))) {
    stop("`x` must be a finite numeric vector with at least two values.",
         call. = FALSE)
  }
  n <- length(x)
  if (!is.numeric(lag) || length(lag) != 1L || is.na(lag) ||
      !is.finite(lag) || lag != floor(lag) || lag < 1 || lag >= n) {
    stop("`lag` must be one integer between 1 and length(x) - 1.",
         call. = FALSE)
  }
  lag <- as.integer(lag)
  data.scale <- max(abs(x))
  if (!is.finite(data.scale) || data.scale <= 0) {
    stop("The series is identically zero; autocorrelations are undefined.",
         call. = FALSE)
  }
  residuals <- x / data.scale
  fitted.center.scaled <- 0
  if (center) {
    fitted.center.scaled <- mean(residuals)
    residuals <- residuals - fitted.center.scaled
  }
  denominator <- sum(residuals^2)
  if (!is.finite(denominator) || denominator <= 0) {
    stop("The centered series has zero or non-finite sum of squares.",
         call. = FALSE)
  }
  autocorrelation <- vapply(seq_len(lag), function(h) {
    sum(residuals[seq_len(n - h)] * residuals[h + seq_len(n - h)]) /
      denominator
  }, numeric(1L))
  names(autocorrelation) <- paste0("lag", seq_len(lag))
  if (type == "Box-Pierce") {
    contributions <- n * autocorrelation^2
  } else {
    contributions <- n * (n + 2) * autocorrelation^2 /
      (n - seq_len(lag))
  }
  statistic <- sum(contributions)
  p.value <- stats::pchisq(statistic, df = lag, lower.tail = FALSE)
  .ch4wn_new_test(
    statistic = stats::setNames(statistic, type),
    p.value = p.value,
    method = paste(type, "test for white noise"),
    data.name = data.name,
    raw.statistic = statistic,
    parameter = stats::setNames(lag, "df"),
    components = list(
      autocorrelations = autocorrelation,
      contributions = contributions,
      sample.size = n,
      truncation.lag = lag,
      fitted.center = fitted.center.scaled * data.scale,
      denominator.scaled = denominator
    ),
    diagnostics = list(
      asymptotic.calibration = "chi-square with lag degrees of freedom",
      centered = center,
      common.data.scale = data.scale,
      no.repair = TRUE
    ),
    call = call
  )
}

#' Feng--Liu--Ma high-dimensional white-noise test
#'
#' Implements the feasible max statistic, the diagonal-deleted U-statistic
#' sum test, and their Fisher combination from Feng, Liu and Ma. The primary
#' paper uses `n` (not `n - h`) in every lagged sample covariance and uses the
#' ordered-pair denominator `n * (n - 1)` in both the sum statistic and the
#' feasible trace estimate.
#'
#' @param x Numeric matrix with time points in rows and coordinates in columns.
#' @param lag Positive lag truncation level, no larger than `n - 2`.
#' @param component Which calibrated result supplies the top-level statistic
#'   and p-value: `"fisher"`, `"sum"`, or `"max"`. All three are returned in
#'   `components` regardless of this choice.
#' @param center `"none"` reproduces the mean-zero primary definition.
#'   `"mean"` subtracts column sample means as an explicit preprocessing step;
#'   the paper's null theorem does not account for estimated means.
#' @param keep_lag Whether to retain lag-specific maxima and sum numerators.
#'
#' @return An `htest` object with raw feasible components and diagnostics.
#' @references Feng, L., Liu, B. and Ma, Y. Testing for High-Dimensional
#'   White Noise. Statistica Sinica. \doi{10.5705/ss.202023.0300}.
#' @examples
#' t <- seq_len(18)
#' x <- cbind(sin(t), cos(t / 2), sin(t / 3 + 0.2))
#' feng_liu_ma_white_noise_test(x, lag = 2)
#' @export
feng_liu_ma_white_noise_test <- function(
    x, lag = 1L, component = c("fisher", "sum", "max"),
    center = c("none", "mean"), keep_lag = FALSE) {
  call <- match.call()
  data.name <- deparse(substitute(x))
  x <- .ch4wn_validate_matrix(x)
  n <- nrow(x)
  p <- ncol(x)
  lag <- .ch4wn_validate_lag(lag, n)
  component <- match.arg(component)
  center <- match.arg(center)
  keep_lag <- .ch4wn_validate_logical(keep_lag, "keep_lag")
  comparison.count <- as.double(lag) * as.double(p)^2
  if (!is.finite(comparison.count) || comparison.count <= 1) {
    stop(paste0(
      "The extreme-value calibration requires lag * p^2 > 1; use ",
      "`white_noise_portmanteau_test()` for the single univariate lag."
    ), call. = FALSE)
  }

  prep <- .ch4wn_stable_residuals(x, center)
  core <- cpp_ch4wn_flm_core(prep$residuals, lag, keep_lag)
  core$marginal_variances_scaled <-
    as.numeric(core$marginal_variances_scaled)
  core$lag_ordered_pair_count <- as.numeric(core$lag_ordered_pair_count)
  if (!is.null(core$lag_maximum_absolute_correlation)) {
    core$lag_maximum_absolute_correlation <-
      as.numeric(core$lag_maximum_absolute_correlation)
  }
  if (!is.null(core$lag_sum_numerator_scaled)) {
    core$lag_sum_numerator_scaled <-
      as.numeric(core$lag_sum_numerator_scaled)
  }
  names(core$marginal_variances_scaled) <- .ch4wn_variable_names(x)
  names(core$maximum_indices) <- c("row.coordinate", "column.coordinate")
  maximum.names <- .ch4wn_variable_names(x)[core$maximum_indices]

  max.centered <- core$maximum_statistic^2 -
    2 * log(comparison.count) + log(log(comparison.count))
  max.tail <- .ch4wn_extreme_tail(max.centered)

  sum.variance <- 2 * lag / (n * (n - 1)) *
    core$trace_estimate_scaled^2
  if (!is.finite(sum.variance) || sum.variance <= 0) {
    stop("The feasible FLM sum variance is not finite and positive.",
         call. = FALSE)
  }
  sum.standard.error <- sqrt(sum.variance)
  sum.z <- core$sum_statistic_scaled / sum.standard.error
  if (!is.finite(sum.z)) {
    stop("The standardized FLM sum statistic is non-finite.", call. = FALSE)
  }
  sum.log.p <- stats::pnorm(sum.z, lower.tail = FALSE, log.p = TRUE)
  sum.p <- exp(sum.log.p)
  fisher <- -2 * (max.tail$log.p.value + sum.log.p)
  fisher.p <- stats::pchisq(fisher, df = 4, lower.tail = FALSE)

  components <- list(
    T.MAX = core$maximum_statistic,
    maximum.absolute.correlation = core$maximum_absolute_correlation,
    signed.correlation.at.maximum = core$signed_correlation_at_maximum,
    maximum.location = c(
      lag = core$maximum_lag,
      row.coordinate = core$maximum_indices[[1L]],
      column.coordinate = core$maximum_indices[[2L]]
    ),
    maximum.coordinate.names = maximum.names,
    maximum.Gumbel.statistic = max.centered,
    p.maximum = max.tail$p.value,
    log.p.maximum = max.tail$log.p.value,
    T.SUM.scaled = core$sum_statistic_scaled,
    sum.numerator.scaled = core$sum_numerator_scaled,
    trace.Sigma.squared.scaled = core$trace_estimate_scaled,
    trace.numerator.scaled = core$trace_numerator_scaled,
    primary.ordered.denominator = core$primary_ordered_denominator,
    sigma.S.squared.scaled = sum.variance,
    sigma.S.scaled = sum.standard.error,
    sum.z = sum.z,
    p.sum = sum.p,
    log.p.sum = sum.log.p,
    Fisher.statistic = fisher,
    p.Fisher = fisher.p,
    marginal.second.moments.scaled = core$marginal_variances_scaled,
    lag.ordered.pair.count = core$lag_ordered_pair_count,
    lag.maximum.absolute.correlation =
      core$lag_maximum_absolute_correlation,
    lag.maximum.indices = core$lag_maximum_indices,
    lag.sum.numerator.scaled = core$lag_sum_numerator_scaled,
    truncation.lag = lag,
    sample.size = n,
    dimension = p
  )
  if (component == "fisher") {
    statistic <- stats::setNames(fisher, "Fisher combination")
    p.value <- fisher.p
    raw <- fisher
    method <- "Feng--Liu--Ma Fisher-combined high-dimensional white-noise test"
  } else if (component == "sum") {
    statistic <- stats::setNames(sum.z, "sum z")
    p.value <- sum.p
    raw <- core$sum_statistic_scaled
    method <- "Feng--Liu--Ma sum-type high-dimensional white-noise test"
  } else {
    statistic <- stats::setNames(max.centered, "maximum Gumbel score")
    p.value <- max.tail$p.value
    raw <- core$maximum_statistic
    method <- "Feng--Liu--Ma max-type high-dimensional white-noise test"
  }
  .ch4wn_new_test(
    statistic = statistic,
    p.value = p.value,
    method = method,
    data.name = data.name,
    raw.statistic = raw,
    components = components,
    diagnostics = list(
      selected.component = component,
      centering = center,
      primary.mean.zero.preprocessing = center == "none",
      common.data.scale = prep$data.scale,
      fourth.power.log.scale = prep$fourth.power.log.scale,
      stable.scaled.raw.components = TRUE,
      sample.autocovariance.divisor = "n at every lag",
      sum.and.trace.denominator = "ordered n * (n - 1)",
      sum.effective.indices = "t,s = 1,...,n-h with t != s",
      maximum.calibration =
        "G(y) = exp{-pi^(-1/2) exp(-y/2)}",
      sum.calibration = "upper standard-normal tail",
      Fisher.calibration = "upper chi-square tail with 4 df",
      asymptotic.applicability = paste(
        "High-dimensional primary-paper conditions are not diagnosed from",
        "one realised data matrix."
      ),
      no.ridge.floor.absolute.value.or.deletion = TRUE
    ),
    call = call
  )
}

#' Zhao--Chen--Wang spatial-sign white-noise test
#'
#' Implements the spatial-sign sum statistic and feasible variance in Zhao,
#' Chen and Wang. The primary null variance is `(H / 2) * tr(Omega^2)^2`.
#'
#' @param x Numeric matrix with time points in rows.
#' @param lag Positive lag truncation level, no larger than `n - 2`.
#' @param center `"none"` is the centered-at-zero primary definition;
#'   `"mean"` is an explicit sample-mean preprocessing extension.
#' @param zero_action `"error"` rejects a zero residual direction.
#'   `"keep"` uses the paper's convention `U(0) = 0` and records the count.
#' @param keep_signs Whether to retain the fitted spatial-sign matrix.
#'
#' @return An `htest` object with the raw statistic, trace estimate, variance,
#'   lag contributions, and zero-direction diagnostics.
#' @references Zhao, P., Chen, D. and Wang, Z. Spatial-sign-based
#'   high-dimensional white noises test.
#'   \doi{10.1080/24754269.2024.2363715}.
#' @examples
#' t <- seq_len(18)
#' x <- cbind(sin(t), cos(t / 2), sin(t / 3 + 0.2))
#' zhao_chen_wang_spatial_sign_white_noise_test(x, lag = 2)
#' @export
zhao_chen_wang_spatial_sign_white_noise_test <- function(
    x, lag = 1L, center = c("none", "mean"),
    zero_action = c("error", "keep"), keep_signs = FALSE) {
  call <- match.call()
  data.name <- deparse(substitute(x))
  x <- .ch4wn_validate_matrix(x)
  n <- nrow(x)
  p <- ncol(x)
  lag <- .ch4wn_validate_lag(lag, n)
  center <- match.arg(center)
  zero_action <- match.arg(zero_action)
  keep_signs <- .ch4wn_validate_logical(keep_signs, "keep_signs")

  prep <- .ch4wn_stable_residuals(x, center)
  core <- cpp_ch4wn_spatial_sign_core(
    prep$residuals, lag, keep_signs
  )
  core$lag_numerators <- as.numeric(core$lag_numerators)
  core$lag_statistics <- as.numeric(core$lag_statistics)
  core$lag_unordered_pair_count <-
    as.numeric(core$lag_unordered_pair_count)
  if (core$zero_sign_count > 0L && zero_action == "error") {
    stop(sprintf(
      paste0(
        "%d residual direction(s) are exactly zero. Set `zero_action = ",
        "\"keep\"` to apply the paper's explicit U(0) = 0 convention."
      ),
      core$zero_sign_count
    ), call. = FALSE)
  }
  if (!is.null(core$signs)) {
    dimnames(core$signs) <- list(
      if (is.null(rownames(x))) paste0("time", seq_len(n)) else rownames(x),
      .ch4wn_variable_names(x)
    )
  }
  variance <- lag / 2 * core$trace_estimate^2
  if (!is.finite(variance) || variance <= 0) {
    stop("The feasible spatial-sign variance is not finite and positive.",
         call. = FALSE)
  }
  standard.error <- sqrt(variance)
  z <- core$statistic / standard.error
  if (!is.finite(z)) {
    stop("The standardized spatial-sign statistic is non-finite.",
         call. = FALSE)
  }
  log.p <- stats::pnorm(z, lower.tail = FALSE, log.p = TRUE)
  p.value <- exp(log.p)
  .ch4wn_new_test(
    statistic = stats::setNames(z, "spatial-sign z"),
    p.value = p.value,
    method = "Zhao--Chen--Wang spatial-sign high-dimensional white-noise test",
    data.name = data.name,
    raw.statistic = core$statistic,
    components = list(
      T.S = core$statistic,
      trace.Omega.squared = core$trace_estimate,
      trace.unordered.sum = core$trace_unordered_sum,
      trace.denominator = core$trace_denominator,
      sigma.S.squared = variance,
      sigma.S = standard.error,
      z = z,
      log.p.value = log.p,
      lag.numerators = core$lag_numerators,
      lag.statistics = core$lag_statistics,
      lag.unordered.pair.count = core$lag_unordered_pair_count,
      signs = core$signs,
      fitted.center = prep$center.original,
      truncation.lag = lag,
      sample.size = n,
      dimension = p
    ),
    diagnostics = list(
      centering = center,
      primary.mean.zero.preprocessing = center == "none",
      common.data.scale = prep$data.scale,
      zero.action = zero_action,
      zero.sign.count = core$zero_sign_count,
      zero.sign.convention = "U(0) = 0",
      primary.variance.factor = "H / 2",
      calibration = "upper standard-normal tail",
      partial.zero.calibration.warning = core$zero_sign_count > 0L,
      asymptotic.applicability = paste(
        "Elliptical null and high-dimensional trace conditions are not",
        "diagnosed from one realised data matrix."
      ),
      no.ridge.floor.absolute.value.or.deletion = TRUE
    ),
    call = call
  )
}

#' Chen--Song--Feng rank-based max white-noise test
#'
#' Implements the two rank procedures for which the primary paper supplies
#' direct, scalable statistics and complete extreme-value calibrations:
#' Spearman's rho and Kendall's tau. The paper does not establish rank-based
#' sum or adaptive Fisher tests; those constructions in the book draft are not
#' used here.
#'
#' @param x Numeric matrix with time points in rows.
#' @param lag Positive lag truncation level, no larger than `n - 2`.
#' @param measure Either `"spearman"` or `"kendall"`.
#' @param keep_lag Whether to retain lag-specific maxima and their locations.
#'
#' @return An `htest` object. Exact ties are rejected because the published
#'   null variance and distribution-free calibration assume continuous margins.
#' @references Chen, D., Song, F. and Feng, L. Rank Based Tests for High
#'   Dimensional White Noise. Statistica Sinica 35, 1323--1347.
#'   \doi{10.5705/ss.202022.0382}.
#' @examples
#' t <- seq_len(18)
#' x <- cbind(sin(t), cos(t / 2), sin(t / 3 + 0.2))
#' chen_song_feng_rank_white_noise_test(x, lag = 2, measure = "spearman")
#' @export
chen_song_feng_rank_white_noise_test <- function(
    x, lag = 1L, measure = c("spearman", "kendall"),
    keep_lag = FALSE) {
  call <- match.call()
  data.name <- deparse(substitute(x))
  x <- .ch4wn_validate_matrix(x)
  n <- nrow(x)
  p <- ncol(x)
  lag <- .ch4wn_validate_lag(lag, n)
  measure <- match.arg(measure)
  keep_lag <- .ch4wn_validate_logical(keep_lag, "keep_lag")
  comparison.count <- as.double(lag) * as.double(p)^2
  if (!is.finite(comparison.count) || comparison.count <= 1) {
    stop(paste0(
      "The extreme-value calibration requires lag * p^2 > 1; use ",
      "`white_noise_portmanteau_test()` for the single univariate lag."
    ), call. = FALSE)
  }
  variable.names <- .ch4wn_variable_names(x)
  duplicated.column <- vapply(
    seq_len(p), function(j) anyDuplicated(x[, j]) > 0L, logical(1L)
  )
  if (any(duplicated.column)) {
    stop(sprintf(
      paste0(
        "The primary continuous-margin rank calibration does not allow ",
        "ties; exact ties occur in: %s."
      ),
      paste(variable.names[duplicated.column], collapse = ", ")
    ), call. = FALSE)
  }

  core <- cpp_ch4wn_rank_max_core(x, lag, measure, keep_lag)
  if (!is.null(core$lag_maximum_standardized_square)) {
    core$lag_maximum_standardized_square <-
      as.numeric(core$lag_maximum_standardized_square)
    core$lag_measure_at_maximum <-
      as.numeric(core$lag_measure_at_maximum)
  }
  names(core$maximum_indices) <- c("row.coordinate", "column.coordinate")
  maximum.names <- variable.names[core$maximum_indices]
  transformed <- core$maximum_standardized_square -
    2 * log(comparison.count) + log(log(comparison.count))
  tail <- .ch4wn_extreme_tail(transformed)
  maximum.effective.n <- n - core$maximum_lag
  if (measure == "spearman") {
    null.variance <- 1 / (maximum.effective.n - 1)
    scaling <- "(n - k) * rho_ij(k)^2"
    label <- "Spearman"
  } else {
    null.variance <- 2 * (2 * maximum.effective.n + 5) /
      (9 * maximum.effective.n * (maximum.effective.n - 1))
    scaling <- paste0(
      "9 * (n-k) * (n-k-1) / {2 * [2 * (n-k) + 5]} * tau_ij(k)^2"
    )
    label <- "Kendall"
  }
  .ch4wn_new_test(
    statistic = stats::setNames(transformed,
                                paste(label, "Gumbel score")),
    p.value = tail$p.value,
    method = paste0(
      "Chen--Song--Feng ", label,
      " rank-max high-dimensional white-noise test"
    ),
    data.name = data.name,
    raw.statistic = core$maximum_standardized_square,
    estimate = stats::setNames(core$measure_at_maximum,
                               tolower(label)),
    components = list(
      maximum.standardized.square = core$maximum_standardized_square,
      measure.at.maximum = core$measure_at_maximum,
      maximum.location = c(
        lag = core$maximum_lag,
        row.coordinate = core$maximum_indices[[1L]],
        column.coordinate = core$maximum_indices[[2L]]
      ),
      maximum.coordinate.names = maximum.names,
      maximum.effective.sample.size = maximum.effective.n,
      primary.null.variance.at.maximum = null.variance,
      primary.scaling = scaling,
      Gumbel.statistic = transformed,
      log.p.value = tail$log.p.value,
      comparison.count = comparison.count,
      lag.maximum.standardized.square =
        core$lag_maximum_standardized_square,
      lag.measure.at.maximum = core$lag_measure_at_maximum,
      lag.maximum.indices = core$lag_maximum_indices,
      truncation.lag = lag,
      sample.size = n,
      dimension = p
    ),
    diagnostics = list(
      measure = measure,
      ties.detected = core$ties_detected,
      continuous.margin.calibration = TRUE,
      calibration = "G(y) = exp{-pi^(-1/2) exp(-y/2)}",
      rank.sum.or.adaptive.test.implemented = FALSE,
      rank.sum.or.adaptive.reason = paste(
        "The primary paper explicitly lists rank sum tests and max/sum",
        "combination theory as future work."
      ),
      review.only.degenerate.methods = c(
        "Hoeffding D", "Blum-Kiefer-Rosenblatt R",
        "Bergsma-Dassios-Yanagimoto tau-star"
      ),
      review.only.reason = paste(
        "The paper gives literal order-5/order-6/order-4 U kernels, but no",
        "legally available scalable author implementation was found; an",
        "exponential or high-polynomial interface is not exposed as practical."
      ),
      asymptotic.applicability = paste(
        "The weak-dependence and dimensional growth conditions are not",
        "diagnosed from one realised data matrix."
      ),
      no.jitter.midrank.or.tie.repair = TRUE
    ),
    call = call
  )
}
