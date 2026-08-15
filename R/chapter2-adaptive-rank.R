#' Zhang--Feng adaptive marginal-rank location tests
#'
#' Implements the one- and two-sample rank tests studied by Zhang and Feng
#' (2024): a maximum of standardized marginal rank scores, the squared-rank
#' sum statistic inherited from Ouyang et al. (2022), and their equal-weight
#' Cauchy combination. Observations are rows and the ordered coordinates are
#' columns.
#'
#' For one sample, let \eqn{U_j} be the sum of the ranks of
#' \eqn{|X_{ij}-\mu_j|} having positive residuals. For two samples,
#' \eqn{U_j} is the usual Mann--Whitney statistic obtained from the pooled
#' ranks of the first sample. If \eqn{e_j} and \eqn{v_j} denote the exact
#' untied-null mean and variance of \eqn{U_j}, the maximum component is
#' \deqn{T_{\max}=\max_j (U_j-e_j)^2/v_j-2\log(p)+\log\log(p).}
#' Its limiting distribution has cdf
#' \eqn{G(y)=\exp\{-\pi^{-1/2}\exp(-y/2)\}}.
#'
#' Write \eqn{M_j=(U_j-e_j)^2}, with exact untied-null mean
#' \eqn{E_0(M_j)} and variance \eqn{\operatorname{var}_0(M_j)}. The sum
#' component is
#' \deqn{T_{\mathrm{sum}}=
#'   \frac{\sqrt p\{p^{-1}\sum_j M_j-E_0(M_j)\}}
#'        {\sqrt{\operatorname{var}_0(M_j)\tau^2}},}
#' and uses an upper-tail standard-normal calibration. The long-run variance
#' \eqn{\tau^2} is never chosen silently. Supply a positive `tau_sq`, or set
#' `tau_method = "ouyang_parzen"` and explicitly supply the lag-window size
#' `lag`. The latter computes
#' \deqn{\widehat\tau^2=1+2\sum_{k=1}^{L-1}w(k/L)\widehat\gamma(k),\qquad
#' \widehat\gamma(k)=\frac{1}{p-k}\sum_{j=1}^{p-k}
#' (R_j-\bar R)(R_{j+k}-\bar R),}
#' where \eqn{R_j=\{M_j-E_0(M_j)\}/
#' \sqrt{\operatorname{var}_0(M_j)}} and \eqn{w} is the Parzen kernel.
#' `lag` is \eqn{L}, so lags 1 through \eqn{L-1} are included. No automatic
#' bandwidth rule, positivity floor, ridge, or other repair is applied.
#'
#' With component p-values \eqn{p_{\max}} and \eqn{p_{\mathrm{sum}}}, the
#' combined test uses
#' \deqn{p_{\mathrm C}=1-F_{\mathrm C}\left\{
#' \tfrac12\cot(\pi p_{\max})+
#' \tfrac12\cot(\pi p_{\mathrm{sum}})\right\},}
#' evaluated with stable log-tail arithmetic. The paper relies on asymptotic
#' independence of the maximum and sum components.
#'
#' Exact zeros or tied absolute residuals are rejected by the one-sample
#' function, and exact pooled ties are rejected by the two-sample function.
#' This is deliberate: the displayed finite-sample moments and asymptotic
#' calibrations are for continuous marginal distributions and Zhang and Feng
#' do not specify a tie correction. The one-sample null additionally assumes
#' coordinatewise symmetry about `mu`. The two-sample construction assumes a
#' common pure-shift model (the two populations otherwise have the same
#' continuous marginal distributions). Both sum calibrations assume that the
#' coordinates, in their supplied order, form a sufficiently weakly dependent
#' (strong-mixing) sequence. Consequently, an Ouyang--Parzen result may change
#' if columns are permuted; column order is part of that estimator's contract.
#'
#' The primary paper contains only the squared-rank sum, maximum, and Cauchy
#' combination above. It does not define the general rank \eqn{L_q} family or
#' a minimum-p adaptive test attributed to it in the accompanying book draft;
#' those procedures are therefore not manufactured here.
#'
#' @param x A numeric matrix or data frame with observations in rows.
#' @param mu For the one-sample test, a finite null location vector with one
#'   entry per column of `x`. `NULL` uses zero.
#' @param y For the two-sample test, a second numeric matrix or data frame with
#'   observations in rows and the same variables as `x`.
#' @param component Which published statistic to return: `"max"` (the
#'   default), `"sum"`, or `"combined"`.
#' @param tau_sq An optional, externally supplied positive long-run variance
#'   for the standardized squared-rank score sequence. It is required for
#'   `"sum"` and `"combined"` unless the explicit Ouyang--Parzen route is
#'   selected. It is not used for `"max"`.
#' @param tau_method `NULL`, `"supplied"`, or `"ouyang_parzen"`. For the
#'   latter, `lag` must be supplied explicitly. There is no default bandwidth.
#' @param lag For `tau_method = "ouyang_parzen"`, the integer lag-window size
#'   \eqn{L} in `1, ..., p`. Lags 1 through \eqn{L-1} are used.
#'
#' @return An object of classes `hd_location_test` and `htest`. Beyond the
#'   standard fields it retains all marginal rank scores, exact null moments,
#'   both component p-values when available, long-run-variance inputs and
#'   autocovariance diagnostics, and the stable Cauchy calculation. For the
#'   combined test, the finite primary `statistic` is the Cauchy angle
#'   \eqn{\arctan(T_{\mathrm C})/\pi}; the formal (possibly infinite)
#'   \eqn{T_{\mathrm C}} is retained in `raw.statistic` and `components`.
#'
#' @references
#' Zhang, J. and Feng, L. (2024). Adaptive rank-based tests for high
#' dimensional mean problems. *Statistics & Probability Letters*, **214**,
#' 110226. \doi{10.1016/j.spl.2024.110226}.
#'
#' Ouyang, M., Xie, Y., and Wang, W. (2022). A new test of high-dimensional
#' mean vector with applications to gene set testing. *Computational
#' Statistics & Data Analysis*, **171**, 107495.
#' \doi{10.1016/j.csda.2022.107495}.
#'
#' @examples
#' x <- matrix(c(
#'   -4, -1,  2,  5,
#'   -2,  3, -5,  1,
#'    1, -4,  3, -6,
#'    3,  5, -1,  2,
#'    6, -2,  7, -3
#' ), nrow = 5, byrow = TRUE)
#' zhang_feng_rank_one_sample_test(x)
#' zhang_feng_rank_one_sample_test(x, component = "sum", tau_sq = 1)
#'
#' y <- x + matrix(rep(c(0.4, -0.2, 0.3, 0.1), each = nrow(x)), ncol = 4)
#' zhang_feng_rank_two_sample_test(x, y)
#'
#' @name zhang_feng_rank_tests
NULL

.zf_component <- function(component) {
  match.arg(component, c("max", "sum", "combined"))
}

.zf_validate_tau_scalar <- function(tau_sq) {
  if (!is.numeric(tau_sq)) {
    stop("`tau_sq` must be one finite, strictly positive number; no " %+%
           "positivity repair is applied.", call. = FALSE)
  }
  tau_sq <- as.numeric(tau_sq)
  if (length(tau_sq) != 1L || is.na(tau_sq) || !is.finite(tau_sq) ||
      tau_sq <= 0) {
    stop("`tau_sq` must be one finite, strictly positive number; no " %+%
           "positivity repair is applied.", call. = FALSE)
  }
  tau_sq
}

.zf_validate_lag <- function(lag, p) {
  if (!is.numeric(lag) || length(lag) != 1L || is.na(lag) ||
      !is.finite(lag) || lag != floor(lag) || lag < 1 || lag > p ||
      lag > .Machine$integer.max) {
    stop(sprintf("`lag` must be an integer in 1, ..., p = %d.", p),
         call. = FALSE)
  }
  as.integer(lag)
}

.zf_resolve_tau <- function(component, tau_sq, tau_method, lag, score) {
  if (identical(component, "max")) {
    if (!is.null(tau_sq) || !is.null(tau_method) || !is.null(lag)) {
      stop("`tau_sq`, `tau_method`, and `lag` are not used when " %+%
             "`component = \"max\"`.", call. = FALSE)
    }
    return(list(value = NULL, source = "not required", details = NULL))
  }

  if (!is.null(tau_method)) {
    if (!is.character(tau_method) || length(tau_method) != 1L ||
        is.na(tau_method)) {
      stop("`tau_method` must be NULL, \"supplied\", or " %+%
             "\"ouyang_parzen\".", call. = FALSE)
    }
    tau_method <- match.arg(tau_method, c("supplied", "ouyang_parzen"))
  }

  if (!is.null(tau_sq)) {
    if (identical(tau_method, "ouyang_parzen")) {
      stop("Choose exactly one long-run-variance route: supply `tau_sq` " %+%
             "or use `tau_method = \"ouyang_parzen\"`.", call. = FALSE)
    }
    if (!is.null(lag)) {
      stop("`lag` is only used with `tau_method = \"ouyang_parzen\"`.",
           call. = FALSE)
    }
    tau_sq <- .zf_validate_tau_scalar(tau_sq)
    return(list(
      value = tau_sq,
      source = "supplied",
      details = list(tau_squared = tau_sq, repair = "none")
    ))
  }

  if (is.null(tau_method)) {
    stop("For a sum or combined test, supply positive `tau_sq` or set " %+%
           "`tau_method = \"ouyang_parzen\"` with an explicit `lag`.",
         call. = FALSE)
  }
  if (identical(tau_method, "supplied")) {
    stop("`tau_method = \"supplied\"` requires a positive `tau_sq`.",
         call. = FALSE)
  }
  if (is.null(lag)) {
    stop("`tau_method = \"ouyang_parzen\"` requires an explicit `lag`; " %+%
           "no automatic bandwidth is used.", call. = FALSE)
  }
  lag <- .zf_validate_lag(lag, length(score))
  details <- cpp_zhang_feng_parzen_tau(score, lag)
  value <- as.numeric(details$tau_squared)
  if (length(value) != 1L || !is.finite(value) || value <= 0) {
    stop(sprintf(
      paste0(
        "The Ouyang--Parzen estimate of `tau_sq` is not strictly positive ",
        "and finite (estimate = %s); no floor, ridge, or absolute-value ",
        "repair is applied."
      ),
      format(value, digits = 17L)
    ), call. = FALSE)
  }
  details$repair <- "none"
  list(value = value, source = "ouyang_parzen", details = details)
}

.zf_gumbel_tail <- function(statistic) {
  log_intensity <- -0.5 * statistic - 0.5 * log(pi)
  log_max <- log(.Machine$double.xmax)
  if (log_intensity > log_max) {
    return(list(
      p.value = 1, log.p.value = 0, log.cdf = -Inf,
      cdf = 0, log.intensity = log_intensity
    ))
  }
  intensity <- exp(log_intensity)
  log_cdf <- -intensity
  if (intensity == 0) {
    p_value <- 0
    log_p_value <- log_intensity
  } else {
    p_value <- -expm1(-intensity)
    log_p_value <- log(p_value)
  }
  list(
    p.value = p_value,
    log.p.value = log_p_value,
    log.cdf = log_cdf,
    cdf = exp(log_cdf),
    log.intensity = log_intensity
  )
}

.zf_normal_tail <- function(statistic) {
  list(
    p.value = stats::pnorm(statistic, lower.tail = FALSE),
    log.p.value = stats::pnorm(
      statistic, lower.tail = FALSE, log.p = TRUE
    ),
    log.cdf = stats::pnorm(statistic, log.p = TRUE),
    cdf = stats::pnorm(statistic)
  )
}

.zf_cotangent_log <- function(p_value, log_p, log_q) {
  # Work with the smaller tail.  This preserves the complement even when the
  # reported probability has rounded to zero or one, while still evaluating
  # the exact cotangent (rather than its small-tail approximation) whenever
  # that tail is representable.
  use_lower_tail <- log_p <= log_q
  log_small <- if (use_lower_tail) log_p else log_q
  small <- exp(log_small)
  if (small == 0) {
    return(c(
      sign = if (use_lower_tail) 1 else -1,
      logabs = -log(pi) - log_small
    ))
  }
  numerator <- cospi(small)
  if (numerator == 0) {
    return(c(sign = 0, logabs = -Inf))
  }
  c(
    sign = if (use_lower_tail) sign(numerator) else -sign(numerator),
    logabs = log(abs(numerator)) - log(sinpi(small))
  )
}

.zf_signed_log_sum <- function(first, second) {
  sign_first <- unname(first[["sign"]])
  sign_second <- unname(second[["sign"]])
  log_first <- unname(first[["logabs"]])
  log_second <- unname(second[["logabs"]])
  if (sign_first == 0) {
    return(c(sign = sign_second, logabs = log_second))
  }
  if (sign_second == 0) {
    return(c(sign = sign_first, logabs = log_first))
  }
  if (sign_first == sign_second) {
    largest <- max(log_first, log_second)
    return(c(
      sign = sign_first,
      logabs = largest + log1p(exp(min(log_first, log_second) - largest))
    ))
  }
  if (log_first == log_second) {
    return(c(sign = 0, logabs = -Inf))
  }
  if (log_first > log_second) {
    return(c(
      sign = sign_first,
      logabs = log_first + log1p(-exp(log_second - log_first))
    ))
  }
  c(
    sign = sign_second,
    logabs = log_second + log1p(-exp(log_first - log_second))
  )
}

.zf_cauchy_combine <- function(max_tail, sum_tail) {
  max_term <- .zf_cotangent_log(
    max_tail$p.value, max_tail$log.p.value, max_tail$log.cdf
  )
  sum_term <- .zf_cotangent_log(
    sum_tail$p.value, sum_tail$log.p.value, sum_tail$log.cdf
  )
  max_term[["logabs"]] <- max_term[["logabs"]] - log(2)
  sum_term[["logabs"]] <- sum_term[["logabs"]] - log(2)
  total <- .zf_signed_log_sum(max_term, sum_term)
  total_sign <- unname(total[["sign"]])
  total_logabs <- unname(total[["logabs"]])
  if (total_sign == 0) {
    angle <- 0
    transform <- 0
  } else if (total_logabs <= log(.Machine$double.xmax)) {
    transform <- total_sign * exp(total_logabs)
    angle <- atan(transform) / pi
  } else {
    transform <- total_sign * Inf
    angle <- total_sign * (0.5 - atan(exp(-total_logabs)) / pi)
  }
  list(
    p.value = 0.5 - angle,
    angle = angle,
    transform = transform,
    transform.sign = total_sign,
    transform.logabs = total_logabs,
    max.term = max_term,
    sum.term = sum_term,
    weights = c(max = 0.5, sum = 0.5)
  )
}

.zf_rank_test_from_scores <- function(scores, component, tau_sq, tau_method,
                                      lag, method, data_name, n, p, call,
                                      variable_names, null_value,
                                      extra_diagnostics = list()) {
  standardized <- as.numeric(scores$standardized_rank_score)
  standardized_squared <- as.numeric(
    scores$standardized_squared_rank_score
  )
  names(standardized) <- names(standardized_squared) <- variable_names

  maximum <- max(abs(standardized))^2 - 2 * log(p) + log(log(p))
  max_tail <- .zf_gumbel_tail(maximum)
  tau <- .zf_resolve_tau(
    component, tau_sq, tau_method, lag, standardized_squared
  )

  sum_statistic <- NULL
  sum_tail <- NULL
  combined <- NULL
  if (!is.null(tau$value)) {
    sum_statistic <- sqrt(p) *
      (as.numeric(scores$mean_squared_rank_score) -
         as.numeric(scores$null_mean_squared_rank_score)) /
      sqrt(as.numeric(scores$null_variance_squared_rank_score) * tau$value)
    if (!is.finite(sum_statistic)) {
      stop("The squared-rank sum statistic is non-finite.", call. = FALSE)
    }
    sum_tail <- .zf_normal_tail(sum_statistic)
    if (identical(component, "combined")) {
      combined <- .zf_cauchy_combine(max_tail, sum_tail)
    }
  }

  if (identical(component, "max")) {
    statistic <- c(`T[max]` = maximum)
    p_value <- max_tail$p.value
    family <- "Zhang--Feng Gumbel limit"
  } else if (identical(component, "sum")) {
    statistic <- c(`T[sum]` = sum_statistic)
    p_value <- sum_tail$p.value
    family <- "standard normal upper tail"
  } else {
    statistic <- c(`Cauchy angle` = combined$angle)
    p_value <- combined$p.value
    family <- "equal-weight Cauchy combination"
  }

  rank_sum <- as.numeric(scores$rank_sum)
  centered <- as.numeric(scores$centered_rank_score)
  squared <- as.numeric(scores$squared_rank_score)
  names(rank_sum) <- names(centered) <- names(squared) <- variable_names
  variance <- c(
    marginal.rank.score = as.numeric(scores$null_variance_rank_sum),
    squared.rank.score =
      as.numeric(scores$null_variance_squared_rank_score)
  )
  if (!is.null(tau$value)) {
    variance <- c(variance, long.run.standardized.square = tau$value)
  }
  raw_statistic <- c(maximum = maximum)
  if (!is.null(sum_statistic)) {
    raw_statistic <- c(raw_statistic, sum = sum_statistic)
  }
  if (!is.null(combined)) {
    raw_statistic <- c(
      raw_statistic,
      cauchy.angle = combined$angle,
      cauchy.transform = combined$transform
    )
  }

  assumptions <- c(
    "continuous marginal distributions (no ties)",
    "ordered-coordinate strong-mixing conditions",
    "high-dimensional asymptotic calibration"
  )
  .new_hd_location_test(
    statistic = statistic,
    parameter = c(dimension = p),
    p.value = p_value,
    method = method,
    data.name = data_name,
    alternative = "two.sided",
    raw.statistic = raw_statistic,
    estimate = standardized,
    null.value = null_value,
    null.distribution = list(
      family = family,
      asymptotic = TRUE,
      upper.tail = TRUE,
      gumbel.cdf = "exp(-pi^(-1/2) * exp(-y/2))",
      assumptions = assumptions
    ),
    variance = variance,
    components = list(
      rank.sum = rank_sum,
      centered.rank.score = centered,
      standardized.rank.score = standardized,
      squared.rank.score = squared,
      standardized.squared.rank.score = standardized_squared,
      mean.squared.rank.score =
        as.numeric(scores$mean_squared_rank_score),
      null.moments = c(
        mean.rank.sum = as.numeric(scores$null_mean_rank_sum),
        variance.rank.sum = as.numeric(scores$null_variance_rank_sum),
        mean.squared.rank.score =
          as.numeric(scores$null_mean_squared_rank_score),
        variance.squared.rank.score =
          as.numeric(scores$null_variance_squared_rank_score)
      ),
      maximum = list(
        statistic = maximum,
        p.value = max_tail$p.value,
        log.p.value = max_tail$log.p.value,
        cdf = max_tail$cdf,
        log.cdf = max_tail$log.cdf,
        log.intensity = max_tail$log.intensity
      ),
      sum = if (is.null(sum_tail)) NULL else list(
        statistic = sum_statistic,
        p.value = sum_tail$p.value,
        log.p.value = sum_tail$log.p.value,
        cdf = sum_tail$cdf,
        log.cdf = sum_tail$log.cdf
      ),
      combined = combined,
      long.run.variance = tau$details
    ),
    diagnostics = c(list(
      selected.component = component,
      tau.source = tau$source,
      tau.squared = tau$value,
      lag.window.size = if (is.null(tau$details)) NULL else
        tau$details$lag_window_size,
      variance.repair = "none",
      tie.correction = "none; exact ties are rejected",
      coordinate.order.sensitive = identical(tau$source, "ouyang_parzen"),
      coordinate.order.contract = paste(
        "columns are an ordered strong-mixing sequence; the explicit",
        "Ouyang--Parzen estimate is not permutation invariant"
      ),
      automatic.bandwidth = FALSE,
      calibration = "asymptotic",
      published.components.only = c("max", "squared-rank sum", "Cauchy"),
      book.general.Lq.minimum.p.implemented = FALSE
    ), extra_diagnostics),
    n = n,
    p = p,
    call = call
  )
}

#' @rdname zhang_feng_rank_tests
#' @export
zhang_feng_rank_one_sample_test <- function(
    x, mu = NULL, component = c("max", "sum", "combined"),
    tau_sq = NULL, tau_method = NULL, lag = NULL) {
  call <- match.call()
  data_name <- deparse1(substitute(x))
  x <- .as_data_matrix(x, "x", min_rows = 2L)
  p <- ncol(x)
  if (p < 2L) {
    stop("The Zhang--Feng maximum calibration requires at least two " %+%
           "variables (`p >= 2`).", call. = FALSE)
  }
  component <- .zf_component(component)
  if (is.null(mu)) {
    mu <- rep.int(0, p)
  } else {
    mu <- .as_location(mu, p, "mu")
  }
  variable_names <- .hotelling_variable_names(x)
  names(mu) <- variable_names
  scores <- cpp_zhang_feng_one_sample_scores(x, mu)
  .zf_rank_test_from_scores(
    scores = scores,
    component = component,
    tau_sq = tau_sq,
    tau_method = tau_method,
    lag = lag,
    method = paste0(
      "Zhang--Feng one-sample adaptive marginal-rank test (", component, ")"
    ),
    data_name = data_name,
    n = c(x = nrow(x)),
    p = p,
    call = call,
    variable_names = variable_names,
    null_value = mu,
    extra_diagnostics = list(
      marginal.null = "coordinatewise symmetry about mu",
      overflow.fallback.columns =
        as.integer(scores$overflow_fallback_columns)
    )
  )
}

#' @rdname zhang_feng_rank_tests
#' @export
zhang_feng_rank_two_sample_test <- function(
    x, y, component = c("max", "sum", "combined"),
    tau_sq = NULL, tau_method = NULL, lag = NULL) {
  call <- match.call()
  x_name <- deparse1(substitute(x))
  y_name <- deparse1(substitute(y))
  x <- .as_data_matrix(x, "x", min_rows = 2L)
  y <- .as_data_matrix(y, "y", min_rows = 2L)
  .check_two_sample_variables(x, y)
  p <- ncol(x)
  if (p < 2L) {
    stop("The Zhang--Feng maximum calibration requires at least two " %+%
           "variables (`p >= 2`).", call. = FALSE)
  }
  component <- .zf_component(component)
  variable_names <- .hotelling_variable_names(x)
  scores <- cpp_zhang_feng_two_sample_scores(x, y)
  .zf_rank_test_from_scores(
    scores = scores,
    component = component,
    tau_sq = tau_sq,
    tau_method = tau_method,
    lag = lag,
    method = paste0(
      "Zhang--Feng two-sample adaptive marginal-rank test (", component, ")"
    ),
    data_name = paste(x_name, "and", y_name),
    n = c(x = nrow(x), y = nrow(y)),
    p = p,
    call = call,
    variable_names = variable_names,
    null_value = stats::setNames(rep.int(0, p), variable_names),
    extra_diagnostics = list(
      marginal.null = "common continuous pure-shift model"
    )
  )
}
