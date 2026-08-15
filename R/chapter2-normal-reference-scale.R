#' Zhang--Zhu--Zhang normal-reference scale-invariant two-sample test
#'
#' Tests equality of two high-dimensional mean vectors when the covariance
#' matrices may differ.  Observations are rows and variables are columns.  If
#' \eqn{S_1} and \eqn{S_2} are the unbiased sample covariance matrices, let
#' \deqn{\widehat\Omega_n=\frac{n_2}{n}S_1+\frac{n_1}{n}S_2,
#' \qquad \widehat D_n=\operatorname{diag}(\widehat\Omega_n).}
#' Notice the crossed sample-size weights.  The scale-invariant statistic is
#' \deqn{T_{n,p}=\frac{n_1n_2}{np}
#' (\bar X_1-\bar X_2)^{\mathsf T}\widehat D_n^{-1}
#' (\bar X_1-\bar X_2).}
#'
#' This is the scale-invariant test of Zhang, Zhu, and Zhang (2023).  It is
#' not the raw-\eqn{L_2} normal-reference test based on
#' \eqn{\|\bar X_1-\bar X_2\|^2}, and it is not the later
#' normal-reference F-type test.  Its feasible reference distribution is the
#' one-parameter Welch--Satterthwaite approximation
#' \eqn{\chi^2_d/d}.
#'
#' Put
#' \eqn{\widehat R_i=\widehat D_n^{-1/2}S_i\widehat D_n^{-1/2}}.  The
#' group-specific bias-corrected squared traces are
#' \deqn{\widehat q_i=
#' \frac{(n_i-1)^2}{(n_i-2)(n_i+1)}\left\{
#' \operatorname{tr}(\widehat R_i^2)-
#' \frac{\operatorname{tr}^2(\widehat R_i)}{n_i-1}\right\}.}
#' The combined estimate and unadjusted degrees of freedom are
#' \deqn{\widehat q=
#' \frac{n_2^2}{n^2}\widehat q_1+
#' \frac{n_1^2}{n^2}\widehat q_2+
#' \frac{2n_1n_2}{n^2}\operatorname{tr}(\widehat R_1\widehat R_2),
#' \qquad \widehat d=\frac{p^2}{\widehat q}.}
#' These corrections apply to the trace estimate; the statistic
#' \eqn{T_{n,p}} itself is not bias-subtracted.
#'
#' With `df_adjustment = "paper"`, the empirical finite-sample rule in the
#' paper is also used.  Define
#' \deqn{c_{n,p}=1+
#' \operatorname{tr}(\widehat R_n^2)/p^{3/2},\qquad
#' \widehat R_n=\frac{n_2}{n}\widehat R_1+
#' \frac{n_1}{n}\widehat R_2.}
#' If \eqn{c_{n,p}\leq 1.2}, the reference degrees of freedom are
#' \eqn{\widehat d/c_{n,p}}; otherwise they remain \eqn{\widehat d}.  The raw
#' squared trace in this rule is deliberately not replaced by the
#' bias-corrected \eqn{\widehat q}.  `df_adjustment = "none"` always uses
#' \eqn{\widehat d}.  Both calibrations are returned regardless of the
#' selected option.
#'
#' The Gaussian oracle mixture has cumulants
#' \eqn{1}, \eqn{2p^{-2}\operatorname{tr}(R_n^2)}, and
#' \eqn{8p^{-3}\operatorname{tr}(R_n^3)}.  The third cumulant and the
#' theoretical quantity
#' \eqn{d^*=\operatorname{tr}^3(R_n^2)/
#' \operatorname{tr}^2(R_n^3)} are not estimated by the paper's feasible
#' calibration and therefore are not substituted here.
#'
#' The C++ kernel translates and rescales each coordinate by a common rule
#' across both groups before computing moments.  It uses a streamed primal
#' covariance identity when \eqn{p\leq n_1+n_2} and an observation-level dual
#' Gram identity otherwise; neither route constructs a persistent
#' \eqn{p\times p} matrix.  No ridge, pseudoinverse, absolute-value repair,
#' trace floor, or degrees-of-freedom clamp is used.  A nonpositive feasible
#' trace estimate is an explicit calibration failure.
#'
#' @param x,y Numeric matrices or data frames with observations in rows and
#'   the same variables in columns.  Each sample must contain at least three
#'   observations.
#' @param alpha Significance level for the upper-tail rejection rule.
#' @param df_adjustment Either `"paper"` for the published conditional
#'   finite-sample degrees-of-freedom adjustment or `"none"` for the
#'   unadjusted Welch--Satterthwaite degrees of freedom.
#'
#' @return An object of class `c("hd_location_test", "htest")`.  The
#'   statistic is \eqn{T_{n,p}} and the p-value is the upper tail of the
#'   selected \eqn{\chi^2_d/d} reference.  `components` contains both
#'   adjusted and unadjusted degrees of freedom, p-values, critical values
#'   and rejection decisions; raw and corrected trace components; diagonal
#'   standardisation; and coordinate contributions.  `diagnostics` records
#'   the crossed weights, trace formulas, primal/dual route, empirical
#'   threshold decision, and no-repair contract.  `log.p.value` in each
#'   calibration remains available when the ordinary tail underflows.
#'
#' @references
#' Zhang, L., Zhu, T., and Zhang, J.-T. (2023). Two-sample
#' Behrens--Fisher problems for high-dimensional data: a normal reference
#' scale-invariant test. *Journal of Applied Statistics*, **50**, 456--476.
#' \doi{10.1080/02664763.2020.1834516}.
#'
#' @examples
#' x <- matrix(c(
#'   0.2, -0.4, 1.1, 0.7,
#'   1.0,  0.3, 0.2, -0.8,
#'  -0.6,  1.2, 0.5, 0.1,
#'   0.4, -0.7, 1.3, 0.9
#' ), ncol = 4, byrow = TRUE)
#' y <- matrix(c(
#'  -0.1, 0.5,  0.4, 1.0,
#'   0.8, 1.1, -0.5, 0.2,
#'  -0.7, 0.2,  1.0, 0.6,
#'   0.3, 0.9,  0.1, -0.4,
#'   1.1, -0.3, 0.8, 0.5
#' ), ncol = 4, byrow = TRUE)
#' zhang_zhu_zhang_two_sample_test(x, y)
#'
#' @export
zhang_zhu_zhang_two_sample_test <- function(
    x, y, alpha = 0.05, df_adjustment = c("paper", "none")) {
  call <- match.call()
  x.name <- deparse1(substitute(x))
  y.name <- deparse1(substitute(y))
  x <- .as_data_matrix(x, "x", min_rows = 3L)
  y <- .as_data_matrix(y, "y", min_rows = 3L)
  .check_two_sample_variables(x, y)

  alpha <- as.numeric(alpha)
  if (length(alpha) != 1L || is.na(alpha) || !is.finite(alpha) ||
      alpha <= 0 || alpha >= 1) {
    stop("`alpha` must be one finite number strictly between zero and one.",
         call. = FALSE)
  }
  df_adjustment <- match.arg(df_adjustment)

  fit <- cpp_zhang_zhu_zhang_two_sample(x, y)
  statistic <- as.numeric(fit$statistic)
  unadjusted <- .zzz23_chisq_calibration(
    statistic, as.numeric(fit$df_unadjusted), alpha
  )
  paper <- .zzz23_chisq_calibration(
    statistic, as.numeric(fit$df_paper), alpha
  )
  selected <- if (identical(df_adjustment, "paper")) paper else unadjusted

  n1 <- nrow(x)
  n2 <- nrow(y)
  p <- ncol(x)
  variable.names <- .hotelling_variable_names(x)
  name.vector <- function(value) {
    stats::setNames(as.numeric(value), variable.names)
  }

  mean.x <- name.vector(fit$mean1)
  mean.y <- name.vector(fit$mean2)
  difference <- name.vector(fit$mean_difference_input)
  null.value <- stats::setNames(numeric(p), variable.names)
  selected.name <- if (identical(df_adjustment, "paper")) {
    "paper-threshold-adjusted"
  } else {
    "unadjusted"
  }

  calibration <- list(
    selected = df_adjustment,
    unadjusted = unadjusted,
    paper = paper
  )

  .new_hd_location_test(
    statistic = c(T.NRSI = statistic),
    parameter = c(df = selected$df),
    p.value = selected$p.value,
    method = paste(
      "Zhang-Zhu-Zhang normal-reference scale-invariant two-sample test",
      "(Welch-Satterthwaite scaled chi-square; not F-type)"
    ),
    data.name = paste(x.name, "and", y.name),
    alternative = "two.sided",
    raw.statistic = c(T.NRSI = statistic),
    estimate = difference,
    null.value = null.value,
    null.distribution = list(
      family = "scaled chi-square",
      parameters = c(df = selected$df, scale = 1 / selected$df),
      exact = FALSE,
      tail = "upper",
      selected.calibration = selected.name,
      assumptions = paste(
        "Zhang-Zhu-Zhang high-dimensional normal-reference assumptions;",
        "independent samples; positive crossed pooled marginal variances;",
        "positive feasible bias-corrected squared-trace estimate"
      )
    ),
    variance = c(
      selected.reference = 2 / selected$df,
      unadjusted.reference = 2 / unadjusted$df,
      paper.reference = 2 / paper$df
    ),
    components = list(
      mean.x = mean.x,
      mean.y = mean.y,
      difference = difference,
      mean.x.standardized = name.vector(fit$mean1_standardized),
      mean.y.standardized = name.vector(fit$mean2_standardized),
      difference.standardized = name.vector(
        fit$mean_difference_standardized
      ),
      difference.input.representable = stats::setNames(
        as.logical(fit$mean_difference_input_representable),
        variable.names
      ),
      T.NRSI = statistic,
      coordinate.contribution = name.vector(
        fit$coordinate_contribution
      ),
      covariance.weights = c(
        group1 = as.numeric(fit$weight_group1_covariance),
        group2 = as.numeric(fit$weight_group2_covariance)
      ),
      variance1.diagonal.standardized = name.vector(
        fit$variance1_diagonal_standardized
      ),
      variance2.diagonal.standardized = name.vector(
        fit$variance2_diagonal_standardized
      ),
      D.hat.diagonal.standardized = name.vector(
        fit$D_hat_diagonal_standardized
      ),
      D.hat.diagonal.input = name.vector(
        fit$D_hat_diagonal_input
      ),
      D.hat.diagonal.input.canonical = name.vector(
        fit$D_hat_diagonal_input_canonical
      ),
      log.D.hat.diagonal.input = name.vector(
        fit$log_D_hat_diagonal_input
      ),
      trace.R1 = as.numeric(fit$trace_R1),
      trace.R2 = as.numeric(fit$trace_R2),
      trace.Rn = as.numeric(fit$trace_Rn),
      trace.R1.squared.raw = as.numeric(
        fit$trace_R1_squared_raw
      ),
      trace.R2.squared.raw = as.numeric(
        fit$trace_R2_squared_raw
      ),
      trace.R1.R2.raw = as.numeric(fit$trace_R1_R2_raw),
      trace.Rn.squared.raw = as.numeric(
        fit$trace_Rn_squared_raw
      ),
      trace.R1.squared.bracket = as.numeric(
        fit$trace_R1_squared_bracket
      ),
      trace.R2.squared.bracket = as.numeric(
        fit$trace_R2_squared_bracket
      ),
      trace.R1.squared.corrected = as.numeric(
        fit$trace_R1_squared_corrected
      ),
      trace.R2.squared.corrected = as.numeric(
        fit$trace_R2_squared_corrected
      ),
      trace.Rn.squared.corrected = as.numeric(
        fit$trace_Rn_squared_corrected
      ),
      trace.correction.factor = c(
        group1 = as.numeric(fit$trace_correction_factor1),
        group2 = as.numeric(fit$trace_correction_factor2)
      ),
      df.hat = unadjusted$df,
      c.np = as.numeric(fit$c_np),
      paper.correction.applied = isTRUE(
        fit$paper_correction_applied
      ),
      df.paper = paper$df,
      df.used = selected$df,
      p.value.unadjusted = unadjusted$p.value,
      log.p.value.unadjusted = unadjusted$log.p.value,
      critical.value.unadjusted = unadjusted$critical.value,
      reject.unadjusted = unadjusted$reject,
      p.value.paper = paper$p.value,
      log.p.value.paper = paper$log.p.value,
      critical.value.paper = paper$critical.value,
      reject.paper = paper$reject,
      calibration = calibration,
      feasible.reference.cumulants = c(
        mean = 1,
        variance.unadjusted = 2 / unadjusted$df,
        third = NA_real_
      ),
      theoretical.third.cumulant.estimated = FALSE,
      n1 = n1,
      n2 = n2,
      p = p
    ),
    diagnostics = list(
      statistic = paste(
        "diagonally standardized nonnegative quadratic statistic;",
        "no statistic bias subtraction"
      ),
      raw.L2.normal.reference = FALSE,
      F.type.reference = FALSE,
      reference = "Welch-Satterthwaite chi-square_df / df",
      tail = "upper",
      df.adjustment = df_adjustment,
      paper.df.rule = paste(
        "df / c_np when c_np <= 1.2, otherwise df;",
        "the statistic is unchanged"
      ),
      paper.correction.applied = isTRUE(
        fit$paper_correction_applied
      ),
      rejection = list(
        alpha = alpha,
        selected = df_adjustment,
        critical.value = selected$critical.value,
        reject = selected$reject,
        unadjusted = unadjusted$reject,
        paper = paper$reject
      ),
      covariance.weights = paste0(
        "group 1 uses n2/n = ",
        format(as.numeric(fit$weight_group1_covariance), digits = 16),
        "; group 2 uses n1/n = ",
        format(as.numeric(fit$weight_group2_covariance), digits = 16)
      ),
      sample.variance.denominators = c(
        group1 = n1 - 1,
        group2 = n2 - 1
      ),
      trace.correction = paste(
        "Eq. (24) whole-bracket multipliers",
        "(n_i-1)^2 / ((n_i-2)(n_i+1))"
      ),
      combined.trace = paste(
        "Eq. (25) crossed squared weights and",
        "2*n1*n2/n^2 cross term"
      ),
      c.np.trace = "raw trace(Rn_hat^2), not bias-corrected qhat",
      trace.computation = as.character(fit$trace_computation),
      trace.route.rule = "primal if p <= n1+n2; otherwise dual",
      constructs.p.by.p.matrix = isTRUE(fit$constructs_p_by_p_matrix),
      minimum.D.hat.standardized = as.numeric(
        fit$minimum_D_hat_standardized
      ),
      maximum.diagonal.identity.error = as.numeric(
        fit$maximum_diagonal_identity_error
      ),
      zero.group.variance.coordinates = c(
        group1 = sum(as.logical(fit$zero_variance_group1)),
        group2 = sum(as.logical(fit$zero_variance_group2))
      ),
      internal.scaling = paste(
        "common per-coordinate anchor and scale across both groups;",
        "long-double compensated moments and Gram sums"
      ),
      column.anchor = name.vector(fit$column_anchor),
      column.scale = name.vector(fit$column_scale),
      column.log.scale = name.vector(fit$column_log_scale),
      subtraction.overflow.fallback.columns = sum(
        as.logical(fit$subtraction_overflow_fallback)
      ),
      regularization = "none",
      pseudoinverse = "none",
      absolute.value.repair = "none",
      trace.floor = "none",
      degrees.of.freedom.clamp = "none"
    ),
    n = c(group1 = n1, group2 = n2),
    p = p,
    call = call
  )
}


.zzz23_chisq_calibration <- function(statistic, df, alpha) {
  statistic <- as.numeric(statistic)
  df <- as.numeric(df)
  if (length(statistic) != 1L || !is.finite(statistic) ||
      statistic < 0 || length(df) != 1L || !is.finite(df) || df <= 0) {
    stop(
      "The scaled chi-square calibration requires finite T >= 0 and df > 0.",
      call. = FALSE
    )
  }

  argument <- df * statistic
  p.value <- stats::pchisq(argument, df = df, lower.tail = FALSE)
  log.p.value <- stats::pchisq(
    argument, df = df, lower.tail = FALSE, log.p = TRUE
  )
  critical.value <- stats::qchisq(1 - alpha, df = df) / df
  if (!is.finite(critical.value) || critical.value < 0 ||
      is.na(p.value) || is.na(log.p.value)) {
    stop(
      "The scaled chi-square tail or critical value is not numerically " %+%
        "defined for the feasible degrees of freedom.",
      call. = FALSE
    )
  }
  list(
    df = df,
    scaled.statistic = argument,
    p.value = p.value,
    log.p.value = log.p.value,
    critical.value = critical.value,
    reject = isTRUE(statistic > critical.value)
  )
}
