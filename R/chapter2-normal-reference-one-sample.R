.zzg22nr_validate_alpha <- function(alpha) {
  alpha <- as.numeric(alpha)
  if (length(alpha) != 1L || is.na(alpha) || !is.finite(alpha) ||
      alpha <= 0 || alpha >= 1) {
    stop("`alpha` must be one finite number strictly between zero and one.",
         call. = FALSE)
  }
  alpha
}


.zzg22nr_restore_power <- function(value, log.scale, power) {
  value <- as.numeric(value)
  log.absolute <- rep.int(-Inf, length(value))
  nonzero <- value != 0
  log.absolute[nonzero] <- log(abs(value[nonzero])) + power * log.scale

  log.maximum <- log(.Machine$double.xmax)
  log.minimum <- log(.Machine$double.xmin) + log(.Machine$double.eps)
  representable <- !nonzero |
    (is.finite(log.absolute) & log.absolute <= log.maximum &
       log.absolute >= log.minimum)
  restored <- rep.int(NA_real_, length(value))
  restored[!nonzero] <- 0
  eligible <- nonzero & representable
  restored[eligible] <- sign(value[eligible]) * exp(log.absolute[eligible])

  list(
    value = restored,
    representable = representable,
    log.absolute = log.absolute
  )
}


.zzg22nr_safe_mean <- function(x) {
  anchor <- max(abs(x))
  if (anchor == 0) {
    return(list(value = numeric(ncol(x)), representable = rep(TRUE, ncol(x))))
  }
  restored <- .zzg22nr_restore_power(
    colMeans(x / anchor), log(anchor), 1
  )
  list(value = restored$value, representable = restored$representable)
}


.zzg22nr_prepare_one_sample <- function(x, mu) {
  anchor <- max(max(abs(x)), max(abs(mu)))
  if (anchor == 0) {
    anchor <- 1
  }
  direct <- suppressWarnings(
    sweep(x, 2L, mu, check.margin = FALSE)
  )
  residual <- sweep(
    x / anchor, 2L, mu / anchor, check.margin = FALSE
  )
  estimate <- .zzg22nr_safe_mean(x)
  list(
    residual = residual,
    log.scale = log(anchor),
    estimate = estimate$value,
    estimate.representable = estimate$representable,
    subtraction.overflow.fallback = any(!is.finite(direct)),
    transformation = "rowwise x_i - mu"
  )
}


.zzg22nr_prepare_paired <- function(x, y, delta) {
  anchor <- max(max(abs(x)), max(abs(y)), max(abs(delta)))
  if (anchor == 0) {
    anchor <- 1
  }
  direct <- suppressWarnings(
    sweep(x - y, 2L, delta, check.margin = FALSE)
  )
  x.unit <- x / anchor
  y.unit <- y / anchor
  residual <- sweep(
    x.unit - y.unit, 2L, delta / anchor, check.margin = FALSE
  )
  estimate <- .zzg22nr_restore_power(
    colMeans(x.unit) - colMeans(y.unit), log(anchor), 1
  )
  list(
    residual = residual,
    log.scale = log(anchor),
    estimate = estimate$value,
    estimate.representable = estimate$representable,
    subtraction.overflow.fallback = any(!is.finite(direct)),
    transformation = "rowwise paired difference x_i - y_i - delta"
  )
}


.zzg22nr_as_contrast_matrix <- function(L, p) {
  if (is.data.frame(L)) {
    L <- data.matrix(L)
  }
  if (is.vector(L) && is.numeric(L)) {
    L <- matrix(L, nrow = 1L)
  }
  if (!is.matrix(L) || !is.numeric(L) || nrow(L) < 1L || ncol(L) != p) {
    stop(sprintf(
      "`L` must be a numeric matrix with at least one row and %d columns.",
      p
    ), call. = FALSE)
  }
  storage.mode(L) <- "double"
  if (anyNA(L) || any(!is.finite(L))) {
    stop("`L` must contain only finite values.", call. = FALSE)
  }
  L
}


.zzg22nr_prepare_linear <- function(x, L, rhs) {
  x.anchor <- max(abs(x))
  L.anchor <- max(abs(L))
  rhs.anchor <- max(abs(rhs))

  if (x.anchor > 0 && L.anchor > 0) {
    log.product.scale <- log(x.anchor) + log(L.anchor)
    product.unit <- (x / x.anchor) %*% t(L / L.anchor)
  } else {
    log.product.scale <- -Inf
    product.unit <- matrix(0, nrow(x), nrow(L))
  }
  log.rhs.scale <- if (rhs.anchor > 0) log(rhs.anchor) else -Inf
  log.scale <- max(log.product.scale, log.rhs.scale)
  if (!is.finite(log.scale)) {
    log.scale <- 0
  }
  product.weight <- if (is.finite(log.product.scale)) {
    exp(log.product.scale - log.scale)
  } else {
    0
  }
  rhs.weight <- if (rhs.anchor > 0) {
    exp(log.rhs.scale - log.scale)
  } else {
    0
  }
  rhs.unit <- if (rhs.anchor > 0) rhs / rhs.anchor else numeric(length(rhs))
  residual <- sweep(
    product.weight * product.unit,
    2L,
    rhs.weight * rhs.unit,
    check.margin = FALSE
  )
  if (anyNA(residual) || any(!is.finite(residual))) {
    stop("The rowwise transform `L x_i - rhs` produced non-finite values.",
         call. = FALSE)
  }

  mean.product.scaled <- product.weight * colMeans(product.unit)
  estimate <- .zzg22nr_restore_power(
    mean.product.scaled, log.scale, 1
  )
  direct <- suppressWarnings(x %*% t(L))
  direct <- suppressWarnings(
    sweep(direct, 2L, rhs, check.margin = FALSE)
  )

  list(
    residual = residual,
    log.scale = log.scale,
    estimate = estimate$value,
    estimate.representable = estimate$representable,
    subtraction.overflow.fallback = any(!is.finite(direct)),
    transformation = "same-iid-unit rowwise transform z_i = L x_i - rhs"
  )
}


.zzg22nr_fit <- function(prepared, alpha, estimate, null.value,
                         variable.names, method, data.name, n, call) {
  fit <- cpp_zhang_zhou_guo_one_sample(prepared$residual)
  total.log.scale <- prepared$log.scale + log(as.numeric(fit$internal_scale))

  restore <- function(value, power) {
    .zzg22nr_restore_power(value, total.log.scale, power)
  }
  T.input <- restore(fit$T_scaled, 2)
  oracle.input <- restore(fit$oracle_Q_scaled, 2)
  trace.S.input <- restore(fit$trace_S_scaled, 2)
  trace.S2.input <- restore(fit$trace_S2_scaled, 4)
  trace.S3.input <- restore(fit$trace_S3_scaled, 6)
  a1.input <- restore(fit$a1_hat_scaled, 2)
  a2.input <- restore(fit$a2_hat_scaled, 4)
  a3.input <- restore(fit$a3_hat_scaled, 6)
  kappa2.input <- restore(fit$kappa2_scaled, 4)
  kappa3.input <- restore(fit$kappa3_scaled, 6)
  beta0.input <- restore(fit$beta0_scaled, 2)
  beta1.input <- restore(fit$beta1_scaled, 2)
  difference.input <- restore(fit$mean_residual_scaled, 1)
  physical.scale <- .zzg22nr_restore_power(1, total.log.scale, 1)

  degrees <- as.numeric(fit$df)
  reference.argument <- as.numeric(fit$reference_argument)
  p.value <- stats::pchisq(
    reference.argument, degrees, lower.tail = FALSE
  )
  log.p.value <- stats::pchisq(
    reference.argument, degrees, lower.tail = FALSE, log.p = TRUE
  )
  critical.reference <- stats::qchisq(
    alpha, degrees, lower.tail = FALSE
  )
  if (!is.finite(critical.reference)) {
    stop("`alpha` gives a non-finite chi-square critical value.", call. = FALSE)
  }
  critical.T.scaled <- as.numeric(fit$beta0_scaled) +
    as.numeric(fit$beta1_scaled) * critical.reference
  critical.T.input <- restore(critical.T.scaled, 2)
  reject <- isTRUE(reference.argument > critical.reference)

  estimate <- stats::setNames(as.numeric(estimate), variable.names)
  null.value <- stats::setNames(as.numeric(null.value), variable.names)
  difference <- stats::setNames(
    as.numeric(difference.input$value), variable.names
  )
  estimate.representable <- stats::setNames(
    as.logical(prepared$estimate.representable), variable.names
  )
  difference.representable <- stats::setNames(
    as.logical(difference.input$representable), variable.names
  )

  .new_hd_location_test(
    statistic = c(chisq.argument = reference.argument),
    parameter = c(df = degrees),
    p.value = p.value,
    method = method,
    data.name = data.name,
    alternative = "two.sided",
    raw.statistic = c(
      T.centered.input = unname(T.input$value),
      T.centered.scaled = as.numeric(fit$T_scaled)
    ),
    estimate = estimate,
    null.value = null.value,
    null.distribution = list(
      family = "shifted-scaled chi-square normal reference",
      parameters = c(
        beta0.scaled = as.numeric(fit$beta0_scaled),
        beta1.scaled = as.numeric(fit$beta1_scaled),
        df = degrees
      ),
      exact = FALSE,
      tail = "upper",
      statistic.transformation = "(T - beta0) / beta1",
      assumptions = paste(
        "independent identically distributed rows; Zhang-Zhou-Guo",
        "high-dimensional conditions; Gaussian normal-reference third",
        "cumulant; positive unbiased trace estimates"
      )
    ),
    variance = c(
      T.scaled = as.numeric(fit$kappa2_scaled),
      T.input = unname(kappa2.input$value)
    ),
    components = list(
      estimate = estimate,
      null.value = null.value,
      mean.null.contrast = difference,
      estimate.input.representable = estimate.representable,
      mean.null.contrast.input.representable = difference.representable,
      oracle.Q.n.scaled = as.numeric(fit$oracle_Q_scaled),
      oracle.Q.n.input = unname(oracle.input$value),
      T.centered.scaled = as.numeric(fit$T_scaled),
      T.centered.input = unname(T.input$value),
      T.by.oracle.minus.trace.scaled =
        as.numeric(fit$T_by_oracle_minus_trace_scaled),
      T.identity.error.scaled = as.numeric(fit$T_identity_error_scaled),
      trace.S.scaled = as.numeric(fit$trace_S_scaled),
      trace.S.input = unname(trace.S.input$value),
      trace.S2.scaled = as.numeric(fit$trace_S2_scaled),
      trace.S2.input = unname(trace.S2.input$value),
      trace.S3.scaled = as.numeric(fit$trace_S3_scaled),
      trace.S3.input = unname(trace.S3.input$value),
      a1.hat.scaled = as.numeric(fit$a1_hat_scaled),
      a1.hat.input = unname(a1.input$value),
      a2.bracket.scaled = as.numeric(fit$a2_bracket_scaled),
      a2.factor = as.numeric(fit$a2_factor),
      a2.hat.scaled = as.numeric(fit$a2_hat_scaled),
      a2.hat.input = unname(a2.input$value),
      a3.bracket.scaled = as.numeric(fit$a3_bracket_scaled),
      a3.factor = as.numeric(fit$a3_factor),
      a3.hat.scaled = as.numeric(fit$a3_hat_scaled),
      a3.hat.input = unname(a3.input$value),
      feasible.reference.cumulants.scaled = c(
        kappa1 = as.numeric(fit$kappa1_scaled),
        kappa2 = as.numeric(fit$kappa2_scaled),
        kappa3 = as.numeric(fit$kappa3_scaled)
      ),
      feasible.reference.cumulants.input = c(
        kappa1 = 0,
        kappa2 = unname(kappa2.input$value),
        kappa3 = unname(kappa3.input$value)
      ),
      motivating.oracle.cumulants.scaled = c(
        kappa1 = as.numeric(fit$a1_hat_scaled),
        kappa2 = 2 * as.numeric(fit$a2_hat_scaled),
        kappa3 = 8 * as.numeric(fit$a3_hat_scaled)
      ),
      beta0.scaled = as.numeric(fit$beta0_scaled),
      beta0.input = unname(beta0.input$value),
      beta1.scaled = as.numeric(fit$beta1_scaled),
      beta1.input = unname(beta1.input$value),
      df = degrees,
      chi.square.argument = reference.argument,
      p.value = p.value,
      log.p.value = log.p.value,
      alpha = alpha,
      critical.chi.square = critical.reference,
      critical.T.scaled = critical.T.scaled,
      critical.T.input = unname(critical.T.input$value),
      reject = reject,
      residual.scale.input = unname(physical.scale$value),
      log.residual.scale.input = total.log.scale,
      n = as.numeric(fit$n),
      p = as.numeric(fit$p)
    ),
    diagnostics = list(
      transformation = prepared$transformation,
      statistic.formula =
        "T = n ||mean(z)||^2 - tr(S) = 2 sum_{i<j} z_i' z_j / (n-1)",
      motivating.oracle.Q.used.as.test.statistic = FALSE,
      centred.feasible.statistic = TRUE,
      trace.a2.formula =
        "v^2 {tr(S^2) - tr(S)^2/v} / {(v-1)(v+2)}",
      trace.a3.formula = paste0(
        "v^4 {tr(S^3) - 3 tr(S)tr(S^2)/v + 2 tr(S)^3/v^2}",
        " / {(v-1)(v+4)(v^2-4)}"
      ),
      trace.a3.denominator.uses.v.plus.4 = TRUE,
      normal.reference = TRUE,
      row.permutation.invariant = TRUE,
      orthogonal.coordinate.invariant = TRUE,
      common.scalar.invariant = TRUE,
      coordinatewise.scale.invariant = FALSE,
      nonnormal.third.cumulant.Upsilon.estimated = FALSE,
      nonnormal.third.cumulant.note = paste(
        "Outside the Gaussian normal reference, the finite-sample third",
        "cumulant also contains 4 n Upsilon/(n-1)^2; it is not estimated."
      ),
      trace.computation = as.character(fit$trace_computation),
      constructs.p.by.p.matrix =
        isTRUE(fit$constructs_p_by_p_matrix),
      pretransform.log.scale = prepared$log.scale,
      internal.residual.scale = as.numeric(fit$internal_scale),
      log.residual.scale.input = total.log.scale,
      residual.scale.input.representable =
        isTRUE(physical.scale$representable),
      subtraction.overflow.fallback =
        isTRUE(prepared$subtraction.overflow.fallback),
      input.representability = list(
        T = isTRUE(T.input$representable),
        oracle.Q = isTRUE(oracle.input$representable),
        trace.S = isTRUE(trace.S.input$representable),
        trace.S2 = isTRUE(trace.S2.input$representable),
        trace.S3 = isTRUE(trace.S3.input$representable),
        a1 = isTRUE(a1.input$representable),
        a2 = isTRUE(a2.input$representable),
        a3 = isTRUE(a3.input$representable),
        beta0 = isTRUE(beta0.input$representable),
        beta1 = isTRUE(beta1.input$representable),
        critical.T = isTRUE(critical.T.input$representable)
      ),
      rejection = list(
        alpha = alpha,
        direction = "upper",
        critical.chi.square = critical.reference,
        reject = reject
      ),
      minimum.sample.size = 4L,
      regularization = "none",
      trace.repair = "none",
      cumulant.repair = "none",
      degrees.of.freedom.clamp = "none"
    ),
    n = n,
    p = length(variable.names),
    call = call
  )
}


#' Zhang--Zhou--Guo normal-reference mean tests
#'
#' Implements the feasible one-sample normal-reference test of Zhang, Zhou,
#' and Guo (2022), plus paired and same-unit linear-hypothesis wrappers.
#' Observations are rows and variables are columns.  For null-centred rows
#' \eqn{Z_i}, their mean \eqn{\bar Z}, and their unbiased sample covariance
#' \eqn{S}, the statistic is
#' \deqn{T=n\|\bar Z\|^2-\operatorname{tr}(S)
#'       =\frac{2}{n-1}\sum_{i<j}Z_i^\mathsf{T}Z_j.}
#' This centring is essential.  The quantity \eqn{n\|\bar Z\|^2} by itself is
#' only the motivating Gaussian oracle quadratic form; it is not the paper's
#' feasible test statistic.  This corrects the conflation in the corresponding
#' short passage of the book draft.
#'
#' Write \eqn{v=n-1}, \eqn{t_r=\operatorname{tr}(S^r)}, and
#' \eqn{a_r=\operatorname{tr}(\Sigma^r)}.  The exact Gaussian/Wishart unbiased
#' trace estimators used here are
#' \deqn{\widehat a_1=t_1,\qquad
#' \widehat a_2=\frac{v^2}{(v-1)(v+2)}
#' \left(t_2-\frac{t_1^2}{v}\right),}
#' \deqn{\widehat a_3=
#' \frac{v^4}{(v-1)(v+4)(v^2-4)}
#' \left(t_3-\frac{3t_1t_2}{v}+\frac{2t_1^3}{v^2}\right).}
#' In particular, the denominator is \eqn{v+4}, as in Appendix (A.33) of the
#' primary paper, not \eqn{v+3}.
#'
#' Under the Gaussian normal reference, the first three cumulants of the
#' centred statistic are
#' \deqn{\kappa_1=0,\qquad
#' \kappa_2=\frac{2n}{n-1}a_2,\qquad
#' \kappa_3=\frac{8n(n-2)}{(n-1)^2}a_3.}
#' Matching these to \eqn{\beta_0+\beta_1\chi_d^2} gives
#' \deqn{\widehat\beta_0=-\frac{n}{n-2}
#' \frac{\widehat a_2^2}{\widehat a_3},\qquad
#' \widehat\beta_1=\frac{n-2}{n-1}
#' \frac{\widehat a_3}{\widehat a_2},\qquad
#' \widehat d=\frac{n(n-1)}{(n-2)^2}
#' \frac{\widehat a_2^3}{\widehat a_3^2}.}
#' The reported `statistic` is the chi-square argument
#' \eqn{(T-\widehat\beta_0)/\widehat\beta_1}; the centred \eqn{T} is retained
#' in `raw.statistic` and `components`.  The p-value is the upper tail of
#' \eqn{\chi^2_{\widehat d}}.
#'
#' The normal-reference approximation deliberately does not estimate the
#' additional non-Gaussian third-cumulant term involving
#' \eqn{\Upsilon=\operatorname{E}\{(Z_1^\mathsf{T}Z_2)^3\}}.  Consequently it
#' is not claimed to be a finite-sample exact calibration outside Gaussian
#' data.  At least four rows are needed.  Non-positive
#' \eqn{\widehat a_2} or \eqn{\widehat a_3} is a calibration failure: no ridge,
#' absolute value, floor, or degrees-of-freedom clamp is applied.
#'
#' `zhang_zhou_guo_paired_test()` applies the one-sample method to the paired
#' rows \eqn{X_i-Y_i-\delta}; the two matrices must therefore describe the same
#' observational units in the same row order.  It is not an independent
#' two-sample test.  `zhang_zhou_guo_linear_hypothesis_test()` tests
#' \eqn{H_0:L\mu=\mathrm{rhs}} by the same-unit row transform
#' \eqn{Z_i=L X_i-\mathrm{rhs}}.  It is not an independent-group MANOVA or a
#' between-group general linear hypothesis procedure.  The Euclidean metric
#' after transformation is intentional, so changing `L` by a non-orthogonal
#' row transformation generally changes the test.
#'
#' The test is invariant to row permutations, orthogonal coordinate changes,
#' and a common nonzero scalar change of units (with the null transformed in
#' the same way).  It is not invariant to arbitrary coordinatewise rescaling.
#'
#' A common safe scale is removed before the compiled calculation.  The
#' chi-square argument, degrees of freedom, and p-value are invariant to that
#' scale.  Input-unit versions of \eqn{T} and the traces are returned when they
#' are representable; otherwise their scaled versions and log scale remain
#' available.  The compiled kernel chooses a primal covariance calculation
#' when the transformed dimension does not exceed \eqn{n}, and an observation-
#' level dual Gram calculation otherwise.
#'
#' @param x A numeric matrix or data frame with observations in rows.  At
#'   least four rows are required.
#' @param mu A finite null mean vector with one value per column of `x`.
#'   `NULL` uses the zero vector.
#' @param alpha A finite significance level strictly between zero and one.
#' @param y For the paired wrapper, a numeric matrix or data frame with exactly
#'   the same dimensions as `x`; row \eqn{i} must be paired with row \eqn{i}.
#' @param delta A finite null paired mean-difference vector.  `NULL` uses zero.
#' @param L A finite numeric contrast matrix with one column per column of
#'   `x`.  A numeric vector is treated as a one-row matrix.  Rows define the
#'   transformed coordinates.
#' @param rhs A finite vector with one value per row of `L`, giving the right-
#'   hand side of \eqn{L\mu=\mathrm{rhs}}.  `NULL` uses zero.
#'
#' @return An object of class `c("hd_location_test", "htest")`.  Besides the
#'   standard fields, `components` contains the centred statistic, motivating
#'   oracle quantity, raw and unbiased trace quantities, all three cumulants,
#'   matched chi-square parameters, critical values, and input-scale
#'   representability flags.  `diagnostics` records the precise formulas,
#'   primal/dual route, wrapper transformation, normal-reference limitation,
#'   safe scaling, and no-repair contract.
#'
#' @references
#' Zhang, J.-T., Zhou, B., and Guo, J. (2022). Testing high-dimensional mean
#' vector with applications. *Statistical Papers*, **63**, 1105--1137.
#' \doi{10.1007/s00362-021-01270-z}.
#'
#' @examples
#' x <- rbind(
#'   c(-1.2, 0.4, 1.1), c(0.3, -0.7, 0.5), c(1.4, 0.8, -0.2),
#'   c(-0.6, 1.5, 0.9), c(0.9, -1.1, 1.3), c(1.7, 0.2, -0.8)
#' )
#' zhang_zhou_guo_one_sample_test(x)
#'
#' y <- x + rbind(
#'   c(0.1, -0.2, 0.3), c(-0.2, 0.1, -0.1), c(0.3, 0.2, -0.2),
#'   c(-0.1, -0.3, 0.2), c(0.2, 0.3, 0.1), c(-0.3, -0.1, -0.3)
#' )
#' zhang_zhou_guo_paired_test(x, y)
#' zhang_zhou_guo_linear_hypothesis_test(x, rbind(c(1, -1, 0)))
#'
#' @name zhang_zhou_guo_tests
NULL


#' @rdname zhang_zhou_guo_tests
#' @export
zhang_zhou_guo_one_sample_test <- function(x, mu = NULL, alpha = 0.05) {
  call <- match.call()
  data.name <- deparse1(substitute(x))
  x <- .as_data_matrix(x, "x", min_rows = 4L)
  p <- ncol(x)
  if (is.null(mu)) {
    mu <- numeric(p)
  } else {
    mu <- .as_location(mu, p, "mu")
  }
  alpha <- .zzg22nr_validate_alpha(alpha)
  variable.names <- .hotelling_variable_names(x)
  prepared <- .zzg22nr_prepare_one_sample(x, mu)

  .zzg22nr_fit(
    prepared = prepared,
    alpha = alpha,
    estimate = prepared$estimate,
    null.value = mu,
    variable.names = variable.names,
    method = paste(
      "Zhang-Zhou-Guo one-sample normal-reference mean test",
      "(centred U-statistic; three-cumulant chi-square)"
    ),
    data.name = data.name,
    n = c(x = nrow(x)),
    call = call
  )
}


#' @rdname zhang_zhou_guo_tests
#' @export
zhang_zhou_guo_paired_test <- function(x, y, delta = NULL, alpha = 0.05) {
  call <- match.call()
  x.name <- deparse1(substitute(x))
  y.name <- deparse1(substitute(y))
  x <- .as_data_matrix(x, "x", min_rows = 4L)
  y <- .as_data_matrix(y, "y", min_rows = 4L)
  .check_two_sample_variables(x, y)
  if (nrow(x) != nrow(y)) {
    stop("Paired Zhang-Zhou-Guo data must have the same number of rows.",
         call. = FALSE)
  }
  p <- ncol(x)
  if (is.null(delta)) {
    delta <- numeric(p)
  } else {
    delta <- .as_location(delta, p, "delta")
  }
  alpha <- .zzg22nr_validate_alpha(alpha)
  variable.names <- .hotelling_variable_names(x)
  prepared <- .zzg22nr_prepare_paired(x, y, delta)

  .zzg22nr_fit(
    prepared = prepared,
    alpha = alpha,
    estimate = prepared$estimate,
    null.value = delta,
    variable.names = variable.names,
    method = paste(
      "Zhang-Zhou-Guo paired-sample normal-reference mean test",
      "(rowwise paired differences; three-cumulant chi-square)"
    ),
    data.name = paste(x.name, "and", y.name),
    n = c(pairs = nrow(x)),
    call = call
  )
}


#' @rdname zhang_zhou_guo_tests
#' @export
zhang_zhou_guo_linear_hypothesis_test <- function(
    x, L, rhs = NULL, alpha = 0.05) {
  call <- match.call()
  x.name <- deparse1(substitute(x))
  L.name <- deparse1(substitute(L))
  x <- .as_data_matrix(x, "x", min_rows = 4L)
  L <- .zzg22nr_as_contrast_matrix(L, ncol(x))
  if (!is.null(colnames(x)) && !is.null(colnames(L)) &&
      !identical(colnames(x), colnames(L))) {
    stop("When both are named, columns of `L` must match columns of `x` in order.",
         call. = FALSE)
  }
  q <- nrow(L)
  if (is.null(rhs)) {
    rhs <- numeric(q)
  } else {
    rhs <- .as_location(rhs, q, "rhs")
  }
  alpha <- .zzg22nr_validate_alpha(alpha)
  variable.names <- if (is.null(rownames(L))) {
    paste0("H", seq_len(q))
  } else {
    rownames(L)
  }
  prepared <- .zzg22nr_prepare_linear(x, L, rhs)

  result <- .zzg22nr_fit(
    prepared = prepared,
    alpha = alpha,
    estimate = prepared$estimate,
    null.value = rhs,
    variable.names = variable.names,
    method = paste(
      "Zhang-Zhou-Guo same-unit linear-hypothesis normal-reference test",
      "(not independent-group MANOVA)"
    ),
    data.name = paste(x.name, "with", L.name),
    n = c(x = nrow(x)),
    call = call
  )
  result$components$L <- L
  result$components$rhs <- stats::setNames(rhs, variable.names)
  result$diagnostics$independent.group.MANOVA <- FALSE
  result$diagnostics$linear.hypothesis.scope <-
    "one rowwise transform per iid observational unit"
  result
}
