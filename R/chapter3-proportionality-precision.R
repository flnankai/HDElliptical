.ch3pp_validate_strict <- function(strict) {
  if (!is.logical(strict) || length(strict) != 1L || is.na(strict)) {
    stop("`strict` must be TRUE or FALSE.", call. = FALSE)
  }
  strict
}


.ch3pp_validate_alpha <- function(alpha) {
  alpha <- as.numeric(alpha)
  if (length(alpha) != 1L || is.na(alpha) || !is.finite(alpha) ||
      alpha <= 0 || alpha >= 1) {
    stop("`alpha` must be one finite number strictly between zero and one.",
         call. = FALSE)
  }
  alpha
}


.ch3pp_warn_or_stop <- function(message, strict) {
  if (isTRUE(strict)) {
    stop(message, call. = FALSE)
  }
  warning(message, call. = FALSE)
  invisible(FALSE)
}


.ch3pp_spatial_fit <- function(x, tol, max_iter, zero_tol) {
  center <- spatial_median(
    x, tol = tol, max_iter = max_iter, zero_tol = zero_tol, warn = FALSE
  )
  scaled <- .center_and_scale(x, as.numeric(center), zero_tol)
  sign.fit <- cpp_spatial_sign(scaled$x, scaled$zero_tol)
  list(
    center = center,
    signs = as.matrix(sign.fit$signs),
    radii.standardized = as.numeric(sign.fit$norms),
    radial.scale = scaled$scale,
    overflow.fallback = scaled$overflow_fallback,
    n.zero = as.numeric(sign.fit$n_zero),
    diagnostics = list(
      converged = isTRUE(attr(center, "converged")),
      iterations = as.integer(attr(center, "iterations")),
      objective = as.numeric(attr(center, "objective")),
      relative.change = as.numeric(attr(center, "relative_change")),
      equation.residual = as.numeric(attr(center, "equation_residual"))
    )
  )
}


.ch3pp_inverse_radial_moments <- function(radii) {
  radii <- as.numeric(radii)
  if (length(radii) == 0L || anyNA(radii) || any(!is.finite(radii)) ||
      any(radii <= 0)) {
    return(NULL)
  }
  # All ratios in Cheng et al. are invariant to a common radial scale.  Using
  # min(R)/R keeps every inverse radius in [0, 1] without a floor or cap.
  normalization <- min(radii)
  inverse.radius <- normalization / radii
  moments <- vapply(1:3, function(k) mean(inverse.radius^k), numeric(1))
  names(moments) <- paste0("m", 1:3, ".scaled")
  list(
    moments = moments,
    inverse.radius.scaled = inverse.radius,
    radius.normalization.standardized = normalization
  )
}


.ch3pp_cheng_delta <- function(radial, n) {
  m1 <- unname(radial$moments[1L])
  m2 <- unname(radial$moments[2L])
  m3 <- unname(radial$moments[3L])
  second.order <- 2 - 2 * m2 / m1^2 + m2^2 / m1^4
  third.order <- -6 * m2^2 / m1^4 +
    2 * m2 * m3 / m1^5 + 8 * m2 / m1^2 - 2 * m3 / m1^3
  list(
    value = second.order / n^2 + third.order / n^3,
    second.order.coefficient = second.order,
    third.order.coefficient = third.order
  )
}


.ch3pp_test_object <- function(statistic, p.value, method, data.name,
                               raw.statistic, alpha, critical.value,
                               rejection, components, diagnostics) {
  structure(
    list(
      statistic = statistic,
      parameter = NULL,
      p.value = p.value,
      method = method,
      data.name = data.name,
      alternative = "two.sided",
      raw.statistic = raw.statistic,
      alpha = alpha,
      critical.value = critical.value,
      rejection = rejection,
      components = components,
      diagnostics = diagnostics
    ),
    class = c("hd_proportionality_test", "htest")
  )
}


#' Cheng--Liu--Peng--Zhang--Zheng two-sample SSCM test
#'
#' Tests equality of two population spatial-sign covariance matrices (SSCMs),
#' equivalently proportionality of the two scatter matrices under elliptical
#' symmetry.  With separately estimated spatial medians and centered signs,
#' the feasible statistic is
#' \deqn{T=p(A+B-2C),}
#' where \eqn{A} and \eqn{B} average squared inner products over ordered
#' within-group pairs and \eqn{C} averages them over all cross pairs.
#'
#' The implementation uses the complete finite-sample bias expression in
#' Cheng et al. (2019), Equation (3.3).  For each group it estimates the
#' inverse radial moments through order three, evaluates the \eqn{n^{-2}} and
#' \eqn{n^{-3}} terms separately, and subtracts \eqn{p\widehat\delta}.
#' Under the null, Remark 3.1 estimates
#' \deqn{q=p^{-1}\mathrm{tr}(\Lambda^2)
#' =\frac{p}{n_1+n_2}\{n_1(A-\delta_1)+n_2(B-\delta_2)\}.}
#' Consequently the *standard-error denominator* and its square are
#' \deqn{\widehat v_0=
#' \frac{2(n_1^{-1}+n_2^{-1})}{p+2}(p\widehat q),\qquad
#' \widehat v_0^2=\{\widehat v_0\}^2.}
#' The reported statistic is \eqn{Z=(T-p\widehat\delta)/\widehat v_0} and
#' the paper's rejection rule is the upper standard-normal tail.
#'
#' A non-positive trace estimate, a zero fitted radius, or an unconverged
#' spatial median invalidates calibration.  No absolute value, variance floor,
#' radius perturbation, or replacement center is used.  With `strict = FALSE`
#' the raw components are returned, but `statistic`, `p.value`, and the
#' rejection decision are `NA` and `diagnostics$calibrated` is `FALSE`.
#' The plug-in bias calibration has the additional regime stated in Remark
#' 3.2 of the primary paper: under bounded eigenvalues it requires
#' \eqn{p=o(n_l^{3/2})} for both groups.  This is stronger than merely using
#' the theorem's \eqn{p=O(n_l^2)} expansion condition.
#'
#' @param x,y Numeric matrices or data frames with observations in rows and
#'   the same variables in columns.  At least two rows are needed per group.
#' @param alpha Nominal level for the upper-tail calibration.
#' @param tol,max_iter Spatial-median convergence controls.
#' @param zero_tol Non-negative tolerance for a zero centered radius.  The
#'   default detects exact zeros only.
#' @param strict If `TRUE`, fail when the feasible calibration is undefined;
#'   otherwise return explicitly uncalibrated raw diagnostics.
#'
#' @return An object of class `c("hd_proportionality_test", "htest")`.
#'   `components` contains all ordered-pair components, bias terms, scaled
#'   inverse radial moments, the trace estimate, the standard-error denominator
#'   and its squared variance, fitted centers, signs, and SSCMs.
#'
#' @references
#' Cheng, G., Liu, B., Peng, L., Zhang, B., and Zheng, S. (2019).
#' Testing the equality of two high-dimensional spatial sign covariance
#' matrices. *Scandinavian Journal of Statistics*, **46**, 257--271.
#' \doi{10.1111/sjos.12350}.
#'
#' @examples
#' x <- matrix(c(-2, 0, 1, 2, 1, -1, 0, 2, -1, -2, 2, 0), 6, 2)
#' y <- matrix(c(-1, 1, 2, -2, 0, 3, 1, -1, 2, 0, -2, 1), 6, 2)
#' cheng_sscm_equality_test(x, y)
#'
#' @export
cheng_sscm_equality_test <- function(
    x, y, alpha = 0.05, tol = 1e-8, max_iter = 500L,
    zero_tol = 0, strict = TRUE) {
  call <- match.call()
  x.name <- deparse1(substitute(x))
  y.name <- deparse1(substitute(y))
  x <- .as_data_matrix(x, "x", min_rows = 2L)
  y <- .as_data_matrix(y, "y", min_rows = 2L)
  .check_two_sample_variables(x, y)
  alpha <- .ch3pp_validate_alpha(alpha)
  controls <- .validate_iteration_controls(tol, max_iter, zero_tol)
  strict <- .ch3pp_validate_strict(strict)

  fit.x <- .ch3pp_spatial_fit(
    x, controls$tol, controls$max_iter, controls$zero_tol
  )
  fit.y <- .ch3pp_spatial_fit(
    y, controls$tol, controls$max_iter, controls$zero_tol
  )
  kernel <- cpp_ch3pp_sscm_components(fit.x$signs, fit.y$signs)
  n1 <- nrow(x)
  n2 <- nrow(y)
  p <- ncol(x)
  T.raw <- p * (kernel$A + kernel$B - 2 * kernel$C)

  failures <- character()
  if (!isTRUE(fit.x$diagnostics$converged)) {
    failures <- c(failures, "the spatial median for `x` did not converge")
  }
  if (!isTRUE(fit.y$diagnostics$converged)) {
    failures <- c(failures, "the spatial median for `y` did not converge")
  }
  if (fit.x$n.zero > 0) {
    failures <- c(failures, "`x` has a zero fitted radius")
  }
  if (fit.y$n.zero > 0) {
    failures <- c(failures, "`y` has a zero fitted radius")
  }

  radial.x <- .ch3pp_inverse_radial_moments(fit.x$radii.standardized)
  radial.y <- .ch3pp_inverse_radial_moments(fit.y$radii.standardized)
  delta.x <- delta.y <- list(
    value = NA_real_, second.order.coefficient = NA_real_,
    third.order.coefficient = NA_real_
  )
  delta.hat <- q.hat <- standard.error <- variance <- centered <- NA_real_
  if (length(failures) == 0L && !is.null(radial.x) && !is.null(radial.y)) {
    delta.x <- .ch3pp_cheng_delta(radial.x, n1)
    delta.y <- .ch3pp_cheng_delta(radial.y, n2)
    delta.hat <- delta.x$value + delta.y$value
    q.hat <- p / (n1 + n2) * (
      n1 * (kernel$A - delta.x$value) +
        n2 * (kernel$B - delta.y$value)
    )
    standard.error <- 2 * (1 / n1 + 1 / n2) / (p + 2) * (p * q.hat)
    variance <- standard.error^2
    centered <- T.raw - p * delta.hat
    if (!is.finite(q.hat) || q.hat <= 0) {
      failures <- c(
        failures,
        "the feasible estimate of p^{-1} tr(Lambda^2) is not positive"
      )
    }
    if (!is.finite(standard.error) || standard.error <= 0 ||
        !is.finite(variance)) {
      failures <- c(failures, "the feasible null standard error is invalid")
    }
  } else if (length(failures) == 0L) {
    failures <- c(failures, "the inverse radial moments are undefined")
  }

  calibrated <- length(failures) == 0L
  if (!calibrated) {
    .ch3pp_warn_or_stop(
      paste0(
        "Cheng SSCM calibration failed: ", paste(unique(failures), collapse = "; "),
        ". No numerical repair or alternative calibration was applied."
      ),
      strict
    )
  }
  z <- if (calibrated) centered / standard.error else NA_real_
  p.value <- if (calibrated) stats::pnorm(z, lower.tail = FALSE) else NA_real_
  critical.value <- stats::qnorm(1 - alpha)
  rejection <- if (calibrated) isTRUE(z > critical.value) else NA

  variable.names <- colnames(x)
  sscm.x <- crossprod(fit.x$signs) / n1
  sscm.y <- crossprod(fit.y$signs) / n2
  if (!is.null(variable.names)) {
    dimnames(sscm.x) <- dimnames(sscm.y) <-
      list(variable.names, variable.names)
    names(fit.x$center) <- names(fit.y$center) <- variable.names
  }
  result <- .ch3pp_test_object(
    statistic = c(Z = z),
    p.value = p.value,
    method = paste(
      "Cheng-Liu-Peng-Zhang-Zheng two-sample SSCM equality test",
      "(feasible median-bias correction)"
    ),
    data.name = paste(x.name, "and", y.name),
    raw.statistic = c(T.SSCM = T.raw),
    alpha = alpha,
    critical.value = c(normal.upper = critical.value),
    rejection = rejection,
    components = list(
      A = as.numeric(kernel$A),
      B = as.numeric(kernel$B),
      C = as.numeric(kernel$C),
      ordered.sum.x = as.numeric(kernel$ordered_sum_x),
      ordered.sum.y = as.numeric(kernel$ordered_sum_y),
      cross.sum = as.numeric(kernel$cross_sum),
      T.SSCM = T.raw,
      delta1.hat = delta.x$value,
      delta2.hat = delta.y$value,
      delta.hat = delta.hat,
      bias.hat = p * delta.hat,
      centered.statistic = centered,
      trace.Lambda2.over.p.hat = q.hat,
      trace.Lambda2.hat = p * q.hat,
      null.standard.error.hat = standard.error,
      null.variance.hat = variance,
      delta1.coefficients = unlist(delta.x[-1L]),
      delta2.coefficients = unlist(delta.y[-1L]),
      radial1 = radial.x,
      radial2 = radial.y,
      center.x = fit.x$center,
      center.y = fit.y$center,
      signs.x = fit.x$signs,
      signs.y = fit.y$signs,
      SSCM.x = sscm.x,
      SSCM.y = sscm.y
    ),
    diagnostics = list(
      calibrated = calibrated,
      failures = unique(failures),
      spatial.median.x = fit.x$diagnostics,
      spatial.median.y = fit.y$diagnostics,
      zero.radii = c(x = fit.x$n.zero, y = fit.y$n.zero),
      radial.scaling = paste(
        "Inverse radii use min(R)/R; all paper moment ratios are unchanged"
      ),
      denominator = "null.standard.error.hat",
      variance = "square of null.standard.error.hat",
      pair.convention = paste(
        "A and B use ordered unequal pairs; C uses every cross pair"
      ),
      tail = "upper normal tail for a two-sided scientific alternative",
      rejection.comparison = "> (complement of the paper's <= acceptance region)",
      plug.in.bias.condition = paste(
        "Primary Remark 3.2: p=o(n_l^(3/2)) under bounded eigenvalues",
        "when delta is estimated"
      ),
      book.erratum = paste(
        "The draft states only p=O(n_l^2) near the feasible test and cites",
        "Theorem 2; the plug-in condition is in Remark 3.2 after Theorem 3.2"
      ),
      no.repair = paste(
        "No ridge, radius perturbation, absolute-value trace repair,",
        "variance floor, replacement center, or alternative calibration"
      ),
      call = call
    )
  )
  result
}


#' Feng--Zhang--Liu spatial-rank proportionality test
#'
#' Implements the high-dimensional spatial-rank test of Feng, Zhang, and Liu
#' (2022).  For each group, the within component averages
#' \eqn{\{U(X_i-X_j)^T U(X_k-X_l)\}^2} over four *ordered, mutually
#' distinct* indices.  The cross component averages over an ordered unequal
#' pair in each group.  If the corresponding unscaled averages are
#' \eqn{A_1,A_2,C_{12}}, then
#' \deqn{T_{HT}=p(A_1+A_2-2C_{12}).}
#'
#' The feasible null variance in the primary paper is implemented literally:
#' \deqn{\widehat\sigma_{0,n}^2=
#' \frac{4(n_1^{-1}+n_2^{-1})^2}{(n_1+n_2)(p+2)^2}
#' p^2(n_1A_1+n_2A_2).}
#' The statistic is calibrated by the upper normal tail.  Pairwise ties use the
#' paper's convention \eqn{U(0)=0}.  A non-positive feasible variance is not
#' floored or replaced.
#'
#' The primary paper normalizes its shape matrix \eqn{\Lambda} to trace one,
#' and the estimand is \eqn{p\,tr\{(S_1-S_2)^2\}} for the SSCM/Kendall
#' functionals.  This is not the trace-\eqn{p} shape-scale formula printed in
#' the current book draft; multiplying \eqn{\Lambda} by \eqn{p} requires the
#' corresponding powers of \eqn{p} in every variance trace.
#'
#' @param x,y Numeric matrices or data frames with at least four observations
#'   per group and the same variables in columns.
#' @param alpha Nominal level for the paper's upper-tail rejection rule.
#' @param zero_tol Non-negative tolerance below which a pairwise difference is
#'   assigned the exact zero spatial sign.
#' @param strict If `TRUE`, fail when the feasible variance is invalid;
#'   otherwise return raw components with no calibrated statistic or p-value.
#'
#' @return An `hd_proportionality_test` object.  Both raw U-statistic averages
#'   and their \eqn{p}-scaled contributions are returned.
#'
#' @references
#' Feng, L., Zhang, X., and Liu, B. (2022). High-dimensional proportionality
#' test of two covariance matrices and its application to gene expression
#' data. *Statistical Theory and Related Fields*, **6**, 161--174.
#' \doi{10.1080/24754269.2021.1984373}.
#'
#' @examples
#' x <- matrix(seq_len(40), 8, 5) + matrix(c(-1, 0, 1, 0), 8, 5)
#' y <- matrix(seq_len(45), 9, 5) + matrix(c(0, 1, -1), 9, 5)
#' feng_spatial_rank_proportionality_test(x, y)
#'
#' @export
feng_spatial_rank_proportionality_test <- function(
    x, y, alpha = 0.05, zero_tol = 0, strict = TRUE) {
  call <- match.call()
  x.name <- deparse1(substitute(x))
  y.name <- deparse1(substitute(y))
  x <- .as_data_matrix(x, "x", min_rows = 4L)
  y <- .as_data_matrix(y, "y", min_rows = 4L)
  .check_two_sample_variables(x, y)
  alpha <- .ch3pp_validate_alpha(alpha)
  zero_tol <- .validate_zero_tol(zero_tol)
  strict <- .ch3pp_validate_strict(strict)

  kernel <- cpp_ch3pp_spatial_rank_components(x, y, zero_tol)
  n1 <- nrow(x)
  n2 <- nrow(y)
  p <- ncol(x)
  A1.raw <- as.numeric(kernel$A1_raw)
  A2.raw <- as.numeric(kernel$A2_raw)
  C12.raw <- as.numeric(kernel$C12_raw)
  T.HT <- p * (A1.raw + A2.raw - 2 * C12.raw)
  variance <- 4 * (1 / n1 + 1 / n2)^2 /
    ((n1 + n2) * (p + 2)^2) * p^2 *
    (n1 * A1.raw + n2 * A2.raw)
  standard.error <- sqrt(variance)
  calibrated <- is.finite(variance) && variance > 0 &&
    is.finite(standard.error) && standard.error > 0
  failures <- if (calibrated) character() else
    "the primary-paper feasible null variance is not finite and positive"
  if (!calibrated) {
    .ch3pp_warn_or_stop(
      paste0(
        "Feng-Zhang-Liu spatial-rank calibration failed: ", failures,
        ". No variance floor or replacement calibration was applied."
      ),
      strict
    )
  }
  z <- if (calibrated) T.HT / standard.error else NA_real_
  p.value <- if (calibrated) stats::pnorm(z, lower.tail = FALSE) else NA_real_
  critical.value <- stats::qnorm(1 - alpha)
  rejection <- if (calibrated) isTRUE(z >= critical.value) else NA

  .ch3pp_test_object(
    statistic = c(Z = z),
    p.value = p.value,
    method = "Feng-Zhang-Liu high-dimensional spatial-rank proportionality test",
    data.name = paste(x.name, "and", y.name),
    raw.statistic = c(T.HT = T.HT),
    alpha = alpha,
    critical.value = c(normal.upper = critical.value),
    rejection = rejection,
    components = list(
      A1.raw = A1.raw,
      A2.raw = A2.raw,
      C12.raw = C12.raw,
      A1 = p * A1.raw,
      A2 = p * A2.raw,
      C12 = p * C12.raw,
      T.HT = T.HT,
      ordered.sum.x = as.numeric(kernel$ordered_sum_x),
      ordered.sum.y = as.numeric(kernel$ordered_sum_y),
      ordered.cross.sum = as.numeric(kernel$ordered_cross_sum),
      sigma0.squared.hat = variance,
      sigma0.hat = standard.error,
      pair.signs.x = as.matrix(kernel$pair_signs_x),
      pair.signs.y = as.matrix(kernel$pair_signs_y),
      zero.pairs = c(
        x = as.numeric(kernel$zero_pairs_x),
        y = as.numeric(kernel$zero_pairs_y)
      )
    ),
    diagnostics = list(
      calibrated = calibrated,
      failures = failures,
      pair.convention = paste(
        "within sums use ordered four-distinct indices (P_n^4);",
        "cross sums use ordered unequal pairs in both groups"
      ),
      variance.formula = paste(
        "4*(1/n1+1/n2)^2/(n1+n2)/(p+2)^2 *",
        "p^2*(n1*A1.raw+n2*A2.raw)"
      ),
      primary.normalization = paste(
        "Lambda has trace one; target is p*tr((S1-S2)^2) for SSCMs"
      ),
      book.erratum = paste(
        "The draft's trace-p Theta substitution omits required powers of p"
      ),
      tail = "upper normal tail for a two-sided scientific alternative",
      rejection.comparison = ">= (as printed in the primary paper)",
      zero.sign = "U(0)=0 for tied pairwise differences",
      no.repair = "No ridge, perturbation, variance floor, or substitution",
      call = call
    )
  )
}


.ch3pp_validate_precision_controls <- function(
    lambda, solver_tol, solver_max_iter, initial_step, max_backtracking) {
  lambda <- as.numeric(lambda)
  if (length(lambda) != 1L || is.na(lambda) || !is.finite(lambda) ||
      lambda <= 0) {
    stop("`lambda` must be one finite positive number.", call. = FALSE)
  }
  solver_tol <- as.numeric(solver_tol)
  if (length(solver_tol) != 1L || is.na(solver_tol) ||
      !is.finite(solver_tol) || solver_tol <= 0) {
    stop("`solver_tol` must be one finite positive number.", call. = FALSE)
  }
  validate.integer <- function(value, name) {
    if (!is.numeric(value) || length(value) != 1L || is.na(value) ||
        !is.finite(value) || value < 1 || value != floor(value) ||
        value > .Machine$integer.max) {
      stop(sprintf("`%s` must be a positive integer.", name), call. = FALSE)
    }
    as.integer(value)
  }
  initial_step <- as.numeric(initial_step)
  if (length(initial_step) != 1L || is.na(initial_step) ||
      !is.finite(initial_step) || initial_step <= 0) {
    stop("`initial_step` must be one finite positive number.", call. = FALSE)
  }
  list(
    lambda = lambda,
    solver_tol = solver_tol,
    solver_max_iter = validate.integer(solver_max_iter, "solver_max_iter"),
    initial_step = initial_step,
    max_backtracking = validate.integer(max_backtracking, "max_backtracking")
  )
}


#' Robust spatial-sign SCLIME and SGLASSO precision estimation
#'
#' Estimates the inverse trace-normalized shape
#' \eqn{V_0=\Lambda_0^{-1}=\{tr(\Sigma_0)/p\}\Sigma_0^{-1}} using the
#' spatial-sign covariance matrix \eqn{\widehat S} centered at the ordinary
#' sample spatial median.  Set \eqn{A=p\widehat S}.  With `method = "sclime"`
#' the paper's problem
#' \deqn{\min_V\|V\|_1\quad\text{subject to}\quad
#' \|AV-I\|_\infty\le\lambda}
#' is solved column by column.  Every column must pass primal feasibility,
#' dual feasibility, l1 stationarity, and relative primal-dual gap checks.
#' The raw solution is then symmetrized exactly as in Lu and Feng (2025): for
#' each transposed pair, retain the entry with smaller absolute value.
#'
#' With `method = "sglasso"`, the function solves
#' \deqn{\min_{V\succ0}\{tr(AV)-\log\det V+\lambda\|V\|_1\}.}
#' Here \eqn{\|V\|_1} is the paper's full elementwise norm, so the diagonal is
#' penalized.  A positive-definite backtracking proximal-gradient algorithm is
#' used and the full subgradient KKT residual must not exceed `solver_tol`.
#'
#' The SCLIME primal-dual iteration is a convergent Chambolle--Pock method;
#' SGLASSO uses the standard convex proximal-gradient majorization inequality.
#' An optimizer stopping because of `solver_max_iter` is not called a valid
#' estimator.  With `strict = FALSE`, a failed run returns `estimate = NULL`
#' plus its last-iterate certificate; it never returns a matrix that looks like
#' a valid paper estimator.  Neither method uses a ridge, eigenvalue floor,
#' pseudoinverse, constraint relaxation, or post-hoc KKT repair.
#'
#' The constants in the paper's theoretical choices of \eqn{\lambda_n} are
#' not observable tuning formulas.  Therefore `lambda` must be supplied rather
#' than silently replacing them by an undocumented default.
#'
#' @param x Numeric matrix or data frame with observations in rows.
#' @param lambda Positive tuning parameter in the primary-paper objective.
#' @param method Either `"sclime"` or `"sglasso"`.
#' @param median_tol,median_max_iter,zero_tol Spatial-median and zero-sign
#'   controls.
#' @param solver_tol Positive tolerance required for every solver certificate.
#' @param solver_max_iter Positive maximum number of primal-dual or proximal
#'   iterations.
#' @param initial_step Initial SGLASSO proximal step; ignored by SCLIME.
#' @param max_backtracking Maximum SGLASSO line-search reductions per update.
#' @param strict If `TRUE`, fail on a spatial-median or solver certificate
#'   failure.  If `FALSE`, return an explicitly invalid diagnostic fit with
#'   `estimate = NULL`.
#'
#' @return An object of class `spatial_sign_precision_fit`.  A successful fit
#'   contains `estimate`, the sample SSCM, \eqn{p\widehat S}, the spatial
#'   median, and complete feasibility/KKT diagnostics.
#'
#' @references
#' Lu, Z. and Feng, L. (2025). Robust sparse precision matrix estimation and
#' its applications. arXiv:2503.03575.
#' \url{https://arxiv.org/abs/2503.03575}.
#'
#' @examples
#' x <- rbind(
#'   c(-2, 0, 1), c(-1, 1, 0), c(0, -1, 2), c(1, 0, -1),
#'   c(2, 1, 1), c(0, 2, -2), c(-1, -2, 0), c(1, -1, 1)
#' )
#' spatial_sign_precision(x, lambda = 0.4, method = "sglasso")
#'
#' @export
spatial_sign_precision <- function(
    x, lambda, method = c("sclime", "sglasso"),
    median_tol = 1e-8, median_max_iter = 500L, zero_tol = 0,
    solver_tol = 1e-7, solver_max_iter = 100000L,
    initial_step = 1, max_backtracking = 100L, strict = TRUE) {
  call <- match.call()
  x.name <- deparse1(substitute(x))
  x <- .as_data_matrix(x, "x", min_rows = 2L)
  method <- match.arg(method)
  median.controls <- .validate_iteration_controls(
    median_tol, median_max_iter, zero_tol
  )
  solver.controls <- .ch3pp_validate_precision_controls(
    lambda, solver_tol, solver_max_iter, initial_step, max_backtracking
  )
  strict <- .ch3pp_validate_strict(strict)
  fit <- .ch3pp_spatial_fit(
    x, median.controls$tol, median.controls$max_iter,
    median.controls$zero_tol
  )
  variable.names <- colnames(x)
  p <- ncol(x)
  sscm.hat <- crossprod(fit$signs) / nrow(x)
  coefficient <- p * sscm.hat
  if (!is.null(variable.names)) {
    dimnames(sscm.hat) <- dimnames(coefficient) <-
      list(variable.names, variable.names)
    names(fit$center) <- variable.names
  }

  failed.fit <- function(stage, message, solver = NULL) {
    .ch3pp_warn_or_stop(message, strict)
    structure(
      list(
        estimate = NULL, valid = FALSE, method = method,
        lambda = solver.controls$lambda, center = fit$center,
        sscm = sscm.hat, coefficient = coefficient,
        data.name = x.name,
        diagnostics = list(
          failure.stage = stage,
          failure = message,
          spatial.median = fit$diagnostics,
          zero.radii = fit$n.zero,
          solver = solver,
          no.repair = paste(
            "No ridge, eigenvalue floor, pseudoinverse, constraint",
            "relaxation, or post-hoc KKT repair"
          ),
          call = call
        )
      ),
      class = "spatial_sign_precision_fit"
    )
  }

  if (!isTRUE(fit$diagnostics$converged)) {
    message <- paste0(
      "Spatial-sign precision estimation failed before optimization: the ",
      "sample spatial median did not converge. No precision estimate was ",
      "constructed."
    )
    return(failed.fit("spatial median", message))
  }
  if (!is.finite(sum(abs(coefficient))) || sum(abs(coefficient)) == 0) {
    message <- paste0(
      "Spatial-sign precision estimation failed before optimization: all ",
      "fitted spatial signs are zero, so the sample SSCM contains no ",
      "directional information. No precision estimate was constructed."
    )
    return(failed.fit("degenerate SSCM", message))
  }

  if (method == "sclime") {
    solver <- tryCatch(
      cpp_ch3pp_sclime(
        coefficient, solver.controls$lambda, solver.controls$solver_tol,
        solver.controls$solver_max_iter
      ),
      error = identity
    )
    if (inherits(solver, "condition")) {
      message <- paste0(
        "Spatial-sign SCLIME optimization failed: ", conditionMessage(solver),
        " No precision estimate was returned and no constraint repair was ",
        "applied."
      )
      return(failed.fit(
        "optimization error", message,
        solver = list(error = conditionMessage(solver))
      ))
    }
    valid <- isTRUE(solver$all_converged) &&
      all(solver$primal_violation <= solver.controls$solver_tol) &&
      all(solver$dual_violation <= solver.controls$solver_tol) &&
      all(solver$stationarity_residual <= solver.controls$solver_tol) &&
      all(solver$relative_gap <= solver.controls$solver_tol) &&
      is.finite(solver$final_feasibility_violation) &&
      solver$final_feasibility_violation <= solver.controls$solver_tol
    certificate <- list(
      all.columns.certified = valid,
      column.converged = as.logical(solver$converged),
      iterations = as.integer(solver$iterations),
      relative.update = as.numeric(solver$relative_update),
      primal.violation = as.numeric(solver$primal_violation),
      stationarity.residual = as.numeric(solver$stationarity_residual),
      dual.violation = as.numeric(solver$dual_violation),
      primal.objective = as.numeric(solver$primal_objective),
      dual.objective = as.numeric(solver$dual_objective),
      duality.gap = as.numeric(solver$duality_gap),
      relative.gap = as.numeric(solver$relative_gap),
      raw.maximum.feasibility.violation = max(solver$primal_violation),
      symmetrized.feasibility.violation = as.numeric(
        solver$final_feasibility_violation
      ),
      symmetrized.minimum.eigenvalue = as.numeric(
        solver$final_minimum_eigenvalue
      ),
      operator.norm = as.numeric(solver$operator_norm),
      primal.dual.step = as.numeric(solver$primal_dual_step)
    )
    last.iterate <- list(
      raw = as.matrix(solver$raw_solution),
      symmetrized = as.matrix(solver$solution),
      dual = as.matrix(solver$dual_solution)
    )
  } else {
    solver <- tryCatch(
      cpp_ch3pp_sglasso(
        coefficient, solver.controls$lambda, solver.controls$solver_tol,
        solver.controls$solver_max_iter, solver.controls$initial_step,
        solver.controls$max_backtracking
      ),
      error = identity
    )
    if (inherits(solver, "condition")) {
      message <- paste0(
        "Spatial-sign SGLASSO optimization failed: ", conditionMessage(solver),
        " No precision estimate was returned and no positive-definite repair ",
        "was applied."
      )
      return(failed.fit(
        "optimization error", message,
        solver = list(error = conditionMessage(solver))
      ))
    }
    valid <- isTRUE(solver$converged) &&
      !isTRUE(solver$backtracking_failed) &&
      is.finite(solver$minimum_eigenvalue) &&
      solver$minimum_eigenvalue > 0 &&
      is.finite(solver$kkt_residual) &&
      solver$kkt_residual <= solver.controls$solver_tol
    certificate <- list(
      kkt.certified = valid,
      converged = isTRUE(solver$converged),
      backtracking.failed = isTRUE(solver$backtracking_failed),
      iterations = as.integer(solver$iterations),
      relative.update = as.numeric(solver$relative_update),
      kkt.residual = as.numeric(solver$kkt_residual),
      objective = as.numeric(solver$objective),
      minimum.eigenvalue = as.numeric(solver$minimum_eigenvalue),
      accepted.step = as.numeric(solver$accepted_step),
      total.backtracking = as.integer(solver$total_backtracking)
    )
    last.iterate <- list(
      symmetrized = as.matrix(solver$solution),
      inverse = as.matrix(solver$inverse)
    )
  }

  if (!is.null(variable.names)) {
    last.iterate <- lapply(last.iterate, function(value) {
      dimnames(value) <- list(variable.names, variable.names)
      value
    })
  }
  if (!valid) {
    message <- paste0(
      "Spatial-sign ", toupper(method), " optimization failed its ",
      "feasibility/KKT certificate. No precision estimate was returned and ",
      "no ridge, tolerance relaxation, or post-hoc repair was applied."
    )
    .ch3pp_warn_or_stop(message, strict)
  }
  estimate <- if (valid) last.iterate$symmetrized else NULL

  structure(
    list(
      estimate = estimate,
      valid = valid,
      method = method,
      lambda = solver.controls$lambda,
      center = fit$center,
      signs = fit$signs,
      sscm = sscm.hat,
      coefficient = coefficient,
      data.name = x.name,
      diagnostics = list(
        failure.stage = if (valid) NULL else "optimization certificate",
        spatial.median = fit$diagnostics,
        zero.radii = fit$n.zero,
        radial.scale = fit$radial.scale,
        overflow.fallback = fit$overflow.fallback,
        solver = certificate,
        last.iterate = if (valid) NULL else last.iterate,
        target = paste(
          "V0 = Lambda0^{-1} = tr(Sigma0)/p * Sigma0^{-1}"
        ),
        scale.identification = "Lambda0 has trace p",
        diagonal.penalty = method == "sglasso",
        sclime.symmetrization = if (method == "sclime")
          "retain the transposed entry with smaller absolute value" else NULL,
        no.repair = paste(
          "No ridge, eigenvalue floor, pseudoinverse, constraint relaxation,",
          "or post-hoc KKT repair"
        ),
        call = call
      )
    ),
    class = "spatial_sign_precision_fit"
  )
}


#' Threshold a spatial-sign precision estimate
#'
#' Applies the support-recovery rule of Lu and Feng (2025),
#' \deqn{\widetilde v_{ij}=\widehat v_{ij}
#' 1\{|\widehat v_{ij}|\ge\tau\}.}
#' The comparison is non-strict and applies to all entries, including the
#' diagonal, exactly as defined in the primary paper.  The paper's theoretical
#' threshold contains unknown constants, so `tau` is deliberately explicit.
#'
#' @param fit A valid object returned by `spatial_sign_precision()`.
#' @param tau One finite non-negative threshold.
#'
#' @return A list containing the thresholded matrix, logical support, sign
#'   matrix, threshold, and source fit.
#'
#' @references
#' Lu, Z. and Feng, L. (2025). Robust sparse precision matrix estimation and
#' its applications. arXiv:2503.03575.
#'
#' @examples
#' x <- rbind(
#'   c(-2, 0), c(-1, 1), c(0, -1), c(1, 0), c(2, 1), c(0, 2)
#' )
#' fit <- spatial_sign_precision(x, 0.4, method = "sglasso")
#' threshold_spatial_sign_precision(fit, tau = 0.1)
#'
#' @export
threshold_spatial_sign_precision <- function(fit, tau) {
  if (!inherits(fit, "spatial_sign_precision_fit")) {
    stop("`fit` must be returned by `spatial_sign_precision()`.",
         call. = FALSE)
  }
  if (!isTRUE(fit$valid) || is.null(fit$estimate)) {
    stop("Cannot threshold an invalid or uncertified precision fit.",
         call. = FALSE)
  }
  tau <- as.numeric(tau)
  if (length(tau) != 1L || is.na(tau) || !is.finite(tau) || tau < 0) {
    stop("`tau` must be one finite non-negative number.", call. = FALSE)
  }
  keep <- abs(fit$estimate) >= tau
  estimate <- fit$estimate
  estimate[!keep] <- 0
  structure(
    list(
      estimate = estimate,
      support = keep & estimate != 0,
      sign = sign(estimate),
      threshold = tau,
      comparison = ">=",
      diagonal.thresholded = TRUE,
      source = fit
    ),
    class = "thresholded_spatial_sign_precision"
  )
}
