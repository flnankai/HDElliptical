.ch4c_alpha_result <- function(statistic, p.value, method, data.name,
                               raw.statistic, components, diagnostics,
                               call, estimate = NULL,
                               parameter = NULL) {
  answer <- list(
    statistic = statistic,
    parameter = parameter,
    p.value = as.numeric(p.value),
    method = method,
    data.name = data.name,
    alternative = "greater",
    null.value = stats::setNames(0, "alpha"),
    raw.statistic = raw.statistic,
    components = components,
    diagnostics = diagnostics,
    call = call
  )
  if (!is.null(estimate)) answer$estimate <- estimate
  answer <- answer[!vapply(answer, is.null, logical(1))]
  class(answer) <- c("hd_alpha_test", "htest")
  answer
}


.ch4c_weighted_formula <- function(directions, radii, h, trace_R2,
                                   weights, weight.label) {
  if (!is.matrix(directions) || !is.numeric(directions) ||
      nrow(directions) < 2L || ncol(directions) < 2L ||
      anyNA(directions) || any(!is.finite(directions))) {
    stop("`directions` must be a finite numeric matrix with at least two rows and columns.",
         call. = FALSE)
  }
  storage.mode(directions) <- "double"
  observations <- nrow(directions)
  if (!is.numeric(radii) || !is.null(dim(radii)) ||
      length(radii) != observations || anyNA(radii) ||
      any(!is.finite(radii)) || any(radii <= 0)) {
    stop("`radii` must contain one strictly positive finite value per direction.",
         call. = FALSE)
  }
  if (!is.numeric(h) || !is.null(dim(h)) ||
      length(h) != observations || anyNA(h) || any(!is.finite(h)) ||
      !(sum(h^2) > 0)) {
    stop("`h` must contain one finite value per direction and have positive squared norm.",
         call. = FALSE)
  }
  if (!is.numeric(trace_R2) || length(trace_R2) != 1L ||
      is.na(trace_R2) || !is.finite(trace_R2) || trace_R2 <= 0) {
    stop("`trace_R2` must be one strictly positive finite number; no absolute-value or floor repair is used.",
         call. = FALSE)
  }
  if (!is.numeric(weights) || !is.null(dim(weights)) ||
      length(weights) != observations || anyNA(weights) ||
      any(!is.finite(weights))) {
    stop("The evaluated radial weights must be finite and have one value per direction.",
         call. = FALSE)
  }
  radii <- as.numeric(radii)
  h <- as.numeric(h)
  weights <- as.numeric(weights)
  trace_R2 <- as.numeric(trace_R2)
  row.square <- rowSums(directions^2)
  if (any(!is.finite(row.square)) ||
      any(abs(row.square - 1) > 1e-8)) {
    stop("Every supplied direction must have unit Euclidean norm; zero or non-unit rows are not repaired.",
         call. = FALSE)
  }
  core <- cpp_ch4_completion_weighted_alpha_q(directions, h, weights)
  log.denominator <- 0.5 * log(2) + log(core$psi2) +
    0.5 * log(trace_R2)
  if (!is.finite(log.denominator) ||
      log.denominator > log(.Machine$double.xmax)) {
    stop("The weighted-alpha studentizing denominator is not representable; no weight cap is applied.",
         call. = FALSE)
  }
  denominator <- exp(log.denominator)
  statistic <- core$Q / denominator
  if (!is.finite(statistic)) {
    stop("The weighted-alpha statistic is not finite; no repair is applied.",
         call. = FALSE)
  }
  list(
    statistic = statistic,
    p.value = stats::pnorm(statistic, lower.tail = FALSE),
    Q = as.numeric(core$Q),
    psi2 = as.numeric(core$psi2),
    trace_R2 = trace_R2,
    denominator = denominator,
    weights = weights,
    directions = directions,
    radii = radii,
    h = h,
    h2 = as.numeric(core$h2),
    weight.label = weight.label
  )
}


#' Oracle weighted spatial-sign alpha statistic
#'
#' Evaluates the formula-complete weighted spatial-sign class in Chapter 4
#' from supplied nuisance-free components. For radial weight function
#' \eqn{K}, the statistic is
#' \deqn{Q_K = \frac{N}{h^T h}\sum_{s \ne t}h_s h_t
#' K(r_s)K(r_t)U_s^T U_t,}
#' divided by
#' \eqn{\{2\hat\psi_{2,K}^2\widehat{\mathrm{tr}(R^2)}\}^{1/2}},
#' where \eqn{\hat\psi_{2,K}=T^{-1}\sum_t K^2(r_t)}.
#'
#' This interface deliberately says `oracle`: it does not estimate directions,
#' radii, factor slopes, diagonal scale, or the trace term, and it does not
#' claim primary-paper feasibility for an arbitrary user function. Use
#' `zhao_chen_zi_inst_alpha_test()` for the primary inverse-norm endpoint.
#'
#' @param directions Numeric \eqn{T} by \eqn{N} matrix of unit spatial signs.
#' @param radii Strictly positive finite length-\eqn{T} radial vector.
#' @param h Finite length-\eqn{T} residualized-intercept vector.
#' @param trace_R2 Strictly positive supplied value of
#'   \eqn{\widehat{\mathrm{tr}(R^2)}}.
#' @param K R function mapping the complete radial vector to a numeric vector
#'   of the same length.
#' @param data_name Optional label used in the returned `htest` object.
#' @return An upper-tail alpha-test object whose diagnostics explicitly record
#'   the oracle scope.
#' @references Zhao, P., Chen, D. and Zi, X. (2022). High-dimensional
#' non-parametric tests for linear asset pricing models. Stat 11, e490.
#' \doi{10.1002/sta4.490}
#' @examples
#' u <- rbind(c(1, 0), c(0, 1), c(-1, 0), c(0, -1))
#' weighted_spatial_sign_alpha_oracle_test(
#'   u, c(1, 2, 3, 4), rep(1, 4), trace_R2 = 3,
#'   K = function(r) 1 / r
#' )
#' @export
weighted_spatial_sign_alpha_oracle_test <- function(
    directions, radii, h, trace_R2, K, data_name = "supplied oracle scores") {
  call <- match.call()
  if (!is.function(K)) {
    stop("`K` must be an R function; no inferred or silently truncated radial weight is used.",
         call. = FALSE)
  }
  if (!is.numeric(radii) || !is.null(dim(radii))) {
    stop("`radii` must be a numeric vector; coercion is not used.",
         call. = FALSE)
  }
  radii.numeric <- as.numeric(radii)
  weights <- tryCatch(
    K(radii.numeric),
    error = function(e) {
      stop(sprintf("`K` failed on the supplied radii: %s", conditionMessage(e)),
           call. = FALSE)
    }
  )
  formula <- .ch4c_weighted_formula(
    directions, radii.numeric, h, trace_R2, weights,
    weight.label = "user-supplied oracle K"
  )
  .ch4c_alpha_result(
    statistic = stats::setNames(formula$statistic, "weighted sign Z"),
    p.value = formula$p.value,
    method = "Book weighted spatial-sign alpha statistic (oracle components)",
    data.name = as.character(data_name)[1L],
    raw.statistic = stats::setNames(formula$Q, "Q.K.alpha"),
    components = formula,
    diagnostics = list(
      scope = "oracle supplied-score formula only",
      primary.feasible.general.K = FALSE,
      nuisance.estimation = "none",
      weight.cap = "none",
      variance.repair = "none"
    ),
    call = call
  )
}


#' Zhao--Chen--Zi inverse-norm sign alpha test
#'
#' Implements the primary inverse-norm endpoint of the weighted spatial-sign
#' alpha class. The spatial directions and diagonal scale come from restricted
#' factor residuals \eqn{Y_t-\hat Bf_t}. Following the displayed Chapter 4
#' formula, the inverse radial weights use unrestricted OLS residuals
#' \eqn{Y_t-\hat\alpha-\hat Bf_t} and the same diagonal. The statistic sets
#' \eqn{K(r)=r^{-1}} in the weighted quadratic form and uses the audited
#' leave-two-out estimate of \eqn{\mathrm{tr}(R^2)}.
#'
#' A zero unrestricted radius makes inverse weighting undefined and is an
#' error. No radius floor, weight cap, absolute-value variance repair, ridge,
#' or generalized inverse is applied. Arbitrary supplied weights belong to the
#' explicitly named oracle interface
#' `weighted_spatial_sign_alpha_oracle_test()`.
#'
#' @inheritParams liu_feng_ma_spatial_sign_alpha_test
#' @param keep_scores Whether to retain the fitted directions, radii, weights,
#'   and residual matrices.
#' @return An upper-tail alpha-test object with literal weighted components and
#'   nuisance-fit diagnostics.
#' @references Zhao, P., Chen, D. and Zi, X. (2022). High-dimensional
#' non-parametric tests for linear asset pricing models. Stat 11, e490.
#' \doi{10.1002/sta4.490}
#' @examples
#' \donttest{
#' f <- cbind(seq(-1, 1, length.out = 12))
#' y <- outer(seq_len(12), 1:3, function(i, j) sin(i + j / 3))
#' zhao_chen_zi_inst_alpha_test(y, f)
#' }
#' @export
zhao_chen_zi_inst_alpha_test <- function(
    returns, factors = NULL, tol = 1e-8, max_iter = 1000L,
    zero_tol = 0, keep_scores = FALSE) {
  call <- match.call()
  data.name <- deparse1(substitute(returns))
  keep_scores <- .ch4_alpha_flag(keep_scores, "keep_scores")
  controls <- .ch4_alpha_controls(tol, max_iter, zero_tol)
  prepared <- .ch4_alpha_prepare(returns, factors, min_assets = 2L)
  core <- .ch4_alpha_lfm_core(prepared, controls, compute_trace = TRUE)
  trace_R2 <- as.numeric(core$trace_R2)
  if (length(trace_R2) != 1L || is.na(trace_R2) ||
      !is.finite(trace_R2) || trace_R2 <= 0) {
    stop("The leave-two-out trace estimate is not strictly positive and finite; no absolute-value or floor repair is used.",
         call. = FALSE)
  }
  ols <- .ch4_alpha_ols(prepared)
  radii <- as.numeric(cpp_ch4_completion_standardized_radii(
    ols$residuals.scaled, as.numeric(core$diagonal)
  ))
  if (any(radii <= controls$zero_tol)) {
    stop("An unrestricted standardized residual radius is zero or below `zero_tol`; inverse-norm weighting is undefined and no floor is applied.",
         call. = FALSE)
  }
  weights <- 1 / radii
  if (any(!is.finite(weights))) {
    stop("An inverse radial weight is not finite; no cap is applied.",
         call. = FALSE)
  }
  formula <- .ch4c_weighted_formula(
    core$directions, radii, as.numeric(core$h), trace_R2, weights,
    weight.label = "K(r) = 1/r"
  )
  canonical.diagonal <- .ch4_alpha_canonical_diagonal(
    core$diagonal, prepared$y.scale
  )
  names(canonical.diagonal$value) <- prepared$asset.names
  names(canonical.diagonal$log) <- prepared$asset.names
  components <- list(
    Q.K.alpha = formula$Q,
    psi2.K.hat = formula$psi2,
    trace.R.squared = formula$trace_R2,
    studentizing.denominator = formula$denominator,
    h = as.numeric(core$h),
    h.squared.norm = formula$h2,
    scale.diagonal = canonical.diagonal$value,
    log.scale.diagonal = canonical.diagonal$log,
    scale.diagonal.scaled.coordinates = as.numeric(core$diagonal),
    diagonal.iterations = as.integer(core$iterations),
    diagonal.equation.residual = as.numeric(core$equation_residual),
    diagonal.relative.update = as.numeric(core$relative_update),
    T = prepared$T,
    N = prepared$N,
    K = prepared$K
  )
  if (keep_scores) {
    components <- c(components, list(
      directions = formula$directions,
      unrestricted.radii = formula$radii,
      inverse.radial.weights = formula$weights,
      restricted.residuals.scaled = core$restricted_residuals,
      unrestricted.residuals.scaled = ols$residuals.scaled
    ))
  }
  .ch4c_alpha_result(
    statistic = stats::setNames(formula$statistic, "INST alpha Z"),
    p.value = formula$p.value,
    method = "Zhao-Chen-Zi inverse-norm spatial-sign alpha test",
    data.name = data.name,
    raw.statistic = stats::setNames(formula$Q, "Q.INST.alpha"),
    components = components,
    diagnostics = list(
      primary.endpoint = "inverse-norm K(r)=1/r",
      general.K.scope = paste(
        "not claimed feasible; use",
        "weighted_spatial_sign_alpha_oracle_test() for supplied scores"
      ),
      direction.residuals = "restricted factor residuals Y - B f",
      radial.residuals = "unrestricted OLS residuals Y - alpha - B f",
      trace.calibration = core$trace_split,
      zero.radius.policy = "error",
      weight.cap = "none",
      variance.repair = "none",
      numerical.repair = "none"
    ),
    call = call,
    estimate = ols$alpha
  )
}


#' Book Gaussian alpha Cauchy benchmark
#'
#' Implements the formula-complete Gaussian benchmark printed in Chapter 4,
#' kept distinct from [gaussian_alpha_combination_test()], which implements
#' the primary Feng--Lan--Liu--Ma Bonferroni rule. With unrestricted OLS
#' residual variances \eqn{\hat\sigma_i^2}, residual correlation \eqn{\hat R},
#' and intercept estimates \eqn{\hat\alpha}, this function uses
#' \deqn{T_{sum} = \frac{T\sum_i\hat\alpha_i^2/\hat\sigma_i^2-N}
#' {\{2\mathrm{tr}(\hat R^2)\}^{1/2}}}
#' and the maximum squared OLS intercept t statistic. Their upper-tail normal
#' and extreme-value p-values are combined by the equal-weight Cauchy rule.
#'
#' The manuscript denotes the trace term abstractly and does not prescribe a
#' unique finite-sample estimator for this Gaussian benchmark. Consequently,
#' `trace_R2` is required rather than silently replacing it by a residual
#' plug-in or a bias-corrected alternative.
#'
#' @inheritParams grs_alpha_test
#' @param trace_R2 Strictly positive supplied value of
#'   \eqn{\widehat{\mathrm{tr}(R^2)}}.
#' @return An upper-tail alpha-test object explicitly labelled as a book
#'   benchmark rather than a primary named procedure.
#' @references High-Dimensional Data Analysis for Elliptically Symmetric
#' Distributions, Chapter 4, equations for the Gaussian alpha benchmark.
#' @examples
#' f <- cbind(seq(-1, 1, length.out = 12))
#' y <- outer(seq_len(12), 1:3, function(i, j) cos(i / 2 + j))
#' book_gaussian_alpha_cauchy_test(y, trace_R2 = 3, factors = f)
#' @export
book_gaussian_alpha_cauchy_test <- function(returns, trace_R2,
                                             factors = NULL) {
  call <- match.call()
  data.name <- deparse1(substitute(returns))
  if (!is.numeric(trace_R2) || length(trace_R2) != 1L ||
      is.na(trace_R2) || !is.finite(trace_R2) || trace_R2 <= 0) {
    stop("`trace_R2` must be one strictly positive supplied finite number; no estimator or repair is selected silently.",
         call. = FALSE)
  }
  trace_R2 <- as.numeric(trace_R2)
  prepared <- .ch4_alpha_prepare(returns, factors, min_assets = 2L)
  ols <- .ch4_alpha_ols(prepared)
  residual.variance.scaled <- as.numeric(ols$native$sse) /
    as.numeric(ols$native$residual_df)
  residual.variance <- residual.variance.scaled * prepared$y.scale^2
  if (any(!is.finite(residual.variance.scaled)) ||
      any(residual.variance.scaled <= 0) ||
      any(!is.finite(residual.variance)) ||
      any(residual.variance <= 0)) {
    stop("Every OLS residual variance must be strictly positive and finite; no floor is used.",
         call. = FALSE)
  }
  quadratic <- prepared$T * sum(ols$alpha^2 / residual.variance)
  sum.statistic <- (quadratic - prepared$N) / sqrt(2 * trace_R2)
  maximum <- max(ols$native$t_squared)
  centered.maximum <- maximum - 2 * log(prepared$N) +
    log(log(prepared$N))
  if (!is.finite(sum.statistic) || !is.finite(centered.maximum)) {
    stop("A Gaussian alpha benchmark component is not finite.",
         call. = FALSE)
  }
  log.p.sum <- stats::pnorm(
    sum.statistic, lower.tail = FALSE, log.p = TRUE
  )
  log.q.sum <- stats::pnorm(
    sum.statistic, lower.tail = TRUE, log.p = TRUE
  )
  max.tail <- .ch4_alpha_gumbel_tail(centered.maximum)
  combined <- .ssmax_cauchy_combine_logtails(
    log.p.sum, log.q.sum,
    max.tail$log.p.value, max.tail$log.cdf
  )
  component.p <- c(
    sum = exp(log.p.sum),
    max = max.tail$p.value
  )
  .ch4c_alpha_result(
    statistic = stats::setNames(combined$statistic, "book Cauchy"),
    p.value = combined$p.value,
    method = "Book Gaussian alpha sum-max Cauchy benchmark",
    data.name = data.name,
    raw.statistic = c(
      sum = sum.statistic,
      maximum = maximum,
      centered.maximum = centered.maximum,
      Cauchy = combined$statistic
    ),
    components = list(
      alpha = ols$alpha,
      residual.variance = stats::setNames(
        residual.variance, prepared$asset.names
      ),
      residual.variance.scaled.coordinates = stats::setNames(
        residual.variance.scaled, prepared$asset.names
      ),
      residual.correlation = ols$native$residual_correlation,
      trace.R.squared = trace_R2,
      quadratic = quadratic,
      sum.statistic = sum.statistic,
      maximum.t.squared = maximum,
      centered.maximum = centered.maximum,
      p.values = component.p,
      log.p.values = c(sum = log.p.sum, max = max.tail$log.p.value),
      Cauchy.statistic.sign = combined$sign,
      log.absolute.Cauchy.statistic = combined$log.absolute,
      T = prepared$T,
      N = prepared$N,
      K = prepared$K
    ),
    diagnostics = list(
      scope = "formula-complete book Gaussian benchmark",
      primary.FLLM.combination = FALSE,
      book.quadratic.factor = "literal T, not h'h",
      trace.estimator = "supplied trace_R2; manuscript estimator left abstract",
      trace.R2.required = TRUE,
      combination = "equal-weight signed-log Cauchy",
      covariance.repair = "none"
    ),
    call = call,
    estimate = ols$alpha
  )
}


#' Fixed-dimensional conditional-factor Wald benchmark
#'
#' Evaluates the Chapter 4 supplied-estimate benchmark
#' \deqn{W=T\hat\delta^T\hat\Omega_\delta^{-1}\hat\delta,}
#' with a chi-squared reference distribution having `length(delta)` degrees
#' of freedom. This function does not estimate spline or kernel nuisances.
#'
#' @param delta Finite estimated conditional-alpha vector.
#' @param covariance Finite, exactly symmetric, strictly positive-definite
#'   covariance matrix for `delta`.
#' @param sample_size Positive integer \eqn{T}.
#' @return A strict Cholesky-based upper-tail Wald `htest` object.
#' @references Li, D. and Yang, L. (2011). Nonparametric tests of conditional
#' factor models. Ang, A. and Kristensen, D. (2012). Testing conditional
#' factor models.
#' @examples
#' conditional_factor_wald_test(c(0.1, -0.2), diag(c(1, 2)), 40)
#' @export
conditional_factor_wald_test <- function(delta, covariance, sample_size) {
  call <- match.call()
  data.name <- deparse1(substitute(delta))
  if (!is.numeric(delta) || !is.null(dim(delta))) {
    stop("`delta` must be a numeric vector; coercion is not used.",
         call. = FALSE)
  }
  delta.names <- names(delta)
  delta <- as.numeric(delta)
  if (!length(delta) || anyNA(delta) || any(!is.finite(delta))) {
    stop("`delta` must be a non-empty finite numeric vector.",
         call. = FALSE)
  }
  if (!is.matrix(covariance) || !is.numeric(covariance) ||
      !identical(dim(covariance), c(length(delta), length(delta))) ||
      anyNA(covariance) || any(!is.finite(covariance))) {
    stop("`covariance` must be a finite square numeric matrix matching `delta`.",
         call. = FALSE)
  }
  storage.mode(covariance) <- "double"
  if (!isTRUE(all(covariance == t(covariance)))) {
    stop("`covariance` must be exactly symmetric; it is not silently symmetrized.",
         call. = FALSE)
  }
  if (!is.numeric(sample_size)) {
    stop("`sample_size` must be numeric; coercion is not used.",
         call. = FALSE)
  }
  sample_size <- .ch4_alpha_positive(
    sample_size, "sample_size", integer = TRUE
  )
  root <- tryCatch(chol(covariance), error = function(e) NULL)
  if (is.null(root)) {
    stop("`covariance` must be strictly positive definite; no ridge, eigenvalue floor, or generalized inverse is used.",
         call. = FALSE)
  }
  solved <- backsolve(root, forwardsolve(t(root), delta))
  statistic <- sample_size * sum(delta * solved)
  if (!is.finite(statistic) || statistic < 0) {
    stop("The conditional-factor Wald statistic is invalid.",
         call. = FALSE)
  }
  df <- length(delta)
  answer <- list(
    statistic = stats::setNames(statistic, "W.T.cond"),
    parameter = stats::setNames(df, "df"),
    p.value = stats::pchisq(statistic, df = df, lower.tail = FALSE),
    method = "Fixed-dimensional conditional-factor Wald benchmark",
    data.name = data.name,
    alternative = "greater",
    estimate = stats::setNames(delta, if (is.null(delta.names)) {
      paste0("delta", seq_along(delta))
    } else {
      delta.names
    }),
    null.value = stats::setNames(rep(0, df), paste0("delta", seq_len(df))),
    components = list(
      covariance = covariance,
      solved.direction = solved,
      sample.size = sample_size
    ),
    diagnostics = list(
      scope = "supplied delta and covariance benchmark only",
      nuisance.estimation = "none",
      linear.solve = "strict Cholesky",
      numerical.repair = "none"
    ),
    call = call
  )
  class(answer) <- c("hd_conditional_alpha_test", "htest")
  answer
}


.ch4c_vector_workload <- function(n, p, q, order, B,
                                  max_kernel_evaluations) {
  max_kernel_evaluations <- .ch4ind_validate_integer(
    max_kernel_evaluations, "max_kernel_evaluations", minimum = 1L
  )
  log.base <- log(p) + log(q) + lchoose(n, order) + lgamma(order + 1)
  log.total <- log.base + log(B + 1)
  if (!is.finite(log.total) ||
      log.total > log(max_kernel_evaluations) + 1e-12) {
    stop(sprintf(
      paste0(
        "The exact high-order U-kernel workload exceeds ",
        "`max_kernel_evaluations` = %d; no approximation or silent ",
        "truncation is used."
      ),
      max_kernel_evaluations
    ), call. = FALSE)
  }
  list(
    maximum = max_kernel_evaluations,
    base = exp(log.base),
    total = exp(log.total),
    log.base = log.base,
    log.total = log.total
  )
}


.ch4c_reject_rank_ties <- function(x, name) {
  tied <- vapply(seq_len(ncol(x)), function(j) anyDuplicated(x[, j]) > 0L,
                 logical(1))
  if (any(tied)) {
    stop(sprintf(
      "`%s` contains an exact tie in column %d; the primary continuous-margin calibration is undefined.",
      name, which(tied)[1L]
    ), call. = FALSE)
  }
  invisible(x)
}


#' Wang--Liu--Feng degenerate rank-U vector-independence test
#'
#' Implements the exact high-order Hoeffding D, Blum--Kiefer--Rosenblatt R,
#' and Bergsma--Dassios--Yanagimoto tau-star examples in the primary paper.
#' Each coordinate-pair U-statistic is computed by enumerating every sample
#' subset and every kernel symmetrization. The maximum, centered sum of
#' squares, intrinsic X-row permutation variance, primary extreme-value tail,
#' and Fisher combination are then evaluated literally.
#'
#' This is an auditable exact implementation, not a scalable algorithm. Its
#' kernel-term workload is
#' \eqn{(B+1)pq {n \choose m}m!} for kernel order \eqn{m=5,6,4}.
#' Computation stops before drawing permutations when this exceeds
#' `max_kernel_evaluations`. No incomplete-U approximation is substituted.
#'
#' The primary tau-star example prints `(n-1)! 4! / n!`, which conflicts with
#' the paper's generic U-statistic definition and its own null second moment.
#' This implementation uses the coherent generic normalization
#' \eqn{1/{n \choose 4}} and records that decision in diagnostics.
#'
#' @param x Numeric \eqn{n} by \eqn{p} matrix.
#' @param y Numeric \eqn{n} by \eqn{q} matrix with matching rows.
#' @param measure One of `"hoeffding_d"`, `"bkr_r"`, or `"tau_star"`.
#' @param component One of `"fisher"`, `"max"`, or `"sum"`.
#' @param B Number of intrinsic X-row permutations; at least two.
#' @param seed Optional non-negative integer. An explicit seed is localized
#'   and preserves the caller's random-number state.
#' @param max_kernel_evaluations Positive integer upper bound on all exact
#'   symmetrized-kernel terms, including the observed statistic and all
#'   permutations.
#' @param keep_estimates Whether to retain the \eqn{p} by \eqn{q} matrix of
#'   observed U-statistics.
#' @param keep_permutation Whether to retain permutation sum statistics and
#'   permutation indices.
#' @return An `htest` object with exact workload and calibration diagnostics.
#' @references Wang, H., Liu, B. and Feng, L. (2026). Testing Independence
#' Between High-Dimensional Random Vectors Using Rank-Based Max-Sum Tests.
#' Scandinavian Journal of Statistics 53, 821--847.
#' \doi{10.1111/sjos.70063}
#' @examples
#' \donttest{
#' x <- cbind(c(1, 4, 2, 7, 3, 6, 5), c(7, 2, 5, 1, 6, 3, 4))
#' y <- cbind(c(2, 7, 4, 1, 6, 3, 5))
#' wang_liu_feng_vector_u_independence_test(
#'   x, y, measure = "tau_star", B = 7, seed = 11
#' )
#' }
#' @export
wang_liu_feng_vector_u_independence_test <- function(
    x, y, measure = c("hoeffding_d", "bkr_r", "tau_star"),
    component = c("fisher", "max", "sum"), B = 199L, seed = NULL,
    max_kernel_evaluations = 50000000L,
    keep_estimates = FALSE, keep_permutation = FALSE) {
  call <- match.call()
  x.name <- deparse1(substitute(x))
  y.name <- deparse1(substitute(y))
  measure <- match.arg(measure)
  component <- match.arg(component)
  order <- c(hoeffding_d = 5L, bkr_r = 6L, tau_star = 4L)[[measure]]
  x <- .ch4ind_validate_matrix(x, "x", min_rows = order, min_cols = 1L)
  y <- .ch4ind_validate_matrix(y, "y", min_rows = order, min_cols = 1L)
  if (nrow(x) != nrow(y)) {
    stop("`x` and `y` must have the same number of rows.", call. = FALSE)
  }
  comparisons <- ncol(x) * ncol(y)
  if (comparisons < 2L) {
    stop("The high-dimensional max calibration requires p * q >= 2.",
         call. = FALSE)
  }
  .ch4c_reject_rank_ties(x, "x")
  .ch4c_reject_rank_ties(y, "y")
  B <- .ch4ind_validate_integer(B, "B", minimum = 2L)
  if (!is.null(seed)) {
    seed <- .ch4ind_validate_integer(seed, "seed", minimum = 0L)
  }
  keep_estimates <- .ch4ind_validate_logical(
    keep_estimates, "keep_estimates"
  )
  keep_permutation <- .ch4ind_validate_logical(
    keep_permutation, "keep_permutation"
  )
  workload <- .ch4c_vector_workload(
    nrow(x), ncol(x), ncol(y), order, B,
    max_kernel_evaluations
  )
  draw.permutations <- function() {
    answer <- replicate(B, sample.int(nrow(x)), simplify = "matrix")
    storage.mode(answer) <- "integer"
    answer
  }
  permutations <- .ch4ind_with_local_seed(seed, draw.permutations())
  core <- cpp_ch4_completion_vector_u_core(
    x, y, match(measure, c("hoeffding_d", "bkr_r", "tau_star")) - 1L,
    permutations, workload$maximum, keep_estimates, keep_permutation
  )
  if (keep_permutation) {
    core$permutation_statistics <- as.numeric(
      core$permutation_statistics
    )
  }
  sum.z <- core$sum_statistic / sqrt(core$permutation_variance)
  log.p.sum <- stats::pnorm(sum.z, lower.tail = FALSE, log.p = TRUE)
  p.sum <- exp(log.p.sum)
  scale.denominator <- c(
    hoeffding_d = 30,
    bkr_r = 90,
    tau_star = 36
  )[[measure]]
  kappa <- 2.467
  max.gumbel <- pi^4 * (nrow(x) - 1) * core$maximum /
    scale.denominator - 2 * log(comparisons) +
    log(log(comparisons)) + pi^4 / 36
  max.tail <- .ch4ind_extreme_tail(
    max.gumbel, log(kappa) - 0.5 * log(pi)
  )
  fisher <- -2 * (log.p.sum + max.tail$log.p.value)
  p.fisher <- stats::pchisq(fisher, df = 4, lower.tail = FALSE)
  x.names <- .ch4ind_names(x, "x")
  y.names <- .ch4ind_names(y, "y")
  names(core$maximum_index) <- c("x", "y")
  core$maximum.names <- c(
    x.names[core$maximum_index[1L]],
    y.names[core$maximum_index[2L]]
  )
  if (keep_estimates) {
    dimnames(core$estimates) <- list(x.names, y.names)
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
    extreme.kappa = kappa,
    fisher = fisher,
    p.fisher = p.fisher,
    workload.limit = workload$maximum
  ))
  selected <- switch(
    component,
    fisher = list(statistic = c(Fisher = fisher), p.value = p.fisher),
    max = list(statistic = c(`max extreme-value` = max.gumbel),
               p.value = max.tail$p.value),
    sum = list(statistic = c(`sum Z` = sum.z), p.value = p.sum)
  )
  measure.label <- c(
    hoeffding_d = "Hoeffding D",
    bkr_r = "Blum-Kiefer-Rosenblatt R",
    tau_star = "Bergsma-Dassios-Yanagimoto tau-star"
  )[[measure]]
  .ch4ind_new_test(
    statistic = selected$statistic,
    p.value = selected$p.value,
    alternative = "greater",
    method = paste0(
      "Wang-Liu-Feng ", measure.label,
      " exact vector-independence test (", component, ")"
    ),
    data.name = paste(x.name, "and", y.name),
    raw.statistic = c(
      maximum.rank.U = core$maximum,
      centered.sum.squares = core$sum_statistic
    ),
    components = components,
    diagnostics = list(
      selected.component = component,
      measure = measure,
      tie.policy = "error; primary continuous-margin calibration",
      sum.calibration = "intrinsic X-row permutation sample variance",
      permutations = B,
      seed = seed,
      explicit.seed.preserves.caller.RNG.state = !is.null(seed),
      max.calibration = paste0(
        "primary degenerate-kernel extreme value; kappa approximately ",
        kappa
      ),
      combination = "Fisher chi-square with 4 degrees of freedom",
      exact.kernel.complexity = "(B+1) * p * q * choose(n,m) * m!",
      scalable.claim = FALSE,
      approximation = "none",
      simulation.adjustment = "not implemented",
      tau.star.normalization = paste(
        "generic primary U-statistic 1/choose(n,4);",
        "the conflicting (n-1)! example typo is not copied"
      ),
      numerical.repair = "none"
    ),
    call = call,
    estimate = c(maximum.rank.U = core$maximum),
    parameter = if (component == "fisher") c(df = 4) else NULL
  )
}






