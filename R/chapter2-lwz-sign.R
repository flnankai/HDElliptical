#' Li--Wang--Zou simpler spatial-sign-based two-sample test
#'
#' Tests equality of two multivariate location vectors using the simplified
#' bias-corrected spatial-sign test (SST) of Li, Wang, and Zou (2016).  For
#' group \eqn{k=1,2}, let \eqn{(\widehat\theta_k,\widehat D_k)} be the
#' full-sample diagonal Hettmansperger--Randles fit satisfying
#' \deqn{n_k^{-1}\sum_i U\{\widehat D_k^{-1/2}
#' (X_{ki}-\widehat\theta_k)\}=0}
#' and
#' \deqn{p n_k^{-1}\mathrm{diag}\left(\sum_i
#' U_{ki}U_{ki}^{\mathsf T}\right)=I_p.}
#' Unlike the earlier Feng--Zou--Wang statistic, this method uses no leave-out
#' fit.  Its uncorrected statistic is
#' \deqn{T_n=-\frac{1}{n_1n_2}\sum_{i=1}^{n_1}\sum_{j=1}^{n_2}
#' U\{\widehat D_1^{-1/2}(X_{1i}-\widehat\theta_2)\}^{\mathsf T}
#' U\{\widehat D_2^{-1/2}(X_{2j}-\widehat\theta_1)\}.}
#'
#' Write
#' \eqn{\widetilde U_{ki}=U\{\widehat D_k^{-1/2}
#' (X_{ki}-\widehat\theta_k)\}},
#' \eqn{\widehat c_k=n_k^{-1}\sum_i
#' \|\widehat D_k^{-1/2}(X_{ki}-\widehat\theta_k)\|^{-1}},
#' \eqn{\widehat c_{nk}=\widehat c_{3-k}/\widehat c_k}, and
#' \eqn{\widehat D_n=\widehat D_1^{1/2}\widehat D_2^{1/2}}.  The feasible
#' bias terms in Proposition 1 are
#' \deqn{\widehat{\operatorname{tr}(A_k)}=
#' \widehat c_{nk}\operatorname{tr}
#' (\widehat D_k\widehat D_n^{-1}),\qquad
#' \widehat\mu_n=\frac{\widehat{\operatorname{tr}(A_1)}}{n_1p}+
#' \frac{\widehat{\operatorname{tr}(A_2)}}{n_2p}.}
#' The within-group squared-trace estimators are
#' \deqn{\widehat{\operatorname{tr}(A_k^2)}=
#' \frac{p^2\widehat c_{nk}^2}{n_k(n_k-1)}
#' \sum_{a=1}^{n_k}\sum_{b\ne a}
#' \{\widetilde U_{kb}^{\mathsf T}\widehat D_k
#' \widehat D_n^{-1}\widetilde U_{ka}\}^2,}
#' and the cross trace is
#' \deqn{\widehat{\operatorname{tr}(A_3^{\mathsf T}A_3)}=
#' \frac{p^2}{n_1n_2}\sum_{i=1}^{n_1}\sum_{j=1}^{n_2}
#' (\widetilde U_{1i}^{\mathsf T}\widetilde U_{2j})^2.}
#' The published feasible variance is
#' \deqn{\widehat\sigma_n^2=
#' \frac{2\widehat{\operatorname{tr}(A_1^2)}}{n_1^2p^2}+
#' \frac{2\widehat{\operatorname{tr}(A_2^2)}}{n_2^2p^2}+
#' \frac{4\widehat{\operatorname{tr}(A_3^{\mathsf T}A_3)}}
#' {n_1n_2p^2}.}
#' In particular, the first two outer denominators are \eqn{n_k^2}; only the
#' ordered-pair trace estimators use \eqn{n_k(n_k-1)}.  The reported statistic
#' is \eqn{Z=(T_n-\widehat\mu_n)/\widehat\sigma_n}.
#'
#' The scientific alternative is two-sided location inequality, but the SST
#' is quadratic in the location difference.  Consequently, the original
#' article rejects for large positive \eqn{Z}, and this function reports the
#' upper standard-normal tail.
#'
#' Each full-sample fit starts at its sample mean and marginal sample
#' variances and iterates the two published estimating equations.  The
#' diagonal equations identify scale only up to a common multiplier; a unit
#' geometric mean is imposed internally.  All bias and variance combinations
#' are invariant to that identification.  `iteration.stable` means that both
#' estimating-equation residuals are at most `tol`.  With `strict = TRUE`, an
#' unstable fit is an error; otherwise the last iterate is returned with a
#' warning and full diagnostics.
#'
#' The spatial-sign convention is \eqn{U(0)=0} for numerator terms.  A zero
#' training radius is nevertheless an error because the recursive location
#' update and \eqn{\widehat c_k} require inverse radii.  No radius perturbation,
#' weight cap, ridge, pseudoinverse, absolute-value repair, numerical floor,
#' or variance substitution is applied.
#'
#' @param x,y Numeric matrices or data frames with observations in rows and
#'   the same variables in columns.  The literal feasible formulas require at
#'   least two observations in each group; every group-specific marginal
#'   sample variance must be positive.
#' @param alpha Significance level for the returned upper-tail rejection rule.
#' @param tol Finite positive tolerance for both full-sample diagonal HR
#'   estimating equations.
#' @param max_iter Positive maximum number of HR updates per group.
#' @param zero_tol Non-negative tolerance below which a standardized radius is
#'   treated as zero.  The default detects exact zeros only.
#' @param strict If `TRUE`, error unless both fits satisfy the estimating
#'   equations within `tol`; if `FALSE`, warn and return the last iterates.
#'
#' @return An object of class `c("hd_location_test", "htest")`.  Its
#'   `components` field contains the uncorrected and bias-corrected statistics,
#'   bias and variance estimates, radial moments, three trace estimates,
#'   ordered/cross-pair kernels, both full-sample fits, signs, radii, and scale
#'   bridges.  `diagnostics` records convergence, denominator conventions,
#'   upper-tail calibration, zero-sign use, numerical scaling, and the
#'   no-repair contract.
#'
#' @references
#' Li, Y., Wang, Z., and Zou, C. (2016). A simpler spatial-sign-based
#' two-sample test for high-dimensional data. *Journal of Multivariate
#' Analysis*, **149**, 192--198.
#' \doi{10.1016/j.jmva.2016.04.004}.
#'
#' @examples
#' set.seed(2016)
#' x <- matrix(stats::rt(32, df = 5), 8, 4)
#' y <- matrix(stats::rt(36, df = 5), 9, 4)
#' li_wang_zou_two_sample_sign_test(x, y, tol = 1e-6)
#'
#' @export
li_wang_zou_two_sample_sign_test <- function(
    x, y, alpha = 0.05, tol = 1e-7, max_iter = 500L,
    zero_tol = 0, strict = TRUE) {
  call <- match.call()
  x.name <- deparse1(substitute(x))
  y.name <- deparse1(substitute(y))
  x <- .as_data_matrix(x, "x", min_rows = 2L)
  y <- .as_data_matrix(y, "y", min_rows = 2L)
  .check_two_sample_variables(x, y)
  alpha <- as.numeric(alpha)
  if (length(alpha) != 1L || is.na(alpha) || !is.finite(alpha) ||
      alpha <= 0 || alpha >= 1) {
    stop("`alpha` must be one finite number strictly between zero and one.",
         call. = FALSE)
  }
  controls <- .validate_iteration_controls(tol, max_iter, zero_tol)
  if (!is.logical(strict) || length(strict) != 1L || is.na(strict)) {
    stop("`strict` must be `TRUE` or `FALSE`.", call. = FALSE)
  }

  fit <- cpp_li_wang_zou_two_sample_sign(
    x, y, controls$tol, controls$max_iter, controls$zero_tol
  )
  if (!isTRUE(fit$all_iteration_stable)) {
    message <- sprintf(
      paste0(
        "Li-Wang-Zou SST did not satisfy the estimating equations for ",
        "%d full-sample fit(s) within `max_iter = %d` updates. No ridge ",
        "or iterate substitution was applied."
      ),
      as.integer(fit$stability_failures), controls$max_iter
    )
    if (isTRUE(strict)) {
      stop(message, call. = FALSE)
    }
    warning(message, call. = FALSE)
  }

  n1 <- nrow(x)
  n2 <- nrow(y)
  p <- ncol(x)
  variable.names <- .hotelling_variable_names(x)
  row.names1 <- rownames(x)
  row.names2 <- rownames(y)
  if (is.null(row.names1)) row.names1 <- paste0("x", seq_len(n1))
  if (is.null(row.names2)) row.names2 <- paste0("y", seq_len(n2))

  name.vector <- function(value) {
    stats::setNames(as.numeric(value), variable.names)
  }
  name.group.vector <- function(value, group) {
    labels <- if (group == 1L) row.names1 else row.names2
    stats::setNames(as.numeric(value), labels)
  }
  name.group.matrix <- function(value, group) {
    value <- as.matrix(value)
    labels <- if (group == 1L) row.names1 else row.names2
    dimnames(value) <- list(labels, variable.names)
    value
  }
  name.square.matrix <- function(value, group) {
    value <- as.matrix(value)
    labels <- if (group == 1L) row.names1 else row.names2
    dimnames(value) <- list(labels, labels)
    value
  }
  name.cross.matrix <- function(value) {
    value <- as.matrix(value)
    dimnames(value) <- list(row.names1, row.names2)
    value
  }

  mean.x <- name.vector(fit$sample_mean1)
  mean.y <- name.vector(fit$sample_mean2)
  difference <- stats::setNames(mean.x - mean.y, variable.names)
  z <- as.numeric(fit$z)
  critical.value <- stats::qnorm(1 - alpha)
  rejection <- isTRUE(z > critical.value)

  full.fit <- list(
    group1 = list(
      location = name.vector(fit$location1),
      location.standardized = name.vector(fit$location1_standardized),
      scale.diagonal = name.vector(fit$diagonal1_input_canonical),
      log.scale.diagonal = name.vector(
        fit$log_diagonal1_input_canonical
      ),
      scale.diagonal.standardized = name.vector(
        fit$diagonal1_standardized
      ),
      diagnostics = fit$fit_diagnostics1
    ),
    group2 = list(
      location = name.vector(fit$location2),
      location.standardized = name.vector(fit$location2_standardized),
      scale.diagonal = name.vector(fit$diagonal2_input_canonical),
      log.scale.diagonal = name.vector(
        fit$log_diagonal2_input_canonical
      ),
      scale.diagonal.standardized = name.vector(
        fit$diagonal2_standardized
      ),
      diagnostics = fit$fit_diagnostics2
    )
  )

  .new_hd_location_test(
    statistic = c(Z = z),
    p.value = stats::pnorm(z, lower.tail = FALSE),
    method = paste(
      "Li-Wang-Zou simpler spatial-sign-based two-sample test",
      "(SST; feasible bias correction)"
    ),
    data.name = paste(x.name, "and", y.name),
    alternative = "two.sided",
    raw.statistic = c(T.SST = as.numeric(fit$T)),
    estimate = difference,
    null.value = stats::setNames(numeric(p), variable.names),
    null.distribution = list(
      family = "normal",
      parameters = c(mean = 0, sd = 1),
      exact = FALSE,
      tail = "upper",
      assumptions = paste(
        "Li-Wang-Zou high-dimensional elliptical assumptions (A1)-(A4);",
        "independent samples; stable full-sample diagonal HR fits"
      )
    ),
    variance = c(estimated = as.numeric(fit$variance)),
    components = list(
      mean.x = mean.x,
      mean.y = mean.y,
      difference = difference,
      T.SST = as.numeric(fit$T),
      bias.hat = as.numeric(fit$bias),
      T.minus.bias = as.numeric(fit$centered),
      sigma2.hat = as.numeric(fit$variance),
      sigma.hat = as.numeric(fit$standard_error),
      c1.hat = as.numeric(fit$c1_hat),
      c2.hat = as.numeric(fit$c2_hat),
      c.n1.hat = as.numeric(fit$c_ratio1_hat),
      c.n2.hat = as.numeric(fit$c_ratio2_hat),
      trace.A1.hat = as.numeric(fit$trace_A1_hat),
      trace.A2.hat = as.numeric(fit$trace_A2_hat),
      trace.A1.squared.hat = as.numeric(fit$trace_A1_squared_hat),
      trace.A2.squared.hat = as.numeric(fit$trace_A2_squared_hat),
      trace.A3tA3.hat = as.numeric(fit$trace_A3tA3_hat),
      variance.term1 = as.numeric(fit$variance_term1),
      variance.term2 = as.numeric(fit$variance_term2),
      variance.term3 = as.numeric(fit$variance_term3),
      numerator.inner.product = name.cross.matrix(
        fit$numerator_inner_product
      ),
      numerator.contribution = name.cross.matrix(
        fit$numerator_contribution
      ),
      trace.A1.ordered.pair.squared = name.square.matrix(
        fit$trace_A1_ordered_pair_squared, 1L
      ),
      trace.A2.ordered.pair.squared = name.square.matrix(
        fit$trace_A2_ordered_pair_squared, 2L
      ),
      trace.A3.cross.pair.squared = name.cross.matrix(
        fit$trace_A3_cross_pair_squared
      ),
      bridge.A1.diagonal.standardized = name.vector(
        fit$bridge1_diagonal_standardized
      ),
      bridge.A2.diagonal.standardized = name.vector(
        fit$bridge2_diagonal_standardized
      ),
      own.direction1 = name.group.matrix(fit$own_direction1, 1L),
      own.direction2 = name.group.matrix(fit$own_direction2, 2L),
      own.radius1 = name.group.vector(fit$own_radius1, 1L),
      own.radius2 = name.group.vector(fit$own_radius2, 2L),
      own.inverse.radius1 = name.group.vector(
        fit$own_inverse_radius1, 1L
      ),
      own.inverse.radius2 = name.group.vector(
        fit$own_inverse_radius2, 2L
      ),
      full.fit = full.fit,
      n1 = as.numeric(fit$n1),
      n2 = as.numeric(fit$n2),
      p = as.numeric(fit$p),
      ordered.pair.count1 = as.numeric(fit$ordered_pair_count1),
      ordered.pair.count2 = as.numeric(fit$ordered_pair_count2),
      cross.pair.count = as.numeric(fit$cross_pair_count)
    ),
    diagnostics = list(
      iteration.stable = isTRUE(fit$all_iteration_stable),
      stability.failures = as.integer(fit$stability_failures),
      convergence.basis = paste(
        "maximum location/diagonal estimating-equation residual;",
        "reported separately for each full-sample fit"
      ),
      full.fit = full.fit,
      score.residual = max(
        as.numeric(fit$fit_diagnostics1$score.residual),
        as.numeric(fit$fit_diagnostics2$score.residual)
      ),
      relative.update = max(
        as.numeric(fit$fit_diagnostics1$relative.update),
        as.numeric(fit$fit_diagnostics2$relative.update)
      ),
      minimum.residual.distance = min(
        as.numeric(fit$fit_diagnostics1$minimum.residual.distance),
        as.numeric(fit$fit_diagnostics2$minimum.residual.distance)
      ),
      fit.structure = paste(
        "exactly two full-sample group fits; no leave-out fit is used"
      ),
      numerator.structure = paste(
        "minus the cross inner product; group-1 observations use the",
        "group-2 full-sample location and group-1 scale, with group",
        "indices reversed for the second sign"
      ),
      bias.correction = paste(
        "T.SST minus trace(A1)/(n1*p) minus trace(A2)/(n2*p)"
      ),
      trace.denominators = paste(
        "within-group ordered distinct pairs use n_k*(n_k-1);",
        "the cross trace uses n1*n2"
      ),
      variance.denominators = paste(
        "published feasible outer factors use n1^2 and n2^2;",
        "not n_k*(n_k-1)"
      ),
      scientific.alternative = "two-sided location inequality",
      calibration.tail = "upper standard-normal tail for quadratic SST",
      rejection = list(
        alpha = alpha,
        critical.value = critical.value,
        rule = "reject when Z > qnorm(1 - alpha)",
        reject = rejection
      ),
      numerator.zero.sign.uses = as.numeric(
        fit$numerator_zero_sign_uses
      ),
      sign.at.zero = paste(
        "U(0) = 0 in numerator terms; zero training radii are errors",
        "because inverse radii are required"
      ),
      radial.moment.scale = paste(
        "c1.hat and c2.hat use the unit-geometric-mean identification",
        "of each D; their ratios and all test quantities are identified"
      ),
      tolerance = controls$tol,
      max.iterations.allowed = controls$max_iter,
      zero.tolerance = controls$zero_tol,
      internal.scaling = paste(
        "one common safe shift and positive coordinatewise scale for both",
        "groups"
      ),
      internal.column.log.scale = stats::setNames(
        as.numeric(fit$column_log_scale), variable.names
      ),
      subtraction.overflow.fallback.columns = as.numeric(
        fit$subtraction_overflow_columns
      ),
      scale.identification = paste(
        "each internal D has unit geometric mean; displayed input D is",
        "canonicalized to maximum one and log ratios are retained"
      ),
      computation = "O(n^2 p) pair aggregation after two full-sample fits",
      regularization = "none",
      radius.perturbation = "none",
      weight.cap = "none",
      pseudoinverse = "none",
      bias.repair = "none",
      variance.repair = "none"
    ),
    n = c(n1 = n1, n2 = n2),
    p = p,
    call = call
  )
}
