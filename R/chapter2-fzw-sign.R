#' Feng--Zou--Wang two-sample multivariate-sign test
#'
#' Tests equality of two multivariate elliptical location vectors with the
#' scalar-invariant spatial-sign procedure of Feng, Zou, and Wang (2016).
#' Observations are rows and variables are columns.  If
#' \eqn{(\widehat\theta_{k,i},\widehat D_{k,i})} is the diagonal HR fit in
#' group \eqn{k} after deleting observation \eqn{i}, the raw statistic is
#' \deqn{R_n=-\frac{1}{n_1n_2}\sum_{i=1}^{n_1}\sum_{j=1}^{n_2}
#' U\{\widehat D_{1,i}^{-1/2}(X_{1i}-\widehat\theta_{2,j})\}^{\mathsf T}
#' U\{\widehat D_{2,j}^{-1/2}(X_{2j}-\widehat\theta_{1,i})\}.}
#' Thus, each sign uses its own group's leave-one-out scale but the other
#' group's leave-one-out location.  This crossed construction removes the
#' high-dimensional location-estimation bias without estimating and
#' subtracting a separate bias term.
#'
#' The p-value uses the feasible null calibration in Proposition 2 of the
#' original article, not the oracle variance for local alternatives.  Let
#' \eqn{\widetilde U_{ki}=U\{\widehat D_{k,i}^{-1/2}
#' (X_{ki}-\widehat\theta_{k,i})\}} and
#' \eqn{\widehat c_k=n_k^{-1}\sum_i
#' \|\widehat D_{k,i}^{-1/2}(X_{ki}-\widehat\theta_{k,i})\|^{-1}}.
#' With \eqn{\widehat D_1,\widehat D_2} denoting the two full-sample
#' diagonal HR fits, the three trace estimators are
#' \deqn{\widehat{\operatorname{tr}(A_1^2)}=
#' \frac{p^2\widehat c_2^2\widehat c_1^{-2}}{n_1(n_1-1)}
#' \sum_{k=1}^{n_1}\sum_{\ell\ne k}
#' (\widetilde U_{1\ell}^{\mathsf T}\widehat D_2^{-1/2}
#' \widehat D_1^{1/2}\widetilde U_{1k})^2,}
#' with the group indices reversed for
#' \eqn{\widehat{\operatorname{tr}(A_2^2)}}, and
#' \deqn{\widehat{\operatorname{tr}(A_3^{\mathsf T}A_3)}=
#' \frac{p^2}{n_1n_2}\sum_{\ell=1}^{n_1}\sum_{k=1}^{n_2}
#' (\widetilde U_{1\ell}^{\mathsf T}\widetilde U_{2k})^2.}
#' These are combined as
#' \deqn{\widehat\sigma_n^2=
#' \frac{2\widehat{\operatorname{tr}(A_1^2)}}{n_1(n_1-1)p^2}+
#' \frac{2\widehat{\operatorname{tr}(A_2^2)}}{n_2(n_2-1)p^2}+
#' \frac{4\widehat{\operatorname{tr}(A_3^{\mathsf T}A_3)}}{n_1n_2p^2}.}
#' The reported statistic is \eqn{Z=R_n/\widehat\sigma_n}; large positive
#' values reject, so the p-value is the upper standard-normal tail.
#'
#' Every full-sample and leave-one-out fit follows the published diagonal HR
#' recursion, initialized by the applicable sample mean and marginal sample
#' variances.  Convergence requires both estimating-equation residuals to be
#' at most `tol`.  A coincident training observation makes the inverse-radius
#' update undefined and is an error.  The spatial-sign convention is
#' \eqn{U(0)=0} in the numerator, but a zero own-sample leave-one-out radius is
#' also an error because its reciprocal is required by
#' \eqn{\widehat c_k}.  Nonconvergence, non-positive scales, and non-positive
#' feasible variance are reported without ridge, absolute-value, floor, or
#' perturbation repairs.
#'
#' Internally, both groups receive one common, safe coordinatewise affine
#' transformation.  This is algebraically neutral under the method's shift
#' and coordinatewise nonzero-scaling invariance and protects extreme finite
#' units.  Returned input-coordinate diagonal scales are canonically
#' normalized to have largest diagonal entry one; log diagonals retain ratios
#' that underflow in the display-scale version.
#'
#' The paper recommends a bootstrap calibration when the dimension is small
#' (it gives \eqn{p\le 50} as a guide) or grows at order \eqn{n^2} or faster.
#' This function implements the article's feasible asymptotic-normal test;
#' diagnostics flag those finite-sample regimes.  It does not silently switch
#' calibration or reproduce the paper's simulation procedure.
#'
#' @param x,y Numeric matrices or data frames with observations in rows and
#'   the same variables in columns.  Each group needs at least three rows,
#'   and each full/leave-one-out fit must have positive marginal variation.
#' @param tol Finite positive tolerance for both diagonal HR estimating
#'   equations.  The default is `1e-7`.
#' @param max_iter Positive integer maximum number of HR updates for every
#'   full-sample and leave-one-out fit.  Nonconvergence is an error.
#'
#' @return An object of class `c("hd_location_test", "htest")`.  In addition
#'   to the test result, `components` contains both full-sample fits, all
#'   leave-one-out fits and directions, inverse radii, bridge diagonals,
#'   ordered/cross pair contributions, three trace estimates, and all three
#'   variance contributions.  `diagnostics` contains iteration counts,
#'   equation residuals, minimum training radii, zero-sign counts, numerical
#'   scaling information, asymptotic-regime flags, and the no-repair policy.
#'
#' @references
#' Feng, L., Zou, C., and Wang, Z. (2016). Multivariate-sign-based
#' high-dimensional tests for the two-sample location problem. *Journal of
#' the American Statistical Association*, **111**(514), 721--735.
#' \doi{10.1080/01621459.2015.1035380}.
#'
#' @examples
#' set.seed(2016)
#' x <- matrix(stats::rt(48, df = 5), 8, 6)
#' y <- matrix(stats::rt(54, df = 5), 9, 6)
#' feng_zou_wang_two_sample_sign_test(x, y, tol = 1e-6)
#'
#' @export
feng_zou_wang_two_sample_sign_test <- function(
    x, y, tol = 1e-7, max_iter = 500L) {
  call <- match.call()
  x.name <- deparse1(substitute(x))
  y.name <- deparse1(substitute(y))
  x <- .as_data_matrix(x, "x", min_rows = 3L)
  y <- .as_data_matrix(y, "y", min_rows = 3L)
  .check_two_sample_variables(x, y)
  controls <- .validate_iteration_controls(tol, max_iter, zero_tol = 0)

  fit <- cpp_feng_zou_wang_two_sample_sign(
    x, y, controls$tol, controls$max_iter
  )
  n1 <- nrow(x)
  n2 <- nrow(y)
  p <- ncol(x)
  variable.names <- .hotelling_variable_names(x)
  row.names1 <- rownames(x)
  row.names2 <- rownames(y)
  if (is.null(row.names1)) {
    row.names1 <- paste0("x", seq_len(n1))
  }
  if (is.null(row.names2)) {
    row.names2 <- paste0("y", seq_len(n2))
  }

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
  z <- unname(fit$z)
  p.value <- stats::pnorm(z, lower.tail = FALSE)

  full.fit <- list(
    group1 = list(
      location = name.vector(fit$full_theta1),
      location.standardized = name.vector(fit$full_theta1_standardised),
      scale.diagonal = name.vector(fit$full_D1_input_canonical),
      log.scale.diagonal = name.vector(fit$full_log_D1_input_canonical),
      scale.diagonal.standardized = name.vector(
        fit$full_D1_standardised
      ),
      iterations = unname(fit$full_iterations1),
      location.equation.residual = unname(
        fit$full_location_residual1
      ),
      scale.equation.residual = unname(fit$full_scale_residual1),
      minimum.training.radius = unname(
        fit$full_minimum_training_radius1
      )
    ),
    group2 = list(
      location = name.vector(fit$full_theta2),
      location.standardized = name.vector(fit$full_theta2_standardised),
      scale.diagonal = name.vector(fit$full_D2_input_canonical),
      log.scale.diagonal = name.vector(fit$full_log_D2_input_canonical),
      scale.diagonal.standardized = name.vector(
        fit$full_D2_standardised
      ),
      iterations = unname(fit$full_iterations2),
      location.equation.residual = unname(
        fit$full_location_residual2
      ),
      scale.equation.residual = unname(fit$full_scale_residual2),
      minimum.training.radius = unname(
        fit$full_minimum_training_radius2
      )
    )
  )

  leaveout.fit <- list(
    group1 = list(
      location = name.group.matrix(fit$leaveout_theta1, 1L),
      location.standardized = name.group.matrix(
        fit$leaveout_theta1_standardised, 1L
      ),
      scale.diagonal = name.group.matrix(
        fit$leaveout_D1_input_canonical, 1L
      ),
      log.scale.diagonal = name.group.matrix(
        fit$leaveout_log_D1_input_canonical, 1L
      ),
      scale.diagonal.standardized = name.group.matrix(
        fit$leaveout_D1_standardised, 1L
      ),
      own.direction = name.group.matrix(fit$own_direction1, 1L),
      own.radius = name.group.vector(fit$own_radius1, 1L),
      own.inverse.radius = name.group.vector(
        fit$own_inverse_radius1, 1L
      ),
      iterations = name.group.vector(fit$leaveout_iterations1, 1L),
      location.equation.residual = name.group.vector(
        fit$leaveout_location_residual1, 1L
      ),
      scale.equation.residual = name.group.vector(
        fit$leaveout_scale_residual1, 1L
      ),
      minimum.training.radius = name.group.vector(
        fit$leaveout_minimum_training_radius1, 1L
      )
    ),
    group2 = list(
      location = name.group.matrix(fit$leaveout_theta2, 2L),
      location.standardized = name.group.matrix(
        fit$leaveout_theta2_standardised, 2L
      ),
      scale.diagonal = name.group.matrix(
        fit$leaveout_D2_input_canonical, 2L
      ),
      log.scale.diagonal = name.group.matrix(
        fit$leaveout_log_D2_input_canonical, 2L
      ),
      scale.diagonal.standardized = name.group.matrix(
        fit$leaveout_D2_standardised, 2L
      ),
      own.direction = name.group.matrix(fit$own_direction2, 2L),
      own.radius = name.group.vector(fit$own_radius2, 2L),
      own.inverse.radius = name.group.vector(
        fit$own_inverse_radius2, 2L
      ),
      iterations = name.group.vector(fit$leaveout_iterations2, 2L),
      location.equation.residual = name.group.vector(
        fit$leaveout_location_residual2, 2L
      ),
      scale.equation.residual = name.group.vector(
        fit$leaveout_scale_residual2, 2L
      ),
      minimum.training.radius = name.group.vector(
        fit$leaveout_minimum_training_radius2, 2L
      )
    )
  )

  all.iterations <- c(
    full.fit$group1$iterations,
    full.fit$group2$iterations,
    unname(leaveout.fit$group1$iterations),
    unname(leaveout.fit$group2$iterations)
  )
  all.location.residuals <- c(
    full.fit$group1$location.equation.residual,
    full.fit$group2$location.equation.residual,
    unname(leaveout.fit$group1$location.equation.residual),
    unname(leaveout.fit$group2$location.equation.residual)
  )
  all.scale.residuals <- c(
    full.fit$group1$scale.equation.residual,
    full.fit$group2$scale.equation.residual,
    unname(leaveout.fit$group1$scale.equation.residual),
    unname(leaveout.fit$group2$scale.equation.residual)
  )
  all.minimum.radii <- c(
    full.fit$group1$minimum.training.radius,
    full.fit$group2$minimum.training.radius,
    unname(leaveout.fit$group1$minimum.training.radius),
    unname(leaveout.fit$group2$minimum.training.radius)
  )
  n.total <- n1 + n2

  .new_hd_location_test(
    statistic = c(Z = z),
    p.value = p.value,
    method = paste(
      "Feng-Zou-Wang two-sample multivariate-sign test",
      "(leave-one-out statistic; feasible asymptotic normal calibration)"
    ),
    data.name = paste(x.name, "and", y.name),
    alternative = "two.sided",
    raw.statistic = c(R.n = unname(fit$R_n)),
    estimate = difference,
    null.value = stats::setNames(numeric(p), variable.names),
    null.distribution = list(
      family = "normal",
      parameters = c(mean = 0, sd = 1),
      exact = FALSE,
      tail = "upper",
      assumptions = paste(
        "Feng-Zou-Wang high-dimensional elliptical conditions (C1)-(C4);",
        "independent samples; converged diagonal HR fits"
      )
    ),
    variance = c(estimated = unname(fit$sigma2_hat)),
    components = list(
      mean.x = mean.x,
      mean.y = mean.y,
      difference = difference,
      R.n = unname(fit$R_n),
      sigma2.hat = unname(fit$sigma2_hat),
      sigma.hat = unname(fit$sigma_hat),
      c1.hat = unname(fit$c1_hat),
      c2.hat = unname(fit$c2_hat),
      trace.A1.squared.hat = unname(fit$trace_A1_squared_hat),
      trace.A2.squared.hat = unname(fit$trace_A2_squared_hat),
      trace.A3tA3.hat = unname(fit$trace_A3tA3_hat),
      variance.term1 = unname(fit$variance_term1),
      variance.term2 = unname(fit$variance_term2),
      variance.term3 = unname(fit$variance_term3),
      numerator.inner.product = name.cross.matrix(
        fit$numerator_inner_product
      ),
      numerator.contribution = name.cross.matrix(
        fit$numerator_contribution
      ),
      trace.A1.ordered.pair.squared = name.square.matrix(
        fit$trace1_pair_squared, 1L
      ),
      trace.A2.ordered.pair.squared = name.square.matrix(
        fit$trace2_pair_squared, 2L
      ),
      trace.A3.cross.pair.squared = name.cross.matrix(
        fit$trace3_pair_squared
      ),
      trace.A1.ordered.pair.squared.sum = unname(
        fit$trace1_pair_squared_sum
      ),
      trace.A2.ordered.pair.squared.sum = unname(
        fit$trace2_pair_squared_sum
      ),
      trace.A3.cross.pair.squared.sum = unname(
        fit$trace3_pair_squared_sum
      ),
      bridge.A1.diagonal.standardized = name.vector(
        fit$bridge1_diagonal
      ),
      bridge.A2.diagonal.standardized = name.vector(
        fit$bridge2_diagonal
      ),
      full.fit = full.fit,
      leaveout.fit = leaveout.fit,
      n1 = unname(fit$n1),
      n2 = unname(fit$n2),
      p = unname(fit$p),
      ordered.pair.count1 = unname(fit$ordered_pair_count1),
      ordered.pair.count2 = unname(fit$ordered_pair_count2),
      cross.pair.count = unname(fit$cross_pair_count)
    ),
    diagnostics = list(
      all.fits.converged = TRUE,
      tolerance = unname(fit$tolerance),
      max.iterations.allowed = unname(fit$max_iterations),
      full.fit = full.fit,
      leaveout.fit = leaveout.fit,
      iteration.summary = c(
        minimum = min(all.iterations),
        median = stats::median(all.iterations),
        mean = mean(all.iterations),
        maximum = max(all.iterations)
      ),
      maximum.location.equation.residual = max(all.location.residuals),
      maximum.scale.equation.residual = max(all.scale.residuals),
      minimum.training.radius = min(all.minimum.radii),
      numerator.zero.signs = name.cross.matrix(
        fit$numerator_zero_signs
      ),
      numerator.zero.sign.uses = unname(fit$numerator_zero_sign_uses),
      sign.at.zero = paste(
        "U(0) = 0 in the numerator; zero training or own leave-one-out",
        "radii are errors because inverse radii are required"
      ),
      numerator.structure = paste(
        "own-group leave-one-out scale and opposite-group leave-one-out",
        "location; minus the cross-sign inner product"
      ),
      variance.structure = paste(
        "Proposition 2: full-sample D bridges, own-group leave-one-out",
        "directions and inverse-radius c estimates"
      ),
      bias.handling = paste(
        "leave-one-out cross-centering makes the null bias asymptotically",
        "negligible; no additive bias estimate is subtracted"
      ),
      calibration = paste(
        "feasible null sigma.hat; the local-alternative oracle variance",
        "is not used in the statistic or p-value"
      ),
      applicability = list(
        p.at.most.50 = p <= 50L,
        p.at.least.total.n.squared = p >= n.total^2,
        total.n = n.total,
        note = paste(
          "The paper recommends bootstrap calibration when p is small",
          "(guide: p <= 50) or grows at order n^2 or faster; this API",
          "reports the feasible asymptotic-normal calibration only."
        )
      ),
      internal.scaling = paste(
        "one common safe shift and positive coordinatewise scale for both",
        "groups"
      ),
      internal.column.anchor = stats::setNames(
        as.numeric(fit$column_anchor), variable.names
      ),
      internal.column.log.scale = stats::setNames(
        as.numeric(fit$column_log_scale), variable.names
      ),
      subtraction.overflow.fallback.columns = unname(
        fit$subtraction_overflow_columns
      ),
      scale.normalization = paste(
        "every internal D uses the coherent unit-geometric-mean",
        "identification; returned input-coordinate D uses maximum diagonal",
        "one for display, with canonical log ratios also retained"
      ),
      computation = paste(
        "literal full-sample and observation-specific leave-one-out",
        "diagonal HR fits"
      ),
      regularization = "none",
      bias.repair = "none",
      variance.repair = "none"
    ),
    n = c(n1 = n1, n2 = n2),
    p = p,
    call = call
  )
}
