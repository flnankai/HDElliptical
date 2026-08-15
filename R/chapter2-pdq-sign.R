.pdq_probability <- function(value, name) {
  value <- as.numeric(value)
  if (length(value) != 1L || is.na(value) || !is.finite(value) ||
      value <= 0 || value >= 1) {
    stop(sprintf("`%s` must be one finite number strictly between zero and one.",
                 name), call. = FALSE)
  }
  value
}

.pdq_positive_integer <- function(value, name) {
  if (!is.numeric(value) || length(value) != 1L || is.na(value) ||
      !is.finite(value) || value < 1 || value != floor(value) ||
      value > .Machine$integer.max) {
    stop(sprintf("`%s` must be a positive integer.", name), call. = FALSE)
  }
  as.integer(value)
}

.pdq_logical_scalar <- function(value, name) {
  if (!is.logical(value) || length(value) != 1L || is.na(value)) {
    stop(sprintf("`%s` must be `TRUE` or `FALSE`.", name), call. = FALSE)
  }
  value
}

.pdq_seed <- function(seed) {
  if (is.null(seed)) return(NULL)
  seed <- as.numeric(seed)
  if (length(seed) != 1L || is.na(seed) || !is.finite(seed) ||
      seed < 0 || seed > 2^32 - 1 || seed != floor(seed)) {
    stop("`seed` must be `NULL` or one integer-valued number in " %+%
           "[0, 2^32 - 1].", call. = FALSE)
  }
  seed
}

#' Feng--Wang pairwise-difference-quantile spatial-sign test
#'
#' Tests equality of two high-dimensional elliptical location vectors using
#' the pairwise-difference-quantile (PDQ) spatial-sign method of Feng and Wang
#' (2026). Observations are rows and variables are columns. The method allows
#' arbitrary within-group correlation and unequal scatter matrices.
#'
#' For group \eqn{k} and coordinate \eqn{j}, let
#' \eqn{M_k=n_k(n_k-1)/2} and form the \eqn{M_k} unordered absolute
#' pairwise differences. The PDQ scale
#' is the exact empirical U-quantile
#' \deqn{q_{kj}=\inf\{t:M_k^{-1}\sum_{i<l}
#' 1(|X_{kij}-X_{klj}|\leq t)\geq a\},}
#' where \eqn{a} is `quantile_prob`. Thus the implementation selects sorted
#' difference number \eqn{\lceil aM_k\rceil}; it does not use an interpolated
#' sample quantile. The diagonal standardizer has entries \eqn{q_{kj}^2}.
#'
#' After computing each full-sample spatial median under its PDQ
#' standardizer, the observed full-sample statistic is the cross-centred
#' quantity
#' \deqn{\widehat R=-{1\over n_1n_2}\sum_{i,j}
#' U\{\widehat D_1^{-1/2}(X_{1i}-\widehat\theta_2)\}^{T}
#' U\{\widehat D_2^{-1/2}(X_{2j}-\widehat\theta_1)\}.}
#' This is not replaced by the fitted-sign quadratic expansion. From fitted
#' within-group signs the function constructs the published matrices
#' \eqn{\widehat K_1,\widehat K_2,\widehat K_3} and subtracts the empirical
#' diagonal term
#' \deqn{\widehat b=\sum_{k=1}^2 n_k^{-2}\sum_i
#' \widehat S_{ki}^{T}\widehat K_k\widehat S_{ki}.}
#' The reported statistic is \eqn{T_{PDQ}=\widehat R-\widehat b}; large values
#' reject the two-sided location null.
#'
#' Calibration fixes all nuisance estimates. For each bootstrap draw,
#' independent Rademacher multipliers are applied to the fitted within-group
#' signs, the published quadratic form \eqn{Q^*} is evaluated, and
#' \eqn{T^*=Q^*-\widehat b} is recorded. The paper's primary rule rejects
#' strictly when \eqn{T_{PDQ}} exceeds the lower empirical
#' \eqn{(1-\mathrm{level})} quantile of the bootstrap values. Because the paper
#' does not prescribe a finite-`B` p-value, the `htest` `p.value` uses the
#' explicitly documented conservative convention
#' \deqn{p_{MC}={1+\#\{T^*\geq T_{PDQ}\}\over B+1}.}
#' Ties enter the upper tail. The primary strict-quantile decision and the
#' auxiliary p-value decision are both returned and can differ on the finite
#' bootstrap grid.
#'
#' The spatial median is computed with a modified Weiszfeld recursion and is
#' accepted when its normalized subgradient residual is at most `tol`.
#' `strict = TRUE` makes non-convergence an error; `strict = FALSE` warns and
#' continues from the last iterate when every published downstream quantity
#' remains defined. A fitted zero residual is always an error because the
#' definition of \eqn{\widehat G_k} requires its inverse radius. Zero
#' cross-centred residuals use the paper's convention \eqn{U(0)=0} and are
#' counted. Zero PDQ scales, singular \eqn{\widehat G_k}, and non-positive
#' diagonal-deleted bootstrap variance fail explicitly. No ridge, scale
#' floor, generalized inverse, absolute-value correction, or random
#' perturbation is used.
#'
#' A common coordinatewise affine preconditioning is used internally to
#' protect calculations with extreme but finite units. It is algebraically
#' neutral under the method's common translation and common nonzero
#' coordinatewise scaling invariances. Actual input-unit PDQ diagonals are
#' returned when representable; logarithms and a canonical maximum-one
#' version retain scale ratios when squaring would overflow or underflow.
#'
#' @param x,y Numeric matrices or data frames containing the two independent
#'   samples. They must have the same variables, at least three rows apiece,
#'   and at least two columns.
#' @param quantile_prob Probability for the coordinatewise PDQ U-quantile.
#'   The default is `0.25`; the paper also mentions `0.5`.
#' @param level Test level strictly between zero and one. This is the paper's
#'   rejection probability (denoted beta there), not `quantile_prob`.
#' @param B Positive number of Rademacher bootstrap draws.
#' @param seed `NULL`, or an integer-valued counter-generator seed in
#'   `[0, 2^32 - 1]`. `NULL` draws and records one seed from R's RNG. An
#'   explicit seed makes the bootstrap reproducible without changing R's RNG
#'   state.
#' @param keep_bootstrap If `TRUE`, retain all bootstrap statistics and the
#'   corresponding Rademacher multiplier matrix.
#' @param tol Finite positive spatial-median subgradient tolerance.
#' @param max_iter Positive maximum number of modified Weiszfeld updates.
#' @param strict Whether spatial-median non-convergence is an error (`TRUE`)
#'   or a warning followed by use of the last defined iterate (`FALSE`).
#'
#' @return An object of class `c("hd_location_test", "htest")`. Its
#'   `components` retain the PDQ quantiles and exact order-statistic ranks,
#'   both spatial-median fits, fitted and cross-centred signs, the
#'   \eqn{\widehat\Omega_k}, \eqn{\widehat G_k}, bridge and K matrices,
#'   observed and fitted quadratic decompositions, diagonal-deletion kernels,
#'   bootstrap variance, tail count, and optional draws. `diagnostics` records
#'   both finite-`B` decisions, convergence, zeros, conditioning, numerical
#'   representation, seed, and the no-repair policy.
#'
#' @references
#' Feng, L. and Wang, H. (2026). High-dimensional two-sample test for
#' elliptical symmetry distribution. arXiv:2605.03265.
#' \doi{10.48550/arXiv.2605.03265}.
#'
#' @examples
#' v <- rbind(c(1, 2, 3), c(2, -3, 1), c(-4, 1, 2), c(3, 4, -2))
#' x <- rbind(v, -v)
#' y <- rbind(1.1 * v, -1.1 * v) + rep(c(0.2, -0.1, 0.15), each = 8)
#' feng_wang_pdq_two_sample_test(x, y, B = 19, seed = 2605)
#'
#' @export
feng_wang_pdq_two_sample_test <- function(
    x, y, quantile_prob = 0.25, level = 0.05, B = 999L,
    seed = NULL, keep_bootstrap = FALSE, tol = 1e-8,
    max_iter = 1000L, strict = TRUE) {
  call <- match.call()
  x.name <- deparse1(substitute(x))
  y.name <- deparse1(substitute(y))
  x <- .as_data_matrix(x, "x", min_rows = 3L)
  y <- .as_data_matrix(y, "y", min_rows = 3L)
  .check_two_sample_variables(x, y)
  if (ncol(x) < 2L) {
    stop("PDQ requires at least two variables because the empirical " %+%
           "G matrices are singular when p = 1.", call. = FALSE)
  }

  quantile_prob <- .pdq_probability(quantile_prob, "quantile_prob")
  level <- .pdq_probability(level, "level")
  B <- .pdq_positive_integer(B, "B")
  seed.requested <- .pdq_seed(seed)
  keep_bootstrap <- .pdq_logical_scalar(keep_bootstrap, "keep_bootstrap")
  strict <- .pdq_logical_scalar(strict, "strict")
  controls <- .validate_iteration_controls(tol, max_iter, zero_tol = 0)
  seed.used <- if (is.null(seed.requested)) {
    as.numeric(sample.int(.Machine$integer.max, 1L) - 1L)
  } else {
    seed.requested
  }

  fit <- cpp_feng_wang_pdq_two_sample(
    x, y, quantile_prob, level, B, seed.used, keep_bootstrap,
    controls$tol, controls$max_iter
  )
  converged <- c(group1 = isTRUE(fit$median_converged1),
                 group2 = isTRUE(fit$median_converged2))
  if (!all(converged)) {
    message <- paste0(
      "The PDQ spatial median did not meet `tol` in ",
      paste(names(converged)[!converged], collapse = " and "),
      " within `max_iter`; no iterate repair was applied."
    )
    if (strict) {
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
    rows <- if (group == 1L) row.names1 else row.names2
    stats::setNames(as.numeric(value), rows)
  }
  name.group.matrix <- function(value, group) {
    rows <- if (group == 1L) row.names1 else row.names2
    value <- as.matrix(value)
    dimnames(value) <- list(rows, variable.names)
    value
  }
  name.variable.matrix <- function(value) {
    value <- as.matrix(value)
    dimnames(value) <- list(variable.names, variable.names)
    value
  }
  name.within.matrix <- function(value, group) {
    rows <- if (group == 1L) row.names1 else row.names2
    value <- as.matrix(value)
    dimnames(value) <- list(rows, rows)
    value
  }
  name.cross.matrix <- function(value) {
    value <- as.matrix(value)
    dimnames(value) <- list(row.names1, row.names2)
    value
  }

  bootstrap.values <- fit$bootstrap
  multiplier.matrix <- fit$multipliers
  if (!is.null(bootstrap.values)) {
    bootstrap.values <- stats::setNames(
      as.numeric(bootstrap.values), paste0("b", seq_len(B))
    )
  }
  if (!is.null(multiplier.matrix)) {
    multiplier.matrix <- as.matrix(multiplier.matrix)
    dimnames(multiplier.matrix) <- list(
      paste0("b", seq_len(B)),
      c(paste0("epsilon.x", seq_len(n1)),
        paste0("epsilon.y", seq_len(n2)))
    )
  }

  observed <- list(
    R.PDQ = as.numeric(fit$R_hat),
    b.hat = as.numeric(fit$b_hat),
    T.PDQ = as.numeric(fit$statistic),
    Q.fitted = as.numeric(fit$Q_hat_fitted),
    Q.fitted.minus.b = as.numeric(fit$fitted_diagonal_deleted),
    identity = "T.PDQ = cross-centred R.PDQ - empirical diagonal b.hat",
    fitted.quadratic.is.not.observed = TRUE
  )
  scales <- list(
    group1 = list(
      pair.count = as.numeric(fit$quantile_pair_count1),
      order.index = as.numeric(fit$quantile_rank1),
      quantile.input = name.vector(fit$quantile_input1),
      log.quantile.input = name.vector(fit$log_quantile_input1),
      quantile.preconditioned = name.vector(
        fit$quantile_preconditioned1
      ),
      quantile.working = name.vector(fit$quantile_working1),
      D.diagonal.input = name.vector(fit$D_input1),
      D.diagonal.input.canonical = name.vector(
        fit$D_input_canonical1
      ),
      D.diagonal.working = name.vector(fit$D_working1)
    ),
    group2 = list(
      pair.count = as.numeric(fit$quantile_pair_count2),
      order.index = as.numeric(fit$quantile_rank2),
      quantile.input = name.vector(fit$quantile_input2),
      log.quantile.input = name.vector(fit$log_quantile_input2),
      quantile.preconditioned = name.vector(
        fit$quantile_preconditioned2
      ),
      quantile.working = name.vector(fit$quantile_working2),
      D.diagonal.input = name.vector(fit$D_input2),
      D.diagonal.input.canonical = name.vector(
        fit$D_input_canonical2
      ),
      D.diagonal.working = name.vector(fit$D_working2)
    )
  )
  medians <- list(
    group1 = list(
      location = name.vector(fit$location1),
      location.preconditioned = name.vector(
        fit$location_preconditioned1
      ),
      location.standardized = name.vector(fit$location_standardized1),
      standardized.data = name.group.matrix(fit$standardized_data1, 1L),
      iterations = as.integer(fit$median_iterations1),
      converged = isTRUE(fit$median_converged1),
      relative.change = as.numeric(fit$median_relative_change1),
      equation.residual = as.numeric(fit$median_equation_residual1),
      objective = as.numeric(fit$median_objective1),
      zero.residual.count = as.numeric(fit$median_zero_count1)
    ),
    group2 = list(
      location = name.vector(fit$location2),
      location.preconditioned = name.vector(
        fit$location_preconditioned2
      ),
      location.standardized = name.vector(fit$location_standardized2),
      standardized.data = name.group.matrix(fit$standardized_data2, 2L),
      iterations = as.integer(fit$median_iterations2),
      converged = isTRUE(fit$median_converged2),
      relative.change = as.numeric(fit$median_relative_change2),
      equation.residual = as.numeric(fit$median_equation_residual2),
      objective = as.numeric(fit$median_objective2),
      zero.residual.count = as.numeric(fit$median_zero_count2)
    )
  )

  fitted.signs1 <- name.group.matrix(fit$fitted_signs1, 1L)
  fitted.signs2 <- name.group.matrix(fit$fitted_signs2, 2L)
  cross.signs1 <- name.group.matrix(fit$cross_signs1, 1L)
  cross.signs2 <- name.group.matrix(fit$cross_signs2, 2L)
  matrices <- list(
    Omega1 = name.variable.matrix(fit$Omega1),
    Omega2 = name.variable.matrix(fit$Omega2),
    G1 = name.variable.matrix(fit$G1),
    G2 = name.variable.matrix(fit$G2),
    G.inverse1 = name.variable.matrix(fit$G_inverse1),
    G.inverse2 = name.variable.matrix(fit$G_inverse2),
    A12.diagonal = name.vector(fit$A12_diagonal),
    A21.diagonal = name.vector(fit$A21_diagonal),
    A12.diagonal.input = name.vector(fit$A12_input_diagonal),
    A21.diagonal.input = name.vector(fit$A21_input_diagonal),
    K1 = name.variable.matrix(fit$K1),
    K2 = name.variable.matrix(fit$K2),
    K3 = name.variable.matrix(fit$K3)
  )
  bootstrap <- list(
    statistics = bootstrap.values,
    multipliers = multiplier.matrix,
    B = B,
    seed = as.numeric(fit$seed),
    exceedances = as.numeric(fit$exceedances),
    p.value = as.numeric(fit$p_value),
    p.value.rule = "(1 + count(T.star >= T.PDQ)) / (B + 1)",
    tie.rule = ">=",
    critical.value = as.numeric(fit$critical_value),
    critical.order.index = as.numeric(fit$critical_rank),
    critical.rule =
      "inf{t: empirical F.star(t) >= 1 - level}",
    primary.rejection.rule = "T.PDQ > critical.value",
    primary.rejected = isTRUE(fit$primary_rejected),
    auxiliary.p.rejected = isTRUE(fit$p_value_rejected),
    decision.disagreement = isTRUE(fit$decision_disagreement),
    conditional.mean.theoretical = 0,
    empirical.mean = as.numeric(fit$bootstrap_mean),
    empirical.variance = as.numeric(fit$bootstrap_empirical_variance),
    variance.diagonal.deleted = as.numeric(
      fit$bootstrap_variance_formula
    ),
    minimum = as.numeric(fit$bootstrap_minimum),
    maximum = as.numeric(fit$bootstrap_maximum)
  )

  p.value <- as.numeric(fit$p_value)
  estimate <- medians$group1$location - medians$group2$location
  names(estimate) <- variable.names
  .new_hd_location_test(
    statistic = c(T.PDQ = as.numeric(fit$statistic)),
    p.value = p.value,
    method = paste(
      "Feng-Wang pairwise-difference-quantile spatial-sign test",
      "(Rademacher wild bootstrap)"
    ),
    data.name = paste(x.name, "and", y.name),
    alternative = "two.sided",
    raw.statistic = c(
      R.PDQ = as.numeric(fit$R_hat),
      b.hat = as.numeric(fit$b_hat),
      T.PDQ = as.numeric(fit$statistic)
    ),
    estimate = estimate,
    null.value = stats::setNames(numeric(p), variable.names),
    null.distribution = list(
      family = "conditional Rademacher wild-bootstrap quadratic form",
      tail = "upper",
      finite.sample.exact = FALSE,
      nuisance.refitted = FALSE,
      primary.rule = "T.PDQ > empirical (1-level) lower quantile",
      auxiliary.p.value.rule =
        "(1 + count(T.star >= T.PDQ)) / (B + 1)",
      assumptions = paste(
        "Independent elliptical samples; growing balanced group sizes;",
        "positive regular PDQ population quantiles and inverse radii;",
        "paper spectral and row-leverage conditions. Arbitrary within-group",
        "correlation and unequal scatter are allowed."
      )
    ),
    variance = c(
      wild.bootstrap.diagonal.deleted =
        as.numeric(fit$bootstrap_variance_formula)
    ),
    components = list(
      scales = scales,
      medians = medians,
      fitted.signs = list(group1 = fitted.signs1, group2 = fitted.signs2),
      fitted.radii = list(
        group1 = name.group.vector(fit$fitted_radii1, 1L),
        group2 = name.group.vector(fit$fitted_radii2, 2L)
      ),
      fitted.inverse.radii = list(
        group1 = name.group.vector(fit$fitted_inverse_radii1, 1L),
        group2 = name.group.vector(fit$fitted_inverse_radii2, 2L)
      ),
      cross.signs = list(group1 = cross.signs1, group2 = cross.signs2),
      cross.inner.products = name.cross.matrix(
        fit$cross_inner_products
      ),
      matrices = matrices,
      observed = observed,
      diagonal.deletion = list(
        within.kernel1 = name.within.matrix(fit$within_kernel1, 1L),
        within.kernel2 = name.within.matrix(fit$within_kernel2, 2L),
        cross.kernel = name.cross.matrix(fit$cross_kernel),
        contribution1 = name.group.vector(
          fit$bias_contributions1, 1L
        ),
        contribution2 = name.group.vector(
          fit$bias_contributions2, 2L
        ),
        bias1 = as.numeric(fit$bias1),
        bias2 = as.numeric(fit$bias2),
        b.hat = as.numeric(fit$b_hat),
        within.denominators = c(group1 = n1^2, group2 = n2^2)
      ),
      bootstrap = bootstrap,
      preconditioning = list(
        x = name.group.matrix(fit$preconditioned_x, 1L),
        y = name.group.matrix(fit$preconditioned_y, 2L),
        magnitude = name.vector(fit$precondition_magnitude),
        midpoint.scaled = name.vector(fit$precondition_midpoint_scaled),
        range.scaled = name.vector(fit$precondition_range_scaled),
        log.unit = name.vector(fit$precondition_log_unit)
      )
    ),
    diagnostics = list(
      quantile.probability = quantile_prob,
      quantile.definition =
        "unordered-pair empirical inverse CDF; no interpolation",
      quantile.order.index = c(
        group1 = as.numeric(fit$quantile_rank1),
        group2 = as.numeric(fit$quantile_rank2)
      ),
      quantile.pair.count = c(
        group1 = as.numeric(fit$quantile_pair_count1),
        group2 = as.numeric(fit$quantile_pair_count2)
      ),
      level = level,
      htest.p.value.is.auxiliary.mc = TRUE,
      htest.p.value.rejected = isTRUE(fit$p_value_rejected),
      primary.paper.rejected = isTRUE(fit$primary_rejected),
      primary.critical.value = as.numeric(fit$critical_value),
      finite.grid.decision.disagreement =
        isTRUE(fit$decision_disagreement),
      primary.tie.rule = "> critical value",
      p.value.tie.rule = ">= observed statistic",
      plus.one.correction = TRUE,
      paper.finite.B.p.value.specified = FALSE,
      B = B,
      minimum.attainable.p = 1 / (B + 1),
      mc.standard.error = sqrt(p.value * (1 - p.value) / (B + 1)),
      seed.requested = seed.requested,
      seed.used = as.numeric(fit$seed),
      counter.generator = "module-local deterministic 32-bit counter hash",
      keep.bootstrap = keep_bootstrap,
      nuisance.refitted.in.bootstrap = FALSE,
      finite.counter.stream.group.swap = paste(
        "the conditional bootstrap law is group-swap invariant, but a",
        "fixed finite counter stream need not give the same Monte Carlo",
        "p-value after swapping labels"
      ),
      median.strict = strict,
      median.tolerance = controls$tol,
      median.max.iterations = controls$max_iter,
      median.converged = converged,
      median.iterations = c(
        group1 = as.integer(fit$median_iterations1),
        group2 = as.integer(fit$median_iterations2)
      ),
      median.equation.residual = c(
        group1 = as.numeric(fit$median_equation_residual1),
        group2 = as.numeric(fit$median_equation_residual2)
      ),
      fitted.zero.residual.count = c(group1 = 0, group2 = 0),
      cross.zero.sign.count = c(
        group1 = as.numeric(fit$cross_zero_count1),
        group2 = as.numeric(fit$cross_zero_count2)
      ),
      U.zero.convention = "U(0) = 0 only where no inverse radius is needed",
      G.eigenvalues = list(
        group1 = as.numeric(fit$G_eigenvalues1),
        group2 = as.numeric(fit$G_eigenvalues2)
      ),
      G.condition.number = c(
        group1 = as.numeric(fit$G_condition1),
        group2 = as.numeric(fit$G_condition2)
      ),
      nonrepresentable.input.quantiles = c(
        group1 = as.numeric(fit$nonrepresentable_quantiles1),
        group2 = as.numeric(fit$nonrepresentable_quantiles2)
      ),
      nonrepresentable.input.diagonals = c(
        group1 = as.numeric(fit$nonrepresentable_diagonals1),
        group2 = as.numeric(fit$nonrepresentable_diagonals2)
      ),
      internal.preconditioning =
        "common coordinatewise midrange/range affine transformation",
      scale.identification =
        "each working PDQ square-root diagonal normalized by its maximum",
      no.repair.policy = paste(
        "no ridge, scale floor, generalized inverse, absolute-value",
        "variance correction, deletion of fitted zero radii, or random",
        "perturbation"
      ),
      asymptotic.calibration = TRUE,
      paper.simulation.reproduced = FALSE
    ),
    n = c(x = n1, y = n2),
    p = p,
    call = call
  )
}
