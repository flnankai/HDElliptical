.wx_positive_integer <- function(value, name,
                                 upper = .Machine$integer.max) {
  if (!is.numeric(value) || length(value) != 1L || is.na(value) ||
      !is.finite(value) || value < 1 || value != floor(value) ||
      value > upper) {
    stop(sprintf("`%s` must be a positive integer no larger than %s.",
                 name, format(upper, scientific = FALSE)), call. = FALSE)
  }
  as.integer(value)
}

.wx_alpha <- function(alpha) {
  alpha <- as.numeric(alpha)
  if (length(alpha) != 1L || is.na(alpha) || !is.finite(alpha) ||
      alpha <= 0 || alpha >= 1) {
    stop("`alpha` must be one finite number strictly between zero and one.",
         call. = FALSE)
  }
  alpha
}

.wx_logical_scalar <- function(value, name) {
  if (!is.logical(value) || length(value) != 1L || is.na(value)) {
    stop(sprintf("`%s` must be `TRUE` or `FALSE`.", name), call. = FALSE)
  }
  value
}

.wx_seed <- function(seed) {
  if (is.null(seed)) return(NULL)
  seed <- as.numeric(seed)
  if (length(seed) != 1L || is.na(seed) || !is.finite(seed) ||
      seed < 0 || seed > 2^32 - 1 || seed != floor(seed)) {
    stop("`seed` must be `NULL` or one integer-valued number in " %+%
           "[0, 2^32 - 1].", call. = FALSE)
  }
  seed
}

.wx_resolve_calibration <- function(n1, n2, calibration, max_exact) {
  m1 <- n1 %/% 2L
  m2 <- n2 %/% 2L
  sign.dimension <- m1 + m2
  reduced.exponent <- sign.dimension - 1L
  can.enumerate <- reduced.exponent <= 30L &&
    2^reduced.exponent <= max_exact
  used <- switch(
    calibration,
    auto = if (can.enumerate) "exact" else "monte_carlo",
    exact = "exact",
    monte_carlo = "monte_carlo"
  )
  if (identical(used, "exact") && !can.enumerate) {
    stop(
      sprintf(
        paste0(
          "Exact Wang--Xu calibration needs 2^(%d) global-sign-reduced ",
          "patterns, exceeding `max_exact = %s`; use ",
          "`calibration = \"monte_carlo\"` or increase `max_exact`."
        ),
        reduced.exponent, format(max_exact, scientific = FALSE)
      ),
      call. = FALSE
    )
  }
  list(
    used = used,
    sign.dimension = sign.dimension,
    reduced.exponent = reduced.exponent,
    can.enumerate = can.enumerate
  )
}

.wx_pair_labels <- function(indices, row.labels) {
  indices <- as.matrix(indices)
  paste(row.labels[indices[, 2L]], row.labels[indices[, 1L]], sep = "-")
}

#' Wang--Xu approximate randomization test
#'
#' Tests equality of two high-dimensional mean vectors under unrestricted
#' and potentially unequal covariance matrices using the approximate
#' randomization calibration of Wang and Xu (2022). Observations are rows and
#' variables are columns.
#'
#' The observed statistic is the full-sample Chen--Qin statistic
#' \deqn{
#' T_{CQ}=\sum_{k=1}^2\frac{2}{n_k(n_k-1)}
#' \sum_{i<j}X_{k,i}^{\mathsf T}X_{k,j}
#' -\frac{2}{n_1n_2}\sum_i\sum_jX_{1,i}^{\mathsf T}X_{2,j}.
#' }
#' The reference sample is not produced by permuting pooled group labels.
#' Within each group the method forms adjacent half-differences
#' \deqn{
#' \widetilde X_{k,i}=(X_{k,2i}-X_{k,2i-1})/2,
#' \qquad i=1,\ldots,m_k,\quad m_k=\lfloor n_k/2\rfloor,
#' }
#' and multiplies each half-difference by an independent Rademacher sign.
#' If a group size is odd, its final row is deliberately unused by the
#' reference distribution, exactly as in the paper. The observed statistic
#' still uses every row.
#'
#' With `calibration = "exact"`, all sign configurations are integrated
#' exactly. Since simultaneous reversal of every sign leaves the statistic
#' unchanged, the implementation evaluates one representative of each
#' global-sign pair. The exact conditional tail is
#' \deqn{2^{-(m_1+m_2)}\sum_e 1\{T_{CQ}(e)\geq T_{CQ}\},}
#' with no plus-one correction. With `calibration = "monte_carlo"`, the
#' paper's practical p-value is used literally:
#' \deqn{
#' \widehat p=\frac{1+\sum_{b=1}^B
#' 1\{T_{CQ}^{(b)}\geq T_{CQ}\}}{B+1}.
#' }
#' Thus ties are always counted in the upper tail. `calibration = "auto"`
#' selects exact enumeration only when its reduced number of configurations
#' does not exceed `max_exact`.
#'
#' Exact enumeration removes Monte Carlo error from the conditional
#' reference distribution; it does not make the overall Behrens--Fisher test
#' finite-sample exact. Its level guarantee is asymptotic under the paper's
#' moment and non-dominating-observation assumptions. Pairing follows the
#' input row order, so arbitrary row reordering can change the finite-sample
#' reference distribution.
#'
#' @param x,y Numeric matrices or data frames containing the two independent
#'   samples. Both must have the same variables and at least four rows.
#' @param alpha Test level strictly between zero and one.
#' @param calibration Reference calculation: `"auto"`, `"exact"`, or
#'   `"monte_carlo"`.
#' @param B Positive number of Monte Carlo sign draws. It is validated but
#'   otherwise ignored under exact calibration.
#' @param seed `NULL`, or an integer-valued counter-generator seed in
#'   `[0, 2^32 - 1]`. Under Monte Carlo calibration, `NULL` draws and records
#'   one seed from R's RNG; an explicit seed makes the result reproducible
#'   without changing R's RNG state. Exact enumeration does not use a seed.
#' @param workers Number of workers. The first implementation deliberately
#'   accepts only `1`; it never claims or silently performs parallel work.
#' @param max_exact Positive maximum number of global-sign-reduced patterns
#'   permitted for exact enumeration.
#' @param keep_randomized If `TRUE`, retain every randomized statistic and
#'   sign pattern. This can require substantial memory.
#'
#' @return An object of class `c("hd_location_test", "htest")`. Raw
#'   components include the observed Chen--Qin decomposition, half-difference
#'   pseudo-samples, pairing maps, scaled quadratic kernel, tail count, and
#'   optional randomized values and signs. Diagnostics distinguish exhaustive
#'   reference calculation from finite-sample exactness and record the
#'   plus-one rule, seed, Monte Carlo resolution, discarded rows, numerical
#'   scaling, and degeneracy.
#'
#' @references
#' Wang, R. and Xu, W. (2022). An approximate randomization test for the
#' high-dimensional two-sample Behrens--Fisher problem under arbitrary
#' covariances. *Biometrika*, 109, 1117--1132.
#' \doi{10.1093/biomet/asac014}.
#'
#' @examples
#' x <- matrix(c(1, 0, 2, 1, 4, 1, 5, 3), 4, 2, byrow = TRUE)
#' y <- matrix(c(0, 2, 1, 4, 3, 5, 4, 7), 4, 2, byrow = TRUE)
#' wang_xu_approx_randomization_test(x, y, calibration = "exact")
#'
#' @export
wang_xu_approx_randomization_test <- function(
    x, y, alpha = 0.05,
    calibration = c("auto", "exact", "monte_carlo"),
    B = 9999L, seed = NULL, workers = 1L,
    max_exact = 1048576L, keep_randomized = FALSE) {
  call <- match.call()
  x.name <- deparse1(substitute(x))
  y.name <- deparse1(substitute(y))
  x <- .as_data_matrix(x, "x", min_rows = 4L)
  y <- .as_data_matrix(y, "y", min_rows = 4L)
  .check_two_sample_variables(x, y)

  alpha <- .wx_alpha(alpha)
  calibration.requested <- match.arg(calibration)
  B <- .wx_positive_integer(B, "B")
  workers <- .wx_positive_integer(workers, "workers")
  if (workers != 1L) {
    stop(
      "`workers > 1` is not implemented in this release; the function " %+%
        "does not silently claim parallel execution.",
      call. = FALSE
    )
  }
  max_exact <- .wx_positive_integer(max_exact, "max_exact")
  keep_randomized <- .wx_logical_scalar(keep_randomized, "keep_randomized")
  seed.requested <- .wx_seed(seed)

  n1 <- nrow(x)
  n2 <- nrow(y)
  p <- ncol(x)
  resolution <- .wx_resolve_calibration(
    n1, n2, calibration.requested, max_exact
  )
  calibration.used <- resolution$used
  if (identical(calibration.used, "monte_carlo")) {
    seed.used <- if (is.null(seed.requested)) {
      as.numeric(sample.int(.Machine$integer.max, 1L) - 1L)
    } else {
      seed.requested
    }
  } else {
    seed.used <- 0
  }

  fit <- cpp_wang_xu_approx_randomization(
    x, y, calibration.used, B, seed.used, max_exact, keep_randomized
  )
  p.value <- as.numeric(fit$p_value)
  variable.names <- .hotelling_variable_names(x)
  row.names.x <- rownames(x)
  row.names.y <- rownames(y)
  if (is.null(row.names.x)) row.names.x <- paste0("x", seq_len(n1))
  if (is.null(row.names.y)) row.names.y <- paste0("y", seq_len(n2))

  named.vector <- function(value) {
    stats::setNames(as.numeric(value), variable.names)
  }
  name.pseudo <- function(value, labels) {
    value <- as.matrix(value)
    dimnames(value) <- list(labels, variable.names)
    value
  }
  pair.labels.x <- .wx_pair_labels(fit$pair_indices_x, row.names.x)
  pair.labels.y <- .wx_pair_labels(fit$pair_indices_y, row.names.y)
  pseudo.x <- name.pseudo(fit$pseudo_x, pair.labels.x)
  pseudo.y <- name.pseudo(fit$pseudo_y, pair.labels.y)
  pseudo.x.scaled <- name.pseudo(fit$pseudo_x_scaled, pair.labels.x)
  pseudo.y.scaled <- name.pseudo(fit$pseudo_y_scaled, pair.labels.y)
  pair.indices.x <- as.matrix(fit$pair_indices_x)
  pair.indices.y <- as.matrix(fit$pair_indices_y)
  rownames(pair.indices.x) <- pair.labels.x
  rownames(pair.indices.y) <- pair.labels.y
  colnames(pair.indices.x) <- colnames(pair.indices.y) <- c("first", "second")

  kernel <- as.matrix(fit$quadratic_kernel_scaled)
  sign.labels <- c(paste0("epsilon.x", seq_len(nrow(pseudo.x))),
                   paste0("epsilon.y", seq_len(nrow(pseudo.y))))
  dimnames(kernel) <- list(sign.labels, sign.labels)
  sign.patterns <- fit$sign_patterns
  if (!is.null(sign.patterns)) {
    sign.patterns <- as.matrix(sign.patterns)
    colnames(sign.patterns) <- sign.labels
  }
  discarded.x <- as.integer(fit$discarded_rows_x)
  discarded.y <- as.integer(fit$discarded_rows_y)
  if (length(discarded.x)) names(discarded.x) <- row.names.x[discarded.x]
  if (length(discarded.y)) names(discarded.y) <- row.names.y[discarded.y]

  plus.one <- isTRUE(fit$plus_one_correction)
  p.rule <- if (plus.one) {
    "(1 + count(T.randomized >= T.observed)) / (B + 1)"
  } else {
    "count(T.randomized >= T.observed) / exhaustive reference size"
  }
  randomization.variance <- as.numeric(fit$randomized_variance)

  .new_hd_location_test(
    statistic = c(T.CQ = as.numeric(fit$observed)),
    p.value = p.value,
    method = paste(
      "Wang-Xu high-dimensional two-sample approximate randomization test",
      sprintf("(%s Rademacher calibration)", calibration.used)
    ),
    data.name = paste(x.name, "and", y.name),
    alternative = "two.sided",
    raw.statistic = c(T.CQ = as.numeric(fit$observed)),
    estimate = named.vector(fit$difference),
    null.value = stats::setNames(numeric(p), variable.names),
    null.distribution = list(
      family = "conditional Rademacher half-difference reference",
      exact = FALSE,
      conditional.reference.exhaustive =
        isTRUE(fit$exact_reference_enumerated),
      tail = "upper",
      tie.rule = ">=",
      p.value.rule = p.rule,
      assumptions = paste(
        "Wang-Xu asymptotics: independent samples, min(n1,n2) increasing,",
        "uniform fourth-moment control, and no dominating observation when",
        "within-group covariances vary; no covariance eigenstructure",
        "restriction"
      )
    ),
    variance = c(randomization = randomization.variance),
    components = list(
      mean.x = named.vector(fit$mean_x),
      mean.y = named.vector(fit$mean_y),
      difference = named.vector(fit$difference),
      observed = list(
        T.CQ = as.numeric(fit$observed),
        mean.difference.squared =
          as.numeric(fit$mean_difference_squared),
        trace.S1 = as.numeric(fit$trace_S1),
        trace.S2 = as.numeric(fit$trace_S2),
        identity = "||mean.x-mean.y||^2 - trace(S1)/n1 - trace(S2)/n2",
        T.CQ.scaled = as.numeric(fit$observed_scaled),
        mean.difference.squared.scaled =
          as.numeric(fit$mean_difference_squared_scaled),
        trace.S1.scaled = as.numeric(fit$trace_S1_scaled),
        trace.S2.scaled = as.numeric(fit$trace_S2_scaled)
      ),
      pseudo.x = pseudo.x,
      pseudo.y = pseudo.y,
      pseudo.x.scaled = pseudo.x.scaled,
      pseudo.y.scaled = pseudo.y.scaled,
      pair.indices.x = pair.indices.x,
      pair.indices.y = pair.indices.y,
      discarded.rows.x = discarded.x,
      discarded.rows.y = discarded.y,
      randomized.quadratic.kernel.scaled = kernel,
      randomized.statistics = fit$randomized_statistics,
      randomized.statistics.scaled = fit$randomized_statistics_scaled,
      sign.patterns = sign.patterns,
      randomized.summary = c(
        mean = as.numeric(fit$randomized_mean),
        variance = randomization.variance,
        minimum = as.numeric(fit$randomized_minimum),
        maximum = as.numeric(fit$randomized_maximum)
      ),
      randomized.summary.scaled = c(
        mean = as.numeric(fit$randomized_mean_scaled),
        variance = as.numeric(fit$randomized_variance_scaled),
        minimum = as.numeric(fit$randomized_minimum_scaled),
        maximum = as.numeric(fit$randomized_maximum_scaled)
      ),
      exceedances = as.numeric(fit$exceedances),
      empirical.tail.without.plus.one = as.numeric(fit$empirical_tail),
      reference.evaluations = as.numeric(fit$reference_evaluations),
      total.sign.configurations = as.numeric(fit$total_sign_configurations),
      log2.total.sign.configurations =
        as.numeric(fit$log2_total_sign_configurations)
    ),
    diagnostics = list(
      alpha = alpha,
      rejected = p.value <= alpha,
      observed.statistic.sample = "all original observations",
      reference.sample =
        "adjacent within-group half-differences in input row order",
      randomization.mechanism = "independent Rademacher sign multipliers",
      pooled.label.permutation = FALSE,
      pairing.order.sensitive = TRUE,
      half.difference.divisor = 2,
      m1 = as.integer(fit$m1),
      m2 = as.integer(fit$m2),
      sign.dimension = as.integer(fit$sign_dimension),
      discarded.rows.x = discarded.x,
      discarded.rows.y = discarded.y,
      calibration.requested = calibration.requested,
      calibration.used = calibration.used,
      conditional.reference.exhaustive =
        isTRUE(fit$exact_reference_enumerated),
      finite.sample.test.exact = FALSE,
      global.sign.symmetry.reduced =
        isTRUE(fit$global_sign_symmetry_reduced),
      reference.evaluations = as.numeric(fit$reference_evaluations),
      total.sign.configurations = as.numeric(fit$total_sign_configurations),
      exceedances = as.numeric(fit$exceedances),
      tie.rule = ">=",
      plus.one.correction = plus.one,
      p.value.rule = p.rule,
      B = if (identical(calibration.used, "monte_carlo")) B else NULL,
      minimum.attainable.p = as.numeric(fit$minimum_attainable_p),
      mc.standard.error = as.numeric(fit$mc_standard_error),
      seed.requested = seed.requested,
      seed.used = if (identical(calibration.used, "monte_carlo"))
        as.numeric(fit$seed) else NULL,
      counter.generator = if (identical(calibration.used, "monte_carlo"))
        "module-local deterministic 32-bit counter hash" else NULL,
      workers = workers,
      parallel.execution = FALSE,
      keep.randomized = keep_randomized,
      randomization.degenerate = isTRUE(fit$randomization_degenerate),
      internal.scaling =
        "one common global scale after common-anchor subtraction",
      internal.scale.factor = as.numeric(fit$global_scale),
      subtraction.overflow.fallback =
        isTRUE(fit$subtraction_overflow_fallback),
      observed.physical.underflow =
        isTRUE(fit$observed_physical_underflow),
      covariance.regularization = "none",
      p.value.repair = "none"
    ),
    n = c(x = n1, y = n2),
    p = p,
    call = call
  )
}
