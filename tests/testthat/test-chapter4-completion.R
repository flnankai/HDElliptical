ch4c_alpha_fixture <- function() {
  time <- seq_len(12)
  factors <- matrix(seq(-0.6, 1.4, length.out = 12), ncol = 1L)
  colnames(factors) <- "market"
  error <- cbind(
    sin(0.7 * time) + 0.2 * cos(1.3 * time),
    1.2 * cos(0.4 * time) + 0.3 * sin(1.1 * time),
    0.8 * sin(0.3 * time + 1) + 0.6 * cos(0.9 * time)
  )
  returns <- sweep(error, 2L, c(0.15, -0.08, 0.05), "+") +
    factors %*% matrix(c(0.6, -0.4, 0.2), nrow = 1L)
  colnames(returns) <- c("A", "B", "C")
  list(returns = returns, factors = factors)
}


ch4c_local_permutations <- function(n, B, seed) {
  had.state <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had.state) old.state <- get(".Random.seed", envir = .GlobalEnv)
  on.exit({
    if (had.state) {
      assign(".Random.seed", old.state, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv,
                      inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  set.seed(seed)
  replicate(B, sample.int(n), simplify = "matrix")
}


ch4c_orderings <- function(values) {
  if (length(values) == 1L) return(matrix(values, nrow = 1L))
  do.call(rbind, lapply(seq_along(values), function(index) {
    cbind(values[index], ch4c_orderings(values[-index]))
  }))
}


ch4c_pair_contrast <- function(values, ordering, anchor) {
  ((values[ordering[1L]] <= values[ordering[anchor]]) -
     (values[ordering[2L]] <= values[ordering[anchor]])) *
    ((values[ordering[3L]] <= values[ordering[anchor]]) -
       (values[ordering[4L]] <= values[ordering[anchor]]))
}


ch4c_separated <- function(values, a, b, c, d) {
  as.integer(
    values[a] < values[c] && values[a] < values[d] &&
      values[b] < values[c] && values[b] < values[d]
  )
}


ch4c_tau_contrast <- function(values, ordering) {
  ch4c_separated(values, ordering[1L], ordering[3L],
                 ordering[2L], ordering[4L]) +
    ch4c_separated(values, ordering[2L], ordering[4L],
                   ordering[1L], ordering[3L]) -
    ch4c_separated(values, ordering[1L], ordering[4L],
                   ordering[2L], ordering[3L]) -
    ch4c_separated(values, ordering[2L], ordering[3L],
                   ordering[1L], ordering[4L])
}


ch4c_kernel <- function(x, y, measure, orderings) {
  terms <- apply(orderings, 1L, function(ordering) {
    if (measure == "hoeffding_d") {
      ch4c_pair_contrast(x, ordering, 5L) *
        ch4c_pair_contrast(y, ordering, 5L)
    } else if (measure == "bkr_r") {
      ch4c_pair_contrast(x, ordering, 5L) *
        ch4c_pair_contrast(y, ordering, 6L)
    } else {
      ch4c_tau_contrast(x, ordering) *
        ch4c_tau_contrast(y, ordering)
    }
  })
  sum(terms) / switch(
    measure,
    hoeffding_d = 16,
    bkr_r = 32,
    tau_star = factorial(4)
  )
}


ch4c_u_matrix <- function(x, y, measure) {
  order <- c(hoeffding_d = 5L, bkr_r = 6L, tau_star = 4L)[[measure]]
  subsets <- combn(nrow(x), order, simplify = FALSE)
  orderings <- ch4c_orderings(seq_len(order))
  answer <- matrix(0, ncol(x), ncol(y))
  for (left in seq_len(ncol(x))) {
    for (right in seq_len(ncol(y))) {
      answer[left, right] <- mean(vapply(subsets, function(index) {
        ch4c_kernel(x[index, left], y[index, right], measure, orderings)
      }, numeric(1)))
    }
  }
  answer
}


ch4c_vector_reference <- function(x, y, measure, permutations) {
  n <- nrow(x)
  estimates <- ch4c_u_matrix(x, y, measure)
  null.second <- switch(
    measure,
    hoeffding_d = 2 * (n^2 + 5 * n - 32) /
      (9 * n * (n - 1) * (n - 3) * (n - 4)),
    bkr_r = 2 * (n^3 - 3 * n^2 - 6 * n + 10) /
      (n * (n - 1) * (n - 2) * (n - 3) * (n - 4)),
    tau_star = 8 * (3 * n^2 + 5 * n - 18) /
      (75 * n * (n - 1) * (n - 2) * (n - 3))
  )
  statistic.from <- function(current) {
    sum(current^2) - length(current) * null.second
  }
  permutation.statistics <- apply(permutations, 2L, function(index) {
    statistic.from(ch4c_u_matrix(x[index, , drop = FALSE], y, measure))
  })
  sum.statistic <- statistic.from(estimates)
  permutation.variance <- stats::var(permutation.statistics)
  sum.z <- sum.statistic / sqrt(permutation.variance)
  p.sum <- stats::pnorm(sum.z, lower.tail = FALSE)
  maximum <- max(estimates)
  comparisons <- length(estimates)
  denominator <- c(hoeffding_d = 30, bkr_r = 90,
                   tau_star = 36)[[measure]]
  max.gumbel <- pi^4 * (n - 1) * maximum / denominator -
    2 * log(comparisons) + log(log(comparisons)) + pi^4 / 36
  p.max <- -expm1(
    -2.467 / sqrt(pi) * exp(-max.gumbel / 2)
  )
  fisher <- -2 * (log(p.sum) + log(p.max))
  list(
    estimates = estimates,
    null.second = null.second,
    maximum = maximum,
    sum.statistic = sum.statistic,
    permutation.statistics = permutation.statistics,
    permutation.variance = permutation.variance,
    sum.z = sum.z,
    p.sum = p.sum,
    max.gumbel = max.gumbel,
    p.max = p.max,
    fisher = fisher,
    p.fisher = stats::pchisq(fisher, 4, lower.tail = FALSE)
  )
}


ch4c_vector_fixture <- function() {
  list(
    x = cbind(
      c(1, 4, 2, 6, 3, 5),
      c(6, 2, 5, 1, 4, 3)
    ),
    y = matrix(c(2, 6, 4, 1, 5, 3), ncol = 1L)
  )
}


test_that("weighted spatial-sign oracle is the literal supplied-score formula", {
  angles <- c(0, pi / 2, pi, 3 * pi / 2, pi / 4)
  directions <- cbind(cos(angles), sin(angles))
  radii <- c(1, 1.5, 2, 2.5, 3)
  h <- c(1, -0.4, 0.7, 1.2, -0.3)
  trace.R2 <- 2.75
  weights <- 1 / radii
  off.diagonal <- 0
  for (first in seq_len(nrow(directions))) {
    for (second in seq_len(nrow(directions))) {
      if (first == second) next
      off.diagonal <- off.diagonal +
        h[first] * h[second] * weights[first] * weights[second] *
        sum(directions[first, ] * directions[second, ])
    }
  }
  Q <- ncol(directions) / sum(h^2) * off.diagonal
  psi2 <- mean(weights^2)
  expected <- Q / sqrt(2 * psi2^2 * trace.R2)
  observed <- weighted_spatial_sign_alpha_oracle_test(
    directions, radii, h, trace.R2, function(value) 1 / value
  )

  expect_s3_class(observed, "hd_alpha_test")
  expect_equal(unname(observed$raw.statistic), Q, tolerance = 2e-14)
  expect_equal(observed$components$psi2, psi2, tolerance = 2e-15)
  expect_equal(unname(observed$statistic), expected, tolerance = 2e-14)
  expect_equal(observed$p.value, pnorm(expected, lower.tail = FALSE),
               tolerance = 2e-15)
  expect_false(observed$diagnostics$primary.feasible.general.K)
  expect_identical(observed$diagnostics$weight.cap, "none")

  non.unit <- directions
  non.unit[1, ] <- 2 * non.unit[1, ]
  expect_error(
    weighted_spatial_sign_alpha_oracle_test(
      non.unit, radii, h, trace.R2, function(value) 1 / value
    ),
    "unit Euclidean norm"
  )
  expect_error(
    weighted_spatial_sign_alpha_oracle_test(
      directions, as.character(radii), h, trace.R2,
      function(value) 1 / value
    ),
    "numeric vector"
  )
  expect_error(
    weighted_spatial_sign_alpha_oracle_test(
      directions, replace(radii, 1L, 0), h, trace.R2,
      function(value) 1 / value
    ),
    "strictly positive"
  )
  expect_error(
    weighted_spatial_sign_alpha_oracle_test(
      directions, radii, h, -trace.R2, function(value) 1 / value
    ),
    "strictly positive"
  )
})


test_that("INST uses unrestricted radii, restricted signs, and literal trace", {
  fixture <- ch4c_alpha_fixture()
  observed <- zhao_chen_zi_inst_alpha_test(
    fixture$returns, fixture$factors, tol = 1e-9,
    max_iter = 3000L, keep_scores = TRUE
  )
  y <- sweep(fixture$returns, 2L,
             apply(abs(fixture$returns), 2L, max), "/")
  f <- sweep(fixture$factors, 2L,
             apply(abs(fixture$factors), 2L, max), "/")
  unrestricted <- lm.fit(cbind(1, f), y)$residuals
  diagonal <- observed$components$scale.diagonal.scaled.coordinates
  standardized <- sweep(unrestricted, 2L, sqrt(diagonal), "/")
  radii <- sqrt(rowSums(standardized^2))
  weights <- 1 / radii
  U <- observed$components$directions
  h <- observed$components$h
  off.diagonal <- 0
  for (first in seq_len(nrow(U))) {
    for (second in seq_len(nrow(U))) {
      if (first == second) next
      off.diagonal <- off.diagonal +
        h[first] * h[second] * weights[first] * weights[second] *
        sum(U[first, ] * U[second, ])
    }
  }
  Q <- ncol(U) / sum(h^2) * off.diagonal
  psi2 <- mean(weights^2)
  denominator <- sqrt(
    2 * psi2^2 * observed$components$trace.R.squared
  )

  expect_equal(observed$components$unrestricted.radii, radii,
               tolerance = 2e-12)
  expect_equal(observed$components$inverse.radial.weights, weights,
               tolerance = 2e-12)
  expect_equal(observed$components$Q.K.alpha, Q, tolerance = 2e-10)
  expect_equal(observed$components$psi2.K.hat, psi2,
               tolerance = 2e-12)
  expect_equal(observed$components$studentizing.denominator, denominator,
               tolerance = 2e-10)
  expect_equal(unname(observed$statistic), Q / denominator,
               tolerance = 2e-10)
  expect_identical(observed$diagnostics$direction.residuals,
                   "restricted factor residuals Y - B f")
  expect_identical(observed$diagnostics$radial.residuals,
                   "unrestricted OLS residuals Y - alpha - B f")
  expect_identical(observed$diagnostics$weight.cap, "none")

  expect_error(
    zhao_chen_zi_inst_alpha_test(
      fixture$returns, fixture$factors, tol = 1e-9,
      max_iter = 3000L, zero_tol = max(radii) + 1
    ),
    "no floor or ridge|inverse-norm weighting is undefined"
  )
})


test_that("book Gaussian Cauchy benchmark matches every displayed formula", {
  fixture <- ch4c_alpha_fixture()
  trace.R2 <- 4.25
  design <- cbind(1, fixture$factors)
  fit <- lm.fit(design, fixture$returns)
  residual.variance <- colSums(fit$residuals^2) /
    (nrow(design) - ncol(design))
  h <- rep(1, nrow(design)) - fixture$factors %*%
    solve(crossprod(fixture$factors),
          crossprod(fixture$factors, rep(1, nrow(design))))
  t.squared <- fit$coefficients[1, ]^2 * sum(h^2) / residual.variance
  quadratic <- nrow(design) *
    sum(fit$coefficients[1, ]^2 / residual.variance)
  sum.statistic <- (quadratic - ncol(fixture$returns)) /
    sqrt(2 * trace.R2)
  maximum <- max(t.squared)
  centered.maximum <- maximum - 2 * log(ncol(fixture$returns)) +
    log(log(ncol(fixture$returns)))
  p.sum <- pnorm(sum.statistic, lower.tail = FALSE)
  p.max <- -expm1(-exp(-centered.maximum / 2) / sqrt(pi))
  cauchy <- 0.5 * tan(pi * (0.5 - p.sum)) +
    0.5 * tan(pi * (0.5 - p.max))
  observed <- book_gaussian_alpha_cauchy_test(
    fixture$returns, trace_R2 = trace.R2, factors = fixture$factors
  )

  expect_equal(observed$components$residual.variance,
               residual.variance, tolerance = 2e-12)
  expect_equal(observed$components$quadratic, quadratic,
               tolerance = 2e-12)
  expect_equal(observed$components$sum.statistic, sum.statistic,
               tolerance = 2e-12)
  expect_equal(observed$components$maximum.t.squared, maximum,
               tolerance = 2e-12)
  expect_equal(observed$components$centered.maximum, centered.maximum,
               tolerance = 2e-12)
  expect_equal(observed$components$p.values,
               c(sum = p.sum, max = p.max), tolerance = 2e-14)
  expect_equal(unname(observed$statistic), cauchy, tolerance = 1e-10)
  expect_equal(observed$p.value, pcauchy(cauchy, lower.tail = FALSE),
               tolerance = 1e-11)
  expect_identical(observed$diagnostics$trace.R2.required, TRUE)
  expect_match(observed$diagnostics$trace.estimator, "supplied")

  expect_error(
    book_gaussian_alpha_cauchy_test(
      fixture$returns, factors = fixture$factors
    ),
    "trace_R2"
  )
  expect_error(
    book_gaussian_alpha_cauchy_test(
      fixture$returns, trace_R2 = -1, factors = fixture$factors
    ),
    "strictly positive supplied"
  )
})


test_that("conditional Wald is a strict supplied-covariance benchmark", {
  delta <- c(first = 0.2, second = -0.1)
  covariance <- matrix(c(2, 0.4, 0.4, 1.5), 2L, 2L)
  sample.size <- 40L
  expected <- sample.size * drop(
    delta %*% solve(covariance, delta)
  )
  observed <- conditional_factor_wald_test(
    delta, covariance, sample.size
  )

  expect_s3_class(observed, "htest")
  expect_equal(unname(observed$statistic), expected, tolerance = 2e-14)
  expect_equal(observed$p.value,
               pchisq(expected, 2, lower.tail = FALSE),
               tolerance = 2e-15)
  expect_identical(names(observed$estimate), names(delta))
  expect_identical(observed$diagnostics$linear.solve,
                   "strict Cholesky")

  asymmetric <- covariance
  asymmetric[1, 2] <- asymmetric[1, 2] + 1e-12
  expect_error(
    conditional_factor_wald_test(delta, asymmetric, sample.size),
    "exactly symmetric"
  )
  expect_error(
    conditional_factor_wald_test(delta, matrix(1, 2L, 2L), sample.size),
    "strictly positive definite"
  )
  expect_error(
    conditional_factor_wald_test(as.character(delta), covariance,
                                 sample.size),
    "numeric vector"
  )
  expect_error(
    conditional_factor_wald_test(delta, covariance, "40"),
    "must be numeric"
  )
})


test_that("all three high-order vector U methods match pure R kernels", {
  fixture <- ch4c_vector_fixture()
  B <- 4L
  seed <- 27L
  permutations <- ch4c_local_permutations(
    nrow(fixture$x), B, seed
  )
  for (measure in c("hoeffding_d", "bkr_r", "tau_star")) {
    reference <- ch4c_vector_reference(
      fixture$x, fixture$y, measure, permutations
    )
    set.seed(710)
    state <- .Random.seed
    observed <- wang_liu_feng_vector_u_independence_test(
      fixture$x, fixture$y, measure = measure, component = "fisher",
      B = B, seed = seed, keep_estimates = TRUE,
      keep_permutation = TRUE
    )

    expect_identical(.Random.seed, state)
    expect_equal(unname(observed$components$permutation.indices),
                 unname(permutations), tolerance = 0)
    expect_equal(unname(observed$components$estimates),
                 reference$estimates, tolerance = 2e-14)
    expect_equal(observed$components$null_second_moment,
                 reference$null.second, tolerance = 2e-15)
    expect_equal(observed$components$maximum,
                 reference$maximum, tolerance = 2e-14)
    expect_equal(observed$components$sum_statistic,
                 reference$sum.statistic, tolerance = 2e-14)
    expect_equal(observed$components$permutation_statistics,
                 reference$permutation.statistics, tolerance = 2e-14)
    expect_equal(observed$components$permutation_variance,
                 reference$permutation.variance, tolerance = 2e-14)
    expect_equal(observed$components$sum.z,
                 reference$sum.z, tolerance = 2e-13)
    expect_equal(observed$components$p.sum,
                 reference$p.sum, tolerance = 2e-14)
    expect_equal(observed$components$max.gumbel,
                 reference$max.gumbel, tolerance = 2e-13)
    expect_equal(observed$components$p.max,
                 reference$p.max, tolerance = 2e-14)
    expect_equal(unname(observed$statistic), reference$fisher,
                 tolerance = 2e-13)
    expect_equal(observed$p.value, reference$p.fisher,
                 tolerance = 2e-14)
    order <- c(hoeffding_d = 5L, bkr_r = 6L,
               tau_star = 4L)[[measure]]
    expected.workload <- (B + 1) * ncol(fixture$x) *
      ncol(fixture$y) * choose(nrow(fixture$x), order) * factorial(order)
    expect_equal(observed$components$total_kernel_evaluations,
                 expected.workload, tolerance = 0)
    expect_identical(observed$diagnostics$scalable.claim, FALSE)
    expect_identical(observed$diagnostics$approximation, "none")
    expect_identical(observed$diagnostics$simulation.adjustment,
                     "not implemented")
  }
})


test_that("vector U guard, ties, RNG, and zero permutation variance are strict", {
  fixture <- ch4c_vector_fixture()
  set.seed(919)
  state <- .Random.seed
  expect_error(
    wang_liu_feng_vector_u_independence_test(
      fixture$x, fixture$y, measure = "tau_star", B = 2L,
      max_kernel_evaluations = 1L
    ),
    "workload exceeds"
  )
  expect_identical(.Random.seed, state)

  tied <- fixture$x
  tied[2L, 1L] <- tied[1L, 1L]
  expect_error(
    wang_liu_feng_vector_u_independence_test(
      tied, fixture$y, measure = "tau_star", B = 2L
    ),
    "exact tie"
  )
  expect_identical(.Random.seed, state)

  small.x <- fixture$x[1:4, , drop = FALSE]
  small.y <- fixture$y[1:4, , drop = FALSE]
  duplicate.seed <- NA_integer_
  for (candidate in 0:10000) {
    draws <- ch4c_local_permutations(4L, 2L, candidate)
    if (identical(draws[, 1L], draws[, 2L])) {
      duplicate.seed <- candidate
      break
    }
  }
  expect_false(is.na(duplicate.seed))
  expect_error(
    wang_liu_feng_vector_u_independence_test(
      small.x, small.y, measure = "tau_star", B = 2L,
      seed = duplicate.seed
    ),
    "variance estimate is not strictly positive"
  )
  expect_identical(.Random.seed, state)
})




