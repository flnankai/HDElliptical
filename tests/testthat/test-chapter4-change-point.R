.cp_test_fixture <- function(n = 20L, p = 6L, changed = FALSE) {
  set.seed(4107L + n + p)
  x <- matrix(stats::rnorm(n * p), n, p)
  x <- sweep(x, 2L, seq(0.7, 1.6, length.out = p), "*")
  x <- sweep(x, 2L, seq(-0.4, 0.5, length.out = p), "+")
  if (changed) {
    index <- seq.int(floor(n / 2) + 1L, n)
    x[index, ] <- sweep(
      x[index, , drop = FALSE], 2L,
      seq(0.15, 0.55, length.out = p), "+"
    )
  }
  x
}

.cp_test_cusum_reference <- function(x) {
  x <- as.matrix(x)
  n <- nrow(x)
  k <- seq_len(n - 1L)
  partial <- apply(x, 2L, cumsum)[k, , drop = FALSE]
  raw <- (partial - outer(k / n, colSums(x))) / sqrt(n)
  list(
    c0 = raw,
    chalf = raw / sqrt((k / n) * (1 - k / n)),
    s2 = colSums(apply(x, 2L, diff)^2) / (2 * (n - 1))
  )
}

.cp_test_dms_reference <- function(x) {
  x <- as.matrix(x)
  n <- nrow(x)
  p <- ncol(x)
  leave_variance <- function(excluded) {
    index <- setdiff(seq.int(2L, n), excluded)
    colSums((x[index, , drop = FALSE] -
      x[index - 1L, , drop = FALSE])^2) / (2 * length(index))
  }
  trace.terms <- vapply(seq_len(n - 3L), function(i) {
    variance <- leave_variance(i + 0:3)
    sum((x[i, ] - x[i + 1L, ]) *
      (x[i + 2L, ] - x[i + 3L, ]) / variance)^2
  }, numeric(1))
  trace.hat <- mean(trace.terms) / 4
  fourth.terms <- vapply(seq_len(n - 2L), function(i) {
    variance <- leave_variance(i + 0:2)
    sum((x[i, ] - x[i + 1L, ]) *
      (x[i + 1L, ] - x[i + 2L, ]) / variance)^2
  }, numeric(1))
  list(
    trace.terms = trace.terms,
    fourth.terms = fourth.terms,
    trace.hat = trace.hat,
    fourth.hat = mean(fourth.terms) - 3 * trace.hat
  )
}

.cp_test_ordered_square_sum <- function(u) {
  gram <- tcrossprod(u)
  diag(gram) <- 0
  sum(gram^2)
}

.cp_test_erht_reference <- function(fit, ridge.index = 1L,
                                    candidate = 1L) {
  local <- fit$components$local
  y <- fit$components$pool.signs
  q <- fit$components$ridge.inverse[[ridge.index]]
  n <- nrow(y)
  delta <- local$delta[candidate, ]
  beta <- local$beta[candidate, ]
  a <- y %*% q %*% t(y) / n
  raw <- local$N.effective[[candidate]] *
    drop(delta %*% q %*% delta)
  kappa <- sum(beta^2 * diag(a))
  offdiag <- outer(beta^2, beta^2) * a^2
  diag(offdiag) <- 0
  sigma2 <- 2 * n * sum(offdiag)
  center <- n * kappa
  variance <- n * sigma2
  list(
    A = a, raw = raw, kappa = kappa, sigma2 = sigma2,
    center = center, variance = variance,
    Z = (raw - center) / sqrt(variance)
  )
}


test_that("CUSUM C++ kernels equal a literal R formula", {
  x <- .cp_test_fixture(13L, 4L)
  reference <- .cp_test_cusum_reference(x)
  result <- ch4_cp_cusum_cpp(x)

  expect_equal(result$cusum0, reference$c0, tolerance = 2e-14)
  expect_equal(result$cusum_half, reference$chalf, tolerance = 3e-14)
  expect_equal(result$difference_variance, reference$s2,
               tolerance = 2e-14)
  expect_equal(result$cusum_half,
               result$cusum0 / sqrt((1:12 / 13) * (1 - 1:12 / 13)),
               tolerance = 2e-14)
})


test_that("classical CUSUM statistic and Kolmogorov calibration are literal", {
  x <- c(-1.2, 0.4, 1.1, -0.8, 0.2, 2.1, 2.5, 1.7)
  scale <- 1.7
  fit <- classical_cusum_test(x, sigma = scale, trim = 2L)
  reference <- .cp_test_cusum_reference(matrix(x, ncol = 1L))
  k <- 2:6
  path <- reference$c0[k, 1L] / scale
  selected <- which.max(abs(path))

  expect_true(fit$valid)
  expect_equal(unname(fit$statistic), max(abs(path)), tolerance = 1e-14)
  expect_equal(unname(fit$estimate[["change.point"]]), k[[selected]])
  expect_equal(fit$path$cusum, path, tolerance = 1e-14)
  expect_equal(fit$p.value, .cp_kolmogorov_survival(max(abs(path))),
               tolerance = 1e-15)
  large.reference <- sum(2 * (-1)^(0:19) *
    exp(-2 * (1:20)^2 * 1.2^2))
  small.cdf.reference <- sqrt(2 * pi) / 0.35 * sum(
    exp(-((2 * (1:10) - 1)^2 * pi^2) / (8 * 0.35^2))
  )
  expect_equal(.cp_kolmogorov_survival(1.2), large.reference,
               tolerance = 2e-15)
  expect_equal(.cp_kolmogorov_survival(0.35), 1 - small.cdf.reference,
               tolerance = 2e-15)
})


test_that("classical CUSUM is translation/scale invariant and seeded calibration restores RNG", {
  x <- .cp_test_fixture(15L, 2L)[, 1L]
  base <- classical_cusum_test(
    x, gamma = 0.5, trim = 3L, calibration = "grid-gaussian",
    calibration_draws = 29L, seed = 99L
  )
  transformed <- classical_cusum_test(
    -3 * x + 10, gamma = 0.5, trim = 3L,
    calibration = "grid-gaussian", calibration_draws = 29L, seed = 99L
  )
  set.seed(700L)
  before <- .Random.seed
  invisible(classical_cusum_test(
    x, gamma = 0.5, trim = 3L, calibration = "grid-gaussian",
    calibration_draws = 9L, seed = 5L
  ))

  expect_equal(unname(base$statistic), unname(transformed$statistic),
               tolerance = 2e-14)
  expect_equal(base$p.value, transformed$p.value)
  expect_equal(base$estimate, transformed$estimate)
  expect_identical(.Random.seed, before)
  expect_error(classical_cusum_test(x, gamma = 0.5),
               "Kolmogorov calibration")
})


test_that("classical CUSUM exposes degeneracy without a scale repair", {
  x <- rep(2, 8)
  expect_error(classical_cusum_test(x), "non-positive")
  expect_warning(invalid <- classical_cusum_test(x, strict = FALSE),
                 "non-positive")
  expect_false(invalid$valid)
  expect_null(invalid$estimate)
  expect_true(is.na(invalid$p.value))
})


test_that("DMS leave-four and leave-three estimators match literal index sets", {
  x <- .cp_test_fixture(18L, 5L)
  reference <- .cp_test_dms_reference(x)
  result <- ch4_cp_dms_moments_cpp(x)

  expect_true(result$valid)
  expect_equal(result$trace_terms, reference$trace.terms,
               tolerance = 3e-13)
  expect_equal(result$fourth_terms, reference$fourth.terms,
               tolerance = 3e-13)
  expect_equal(result$trace_hat, reference$trace.hat,
               tolerance = 3e-13)
  expect_equal(result$fourth_hat, reference$fourth.hat,
               tolerance = 3e-13)
  expect_match(result$leaveout_definition, "no bridge")

  first.variance <- colSums((x[5:18, ] - x[4:17, ])^2) / (2 * 14)
  first.term <- sum((x[1, ] - x[2, ]) *
    (x[3, ] - x[4, ]) / first.variance)^2
  expect_equal(result$trace_terms[[1L]], first.term, tolerance = 3e-13)
})


test_that("DMS dense, variance, max and Fisher formulas are self-consistent", {
  x <- .cp_test_fixture(24L, 6L, changed = TRUE)
  fit <- wang_feng_dms_test(x, gamma = 0)
  cusum <- .cp_test_cusum_reference(x)
  moments <- .cp_test_dms_reference(x)
  standardized.half <- sweep(cusum$chalf, 2L, sqrt(cusum$s2), "/")
  dense.by.k <- rowSums(standardized.half^2)
  s <- sum(dense.by.k)
  variance <- ((2 * pi^2 - 18) / 3) * 24^2 * moments$trace.hat +
    ((15 - pi^2) / 3) * 24 * (moments$fourth.hat - 6^2)
  max.path <- sweep(cusum$c0, 2L, sqrt(cusum$s2), "/")
  m <- max(abs(max.path))
  pivot <- 2 * m^2 - log(12)
  p.max <- 1 - exp(-exp(-pivot))
  z <- (s - 26 * 6) / sqrt(variance)
  p.sum <- stats::pnorm(z, lower.tail = FALSE)

  expect_true(fit$valid)
  expect_equal(fit$components$sum$statistic, s, tolerance = 3e-12)
  expect_equal(fit$components$sum$by.location, dense.by.k,
               tolerance = 3e-12)
  expect_equal(fit$components$sum$center, (24 + 2) * 6)
  expect_equal(fit$components$trace.hat, moments$trace.hat,
               tolerance = 3e-12)
  expect_equal(fit$components$fourth.moment.hat, moments$fourth.hat,
               tolerance = 3e-12)
  expect_equal(fit$components$variance.hat, variance, tolerance = 3e-12)
  expect_equal(fit$components$variance.sd^2, variance, tolerance = 3e-12)
  expect_equal(fit$components$max$statistic, m, tolerance = 3e-13)
  expect_equal(fit$components$max$pivot, pivot, tolerance = 3e-13)
  expect_lte(abs(fit$components$max$p.value - p.max), 2e-15)
  expect_equal(fit$components$sum$z, z, tolerance = 3e-13)
  expect_lte(abs(fit$components$sum$p.value - p.sum), 2e-15)
  expect_equal(unname(fit$statistic), -2 * log(p.max * p.sum),
               tolerance = 3e-13)
  expect_lte(
    abs(fit$p.value -
          stats::pchisq(unname(fit$statistic), 4, lower.tail = FALSE)),
    2e-15
  )
})


test_that("weighted DMS locks the primary h, A/D, trimming and right tail", {
  x <- .cp_test_fixture(24L, 6L, changed = TRUE)
  lambda <- 4L
  fit <- wang_feng_dms_test(x, gamma = 0.5, lambda = lambda)
  reference <- .cp_test_cusum_reference(x)
  path <- sweep(reference$chalf, 2L, sqrt(reference$s2), "/")
  m <- max(abs(path[lambda:(24 - lambda), , drop = FALSE]))
  h <- ((lambda / 24)^(-1) - 1)^2
  argument <- 6 * log(h)
  pivot <- sqrt(2 * log(argument)) * m -
    (2 * log(argument) + 0.5 * log(log(argument)) - 0.5 * log(pi))

  expect_equal(fit$components$max$statistic, m, tolerance = 3e-13)
  expect_equal(fit$components$max$pivot, pivot, tolerance = 3e-13)
  expect_equal(fit$components$max$p.value,
               1 - exp(-exp(-pivot)), tolerance = 2e-15)
  expect_equal(fit$components$max$lambda, lambda)
  expect_error(wang_feng_dms_test(x, gamma = 0.5), "lambda")
})


test_that("DMS is common-translation and coordinate-scale invariant", {
  x <- .cp_test_fixture(24L, 6L, changed = TRUE)
  base <- wang_feng_dms_test(x, gamma = 0)
  translated <- wang_feng_dms_test(
    sweep(x, 2L, seq(10, 60, by = 10), "+"), gamma = 0
  )
  scaled <- wang_feng_dms_test(
    sweep(x, 2L, c(0.5, 2, 3, 0.75, 1.4, 4), "*"), gamma = 0
  )

  expect_equal(unname(base$statistic), unname(translated$statistic),
               tolerance = 2e-11)
  expect_equal(unname(base$statistic), unname(scaled$statistic),
               tolerance = 2e-11)
  expect_equal(base$components$variance.hat,
               translated$components$variance.hat, tolerance = 2e-10)
  expect_equal(base$components$variance.hat,
               scaled$components$variance.hat, tolerance = 2e-10)
  expect_equal(base$estimate, translated$estimate)
  expect_equal(base$estimate, scaled$estimate)

  degenerate <- cbind(seq_len(8), rep(1, 8))
  expect_warning(invalid <- wang_feng_dms_test(
    degenerate, strict = FALSE
  ), "non-positive")
  expect_false(invalid$valid)
  expect_null(invalid$estimate)
})


test_that("signed-log Cauchy aggregation handles extreme and reverse endpoints", {
  ordinary <- .cp_cauchy_combine(c(0.2, 0.7), c(0.3, 0.7))
  direct <- sum(c(0.3, 0.7) * tan(pi * (0.5 - c(0.2, 0.7))))
  tiny <- .cp_cauchy_combine(c(1e-300, 0.3))
  cancelling <- .cp_cauchy_combine(c(0.2, 0.8))

  expect_equal(ordinary$statistic, direct, tolerance = 2e-15)
  expect_equal(ordinary$p.value, atan2(1, direct) / pi,
               tolerance = 2e-15)
  expect_true(tiny$log.absolute.statistic > 680)
  expect_true(tiny$p.value > 0)
  expect_equal(cancelling$statistic, 0, tolerance = 2e-15)
  expect_equal(cancelling$p.value, 0.5, tolerance = 2e-15)
  expect_equal(.cp_cauchy_combine(c(0, 0.4))$p.value, 0)
  expect_equal(.cp_cauchy_combine(c(1, 0.6))$p.value, 1)
  expect_error(.cp_cauchy_combine(c(0, 1)), "indeterminate")
})


test_that("scaled HR one-step map and diagnostics are literal", {
  x <- .cp_test_fixture(10L, 4L)
  theta <- colMeans(x)
  diagonal <- apply(x, 2L, stats::var)
  residual <- sweep(sweep(x, 2L, theta, "-"), 2L,
                    sqrt(diagonal), "/")
  radius <- sqrt(rowSums(residual^2))
  u <- residual / radius
  theta.next <- theta + sqrt(diagonal) * colSums(u) / sum(1 / radius)
  diagonal.next <- 4 * diagonal * colMeans(u^2)
  fit <- ch4_cp_scaled_hr_cpp(x, tol = 1e-15, max_iter = 1L)

  expect_false(fit$valid)
  expect_false(fit$stable)
  expect_equal(as.numeric(fit$location), theta.next, tolerance = 3e-14)
  expect_equal(as.numeric(fit$diagonal), diagonal.next,
               tolerance = 3e-14)
  expect_equal(fit$iterations, 1L)
  expect_match(fit$scale_identification, "relative")
  expect_true(is.finite(fit$score_l2))
  expect_true(is.finite(fit$diagonal_residual))
})


test_that("spatial-sign feasible radial, trace and L2 formulas are literal", {
  x <- .cp_test_fixture(30L, 3L, changed = TRUE)
  lambda <- 6L
  m <- 6L
  fit <- spatial_sign_change_point_test(
    x, lambda = lambda, endpoint_fraction = 0.2,
    variant = "unweighted", fv_draws = 39L, seed = 12L
  )
  estimation <- fit$components$estimation
  diagonal <- fit$components$common.diagonal
  first.residual <- sweep(x[1:m, , drop = FALSE], 2L,
                          estimation$first.endpoint.location, "-")
  first.residual <- sweep(first.residual, 2L, sqrt(diagonal), "/")
  last.residual <- sweep(x[(30 - m + 1):30, , drop = FALSE], 2L,
                         estimation$last.endpoint.location, "-")
  last.residual <- sweep(last.residual, 2L, sqrt(diagonal), "/")
  first.radius <- sqrt(rowSums(first.residual^2))
  last.radius <- sqrt(rowSums(last.residual^2))
  first.sign <- first.residual / first.radius
  last.sign <- last.residual / last.radius
  zeta <- 0.5 * (mean(1 / first.radius) + mean(1 / last.radius))
  ordered <- c(
    first = .cp_test_ordered_square_sum(first.sign),
    last = .cp_test_ordered_square_sum(last.sign)
  )
  trace <- 3^2 * sum(ordered) / (2 * m * (m - 1))
  signs <- estimation$full.sample.signs
  partial <- apply(signs, 2L, cumsum)
  total <- colSums(signs)
  l2 <- vapply(lambda:(30 - lambda), function(k) {
    u <- k / 30
    transformed <- sqrt(3 / 30) * (partial[k, ] - u * total)
    sum(transformed^2) - u * (1 - u) * 3
  }, numeric(1))
  s <- (1 - 30^(-0.5)) * max(l2)

  expect_true(fit$valid)
  expect_equal(estimation$first.endpoint.signs, first.sign,
               tolerance = 3e-12)
  expect_equal(estimation$last.endpoint.signs, last.sign,
               tolerance = 3e-12)
  expect_equal(fit$components$zeta1.hat, zeta, tolerance = 3e-12)
  expect_equal(fit$components$trace.ordered.sums, ordered,
               tolerance = 3e-12)
  expect_equal(fit$components$trace.R2.hat, trace, tolerance = 3e-12)
  expect_equal(fit$components$max.L2$path, l2, tolerance = 3e-12)
  expect_equal(fit$components$max.L2$statistic, s, tolerance = 3e-12)
  expect_equal(fit$components$max.L2$pivot, s / sqrt(2 * trace),
               tolerance = 3e-12)
  expect_equal(fit$components$finite.sample.factor, 1 - 30^(-0.5))
  expect_match(fit$diagnostics$trace.denominator, "2\\*m\\*\\(m-1\\)")
})


test_that("spatial-sign max-Linf and weighted L2 pivots use primary constants", {
  x <- .cp_test_fixture(30L, 3L, changed = TRUE)
  unweighted <- spatial_sign_change_point_test(
    x, lambda = 6L, endpoint_fraction = 0.2,
    variant = "unweighted", fv_draws = 19L, seed = 2L
  )
  weighted <- spatial_sign_change_point_test(
    x, lambda = 6L, endpoint_fraction = 0.2,
    variant = "weighted", fv_draws = 19L, seed = 2L
  )
  p <- ncol(x)
  zeta <- unweighted$components$zeta1.hat
  m0 <- unweighted$components$max.Linf$statistic
  expected.max0 <- 2 * p * zeta^2 * m0^2 - log(2 * p)
  h <- ((6 / 30)^(-1) - 1)^2
  expected.max.half <- .cp_A(p * log(h)) * sqrt(p) *
    weighted$components$zeta1.hat *
    weighted$components$max.Linf$statistic - .cp_D(p * log(h))
  argument <- log(30^2 / 6^2)
  expected.l2.half <- .cp_A(argument) *
    abs(weighted$components$max.L2$statistic) /
    sqrt(2 * weighted$components$trace.R2.hat) - .cp_D(argument)

  expect_equal(unweighted$components$max.Linf$pivot, expected.max0,
               tolerance = 3e-12)
  expect_equal(weighted$components$max.Linf$pivot, expected.max.half,
               tolerance = 3e-12)
  expect_equal(weighted$components$max.L2$pivot, expected.l2.half,
               tolerance = 3e-12)
  expect_equal(weighted$components$max.L2$p.value,
               .cp_gumbel_survival(expected.l2.half, factor = 2),
               tolerance = 2e-15)
  expect_equal(unweighted$p.value,
               stats::pchisq(unname(unweighted$statistic), 4,
                              lower.tail = FALSE), tolerance = 2e-15)
})


test_that("spatial-sign scan is translation and diagonal-scale invariant", {
  x <- .cp_test_fixture(30L, 3L, changed = TRUE)
  arguments <- list(
    lambda = 6L, endpoint_fraction = 0.2,
    variant = "unweighted", fv_draws = 19L, seed = 7L
  )
  base <- do.call(spatial_sign_change_point_test, c(list(x = x), arguments))
  translated <- do.call(
    spatial_sign_change_point_test,
    c(list(x = sweep(x, 2L, c(10, -20, 30), "+")), arguments)
  )
  scaled <- do.call(
    spatial_sign_change_point_test,
    c(list(x = sweep(x, 2L, c(2, 0.5, 3), "*")), arguments)
  )

  expect_equal(unname(base$statistic), unname(translated$statistic),
               tolerance = 3e-8)
  expect_equal(unname(base$statistic), unname(scaled$statistic),
               tolerance = 3e-8)
  expect_equal(base$components$max.Linf$pivot,
               translated$components$max.Linf$pivot, tolerance = 3e-8)
  expect_equal(base$components$max.Linf$pivot,
               scaled$components$max.Linf$pivot, tolerance = 3e-8)
  expect_equal(base$components$max.L2$pivot,
               translated$components$max.L2$pivot, tolerance = 3e-8)
  expect_equal(base$components$max.L2$pivot,
               scaled$components$max.L2$pivot, tolerance = 3e-8)
  expect_equal(base$estimate, translated$estimate)
  expect_equal(base$estimate, scaled$estimate)
})


test_that("spatial-sign iteration and degeneracy failures are explicit", {
  x <- .cp_test_fixture(16L, 4L)
  expect_warning(nonconverged <- spatial_sign_change_point_test(
    x, lambda = 4L, endpoint_fraction = 0.25,
    fv_draws = 9L, max_iter = 1L, strict = FALSE
  ), "did not stabilize")
  expect_false(nonconverged$valid)
  expect_null(nonconverged$estimate)

  degenerate <- matrix(rep(1:16, 4L), 16L, 4L)
  expect_warning(invalid <- spatial_sign_change_point_test(
    degenerate, lambda = 4L, endpoint_fraction = 0.25,
    fv_draws = 9L, strict = FALSE
  ), "initializer|zero")
  expect_false(invalid$valid)
  expect_null(invalid$estimate)
})


test_that("modified Weiszfeld median certifies coincident observations", {
  x <- rbind(c(0, 0), c(0, 0), c(1, 0), c(-1, 0))
  fit <- ch4_cp_spatial_median_cpp(x, tol = 1e-12, max_iter = 100L)

  expect_true(fit$valid)
  expect_true(fit$stable)
  expect_equal(as.numeric(fit$location), c(0, 0), tolerance = 1e-12)
  expect_gte(fit$coincident_observations, 2L)
  expect_equal(fit$score_residual, 0, tolerance = 1e-14)
  expect_match(fit$convergence_basis, "Weiszfeld")
})


test_that("ERHT local raw, center, ordered variance and Z are literal", {
  x <- .cp_test_fixture(14L, 4L, changed = TRUE)
  fit <- erht_change_point_test(
    x, ridge = c(0.4, 0.9), epsilon = 0.25,
    calibration = "none"
  )
  reference <- .cp_test_erht_reference(fit, 2L, 3L)
  local <- fit$components$local

  expect_true(fit$valid)
  expect_equal(local$raw[3L, 2L], reference$raw, tolerance = 3e-12)
  expect_equal(local$kappa[3L, 2L], reference$kappa,
               tolerance = 3e-12)
  expect_equal(local$sigma2[3L, 2L], reference$sigma2,
               tolerance = 3e-12)
  expect_equal(local$center[3L, 2L], reference$center,
               tolerance = 3e-12)
  expect_equal(local$variance[3L, 2L], reference$variance,
               tolerance = 3e-12)
  expect_equal(local$Z[3L, 2L], reference$Z, tolerance = 3e-12)
  expect_equal(local$center[3L, 2L], nrow(x) * local$kappa[3L, 2L],
               tolerance = 2e-15)
  expect_equal(local$variance[3L, 2L], nrow(x) * local$sigma2[3L, 2L],
               tolerance = 2e-15)
  expect_equal(fit$components$spatial.sign.covariance,
               crossprod(fit$components$pool.signs) / nrow(x),
               tolerance = 3e-14)
  expect_equal(rowSums(fit$components$pool.signs^2), rep(4, nrow(x)),
               tolerance = 3e-10)
})


test_that("ERHT adjacent candidates and temporal covariance match formulas", {
  single <- .cp_erht_candidates(12L, 0.25, "single")
  multiple <- .cp_erht_candidates(12L, 0.25, "multiple")
  covariance <- .cp_erht_temporal_covariance(single)
  phi <- function(row) {
    grid <- (seq_len(12L) - 0.5) / 12
    as.numeric(grid >= row$t2 & grid < row$t3) /
      (row$t3 - row$t2) -
      as.numeric(grid >= row$t1 & grid < row$t2) /
      (row$t2 - row$t1)
  }
  phi.matrix <- t(vapply(seq_len(nrow(single)), function(i) {
    phi(single[i, ])
  }, numeric(12)))
  inner <- tcrossprod(phi.matrix) / 12
  reference <- inner^2 / sqrt(outer(diag(inner)^2, diag(inner)^2))

  expect_equal(single$q1, rep(0L, 7L))
  expect_equal(single$q2, 3:9)
  expect_equal(single$q3, rep(12L, 7L))
  expect_equal(nrow(multiple), choose(5L, 3L))
  expect_true(all(multiple$q2 - multiple$q1 >= 3L))
  expect_true(all(multiple$q3 - multiple$q2 >= 3L))
  expect_equal(covariance, reference, tolerance = 3e-14)
  expect_equal(diag(covariance), rep(1, nrow(single)))
})


test_that("ERHT Gaussian calibration is fixed-seed reproducible and explicit", {
  x <- .cp_test_fixture(14L, 4L, changed = TRUE)
  set.seed(912L)
  before <- .Random.seed
  first <- erht_change_point_test(
    x, ridge = c(0.4, 0.9), epsilon = 0.25,
    calibration_draws = 29L, seed = 44L, keep_calibration = TRUE
  )
  after <- .Random.seed
  second <- erht_change_point_test(
    x, ridge = c(0.4, 0.9), epsilon = 0.25,
    calibration_draws = 29L, seed = 44L, keep_calibration = TRUE
  )

  expect_identical(after, before)
  expect_equal(first$p.value, second$p.value)
  expect_equal(first$components$ridge.results,
               second$components$ridge.results)
  expect_equal(first$diagnostics$calibration.reference,
               second$diagnostics$calibration.reference)
  expect_false(first$diagnostics$calibration.exact)
  expect_match(first$components$combined$exact.scope, "r_E")
  expect_equal(nrow(first$components$ridge.results), 2L)
  expect_equal(length(first$diagnostics$calibration.reference), 29L)
  expect_true(all(first$components$ridge.results$p.value >= 1 / 30))
})


test_that("ERHT is translation, global-scale and signed-permutation invariant", {
  x <- .cp_test_fixture(14L, 4L, changed = TRUE)
  base <- erht_change_point_test(
    x, ridge = 0.6, epsilon = 0.25,
    calibration_draws = 19L, seed = 8L
  )
  translated <- erht_change_point_test(
    sweep(x, 2L, c(20, -10, 5, 30), "+"), ridge = 0.6,
    epsilon = 0.25, calibration_draws = 19L, seed = 8L
  )
  scaled <- erht_change_point_test(
    7 * x, ridge = 0.6, epsilon = 0.25,
    calibration_draws = 19L, seed = 8L
  )
  transformed.x <- sweep(
    x[, c(3, 1, 4, 2), drop = FALSE],
    2L, c(-1, 1, -1, 1), "*"
  )
  transformed <- erht_change_point_test(
    transformed.x, ridge = 0.6, epsilon = 0.25,
    calibration_draws = 19L, seed = 8L
  )

  expect_equal(unname(base$statistic), unname(translated$statistic),
               tolerance = 2e-8)
  expect_equal(unname(base$statistic), unname(scaled$statistic),
               tolerance = 2e-8)
  expect_equal(unname(base$statistic), unname(transformed$statistic),
               tolerance = 2e-8)
  expect_equal(base$p.value, translated$p.value)
  expect_equal(base$p.value, scaled$p.value)
  expect_equal(base$p.value, transformed$p.value)
  expect_equal(base$estimate, translated$estimate)
  expect_equal(base$estimate, scaled$estimate)
  expect_equal(base$estimate, transformed$estimate)
})


test_that("ERHT multiple scan works and no-calibration contract is honest", {
  x <- .cp_test_fixture(16L, 4L, changed = TRUE)
  fit <- erht_change_point_test(
    x, ridge = 0.5, scan = "multiple", epsilon = 0.25,
    calibration = "none"
  )

  expect_true(fit$valid)
  expect_true(nrow(fit$components$local$candidates) > 1L)
  expect_true(is.finite(unname(fit$statistic)))
  expect_true(is.na(fit$p.value))
  expect_true(is.na(fit$components$reject))
  expect_true(is.na(fit$diagnostics$calibration.exact))
  expect_equal(fit$diagnostics$resolved.ridge, 0.5)
  expect_match(fit$diagnostics$ridge.parameterization, "Actual rho")
})


test_that("ERHT nonconvergence and coincident residual failures return no estimate", {
  x <- .cp_test_fixture(12L, 3L)
  expect_warning(nonconverged <- erht_change_point_test(
    x, ridge = 0.5, epsilon = 0.25, calibration = "none",
    max_iter = 1L, strict = FALSE
  ), "did not stabilize")
  expect_false(nonconverged$valid)
  expect_null(nonconverged$estimate)

  identical.rows <- matrix(1, 12L, 3L)
  expect_warning(invalid <- erht_change_point_test(
    identical.rows, ridge = 0.5, epsilon = 0.25,
    calibration = "none", strict = FALSE
  ), "inverse-distance|variance|residual")
  expect_false(invalid$valid)
  expect_null(invalid$estimate)
  expect_error(erht_change_point_test(x, ridge = 0), "positive finite")
})


test_that("ERHT moments use ordered off-diagonal pairs without subtraction", {
  a <- matrix(c(
    1.2, 0.3, -0.2,
    0.3, 0.9, 0.4,
    -0.2, 0.4, 1.1
  ), 3L, 3L, byrow = TRUE)
  beta <- c(-0.5, 0.2, 0.7)
  fit <- ch4_cp_erht_moments_cpp(a, beta)
  reference.kappa <- sum(beta^2 * diag(a))
  terms <- outer(beta^2, beta^2) * a^2
  diag(terms) <- 0
  offdiag <- sum(terms)

  expect_true(fit$valid)
  expect_equal(fit$kappa, reference.kappa, tolerance = 2e-15)
  expect_equal(fit$ordered_offdiagonal_sum, offdiag, tolerance = 2e-15)
  expect_equal(fit$sigma2, 2 * 3 * offdiag, tolerance = 2e-15)

  diagonal <- ch4_cp_erht_moments_cpp(diag(3), beta)
  expect_false(diagonal$valid)
  expect_match(diagonal$failure, "non-positive")
})


test_that("WBS triple collection equals the primary integer constraints", {
  triples <- .cp_wbs_candidate_triples(8L, 0.25, 1L)
  literal <- do.call(rbind, lapply(2:6, function(q2) {
    rows <- list()
    id <- 1L
    for (q1 in 0:(q2 - 2L)) for (q3 in (q2 + 2L):8L) {
      rows[[id]] <- c(q1 = q1, q2 = q2, q3 = q3)
      id <- id + 1L
    }
    do.call(rbind, rows)
  }))

  expect_equal(as.matrix(triples[, c("q1", "q2", "q3")]), literal)
  expect_true(all(triples$q2 - triples$q1 >= 2L))
  expect_true(all(triples$q3 - triples$q2 >= 2L))
  expect_equal(nrow(.cp_wbs_candidate_triples(8L, 0.25, 2L)), 22L)
})


test_that("WBS applies strict threshold, refinement and supplied tuning", {
  x <- .cp_test_fixture(12L, 3L, changed = TRUE)
  controls <- .cp_iteration_controls(1e-8, 1000L, 0)
  score <- .cp_wbs_interval_score(
    x, 1L, 12L, ridge = 0.5, epsilon = 0.25,
    triple_step = 1L, controls = controls
  )
  intervals <- matrix(c(1L, 12L), ncol = 2L)
  exact <- erht_wbs(
    x, ridge = 0.5, threshold = score$score,
    min_interval = 8L, refinement_radius = 3L,
    deletion_radius = 1L, intervals = intervals,
    epsilon = 0.25, max_changes = 1L
  )
  below <- erht_wbs(
    x, ridge = 0.5, threshold = score$score - 1e-10,
    min_interval = 8L, refinement_radius = 3L,
    deletion_radius = 1L, intervals = intervals,
    epsilon = 0.25, max_changes = 1L
  )

  expect_true(exact$valid)
  expect_length(exact$estimate, 0L)
  expect_true(below$valid)
  expect_length(below$estimate, 1L)
  expect_equal(below$tuning$ridge, 0.5)
  expect_equal(below$tuning$threshold, score$score - 1e-10)
  expect_equal(below$tuning$refinement.radius, 3L)
  expect_equal(below$tuning$deletion.radius, 1L)
  expect_true(below$tuning$exhaustive.triples)
  expect_match(below$diagnostics$threshold.rule, "strict")
  expect_match(below$diagnostics$threshold.source, "user supplied")
})


test_that("WBS random intervals are reproducible and failures are not partial fits", {
  x <- .cp_test_fixture(12L, 3L, changed = TRUE)
  first <- erht_wbs(
    x, ridge = 0.5, threshold = 100,
    min_interval = 8L, refinement_radius = 3L,
    deletion_radius = 1L, M = 9L, seed = 3L,
    epsilon = 0.25
  )
  second <- erht_wbs(
    x, ridge = 0.5, threshold = 100,
    min_interval = 8L, refinement_radius = 3L,
    deletion_radius = 1L, M = 9L, seed = 3L,
    epsilon = 0.25
  )

  expect_equal(first$intervals, second$intervals)
  expect_equal(first$estimate, second$estimate)
  expect_equal(first$diagnostics$interval.source,
               "uniform unordered endpoint pairs")
  expect_true(all(first$intervals[, 1L] < first$intervals[, 2L]))
  expect_error(erht_wbs(
    x, ridge = 0.5, threshold = 1, min_interval = 8L,
    refinement_radius = 3L, deletion_radius = 1L,
    epsilon = 0.25
  ), "Supply either")

  expect_warning(invalid <- erht_wbs(
    matrix(1, 12L, 3L), ridge = 0.5, threshold = 1,
    min_interval = 8L, refinement_radius = 3L,
    deletion_radius = 1L, intervals = matrix(c(1, 12), ncol = 2),
    epsilon = 0.25, strict = FALSE
  ), "inverse-distance|variance|residual")
  expect_false(invalid$valid)
  expect_null(invalid$estimate)
})


test_that("module helpers reject malformed inputs and apply no hidden repair", {
  x <- .cp_test_fixture(10L, 3L)
  expect_error(classical_cusum_test(cbind(x[, 1], x[, 2])), "univariate")
  expect_error(wang_feng_dms_test(x[1:5, ]), "at least 6")
  expect_error(spatial_sign_change_point_test(
    x, lambda = 6L, endpoint_fraction = 0.2, fv_draws = 9L
  ), "no spatial-sign scan")
  expect_error(erht_change_point_test(
    x, ridge = NA_real_, calibration = "none"
  ), "positive finite")
  expect_error(erht_change_point_test(
    x, ridge = 0.5, epsilon = 0.5, scan = "multiple",
    calibration = "none"
  ), "permitted finite range")
  expect_error(.cp_psd_root(matrix(c(1, 2, 2, 1), 2L)),
               "not positive semidefinite")
  expect_equal(.cp_row_norms(matrix(c(1e200, 1e200), nrow = 1L)),
               sqrt(2) * 1e200, tolerance = 2e-15)
})
