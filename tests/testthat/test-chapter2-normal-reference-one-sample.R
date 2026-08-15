zzg22nr_fixture <- function() {
  rbind(
    c(-1.2, 0.4, 1.1),
    c(0.3, -0.7, 0.5),
    c(1.4, 0.8, -0.2),
    c(-0.6, 1.5, 0.9),
    c(0.9, -1.1, 1.3),
    c(1.7, 0.2, -0.8)
  )
}


zzg22nr_reference <- function(x, mu = numeric(ncol(x)), alpha = 0.05) {
  z <- sweep(x, 2L, mu, check.margin = FALSE)
  n <- nrow(z)
  v <- n - 1
  S <- stats::cov(z)
  if (ncol(z) == 1L) {
    S <- matrix(S, 1L, 1L)
  }
  t1 <- sum(diag(S))
  t2 <- sum(S^2)
  t3 <- sum(diag(S %*% S %*% S))
  oracle <- n * sum(colMeans(z)^2)
  T <- oracle - t1
  T.u <- 2 / v * sum(vapply(
    seq_len(n - 1L),
    function(i) sum(z[i, , drop = TRUE] *
                      colSums(z[(i + 1L):n, , drop = FALSE])),
    numeric(1L)
  ))
  a1 <- t1
  a2.bracket <- t2 - t1^2 / v
  a2.factor <- v^2 / ((v - 1) * (v + 2))
  a2 <- a2.factor * a2.bracket
  a3.bracket <- t3 - 3 * t1 * t2 / v + 2 * t1^3 / v^2
  a3.factor <- v^4 / ((v - 1) * (v + 4) * (v^2 - 4))
  a3 <- a3.factor * a3.bracket
  kappa2 <- 2 * n / v * a2
  kappa3 <- 8 * n * (n - 2) / v^2 * a3
  beta0 <- -n / (n - 2) * a2^2 / a3
  beta1 <- (n - 2) / v * a3 / a2
  df <- n * v / (n - 2)^2 * a2^3 / a3^2
  argument <- (T - beta0) / beta1
  critical.argument <- stats::qchisq(alpha, df, lower.tail = FALSE)
  list(
    z = z, n = n, v = v, S = S, t1 = t1, t2 = t2, t3 = t3,
    oracle = oracle, T = T, T.u = T.u, a1 = a1,
    a2.bracket = a2.bracket, a2.factor = a2.factor, a2 = a2,
    a3.bracket = a3.bracket, a3.factor = a3.factor, a3 = a3,
    kappa2 = kappa2, kappa3 = kappa3,
    beta0 = beta0, beta1 = beta1, df = df, argument = argument,
    p.value = stats::pchisq(argument, df, lower.tail = FALSE),
    log.p.value = stats::pchisq(
      argument, df, lower.tail = FALSE, log.p = TRUE
    ),
    critical.argument = critical.argument,
    critical.T = beta0 + beta1 * critical.argument
  )
}


test_that("Zhang-Zhou-Guo matches every feasible one-sample formula", {
  x <- zzg22nr_fixture()
  mu <- c(0.15, -0.2, 0.1)
  reference <- zzg22nr_reference(x, mu, alpha = 0.08)
  result <- zhang_zhou_guo_one_sample_test(x, mu, alpha = 0.08)
  component <- result$components

  expect_s3_class(result, "hd_location_test")
  expect_s3_class(result, "htest")
  expect_equal(reference$T, reference$T.u, tolerance = 2e-15)
  expect_equal(component$oracle.Q.n.input, reference$oracle,
               tolerance = 3e-14)
  expect_equal(component$T.centered.input, reference$T,
               tolerance = 3e-14)
  expect_equal(component$T.by.oracle.minus.trace.scaled,
               component$T.centered.scaled, tolerance = 3e-15)
  expect_lt(component$T.identity.error.scaled, 2e-15)
  expect_gt(abs(component$oracle.Q.n.input - component$T.centered.input),
            0.1)
  expect_equal(component$trace.S.input, reference$t1,
               tolerance = 3e-14)
  expect_equal(component$trace.S2.input, reference$t2,
               tolerance = 5e-14)
  expect_equal(component$trace.S3.input, reference$t3,
               tolerance = 7e-14)
  expect_equal(component$a1.hat.input, reference$a1,
               tolerance = 3e-14)
  expect_equal(component$a2.bracket.scaled /
                 component$trace.S2.scaled,
               reference$a2.bracket / reference$t2,
               tolerance = 4e-14)
  expect_equal(component$a2.factor, reference$a2.factor,
               tolerance = 2e-15)
  expect_equal(component$a2.hat.input, reference$a2,
               tolerance = 8e-14)
  expect_equal(component$a3.bracket.scaled /
                 component$trace.S3.scaled,
               reference$a3.bracket / reference$t3,
               tolerance = 8e-14)
  expect_equal(component$a3.factor, reference$a3.factor,
               tolerance = 2e-15)
  expect_equal(component$a3.hat.input, reference$a3,
               tolerance = 2e-13)
  expect_equal(component$feasible.reference.cumulants.input[["kappa1"]],
               0, tolerance = 0)
  expect_equal(component$feasible.reference.cumulants.input[["kappa2"]],
               reference$kappa2, tolerance = 1e-13)
  expect_equal(component$feasible.reference.cumulants.input[["kappa3"]],
               reference$kappa3, tolerance = 3e-13)
  expect_equal(component$beta0.input, reference$beta0,
               tolerance = 2e-13)
  expect_equal(component$beta1.input, reference$beta1,
               tolerance = 2e-13)
  expect_equal(component$df, reference$df, tolerance = 3e-13)
  expect_equal(unname(result$statistic), reference$argument,
               tolerance = 2e-13)
  expect_equal(unname(result$parameter), reference$df,
               tolerance = 3e-13)
  expect_equal(result$p.value, reference$p.value, tolerance = 2e-15)
  expect_equal(component$log.p.value, reference$log.p.value,
               tolerance = 2e-15)
  expect_equal(component$critical.chi.square, reference$critical.argument,
               tolerance = 2e-14)
  expect_equal(component$critical.T.input, reference$critical.T,
               tolerance = 3e-13)
  expect_identical(component$reject,
                   isTRUE(reference$argument > reference$critical.argument))
  expect_identical(result$null.distribution$tail, "upper")
  expect_false(result$diagnostics$motivating.oracle.Q.used.as.test.statistic)
  expect_true(result$diagnostics$centred.feasible.statistic)
  expect_true(result$diagnostics$trace.a3.denominator.uses.v.plus.4)
  expect_false(result$diagnostics$nonnormal.third.cumulant.Upsilon.estimated)
  expect_true(result$diagnostics$orthogonal.coordinate.invariant)
  expect_false(result$diagnostics$coordinatewise.scale.invariant)
  expect_identical(result$diagnostics$regularization, "none")
  expect_identical(result$diagnostics$trace.repair, "none")
})


test_that("the appendix v+4 third-trace denominator is locked", {
  x <- zzg22nr_fixture()
  reference <- zzg22nr_reference(x)
  result <- zhang_zhou_guo_one_sample_test(x)
  v <- nrow(x) - 1
  wrong.factor <- v^4 / ((v - 1) * (v + 3) * (v^2 - 4))
  wrong.a3 <- wrong.factor * reference$a3.bracket

  expect_equal(result$components$a3.factor,
               v^4 / ((v - 1) * (v + 4) * (v^2 - 4)),
               tolerance = 0)
  expect_gt(abs(result$components$a3.hat.input - wrong.a3), 1e-3)
  expect_equal(
    result$components$beta0.scaled +
      result$components$beta1.scaled * result$components$df,
    0,
    tolerance = 3e-15
  )
})


test_that("row, null-translation, orthogonal, and common-scale properties hold", {
  x <- zzg22nr_fixture()
  mu <- c(0.15, -0.2, 0.1)
  baseline <- zhang_zhou_guo_one_sample_test(x, mu)

  row.result <- zhang_zhou_guo_one_sample_test(
    x[c(6, 2, 5, 1, 4, 3), , drop = FALSE], mu
  )
  expect_equal(row.result$statistic, baseline$statistic,
               tolerance = 5e-14)
  expect_equal(row.result$p.value, baseline$p.value, tolerance = 2e-15)

  shift <- c(3.2, -1.7, 5.1)
  translated <- zhang_zhou_guo_one_sample_test(
    sweep(x, 2L, shift, "+"), mu + shift
  )
  expect_equal(translated$statistic, baseline$statistic,
               tolerance = 2e-13)
  expect_equal(translated$p.value, baseline$p.value, tolerance = 2e-14)

  angle <- 0.37
  Q <- matrix(c(
    cos(angle), -sin(angle), 0,
    sin(angle), cos(angle), 0,
    0, 0, 1
  ), 3L, 3L, byrow = TRUE)
  rotated <- zhang_zhou_guo_one_sample_test(x %*% Q, drop(mu %*% Q))
  expect_equal(rotated$statistic, baseline$statistic,
               tolerance = 8e-14)
  expect_equal(rotated$p.value, baseline$p.value, tolerance = 5e-15)

  for (multiplier in c(-7, 1e-200, 1e200)) {
    scaled <- zhang_zhou_guo_one_sample_test(x * multiplier, mu * multiplier)
    expect_equal(scaled$statistic, baseline$statistic,
                 tolerance = 2e-13)
    expect_equal(scaled$p.value, baseline$p.value, tolerance = 2e-14)
    expect_equal(scaled$components$df, baseline$components$df,
                 tolerance = 3e-13)
  }
  huge <- zhang_zhou_guo_one_sample_test(x * 1e200, mu * 1e200)
  tiny <- zhang_zhou_guo_one_sample_test(x * 1e-200, mu * 1e-200)
  expect_true(all(is.finite(c(
    huge$components$T.centered.scaled,
    huge$components$a2.hat.scaled,
    huge$components$a3.hat.scaled,
    tiny$components$T.centered.scaled,
    tiny$components$a2.hat.scaled,
    tiny$components$a3.hat.scaled
  ))))
  expect_true(is.na(huge$components$T.centered.input))
  expect_true(is.na(tiny$components$T.centered.input))
  expect_false(huge$diagnostics$input.representability$T)
  expect_false(tiny$diagnostics$input.representability$T)
})


test_that("minimum n, p=1, subnormal values, and rank deficiency work", {
  x <- matrix(c(-2, -0.5, 1, 3), ncol = 1L)
  reference <- zzg22nr_reference(x)
  result <- zhang_zhou_guo_one_sample_test(x)
  expect_equal(result$components$T.centered.input, reference$T,
               tolerance = 2e-14)
  expect_equal(result$components$a2.hat.input, reference$a2,
               tolerance = 3e-14)
  expect_equal(result$components$a3.hat.input, reference$a3,
               tolerance = 5e-14)
  expect_equal(result$p.value, reference$p.value, tolerance = 2e-15)

  base <- matrix(c(-3, -1, 1, 2, 4, 6), ncol = 1L)
  ordinary <- zhang_zhou_guo_one_sample_test(base)
  smallest <- 2^-1074
  subnormal <- zhang_zhou_guo_one_sample_test(base * smallest)
  expect_equal(subnormal$statistic, ordinary$statistic,
               tolerance = 3e-13)
  expect_equal(subnormal$p.value, ordinary$p.value, tolerance = 3e-14)
  expect_true(is.finite(subnormal$components$T.centered.scaled))
  expect_true(is.na(subnormal$components$T.centered.input))

  rank.deficient <- cbind(
    zzg22nr_fixture(),
    zzg22nr_fixture()[, 1] + zzg22nr_fixture()[, 2],
    2 * zzg22nr_fixture()[, 3]
  )
  allowed <- zhang_zhou_guo_one_sample_test(rank.deficient)
  expect_true(is.finite(allowed$statistic))
  expect_true(is.finite(allowed$p.value))
  expect_identical(allowed$diagnostics$trace.computation, "primal")

  dual.data <- outer(seq_len(6), seq_len(11), function(i, j) {
    sin(i * j / 5) + (i - 3.5) * (j + 1)^2 / 200
  })
  dual.reference <- zzg22nr_reference(dual.data)
  dual <- zhang_zhou_guo_one_sample_test(dual.data)
  expect_identical(dual$diagnostics$trace.computation, "dual")
  expect_false(dual$diagnostics$constructs.p.by.p.matrix)
  expect_equal(dual$components$T.centered.input, dual.reference$T,
               tolerance = 2e-12)
  expect_equal(dual$components$trace.S2.input, dual.reference$t2,
               tolerance = 5e-12)
  expect_equal(dual$components$trace.S3.input, dual.reference$t3,
               tolerance = 2e-11)
  expect_equal(dual$components$a2.hat.input, dual.reference$a2,
               tolerance = 8e-12)
  expect_equal(dual$components$a3.hat.input, dual.reference$a3,
               tolerance = 3e-11)
  expect_equal(dual$p.value, dual.reference$p.value, tolerance = 2e-13)
  expect_true(is.finite(dual$p.value))
})


test_that("overflow-safe null subtraction remains operational", {
  x <- matrix(c(
    1.7e308, 1.3e308, 0.9e308, 0.2e308,
    -0.4e308, -1.1e308
  ), ncol = 1L)
  result <- zhang_zhou_guo_one_sample_test(x, mu = -1.7e308)

  expect_true(result$diagnostics$subtraction.overflow.fallback)
  expect_true(is.finite(result$statistic))
  expect_true(is.finite(result$p.value))
  expect_true(is.finite(result$components$T.centered.scaled))
  expect_true(is.na(result$components$T.centered.input))
})


test_that("paired wrapper is exactly the rowwise one-sample transformation", {
  x <- zzg22nr_fixture()
  perturbation <- rbind(
    c(0.2, -0.1, 0.3), c(-0.4, 0.2, -0.2),
    c(0.1, 0.5, -0.1), c(-0.3, -0.2, 0.4),
    c(0.5, 0.1, 0.2), c(-0.2, -0.4, -0.3)
  )
  y <- x - perturbation
  delta <- c(0.05, -0.1, 0.02)
  paired <- zhang_zhou_guo_paired_test(x, y, delta)
  direct <- zhang_zhou_guo_one_sample_test(x - y, delta)

  expect_equal(paired$statistic, direct$statistic, tolerance = 8e-14)
  expect_equal(paired$p.value, direct$p.value, tolerance = 5e-15)
  expect_equal(paired$components$T.centered.input,
               direct$components$T.centered.input, tolerance = 4e-14)
  expect_equal(unname(paired$estimate), unname(colMeans(x - y)),
               tolerance = 2e-14)
  expect_match(paired$diagnostics$transformation, "paired difference")

  shift <- c(1.2, -3.4, 5.6)
  translated <- zhang_zhou_guo_paired_test(
    sweep(x, 2L, shift, "+"), sweep(y, 2L, shift, "+"), delta
  )
  expect_equal(translated$statistic, paired$statistic,
               tolerance = 2e-13)
  expect_equal(translated$p.value, paired$p.value, tolerance = 2e-14)
})


test_that("linear wrapper is exactly the documented same-unit transform", {
  x <- zzg22nr_fixture()
  colnames(x) <- c("a", "b", "c")
  L <- rbind(sum = c(1, 0.5, -0.2), difference = c(0, 1, -1))
  colnames(L) <- colnames(x)
  rhs <- c(0.1, -0.2)
  linear <- zhang_zhou_guo_linear_hypothesis_test(x, L, rhs)
  z <- x %*% t(L)
  direct <- zhang_zhou_guo_one_sample_test(z, rhs)

  expect_equal(linear$statistic, direct$statistic, tolerance = 1e-13)
  expect_equal(linear$p.value, direct$p.value, tolerance = 8e-15)
  expect_equal(linear$components$T.centered.input,
               direct$components$T.centered.input, tolerance = 7e-14)
  expect_equal(unname(linear$estimate), unname(drop(L %*% colMeans(x))),
               tolerance = 3e-14)
  expect_identical(names(linear$estimate), rownames(L))
  expect_false(linear$diagnostics$independent.group.MANOVA)
  expect_match(linear$diagnostics$linear.hypothesis.scope, "iid")

  overflow.safe <- zhang_zhou_guo_linear_hypothesis_test(
    x * 1e200, L * 1e200, rhs = numeric(nrow(L))
  )
  zero.rhs <- zhang_zhou_guo_linear_hypothesis_test(
    x, L, rhs = numeric(nrow(L))
  )
  expect_equal(overflow.safe$statistic, zero.rhs$statistic,
               tolerance = 3e-13)
  expect_equal(overflow.safe$p.value, zero.rhs$p.value,
               tolerance = 3e-14)
  expect_true(overflow.safe$diagnostics$subtraction.overflow.fallback)
})


test_that("invalid inputs and trace degeneracy fail without repair", {
  x <- zzg22nr_fixture()
  expect_error(zhang_zhou_guo_one_sample_test(x[1:3, ]),
               "at least 4 row")
  expect_error(zhang_zhou_guo_one_sample_test(x, mu = c(0, 1)),
               "length 3")
  expect_error(zhang_zhou_guo_one_sample_test(x, alpha = 0),
               "strictly between")
  expect_error(zhang_zhou_guo_one_sample_test(x, alpha = 1),
               "strictly between")
  expect_error(zhang_zhou_guo_one_sample_test(x, alpha = c(0.1, 0.2)),
               "one finite")

  expect_error(
    zhang_zhou_guo_one_sample_test(matrix(1, 6, 3)),
    "positive sample variation|trace\\(Sigma"
  )

  helmert <- stats::contr.helmert(6)
  helmert <- sweep(helmert, 2L, sqrt(colSums(helmert^2)), "/")
  expect_error(
    zhang_zhou_guo_one_sample_test(helmert),
    "trace\\(Sigma\\^2\\) estimate"
  )
  expect_error(
    zhang_zhou_guo_one_sample_test(helmert[, 1:4, drop = FALSE]),
    "trace\\(Sigma\\^3\\) estimate"
  )

  expect_error(zhang_zhou_guo_paired_test(x, x[-1, ]),
               "same number of rows|at least 4 row")
  expect_error(zhang_zhou_guo_paired_test(x, x[, 1:2]),
               "same number of columns")
  expect_error(zhang_zhou_guo_paired_test(x, x, delta = c(0, 1)),
               "length 3")

  expect_error(zhang_zhou_guo_linear_hypothesis_test(x, matrix(1, 2, 2)),
               "3 columns")
  expect_error(zhang_zhou_guo_linear_hypothesis_test(
    x, matrix(NA_real_, 1, 3)
  ),
               "finite")
  expect_error(zhang_zhou_guo_linear_hypothesis_test(
    x, matrix(1, 2, 3), rhs = 0
  ), "length 2")

  named.x <- x
  colnames(named.x) <- c("a", "b", "c")
  named.L <- matrix(1, 1, 3, dimnames = list(NULL, c("b", "a", "c")))
  expect_error(zhang_zhou_guo_linear_hypothesis_test(named.x, named.L),
               "match columns")
  expect_error(zhang_zhou_guo_linear_hypothesis_test(
    x, matrix(0, 2, 3), rhs = c(0, 0)
  ), "nonzero|null-centred|positive sample variation")
})
