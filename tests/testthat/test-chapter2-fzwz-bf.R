fzwz_variance <- function(z) {
  apply(z, 2L, stats::var)
}

fzwz_third <- function(z) {
  centered <- sweep(z, 2L, colMeans(z), "-")
  colMeans(centered^3)
}

fzwz_within_reference <- function(target, other, gamma,
                                  target_is_group1) {
  n <- nrow(target)
  result <- 0
  for (i1 in seq_len(n)) {
    for (i2 in setdiff(seq_len(n), i1)) {
      for (i3 in setdiff(seq_len(n), c(i1, i2))) {
        for (i4 in setdiff(seq_len(n), c(i1, i2, i3))) {
          keep <- -c(i1, i2, i3, i4)
          leave.variance <- fzwz_variance(target[keep, , drop = FALSE])
          denominator <- if (target_is_group1) {
            leave.variance + gamma * fzwz_variance(other)
          } else {
            fzwz_variance(other) + gamma * leave.variance
          }
          first <- sum(
            (target[i1, ] - target[i2, ]) *
              (target[i3, ] - target[i4, ]) / denominator
          )
          second <- sum(
            (target[i3, ] - target[i2, ]) *
              (target[i1, ] - target[i4, ]) / denominator
          )
          result <- result + first * second
        }
      }
    }
  }
  result / (2 * n * (n - 1) * (n - 2) * (n - 3))
}

fzwz_cross_reference <- function(x, y, gamma) {
  n1 <- nrow(x)
  n2 <- nrow(y)
  result <- 0
  for (i1 in seq_len(n1)) {
    for (i2 in setdiff(seq_len(n1), i1)) {
      variance1 <- fzwz_variance(x[-c(i1, i2), , drop = FALSE])
      for (j1 in seq_len(n2)) {
        for (j2 in setdiff(seq_len(n2), j1)) {
          variance2 <- fzwz_variance(y[-c(j1, j2), , drop = FALSE])
          denominator <- variance1 + gamma * variance2
          inner <- sum(
            (x[i1, ] - x[i2, ]) * (y[j1, ] - y[j2, ]) /
              denominator
          )
          result <- result + inner^2
        }
      }
    }
  }
  result / (4 * n1 * (n1 - 1) * n2 * (n2 - 1))
}

fzwz_reference <- function(x, y) {
  n1 <- nrow(x)
  n2 <- nrow(y)
  gamma <- n1 / n2
  variance1 <- fzwz_variance(x)
  variance2 <- fzwz_variance(y)
  third1 <- fzwz_third(x)
  third2 <- fzwz_third(y)
  denominator <- variance1 + gamma * variance2
  difference <- colMeans(x) - colMeans(y)
  A <- difference^2 - variance1 / n1 - variance2 / n2
  T.BF <- sum(A / denominator)
  b1 <- sum(
    2 * variance1^2 / (n1 * (n1 - 1) * denominator^2) +
      2 * gamma * variance2^2 /
        (n2 * (n2 - 1) * denominator^2)
  )
  b2 <- sum(
    2 * (third1 / n1 - gamma * third2 / n2)^2 /
      denominator^3
  )
  trace1 <- fzwz_within_reference(x, y, gamma, TRUE)
  trace2 <- fzwz_within_reference(y, x, gamma, FALSE)
  trace12 <- fzwz_cross_reference(x, y, gamma)
  coefficients <- c(
    group1 = 2 / (n1 * (n1 - 1)),
    group2 = 2 / (n2 * (n2 - 1)),
    cross = 4 / (n1 * n2)
  )
  null.variance <- coefficients[[1L]] * trace1 +
    coefficients[[2L]] * trace2 + coefficients[[3L]] * trace12
  z <- (T.BF - b1 - b2) / sqrt(null.variance)

  list(
    mean.x = colMeans(x),
    mean.y = colMeans(y),
    difference = difference,
    variance1 = variance1,
    variance2 = variance2,
    D = denominator,
    A = A,
    T.BF = T.BF,
    Q3 = n1 * T.BF,
    b1 = b1,
    b2 = b2,
    null.centering = b1 + b2,
    trace1 = trace1,
    trace2 = trace2,
    trace12 = trace12,
    coefficients = coefficients,
    variance = null.variance,
    z = z
  )
}


test_that("FZWZ BF matches the literal published leave-out formulas", {
  x <- cbind(
    c(-1.2, -0.1, 0.4, 1.3, 2.1, -0.7),
    c(0.2, 1.7, -1.1, 0.8, 2.4, -0.4),
    c(2.0, -0.5, 0.7, 1.4, -1.6, 0.1)
  )
  y <- cbind(
    c(-0.8, 0.6, 1.7, -1.1, 0.2, 2.5, -0.3),
    c(1.3, -0.6, 0.4, 2.0, -1.4, 0.9, 0.1),
    c(-1.2, 0.3, 1.1, -0.4, 2.2, 0.8, -0.9)
  )
  reference <- fzwz_reference(x, y)
  result <- feng_zou_wang_zhu_two_sample_test(x, y)

  expect_s3_class(result, "hd_location_test")
  expect_s3_class(result, "htest")
  expect_equal(unname(result$statistic), reference$z, tolerance = 3e-12)
  expect_equal(result$p.value,
               stats::pnorm(reference$z, lower.tail = FALSE),
               tolerance = 1e-14)
  expect_equal(result$components$T.BF, reference$T.BF,
               tolerance = 2e-13)
  expect_equal(result$components$Q3, reference$Q3, tolerance = 2e-13)
  expect_equal(result$components$b1, reference$b1, tolerance = 2e-13)
  expect_equal(result$components$b2, reference$b2, tolerance = 2e-13)
  expect_equal(result$components$null.centering,
               reference$null.centering, tolerance = 2e-13)
  expect_equal(result$components$trace1, reference$trace1,
               tolerance = 2e-12)
  expect_equal(result$components$trace2, reference$trace2,
               tolerance = 2e-12)
  expect_equal(result$components$trace12, reference$trace12,
               tolerance = 2e-12)
  expect_equal(result$components$variance.coefficients,
               reference$coefficients, tolerance = 1e-15)
  expect_equal(result$variance[[1L]], reference$variance,
               tolerance = 2e-12)
  expect_equal(unname(result$components$variance1.diagonal),
               reference$variance1, tolerance = 2e-13)
  expect_equal(unname(result$components$variance2.diagonal),
               reference$variance2, tolerance = 2e-13)
  expect_equal(unname(result$components$D.hat.diagonal), reference$D,
               tolerance = 2e-13)
  expect_equal(unname(result$components$A.coordinate), reference$A,
               tolerance = 2e-13)
})


test_that("FZWZ BF uses one gamma in the group-2 b1 correction", {
  set.seed(2421)
  x <- matrix(stats::rnorm(21, sd = 0.7), 7, 3)
  y <- matrix(stats::rnorm(30, sd = 2.3), 10, 3)
  result <- feng_zou_wang_zhu_two_sample_test(x, y)
  n1 <- nrow(x)
  n2 <- nrow(y)
  gamma <- n1 / n2
  variance1 <- apply(x, 2L, stats::var)
  variance2 <- apply(y, 2L, stats::var)
  denominator <- variance1 + gamma * variance2
  published <- sum(
    2 * variance1^2 / (n1 * (n1 - 1) * denominator^2) +
      2 * gamma * variance2^2 /
        (n2 * (n2 - 1) * denominator^2)
  )
  extra.gamma <- sum(
    2 * variance1^2 / (n1 * (n1 - 1) * denominator^2) +
      2 * gamma^2 * variance2^2 /
        (n2 * (n2 - 1) * denominator^2)
  )

  expect_equal(result$components$b1, published, tolerance = 2e-13)
  expect_gt(abs(result$components$b1 - extra.gamma), 1e-4)
  expect_identical(result$diagnostics$b1.group2.gamma.power, 1)
})


test_that("FZWZ BF exposes the published permutation denominators", {
  set.seed(2422)
  x <- matrix(stats::rnorm(24), 6, 4)
  y <- matrix(stats::rnorm(32), 8, 4)
  result <- feng_zou_wang_zhu_two_sample_test(x, y)

  expect_equal(result$components$P4.n1, 6 * 5 * 4 * 3)
  expect_equal(result$components$P4.n2, 8 * 7 * 6 * 5)
  expect_equal(result$components$P2.n1, 6 * 5)
  expect_equal(result$components$P2.n2, 8 * 7)
  expect_equal(result$components$within1.combinations, choose(6, 4))
  expect_equal(result$components$within2.combinations, choose(8, 4))
  expect_equal(result$components$cross.pair.combinations,
               choose(6, 2) * choose(8, 2))
  expect_equal(
    result$components$variance.coefficients,
    c(group1 = 2 / (6 * 5), group2 = 2 / (8 * 7),
      cross = 4 / (6 * 8)),
    tolerance = 1e-15
  )
  expect_match(result$diagnostics$within.trace.estimator,
               "2 \\* P_n\\^4")
  expect_match(result$diagnostics$cross.trace.estimator,
               "4 \\* P_n1\\^2 \\* P_n2\\^2")
})


test_that("FZWZ BF uses the dual high-dimensional path without p by p storage", {
  set.seed(2423)
  x <- matrix(stats::rnorm(180), 6, 30)
  y <- matrix(stats::rnorm(210, 0.1), 7, 30)
  result <- feng_zou_wang_zhu_two_sample_test(x, y)

  expect_true(is.finite(unname(result$statistic)))
  expect_true(is.finite(result$p.value))
  expect_identical(result$diagnostics$trace.computation,
                   "local 4 by 4 dual Gram symmetry reduction")
  expect_identical(result$diagnostics$constructs.p.by.p.matrix, FALSE)
  expect_equal(result$p, 30L)
})


test_that("FZWZ BF handles unequal sizes and exchanging groups", {
  set.seed(2431)
  x <- matrix(stats::rnorm(28), 7, 4)
  y <- matrix(stats::rnorm(44, 0.2, 1.4), 11, 4)
  baseline <- feng_zou_wang_zhu_two_sample_test(x, y)
  swapped <- feng_zou_wang_zhu_two_sample_test(y, x)
  gamma <- nrow(x) / nrow(y)

  expect_equal(swapped$statistic, baseline$statistic, tolerance = 5e-12)
  expect_equal(swapped$p.value, baseline$p.value, tolerance = 5e-14)
  expect_equal(swapped$components$Q3, baseline$components$Q3,
               tolerance = 5e-12)
  expect_equal(swapped$components$T.BF,
               gamma * baseline$components$T.BF, tolerance = 5e-12)
  expect_equal(swapped$components$null.centering,
               gamma * baseline$components$null.centering,
               tolerance = 5e-12)
  expect_equal(swapped$variance[[1L]],
               gamma^2 * baseline$variance[[1L]], tolerance = 8e-12)
  expect_equal(swapped$components$trace1,
               gamma^2 * baseline$components$trace2,
               tolerance = 8e-12)
  expect_equal(swapped$components$trace2,
               gamma^2 * baseline$components$trace1,
               tolerance = 8e-12)
  expect_equal(swapped$components$trace12,
               gamma^2 * baseline$components$trace12,
               tolerance = 8e-12)
  expect_equal(swapped$estimate, -baseline$estimate, tolerance = 2e-14)
})


test_that("FZWZ BF is invariant to translation and diagonal units", {
  set.seed(2432)
  x <- matrix(stats::rnorm(40), 8, 5)
  y <- matrix(stats::rnorm(50, 0.15, 1.3), 10, 5)
  baseline <- feng_zou_wang_zhu_two_sample_test(x, y)

  shift <- c(30, -17, 4.5, 100, -0.25)
  translated <- feng_zou_wang_zhu_two_sample_test(
    sweep(x, 2L, shift, "+"), sweep(y, 2L, shift, "+")
  )
  expect_equal(translated$statistic, baseline$statistic, tolerance = 8e-12)
  expect_equal(translated$components$T.BF, baseline$components$T.BF,
               tolerance = 8e-12)
  expect_equal(translated$components$null.centering,
               baseline$components$null.centering, tolerance = 8e-12)
  expect_equal(translated$variance, baseline$variance, tolerance = 8e-12)
  expect_equal(translated$estimate, baseline$estimate, tolerance = 8e-14)

  units <- c(-4, 0.2, 7.5, -0.125, 3)
  rescaled <- feng_zou_wang_zhu_two_sample_test(
    sweep(x, 2L, units, "*"), sweep(y, 2L, units, "*")
  )
  expect_equal(rescaled$statistic, baseline$statistic, tolerance = 8e-12)
  expect_equal(rescaled$components$T.BF, baseline$components$T.BF,
               tolerance = 8e-12)
  expect_equal(rescaled$components$b1, baseline$components$b1,
               tolerance = 8e-12)
  expect_equal(rescaled$components$b2, baseline$components$b2,
               tolerance = 8e-12)
  expect_equal(rescaled$components$trace1, baseline$components$trace1,
               tolerance = 8e-12)
  expect_equal(rescaled$components$trace2, baseline$components$trace2,
               tolerance = 8e-12)
  expect_equal(rescaled$components$trace12, baseline$components$trace12,
               tolerance = 8e-12)

  set.seed(2433)
  reordered <- feng_zou_wang_zhu_two_sample_test(
    x[sample.int(nrow(x)), , drop = FALSE],
    y[sample.int(nrow(y)), , drop = FALSE]
  )
  expect_equal(reordered$statistic, baseline$statistic, tolerance = 8e-12)
  expect_equal(reordered$components$T.BF, baseline$components$T.BF,
               tolerance = 8e-12)
  expect_equal(reordered$components$trace1, baseline$components$trace1,
               tolerance = 8e-12)
  expect_equal(reordered$components$trace2, baseline$components$trace2,
               tolerance = 8e-12)
  expect_equal(reordered$components$trace12, baseline$components$trace12,
               tolerance = 8e-12)
})


test_that("FZWZ BF supports p equals one and the minimum sample sizes", {
  x <- matrix(c(-2.1, -0.4, 0.3, 1.2, 2.7, 3.4), ncol = 1)
  y <- matrix(c(-1.7, -0.8, 0.9, 1.8, 2.2, 4.1), ncol = 1)
  reference <- fzwz_reference(x, y)
  result <- feng_zou_wang_zhu_two_sample_test(x, y)

  expect_true(is.finite(unname(result$statistic)))
  expect_true(is.finite(result$p.value))
  expect_equal(unname(result$statistic), reference$z, tolerance = 3e-12)
  expect_equal(result$components$trace1, reference$trace1,
               tolerance = 3e-12)
  expect_equal(result$components$trace2, reference$trace2,
               tolerance = 3e-12)
  expect_equal(result$components$trace12, reference$trace12,
               tolerance = 3e-12)
  expect_equal(result$n, c(group1 = 6, group2 = 6))
  expect_equal(result$p, 1L)
})


test_that("FZWZ BF remains stable under extreme finite measurement units", {
  set.seed(2441)
  x <- matrix(stats::rnorm(32), 8, 4)
  y <- matrix(stats::rnorm(36, 0.1), 9, 4)
  baseline <- feng_zou_wang_zhu_two_sample_test(x, y)

  huge <- feng_zou_wang_zhu_two_sample_test(1e150 * x, 1e150 * y)
  tiny <- feng_zou_wang_zhu_two_sample_test(1e-150 * x, 1e-150 * y)
  expect_equal(huge$statistic, baseline$statistic, tolerance = 1e-11)
  expect_equal(tiny$statistic, baseline$statistic, tolerance = 1e-11)
  expect_equal(huge$components$T.BF, baseline$components$T.BF,
               tolerance = 1e-11)
  expect_equal(tiny$variance, baseline$variance, tolerance = 1e-11)

  units <- c(-3e150, 2e-150, -7e75, 0.125)
  mixed <- feng_zou_wang_zhu_two_sample_test(
    sweep(x, 2L, units, "*"), sweep(y, 2L, units, "*")
  )
  expect_true(is.finite(unname(mixed$statistic)))
  expect_true(is.finite(mixed$p.value))
  expect_equal(mixed$statistic, baseline$statistic, tolerance = 1e-11)
  expect_equal(mixed$components$null.centering,
               baseline$components$null.centering, tolerance = 1e-11)
})


test_that("FZWZ BF rejects malformed and degenerate inputs without repair", {
  expect_error(
    feng_zou_wang_zhu_two_sample_test(
      matrix(stats::rnorm(15), 5, 3), matrix(stats::rnorm(18), 6, 3)
    ),
    "at least 6 row"
  )
  expect_error(
    feng_zou_wang_zhu_two_sample_test(
      matrix(stats::rnorm(18), 6, 3), matrix(stats::rnorm(24), 6, 4)
    ),
    "same number of columns"
  )

  named.x <- matrix(stats::rnorm(18), 6, 3,
                    dimnames = list(NULL, c("a", "b", "c")))
  named.y <- matrix(stats::rnorm(21), 7, 3,
                    dimnames = list(NULL, c("b", "a", "c")))
  expect_error(
    feng_zou_wang_zhu_two_sample_test(named.x, named.y),
    "same names in the same order"
  )
  named.x[1, 1] <- Inf
  expect_error(
    feng_zou_wang_zhu_two_sample_test(named.x, unname(named.y)),
    "finite values"
  )

  expect_error(
    feng_zou_wang_zhu_two_sample_test(
      matrix(1, 6, 2), matrix(-3, 7, 2)
    ),
    "full-sample combined marginal variance"
  )

  # The full variance is positive, but omitting the four nonzero observations
  # from group 1 leaves two identical values while group 2 is constant.
  expect_error(
    feng_zou_wang_zhu_two_sample_test(
      matrix(c(0, 0, 1, 2, 3, 4), ncol = 1),
      matrix(2, 6, 1)
    ),
    "leave-four-out combined marginal variance"
  )
})


test_that("FZWZ BF preserves names and documents its numerical contract", {
  set.seed(2451)
  variables <- c("height", "width", "depth")
  x <- matrix(stats::rnorm(21), 7, 3,
              dimnames = list(NULL, variables))
  y <- matrix(stats::rnorm(24), 8, 3,
              dimnames = list(NULL, variables))
  result <- feng_zou_wang_zhu_two_sample_test(x, y)

  expect_identical(names(result$estimate), variables)
  expect_identical(names(result$components$D.hat.diagonal), variables)
  expect_identical(names(result$components$A.coordinate), variables)
  expect_identical(names(result$diagnostics$internal.scale.factors),
                   variables)
  expect_identical(result$diagnostics$constructs.p.by.p.matrix, FALSE)
  expect_identical(result$diagnostics$regularization, "none")
  expect_identical(result$diagnostics$variance.repair, "none")
  expect_identical(result$null.distribution$family, "normal")
  expect_identical(result$null.distribution$tail, "upper")
  expect_false(result$null.distribution$exact)
  expect_equal(result$n, c(group1 = 7, group2 = 8))
  expect_equal(result$p, 3L)
  expect_equal(result$components$Q3,
               nrow(x) * result$components$T.BF,
               tolerance = 1e-14)
  expect_equal(
    unname(result$statistic),
    (result$components$T.BF - result$components$null.centering) /
      sqrt(result$variance[[1L]]),
    tolerance = 1e-14
  )
})
