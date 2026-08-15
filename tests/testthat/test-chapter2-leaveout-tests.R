pa_reference <- function(x, mu = numeric(ncol(x))) {
  z <- sweep(x, 2L, mu, "-")
  n <- nrow(z)
  t1 <- 0
  t2 <- 0
  min_variance <- Inf
  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      if (i == j) {
        next
      }
      keep <- setdiff(seq_len(n), c(i, j))
      mean.minus <- colMeans(z[keep, , drop = FALSE])
      variance.minus <- apply(z[keep, , drop = FALSE], 2L, stats::var)
      min_variance <- min(min_variance, variance.minus)
      t1 <- t1 + sum(z[i, ] * z[j, ] / variance.minus)
      left <- sum(z[i, ] * (z[j, ] - mean.minus) / variance.minus)
      right <- sum(z[j, ] * (z[i, ] - mean.minus) / variance.minus)
      t2 <- t2 + left * right
    }
  }
  correction <- (n - 5) / (n - 3)
  u <- correction * t1 / (n * (n - 1))
  variance <- correction^2 * 2 * t2 / (n * (n - 1))^2
  list(
    T1 = t1,
    T2 = t2,
    correction = correction,
    U = u,
    variance = variance,
    z = u / sqrt(variance),
    min.variance = min_variance
  )
}


cq_reference <- function(x, y) {
  n1 <- nrow(x)
  n2 <- nrow(y)
  within.u <- function(z) {
    n <- nrow(z)
    total <- 0
    for (i in seq_len(n)) {
      for (j in seq_len(n)) {
        if (i != j) {
          total <- total + sum(z[i, ] * z[j, ])
        }
      }
    }
    total / (n * (n - 1))
  }
  cross.u <- 0
  for (i in seq_len(n1)) {
    for (j in seq_len(n2)) {
      cross.u <- cross.u + sum(x[i, ] * y[j, ])
    }
  }
  t.direct <- within.u(x) + within.u(y) - 2 * cross.u / (n1 * n2)

  within.trace <- function(z) {
    n <- nrow(z)
    total <- 0
    for (i in seq_len(n)) {
      for (j in seq_len(n)) {
        if (i == j) {
          next
        }
        mean.minus <- colMeans(z[-c(i, j), , drop = FALSE])
        total <- total +
          sum(z[i, ] * (z[j, ] - mean.minus)) *
          sum(z[j, ] * (z[i, ] - mean.minus))
      }
    }
    total / (n * (n - 1))
  }

  a12.total <- 0
  for (i in seq_len(n1)) {
    mean.x.minus <- colMeans(x[-i, , drop = FALSE])
    for (j in seq_len(n2)) {
      mean.y.minus <- colMeans(y[-j, , drop = FALSE])
      a12.total <- a12.total +
        sum(x[i, ] * (y[j, ] - mean.y.minus)) *
        sum(y[j, ] * (x[i, ] - mean.x.minus))
    }
  }
  a1 <- within.trace(x)
  a2 <- within.trace(y)
  a12 <- a12.total / (n1 * n2)
  variance <- 2 * a1 / (n1 * (n1 - 1)) +
    2 * a2 / (n2 * (n2 - 1)) + 4 * a12 / (n1 * n2)
  list(
    T = t.direct,
    T.fast = sum((colMeans(x) - colMeans(y))^2) -
      sum(diag(stats::cov(x))) / n1 - sum(diag(stats::cov(y))) / n2,
    A1 = a1,
    A2 = a2,
    A12 = a12,
    variance = variance,
    z = t.direct / sqrt(variance)
  )
}


test_that("Park-Ayyala matches a literal leave-two-out reference", {
  set.seed(2101)
  x <- matrix(stats::rnorm(24), 8, 3)
  x <- sweep(x, 2L, c(0.4, -0.2, 0.7), "+")
  mu <- c(0.1, -0.3, 0.2)
  reference <- pa_reference(x, mu)
  result <- park_ayyala_one_sample_test(x, mu)

  expect_s3_class(result, "hd_location_test")
  expect_equal(unname(result$statistic), reference$z, tolerance = 2e-12)
  expect_equal(result$components$ordered.pair.sum, reference$T1,
               tolerance = 2e-11)
  expect_equal(result$components$variance.pair.sum, reference$T2,
               tolerance = 2e-11)
  expect_equal(result$components$U.PA, reference$U, tolerance = 2e-12)
  expect_equal(result$components$finite.sample.correction,
               reference$correction, tolerance = 1e-14)
  expect_equal(result$variance[[1L]], reference$variance, tolerance = 2e-11)
  expect_equal(
    result$variance[[1L]],
    result$components$finite.sample.correction^2 * 2 *
      result$components$variance.pair.sum /
      (nrow(x) * (nrow(x) - 1))^2,
    tolerance = 1e-14
  )
  expect_equal(
    unname(result$statistic),
    reference$T1 / sqrt(2 * reference$T2),
    tolerance = 2e-12
  )
  expect_equal(result$diagnostics$minimum.scaled.leaveout.variance,
               pa_reference(
                 {
                   residuals <- sweep(x, 2L, mu, "-")
                   sweep(residuals, 2L, apply(abs(residuals), 2L, max), "/")
                 },
                 numeric(ncol(x))
               )$min.variance,
               tolerance = 2e-12)
  expect_equal(result$p.value,
               stats::pnorm(reference$z, lower.tail = FALSE),
               tolerance = 1e-14)
})


test_that("Park-Ayyala has its stated invariances and strict failures", {
  set.seed(2102)
  x <- matrix(stats::rnorm(27), 9, 3)
  mu <- c(0.2, -0.1, 0.3)
  baseline <- park_ayyala_one_sample_test(x, mu)

  scales <- c(-4, 0.25, 7)
  rescaled <- park_ayyala_one_sample_test(
    sweep(x, 2L, scales, "*"), mu * scales
  )
  expect_equal(rescaled$statistic, baseline$statistic, tolerance = 2e-11)

  shift <- c(10, -20, 5)
  translated <- park_ayyala_one_sample_test(
    sweep(x, 2L, shift, "+"), mu + shift
  )
  expect_equal(translated$statistic, baseline$statistic, tolerance = 2e-11)

  huge.shift <- rep(1e10, ncol(x))
  shifted.x <- sweep(x, 2L, huge.shift, "+")
  shifted.mu <- mu + huge.shift
  shifted.actual.residuals <- sweep(shifted.x, 2L, shifted.mu, "-")
  stable.shifted <- park_ayyala_one_sample_test(
    shifted.actual.residuals, numeric(ncol(x))
  )
  huge.translated <- park_ayyala_one_sample_test(shifted.x, shifted.mu)
  expect_equal(huge.translated$statistic, stable.shifted$statistic,
               tolerance = 2e-12)
  expect_equal(unname(huge.translated$components$difference),
               colMeans(shifted.actual.residuals), tolerance = 1e-12)
  expect_equal(
    park_ayyala_one_sample_test(x[sample(nrow(x)), ], mu)$statistic,
    baseline$statistic,
    tolerance = 2e-11
  )

  expect_equal(
    park_ayyala_one_sample_test(1e150 * x, 1e150 * mu)$statistic,
    baseline$statistic,
    tolerance = 2e-11
  )
  expect_equal(
    park_ayyala_one_sample_test(1e-150 * x, 1e-150 * mu)$statistic,
    baseline$statistic,
    tolerance = 2e-11
  )

  # Finite operands with opposite signs can have an unrepresentable raw
  # difference. The kernel must use its pre-scaled subtraction fallback.
  extreme.pattern <- c(-1, -0.8, -0.3, 0.2, 0.7, 1, -0.6, 0.5)
  extreme.x <- cbind(1.7e308 * extreme.pattern, rev(extreme.pattern))
  extreme.mu <- c(-1.7e308, 0)
  extreme <- park_ayyala_one_sample_test(extreme.x, extreme.mu)
  expect_true(is.finite(unname(extreme$statistic)))
  expect_true(is.finite(extreme$components$variance.pair.sum))

  expect_error(
    park_ayyala_one_sample_test(matrix(stats::rnorm(15), 5, 3)),
    "at least 6 row"
  )
  expect_error(
    park_ayyala_one_sample_test(cbind(seq_len(6), rep(0, 6))),
    "leave-two-out marginal variance"
  )
})


test_that("Chen-Qin matches direct U-statistic and leave-out references", {
  set.seed(2201)
  x <- matrix(stats::rnorm(21), 7, 3)
  y <- matrix(stats::rnorm(24, 0.25), 8, 3)
  reference <- cq_reference(x, y)
  result <- chen_qin_two_sample_test(x, y)

  expect_s3_class(result, "hd_location_test")
  expect_equal(reference$T, reference$T.fast, tolerance = 2e-14)
  expect_equal(result$components$T.CQ, reference$T, tolerance = 2e-13)
  expect_equal(result$components$A1, reference$A1, tolerance = 2e-12)
  expect_equal(result$components$A2, reference$A2, tolerance = 2e-12)
  expect_equal(result$components$A12, reference$A12, tolerance = 2e-12)
  expect_equal(result$variance[[1L]], reference$variance, tolerance = 2e-12)
  expect_equal(unname(result$statistic), reference$z, tolerance = 2e-12)
  expect_equal(
    unname(result$statistic),
    result$components$T.CQ.scaled /
      sqrt(result$components$denominator.variance.scaled),
    tolerance = 2e-14
  )
  expect_equal(result$p.value,
               stats::pnorm(reference$z, lower.tail = FALSE),
               tolerance = 1e-14)
})


test_that("Chen-Qin invariances and documented translation behavior hold", {
  set.seed(2202)
  x <- matrix(stats::rnorm(28), 7, 4)
  y <- matrix(stats::rnorm(32, 0.2), 8, 4)
  baseline <- chen_qin_two_sample_test(x, y)
  swapped <- chen_qin_two_sample_test(y, x)
  expect_equal(swapped$statistic, baseline$statistic, tolerance = 2e-11)
  expect_equal(swapped$components$T.CQ, baseline$components$T.CQ,
               tolerance = 2e-12)

  q <- qr.Q(qr(matrix(c(
    1, 2, 0, 1,
    -1, 1, 2, 0,
    0, 1, 1, -1,
    2, 0, 1, 1
  ), 4, 4)))
  rotated <- chen_qin_two_sample_test(x %*% q, y %*% q)
  expect_equal(rotated$statistic, baseline$statistic, tolerance = 3e-11)

  shift <- c(20, -7, 3, 11)
  translated <- chen_qin_two_sample_test(
    sweep(x, 2L, shift, "+"), sweep(y, 2L, shift, "+")
  )
  expect_equal(translated$components$T.CQ, baseline$components$T.CQ,
               tolerance = 2e-11)
  expect_false(translated$diagnostics$finite.sample.variance.translation.invariant)

  huge.shift <- rep(1e10, ncol(x))
  huge.translated <- chen_qin_two_sample_test(
    sweep(x, 2L, huge.shift, "+"), sweep(y, 2L, huge.shift, "+")
  )
  stable.shifted.reference <- cq_reference(
    sweep(sweep(x, 2L, huge.shift, "+"), 2L, huge.shift, "-"),
    sweep(sweep(y, 2L, huge.shift, "+"), 2L, huge.shift, "-")
  )
  expect_equal(huge.translated$components$T.CQ,
               stable.shifted.reference$T, tolerance = 2e-5)
  expect_equal(huge.translated$components$T.CQ,
               baseline$components$T.CQ, tolerance = 2e-5)
  expect_equal(huge.translated$components$A12,
               baseline$components$A12, tolerance = 2e-5)
  # Compare against a stable group-mean difference formed after removing the
  # shared anchor; the unequal sample sizes preclude rowwise cancellation.
  actual.shifted.difference <-
    colMeans(sweep(x, 2L, huge.shift, "+") - huge.shift) -
    colMeans(sweep(y, 2L, huge.shift, "+") - huge.shift)
  expect_equal(unname(huge.translated$estimate),
               actual.shifted.difference, tolerance = 1e-12)

  huge <- chen_qin_two_sample_test(1e150 * x, 1e150 * y)
  tiny <- chen_qin_two_sample_test(1e-150 * x, 1e-150 * y)
  expect_equal(huge$statistic, baseline$statistic, tolerance = 3e-11)
  expect_equal(tiny$statistic, baseline$statistic, tolerance = 3e-11)
  for (result in list(huge, tiny)) {
    scaled <- unlist(result$components[c(
      "T.CQ.scaled", "A1.scaled", "A2.scaled", "A12.scaled",
      "denominator.variance.scaled"
    )])
    expect_true(all(is.finite(scaled)))
    expect_equal(
      unname(result$statistic),
      result$components$T.CQ.scaled /
        sqrt(result$components$denominator.variance.scaled),
      tolerance = 2e-14
    )
  }
  expect_true(is.infinite(huge$components$denominator.variance))
  expect_identical(tiny$components$denominator.variance, 0)
})


test_that("Chen-Qin rejects malformed and degenerate inputs", {
  expect_error(
    chen_qin_two_sample_test(matrix(1:4, 2), matrix(1:6, 3)),
    "at least 3 row"
  )
  expect_error(
    chen_qin_two_sample_test(matrix(1:9, 3, 3), matrix(1:8, 4, 2)),
    "same number of columns"
  )
  x <- matrix(rep(c(1, 2), each = 3), 3, 2)
  y <- matrix(rep(c(3, 4), each = 3), 3, 2)
  expect_error(
    chen_qin_two_sample_test(x, y),
    "strictly positive original leave-out variance"
  )

  x <- matrix(stats::rnorm(18), 6, 3,
              dimnames = list(NULL, c("a", "b", "c")))
  y <- matrix(stats::rnorm(18), 6, 3,
              dimnames = list(NULL, c("b", "a", "c")))
  expect_error(chen_qin_two_sample_test(x, y), "same names in the same order")
})
