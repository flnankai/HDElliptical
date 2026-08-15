.reference_unit_rows <- function(x, zero_tol = 0) {
  norms <- sqrt(rowSums(x^2))
  result <- matrix(0, nrow(x), ncol(x))
  nonzero <- norms > zero_tol
  result[nonzero, ] <- x[nonzero, , drop = FALSE] / norms[nonzero]
  result
}

.reference_signed_ranks <- function(centered, zero_tol = 0) {
  n <- nrow(centered)
  result <- matrix(0, n, ncol(centered))
  for (i in seq_len(n)) {
    pair_sums <- sweep(centered, 2L, centered[i, ], "+")
    result[i, ] <- colMeans(.reference_unit_rows(pair_sums, zero_tol))
  }
  result
}

.reference_pooled_ranks <- function(x, y, zero_tol = 0) {
  pooled <- rbind(x, y)
  total <- nrow(pooled)
  result <- matrix(0, total, ncol(pooled))
  for (i in seq_len(total)) {
    differences <- sweep(pooled, 2L, pooled[i, ], function(a, b) b - a)
    result[i, ] <- colMeans(.reference_unit_rows(differences, zero_tol))
  }
  result
}


test_that("spatial sign test matches an independent R calculation", {
  x <- rbind(
    c(2.0, 1.0), c(-0.5, 1.5), c(1.2, -0.8), c(3.0, 0.2),
    c(-1.0, -2.0), c(0.4, 2.2), c(1.7, -1.3), c(-2.1, 0.6)
  )
  colnames(x) <- c("first", "second")
  mu <- c(0.25, -0.1)
  signs <- .reference_unit_rows(sweep(x, 2L, mu, "-"))
  mean_sign <- colMeans(signs)
  second_moment <- crossprod(signs) / nrow(x)
  reference <- nrow(x) * drop(
    mean_sign %*% solve(second_moment, mean_sign)
  )

  observed <- spatial_sign_test(x, mu)
  expect_equal(unname(observed$statistic), reference, tolerance = 1e-12)
  expect_equal(unname(observed$components$mean.sign), mean_sign,
               tolerance = 1e-12)
  expect_equal(unname(observed$variance), second_moment, tolerance = 1e-12)
  expect_equal(observed$p.value,
               pchisq(reference, df = ncol(x), lower.tail = FALSE),
               tolerance = 1e-14)
  expect_identical(dimnames(observed$variance), list(colnames(x), colnames(x)))
})

test_that("spatial sign test has the expected similarity invariances", {
  set.seed(2101)
  x <- matrix(rnorm(90, 0.2), 30, 3)
  mu <- c(-0.1, 0.25, 0.05)
  rotation <- qr.Q(qr(matrix(c(
    1, 2, 3, -2, 1, 0.5, 0.2, -1, 2
  ), 3, 3)))
  shift <- c(7, -4, 2)
  scale <- 3.5
  transformed_x <- scale * (x %*% rotation) +
    matrix(shift, nrow(x), 3, byrow = TRUE)
  transformed_mu <- scale * drop(mu %*% rotation) + shift

  observed <- spatial_sign_test(x, mu)
  transformed <- spatial_sign_test(transformed_x, transformed_mu)
  expect_equal(transformed$statistic, observed$statistic, tolerance = 1e-11)
  expect_equal(transformed$p.value, observed$p.value, tolerance = 1e-12)
})

test_that("spatial sign test exposes zero handling and rejects degeneracy", {
  x <- rbind(c(0, 0), c(1, 0), c(0, 1), c(-1, 1), c(2, -1), c(-2, -1))
  observed <- spatial_sign_test(x)

  expect_s3_class(observed, "hd_location_test")
  expect_s3_class(observed, "htest")
  expect_identical(observed$null.distribution$family, "chi-squared")
  expect_false(observed$null.distribution$exact)
  expect_equal(observed$diagnostics$n.zero.signs, 1)
  expect_output(print(observed), "spatial-sign")

  expect_error(
    spatial_sign_test(matrix(seq_len(9), 3, 3)),
    "requires `n > p`",
    fixed = TRUE
  )
  expect_error(
    spatial_sign_test(matrix(rep(c(1, 0), 6), 6, 2, byrow = TRUE)),
    "not positive definite|ill-conditioned"
  )
  expect_error(spatial_sign_test(x, zero_tol = -1), "non-negative")
})


test_that("signed-rank test matches the Hoeffding-scaled R reference", {
  x <- rbind(
    c(2.2, 1.0), c(-0.4, 1.7), c(1.3, -0.9), c(3.1, 0.3),
    c(-1.2, -2.1), c(0.6, 2.4), c(1.8, -1.1), c(-2.0, 0.8),
    c(0.9, 0.4)
  )
  mu <- c(0.2, -0.15)
  ranks <- .reference_signed_ranks(sweep(x, 2L, mu, "-"))
  mean_rank <- colMeans(ranks)
  second_moment <- crossprod(ranks) / nrow(x)
  unscaled_quadratic <- nrow(x) * drop(
    mean_rank %*% solve(second_moment, mean_rank)
  )
  reference <- unscaled_quadratic / 4

  observed <- spatial_signed_rank_test(x, mu)
  expect_equal(unname(observed$statistic), reference, tolerance = 1e-11)
  expect_equal(unname(observed$components$mean.signed.rank), mean_rank,
               tolerance = 1e-12)
  expect_equal(unname(observed$components$rank.second.moment), second_moment,
               tolerance = 1e-12)
  expect_equal(unname(observed$variance), 4 * second_moment,
               tolerance = 1e-12)
  expect_equal(observed$components$quadratic.multiplier, nrow(x) / 4)
  expect_equal(4 * unname(observed$statistic), unscaled_quadratic,
               tolerance = 1e-11)
  expect_identical(observed$diagnostics$hoeffding.factor, 4)
})

test_that("signed-rank test is translation, scale, and orthogonally invariant", {
  set.seed(2102)
  x <- matrix(rnorm(100, 0.35), 25, 4)
  mu <- c(0.1, -0.2, 0.15, 0.3)
  rotation <- qr.Q(qr(matrix(rnorm(16), 4, 4)))
  shift <- c(10, -3, 5, 1)
  scale <- 0.4
  transformed_x <- scale * (x %*% rotation) +
    matrix(shift, nrow(x), 4, byrow = TRUE)
  transformed_mu <- scale * drop(mu %*% rotation) + shift

  observed <- spatial_signed_rank_test(x, mu)
  transformed <- spatial_signed_rank_test(transformed_x, transformed_mu)
  expect_equal(transformed$statistic, observed$statistic, tolerance = 1e-10)
  expect_equal(transformed$p.value, observed$p.value, tolerance = 1e-11)
})

test_that("signed-rank test counts zero pair sums and rejects singular moments", {
  x <- rbind(c(1, 0), c(-1, 0), c(0, 1), c(0, -2), c(2, 1), c(-1, 2))
  observed <- spatial_signed_rank_test(x)
  expect_gte(observed$diagnostics$n.zero.ordered.pairs, 2)
  expect_equal(observed$diagnostics$n.ordered.pairs, nrow(x)^2)
  expect_output(print(observed), "signed-rank")

  collinear <- cbind(seq_len(7), 2 * seq_len(7))
  expect_error(
    spatial_signed_rank_test(collinear),
    "not positive definite|ill-conditioned"
  )
  expect_error(
    spatial_signed_rank_test(matrix(rnorm(16), 4, 4)),
    "requires `n > p`",
    fixed = TRUE
  )
  expect_error(spatial_signed_rank_test(x, tol = 0), "strictly between")
})


test_that("pooled spatial-rank test matches an independent R calculation", {
  x <- rbind(c(1, 2), c(2, 0), c(0, 1), c(3, 4), c(2, 3), c(-1, 2))
  y <- rbind(
    c(-1, 0), c(0, -2), c(1, 1), c(-2, 2), c(2, -1), c(0, 3), c(-3, -1)
  )
  ranks <- .reference_pooled_ranks(x, y)
  n1 <- nrow(x)
  n2 <- nrow(y)
  total <- n1 + n2
  mean_x <- colMeans(ranks[seq_len(n1), , drop = FALSE])
  mean_y <- colMeans(ranks[n1 + seq_len(n2), , drop = FALSE])
  difference <- mean_x - mean_y
  rank_scatter <- crossprod(ranks)
  covariance <- rank_scatter / (total - 1)
  old_second_moment <- rank_scatter / total
  reference <- n1 * n2 / total * drop(
    difference %*% solve(covariance, difference)
  )
  old_reference <- n1 * n2 / total * drop(
    difference %*% solve(old_second_moment, difference)
  )

  observed <- spatial_rank_test(x, y)
  expect_equal(unname(observed$statistic), reference, tolerance = 1e-11)
  expect_equal(unname(observed$components$mean.rank.x), mean_x,
               tolerance = 1e-12)
  expect_equal(unname(observed$components$mean.rank.y), mean_y,
               tolerance = 1e-12)
  expect_equal(unname(observed$components$difference), difference,
               tolerance = 1e-12)
  expect_equal(unname(observed$variance), covariance, tolerance = 1e-12)
  expect_equal(observed$diagnostics$n.zero.ordered.pairs, total)
  expect_equal(observed$diagnostics$covariance.denominator, total - 1)
  expect_equal(reference, (total - 1) / total * old_reference,
               tolerance = 1e-12)
  expect_equal(unname(observed$statistic),
               (total - 1) / total * old_reference,
               tolerance = 1e-11)
})

test_that("pooled spatial-rank test is similarity invariant and group symmetric", {
  set.seed(2103)
  x <- matrix(rnorm(72, 0.4), 24, 3)
  y <- matrix(rnorm(81, -0.1), 27, 3)
  rotation <- qr.Q(qr(matrix(rnorm(9), 3, 3)))
  shift <- c(8, -5, 2)
  transformed_x <- 2.75 * (x %*% rotation) +
    matrix(shift, nrow(x), 3, byrow = TRUE)
  transformed_y <- 2.75 * (y %*% rotation) +
    matrix(shift, nrow(y), 3, byrow = TRUE)

  observed <- spatial_rank_test(x, y)
  transformed <- spatial_rank_test(transformed_x, transformed_y)
  swapped <- spatial_rank_test(y, x)
  expect_equal(transformed$statistic, observed$statistic, tolerance = 1e-10)
  expect_equal(transformed$p.value, observed$p.value, tolerance = 1e-11)
  expect_equal(swapped$statistic, observed$statistic, tolerance = 1e-11)
  expect_equal(swapped$estimate, -observed$estimate, tolerance = 1e-12)
})

test_that("pooled spatial-rank test validates ties, variables, and rank", {
  x <- rbind(c(0, 0), c(0, 0), c(1, 0), c(0, 1))
  y <- rbind(c(-1, 0), c(0, -1), c(2, 1), c(-1, 2))
  observed <- spatial_rank_test(x, y)
  expect_gt(observed$diagnostics$n.zero.ordered.pairs, nrow(x) + nrow(y))
  expect_output(print(observed), "pooled spatial-rank")

  collinear_x <- cbind(1:5, 2 * (1:5))
  collinear_y <- cbind(6:10, 2 * (6:10))
  expect_error(
    spatial_rank_test(collinear_x, collinear_y),
    "not positive definite|ill-conditioned"
  )
  named_x <- matrix(rnorm(18), 6, 3,
                    dimnames = list(NULL, c("a", "b", "c")))
  named_y <- matrix(rnorm(21), 7, 3,
                    dimnames = list(NULL, c("a", "c", "b")))
  expect_error(spatial_rank_test(named_x, named_y), "same names")
  expect_error(
    spatial_rank_test(matrix(1:8, 2, 4), matrix(9:16, 2, 4)),
    "requires `n1 + n2 > p`",
    fixed = TRUE
  )
  expect_error(spatial_rank_test(x, y, zero_tol = Inf), "finite")
})

test_that("spatial direction kernels stabilize huge and subnormal values", {
  huge <- 1.7e308
  tiny <- .Machine$double.xmin * .Machine$double.eps
  expect_gt(tiny, 0)

  one_sample <- matrix(c(huge, -huge, tiny, -tiny, 1, -2), ncol = 1L)
  sign_result <- spatial_sign_test(one_sample)
  signed_rank_result <- spatial_signed_rank_test(one_sample)

  expect_true(is.finite(sign_result$statistic))
  expect_equal(sign_result$diagnostics$n.zero.signs, 0)
  expect_true(is.finite(signed_rank_result$statistic))
  # Only the two opposite-value pairs have zero sums: two ordered pairs each.
  expect_equal(signed_rank_result$diagnostics$n.zero.ordered.pairs, 4)

  x <- matrix(c(huge, tiny, 1), ncol = 1L)
  y <- matrix(c(-huge, -tiny, -2), ncol = 1L)
  pooled_result <- spatial_rank_test(x, y)
  expect_true(is.finite(pooled_result$statistic))
  # All observations are distinct, so only the six self-differences are zero.
  expect_equal(pooled_result$diagnostics$n.zero.ordered.pairs, 6)
})
