test_that("one-sample Hotelling agrees with the univariate t test", {
  x <- matrix(c(1.2, 2.1, 0.7, 1.8, 2.5, -0.1, 1.4, 0.9), ncol = 1L)
  mu <- 0.8

  observed <- hotelling_one_sample_test(x, mu = mu)
  reference <- t.test(x[, 1L], mu = mu)

  expect_equal(unname(observed$raw.statistic),
               unname(reference$statistic)^2, tolerance = 1e-12)
  expect_equal(unname(observed$statistic),
               unname(reference$statistic)^2, tolerance = 1e-12)
  expect_equal(observed$p.value, reference$p.value, tolerance = 1e-12)
  expect_equal(observed$parameter, c(df1 = 1, df2 = nrow(x) - 1))
})

test_that("one-sample Hotelling matches its direct formula and is affine invariant", {
  x <- rbind(
    c(2, 1), c(0, 1), c(1, 3), c(4, 2), c(3, -1), c(-1, 0)
  )
  colnames(x) <- c("first", "second")
  mu <- c(0.5, 0.25)

  observed <- hotelling_one_sample_test(x, mu)
  difference <- colMeans(x) - mu
  direct_t2 <- nrow(x) * drop(difference %*% solve(cov(x), difference))
  direct_f <- (nrow(x) - ncol(x)) /
    (ncol(x) * (nrow(x) - 1)) * direct_t2

  expect_equal(unname(observed$raw.statistic), direct_t2, tolerance = 1e-12)
  expect_equal(unname(observed$statistic), direct_f, tolerance = 1e-12)
  expect_equal(observed$variance, cov(x), tolerance = 1e-12)
  expect_identical(names(observed$estimate), colnames(x))
  expect_identical(dimnames(observed$variance), list(colnames(x), colnames(x)))

  transform <- matrix(c(2, 0.4, -0.7, 1.3), 2, 2, byrow = TRUE)
  transformed <- hotelling_one_sample_test(
    x %*% t(transform), drop(transform %*% mu)
  )
  expect_equal(transformed$raw.statistic, observed$raw.statistic,
               tolerance = 1e-11)
  expect_equal(transformed$p.value, observed$p.value, tolerance = 1e-11)
})

test_that("one-sample Hotelling returns the shared htest-compatible contract", {
  set.seed(2031)
  observed <- hotelling_one_sample_test(matrix(rnorm(30), 10, 3))

  expect_s3_class(observed, "hd_location_test")
  expect_s3_class(observed, "htest")
  expect_named(observed, c(
    "statistic", "parameter", "p.value", "method", "data.name",
    "alternative", "raw.statistic", "estimate", "null.value",
    "null.distribution", "variance", "components", "diagnostics", "n",
    "p", "call"
  ))
  expect_identical(names(observed$statistic), "F")
  expect_identical(observed$null.distribution$family, "F")
  expect_true(observed$null.distribution$exact)
  expect_identical(observed$diagnostics$solver, "Cholesky")
  expect_identical(names(observed$n), "x")
  expect_output(print(observed), "Hotelling")
})

test_that("one-sample Hotelling rejects invalid dimensions and covariance matrices", {
  expect_error(
    hotelling_one_sample_test(matrix(seq_len(9), 3, 3)),
    "requires `n > p`",
    fixed = TRUE
  )
  singular <- cbind(seq_len(6), 2 * seq_len(6))
  expect_error(
    hotelling_one_sample_test(singular),
    "not positive definite|ill-conditioned"
  )
  expect_error(
    hotelling_one_sample_test(matrix(rnorm(20), 10, 2), tol = 1),
    "strictly between zero and one"
  )
  expect_error(
    hotelling_one_sample_test(matrix(c(1, NA_real_, 2, 3), 2, 2)),
    "finite"
  )
})

test_that("two-sample Hotelling agrees with the pooled univariate t test", {
  x <- matrix(c(1.1, 0.8, 1.7, 2.0, 1.4, 0.9), ncol = 1L)
  y <- matrix(c(0.2, 0.7, 1.0, 0.4, 0.5, 1.2, 0.1), ncol = 1L)

  observed <- hotelling_two_sample_test(x, y)
  reference <- t.test(x[, 1L], y[, 1L], var.equal = TRUE)

  expect_equal(unname(observed$raw.statistic),
               unname(reference$statistic)^2, tolerance = 1e-12)
  expect_equal(unname(observed$statistic),
               unname(reference$statistic)^2, tolerance = 1e-12)
  expect_equal(observed$p.value, reference$p.value, tolerance = 1e-12)
})

test_that("two-sample Hotelling matches its direct pooled formula", {
  x <- rbind(c(1, 2), c(2, 0), c(0, 1), c(3, 4), c(2, 3))
  y <- rbind(c(-1, 0), c(0, -2), c(1, 1), c(-2, 2), c(2, -1), c(0, 3))
  colnames(x) <- colnames(y) <- c("a", "b")
  n1 <- nrow(x)
  n2 <- nrow(y)
  total <- n1 + n2
  difference <- colMeans(x) - colMeans(y)
  pooled <- ((n1 - 1) * cov(x) + (n2 - 1) * cov(y)) / (total - 2)
  direct_t2 <- n1 * n2 / total *
    drop(difference %*% solve(pooled, difference))
  direct_f <- (total - ncol(x) - 1) /
    (ncol(x) * (total - 2)) * direct_t2

  observed <- hotelling_two_sample_test(x, y)
  expect_equal(unname(observed$raw.statistic), direct_t2, tolerance = 1e-12)
  expect_equal(unname(observed$statistic), direct_f, tolerance = 1e-12)
  expect_equal(observed$variance, pooled, tolerance = 1e-12)
  expect_equal(observed$components$difference, difference, tolerance = 1e-12)
  expect_equal(observed$n, c(x = n1, y = n2))
})

test_that("two-sample Hotelling is affine invariant and symmetric in the samples", {
  x <- rbind(c(1, 2), c(2, 0), c(0, 1), c(3, 4), c(2, 3))
  y <- rbind(c(-1, 0), c(0, -2), c(1, 1), c(-2, 2), c(2, -1), c(0, 3))
  transform <- matrix(c(1.5, -0.2, 0.8, 2.0), 2, 2, byrow = TRUE)

  observed <- hotelling_two_sample_test(x, y)
  transformed <- hotelling_two_sample_test(x %*% t(transform),
                                            y %*% t(transform))
  swapped <- hotelling_two_sample_test(y, x)

  expect_equal(transformed$raw.statistic, observed$raw.statistic,
               tolerance = 1e-11)
  expect_equal(transformed$p.value, observed$p.value, tolerance = 1e-11)
  expect_equal(swapped$raw.statistic, observed$raw.statistic,
               tolerance = 1e-12)
  expect_equal(swapped$p.value, observed$p.value, tolerance = 1e-12)
  expect_equal(swapped$estimate, -observed$estimate, tolerance = 1e-12)
})

test_that("two-sample Hotelling validates variables, degrees of freedom, and rank", {
  expect_error(
    hotelling_two_sample_test(matrix(1:8, 2, 4), matrix(9:16, 2, 4)),
    "requires `n1 + n2 > p + 1`",
    fixed = TRUE
  )
  x <- matrix(rnorm(15), 5, 3, dimnames = list(NULL, c("a", "b", "c")))
  y <- matrix(rnorm(18), 6, 3, dimnames = list(NULL, c("a", "c", "b")))
  expect_error(
    hotelling_two_sample_test(x, y),
    "same names in the same order"
  )
  singular_x <- cbind(1:4, 2 * (1:4))
  singular_y <- cbind(5:8, 2 * (5:8) + 1)
  expect_error(
    hotelling_two_sample_test(singular_x, singular_y),
    "not positive definite|ill-conditioned"
  )
  expect_error(
    hotelling_two_sample_test(matrix(rnorm(10), 5, 2),
                              matrix(rnorm(15), 5, 3)),
    "same number of columns"
  )
})
