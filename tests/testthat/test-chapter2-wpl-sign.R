wpl_sign_matrix_reference <- function(x, mu) {
  residual <- sweep(x, 2L, mu, "-")
  norms <- sqrt(rowSums(residual^2))
  z <- matrix(0, nrow(x), ncol(x))
  nonzero <- norms > 0
  z[nonzero, ] <- residual[nonzero, , drop = FALSE] / norms[nonzero]
  z
}


wpl_reference <- function(x, mu = numeric(ncol(x))) {
  z <- wpl_sign_matrix_reference(x, mu)
  n <- nrow(z)
  T.raw <- 0
  for (j in seq_len(n - 1L)) {
    for (k in seq.int(j + 1L, n)) {
      T.raw <- T.raw + sum(z[j, ] * z[k, ])
    }
  }

  cross.validation.numerator <- 0
  for (j in seq_len(n)) {
    for (k in seq_len(n)) {
      if (j == k) {
        next
      }
      leaveout.mean <- colMeans(z[-c(j, k), , drop = FALSE])
      cross.validation.numerator <- cross.validation.numerator +
        sum(z[j, ] * (z[k, ] - leaveout.mean)) *
        sum(z[k, ] * (z[j, ] - leaveout.mean))
    }
  }
  trace.B2.hat <- cross.validation.numerator / (n * (n - 1))
  variance.hat <- n * (n - 1) * trace.B2.hat / 2

  list(
    z = z,
    T.raw = T.raw,
    cross.validation.numerator = cross.validation.numerator,
    trace.B2.hat = trace.B2.hat,
    variance.hat = variance.hat,
    statistic = T.raw / sqrt(variance.hat)
  )
}


wpl_equation8_reference <- function(z) {
  n <- nrow(z)
  sign.crossproduct <- crossprod(z)
  z.star <- colSums(z) / (n - 2)
  -n / (n - 2)^2 +
    (n - 1) / (n * (n - 2)^2) * sum(sign.crossproduct^2) +
    (1 - 2 * n) / (n * (n - 1)) *
      drop(crossprod(z.star, sign.crossproduct %*% z.star)) +
    2 / n * sum(z.star^2) +
    (n - 2)^2 / (n * (n - 1)) * sum(z.star^2)^2
}


test_that("WPL matches the literal equation (7) reference", {
  set.seed(2601)
  x <- matrix(stats::rt(63, df = 5), 9, 7)
  mu <- c(-0.2, 0.1, 0, 0.3, -0.1, 0.25, -0.35)
  reference <- wpl_reference(x, mu)
  result <- wang_peng_li_one_sample_test(x, mu)

  expect_s3_class(result, "hd_location_test")
  expect_s3_class(result, "htest")
  expect_equal(unname(result$statistic), reference$statistic,
               tolerance = 3e-14)
  expect_equal(unname(result$raw.statistic), reference$T.raw,
               tolerance = 3e-14)
  expect_equal(result$components$T.WPL, reference$T.raw,
               tolerance = 3e-14)
  expect_equal(result$components$cross.validation.numerator,
               reference$cross.validation.numerator, tolerance = 5e-14)
  expect_equal(result$components$trace.B2.hat, reference$trace.B2.hat,
               tolerance = 5e-14)
  expect_equal(result$components$variance.hat, reference$variance.hat,
               tolerance = 5e-14)
  expect_equal(unname(result$variance), reference$variance.hat,
               tolerance = 5e-14)
  expect_equal(result$components$standard.error,
               sqrt(reference$variance.hat), tolerance = 5e-14)
  expect_equal(result$p.value,
               stats::pnorm(reference$statistic, lower.tail = FALSE),
               tolerance = 1e-15)
  expect_equal(unname(result$components$direction.sum),
               colSums(reference$z), tolerance = 3e-14)
  expect_equal(unname(result$components$direction.mean),
               colMeans(reference$z), tolerance = 3e-14)
  expect_equal(result$components$sum.direction.norm.squared, nrow(x),
               tolerance = 3e-14)
  expect_equal(result$components$pair.count, choose(nrow(x), 2))
  expect_equal(result$components$ordered.pair.count,
               nrow(x) * (nrow(x) - 1))
  expect_equal(result$null.distribution$tail, "upper")
  expect_false(result$null.distribution$exact)
  expect_identical(result$alternative, "two.sided")
  expect_identical(result$diagnostics$variance.repair, "none")
  expect_false(result$diagnostics$equation.8.shortcut.used)
})


test_that("literal estimator agrees with equation (8) when every sign is unit", {
  set.seed(2602)
  x <- matrix(stats::rnorm(96), 12, 8)
  result <- wang_peng_li_one_sample_test(x)
  signs <- wpl_sign_matrix_reference(x, numeric(ncol(x)))

  expect_equal(result$components$trace.B2.hat,
               wpl_equation8_reference(signs), tolerance = 8e-14)
  expect_equal(rowSums(signs^2), rep(1, nrow(x)), tolerance = 2e-15)
})


test_that("WPL is permutation, orthogonal, translation, and common-scale invariant", {
  set.seed(2603)
  x <- matrix(stats::rnorm(108), 12, 9)
  mu <- stats::rnorm(9, sd = 0.2)
  baseline <- wang_peng_li_one_sample_test(x, mu)

  permuted <- wang_peng_li_one_sample_test(x[sample.int(nrow(x)), ], mu)
  expect_equal(permuted$statistic, baseline$statistic, tolerance = 5e-14)
  expect_equal(permuted$raw.statistic, baseline$raw.statistic,
               tolerance = 5e-14)
  expect_equal(permuted$components$trace.B2.hat,
               baseline$components$trace.B2.hat, tolerance = 8e-14)

  q <- qr.Q(qr(matrix(stats::rnorm(81), 9, 9)))
  rotated <- wang_peng_li_one_sample_test(
    x %*% q, drop(mu %*% q)
  )
  expect_equal(rotated$statistic, baseline$statistic, tolerance = 2e-13)
  expect_equal(rotated$raw.statistic, baseline$raw.statistic,
               tolerance = 2e-13)
  expect_equal(rotated$components$trace.B2.hat,
               baseline$components$trace.B2.hat, tolerance = 3e-13)

  shift <- seq(-1.5, 1.5, length.out = ncol(x))
  translated <- wang_peng_li_one_sample_test(
    sweep(x, 2L, shift, "+"), mu + shift
  )
  expect_equal(translated$statistic, baseline$statistic,
               tolerance = 8e-14)
  expect_equal(translated$raw.statistic, baseline$raw.statistic,
               tolerance = 8e-14)
  expect_equal(translated$components$trace.B2.hat,
               baseline$components$trace.B2.hat, tolerance = 1e-13)

  rescaled <- wang_peng_li_one_sample_test(3.75 * x, 3.75 * mu)
  expect_equal(rescaled$statistic, baseline$statistic, tolerance = 5e-14)
  expect_equal(rescaled$raw.statistic, baseline$raw.statistic,
               tolerance = 5e-14)
  expect_equal(rescaled$components$trace.B2.hat,
               baseline$components$trace.B2.hat, tolerance = 8e-14)
})


test_that("WPL is generally not invariant to unequal coordinate scaling", {
  x <- rbind(
    c(1.2, -0.4, 0.3, 0.8),
    c(-0.7, 1.1, 0.5, -0.2),
    c(0.4, 0.6, -1.3, 0.7),
    c(1.5, 0.2, 0.9, -1.0),
    c(-1.1, -0.8, 0.4, 1.2),
    c(0.3, -1.4, 1.1, 0.5),
    c(0.9, 1.3, -0.6, -0.7),
    c(-0.5, 0.7, 1.4, 0.1),
    c(0.6, -0.9, -0.8, 1.5),
    c(-1.3, 0.5, 0.7, -0.4)
  )
  mu <- c(0.1, -0.2, 0.05, 0.15)
  baseline <- wang_peng_li_one_sample_test(x, mu)
  coordinate.scale <- c(0.08, 0.7, 4, 13)
  changed.units <- wang_peng_li_one_sample_test(
    sweep(x, 2L, coordinate.scale, "*"), mu * coordinate.scale
  )

  expect_gt(abs(unname(changed.units$statistic - baseline$statistic)), 1e-3)
  expect_gt(abs(changed.units$components$T.WPL -
                  baseline$components$T.WPL), 1e-3)
  expect_gt(abs(changed.units$components$trace.B2.hat -
                  baseline$components$trace.B2.hat), 1e-4)
})


test_that("U(0) is retained literally and reported", {
  set.seed(2604)
  mu <- c(-0.3, 0.1, 0.25, -0.2, 0.4)
  x <- matrix(stats::rnorm(55), 11, 5)
  x[4, ] <- mu
  reference <- wpl_reference(x, mu)
  result <- wang_peng_li_one_sample_test(x, mu)

  expect_equal(unname(result$statistic), reference$statistic,
               tolerance = 5e-14)
  expect_equal(result$components$trace.B2.hat, reference$trace.B2.hat,
               tolerance = 8e-14)
  expect_equal(result$diagnostics$zero.signs, 1)
  expect_equal(result$diagnostics$nonzero.signs, nrow(x) - 1)
  expect_equal(result$diagnostics$zero.sign.proportion, 1 / nrow(x))
  expect_identical(result$diagnostics$sign.at.zero, "U(0) = 0")
  expect_equal(result$components$sum.direction.norm.squared,
               nrow(x) - 1, tolerance = 3e-14)

  # Equation (8) assumes every sign has squared norm one and is not the same
  # estimator for this documented U(0)=0 sample.
  expect_gt(abs(result$components$trace.B2.hat -
                  wpl_equation8_reference(reference$z)), 1e-5)
})


test_that("WPL remains finite at extreme common scales", {
  set.seed(2605)
  x <- matrix(stats::runif(96, -0.8, 0.8), 12, 8)
  mu <- seq(-0.18, 0.17, length.out = 8)
  baseline <- wang_peng_li_one_sample_test(x, mu)

  for (scale in c(1e-300, 1e300)) {
    extreme <- wang_peng_li_one_sample_test(x * scale, mu * scale)
    expect_true(is.finite(extreme$statistic))
    expect_true(is.finite(extreme$p.value))
    expect_true(is.finite(extreme$components$trace.B2.hat))
    expect_equal(extreme$statistic, baseline$statistic, tolerance = 3e-13)
    expect_equal(extreme$raw.statistic, baseline$raw.statistic,
                 tolerance = 3e-13)
    expect_equal(extreme$components$trace.B2.hat,
                 baseline$components$trace.B2.hat, tolerance = 4e-13)
  }
})


test_that("overflowing direct subtraction uses an equivalent safe direction", {
  set.seed(2606)
  n <- 12
  safe.mu <- c(-1.05, 0.35, -0.2, 0.55, -0.4, 0.1)
  safe.x <- cbind(
    seq(0.9, 1.4, length.out = n),
    matrix(stats::runif(n * 5, -1.25, 1.25), n, 5)
  )
  baseline <- wang_peng_li_one_sample_test(safe.x, safe.mu)
  huge <- wang_peng_li_one_sample_test(safe.x * 1e308, safe.mu * 1e308)

  expect_equal(huge$diagnostics$subtraction.overflow.fallback.rows, n)
  expect_true(is.finite(huge$statistic))
  expect_equal(huge$statistic, baseline$statistic, tolerance = 8e-13)
  expect_equal(huge$raw.statistic, baseline$raw.statistic,
               tolerance = 8e-13)
  expect_equal(huge$components$trace.B2.hat,
               baseline$components$trace.B2.hat, tolerance = 1e-12)
})


test_that("WPL rejects invalid and degenerate calibrations without repair", {
  expect_error(
    wang_peng_li_one_sample_test(matrix(1:8, 2, 4)),
    "at least 3 row"
  )
  expect_error(
    wang_peng_li_one_sample_test(matrix(1:12, 3, 4), mu = 1:3),
    "length 4"
  )
  bad <- matrix(stats::rnorm(20), 5, 4)
  bad[2, 3] <- NA_real_
  expect_error(wang_peng_li_one_sample_test(bad), "finite")
  expect_error(
    wang_peng_li_one_sample_test(matrix(0, 8, 5)),
    "strictly positive.*8 zero spatial sign"
  )
  same.direction <- matrix(rep(c(1, 2, 3, 4), each = 9), 9, 4)
  expect_error(
    wang_peng_li_one_sample_test(same.direction),
    "strictly positive"
  )
})


test_that("WPL preserves variable names and raw reconstruction identities", {
  set.seed(2607)
  variables <- c("gene_a", "gene_b", "gene_c", "gene_d", "gene_e")
  x <- matrix(stats::rnorm(60), 12, 5,
              dimnames = list(NULL, variables))
  mu <- stats::setNames(c(0.1, -0.2, 0, 0.15, -0.05), variables)
  result <- wang_peng_li_one_sample_test(as.data.frame(x), mu)

  expect_identical(names(result$estimate), variables)
  expect_identical(names(result$null.value), variables)
  expect_identical(names(result$components$direction.sum), variables)
  expect_identical(names(result$components$direction.mean), variables)
  expect_equal(result$components$direction.mean,
               result$components$direction.sum / nrow(x), tolerance = 1e-15)
  expect_equal(result$components$trace.B2.hat,
               result$components$cross.validation.numerator /
                 result$components$ordered.pair.count,
               tolerance = 1e-15)
  expect_equal(result$components$variance.hat,
               result$components$pair.count *
                 result$components$trace.B2.hat,
               tolerance = 1e-15)
  expect_equal(unname(result$statistic),
               result$components$T.WPL / result$components$standard.error,
               tolerance = 1e-15)
})
