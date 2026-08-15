sd_reference <- function(x, mu = numeric(ncol(x))) {
  N <- nrow(x)
  p <- ncol(x)
  nu <- N - 1
  sample.mean <- colMeans(x)
  difference <- sample.mean - mu
  S <- stats::cov(x)
  variance.diagonal <- diag(S)
  inverse.sqrt.D <- diag(1 / sqrt(variance.diagonal), nrow = p)
  R <- inverse.sqrt.D %*% S %*% inverse.sqrt.D
  trace.R2 <- sum(R^2)
  A <- N * sum(difference^2 / variance.diagonal)
  correction <- 1 + trace.R2 / p^(3 / 2)
  null.center <- nu * p / (nu - 2)
  variance <- 2 * (trace.R2 - p^2 / nu) * correction
  z <- (A - null.center) / sqrt(variance)

  list(
    z = z,
    A = A,
    sample.mean = sample.mean,
    difference = difference,
    variance.diagonal = variance.diagonal,
    N = N,
    nu = nu,
    p = p,
    trace.R2 = trace.R2,
    c = correction,
    null.center = null.center,
    variance = variance
  )
}

bs_reference <- function(x, y) {
  n1 <- nrow(x)
  n2 <- nrow(y)
  N <- n1 + n2
  nu <- N - 2
  p <- ncol(x)
  mean.x <- colMeans(x)
  mean.y <- colMeans(y)
  difference <- mean.x - mean.y
  Sp <- ((n1 - 1) * stats::cov(x) + (n2 - 1) * stats::cov(y)) / nu
  trace.Sp <- sum(diag(Sp))
  trace.Sp2 <- sum(Sp^2)
  a <- N / (n1 * n2)
  xi <- sum(difference^2) - a * trace.Sp
  B <- trace.Sp2 - trace.Sp^2 / nu
  trace.Sigma2.hat <- nu^2 / (N * (N - 3)) * B
  variance <- 2 * a^2 * (N - 1) / nu * trace.Sigma2.hat
  z <- xi / sqrt(variance)

  list(
    z = z,
    xi = xi,
    mean.x = mean.x,
    mean.y = mean.y,
    difference = difference,
    n1 = n1,
    n2 = n2,
    N = N,
    nu = nu,
    p = p,
    a = a,
    trace.Sp = trace.Sp,
    trace.Sp2 = trace.Sp2,
    B = B,
    trace.Sigma2.hat = trace.Sigma2.hat,
    variance = variance
  )
}

test_that("Srivastava-Du agrees with a direct covariance-matrix calculation", {
  x <- matrix(
    c(
      2, 1, 4, 3,
      1, 3, 2, 5,
      4, 2, 1, 2,
      3, 5, 3, 1,
      5, 4, 6, 4,
      6, 2, 5, 7,
      2, 6, 4, 6
    ),
    ncol = 4,
    byrow = TRUE,
    dimnames = list(NULL, paste0("v", 1:4))
  )
  mu <- c(2, 2.5, 3, 3.5)
  ref <- sd_reference(x, mu)
  fit <- srivastava_du_one_sample_test(x, mu)

  expect_s3_class(fit, "hd_location_test")
  expect_s3_class(fit, "htest")
  expect_equal(unname(fit$statistic), ref$z, tolerance = 1e-12)
  expect_equal(unname(fit$raw.statistic), ref$A, tolerance = 1e-12)
  expect_equal(fit$p.value, stats::pnorm(ref$z, lower.tail = FALSE))
  expect_equal(unname(fit$estimate), unname(ref$sample.mean))
  expect_equal(unname(fit$null.value), mu)
  expect_equal(unname(fit$variance), ref$variance, tolerance = 1e-12)
  expect_equal(fit$components$difference, setNames(ref$difference, colnames(x)))
  expect_equal(
    fit$components$variance.diagonal,
    setNames(ref$variance.diagonal, colnames(x)),
    tolerance = 1e-12
  )
  expect_equal(fit$components$trace.R2, ref$trace.R2, tolerance = 1e-12)
  expect_equal(fit$components$c, ref$c, tolerance = 1e-12)
  expect_equal(fit$components$null.center, ref$null.center)
  expect_equal(fit$null.distribution$family, "normal")
  expect_equal(fit$null.distribution$tail, "upper")
  expect_identical(fit$alternative, "two.sided")
  expect_identical(fit$diagnostics$gram.type, "primal")
  expect_equal(fit$diagnostics$gram.dimension, ncol(x))
  expect_identical(fit$diagnostics$constructs.p.by.p.matrix, TRUE)
})

test_that("Srivastava-Du is invariant to coordinate units and signed permutation", {
  set.seed(2041)
  x <- matrix(stats::rnorm(11 * 7), nrow = 11)
  mu <- seq(-0.3, 0.3, length.out = 7)
  baseline <- srivastava_du_one_sample_test(x, mu)

  scale <- c(-3, 0.25, 2, -0.75, 4, 1.5, -2.5)
  shift <- seq(2, 8, length.out = 7)
  transformed <- sweep(sweep(x, 2, scale, `*`), 2, shift, `+`)
  transformed.mu <- mu * scale + shift
  changed.units <- srivastava_du_one_sample_test(transformed, transformed.mu)
  expect_equal(changed.units$statistic, baseline$statistic, tolerance = 2e-12)
  expect_equal(changed.units$p.value, baseline$p.value, tolerance = 2e-12)
  expect_equal(
    unname(changed.units$estimate),
    unname(baseline$estimate) * scale + shift,
    tolerance = 2e-12
  )
  expect_equal(
    unname(changed.units$components$difference),
    unname(baseline$components$difference) * scale,
    tolerance = 2e-12
  )
  expect_equal(
    unname(changed.units$components$variance.diagonal),
    unname(baseline$components$variance.diagonal) * scale^2,
    tolerance = 2e-12
  )

  permutation <- c(5, 1, 7, 3, 2, 6, 4)
  permuted <- srivastava_du_one_sample_test(
    x[, permutation, drop = FALSE],
    mu[permutation]
  )
  expect_equal(permuted$statistic, baseline$statistic, tolerance = 2e-12)
  expect_equal(permuted$p.value, baseline$p.value, tolerance = 2e-12)
})

test_that("Srivastava-Du uses its sample-space Gram representation when p is large", {
  set.seed(2042)
  x <- matrix(stats::rnorm(8 * 1500), nrow = 8)
  fit <- srivastava_du_one_sample_test(x)

  expect_true(is.finite(unname(fit$statistic)))
  expect_identical(fit$diagnostics$trace.computation, "dual Gram matrix")
  expect_identical(fit$diagnostics$gram.type, "dual")
  expect_identical(fit$diagnostics$constructs.p.by.p.matrix, FALSE)
  expect_equal(fit$diagnostics$gram.dimension, nrow(x))
  expect_false(any(vapply(fit$components, is.matrix, logical(1))))
})

test_that("Srivastava-Du remains finite under extreme column scales", {
  set.seed(2044)
  x <- matrix(stats::rnorm(12 * 6), nrow = 12)
  mu <- seq(-0.25, 0.25, length.out = 6)
  baseline <- srivastava_du_one_sample_test(x, mu)
  factors <- c(1e-150, -1e150, 3e-145, -2e145, 5e-130, -4e130)
  extreme.x <- sweep(x, 2, factors, `*`)
  extreme.mu <- mu * factors
  extreme <- srivastava_du_one_sample_test(extreme.x, extreme.mu)

  expect_true(all(is.finite(c(extreme$statistic, extreme$p.value))))
  expect_equal(extreme$statistic, baseline$statistic, tolerance = 2e-12)
  expect_equal(extreme$p.value, baseline$p.value, tolerance = 2e-12)
  expect_equal(
    unname(extreme$estimate) / factors,
    unname(baseline$estimate),
    tolerance = 2e-12
  )
  expect_equal(
    unname(extreme$components$difference) / factors,
    unname(baseline$components$difference),
    tolerance = 2e-12
  )
  expect_equal(
    unname(extreme$components$variance.diagonal) / factors^2,
    unname(baseline$components$variance.diagonal),
    tolerance = 2e-12
  )
  expect_equal(extreme$raw.statistic, baseline$raw.statistic, tolerance = 2e-12)
})

test_that("Srivastava-Du rejects undersized or degenerate inputs without repair", {
  expect_error(
    srivastava_du_one_sample_test(matrix(1:9, nrow = 3)),
    "at least 4 row"
  )

  constant.column <- cbind(1:6, rep(2, 6), c(1, 3, 2, 6, 5, 4))
  expect_error(
    srivastava_du_one_sample_test(constant.column),
    "strictly positive"
  )

  nonfinite <- matrix(stats::rnorm(24), nrow = 6)
  nonfinite[2, 3] <- Inf
  expect_error(srivastava_du_one_sample_test(nonfinite), "finite")
  expect_error(
    srivastava_du_one_sample_test(matrix(stats::rnorm(30), nrow = 6), 1:4),
    "length 5"
  )
})

test_that("Bai-Saranadasa agrees with its original pooled-covariance formula", {
  x <- matrix(
    c(
      2, 1, 4,
      1, 3, 2,
      4, 2, 1,
      3, 5, 3,
      5, 4, 6,
      6, 2, 5
    ),
    ncol = 3,
    byrow = TRUE,
    dimnames = list(NULL, c("a", "b", "c"))
  )
  y <- matrix(
    c(
      0, 2, 3,
      2, 1, 1,
      1, 4, 2,
      3, 3, 5,
      2, 5, 4
    ),
    ncol = 3,
    byrow = TRUE,
    dimnames = list(NULL, c("a", "b", "c"))
  )
  ref <- bs_reference(x, y)
  fit <- bai_saranadasa_two_sample_test(x, y)

  expect_s3_class(fit, "hd_location_test")
  expect_s3_class(fit, "htest")
  expect_equal(unname(fit$statistic), ref$z, tolerance = 1e-12)
  expect_equal(unname(fit$raw.statistic), ref$xi, tolerance = 1e-12)
  expect_equal(fit$p.value, stats::pnorm(ref$z, lower.tail = FALSE))
  expect_equal(
    unname(fit$estimate),
    unname(ref$difference),
    tolerance = 1e-12
  )
  expect_equal(fit$components$trace.Sp, ref$trace.Sp, tolerance = 1e-12)
  expect_equal(fit$components$trace.Sp2, ref$trace.Sp2, tolerance = 1e-12)
  expect_equal(fit$components$B, ref$B, tolerance = 1e-12)
  expect_equal(
    fit$components$trace.Sigma2.hat,
    ref$trace.Sigma2.hat,
    tolerance = 1e-12
  )
  expect_equal(unname(fit$variance), ref$variance, tolerance = 1e-12)
  expect_equal(fit$diagnostics$covariance.model, "common")
  expect_equal(fit$diagnostics$pooled.covariance.denominator, ref$nu)
  expect_identical(fit$diagnostics$gram.type, "primal")
  expect_equal(fit$diagnostics$gram.dimension, ncol(x))
  expect_identical(fit$diagnostics$constructs.p.by.p.matrix, TRUE)
  expect_equal(fit$null.distribution$tail, "upper")
})

test_that("Bai-Saranadasa respects orthogonal, global-scale, and group symmetries", {
  set.seed(2043)
  x <- matrix(stats::rnorm(10 * 6), nrow = 10)
  y <- matrix(stats::rnorm(13 * 6, mean = 0.2), nrow = 13)
  baseline <- bai_saranadasa_two_sample_test(x, y)

  q <- qr.Q(qr(matrix(stats::rnorm(36), nrow = 6)))
  shift <- seq(-2, 3, length.out = 6)
  x.transformed <- sweep(2.75 * x %*% q, 2, shift, `+`)
  y.transformed <- sweep(2.75 * y %*% q, 2, shift, `+`)
  transformed <- bai_saranadasa_two_sample_test(x.transformed, y.transformed)
  expect_equal(transformed$statistic, baseline$statistic, tolerance = 3e-12)
  expect_equal(transformed$p.value, baseline$p.value, tolerance = 3e-12)
  expect_equal(
    unname(transformed$raw.statistic),
    unname(baseline$raw.statistic) * 2.75^2,
    tolerance = 3e-12
  )
  expect_equal(
    transformed$components$trace.Sp,
    baseline$components$trace.Sp * 2.75^2,
    tolerance = 3e-12
  )
  expect_equal(
    transformed$components$trace.Sp2,
    baseline$components$trace.Sp2 * 2.75^4,
    tolerance = 3e-12
  )
  expect_equal(
    transformed$components$B,
    baseline$components$B * 2.75^4,
    tolerance = 3e-12
  )
  expect_equal(
    transformed$components$trace.Sigma2.hat,
    baseline$components$trace.Sigma2.hat * 2.75^4,
    tolerance = 3e-12
  )
  expect_equal(
    unname(transformed$variance),
    unname(baseline$variance) * 2.75^4,
    tolerance = 3e-12
  )

  swapped <- bai_saranadasa_two_sample_test(y, x)
  expect_equal(swapped$statistic, baseline$statistic, tolerance = 2e-12)
  expect_equal(swapped$p.value, baseline$p.value, tolerance = 2e-12)
  expect_equal(unname(swapped$estimate), -unname(baseline$estimate))
  expect_equal(unname(swapped$n), rev(unname(baseline$n)))
})

test_that("Bai-Saranadasa selects a dual Gram matrix when p exceeds total n", {
  set.seed(2045)
  x <- matrix(stats::rnorm(4 * 300), nrow = 4)
  y <- matrix(stats::rnorm(5 * 300), nrow = 5)
  fit <- bai_saranadasa_two_sample_test(x, y)

  expect_true(is.finite(unname(fit$statistic)))
  expect_identical(fit$diagnostics$trace.computation, "dual Gram matrix")
  expect_identical(fit$diagnostics$gram.type, "dual")
  expect_equal(fit$diagnostics$gram.dimension, nrow(x) + nrow(y))
  expect_identical(fit$diagnostics$constructs.p.by.p.matrix, FALSE)
})

test_that("Bai-Saranadasa standardisation survives extreme global scales", {
  set.seed(2046)
  x <- matrix(stats::rnorm(10 * 5), nrow = 10)
  y <- matrix(stats::rnorm(12 * 5, mean = 0.35), nrow = 12)
  baseline <- bai_saranadasa_two_sample_test(x, y)

  for (scale in c(1e-150, 1e150)) {
    extreme <- bai_saranadasa_two_sample_test(x * scale, y * scale)
    expect_true(all(is.finite(c(extreme$statistic, extreme$p.value))))
    expect_true(is.finite(extreme$diagnostics$internal.scaled.xi))
    expect_true(
      is.finite(extreme$diagnostics$internal.scaled.denominator.variance) &&
        extreme$diagnostics$internal.scaled.denominator.variance > 0
    )
    expect_equal(extreme$statistic, baseline$statistic, tolerance = 3e-12)
    expect_equal(extreme$p.value, baseline$p.value, tolerance = 3e-12)
    expect_equal(
      unname(extreme$estimate) / scale,
      unname(baseline$estimate),
      tolerance = 3e-12
    )
    expect_equal(
      extreme$components$trace.Sp / scale^2,
      baseline$components$trace.Sp,
      tolerance = 3e-12
    )
    expect_equal(
      unname(extreme$raw.statistic) / scale^2,
      unname(baseline$raw.statistic),
      tolerance = 3e-12
    )
  }

  huge <- bai_saranadasa_two_sample_test(x * 1e150, y * 1e150)
  tiny <- bai_saranadasa_two_sample_test(x * 1e-150, y * 1e-150)
  expect_true(is.infinite(unname(huge$variance)))
  expect_equal(unname(tiny$variance), 0)
  expect_true(all(is.finite(c(huge$statistic, huge$p.value))))
  expect_true(all(is.finite(c(tiny$statistic, tiny$p.value))))
})

test_that("Bai-Saranadasa keeps the standard location-test contract", {
  x <- matrix(c(0, 1, 2, 4), ncol = 1)
  y <- matrix(c(-1, 2, 3), ncol = 1)
  fit <- bai_saranadasa_two_sample_test(x, y)

  expect_named(
    fit,
    c(
      "statistic", "parameter", "p.value", "estimate",
      "null.value", "alternative", "method", "data.name", "call",
      "raw.statistic", "null.distribution", "variance", "components",
      "diagnostics", "n", "p"
    ),
    ignore.order = TRUE
  )
  expect_named(fit$statistic, "Z")
  expect_named(fit$raw.statistic, "xi")
  expect_named(fit$variance, "denominator")
  expect_true(is.finite(fit$p.value) && fit$p.value >= 0 && fit$p.value <= 1)
  expect_identical(fit$diagnostics$regularization, "none")
  expect_identical(fit$diagnostics$variance.repair, "none")
})

test_that("Bai-Saranadasa rejects incompatible and degenerate inputs", {
  x <- matrix(stats::rnorm(30), nrow = 6)
  y <- matrix(stats::rnorm(35), nrow = 7)
  expect_error(
    bai_saranadasa_two_sample_test(x[, 1:4], y),
    "same number of columns"
  )
  expect_error(
    bai_saranadasa_two_sample_test(x[1, , drop = FALSE], y),
    "at least 2 row"
  )

  colnames(x) <- paste0("x", seq_len(ncol(x)))
  colnames(y) <- paste0("y", seq_len(ncol(y)))
  expect_error(bai_saranadasa_two_sample_test(x, y), "same names")

  constant.x <- matrix(1, nrow = 4, ncol = 3)
  constant.y <- matrix(1, nrow = 5, ncol = 3)
  expect_error(
    bai_saranadasa_two_sample_test(constant.x, constant.y),
    "strictly positive"
  )

  nonfinite <- matrix(stats::rnorm(30), nrow = 6)
  nonfinite[1, 1] <- NA_real_
  expect_error(bai_saranadasa_two_sample_test(nonfinite, y), "finite")
})
