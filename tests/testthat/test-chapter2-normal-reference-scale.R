zzz23_fixture <- function() {
  list(
    x = rbind(
      c(0.7, -1.1, 0.2, 1.3),
      c(-0.4, 0.8, 1.2, -0.6),
      c(1.1, 0.3, -0.9, 0.5),
      c(-1.2, -0.5, 0.7, 1.0),
      c(0.2, 1.4, -0.4, -1.1),
      c(1.5, -0.2, 0.5, 0.1),
      c(-0.8, 0.6, -1.3, 0.9)
    ),
    y = rbind(
      c(-0.3, 0.7, -0.1, 1.1),
      c(1.2, -0.6, 0.8, -0.4),
      c(-1.1, 0.2, 1.3, 0.6),
      c(0.5, 1.1, -0.7, -0.9),
      c(0.9, -1.3, 0.4, 0.2),
      c(-0.6, -0.4, -1.2, 1.4),
      c(1.4, 0.5, -0.3, -0.8),
      c(-0.9, 1.3, 0.9, -0.2)
    )
  )
}


zzz23_reference <- function(x, y, alpha = 0.05) {
  n1 <- nrow(x)
  n2 <- nrow(y)
  n <- n1 + n2
  p <- ncol(x)
  weight1 <- n2 / n
  weight2 <- n1 / n
  S1 <- stats::cov(x)
  S2 <- stats::cov(y)
  if (p == 1L) {
    S1 <- matrix(S1, 1L, 1L)
    S2 <- matrix(S2, 1L, 1L)
  }
  omega <- weight1 * S1 + weight2 * S2
  D <- diag(omega)
  inverse.root <- diag(1 / sqrt(D), nrow = p, ncol = p)
  R1 <- inverse.root %*% S1 %*% inverse.root
  R2 <- inverse.root %*% S2 %*% inverse.root
  Rn <- weight1 * R1 + weight2 * R2
  difference <- colMeans(x) - colMeans(y)
  coordinate <- n1 * n2 / (n * p) * difference^2 / D
  T <- sum(coordinate)
  trace1 <- sum(diag(R1))
  trace2 <- sum(diag(R2))
  square1 <- sum(R1^2)
  square2 <- sum(R2^2)
  cross <- sum(R1 * R2)
  factor1 <- (n1 - 1)^2 / ((n1 - 2) * (n1 + 1))
  factor2 <- (n2 - 1)^2 / ((n2 - 2) * (n2 + 1))
  bracket1 <- square1 - trace1^2 / (n1 - 1)
  bracket2 <- square2 - trace2^2 / (n2 - 1)
  corrected1 <- factor1 * bracket1
  corrected2 <- factor2 * bracket2
  qhat <- weight1^2 * corrected1 + weight2^2 * corrected2 +
    2 * weight1 * weight2 * cross
  raw <- sum(Rn^2)
  df <- p^2 / qhat
  c.np <- 1 + raw / p^(3 / 2)
  df.paper <- if (c.np <= 1.2) df / c.np else df
  calibrate <- function(degrees) {
    list(
      p.value = stats::pchisq(
        degrees * T, degrees, lower.tail = FALSE
      ),
      log.p.value = stats::pchisq(
        degrees * T, degrees, lower.tail = FALSE, log.p = TRUE
      ),
      critical = stats::qchisq(1 - alpha, degrees) / degrees
    )
  }
  list(
    n1 = n1, n2 = n2, n = n, p = p,
    weight1 = weight1, weight2 = weight2,
    S1 = S1, S2 = S2, omega = omega, D = D,
    R1 = R1, R2 = R2, Rn = Rn,
    difference = difference, coordinate = coordinate, T = T,
    trace1 = trace1, trace2 = trace2,
    square1 = square1, square2 = square2, cross = cross,
    factor1 = factor1, factor2 = factor2,
    bracket1 = bracket1, bracket2 = bracket2,
    corrected1 = corrected1, corrected2 = corrected2,
    qhat = qhat, raw = raw, df = df, c.np = c.np,
    df.paper = df.paper,
    unadjusted = calibrate(df),
    paper = calibrate(df.paper)
  )
}


zzz23_fourier_design <- function(n, p, phase = 0) {
  stopifnot(p %% 2L == 0L)
  index <- seq.int(0, n - 1) + phase
  answer <- matrix(NA_real_, n, p)
  for (frequency in seq_len(p / 2L)) {
    angle <- 2 * pi * frequency * index / n
    answer[, 2L * frequency - 1L] <- sin(angle)
    answer[, 2L * frequency] <- cos(angle)
  }
  answer
}


test_that("Zhang-Zhu-Zhang matches the full-matrix feasible formulas", {
  fixture <- zzz23_fixture()
  reference <- zzz23_reference(fixture$x, fixture$y, alpha = 0.1)
  result <- zhang_zhu_zhang_two_sample_test(
    fixture$x, fixture$y, alpha = 0.1
  )
  component <- result$components

  expect_s3_class(result, "hd_location_test")
  expect_s3_class(result, "htest")
  expect_equal(unname(result$statistic), reference$T,
               tolerance = 3e-13)
  expect_equal(unname(result$raw.statistic), reference$T,
               tolerance = 3e-13)
  expect_equal(component$T.NRSI, reference$T, tolerance = 3e-13)
  expect_equal(unname(component$coordinate.contribution),
               reference$coordinate, tolerance = 3e-13)
  expect_equal(sum(component$coordinate.contribution), component$T.NRSI,
               tolerance = 3e-15)
  expect_equal(unname(component$D.hat.diagonal.input), reference$D,
               tolerance = 4e-13)
  expect_equal(component$trace.R1, reference$trace1,
               tolerance = 4e-13)
  expect_equal(component$trace.R2, reference$trace2,
               tolerance = 4e-13)
  expect_equal(component$trace.Rn, reference$p, tolerance = 4e-13)
  expect_equal(component$trace.R1.squared.raw, reference$square1,
               tolerance = 6e-13)
  expect_equal(component$trace.R2.squared.raw, reference$square2,
               tolerance = 6e-13)
  expect_equal(component$trace.R1.R2.raw, reference$cross,
               tolerance = 6e-13)
  expect_equal(component$trace.Rn.squared.raw, reference$raw,
               tolerance = 6e-13)
  expect_equal(component$trace.R1.squared.bracket, reference$bracket1,
               tolerance = 6e-13)
  expect_equal(component$trace.R2.squared.bracket, reference$bracket2,
               tolerance = 6e-13)
  expect_equal(component$trace.R1.squared.corrected,
               reference$corrected1, tolerance = 7e-13)
  expect_equal(component$trace.R2.squared.corrected,
               reference$corrected2, tolerance = 7e-13)
  expect_equal(component$trace.Rn.squared.corrected, reference$qhat,
               tolerance = 7e-13)
  expect_equal(component$df.hat, reference$df, tolerance = 7e-13)
  expect_equal(component$c.np, reference$c.np, tolerance = 5e-14)
  expect_equal(component$df.paper, reference$df.paper,
               tolerance = 7e-13)
  expect_equal(component$p.value.unadjusted,
               reference$unadjusted$p.value, tolerance = 2e-15)
  expect_equal(component$log.p.value.unadjusted,
               reference$unadjusted$log.p.value, tolerance = 2e-15)
  expect_equal(component$critical.value.unadjusted,
               reference$unadjusted$critical, tolerance = 2e-15)
  expect_equal(component$p.value.paper, reference$paper$p.value,
               tolerance = 2e-15)
  expect_equal(component$critical.value.paper,
               reference$paper$critical, tolerance = 2e-15)
  expect_equal(result$p.value, reference$paper$p.value,
               tolerance = 2e-15)
  expect_equal(unname(result$parameter), reference$df.paper,
               tolerance = 7e-13)
  expect_identical(result$null.distribution$family,
                   "scaled chi-square")
  expect_identical(result$null.distribution$tail, "upper")
  expect_false(result$diagnostics$raw.L2.normal.reference)
  expect_false(result$diagnostics$F.type.reference)
  expect_identical(result$diagnostics$regularization, "none")
  expect_identical(result$diagnostics$trace.floor, "none")
  expect_identical(result$diagnostics$degrees.of.freedom.clamp, "none")
  expect_identical(result$diagnostics$trace.computation, "primal")
})


test_that("crossed weights, Eq. 24, Eq. 25, c_np, and right tail are locked", {
  fixture <- zzz23_fixture()
  result <- zhang_zhu_zhang_two_sample_test(
    fixture$x, fixture$y, alpha = 0.07, df_adjustment = "none"
  )
  component <- result$components
  n1 <- nrow(fixture$x)
  n2 <- nrow(fixture$y)
  n <- n1 + n2
  p <- ncol(fixture$x)

  expect_equal(component$covariance.weights,
               c(group1 = n2 / n, group2 = n1 / n), tolerance = 0)
  expect_equal(component$trace.correction.factor,
               c(group1 = (n1 - 1)^2 / ((n1 - 2) * (n1 + 1)),
                 group2 = (n2 - 1)^2 / ((n2 - 2) * (n2 + 1))),
               tolerance = 2e-15)
  expect_equal(
    component$trace.R1.squared.corrected,
    component$trace.correction.factor[["group1"]] *
      (component$trace.R1.squared.raw -
         component$trace.R1^2 / (n1 - 1)),
    tolerance = 3e-15
  )
  expect_equal(
    component$trace.R2.squared.corrected,
    component$trace.correction.factor[["group2"]] *
      (component$trace.R2.squared.raw -
         component$trace.R2^2 / (n2 - 1)),
    tolerance = 3e-15
  )
  reconstructed <- (n2 / n)^2 *
    component$trace.R1.squared.corrected +
    (n1 / n)^2 * component$trace.R2.squared.corrected +
    2 * n1 * n2 / n^2 * component$trace.R1.R2.raw
  expect_equal(component$trace.Rn.squared.corrected, reconstructed,
               tolerance = 4e-15)
  expect_equal(component$df.hat,
               p^2 / component$trace.Rn.squared.corrected,
               tolerance = 3e-15)
  expect_equal(component$c.np,
               1 + component$trace.Rn.squared.raw / p^(3 / 2),
               tolerance = 3e-15)
  expect_gt(abs(
    component$c.np -
      (1 + component$trace.Rn.squared.corrected / p^(3 / 2))
  ), 1e-4)
  expect_equal(result$p.value,
               stats::pchisq(component$df.hat * component$T.NRSI,
                             component$df.hat, lower.tail = FALSE),
               tolerance = 0)
  expect_equal(component$critical.value.unadjusted,
               stats::qchisq(0.93, component$df.hat) /
                 component$df.hat,
               tolerance = 2e-15)
  expect_identical(
    result$diagnostics$rejection$reject,
    isTRUE(component$T.NRSI > component$critical.value.unadjusted)
  )
  expect_equal(unname(result$variance[["unadjusted.reference"]]),
               2 / component$df.hat, tolerance = 0)

  wrong.omega <- n1 / n * stats::cov(fixture$x) +
    n2 / n * stats::cov(fixture$y)
  wrong.T <- n1 * n2 / (n * p) * sum(
    (colMeans(fixture$x) - colMeans(fixture$y))^2 /
      diag(wrong.omega)
  )
  expect_gt(abs(component$T.NRSI - wrong.T), 1e-4)
})


test_that("paper and none df choices implement both threshold branches", {
  high.c <- zzz23_fixture()
  high.paper <- zhang_zhu_zhang_two_sample_test(
    high.c$x, high.c$y, df_adjustment = "paper"
  )
  high.none <- zhang_zhu_zhang_two_sample_test(
    high.c$x, high.c$y, df_adjustment = "none"
  )
  expect_gt(high.paper$components$c.np, 1.2)
  expect_false(high.paper$components$paper.correction.applied)
  expect_equal(high.paper$components$df.paper,
               high.paper$components$df.hat, tolerance = 0)
  expect_equal(high.paper$p.value, high.none$p.value, tolerance = 0)

  p <- 36L
  x <- zzz23_fourier_design(80L, p, phase = 0)
  y <- zzz23_fourier_design(83L, p, phase = 0.37)
  y <- sweep(y, 2L, seq_len(p) / 500, "+")
  low.paper <- zhang_zhu_zhang_two_sample_test(
    x, y, df_adjustment = "paper"
  )
  low.none <- zhang_zhu_zhang_two_sample_test(
    x, y, df_adjustment = "none"
  )

  expect_lte(low.paper$components$c.np, 1.2)
  expect_true(low.paper$components$paper.correction.applied)
  expect_equal(low.paper$components$df.paper,
               low.paper$components$df.hat / low.paper$components$c.np,
               tolerance = 3e-15)
  expect_equal(unname(low.paper$parameter),
               low.paper$components$df.paper, tolerance = 0)
  expect_equal(unname(low.none$parameter), low.none$components$df.hat,
               tolerance = 0)
  expect_equal(low.paper$statistic, low.none$statistic, tolerance = 0)
  expect_equal(low.paper$p.value, low.paper$components$p.value.paper,
               tolerance = 0)
  expect_equal(low.none$p.value,
               low.none$components$p.value.unadjusted, tolerance = 0)
  expect_false(isTRUE(all.equal(low.paper$p.value, low.none$p.value,
                                tolerance = 1e-12)))
})


test_that("streamed primal and dual traces agree with matrix references", {
  fixture <- zzz23_fixture()
  primal <- zhang_zhu_zhang_two_sample_test(fixture$x, fixture$y)
  expect_identical(primal$diagnostics$trace.computation, "primal")
  expect_false(primal$diagnostics$constructs.p.by.p.matrix)

  n1 <- 5L
  n2 <- 7L
  p <- 20L
  x <- outer(seq_len(n1), seq_len(p), function(i, j) {
    sin((i + 0.3) * j / 4) + cos((2 * i - j) / 7)
  })
  y <- outer(seq_len(n2), seq_len(p), function(i, j) {
    cos((i + 0.7) * j / 5) - sin((i + 2 * j) / 9) + j / 100
  })
  reference <- zzz23_reference(x, y)
  dual <- zhang_zhu_zhang_two_sample_test(x, y)

  expect_identical(dual$diagnostics$trace.computation, "dual")
  expect_false(dual$diagnostics$constructs.p.by.p.matrix)
  expect_equal(unname(dual$statistic), reference$T,
               tolerance = 2e-11)
  expect_equal(dual$components$trace.R1.squared.raw,
               reference$square1, tolerance = 3e-10)
  expect_equal(dual$components$trace.R2.squared.raw,
               reference$square2, tolerance = 3e-10)
  expect_equal(dual$components$trace.R1.R2.raw,
               reference$cross, tolerance = 3e-10)
  expect_equal(dual$components$trace.Rn.squared.corrected,
               reference$qhat, tolerance = 4e-10)
  expect_equal(dual$components$df.hat, reference$df,
               tolerance = 4e-10)
})


test_that("group, row, coordinate, translation, and scale invariance hold", {
  fixture <- zzz23_fixture()
  baseline <- zhang_zhu_zhang_two_sample_test(fixture$x, fixture$y)

  swapped <- zhang_zhu_zhang_two_sample_test(fixture$y, fixture$x)
  expect_equal(swapped$statistic, baseline$statistic,
               tolerance = 3e-13)
  expect_equal(swapped$p.value, baseline$p.value, tolerance = 2e-15)
  expect_equal(swapped$components$trace.R1.squared.corrected,
               baseline$components$trace.R2.squared.corrected,
               tolerance = 7e-13)

  row.result <- zhang_zhu_zhang_two_sample_test(
    fixture$x[c(7, 2, 5, 1, 6, 3, 4), , drop = FALSE],
    fixture$y[c(4, 8, 1, 6, 2, 7, 3, 5), , drop = FALSE]
  )
  expect_equal(row.result$statistic, baseline$statistic,
               tolerance = 3e-13)
  expect_equal(row.result$p.value, baseline$p.value, tolerance = 2e-15)

  permutation <- c(4, 2, 1, 3)
  column.result <- zhang_zhu_zhang_two_sample_test(
    fixture$x[, permutation, drop = FALSE],
    fixture$y[, permutation, drop = FALSE]
  )
  expect_equal(column.result$statistic, baseline$statistic,
               tolerance = 3e-13)
  expect_equal(column.result$p.value, baseline$p.value, tolerance = 2e-15)

  multiplier <- c(1e-150, -2e150, 3e-120, -4e120)
  shift.multiplier <- c(1e8, -2e8, 3e8, -4e8)
  x.scaled <- sweep(fixture$x, 2L, multiplier, "*")
  y.scaled <- sweep(fixture$y, 2L, multiplier, "*")
  common.shift <- multiplier * shift.multiplier
  x.scaled <- sweep(x.scaled, 2L, common.shift, "+")
  y.scaled <- sweep(y.scaled, 2L, common.shift, "+")
  scaled <- zhang_zhu_zhang_two_sample_test(x.scaled, y.scaled)

  expect_equal(scaled$statistic, baseline$statistic,
               tolerance = 3e-7)
  expect_equal(scaled$p.value, baseline$p.value, tolerance = 4e-8)
  expect_equal(scaled$components$df.hat,
               baseline$components$df.hat, tolerance = 8e-8)
  expected.log.D <- baseline$components$log.D.hat.diagonal.input +
    2 * log(abs(multiplier))
  expect_equal(
    unname(scaled$components$log.D.hat.diagonal.input -
             expected.log.D[[1L]]),
    unname(expected.log.D - expected.log.D[[1L]]),
    tolerance = 2e-7
  )
})


test_that("minimum sample size, p equals one, and zero group variance work", {
  x <- matrix(c(-1, 0, 1), ncol = 1L)
  y <- matrix(c(-0.5, 0.2, 0.9, 1.4), ncol = 1L)
  reference <- zzz23_reference(x, y)
  result <- zhang_zhu_zhang_two_sample_test(x, y)
  expect_equal(unname(result$statistic), reference$T,
               tolerance = 2e-14)
  expect_equal(result$components$trace.Rn.squared.corrected,
               reference$qhat, tolerance = 2e-14)
  expect_equal(result$components$df.hat, reference$df,
               tolerance = 2e-14)

  group1 <- cbind(rep(2, 5), c(-2, -1, 0, 1, 2))
  group2 <- cbind(c(-1, 0, 1, 2, 3, 4), c(2, -1, 3, 0, 1, -2))
  allowed <- zhang_zhu_zhang_two_sample_test(group1, group2)
  expect_equal(
    allowed$diagnostics$zero.group.variance.coordinates[["group1"]],
    1
  )
  expect_true(is.finite(allowed$statistic))
  expect_true(is.finite(allowed$p.value))
})


test_that("tails retain log probabilities after ordinary underflow", {
  pattern <- rbind(
    c(-3, -1, 2), c(-2, 1, -1), c(-1, 3, 0),
    c(1, -2, 3), c(2, 0, -3), c(3, -1, 1)
  )
  x <- pattern * 1e-100
  y <- matrix(rep(c(1, 2, 3), each = 7L), nrow = 7L)
  result <- zhang_zhu_zhang_two_sample_test(x, y)

  expect_identical(result$p.value, 0)
  expect_true(is.finite(result$components$log.p.value.paper))
  expect_lt(result$components$log.p.value.paper, -1e100)
  expect_identical(result$diagnostics$regularization, "none")
})


test_that("invalid inputs and degenerate calibration fail without repair", {
  fixture <- zzz23_fixture()
  expect_error(
    zhang_zhu_zhang_two_sample_test(fixture$x[1:2, ], fixture$y),
    "at least 3 row"
  )
  expect_error(
    zhang_zhu_zhang_two_sample_test(fixture$x, fixture$y[, 1:3]),
    "same number of columns"
  )
  expect_error(
    zhang_zhu_zhang_two_sample_test(fixture$x, fixture$y, alpha = 0),
    "strictly between"
  )
  expect_error(
    zhang_zhu_zhang_two_sample_test(fixture$x, fixture$y, alpha = 1),
    "strictly between"
  )
  expect_error(
    zhang_zhu_zhang_two_sample_test(
      fixture$x, fixture$y, df_adjustment = "always"
    ),
    "arg"
  )

  constant <- fixture$x
  constant[, 2] <- 4
  constant.y <- fixture$y
  constant.y[, 2] <- 4
  expect_error(
    zhang_zhu_zhang_two_sample_test(constant, constant.y),
    "positive variation|pooled marginal variance"
  )

  orthogonal <- cbind(
    c(1, -1, 1, -1),
    c(1, 1, -1, -1),
    c(1, -1, -1, 1)
  )
  disjoint1 <- cbind(orthogonal, matrix(0, 4, 3))
  disjoint2 <- cbind(matrix(0, 4, 3), orthogonal)
  expect_error(
    zhang_zhu_zhang_two_sample_test(disjoint1, disjoint2),
    "bias-corrected trace estimate"
  )
})
