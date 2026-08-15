skk_reference <- function(x, y) {
  n1 <- nrow(x)
  n2 <- nrow(y)
  p <- ncol(x)
  mean.x <- colMeans(x)
  mean.y <- colMeans(y)
  centered.x <- sweep(x, 2L, mean.x, "-")
  centered.y <- sweep(y, 2L, mean.y, "-")
  s1 <- crossprod(centered.x) / (n1 - 1)
  s2 <- crossprod(centered.y) / (n2 - 1)
  d <- diag(s1) / n1 + diag(s2) / n2
  difference <- mean.x - mean.y
  Q <- sum(difference^2 / d)
  q <- (Q - p) / sqrt(p)

  root.outer <- outer(sqrt(d), sqrt(d))
  a1 <- s1 / root.outer
  a2 <- s2 / root.outer
  trace.a1 <- sum(diag(a1))
  trace.a2 <- sum(diag(a2))
  trace.a1.squared <- sum(a1 * t(a1))
  trace.a2.squared <- sum(a2 * t(a2))
  trace.a1.a2 <- sum(a1 * t(a2))
  F1 <- (trace.a1.squared - trace.a1^2 / (n1 - 1)) / p
  F2 <- (trace.a2.squared - trace.a2^2 / (n2 - 1)) / p
  G <- trace.a1.a2 / p
  variance.q <- 2 * F1 / n1^2 + 2 * F2 / n2^2 +
    4 * G / (n1 * n2)
  r.hat <- a1 / n1 + a2 / n2
  trace.r.hat2 <- sum(r.hat * t(r.hat))
  c.hat <- 1 + trace.r.hat2 / p^(3 / 2)
  denominator.variance <- variance.q * c.hat

  list(
    mean.x = mean.x,
    mean.y = mean.y,
    difference = difference,
    S1 = s1,
    S2 = s2,
    D = d,
    A1 = a1,
    A2 = a2,
    Q = Q,
    q = q,
    trace.A1 = trace.a1,
    trace.A2 = trace.a2,
    trace.A1.squared = trace.a1.squared,
    trace.A2.squared = trace.a2.squared,
    trace.A1.A2 = trace.a1.a2,
    F1 = F1,
    F2 = F2,
    G = G,
    variance.q = variance.q,
    trace.R.hat2 = trace.r.hat2,
    c.hat = c.hat,
    denominator.variance = denominator.variance,
    z = q / sqrt(denominator.variance)
  )
}


test_that("SKK matches a literal corrected-formula reference", {
  set.seed(2301)
  x <- matrix(stats::rnorm(32), 8, 4)
  y <- matrix(stats::rnorm(36, 0.25), 9, 4)
  x <- sweep(x, 2L, c(0.5, 2, 4, 0.25), "*")
  y <- sweep(y, 2L, c(0.5, 2, 4, 0.25), "*")
  reference <- skk_reference(x, y)
  result <- srivastava_katayama_kano_two_sample_test(x, y)

  expect_s3_class(result, "hd_location_test")
  expect_s3_class(result, "htest")
  expect_equal(unname(result$statistic), reference$z, tolerance = 3e-12)
  expect_equal(result$p.value,
               stats::pnorm(reference$z, lower.tail = FALSE),
               tolerance = 1e-14)
  expect_equal(result$components$Q.raw, reference$Q, tolerance = 2e-13)
  expect_equal(result$components$q.scaled, reference$q, tolerance = 2e-13)
  expect_equal(result$components$F1, reference$F1, tolerance = 2e-12)
  expect_equal(result$components$F2, reference$F2, tolerance = 2e-12)
  expect_equal(result$components$G, reference$G, tolerance = 2e-12)
  expect_equal(result$components$q.variance, reference$variance.q,
               tolerance = 2e-12)
  expect_equal(result$components$trace.R.hat2, reference$trace.R.hat2,
               tolerance = 2e-12)
  expect_equal(result$components$c.hat, reference$c.hat,
               tolerance = 2e-12)
  expect_equal(result$variance[[1L]], reference$denominator.variance,
               tolerance = 2e-12)
  expect_equal(unname(result$components$D.hat.diagonal), reference$D,
               tolerance = 2e-13)
  expect_equal(unname(result$components$variance1.diagonal),
               unname(diag(reference$S1)), tolerance = 2e-13)
  expect_equal(unname(result$components$variance2.diagonal),
               unname(diag(reference$S2)), tolerance = 2e-13)
})


test_that("SKK uses sqrt(p), n_k - 1, and sample R-hat corrections", {
  set.seed(2302)
  x <- matrix(stats::rnorm(35), 7, 5)
  y <- matrix(stats::rnorm(45, 0.15), 9, 5)
  reference <- skk_reference(x, y)
  result <- srivastava_katayama_kano_two_sample_test(x, y)

  expect_equal(result$components$q.scaled,
               (result$components$Q.raw - 5) / sqrt(5),
               tolerance = 1e-14)
  expect_gt(abs(result$components$q.scaled -
                  (result$components$Q.raw - 5)), 1e-5)

  wrong.F1 <- (reference$trace.A1.squared -
                 reference$trace.A1^2 / nrow(x)) / ncol(x)
  wrong.F2 <- (reference$trace.A2.squared -
                 reference$trace.A2^2 / nrow(y)) / ncol(y)
  expect_equal(result$components$F1, reference$F1, tolerance = 1e-13)
  expect_equal(result$components$F2, reference$F2, tolerance = 1e-13)
  expect_gt(abs(result$components$F1 - wrong.F1), 1e-5)
  expect_gt(abs(result$components$F2 - wrong.F2), 1e-5)

  expect_equal(result$components$trace.R.hat2,
               sum((reference$A1 / nrow(x) + reference$A2 / nrow(y))^2),
               tolerance = 2e-13)
  expect_equal(result$components$c.hat,
               1 + result$components$trace.R.hat2 / ncol(x)^(3 / 2),
               tolerance = 1e-14)
  expect_equal(
    result$diagnostics$F.bias.correction.denominators,
    c(group1 = nrow(x) - 1, group2 = nrow(y) - 1)
  )
})


test_that("SKK agrees with the reference in primal and dual dimensions", {
  dimensions <- list(
    c(n1 = 7, n2 = 8, p = 3),
    c(n1 = 5, n2 = 6, p = 24)
  )
  expected.gram <- c("primal", "dual")
  expected.dimension <- c(3, 11)

  for (i in seq_along(dimensions)) {
    set.seed(2310 + i)
    d <- dimensions[[i]]
    x <- matrix(stats::rnorm(d[["n1"]] * d[["p"]]),
                d[["n1"]], d[["p"]])
    y <- matrix(stats::rnorm(d[["n2"]] * d[["p"]], 0.1),
                d[["n2"]], d[["p"]])
    result <- srivastava_katayama_kano_two_sample_test(x, y)
    reference <- skk_reference(x, y)

    expect_equal(unname(result$statistic), reference$z, tolerance = 8e-12)
    expect_equal(result$components$F1, reference$F1, tolerance = 8e-12)
    expect_equal(result$components$F2, reference$F2, tolerance = 8e-12)
    expect_equal(result$components$G, reference$G, tolerance = 8e-12)
    expect_identical(result$diagnostics$gram.type, expected.gram[[i]])
    expect_equal(result$diagnostics$gram.dimension,
                 expected.dimension[[i]])
    expect_identical(
      result$diagnostics$constructs.p.by.p.matrix,
      identical(expected.gram[[i]], "primal")
    )
  }
})


test_that("SKK is invariant to diagonal units and exchanging groups", {
  set.seed(2321)
  x <- matrix(stats::rnorm(40), 8, 5)
  y <- matrix(stats::rnorm(50, 0.2), 10, 5)
  baseline <- srivastava_katayama_kano_two_sample_test(x, y)

  units <- c(-4, 0.25, 7, -2.5, 0.1)
  rescaled <- srivastava_katayama_kano_two_sample_test(
    sweep(x, 2L, units, "*"), sweep(y, 2L, units, "*")
  )
  expect_equal(rescaled$statistic, baseline$statistic, tolerance = 3e-12)
  expect_equal(rescaled$components$Q.raw, baseline$components$Q.raw,
               tolerance = 3e-12)
  expect_equal(rescaled$components$q.variance,
               baseline$components$q.variance, tolerance = 3e-12)
  expect_equal(rescaled$components$c.hat, baseline$components$c.hat,
               tolerance = 3e-12)

  swapped <- srivastava_katayama_kano_two_sample_test(y, x)
  expect_equal(swapped$statistic, baseline$statistic, tolerance = 3e-12)
  expect_equal(swapped$components$Q.raw, baseline$components$Q.raw,
               tolerance = 3e-12)
  expect_equal(swapped$components$F1, baseline$components$F2,
               tolerance = 3e-12)
  expect_equal(swapped$components$F2, baseline$components$F1,
               tolerance = 3e-12)
  expect_equal(swapped$components$G, baseline$components$G,
               tolerance = 3e-12)
  expect_equal(swapped$components$difference,
               -baseline$components$difference,
               tolerance = 1e-14)

  shift <- c(10, -7, 3, 20, -12)
  translated <- srivastava_katayama_kano_two_sample_test(
    sweep(x, 2L, shift, "+"), sweep(y, 2L, shift, "+")
  )
  expect_equal(translated$statistic, baseline$statistic, tolerance = 8e-12)
})


test_that("SKK internal scaling handles extreme finite units", {
  set.seed(2322)
  x <- matrix(stats::rnorm(45), 9, 5)
  y <- matrix(stats::rnorm(55, 0.1), 11, 5)
  baseline <- srivastava_katayama_kano_two_sample_test(x, y)

  expect_equal(
    srivastava_katayama_kano_two_sample_test(1e150 * x, 1e150 * y)$statistic,
    baseline$statistic,
    tolerance = 8e-12
  )
  expect_equal(
    srivastava_katayama_kano_two_sample_test(1e-150 * x, 1e-150 * y)$statistic,
    baseline$statistic,
    tolerance = 8e-12
  )

  extreme.units <- c(-3e150, 2e-150, 7, -4e75, 0.125)
  extreme <- srivastava_katayama_kano_two_sample_test(
    sweep(x, 2L, extreme.units, "*"),
    sweep(y, 2L, extreme.units, "*")
  )
  expect_true(is.finite(unname(extreme$statistic)))
  expect_true(is.finite(extreme$p.value))
  expect_equal(extreme$statistic, baseline$statistic, tolerance = 8e-12)
})


test_that("SKK rejects malformed and degenerate inputs without repair", {
  expect_error(
    srivastava_katayama_kano_two_sample_test(
      matrix(1:3, 1, 3), matrix(stats::rnorm(9), 3, 3)
    ),
    "at least 2 row"
  )
  expect_error(
    srivastava_katayama_kano_two_sample_test(
      matrix(stats::rnorm(12), 4, 3), matrix(stats::rnorm(8), 4, 2)
    ),
    "same number of columns"
  )

  x <- cbind(stats::rnorm(6), rep(2, 6))
  y <- cbind(stats::rnorm(7), rep(-3, 7))
  expect_error(
    srivastava_katayama_kano_two_sample_test(x, y),
    "combined marginal variance"
  )

  x.orthogonal <- rbind(c(-1, 0), c(1, 0))
  y.orthogonal <- rbind(c(0, -1), c(0, 1))
  expect_error(
    srivastava_katayama_kano_two_sample_test(x.orthogonal, y.orthogonal),
    "strictly positive variance estimate for q"
  )

  x.named <- matrix(stats::rnorm(18), 6, 3,
                    dimnames = list(NULL, c("a", "b", "c")))
  y.named <- matrix(stats::rnorm(21), 7, 3,
                    dimnames = list(NULL, c("b", "a", "c")))
  expect_error(
    srivastava_katayama_kano_two_sample_test(x.named, y.named),
    "same names in the same order"
  )
  x.named[1, 1] <- Inf
  expect_error(
    srivastava_katayama_kano_two_sample_test(x.named, unname(y.named)),
    "finite values"
  )
})


test_that("SKK preserves names and reports its numerical policy", {
  set.seed(2330)
  variables <- c("height", "width", "depth")
  x <- matrix(stats::rnorm(21), 7, 3,
              dimnames = list(NULL, variables))
  y <- matrix(stats::rnorm(24), 8, 3,
              dimnames = list(NULL, variables))
  result <- srivastava_katayama_kano_two_sample_test(x, y)

  expect_identical(names(result$estimate), variables)
  expect_identical(names(result$components$D.hat.diagonal), variables)
  expect_identical(names(result$diagnostics$internal.scale.factors), variables)
  expect_identical(result$diagnostics$regularization, "none")
  expect_identical(result$diagnostics$variance.repair, "none")
  expect_identical(result$null.distribution$tail, "upper")
  expect_equal(result$n, c(group1 = 7, group2 = 8))
  expect_equal(result$p, 3L)
})
