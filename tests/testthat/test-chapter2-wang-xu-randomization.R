wx_ref_cq <- function(x, y) {
  n1 <- nrow(x)
  n2 <- nrow(y)
  within.x <- 0
  within.y <- 0
  cross <- 0
  for (i in seq_len(n1)) {
    for (j in seq_len(n1)) {
      if (i != j) within.x <- within.x + sum(x[i, ] * x[j, ])
    }
  }
  for (i in seq_len(n2)) {
    for (j in seq_len(n2)) {
      if (i != j) within.y <- within.y + sum(y[i, ] * y[j, ])
    }
  }
  for (i in seq_len(n1)) {
    for (j in seq_len(n2)) {
      cross <- cross + sum(x[i, ] * y[j, ])
    }
  }
  terms <- c(
    within.x = within.x / (n1 * (n1 - 1)),
    within.y = within.y / (n2 * (n2 - 1)),
    cross = -2 * cross / (n1 * n2)
  )
  c(T.CQ = sum(terms), terms)
}

wx_ref_half_differences <- function(x) {
  m <- nrow(x) %/% 2L
  answer <- matrix(0, m, ncol(x))
  for (i in seq_len(m)) {
    answer[i, ] <- (x[2L * i, ] - x[2L * i - 1L, ]) / 2
  }
  answer
}

wx_ref_randomized <- function(pseudo.x, pseudo.y, signs) {
  m1 <- nrow(pseudo.x)
  m2 <- nrow(pseudo.y)
  s1 <- signs[seq_len(m1)]
  s2 <- signs[m1 + seq_len(m2)]
  within.x <- 0
  within.y <- 0
  cross <- 0
  for (i in seq_len(m1 - 1L)) {
    for (j in seq.int(i + 1L, m1)) {
      within.x <- within.x +
        2 * s1[i] * s1[j] * sum(pseudo.x[i, ] * pseudo.x[j, ])
    }
  }
  for (i in seq_len(m2 - 1L)) {
    for (j in seq.int(i + 1L, m2)) {
      within.y <- within.y +
        2 * s2[i] * s2[j] * sum(pseudo.y[i, ] * pseudo.y[j, ])
    }
  }
  for (i in seq_len(m1)) {
    for (j in seq_len(m2)) {
      cross <- cross -
        2 * s1[i] * s2[j] * sum(pseudo.x[i, ] * pseudo.y[j, ])
    }
  }
  within.x / (m1 * (m1 - 1)) +
    within.y / (m2 * (m2 - 1)) + cross / (m1 * m2)
}

wx_ref_all_signs <- function(dimension) {
  as.matrix(expand.grid(rep(list(c(-1L, 1L)), dimension)))
}

wx_ref_u32_add <- function(a, b) (a + b) %% 2^32

wx_ref_u32_multiply <- function(a, b) {
  a0 <- a %% 2^16
  a1 <- floor(a / 2^16)
  b0 <- b %% 2^16
  b1 <- floor(b / 2^16)
  low <- a0 * b0
  middle <- (a0 * b1 + a1 * b0) %% 2^16
  (low + middle * 2^16) %% 2^32
}

wx_ref_u32_xor <- function(a, b) {
  powers <- 2^(0:31)
  bits.a <- floor(a / powers) %% 2
  bits.b <- floor(b / powers) %% 2
  sum(((bits.a + bits.b) %% 2) * powers)
}

wx_ref_mix32 <- function(value) {
  value <- wx_ref_u32_xor(value, floor(value / 2^16))
  value <- wx_ref_u32_multiply(value, 2146121005)
  value <- wx_ref_u32_xor(value, floor(value / 2^15))
  value <- wx_ref_u32_multiply(value, 2221713035)
  wx_ref_u32_xor(value, floor(value / 2^16))
}

wx_ref_counter_word <- function(counter, seed) {
  low <- counter %% 2^32
  high <- floor(counter / 2^32)
  high.key <- wx_ref_mix32(wx_ref_u32_add(high, 2654435769))
  keyed <- wx_ref_u32_add(
    wx_ref_u32_xor(seed, high.key),
    wx_ref_u32_multiply(2654435769, wx_ref_u32_add(low, 1))
  )
  wx_ref_mix32(keyed)
}

wx_ref_counter_signs <- function(B, dimension, seed) {
  answer <- matrix(-1L, B, dimension)
  for (b in seq_len(B)) {
    for (j in seq_len(dimension)) {
      counter <- (b - 1) * dimension + j - 1
      answer[b, j] <- if (wx_ref_counter_word(counter, seed) %% 2 == 1) 1L else -1L
    }
  }
  answer
}

wx_fixture <- function() {
  list(
    x = matrix(c(
      1, 2, -1,
      2, 0, 3,
      4, 1, 2,
      -1, 3, 0,
      5, -2, 1
    ), 5, 3, byrow = TRUE),
    y = matrix(c(
      0, 1, 2,
      3, -1, 1,
      2, 4, -2,
      -2, 2, 3,
      1, -3, 4,
      4, 0, -1
    ), 6, 3, byrow = TRUE)
  )
}

test_that("observed CQ and adjacent half-differences match literal formulas", {
  data <- wx_fixture()
  rownames(data$x) <- paste0("A", seq_len(nrow(data$x)))
  rownames(data$y) <- paste0("B", seq_len(nrow(data$y)))
  colnames(data$x) <- colnames(data$y) <- c("u", "v", "w")
  result <- wang_xu_approx_randomization_test(
    data$x, data$y, calibration = "exact", keep_randomized = TRUE
  )

  literal <- wx_ref_cq(data$x, data$y)
  centered.identity <- sum((colMeans(data$x) - colMeans(data$y))^2) -
    sum(diag(stats::cov(data$x))) / nrow(data$x) -
    sum(diag(stats::cov(data$y))) / nrow(data$y)
  expect_s3_class(result, "htest")
  expect_s3_class(result, "hd_location_test")
  expect_equal(unname(result$statistic), unname(literal[["T.CQ"]]),
               tolerance = 2e-14)
  expect_equal(unname(result$raw.statistic), unname(literal[["T.CQ"]]),
               tolerance = 2e-14)
  expect_equal(unname(result$statistic), centered.identity,
               tolerance = 2e-14)
  expect_equal(result$estimate, colMeans(data$x) - colMeans(data$y),
               tolerance = 2e-14, ignore_attr = TRUE)
  expect_equal(result$components$pseudo.x,
               wx_ref_half_differences(data$x),
               tolerance = 2e-14, ignore_attr = TRUE)
  expect_equal(result$components$pseudo.y,
               wx_ref_half_differences(data$y),
               tolerance = 2e-14, ignore_attr = TRUE)
  expect_equal(unname(result$components$pair.indices.x),
               matrix(c(1L, 2L, 3L, 4L), 2, 2, byrow = TRUE))
  expect_equal(unname(result$components$pair.indices.y),
               matrix(1:6, 3, 2, byrow = TRUE))
  expect_identical(unname(result$components$discarded.rows.x), 5L)
  expect_length(result$components$discarded.rows.y, 0L)
  expect_identical(result$diagnostics$observed.statistic.sample,
                   "all original observations")
  expect_false(result$diagnostics$pooled.label.permutation)
  expect_true(result$diagnostics$pairing.order.sensitive)
  expect_identical(result$alternative, "two.sided")
})

test_that("exact calibration exhausts the literal Rademacher law", {
  data <- wx_fixture()
  x <- data$x[1:4, , drop = FALSE]
  y <- data$y[1:4, , drop = FALSE]
  result <- wang_xu_approx_randomization_test(
    x, y, calibration = "exact", seed = 123,
    max_exact = 8L, keep_randomized = TRUE
  )
  pseudo.x <- wx_ref_half_differences(x)
  pseudo.y <- wx_ref_half_differences(y)
  all.signs <- wx_ref_all_signs(4L)
  all.statistics <- apply(
    all.signs, 1L, wx_ref_randomized,
    pseudo.x = pseudo.x, pseudo.y = pseudo.y
  )
  expected.p <- mean(all.statistics >= unname(wx_ref_cq(x, y)[["T.CQ"]]))

  expect_equal(result$p.value, expected.p, tolerance = 0)
  expect_equal(sort(rep(result$components$randomized.statistics, each = 2L)),
               sort(all.statistics), tolerance = 2e-14)
  expect_equal(
    vapply(seq_len(nrow(result$components$sign.patterns)), function(i) {
      wx_ref_randomized(
        pseudo.x, pseudo.y, result$components$sign.patterns[i, ]
      )
    }, numeric(1)),
    result$components$randomized.statistics,
    tolerance = 2e-14
  )
  expect_equal(
    vapply(seq_len(nrow(result$components$sign.patterns)), function(i) {
      signs <- result$components$sign.patterns[i, ]
      as.numeric(crossprod(
        signs,
        result$components$randomized.quadratic.kernel.scaled %*% signs
      ))
    }, numeric(1)),
    result$components$randomized.statistics.scaled,
    tolerance = 2e-14
  )
  expect_true(all(result$components$sign.patterns[, 1L] == 1L))
  expect_equal(result$components$reference.evaluations, 8)
  expect_equal(result$components$total.sign.configurations, 16)
  expect_false(result$diagnostics$plus.one.correction)
  expect_true(result$diagnostics$conditional.reference.exhaustive)
  expect_true(result$diagnostics$global.sign.symmetry.reduced)
  expect_false(result$diagnostics$finite.sample.test.exact)
  expect_identical(result$diagnostics$tie.rule, ">=")

  set.seed(991)
  state <- .Random.seed
  invisible(wang_xu_approx_randomization_test(
    x, y, calibration = "exact", seed = 4294967295
  ))
  expect_identical(.Random.seed, state)
})

test_that("Monte Carlo signs match the independent counter reference", {
  data <- wx_fixture()
  x <- data$x[1:4, , drop = FALSE]
  y <- data$y[1:4, , drop = FALSE]
  B <- 19L
  seed <- 314159
  result <- wang_xu_approx_randomization_test(
    x, y, calibration = "monte_carlo", B = B, seed = seed,
    keep_randomized = TRUE
  )
  expected.signs <- wx_ref_counter_signs(B, 4L, seed)
  pseudo.x <- wx_ref_half_differences(x)
  pseudo.y <- wx_ref_half_differences(y)
  expected.statistics <- vapply(seq_len(B), function(i) {
    wx_ref_randomized(pseudo.x, pseudo.y, expected.signs[i, ])
  }, numeric(1))
  observed <- unname(wx_ref_cq(x, y)[["T.CQ"]])
  exceedances <- sum(expected.statistics >= observed)

  expect_equal(unname(result$components$sign.patterns), expected.signs)
  expect_equal(result$components$randomized.statistics,
               expected.statistics, tolerance = 2e-14)
  expect_equal(result$components$exceedances, exceedances)
  expect_equal(result$p.value, (1 + exceedances) / (B + 1), tolerance = 0)
  expect_true(result$diagnostics$plus.one.correction)
  expect_false(result$diagnostics$conditional.reference.exhaustive)
  expect_equal(result$diagnostics$minimum.attainable.p, 1 / (B + 1))
  expect_equal(result$diagnostics$seed.used, seed)
  expect_equal(result$diagnostics$workers, 1L)
  expect_false(result$diagnostics$parallel.execution)

  repeat.result <- wang_xu_approx_randomization_test(
    x, y, calibration = "monte_carlo", B = B, seed = seed,
    keep_randomized = TRUE
  )
  expect_identical(result$components$sign.patterns,
                   repeat.result$components$sign.patterns)
  expect_identical(result$p.value, repeat.result$p.value)

  set.seed(992)
  state <- .Random.seed
  invisible(wang_xu_approx_randomization_test(
    x, y, calibration = "monte_carlo", B = 5L, seed = seed
  ))
  expect_identical(.Random.seed, state)

  set.seed(993)
  first <- wang_xu_approx_randomization_test(
    x, y, calibration = "monte_carlo", B = 7L, seed = NULL,
    keep_randomized = TRUE
  )
  set.seed(993)
  second <- wang_xu_approx_randomization_test(
    x, y, calibration = "monte_carlo", B = 7L, seed = NULL,
    keep_randomized = TRUE
  )
  expect_identical(first$diagnostics$seed.used, second$diagnostics$seed.used)
  expect_identical(first$components$sign.patterns,
                   second$components$sign.patterns)
})

test_that("auto selection and exact caps are explicit", {
  data <- wx_fixture()
  x <- data$x[1:4, , drop = FALSE]
  y <- data$y[1:4, , drop = FALSE]
  exact <- wang_xu_approx_randomization_test(
    x, y, calibration = "auto", max_exact = 8L
  )
  monte <- wang_xu_approx_randomization_test(
    x, y, calibration = "auto", max_exact = 7L,
    B = 11L, seed = 7
  )
  expect_identical(exact$diagnostics$calibration.used, "exact")
  expect_identical(monte$diagnostics$calibration.used, "monte_carlo")
  expect_equal(exact$diagnostics$reference.evaluations, 8)
  expect_equal(monte$diagnostics$reference.evaluations, 11)
  expect_error(
    wang_xu_approx_randomization_test(
      x, y, calibration = "exact", max_exact = 7L
    ),
    "max_exact"
  )
  compact <- wang_xu_approx_randomization_test(
    x, y, calibration = "exact", keep_randomized = FALSE
  )
  expect_null(compact$components$randomized.statistics)
  expect_null(compact$components$randomized.statistics.scaled)
  expect_null(compact$components$sign.patterns)
})

test_that("the reference respects geometric invariances but pairing is ordered", {
  data <- wx_fixture()
  x <- data$x
  y <- data$y
  baseline <- wang_xu_approx_randomization_test(
    x, y, calibration = "exact", keep_randomized = TRUE
  )
  translation <- c(1e12, -2e12, 3e12)
  translated <- wang_xu_approx_randomization_test(
    sweep(x, 2L, translation, "+"),
    sweep(y, 2L, translation, "+"),
    calibration = "exact", keep_randomized = TRUE
  )
  theta <- 0.37
  rotation <- matrix(c(
    cos(theta), -sin(theta), 0,
    sin(theta), cos(theta), 0,
    0, 0, 1
  ), 3, 3, byrow = TRUE)
  rotated <- wang_xu_approx_randomization_test(
    x %*% rotation, y %*% rotation,
    calibration = "exact", keep_randomized = TRUE
  )
  scaled <- wang_xu_approx_randomization_test(
    -3.5 * x, -3.5 * y,
    calibration = "exact", keep_randomized = TRUE
  )
  swapped <- wang_xu_approx_randomization_test(
    y, x, calibration = "exact", keep_randomized = TRUE
  )

  expect_equal(translated$p.value, baseline$p.value, tolerance = 0)
  expect_equal(unname(translated$statistic), unname(baseline$statistic),
               tolerance = 2e-3)
  expect_equal(rotated$p.value, baseline$p.value, tolerance = 0)
  expect_equal(unname(rotated$statistic), unname(baseline$statistic),
               tolerance = 2e-13)
  expect_equal(scaled$p.value, baseline$p.value, tolerance = 0)
  expect_equal(unname(scaled$statistic),
               3.5^2 * unname(baseline$statistic), tolerance = 2e-13)
  expect_equal(swapped$p.value, baseline$p.value, tolerance = 0)
  expect_equal(unname(swapped$statistic), unname(baseline$statistic),
               tolerance = 2e-14)

  reverse.within.pairs <- c(2L, 1L, 4L, 3L, 5L)
  reversed <- wang_xu_approx_randomization_test(
    x[reverse.within.pairs, , drop = FALSE], y,
    calibration = "exact", keep_randomized = TRUE
  )
  expect_equal(reversed$p.value, baseline$p.value, tolerance = 0)
  expect_equal(sort(reversed$components$randomized.statistics),
               sort(baseline$components$randomized.statistics),
               tolerance = 2e-14)

  x.order <- matrix(c(0, 10, 11, 12, 13, 100), ncol = 1)
  y.order <- matrix(c(0, 1, 2, 3, 4, 5), ncol = 1)
  ordered <- wang_xu_approx_randomization_test(
    x.order, y.order, calibration = "exact", keep_randomized = TRUE
  )
  reordered <- wang_xu_approx_randomization_test(
    x.order[c(1, 3, 5, 2, 4, 6), , drop = FALSE], y.order,
    calibration = "exact", keep_randomized = TRUE
  )
  expect_equal(unname(ordered$statistic), unname(reordered$statistic),
               tolerance = 2e-14)
  expect_false(isTRUE(all.equal(
    sort(ordered$components$randomized.statistics),
    sort(reordered$components$randomized.statistics), tolerance = 1e-14
  )))
})

test_that("ties, point-mass references, and plus-one endpoints are literal", {
  x0 <- matrix(rep(c(2, -1), each = 4L), 4, 2)
  y0 <- x0
  exact.tie <- wang_xu_approx_randomization_test(
    x0, y0, calibration = "exact", keep_randomized = TRUE
  )
  mc.tie <- wang_xu_approx_randomization_test(
    x0, y0, calibration = "monte_carlo", B = 7L, seed = 1,
    keep_randomized = TRUE
  )
  expect_equal(unname(exact.tie$statistic), 0, tolerance = 0)
  expect_true(all(exact.tie$components$randomized.statistics == 0))
  expect_equal(exact.tie$p.value, 1, tolerance = 0)
  expect_equal(mc.tie$p.value, 1, tolerance = 0)
  expect_true(exact.tie$diagnostics$randomization.degenerate)
  expect_true(mc.tie$diagnostics$randomization.degenerate)
  expect_equal(mc.tie$components$exceedances, 7)

  y.shifted <- y0
  y.shifted[, 1L] <- y.shifted[, 1L] + 1
  exact.none <- wang_xu_approx_randomization_test(
    x0, y.shifted, calibration = "exact", keep_randomized = TRUE
  )
  mc.none <- wang_xu_approx_randomization_test(
    x0, y.shifted, calibration = "monte_carlo", B = 7L, seed = 1,
    keep_randomized = TRUE
  )
  expect_gt(unname(exact.none$statistic), 0)
  expect_true(all(exact.none$components$randomized.statistics == 0))
  expect_equal(exact.none$p.value, 0, tolerance = 0)
  expect_equal(mc.none$components$exceedances, 0)
  expect_equal(mc.none$p.value, 1 / 8, tolerance = 0)
  expect_equal(mc.none$diagnostics$minimum.attainable.p, 1 / 8)
  expect_identical(mc.none$diagnostics$p.value.repair, "none")
})

test_that("univariate, extreme-scale, and no-repair contracts are explicit", {
  x <- matrix(c(-3, -1, 2, 5), ncol = 1)
  y <- matrix(c(-2, 0, 4, 7), ncol = 1)
  base <- wang_xu_approx_randomization_test(
    x, y, calibration = "exact", keep_randomized = TRUE
  )
  large <- wang_xu_approx_randomization_test(
    x * 1e145, y * 1e145,
    calibration = "exact", keep_randomized = TRUE
  )
  small <- wang_xu_approx_randomization_test(
    x * 1e-145, y * 1e-145,
    calibration = "exact", keep_randomized = TRUE
  )
  expect_equal(base$p, 1L)
  expect_equal(large$p.value, base$p.value, tolerance = 0)
  expect_equal(small$p.value, base$p.value, tolerance = 0)
  expect_equal(unname(large$statistic) / 1e290,
               unname(base$statistic), tolerance = 2e-14)
  expect_equal(unname(small$statistic) / 1e-290,
               unname(base$statistic), tolerance = 2e-14)
  expect_identical(base$diagnostics$covariance.regularization, "none")
  expect_identical(base$diagnostics$p.value.repair, "none")
  expect_true(all(is.finite(base$components$randomized.statistics.scaled)))

  extremes.x <- matrix(c(
    .Machine$double.xmax, 0,
    -.Machine$double.xmax, 0,
    .Machine$double.xmax, 1,
    -.Machine$double.xmax, -1
  ), 4, 2, byrow = TRUE)
  extremes.y <- -extremes.x
  expect_error(
    wang_xu_approx_randomization_test(
      extremes.x, extremes.y, calibration = "exact"
    ),
    "not representable"
  )
})

test_that("invalid controls and sample boundaries fail clearly", {
  data <- wx_fixture()
  x <- data$x[1:4, , drop = FALSE]
  y <- data$y[1:4, , drop = FALSE]
  expect_error(wang_xu_approx_randomization_test(x[1:3, ], y),
               "at least 4 row")
  expect_error(wang_xu_approx_randomization_test(x, y[1:3, ]),
               "at least 4 row")
  expect_error(wang_xu_approx_randomization_test(x, y[, 1:2]),
               "same number of columns")
  bad <- x
  bad[1, 1] <- NA_real_
  expect_error(wang_xu_approx_randomization_test(bad, y),
               "only finite")
  expect_error(wang_xu_approx_randomization_test(x, y, alpha = 0),
               "strictly between")
  expect_error(wang_xu_approx_randomization_test(x, y, B = 0),
               "positive integer")
  expect_error(wang_xu_approx_randomization_test(x, y, seed = -1),
               "2\\^32")
  expect_error(wang_xu_approx_randomization_test(x, y, seed = 2^32),
               "2\\^32")
  expect_error(wang_xu_approx_randomization_test(x, y, seed = 1.5),
               "integer-valued")
  expect_error(wang_xu_approx_randomization_test(x, y, workers = 2),
               "not implemented")
  expect_error(wang_xu_approx_randomization_test(x, y, max_exact = 0),
               "positive integer")
  expect_error(wang_xu_approx_randomization_test(
    x, y, keep_randomized = NA
  ), "TRUE.*FALSE")
  expect_error(wang_xu_approx_randomization_test(
    x, y, calibration = "labels"
  ), "arg")

  named.x <- x
  named.y <- y
  colnames(named.x) <- c("a", "b", "c")
  colnames(named.y) <- c("a", "c", "b")
  expect_error(wang_xu_approx_randomization_test(named.x, named.y),
               "same names")
})
