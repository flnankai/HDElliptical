zf_one_rank_reference <- function(x, mu = rep(0, ncol(x))) {
  residual <- sweep(x, 2L, mu, "-")
  vapply(seq_len(ncol(x)), function(j) {
    sum(rank(abs(residual[, j]), ties.method = "average") *
          (residual[, j] > 0))
  }, numeric(1L))
}

zf_two_rank_reference <- function(x, y) {
  n <- nrow(x)
  vapply(seq_len(ncol(x)), function(j) {
    pooled.rank <- rank(c(x[, j], y[, j]), ties.method = "average")
    sum(pooled.rank[seq_len(n)]) - n * (n + 1) / 2
  }, numeric(1L))
}

zf_one_moment_reference <- function(n) {
  mean.u <- n * (n + 1) / 4
  variance.u <- n * (n + 1) * (2 * n + 1) / 24
  variance.square <-
    (6 * n + 5 * n^2 - 30 * n^3 - 25 * n^4 +
       24 * n^5 + 20 * n^6) / 1440
  c(
    mean.rank.sum = mean.u,
    variance.rank.sum = variance.u,
    mean.squared.rank.score = variance.u,
    variance.squared.rank.score = variance.square
  )
}

zf_two_moment_reference <- function(n, m) {
  total <- n + m
  product <- n * m
  mean.u <- product / 2
  variance.u <- product * (total + 1) / 12
  variance.square <-
    (product * (5 * total + 8) - 3 * total * (total + 1)) *
    (total + 1) * product / 360
  c(
    mean.rank.sum = mean.u,
    variance.rank.sum = variance.u,
    mean.squared.rank.score = variance.u,
    variance.squared.rank.score = variance.square
  )
}

zf_statistic_reference <- function(rank.sum, moments, tau.squared) {
  p <- length(rank.sum)
  centered <- rank.sum - moments[["mean.rank.sum"]]
  marginal <- centered / sqrt(moments[["variance.rank.sum"]])
  square <- centered^2
  maximum <- max(marginal^2) - 2 * log(p) + log(log(p))
  max.p <- -expm1(-exp(-maximum / 2) / sqrt(pi))
  sum.statistic <- sqrt(p) *
    (mean(square) - moments[["mean.squared.rank.score"]]) /
    sqrt(moments[["variance.squared.rank.score"]] * tau.squared)
  sum.p <- stats::pnorm(sum.statistic, lower.tail = FALSE)
  cauchy <- 0.5 / tan(pi * max.p) + 0.5 / tan(pi * sum.p)
  list(
    rank.sum = rank.sum,
    centered = centered,
    marginal = marginal,
    square = square,
    maximum = maximum,
    max.p = max.p,
    sum = sum.statistic,
    sum.p = sum.p,
    cauchy = cauchy,
    combined.p = stats::pcauchy(cauchy, lower.tail = FALSE)
  )
}

zf_parzen_reference <- function(score, lag) {
  centered <- score - mean(score)
  lags <- if (lag == 1L) integer() else seq_len(lag - 1L)
  gamma <- vapply(lags, function(k) {
    sum(head(centered, -k) * tail(centered, -k)) / (length(score) - k)
  }, numeric(1L))
  ratio <- lags / lag
  weight <- ifelse(
    ratio < 0.5,
    1 - 6 * ratio^2 + 6 * ratio^3,
    2 * (1 - ratio)^3
  )
  list(
    gamma = gamma,
    weight = weight,
    contribution = gamma * weight,
    tau.squared = 1 + 2 * sum(gamma * weight)
  )
}

zf_example_one_sample <- function() {
  matrix(c(
    -4, -1,  2,  5,
    -2,  3, -5,  1,
     1, -4,  3, -6,
     3,  5, -1,  2,
     6, -2,  7, -3
  ), nrow = 5, byrow = TRUE)
}

zf_example_two_sample <- function() {
  matrix(c(
    -3.5, -0.6,  2.7,  5.6,
    -1.4,  3.6, -4.4,  1.7,
     1.8, -3.3,  4.1, -5.2,
     3.7,  5.8, -0.2,  2.9,
     6.5, -1.1,  7.9, -2.1,
     8.2,  7.1,  9.3,  4.4
  ), nrow = 6, byrow = TRUE)
}

zf_signed_rank_matrix <- function(targets, n = 5L) {
  vapply(targets, function(target) {
    remaining <- target
    positive <- logical(n)
    for (rank in n:1L) {
      if (rank <= remaining) {
        positive[[rank]] <- TRUE
        remaining <- remaining - rank
      }
    }
    stopifnot(remaining == 0)
    ifelse(positive, 1, -1) * seq_len(n)
  }, numeric(n))
}


test_that("one-sample components match literal signed-rank formulas", {
  x <- zf_example_one_sample()
  colnames(x) <- paste0("X", seq_len(ncol(x)))
  tau.squared <- 1.7
  reference <- zf_statistic_reference(
    zf_one_rank_reference(x), zf_one_moment_reference(nrow(x)), tau.squared
  )
  maximum <- zhang_feng_rank_one_sample_test(x)
  summed <- zhang_feng_rank_one_sample_test(
    x, component = "sum", tau_sq = tau.squared
  )
  combined <- zhang_feng_rank_one_sample_test(
    x, component = "combined", tau_sq = tau.squared
  )

  expect_equal(unname(maximum$components$rank.sum), reference$rank.sum)
  expect_equal(unname(maximum$components$centered.rank.score),
               reference$centered)
  expect_equal(unname(maximum$components$standardized.rank.score),
               reference$marginal, tolerance = 2e-15)
  expect_equal(unname(maximum$components$squared.rank.score),
               reference$square)
  expect_equal(maximum$components$null.moments,
               zf_one_moment_reference(nrow(x)))
  expect_equal(unname(maximum$statistic), reference$maximum,
               tolerance = 2e-15)
  expect_equal(maximum$p.value, reference$max.p, tolerance = 2e-15)
  expect_equal(unname(summed$statistic), reference$sum, tolerance = 2e-15)
  expect_equal(summed$p.value, reference$sum.p, tolerance = 2e-15)
  expect_equal(combined$components$combined$transform,
               reference$cauchy, tolerance = 2e-12)
  expect_equal(combined$p.value, reference$combined.p, tolerance = 2e-15)
  expect_identical(combined$diagnostics$tau.source, "supplied")
  expect_identical(combined$diagnostics$variance.repair, "none")
})


test_that("two-sample components match literal WMW formulas for unequal n", {
  x <- zf_example_one_sample()
  y <- zf_example_two_sample()
  tau.squared <- 0.85
  reference <- zf_statistic_reference(
    zf_two_rank_reference(x, y),
    zf_two_moment_reference(nrow(x), nrow(y)), tau.squared
  )
  maximum <- zhang_feng_rank_two_sample_test(x, y)
  summed <- zhang_feng_rank_two_sample_test(
    x, y, component = "sum", tau_sq = tau.squared
  )
  combined <- zhang_feng_rank_two_sample_test(
    x, y, component = "combined", tau_sq = tau.squared
  )

  expect_equal(unname(maximum$components$rank.sum), reference$rank.sum)
  expect_equal(maximum$components$null.moments,
               zf_two_moment_reference(nrow(x), nrow(y)))
  expect_equal(unname(maximum$statistic), reference$maximum,
               tolerance = 2e-15)
  expect_equal(maximum$p.value, reference$max.p, tolerance = 2e-15)
  expect_equal(unname(summed$statistic), reference$sum, tolerance = 2e-15)
  expect_equal(summed$p.value, reference$sum.p, tolerance = 2e-15)
  expect_equal(combined$components$combined$transform,
               reference$cauchy, tolerance = 2e-12)
  expect_equal(combined$p.value, reference$combined.p, tolerance = 2e-15)
  expect_equal(unname(maximum$n), c(5, 6))
})


test_that("published squared-score null moments agree with exact enumeration", {
  n <- 6L
  signed.rank.sum <- vapply(0:(2^n - 1L), function(mask) {
    included <- as.logical(intToBits(mask)[seq_len(n)])
    sum(seq_len(n)[included])
  }, numeric(1L))
  one.moments <- zf_one_moment_reference(n)
  one.square <- (signed.rank.sum - one.moments[["mean.rank.sum"]])^2
  expect_equal(mean(signed.rank.sum), one.moments[["mean.rank.sum"]])
  expect_equal(mean((signed.rank.sum - mean(signed.rank.sum))^2),
               one.moments[["variance.rank.sum"]])
  expect_equal(mean(one.square),
               one.moments[["mean.squared.rank.score"]])
  expect_equal(mean((one.square - mean(one.square))^2),
               one.moments[["variance.squared.rank.score"]])

  n1 <- 3L
  n2 <- 4L
  combinations <- utils::combn(n1 + n2, n1)
  mann.whitney <- colSums(combinations) - n1 * (n1 + 1) / 2
  two.moments <- zf_two_moment_reference(n1, n2)
  two.square <- (mann.whitney - two.moments[["mean.rank.sum"]])^2
  expect_equal(mean(mann.whitney), two.moments[["mean.rank.sum"]])
  expect_equal(mean((mann.whitney - mean(mann.whitney))^2),
               two.moments[["variance.rank.sum"]])
  expect_equal(mean(two.square),
               two.moments[["mean.squared.rank.score"]])
  expect_equal(mean((two.square - mean(two.square))^2),
               two.moments[["variance.squared.rank.score"]])
})


test_that("explicit Ouyang--Parzen estimator uses p-k denominators", {
  x <- zf_signed_rank_matrix(c(0, 1, 3, 6, 10, 15, 8, 4))
  preliminary <- zhang_feng_rank_one_sample_test(x)
  score <- unname(
    preliminary$components$standardized.squared.rank.score
  )
  reference <- zf_parzen_reference(score, lag = 5L)
  result <- zhang_feng_rank_one_sample_test(
    x, component = "sum", tau_method = "ouyang_parzen", lag = 5L
  )
  details <- result$components$long.run.variance

  expect_equal(details$lags, 1:4)
  expect_equal(details$autocovariance, reference$gamma,
               tolerance = 2e-15)
  expect_equal(details$weight, reference$weight, tolerance = 2e-15)
  expect_equal(details$weighted_autocovariance, reference$contribution,
               tolerance = 2e-15)
  expect_equal(details$tau_squared, reference$tau.squared,
               tolerance = 2e-15)
  expect_identical(details$denominator, "p - k")
  expect_identical(details$repair, "none")
  expect_true(result$diagnostics$coordinate.order.sensitive)

  lag.one <- zhang_feng_rank_one_sample_test(
    x, component = "sum", tau_method = "ouyang_parzen", lag = 1L
  )
  expect_equal(lag.one$components$long.run.variance$tau_squared, 1)
  expect_length(lag.one$components$long.run.variance$lags, 0L)
})


test_that("rank tests have their stated invariances and order contract", {
  x <- zf_example_one_sample()
  y <- zf_example_two_sample()
  shift <- c(16, -32, 8, 64)
  scale <- c(-3, 0.25, 7, -0.5)
  one <- zhang_feng_rank_one_sample_test(
    x, component = "combined", tau_sq = 1.2
  )
  one.shift <- zhang_feng_rank_one_sample_test(
    sweep(x, 2L, shift, "+"), mu = shift,
    component = "combined", tau_sq = 1.2
  )
  one.scale <- zhang_feng_rank_one_sample_test(
    sweep(x, 2L, scale, "*"),
    component = "combined", tau_sq = 1.2
  )
  one.rows <- zhang_feng_rank_one_sample_test(
    x[c(5, 2, 4, 1, 3), ], component = "combined", tau_sq = 1.2
  )
  expect_equal(one.shift$p.value, one$p.value, tolerance = 2e-15)
  expect_equal(one.scale$p.value, one$p.value, tolerance = 2e-15)
  expect_equal(one.rows$p.value, one$p.value, tolerance = 2e-15)

  two <- zhang_feng_rank_two_sample_test(
    x, y, component = "combined", tau_sq = 1.2
  )
  two.shift <- zhang_feng_rank_two_sample_test(
    sweep(x, 2L, shift, "+"), sweep(y, 2L, shift, "+"),
    component = "combined", tau_sq = 1.2
  )
  two.scale <- zhang_feng_rank_two_sample_test(
    sweep(x, 2L, scale, "*"), sweep(y, 2L, scale, "*"),
    component = "combined", tau_sq = 1.2
  )
  two.exchange <- zhang_feng_rank_two_sample_test(
    y, x, component = "combined", tau_sq = 1.2
  )
  expect_equal(two.shift$p.value, two$p.value, tolerance = 2e-15)
  expect_equal(two.scale$p.value, two$p.value, tolerance = 2e-15)
  expect_equal(two.exchange$p.value, two$p.value, tolerance = 2e-15)

  ordered <- zf_signed_rank_matrix(c(0, 1, 3, 6, 10, 15, 8, 4))
  permuted <- ordered[, c(1, 4, 8, 2, 6, 3, 7, 5)]
  supplied <- zhang_feng_rank_one_sample_test(
    ordered, component = "sum", tau_sq = 1
  )
  supplied.permuted <- zhang_feng_rank_one_sample_test(
    permuted, component = "sum", tau_sq = 1
  )
  parzen <- zhang_feng_rank_one_sample_test(
    ordered, component = "sum", tau_method = "ouyang_parzen", lag = 5
  )
  parzen.permuted <- zhang_feng_rank_one_sample_test(
    permuted, component = "sum", tau_method = "ouyang_parzen", lag = 5
  )
  expect_equal(supplied.permuted$p.value, supplied$p.value,
               tolerance = 2e-15)
  expect_false(isTRUE(all.equal(
    parzen.permuted$diagnostics$tau.squared,
    parzen$diagnostics$tau.squared,
    tolerance = 1e-12
  )))
})


test_that("rank kernels preserve extreme scales and overflowed residual order", {
  x <- zf_example_one_sample()
  ordinary <- zhang_feng_rank_one_sample_test(x)
  huge <- zhang_feng_rank_one_sample_test(x * 1e300)
  tiny <- zhang_feng_rank_one_sample_test(x * 1e-300)
  expect_equal(huge$statistic, ordinary$statistic, tolerance = 0)
  expect_equal(tiny$statistic, ordinary$statistic, tolerance = 0)

  extreme.x <- cbind(
    c(1.7e308, 1.6e308, 1.5e308),
    c(-1.7e308, -1.5e308, -1.3e308)
  )
  extreme.mu <- c(-1.7e308, 1.7e308)
  overflow <- zhang_feng_rank_one_sample_test(extreme.x, extreme.mu)
  scaled <- zhang_feng_rank_one_sample_test(
    extreme.x / 1.7e308, extreme.mu / 1.7e308
  )
  expect_equal(overflow$statistic, scaled$statistic, tolerance = 0)
  expect_equal(overflow$p.value, scaled$p.value, tolerance = 0)
  expect_identical(overflow$diagnostics$overflow.fallback.columns, 2L)

  y <- zf_example_two_sample()
  two <- zhang_feng_rank_two_sample_test(x, y)
  two.huge <- zhang_feng_rank_two_sample_test(x * 1e300, y * 1e300)
  two.tiny <- zhang_feng_rank_two_sample_test(x * 1e-300, y * 1e-300)
  expect_equal(two.huge$statistic, two$statistic, tolerance = 0)
  expect_equal(two.tiny$statistic, two$statistic, tolerance = 0)
})


test_that("stable Gumbel and Cauchy tails remain finite at underflow", {
  max.tail <- HDElliptical:::.zf_gumbel_tail(2000)
  sum.tail <- HDElliptical:::.zf_normal_tail(40)
  combined <- HDElliptical:::.zf_cauchy_combine(max.tail, sum.tail)
  expect_identical(max.tail$p.value, 0)
  expect_true(is.finite(max.tail$log.p.value))
  expect_identical(sum.tail$p.value, 0)
  expect_true(is.finite(sum.tail$log.p.value))
  expect_true(is.finite(combined$angle))
  expect_true(is.finite(combined$p.value))
  expect_gte(combined$p.value, 0)
  expect_lte(combined$p.value, 1)
  expect_true(combined$transform.sign > 0)

  moderate.max <- HDElliptical:::.zf_gumbel_tail(1.25)
  moderate.sum <- HDElliptical:::.zf_normal_tail(-0.35)
  moderate <- HDElliptical:::.zf_cauchy_combine(
    moderate.max, moderate.sum
  )
  direct <- 0.5 / tan(pi * moderate.max$p.value) +
    0.5 / tan(pi * moderate.sum$p.value)
  expect_equal(moderate$transform, direct, tolerance = 2e-15)
  expect_equal(moderate$p.value,
               stats::pcauchy(direct, lower.tail = FALSE),
               tolerance = 2e-15)
})


test_that("return objects expose the htest and diagnostic contract", {
  x <- zf_example_one_sample()
  colnames(x) <- c("a", "b", "c", "d")
  result <- zhang_feng_rank_one_sample_test(
    x, component = "combined", tau_sq = 1
  )
  expect_s3_class(result, "hd_location_test")
  expect_s3_class(result, "htest")
  expect_identical(names(result$estimate), colnames(x))
  expect_identical(names(result$null.value), colnames(x))
  expect_identical(result$diagnostics$selected.component, "combined")
  expect_false(result$diagnostics$automatic.bandwidth)
  expect_false(result$diagnostics$coordinate.order.sensitive)
  expect_false(result$diagnostics$book.general.Lq.minimum.p.implemented)
  expect_equal(result$components$combined$weights,
               c(max = 0.5, sum = 0.5))
  expect_match(result$null.distribution$family, "Cauchy")
  expect_match(result$diagnostics$marginal.null, "symmetry")

  two <- zhang_feng_rank_two_sample_test(
    x, zf_example_two_sample(), component = "sum", tau_sq = 1
  )
  expect_identical(unname(two$null.value), rep(0, ncol(x)))
  expect_match(two$diagnostics$marginal.null, "pure-shift")
})


test_that("continuous-data and explicit-tau contracts fail loudly", {
  x <- zf_example_one_sample()
  y <- zf_example_two_sample()
  zero <- x
  zero[1, 1] <- 0
  tied.absolute <- x
  tied.absolute[2, 1] <- -tied.absolute[1, 1]
  tied.pooled <- y
  tied.pooled[1, 1] <- x[1, 1]

  expect_error(zhang_feng_rank_one_sample_test(zero), "zero residual")
  expect_error(zhang_feng_rank_one_sample_test(tied.absolute),
               "absolute residuals")
  expect_error(zhang_feng_rank_two_sample_test(x, tied.pooled),
               "pooled ties|pooled values")
  expect_error(zhang_feng_rank_one_sample_test(x[, 1, drop = FALSE]),
               "p >= 2")
  expect_error(zhang_feng_rank_two_sample_test(
    x[, 1, drop = FALSE], y[, 1, drop = FALSE]
  ), "p >= 2")
  expect_error(zhang_feng_rank_two_sample_test(x, y[, -1]),
               "same number")
  x.named <- x
  y.named <- y
  colnames(x.named) <- letters[1:4]
  colnames(y.named) <- letters[4:1]
  expect_error(zhang_feng_rank_two_sample_test(x.named, y.named),
               "same names")
  expect_error(zhang_feng_rank_one_sample_test(x, component = "other"))
  expect_error(zhang_feng_rank_one_sample_test(x, tau_sq = 1),
               "not used")
  expect_error(zhang_feng_rank_one_sample_test(x, component = "sum"),
               "supply positive")
  expect_error(zhang_feng_rank_one_sample_test(
    x, component = "sum", tau_method = "supplied"
  ), "requires a positive")
  expect_error(zhang_feng_rank_one_sample_test(
    x, component = "sum", tau_sq = 1, lag = 2
  ), "only used")
  expect_error(zhang_feng_rank_one_sample_test(
    x, component = "sum", tau_sq = 1,
    tau_method = "ouyang_parzen", lag = 2
  ), "exactly one")
  for (bad in list(0, -1, Inf, NA_real_, c(1, 2), "1", TRUE)) {
    expect_error(zhang_feng_rank_one_sample_test(
      x, component = "sum", tau_sq = bad
    ), "strictly positive")
  }
  expect_error(zhang_feng_rank_one_sample_test(
    x, component = "sum", tau_method = "ouyang_parzen"
  ), "explicit `lag`")
  for (bad in list(0, 1.5, Inf, NA_real_, 5, c(1, 2))) {
    expect_error(zhang_feng_rank_one_sample_test(
      x, component = "sum", tau_method = "ouyang_parzen", lag = bad
    ), "integer in")
  }
  expect_error(zhang_feng_rank_one_sample_test(
    x, component = "sum", tau_method = "unknown", lag = 2
  ), "arg")

  nonpositive <- zf_signed_rank_matrix(c(0, 1, 3, 6, 10, 15, 8, 4))
  expect_error(zhang_feng_rank_one_sample_test(
    nonpositive, component = "sum",
    tau_method = "ouyang_parzen", lag = 7
  ), "not strictly positive")
})
