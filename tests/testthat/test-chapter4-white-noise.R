ch4wn_fixture <- function() {
  structure(
    matrix(c(
      -2.3,  0.4,  1.1,
      -1.4,  1.8, -0.7,
      -0.8, -1.9,  0.3,
       0.2,  0.9,  2.2,
       1.0, -0.6, -1.6,
       1.7,  1.3,  0.5,
       2.5, -1.2,  1.4,
       3.1,  0.1, -0.9,
       3.8,  2.4,  0.8
    ), ncol = 3L, byrow = TRUE),
    dimnames = list(paste0("t", 1:9), c("a", "b", "c"))
  )
}

ch4wn_ev_reference <- function(value) {
  -expm1(-exp(-value / 2) / sqrt(pi))
}

ch4wn_flm_reference <- function(x, lag) {
  x <- unname(x) / max(abs(x))
  n <- nrow(x)
  p <- ncol(x)
  second <- colSums(x^2) / n
  max.correlation <- -Inf
  max.lag <- NA_integer_
  max.index <- c(NA_integer_, NA_integer_)
  for (h in seq_len(lag)) {
    m <- n - h
    covariance <- crossprod(x[seq_len(m), , drop = FALSE],
                            x[h + seq_len(m), , drop = FALSE]) / n
    correlation <- covariance / sqrt(outer(second, second))
    location <- which(abs(correlation) == max(abs(correlation)),
                      arr.ind = TRUE)[1L, ]
    if (abs(correlation[location[1L], location[2L]]) > max.correlation) {
      max.correlation <- abs(correlation[location[1L], location[2L]])
      max.lag <- h
      max.index <- location
    }
  }
  gram <- tcrossprod(x)
  trace.numerator <- sum(gram[row(gram) != col(gram)]^2)
  sum.numerator <- 0
  lag.numerator <- numeric(lag)
  for (h in seq_len(lag)) {
    m <- n - h
    current <- 0
    for (t in seq_len(m)) for (s in seq_len(m)) {
      if (s != t) current <- current + gram[t, s] * gram[t + h, s + h]
    }
    lag.numerator[h] <- current
    sum.numerator <- sum.numerator + current
  }
  denominator <- n * (n - 1)
  trace.hat <- trace.numerator / denominator
  sum.statistic <- sum.numerator / denominator
  sigma2 <- 2 * lag / denominator * trace.hat^2
  sum.z <- sum.statistic / sqrt(sigma2)
  maximum <- sqrt(n) * max.correlation
  comparisons <- lag * p^2
  max.gumbel <- maximum^2 - 2 * log(comparisons) + log(log(comparisons))
  p.max <- ch4wn_ev_reference(max.gumbel)
  p.sum <- pnorm(sum.z, lower.tail = FALSE)
  fisher <- -2 * (log(p.max) + pnorm(sum.z, lower.tail = FALSE,
                                     log.p = TRUE))
  list(
    second = second,
    maximum = maximum,
    max.correlation = max.correlation,
    max.lag = max.lag,
    max.index = max.index,
    max.gumbel = max.gumbel,
    p.max = p.max,
    trace.numerator = trace.numerator,
    trace.hat = trace.hat,
    sum.numerator = sum.numerator,
    lag.numerator = lag.numerator,
    sum.statistic = sum.statistic,
    sigma2 = sigma2,
    sum.z = sum.z,
    p.sum = p.sum,
    fisher = fisher,
    p.fisher = pchisq(fisher, 4, lower.tail = FALSE),
    denominator = denominator
  )
}

ch4wn_unit <- function(value) {
  scale <- max(abs(value))
  if (scale == 0) return(numeric(length(value)))
  value <- value / scale
  value / sqrt(sum(value^2))
}

ch4wn_sign_reference <- function(x, lag) {
  x <- unname(x) / max(abs(x))
  n <- nrow(x)
  signs <- t(vapply(seq_len(n), function(i) ch4wn_unit(x[i, ]),
                    numeric(ncol(x))))
  gram <- tcrossprod(signs)
  trace.sum <- sum(gram[upper.tri(gram)]^2)
  trace.hat <- 2 * trace.sum / (n * (n - 1))
  lag.numerator <- numeric(lag)
  lag.statistic <- numeric(lag)
  for (h in seq_len(lag)) {
    for (s in (h + 1L):(n - 1L)) for (tt in (s + 1L):n) {
      lag.numerator[h] <- lag.numerator[h] +
        gram[s - h, tt - h] * gram[s, tt]
    }
    lag.statistic[h] <- lag.numerator[h] / (n - h)
  }
  statistic <- sum(lag.statistic)
  sigma2 <- lag / 2 * trace.hat^2
  z <- statistic / sqrt(sigma2)
  list(
    signs = signs,
    trace.sum = trace.sum,
    trace.hat = trace.hat,
    lag.numerator = lag.numerator,
    lag.statistic = lag.statistic,
    statistic = statistic,
    sigma2 = sigma2,
    z = z,
    p.value = pnorm(z, lower.tail = FALSE)
  )
}

ch4wn_rank_reference <- function(x, lag, measure) {
  x <- unname(x)
  n <- nrow(x)
  p <- ncol(x)
  maximum <- -Inf
  maximum.value <- NA_real_
  maximum.location <- c(NA_integer_, NA_integer_, NA_integer_)
  lag.maximum <- numeric(lag)
  for (h in seq_len(lag)) {
    m <- n - h
    current.maximum <- -Inf
    for (i in seq_len(p)) for (j in seq_len(p)) {
      first <- x[seq_len(m), i]
      second <- x[h + seq_len(m), j]
      if (measure == "spearman") {
        value <- cor(rank(first), rank(second))
        score <- m * value^2
      } else {
        concordance <- 0
        for (s in seq_len(m - 1L)) for (tt in (s + 1L):m) {
          concordance <- concordance +
            sign(first[tt] - first[s]) * sign(second[tt] - second[s])
        }
        value <- 2 * concordance / (m * (m - 1))
        score <- 9 * m * (m - 1) / (2 * (2 * m + 5)) * value^2
      }
      current.maximum <- max(current.maximum, score)
      if (score > maximum) {
        maximum <- score
        maximum.value <- value
        maximum.location <- c(h, i, j)
      }
    }
    lag.maximum[h] <- current.maximum
  }
  comparisons <- lag * p^2
  transformed <- maximum - 2 * log(comparisons) + log(log(comparisons))
  list(
    maximum = maximum,
    value = maximum.value,
    location = maximum.location,
    lag.maximum = lag.maximum,
    transformed = transformed,
    p.value = ch4wn_ev_reference(transformed)
  )
}

test_that("classical portmanteau statistics match literal formulas", {
  x <- c(2.1, -0.7, 1.3, 0.2, -1.9, 0.8, 2.7, -0.4)
  residuals <- x - mean(x)
  n <- length(x)
  lag <- 3L
  rho <- vapply(seq_len(lag), function(h) {
    sum(residuals[seq_len(n - h)] * residuals[h + seq_len(n - h)]) /
      sum(residuals^2)
  }, numeric(1L))
  bp.ref <- n * sum(rho^2)
  lb.ref <- n * (n + 2) * sum(rho^2 / (n - seq_len(lag)))

  bp <- white_noise_portmanteau_test(x, lag, "Box-Pierce")
  lb <- white_noise_portmanteau_test(x, lag, "Ljung-Box")
  expect_s3_class(bp, "htest")
  expect_s3_class(lb, "hd_white_noise_test")
  expect_equal(unname(bp$statistic), bp.ref, tolerance = 2e-15)
  expect_equal(unname(lb$statistic), lb.ref, tolerance = 2e-15)
  expect_equal(unname(bp$components$autocorrelations), rho,
               tolerance = 2e-15)
  expect_equal(bp$p.value, pchisq(bp.ref, lag, lower.tail = FALSE),
               tolerance = 1e-12)
  expect_equal(lb$p.value, pchisq(lb.ref, lag, lower.tail = FALSE),
               tolerance = 1e-12)
  expect_identical(unname(bp$parameter), lag)
})

test_that("portmanteau preprocessing invariance and degeneracy are explicit", {
  x <- c(-2, 0.3, 1.1, -0.8, 2.4, 0.7)
  base <- white_noise_portmanteau_test(x, lag = 2)
  shifted <- white_noise_portmanteau_test(x + 100, lag = 2)
  scaled <- white_noise_portmanteau_test(x * -1e250, lag = 2)
  expect_equal(base$p.value, shifted$p.value, tolerance = 2e-14)
  expect_equal(base$p.value, scaled$p.value, tolerance = 2e-14)
  expect_error(white_noise_portmanteau_test(rep(1, 6)), "zero")
  expect_error(white_noise_portmanteau_test(x, lag = 6), "between")
  expect_error(white_noise_portmanteau_test(c(x, NA_real_)), "finite")
})

test_that("FLM primary n-divisor and ordered-pair formulas are locked", {
  x <- ch4wn_fixture()
  lag <- 2L
  ref <- ch4wn_flm_reference(x, lag)
  result <- feng_liu_ma_white_noise_test(x, lag, keep_lag = TRUE)

  expect_s3_class(result, "htest")
  expect_s3_class(result, "hd_white_noise_test")
  expect_equal(result$components$T.MAX, ref$maximum, tolerance = 3e-13)
  expect_equal(result$components$maximum.absolute.correlation,
               ref$max.correlation, tolerance = 3e-13)
  expect_equal(result$components$maximum.Gumbel.statistic,
               ref$max.gumbel, tolerance = 4e-13)
  expect_equal(result$components$p.maximum, ref$p.max, tolerance = 3e-15)
  expect_equal(result$components$trace.numerator.scaled,
               ref$trace.numerator, tolerance = 3e-12)
  expect_equal(result$components$trace.Sigma.squared.scaled,
               ref$trace.hat, tolerance = 3e-14)
  expect_equal(result$components$sum.numerator.scaled,
               ref$sum.numerator, tolerance = 4e-12)
  expect_equal(result$components$T.SUM.scaled,
               ref$sum.statistic, tolerance = 4e-14)
  expect_equal(result$components$sigma.S.squared.scaled,
               ref$sigma2, tolerance = 4e-15)
  expect_equal(result$components$sum.z, ref$sum.z, tolerance = 4e-13)
  expect_equal(result$components$p.sum, ref$p.sum, tolerance = 3e-15)
  expect_equal(result$components$Fisher.statistic,
               ref$fisher, tolerance = 5e-13)
  expect_equal(result$p.value, ref$p.fisher, tolerance = 3e-15)
  expect_equal(result$components$primary.ordered.denominator,
               nrow(x) * (nrow(x) - 1), tolerance = 0)
  expect_equal(result$components$lag.sum.numerator.scaled,
               ref$lag.numerator, tolerance = 3e-12)
  expect_equal(result$components$lag.ordered.pair.count,
               (nrow(x) - seq_len(lag)) *
                 (nrow(x) - seq_len(lag) - 1), tolerance = 0)
  expect_identical(result$diagnostics$sample.autocovariance.divisor,
                   "n at every lag")
  expect_identical(result$diagnostics$sum.and.trace.denominator,
                   "ordered n * (n - 1)")
})

test_that("FLM component selection changes only the htest view", {
  x <- ch4wn_fixture()
  fisher <- feng_liu_ma_white_noise_test(x, 2, "fisher")
  sum <- feng_liu_ma_white_noise_test(x, 2, "sum")
  maximum <- feng_liu_ma_white_noise_test(x, 2, "max")
  expect_equal(fisher$p.value, fisher$components$p.Fisher, tolerance = 0)
  expect_equal(sum$p.value, sum$components$p.sum, tolerance = 0)
  expect_equal(maximum$p.value, maximum$components$p.maximum, tolerance = 0)
  expect_equal(unname(sum$statistic), sum$components$sum.z, tolerance = 0)
  expect_equal(unname(maximum$statistic),
               maximum$components$maximum.Gumbel.statistic, tolerance = 0)
  expect_equal(fisher$components$T.SUM.scaled,
               maximum$components$T.SUM.scaled, tolerance = 0)
  expect_equal(sum$components$T.MAX, fisher$components$T.MAX, tolerance = 0)
})

test_that("FLM invariances and non-invariance match the statistic", {
  x <- ch4wn_fixture()
  base <- feng_liu_ma_white_noise_test(x, 2)
  common <- feng_liu_ma_white_noise_test(x * -1e250, 2)
  permuted <- feng_liu_ma_white_noise_test(x[, c(3, 1, 2)], 2)
  centered <- feng_liu_ma_white_noise_test(x, 2, center = "mean")
  translated <- feng_liu_ma_white_noise_test(
    sweep(x, 2, c(100, -70, 25), "+"), 2, center = "mean"
  )
  expect_equal(base$p.value, common$p.value, tolerance = 2e-13)
  expect_equal(base$p.value, permuted$p.value, tolerance = 2e-13)
  expect_equal(centered$p.value, translated$p.value, tolerance = 2e-13)
  expect_false(centered$diagnostics$primary.mean.zero.preprocessing)

  scaled.columns <- sweep(x, 2, c(1, 6, 0.2), "*")
  max.base <- feng_liu_ma_white_noise_test(x, 2, "max")
  max.scaled <- feng_liu_ma_white_noise_test(scaled.columns, 2, "max")
  sum.base <- feng_liu_ma_white_noise_test(x, 2, "sum")
  sum.scaled <- feng_liu_ma_white_noise_test(scaled.columns, 2, "sum")
  expect_equal(max.base$p.value, max.scaled$p.value, tolerance = 2e-13)
  expect_gt(abs(sum.base$components$sum.z - sum.scaled$components$sum.z),
            1e-4)
})

test_that("FLM boundary failures do not repair invalid variances", {
  x <- ch4wn_fixture()
  constant <- x
  constant[, 2] <- 1
  expect_error(feng_liu_ma_white_noise_test(constant, center = "mean"),
               "strictly positive")
  expect_error(feng_liu_ma_white_noise_test(matrix(1:6, ncol = 1), lag = 1),
               "lag \\* p\\^2")
  expect_error(feng_liu_ma_white_noise_test(x, lag = nrow(x) - 1),
               "n - 2")
  expect_error(feng_liu_ma_white_noise_test(x, keep_lag = NA),
               "TRUE or FALSE")
  expect_error(feng_liu_ma_white_noise_test(x + Inf), "finite")
})

test_that("spatial-sign statistic and H-over-two variance match literal sums", {
  x <- ch4wn_fixture()
  lag <- 2L
  ref <- ch4wn_sign_reference(x, lag)
  result <- zhao_chen_wang_spatial_sign_white_noise_test(
    x, lag, keep_signs = TRUE
  )
  expect_s3_class(result, "htest")
  expect_equal(unname(result$components$signs), ref$signs,
               tolerance = 3e-15)
  expect_equal(result$components$trace.unordered.sum,
               ref$trace.sum, tolerance = 3e-13)
  expect_equal(result$components$trace.Omega.squared,
               ref$trace.hat, tolerance = 3e-15)
  expect_equal(result$components$lag.numerators,
               ref$lag.numerator, tolerance = 4e-13)
  expect_equal(result$components$lag.statistics,
               ref$lag.statistic, tolerance = 4e-14)
  expect_equal(result$components$T.S, ref$statistic, tolerance = 5e-14)
  expect_equal(result$components$sigma.S.squared,
               ref$sigma2, tolerance = 4e-15)
  expect_equal(unname(result$statistic), ref$z, tolerance = 4e-13)
  expect_equal(result$p.value, ref$p.value, tolerance = 3e-15)
  expect_identical(result$diagnostics$primary.variance.factor, "H / 2")
  expect_false(isTRUE(all.equal(
    result$components$sigma.S.squared,
    lag^2 * result$components$trace.Omega.squared^2
  )))
})

test_that("spatial-sign geometric invariances and scaling counterexample hold", {
  x <- ch4wn_fixture()
  q <- qr.Q(qr(matrix(c(1, 2, 3, -2, 1, 1, 1, -3, 2), 3, 3)))
  base <- zhao_chen_wang_spatial_sign_white_noise_test(x, 2)
  rotated <- zhao_chen_wang_spatial_sign_white_noise_test(x %*% q, 2)
  permuted <- zhao_chen_wang_spatial_sign_white_noise_test(
    x[, c(3, 1, 2)], 2
  )
  common <- zhao_chen_wang_spatial_sign_white_noise_test(x * -1e250, 2)
  centered <- zhao_chen_wang_spatial_sign_white_noise_test(
    x, 2, center = "mean"
  )
  translated <- zhao_chen_wang_spatial_sign_white_noise_test(
    sweep(x, 2, c(40, -80, 120), "+"), 2, center = "mean"
  )
  expect_equal(base$p.value, rotated$p.value, tolerance = 3e-13)
  expect_equal(base$p.value, permuted$p.value, tolerance = 3e-13)
  expect_equal(base$p.value, common$p.value, tolerance = 3e-13)
  expect_equal(centered$p.value, translated$p.value, tolerance = 3e-13)

  coordinate.scaled <- zhao_chen_wang_spatial_sign_white_noise_test(
    sweep(x, 2, c(1, 9, 0.15), "*"), 2
  )
  expect_gt(abs(base$components$z - coordinate.scaled$components$z), 1e-4)
})

test_that("spatial-sign zeros and degenerate inputs have explicit contracts", {
  x <- ch4wn_fixture()
  x[1, ] <- 0
  expect_error(zhao_chen_wang_spatial_sign_white_noise_test(x, 1),
               "exactly zero")
  kept <- zhao_chen_wang_spatial_sign_white_noise_test(
    x, 1, zero_action = "keep", keep_signs = TRUE
  )
  expect_identical(kept$diagnostics$zero.sign.count, 1L)
  expect_true(kept$diagnostics$partial.zero.calibration.warning)
  expect_equal(unname(kept$components$signs[1, ]), numeric(ncol(x)),
               tolerance = 0)
  expect_error(zhao_chen_wang_spatial_sign_white_noise_test(
    matrix(0, 5, 2), 1, zero_action = "keep"
  ), "identically zero")
  expect_error(zhao_chen_wang_spatial_sign_white_noise_test(x, lag = 8),
               "n - 2")
})

test_that("Spearman and Kendall rank-max tests match literal formulas", {
  x <- ch4wn_fixture()
  for (measure in c("spearman", "kendall")) {
    ref <- ch4wn_rank_reference(x, 2, measure)
    result <- chen_song_feng_rank_white_noise_test(
      x, 2, measure, keep_lag = TRUE
    )
    expect_s3_class(result, "htest")
    expect_equal(result$components$maximum.standardized.square,
                 ref$maximum, tolerance = 4e-13, info = measure)
    expect_equal(result$components$measure.at.maximum,
                 ref$value, tolerance = 3e-14, info = measure)
    expect_equal(unname(result$components$maximum.location),
                 ref$location, tolerance = 0, info = measure)
    expect_equal(result$components$lag.maximum.standardized.square,
                 ref$lag.maximum, tolerance = 4e-13, info = measure)
    expect_equal(result$components$Gumbel.statistic,
                 ref$transformed, tolerance = 4e-13, info = measure)
    expect_equal(result$p.value, ref$p.value, tolerance = 3e-15,
                 info = measure)
  }
})

test_that("published Spearman and Kendall normalizations are distinct and locked", {
  x <- ch4wn_fixture()
  spearman <- chen_song_feng_rank_white_noise_test(x, 2, "spearman")
  kendall <- chen_song_feng_rank_white_noise_test(x, 2, "kendall")
  ms <- nrow(x) - spearman$components$maximum.location[["lag"]]
  mk <- nrow(x) - kendall$components$maximum.location[["lag"]]
  expect_equal(spearman$components$primary.null.variance.at.maximum,
               1 / (ms - 1), tolerance = 0)
  expect_equal(kendall$components$primary.null.variance.at.maximum,
               2 * (2 * mk + 5) / (9 * mk * (mk - 1)), tolerance = 0)
  expect_identical(spearman$components$primary.scaling,
                   "(n - k) * rho_ij(k)^2")
  expect_match(kendall$components$primary.scaling, "2 \\* \\[2")
})

test_that("rank-max tests have coordinatewise monotone invariance", {
  x <- ch4wn_fixture()
  transformed <- cbind(
    exp(x[, 1] / 4),
    -x[, 2]^3,
    10 + 7 * x[, 3]
  )
  permuted <- x[, c(3, 1, 2)]
  reversed <- x[nrow(x):1, , drop = FALSE]
  for (measure in c("spearman", "kendall")) {
    base <- chen_song_feng_rank_white_noise_test(x, 2, measure)
    monotone <- chen_song_feng_rank_white_noise_test(
      transformed, 2, measure
    )
    coordinate <- chen_song_feng_rank_white_noise_test(
      permuted, 2, measure
    )
    time.reversed <- chen_song_feng_rank_white_noise_test(
      reversed, 2, measure
    )
    extreme <- chen_song_feng_rank_white_noise_test(
      x * 1e250, 2, measure
    )
    expect_equal(base$p.value, monotone$p.value, tolerance = 3e-13,
                 info = measure)
    expect_equal(base$p.value, coordinate$p.value, tolerance = 3e-13,
                 info = measure)
    expect_equal(base$p.value, time.reversed$p.value, tolerance = 3e-13,
                 info = measure)
    expect_equal(base$p.value, extreme$p.value, tolerance = 3e-13,
                 info = measure)
  }
})

test_that("rank scope and continuous-margin boundary are explicit", {
  x <- ch4wn_fixture()
  result <- chen_song_feng_rank_white_noise_test(x, 1, "spearman")
  expect_false(result$diagnostics$rank.sum.or.adaptive.test.implemented)
  expect_match(result$diagnostics$rank.sum.or.adaptive.reason, "future work")
  expect_length(result$diagnostics$review.only.degenerate.methods, 3L)
  expect_true(result$diagnostics$continuous.margin.calibration)
  tied <- x
  tied[2, 1] <- tied[1, 1]
  expect_error(chen_song_feng_rank_white_noise_test(tied, 1), "ties")
  expect_error(chen_song_feng_rank_white_noise_test(
    matrix(c(3, 1, 2, 4), ncol = 1), 1
  ), "lag \\* p\\^2")
  expect_error(chen_song_feng_rank_white_noise_test(x, lag = 8), "n - 2")
  expect_error(chen_song_feng_rank_white_noise_test(x, keep_lag = NA),
               "TRUE or FALSE")
})
