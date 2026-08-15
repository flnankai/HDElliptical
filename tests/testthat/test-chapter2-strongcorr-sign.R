zfsc_fixture <- function() {
  rbind(
    c(0.8, -1.1, 0.3),
    c(-0.5, 0.7, 1.4),
    c(1.6, 0.2, -0.8),
    c(-1.3, -0.4, 0.9),
    c(0.1, 1.5, -0.2),
    c(1.0, -0.7, 0.5),
    c(-0.9, 0.4, -1.2)
  )
}

zfsc_ref_sign <- function(x, center) {
  residual <- sweep(x, 2L, center, "-")
  radius <- sqrt(rowSums(residual^2))
  answer <- matrix(0, nrow(x), ncol(x))
  nonzero <- radius > 0
  answer[nonzero, ] <- residual[nonzero, , drop = FALSE] / radius[nonzero]
  answer
}

zfsc_ref_pair_sum <- function(signs) {
  answer <- 0
  if (nrow(signs) < 2L) return(answer)
  for (i in seq_len(nrow(signs) - 1L)) {
    for (j in seq.int(i + 1L, nrow(signs))) {
      answer <- answer + sum(signs[i, ] * signs[j, ])
    }
  }
  answer
}

zfsc_ref_bootstrap <- function(signs, multiplier) {
  apply(multiplier, 1L, function(draw) {
    weighted <- signs * draw
    zfsc_ref_pair_sum(weighted)
  })
}

zfsc_ref_rademacher <- function(B, n, seed) {
  set.seed(seed)
  matrix(
    ifelse(stats::runif(B * n) < 0.5, -1, 1),
    nrow = B, ncol = n, byrow = TRUE
  )
}

zfsc_ref_gaussian <- function(B, n, seed) {
  set.seed(seed)
  matrix(stats::rnorm(B * n), nrow = B, ncol = n, byrow = TRUE)
}


test_that("Rademacher calibration matches literal pair formulas", {
  x <- zfsc_fixture()
  colnames(x) <- c("a", "b", "c")
  rownames(x) <- paste0("row", seq_len(nrow(x)))
  mu <- c(a = 0.25, b = -0.15, c = 0.1)
  B <- 29L
  alpha <- 0.2
  seed <- 2601L
  result <- zhao_feng_strongcorr_sign_test(
    x, mu, alpha = alpha, multiplier = "rademacher",
    B = B, seed = seed, keep_bootstrap = TRUE
  )

  fitted <- spatial_median(x, tol = 1e-8, max_iter = 1000L, warn = FALSE)
  observed.signs <- zfsc_ref_sign(x, mu)
  fitted.signs <- zfsc_ref_sign(x, fitted)
  observed <- zfsc_ref_pair_sum(observed.signs)
  draws <- zfsc_ref_rademacher(B, nrow(x), seed)
  bootstrap <- zfsc_ref_bootstrap(fitted.signs, draws)
  critical.index <- ceiling((1 - alpha) * B)
  critical <- sort(bootstrap)[critical.index]
  exceedances <- sum(bootstrap >= observed)

  expect_s3_class(result, "hd_location_test")
  expect_s3_class(result, "htest")
  expect_equal(unname(result$statistic), observed, tolerance = 3e-14)
  expect_equal(unname(result$raw.statistic), observed, tolerance = 3e-14)
  expect_equal(result$estimate, fitted, tolerance = 3e-14,
               ignore_attr = TRUE)
  expect_equal(
    unname(result$components$observed.signs.null.centered),
    observed.signs, tolerance = 3e-14
  )
  expect_equal(
    unname(result$components$bootstrap.signs.median.centered),
    fitted.signs, tolerance = 3e-14
  )
  expect_equal(result$components$bootstrap.statistics.raw,
               bootstrap, tolerance = 6e-14)
  expect_equal(result$components$critical.order.index,
               critical.index, tolerance = 0)
  expect_equal(result$components$critical.value.raw,
               critical, tolerance = 6e-14)
  expect_identical(result$components$reject.paper.critical,
                   observed > critical)
  expect_equal(result$components$exceedances.including.ties,
               exceedances, tolerance = 0)
  expect_equal(result$p.value, (1 + exceedances) / (B + 1), tolerance = 0)
  expect_identical(result$components$reject.plus.one.p.value,
                   result$p.value <= alpha)
  expect_equal(
    unname(result$components$bootstrap.summary.raw),
    unname(c(
      mean = mean(bootstrap),
      variance.population = mean((bootstrap - mean(bootstrap))^2),
      minimum = min(bootstrap), maximum = max(bootstrap)
    )),
    tolerance = 8e-14
  )
  expect_identical(result$alternative, "two.sided")
  expect_identical(result$null.distribution$tail, "upper")
})


test_that("null and fitted centering roles cannot be interchanged", {
  x <- zfsc_fixture()
  mu <- c(0.7, -0.6, 0.4)
  result <- zhao_feng_strongcorr_sign_test(
    x, mu, B = 17L, seed = 81L, keep_bootstrap = TRUE
  )
  null.statistic <- zfsc_ref_pair_sum(zfsc_ref_sign(x, mu))
  fitted.statistic <- zfsc_ref_pair_sum(
    zfsc_ref_sign(x, result$components$fitted.spatial.median)
  )
  wrong.draws <- zfsc_ref_rademacher(17L, nrow(x), 81L)
  wrong.bootstrap <- zfsc_ref_bootstrap(
    zfsc_ref_sign(x, mu), wrong.draws
  )

  expect_equal(unname(result$statistic), null.statistic, tolerance = 3e-14)
  expect_equal(result$components$fitted.raw.pair.sum,
               fitted.statistic, tolerance = 3e-14)
  expect_gt(abs(null.statistic - fitted.statistic), 0.1)
  expect_false(isTRUE(all.equal(
    result$components$bootstrap.statistics.raw,
    wrong.bootstrap, tolerance = 1e-12
  )))
  expect_match(result$diagnostics$centering$observed, "null location")
  expect_match(result$diagnostics$centering$bootstrap,
               "ordinary spatial median")
})


test_that("Gaussian bootstrap uses the variable squared-multiplier diagonal", {
  x <- zfsc_fixture()
  mu <- c(-0.1, 0.2, -0.3)
  B <- 23L
  seed <- 902L
  result <- zhao_feng_strongcorr_sign_test(
    x, mu, alpha = 0.13, multiplier = "gaussian",
    B = B, seed = seed, keep_bootstrap = TRUE
  )
  fitted.signs <- unname(result$components$bootstrap.signs.median.centered)
  draws <- zfsc_ref_gaussian(B, nrow(x), seed)
  direct <- zfsc_ref_bootstrap(fitted.signs, draws)
  weighted.sums <- draws %*% fitted.signs
  variable.diagonal <- as.numeric(
    (draws^2) %*% rowSums(fitted.signs^2)
  )
  identity <- 0.5 * (rowSums(weighted.sums^2) - variable.diagonal)
  wrong.fixed.diagonal <- 0.5 * (
    rowSums(weighted.sums^2) - sum(rowSums(fitted.signs^2))
  )

  expect_equal(result$components$bootstrap.statistics.raw,
               direct, tolerance = 8e-14)
  expect_equal(result$components$bootstrap.statistics.raw,
               identity, tolerance = 8e-14)
  expect_gt(max(abs(direct - wrong.fixed.diagonal)), 0.1)
  expect_identical(result$components$multiplier, "gaussian")
})


test_that("tau and pair normalizations cancel without a plug-in", {
  x <- zfsc_fixture()
  result <- zhao_feng_strongcorr_sign_test(
    x, B = 31L, seed = 44L, keep_bootstrap = TRUE
  )
  component <- result$components
  root.pairs <- sqrt(choose(nrow(x), 2))
  arbitrary.tau <- 0.037

  expect_equal(component$root.pair.count, root.pairs, tolerance = 0)
  expect_equal(
    component$observed.pair.normalized.tau.free,
    component$observed.raw / root.pairs, tolerance = 3e-15
  )
  expect_equal(
    component$bootstrap.statistics.pair.normalized.tau.free,
    component$bootstrap.statistics.raw / root.pairs, tolerance = 3e-15
  )
  expect_equal(
    component$critical.value.pair.normalized.tau.free,
    component$critical.value.raw / root.pairs, tolerance = 3e-15
  )
  expect_identical(
    component$observed.raw > component$critical.value.raw,
    component$observed.raw / sqrt(arbitrary.tau) / root.pairs >
      component$critical.value.raw / sqrt(arbitrary.tau) / root.pairs
  )
  expect_identical(component$tau.estimated, FALSE)
  expect_identical(component$tau.factor.used, FALSE)
  expect_identical(result$diagnostics$scale$tau.estimated, FALSE)
  expect_identical(result$diagnostics$regularization, "none")
  expect_identical(result$diagnostics$numerical.floor, "none")
  expect_identical(result$diagnostics$perturbation, "none")
})


test_that("paper critical and auxiliary plus-one decisions are both literal", {
  x <- zfsc_fixture()
  B <- 9L
  alpha <- 0.3
  result <- zhao_feng_strongcorr_sign_test(
    x, alpha = alpha, B = B, seed = 991L, keep_bootstrap = TRUE
  )
  bootstrap <- result$components$bootstrap.statistics.raw
  observed <- unname(result$statistic)
  type1 <- as.numeric(stats::quantile(
    bootstrap, probs = 1 - alpha, type = 1, names = FALSE
  ))

  expect_equal(result$components$critical.order.index,
               ceiling((1 - alpha) * B), tolerance = 0)
  expect_equal(result$components$critical.value.raw,
               type1, tolerance = 0)
  expect_identical(result$components$reject.paper.critical,
                   observed > type1)
  expect_equal(
    result$p.value,
    (1 + sum(bootstrap >= observed)) / (B + 1), tolerance = 0
  )
  expect_match(result$diagnostics$rejection$paper.rule, "strict")
  expect_match(result$diagnostics$rejection$auxiliary.p.value.rule, ">=")
  expect_true(result$diagnostics$rejection$finite.B.decisions.can.differ)

  tied <- matrix(rep(c(2, -1, 3), each = 5L), 5, 3)
  tied.result <- zhao_feng_strongcorr_sign_test(
    tied, mu = c(2, -1, 3), alpha = 0.2,
    B = 7L, seed = 1L, keep_bootstrap = TRUE
  )
  expect_true(all(tied.result$components$bootstrap.statistics.raw == 0))
  expect_equal(unname(tied.result$statistic), 0, tolerance = 0)
  expect_equal(tied.result$components$critical.value.raw, 0, tolerance = 0)
  expect_false(tied.result$components$reject.paper.critical)
  expect_equal(tied.result$components$exceedances.including.ties, 7,
               tolerance = 0)
  expect_equal(tied.result$p.value, 1, tolerance = 0)
  expect_true(tied.result$diagnostics$randomization$degenerate)
})


test_that("explicit and generated seeds obey the recorded RNG contract", {
  x <- zfsc_fixture()
  set.seed(727)
  state <- .Random.seed
  first <- zhao_feng_strongcorr_sign_test(
    x, B = 13L, seed = 31415L, keep_bootstrap = TRUE
  )
  expect_identical(.Random.seed, state)
  second <- zhao_feng_strongcorr_sign_test(
    x, B = 13L, seed = 31415L, keep_bootstrap = TRUE
  )
  expect_identical(first$components$bootstrap.statistics.raw,
                   second$components$bootstrap.statistics.raw)
  expect_identical(first$p.value, second$p.value)
  expect_equal(first$components$seed.used, 31415L, tolerance = 0)

  set.seed(812)
  expected.seed <- as.integer(
    sample.int(.Machine$integer.max, 1L) - 1L
  )
  expected.state <- .Random.seed
  set.seed(812)
  generated <- zhao_feng_strongcorr_sign_test(
    x, B = 11L, seed = NULL, keep_bootstrap = TRUE
  )
  expect_equal(generated$components$seed.used,
               expected.seed, tolerance = 0)
  expect_identical(.Random.seed, expected.state)

  set.seed(812)
  repeated <- zhao_feng_strongcorr_sign_test(
    x, B = 11L, seed = NULL, keep_bootstrap = TRUE
  )
  expect_identical(
    generated$components$bootstrap.statistics.raw,
    repeated$components$bootstrap.statistics.raw
  )
})


test_that("translation orthogonal and common-scale invariances hold", {
  x <- zfsc_fixture()
  mu <- c(0.2, -0.1, 0.3)
  baseline <- zhao_feng_strongcorr_sign_test(
    x, mu, B = 37L, seed = 606L, keep_bootstrap = TRUE,
    tol = 1e-12
  )
  shift <- c(200, -350, 125)
  translated <- zhao_feng_strongcorr_sign_test(
    sweep(x, 2L, shift, "+"), mu + shift,
    B = 37L, seed = 606L, keep_bootstrap = TRUE, tol = 1e-12
  )
  theta <- 0.43
  rotation <- matrix(c(
    cos(theta), -sin(theta), 0,
    sin(theta), cos(theta), 0,
    0, 0, 1
  ), 3, 3, byrow = TRUE)
  rotated <- zhao_feng_strongcorr_sign_test(
    x %*% rotation, as.numeric(mu %*% rotation),
    B = 37L, seed = 606L, keep_bootstrap = TRUE, tol = 1e-12
  )
  scaled <- zhao_feng_strongcorr_sign_test(
    -4.25 * x, -4.25 * mu,
    B = 37L, seed = 606L, keep_bootstrap = TRUE, tol = 1e-12
  )

  for (candidate in list(translated, rotated, scaled)) {
    expect_equal(unname(candidate$statistic), unname(baseline$statistic),
                 tolerance = 2e-10)
    expect_equal(candidate$components$bootstrap.statistics.raw,
                 baseline$components$bootstrap.statistics.raw,
                 tolerance = 2e-10)
    expect_equal(candidate$p.value, baseline$p.value, tolerance = 0)
  }
})


test_that("zero signs use the generalized diagonal identity", {
  x <- rbind(
    c(0, 0), c(1, 0), c(0, 1), c(-1, 0), c(0, -2), c(2, 1)
  )
  mu <- c(0, 0)
  result <- zhao_feng_strongcorr_sign_test(
    x, mu, B = 19L, seed = 5L, keep_bootstrap = TRUE
  )
  signs <- zfsc_ref_sign(x, mu)
  generalized <- 0.5 * (
    sum(colSums(signs)^2) - sum(rowSums(signs^2))
  )
  wrong.unit.diagonal <- 0.5 * (sum(colSums(signs)^2) - nrow(x))

  expect_equal(unname(result$statistic), generalized, tolerance = 3e-14)
  expect_equal(result$components$observed.sign.diagonal,
               sum(rowSums(signs^2)), tolerance = 3e-14)
  expect_equal(result$diagnostics$observed.zero.residuals, 1,
               tolerance = 0)
  expect_equal(generalized - wrong.unit.diagonal, 0.5, tolerance = 3e-14)
  expect_match(result$diagnostics$pair.sum.identity, "zero signs")
})


test_that("strict and non-strict median contracts expose convergence", {
  x <- zfsc_fixture()
  expect_error(
    zhao_feng_strongcorr_sign_test(
      x, B = 3L, seed = 1L, tol = 1e-16, max_iter = 1L
    ),
    "did not converge"
  )
  expect_warning(
    last <- zhao_feng_strongcorr_sign_test(
      x, B = 3L, seed = 1L, tol = 1e-16,
      max_iter = 1L, strict = FALSE
    ),
    "last finite iterate"
  )
  expect_false(last$diagnostics$spatial.median$converged)
  expect_identical(last$diagnostics$spatial.median$strict, FALSE)
  expect_true(is.finite(last$diagnostics$spatial.median$equation.residual))
  expect_true(is.finite(last$diagnostics$spatial.median$relative.change))
  expect_true(is.finite(last$statistic))
  expect_true(is.finite(last$p.value))
})


test_that("extreme residual subtraction is finite and never repaired", {
  maximum <- 0.95 * .Machine$double.xmax
  x <- rbind(
    c(maximum, 0), c(-maximum, 0),
    c(maximum, 1), c(-maximum, -1)
  )
  direct <- cpp_zhao_feng_strongcorr_sign_bootstrap(
    x, c(-maximum, 0), c(0, 0), "rademacher",
    7L, 0.05, TRUE
  )
  expect_true(is.infinite(x[1, 1] - (-maximum)))
  expect_true(is.finite(direct$observed_raw))
  expect_true(all(is.finite(direct$observed_signs)))
  expect_true(direct$observed_overflow_fallback_rows %in% c(0, 2))
  expect_true(all(is.finite(direct$bootstrap_statistics_raw)))

  ordinary <- zfsc_fixture()
  mu <- c(0.2, -0.1, 0.3)
  baseline <- zhao_feng_strongcorr_sign_test(
    ordinary, mu, B = 17L, seed = 90L, keep_bootstrap = TRUE
  )
  large <- zhao_feng_strongcorr_sign_test(
    ordinary * 1e150, mu * 1e150,
    B = 17L, seed = 90L, keep_bootstrap = TRUE
  )
  small <- zhao_feng_strongcorr_sign_test(
    ordinary * 1e-150, mu * 1e-150,
    B = 17L, seed = 90L, keep_bootstrap = TRUE
  )
  expect_equal(unname(large$statistic), unname(baseline$statistic),
               tolerance = 3e-12)
  expect_equal(unname(small$statistic), unname(baseline$statistic),
               tolerance = 3e-12)
  expect_equal(large$components$bootstrap.statistics.raw,
               baseline$components$bootstrap.statistics.raw,
               tolerance = 4e-12)
  expect_equal(small$components$bootstrap.statistics.raw,
               baseline$components$bootstrap.statistics.raw,
               tolerance = 4e-12)
  expect_identical(baseline$diagnostics$absolute.value.repair, "none")
  expect_identical(baseline$diagnostics$pseudoinverse, "none")
})


test_that("compact output and invalid boundaries are explicit", {
  x <- zfsc_fixture()
  compact <- zhao_feng_strongcorr_sign_test(
    x, B = 5L, seed = 1L, keep_bootstrap = FALSE
  )
  expect_null(compact$components$bootstrap.statistics.raw)
  expect_null(
    compact$components$bootstrap.statistics.pair.normalized.tau.free
  )
  expect_error(zhao_feng_strongcorr_sign_test(x[1, , drop = FALSE]),
               "at least 2 row")
  bad <- x
  bad[1, 1] <- NA_real_
  expect_error(zhao_feng_strongcorr_sign_test(bad), "only finite")
  expect_error(zhao_feng_strongcorr_sign_test(x, mu = 1:2), "length 3")
  expect_error(zhao_feng_strongcorr_sign_test(x, alpha = 0),
               "strictly between")
  expect_error(zhao_feng_strongcorr_sign_test(x, alpha = 1),
               "strictly between")
  expect_error(zhao_feng_strongcorr_sign_test(x, multiplier = "mammen"),
               "arg")
  expect_error(zhao_feng_strongcorr_sign_test(x, B = 0),
               "positive integer")
  expect_error(zhao_feng_strongcorr_sign_test(x, B = 1.5),
               "positive integer")
  expect_error(zhao_feng_strongcorr_sign_test(x, seed = -1),
               "zero through")
  expect_error(zhao_feng_strongcorr_sign_test(
    x, seed = .Machine$integer.max + 1
  ), "zero through")
  expect_error(zhao_feng_strongcorr_sign_test(x, seed = 1.5),
               "zero through")
  expect_error(zhao_feng_strongcorr_sign_test(x, keep_bootstrap = NA),
               "TRUE or FALSE")
  expect_error(zhao_feng_strongcorr_sign_test(x, strict = NA),
               "TRUE or FALSE")
  expect_error(zhao_feng_strongcorr_sign_test(x, tol = 0),
               "positive")
  expect_error(zhao_feng_strongcorr_sign_test(x, max_iter = 0),
               "positive integer")
})
