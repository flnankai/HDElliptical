clx_reference <- function(x, y, precision, oracle) {
  n1 <- nrow(x)
  n2 <- nrow(y)
  total <- n1 + n2
  effective.n <- n1 * n2 / total
  mean.x <- colMeans(x)
  mean.y <- colMeans(y)
  difference <- mean.x - mean.y
  transformed.difference <- as.numeric(precision %*% difference)

  if (oracle) {
    variance.x <- variance.y <- NULL
    denominator <- diag(precision)
  } else {
    transformed.x <- x %*% t(precision)
    transformed.y <- y %*% t(precision)
    variance.x <- colMeans(sweep(transformed.x, 2, colMeans(transformed.x))^2)
    variance.y <- colMeans(sweep(transformed.y, 2, colMeans(transformed.y))^2)
    denominator <- (n1 * variance.x + n2 * variance.y) / total
  }
  scores <- sqrt(effective.n) * transformed.difference / sqrt(denominator)
  squared.scores <- scores^2
  maximum <- max(squared.scores)
  centered.maximum <- maximum - 2 * log(ncol(x)) + log(log(ncol(x)))

  list(
    M = maximum,
    G = centered.maximum,
    argmax = which.max(squared.scores),
    scores = scores,
    squared.scores = squared.scores,
    denominator = denominator,
    variance.x = variance.x,
    variance.y = variance.y,
    mean.x = mean.x,
    mean.y = mean.y,
    difference = difference,
    transformed.difference = transformed.difference,
    effective.n = effective.n
  )
}

clx_adaptive_reference <- function(x, y, delta, eigen.floor) {
  n1 <- nrow(x)
  n2 <- nrow(y)
  total <- n1 + n2
  p <- ncol(x)
  centered.x <- sweep(x, 2, colMeans(x))
  centered.y <- sweep(y, 2, colMeans(y))
  pooled <- (crossprod(centered.x) + crossprod(centered.y)) / total
  theta <- lambda <- matrix(0, p, p)
  for (i in seq_len(p)) {
    for (j in seq_len(p)) {
      products <- c(
        centered.x[, i] * centered.x[, j],
        centered.y[, i] * centered.y[, j]
      )
      theta[i, j] <- sum((products - pooled[i, j])^2) / total
      lambda[i, j] <- delta * sqrt(theta[i, j] * log(p) / total)
    }
  }
  thresholded <- pooled * (abs(pooled) >= lambda)
  decomposition <- eigen(thresholded, symmetric = TRUE)
  scale <- max(abs(decomposition$values), abs(diag(pooled)))
  absolute.floor <- eigen.floor * scale
  adjusted.values <- decomposition$values
  adjusted.values[adjusted.values < absolute.floor] <- absolute.floor
  adjusted <- decomposition$vectors %*% (adjusted.values *
    t(decomposition$vectors))
  precision <- decomposition$vectors %*% ((1 / adjusted.values) *
    t(decomposition$vectors))

  list(
    pooled = pooled,
    theta = theta,
    lambda = lambda,
    thresholded = thresholded,
    adjusted = (adjusted + t(adjusted)) / 2,
    precision = (precision + t(precision)) / 2,
    eigenvalues.before = decomposition$values,
    eigenvalues.after = adjusted.values,
    absolute.floor = absolute.floor
  )
}

clx_stable_reference <- function(x, y, precision, oracle) {
  n1 <- nrow(x)
  n2 <- nrow(y)
  total <- n1 + n2
  anchor <- x[1, ]
  delta.x <- sweep(x, 2, anchor, `-`)
  delta.y <- sweep(y, 2, anchor, `-`)
  mean.delta.x <- colMeans(delta.x)
  mean.delta.y <- colMeans(delta.y)
  centered.x <- sweep(delta.x, 2, mean.delta.x, `-`)
  centered.y <- sweep(delta.y, 2, mean.delta.y, `-`)
  difference <- mean.delta.x - mean.delta.y
  transformed.difference <- as.numeric(precision %*% difference)

  if (oracle) {
    variance.x <- variance.y <- NULL
    denominator <- diag(precision)
  } else {
    transformed.x <- centered.x %*% t(precision)
    transformed.y <- centered.y %*% t(precision)
    variance.x <- colMeans(transformed.x^2)
    variance.y <- colMeans(transformed.y^2)
    denominator <- (n1 * variance.x + n2 * variance.y) / total
  }
  squared.scores <- (n1 * n2 / total) *
    transformed.difference^2 / denominator
  maximum <- max(squared.scores)

  list(
    M = maximum,
    G = maximum - 2 * log(ncol(x)) + log(log(ncol(x))),
    difference = difference,
    denominator = denominator,
    variance.x = variance.x,
    variance.y = variance.y
  )
}

clx_adaptive_stable_reference <- function(x, y, delta, eigen.floor) {
  total <- nrow(x) + nrow(y)
  p <- ncol(x)
  center.relative.to.first.row <- function(z) {
    delta.z <- sweep(z, 2, z[1, ], `-`)
    sweep(delta.z, 2, colMeans(delta.z), `-`)
  }
  centered.x <- center.relative.to.first.row(x)
  centered.y <- center.relative.to.first.row(y)
  internal.scale <- max(abs(centered.x), abs(centered.y))
  if (internal.scale == 0) {
    internal.scale <- 1
  }
  scaled.x <- centered.x / internal.scale
  scaled.y <- centered.y / internal.scale
  pooled.scaled <- (crossprod(scaled.x) + crossprod(scaled.y)) / total
  theta.scaled <- lambda.scaled <- matrix(0, p, p)
  for (i in seq_len(p)) {
    for (j in seq_len(p)) {
      products <- c(
        scaled.x[, i] * scaled.x[, j],
        scaled.y[, i] * scaled.y[, j]
      )
      theta.scaled[i, j] <-
        sum((products - pooled.scaled[i, j])^2) / total
      lambda.scaled[i, j] <- delta *
        sqrt(theta.scaled[i, j] * log(p) / total)
    }
  }
  thresholded.scaled <- pooled.scaled *
    (abs(pooled.scaled) >= lambda.scaled)
  decomposition <- eigen(thresholded.scaled, symmetric = TRUE)
  spectral.scale <- max(
    abs(decomposition$values), abs(diag(pooled.scaled))
  )
  absolute.floor <- eigen.floor * spectral.scale
  adjusted.values <- pmax(decomposition$values, absolute.floor)
  precision.scaled <- decomposition$vectors %*%
    ((1 / adjusted.values) * t(decomposition$vectors))
  precision.scaled <- (precision.scaled + t(precision.scaled)) / 2

  list(
    precision = precision.scaled / internal.scale^2,
    pooled.scaled = pooled.scaled,
    theta.scaled = theta.scaled,
    lambda.scaled = lambda.scaled,
    internal.scale = internal.scale
  )
}

test_that("oracle CLX agrees with a direct equation (2) calculation", {
  x <- matrix(
    c(
      2, 1, 3, 0,
      1, 4, 2, 2,
      3, 2, 5, 1,
      4, 3, 1, 3,
      2, 5, 4, 2,
      5, 2, 3, 4
    ),
    byrow = TRUE,
    ncol = 4,
    dimnames = list(NULL, c("geneA", "geneB", "geneC", "geneD"))
  )
  y <- matrix(
    c(
      0, 2, 1, 3,
      2, 1, 3, 1,
      1, 3, 2, 4,
      3, 0, 4, 2,
      2, 4, 0, 5
    ),
    byrow = TRUE,
    ncol = 4,
    dimnames = list(NULL, colnames(x))
  )
  precision <- matrix(
    c(
      1.5, 0.2, 0.0, 0.1,
      0.2, 1.2, 0.15, 0.0,
      0.0, 0.15, 1.0, 0.25,
      0.1, 0.0, 0.25, 1.4
    ),
    4, 4, byrow = TRUE
  )
  reference <- clx_reference(x, y, precision, oracle = TRUE)
  fit <- cai_liu_xia_two_sample_test(
    x, y, precision = precision, precision_source = "oracle"
  )

  expect_s3_class(fit, "hd_location_test")
  expect_s3_class(fit, "htest")
  expect_equal(unname(fit$raw.statistic), reference$M, tolerance = 1e-13)
  expect_equal(unname(fit$statistic), reference$G, tolerance = 1e-13)
  expect_equal(
    fit$p.value,
    1 - exp(-pi^(-0.5) * exp(-reference$G / 2)),
    tolerance = 1e-15
  )
  expect_equal(
    unname(fit$components$transformed.difference),
    reference$transformed.difference,
    tolerance = 1e-13
  )
  expect_equal(
    unname(fit$components$standardized.scores),
    reference$scores,
    tolerance = 1e-13
  )
  expect_equal(unname(fit$variance), reference$denominator)
  expect_equal(names(fit$variance), colnames(x))
  expect_equal(fit$components$maximum.coordinate$index, reference$argmax)
  expect_identical(
    fit$components$maximum.coordinate$name,
    colnames(x)[[reference$argmax]]
  )
  expect_identical(fit$diagnostics$precision.source, "oracle")
  expect_identical(fit$diagnostics$statistic.path, "oracle")
  expect_match(fit$diagnostics$denominator.source, "known population")
  expect_false("transformed.variance.x" %in% names(fit$components))
  expect_identical(dimnames(fit$components$precision),
                   list(colnames(x), colnames(x)))
})

test_that("supplied feasible CLX uses transformed empirical variances", {
  set.seed(2201)
  x <- matrix(stats::rnorm(14 * 5), 14, 5)
  y <- matrix(stats::rnorm(17 * 5, 0.25), 17, 5)
  base <- matrix(c(
    1, .2, 0, 0, .1,
    .2, 1.4, .1, 0, 0,
    0, .1, 1.2, .15, 0,
    0, 0, .15, 1.1, .2,
    .1, 0, 0, .2, 1.3
  ), 5, 5, byrow = TRUE)
  reference <- clx_reference(x, y, base, oracle = FALSE)
  feasible <- cai_liu_xia_two_sample_test(
    x, y, precision = base, precision_source = "feasible"
  )
  oracle <- cai_liu_xia_two_sample_test(
    x, y, precision = base, precision_source = "oracle"
  )

  expect_equal(feasible$raw.statistic, c(M = reference$M),
               tolerance = 2e-13)
  expect_equal(feasible$statistic, c(G = reference$G), tolerance = 2e-13)
  expect_equal(unname(feasible$components$transformed.variance.x),
               reference$variance.x, tolerance = 2e-13)
  expect_equal(unname(feasible$components$transformed.variance.y),
               reference$variance.y, tolerance = 2e-13)
  expect_equal(unname(feasible$components$denominator),
               reference$denominator, tolerance = 2e-13)
  expect_false(isTRUE(all.equal(reference$denominator, diag(base))))
  expect_false(isTRUE(all.equal(feasible$raw.statistic,
                                oracle$raw.statistic)))
  expect_identical(feasible$diagnostics$statistic.path, "feasible")
  expect_equal(
    feasible$diagnostics$transformed.variance.divisors,
    c(x = nrow(x), y = nrow(y))
  )
  expect_match(feasible$diagnostics$denominator.source, "divisors n1 and n2")
})

test_that("adaptive covariance, theta, and lambda reproduce the paper", {
  set.seed(2202)
  x <- matrix(stats::rnorm(55 * 4), 55, 4)
  y <- matrix(stats::rnorm(47 * 4, 0.15), 47, 4)
  delta <- 0.65
  eigen.floor <- 1e-10
  reference <- clx_adaptive_reference(x, y, delta, eigen.floor)
  backend <- cpp_clx_adaptive_precision(x, y, delta, eigen.floor)

  expect_equal(backend$pooled_covariance, reference$pooled,
               tolerance = 2e-13)
  expect_equal(backend$theta, reference$theta, tolerance = 2e-13)
  expect_equal(backend$thresholds, reference$lambda, tolerance = 2e-13)
  expect_equal(backend$thresholded_covariance, reference$thresholded,
               tolerance = 2e-13)
  expect_equal(backend$adjusted_covariance, reference$adjusted,
               tolerance = 2e-12)
  expect_equal(backend$precision, reference$precision, tolerance = 2e-12)
  expect_equal(backend$sample_size, nrow(x) + nrow(y))

  fit <- cai_liu_xia_two_sample_test(
    x, y, precision_source = "adaptive",
    threshold_delta = delta, eigen_floor = eigen.floor
  )
  statistic.reference <- clx_reference(
    x, y, reference$precision, oracle = FALSE
  )
  expect_equal(unname(fit$components$precision), reference$precision,
               tolerance = 2e-12)
  expect_equal(unname(fit$raw.statistic), statistic.reference$M,
               tolerance = 3e-12)
  expect_equal(unname(fit$components$denominator),
               statistic.reference$denominator, tolerance = 3e-12)
  expect_false(isTRUE(all.equal(
    unname(fit$components$denominator),
    diag(fit$components$precision)
  )))
  adaptive <- fit$diagnostics$adaptive.thresholding
  expect_equal(adaptive$pooled.covariance.denominator, nrow(x) + nrow(y))
  expect_equal(adaptive$theta.denominator, nrow(x) + nrow(y))
  expect_equal(adaptive$threshold.sample.size, nrow(x) + nrow(y))
  expect_equal(adaptive$threshold.delta, delta)
  expect_match(adaptive$threshold.formula, "n1 \\+ n2")
  expect_identical(fit$diagnostics$statistic.path, "feasible")
})

test_that("adaptive eigen-floor adjustment is explicit and scale equivariant", {
  x <- rbind(
    c(1, 2, 3, 5, 8),
    c(-1, -2, -3, -5, -8)
  )
  y <- rbind(
    c(2, -1, 4, -3, 6),
    c(-2, 1, -4, 3, -6)
  )

  expect_error(
    cai_liu_xia_two_sample_test(
      x, y, threshold_delta = 1e-12, eigen_floor = 0
    ),
    "not positive definite"
  )
  fit <- cai_liu_xia_two_sample_test(
    x, y, threshold_delta = 1e-12, eigen_floor = 1e-6
  )
  diagnostics <- fit$diagnostics$adaptive.thresholding
  expect_true(diagnostics$eigen.floor.applied)
  expect_gt(diagnostics$adjusted.eigenvalues, 0)
  expect_equal(
    diagnostics$minimum.eigenvalue.after,
    diagnostics$eigen.floor.absolute,
    tolerance = 1e-12
  )
  expect_true(fit$diagnostics$precision$positive.definite)

  multiplier <- 7.5
  scaled <- cai_liu_xia_two_sample_test(
    multiplier * x, multiplier * y,
    threshold_delta = 1e-12, eigen_floor = 1e-6
  )
  expect_equal(scaled$statistic, fit$statistic, tolerance = 5e-8)
  expect_equal(scaled$p.value, fit$p.value, tolerance = 5e-8)
  expect_equal(
    scaled$components$precision * multiplier^2,
    fit$components$precision,
    tolerance = 5e-8
  )
  expect_equal(
    scaled$diagnostics$adaptive.thresholding$eigen.floor.absolute /
      multiplier^2,
    diagnostics$eigen.floor.absolute,
    tolerance = 1e-12
  )
})

test_that("CLX extreme-value CDF, quantile, and critical values invert", {
  probabilities <- c(1e-12, 0.01, 0.25, 0.5, 0.95, 1 - 1e-12)
  quantiles <- .clx_extreme_value_quantile(probabilities)
  expect_equal(
    .clx_extreme_value_cdf(quantiles),
    probabilities,
    tolerance = 2e-15
  )
  expect_equal(
    .clx_extreme_value_survival(quantiles),
    1 - probabilities,
    tolerance = 2e-15
  )

  for (alpha in c(1e-6, 0.01, 0.05, 0.25)) {
    critical <- .clx_critical_value(alpha, p = 250)
    centered <- critical - 2 * log(250) + log(log(250))
    expect_lt(
      abs(.clx_extreme_value_survival(centered) - alpha),
      1e-15
    )
  }
  expect_equal(.clx_extreme_value_cdf(c(-Inf, Inf)), c(0, 1))
  expect_equal(.clx_extreme_value_survival(c(-Inf, Inf)), c(1, 0))
  expect_equal(.clx_extreme_value_survival(c(-1e4, 1e4)), c(1, 0))
})

test_that("oracle and feasible paths respect CLX invariances", {
  set.seed(2203)
  x <- matrix(stats::rnorm(18 * 6), 18, 6)
  y <- matrix(stats::rnorm(21 * 6, 0.2), 21, 6)
  generator <- matrix(stats::rnorm(36), 6, 6)
  precision <- crossprod(generator) + diag(6)
  baseline <- cai_liu_xia_two_sample_test(
    x, y, precision = precision, precision_source = "oracle"
  )

  shift <- seq(-3, 2, length.out = 6)
  translated <- cai_liu_xia_two_sample_test(
    sweep(x, 2, shift, `+`), sweep(y, 2, shift, `+`),
    precision = precision, precision_source = "oracle"
  )
  expect_equal(translated$statistic, baseline$statistic, tolerance = 2e-12)
  expect_equal(translated$p.value, baseline$p.value, tolerance = 2e-12)

  swapped <- cai_liu_xia_two_sample_test(
    y, x, precision = precision, precision_source = "oracle"
  )
  expect_equal(swapped$statistic, baseline$statistic, tolerance = 2e-12)
  expect_equal(swapped$p.value, baseline$p.value, tolerance = 2e-12)
  expect_equal(unname(swapped$estimate), -unname(baseline$estimate),
               tolerance = 2e-12)

  scales <- c(-2, 0.5, 3, -1.5, 4, 0.75)
  inverse.scale <- diag(1 / scales)
  transformed.precision <- inverse.scale %*% precision %*% inverse.scale
  changed.units <- cai_liu_xia_two_sample_test(
    sweep(x, 2, scales, `*`), sweep(y, 2, scales, `*`),
    precision = transformed.precision, precision_source = "oracle"
  )
  expect_equal(changed.units$statistic, baseline$statistic, tolerance = 4e-12)
  expect_equal(changed.units$p.value, baseline$p.value, tolerance = 4e-12)

  permutation <- c(4, 1, 6, 2, 5, 3)
  permuted <- cai_liu_xia_two_sample_test(
    x[, permutation, drop = FALSE], y[, permutation, drop = FALSE],
    precision = precision[permutation, permutation],
    precision_source = "feasible"
  )
  feasible <- cai_liu_xia_two_sample_test(
    x, y, precision = precision, precision_source = "feasible"
  )
  expect_equal(permuted$statistic, feasible$statistic, tolerance = 4e-12)
  expect_equal(permuted$p.value, feasible$p.value, tolerance = 4e-12)

  row.permuted <- cai_liu_xia_two_sample_test(
    x[sample(nrow(x)), , drop = FALSE],
    y[sample(nrow(y)), , drop = FALSE],
    precision = precision, precision_source = "feasible"
  )
  expect_equal(row.permuted$statistic, feasible$statistic, tolerance = 4e-12)
})

test_that("adaptive CLX respects group, location, and variable symmetries", {
  set.seed(2204)
  x <- matrix(stats::rnorm(35 * 5), 35, 5)
  y <- matrix(stats::rnorm(42 * 5, 0.1), 42, 5)
  baseline <- cai_liu_xia_two_sample_test(x, y, threshold_delta = 1.25)

  shifted <- seq(-1, 1, length.out = 5)
  translated <- cai_liu_xia_two_sample_test(
    sweep(x, 2, shifted, `+`), sweep(y, 2, shifted, `+`),
    threshold_delta = 1.25
  )
  expect_equal(translated$statistic, baseline$statistic, tolerance = 1e-11)

  swapped <- cai_liu_xia_two_sample_test(y, x, threshold_delta = 1.25)
  expect_equal(swapped$statistic, baseline$statistic, tolerance = 1e-11)
  expect_equal(swapped$p.value, baseline$p.value, tolerance = 1e-11)

  permutation <- c(3, 5, 1, 4, 2)
  permuted <- cai_liu_xia_two_sample_test(
    x[, permutation, drop = FALSE], y[, permutation, drop = FALSE],
    threshold_delta = 1.25
  )
  expect_equal(permuted$statistic, baseline$statistic, tolerance = 1e-11)
  expect_equal(permuted$p.value, baseline$p.value, tolerance = 1e-11)
})

test_that("CLX uses stable centering for large common translations", {
  set.seed(2210)
  x <- matrix(stats::rnorm(24 * 5), 24, 5)
  y <- matrix(stats::rnorm(27 * 5, 0.15), 27, 5)
  generator <- matrix(stats::rnorm(25), 5, 5)
  precision <- crossprod(generator) + diag(5)
  shift <- c(1, -1, 1, -1, 1) * 1e14
  shifted.x <- sweep(x, 2, shift, `+`)
  shifted.y <- sweep(y, 2, shift, `+`)

  for (source in c("oracle", "feasible")) {
    oracle <- identical(source, "oracle")
    reference <- clx_stable_reference(
      shifted.x, shifted.y, precision, oracle
    )
    fit <- cai_liu_xia_two_sample_test(
      shifted.x, shifted.y,
      precision = precision, precision_source = source
    )
    expect_equal(unname(fit$raw.statistic), reference$M,
                 tolerance = 5e-12)
    expect_equal(unname(fit$statistic), reference$G,
                 tolerance = 5e-12)
    expect_equal(unname(fit$estimate), reference$difference,
                 tolerance = 5e-13)
    expect_equal(unname(fit$variance), reference$denominator,
                 tolerance = 5e-12)
  }

  delta <- 0.8
  eigen.floor <- sqrt(.Machine$double.eps)
  adaptive.reference <- clx_adaptive_stable_reference(
    shifted.x, shifted.y, delta, eigen.floor
  )
  statistic.reference <- clx_stable_reference(
    shifted.x, shifted.y, adaptive.reference$precision, oracle = FALSE
  )
  adaptive <- cai_liu_xia_two_sample_test(
    shifted.x, shifted.y, threshold_delta = delta,
    eigen_floor = eigen.floor
  )
  backend <- cpp_clx_adaptive_precision(
    shifted.x, shifted.y, delta, eigen.floor
  )

  expect_equal(unname(adaptive$raw.statistic), statistic.reference$M,
               tolerance = 2e-11)
  expect_equal(unname(adaptive$statistic), statistic.reference$G,
               tolerance = 2e-11)
  expect_equal(backend$pooled_covariance_scaled,
               adaptive.reference$pooled.scaled, tolerance = 3e-13)
  expect_equal(backend$theta_scaled,
               adaptive.reference$theta.scaled, tolerance = 3e-13)
  expect_equal(backend$thresholds_scaled,
               adaptive.reference$lambda.scaled, tolerance = 3e-13)
})

test_that("CLX remains stable under representable extreme global scales", {
  set.seed(2211)
  x <- matrix(stats::rnorm(26 * 4), 26, 4)
  y <- matrix(stats::rnorm(29 * 4, 0.2), 29, 4)
  generator <- matrix(stats::rnorm(16), 4, 4)
  precision <- crossprod(generator) + diag(4)
  delta <- 0.75
  eigen.floor <- sqrt(.Machine$double.eps)

  baseline.oracle <- cai_liu_xia_two_sample_test(
    x, y, precision = precision, precision_source = "oracle"
  )
  baseline.feasible <- cai_liu_xia_two_sample_test(
    x, y, precision = precision, precision_source = "feasible"
  )
  baseline.adaptive <- cai_liu_xia_two_sample_test(
    x, y, threshold_delta = delta, eigen_floor = eigen.floor
  )

  for (multiplier in c(1e-150, 1e150)) {
    changed.precision <- precision / multiplier^2
    scaled.oracle <- cai_liu_xia_two_sample_test(
      multiplier * x, multiplier * y,
      precision = changed.precision, precision_source = "oracle"
    )
    scaled.feasible <- cai_liu_xia_two_sample_test(
      multiplier * x, multiplier * y,
      precision = changed.precision, precision_source = "feasible"
    )
    scaled.adaptive <- cai_liu_xia_two_sample_test(
      multiplier * x, multiplier * y,
      threshold_delta = delta, eigen_floor = eigen.floor
    )
    backend <- cpp_clx_adaptive_precision(
      multiplier * x, multiplier * y, delta, eigen.floor
    )

    expect_equal(scaled.oracle$statistic, baseline.oracle$statistic,
                 tolerance = 2e-11)
    expect_equal(scaled.feasible$statistic, baseline.feasible$statistic,
                 tolerance = 2e-11)
    expect_equal(scaled.adaptive$statistic, baseline.adaptive$statistic,
                 tolerance = 3e-10)
    expect_true(all(is.finite(backend$pooled_covariance_scaled)))
    expect_true(all(is.finite(backend$theta_scaled)))
    expect_true(all(is.finite(backend$thresholds_scaled)))
    expect_true(all(is.finite(backend$precision_scaled)))
  }
})

test_that("feasible CLX accepts diagnosed indefinite estimates", {
  set.seed(2205)
  x <- matrix(stats::rnorm(30), 15, 2)
  y <- matrix(stats::rnorm(34), 17, 2)
  indefinite <- matrix(c(1, 2, 2, 1), 2, 2)
  fit <- cai_liu_xia_two_sample_test(
    x, y, precision = indefinite, precision_source = "feasible"
  )

  expect_true(is.finite(unname(fit$statistic)))
  expect_false(fit$diagnostics$precision$positive.definite)
  expect_lt(fit$diagnostics$precision$minimum.eigenvalue, 0)
  expect_error(
    cai_liu_xia_two_sample_test(
      x, y, precision = indefinite, precision_source = "oracle"
    ),
    "positive definite"
  )
})

test_that("CLX keeps the standard high-dimensional location-test contract", {
  set.seed(2206)
  x <- matrix(stats::rnorm(48), 12, 4)
  y <- matrix(stats::rnorm(60), 15, 4)
  fit <- cai_liu_xia_two_sample_test(
    x, y, precision = diag(4), precision_source = "oracle"
  )

  expect_named(
    fit,
    c(
      "statistic", "parameter", "p.value", "estimate", "null.value",
      "alternative", "method", "data.name", "call", "raw.statistic",
      "null.distribution", "variance", "components", "diagnostics", "n",
      "p"
    ),
    ignore.order = TRUE
  )
  expect_named(fit$statistic, "G")
  expect_named(fit$raw.statistic, "M")
  expect_identical(fit$alternative, "two.sided")
  expect_identical(fit$null.distribution$tail, "upper")
  expect_false(fit$null.distribution$exact)
  expect_true(is.finite(fit$p.value) && fit$p.value >= 0 && fit$p.value <= 1)
  expect_equal(fit$components$G, unname(fit$statistic))
  expect_equal(fit$components$M, unname(fit$raw.statistic))
})

test_that("CLX rejects ambiguous, incompatible, and degenerate inputs", {
  set.seed(2207)
  x <- matrix(stats::rnorm(30), 6, 5)
  y <- matrix(stats::rnorm(35), 7, 5)

  expect_error(
    cai_liu_xia_two_sample_test(x[, 1, drop = FALSE], y[, 1, drop = FALSE]),
    "p >= 2"
  )
  expect_error(cai_liu_xia_two_sample_test(x[1, , drop = FALSE], y),
               "at least 2 row")
  expect_error(cai_liu_xia_two_sample_test(x[, 1:4], y),
               "same number of columns")

  named.x <- x
  named.y <- y
  colnames(named.x) <- paste0("x", 1:5)
  colnames(named.y) <- paste0("y", 1:5)
  expect_error(cai_liu_xia_two_sample_test(named.x, named.y), "same names")

  expect_error(
    cai_liu_xia_two_sample_test(
      x, y, precision = diag(5), precision_source = "adaptive"
    ),
    "must be `NULL`"
  )
  expect_error(
    cai_liu_xia_two_sample_test(x, y, precision_source = "oracle"),
    "is required"
  )
  expect_error(
    cai_liu_xia_two_sample_test(x, y, precision_source = "feasible"),
    "is required"
  )
  expect_error(
    cai_liu_xia_two_sample_test(
      x, y, precision = diag(4), precision_source = "oracle"
    ),
    "5 by 5"
  )

  nonfinite.precision <- diag(5)
  nonfinite.precision[1, 1] <- Inf
  expect_error(
    cai_liu_xia_two_sample_test(
      x, y, precision = nonfinite.precision, precision_source = "oracle"
    ),
    "finite"
  )
  asymmetric <- diag(5)
  asymmetric[1, 2] <- 0.5
  expect_error(
    cai_liu_xia_two_sample_test(
      x, y, precision = asymmetric, precision_source = "feasible"
    ),
    "symmetric"
  )

  expect_error(cai_liu_xia_two_sample_test(x, y, threshold_delta = 0),
               "threshold_delta")
  expect_error(cai_liu_xia_two_sample_test(x, y, threshold_delta = Inf),
               "threshold_delta")
  expect_error(cai_liu_xia_two_sample_test(x, y, eigen_floor = -1),
               "eigen_floor")
  expect_error(cai_liu_xia_two_sample_test(x, y, eigen_floor = 1),
               "eigen_floor")

  zero.precision <- matrix(0, 5, 5)
  expect_error(
    cai_liu_xia_two_sample_test(
      x, y, precision = zero.precision, precision_source = "feasible"
    ),
    "strictly positive"
  )
  constant.x <- matrix(1, 4, 3)
  constant.y <- matrix(2, 5, 3)
  expect_error(cai_liu_xia_two_sample_test(constant.x, constant.y),
               "positive numerical scale")

  nonfinite.x <- x
  nonfinite.x[1, 1] <- NA_real_
  expect_error(cai_liu_xia_two_sample_test(nonfinite.x, y), "finite")
  expect_error(.clx_extreme_value_quantile(c(0.5, 1.1)), "\\[0, 1\\]")
  expect_error(.clx_critical_value(0, 10), "alpha")
  expect_error(.clx_critical_value(0.05, 1), "at least 2")
})

test_that("near-symmetric supplied precision is diagnosed and symmetrized", {
  set.seed(2208)
  x <- matrix(stats::rnorm(60), 15, 4)
  y <- matrix(stats::rnorm(68), 17, 4)
  precision <- diag(4)
  precision[1, 2] <- 5e-10
  precision[2, 1] <- 0
  fit <- cai_liu_xia_two_sample_test(
    x, y, precision = precision, precision_source = "feasible"
  )

  expect_true(fit$diagnostics$precision$symmetrized)
  expect_equal(
    fit$components$precision,
    t(fit$components$precision),
    tolerance = 0
  )
  expect_equal(fit$components$precision[1, 2], 2.5e-10)
})
