yzf_test_fixture <- function() {
  set.seed(1)
  matrix(
    stats::runif(80, -1, 1), 8, 10,
    dimnames = list(NULL, paste0("v", seq_len(10)))
  )
}


yzf_test_prepare <- function(x, center) {
  residual <- sweep(x, 2L, center, "-")
  scale <- apply(abs(residual), 2L, max)
  if (any(!is.finite(scale)) || any(scale <= 0)) {
    stop("invalid pure-R residual scale")
  }
  list(data = sweep(residual, 2L, scale, "/"), scale = scale)
}


yzf_test_score <- function(data, location, diagonal) {
  epsilon <- sweep(data, 2L, location, "-")
  epsilon <- sweep(epsilon, 2L, sqrt(diagonal), "/")
  radius <- sqrt(rowSums(epsilon^2))
  if (any(!is.finite(radius)) || any(radius <= 0)) {
    stop("singular pure-R training radius")
  }
  list(
    epsilon = epsilon,
    radius = radius,
    direction = epsilon / radius
  )
}


yzf_test_fit <- function(data, leave = integer(), m = 0,
                         tol = 1e-7, max_iter = 2000L) {
  retained <- data[setdiff(seq_len(nrow(data)), leave), , drop = FALSE]
  location <- colMeans(retained)
  diagonal <- apply(retained, 2L, stats::var)
  initial.diagonal <- diagonal
  stable <- FALSE
  relative.update <- Inf
  for (iteration in seq_len(max_iter)) {
    score <- yzf_test_score(retained, location, diagonal)
    numerator <- colSums(score$direction * score$radius^m)
    denominator <- sum(score$radius^(m - 1))
    next.location <- location + sqrt(diagonal) * numerator / denominator
    next.diagonal <- ncol(data) * diagonal *
      colMeans(score$direction^2)
    location.update <- max(abs(
      (next.location - location) / sqrt(initial.diagonal)
    ))
    diagonal.update <- max(abs(log(next.diagonal / diagonal)))
    relative.update <- max(location.update, diagonal.update)
    location <- next.location
    diagonal <- next.diagonal
    if (relative.update <= tol) {
      stable <- TRUE
      break
    }
  }
  score <- yzf_test_score(retained, location, diagonal)
  weighted.mean.direction <- colSums(
    score$direction * score$radius^m
  ) / sum(score$radius^m)
  list(
    location = location,
    diagonal = diagonal,
    score = score,
    iterations = iteration,
    stable = stable,
    relative.update = relative.update,
    location.score = sqrt(sum(weighted.mean.direction^2)),
    diagonal.score = max(abs(
      ncol(data) * colMeans(score$direction^2) - 1
    ))
  )
}


yzf_test_full_reference <- function(x, mu, m, tol = 1e-7,
                                    max_iter = 2000L) {
  prepared <- yzf_test_prepare(x, mu)
  fit <- yzf_test_fit(prepared$data, m = m, tol = tol,
                      max_iter = max_iter)
  score <- fit$score
  zeta.first <- mean(score$radius^(m - 1))
  zeta.second <- mean(score$radius^(2 * m))
  standardized.location <- fit$location / sqrt(fit$diagonal)
  M <- nrow(x) * max(abs(standardized.location))^2 * ncol(x) *
    (1 - nrow(x)^(-0.5)) * zeta.first^2 / zeta.second
  list(
    prepared = prepared,
    fit = fit,
    location = mu + prepared$scale * fit$location,
    zeta.first = zeta.first,
    zeta.second = zeta.second,
    moment.ratio = zeta.first^2 / zeta.second,
    M = M,
    T = M - 2 * log(ncol(x)) + log(log(ncol(x)))
  )
}


yzf_test_direction <- function(value, location, diagonal) {
  epsilon <- (value - location) / sqrt(diagonal)
  radius <- sqrt(sum(epsilon^2))
  if (radius == 0) {
    return(list(direction = numeric(length(value)), radius = 0))
  }
  list(direction = epsilon / radius, radius = radius)
}


yzf_test_sum_reference <- function(x, mu, m, tol = 1e-7,
                                   max_iter = 2000L) {
  prepared <- yzf_test_prepare(x, mu)
  data <- prepared$data
  n <- nrow(data)
  p <- ncol(data)
  pairs <- t(utils::combn(n, 2L))
  pair.count <- nrow(pairs)
  statistic.kernel <- variance.kernel <- wrong.kernel <-
    trace.inner.squared <- numeric(pair.count)
  test.inner <- numeric(pair.count)
  variance.factor <- matrix(NA_real_, pair.count, 2L)
  endpoint.radius <- endpoint.weight <- matrix(NA_real_, pair.count, 2L)
  sign.mean <- matrix(NA_real_, pair.count, p)
  location <- diagonal <- matrix(NA_real_, pair.count, p)

  for (pair in seq_len(pair.count)) {
    i <- pairs[pair, 1L]
    j <- pairs[pair, 2L]
    fit <- yzf_test_fit(
      data, c(i, j), m = 0, tol = tol, max_iter = max_iter
    )
    location[pair, ] <- fit$location
    diagonal[pair, ] <- fit$diagonal
    left <- yzf_test_direction(data[i, ], numeric(p), fit$diagonal)
    right <- yzf_test_direction(data[j, ], numeric(p), fit$diagonal)
    endpoint.radius[pair, ] <- c(left$radius, right$radius)
    endpoint.weight[pair, ] <- endpoint.radius[pair, ]^m
    inner <- sum(left$direction * right$direction)
    test.inner[pair] <- inner
    statistic.kernel[pair] <- prod(endpoint.weight[pair, ]) * inner

    keep <- setdiff(seq_len(n), c(i, j))
    retained <- lapply(
      keep,
      function(k) yzf_test_direction(
        data[k, ], numeric(p), fit$diagonal
      )
    )
    retained.direction <- do.call(
      rbind, lapply(retained, `[[`, "direction")
    )
    retained.radius <- vapply(retained, `[[`, numeric(1L), "radius")
    sign.mean[pair, ] <- colMeans(retained.direction)
    first <- sum((left$direction - sign.mean[pair, ]) *
                   right$direction)
    second <- sum((right$direction - sign.mean[pair, ]) *
                    left$direction)
    variance.factor[pair, ] <- c(first, second)
    variance.kernel[pair] <- prod(endpoint.weight[pair, ])^2 *
      first * second

    weighted.mean <- colSums(
      retained.direction * retained.radius^m
    ) / sum(retained.radius^m)
    wrong.first <- sum((left$direction - weighted.mean) *
                         right$direction)
    wrong.second <- sum((right$direction - weighted.mean) *
                          left$direction)
    wrong.kernel[pair] <- prod(endpoint.weight[pair, ])^2 *
      wrong.first * wrong.second

    trace.left <- yzf_test_direction(data[i, ], fit$location, fit$diagonal)
    trace.right <- yzf_test_direction(data[j, ], fit$location, fit$diagonal)
    trace.inner.squared[pair] <- sum(
      trace.left$direction * trace.right$direction
    )^2
  }

  T <- mean(statistic.kernel)
  ordered.kernel.sum <- 2 * sum(variance.kernel)
  variance <- 2 * ordered.kernel.sum / n^4
  wrong.variance <- 4 * sum(wrong.kernel) / n^4
  zeta.2m <- mean(endpoint.weight^2)
  trace.R2 <- p^2 * mean(trace.inner.squared)
  list(
    pairs = pairs,
    T = T,
    variance = variance,
    z = T / sqrt(variance),
    ordered.kernel.sum = ordered.kernel.sum,
    wrong.variance = wrong.variance,
    zeta.2m = zeta.2m,
    trace.R2 = trace.R2,
    statistic.kernel = statistic.kernel,
    test.inner = test.inner,
    variance.factor = variance.factor,
    variance.kernel = variance.kernel,
    endpoint.radius = endpoint.radius,
    endpoint.weight = endpoint.weight,
    sign.mean = sign.mean,
    location = location,
    diagonal = diagonal
  )
}


test_that("weighted scaled median matches all corrected estimating equations", {
  x <- yzf_test_fixture()
  center <- colMeans(x)
  prepared <- yzf_test_prepare(x, center)

  for (m in c(-1, 0, 0.5, 1)) {
    reference <- yzf_test_fit(prepared$data, m = m)
    result <- weighted_scaled_spatial_median(
      x, m = m, tol = 1e-7, max_iter = 2000L
    )
    expect_s3_class(result, "weighted_scaled_spatial_median")
    expect_equal(
      unname(result$weighted.location.standardized),
      unname(reference$location),
      tolerance = 3e-9
    )
    expect_equal(
      unname(result$scale.diagonal.standardized),
      unname(reference$diagonal),
      tolerance = 5e-9
    )
    expect_lt(abs(
      result$diagnostics$location.score.residual -
        reference$location.score
    ), 3e-12)
    expect_lt(abs(
      result$diagnostics$diagonal.score.residual -
        reference$diagonal.score
    ), 3e-12)
    expect_true(result$diagnostics$iteration.stable)
  }

  mean.fit <- weighted_scaled_spatial_median(
    x, m = 1, tol = 1e-8, max_iter = 2000L
  )
  expect_equal(
    unname(mean.fit$weighted.location), unname(colMeans(x)),
    tolerance = 3e-15
  )
  expect_match(mean.fit$diagnostics$corrected.location.update,
               "D\\^\\(1/2\\)")
  expect_match(mean.fit$diagnostics$corrected.scale.update,
               "D\\^\\(1/2\\)")
  expect_identical(mean.fit$diagnostics$regularization, "none")
})


test_that("weighted max uses literal moments, finite correction, and Gumbel law", {
  x <- yzf_test_fixture()
  mu <- rep(0, ncol(x))

  for (m in c(-1, 0, 1)) {
    reference <- yzf_test_full_reference(x, mu, m)
    result <- yan_zhao_feng_weighted_max_test(
      x, mu, m = m, tol = 1e-7, max_iter = 2000L
    )
    component <- result$components
    expected.p <- -expm1(-pi^(-0.5) * exp(-reference$T / 2))
    expect_s3_class(result, "hd_location_test")
    expect_s3_class(result, "htest")
    expect_equal(component$M.MAX, reference$M, tolerance = 8e-10)
    expect_equal(component$T.MAX, reference$T, tolerance = 8e-10)
    expect_equal(component$zeta.m.minus.1.hat,
                 reference$zeta.first, tolerance = 3e-10)
    expect_equal(component$zeta.2m.hat,
                 reference$zeta.second, tolerance = 3e-10)
    expect_equal(component$moment.ratio.hat,
                 reference$moment.ratio, tolerance = 3e-10)
    expect_equal(component$finite.sample.correction,
                 1 - nrow(x)^(-0.5), tolerance = 0)
    expect_equal(result$p.value, expected.p, tolerance = 3e-14)
  }

  inverse <- yan_zhao_feng_weighted_max_test(
    x, m = -1, tol = 1e-7, max_iter = 2000L
  )$components
  sign <- yan_zhao_feng_weighted_max_test(
    x, m = 0, tol = 1e-7, max_iter = 2000L
  )$components
  mean <- yan_zhao_feng_weighted_max_test(
    x, m = 1, tol = 1e-7, max_iter = 2000L
  )$components
  expect_equal(inverse$moment.ratio.hat,
               inverse$zeta.2m.hat, tolerance = 3e-14)
  expect_equal(sign$moment.ratio.hat,
               sign$zeta.m.minus.1.hat^2, tolerance = 3e-14)
  expect_equal(mean$moment.ratio.hat,
               1 / mean$zeta.2m.hat, tolerance = 3e-14)
})


test_that("general-m sum and direct variance match literal pure-R formulas", {
  x <- yzf_test_fixture()
  m <- 0.5
  reference <- yzf_test_sum_reference(x, numeric(ncol(x)), m)
  result <- yan_zhao_feng_weighted_maxsum_test(
    x, m = m, tol = 1e-7, max_iter = 2000L
  )
  component <- result$components

  expect_equal(component$T.SUM, reference$T, tolerance = 8e-10)
  expect_equal(component$sigma2.direct.hat,
               reference$variance, tolerance = 8e-10)
  expect_equal(component$Z.SUM, reference$z, tolerance = 8e-10)
  expect_equal(
    unname(component$pair.test.inner.product),
    reference$test.inner, tolerance = 8e-10
  )
  expect_equal(
    unname(component$pair.statistic.kernel),
    reference$statistic.kernel, tolerance = 8e-10
  )
  expect_equal(
    unname(component$pair.variance.factor),
    reference$variance.factor, tolerance = 8e-10
  )
  expect_equal(
    unname(component$pair.variance.kernel),
    reference$variance.kernel, tolerance = 8e-10
  )
  expect_equal(
    unname(component$endpoint.radius),
    reference$endpoint.radius, tolerance = 8e-10
  )
  expect_equal(
    unname(component$endpoint.radial.weight),
    reference$endpoint.weight, tolerance = 8e-10
  )
  expect_equal(component$ordered.variance.kernel.sum,
               2 * sum(component$pair.variance.kernel), tolerance = 2e-15)
  expect_equal(component$sigma2.direct.hat,
               2 * component$ordered.variance.kernel.sum / nrow(x)^4,
               tolerance = 2e-15)
  expect_equal(component$zeta.2m.crossfit.hat,
               reference$zeta.2m, tolerance = 8e-10)
  expect_equal(component$trace.R2.hat,
               reference$trace.R2, tolerance = 8e-9)
  expect_gt(abs(reference$wrong.variance - reference$variance), 1e-5)
  expect_match(result$diagnostics$variance.sign.center, "unweighted")
  expect_match(result$diagnostics$primary.sum.variance,
               "2\\*n\\^\\(-4\\)")

  ordinary.cauchy <- 0.5 * tan(pi * (0.5 - component$p.MAX)) +
    0.5 * tan(pi * (0.5 - component$p.SUM))
  expect_equal(component$Cauchy.statistic,
               ordinary.cauchy, tolerance = 2e-14)
  expect_equal(result$p.value,
               stats::pcauchy(ordinary.cauchy, lower.tail = FALSE),
               tolerance = 2e-15)
})


test_that("m=-1 sum components numerically cross-lock to INST", {
  x <- yzf_test_fixture()
  weighted <- yan_zhao_feng_weighted_maxsum_test(
    x, m = -1, tol = 1e-7, max_iter = 2000L
  )
  inst <- inst_one_sample_test(
    x, tol = 1e-7, max_iter = 2000L
  )
  component <- weighted$components
  reference <- inst$components
  comparison_tolerance <- 128 * .Machine$double.eps

  expect_equal(component$T.SUM, reference$T.INST,
               tolerance = comparison_tolerance)
  expect_equal(component$sigma2.direct.hat,
               reference$sigma2.direct.hat,
               tolerance = comparison_tolerance)
  expect_equal(component$sigma.direct.hat,
               reference$sigma.direct.hat,
               tolerance = comparison_tolerance)
  expect_equal(component$Z.SUM, unname(inst$statistic),
               tolerance = comparison_tolerance)
  expect_equal(component$p.SUM, inst$p.value,
               tolerance = comparison_tolerance)
  expect_equal(component$ordered.variance.kernel.sum,
               reference$ordered.variance.kernel.sum,
               tolerance = comparison_tolerance)
  expect_equal(component$zeta.2m.crossfit.hat,
               reference$nu2.IN.hat,
               tolerance = comparison_tolerance)
  expect_equal(component$trace.R2.hat, reference$trace.R2.hat,
               tolerance = comparison_tolerance)
  expect_equal(
    component$sigma2.oracle.pair.adjusted.diagnostic,
    reference$sigma2.oracle.factorized.diagnostic,
    tolerance = comparison_tolerance
  )
  expect_equal(
    component$trace.weighted.score.covariance.squared.hat,
    reference$trace.weighted.score.covariance.squared.hat,
    tolerance = comparison_tolerance
  )
  expect_equal(unname(component$pair.test.inner.product),
               unname(reference$pair.test.inner.product),
               tolerance = comparison_tolerance)
  expect_equal(unname(component$pair.statistic.kernel),
               unname(reference$pair.statistic.kernel),
               tolerance = comparison_tolerance)
  expect_equal(unname(component$pair.variance.factor),
               unname(reference$pair.variance.factor),
               tolerance = comparison_tolerance)
  expect_equal(unname(component$pair.variance.kernel),
               unname(reference$pair.variance.kernel),
               tolerance = comparison_tolerance)
  expect_equal(unname(component$endpoint.radius),
               unname(reference$endpoint.radius),
               tolerance = comparison_tolerance)
  expect_equal(unname(component$endpoint.radial.weight),
               unname(reference$endpoint.inverse.radius),
               tolerance = comparison_tolerance)
  expect_equal(unname(component$leaveout.location.standardized),
               unname(reference$leaveout.location.standardized),
               tolerance = comparison_tolerance)
  expect_equal(unname(component$leaveout.scale.diagonal.standardized),
               unname(reference$leaveout.scale.diagonal.standardized),
               tolerance = comparison_tolerance)
})


test_that("weighted max-sum has exactly its stated invariances", {
  x <- yzf_test_fixture()
  mu <- seq(-0.09, 0.09, length.out = ncol(x))
  baseline <- yan_zhao_feng_weighted_maxsum_test(
    x, mu, m = 0.5, tol = 1e-7, max_iter = 2000L
  )
  rows <- c(8, 2, 6, 1, 5, 3, 7, 4)
  row.result <- yan_zhao_feng_weighted_maxsum_test(
    x[rows, , drop = FALSE], mu, m = 0.5,
    tol = 1e-7, max_iter = 2000L
  )
  expect_equal(row.result$components$T.MAX,
               baseline$components$T.MAX, tolerance = 3e-8)
  expect_equal(row.result$components$T.SUM,
               baseline$components$T.SUM, tolerance = 3e-8)
  expect_equal(row.result$p.value, baseline$p.value, tolerance = 3e-8)

  multiplier <- c(-3, 0.25, 2, -5, 7, -0.4, 1.5, 9, -2, 0.1)
  scaled <- sweep(x, 2L, multiplier, "*")
  scaled.result <- yan_zhao_feng_weighted_maxsum_test(
    scaled, mu * multiplier, m = 0.5,
    tol = 1e-7, max_iter = 2000L
  )
  expect_equal(scaled.result$components$T.MAX,
               baseline$components$T.MAX, tolerance = 3e-9)
  expect_equal(scaled.result$components$Z.SUM,
               baseline$components$Z.SUM, tolerance = 3e-9)
  expect_equal(scaled.result$p.value, baseline$p.value, tolerance = 3e-9)

  shift <- seq(-5, 4, length.out = ncol(x))
  translated <- sweep(x, 2L, shift, "+")
  translated.result <- yan_zhao_feng_weighted_maxsum_test(
    translated, mu + shift, m = 0.5,
    tol = 1e-7, max_iter = 2000L
  )
  expect_equal(translated.result$components$T.MAX,
               baseline$components$T.MAX, tolerance = 3e-8)
  expect_equal(translated.result$components$Z.SUM,
               baseline$components$Z.SUM, tolerance = 3e-8)
  expect_equal(translated.result$p.value,
               baseline$p.value, tolerance = 3e-8)

  mixed <- x
  mixed[, 1L] <- x[, 1L] + 0.7 * x[, 2L]
  mixed.result <- yan_zhao_feng_weighted_maxsum_test(
    mixed, replace(mu, 1L, mu[1L] + 0.7 * mu[2L]), m = 0.5,
    tol = 1e-7, max_iter = 2000L
  )
  expect_gt(abs(mixed.result$components$T.MAX -
                  baseline$components$T.MAX), 1e-4)
})


test_that("Gumbel and signed-log Cauchy tails remain stable at extremes", {
  ordinary <- .yzf_gumbel_tail(0)
  expect_equal(ordinary$p.value,
               -expm1(-pi^(-0.5)), tolerance = 2e-15)
  expect_equal(ordinary$log.cdf, -pi^(-0.5), tolerance = 2e-15)
  upper <- .yzf_gumbel_tail(2000)
  lower <- .yzf_gumbel_tail(-2000)
  expect_identical(upper$p.value, 0)
  expect_true(is.finite(upper$log.p.value))
  expect_identical(lower$p.value, 1)
  expect_identical(lower$log.cdf, -Inf)

  p1 <- 0.031
  p2 <- 0.42
  stable <- .yzf_cauchy_combine(
    p1, log(p1), log1p(-p1), p2, log(p2), log1p(-p2)
  )
  raw <- 0.5 * tan(pi * (0.5 - p1)) +
    0.5 * tan(pi * (0.5 - p2))
  expect_equal(stable$statistic, raw, tolerance = 2e-14)
  expect_equal(stable$p.value,
               stats::pcauchy(raw, lower.tail = FALSE), tolerance = 2e-15)

  tiny <- .yzf_cauchy_combine(
    0, -1000, 0, 0, -1000, 0
  )
  expect_identical(tiny$p.value, 0)
  expect_equal(tiny$log.absolute.statistic,
               1000 - log(pi), tolerance = 3e-13)
  expect_error(
    .yzf_cauchy_combine(0, -Inf, 0, 1, 0, -Inf),
    "indeterminate"
  )
})


test_that("strict controls and all no-repair failure contracts are explicit", {
  x <- yzf_test_fixture()
  expect_error(
    weighted_scaled_spatial_median(
      x, m = 0.5, max_iter = 1L, strict = TRUE
    ),
    "did not stabilize"
  )
  expect_warning(
    loose.fit <- weighted_scaled_spatial_median(
      x, m = 0.5, max_iter = 1L, strict = FALSE
    ),
    "did not stabilize"
  )
  expect_false(loose.fit$diagnostics$iteration.stable)

  expect_error(
    yan_zhao_feng_weighted_maxsum_test(
      x, m = 0.5, max_iter = 1L, strict = TRUE
    ),
    "did not stabilize"
  )
  expect_warning(
    loose.test <- yan_zhao_feng_weighted_maxsum_test(
      x, m = 0.5, max_iter = 1L, strict = FALSE
    ),
    "did not stabilize"
  )
  expect_false(loose.test$diagnostics$iteration.stable)
  expect_false(loose.test$diagnostics$pair.iteration.stable)

  zero.at.mean <- rbind(c(-1, -2), c(0, 0), c(1, 2), c(2, -1))
  expect_error(
    yan_zhao_feng_weighted_max_test(zero.at.mean, m = -1),
    "r\\^\\(m-1\\).*singular.*no perturbation"
  )
  expect_error(
    yan_zhao_feng_weighted_max_test(
      zero.at.mean, m = 0, zero_tol = 100
    ),
    "at or below `zero_tol`"
  )
  expect_error(
    yan_zhao_feng_weighted_max_test(cbind(1:5, rep(2, 5)), m = 1),
    "constant.*No ridge"
  )

  negative.variance <- matrix(
    c(0.906749379821122, 0.50453615700826, 0.810187134891748,
      0.946549799758941, -0.183502823580056,
      -0.240321868099272, 0.92634824058041, -0.863568638917059,
      -0.735817824490368, 0.98909410322085),
    nrow = 5L, ncol = 2L
  )
  expect_error(
    yan_zhao_feng_weighted_maxsum_test(
      negative.variance, m = 1, tol = 1e-6, max_iter = 2000L
    ),
    "direct feasible variance.*strictly positive.*no absolute value"
  )
  expect_error(weighted_scaled_spatial_median(x, m = 1.01),
               "no greater than one")
  expect_error(weighted_scaled_spatial_median(x, m = Inf), "finite")
  expect_error(yan_zhao_feng_weighted_max_test(matrix(1:6, ncol = 1)),
               "at least two variables")
  expect_error(yan_zhao_feng_weighted_maxsum_test(x[1:3, ]),
               "at least 4 row")
  expect_identical(loose.test$diagnostics$regularization, "none")
  expect_identical(loose.test$diagnostics$weight.cap, "none")
  expect_identical(loose.test$diagnostics$variance.repair, "none")
})


test_that("extreme coordinate units are neutral and metadata are complete", {
  x <- yzf_test_fixture()
  baseline <- yan_zhao_feng_weighted_maxsum_test(
    x, m = 1, tol = 1e-7, max_iter = 2000L
  )
  multiplier <- c(
    1e150, -1e-150, 1e100, -1e-100, 3, -0.2, 7, -9, 1e50, -1e50
  )
  extreme <- sweep(x, 2L, multiplier, "*")
  result <- yan_zhao_feng_weighted_maxsum_test(
    extreme, m = 1, tol = 1e-7, max_iter = 2000L
  )
  expect_true(all(is.finite(result$statistic)))
  expect_true(is.finite(result$p.value))
  expect_equal(result$components$T.MAX,
               baseline$components$T.MAX, tolerance = 5e-9)
  expect_equal(result$components$Z.SUM,
               baseline$components$Z.SUM, tolerance = 5e-9)
  expect_equal(result$p.value, baseline$p.value, tolerance = 5e-9)
  expect_identical(names(result$estimate), colnames(x))
  expect_identical(names(result$null.value), colnames(x))
  expect_identical(colnames(result$components$leaveout.location),
                   colnames(x))
  expect_identical(result$diagnostics$primary.sum.calibration,
                   "direct feasible variance only")
  expect_identical(result$diagnostics$oracle.factorized.variance.role,
                   "diagnostic only")
  expect_match(result$diagnostics$applicability, "asymptotic")
  expect_identical(result$alternative, "two.sided")
})
