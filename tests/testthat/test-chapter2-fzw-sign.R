fzw_direction_reference <- function(value) {
  radius <- sqrt(sum(value^2))
  if (radius == 0) {
    return(rep(0, length(value)))
  }
  value / radius
}


fzw_prepare_reference <- function(x, y) {
  anchor <- x[1L, ]
  x.centered <- sweep(x, 2L, anchor, "-")
  y.centered <- sweep(y, 2L, anchor, "-")
  scale <- apply(abs(rbind(x.centered, y.centered)), 2L, max)
  list(
    x = sweep(x.centered, 2L, scale, "/"),
    y = sweep(y.centered, 2L, scale, "/"),
    anchor = anchor,
    scale = scale
  )
}


fzw_score_reference <- function(data, keep, theta, diagonal) {
  epsilon <- sweep(data[keep, , drop = FALSE], 2L, theta, "-")
  epsilon <- sweep(epsilon, 2L, sqrt(diagonal), "/")
  radii <- sqrt(rowSums(epsilon^2))
  if (any(radii == 0)) {
    stop("coincident training observation")
  }
  signs <- epsilon / radii
  p <- ncol(data)
  list(
    sign.sum = colSums(signs),
    sign.square.sum = colSums(signs^2),
    inverse.radius.sum = sum(1 / radii),
    minimum.radius = min(radii),
    location.residual = max(abs(colMeans(signs))),
    scale.residual = max(abs(p * colMeans(signs^2) - 1))
  )
}


fzw_fit_reference <- function(data, leave = integer(), tol = 1e-7,
                              max_iter = 500L) {
  keep <- setdiff(seq_len(nrow(data)), leave)
  theta <- colMeans(data[keep, , drop = FALSE])
  diagonal <- apply(data[keep, , drop = FALSE], 2L, stats::var)
  diagonal <- diagonal / exp(mean(log(diagonal)))
  p <- ncol(data)
  m <- length(keep)

  for (iteration in 0:max_iter) {
    score <- fzw_score_reference(data, keep, theta, diagonal)
    if (score$location.residual <= tol && score$scale.residual <= tol) {
      return(list(
        theta = theta,
        diagonal = diagonal,
        iterations = iteration,
        location.residual = score$location.residual,
        scale.residual = score$scale.residual,
        minimum.training.radius = score$minimum.radius
      ))
    }
    if (iteration == max_iter) {
      stop("reference did not converge")
    }
    theta <- theta + sqrt(diagonal) * score$sign.sum /
      score$inverse.radius.sum
    diagonal <- diagonal * p * score$sign.square.sum / m
    diagonal <- diagonal / exp(mean(log(diagonal)))
  }
  stop("unreachable")
}


fzw_group_reference <- function(data, tol = 1e-7, max_iter = 500L) {
  n <- nrow(data)
  p <- ncol(data)
  full <- fzw_fit_reference(data, tol = tol, max_iter = max_iter)
  theta <- diagonal <- matrix(NA_real_, n, p)
  iterations <- location.residual <- scale.residual <-
    minimum.training.radius <- numeric(n)
  own.direction <- matrix(NA_real_, n, p)
  own.radius <- numeric(n)
  for (i in seq_len(n)) {
    fit <- fzw_fit_reference(data, i, tol, max_iter)
    theta[i, ] <- fit$theta
    diagonal[i, ] <- fit$diagonal
    iterations[i] <- fit$iterations
    location.residual[i] <- fit$location.residual
    scale.residual[i] <- fit$scale.residual
    minimum.training.radius[i] <- fit$minimum.training.radius
    epsilon <- (data[i, ] - fit$theta) / sqrt(fit$diagonal)
    own.radius[i] <- sqrt(sum(epsilon^2))
    own.direction[i, ] <- fzw_direction_reference(epsilon)
  }
  list(
    full = full,
    theta = theta,
    diagonal = diagonal,
    iterations = iterations,
    location.residual = location.residual,
    scale.residual = scale.residual,
    minimum.training.radius = minimum.training.radius,
    own.direction = own.direction,
    own.radius = own.radius,
    own.inverse.radius = 1 / own.radius,
    c.hat = mean(1 / own.radius)
  )
}


fzw_reference <- function(x, y, tol = 1e-7, max_iter = 500L) {
  prepared <- fzw_prepare_reference(x, y)
  x.data <- prepared$x
  y.data <- prepared$y
  n1 <- nrow(x)
  n2 <- nrow(y)
  p <- ncol(x)
  group1 <- fzw_group_reference(x.data, tol, max_iter)
  group2 <- fzw_group_reference(y.data, tol, max_iter)

  numerator.inner <- numerator.contribution <- matrix(NA_real_, n1, n2)
  for (i in seq_len(n1)) {
    for (j in seq_len(n2)) {
      sign1 <- fzw_direction_reference(
        (x.data[i, ] - group2$theta[j, ]) /
          sqrt(group1$diagonal[i, ])
      )
      sign2 <- fzw_direction_reference(
        (y.data[j, ] - group1$theta[i, ]) /
          sqrt(group2$diagonal[j, ])
      )
      numerator.inner[i, j] <- sum(sign1 * sign2)
      numerator.contribution[i, j] <- -numerator.inner[i, j]
    }
  }
  R.n <- mean(numerator.contribution)

  bridge1 <- sqrt(group1$full$diagonal / group2$full$diagonal)
  bridge2 <- 1 / bridge1
  trace1.pairs <- matrix(0, n1, n1)
  trace2.pairs <- matrix(0, n2, n2)
  trace3.pairs <- matrix(NA_real_, n1, n2)
  for (k in seq_len(n1)) {
    for (l in setdiff(seq_len(n1), k)) {
      trace1.pairs[l, k] <- sum(
        group1$own.direction[l, ] * bridge1 *
          group1$own.direction[k, ]
      )^2
    }
  }
  for (k in seq_len(n2)) {
    for (l in setdiff(seq_len(n2), k)) {
      trace2.pairs[l, k] <- sum(
        group2$own.direction[l, ] * bridge2 *
          group2$own.direction[k, ]
      )^2
    }
  }
  for (l in seq_len(n1)) {
    for (k in seq_len(n2)) {
      trace3.pairs[l, k] <- sum(
        group1$own.direction[l, ] * group2$own.direction[k, ]
      )^2
    }
  }
  trace.A1 <- p^2 * group2$c.hat^2 / group1$c.hat^2 *
    sum(trace1.pairs) / (n1 * (n1 - 1))
  trace.A2 <- p^2 * group1$c.hat^2 / group2$c.hat^2 *
    sum(trace2.pairs) / (n2 * (n2 - 1))
  trace.A3 <- p^2 * sum(trace3.pairs) / (n1 * n2)
  variance.term1 <- 2 * trace.A1 / (n1 * (n1 - 1) * p^2)
  variance.term2 <- 2 * trace.A2 / (n2 * (n2 - 1) * p^2)
  variance.term3 <- 4 * trace.A3 / (n1 * n2 * p^2)
  sigma2 <- variance.term1 + variance.term2 + variance.term3

  list(
    prepared = prepared,
    group1 = group1,
    group2 = group2,
    numerator.inner = numerator.inner,
    numerator.contribution = numerator.contribution,
    R.n = R.n,
    bridge1 = bridge1,
    bridge2 = bridge2,
    trace1.pairs = trace1.pairs,
    trace2.pairs = trace2.pairs,
    trace3.pairs = trace3.pairs,
    trace.A1 = trace.A1,
    trace.A2 = trace.A2,
    trace.A3 = trace.A3,
    variance.term1 = variance.term1,
    variance.term2 = variance.term2,
    variance.term3 = variance.term3,
    sigma2 = sigma2,
    statistic = R.n / sqrt(sigma2)
  )
}


test_that("Feng-Zou-Wang matches the literal feasible pure-R formulas", {
  x <- rbind(
    c(0.7, -1.1, 0.2, 1.3),
    c(-0.4, 0.8, 1.2, -0.6),
    c(1.1, 0.3, -0.9, 0.5),
    c(-1.2, -0.5, 0.7, 1.0),
    c(0.2, 1.4, -0.4, -1.1),
    c(1.5, -0.2, 0.5, 0.1),
    c(-0.8, 0.6, -1.3, 0.9)
  )
  y <- rbind(
    c(-0.3, 0.7, -0.1, 1.1),
    c(1.2, -0.6, 0.8, -0.4),
    c(-1.1, 0.2, 1.3, 0.6),
    c(0.5, 1.1, -0.7, -0.9),
    c(0.9, -1.3, 0.4, 0.2),
    c(-0.6, -0.4, -1.2, 1.4),
    c(1.4, 0.5, -0.3, -0.8),
    c(-0.9, 1.3, 0.9, -0.2)
  )
  reference <- fzw_reference(x, y, tol = 1e-7)
  result <- feng_zou_wang_two_sample_sign_test(x, y, tol = 1e-7)
  components <- result$components

  expect_s3_class(result, "hd_location_test")
  expect_s3_class(result, "htest")
  expect_equal(unname(result$statistic), reference$statistic,
               tolerance = 2e-10)
  expect_equal(unname(result$raw.statistic), reference$R.n,
               tolerance = 2e-12)
  expect_equal(components$R.n, reference$R.n, tolerance = 2e-12)
  expect_equal(components$sigma2.hat, reference$sigma2,
               tolerance = 3e-11)
  expect_equal(unname(result$variance), reference$sigma2,
               tolerance = 3e-11)
  expect_equal(components$c1.hat, reference$group1$c.hat,
               tolerance = 2e-11)
  expect_equal(components$c2.hat, reference$group2$c.hat,
               tolerance = 2e-11)
  expect_equal(components$trace.A1.squared.hat, reference$trace.A1,
               tolerance = 2e-9)
  expect_equal(components$trace.A2.squared.hat, reference$trace.A2,
               tolerance = 2e-9)
  expect_equal(components$trace.A3tA3.hat, reference$trace.A3,
               tolerance = 2e-9)
  expect_equal(components$variance.term1, reference$variance.term1,
               tolerance = 2e-11)
  expect_equal(components$variance.term2, reference$variance.term2,
               tolerance = 2e-11)
  expect_equal(components$variance.term3, reference$variance.term3,
               tolerance = 2e-11)
  expect_equal(unname(components$numerator.inner.product),
               reference$numerator.inner, tolerance = 2e-11)
  expect_equal(unname(components$numerator.contribution),
               reference$numerator.contribution, tolerance = 2e-11)
  expect_equal(unname(components$trace.A1.ordered.pair.squared),
               reference$trace1.pairs, tolerance = 2e-11)
  expect_equal(unname(components$trace.A2.ordered.pair.squared),
               reference$trace2.pairs, tolerance = 2e-11)
  expect_equal(unname(components$trace.A3.cross.pair.squared),
               reference$trace3.pairs, tolerance = 2e-11)
  expect_equal(unname(components$bridge.A1.diagonal.standardized),
               reference$bridge1, tolerance = 2e-11)
  expect_equal(unname(components$bridge.A2.diagonal.standardized),
               reference$bridge2, tolerance = 2e-11)
  expect_equal(
    unname(components$leaveout.fit$group1$location.standardized),
    reference$group1$theta, tolerance = 3e-10
  )
  expect_equal(
    unname(components$leaveout.fit$group2$location.standardized),
    reference$group2$theta, tolerance = 3e-10
  )
  expect_equal(
    unname(components$leaveout.fit$group1$scale.diagonal.standardized),
    reference$group1$diagonal, tolerance = 3e-10
  )
  expect_equal(
    unname(components$leaveout.fit$group2$scale.diagonal.standardized),
    reference$group2$diagonal, tolerance = 3e-10
  )
  expect_equal(
    unname(components$leaveout.fit$group1$own.direction),
    reference$group1$own.direction, tolerance = 3e-10
  )
  expect_equal(
    unname(components$leaveout.fit$group2$own.direction),
    reference$group2$own.direction, tolerance = 3e-10
  )
  expect_equal(result$p.value,
               stats::pnorm(reference$statistic, lower.tail = FALSE),
               tolerance = 1e-15)
  expect_identical(result$null.distribution$tail, "upper")
  expect_false(result$null.distribution$exact)
})


test_that("cross-centering, inverse radii, and ordered denominators are locked", {
  set.seed(2811)
  x <- matrix(stats::runif(42, -1, 1), 7, 6)
  y <- matrix(stats::runif(48, -1, 1), 8, 6)
  result <- feng_zou_wang_two_sample_sign_test(x, y, tol = 1e-6)
  cmp <- result$components
  g1 <- cmp$leaveout.fit$group1
  g2 <- cmp$leaveout.fit$group2
  prepared <- fzw_prepare_reference(x, y)

  rebuilt <- wrong.own.center <- matrix(NA_real_, nrow(x), nrow(y))
  for (i in seq_len(nrow(x))) {
    for (j in seq_len(nrow(y))) {
      u1 <- fzw_direction_reference(
        (prepared$x[i, ] - unname(g2$location.standardized[j, ])) /
          sqrt(unname(g1$scale.diagonal.standardized[i, ]))
      )
      u2 <- fzw_direction_reference(
        (prepared$y[j, ] - unname(g1$location.standardized[i, ])) /
          sqrt(unname(g2$scale.diagonal.standardized[j, ]))
      )
      rebuilt[i, j] <- -sum(u1 * u2)
      wrong.u1 <- fzw_direction_reference(
        (prepared$x[i, ] - unname(g1$location.standardized[i, ])) /
          sqrt(unname(g1$scale.diagonal.standardized[i, ]))
      )
      wrong.u2 <- fzw_direction_reference(
        (prepared$y[j, ] - unname(g2$location.standardized[j, ])) /
          sqrt(unname(g2$scale.diagonal.standardized[j, ]))
      )
      wrong.own.center[i, j] <- -sum(wrong.u1 * wrong.u2)
    }
  }
  expect_equal(unname(cmp$numerator.contribution), rebuilt,
               tolerance = 3e-14)
  expect_equal(cmp$R.n, mean(rebuilt), tolerance = 3e-15)
  expect_gt(max(abs(wrong.own.center - rebuilt)), 1e-3)
  expect_match(result$diagnostics$numerator.structure,
               "opposite-group leave-one-out")

  expect_equal(cmp$c1.hat, mean(1 / unname(g1$own.radius)),
               tolerance = 2e-15)
  expect_equal(cmp$c2.hat, mean(1 / unname(g2$own.radius)),
               tolerance = 2e-15)
  expect_gt(abs(cmp$c1.hat - mean(unname(g1$own.radius))), 1e-3)
  expect_equal(cmp$ordered.pair.count1, nrow(x) * (nrow(x) - 1))
  expect_equal(cmp$ordered.pair.count2, nrow(y) * (nrow(y) - 1))
  expect_equal(cmp$cross.pair.count, nrow(x) * nrow(y))
  expect_equal(unname(diag(cmp$trace.A1.ordered.pair.squared)),
               rep(0, nrow(x)))
  expect_equal(unname(diag(cmp$trace.A2.ordered.pair.squared)),
               rep(0, nrow(y)))

  p <- ncol(x)
  expect_equal(
    cmp$trace.A1.squared.hat,
    p^2 * cmp$c2.hat^2 / cmp$c1.hat^2 *
      sum(cmp$trace.A1.ordered.pair.squared) /
      cmp$ordered.pair.count1,
    tolerance = 3e-14
  )
  expect_equal(
    cmp$trace.A2.squared.hat,
    p^2 * cmp$c1.hat^2 / cmp$c2.hat^2 *
      sum(cmp$trace.A2.ordered.pair.squared) /
      cmp$ordered.pair.count2,
    tolerance = 3e-14
  )
  expect_equal(
    cmp$trace.A3tA3.hat,
    p^2 * sum(cmp$trace.A3.cross.pair.squared) / cmp$cross.pair.count,
    tolerance = 3e-14
  )
  expect_equal(cmp$sigma2.hat,
               cmp$variance.term1 + cmp$variance.term2 +
                 cmp$variance.term3,
               tolerance = 3e-15)
  expect_match(result$diagnostics$calibration, "local-alternative oracle")
  expect_identical(result$diagnostics$regularization, "none")
  expect_identical(result$diagnostics$variance.repair, "none")
})


test_that("Feng-Zou-Wang is permutation, group-swap, and scalar invariant", {
  set.seed(2812)
  x <- matrix(stats::runif(48, -0.9, 0.9), 8, 6)
  y <- matrix(stats::runif(54, -0.9, 0.9), 9, 6)
  baseline <- feng_zou_wang_two_sample_sign_test(x, y, tol = 1e-6)
  permuted <- feng_zou_wang_two_sample_sign_test(
    x[c(5, 1, 8, 2, 7, 3, 6, 4), ],
    y[c(3, 8, 1, 9, 4, 7, 2, 6, 5), ],
    tol = 1e-6
  )
  swapped <- feng_zou_wang_two_sample_sign_test(y, x, tol = 1e-6)

  expect_equal(permuted$statistic, baseline$statistic, tolerance = 3e-8)
  expect_equal(permuted$raw.statistic, baseline$raw.statistic,
               tolerance = 3e-9)
  expect_equal(permuted$variance, baseline$variance, tolerance = 3e-9)
  expect_equal(swapped$statistic, baseline$statistic, tolerance = 3e-8)
  expect_equal(swapped$raw.statistic, baseline$raw.statistic,
               tolerance = 3e-9)
  expect_equal(swapped$variance, baseline$variance, tolerance = 3e-9)
  expect_equal(swapped$components$trace.A1.squared.hat,
               baseline$components$trace.A2.squared.hat,
               tolerance = 3e-8)
  expect_equal(swapped$components$trace.A2.squared.hat,
               baseline$components$trace.A1.squared.hat,
               tolerance = 3e-8)

  coordinate.scale <- c(-0.08, 0.3, -1.7, 4, 11, -0.6)
  shift <- seq(-0.8, 0.9, length.out = 6)
  transformed <- feng_zou_wang_two_sample_sign_test(
    sweep(sweep(x, 2L, coordinate.scale, "*"), 2L, shift, "+"),
    sweep(sweep(y, 2L, coordinate.scale, "*"), 2L, shift, "+"),
    tol = 1e-6
  )
  expect_equal(transformed$statistic, baseline$statistic,
               tolerance = 4e-8)
  expect_equal(transformed$raw.statistic, baseline$raw.statistic,
               tolerance = 4e-9)
  expect_equal(transformed$variance, baseline$variance, tolerance = 4e-9)
  expect_equal(transformed$components$trace.A1.squared.hat,
               baseline$components$trace.A1.squared.hat,
               tolerance = 5e-8)
  expect_equal(transformed$components$trace.A2.squared.hat,
               baseline$components$trace.A2.squared.hat,
               tolerance = 5e-8)
  expect_equal(transformed$components$trace.A3tA3.hat,
               baseline$components$trace.A3tA3.hat,
               tolerance = 5e-8)
})


test_that("the scalar-invariant test is generally not rotation invariant", {
  set.seed(2813)
  x <- matrix(stats::runif(48, -1, 1), 8, 6)
  y <- matrix(stats::runif(54, -1, 1), 9, 6)
  baseline <- feng_zou_wang_two_sample_sign_test(x, y, tol = 1e-6)
  q <- qr.Q(qr(matrix(stats::rnorm(36), 6, 6)))
  rotated <- feng_zou_wang_two_sample_sign_test(
    x %*% q, y %*% q, tol = 1e-6
  )

  expect_gt(abs(unname(rotated$statistic - baseline$statistic)), 1e-3)
  expect_gt(abs(unname(rotated$raw.statistic - baseline$raw.statistic)),
            1e-4)
})


test_that("Feng-Zou-Wang handles extreme finite coordinate units", {
  set.seed(2814)
  x <- matrix(stats::runif(48, -0.7, 0.7), 8, 6)
  y <- matrix(stats::runif(54, -0.7, 0.7), 9, 6)
  baseline <- feng_zou_wang_two_sample_sign_test(x, y, tol = 1e-6)
  coordinate.scale <- c(1e-300, -1e300, 1e-200, -1e200, 1e-100, 1e100)
  extreme <- feng_zou_wang_two_sample_sign_test(
    sweep(x, 2L, coordinate.scale, "*"),
    sweep(y, 2L, coordinate.scale, "*"),
    tol = 1e-6
  )
  expect_true(is.finite(extreme$statistic))
  expect_true(is.finite(extreme$p.value))
  expect_equal(extreme$statistic, baseline$statistic, tolerance = 5e-8)
  expect_equal(extreme$raw.statistic, baseline$raw.statistic,
               tolerance = 5e-9)
  expect_equal(extreme$variance, baseline$variance, tolerance = 5e-9)
  expect_true(all(is.finite(
    extreme$components$leaveout.fit$group1$log.scale.diagonal
  )))
  expect_true(all(is.finite(
    extreme$components$leaveout.fit$group2$log.scale.diagonal
  )))

  set.seed(2815)
  x.safe <- matrix(stats::runif(48, -1.25, 1.25), 8, 6)
  y.safe <- matrix(stats::runif(54, -1.25, 1.25), 9, 6)
  x.safe[1, 1] <- -1.15
  y.safe[1, 1] <- 1.2
  safe <- feng_zou_wang_two_sample_sign_test(x.safe, y.safe, tol = 1e-6)
  huge <- feng_zou_wang_two_sample_sign_test(
    x.safe * 1e308, y.safe * 1e308, tol = 1e-6
  )
  expect_gte(huge$diagnostics$subtraction.overflow.fallback.columns, 1)
  expect_true(all(is.finite(
    huge$components$leaveout.fit$group1$location
  )))
  expect_true(all(is.finite(
    huge$components$leaveout.fit$group2$location
  )))
  expect_equal(huge$statistic, safe$statistic, tolerance = 8e-8)
  expect_equal(huge$raw.statistic, safe$raw.statistic, tolerance = 8e-9)
  expect_equal(huge$variance, safe$variance, tolerance = 8e-9)
})


test_that("U(0)=0 is used only where no inverse-radius divisor is needed", {
  a <- c(1.0, 0.2, -0.4, 0.7)
  b <- c(-0.3, 1.1, 0.5, -0.2)
  c.vector <- c(0.6, -0.7, 1.2, 0.4)
  d <- c(0.4, 0.8, -1.0, 0.3)
  y <- rbind(
    a, -a, b, -b, c.vector, -c.vector, d, -d,
    c(0.5, -0.2, 0.9, 1.1)
  )
  set.seed(1)
  x <- matrix(stats::runif(40, -1, 1), 10, 4)
  x[1, ] <- 0
  result <- feng_zou_wang_two_sample_sign_test(x, y, tol = 1e-7)

  expect_gt(result$diagnostics$numerator.zero.sign.uses, 0)
  expect_equal(result$components$numerator.contribution[1, 9], 0,
               tolerance = 1e-13)
  expect_true(all(
    result$components$leaveout.fit$group1$own.radius > 0
  ))
  expect_true(all(
    result$components$leaveout.fit$group2$own.radius > 0
  ))
  expect_match(result$diagnostics$sign.at.zero, "U\\(0\\) = 0")
  expect_match(result$diagnostics$sign.at.zero, "inverse radii are required")
})


test_that("fit diagnostics and asymptotic-calibration flags are explicit", {
  set.seed(2811)
  x <- matrix(stats::runif(42, -1, 1), 7, 6)
  y <- matrix(stats::runif(48, -1, 1), 8, 6)
  result <- feng_zou_wang_two_sample_sign_test(x, y, tol = 1e-6)

  expect_true(result$diagnostics$all.fits.converged)
  expect_lte(result$diagnostics$maximum.location.equation.residual, 1e-6)
  expect_lte(result$diagnostics$maximum.scale.equation.residual, 1e-6)
  expect_gt(result$diagnostics$minimum.training.radius, 0)
  expect_length(result$diagnostics$leaveout.fit$group1$iterations, nrow(x))
  expect_length(result$diagnostics$leaveout.fit$group2$iterations, nrow(y))
  expect_true(result$diagnostics$applicability$p.at.most.50)
  expect_false(
    result$diagnostics$applicability$p.at.least.total.n.squared
  )
  expect_match(result$diagnostics$applicability$note, "bootstrap")
  expect_match(result$diagnostics$bias.handling, "no additive bias")
  expect_identical(result$diagnostics$bias.repair, "none")
})


test_that("Feng-Zou-Wang rejects invalid, degenerate, and unconverged fits", {
  expect_error(
    feng_zou_wang_two_sample_sign_test(
      matrix(1:8, 2, 4), matrix(1:12, 3, 4)
    ),
    "at least 3 row"
  )
  expect_error(
    feng_zou_wang_two_sample_sign_test(
      matrix(stats::rnorm(24), 6, 4),
      matrix(stats::rnorm(30), 6, 5)
    ),
    "same number of columns"
  )
  expect_error(
    feng_zou_wang_two_sample_sign_test(
      matrix(stats::rnorm(24), 6, 4),
      matrix(stats::rnorm(24), 6, 4), tol = 0
    ),
    "positive"
  )
  expect_error(
    feng_zou_wang_two_sample_sign_test(
      matrix(stats::rnorm(24), 6, 4),
      matrix(stats::rnorm(24), 6, 4), max_iter = 2.5
    ),
    "positive integer"
  )

  set.seed(2817)
  x <- matrix(stats::runif(30, -1, 1), 6, 5)
  y <- matrix(stats::runif(35, -1, 1), 7, 5)
  x[, 3] <- 1
  expect_error(
    feng_zou_wang_two_sample_sign_test(x, y),
    "degenerate.*No ridge"
  )

  coincident <- rbind(
    c(-1, -2, -3),
    c(0, 0, 0),
    c(1, 2, 3),
    c(2, -1, 1),
    c(-2, 1, -1)
  )
  other <- rbind(
    c(0.2, -0.4, 0.7),
    c(-0.8, 1.1, -0.2),
    c(1.3, 0.5, -0.9),
    c(-0.3, -1.2, 1.0),
    c(0.9, -0.1, 0.4)
  )
  expect_error(
    feng_zou_wang_two_sample_sign_test(coincident[1:3, ], other),
    "training observation exactly"
  )
  expect_error(
    feng_zou_wang_two_sample_sign_test(
      matrix(stats::runif(30), 6, 5),
      matrix(stats::runif(35), 7, 5),
      tol = 1e-14, max_iter = 1L
    ),
    "did not converge"
  )
})
