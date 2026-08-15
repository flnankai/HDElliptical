lwz_direction_reference <- function(value, location, diagonal,
                                    zero_tol = 0, allow_zero = FALSE) {
  epsilon <- (value - location) / sqrt(diagonal)
  radius <- sqrt(sum(epsilon^2))
  if (!is.finite(radius) || radius <= zero_tol) {
    if (allow_zero) {
      return(list(direction = numeric(length(value)), radius = 0))
    }
    stop("singular reference radius")
  }
  list(direction = epsilon / radius, radius = radius)
}


lwz_prepare_reference <- function(x, y) {
  anchor <- x[1L, ]
  x.centered <- sweep(x, 2L, anchor, "-")
  y.centered <- sweep(y, 2L, anchor, "-")
  scale <- apply(abs(rbind(x.centered, y.centered)), 2L, max)
  if (any(!is.finite(scale)) || any(scale <= 0)) {
    stop("invalid pooled reference scale")
  }
  list(
    x = sweep(x.centered, 2L, scale, "/"),
    y = sweep(y.centered, 2L, scale, "/"),
    anchor = anchor,
    scale = scale
  )
}


lwz_score_reference <- function(data, location, diagonal, zero_tol = 0) {
  evaluated <- lapply(
    seq_len(nrow(data)),
    function(i) lwz_direction_reference(
      data[i, ], location, diagonal, zero_tol
    )
  )
  directions <- t(vapply(
    evaluated, function(value) value$direction, numeric(ncol(data))
  ))
  radii <- vapply(evaluated, function(value) value$radius, numeric(1L))
  list(
    directions = directions,
    radii = radii,
    sign.sum = colSums(directions),
    sign.square.sum = colSums(directions^2),
    inverse.radius.sum = sum(1 / radii),
    location.residual = max(abs(colMeans(directions))),
    diagonal.residual = max(abs(
      ncol(data) * colMeans(directions^2) - 1
    ))
  )
}


lwz_fit_reference <- function(data, tol = 1e-7, max_iter = 2000L,
                              zero_tol = 0) {
  n <- nrow(data)
  p <- ncol(data)
  location <- colMeans(data)
  diagonal <- apply(data, 2L, stats::var)
  if (any(!is.finite(diagonal)) || any(diagonal <= 0)) {
    stop("non-positive reference marginal variance")
  }
  diagonal <- diagonal / exp(mean(log(diagonal)))
  initial.diagonal <- diagonal
  stable <- FALSE
  relative.update <- Inf
  location.update <- Inf
  diagonal.update <- Inf
  minimum.radius <- Inf

  for (iteration in 0:max_iter) {
    score <- lwz_score_reference(data, location, diagonal, zero_tol)
    minimum.radius <- min(minimum.radius, score$radii)
    score.residual <- max(
      score$location.residual, score$diagonal.residual
    )
    if (score.residual <= tol) {
      stable <- TRUE
      break
    }
    if (iteration == max_iter) break
    next.location <- location + sqrt(diagonal) *
      score$sign.sum / score$inverse.radius.sum
    next.diagonal <- diagonal * p * score$sign.square.sum / n
    next.diagonal <- next.diagonal / exp(mean(log(next.diagonal)))
    location.update <- max(abs(
      (next.location - location) / sqrt(initial.diagonal)
    ))
    diagonal.update <- max(abs(log(next.diagonal / diagonal)))
    relative.update <- max(location.update, diagonal.update)
    location <- next.location
    diagonal <- next.diagonal
  }
  final <- lwz_score_reference(data, location, diagonal, zero_tol)
  list(
    location = location,
    diagonal = diagonal,
    iterations = iteration,
    iteration.stable = stable,
    relative.update = relative.update,
    location.update = location.update,
    diagonal.update = diagonal.update,
    location.score = final$location.residual,
    diagonal.score = final$diagonal.residual,
    score.residual = max(
      final$location.residual, final$diagonal.residual
    ),
    minimum.residual.distance = minimum.radius,
    directions = final$directions,
    radii = final$radii,
    inverse.radii = 1 / final$radii,
    c.hat = mean(1 / final$radii)
  )
}


lwz_reference <- function(x, y, tol = 1e-7, max_iter = 2000L,
                          zero_tol = 0) {
  prepared <- lwz_prepare_reference(x, y)
  n1 <- nrow(x)
  n2 <- nrow(y)
  p <- ncol(x)
  fit1 <- lwz_fit_reference(
    prepared$x, tol, max_iter, zero_tol
  )
  fit2 <- lwz_fit_reference(
    prepared$y, tol, max_iter, zero_tol
  )
  c.ratio1 <- fit2$c.hat / fit1$c.hat
  c.ratio2 <- fit1$c.hat / fit2$c.hat
  bridge1 <- sqrt(fit1$diagonal / fit2$diagonal)
  bridge2 <- 1 / bridge1
  trace.A1 <- c.ratio1 * sum(bridge1)
  trace.A2 <- c.ratio2 * sum(bridge2)

  trace1.pairs <- matrix(0, n1, n1)
  trace2.pairs <- matrix(0, n2, n2)
  trace3.pairs <- matrix(NA_real_, n1, n2)
  for (k in seq_len(n1)) {
    for (ell in setdiff(seq_len(n1), k)) {
      trace1.pairs[ell, k] <- sum(
        fit1$directions[ell, ] * bridge1 * fit1$directions[k, ]
      )^2
    }
  }
  for (k in seq_len(n2)) {
    for (ell in setdiff(seq_len(n2), k)) {
      trace2.pairs[ell, k] <- sum(
        fit2$directions[ell, ] * bridge2 * fit2$directions[k, ]
      )^2
    }
  }
  for (i in seq_len(n1)) {
    for (j in seq_len(n2)) {
      trace3.pairs[i, j] <- sum(
        fit1$directions[i, ] * fit2$directions[j, ]
      )^2
    }
  }
  trace.A1.squared <- p^2 * c.ratio1^2 * sum(trace1.pairs) /
    (n1 * (n1 - 1))
  trace.A2.squared <- p^2 * c.ratio2^2 * sum(trace2.pairs) /
    (n2 * (n2 - 1))
  trace.A3 <- p^2 * sum(trace3.pairs) / (n1 * n2)

  numerator.inner <- numerator.contribution <- matrix(
    NA_real_, n1, n2
  )
  wrong.own.center <- matrix(NA_real_, n1, n2)
  for (i in seq_len(n1)) {
    for (j in seq_len(n2)) {
      sign1 <- lwz_direction_reference(
        prepared$x[i, ], fit2$location, fit1$diagonal,
        zero_tol, allow_zero = TRUE
      )$direction
      sign2 <- lwz_direction_reference(
        prepared$y[j, ], fit1$location, fit2$diagonal,
        zero_tol, allow_zero = TRUE
      )$direction
      numerator.inner[i, j] <- sum(sign1 * sign2)
      numerator.contribution[i, j] <- -numerator.inner[i, j]
      wrong.own.center[i, j] <- -sum(
        fit1$directions[i, ] * fit2$directions[j, ]
      )
    }
  }
  T <- mean(numerator.contribution)
  bias <- trace.A1 / (n1 * p) + trace.A2 / (n2 * p)
  variance.term1 <- 2 * trace.A1.squared / (n1^2 * p^2)
  variance.term2 <- 2 * trace.A2.squared / (n2^2 * p^2)
  variance.term3 <- 4 * trace.A3 / (n1 * n2 * p^2)
  variance <- variance.term1 + variance.term2 + variance.term3
  list(
    prepared = prepared,
    fit1 = fit1,
    fit2 = fit2,
    c.ratio1 = c.ratio1,
    c.ratio2 = c.ratio2,
    bridge1 = bridge1,
    bridge2 = bridge2,
    trace.A1 = trace.A1,
    trace.A2 = trace.A2,
    trace.A1.squared = trace.A1.squared,
    trace.A2.squared = trace.A2.squared,
    trace.A3 = trace.A3,
    trace1.pairs = trace1.pairs,
    trace2.pairs = trace2.pairs,
    trace3.pairs = trace3.pairs,
    numerator.inner = numerator.inner,
    numerator.contribution = numerator.contribution,
    wrong.own.center = wrong.own.center,
    T = T,
    bias = bias,
    centered = T - bias,
    variance.term1 = variance.term1,
    variance.term2 = variance.term2,
    variance.term3 = variance.term3,
    variance = variance,
    z = (T - bias) / sqrt(variance)
  )
}


lwz_fixture <- function() {
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


test_that("Li-Wang-Zou SST matches the literal feasible formulas", {
  fixture <- lwz_fixture()
  reference <- lwz_reference(
    fixture$x, fixture$y, tol = 1e-7
  )
  result <- li_wang_zou_two_sample_sign_test(
    fixture$x, fixture$y, tol = 1e-7, max_iter = 2000L
  )
  cmp <- result$components

  expect_s3_class(result, "hd_location_test")
  expect_s3_class(result, "htest")
  expect_equal(unname(result$statistic), reference$z,
               tolerance = 3e-10)
  expect_equal(unname(result$raw.statistic), reference$T,
               tolerance = 3e-12)
  expect_equal(cmp$T.SST, reference$T, tolerance = 3e-12)
  expect_equal(cmp$bias.hat, reference$bias, tolerance = 3e-11)
  expect_equal(cmp$T.minus.bias, reference$centered,
               tolerance = 3e-11)
  expect_equal(cmp$sigma2.hat, reference$variance,
               tolerance = 4e-11)
  expect_equal(unname(result$variance), reference$variance,
               tolerance = 4e-11)
  expect_equal(cmp$c1.hat, reference$fit1$c.hat,
               tolerance = 3e-11)
  expect_equal(cmp$c2.hat, reference$fit2$c.hat,
               tolerance = 3e-11)
  expect_equal(cmp$c.n1.hat, reference$c.ratio1,
               tolerance = 3e-11)
  expect_equal(cmp$c.n2.hat, reference$c.ratio2,
               tolerance = 3e-11)
  expect_equal(cmp$trace.A1.hat, reference$trace.A1,
               tolerance = 3e-10)
  expect_equal(cmp$trace.A2.hat, reference$trace.A2,
               tolerance = 3e-10)
  expect_equal(cmp$trace.A1.squared.hat,
               reference$trace.A1.squared, tolerance = 3e-9)
  expect_equal(cmp$trace.A2.squared.hat,
               reference$trace.A2.squared, tolerance = 3e-9)
  expect_equal(cmp$trace.A3tA3.hat, reference$trace.A3,
               tolerance = 3e-9)
  expect_equal(cmp$variance.term1, reference$variance.term1,
               tolerance = 4e-11)
  expect_equal(cmp$variance.term2, reference$variance.term2,
               tolerance = 4e-11)
  expect_equal(cmp$variance.term3, reference$variance.term3,
               tolerance = 4e-11)
  expect_equal(unname(cmp$numerator.inner.product),
               reference$numerator.inner, tolerance = 3e-11)
  expect_equal(unname(cmp$numerator.contribution),
               reference$numerator.contribution, tolerance = 3e-11)
  expect_equal(unname(cmp$trace.A1.ordered.pair.squared),
               reference$trace1.pairs, tolerance = 3e-10)
  expect_equal(unname(cmp$trace.A2.ordered.pair.squared),
               reference$trace2.pairs, tolerance = 3e-10)
  expect_equal(unname(cmp$trace.A3.cross.pair.squared),
               reference$trace3.pairs, tolerance = 3e-10)
  expect_equal(unname(cmp$bridge.A1.diagonal.standardized),
               reference$bridge1, tolerance = 3e-11)
  expect_equal(unname(cmp$bridge.A2.diagonal.standardized),
               reference$bridge2, tolerance = 3e-11)
  expect_equal(unname(cmp$full.fit$group1$location.standardized),
               reference$fit1$location, tolerance = 3e-10)
  expect_equal(unname(cmp$full.fit$group2$location.standardized),
               reference$fit2$location, tolerance = 3e-10)
  expect_equal(unname(cmp$own.direction1),
               reference$fit1$directions, tolerance = 3e-10)
  expect_equal(unname(cmp$own.direction2),
               reference$fit2$directions, tolerance = 3e-10)
  expect_equal(result$p.value,
               stats::pnorm(reference$z, lower.tail = FALSE),
               tolerance = 1e-15)
})


test_that("bias, ordered traces, n-squared variance, and upper tail are locked", {
  fixture <- lwz_fixture()
  result <- li_wang_zou_two_sample_sign_test(
    fixture$x, fixture$y, alpha = 0.1,
    tol = 1e-7, max_iter = 2000L
  )
  cmp <- result$components
  n1 <- nrow(fixture$x)
  n2 <- nrow(fixture$y)
  p <- ncol(fixture$x)

  expect_equal(cmp$ordered.pair.count1, n1 * (n1 - 1))
  expect_equal(cmp$ordered.pair.count2, n2 * (n2 - 1))
  expect_equal(cmp$cross.pair.count, n1 * n2)
  expect_equal(unname(diag(cmp$trace.A1.ordered.pair.squared)),
               numeric(n1))
  expect_equal(unname(diag(cmp$trace.A2.ordered.pair.squared)),
               numeric(n2))
  expect_equal(
    cmp$trace.A1.squared.hat,
    p^2 * cmp$c.n1.hat^2 *
      sum(cmp$trace.A1.ordered.pair.squared) / (n1 * (n1 - 1)),
    tolerance = 3e-14
  )
  expect_equal(
    cmp$trace.A2.squared.hat,
    p^2 * cmp$c.n2.hat^2 *
      sum(cmp$trace.A2.ordered.pair.squared) / (n2 * (n2 - 1)),
    tolerance = 3e-14
  )
  expect_equal(
    cmp$trace.A3tA3.hat,
    p^2 * sum(cmp$trace.A3.cross.pair.squared) / (n1 * n2),
    tolerance = 3e-14
  )
  expect_equal(
    cmp$bias.hat,
    cmp$trace.A1.hat / (n1 * p) + cmp$trace.A2.hat / (n2 * p),
    tolerance = 3e-15
  )
  expect_equal(cmp$c.n1.hat, cmp$c2.hat / cmp$c1.hat,
               tolerance = 2e-15)
  expect_equal(cmp$c.n2.hat, cmp$c1.hat / cmp$c2.hat,
               tolerance = 2e-15)
  expect_equal(
    cmp$trace.A1.hat,
    cmp$c.n1.hat * sum(cmp$bridge.A1.diagonal.standardized),
    tolerance = 3e-15
  )
  expect_equal(
    cmp$trace.A2.hat,
    cmp$c.n2.hat * sum(cmp$bridge.A2.diagonal.standardized),
    tolerance = 3e-15
  )
  expect_equal(cmp$T.minus.bias, cmp$T.SST - cmp$bias.hat,
               tolerance = 3e-15)
  expect_equal(
    cmp$variance.term1,
    2 * cmp$trace.A1.squared.hat / (n1^2 * p^2),
    tolerance = 3e-15
  )
  expect_equal(
    cmp$variance.term2,
    2 * cmp$trace.A2.squared.hat / (n2^2 * p^2),
    tolerance = 3e-15
  )
  wrong1 <- 2 * cmp$trace.A1.squared.hat /
    (n1 * (n1 - 1) * p^2)
  wrong2 <- 2 * cmp$trace.A2.squared.hat /
    (n2 * (n2 - 1) * p^2)
  expect_gt(abs(cmp$variance.term1 - wrong1), 1e-6)
  expect_gt(abs(cmp$variance.term2 - wrong2), 1e-6)
  expect_equal(cmp$sigma2.hat,
               cmp$variance.term1 + cmp$variance.term2 +
                 cmp$variance.term3,
               tolerance = 3e-15)
  expect_equal(unname(result$statistic),
               cmp$T.minus.bias / cmp$sigma.hat,
               tolerance = 3e-15)
  expect_equal(result$p.value,
               stats::pnorm(unname(result$statistic), lower.tail = FALSE),
               tolerance = 0)
  expect_identical(result$null.distribution$tail, "upper")
  expect_identical(result$alternative, "two.sided")
  expect_identical(
    result$diagnostics$rejection$reject,
    isTRUE(unname(result$statistic) > stats::qnorm(0.9))
  )
  expect_match(result$diagnostics$variance.denominators, "n1\\^2")
  expect_match(result$diagnostics$trace.denominators, "ordered distinct")
})


test_that("SST uses full-sample crossed locations and no leave-out fit", {
  fixture <- lwz_fixture()
  reference <- lwz_reference(fixture$x, fixture$y, tol = 1e-7)
  result <- li_wang_zou_two_sample_sign_test(
    fixture$x, fixture$y, tol = 1e-7, max_iter = 2000L
  )
  cmp <- result$components

  expect_equal(unname(cmp$numerator.contribution),
               reference$numerator.contribution, tolerance = 3e-11)
  expect_gt(max(abs(reference$wrong.own.center -
                      reference$numerator.contribution)), 1e-3)
  expect_match(result$diagnostics$fit.structure,
               "exactly two full-sample.*no leave-out")
  expect_match(result$diagnostics$numerator.structure,
               "group-2 full-sample location")
  expect_false("leaveout.fit" %in% names(cmp))
  expect_length(cmp$full.fit, 2L)
  expect_equal(cmp$c1.hat, mean(cmp$own.inverse.radius1),
               tolerance = 2e-15)
  expect_equal(cmp$c2.hat, mean(cmp$own.inverse.radius2),
               tolerance = 2e-15)
  expect_identical(result$diagnostics$regularization, "none")
  expect_identical(result$diagnostics$pseudoinverse, "none")
  expect_identical(result$diagnostics$bias.repair, "none")
  expect_identical(result$diagnostics$variance.repair, "none")
})


test_that("SST is permutation, group-swap, and diagonal-scale invariant", {
  fixture <- lwz_fixture()
  x <- fixture$x
  y <- fixture$y
  baseline <- li_wang_zou_two_sample_sign_test(
    x, y, tol = 1e-7, max_iter = 2000L
  )
  permuted <- li_wang_zou_two_sample_sign_test(
    x[c(4, 1, 7, 2, 6, 3, 5), ],
    y[c(3, 8, 1, 7, 2, 6, 4, 5), ],
    tol = 1e-7, max_iter = 2000L
  )
  swapped <- li_wang_zou_two_sample_sign_test(
    y, x, tol = 1e-7, max_iter = 2000L
  )

  expect_equal(permuted$statistic, baseline$statistic,
               tolerance = 4e-8)
  expect_equal(permuted$raw.statistic, baseline$raw.statistic,
               tolerance = 4e-9)
  expect_equal(permuted$components$bias.hat,
               baseline$components$bias.hat, tolerance = 4e-9)
  expect_equal(permuted$variance, baseline$variance,
               tolerance = 4e-9)
  expect_equal(swapped$statistic, baseline$statistic,
               tolerance = 4e-8)
  expect_equal(swapped$raw.statistic, baseline$raw.statistic,
               tolerance = 4e-9)
  expect_equal(swapped$components$bias.hat,
               baseline$components$bias.hat, tolerance = 4e-9)
  expect_equal(swapped$variance, baseline$variance,
               tolerance = 4e-9)
  expect_equal(swapped$components$trace.A1.hat,
               baseline$components$trace.A2.hat, tolerance = 4e-8)
  expect_equal(swapped$components$trace.A2.hat,
               baseline$components$trace.A1.hat, tolerance = 4e-8)

  coordinate.scale <- c(-0.08, 0.4, -3.5, 120)
  shift <- c(1e6, -2e6, 3e6, -4e6)
  transformed <- li_wang_zou_two_sample_sign_test(
    sweep(sweep(x, 2L, coordinate.scale, "*"), 2L, shift, "+"),
    sweep(sweep(y, 2L, coordinate.scale, "*"), 2L, shift, "+"),
    tol = 1e-7, max_iter = 2000L
  )
  expect_equal(transformed$statistic, baseline$statistic,
               tolerance = 3e-6)
  expect_equal(transformed$raw.statistic, baseline$raw.statistic,
               tolerance = 3e-7)
  expect_equal(transformed$components$bias.hat,
               baseline$components$bias.hat, tolerance = 3e-7)
  expect_equal(transformed$variance, baseline$variance,
               tolerance = 3e-7)
})


test_that("fit convergence and singular contracts are explicit", {
  fixture <- lwz_fixture()
  stable <- li_wang_zou_two_sample_sign_test(
    fixture$x, fixture$y, tol = 1e-7, max_iter = 2000L
  )
  expect_true(stable$diagnostics$iteration.stable)
  expect_lte(stable$diagnostics$score.residual, 1e-7)
  expect_gt(stable$diagnostics$minimum.residual.distance, 0)
  expect_true(stable$components$full.fit$group1$diagnostics$iteration.stable)
  expect_true(stable$components$full.fit$group2$diagnostics$iteration.stable)

  expect_error(
    li_wang_zou_two_sample_sign_test(
      fixture$x, fixture$y, max_iter = 1L, strict = TRUE
    ),
    "did not satisfy the estimating equations"
  )
  expect_warning(
    unstable <- li_wang_zou_two_sample_sign_test(
      fixture$x, fixture$y, max_iter = 1L, strict = FALSE
    ),
    "did not satisfy the estimating equations"
  )
  expect_false(unstable$diagnostics$iteration.stable)
  expect_gt(unstable$diagnostics$stability.failures, 0L)

  minimum <- li_wang_zou_two_sample_sign_test(
    matrix(c(-1, 1), ncol = 1L),
    matrix(c(-2, 2), ncol = 1L)
  )
  expect_true(is.finite(minimum$statistic))
  expect_equal(minimum$n, c(n1 = 2L, n2 = 2L))
  expect_equal(minimum$diagnostics$relative.update, 0)

  expect_error(
    li_wang_zou_two_sample_sign_test(
      matrix(1, nrow = 1L), matrix(c(-1, 1), ncol = 1L)
    ),
    "at least 2 row"
  )
  expect_error(
    li_wang_zou_two_sample_sign_test(
      matrix(stats::rnorm(12), 4, 3), matrix(stats::rnorm(16), 4, 4)
    ),
    "same number of columns"
  )
  constant <- matrix(stats::rnorm(24), 8, 3)
  constant[, 2] <- 1
  expect_error(
    li_wang_zou_two_sample_sign_test(
      constant, matrix(stats::rnorm(27), 9, 3)
    ),
    "group-specific marginal sample variance.*No ridge"
  )
  expect_error(
    li_wang_zou_two_sample_sign_test(
      matrix(c(-1, 0, 1), ncol = 1L),
      matrix(c(-2, 2), ncol = 1L)
    ),
    "training observation.*no perturbation"
  )
  expect_error(
    li_wang_zou_two_sample_sign_test(
      fixture$x, fixture$y, zero_tol = 100
    ),
    "at or below `zero_tol`"
  )
  expect_error(
    li_wang_zou_two_sample_sign_test(
      fixture$x, fixture$y, tol = 0
    ),
    "finite positive"
  )
  expect_error(
    li_wang_zou_two_sample_sign_test(
      fixture$x, fixture$y, max_iter = 0
    ),
    "positive integer"
  )
  expect_error(
    li_wang_zou_two_sample_sign_test(
      fixture$x, fixture$y, strict = NA
    ),
    "TRUE.*FALSE"
  )
})


test_that("safe scaling handles extreme finite coordinate units", {
  fixture <- lwz_fixture()
  baseline <- li_wang_zou_two_sample_sign_test(
    fixture$x, fixture$y, tol = 1e-7, max_iter = 2000L
  )
  huge <- li_wang_zou_two_sample_sign_test(
    fixture$x * 1e308, fixture$y * 1e308,
    tol = 1e-7, max_iter = 2000L
  )
  expect_gte(huge$diagnostics$subtraction.overflow.fallback.columns, 1)
  expect_true(is.finite(huge$statistic))
  expect_true(is.finite(huge$p.value))
  expect_true(all(is.finite(
    huge$components$full.fit$group1$location
  )))
  expect_true(all(is.finite(
    huge$components$full.fit$group2$location
  )))
  expect_equal(huge$statistic, baseline$statistic, tolerance = 5e-8)
  expect_equal(huge$raw.statistic, baseline$raw.statistic,
               tolerance = 5e-9)
  expect_equal(huge$components$bias.hat,
               baseline$components$bias.hat, tolerance = 5e-9)
  expect_equal(huge$variance, baseline$variance, tolerance = 5e-9)
})


test_that("SST result is self-consistent and preserves variable names", {
  fixture <- lwz_fixture()
  variables <- c("gene_a", "gene_b", "gene_c", "gene_d")
  colnames(fixture$x) <- variables
  colnames(fixture$y) <- variables
  result <- li_wang_zou_two_sample_sign_test(
    as.data.frame(fixture$x), as.data.frame(fixture$y),
    alpha = 0.1, tol = 1e-7, max_iter = 2000L
  )
  cmp <- result$components

  expect_identical(names(result$estimate), variables)
  expect_identical(names(result$null.value), variables)
  expect_identical(names(cmp$full.fit$group1$location), variables)
  expect_identical(names(cmp$full.fit$group2$scale.diagonal), variables)
  expect_identical(colnames(cmp$own.direction1), variables)
  expect_identical(colnames(cmp$own.direction2), variables)
  expect_equal(unname(result$statistic),
               cmp$T.minus.bias / cmp$sigma.hat,
               tolerance = 2e-15)
  expect_equal(unname(result$variance), cmp$sigma2.hat, tolerance = 0)
  expect_identical(result$diagnostics$scientific.alternative,
                   "two-sided location inequality")
  expect_identical(result$diagnostics$calibration.tail,
                   "upper standard-normal tail for quadratic SST")
  expect_identical(result$diagnostics$radius.perturbation, "none")
  expect_identical(result$diagnostics$weight.cap, "none")
  expect_identical(result$diagnostics$variance.repair, "none")
})
