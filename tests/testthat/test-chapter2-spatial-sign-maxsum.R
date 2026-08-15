ssmax_fixture <- function() {
  rbind(
    c(0.7, -1.1, 0.2, 1.3),
    c(-0.4, 0.8, 1.2, -0.6),
    c(1.1, 0.3, -0.9, 0.5),
    c(-1.2, -0.5, 0.7, 1.0),
    c(0.2, 1.4, -0.4, -1.1),
    c(1.5, -0.2, 0.5, 0.1),
    c(-0.8, 0.6, -1.3, 0.9)
  )
}


ssmax_normalize_diagonal_reference <- function(diagonal) {
  exp(log(diagonal) - mean(log(diagonal)))
}


ssmax_score_reference <- function(data, theta, diagonal) {
  residual <- sweep(data, 2L, theta, "-")
  residual <- sweep(residual, 2L, sqrt(diagonal), "/")
  radius <- sqrt(rowSums(residual^2))
  if (any(radius == 0)) {
    stop("coincident residual in reference")
  }
  direction <- residual / radius
  list(
    direction = direction,
    radius = radius,
    inverse.radius = 1 / radius,
    location.residual = max(abs(colMeans(direction))),
    diagonal.residual = max(abs(
      ncol(data) * colMeans(direction^2) - 1
    )),
    direction.sum = colSums(direction),
    direction.square.sum = colSums(direction^2)
  )
}


ssmax_fit_reference <- function(x, origin, tol = 1e-9,
                                max_iter = 1000L) {
  residual <- sweep(x, 2L, origin, "-")
  scale <- apply(abs(residual), 2L, max)
  data <- sweep(residual, 2L, scale, "/")
  theta <- colMeans(data)
  diagonal <- ssmax_normalize_diagonal_reference(
    apply(data, 2L, stats::var)
  )
  initial.diagonal <- diagonal
  relative.update <- location.update <- diagonal.update <- Inf

  for (iteration in seq_len(max_iter)) {
    score <- ssmax_score_reference(data, theta, diagonal)
    theta.new <- theta + sqrt(diagonal) * score$direction.sum /
      sum(score$inverse.radius)
    diagonal.new <- diagonal * ncol(data) *
      score$direction.square.sum / nrow(data)
    diagonal.new <- ssmax_normalize_diagonal_reference(diagonal.new)
    location.update <- max(abs(
      (theta.new - theta) / sqrt(initial.diagonal)
    ))
    diagonal.update <- max(abs(log(diagonal.new / diagonal)))
    relative.update <- max(location.update, diagonal.update)
    theta <- theta.new
    diagonal <- diagonal.new

    if (relative.update <= tol) {
      score <- ssmax_score_reference(data, theta, diagonal)
      score.residual <- max(
        score$location.residual, score$diagonal.residual
      )
      raw.diagonal <- scale^2 * diagonal
      return(list(
        origin = origin,
        scale = scale,
        data = data,
        theta = theta,
        diagonal = diagonal,
        location = origin + scale * theta,
        diagonal.input = raw.diagonal / max(raw.diagonal),
        score = score,
        score.residual = score.residual,
        iterations = iteration,
        relative.update = relative.update,
        location.update = location.update,
        diagonal.update = diagonal.update
      ))
    }
    if (iteration == max_iter) {
      stop("reference fit did not converge")
    }
  }
  stop("unreachable")
}


test_that("scaled spatial median matches the literal full-sample HR recursion", {
  x <- ssmax_fixture()
  colnames(x) <- paste0("feature", seq_len(ncol(x)))
  reference <- ssmax_fit_reference(x, colMeans(x), tol = 1e-9)
  result <- scaled_spatial_median(x, tol = 1e-9)

  expect_s3_class(result, "scaled_spatial_median")
  expect_equal(unname(result$location), unname(reference$location),
               tolerance = 2e-12)
  expect_equal(unname(result$scale.diagonal),
               unname(reference$diagonal.input),
               tolerance = 3e-12)
  expect_equal(unname(result$scale.diagonal.standardized),
               unname(reference$diagonal), tolerance = 3e-12)
  expect_equal(unname(result$directions),
               unname(reference$score$direction),
               tolerance = 3e-12)
  expect_equal(unname(result$radii.standardized),
               unname(reference$score$radius),
               tolerance = 3e-12)
  expect_equal(unname(result$inverse.radii.standardized),
               unname(reference$score$inverse.radius), tolerance = 3e-12)
  expect_equal(result$inverse.radius.moment.standardized,
               mean(reference$score$inverse.radius), tolerance = 3e-12)
  expect_lt(abs(
    unname(result$diagnostics$location.equation.residual) -
      unname(reference$score$location.residual)
  ), 1e-14)
  expect_lt(abs(
    unname(result$diagnostics$diagonal.equation.residual) -
      unname(reference$score$diagonal.residual)
  ), 1e-14)
  expect_lt(abs(
    unname(result$diagnostics$score.residual) -
      unname(reference$score.residual)
  ), 1e-14)
  expect_equal(result$diagnostics$iterations, reference$iterations)
  expect_lt(abs(
    unname(result$diagnostics$relative.update) -
      unname(reference$relative.update)
  ), 1e-14)
  expect_equal(prod(result$scale.diagonal.standardized), 1,
               tolerance = 2e-14)
  expect_equal(max(result$scale.diagonal), 1, tolerance = 2e-14)
  expect_equal(unname(result$log.scale.diagonal),
               log(unname(result$scale.diagonal)), tolerance = 2e-14)
  expect_true(result$diagnostics$iteration.stable)
  expect_identical(result$diagnostics$regularization, "none")
  expect_identical(result$diagnostics$numerical.floor, "none")
  expect_identical(result$diagnostics$perturbation, "none")
  expect_identical(result$diagnostics$pseudoinverse, "none")
})


test_that("max statistic locks every primary-paper factor and Gumbel tail", {
  x <- ssmax_fixture()
  mu <- c(0.1, -0.15, 0.05, 0.2)
  reference <- ssmax_fit_reference(x, mu, tol = 1e-9)
  result <- spatial_sign_max_test(x, mu, alpha = 0.05, tol = 1e-9)
  component <- result$components
  n <- nrow(x)
  p <- ncol(x)
  zeta1 <- mean(reference$score$inverse.radius)
  standardized.location <- reference$theta / sqrt(reference$diagonal)
  multiplier <- n * p * zeta1^2 * (1 - n^(-0.5))
  coordinate <- multiplier * standardized.location^2
  T.max <- max(coordinate)
  centered <- T.max - 2 * log(p) + log(log(p))
  p.max <- -expm1(-pi^(-0.5) * exp(-centered / 2))

  expect_s3_class(result, "hd_location_test")
  expect_s3_class(result, "htest")
  expect_equal(unname(result$raw.statistic), T.max, tolerance = 4e-12)
  expect_equal(component$T.MAX, T.max, tolerance = 4e-12)
  expect_equal(unname(component$coordinate.statistic), coordinate,
               tolerance = 4e-12)
  expect_equal(component$zeta1.hat.standardized, zeta1,
               tolerance = 3e-12)
  expect_equal(component$finite.sample.factor, 1 - n^(-0.5),
               tolerance = 1e-15)
  expect_equal(
    component$T.MAX / max(standardized.location^2),
    n * p * component$zeta1.hat.standardized^2 *
      component$finite.sample.factor,
    tolerance = 3e-12
  )
  expect_equal(unname(result$statistic), centered, tolerance = 4e-12)
  expect_equal(result$p.value, p.max, tolerance = 2e-15)
  expect_equal(component$log.p.MAX, log(p.max), tolerance = 2e-15)
  expect_equal(component$log.one.minus.p.MAX,
               -pi^(-0.5) * exp(-centered / 2), tolerance = 2e-15)
  expect_equal(.ssmax_gumbel_tail(0)$p.value,
               0.4311790581359798, tolerance = 2e-15)
  expect_equal(.ssmax_gumbel_quantile(0.95),
               4.7956606122349275, tolerance = 2e-15)
  expect_equal(
    component$critical.value.raw,
    2 * log(p) - log(log(p)) + .ssmax_gumbel_quantile(0.95),
    tolerance = 2e-15
  )
  expect_identical(result$null.distribution$tail, "upper")
  expect_identical(result$alternative, "two.sided")
  expect_identical(result$diagnostics$regularization, "none")
})


test_that("max test is permutation, translation, and diagonal-scale invariant", {
  x <- ssmax_fixture()
  mu <- c(0.1, -0.15, 0.05, 0.2)
  baseline <- spatial_sign_max_test(x, mu, tol = 1e-8)

  row.permutation <- c(7, 2, 5, 1, 6, 3, 4)
  row.result <- spatial_sign_max_test(
    x[row.permutation, , drop = FALSE], mu, tol = 1e-8
  )
  expect_equal(row.result$raw.statistic, baseline$raw.statistic,
               tolerance = 2e-12)
  expect_equal(row.result$p.value, baseline$p.value,
               tolerance = 2e-13)

  column.permutation <- c(4, 2, 1, 3)
  column.result <- spatial_sign_max_test(
    x[, column.permutation, drop = FALSE], mu[column.permutation],
    tol = 1e-8
  )
  expect_equal(column.result$raw.statistic, baseline$raw.statistic,
               tolerance = 2e-12)
  expect_equal(column.result$p.value, baseline$p.value,
               tolerance = 2e-13)

  multiplier <- c(-3, 0.25, 2, -5)
  scaled <- sweep(x, 2L, multiplier, "*")
  scaled.result <- spatial_sign_max_test(
    scaled, mu * multiplier, tol = 1e-8
  )
  expected.diagonal <- baseline$components$scale.diagonal * multiplier^2
  expected.diagonal <- expected.diagonal / max(expected.diagonal)
  expect_equal(scaled.result$raw.statistic, baseline$raw.statistic,
               tolerance = 5e-12)
  expect_equal(scaled.result$p.value, baseline$p.value,
               tolerance = 3e-13)
  expect_equal(unname(scaled.result$components$scale.diagonal),
               unname(expected.diagonal), tolerance = 5e-12)

  shift <- c(1e8, -2e8, 3e8, -4e8)
  translated <- sweep(x, 2L, shift, "+")
  translated.result <- spatial_sign_max_test(
    translated, mu + shift, tol = 1e-8
  )
  expect_equal(translated.result$raw.statistic,
               baseline$raw.statistic, tolerance = 2e-7)
  expect_equal(translated.result$p.value, baseline$p.value,
               tolerance = 5e-8)
  expect_equal(
    unname(translated.result$components$location.minus.null),
    unname(baseline$components$location.minus.null),
    tolerance = 3e-7
  )
})


test_that("stable Cauchy algebra agrees with ordinary formulas and endpoints", {
  for (u in c(0.1, 0.5, 0.9, 1e-100)) {
    equal <- .ssmax_cauchy_combine_logtails(
      log(u), log1p(-u), log(u), log1p(-u)
    )
    expect_equal(equal$p.value, u, tolerance = 3e-15)
    expect_equal(equal$log.p.value, log(u), tolerance = 3e-14)
  }

  complementary <- .ssmax_cauchy_combine_logtails(
    log(0.13), log(0.87), log(0.87), log(0.13)
  )
  expect_equal(complementary$p.value, 0.5, tolerance = 2e-15)
  expect_equal(complementary$statistic, 0, tolerance = 2e-15)

  p1 <- 0.031
  p2 <- 0.42
  ordinary <- atan2(
    2 * sin(pi * p1) * sin(pi * p2),
    sin(pi * (p1 + p2))
  ) / pi
  stable <- .ssmax_cauchy_combine_logtails(
    log(p1), log1p(-p1), log(p2), log1p(-p2)
  )
  expect_equal(stable$p.value, ordinary, tolerance = 2e-15)

  log.extreme <- .ssmax_cauchy_combine_logtails(
    -1000, 0, -1000, 0
  )
  expect_equal(log.extreme$p.value, 0)
  expect_equal(log.extreme$log.p.value, -1000, tolerance = 2e-13)
  expect_true(is.finite(log.extreme$log.absolute))

  zero <- .ssmax_cauchy_combine_logtails(-Inf, 0, -Inf, 0)
  one <- .ssmax_cauchy_combine_logtails(0, -Inf, 0, -Inf)
  expect_identical(zero$p.value, 0)
  expect_identical(one$p.value, 1)
  expect_error(
    .ssmax_cauchy_combine_logtails(-Inf, 0, 0, -Inf),
    "indeterminate"
  )
  expect_error(
    .ssmax_cauchy_combine_logtails(0, -Inf, -Inf, 0),
    "indeterminate"
  )
})


test_that("max-sum directly uses the feasible Feng-Sun component", {
  x <- ssmax_fixture()
  mu <- c(0.1, -0.15, 0.05, 0.2)
  max.reference <- spatial_sign_max_test(x, mu, tol = 1e-6)
  sum.reference <- feng_sun_one_sample_test(x, mu, tol = 1e-6)
  result <- spatial_sign_maxsum_test(x, mu, tol = 1e-6)
  component <- result$components

  expect_s3_class(result, "hd_location_test")
  expect_equal(component$max.test$raw.statistic,
               max.reference$raw.statistic, tolerance = 2e-13)
  expect_equal(component$max.test$p.value,
               max.reference$p.value, tolerance = 2e-15)
  expect_equal(component$sum.test$raw.statistic,
               sum.reference$raw.statistic, tolerance = 2e-13)
  expect_equal(component$sum.test$statistic,
               sum.reference$statistic, tolerance = 2e-13)
  expect_equal(component$p.SUM,
               stats::pnorm(unname(sum.reference$statistic),
                            lower.tail = FALSE), tolerance = 2e-15)

  ordinary.statistic <-
    0.5 * tan(pi * (0.5 - component$p.MAX)) +
    0.5 * tan(pi * (0.5 - component$p.SUM))
  ordinary.p <- stats::pcauchy(
    ordinary.statistic, lower.tail = FALSE
  )
  expect_equal(component$Cauchy.statistic, ordinary.statistic,
               tolerance = 2e-14)
  expect_equal(result$p.value, ordinary.p, tolerance = 2e-15)
  expect_equal(unname(result$statistic),
               atan(ordinary.statistic), tolerance = 2e-15)
  expect_equal(component$log.p.Cauchy, log(result$p.value),
               tolerance = 2e-15)
  expect_identical(result$diagnostics$probability.clipping, "none")
  expect_true(result$diagnostics$sum.leaveout.fits.converged)
})


test_that("strict and non-strict contracts expose the last finite iterate", {
  x <- ssmax_fixture()
  expect_error(
    scaled_spatial_median(x, tol = 1e-15, max_iter = 1L),
    "did not stabilize"
  )
  expect_warning(
    last <- scaled_spatial_median(
      x, tol = 1e-15, max_iter = 1L, strict = FALSE
    ),
    "last finite iterate"
  )
  expect_false(last$diagnostics$iteration.stable)
  expect_identical(last$diagnostics$strict, FALSE)
  expect_true(is.finite(last$diagnostics$relative.update))
  expect_true(is.finite(last$diagnostics$score.residual))

  expect_error(
    spatial_sign_max_test(x, tol = 1e-15, max_iter = 1L),
    "did not stabilize"
  )
  expect_warning(
    max.last <- spatial_sign_max_test(
      x, tol = 1e-15, max_iter = 1L, strict = FALSE
    ),
    "last finite iterate"
  )
  expect_false(max.last$diagnostics$iteration.stable)
  expect_true(is.finite(max.last$statistic))
  expect_true(is.finite(max.last$p.value))

  slow.score <- matrix(
    c(-2, -1, 0, 1, 2, 3, -1, 2, 1, -2, 3, 0), 6, 2
  )
  update.stable <- scaled_spatial_median(slow.score, tol = 1e-6)
  expect_true(update.stable$diagnostics$iteration.stable)
  expect_lte(update.stable$diagnostics$relative.update, 1e-6)
  expect_gt(update.stable$diagnostics$score.residual, 1e-6)
  expect_match(update.stable$diagnostics$score.residual.role,
               "separately from iterate stability")
})


test_that("degeneracy and input boundaries fail without numerical repairs", {
  expect_error(scaled_spatial_median(matrix(1:4, 1, 4)),
               "at least 2 row")
  expect_error(spatial_sign_max_test(matrix(1:8, 4, 2), mu = 1:3),
               "length 2")
  expect_error(spatial_sign_max_test(matrix(1:8, 4, 2), alpha = 0),
               "strictly between")
  expect_error(spatial_sign_max_test(matrix(1:8, 4, 2), alpha = 1),
               "strictly between")
  expect_error(spatial_sign_max_test(matrix(1:8, 4, 2), strict = NA),
               "TRUE or FALSE")
  expect_error(spatial_sign_max_test(matrix(1:8, 4, 2), tol = 0),
               "positive")
  expect_error(spatial_sign_max_test(matrix(1:8, 4, 2), max_iter = 0),
               "positive integer")
  expect_error(spatial_sign_max_test(matrix(1:8, 4, 2), zero_tol = -1),
               "non-negative")
  expect_error(spatial_sign_max_test(matrix(1:5, 5, 1)),
               "at least two variables")
  expect_error(spatial_sign_maxsum_test(matrix(1:9, 3, 3)),
               "at least 4 row")

  degenerate <- cbind(seq_len(6), rep(1, 6), seq_len(6)^2)
  expect_error(scaled_spatial_median(degenerate), "positive variation")

  coincident <- rbind(
    c(0, 0), c(1, 1), c(-1, -1), c(1, -1), c(-1, 1)
  )
  expect_error(scaled_spatial_median(coincident), "undefined")
  expect_error(scaled_spatial_median(ssmax_fixture(), zero_tol = 10),
               "zero_tol")
})


test_that("overflow-safe target subtraction remains finite and diagnosed", {
  scale <- 1.6e308
  x <- rbind(
    c(scale, scale), c(-scale, scale),
    c(scale, -scale), c(-scale, -scale),
    c(0.5 * scale, 0.2 * scale),
    c(-0.4 * scale, -0.3 * scale)
  )
  result <- spatial_sign_max_test(
    x, c(scale, -scale), tol = 1e-6
  )

  expect_s3_class(result, "hd_location_test")
  expect_true(is.finite(result$statistic))
  expect_true(is.finite(result$p.value))
  expect_true(all(is.finite(result$estimate)))
  expect_equal(
    result$diagnostics$subtraction.overflow.fallback.columns, 2
  )
  expect_identical(result$diagnostics$regularization, "none")
})
