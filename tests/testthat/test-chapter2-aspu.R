aspu_double_factorial_reference <- function(value) {
  if (value <= 0) {
    return(1)
  }
  prod(seq.int(value, 1, by = -2))
}

aspu_normal_moment_reference <- function(power, variance = 1) {
  if (power %% 2L) {
    0
  } else {
    aspu_double_factorial_reference(power - 1L) * variance^(power / 2)
  }
}

aspu_bivariate_moment_reference <- function(first, second, variance.first,
                                             variance.second, covariance) {
  if ((first + second) %% 2L) {
    return(0)
  }
  answer <- 0
  for (cross in 0:min(first, second)) {
    if ((first - cross) %% 2L || (second - cross) %% 2L) {
      next
    }
    answer <- answer + choose(first, cross) * choose(second, cross) *
      factorial(cross) *
      aspu_double_factorial_reference(first - cross - 1L) *
      aspu_double_factorial_reference(second - cross - 1L) *
      covariance^cross * variance.first^((first - cross) / 2) *
      variance.second^((second - cross) / 2)
  }
  answer
}

aspu_power_reference <- function(scores, covariance, powers) {
  observed <- vapply(
    powers, function(power) sum(scores^power), numeric(1L)
  )
  marginal <- matrix(
    unlist(lapply(powers, function(power) {
      vapply(diag(covariance), function(variance) {
        aspu_normal_moment_reference(power, variance)
      }, numeric(1L))
    }), use.names = FALSE),
    nrow = length(scores), ncol = length(powers)
  )
  null.mean <- colSums(marginal)
  null.covariance <- matrix(0, length(powers), length(powers))
  for (first in seq_along(powers)) {
    for (second in seq_along(powers)) {
      for (j in seq_along(scores)) {
        for (k in seq_along(scores)) {
          null.covariance[first, second] <-
            null.covariance[first, second] +
            aspu_bivariate_moment_reference(
              powers[[first]], powers[[second]],
              covariance[j, j], covariance[k, k], covariance[j, k]
            ) - marginal[j, first] * marginal[k, second]
        }
      }
    }
  }
  null.variance <- diag(null.covariance)
  standardized <- (observed - null.mean) / sqrt(null.variance)
  list(
    observed = observed,
    null.mean = null.mean,
    null.variance = null.variance,
    null.covariance = null.covariance,
    null.correlation = stats::cov2cor(null.covariance),
    standardized = standardized
  )
}

aspu_group_reference <- function(z, correlation, tail, steps = 128L) {
  if (length(z) == 1L) {
    if (identical(tail, "two.sided")) {
      return(2 * stats::pnorm(abs(z), lower.tail = FALSE))
    }
    return(stats::pnorm(z, lower.tail = FALSE))
  }
  threshold <- if (identical(tail, "two.sided")) max(abs(z)) else max(z)
  lower <- if (identical(tail, "two.sided")) {
    rep(-threshold, length(z))
  } else {
    rep(-Inf, length(z))
  }
  1 - as.numeric(mvtnorm::pmvnorm(
    lower = lower, upper = rep(threshold, length(z)),
    corr = correlation,
    algorithm = mvtnorm::Miwa(steps = steps)
  ))
}


test_that("book-studentized aSPU matches literal standardisation and moments", {
  x <- cbind(
    c(-1.2, 0.3, 1.7, -0.8, 2.1, 0.5, -0.2),
    c(0.4, -1.1, 0.8, 2.0, -0.5, 1.4, 0.1),
    c(1.5, -0.7, 0.2, 1.1, -1.4, 0.6, 2.2),
    c(-0.2, 1.3, -1.7, 0.9, 0.5, 2.4, -0.8)
  )
  y <- cbind(
    c(-0.6, 1.1, -1.4, 0.8, 2.5, -0.1, 1.7, 0.3, -0.9),
    c(1.8, -0.3, 0.7, -1.2, 0.5, 2.1, -0.8, 1.0, 0.2),
    c(-1.0, 0.4, 1.6, -0.2, 0.9, 2.3, -0.6, 1.2, 0.1),
    c(0.7, -1.5, 0.3, 1.9, -0.4, 0.8, 2.0, -0.7, 1.1)
  )
  powers <- c(1:6, Inf)
  finite <- powers[is.finite(powers)]
  n1 <- nrow(x)
  n2 <- nrow(y)
  pooled <- ((n1 - 1) * stats::cov(x) +
               (n2 - 1) * stats::cov(y)) / (n1 + n2 - 2)
  standard.error <- sqrt(diag(pooled) * (1 / n1 + 1 / n2))
  w <- (colMeans(x) - colMeans(y)) / standard.error
  correlation <- stats::cov2cor(pooled)
  reference <- aspu_power_reference(w, correlation, finite)
  odd <- finite %% 2L == 1L
  odd.p <- aspu_group_reference(
    reference$standardized[odd],
    reference$null.correlation[odd, odd, drop = FALSE], "two.sided"
  )
  even.p <- aspu_group_reference(
    reference$standardized[!odd],
    reference$null.correlation[!odd, !odd, drop = FALSE], "upper"
  )
  maximum <- max(abs(w))
  g <- maximum^2 - 2 * log(ncol(x)) + log(log(ncol(x)))
  infinity.p <- 1 - exp(-exp(-g / 2) / sqrt(pi))
  expected.p <- 1 - (1 - min(odd.p, even.p, infinity.p))^3
  result <- xu_lin_wei_pan_aspu_test(
    x, y, powers = powers, score_scale = "book_studentized"
  )

  expect_equal(unname(result$components$standard.error), standard.error,
               tolerance = 2e-14)
  expect_equal(unname(result$components$W), w, tolerance = 2e-14)
  expect_equal(unname(result$components$coordinate.correlation),
               unname(correlation), tolerance = 2e-14)
  expect_equal(unname(result$components$finite$statistic),
               reference$observed, tolerance = 2e-13)
  expect_equal(unname(result$components$finite$null.mean),
               reference$null.mean, tolerance = 2e-13)
  expect_equal(unname(result$components$finite$null.variance),
               reference$null.variance, tolerance = 3e-11)
  expect_equal(unname(result$components$finite$null.covariance),
               unname(reference$null.covariance), tolerance = 3e-11)
  expect_equal(unname(result$components$finite$null.correlation),
               unname(reference$null.correlation), tolerance = 3e-13)
  expect_equal(unname(result$components$finite$standardized.statistic),
               reference$standardized, tolerance = 3e-13)
  expect_equal(result$components$groups$odd$p.value, odd.p,
               tolerance = 2e-14)
  expect_equal(result$components$groups$even$p.value, even.p,
               tolerance = 2e-14)
  expect_equal(result$components$maximum$M, maximum, tolerance = 2e-14)
  expect_equal(result$components$maximum$G, g, tolerance = 2e-14)
  expect_equal(result$components$maximum$p.value, infinity.p,
               tolerance = 2e-14)
  expect_equal(result$p.value, expected.p, tolerance = 3e-14)
})


test_that("default aSPU matches the primary paper raw-coordinate formulas", {
  x <- cbind(
    c(-1.4, -0.2, 0.7, 1.5, 2.1, 0.3),
    c(0.2, 1.8, -0.9, 0.6, 2.4, -0.5),
    c(1.7, -0.6, 0.4, 1.2, -1.1, 2.0)
  )
  y <- cbind(
    c(-0.8, 0.5, 1.9, -1.2, 0.1, 2.6, -0.4, 0.9),
    c(1.1, -0.7, 0.3, 2.0, -1.3, 0.8, 1.5, -0.1),
    c(-1.0, 0.2, 1.4, -0.5, 2.2, 0.7, -0.8, 1.1)
  )
  finite <- 1:4
  n1 <- nrow(x)
  n2 <- nrow(y)
  pooled <- ((n1 - 1) * stats::cov(x) +
               (n2 - 1) * stats::cov(y)) / (n1 + n2 - 2)
  mean.covariance <- pooled * (1 / n1 + 1 / n2)
  difference <- colMeans(x) - colMeans(y)
  reference <- aspu_power_reference(
    difference, mean.covariance, finite
  )
  odd <- finite %% 2L == 1L
  odd.p <- aspu_group_reference(
    reference$standardized[odd],
    reference$null.correlation[odd, odd, drop = FALSE], "two.sided"
  )
  even.p <- aspu_group_reference(
    reference$standardized[!odd],
    reference$null.correlation[!odd, !odd, drop = FALSE], "upper"
  )
  result <- xu_lin_wei_pan_aspu_test(x, y, powers = finite)
  normalizer <- max(sqrt(diag(mean.covariance)))

  expect_identical(result$diagnostics$score.scale, "paper_raw")
  expect_equal(unname(result$components$finite$score), difference,
               tolerance = 2e-14)
  expect_equal(result$components$finite$statistic,
               stats::setNames(reference$observed, as.character(finite)),
               tolerance = 3e-13)
  expect_equal(result$components$finite$null.mean,
               stats::setNames(reference$null.mean, as.character(finite)),
               tolerance = 3e-13)
  expect_equal(result$components$finite$null.variance,
               stats::setNames(reference$null.variance, as.character(finite)),
               tolerance = 3e-12)
  expect_equal(unname(result$components$finite$null.covariance),
               unname(reference$null.covariance), tolerance = 3e-12)
  expect_equal(unname(result$components$finite$null.correlation),
               unname(reference$null.correlation), tolerance = 4e-13)
  expect_equal(unname(result$components$finite$standardized.statistic),
               reference$standardized, tolerance = 4e-13)
  expect_equal(
    unname(result$components$finite$null.coordinate.covariance),
    unname(mean.covariance), tolerance = 3e-14
  )
  expect_equal(result$components$finite$common.normalization.scale,
               normalizer, tolerance = 2e-14)
  expect_equal(unname(result$components$finite$normalized.score),
               difference / normalizer, tolerance = 3e-14)
  expect_equal(result$components$groups$odd$p.value, odd.p,
               tolerance = 3e-14)
  expect_equal(result$components$groups$even$p.value, even.p,
               tolerance = 3e-14)
  expect_equal(result$p.value, 1 - (1 - min(odd.p, even.p))^2,
               tolerance = 5e-13)
})


test_that("odd, even, and maximum families use their published tails", {
  set.seed(2811)
  x <- matrix(stats::rnorm(50), 10, 5)
  y <- matrix(stats::rnorm(65, 0.2), 13, 5)
  odd <- xu_lin_wei_pan_aspu_test(x, y, powers = c(1, 3, 5))
  mixed <- xu_lin_wei_pan_aspu_test(x, y, powers = c(1, 3, 2, 4))
  maximum <- xu_lin_wei_pan_aspu_test(x, y, powers = Inf)

  expect_equal(odd$components$combination.exponent, 1L)
  expect_equal(odd$p.value, odd$components$groups$odd$p.value,
               tolerance = 1e-15)
  expect_identical(odd$components$groups$odd$tail,
                   "joint two-sided normal")
  expect_equal(mixed$components$combination.exponent, 2L)
  expect_equal(
    mixed$p.value,
    1 - (1 - min(mixed$components$group.p.value))^2,
    tolerance = 5e-13
  )
  expect_identical(mixed$components$groups$even$tail,
                   "joint upper normal")
  expect_equal(maximum$components$combination.exponent, 1L)
  expect_equal(maximum$p.value,
               maximum$components$groups$infinite$p.value,
               tolerance = 1e-15)
  expect_identical(maximum$diagnostics$combination.groups, "infinite")
  expect_equal(odd$components$finite$p.value[["SPU_1"]],
               2 * stats::pnorm(
                 abs(odd$components$finite$standardized.statistic[["1"]]),
                 lower.tail = FALSE
               ), tolerance = 1e-15)
})


test_that("common, unequal, and supplied correlation paths are explicit", {
  set.seed(2812)
  variables <- c("a", "b", "c", "d")
  x <- matrix(stats::rnorm(40), 10, 4,
              dimnames = list(NULL, variables))
  y <- matrix(stats::rnorm(56, 0.1, 1.5), 14, 4,
              dimnames = list(NULL, variables))
  n1 <- nrow(x)
  n2 <- nrow(y)
  mean.covariance <- stats::cov(x) / n1 + stats::cov(y) / n2
  unequal <- xu_lin_wei_pan_aspu_test(
    x, y, powers = c(1, 2), correlation_source = "unequal"
  )
  expect_equal(
    unname(unequal$components$standard.error),
    unname(sqrt(diag(mean.covariance))), tolerance = 2e-14
  )
  expect_equal(
    unname(unequal$components$coordinate.correlation),
    unname(stats::cov2cor(mean.covariance)), tolerance = 2e-14
  )

  supplied.correlation <- stats::setNames(numeric(), character())
  supplied.correlation <- toeplitz(0.45^(0:3))
  dimnames(supplied.correlation) <- list(variables, variables)
  supplied.se <- stats::setNames(c(0.7, 1.1, 0.9, 1.4), variables)
  supplied <- xu_lin_wei_pan_aspu_test(
    x, y, powers = c(1, 2), correlation_source = "supplied",
    supplied_correlation = supplied.correlation,
    standard_errors = supplied.se
  )
  expect_equal(unname(supplied$components$standard.error),
               unname(supplied.se), tolerance = 2e-14)
  expect_equal(unname(supplied$components$W),
               unname((colMeans(x) - colMeans(y)) / supplied.se),
               tolerance = 2e-14)
  expect_equal(supplied$components$coordinate.correlation,
               supplied.correlation, tolerance = 2e-14)
  expect_identical(supplied$diagnostics$standard.error.source, "supplied")

  estimated.se <- xu_lin_wei_pan_aspu_test(
    x, y, powers = c(1, 2), correlation_source = "supplied",
    supplied_correlation = supplied.correlation
  )
  expect_equal(unname(estimated.se$components$standard.error),
               unname(sqrt(diag(mean.covariance))), tolerance = 2e-14)
  expect_identical(estimated.se$diagnostics$standard.error.source,
                   "unequal sample variances")
})


test_that("fixed hard banding follows the paper covariance rule", {
  set.seed(2813)
  x <- matrix(stats::rnorm(240), 60, 4)
  y <- matrix(stats::rnorm(280, 0.1), 70, 4)
  pooled <- ((nrow(x) - 1) * stats::cov(x) +
               (nrow(y) - 1) * stats::cov(y)) /
    (nrow(x) + nrow(y) - 2)
  banded <- pooled
  banded[abs(row(banded) - col(banded)) > 1] <- 0
  common <- xu_lin_wei_pan_aspu_test(
    x, y, powers = c(1, 2), bandwidth = 1
  )
  unequal <- xu_lin_wei_pan_aspu_test(
    x, y, powers = c(1, 2), correlation_source = "unequal",
    bandwidth = c(0, 0)
  )

  expect_equal(unname(common$components$coordinate.correlation),
               unname(stats::cov2cor(banded)), tolerance = 2e-14)
  expect_equal(unname(unequal$components$coordinate.correlation), diag(4),
               tolerance = 2e-14)
  expect_equal(common$diagnostics$bandwidth, c(first = 1L, second = 1L))
  expect_equal(unequal$diagnostics$bandwidth,
               c(first = 0L, second = 0L))
  expect_match(common$diagnostics$banding.rule, "hard zero")
})


test_that("aSPU is stable under exchange, translation, and positive units", {
  set.seed(2814)
  raw.x <- matrix(stats::rnorm(48), 12, 4)
  raw.y <- matrix(stats::rnorm(60, 0.15, 1.2), 15, 4)
  scale <- apply(rbind(raw.x, raw.y), 2L, function(z) 4 * max(abs(z)))
  x <- sweep(raw.x, 2L, scale, "/")
  y <- sweep(raw.y, 2L, scale, "/")
  baseline <- xu_lin_wei_pan_aspu_test(
    x, y, score_scale = "book_studentized",
    correlation_source = "unequal"
  )
  swapped <- xu_lin_wei_pan_aspu_test(
    y, x, score_scale = "book_studentized",
    correlation_source = "unequal"
  )
  shift <- c(13, -7, 2.5, 40)
  translated <- xu_lin_wei_pan_aspu_test(
    sweep(x, 2L, shift, "+"), sweep(y, 2L, shift, "+"),
    score_scale = "book_studentized",
    correlation_source = "unequal"
  )
  units <- c(1.7e308, 1e-308, 8e250, 3e-200)
  rescaled <- xu_lin_wei_pan_aspu_test(
    sweep(x, 2L, units, "*"), sweep(y, 2L, units, "*"),
    score_scale = "book_studentized",
    correlation_source = "unequal"
  )
  permuted <- xu_lin_wei_pan_aspu_test(
    x[c(8, 1, 12, 4, 2, 10, 6, 9, 3, 11, 5, 7), , drop = FALSE],
    y[c(15, 2, 9, 1, 12, 5, 14, 7, 3, 11, 4, 13, 8, 6, 10), ,
      drop = FALSE],
    score_scale = "book_studentized",
    correlation_source = "unequal"
  )

  expect_lte(abs(swapped$p.value - baseline$p.value), 3e-13)
  expect_equal(swapped$components$W, -baseline$components$W,
               tolerance = 3e-14)
  expect_lte(abs(swapped$components$groups$odd$p.value -
                   baseline$components$groups$odd$p.value), 3e-13)
  expect_lte(abs(swapped$components$groups$even$p.value -
                   baseline$components$groups$even$p.value), 3e-13)
  expect_lte(abs(translated$p.value - baseline$p.value), 5e-13)
  expect_equal(translated$components$W, baseline$components$W,
               tolerance = 5e-13)
  expect_lte(abs(rescaled$p.value - baseline$p.value), 5e-13)
  expect_equal(rescaled$components$W, baseline$components$W,
               tolerance = 5e-13)
  expect_equal(rescaled$components$coordinate.correlation,
               baseline$components$coordinate.correlation,
               tolerance = 5e-13)
  expect_lte(abs(permuted$p.value - baseline$p.value), 5e-13)
  expect_equal(permuted$components$finite$null.covariance,
               baseline$components$finite$null.covariance,
               tolerance = 5e-12)
})


test_that("paper-raw and book-studentized paths have distinct scale contracts", {
  set.seed(28141)
  x <- matrix(stats::rnorm(48), 12, 4)
  y <- matrix(stats::rnorm(60, 0.25, 1.3), 15, 4)
  units <- c(12, 0.15, 3.5, 0.4)
  scaled.x <- sweep(x, 2L, units, "*")
  scaled.y <- sweep(y, 2L, units, "*")
  paper <- xu_lin_wei_pan_aspu_test(x, y, powers = 1:4)
  paper.scaled <- xu_lin_wei_pan_aspu_test(
    scaled.x, scaled.y, powers = 1:4
  )
  book <- xu_lin_wei_pan_aspu_test(
    x, y, powers = 1:4, score_scale = "book_studentized"
  )
  book.scaled <- xu_lin_wei_pan_aspu_test(
    scaled.x, scaled.y, powers = 1:4,
    score_scale = "book_studentized"
  )
  paper.swapped <- xu_lin_wei_pan_aspu_test(y, x, powers = 1:4)
  difference.scaled <- (colMeans(x) - colMeans(y)) * units

  expect_equal(
    unname(paper.scaled$components$finite$statistic),
    vapply(1:4, function(power) sum(difference.scaled^power), numeric(1L)),
    tolerance = 5e-13
  )
  expect_gt(abs(paper.scaled$p.value - paper$p.value), 1e-5)
  expect_equal(book.scaled$p.value, book$p.value, tolerance = 5e-13)
  expect_equal(book.scaled$components$finite$standardized.statistic,
               book$components$finite$standardized.statistic,
               tolerance = 5e-13)
  expect_equal(paper.swapped$p.value, paper$p.value, tolerance = 5e-13)
  expect_equal(paper.swapped$components$groups$odd$p.value,
               paper$components$groups$odd$p.value, tolerance = 5e-13)

  huge <- xu_lin_wei_pan_aspu_test(1e150 * x, 1e150 * y, powers = c(1, 2))
  tiny <- xu_lin_wei_pan_aspu_test(1e-150 * x, 1e-150 * y,
                                  powers = c(1, 2))
  expect_equal(huge$p.value,
               xu_lin_wei_pan_aspu_test(x, y, powers = c(1, 2))$p.value,
               tolerance = 5e-13)
  expect_equal(tiny$p.value,
               xu_lin_wei_pan_aspu_test(x, y, powers = c(1, 2))$p.value,
               tolerance = 5e-13)
  expect_true(is.infinite(huge$components$finite$null.variance[["2"]]))
  expect_identical(tiny$components$finite$null.variance[["2"]], 0)
})


test_that("PSD adjustment is opt-in and fully diagnosed", {
  set.seed(2815)
  x <- matrix(stats::rnorm(30), 10, 3)
  y <- matrix(stats::rnorm(36), 12, 3)
  indefinite <- matrix(
    c(1, 0.9, 0.9, 0.9, 1, -0.9, 0.9, -0.9, 1), 3, 3
  )
  expect_lt(min(eigen(indefinite, symmetric = TRUE)$values), 0)
  expect_error(
    xu_lin_wei_pan_aspu_test(
      x, y, powers = c(1, 2), correlation_source = "supplied",
      supplied_correlation = indefinite
    ),
    "not positive semidefinite"
  )
  adjusted <- xu_lin_wei_pan_aspu_test(
    x, y, powers = c(1, 2), correlation_source = "supplied",
    supplied_correlation = indefinite, psd_adjust = "eigen_clip"
  )

  expect_true(adjusted$diagnostics$coordinate.psd$adjusted)
  expect_gt(adjusted$diagnostics$coordinate.psd$frobenius.adjustment, 0)
  expect_lt(adjusted$diagnostics$coordinate.psd$minimum.eigenvalue.before, 0)
  expect_gt(adjusted$diagnostics$coordinate.psd$minimum.eigenvalue.after, 0)
  expect_true(is.finite(adjusted$p.value))
  expect_true(all(diag(adjusted$components$coordinate.correlation) == 1))
  expect_identical(adjusted$diagnostics$variance.repair, "none")
})


test_that("finite powers support p one and minimum group sizes", {
  x <- matrix(c(-2, 3), ncol = 1)
  y <- matrix(c(-1, 4), ncol = 1)
  result <- xu_lin_wei_pan_aspu_test(
    x, y, powers = c(1, 2), score_scale = "book_studentized"
  )
  pooled <- (stats::var(x[, 1]) + stats::var(y[, 1])) / 2
  w <- (mean(x) - mean(y)) / sqrt(pooled)

  expect_s3_class(result, "hd_location_test")
  expect_s3_class(result, "htest")
  expect_equal(result$p, 1L)
  expect_equal(result$n, c(x = 2, y = 2))
  expect_equal(unname(result$components$W), w, tolerance = 2e-14)
  expect_equal(result$components$finite$statistic,
               c(`1` = w, `2` = w^2), tolerance = 2e-14)
  expect_equal(result$components$combination.exponent, 2L)
  expect_error(xu_lin_wei_pan_aspu_test(x, y, powers = Inf),
               "requires `p >= 2`")
})


test_that("aSPU preserves the htest and diagnostic contract", {
  set.seed(2816)
  variables <- c("height", "width", "depth")
  x <- matrix(stats::rnorm(30), 10, 3,
              dimnames = list(NULL, variables))
  y <- matrix(stats::rnorm(36), 12, 3,
              dimnames = list(NULL, variables))
  result <- xu_lin_wei_pan_aspu_test(x, y)

  expect_named(result$statistic, "minimum group p")
  expect_identical(names(result$estimate), variables)
  expect_identical(names(result$components$W), variables)
  expect_identical(dimnames(result$components$coordinate.correlation),
                   list(variables, variables))
  expect_identical(result$diagnostics$calibration,
                   "analytical; no permutation or bootstrap")
  expect_match(result$diagnostics$finite.power.coordinate.definition,
               "unstandardised.*primary paper")
  expect_identical(result$diagnostics$score.scale, "paper_raw")
  expect_identical(result$diagnostics$integration$odd$algorithm,
                   "mvtnorm::Miwa")
  expect_identical(result$diagnostics$integration$even$algorithm,
                   "mvtnorm::Miwa")
  expect_true(result$diagnostics$group.asymptotic.independence)
  expect_equal(result$components$combination.exponent, 3L)
  expect_identical(result$null.distribution$exact, FALSE)
  expect_equal(result$parameter, c(p = 3, groups = 3))
  expect_true(result$p.value >= 0 && result$p.value <= 1)
})


test_that("aSPU rejects malformed and degenerate inputs without repair", {
  set.seed(2817)
  x <- matrix(stats::rnorm(30), 10, 3)
  y <- matrix(stats::rnorm(36), 12, 3)

  expect_error(xu_lin_wei_pan_aspu_test(x[1, , drop = FALSE], y),
               "at least 2 row")
  expect_error(xu_lin_wei_pan_aspu_test(x, y[, 1:2]),
               "same number of columns")
  named.x <- x
  named.y <- y
  colnames(named.x) <- c("a", "b", "c")
  colnames(named.y) <- c("b", "a", "c")
  expect_error(xu_lin_wei_pan_aspu_test(named.x, named.y), "same names")
  x.bad <- x
  x.bad[1, 1] <- Inf
  expect_error(xu_lin_wei_pan_aspu_test(x.bad, y), "finite values")

  expect_error(xu_lin_wei_pan_aspu_test(x, y, powers = numeric()),
               "non-empty")
  expect_error(xu_lin_wei_pan_aspu_test(x, y, powers = c(1, 1)),
               "duplicates")
  expect_error(xu_lin_wei_pan_aspu_test(x, y, powers = 1.5),
               "positive integers")
  expect_error(xu_lin_wei_pan_aspu_test(x, y, powers = 0),
               "positive integers")
  expect_error(xu_lin_wei_pan_aspu_test(x, y, powers = -Inf),
               "only allowed")
  expect_error(xu_lin_wei_pan_aspu_test(x, y, powers = NA_real_),
               "without missing")
  expect_error(xu_lin_wei_pan_aspu_test(x, y, score_scale = "zscore"),
               "paper_raw")

  expect_error(xu_lin_wei_pan_aspu_test(
    x, y, correlation_source = "supplied"
  ), "supplied_correlation.*required")
  expect_error(xu_lin_wei_pan_aspu_test(
    x, y, supplied_correlation = diag(3)
  ), "only available")
  expect_error(xu_lin_wei_pan_aspu_test(x, y, bandwidth = -1),
               "0, \\.\\.\\., p - 1")
  expect_error(xu_lin_wei_pan_aspu_test(x, y, bandwidth = c(0, 1)),
               "one value")
  expect_error(xu_lin_wei_pan_aspu_test(
    x, y, correlation_source = "supplied",
    supplied_correlation = diag(3), bandwidth = 0
  ), "must be `NULL`")

  expect_error(xu_lin_wei_pan_aspu_test(
    x, y, correlation_source = "supplied",
    supplied_correlation = diag(2)
  ), "3 by 3")
  nonsymmetric <- diag(3)
  nonsymmetric[1, 2] <- 0.4
  expect_error(xu_lin_wei_pan_aspu_test(
    x, y, correlation_source = "supplied",
    supplied_correlation = nonsymmetric
  ), "symmetric")
  bad.diagonal <- diag(c(1, 0.9, 1))
  expect_error(xu_lin_wei_pan_aspu_test(
    x, y, correlation_source = "supplied",
    supplied_correlation = bad.diagonal
  ), "unit diagonal")
  expect_error(xu_lin_wei_pan_aspu_test(
    x, y, correlation_source = "supplied",
    supplied_correlation = diag(3), standard_errors = c(1, 0, 1)
  ), "strictly positive")
  named.correlation <- diag(3)
  dimnames(named.correlation) <- list(c("a", "b", "c"),
                                      c("a", "b", "c"))
  expect_error(xu_lin_wei_pan_aspu_test(
    x, y, correlation_source = "supplied",
    supplied_correlation = named.correlation
  ), "Dimnames.*must match")

  expect_error(xu_lin_wei_pan_aspu_test(x, y, psd_tol = 0),
               "strictly between")
  expect_error(xu_lin_wei_pan_aspu_test(x, y, miwa_steps = 1.5),
               "must be one integer")
  expect_error(xu_lin_wei_pan_aspu_test(
    x, y, powers = seq(1, 41, by = 2)
  ), "at most 20")

  x.constant <- cbind(seq_len(10), 1, stats::rnorm(10))
  y.constant <- cbind(seq_len(12) + 0.2, 1, stats::rnorm(12))
  expect_error(xu_lin_wei_pan_aspu_test(x.constant, y.constant),
               "marginal variance.*strictly positive|standard error.*strictly positive")
})
