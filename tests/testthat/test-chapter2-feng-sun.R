fs_direction_reference <- function(value) {
  radius <- sqrt(sum(value^2))
  if (radius == 0) {
    return(rep(0, length(value)))
  }
  value / radius
}


fs_prepare_reference <- function(x, mu) {
  residual <- sweep(x, 2L, mu, "-")
  scale <- apply(abs(residual), 2L, max)
  list(
    residual = sweep(residual, 2L, scale, "/"),
    scale = scale
  )
}


fs_score_reference <- function(data, keep, theta, diagonal) {
  epsilon <- sweep(data[keep, , drop = FALSE], 2L, theta, "-")
  epsilon <- sweep(epsilon, 2L, sqrt(diagonal), "/")
  radii <- sqrt(rowSums(epsilon^2))
  if (any(radii == 0)) {
    stop("coincident training observation")
  }
  signs <- epsilon / radii
  p <- ncol(data)
  m <- length(keep)
  mean.sign <- colMeans(signs)
  scale.equation <- p * colMeans(signs^2)
  list(
    sign.sum = colSums(signs),
    sign.square.sum = colSums(signs^2),
    inverse.radius.sum = sum(1 / radii),
    location.residual = max(abs(mean.sign)),
    scale.residual = max(abs(scale.equation - 1))
  )
}


fs_fit_reference <- function(data, leave, tol, max_iter) {
  keep <- setdiff(seq_len(nrow(data)), leave)
  fit.data <- data[keep, , drop = FALSE]
  theta <- colMeans(fit.data)
  diagonal <- apply(fit.data, 2L, stats::var)
  diagonal <- diagonal / max(diagonal)
  p <- ncol(data)
  m <- length(keep)

  for (iteration in 0:max_iter) {
    score <- fs_score_reference(data, keep, theta, diagonal)
    if (score$location.residual <= tol && score$scale.residual <= tol) {
      return(list(
        theta = theta,
        diagonal = diagonal,
        iterations = iteration,
        location.residual = score$location.residual,
        scale.residual = score$scale.residual
      ))
    }
    if (iteration == max_iter) {
      stop("reference did not converge")
    }
    theta.new <- theta + sqrt(diagonal) * score$sign.sum /
      score$inverse.radius.sum
    diagonal.new <- diagonal * p * score$sign.square.sum / m
    diagonal.new <- diagonal.new / max(diagonal.new)
    theta <- theta.new
    diagonal <- diagonal.new
  }
  stop("unreachable")
}


fs_reference <- function(x, mu = numeric(ncol(x)), tol = 1e-7,
                         max_iter = 500L) {
  prepared <- fs_prepare_reference(x, mu)
  data <- prepared$residual
  n <- nrow(data)
  p <- ncol(data)
  pairs <- utils::combn(n, 2L)
  pair.count <- ncol(pairs)
  theta.standardized <- matrix(NA_real_, pair.count, p)
  theta.raw <- matrix(NA_real_, pair.count, p)
  diagonal.standardized <- matrix(NA_real_, pair.count, p)
  diagonal.raw <- matrix(NA_real_, pair.count, p)
  test.inner <- variance.inner.squared <- numeric(pair.count)
  iterations <- location.residual <- scale.residual <- numeric(pair.count)

  for (index in seq_len(pair.count)) {
    i <- pairs[1L, index]
    j <- pairs[2L, index]
    fit <- fs_fit_reference(data, c(i, j), tol, max_iter)
    theta.standardized[index, ] <- fit$theta
    theta.raw[index, ] <- mu + prepared$scale * fit$theta
    diagonal.standardized[index, ] <- fit$diagonal
    raw.diagonal <- prepared$scale^2 * fit$diagonal
    diagonal.raw[index, ] <- raw.diagonal / max(raw.diagonal)
    iterations[index] <- fit$iterations
    location.residual[index] <- fit$location.residual
    scale.residual[index] <- fit$scale.residual

    test.i <- fs_direction_reference(data[i, ] / sqrt(fit$diagonal))
    test.j <- fs_direction_reference(data[j, ] / sqrt(fit$diagonal))
    test.inner[index] <- sum(test.i * test.j)
    variance.i <- fs_direction_reference(
      (data[i, ] - fit$theta) / sqrt(fit$diagonal)
    )
    variance.j <- fs_direction_reference(
      (data[j, ] - fit$theta) / sqrt(fit$diagonal)
    )
    variance.inner.squared[index] <- sum(variance.i * variance.j)^2
  }

  test.sum <- sum(test.inner)
  variance.square.sum <- sum(variance.inner.squared)
  T.SS <- test.sum / pair.count
  trace.R2.hat <- p^2 * variance.square.sum / pair.count
  sigma2.hat <- variance.square.sum / pair.count^2
  list(
    prepared = prepared,
    pairs = t(pairs),
    theta.standardized = theta.standardized,
    theta.raw = theta.raw,
    diagonal.standardized = diagonal.standardized,
    diagonal.raw = diagonal.raw,
    test.inner = test.inner,
    variance.inner.squared = variance.inner.squared,
    iterations = iterations,
    location.residual = location.residual,
    scale.residual = scale.residual,
    test.sum = test.sum,
    variance.square.sum = variance.square.sum,
    T.SS = T.SS,
    trace.R2.hat = trace.R2.hat,
    sigma2.hat = sigma2.hat,
    statistic = T.SS / sqrt(sigma2.hat)
  )
}


test_that("Feng-Sun matches a literal feasible pure-R reference", {
  x <- rbind(
    c(0.7, -1.1, 0.2, 1.3),
    c(-0.4, 0.8, 1.2, -0.6),
    c(1.1, 0.3, -0.9, 0.5),
    c(-1.2, -0.5, 0.7, 1.0),
    c(0.2, 1.4, -0.4, -1.1),
    c(1.5, -0.2, 0.5, 0.1),
    c(-0.8, 0.6, -1.3, 0.9)
  )
  mu <- c(0.1, -0.15, 0.05, 0.2)
  reference <- fs_reference(x, mu, tol = 1e-7)
  result <- feng_sun_one_sample_test(x, mu, tol = 1e-7)

  expect_s3_class(result, "hd_location_test")
  expect_s3_class(result, "htest")
  expect_equal(unname(result$statistic), reference$statistic,
               tolerance = 2e-12)
  expect_equal(unname(result$raw.statistic), reference$T.SS,
               tolerance = 2e-13)
  expect_equal(result$components$T.SS, reference$T.SS,
               tolerance = 2e-13)
  expect_equal(result$components$test.inner.product.sum,
               reference$test.sum, tolerance = 3e-13)
  expect_equal(result$components$variance.inner.product.squared.sum,
               reference$variance.square.sum, tolerance = 3e-13)
  expect_equal(result$components$trace.R2.hat, reference$trace.R2.hat,
               tolerance = 4e-12)
  expect_equal(result$components$sigma2.hat, reference$sigma2.hat,
               tolerance = 3e-13)
  expect_equal(unname(result$variance), reference$sigma2.hat,
               tolerance = 3e-13)
  expect_equal(result$components$sigma.hat, sqrt(reference$sigma2.hat),
               tolerance = 3e-13)
  expect_equal(unname(result$components$pair.test.inner.product),
               reference$test.inner, tolerance = 3e-13)
  expect_equal(
    unname(result$components$pair.variance.inner.product.squared),
    reference$variance.inner.squared,
    tolerance = 3e-13
  )
  expect_equal(unname(result$components$leaveout.location.standardized),
               reference$theta.standardized, tolerance = 8e-12)
  expect_equal(unname(result$components$leaveout.location),
               reference$theta.raw, tolerance = 8e-12)
  expect_equal(
    unname(result$components$leaveout.scale.diagonal.standardized),
    reference$diagonal.standardized,
    tolerance = 1e-11
  )
  expect_equal(unname(result$components$leaveout.scale.diagonal),
               reference$diagonal.raw, tolerance = 1e-11)
  expect_equal(unname(result$diagnostics$pair.iterations),
               reference$iterations)
  location_residual_error <- max(abs(
    unname(result$diagnostics$pair.location.equation.residual) -
      reference$location.residual
  ))
  scale_residual_error <- max(abs(
    unname(result$diagnostics$pair.scale.equation.residual) -
      reference$scale.residual
  ))
  expect_true(is.finite(location_residual_error))
  expect_true(is.finite(scale_residual_error))
  expect_lt(location_residual_error, 5e-12)
  expect_lt(scale_residual_error, 5e-12)
  expect_equal(result$p.value,
               stats::pnorm(reference$statistic, lower.tail = FALSE),
               tolerance = 1e-15)
  expect_identical(result$null.distribution$tail, "upper")
  expect_false(result$null.distribution$exact)
  expect_identical(result$diagnostics$regularization, "none")
  expect_identical(result$diagnostics$variance.repair, "none")
})


test_that("ordered-pair denominators and theta roles are locked", {
  set.seed(2611)
  x <- matrix(stats::runif(48, -1, 1), 8, 6)
  mu <- seq(-0.12, 0.13, length.out = 6)
  result <- feng_sun_one_sample_test(x, mu, tol = 1e-6)
  components <- result$components
  prepared <- fs_prepare_reference(x, mu)
  pair.count <- choose(nrow(x), 2)
  ordered.count <- nrow(x) * (nrow(x) - 1)

  expect_equal(components$pair.count, pair.count)
  expect_equal(components$ordered.pair.count, ordered.count)
  expect_equal(components$T.SS,
               sum(components$pair.test.inner.product) / pair.count,
               tolerance = 1e-15)
  expect_equal(
    components$trace.R2.hat,
    ncol(x)^2 / ordered.count *
      2 * sum(components$pair.variance.inner.product.squared),
    tolerance = 2e-15
  )
  expect_equal(
    components$sigma2.hat,
    2 * components$trace.R2.hat /
      (ordered.count * ncol(x)^2),
    tolerance = 2e-15
  )

  test.rebuilt <- variance.rebuilt <- wrong.centered.test <-
    numeric(pair.count)
  for (index in seq_len(pair.count)) {
    i <- components$pair.index$i[index]
    j <- components$pair.index$j[index]
    diagonal <- components$leaveout.scale.diagonal.standardized[index, ]
    theta <- components$leaveout.location.standardized[index, ]
    test.i <- fs_direction_reference(
      prepared$residual[i, ] / sqrt(diagonal)
    )
    test.j <- fs_direction_reference(
      prepared$residual[j, ] / sqrt(diagonal)
    )
    test.rebuilt[index] <- sum(test.i * test.j)

    variance.i <- fs_direction_reference(
      (prepared$residual[i, ] - theta) / sqrt(diagonal)
    )
    variance.j <- fs_direction_reference(
      (prepared$residual[j, ] - theta) / sqrt(diagonal)
    )
    variance.rebuilt[index] <- sum(variance.i * variance.j)^2
    wrong.centered.test[index] <- sum(variance.i * variance.j)
  }
  expect_equal(test.rebuilt,
               unname(components$pair.test.inner.product),
               tolerance = 3e-15)
  expect_equal(variance.rebuilt,
               unname(components$pair.variance.inner.product.squared),
               tolerance = 3e-15)
  expect_gt(max(abs(wrong.centered.test - test.rebuilt)), 1e-3)
  expect_match(result$diagnostics$test.location.center,
               "theta is not used")
  expect_match(result$diagnostics$variance.location.center,
               "leave-two-out theta")
})


test_that("Feng-Sun is row-permutation and coordinatewise-scale invariant", {
  set.seed(2612)
  x <- matrix(stats::runif(63, -0.9, 0.9), 9, 7)
  mu <- stats::runif(7, -0.15, 0.15)
  baseline <- feng_sun_one_sample_test(x, mu, tol = 1e-6)

  permuted <- feng_sun_one_sample_test(
    x[c(4, 1, 8, 2, 9, 3, 7, 5, 6), ], mu, tol = 1e-6
  )
  expect_equal(permuted$statistic, baseline$statistic, tolerance = 2e-8)
  expect_equal(permuted$raw.statistic, baseline$raw.statistic,
               tolerance = 2e-9)
  expect_equal(permuted$components$trace.R2.hat,
               baseline$components$trace.R2.hat, tolerance = 3e-8)

  coordinate.scale <- c(-0.08, 0.3, -1.7, 4, 11, -0.6, 2.2)
  shift <- seq(-0.8, 0.9, length.out = 7)
  transformed <- feng_sun_one_sample_test(
    sweep(sweep(x, 2L, coordinate.scale, "*"), 2L, shift, "+"),
    mu * coordinate.scale + shift,
    tol = 1e-6
  )
  expect_equal(transformed$statistic, baseline$statistic,
               tolerance = 5e-9)
  expect_equal(transformed$raw.statistic, baseline$raw.statistic,
               tolerance = 5e-10)
  expect_equal(transformed$components$trace.R2.hat,
               baseline$components$trace.R2.hat, tolerance = 8e-9)
  expect_equal(
    transformed$components$leaveout.scale.diagonal.standardized,
    baseline$components$leaveout.scale.diagonal.standardized,
    tolerance = 2e-9
  )
})


test_that("Feng-Sun is generally not invariant under a dense rotation", {
  set.seed(2613)
  x <- matrix(stats::runif(56, -1, 1), 8, 7)
  mu <- stats::runif(7, -0.1, 0.1)
  baseline <- feng_sun_one_sample_test(x, mu, tol = 1e-6)
  q <- qr.Q(qr(matrix(stats::rnorm(49), 7, 7)))
  rotated <- feng_sun_one_sample_test(
    x %*% q, drop(mu %*% q), tol = 1e-6
  )

  expect_gt(abs(unname(rotated$statistic - baseline$statistic)), 1e-3)
  expect_gt(abs(rotated$components$T.SS - baseline$components$T.SS),
            1e-4)
})


test_that("Feng-Sun handles extreme coordinate units and subtraction overflow", {
  set.seed(2614)
  x <- matrix(stats::runif(54, -0.7, 0.7), 9, 6)
  mu <- seq(-0.12, 0.13, length.out = 6)
  baseline <- feng_sun_one_sample_test(x, mu, tol = 1e-6)
  coordinate.scale <- c(1e-300, 1e300, 1e-200, 1e200, 1e-100, 1e100)
  extreme <- feng_sun_one_sample_test(
    sweep(x, 2L, coordinate.scale, "*"),
    mu * coordinate.scale,
    tol = 1e-6
  )

  expect_true(is.finite(extreme$statistic))
  expect_true(is.finite(extreme$p.value))
  expect_equal(extreme$statistic, baseline$statistic, tolerance = 3e-8)
  expect_equal(extreme$raw.statistic, baseline$raw.statistic,
               tolerance = 3e-9)
  expect_equal(extreme$components$trace.R2.hat,
               baseline$components$trace.R2.hat, tolerance = 5e-8)

  set.seed(2701)
  safe.mu <- c(-1.05, 0.2, -0.3, 0.4, 0.1, -0.2)
  safe.x <- matrix(stats::runif(54, -1.25, 1.25), 9, 6)
  safe.x[1, 1] <- 1.2
  safe <- feng_sun_one_sample_test(safe.x, safe.mu, tol = 1e-6)
  huge <- feng_sun_one_sample_test(
    safe.x * 1e308, safe.mu * 1e308, tol = 1e-6
  )
  expect_equal(huge$diagnostics$subtraction.overflow.fallback.columns, 1)
  expect_true(all(is.finite(huge$components$leaveout.location)))
  expect_equal(huge$statistic, safe$statistic, tolerance = 4e-8)
  expect_equal(huge$raw.statistic, safe$raw.statistic, tolerance = 4e-9)
  expect_equal(huge$components$trace.R2.hat,
               safe$components$trace.R2.hat, tolerance = 8e-8)
})


test_that("Feng-Sun reports exact null signs without perturbing them", {
  set.seed(2615)
  mu <- c(0.1, -0.2, 0.05, 0.15, -0.1)
  x <- matrix(stats::runif(45, -1, 1), 9, 5)
  x[3, ] <- mu
  result <- feng_sun_one_sample_test(x, mu, tol = 1e-6)
  reference <- fs_reference(x, mu, tol = 1e-6)

  expect_equal(unname(result$statistic), reference$statistic,
               tolerance = 3e-11)
  expect_equal(result$diagnostics$test.zero.sign.pair.uses, nrow(x) - 1)
  expect_gt(result$diagnostics$variance.zero.sign.pair.uses, -1)
  expect_match(result$diagnostics$sign.at.zero, "U\\(0\\) = 0")
  expect_identical(result$diagnostics$variance.repair, "none")
})


test_that("Feng-Sun rejects invalid, degenerate, or unconverged fits", {
  expect_error(
    feng_sun_one_sample_test(matrix(1:9, 3, 3)),
    "at least 4 row"
  )
  expect_error(
    feng_sun_one_sample_test(matrix(1:20, 5, 4), mu = 1:3),
    "length 4"
  )
  expect_error(
    feng_sun_one_sample_test(matrix(stats::rnorm(24), 6, 4), tol = 0),
    "positive"
  )
  expect_error(
    feng_sun_one_sample_test(matrix(stats::rnorm(24), 6, 4),
                             max_iter = 2.5),
    "positive integer"
  )

  constant <- matrix(stats::rnorm(30), 6, 5)
  constant[, 2] <- 4
  expect_error(feng_sun_one_sample_test(constant), "positive variation")

  leaveout.constant <- rbind(
    c(1, -1, 0.5),
    c(2, 0.2, -0.7),
    c(0, 0.8, 1.1),
    c(0, -0.4, 0.3)
  )
  expect_error(
    feng_sun_one_sample_test(leaveout.constant),
    "leave-two-out initial marginal variance"
  )

  coincident <- rbind(
    c(0.7, -0.5),
    c(-0.8, 0.9),
    c(-1, -2),
    c(0, 0),
    c(1, 2)
  )
  expect_error(
    feng_sun_one_sample_test(coincident),
    "training observation exactly at its current location"
  )

  set.seed(2616)
  expect_error(
    feng_sun_one_sample_test(
      matrix(stats::runif(42, -1, 1), 7, 6), max_iter = 1L
    ),
    "did not converge"
  )
})


test_that("Feng-Sun preserves names and exposes convergence diagnostics", {
  set.seed(2617)
  variables <- c("a", "b", "c", "d", "e")
  x <- matrix(stats::runif(40, -1, 1), 8, 5,
              dimnames = list(NULL, variables))
  mu <- stats::setNames(c(0.1, -0.1, 0.05, 0, 0.15), variables)
  result <- feng_sun_one_sample_test(as.data.frame(x), mu, tol = 1e-6)

  expect_identical(names(result$estimate), variables)
  expect_identical(names(result$null.value), variables)
  expect_identical(colnames(result$components$leaveout.location), variables)
  expect_identical(
    colnames(result$components$leaveout.scale.diagonal), variables
  )
  expect_identical(rownames(result$components$leaveout.location),
                   rownames(result$components$pair.index))
  expect_true(result$diagnostics$all.leaveout.fits.converged)
  expect_lte(result$diagnostics$maximum.location.equation.residual, 1e-6)
  expect_lte(result$diagnostics$maximum.scale.equation.residual, 1e-6)
  expect_equal(unname(apply(
    result$components$leaveout.scale.diagonal, 1L, max
  )),
               rep(1, choose(nrow(x), 2)), tolerance = 1e-15)
  expect_equal(unname(apply(
    result$components$leaveout.log.scale.diagonal, 1L, max
  )),
               rep(0, choose(nrow(x), 2)), tolerance = 1e-15)
  expect_equal(unname(result$statistic),
               result$components$T.SS / result$components$sigma.hat,
               tolerance = 1e-15)
})
