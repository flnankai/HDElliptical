inst_reference_direction <- function(value, location, diagonal,
                                     zero_tol = 0) {
  epsilon <- (value - location) / sqrt(diagonal)
  radius <- sqrt(sum(epsilon^2))
  if (!is.finite(radius) || radius <= zero_tol) {
    stop("singular reference radius")
  }
  list(direction = epsilon / radius, radius = radius)
}


inst_reference_prepare <- function(x, mu) {
  residual <- sweep(x, 2L, mu, "-")
  scale <- apply(abs(residual), 2L, max)
  if (any(!is.finite(scale)) || any(scale <= 0)) {
    stop("invalid reference residual scale")
  }
  list(
    residual = sweep(residual, 2L, scale, "/"),
    scale = scale
  )
}


inst_reference_fit <- function(data, leave, tol = 1e-8,
                               max_iter = 1000L, zero_tol = 0) {
  keep <- setdiff(seq_len(nrow(data)), leave)
  retained <- data[keep, , drop = FALSE]
  location <- colMeans(retained)
  diagonal <- apply(retained, 2L, stats::var)
  if (any(!is.finite(diagonal)) || any(diagonal <= 0)) {
    stop("non-positive reference marginal variance")
  }
  initial.diagonal <- diagonal
  stable <- FALSE
  relative.update <- Inf
  location.update <- Inf
  diagonal.update <- Inf

  score <- function(location, diagonal) {
    epsilon <- sweep(retained, 2L, location, "-")
    epsilon <- sweep(epsilon, 2L, sqrt(diagonal), "/")
    radii <- sqrt(rowSums(epsilon^2))
    if (any(radii <= zero_tol)) stop("singular training radius")
    directions <- epsilon / radii
    list(
      directions = directions,
      radii = radii,
      location.residual = sqrt(sum(colMeans(directions)^2)),
      diagonal.residual = max(abs(
        ncol(data) * colMeans(directions^2) - 1
      ))
    )
  }

  for (update in seq_len(max_iter)) {
    current <- score(location, diagonal)
    minimum.radius <- min(current$radii)
    standardized.step <- colSums(current$directions) * minimum.radius /
      sum(minimum.radius / current$radii)
    next.location <- location + sqrt(diagonal) * standardized.step
    next.diagonal <- ncol(data) * diagonal *
      colMeans(current$directions^2)
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
  final <- score(location, diagonal)
  list(
    location = location,
    diagonal = diagonal,
    iterations = update,
    iteration.stable = stable,
    relative.update = relative.update,
    location.update = location.update,
    diagonal.update = diagonal.update,
    location.score = final$location.residual,
    diagonal.score = final$diagonal.residual,
    score.residual = max(
      final$location.residual, final$diagonal.residual
    ),
    minimum.residual.distance = min(final$radii)
  )
}


inst_reference <- function(x, mu = numeric(ncol(x)), tol = 1e-8,
                           max_iter = 1000L, zero_tol = 0) {
  prepared <- inst_reference_prepare(x, mu)
  data <- prepared$residual
  n <- nrow(data)
  p <- ncol(data)
  pairs <- t(utils::combn(n, 2L))
  pair.count <- nrow(pairs)
  fits <- vector("list", pair.count)
  locations <- diagonals <- matrix(NA_real_, pair.count, p)
  radii <- inverse.radii <- matrix(NA_real_, pair.count, 2L)
  test.inner <- statistic.kernel <- numeric(pair.count)
  sign.mean <- matrix(NA_real_, pair.count, p)
  variance.factor <- matrix(NA_real_, pair.count, 2L)
  variance.kernel <- trace.inner.squared <- numeric(pair.count)
  wrong.theta.kernel <- numeric(pair.count)

  for (pair in seq_len(pair.count)) {
    i <- pairs[pair, 1L]
    j <- pairs[pair, 2L]
    fits[[pair]] <- inst_reference_fit(
      data, c(i, j), tol, max_iter, zero_tol
    )
    location <- fits[[pair]]$location
    diagonal <- fits[[pair]]$diagonal
    locations[pair, ] <- location
    diagonals[pair, ] <- diagonal

    endpoint.i <- inst_reference_direction(
      data[i, ], numeric(p), diagonal, zero_tol
    )
    endpoint.j <- inst_reference_direction(
      data[j, ], numeric(p), diagonal, zero_tol
    )
    radii[pair, ] <- c(endpoint.i$radius, endpoint.j$radius)
    inverse.radii[pair, ] <- 1 / radii[pair, ]
    test.inner[pair] <- sum(
      endpoint.i$direction * endpoint.j$direction
    )
    statistic.kernel[pair] <- test.inner[pair] /
      prod(radii[pair, ])

    keep <- setdiff(seq_len(n), c(i, j))
    leaveout.signs <- t(vapply(
      keep,
      function(k) inst_reference_direction(
        data[k, ], numeric(p), diagonal, zero_tol
      )$direction,
      numeric(p)
    ))
    sign.mean[pair, ] <- colMeans(leaveout.signs)
    variance.factor[pair, 1L] <- sum(
      (endpoint.i$direction - sign.mean[pair, ]) *
        endpoint.j$direction
    )
    variance.factor[pair, 2L] <- sum(
      (endpoint.j$direction - sign.mean[pair, ]) *
        endpoint.i$direction
    )
    variance.kernel[pair] <- prod(inverse.radii[pair, ]^2) *
      prod(variance.factor[pair, ])

    trace.i <- inst_reference_direction(
      data[i, ], location, diagonal
    )$direction
    trace.j <- inst_reference_direction(
      data[j, ], location, diagonal
    )$direction
    trace.inner.squared[pair] <- sum(trace.i * trace.j)^2
    wrong.theta.kernel[pair] <- sum(trace.i * trace.j) /
      prod(radii[pair, ])
  }

  ordered.variance.kernel.sum <- 2 * sum(variance.kernel)
  variance <- 2 * ordered.variance.kernel.sum / n^4
  T <- mean(statistic.kernel)
  nu2 <- mean(inverse.radii^2)
  trace.R2 <- p^2 * mean(trace.inner.squared)
  oracle.variance <- 2 * nu2^2 * trace.R2 /
    (n * (n - 1) * p^2)

  list(
    prepared = prepared,
    pairs = pairs,
    fits = fits,
    locations = locations,
    diagonals = diagonals,
    radii = radii,
    inverse.radii = inverse.radii,
    test.inner = test.inner,
    statistic.kernel = statistic.kernel,
    sign.mean = sign.mean,
    variance.factor = variance.factor,
    variance.kernel = variance.kernel,
    trace.inner.squared = trace.inner.squared,
    wrong.theta.kernel = wrong.theta.kernel,
    ordered.variance.kernel.sum = ordered.variance.kernel.sum,
    T = T,
    variance = variance,
    z = T / sqrt(variance),
    inverse.radius.mean = mean(inverse.radii),
    nu2 = nu2,
    trace.R2 = trace.R2,
    oracle.variance = oracle.variance,
    weighted.trace = variance * n * (n - 1) / 2
  )
}


test_that("INST matches a literal primary-paper feasible reference", {
  x <- rbind(
    c(0.8, -1.2, 0.3),
    c(-0.5, 0.9, 1.1),
    c(1.2, 0.4, -0.8),
    c(-1.1, -0.6, 0.7),
    c(0.3, 1.3, -0.5),
    c(1.4, -0.3, 0.6),
    c(-0.7, 0.5, -1.2),
    c(0.1, -0.8, 1.4)
  )
  mu <- c(0.12, -0.08, 0.17)
  tolerance <- 1e-7
  reference <- inst_reference(x, mu, tolerance, 2000L)
  result <- inst_one_sample_test(
    x, mu, tol = tolerance, max_iter = 2000L
  )

  expect_s3_class(result, "hd_location_test")
  expect_s3_class(result, "htest")
  expect_equal(unname(result$statistic), reference$z,
               tolerance = 2e-10)
  expect_equal(result$components$T.INST, reference$T,
               tolerance = 2e-11)
  expect_equal(result$components$sigma2.direct.hat,
               reference$variance, tolerance = 3e-11)
  expect_equal(result$components$sigma.direct.hat,
               sqrt(reference$variance), tolerance = 3e-11)
  expect_equal(result$components$ordered.variance.kernel.sum,
               reference$ordered.variance.kernel.sum,
               tolerance = 4e-10)
  expect_equal(result$components$inverse.radius.mean.hat,
               reference$inverse.radius.mean, tolerance = 3e-11)
  expect_equal(result$components$nu2.IN.hat, reference$nu2,
               tolerance = 3e-11)
  expect_identical(result$components$nu2.IN.hat,
                   result$components$c0.IN.hat)
  expect_equal(result$components$trace.R2.hat,
               reference$trace.R2, tolerance = 3e-10)
  expect_equal(
    result$components$sigma2.oracle.factorized.diagnostic,
    reference$oracle.variance,
    tolerance = 5e-11
  )
  expect_equal(
    result$components$trace.weighted.score.covariance.squared.hat,
    reference$weighted.trace,
    tolerance = 5e-11
  )
  expect_equal(unname(result$components$pair.test.inner.product),
               reference$test.inner, tolerance = 3e-11)
  expect_equal(unname(result$components$pair.statistic.kernel),
               reference$statistic.kernel, tolerance = 3e-11)
  expect_equal(unname(result$components$pair.variance.factor),
               reference$variance.factor, tolerance = 4e-11)
  expect_equal(unname(result$components$pair.variance.kernel),
               reference$variance.kernel, tolerance = 6e-10)
  expect_equal(
    unname(result$components$pair.trace.inner.product.squared),
    reference$trace.inner.squared,
    tolerance = 4e-10
  )
  expect_equal(unname(result$components$endpoint.radius),
               reference$radii, tolerance = 3e-10)
  expect_equal(unname(result$components$endpoint.inverse.radius),
               reference$inverse.radii, tolerance = 3e-10)
  expect_equal(
    unname(result$components$leaveout.location.standardized),
    reference$locations,
    tolerance = 2e-9
  )
  expect_equal(
    unname(result$components$leaveout.scale.diagonal.standardized),
    reference$diagonals,
    tolerance = 2e-9
  )
  expect_equal(result$p.value,
               stats::pnorm(reference$z, lower.tail = FALSE),
               tolerance = 1e-15)
  expect_identical(result$null.distribution$tail, "upper")
  expect_false(result$null.distribution$exact)
})


test_that("published n^-4 factor and location roles are locked", {
  set.seed(28102)
  x <- matrix(stats::runif(27, -1, 1), 9, 3)
  mu <- c(-0.13, 0.09, 0.18)
  result <- inst_one_sample_test(x, mu, tol = 1e-7, max_iter = 2000L)
  components <- result$components
  n <- nrow(x)
  p <- ncol(x)

  expect_equal(components$pair.count, choose(n, 2))
  expect_equal(components$ordered.pair.count, n * (n - 1))
  expect_equal(
    components$T.INST,
    mean(components$pair.statistic.kernel),
    tolerance = 2e-15
  )
  expect_equal(
    components$ordered.variance.kernel.sum,
    2 * sum(components$pair.variance.kernel),
    tolerance = 3e-15
  )
  expect_equal(
    components$sigma2.direct.hat,
    2 * components$ordered.variance.kernel.sum / n^4,
    tolerance = 2e-15
  )
  wrong.pair.denominator <-
    2 * components$ordered.variance.kernel.sum /
    (n^2 * (n - 1)^2)
  expect_gt(abs(components$sigma2.direct.hat - wrong.pair.denominator),
            1e-6)
  expect_equal(
    components$nu2.IN.hat,
    mean(components$endpoint.inverse.radius^2),
    tolerance = 2e-15
  )
  expect_equal(
    components$trace.R2.hat,
    p^2 * mean(components$pair.trace.inner.product.squared),
    tolerance = 3e-15
  )
  expect_equal(
    components$sigma2.oracle.factorized.diagnostic,
    2 * components$nu2.IN.hat^2 * components$trace.R2.hat /
      (n * (n - 1) * p^2),
    tolerance = 3e-15
  )
  expect_equal(
    components$trace.weighted.score.covariance.squared.hat,
    components$sigma2.direct.hat * n * (n - 1) / 2,
    tolerance = 2e-15
  )

  reference <- inst_reference(x, mu, tol = 1e-7, max_iter = 2000L)
  expect_gt(max(abs(reference$wrong.theta.kernel -
                      reference$statistic.kernel)), 1e-3)
  expect_match(result$diagnostics$endpoint.center,
               "null location mu only")
  expect_match(result$diagnostics$endpoint.center,
               "only in estimating D_ij")
  expect_match(result$diagnostics$variance.sign.center,
               "null-centred")
  expect_match(result$diagnostics$primary.variance.estimator,
               "2 \\* n\\^\\(-4\\)")
  expect_identical(result$diagnostics$primary.calibration,
                   "direct feasible variance only")
  expect_identical(result$diagnostics$oracle.factorized.variance.role,
                   "diagnostic only")
})


test_that("INST is permutation, diagonal-scale, and common-translation invariant", {
  set.seed(28103)
  x <- matrix(stats::runif(30, -0.9, 0.9), 10, 3)
  mu <- c(-0.11, 0.07, 0.16)
  baseline <- inst_one_sample_test(
    x, mu, tol = 1e-7, max_iter = 2000L
  )

  permuted <- inst_one_sample_test(
    x[c(4, 1, 9, 2, 10, 3, 8, 5, 7, 6), ], mu,
    tol = 1e-7, max_iter = 2000L
  )
  expect_equal(permuted$statistic, baseline$statistic,
               tolerance = 3e-8)
  expect_equal(permuted$raw.statistic, baseline$raw.statistic,
               tolerance = 3e-9)
  expect_equal(permuted$variance, baseline$variance,
               tolerance = 5e-9)
  expect_equal(permuted$components$nu2.IN.hat,
               baseline$components$nu2.IN.hat, tolerance = 4e-9)
  expect_equal(permuted$components$trace.R2.hat,
               baseline$components$trace.R2.hat, tolerance = 5e-8)

  coordinate.scale <- c(-0.07, 3.5, -120)
  # Keep the translated observations exactly informative at double precision;
  # the separate overflow test below exercises much larger operands.
  shift <- c(1e6, -2e6, 3e6)
  transformed <- inst_one_sample_test(
    sweep(sweep(x, 2L, coordinate.scale, "*"), 2L, shift, "+"),
    mu * coordinate.scale + shift,
    tol = 1e-7, max_iter = 2000L
  )
  expect_equal(transformed$statistic, baseline$statistic,
               tolerance = 2e-6)
  expect_equal(transformed$raw.statistic, baseline$raw.statistic,
               tolerance = 2e-7)
  expect_equal(transformed$variance, baseline$variance,
               tolerance = 3e-7)
  expect_equal(transformed$components$nu2.IN.hat,
               baseline$components$nu2.IN.hat, tolerance = 2e-7)
  expect_equal(transformed$components$trace.R2.hat,
               baseline$components$trace.R2.hat, tolerance = 2e-6)
})


test_that("INST exposes distinct iteration stability and score diagnostics", {
  set.seed(28105)
  x <- matrix(stats::runif(27, -0.8, 0.8), 9, 3)
  stable <- inst_one_sample_test(
    x, tol = 1e-7, max_iter = 2000L
  )
  expect_true(stable$diagnostics$iteration.stable)
  expect_equal(stable$diagnostics$stability.failures, 0L)
  expect_lte(stable$diagnostics$relative.update, 1e-7)
  expect_true(is.finite(stable$diagnostics$score.residual))
  expect_true(is.finite(stable$diagnostics$minimum.residual.distance))
  expect_identical(
    stable$diagnostics$convergence.basis,
    "relative location and log-diagonal iterate update"
  )
  expect_true(all(stable$diagnostics$pair.fit$iteration.stable))
  expect_equal(
    stable$diagnostics$relative.update,
    max(stable$diagnostics$pair.fit$relative.update),
    tolerance = 0
  )
  expect_equal(
    stable$diagnostics$score.residual,
    max(stable$diagnostics$pair.fit$score.residual),
    tolerance = 0
  )

  expect_error(
    inst_one_sample_test(x, max_iter = 1L, strict = TRUE),
    "did not become iteration-stable"
  )
  expect_warning(
    unstable <- inst_one_sample_test(x, max_iter = 1L, strict = FALSE),
    "did not become iteration-stable"
  )
  expect_false(unstable$diagnostics$iteration.stable)
  expect_gt(unstable$diagnostics$stability.failures, 0L)
  expect_true(any(!unstable$diagnostics$pair.fit$iteration.stable))
  expect_identical(unstable$diagnostics$regularization, "none")
  expect_identical(unstable$diagnostics$variance.repair, "none")
})


test_that("INST rejects singular and invalid inputs without repair", {
  expect_error(
    inst_one_sample_test(matrix(1:9, 3, 3)),
    "at least 4 row"
  )
  expect_error(
    inst_one_sample_test(matrix(1:12, 4, 3), mu = 1:2),
    "length 3"
  )
  bad <- matrix(stats::rnorm(24), 8, 3)
  bad[2, 1] <- NA_real_
  expect_error(inst_one_sample_test(bad), "finite")
  expect_error(inst_one_sample_test(matrix(stats::rnorm(24), 8, 3),
                                    alpha = 1),
               "strictly between")
  expect_error(inst_one_sample_test(matrix(stats::rnorm(24), 8, 3),
                                    strict = NA),
               "TRUE.*FALSE")
  expect_error(inst_one_sample_test(matrix(stats::rnorm(24), 8, 3),
                                    tol = 0),
               "finite positive")
  expect_error(inst_one_sample_test(matrix(stats::rnorm(24), 8, 3),
                                    max_iter = 0),
               "positive integer")
  expect_error(inst_one_sample_test(matrix(stats::rnorm(24), 8, 3),
                                    zero_tol = -1),
               "non-negative")

  constant.coordinate <- cbind(
    matrix(stats::rnorm(16), 8, 2),
    rep(2, 8)
  )
  expect_error(
    inst_one_sample_test(constant.coordinate, mu = c(0, 0, 0)),
    "positive variation.*No ridge"
  )

  mu <- c(0.1, -0.2, 0.3)
  at.null <- matrix(stats::runif(27, -1, 1), 9, 3)
  at.null[4, ] <- mu
  expect_error(
    inst_one_sample_test(at.null, mu),
    "inverse-norm weighting is singular.*no perturbation"
  )
  expect_error(
    inst_one_sample_test(matrix(stats::runif(27, -1, 1), 9, 3),
                         zero_tol = 100),
    "at or below `zero_tol`"
  )

  # In one dimension the leave-two-out spatial median can equal a retained
  # observation exactly.  The recursion must stop instead of perturbing it.
  expect_error(
    inst_one_sample_test(matrix(c(-2, 0, 2, 5, 9), ncol = 1L)),
    "retained observation.*no perturbation"
  )

  # This deterministic fixture has a negative direct S.3 estimate.  The
  # public function must expose the estimator's failure, not take abs() or
  # replace it by an oracle-factorized diagnostic.
  negative.variance <- matrix(
    c(0.626492443494499, 0.298304968979210,
      0.0474627600051463, 0.482336722780019),
    ncol = 1L
  )
  expect_error(
    inst_one_sample_test(
      negative.variance, tol = 1e-7, max_iter = 2000L
    ),
    "direct feasible variance.*strictly positive.*no absolute value"
  )
})


test_that("overflow-safe subtraction preserves the feasible calculation", {
  x <- rbind(
    c(0.8, -1.2, 0.3),
    c(-0.5, 0.9, 1.1),
    c(1.2, 0.4, -0.8),
    c(-1.1, -0.6, 0.7),
    c(0.3, 1.3, -0.5),
    c(1.4, -0.3, 0.6),
    c(-0.7, 0.5, -1.2),
    c(0.1, -0.8, 1.4)
  )
  mu <- c(-0.5, -0.08, 0.17)
  baseline <- inst_one_sample_test(
    x, mu, tol = 1e-7, max_iter = 2000L
  )
  huge <- inst_one_sample_test(
    x * 1e308, mu * 1e308,
    tol = 1e-7, max_iter = 2000L
  )

  expect_equal(huge$diagnostics$subtraction.overflow.fallback.columns, 1)
  expect_true(is.finite(huge$statistic))
  expect_true(is.finite(huge$p.value))
  expect_true(all(is.finite(huge$components$leaveout.location)))
  expect_equal(huge$statistic, baseline$statistic, tolerance = 3e-7)
  expect_equal(huge$raw.statistic, baseline$raw.statistic,
               tolerance = 3e-8)
  expect_equal(huge$variance, baseline$variance, tolerance = 5e-8)
  expect_equal(huge$components$nu2.IN.hat,
               baseline$components$nu2.IN.hat, tolerance = 4e-8)
  expect_equal(huge$components$trace.R2.hat,
               baseline$components$trace.R2.hat, tolerance = 4e-7)
})


test_that("INST result is self-consistent and preserves variable names", {
  set.seed(28106)
  variables <- c("gene_a", "gene_b", "gene_c")
  x <- matrix(stats::runif(30, -1, 1), 10, 3,
              dimnames = list(NULL, variables))
  mu <- stats::setNames(c(0.05, -0.12, 0.18), variables)
  result <- inst_one_sample_test(
    as.data.frame(x), mu, alpha = 0.1,
    tol = 1e-7, max_iter = 2000L
  )

  expect_identical(names(result$estimate), variables)
  expect_identical(names(result$null.value), variables)
  expect_identical(
    colnames(result$components$leaveout.location), variables
  )
  expect_identical(
    colnames(result$components$leaveout.scale.diagonal), variables
  )
  expect_equal(
    unname(result$statistic),
    result$components$T.INST / result$components$sigma.direct.hat,
    tolerance = 2e-15
  )
  expect_equal(unname(result$variance),
               result$components$sigma2.direct.hat, tolerance = 0)
  expect_equal(result$p.value,
               stats::pnorm(unname(result$statistic), lower.tail = FALSE),
               tolerance = 0)
  expect_identical(
    result$diagnostics$rejection$reject,
    isTRUE(unname(result$statistic) > stats::qnorm(0.9))
  )
  expect_identical(result$diagnostics$radius.perturbation, "none")
  expect_identical(result$diagnostics$weight.cap, "none")
  expect_identical(result$diagnostics$variance.repair, "none")
})
