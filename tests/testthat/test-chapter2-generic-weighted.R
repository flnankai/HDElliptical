ch2gw_asymmetric_fixture <- function() {
  x <- outer(
    seq_len(15), seq_len(4),
    function(i, j) sin(0.37 * i + 0.71 * j) +
      0.3 * cos(0.11 * i * j) + 0.02 * i * j
  )
  colnames(x) <- paste0("v", seq_len(ncol(x)))
  rownames(x) <- paste0("row", seq_len(nrow(x)))
  x
}


ch2gw_cube_fixture <- function() {
  as.matrix(expand.grid(
    first = c(-1, 1), second = c(-1, 1), third = c(-1, 1)
  ))
}


ch2gw_weights_reference <- function(radii, kind, power = NULL,
                                    callback = NULL) {
  switch(
    kind,
    constant = rep(1, length(radii)),
    inverse_norm = 1 / radii,
    power = radii^power,
    callback = vapply(radii, callback, numeric(1))
  )
}


ch2gw_literal_fit <- function(x, kind, power = NULL, callback = NULL,
                              tol = 1e-8, max_iter = 1000L) {
  location <- colMeans(x)
  diagonal <- apply(x, 2L, var)
  initial.diagonal <- diagonal
  location.history <- list(location)
  diagonal.history <- list(diagonal)
  for (iteration in seq_len(max_iter)) {
    residual <- sweep(x, 2L, location, "-")
    standardized <- sweep(residual, 2L, sqrt(diagonal), "/")
    radii <- sqrt(rowSums(standardized^2))
    directions <- standardized / radii
    weights <- ch2gw_weights_reference(
      radii, kind, power = power, callback = callback
    )
    denominator <- sum(weights / radii)
    numerator <- colSums(directions * weights)
    next.location <- location + sqrt(diagonal) * numerator / denominator
    next.diagonal <- ncol(x) * diagonal * colMeans(directions^2)
    location.update <- max(
      abs(next.location - location) / sqrt(initial.diagonal)
    )
    diagonal.update <- max(abs(log(next.diagonal / diagonal)))
    relative.update <- max(location.update, diagonal.update)
    location <- next.location
    diagonal <- next.diagonal
    location.history[[length(location.history) + 1L]] <- location
    diagonal.history[[length(diagonal.history) + 1L]] <- diagonal
    if (relative.update <= tol) break
  }
  residual <- sweep(x, 2L, location, "-")
  standardized <- sweep(residual, 2L, sqrt(diagonal), "/")
  radii <- sqrt(rowSums(standardized^2))
  directions <- standardized / radii
  weights <- ch2gw_weights_reference(
    radii, kind, power = power, callback = callback
  )
  list(
    location = location,
    diagonal = diagonal,
    directions = directions,
    radii = radii,
    weights = weights,
    numerator = colSums(directions * weights),
    denominator = sum(weights / radii),
    diagonal.equation = ncol(x) * colMeans(directions^2),
    iterations = iteration,
    relative.update = relative.update,
    location.history = do.call(rbind, location.history),
    diagonal.history = do.call(rbind, diagonal.history)
  )
}


ch2gw_oracle_reference <- function(x, theta, diagonal, weights.fun) {
  residual <- sweep(x, 2L, theta, "-")
  standardized <- sweep(residual, 2L, sqrt(diagonal), "/")
  radii <- sqrt(rowSums(standardized^2))
  directions <- standardized / radii
  weights <- vapply(radii, weights.fun, numeric(1))
  weighted <- directions * weights
  pair.total <- 0
  for (first in seq_len(nrow(x) - 1L)) {
    for (second in seq.int(first + 1L, nrow(x))) {
      pair.total <- pair.total + sum(weighted[first, ] * weighted[second, ])
    }
  }
  list(
    score = 2 * pair.total / (nrow(x) * (nrow(x) - 1)),
    directions = directions,
    radii = radii,
    weights = weights,
    empirical.nu2 = mean(weights^2),
    weighted.mean = colMeans(weighted)
  )
}


test_that("generic weighted HR recursion matches every literal update", {
  x <- ch2gw_asymmetric_fixture()
  power <- 0.5
  reference <- ch2gw_literal_fit(
    x, kind = "power", power = power, tol = 1e-9, max_iter = 2000L
  )
  observed <- generic_weighted_hr_location(
    x, K = "power", power = power, tol = 1e-9,
    max_iter = 2000L, keep_history = TRUE
  )

  expect_s3_class(observed, "generic_weighted_hr_location")
  expect_equal(unname(observed$location), unname(reference$location),
               tolerance = 2e-12)
  expect_equal(unname(observed$scale.diagonal), unname(reference$diagonal),
               tolerance = 2e-12)
  expect_equal(unname(observed$directions), unname(reference$directions),
               tolerance = 2e-12)
  expect_equal(unname(observed$radii), unname(reference$radii),
               tolerance = 2e-12)
  expect_equal(unname(observed$weights), unname(reference$weights),
               tolerance = 2e-12)
  expect_lte(
    max(abs(
      unname(observed$weighted.location.numerator) -
        unname(reference$numerator)
    )),
    2e-12
  )
  expect_equal(observed$weighted.location.denominator,
               reference$denominator, tolerance = 2e-12)
  expect_equal(unname(observed$diagonal.equation),
               unname(reference$diagonal.equation), tolerance = 2e-12)
  expect_equal(observed$diagnostics$iterations,
               reference$iterations, tolerance = 0)
  expect_lte(
    abs(
      observed$diagnostics$relative.update -
        reference$relative.update
    ),
    2e-12
  )
  expect_equal(unname(observed$location.history),
               unname(reference$location.history), tolerance = 2e-12)
  expect_equal(unname(observed$scale.diagonal.history),
               unname(reference$diagonal.history), tolerance = 2e-12)
  expect_identical(observed$diagnostics$diagonal.update,
                   "unweighted HR equation")
  expect_identical(observed$diagnostics$numerical.repair, "none")
})


test_that("constant, inverse-norm, and power special cases are locked", {
  x <- ch2gw_asymmetric_fixture()
  constant <- generic_weighted_hr_location(
    x, K = "constant", tol = 1e-8, max_iter = 2000L
  )
  power.zero <- generic_weighted_hr_location(
    x, K = "power", power = 0, tol = 1e-8, max_iter = 2000L
  )
  callback.constant <- generic_weighted_hr_location(
    x, K = function(radius) 1, tol = 1e-8, max_iter = 2000L
  )
  expect_equal(constant$location, power.zero$location, tolerance = 0)
  expect_equal(constant$scale.diagonal, power.zero$scale.diagonal,
               tolerance = 0)
  expect_equal(constant$location, callback.constant$location,
               tolerance = 0)
  expect_equal(constant$scale.diagonal, callback.constant$scale.diagonal,
               tolerance = 0)

  cube <- ch2gw_cube_fixture()
  inverse <- generic_weighted_hr_location(
    cube, K = "inverse_norm", tol = 1e-12, max_iter = 20L
  )
  power.minus.one <- generic_weighted_hr_location(
    cube, K = "power", power = -1, tol = 1e-12, max_iter = 20L
  )
  callback.inverse <- generic_weighted_hr_location(
    cube, K = function(radius) 1 / radius,
    tol = 1e-12, max_iter = 20L
  )
  expect_equal(unname(inverse$location), numeric(3), tolerance = 2e-16)
  expect_equal(inverse$location, power.minus.one$location, tolerance = 0)
  expect_equal(inverse$scale.diagonal,
               power.minus.one$scale.diagonal, tolerance = 0)
  expect_equal(inverse$location, callback.inverse$location, tolerance = 0)
  expect_equal(inverse$scale.diagonal,
               callback.inverse$scale.diagonal, tolerance = 0)
})


test_that("custom K receives scalars and preserves applicable equivariances", {
  x <- ch2gw_asymmetric_fixture()
  callback.lengths <- integer()
  weight <- function(radius) {
    callback.lengths <<- c(callback.lengths, length(radius))
    1 + 0.1 * radius
  }
  base <- generic_weighted_hr_location(
    x, K = weight, tol = 1e-8, max_iter = 2000L
  )
  expect_true(length(callback.lengths) > nrow(x))
  expect_true(all(callback.lengths == 1L))
  expect_identical(base$diagnostics$callback.evaluation,
                   "one scalar radius per R call")

  shift <- c(1.2, -0.8, 2.1, 0.4)
  translated <- generic_weighted_hr_location(
    sweep(x, 2L, shift, "+"), K = function(radius) 1 + 0.1 * radius,
    tol = 1e-8, max_iter = 2000L
  )
  expect_equal(unname(translated$location),
               unname(base$location) + shift, tolerance = 3e-12)
  expect_equal(unname(translated$scale.diagonal),
               unname(base$scale.diagonal), tolerance = 3e-12)
  expect_equal(unname(translated$radii), unname(base$radii),
               tolerance = 3e-12)

  multiplier <- 3.7
  scaled <- generic_weighted_hr_location(
    multiplier * x, K = function(radius) 1 + 0.1 * radius,
    tol = 1e-8, max_iter = 2000L
  )
  expect_equal(unname(scaled$location),
               multiplier * unname(base$location), tolerance = 4e-12)
  expect_equal(unname(scaled$scale.diagonal),
               multiplier^2 * unname(base$scale.diagonal),
               tolerance = 4e-12)
  expect_equal(unname(scaled$radii), unname(base$radii),
               tolerance = 4e-12)

  ordering <- c(3L, 1L, 4L, 2L)
  signs <- c(-1, 1, -1, 1)
  transformed.x <- sweep(x[, ordering, drop = FALSE], 2L, signs, "*")
  transformed <- generic_weighted_hr_location(
    transformed.x, K = function(radius) 1 + 0.1 * radius,
    tol = 1e-8, max_iter = 2000L
  )
  expect_equal(unname(transformed$location),
               unname(base$location)[ordering] * signs,
               tolerance = 4e-12)
  expect_equal(unname(transformed$scale.diagonal),
               unname(base$scale.diagonal)[ordering],
               tolerance = 4e-12)
  expect_equal(unname(transformed$directions),
               sweep(unname(base$directions)[, ordering, drop = FALSE],
                     2L, signs, "*"), tolerance = 4e-12)
  expect_equal(unname(transformed$radii), unname(base$radii),
               tolerance = 4e-12)
})


test_that("weighted HR callback and no-repair failures are strict", {
  x <- ch2gw_asymmetric_fixture()
  expect_error(
    generic_weighted_hr_location(
      x, K = function(radius) c(1, 2), max_iter = 10L
    ),
    "one finite numeric scalar"
  )
  expect_error(
    generic_weighted_hr_location(
      x, K = function(radius) NA_real_, max_iter = 10L
    ),
    "one finite numeric scalar"
  )
  expect_error(
    generic_weighted_hr_location(
      x, K = function(radius) "bad", max_iter = 10L
    ),
    "one finite numeric scalar"
  )
  expect_error(
    generic_weighted_hr_location(
      x, K = function(radius) stop("callback sentinel"), max_iter = 10L
    ),
    "callback sentinel"
  )
  expect_error(
    generic_weighted_hr_location(
      x, K = function(radius) -1, max_iter = 10L
    ),
    "denominator.*not strictly positive"
  )
  expect_error(
    generic_weighted_hr_location(
      x, K = function(radius) 0, max_iter = 10L
    ),
    "denominator.*not strictly positive"
  )
  expect_error(
    generic_weighted_hr_location(
      x, K = "constant", tol = 1e-16, max_iter = 1L
    ),
    "did not converge"
  )
  expect_error(
    generic_weighted_hr_location(
      x, K = "constant", initial_diagonal = c(1, 1, 0, 1)
    ),
    "strictly positive"
  )
  expect_error(
    generic_weighted_hr_location(x, K = "constant", power = 0),
    "only used"
  )
  nonfinite <- x
  nonfinite[1, 1] <- Inf
  expect_error(generic_weighted_hr_location(nonfinite), "finite values")

  coincident <- rbind(
    c(-1, -1), c(0, 0), c(1, 1), c(2, -2), c(-2, 2)
  )
  expect_equal(colMeans(coincident), c(0, 0), tolerance = 0)
  expect_error(
    generic_weighted_hr_location(coincident, K = "constant"),
    "radius is zero"
  )
})


test_that("oracle weighted-sign score matches the literal pairwise formula", {
  x <- ch2gw_asymmetric_fixture()
  theta <- c(0.2, -0.15, 0.4, 0.1)
  diagonal <- c(1.2, 0.8, 1.7, 0.6)
  weight <- function(radius) 1 + 0.2 * radius
  reference <- ch2gw_oracle_reference(x, theta, diagonal, weight)
  observed <- oracle_weighted_sign_sum_test(
    x, theta, diagonal, K = weight, keep_scores = TRUE
  )

  expect_s3_class(observed, "oracle_weighted_sign_sum_test")
  expect_false(observed$calibrated)
  expect_false("p.value" %in% names(observed))
  expect_equal(unname(observed$statistic), reference$score,
               tolerance = 2e-14)
  expect_equal(unname(observed$score), reference$score,
               tolerance = 2e-14)
  expect_equal(unname(observed$directions), unname(reference$directions),
               tolerance = 2e-14)
  expect_equal(unname(observed$radii), unname(reference$radii),
               tolerance = 2e-14)
  expect_equal(unname(observed$weights), unname(reference$weights),
               tolerance = 2e-14)
  expect_equal(observed$empirical.nu2, reference$empirical.nu2,
               tolerance = 2e-14)
  expect_equal(unname(observed$weighted.score.mean),
               unname(reference$weighted.mean), tolerance = 2e-14)
  expect_identical(observed$diagnostics$calibration, "none")
  expect_false(observed$diagnostics$feasible.leaveout.or.plugin)
})


test_that("oracle calibration is available only from complete supplied inputs", {
  x <- ch2gw_asymmetric_fixture()
  theta <- c(0.2, -0.15, 0.4, 0.1)
  diagonal <- c(1.2, 0.8, 1.7, 0.6)
  raw <- oracle_weighted_sign_sum_test(
    x, theta, diagonal, K = "constant"
  )
  supplied.sd <- 0.37
  direct <- oracle_weighted_sign_sum_test(
    x, theta, diagonal, K = "constant", null_sd = supplied.sd
  )
  expect_true(direct$calibrated)
  expect_equal(unname(direct$statistic),
               unname(raw$score) / supplied.sd, tolerance = 2e-15)
  expect_equal(direct$p.value,
               pnorm(unname(raw$score) / supplied.sd,
                     lower.tail = FALSE), tolerance = 2e-15)
  expect_identical(direct$diagnostics$calibration,
                   "supplied null_sd")

  nu2 <- 1.8
  trace.R2 <- 5.4
  expected.sd <- sqrt(
    2 * nu2^2 * trace.R2 /
      (nrow(x) * (nrow(x) - 1) * ncol(x)^2)
  )
  moments <- oracle_weighted_sign_sum_test(
    x, theta, diagonal, K = "constant",
    nu2 = nu2, trace_R2 = trace.R2
  )
  expect_true(moments$calibrated)
  expect_equal(moments$null.sd, expected.sd, tolerance = 2e-15)
  expect_equal(unname(moments$statistic),
               unname(raw$score) / expected.sd, tolerance = 2e-14)
  expect_equal(moments$p.value,
               pnorm(unname(raw$score) / expected.sd,
                     lower.tail = FALSE), tolerance = 2e-15)
  expect_identical(moments$diagnostics$calibration,
                   "supplied nu2 and trace_R2")

  expect_error(
    oracle_weighted_sign_sum_test(
      x, theta, diagonal, nu2 = nu2
    ),
    "must be supplied together"
  )
  expect_error(
    oracle_weighted_sign_sum_test(
      x, theta, diagonal, trace_R2 = trace.R2
    ),
    "must be supplied together"
  )
  expect_error(
    oracle_weighted_sign_sum_test(
      x, theta, diagonal, null_sd = supplied.sd,
      nu2 = nu2, trace_R2 = trace.R2
    ),
    "not conflicting"
  )
  expect_error(
    oracle_weighted_sign_sum_test(
      x, theta, diagonal, null_sd = 0
    ),
    "strictly positive"
  )
})


test_that("oracle built-ins and reference-geometry invariances are exact", {
  cube <- ch2gw_cube_fixture()
  theta <- c(0, 0, 0)
  diagonal <- c(1.1, 0.9, 1.4)
  constant <- oracle_weighted_sign_sum_test(
    cube, theta, diagonal, K = "constant"
  )
  power.zero <- oracle_weighted_sign_sum_test(
    cube, theta, diagonal, K = "power", power = 0
  )
  inverse <- oracle_weighted_sign_sum_test(
    cube, theta, diagonal, K = "inverse_norm"
  )
  power.minus.one <- oracle_weighted_sign_sum_test(
    cube, theta, diagonal, K = "power", power = -1
  )
  expect_equal(constant$score, power.zero$score, tolerance = 0)
  expect_equal(inverse$score, power.minus.one$score, tolerance = 0)

  x <- ch2gw_asymmetric_fixture()
  theta <- c(0.2, -0.15, 0.4, 0.1)
  diagonal <- c(1.2, 0.8, 1.7, 0.6)
  callback <- function(radius) 1 + 0.2 * radius
  base <- oracle_weighted_sign_sum_test(
    x, theta, diagonal, K = callback, keep_scores = TRUE
  )
  shift <- c(1.1, -0.7, 0.3, 2)
  multiplier <- 2.6
  shifted.scaled <- oracle_weighted_sign_sum_test(
    multiplier * sweep(x, 2L, shift, "+"),
    multiplier * (theta + shift), multiplier^2 * diagonal,
    K = callback, keep_scores = TRUE
  )
  expect_equal(shifted.scaled$score, base$score, tolerance = 3e-14)
  expect_equal(unname(shifted.scaled$radii), unname(base$radii),
               tolerance = 3e-14)
  expect_equal(unname(shifted.scaled$directions),
               unname(base$directions), tolerance = 3e-14)

  ordering <- c(4L, 2L, 1L, 3L)
  signs <- c(-1, 1, -1, 1)
  transformed <- oracle_weighted_sign_sum_test(
    sweep(x[, ordering, drop = FALSE], 2L, signs, "*"),
    theta[ordering] * signs, diagonal[ordering],
    K = callback, keep_scores = TRUE
  )
  expect_equal(transformed$score, base$score, tolerance = 3e-14)
  expect_equal(unname(transformed$radii), unname(base$radii),
               tolerance = 3e-14)
  expect_equal(unname(transformed$directions),
               sweep(unname(base$directions)[, ordering, drop = FALSE],
                     2L, signs, "*"), tolerance = 3e-14)
})


test_that("oracle geometry and callback failures never receive repairs", {
  x <- ch2gw_asymmetric_fixture()
  theta <- c(0.2, -0.15, 0.4, 0.1)
  diagonal <- c(1.2, 0.8, 1.7, 0.6)
  expect_error(
    oracle_weighted_sign_sum_test(
      x, theta, diagonal, K = function(radius) Inf
    ),
    "one finite numeric scalar"
  )
  expect_error(
    oracle_weighted_sign_sum_test(
      x, theta, diagonal, K = function(radius) 0
    ),
    "second radial-weight moment is zero"
  )
  expect_error(
    oracle_weighted_sign_sum_test(
      x, x[1, ], diagonal, K = "constant"
    ),
    "radius is zero"
  )
  expect_error(
    oracle_weighted_sign_sum_test(
      x, theta, replace(diagonal, 2L, -1), K = "constant"
    ),
    "strictly positive"
  )
  expect_error(
    oracle_weighted_sign_sum_test(
      x, theta, diagonal, K = "inverse_norm", power = -1
    ),
    "only used"
  )
})








