test_that("spatial signs implement the zero convention and scale invariance", {
  expect_equal(as.numeric(spatial_sign(c(3, 4))), c(0.6, 0.8),
               tolerance = 1e-14)
  signs <- spatial_sign(rbind(c(0, 0), c(3, 4), c(-6, -8)))
  expect_equal(unname(signs), rbind(c(0, 0), c(0.6, 0.8), c(-0.6, -0.8)),
               ignore_attr = TRUE)
  expect_equal(attr(signs, "n_zero"), 1)
  expect_equal(
    spatial_sign(matrix(c(1, 2), 1), center = c(1, 2)),
    matrix(c(0, 0), 1),
    ignore_attr = TRUE
  )
  huge <- spatial_sign(c(1e200, -1e200))
  expect_true(all(is.finite(huge)))
  expect_equal(as.numeric(huge), c(1, -1) / sqrt(2), tolerance = 1e-14)
  expect_identical(attr(huge, "norms"), attr(huge, "norm"))
  near_limit <- spatial_sign(c(1.7e308, 1.7e308))
  expect_equal(as.numeric(near_limit), rep(1 / sqrt(2), 2),
               tolerance = 1e-14)
  expect_true(is.infinite(attr(near_limit, "norm")))
  subnormal <- spatial_sign(rep(5e-324, 2))
  expect_equal(as.numeric(subnormal), rep(1 / sqrt(2), 2),
               tolerance = 1e-14)
})


test_that("modified Weiszfeld iteration returns a spatial median", {
  x <- rbind(c(-1, 0), c(1, 0), c(0, -1), c(0, 1))
  estimate <- spatial_median(x)
  expect_equal(as.numeric(estimate), c(0, 0), tolerance = 1e-10)
  expect_true(isTRUE(attr(estimate, "converged")))
  expect_lte(attr(estimate, "equation_residual"), 1e-8)

  shift <- c(10, -7)
  shifted <- spatial_median(sweep(x, 2, shift, "+"))
  expect_equal(as.numeric(shifted), as.numeric(estimate + shift),
               tolerance = 1e-8)

  set.seed(1234)
  heavy <- matrix(stats::rt(120, df = 2), 40, 3)
  huge_shift <- c(1e9, -2e9, 3e9)
  baseline <- spatial_median(heavy, tol = 1e-9)
  translated <- spatial_median(sweep(heavy, 2, huge_shift, "+"), tol = 1e-9)
  expect_equal(as.numeric(translated - huge_shift), as.numeric(baseline),
               tolerance = 2e-6)
  expect_lte(attr(translated, "equation_residual"), 1e-7)

  one_dimensional <- spatial_median(matrix(c(-2, 0, 1, 100), ncol = 1))
  expect_gte(one_dimensional, 0)
  expect_lte(one_dimensional, 1)

  near_observation <- rbind(c(0, 0), c(1e-10, 1), c(1, 1e-10))
  near_fit <- spatial_median(near_observation, tol = 1e-10, max_iter = 5000)
  expect_true(isTRUE(attr(near_fit, "converged")))
  expect_lte(attr(near_fit, "equation_residual"), 1e-10)
  expect_equal(
    as.numeric(near_fit), rep((3 - sqrt(3)) / 6, 2), tolerance = 2e-8
  )

  subnormal_start <- rbind(
    c(0, 0), c(1e-320, 0), c(1, 0), c(0, 1)
  )
  subnormal_fit <- spatial_median(
    subnormal_start, initial = c(0, 0), tol = 1e-9,
    max_iter = 10000, warn = FALSE
  )
  expect_true(isTRUE(attr(subnormal_fit, "converged")))
  expect_lte(attr(subnormal_fit, "equation_residual"), 1e-9)
  expect_equal(as.numeric(subnormal_fit), c(0, 0), tolerance = 1e-300)
})


test_that("spatial ranks agree with a direct R implementation", {
  x <- rbind(c(0, 0), c(1, 0), c(0, 2))
  direct_sign <- function(z) {
    length <- sqrt(sum(z^2))
    if (length == 0) z else z / length
  }
  expected <- t(vapply(seq_len(nrow(x)), function(i) {
    rowMeans(vapply(seq_len(nrow(x)), function(j) {
      direct_sign(x[i, ] - x[j, ])
    }, numeric(ncol(x))))
  }, numeric(ncol(x))))
  expect_equal(spatial_rank(x), expected, tolerance = 1e-14,
               ignore_attr = TRUE)

  extreme <- matrix(c(-1e308, 1e308), ncol = 1)
  expect_equal(as.numeric(spatial_rank(extreme)), c(-0.5, 0.5),
               tolerance = 1e-14)

  smallest <- 5e-324
  mixed_scale <- rbind(c(0, 0), c(smallest, smallest), c(1, 0))
  stable_sign <- function(left, right) {
    difference <- left - right
    scale <- max(abs(difference))
    if (scale == 0) {
      return(numeric(length(left)))
    }
    scaled <- difference / scale
    scaled / sqrt(sum(scaled^2))
  }
  expected_mixed <- t(vapply(seq_len(nrow(mixed_scale)), function(i) {
    rowMeans(vapply(seq_len(nrow(mixed_scale)), function(j) {
      stable_sign(mixed_scale[i, ], mixed_scale[j, ])
    }, numeric(ncol(mixed_scale))))
  }, numeric(ncol(mixed_scale))))
  expect_equal(spatial_rank(mixed_scale), expected_mixed,
               tolerance = 1e-14, ignore_attr = TRUE)

  epsilon <- .Machine$double.eps
  large_left <- c(1.23456789012345e308, 1.65432109876543e308)
  large_right <- large_left * (1 - epsilon * c(4, 2))
  finite_difference <- large_left - large_right
  expected_direction <- finite_difference / max(abs(finite_difference))
  expected_direction <- expected_direction /
    sqrt(sum(expected_direction^2))
  large_rank <- spatial_rank(rbind(large_left, large_right))
  expect_equal(as.numeric(large_rank[1, ]), expected_direction / 2,
               tolerance = 1e-14)
  expect_equal(as.numeric(large_rank[2, ]), -expected_direction / 2,
               tolerance = 1e-14)
})


test_that("SSCM has the stated trace and orthogonal equivariance", {
  x <- rbind(c(1, 2), c(-2, 1), c(3, -1), c(-1, -3))
  estimate <- sscm(x, center = "none")
  expect_equal(sum(diag(estimate)), 1, tolerance = 1e-14)
  expect_equal(estimate, t(estimate), tolerance = 1e-14,
               ignore_attr = TRUE)

  angle <- 0.37
  rotation <- matrix(c(cos(angle), -sin(angle), sin(angle), cos(angle)), 2)
  rotated <- sscm(x %*% rotation, center = "none")
  expected <- t(rotation) %*% unclass(estimate) %*% rotation
  expect_equal(unclass(rotated), expected, tolerance = 1e-12,
               ignore_attr = TRUE)
})


test_that("spatial Kendall matrix is the exact pairwise U-statistic", {
  x <- rbind(c(0, 0), c(1, 0), c(0, 2), c(-1, 1))
  estimate <- spatial_kendall(x)
  pairs <- combn(nrow(x), 2)
  expected <- matrix(0, ncol(x), ncol(x))
  for (j in seq_len(ncol(pairs))) {
    difference <- x[pairs[1, j], ] - x[pairs[2, j], ]
    direction <- difference / sqrt(sum(difference^2))
    expected <- expected + tcrossprod(direction)
  }
  expected <- expected / ncol(pairs)
  expect_equal(unclass(estimate), expected, tolerance = 1e-14,
               ignore_attr = TRUE)
  expect_equal(sum(diag(estimate)), 1, tolerance = 1e-14)
  expect_equal(attr(estimate, "n_pairs"), 6)
  expect_equal(
    unclass(spatial_kendall(sweep(x, 2, c(4, -9), "+"))),
    unclass(estimate), tolerance = 1e-14, ignore_attr = TRUE
  )

  extreme <- matrix(c(-1e308, 0, 1e308), ncol = 1)
  extreme_estimate <- spatial_kendall(extreme)
  expect_equal(unclass(extreme_estimate), matrix(1), tolerance = 1e-14,
               ignore_attr = TRUE)
  expect_equal(attr(extreme_estimate, "n_zero_pairs"), 0)

  smallest <- 5e-324
  mixed_scale <- rbind(c(0, 0), c(smallest, smallest), c(1, 0))
  mixed_estimate <- spatial_kendall(mixed_scale)
  mixed_expected <- matrix(c(5 / 6, 1 / 6, 1 / 6, 1 / 6), 2)
  expect_equal(unclass(mixed_estimate), mixed_expected, tolerance = 1e-14,
               ignore_attr = TRUE)
  expect_equal(sum(diag(mixed_estimate)), 1, tolerance = 1e-14)
  expect_equal(attr(mixed_estimate, "n_zero_pairs"), 0)

  epsilon <- .Machine$double.eps
  large_left <- c(1.23456789012345e308, 1.65432109876543e308)
  large_right <- large_left * (1 - epsilon * c(4, 2))
  finite_difference <- large_left - large_right
  expected_direction <- finite_difference / max(abs(finite_difference))
  expected_direction <- expected_direction /
    sqrt(sum(expected_direction^2))
  expect_equal(
    unclass(spatial_kendall(rbind(large_left, large_right))),
    tcrossprod(expected_direction), tolerance = 1e-14, ignore_attr = TRUE
  )
})


test_that("rank covariance is the average rank outer product", {
  set.seed(11)
  x <- matrix(rnorm(24), 8, 3)
  ranks <- spatial_rank(x)
  expect_equal(spatial_rank_covariance(x), crossprod(ranks) / nrow(x),
               tolerance = 1e-14, ignore_attr = TRUE)
})


test_that("shape normalizations fix trace or determinant", {
  shape <- matrix(c(4, 1, 1, 2), 2)
  expect_equal(sum(diag(normalize_shape(shape))), 2, tolerance = 1e-14)
  expect_equal(det(normalize_shape(shape, "determinant")), 1,
               tolerance = 1e-12)
  expect_error(normalize_shape(matrix(c(1, 2, 2, 1), 2)),
               "positive definite")
  extreme <- normalize_shape(diag(c(1e308, 2)))
  expect_true(all(is.finite(extreme)))
  expect_equal(sum(diag(extreme)), 2, tolerance = 1e-14)
})


test_that("Tyler iteration solves its equation and has the stated equivariance", {
  set.seed(12)
  shape <- matrix(c(2, 0.5, 0.5, 1), 2)
  x <- relliptical(
    160, location = c(0, 0), shape = shape,
    radial = function(n) abs(stats::rt(n, df = 2))
  )
  estimate <- tyler_shape(x, center = "none", tol = 1e-10)
  expect_true(isTRUE(attr(estimate, "converged")))
  expect_lte(attr(estimate, "equation_residual"), 1e-9)
  expect_equal(sum(diag(estimate)), 2, tolerance = 1e-12)
  expect_equal(
    unclass(tyler_shape(7 * x, center = "none", tol = 1e-10)),
    unclass(estimate), tolerance = 1e-8, ignore_attr = TRUE
  )
  expect_equal(
    unclass(tyler_shape(1e-12 * x, center = "none", tol = 1e-10)),
    unclass(estimate), tolerance = 1e-8, ignore_attr = TRUE
  )

  transform <- matrix(c(1.3, 0.2, -0.4, 0.9), 2)
  transformed <- tyler_shape(x %*% t(transform), center = "none",
                             tol = 1e-10)
  expected <- normalize_shape(transform %*% unclass(estimate) %*% t(transform))
  expect_equal(unclass(transformed), expected, tolerance = 2e-7,
               ignore_attr = TRUE)

  shifted_x <- sweep(x, 2, c(3, -5), "+")
  default_fit <- tyler_shape(shifted_x, tol = 1e-10)
  affine_x <- sweep(shifted_x %*% t(transform), 2, c(-4, 7), "+")
  affine_fit <- tyler_shape(affine_x, tol = 1e-10)
  default_expected <- normalize_shape(
    transform %*% unclass(default_fit) %*% t(transform)
  )
  expect_equal(unclass(affine_fit), default_expected, tolerance = 2e-7,
               ignore_attr = TRUE)
})


test_that("HR estimates satisfy the joint estimating equations", {
  set.seed(13)
  x <- relliptical(
    180, location = c(2, -1), shape = matrix(c(2, 0.4, 0.4, 1), 2),
    radial = function(n) sqrt(2 * stats::rf(n, 2, 5))
  )
  estimate <- hr_estimator(x, tol = 1e-9, max_iter = 1000)
  expect_s3_class(estimate, "hd_hr")
  expect_true(isTRUE(estimate$converged))
  expect_lte(estimate$location_equation_residual, 1e-6)
  expect_lte(estimate$shape_equation_residual, 1e-5)
  expect_equal(sum(diag(estimate$shape)), 2, tolerance = 1e-12)

  huge_shift <- c(1e8, -2e8)
  translated <- hr_estimator(
    sweep(x, 2, huge_shift, "+"), tol = 1e-9, max_iter = 1000
  )
  expect_equal(translated$location - huge_shift, estimate$location,
               tolerance = 2e-6, ignore_attr = TRUE)
  expect_equal(translated$shape, estimate$shape, tolerance = 2e-6,
               ignore_attr = TRUE)
  expect_lte(translated$location_equation_residual, 1e-6)

  tiny <- hr_estimator(1e-12 * x, tol = 1e-9, max_iter = 1000)
  expect_equal(tiny$location / 1e-12, estimate$location,
               tolerance = 2e-6, ignore_attr = TRUE)
  expect_equal(tiny$shape, estimate$shape, tolerance = 2e-6,
               ignore_attr = TRUE)

  false_stop_data <- rbind(
    c(1e-10, 0), c(1, 0), c(2, 0), c(0, 1), c(0, 2), c(0, -1)
  )
  guarded <- hr_estimator(
    false_stop_data, initial_location = c(0, 0), initial_shape = diag(2),
    tol = 1e-8, max_iter = 2, warn = FALSE
  )
  expect_false(isTRUE(guarded$converged))
  expect_gt(max(guarded$location_equation_residual,
                guarded$shape_equation_residual), 1e-8)
  expect_equal(
    guarded$equation_residual,
    max(guarded$location_equation_residual,
        guarded$shape_equation_residual),
    tolerance = 1e-14
  )
  expect_lte(guarded$best_iteration, guarded$iterations)

  one_dimensional <- hr_estimator(matrix(c(-2, 0, 1, 100), ncol = 1))
  expect_equal(one_dimensional$location, 0.5, ignore_attr = TRUE)
  expect_equal(one_dimensional$shape, matrix(1), ignore_attr = TRUE)
  expect_true(one_dimensional$converged)
  expect_true(is.na(one_dimensional$shape_equation_residual))
  expect_equal(one_dimensional$best_iteration, 0L)
  expect_equal(one_dimensional$equation_residual, 0)

  tied_one_dimensional <- hr_estimator(matrix(c(0, 0, 1), ncol = 1))
  expect_equal(tied_one_dimensional$location, 0, ignore_attr = TRUE)
  expect_equal(tied_one_dimensional$location_equation_residual, 0)
  expect_equal(tied_one_dimensional$equation_residual, 0)

  expect_error(
    hr_estimator(matrix(c(1, 2), ncol = 1), tol = -1),
    "`tol` must be a finite positive number"
  )
  expect_error(
    hr_estimator(matrix(c(1, 2), ncol = 1), max_iter = 0),
    "`max_iter` must be a positive integer"
  )
  expect_error(
    hr_estimator(matrix(c(1, 2), ncol = 1), max_iter = 1.5),
    "`max_iter` must be a positive integer"
  )
  expect_error(
    hr_estimator(matrix(c(1, 2), ncol = 1), zero_tol = -1),
    "`zero_tol` must be a finite non-negative number"
  )

  coincident_start <- rbind(c(0, 0), c(1, 0), c(-0.5, 0.1))
  expect_equal(as.numeric(spatial_median(coincident_start)), c(0, 0),
               tolerance = 1e-12)
  fallback_fit <- hr_estimator(
    coincident_start, tol = 1e-9, max_iter = 10000, warn = FALSE
  )
  expect_true(fallback_fit$converged)
  expect_lte(fallback_fit$equation_residual, 1e-9)
})


test_that("ACG likelihood is invariant to row and shape scaling", {
  set.seed(14)
  x <- matrix(rnorm(30), 10, 3)
  shape <- crossprod(matrix(rnorm(9), 3, 3)) + diag(3)
  expect_equal(acg_loglik(x, shape), acg_loglik(5 * x, 9 * shape),
               tolerance = 1e-12)
  expect_equal(acg_loglik(x, shape), acg_loglik(1e200 * x, shape),
               tolerance = 1e-12)
  expect_true(is.finite(acg_loglik(x, diag(c(1e8, 1, 1e-8)))))
  expect_equal(acg_loglik(matrix(c(1.7e308, 1.7e308), 1), diag(2)),
               0, tolerance = 1e-14)
  expect_error(acg_loglik(rbind(c(0, 0), c(1, 1)), diag(2)), "nonzero")
})


test_that("elliptical and spherical generators honor dimensions and radii", {
  set.seed(15)
  x <- relliptical(20, location = c(a = 1, b = -2),
                   shape = diag(c(1, 0)), radial = rep(2, 20))
  expect_equal(dim(x), c(20L, 2L))
  expect_equal(colnames(x), c("a", "b"))
  expect_equal(x[, 2], rep(-2, 20), tolerance = 1e-14)
  expect_equal(dim(rspherical(7, 4)), c(7L, 4L))
  extreme <- relliptical(3, shape = diag(c(1e308, 2)), radial = rep(1, 3))
  expect_true(all(is.finite(extreme)))
  degenerate <- relliptical(
    4, location = c(1, -2), shape = matrix(0, 2, 2)
  )
  expect_equal(degenerate, matrix(rep(c(1, -2), each = 4), 4, 2),
               ignore_attr = TRUE)
  expect_error(relliptical(2.9), "positive integer")
  expect_error(rspherical(2, 2.9), "positive integer")
})


test_that("input validation fails early and clearly", {
  expect_error(spatial_sign(c(1, NA_real_)), "finite")
  expect_error(spatial_median(matrix(1:4, 2), weights = c(1, -1)),
               "non-negative")
  expect_error(spatial_rank(matrix(1:4, 2), matrix(1:6, 2)),
               "same number of columns")
  expect_error(spatial_kendall(matrix(1:3, 1)), "at least 2 row")
  expect_error(tyler_shape(matrix(rnorm(9), 3, 3)), "n > p")
  expect_error(relliptical(2, shape = matrix(c(1, 2, 2, 1), 2)),
               "positive semidefinite")
})
