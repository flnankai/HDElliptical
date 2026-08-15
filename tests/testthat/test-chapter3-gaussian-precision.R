ch3gp_covariance_reference <- function(x, divisor = c("n", "n-1"),
                                       center = TRUE) {
  divisor <- match.arg(divisor)
  z <- if (center) sweep(x, 2L, colMeans(x), "-") else x
  crossprod(z) / if (divisor == "n") nrow(x) else nrow(x) - 1
}


ch3gp_soft_reference <- function(x, lambda) {
  sign(x) * pmax(abs(x) - lambda, 0)
}


ch3gp_ec2_two_reference <- function(x, lambda, tau,
                                    divisor = c("n", "n-1")) {
  covariance <- ch3gp_covariance_reference(x, match.arg(divisor))
  scale <- sqrt(diag(covariance))
  correlation <- covariance / outer(scale, scale)
  value <- ch3gp_soft_reference(correlation[1L, 2L], lambda)
  value <- sign(value) * min(abs(value), 1 - tau)
  fitted.correlation <- matrix(c(1, value, value, 1), 2L, 2L)
  list(
    correlation = fitted.correlation,
    covariance = fitted.correlation * outer(scale, scale),
    sample.covariance = covariance
  )
}


ch3gp_glasso_two_reference <- function(x, lambda,
                                       divisor = c("n", "n-1"),
                                       center = TRUE) {
  scatter <- ch3gp_covariance_reference(
    x, match.arg(divisor), center = center
  )
  dual.covariance <- scatter
  dual.covariance[1L, 2L] <- dual.covariance[2L, 1L] <-
    ch3gp_soft_reference(scatter[1L, 2L], lambda)
  list(
    scatter = scatter,
    dual.covariance = dual.covariance,
    precision = solve(dual.covariance)
  )
}


ch3gp_glasso_kkt_reference <- function(scatter, precision, lambda) {
  gradient <- scatter - solve(precision)
  diagonal <- max(abs(diag(gradient)))
  off.diagonal <- 0
  for (i in seq_len(nrow(precision))) {
    for (j in seq_len(i - 1L)) {
      if (precision[i, j] > 0) {
        component <- abs(gradient[i, j] + lambda)
      } else if (precision[i, j] < 0) {
        component <- abs(gradient[i, j] - lambda)
      } else {
        component <- max(0, abs(gradient[i, j]) - lambda)
      }
      off.diagonal <- max(off.diagonal, component)
    }
  }
  c(diagonal = diagonal, off.diagonal = off.diagonal)
}


ch3gp_smoke <- rbind(
  c(-2, 0), c(-1, -1), c(1, 1), c(2, 0)
)
colnames(ch3gp_smoke) <- c("first", "second")


test_that("EC2 matches the exact bivariate active-constraint formula", {
  x <- rbind(
    c(-2, -1.9), c(-1, -1.2), c(1, 1.1), c(2, 2)
  )
  reference <- ch3gp_ec2_two_reference(
    x, lambda = 0.1, tau = 0.3
  )
  fit <- ec2_covariance(
    x, lambda = 0.1, tau = 0.3, solver_tol = 1e-8
  )

  expect_true(fit$valid)
  expect_equal(fit$sample.covariance, reference$sample.covariance,
               tolerance = 1e-14, ignore_attr = TRUE)
  expect_equal(fit$correlation, reference$correlation,
               tolerance = 2e-7, ignore_attr = TRUE)
  expect_equal(fit$estimate, reference$covariance,
               tolerance = 2e-7, ignore_attr = TRUE)
  expect_equal(unname(diag(fit$correlation)), rep(1, 2))
  expect_gte(min(eigen(fit$correlation, symmetric = TRUE,
                       only.values = TRUE)$values), 0.3 - 1e-8)
  expect_false(fit$diagnostics$solver$soft.threshold.feasible)
})


test_that("EC2 returns the exact soft-threshold solution when feasible", {
  x <- rbind(
    c(-2, -1.9), c(-1, -1.2), c(1, 1.1), c(2, 2)
  )
  reference <- ch3gp_ec2_two_reference(
    x, lambda = 0.6, tau = 0.1
  )
  fit <- ec2_covariance(x, lambda = 0.6, tau = 0.1)

  expect_true(fit$valid)
  expect_equal(fit$correlation, reference$correlation,
               tolerance = 1e-14, ignore_attr = TRUE)
  expect_identical(fit$diagnostics$solver$iterations, 0L)
  expect_true(fit$diagnostics$solver$soft.threshold.feasible)
  expect_lte(fit$diagnostics$solver$kkt.maximum, 1e-12)
})


test_that("EC2 certificates reconstruct its convex KKT system", {
  x <- rbind(
    c(-2, -1.9), c(-1, -1.2), c(1, 1.1), c(2, 2)
  )
  tolerance <- 1e-7
  fit <- ec2_covariance(
    x, lambda = 0, tau = 0.4, solver_tol = tolerance
  )
  certificate <- fit$diagnostics$solver

  expect_true(fit$valid)
  expect_lte(certificate$primal.equality.residual, tolerance)
  expect_lte(certificate$sparse.fixed.point.residual, tolerance)
  expect_lte(certificate$spectral.fixed.point.residual, tolerance)
  expect_lte(certificate$diagonal.residual, tolerance)
  expect_lte(certificate$symmetry.residual, tolerance)
  expect_lte(certificate$eigenvalue.feasibility.violation, tolerance)
  expect_lte(certificate$off.diagonal.stationarity.residual, tolerance)
  expect_lte(certificate$spectral.dual.violation, tolerance)
  expect_lte(certificate$spectral.complementarity.residual, tolerance)
  expect_lte(certificate$kkt.maximum, tolerance)
})


test_that("EC2 makes the finite covariance divisor explicit", {
  n <- nrow(ch3gp_smoke)
  fit.n <- ec2_covariance(
    ch3gp_smoke, lambda = 0.4, tau = 0.1, divisor = "n"
  )
  fit.n1 <- ec2_covariance(
    ch3gp_smoke, lambda = 0.4, tau = 0.1, divisor = "n-1"
  )

  expect_equal(fit.n$correlation, fit.n1$correlation,
               tolerance = 1e-14)
  expect_equal(fit.n1$estimate, fit.n$estimate * n / (n - 1),
               tolerance = 1e-13)
  expect_identical(fit.n$diagnostics$covariance.divisor, n)
  expect_identical(fit.n1$diagnostics$covariance.divisor, n - 1)
})


test_that("EC2 is translation, signed-scale, and permutation equivariant", {
  base <- ec2_covariance(
    ch3gp_smoke, lambda = 0.2, tau = 0.1
  )
  shifted <- ec2_covariance(
    sweep(ch3gp_smoke, 2L, c(11, -7), "+"),
    lambda = 0.2, tau = 0.1
  )
  diagonal <- diag(c(-2, 3))
  scaled <- ec2_covariance(
    ch3gp_smoke %*% diagonal, lambda = 0.2, tau = 0.1
  )
  permutation <- c(2L, 1L)
  permuted <- ec2_covariance(
    ch3gp_smoke[, permutation], lambda = 0.2, tau = 0.1
  )

  expect_equal(shifted$estimate, base$estimate, tolerance = 1e-12)
  expect_equal(scaled$estimate, diagonal %*% base$estimate %*% diagonal,
               tolerance = 2e-7, ignore_attr = TRUE)
  expect_equal(permuted$estimate,
               base$estimate[permutation, permutation],
               tolerance = 1e-12, ignore_attr = TRUE)
})


test_that("EC2 rejects unspecified branches and degenerate inputs", {
  constant <- cbind(ch3gp_smoke[, 1L], 1)

  expect_error(
    ec2_covariance(ch3gp_smoke, 0.2, 0.1, penalty = "adaptive"),
    "review-only"
  )
  expect_error(ec2_covariance(ch3gp_smoke, 0.2, 0), "greater than")
  expect_error(ec2_covariance(ch3gp_smoke, 0.2, 1.1), "at most")
  expect_error(ec2_covariance(constant, 0.2, 0.1),
               "positive empirical variance")
  expect_error(ec2_covariance(ch3gp_smoke, -0.1, 0.1), "at least")
})


test_that("EC2 strict false exposes an uncertified last iterate", {
  x <- rbind(
    c(-2, -1.9), c(-1, -1.2), c(1, 1.1), c(2, 2)
  )
  expect_warning(
    fit <- ec2_covariance(
      x, lambda = 0, tau = 0.5, solver_tol = 1e-14,
      solver_max_iter = 1L, strict = FALSE
    ),
    "without all feasibility"
  )
  expect_false(fit$valid)
  expect_null(fit$estimate)
  expect_null(fit$correlation)
  expect_true(is.list(fit$diagnostics$last.iterate))
  expect_error(
    ec2_covariance(
      x, lambda = 0, tau = 0.5, solver_tol = 1e-14,
      solver_max_iter = 1L, strict = TRUE
    ),
    "without all feasibility"
  )
})


test_that("graphical lasso matches the exact bivariate dual formula", {
  lambda <- 0.2
  reference <- ch3gp_glasso_two_reference(ch3gp_smoke, lambda)
  fit <- gaussian_graphical_lasso(
    ch3gp_smoke, lambda, solver_tol = 1e-8
  )

  expect_true(fit$valid)
  expect_equal(fit$sample.covariance, reference$scatter,
               tolerance = 1e-14, ignore_attr = TRUE)
  expect_equal(fit$estimate, reference$precision,
               tolerance = 3e-7, ignore_attr = TRUE)
  expect_equal(solve(fit$estimate), reference$dual.covariance,
               tolerance = 3e-7, ignore_attr = TRUE)
  expect_gt(min(eigen(fit$estimate, symmetric = TRUE,
                      only.values = TRUE)$values), 0)
})


test_that("graphical lasso objective and KKT are reconstructed directly", {
  lambda <- 0.2
  tolerance <- 1e-7
  fit <- gaussian_graphical_lasso(
    ch3gp_smoke, lambda, solver_tol = tolerance
  )
  scatter <- ch3gp_covariance_reference(ch3gp_smoke)
  kkt <- ch3gp_glasso_kkt_reference(scatter, fit$estimate, lambda)
  logdet <- determinant(fit$estimate, logarithm = TRUE)
  objective <- sum(scatter * t(fit$estimate)) -
    as.numeric(logdet$modulus) +
    lambda * sum(abs(fit$estimate[row(fit$estimate) !=
                                    col(fit$estimate)]))
  history <- fit$diagnostics$solver$objective.history

  expect_lte(
    abs(fit$diagnostics$solver$diagonal.kkt.residual -
          unname(kkt["diagonal"])),
    2e-12
  )
  expect_lte(
    abs(fit$diagnostics$solver$off.diagonal.kkt.residual -
          unname(kkt["off.diagonal"])),
    2e-12
  )
  expect_lte(max(kkt), tolerance)
  expect_equal(fit$diagnostics$solver$objective, objective,
               tolerance = 1e-12)
  expect_true(all(diff(history) <=
                  1e-12 * pmax(1, abs(head(history, -1L)))))
})


test_that("zero-penalty graphical lasso is the exact SPD inverse", {
  fit <- gaussian_graphical_lasso(ch3gp_smoke, lambda = 0)

  expect_true(fit$valid)
  expect_identical(fit$diagnostics$solver$iterations, 0L)
  expect_equal(fit$estimate, solve(fit$sample.covariance),
               tolerance = 1e-14)
  expect_lte(fit$diagnostics$solver$kkt.maximum, 1e-13)
})


test_that("large graphical penalty zeros only off-diagonal entries", {
  scatter <- ch3gp_covariance_reference(ch3gp_smoke)
  lambda <- abs(scatter[1L, 2L]) + 0.1
  fit <- gaussian_graphical_lasso(ch3gp_smoke, lambda)

  expect_true(fit$valid)
  expect_equal(diag(fit$estimate), 1 / diag(scatter),
               tolerance = 1e-14)
  expect_equal(fit$estimate[1L, 2L], 0)
  expect_false(fit$diagnostics$diagonal.penalty)
})


test_that("graphical lasso respects translation, units, and permutation", {
  lambda <- 0.2
  base <- gaussian_graphical_lasso(ch3gp_smoke, lambda)
  shifted <- gaussian_graphical_lasso(
    sweep(ch3gp_smoke, 2L, c(9, -4), "+"), lambda
  )
  factor <- 3
  scaled <- gaussian_graphical_lasso(
    factor * ch3gp_smoke, factor^2 * lambda
  )
  permutation <- c(2L, 1L)
  permuted <- gaussian_graphical_lasso(
    ch3gp_smoke[, permutation], lambda
  )

  expect_equal(shifted$estimate, base$estimate, tolerance = 1e-11)
  expect_equal(scaled$estimate, base$estimate / factor^2,
               tolerance = 3e-7)
  expect_equal(permuted$estimate,
               base$estimate[permutation, permutation],
               tolerance = 3e-7, ignore_attr = TRUE)
})


test_that("graphical lasso does not repair rank or variance failures", {
  singular <- rbind(c(-1, -2, -3), c(0, 0, 0), c(1, 2, 3))
  constant <- cbind(ch3gp_smoke[, 1L], 1)

  expect_error(
    gaussian_graphical_lasso(singular, lambda = 0),
    "strictly positive definite"
  )
  expect_warning(
    invalid <- gaussian_graphical_lasso(
      singular, lambda = 0, strict = FALSE
    ),
    "optimization failed"
  )
  expect_false(invalid$valid)
  expect_null(invalid$estimate)
  expect_error(
    gaussian_graphical_lasso(constant, lambda = 0.2),
    "positive empirical variance"
  )
})


test_that("graphical lasso strict false never promotes an early iterate", {
  expect_warning(
    fit <- gaussian_graphical_lasso(
      ch3gp_smoke, lambda = 0.2, solver_tol = 1e-14,
      solver_max_iter = 1L, strict = FALSE
    ),
    "without every SPD"
  )
  expect_false(fit$valid)
  expect_null(fit$estimate)
  expect_true(is.matrix(fit$diagnostics$last.iterate$precision))
  expect_error(
    gaussian_graphical_lasso(
      ch3gp_smoke, lambda = 0.2, solver_tol = 1e-14,
      solver_max_iter = 1L
    ),
    "without every SPD"
  )
})


test_that("CLIME zero penalty is the exact inverse when it exists", {
  fit <- clime_precision(ch3gp_smoke, lambda = 0)

  expect_true(fit$valid)
  expect_equal(fit$estimate, solve(fit$sample.covariance),
               tolerance = 1e-14)
  expect_equal(fit$diagnostics$raw.solution, fit$estimate,
               tolerance = 1e-14)
  expect_true(all(fit$diagnostics$solver$iterations == 0L))
  expect_lte(max(fit$diagnostics$solver$primal.violation), 1e-13)
  expect_lte(max(fit$diagnostics$solver$relative.gap), 1e-13)
})


test_that("CLIME matches the exact diagonal column programs", {
  x <- rbind(c(-1, 0), c(1, 0), c(0, -2), c(0, 2))
  lambda <- 0.2
  scatter <- ch3gp_covariance_reference(x)
  reference <- diag((1 - lambda) / diag(scatter))
  fit <- clime_precision(x, lambda, solver_tol = 1e-7)

  expect_true(fit$valid)
  expect_equal(fit$diagnostics$raw.solution, reference,
               tolerance = 3e-7, ignore_attr = TRUE)
  expect_equal(fit$estimate, reference,
               tolerance = 3e-7, ignore_attr = TRUE)
})


test_that("CLIME primal-dual certificates reconstruct every raw column", {
  lambda <- 0.2
  tolerance <- 1e-6
  fit <- clime_precision(
    ch3gp_smoke, lambda, solver_tol = tolerance
  )
  scatter <- ch3gp_covariance_reference(ch3gp_smoke)
  raw <- fit$diagnostics$raw.solution
  dual <- fit$diagnostics$dual.solution
  identity <- diag(ncol(raw))

  expect_true(fit$valid)
  expect_lte(max(abs(scatter %*% raw - identity)), lambda + tolerance)
  for (j in seq_len(ncol(raw))) {
    stationarity <- t(scatter) %*% dual[, j]
    nonzero <- raw[, j] != 0
    expect_lte(max(c(0, abs(stationarity[nonzero] +
                            sign(raw[nonzero, j])))), tolerance)
    expect_lte(max(c(0, abs(stationarity[!nonzero]) - 1)), tolerance)
    primal.objective <- sum(abs(raw[, j]))
    dual.objective <- -dual[j, j] - lambda * sum(abs(dual[, j]))
    expect_equal(fit$diagnostics$solver$primal.objective[j],
                 primal.objective, tolerance = 1e-12)
    expect_equal(fit$diagnostics$solver$dual.objective[j],
                 dual.objective, tolerance = 1e-12)
    expect_lte(abs(primal.objective - dual.objective) /
                 max(1, abs(primal.objective), abs(dual.objective)),
               tolerance)
  }
})


test_that("CLIME uses deterministic smaller-absolute symmetrization", {
  fit <- clime_precision(ch3gp_smoke, 0.2, solver_tol = 1e-6)
  raw <- fit$diagnostics$raw.solution
  expected <- raw
  for (i in seq_len(nrow(raw))) {
    for (j in i:ncol(raw)) {
      value <- if (abs(raw[i, j]) <= abs(raw[j, i])) {
        raw[i, j]
      } else {
        raw[j, i]
      }
      expected[i, j] <- expected[j, i] <- value
    }
  }

  expect_equal(fit$estimate, expected, tolerance = 1e-14)
  expect_equal(fit$estimate, t(fit$estimate), tolerance = 0)
  expect_match(fit$diagnostics$symmetrization, "ties retain")
})


test_that("CLIME reports rather than repairs lack of positive definiteness", {
  fit <- clime_precision(ch3gp_smoke, lambda = 1)

  expect_true(fit$valid)
  expect_equal(fit$estimate, matrix(0, 2L, 2L), ignore_attr = TRUE)
  expect_false(fit$diagnostics$symmetrized.positive.definite)
  expect_equal(fit$diagnostics$symmetrized.minimum.eigenvalue, 0)
  expect_match(fit$diagnostics$positive.definiteness.repair,
               "not guaranteed")
})


test_that("CLIME makes its divisor and measurement units explicit", {
  n <- nrow(ch3gp_smoke)
  fit.n <- clime_precision(ch3gp_smoke, 0, divisor = "n")
  fit.n1 <- clime_precision(ch3gp_smoke, 0, divisor = "n-1")
  factor <- 3
  scaled <- clime_precision(factor * ch3gp_smoke, 0.2)
  base <- clime_precision(ch3gp_smoke, 0.2)

  expect_equal(fit.n1$estimate, fit.n$estimate * (n - 1) / n,
               tolerance = 1e-14)
  expect_equal(scaled$estimate, base$estimate / factor^2,
               tolerance = 3e-7)
  expect_identical(fit.n$diagnostics$covariance.divisor, n)
  expect_identical(fit.n1$diagnostics$covariance.divisor, n - 1)
})


test_that("CLIME respects translation and variable permutation", {
  lambda <- 0.2
  base <- clime_precision(ch3gp_smoke, lambda, solver_tol = 1e-6)
  shifted <- clime_precision(
    sweep(ch3gp_smoke, 2L, c(9, -4), "+"),
    lambda, solver_tol = 1e-6
  )
  permutation <- c(2L, 1L)
  permuted <- clime_precision(
    ch3gp_smoke[, permutation], lambda, solver_tol = 1e-6
  )

  expect_equal(shifted$estimate, base$estimate, tolerance = 3e-7)
  expect_equal(permuted$estimate,
               base$estimate[permutation, permutation],
               tolerance = 3e-7, ignore_attr = TRUE)
})


test_that("CLIME strict false returns invalid rather than a pseudo-fit", {
  singular <- rbind(c(-1, -2, -3), c(0, 0, 0), c(1, 2, 3))
  expect_warning(
    fit <- clime_precision(
      singular, lambda = 0, solver_tol = 1e-14,
      solver_max_iter = 1L, strict = FALSE
    ),
    "without every raw-column"
  )
  expect_false(fit$valid)
  expect_null(fit$estimate)
  expect_true(is.matrix(fit$diagnostics$last.symmetrized.iterate))
  expect_error(
    clime_precision(
      singular, lambda = 0, solver_tol = 1e-14,
      solver_max_iter = 1L, strict = TRUE
    ),
    "without every raw-column"
  )
})


test_that("shared Gaussian precision input contracts reject invalid values", {
  constant <- cbind(ch3gp_smoke[, 1L], 1)
  functions <- list(
    function() gaussian_graphical_lasso(constant, 0.2),
    function() clime_precision(constant, 0.2)
  )
  for (fun in functions) {
    expect_error(fun(), "positive empirical variance")
  }
  expect_error(gaussian_graphical_lasso(ch3gp_smoke, -1), "at least")
  expect_error(clime_precision(ch3gp_smoke, -1), "at least")
  expect_error(gaussian_graphical_lasso(ch3gp_smoke, 0.2,
                                        solver_max_iter = 1.5),
               "positive integer")
  expect_error(clime_precision(ch3gp_smoke, 0.2, strict = NA),
               "TRUE or FALSE")
})


test_that("center false uses the displayed uncentered divisor-n scatter", {
  x <- rbind(c(1, 0), c(2, 1), c(0, 2), c(1, 3))
  scatter <- crossprod(x) / nrow(x)
  glasso <- gaussian_graphical_lasso(
    x, lambda = 0, center = FALSE
  )
  clime <- clime_precision(x, lambda = 0, center = FALSE)

  expect_equal(glasso$sample.covariance, scatter,
               tolerance = 1e-14, ignore_attr = TRUE)
  expect_equal(glasso$estimate, solve(scatter),
               tolerance = 1e-14, ignore_attr = TRUE)
  expect_equal(clime$estimate, solve(scatter),
               tolerance = 1e-14, ignore_attr = TRUE)
})


test_that("Gaussian covariance and precision solvers do not touch R RNG", {
  set.seed(3107)
  before <- .Random.seed
  invisible(ec2_covariance(ch3gp_smoke, 0.2, 0.1))
  expect_identical(.Random.seed, before)
  invisible(gaussian_graphical_lasso(ch3gp_smoke, 0.2))
  expect_identical(.Random.seed, before)
  invisible(clime_precision(ch3gp_smoke, 0.2, solver_tol = 1e-6))
  expect_identical(.Random.seed, before)
})
