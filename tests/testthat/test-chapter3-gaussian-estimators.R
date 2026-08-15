ch3g_estimator_fixture <- function() {
  matrix(c(
    -2.0,  0.2,  1.0, -0.4,
    -1.3,  1.5, -0.7,  0.8,
    -0.6, -1.2,  0.4,  1.4,
     0.1,  0.7,  1.9, -1.1,
     0.7, -0.5, -1.4,  0.2,
     1.4,  1.0,  0.6,  1.7,
     2.1, -0.9,  1.3, -0.8,
     2.8,  0.3, -0.9,  0.5
  ), ncol = 4L, byrow = TRUE)
}


ch3g_threshold_reference <- function(z, lambda, rule, shape) {
  magnitude <- abs(z)
  sign.z <- sign(z)
  result <- matrix(0, nrow(z), ncol(z))
  active <- magnitude > lambda
  if (rule == "hard") {
    result[active] <- z[active]
  } else if (rule == "soft") {
    result[active] <- sign.z[active] *
      (magnitude[active] - lambda[active])
  } else if (rule == "scad") {
    first <- active & magnitude <= 2 * lambda
    middle <- magnitude > 2 * lambda & magnitude <= shape * lambda
    last <- magnitude > shape * lambda
    result[first] <- sign.z[first] * (magnitude[first] - lambda[first])
    result[middle] <- ((shape - 1) * z[middle] -
      sign.z[middle] * shape * lambda[middle]) / (shape - 2)
    result[last] <- z[last]
  } else {
    result[active] <- z[active] *
      (1 - (lambda[active] / magnitude[active])^(shape + 1))
  }
  result
}


test_that("Bickel-Levina uses the primary divisor-n covariance and threshold", {
  x <- ch3g_estimator_fixture()
  z <- sweep(x, 2L, colMeans(x), "-")
  S.n <- crossprod(z) / nrow(x)
  lambda <- 0.28
  reference <- ifelse(abs(S.n) > lambda, S.n, 0)
  fit <- bickel_levina_covariance_threshold(x, threshold = lambda)

  expect_s3_class(fit, "hd_covariance_estimator")
  expect_equal(fit$sample.covariance, S.n, tolerance = 2e-14)
  expect_equal(fit$estimate, reference, tolerance = 2e-14)
  expect_equal(fit$threshold, matrix(lambda, 4L, 4L), tolerance = 0)
  expect_identical(fit$diagnostics$divisor.convention, "n")
  expect_equal(fit$diagnostics$covariance.divisor, nrow(x))

  constant <- 0.73
  automatic <- bickel_levina_covariance_threshold(x, constant = constant)
  automatic.lambda <- constant * sqrt(log(ncol(x)) / nrow(x))
  expect_equal(automatic$threshold,
               matrix(automatic.lambda, 4L, 4L), tolerance = 2e-15)
  expect_equal(automatic$estimate,
               ifelse(abs(S.n) > automatic.lambda, S.n, 0),
               tolerance = 2e-14)

  n.minus.one <- bickel_levina_covariance_threshold(
    x, threshold = lambda, divisor = "n-1"
  )
  S.unbiased <- crossprod(z) / (nrow(x) - 1)
  expect_equal(n.minus.one$sample.covariance, S.unbiased,
               tolerance = 2e-14)
  expect_identical(n.minus.one$diagnostics$divisor.convention, "n-1")
  expect_false(isTRUE(all.equal(fit$sample.covariance,
                                n.minus.one$sample.covariance)))
})


test_that("Rothman-Levina-Zhu locks hard, soft, SCAD, and adaptive lasso", {
  x <- ch3g_estimator_fixture()
  S <- crossprod(sweep(x, 2L, colMeans(x), "-")) / nrow(x)
  lambda <- matrix(0.31, ncol(x), ncol(x))
  settings <- list(
    hard = 1,
    soft = 1,
    scad = 3.7,
    adaptive_lasso = 1.6
  )
  for (rule in names(settings)) {
    shape <- settings[[rule]]
    fit <- rothman_levina_zhu_covariance_threshold(
      x, threshold = lambda[1L], rule = rule,
      scad_a = if (rule == "scad") shape else 3.7,
      adaptive_eta = if (rule == "adaptive_lasso") shape else 1
    )
    reference <- ch3g_threshold_reference(S, lambda, rule, shape)
    expect_equal(fit$estimate, reference, tolerance = 3e-14)
    expect_identical(fit$components$rule, rule)
  }
  adaptive <- rothman_levina_zhu_covariance_threshold(
    x, threshold = 0.31, rule = "adaptive_lasso", adaptive_eta = 1.6
  )
  expect_equal(
    adaptive$estimate,
    sign(S) * pmax(abs(S) - 0.31^(2.6) * abs(S)^(-1.6), 0),
    tolerance = 3e-14
  )
  expect_identical(adaptive$diagnostics$book.MCP.attribution,
                   "not used; primary method is adaptive lasso")
})


test_that("Cai-Liu variability and thresholds match literal divisor-n loops", {
  x <- ch3g_estimator_fixture()
  n <- nrow(x)
  p <- ncol(x)
  z <- sweep(x, 2L, colMeans(x), "-")
  S <- crossprod(z) / n
  theta <- matrix(0, p, p)
  for (i in seq_len(p)) for (j in seq_len(p)) {
    theta[i, j] <- mean((z[, i] * z[, j] - S[i, j])^2)
  }
  delta <- 1.35
  lambda <- delta * sqrt(theta * log(p) / n)
  fit <- cai_liu_adaptive_covariance_threshold(x, delta = delta)

  expect_equal(fit$sample.covariance, S, tolerance = 3e-14)
  expect_equal(fit$components$theta, theta, tolerance = 4e-14)
  expect_equal(fit$threshold, lambda, tolerance = 4e-14)
  expect_equal(fit$estimate, ifelse(abs(S) > lambda, S, 0),
               tolerance = 4e-14)
  expect_true(fit$diagnostics$primary.finite.sample.formula)
  expect_equal(fit$diagnostics$theta.outer.divisor, n)

  mixed <- cai_liu_adaptive_covariance_threshold(
    x, delta = delta, divisor = "n-1"
  )
  S.mixed <- crossprod(z) / (n - 1)
  theta.mixed <- matrix(0, p, p)
  for (i in seq_len(p)) for (j in seq_len(p)) {
    theta.mixed[i, j] <- mean((z[, i] * z[, j] - S.mixed[i, j])^2)
  }
  expect_equal(mixed$sample.covariance, S.mixed, tolerance = 3e-14)
  expect_equal(mixed$components$theta, theta.mixed, tolerance = 4e-14)
  expect_false(mixed$diagnostics$primary.finite.sample.formula)
})


test_that("POET equals PCA removal plus thresholded orthogonal complement", {
  x <- ch3g_estimator_fixture()
  n <- nrow(x)
  p <- ncol(x)
  z <- sweep(x, 2L, colMeans(x), "-")
  decomposition <- svd(z, nu = 0L, nv = min(n - 1L, p))
  vector <- decomposition$v[, 1L, drop = FALSE]
  value <- decomposition$d[1L]^2 / n
  low.rank <- tcrossprod(vector * sqrt(value))
  residual.data <- z - z %*% vector %*% t(vector)
  residual.covariance <- crossprod(residual.data) / n
  tau <- 0.37
  thresholds <- tau * sqrt(outer(diag(residual.covariance),
                                 diag(residual.covariance)))
  thresholded <- ch3g_threshold_reference(
    residual.covariance, thresholds, "hard", 1
  )
  diag(thresholded) <- diag(residual.covariance)
  fit <- poet_covariance(
    x, factors = 1, threshold = tau,
    threshold_type = "correlation", rule = "hard"
  )

  expect_equal(fit$sample.covariance, crossprod(z) / n,
               tolerance = 4e-14)
  expect_equal(fit$components$low.rank, low.rank, tolerance = 5e-14)
  expect_equal(fit$components$orthogonal.complement,
               residual.covariance, tolerance = 5e-14)
  expect_equal(fit$threshold, thresholds, tolerance = 5e-14)
  expect_equal(fit$components$thresholded.orthogonal.complement,
               thresholded, tolerance = 5e-14)
  expect_equal(fit$estimate, low.rank + thresholded, tolerance = 6e-14)
  expect_equal(fit$components$leading.eigenvalues, value,
               tolerance = 3e-14)
  expect_identical(fit$diagnostics$factor.count.source, "supplied")
  expect_identical(fit$diagnostics$positive.definiteness.repair, "none")

  variance.fit <- poet_covariance(
    x, factors = 0, threshold = 0.9, threshold_type = "variance"
  )
  theta <- matrix(0, p, p)
  S <- crossprod(z) / n
  for (i in seq_len(p)) for (j in seq_len(p)) {
    theta[i, j] <- mean((z[, i] * z[, j] - S[i, j])^2)
  }
  omega <- 1 / sqrt(p) + sqrt(log(p) / n)
  expect_equal(variance.fit$components$theta, theta, tolerance = 4e-14)
  expect_equal(variance.fit$threshold, 0.9 * sqrt(theta) * omega,
               tolerance = 5e-14)
})


test_that("covariance estimators have their stated scale and permutation behavior", {
  x <- ch3g_estimator_fixture()
  permutation <- c(4, 2, 1, 3)
  scale <- 1e120
  lambda <- 0.22

  base <- rothman_levina_zhu_covariance_threshold(
    x, threshold = lambda, rule = "soft"
  )
  large <- rothman_levina_zhu_covariance_threshold(
    scale * x, threshold = scale^2 * lambda, rule = "soft"
  )
  expect_equal(large$estimate / scale^2, base$estimate, tolerance = 3e-13)
  translated <- rothman_levina_zhu_covariance_threshold(
    sweep(x, 2L, c(7, -9, 2, 5), "+"), threshold = lambda, rule = "soft"
  )
  expect_equal(translated$estimate, base$estimate, tolerance = 3e-14)
  permuted <- rothman_levina_zhu_covariance_threshold(
    x[, permutation], threshold = lambda, rule = "soft"
  )
  expect_equal(permuted$estimate,
               base$estimate[permutation, permutation], tolerance = 3e-14)

  cai <- cai_liu_adaptive_covariance_threshold(x, delta = 1.1)
  cai.large <- cai_liu_adaptive_covariance_threshold(scale * x, delta = 1.1)
  expect_equal(cai.large$estimate / scale^2, cai$estimate,
               tolerance = 4e-13)
  poet <- poet_covariance(x, factors = 1, threshold = 0.4)
  poet.small <- poet_covariance(1e-120 * x, factors = 1, threshold = 0.4)
  expect_equal(poet.small$estimate / 1e-240, poet$estimate,
               tolerance = 8e-13)
})


test_that("covariance estimators enforce tuning and solver boundaries", {
  x <- ch3g_estimator_fixture()
  expect_error(bickel_levina_covariance_threshold(x, threshold = -1),
               "non-negative")
  expect_error(bickel_levina_covariance_threshold(x, divisor = "unbiased"))
  expect_error(rothman_levina_zhu_covariance_threshold(x, rule = "mcp"))
  expect_error(rothman_levina_zhu_covariance_threshold(
    x, threshold = 0.2, rule = "scad", scad_a = 2
  ), "greater than 2")
  expect_error(rothman_levina_zhu_covariance_threshold(
    x, threshold = 0.2, rule = "adaptive_lasso", adaptive_eta = -0.1
  ), "non-negative")
  expect_error(cai_liu_adaptive_covariance_threshold(x, delta = Inf),
               "finite")
  expect_error(poet_covariance(x, factors = 1.5, threshold = 0.2),
               "integer")
  expect_error(poet_covariance(x, factors = 5, threshold = 0.2),
               "rank bound")
  expect_error(poet_covariance(x, factors = 1, threshold = -0.2),
               "non-negative")
  expect_error(poet_covariance(matrix(0, 5L, 3L), factors = 1,
                               threshold = 0.2), "nonzero")
})
