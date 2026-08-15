erht_fixture <- function() {
  rbind(
    c(0.7, -1.1, 0.2, 1.3, -0.4),
    c(-0.4, 0.8, 1.2, -0.6, 0.9),
    c(1.1, 0.3, -0.9, 0.5, -1.2),
    c(-1.2, -0.5, 0.7, 1.0, 0.3),
    c(0.2, 1.4, -0.4, -1.1, 0.6),
    c(1.5, -0.2, 0.5, 0.1, 1.1),
    c(-0.8, 0.6, -1.3, 0.9, -0.7),
    c(0.5, -1.4, 0.9, -0.2, 1.4),
    c(-1.5, 0.1, 0.4, 1.2, -0.1),
    c(0.9, 1.0, -0.6, -0.8, 0.5),
    c(-0.1, -0.9, 1.4, 0.4, -1.0)
  )
}


erht_reference <- function(x, mu, theta, rho) {
  n <- nrow(x)
  p <- ncol(x)
  residual <- sweep(x, 2L, theta)
  radii <- sqrt(rowSums(residual^2))
  U <- residual / radii
  R <- p * crossprod(U) / n
  Q <- solve(R + rho * diag(p))
  difference <- theta - mu
  T <- as.numeric(n * crossprod(difference, Q %*% difference))
  A <- p * U %*% Q %*% t(U) / n
  w <- sqrt(p) / radii
  kappa <- mean(diag(A))
  e <- mean(w)
  tt <- mean(w^2)
  b1 <- mean(diag(A) * w)
  b2 <- mean(diag(A) * w^2)
  D <- (e - b1)^2 + kappa * (tt - b2)
  mu.fun <- kappa / D
  A2 <- A^2
  diag(A2) <- 0
  psi <- function(a, b) {
    sum(A2 * outer(w^a, w^b)) / n
  }
  psi.values <- c(
    psi00 = psi(0, 0), psi01 = psi(0, 1),
    psi02 = psi(0, 2), psi11 = psi(1, 1),
    psi12 = psi(1, 2), psi22 = psi(2, 2)
  )
  Gamma <- matrix(c(
    2 * psi.values["psi00"], 2 * psi.values["psi01"],
    2 * psi.values["psi11"],
    2 * psi.values["psi01"],
    psi.values["psi02"] + psi.values["psi11"],
    2 * psi.values["psi12"],
    2 * psi.values["psi11"], 2 * psi.values["psi12"],
    2 * psi.values["psi22"]
  ), 3L, 3L, byrow = TRUE)
  g <- c((e - b1)^2, 2 * kappa * (e - b1), kappa^2) / D^2
  sigma2 <- as.numeric(crossprod(g, Gamma %*% g))
  center <- n * mu.fun
  variance <- n * sigma2
  list(
    T = T, A = A, R = R, Q = Q, radii = radii, w = w,
    kappa = kappa, e = e, t = tt, b1 = b1, b2 = b2,
    D = D, mu = mu.fun, psi = psi.values, Gamma = Gamma, g = g,
    sigma2 = sigma2, center = center, variance = variance,
    Z = (T - center) / sqrt(variance)
  )
}


erht_svd_reference <- function(residual, difference, rho) {
  n <- nrow(residual)
  p <- ncol(residual)
  radii <- sqrt(rowSums(residual^2))
  signs <- residual / radii
  decomposition <- svd(signs, nu = min(n, p), nv = min(n, p))
  eigenvalues <- (p / n) * decomposition$d^2
  coefficients <- as.numeric(crossprod(decomposition$v, difference))
  null.part <- difference - as.numeric(decomposition$v %*% coefficients)
  T <- n * (
    sum(coefficients^2 / (eigenvalues + rho)) +
      sum(null.part^2) / rho
  )

  a <- eigenvalues / (eigenvalues + rho)
  b <- rho / (eigenvalues + rho)
  if (n <= p) {
    B <- sweep(decomposition$u, 2L, b, "*") %*% t(decomposition$u)
  } else {
    B <- diag(n) -
      sweep(decomposition$u, 2L, a, "*") %*% t(decomposition$u)
  }
  w <- sqrt(p) / radii
  kappa <- sum(a) / n
  e.minus.b1 <- mean(diag(B) * w)
  t.minus.b2 <- mean(diag(B) * w^2)
  D <- e.minus.b1^2 + kappa * t.minus.b2
  mu.fun <- kappa / D
  B2 <- B^2
  diag(B2) <- 0
  psi <- function(left, right) {
    sum(B2 * outer(w^left, w^right)) / n
  }
  psi.values <- c(
    psi00 = psi(0, 0), psi01 = psi(0, 1),
    psi02 = psi(0, 2), psi11 = psi(1, 1),
    psi12 = psi(1, 2), psi22 = psi(2, 2)
  )
  Gamma <- matrix(c(
    2 * psi.values["psi00"], 2 * psi.values["psi01"],
    2 * psi.values["psi11"],
    2 * psi.values["psi01"],
    psi.values["psi02"] + psi.values["psi11"],
    2 * psi.values["psi12"],
    2 * psi.values["psi11"], 2 * psi.values["psi12"],
    2 * psi.values["psi22"]
  ), 3L, 3L, byrow = TRUE)
  gradient <- c(
    e.minus.b1^2,
    2 * kappa * e.minus.b1,
    kappa^2
  ) / D^2
  sigma2 <- as.numeric(crossprod(gradient, Gamma %*% gradient))
  center <- n * mu.fun
  variance <- n * sigma2
  list(
    T = T, center = center, variance = variance,
    Z = (T - center) / sqrt(variance),
    kappa = kappa, D = D, psi = psi.values,
    g = gradient, sigma2 = sigma2
  )
}


test_that("fixed-ridge ERHT matches the literal feasible equations", {
  x <- erht_fixture()
  mu <- c(0.1, -0.2, 0.05, 0.15, -0.1)
  fit <- elliptical_regularized_hotelling_test(x, mu, rho = 0.4)
  ref <- erht_reference(x, mu, fit$estimate, 0.4)

  expect_s3_class(fit, "htest")
  expect_s3_class(fit, "hd_location_test")
  expect_equal(unname(fit$statistic), ref$Z, tolerance = 2e-10)
  expect_equal(unname(fit$raw.statistic), ref$T, tolerance = 2e-10)
  expect_equal(fit$components$marginal$center, ref$center,
               tolerance = 2e-10)
  expect_equal(fit$components$marginal$variance, ref$variance,
               tolerance = 2e-10)
  expect_equal(fit$components$e, ref$e, tolerance = 2e-10)
  expect_equal(fit$components$t, ref$t, tolerance = 2e-10)
  expect_equal(unname(fit$components$kappa), ref$kappa,
               tolerance = 2e-10)
  expect_equal(fit$components$b1, ref$b1, tolerance = 2e-10)
  expect_equal(fit$components$b2, ref$b2, tolerance = 2e-10)
  expect_equal(fit$components$D, ref$D, tolerance = 2e-10)
  expect_equal(fit$components$mu.functional, ref$mu, tolerance = 2e-10)
  expect_equal(fit$components$sigma.D.squared, ref$sigma2,
               tolerance = 2e-10)
  expect_equal(as.numeric(fit$components$psi), unname(ref$psi),
               tolerance = 3e-10)
  expect_equal(as.numeric(fit$components$g), unname(ref$g),
               tolerance = 3e-10)
  expect_equal(fit$p.value, pnorm(ref$Z, lower.tail = FALSE),
               tolerance = 1e-14)
  expect_identical(fit$alternative, "two.sided")
  expect_false(fit$diagnostics$extra.regularization)
  expect_match(fit$diagnostics$bartlett.center.correction,
               "does not specify")
})


test_that("ERHT Cauchy aggregation matches all fixed-ridge margins", {
  x <- erht_fixture()
  rho <- c(0.2, 0.5, 1.0)
  weights <- c(0.2, 0.3, 0.5)
  fit <- elliptical_regularized_hotelling_cauchy_test(
    x, rho = rho, weights = weights
  )
  individual <- lapply(rho, function(value) {
    elliptical_regularized_hotelling_test(x, rho = value)
  })
  z <- vapply(individual, function(object) unname(object$statistic), numeric(1))
  p.value <- pnorm(z, lower.tail = FALSE)
  direct <- sum(weights * tan(pi * (0.5 - p.value)))
  direct.p <- 0.5 - atan(direct) / pi

  expect_equal(fit$components$marginal$Z, z, tolerance = 2e-11)
  expect_equal(fit$components$marginal$p.value, p.value,
               tolerance = 1e-14)
  expect_equal(fit$components$cauchy[["statistic"]], direct,
               tolerance = 2e-12)
  expect_equal(fit$p.value, direct.p, tolerance = 2e-14)
  expect_equal(unname(fit$components$cauchy.weights), weights)
  expect_false(fit$null.distribution$fixed.level.exact)
  expect_equal(fit$parameter[["K"]], 3)
})


test_that("ERHT is translation, scalar, orthogonal, and row invariant", {
  x <- erht_fixture()
  mu <- c(0.1, -0.2, 0.05, 0.15, -0.1)
  base <- elliptical_regularized_hotelling_test(x, mu, rho = 0.6)
  shift <- c(1e8, -2e8, 3e8, -4e8, 5e8)
  translated <- elliptical_regularized_hotelling_test(
    sweep(x, 2L, shift, "+"), mu + shift, rho = 0.6
  )
  scaled <- elliptical_regularized_hotelling_test(
    -3.25 * x, -3.25 * mu, rho = 0.6
  )
  set.seed(22)
  q <- qr.Q(qr(matrix(rnorm(25), 5, 5)))
  rotated <- elliptical_regularized_hotelling_test(
    x %*% q, as.numeric(mu %*% q), rho = 0.6
  )
  permuted <- elliptical_regularized_hotelling_test(
    x[c(7, 2, 11, 1, 8, 5, 3, 10, 4, 9, 6), ], mu, rho = 0.6
  )

  expect_equal(translated$statistic, base$statistic, tolerance = 5e-7)
  expect_equal(scaled$statistic, base$statistic, tolerance = 2e-9)
  expect_equal(rotated$statistic, base$statistic, tolerance = 2e-8)
  expect_equal(permuted$statistic, base$statistic, tolerance = 2e-9)
  expect_equal(unname(scaled$raw.statistic),
               3.25^2 * unname(base$raw.statistic), tolerance = 2e-8)
  expect_equal(scaled$components$marginal$center,
               3.25^2 * base$components$marginal$center,
               tolerance = 2e-8)
  expect_equal(scaled$components$marginal$variance,
               3.25^4 * base$components$marginal$variance,
               tolerance = 2e-7)
})


test_that("ERHT companion route works when p exceeds n", {
  base <- erht_fixture()
  x <- cbind(base, base^2, sin(base), cos(base), base[, 1:3]^3)
  fit <- elliptical_regularized_hotelling_test(
    x, rho = 0.7, keep_companion = TRUE
  )
  expect_true(is.finite(unname(fit$statistic)))
  expect_equal(length(fit$components$companion$gram.eigenvalues), nrow(x))
  expect_equal(dim(fit$components$companion$matrices[[1L]]),
               c(nrow(x), nrow(x)))
  expect_equal(sum(fit$components$companion$eigenvalues[1, ]),
               unname(fit$components$kappa) * nrow(x), tolerance = 2e-10)
})


test_that("ERHT ridge quadratic remains stable for small positive rho", {
  set.seed(7)
  x <- matrix(rnorm(60), 20, 3)
  difference <- c(0.7, -0.2, 0.4)
  residual <- x
  radii <- sqrt(rowSums(residual^2))
  signs <- residual / radii
  shape <- 3 * crossprod(signs) / 20

  native <- cpp_erht_grid(x, difference, 1e-16, FALSE)
  direct <- as.numeric(
    20 * crossprod(
      difference,
      solve(shape + 1e-16 * diag(3), difference)
    )
  )
  expect_equal(
    as.numeric(native$raw_statistic_scaled), direct, tolerance = 2e-12
  )
  expect_match(native$quadratic_route, "row-space/null-space")
  expect_false(native$quadratic_roundoff_clipped)
})


test_that("ERHT ridge labels distinguish adjacent doubles", {
  rho <- c(0.1, 0.1 + .Machine$double.eps)
  fit <- elliptical_regularized_hotelling_cauchy_test(
    erht_fixture(), rho = rho
  )
  expect_length(unique(rownames(fit$components$marginal)), 2L)
  expect_length(unique(names(fit$components$cauchy.weights)), 2L)
  expect_length(unique(names(fit$components$kappa)), 2L)
})


test_that("ERHT normalizes finite positive Cauchy weights stably", {
  x <- erht_fixture()
  fit <- elliptical_regularized_hotelling_cauchy_test(
    x, rho = c(0.2, 0.6),
    weights = rep(.Machine$double.xmax, 2L)
  )
  expect_equal(unname(fit$components$cauchy.weights), c(0.5, 0.5))
  expect_error(
    elliptical_regularized_hotelling_cauchy_test(
      x, rho = c(0.2, 0.6),
      weights = c(1e308, .Machine$double.xmin * .Machine$double.eps)
    ),
    "representable double-precision range"
  )
  expect_error(
    elliptical_regularized_hotelling_cauchy_test(
      x, rho = c(0.2, 0.6, 1),
      weights = c(1, 1, .Machine$double.xmin * .Machine$double.eps)
    ),
    "representable double-precision range"
  )
})


test_that("ERHT retains small positive spectral directions", {
  set.seed(17)
  base <- matrix(rnorm(90), 30, 3)
  x <- sweep(base, 2L, c(1, 1e-8, 1e-8), "*")
  radii <- sqrt(rowSums(x^2))
  signs <- x / radii
  shape <- 3 * crossprod(signs) / 30
  difference <- c(0, 1, 0)

  for (rho in c(1e-16, 1e-14)) {
    native <- cpp_erht_grid(x, difference, rho, FALSE)
    reference <- erht_reference(
      x, mu = -difference, theta = numeric(3), rho = rho
    )
    expect_equal(
      as.numeric(native$raw_statistic_scaled), reference$T,
      tolerance = 3e-11
    )
    expect_equal(
      as.numeric(native$center_scaled), reference$center,
      tolerance = 3e-10
    )
    expect_equal(
      as.numeric(native$variance_scaled), reference$variance,
      tolerance = 3e-9
    )
    expect_equal(
      as.numeric(native$z), reference$Z, tolerance = 3e-9
    )
  }
})


test_that("ERHT complement calibration is stable when p exceeds n", {
  set.seed(9)
  n <- 4L
  p <- 19L
  x <- matrix(rnorm(n * p), n, p)
  difference <- rnorm(p)

  for (rho in c(1e-16, 1e-18)) {
    native <- cpp_erht_grid(x, difference, rho, FALSE)
    reference <- erht_svd_reference(x, difference, rho)
    expect_equal(
      as.numeric(native$raw_statistic_scaled), reference$T,
      tolerance = 4e-12
    )
    expect_equal(
      as.numeric(native$center_scaled), reference$center,
      tolerance = 4e-11
    )
    expect_equal(
      as.numeric(native$variance_scaled), reference$variance,
      tolerance = 4e-10
    )
    expect_equal(
      as.numeric(native$z), reference$Z, tolerance = 4e-10
    )
  }
})


test_that("ERHT preserves calibration under extreme common scales", {
  x <- erht_fixture()
  base <- elliptical_regularized_hotelling_cauchy_test(
    x, rho = c(0.2, 0.6, 1)
  )
  tiny <- elliptical_regularized_hotelling_cauchy_test(
    1e-150 * x, rho = c(0.2, 0.6, 1)
  )
  huge <- elliptical_regularized_hotelling_cauchy_test(
    1e150 * x, rho = c(0.2, 0.6, 1)
  )
  expect_equal(tiny$components$marginal$Z,
               base$components$marginal$Z, tolerance = 3e-8)
  expect_equal(huge$components$marginal$Z,
               base$components$marginal$Z, tolerance = 3e-8)
  expect_equal(tiny$p.value, base$p.value, tolerance = 2e-9)
  expect_equal(huge$p.value, base$p.value, tolerance = 2e-9)
  expect_true(all(is.finite(tiny$components$scaled$variance)))
  expect_true(all(is.finite(huge$components$scaled$variance)))
})


test_that("ERHT separates residual and distant-null numerical scales", {
  x <- erht_fixture()
  mu <- rep(1e78, ncol(x))
  fit <- elliptical_regularized_hotelling_test(x, mu, rho = 0.4)
  ref <- erht_reference(x, mu, fit$estimate, 0.4)

  expect_true(is.finite(unname(fit$raw.statistic)))
  expect_true(is.finite(unname(fit$statistic)))
  expect_equal(unname(fit$raw.statistic), ref$T, tolerance = 3e-10)
  expect_equal(unname(fit$statistic), ref$Z, tolerance = 3e-10)
  expect_gt(
    fit$diagnostics$internal.null.displacement.scale,
    1e70 * fit$diagnostics$internal.residual.scale
  )
  expect_true(all(is.finite(fit$components$scaled$psi)))
  expect_true(all(is.finite(fit$components$scaled$variance)))
  expect_length(fit$diagnostics$quadratic.roundoff.clipped, 1L)
})


test_that("ERHT stable Cauchy aggregation handles extreme log tails", {
  same <- .erht_cauchy_grid(
    list(
      p.value = c(0, 0),
      log.p.value = c(-1000, -1000),
      log.one.minus.p.value = c(0, 0)
    ),
    c(0.5, 0.5)
  )
  expect_equal(same$log.absolute.statistic, 1000 - log(pi),
               tolerance = 1e-12)
  expect_equal(same$log.p.value, -1000, tolerance = 1e-10)

  cancelling <- .erht_cauchy_grid(
    list(
      p.value = c(0.2, 0.8),
      log.p.value = log(c(0.2, 0.8)),
      log.one.minus.p.value = log(c(0.8, 0.2))
    ),
    c(0.5, 0.5)
  )
  expect_equal(cancelling$statistic, 0, tolerance = 1e-15)
  expect_equal(cancelling$p.value, 0.5, tolerance = 1e-15)
})


test_that("ERHT rejects undefined and unsupported inputs without repair", {
  x <- erht_fixture()
  expect_error(
    elliptical_regularized_hotelling_test(x[, 1, drop = FALSE]),
    "at least two variables"
  )
  expect_error(
    elliptical_regularized_hotelling_test(x[1:2, ]),
    "at least 3 row"
  )
  expect_error(elliptical_regularized_hotelling_test(x, rho = 0),
               "strictly positive")
  expect_error(elliptical_regularized_hotelling_test(x, rho = c(0.2, 0.4)),
               "exactly one")
  expect_error(
    elliptical_regularized_hotelling_cauchy_test(x, rho = 0.5),
    "at least two"
  )
  expect_error(
    elliptical_regularized_hotelling_cauchy_test(
      x, rho = c(0.5, 0.2)
    ),
    "strictly increasing"
  )
  expect_error(
    elliptical_regularized_hotelling_cauchy_test(
      x, rho = c(0.2, 0.5), weights = c(1, 0)
    ),
    "strictly positive"
  )
  coincident <- rbind(c(0, 0), c(1, 0), c(-0.5, 0.1))
  expect_error(
    elliptical_regularized_hotelling_test(coincident, max_iter = 10000L),
    "inverse-distance calibration is undefined"
  )
  expect_error(
    elliptical_regularized_hotelling_test(x, max_iter = 1L),
    "did not satisfy"
  )
})
