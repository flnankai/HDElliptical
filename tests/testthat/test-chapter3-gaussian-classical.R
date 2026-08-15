ch3g_classical_fixture <- function() {
  matrix(c(
    -2.0,  0.4,  1.1,
    -1.2,  1.8, -0.7,
    -0.5, -1.5,  0.3,
     0.2,  0.9,  2.1,
     0.8, -0.6, -1.3,
     1.5,  1.1,  0.5,
     2.2, -1.0,  1.4,
     2.9,  0.1, -0.9
  ), ncol = 3L, byrow = TRUE)
}


ch3g_hp_reference <- function(x, shape0, center, score) {
  n <- nrow(x)
  p <- ncol(x)
  residual <- sweep(x, 2L, center, "-")
  z <- t(backsolve(chol(shape0), t(residual), transpose = TRUE))
  radii <- sqrt(rowSums(z^2))
  directions <- z / radii
  ranks <- rank(radii, ties.method = "first")
  u <- ranks / (n + 1)
  scores <- switch(
    score,
    sign = rep(1, n),
    wilcoxon = u,
    spearman = u^2,
    vdw = qchisq(u, df = p)
  )
  second <- switch(
    score, sign = 1, wilcoxon = 1 / 3, spearman = 1 / 5,
    vdw = p * (p + 2)
  )
  scatter <- crossprod(directions, directions * scores) / n
  contrast <- sum(scatter^2) - sum(diag(scatter))^2 / p
  list(
    statistic = n * p * (p + 2) * contrast / (2 * second),
    scatter = scatter,
    radii = radii,
    directions = directions,
    scores = scores,
    second = second
  )
}


test_that("classical covariance LRT matches the Gaussian likelihood formula", {
  x <- ch3g_classical_fixture()
  sigma0 <- matrix(c(1.7, 0.2, -0.1, 0.2, 1.2, 0.3,
                     -0.1, 0.3, 1.5), 3L, 3L)
  result <- gaussian_covariance_lrt(x, sigma0)
  residual <- sweep(x, 2L, colMeans(x), "-")
  S <- crossprod(residual) / nrow(x)
  relative <- t(backsolve(chol(sigma0), t(residual), transpose = TRUE))
  values <- eigen(crossprod(relative) / nrow(x), symmetric = TRUE,
                  only.values = TRUE)$values
  statistic <- nrow(x) * (sum(values) - sum(log(values)) - ncol(x))

  expect_s3_class(result, "htest")
  expect_s3_class(result, "hd_covariance_test")
  expect_equal(unname(result$statistic), statistic, tolerance = 2e-12)
  expect_equal(result$p.value, pchisq(statistic, 6, lower.tail = FALSE),
               tolerance = 2e-15)
  expect_equal(result$components$relative.eigenvalues, values,
               tolerance = 2e-13)
  expect_equal(result$diagnostics$covariance.divisor, nrow(x))
  expect_equal(
    unname(gaussian_covariance_lrt(x + rep(c(4, -3, 8), each = nrow(x)),
                                  sigma0)$statistic),
    unname(result$statistic), tolerance = 2e-12
  )
  expect_equal(
    unname(gaussian_covariance_lrt(1e120 * x, 1e240 * sigma0)$statistic),
    unname(result$statistic), tolerance = 2e-11
  )

  known <- gaussian_covariance_lrt(x, sigma0, center = FALSE)
  S0 <- crossprod(x) / nrow(x)
  relative0 <- solve(sigma0, S0)
  direct0 <- nrow(x) * (sum(diag(relative0)) -
    as.numeric(determinant(relative0, logarithm = TRUE)$modulus) - ncol(x))
  expect_equal(unname(known$statistic), direct0, tolerance = 2e-12)
})


test_that("Mauchly, John, and Nagao lock their corrected fixed-p scalings", {
  x <- ch3g_classical_fixture()
  n <- nrow(x)
  p <- ncol(x)
  m <- n - 1
  S <- stats::cov(x)
  values <- eigen(S, symmetric = TRUE, only.values = TRUE)$values

  mauchly <- mauchly_sphericity_test(x)
  log.V <- sum(log(values)) - p * log(mean(values))
  rho <- 1 - (2 * p^2 + p + 2) / (6 * p * m)
  expect_equal(unname(mauchly$statistic), -rho * m * log.V,
               tolerance = 2e-13)
  expect_equal(mauchly$components$rho, rho, tolerance = 0)
  expect_identical(mauchly$diagnostics$correction.denominator,
                   "6 * p * residual.df")

  john <- john_sphericity_test(x)
  U <- p * sum(S^2) / sum(diag(S))^2 - 1
  expect_equal(john$raw.statistic[["U"]], U, tolerance = 2e-14)
  expect_equal(unname(john$statistic), m * p * U / 2,
               tolerance = 2e-13)
  expect_identical(john$diagnostics$multiplier, "residual.df * p / 2")

  nagao <- nagao_identity_test(x)
  distance <- sum((S - diag(p))^2)
  expect_equal(unname(nagao$statistic), m * distance / 2,
               tolerance = 2e-13)
  expect_equal(nagao$components$frobenius.distance.squared, distance,
               tolerance = 2e-14)

  Q <- qr.Q(qr(matrix(c(1, 2, 3, -2, 1, 1, 1, -1, 2), 3L)))
  shifted <- sweep(x %*% Q, 2L, c(9, -5, 2), "+")
  for (fun in list(mauchly_sphericity_test, john_sphericity_test)) {
    expect_equal(unname(fun(shifted)$statistic),
                 unname(fun(x)$statistic), tolerance = 3e-12)
    expect_equal(unname(fun(1e-120 * x)$statistic),
                 unname(fun(x)$statistic), tolerance = 3e-12)
  }
  expect_equal(unname(nagao_identity_test(shifted)$statistic),
               unname(nagao$statistic), tolerance = 3e-12)
})


test_that("Hallin-Paindaveine signed ranks match the double-sum quadratic", {
  x <- ch3g_classical_fixture()
  center <- c(0.15, -0.25, 0.35)
  shape0 <- matrix(c(2.0, 0.3, -0.2, 0.3, 1.4, 0.1,
                     -0.2, 0.1, 0.9), 3L, 3L)

  for (score in c("sign", "wilcoxon", "spearman", "vdw")) {
    reference <- ch3g_hp_reference(x, shape0, center, score)
    result <- hallin_paindaveine_shape_test(
      x, shape0 = shape0, center = center, score = score
    )
    expect_equal(unname(result$statistic), reference$statistic,
                 tolerance = 3e-12)
    expect_equal(result$components$weighted.sign.scatter,
                 reference$scatter, tolerance = 3e-14)
    expect_equal(result$components$scores, reference$scores,
                 tolerance = 3e-14)
    expect_equal(result$components$score.second.moment, reference$second,
                 tolerance = 0)
    expect_equal(exp(result$components$log.radii) /
                   max(exp(result$components$log.radii)),
                 reference$radii / max(reference$radii), tolerance = 3e-14)
  }

  sign.fit <- hallin_paindaveine_shape_test(
    x, center = center, score = "sign"
  )
  spatial <- spatial_sign_sphericity_test(x, center = center)
  omega <- sign.fit$components$weighted.sign.scatter
  q.s <- ncol(x) * (sum(omega^2) - 1 / ncol(x))
  expect_equal(unname(spatial$statistic), unname(sign.fit$statistic),
               tolerance = 2e-14)
  expect_equal(spatial$raw.statistic[["Q.S"]], q.s, tolerance = 2e-14)

  A <- matrix(c(2, 0.3, -0.2, -0.4, 1.3, 0.1, 0.2, -0.5, 0.8), 3L)
  transformed <- hallin_paindaveine_shape_test(
    x %*% A + rep(c(5, -7, 3), each = nrow(x)),
    shape0 = t(A) %*% shape0 %*% A,
    center = as.numeric(center %*% A + c(5, -7, 3)),
    score = "wilcoxon"
  )
  original <- hallin_paindaveine_shape_test(
    x, shape0 = shape0, center = center, score = "wilcoxon"
  )
  expect_equal(unname(transformed$statistic), unname(original$statistic),
               tolerance = 2e-11)
  expect_equal(
    unname(hallin_paindaveine_shape_test(
      x, shape0 = 11 * shape0, center = center, score = "vdw"
    )$statistic),
    unname(hallin_paindaveine_shape_test(
      x, shape0 = shape0, center = center, score = "vdw"
    )$statistic), tolerance = 2e-12
  )
})


test_that("classical covariance tests reject undefined boundary inputs", {
  x <- ch3g_classical_fixture()
  expect_error(gaussian_covariance_lrt(x, diag(c(1, 1, 0))),
               "positive definite")
  expect_error(gaussian_covariance_lrt(x[, 1:2, drop = FALSE], diag(3)),
               "dimension")
  expect_error(gaussian_covariance_lrt(matrix(1:9, 3L, 3L), diag(3)),
               "residual degrees")
  expect_error(mauchly_sphericity_test(cbind(1:8, 2 * (1:8))),
               "positive definite")
  expect_error(john_sphericity_test(matrix(1, 5L, 2L)), "non-positive")
  expect_error(nagao_identity_test(x, center = NA), "TRUE or FALSE")
  expect_error(hallin_paindaveine_shape_test(x, shape0 = diag(c(1, -1, 1)),
                                             center = c(0, 0, 0)),
               "positive definite")
  expect_error(hallin_paindaveine_shape_test(
    rbind(c(1, 0), c(-1, 0), c(0, 1), c(0, -1)),
    center = c(0, 0), score = "wilcoxon"
  ), "distinct radii")
  expect_error(hallin_paindaveine_shape_test(
    rbind(c(0, 0), c(1, 0), c(0, 1)), center = c(0, 0), score = "sign"
  ), "exactly zero")
  expect_error(spatial_sign_sphericity_test(matrix(1:5, ncol = 1L), center = 3),
               "two variables")
})
