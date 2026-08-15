osh_fixture <- function() {
  structure(
    matrix(c(
      -2.2,  0.4,  1.1,
      -1.4,  1.8, -0.7,
      -0.8, -1.9,  0.3,
       0.2,  0.9,  2.1,
       1.0, -0.6, -1.4,
       1.7,  1.3,  0.5,
       2.5, -1.2,  1.4,
       3.1,  0.1, -0.9,
      -0.3,  2.4,  1.8
    ), ncol = 3L, byrow = TRUE),
    dimnames = list(paste0("r", 1:9), c("a", "b", "c"))
  )
}

osh_pool_fixture <- function() {
  x <- matrix(c(
    -2.0,  0.0,
    -1.0,  2.0,
     0.0, -2.0,
     1.0,  1.0,
     2.0, -1.0,
     3.0,  2.0,
    -1.0, -2.0,
     2.0,  0.0
  ), ncol = 2L, byrow = TRUE,
  dimnames = list(NULL, c("u", "v")))
  y <- matrix(c(
    -1.5,  1.1,
    -0.2,  2.8,
     0.7, -1.6,
     1.8,  1.7,
     2.9, -0.2,
    -2.2, -1.0,
     0.5,  2.1,
     3.2,  1.0
  ), ncol = 2L, byrow = TRUE,
  dimnames = list(NULL, c("u", "v")))
  list(first = x, second = y)
}

osh_unit <- function(value) {
  scale <- max(abs(value))
  if (scale == 0) return(numeric(length(value)))
  scaled <- value / scale
  scaled / sqrt(sum(scaled^2))
}

osh_signs <- function(x, center) {
  t(vapply(seq_len(nrow(x)), function(index) {
    osh_unit(x[index, ] - center)
  }, numeric(ncol(x))))
}

osh_kurtosis <- function(x) {
  n <- nrow(x)
  p <- ncol(x)
  centered <- sweep(x, 2L, colMeans(x))
  g2 <- colMeans(centered^4) / colMeans(centered^2)^2 - 3
  corrected <- (n - 1) / ((n - 2) * (n - 3)) * ((n + 1) * g2 + 6)
  raw <- mean(corrected) / 3
  list(
    g2 = g2,
    corrected = corrected,
    raw = raw,
    kappa = max(-2 / (p + 2), raw)
  )
}

osh_radial <- function(radii) {
  n <- length(radii)
  m1 <- mean(radii^(-1))
  ratio2 <- mean(radii^(-2)) / m1^2
  ratio3 <- mean(radii^(-3)) / m1^3
  delta <- (2 - 2 * ratio2 + ratio2^2) / n^2 +
    (8 * ratio2 - 6 * ratio2^2 + 2 * ratio2 * ratio3 -
       2 * ratio3) / n^3
  c(ratio2 = ratio2, ratio3 = ratio3, delta = delta)
}

osh_basic_map <- function(lambda, p) {
  if (lambda == 0) return(0)
  0.5 * stats::integrate(function(t) {
    lambda / (1 - t + t * lambda) * (1 - t)^(p / 2 - 1)
  }, lower = 0, upper = 1, rel.tol = 1e-12,
  subdivisions = 1000L)$value
}

osh_active_qp <- function(B, c, lower, convex, tolerance = 1e-8) {
  dimension <- length(c)
  candidates <- list()
  objective <- numeric()
  for (mask in 0:(2^dimension - 1L)) {
    free <- which(as.logical(intToBits(mask)[seq_len(dimension)]))
    active <- setdiff(seq_len(dimension), free)
    value <- lower
    if (convex && length(free) == 0L) {
      if (abs(sum(lower) - 1) > tolerance) next
    } else if (length(free) > 0L) {
      rhs <- c[free]
      if (length(active) > 0L) {
        rhs <- rhs - B[free, active, drop = FALSE] %*% lower[active]
      }
      if (convex) {
        system <- rbind(
          cbind(B[free, free, drop = FALSE], 1),
          c(rep(1, length(free)), 0)
        )
        solution <- tryCatch(
          solve(system, c(rhs, 1 - sum(lower[active]))),
          error = function(condition) NULL
        )
        if (is.null(solution)) next
        value[free] <- solution[seq_along(free)]
      } else {
        solution <- tryCatch(
          solve(B[free, free, drop = FALSE], rhs),
          error = function(condition) NULL
        )
        if (is.null(solution)) next
        value[free] <- solution
      }
    }
    if (any(value < lower - tolerance)) next
    if (convex && abs(sum(value) - 1) > tolerance) next
    gradient <- as.numeric(B %*% value - c)
    if (convex) {
      level <- if (length(free)) mean(gradient[free]) else min(gradient)
      if (length(free) &&
          max(abs(gradient[free] - level)) > tolerance) next
      if (length(active) && any(gradient[active] < level - tolerance)) next
    } else {
      if (length(free) && max(abs(gradient[free])) > tolerance) next
      if (length(active) && any(gradient[active] < -tolerance)) next
    }
    candidates[[length(candidates) + 1L]] <- value
    objective <- c(
      objective,
      0.5 * drop(crossprod(value, B %*% value)) - sum(c * value)
    )
  }
  if (!length(candidates)) stop("No active-set QP reference was feasible.")
  candidates[[which.min(objective)]]
}


test_that("Ell1 follows the published SSCM formula literally", {
  x <- osh_fixture()
  result <- ollila_raninen_shrinkage_covariance(
    x, sphericity = "ell1", tol = 1e-11, keep_signs = TRUE
  )
  n <- nrow(x)
  p <- ncol(x)
  S <- stats::cov(x)
  eta <- sum(diag(S)) / p
  signs <- osh_signs(unname(x), unname(result$components$fitted.location))
  sscm <- crossprod(signs) / n
  gamma.raw <- n / (n - 1) * (p * sum(sscm^2) - p / n)
  gamma <- min(p, max(1, gamma.raw))
  kappa <- osh_kurtosis(unname(x))$kappa
  beta <- (gamma - 1) / (
    gamma - 1 + kappa * (2 * gamma + p) / n +
      (gamma + p) / (n - 1)
  )
  expected <- beta * S + (1 - beta) * eta * diag(p)

  expect_s3_class(result, "hd_covariance_estimator")
  expect_equal(unname(result$sample.covariance), unname(S),
               tolerance = 2e-14)
  expect_equal(result$components$eta, eta, tolerance = 2e-14)
  expect_equal(unname(result$components$fitted.signs), signs,
               tolerance = 2e-13)
  expect_equal(unname(result$components$sscm), sscm, tolerance = 2e-13)
  expect_equal(result$components$ell1$raw, gamma.raw, tolerance = 2e-13)
  expect_equal(result$components$ell1$projected, gamma,
               tolerance = 2e-13)
  expect_equal(result$components$scm.data.weight, beta,
               tolerance = 2e-13)
  expect_equal(result$components$shrinkage.intensity, 1 - beta,
               tolerance = 2e-13)
  expect_equal(unname(result$estimate), unname(expected),
               tolerance = 3e-13)
  expect_equal(result$diagnostics$covariance.divisor, n - 1,
               tolerance = 0)
  expect_identical(result$diagnostics$ridge, "none")
  expect_identical(result$diagnostics$deletion, "none")
})


test_that("Ell2, Ell3, and corrected kurtosis match literal formulas", {
  x <- osh_fixture()
  n <- nrow(x)
  p <- ncol(x)
  S <- stats::cov(x)
  kurtosis <- osh_kurtosis(unname(x))
  a.n <- n / (n + kurtosis$kappa) *
    (n / (n - 1) + kurtosis$kappa)
  b.n <- (kurtosis$kappa + n) * (n - 1)^2 /
    ((n - 2) * (3 * kurtosis$kappa * (n - 1) + n * (n + 1)))
  gamma.raw <- b.n * (
    p * sum(S^2) / sum(diag(S))^2 - a.n * p / n
  )
  gamma <- min(p, max(1, gamma.raw))
  ell2 <- ollila_raninen_shrinkage_covariance(x, "ell2")
  ell3 <- ollila_raninen_shrinkage_covariance(x, "ell3", tol = 1e-11)

  expect_equal(ell2$components$marginal.g2, kurtosis$g2,
               tolerance = 2e-14)
  expect_equal(ell2$components$corrected.marginal.kurtosis,
               kurtosis$corrected, tolerance = 2e-14)
  expect_equal(ell2$components$kappa, kurtosis$kappa,
               tolerance = 2e-14)
  expect_equal(ell2$diagnostics$kappa.raw, kurtosis$raw,
               tolerance = 2e-14)
  expect_equal(ell2$components$ell2$a.n, a.n, tolerance = 2e-14)
  expect_equal(ell2$components$ell2$b.n, b.n, tolerance = 2e-14)
  expect_equal(ell2$components$ell2$raw, gamma.raw,
               tolerance = 3e-14)
  expect_equal(ell2$components$ell2$projected, gamma,
               tolerance = 3e-14)
  expect_identical(
    ell3$components$sphericity.selected,
    if (ell3$components$ell1$projected < gamma) "ell1" else "ell2"
  )
  expect_equal(
    ell3$components$sphericity,
    min(ell3$components$ell1$projected, gamma), tolerance = 1e-14
  )

  projected <- .os_kurtosis_from_g2(rep(-3, p), n = 10L, p = p)
  expect_equal(projected$estimate, -2 / (p + 2), tolerance = 0)
  expect_true(projected$below.lower.before.projection)
  expect_identical(projected$official.code.boundary.repair,
                   "not used (no 0.99 multiplier)")

  one <- matrix(c(-3, -1, 2, 4, 7), ncol = 1L)
  tie <- ollila_raninen_shrinkage_covariance(one, "ell3")
  expect_equal(tie$components$sphericity, 1, tolerance = 0)
  expect_identical(tie$components$sphericity.selected, "ell2")
  expect_true(tie$diagnostics$p.one.sphericity.known)
})


test_that("shrinkage covariance has its exact finite-sample equivariances", {
  x <- osh_fixture()
  shift <- c(1e4, -7, 13)
  signed.permutation <- matrix(c(
    0, 0, -1,
    1, 0,  0,
    0, 1,  0
  ), 3L, 3L, byrow = TRUE)
  fit <- ollila_raninen_shrinkage_covariance(x, "ell1", tol = 1e-11)
  translated <- ollila_raninen_shrinkage_covariance(
    sweep(x, 2L, shift, "+"), "ell1", tol = 1e-11
  )
  transformed.x <- unname(x) %*% signed.permutation
  colnames(transformed.x) <- colnames(x)
  transformed <- ollila_raninen_shrinkage_covariance(
    transformed.x, "ell1", tol = 1e-11
  )
  scaled <- ollila_raninen_shrinkage_covariance(
    x * 1e100, "ell1", tol = 1e-11
  )
  tiny <- ollila_raninen_shrinkage_covariance(
    x * 1e-100, "ell1", tol = 1e-11
  )

  expect_equal(translated$estimate, fit$estimate, tolerance = 2e-10)
  expect_equal(
    unname(transformed$estimate),
    t(signed.permutation) %*% unname(fit$estimate) %*% signed.permutation,
    tolerance = 3e-12
  )
  expect_equal(unname(scaled$estimate) / 1e200,
               unname(fit$estimate), tolerance = 3e-12)
  expect_equal(unname(tiny$estimate) / 1e-200,
               unname(fit$estimate), tolerance = 3e-12)
  expect_equal(scaled$components$scm.data.weight,
               fit$components$scm.data.weight, tolerance = 2e-13)
  expect_equal(tiny$components$scm.data.weight,
               fit$components$scm.data.weight, tolerance = 2e-13)
})


test_that("RSSCM matches the literal data weight and handles a equals one", {
  x <- osh_fixture()
  center <- colMeans(x)
  signs <- osh_signs(unname(x), center)
  sscm <- crossprod(signs) / nrow(x)
  V <- ncol(x) * sscm
  a <- sum(V^2) / ncol(x)
  alpha.raw <- (nrow(x) / (nrow(x) - 1) *
                  (a - ncol(x) / nrow(x)) - 1) / (a - 1)
  alpha <- min(1, max(0, alpha.raw))
  expected <- alpha * V + (1 - alpha) * diag(ncol(x))
  result <- regularized_spatial_sign_covariance(
    x, center = "mean", keep_signs = TRUE
  )

  expect_s3_class(result, "hd_covariance_estimator")
  expect_equal(unname(result$components$fitted.signs), unname(signs),
               tolerance = 2e-13)
  expect_equal(unname(result$components$sscm), unname(sscm),
               tolerance = 2e-13)
  expect_equal(result$components$a, a, tolerance = 2e-13)
  expect_equal(result$components$alpha.raw, alpha.raw,
               tolerance = 2e-13)
  expect_equal(result$components$alpha.data, alpha, tolerance = 2e-13)
  expect_equal(unname(result$estimate), unname(expected), tolerance = 3e-13)
  expect_equal(sum(diag(result$estimate)), ncol(x), tolerance = 2e-14)

  axes <- matrix(c(1, 0, -1, 0, 0, 1, 0, -1),
                 ncol = 2L, byrow = TRUE)
  signal.free <- regularized_spatial_sign_covariance(
    axes, center = "none"
  )
  expect_equal(signal.free$components$a, 1, tolerance = 0)
  expect_true(is.na(signal.free$components$alpha.raw))
  expect_equal(signal.free$components$alpha.data, 0, tolerance = 0)
  expect_equal(unname(signal.free$estimate), diag(2), tolerance = 0)
  expect_true(signal.free$diagnostics$signal.free)
})


test_that("BASIC inversion agrees with direct integration", {
  x <- osh_fixture()
  result <- basic_shape(
    x, center = "mean", keep_signs = TRUE,
    quadrature_order = 48L, inversion_tol = 1e-11,
    integration_tol = 1e-10
  )
  p <- ncol(x)
  mapped.reference <- vapply(
    as.numeric(result$components$lambda.raw), osh_basic_map,
    numeric(1), p = p
  )
  expected <- result$components$eigenvectors %*%
    (as.numeric(result$components$lambda.normalized) *
       t(result$components$eigenvectors))

  expect_s3_class(result, "hd_shape_estimator")
  expect_equal(mapped.reference,
               as.numeric(result$components$sscm.eigenvalues),
               tolerance = 5e-10)
  expect_equal(as.numeric(result$components$delta.mapped),
               as.numeric(result$components$sscm.eigenvalues),
               tolerance = 2e-11)
  expect_equal(unname(result$estimate), unname(expected),
               tolerance = 2e-13)
  expect_equal(sum(diag(result$estimate)), p, tolerance = 2e-13)
  expect_true(result$diagnostics$inversion.converged)
  expect_lte(result$diagnostics$maximum.integration.error, 1e-10)
  expect_lte(result$diagnostics$maximum.inversion.error, 1e-11)
  expect_identical(result$diagnostics$extrapolation, "none")
  expect_identical(result$diagnostics$numerical.floor, "none")

  rank.deficient <- rbind(
    c(1, 2, 3, 4),
    c(-2, 1, -1, 2)
  )
  boundary <- basic_shape(
    rank.deficient, center = "none", quadrature_order = 48L
  )
  expect_equal(sum(boundary$diagnostics$boundary.zero), 2)
  expect_equal(sum(diag(boundary$estimate)), 4, tolerance = 2e-13)
  expect_equal(sum(eigen(boundary$estimate, symmetric = TRUE,
                         only.values = TRUE)$values < -1e-12), 0)

  one <- matrix(c(-2, -1, 1, 3), ncol = 1L)
  expect_equal(unname(basic_shape(one)$estimate), matrix(1), tolerance = 0)
  expect_true(basic_shape(one)$diagnostics$p.one.special.case)
  expect_error(
    cpp_ollila_basic_inverse(c(1.1, 0), 2L, 16L,
                             1e-10, 1e-10, 100L, 100L),
    "lie in"
  )
  expect_error(
    cpp_ollila_basic_inverse(c(0.7, 0.3), 2L, 16L,
                             1e-14, 1e-10, 1L, 100L),
    "did not converge"
  )

  high.p.delta <- c(.15, rep(.85 / 24, 24))
  high.p <- cpp_ollila_basic_inverse(
    high.p.delta, 25L, 64L, 1e-11, 1e-11, 100L, 100L
  )
  expect_equal(as.numeric(high.p$delta_mapped), high.p.delta,
               tolerance = 2e-11)
  expect_equal(sum(high.p$lambda_normalized), 25, tolerance = 2e-13)
  expect_lte(high.p$max_integration_error, 1e-11)
})


test_that("BASICS applies BASIC to RSSCM eigenvalues divided by p", {
  x <- osh_fixture()
  result <- basics_shape(
    x, center = "mean", quadrature_order = 48L,
    inversion_tol = 1e-11
  )
  delta <- result$components$alpha.data *
    eigen(result$components$sscm, symmetric = TRUE,
          only.values = TRUE)$values +
    (1 - result$components$alpha.data) / ncol(x)
  mapped.reference <- vapply(
    as.numeric(result$components$lambda.raw), osh_basic_map,
    numeric(1), p = ncol(x)
  )

  expect_equal(sort(result$components$rsscm.eigenvalues.divided.by.p),
               sort(delta), tolerance = 2e-13)
  expect_equal(mapped.reference,
               as.numeric(result$components$rsscm.eigenvalues.divided.by.p),
               tolerance = 5e-10)
  expect_equal(sum(diag(result$estimate)), ncol(x), tolerance = 2e-13)
  expect_true(result$diagnostics$inversion.converged)
  expect_identical(result$diagnostics$extrapolation, "none")

  axes <- matrix(c(1, 0, -1, 0, 0, 1, 0, -1),
                 ncol = 2L, byrow = TRUE)
  signal.free <- basics_shape(axes, center = "none")
  expect_equal(unname(signal.free$estimate), diag(2), tolerance = 2e-10)
  expect_equal(signal.free$components$alpha.data, 0, tolerance = 0)
  expect_true(signal.free$diagnostics$signal.free)
  one <- matrix(c(-2, -1, 1, 3), ncol = 1L)
  expect_equal(unname(basics_shape(one)$estimate), matrix(1), tolerance = 0)
})


test_that("spatial-sign shape estimators have the stated equivariances", {
  x <- osh_fixture()
  theta <- 0.37
  Q <- matrix(c(
    cos(theta), -sin(theta), 0,
    sin(theta),  cos(theta), 0,
    0,           0,          1
  ), 3L, 3L, byrow = TRUE)
  shifted <- sweep(x, 2L, c(7, -9, 4), "+")
  rotated <- unname(x) %*% Q
  colnames(rotated) <- colnames(x)
  methods <- list(
    rsscm = function(z) regularized_spatial_sign_covariance(
      z, center = "mean"
    ),
    basic = function(z) basic_shape(
      z, center = "mean", quadrature_order = 40L
    ),
    basics = function(z) basics_shape(
      z, center = "mean", quadrature_order = 40L
    )
  )
  for (method in methods) {
    fit <- method(x)
    translated <- method(shifted)
    orthogonal <- method(rotated)
    huge <- method(x * 1e250)
    tiny <- method(x * 1e-250)
    expect_equal(translated$estimate, fit$estimate, tolerance = 2e-12)
    expect_equal(
      unname(orthogonal$estimate),
      t(Q) %*% unname(fit$estimate) %*% Q,
      tolerance = 3e-9
    )
    expect_equal(unname(huge$estimate), unname(fit$estimate),
                 tolerance = 3e-10)
    expect_equal(unname(tiny$estimate), unname(fit$estimate),
                 tolerance = 3e-10)
  }
})


test_that("pooling components match literal real-valued formulas", {
  data <- osh_pool_fixture()
  result <- linear_pool_covariance(
    data, method = "linear", identity = FALSE,
    tol = 1e-11, solver_tol = 1e-12, solver_max_iter = 50000L
  )
  K <- length(data)
  p <- ncol(data[[1]])
  eta <- gamma <- gamma0 <- kappa <- delta.diag <- numeric(K)
  sscm <- vector("list", K)
  covariance <- vector("list", K)
  for (index in seq_len(K)) {
    x <- data[[index]]
    n <- nrow(x)
    covariance[[index]] <- stats::cov(x)
    eta[index] <- sum(diag(covariance[[index]])) / p
    signs <- osh_signs(x, result$spatial.centers[[index]])
    radii <- sqrt(rowSums(sweep(x, 2L,
                                result$spatial.centers[[index]])^2))
    sscm[[index]] <- crossprod(signs) / n
    gamma0[index] <- p * n / (n - 1) *
      (sum(sscm[[index]]^2) - 1 / n)
    radial <- osh_radial(radii)
    gamma.raw <- gamma0[index] - p * radial[["delta"]]
    gamma[index] <- min(p, max(1, gamma.raw))
    kappa[index] <- osh_kurtosis(x)$kappa
    tr.sigma2 <- p * eta[index]^2 * gamma[index]
    tau1 <- 1 / (n - 1) + kappa[index] / n
    tau2 <- kappa[index] / n
    delta.diag[index] <- (
      tau1 * sum(diag(covariance[[index]]))^2 +
        (tau1 + tau2) * tr.sigma2
    ) / p
    expect_equal(result$components$radial.bias[[index]]$ratio2,
                 radial[["ratio2"]], tolerance = 3e-11)
    expect_equal(result$components$radial.bias[[index]]$ratio3,
                 radial[["ratio3"]], tolerance = 3e-11)
    expect_equal(result$components$radial.bias[[index]]$delta,
                 radial[["delta"]], tolerance = 3e-11)
  }
  C <- matrix(0, K, K)
  diag(C) <- eta^2 * gamma
  C[1, 2] <- C[2, 1] <- sum(sscm[[1]] * sscm[[2]]) *
    sum(diag(covariance[[1]])) * sum(diag(covariance[[2]])) / p
  Delta <- diag(delta.diag, K)

  expect_s3_class(result, "linear_pool_covariance")
  expect_equal(unname(result$sample.covariances[[1]]),
               unname(covariance[[1]]),
               tolerance = 3e-14)
  expect_equal(unname(result$sample.covariances[[2]]),
               unname(covariance[[2]]),
               tolerance = 3e-14)
  expect_equal(unname(result$spatial.sign.covariances[[1]]),
               unname(sscm[[1]]),
               tolerance = 3e-12)
  expect_equal(unname(result$spatial.sign.covariances[[2]]),
               unname(sscm[[2]]),
               tolerance = 3e-12)
  expect_equal(unname(result$eta), eta, tolerance = 3e-14)
  expect_equal(unname(result$gamma0.raw), gamma0, tolerance = 3e-12)
  expect_equal(unname(result$gamma), gamma, tolerance = 3e-12)
  expect_equal(unname(result$kappa), kappa, tolerance = 3e-14)
  expect_equal(unname(result$C), C, tolerance = 3e-11)
  expect_equal(unname(result$Delta), Delta, tolerance = 3e-11)
  expect_length(result$diagnostics$spatial.median, K)
  expect_true(all(vapply(
    result$diagnostics$spatial.median,
    function(value) isTRUE(value$converged), logical(1)
  )))
  expect_length(
    result$diagnostics$residual.radius.scale.relative.to.global, K
  )
  for (target in seq_len(K)) {
    reference <- osh_active_qp(
      result$components$objective.C.scaled +
        result$components$objective.Delta.scaled,
      result$components$objective.C.scaled[, target],
      numeric(K), convex = FALSE
    )
    expect_equal(unname(result$coefficients[, target]), reference,
                 tolerance = 2e-8)
    expected <- Reduce(`+`, Map(`*`, result$sample.covariances,
                                result$coefficients[, target]))
    expect_equal(unname(result$estimates[[target]]), unname(expected),
                 tolerance = 3e-13)
    expect_true(result$solver[[target]]$converged)
    expect_lte(result$solver[[target]]$feasibility_residual, 1e-12)
    expect_identical(result$solver[[target]]$ridge, "none")
  }
})


test_that("convex pooling and identity lower bounds match an active-set QP", {
  data <- osh_pool_fixture()
  result <- linear_pool_covariance(
    data, method = "convex", identity = TRUE, identity_lower = c(.1, .2),
    tol = 1e-11, solver_tol = 1e-12, solver_max_iter = 50000L
  )
  B <- result$components$objective.C.scaled +
    result$components$objective.Delta.scaled
  for (target in seq_along(data)) {
    lower <- c(0, 0, c(.1, .2)[target])
    reference <- osh_active_qp(
      B, result$components$objective.C.scaled[, target], lower,
      convex = TRUE
    )
    expect_equal(unname(result$coefficients[, target]), reference,
                 tolerance = 3e-8)
    expect_equal(sum(result$coefficients[, target]), 1, tolerance = 2e-12)
    expect_gte(result$coefficients["identity", target], lower[3] - 1e-12)
    expect_true(result$solver[[target]]$converged)
  }
  expect_equal(unname(diag(result$C.extended)[3]), 1, tolerance = 0)
  expect_equal(result$Delta.extended[3, 3], 0, tolerance = 0)
  expect_identical(result$diagnostics$identity.default.difference,
                   "default is zero; official MATLAB convenience default 1e-8 is not imposed")

  singular <- cpp_ollila_pool_qp(
    matrix(0, 3, 3), numeric(3), numeric(3), TRUE, 1e-12, 10L
  )
  expect_true(singular$diagnostics$converged)
  expect_equal(sum(singular$solution), 1, tolerance = 2e-15)
  expect_equal(as.numeric(singular$solution), rep(1 / 3, 3),
               tolerance = 2e-15)
  expect_identical(singular$diagnostics$regularization, "none")
  expect_error(
    cpp_ollila_pool_qp(-diag(2), c(1, 1), c(0, 0), FALSE,
                       1e-10, 100L),
    "not positive semidefinite"
  )
})


test_that("pooling has the appropriate target-dependent equivariances", {
  data <- osh_pool_fixture()
  fit <- linear_pool_covariance(
    data, identity = FALSE, tol = 1e-11, solver_tol = 1e-11
  )
  shifted.data <- Map(function(value, shift) {
    sweep(value, 2L, shift, "+")
  }, data, list(c(13, -7), c(-5, 11)))
  shifted <- linear_pool_covariance(
    shifted.data, identity = FALSE, tol = 1e-11, solver_tol = 1e-11
  )
  scaled <- linear_pool_covariance(
    lapply(data, `*`, 1e40), identity = FALSE,
    tol = 1e-11, solver_tol = 1e-11
  )
  signed.permutation <- matrix(c(0, -1, 1, 0), 2L, 2L, byrow = TRUE)
  transformed.data <- lapply(data, function(value) {
    answer <- unname(value) %*% signed.permutation
    colnames(answer) <- colnames(value)
    answer
  })
  transformed <- linear_pool_covariance(
    transformed.data, identity = FALSE, tol = 1e-11,
    solver_tol = 1e-11
  )

  expect_equal(shifted$coefficients, fit$coefficients, tolerance = 3e-9)
  expect_equal(scaled$coefficients, fit$coefficients, tolerance = 3e-9)
  expect_equal(transformed$coefficients, fit$coefficients, tolerance = 3e-9)
  for (target in seq_along(data)) {
    expect_equal(shifted$estimates[[target]], fit$estimates[[target]],
                 tolerance = 3e-9)
    expect_equal(unname(scaled$estimates[[target]]) / 1e80,
                 unname(fit$estimates[[target]]), tolerance = 3e-9)
    expect_equal(
      unname(transformed$estimates[[target]]),
      t(signed.permutation) %*% unname(fit$estimates[[target]]) %*%
        signed.permutation,
      tolerance = 3e-9
    )
  }
  expect_error(
    linear_pool_covariance(lapply(data, `*`, 1e100), identity = FALSE),
    "fourth-order scale"
  )
  expect_error(
    linear_pool_covariance(lapply(data, `*`, 1e-100), identity = FALSE),
    "fourth-order scale"
  )
})


test_that("all undefined and nonconverged cases are explicit", {
  x <- osh_fixture()
  zero.variable <- cbind(x[, 1:2], constant = 1)
  zero.residual <- rbind(
    c(-2, 0), c(2, 0), c(0, -2), c(0, 2), c(0, 0)
  )

  expect_error(ollila_raninen_shrinkage_covariance(x[1:3, ], "ell2"),
               "at least")
  expect_error(ollila_raninen_shrinkage_covariance(zero.variable, "ell2"),
               "zero-variance")
  expect_error(
    regularized_spatial_sign_covariance(zero.residual, center = c(0, 0)),
    "radius is exactly zero"
  )
  expect_error(basic_shape(zero.residual, center = c(0, 0)),
               "radius is exactly zero")
  expect_error(basics_shape(zero.residual, center = c(0, 0)),
               "radius is exactly zero")
  expect_error(
    regularized_spatial_sign_covariance(x, max_iter = 1L, strict = TRUE),
    "did not converge"
  )
  expect_warning(
    nonstrict <- regularized_spatial_sign_covariance(
      x, max_iter = 1L, strict = FALSE
    ),
    "last finite iterate"
  )
  expect_false(nonstrict$diagnostics$center.details$converged)

  data <- osh_pool_fixture()
  expect_error(linear_pool_covariance(list(data[[1]][1:3, ])),
               "at least")
  renamed <- data
  colnames(renamed[[2]]) <- c("v", "u")
  expect_error(linear_pool_covariance(renamed), "Column names")
  expect_error(linear_pool_covariance(data, identity = FALSE,
                                      identity_lower = .1),
               "must be zero")
  expect_error(linear_pool_covariance(data, method = "convex",
                                      identity_lower = 1.1),
               "at most one")
  expect_error(
    linear_pool_covariance(data, identity = FALSE, solver_tol = 1e-16,
                           solver_max_iter = 1L, strict = TRUE),
    "did not converge"
  )
  warning.messages <- character()
  uncalibrated <- withCallingHandlers(
    linear_pool_covariance(
      data, identity = FALSE, solver_tol = 1e-16,
      solver_max_iter = 1L, strict = FALSE
    ),
    warning = function(condition) {
      warning.messages <<- c(warning.messages, conditionMessage(condition))
      invokeRestart("muffleWarning")
    }
  )
  expect_length(warning.messages, 2L)
  expect_true(all(grepl("uncalibrated", warning.messages, fixed = TRUE)))
  expect_false(uncalibrated$diagnostics$calibrated)
  expect_true(all(vapply(
    uncalibrated$solver, function(value) !value$converged, logical(1)
  )))
  expect_identical(uncalibrated$diagnostics$ridge, "none")
  expect_identical(uncalibrated$diagnostics$pseudoinverse, "none")
})
