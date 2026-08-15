ch3pp_fixture_x <- function(n = 12L) {
  i <- seq_len(n)
  cbind(
    sin(0.7 * i) + 0.08 * i,
    cos(1.1 * i) - 0.03 * i,
    ((7 * i) %% 13 - 6) / 4,
    sin(0.3 * i) + cos(0.9 * i),
    ((5 * i + 2) %% 17 - 8) / 5
  )
}


ch3pp_fixture_y <- function(n = 15L) {
  i <- seq_len(n)
  base <- cbind(
    cos(0.6 * i) + 0.04 * i,
    sin(1.2 * i) - 0.02 * i,
    ((9 * i + 1) %% 19 - 9) / 5,
    cos(0.35 * i) - sin(0.8 * i),
    ((4 * i + 3) %% 23 - 11) / 6
  )
  sweep(base, 2L, c(0.15, -0.2, 0.1, 0.05, -0.1), "+")
}


ch3pp_sign_reference <- function(z, zero_tol = 0) {
  z <- as.numeric(z)
  scale <- max(abs(z))
  if (scale == 0) return(numeric(length(z)))
  standardized <- z / scale
  norm.standardized <- sqrt(sum(standardized^2))
  radius <- scale * norm.standardized
  if (radius <= zero_tol) return(numeric(length(z)))
  standardized / norm.standardized
}


ch3pp_centered_signs_reference <- function(x, center) {
  t(vapply(
    seq_len(nrow(x)),
    function(i) ch3pp_sign_reference(x[i, ] - center),
    numeric(ncol(x))
  ))
}


ch3pp_cheng_delta_reference <- function(radii, n) {
  inverse <- min(radii) / radii
  m1 <- mean(inverse)
  m2 <- mean(inverse^2)
  m3 <- mean(inverse^3)
  d2 <- 2 - 2 * m2 / m1^2 + m2^2 / m1^4
  d3 <- -6 * m2^2 / m1^4 + 2 * m2 * m3 / m1^5 +
    8 * m2 / m1^2 - 2 * m3 / m1^3
  list(value = d2 / n^2 + d3 / n^3, d2 = d2, d3 = d3)
}


ch3pp_cheng_reference <- function(x, y, center.x, center.y) {
  ux <- ch3pp_centered_signs_reference(x, center.x)
  uy <- ch3pp_centered_signs_reference(y, center.y)
  gx <- tcrossprod(ux)
  gy <- tcrossprod(uy)
  gxy <- ux %*% t(uy)
  n1 <- nrow(x)
  n2 <- nrow(y)
  p <- ncol(x)
  A <- (sum(gx^2) - sum(diag(gx)^2)) / (n1 * (n1 - 1))
  B <- (sum(gy^2) - sum(diag(gy)^2)) / (n2 * (n2 - 1))
  C <- sum(gxy^2) / (n1 * n2)
  rx <- sqrt(rowSums(sweep(x, 2L, center.x, "-")^2))
  ry <- sqrt(rowSums(sweep(y, 2L, center.y, "-")^2))
  d1 <- ch3pp_cheng_delta_reference(rx, n1)
  d2 <- ch3pp_cheng_delta_reference(ry, n2)
  delta <- d1$value + d2$value
  T <- p * (A + B - 2 * C)
  q <- p / (n1 + n2) *
    (n1 * (A - d1$value) + n2 * (B - d2$value))
  standard.error <- 2 * (1 / n1 + 1 / n2) / (p + 2) * (p * q)
  list(
    ux = ux, uy = uy, A = A, B = B, C = C, T = T,
    d1 = d1, d2 = d2, delta = delta, q = q,
    standard.error = standard.error,
    variance = standard.error^2,
    z = (T - p * delta) / standard.error
  )
}


ch3pp_rank_reference <- function(x, y, zero_tol = 0) {
  n1 <- nrow(x)
  n2 <- nrow(y)
  p <- ncol(x)
  ux <- vector("list", n1 * n1)
  uy <- vector("list", n2 * n2)
  index.x <- function(i, j) (i - 1L) * n1 + j
  index.y <- function(i, j) (i - 1L) * n2 + j
  for (i in seq_len(n1)) for (j in seq_len(n1)) {
    ux[[index.x(i, j)]] <- ch3pp_sign_reference(x[i, ] - x[j, ], zero_tol)
  }
  for (i in seq_len(n2)) for (j in seq_len(n2)) {
    uy[[index.y(i, j)]] <- ch3pp_sign_reference(y[i, ] - y[j, ], zero_tol)
  }
  sum.x <- 0
  for (i in seq_len(n1)) for (j in seq_len(n1)) {
    if (j == i) next
    for (k in seq_len(n1)) for (l in seq_len(n1)) {
      if (length(unique(c(i, j, k, l))) != 4L) next
      sum.x <- sum.x + sum(
        ux[[index.x(i, j)]] * ux[[index.x(k, l)]]
      )^2
    }
  }
  sum.y <- 0
  for (i in seq_len(n2)) for (j in seq_len(n2)) {
    if (j == i) next
    for (k in seq_len(n2)) for (l in seq_len(n2)) {
      if (length(unique(c(i, j, k, l))) != 4L) next
      sum.y <- sum.y + sum(
        uy[[index.y(i, j)]] * uy[[index.y(k, l)]]
      )^2
    }
  }
  sum.cross <- 0
  for (i in seq_len(n1)) for (j in seq_len(n1)) {
    if (j == i) next
    for (k in seq_len(n2)) for (l in seq_len(n2)) {
      if (l == k) next
      sum.cross <- sum.cross + sum(
        ux[[index.x(i, j)]] * uy[[index.y(k, l)]]
      )^2
    }
  }
  denom.x <- n1 * (n1 - 1) * (n1 - 2) * (n1 - 3)
  denom.y <- n2 * (n2 - 1) * (n2 - 2) * (n2 - 3)
  denom.cross <- n1 * (n1 - 1) * n2 * (n2 - 1)
  A1 <- sum.x / denom.x
  A2 <- sum.y / denom.y
  C12 <- sum.cross / denom.cross
  T <- p * (A1 + A2 - 2 * C12)
  variance <- 4 * (1 / n1 + 1 / n2)^2 /
    ((n1 + n2) * (p + 2)^2) * p^2 * (n1 * A1 + n2 * A2)
  list(
    sum.x = sum.x, sum.y = sum.y, sum.cross = sum.cross,
    A1 = A1, A2 = A2, C12 = C12, T = T,
    variance = variance, standard.error = sqrt(variance),
    z = T / sqrt(variance)
  )
}


test_that("Cheng feasible statistic matches a literal unequal-n reference", {
  x <- ch3pp_fixture_x(12)
  y <- ch3pp_fixture_y(15)
  result <- cheng_sscm_equality_test(x, y, tol = 1e-10)
  reference <- ch3pp_cheng_reference(
    x, y, result$components$center.x, result$components$center.y
  )

  expect_s3_class(result, "hd_proportionality_test")
  expect_true(result$diagnostics$calibrated)
  expect_equal(result$components$signs.x, reference$ux, tolerance = 2e-12)
  expect_equal(result$components$signs.y, reference$uy, tolerance = 2e-12)
  expect_equal(result$components$A, reference$A, tolerance = 2e-13)
  expect_equal(result$components$B, reference$B, tolerance = 2e-13)
  expect_equal(result$components$C, reference$C, tolerance = 2e-13)
  expect_equal(result$components$T.SSCM, reference$T, tolerance = 2e-12)
  expect_equal(result$components$delta1.hat, reference$d1$value,
               tolerance = 2e-12)
  expect_equal(result$components$delta2.hat, reference$d2$value,
               tolerance = 2e-12)
  expect_equal(result$components$delta.hat, reference$delta,
               tolerance = 2e-12)
  expect_equal(result$components$trace.Lambda2.over.p.hat, reference$q,
               tolerance = 2e-12)
  expect_equal(result$components$null.standard.error.hat,
               reference$standard.error, tolerance = 2e-12)
  expect_equal(result$components$null.variance.hat,
               reference$variance, tolerance = 2e-12)
  expect_equal(unname(result$statistic), reference$z, tolerance = 2e-12)
  expect_equal(result$p.value, pnorm(reference$z, lower.tail = FALSE),
               tolerance = 2e-14)
})


test_that("Cheng standard-error and squared-variance fields cannot be confused", {
  result <- cheng_sscm_equality_test(
    ch3pp_fixture_x(13), ch3pp_fixture_y(16), tol = 1e-10
  )
  cmp <- result$components
  n1 <- 13
  n2 <- 16
  p <- 5
  expected.denominator <- 2 * (1 / n1 + 1 / n2) / (p + 2) *
    cmp$trace.Lambda2.hat

  expect_equal(cmp$null.standard.error.hat, expected.denominator,
               tolerance = 1e-14)
  expect_equal(cmp$null.variance.hat, expected.denominator^2,
               tolerance = 1e-14)
  expect_equal(unname(result$statistic),
               cmp$centered.statistic / cmp$null.standard.error.hat,
               tolerance = 1e-14)
  expect_identical(result$diagnostics$denominator,
                   "null.standard.error.hat")
  expect_match(result$diagnostics$variance, "square")
  expect_match(result$diagnostics$plug.in.bias.condition, "3/2")
  expect_match(result$diagnostics$book.erratum, "Remark 3.2")
  expect_match(result$diagnostics$rejection.comparison, "^>")
})


test_that("Cheng test has its exact transformation invariances", {
  x <- ch3pp_fixture_x(13)
  y <- ch3pp_fixture_y(16)
  base <- cheng_sscm_equality_test(x, y, tol = 1e-10)
  swapped <- cheng_sscm_equality_test(y, x, tol = 1e-10)
  shifted.scaled <- cheng_sscm_equality_test(
    7 * sweep(x, 2L, c(10, -8, 4, 3, -2), "+"),
    0.3 * sweep(y, 2L, c(-5, 9, 2, -7, 1), "+"),
    tol = 1e-10
  )
  q <- qr.Q(qr(matrix(c(
    1, 2, 0, 1, -1,
    0, 1, 1, -2, 1,
    2, 0, 1, 1, 0,
    1, -1, 2, 0, 1,
    0, 2, -1, 1, 2
  ), 5, 5)))
  rotated <- cheng_sscm_equality_test(x %*% q, y %*% q, tol = 1e-10)

  for (candidate in list(swapped, shifted.scaled, rotated)) {
    expect_equal(candidate$components$T.SSCM, base$components$T.SSCM,
                 tolerance = 2e-10)
    expect_equal(candidate$components$delta.hat, base$components$delta.hat,
                 tolerance = 2e-10)
    expect_equal(candidate$components$null.variance.hat,
                 base$components$null.variance.hat, tolerance = 2e-10)
    expect_equal(unname(candidate$statistic), unname(base$statistic),
                 tolerance = 2e-10)
    expect_equal(candidate$p.value, base$p.value, tolerance = 2e-11)
  }
})


test_that("Cheng zero-radius and nonconvergence contracts never fake calibration", {
  x.zero <- rbind(
    c(-1, 0), c(0, -1), c(0, 0), c(0, 1), c(1, 0)
  )
  y <- cbind(sin(seq_len(8)), cos(0.7 * seq_len(8)))
  expect_error(
    cheng_sscm_equality_test(x.zero, y),
    "zero fitted radius"
  )
  uncalibrated <- NULL
  expect_warning(
    uncalibrated <- cheng_sscm_equality_test(x.zero, y, strict = FALSE),
    "No numerical repair"
  )
  expect_false(uncalibrated$diagnostics$calibrated)
  expect_true(is.na(uncalibrated$p.value))
  expect_true(is.na(unname(uncalibrated$statistic)))
  expect_true(is.finite(uncalibrated$components$T.SSCM))
  expect_error(
    cheng_sscm_equality_test(
      ch3pp_fixture_x(12), ch3pp_fixture_y(15), max_iter = 1L
    ),
    "did not converge"
  )
})


test_that("Feng spatial-rank statistic matches ordered four-index formulas", {
  x <- ch3pp_fixture_x(5)[, 1:3]
  y <- ch3pp_fixture_y(6)[, 1:3]
  result <- feng_spatial_rank_proportionality_test(x, y)
  reference <- ch3pp_rank_reference(x, y)

  expect_true(result$diagnostics$calibrated)
  expect_equal(result$components$ordered.sum.x, reference$sum.x,
               tolerance = 3e-12)
  expect_equal(result$components$ordered.sum.y, reference$sum.y,
               tolerance = 3e-12)
  expect_equal(result$components$ordered.cross.sum, reference$sum.cross,
               tolerance = 3e-12)
  expect_equal(result$components$A1.raw, reference$A1, tolerance = 3e-13)
  expect_equal(result$components$A2.raw, reference$A2, tolerance = 3e-13)
  expect_equal(result$components$C12.raw, reference$C12, tolerance = 3e-13)
  expect_equal(result$components$A1, 3 * reference$A1, tolerance = 3e-13)
  expect_equal(result$components$T.HT, reference$T, tolerance = 3e-12)
  expect_equal(result$components$sigma0.squared.hat, reference$variance,
               tolerance = 3e-13)
  expect_equal(result$components$sigma0.hat, reference$standard.error,
               tolerance = 3e-13)
  expect_equal(unname(result$statistic), reference$z, tolerance = 3e-12)
})


test_that("Feng feasible variance locks the pooled and p scaling factors", {
  x <- ch3pp_fixture_x(6)[, 1:4]
  y <- ch3pp_fixture_y(8)[, 1:4]
  result <- feng_spatial_rank_proportionality_test(x, y)
  cmp <- result$components
  expected <- 4 * (1 / 6 + 1 / 8)^2 / (6 + 8) / (4 + 2)^2 *
    4^2 * (6 * cmp$A1.raw + 8 * cmp$A2.raw)
  wrong.without.pool <- 4 * (1 / 6 + 1 / 8)^2 / (4 + 2)^2 *
    4^2 * (6 * cmp$A1.raw + 8 * cmp$A2.raw)

  expect_equal(cmp$sigma0.squared.hat, expected, tolerance = 1e-14)
  expect_false(isTRUE(all.equal(cmp$sigma0.squared.hat, wrong.without.pool)))
  expect_equal(cmp$T.HT, cmp$A1 + cmp$A2 - 2 * cmp$C12,
               tolerance = 1e-14)
  expect_match(result$diagnostics$primary.normalization, "trace one")
  expect_match(result$diagnostics$book.erratum, "omits required powers")
  expect_match(result$diagnostics$rejection.comparison, ">=")
})


test_that("Feng test is invariant to group order, location, scalar, and rotation", {
  x <- ch3pp_fixture_x(6)[, 1:4]
  y <- ch3pp_fixture_y(7)[, 1:4]
  base <- feng_spatial_rank_proportionality_test(x, y)
  swapped <- feng_spatial_rank_proportionality_test(y, x)
  transformed <- feng_spatial_rank_proportionality_test(
    9 * sweep(x, 2L, c(1e8, -2e8, 3e8, -4e8), "+"),
    0.2 * sweep(y, 2L, c(-3e8, 2e8, -1e8, 4e8), "+")
  )
  q <- qr.Q(qr(matrix(c(
    1, 2, 0, -1,
    0, 1, 2, 1,
    2, 0, 1, 1,
    1, -1, 0, 2
  ), 4, 4)))
  rotated <- feng_spatial_rank_proportionality_test(x %*% q, y %*% q)

  for (candidate in list(swapped, transformed, rotated)) {
    expect_equal(candidate$components$T.HT, base$components$T.HT,
                 tolerance = 2e-7)
    expect_equal(candidate$components$sigma0.squared.hat,
                 base$components$sigma0.squared.hat, tolerance = 3e-8)
    expect_equal(unname(candidate$statistic), unname(base$statistic),
                 tolerance = 2e-7)
    expect_equal(candidate$p.value, base$p.value, tolerance = 1e-6)
  }
})


test_that("Feng ties, minimum n, and zero variance have explicit contracts", {
  x <- rbind(
    c(0, 0), c(0, 0), c(1, 0), c(0, 1), c(1, 1)
  )
  y <- rbind(
    c(-1, 0), c(0, -1), c(1, 0), c(0, 1), c(2, 1), c(-1, 2)
  )
  tied <- feng_spatial_rank_proportionality_test(x, y)
  expect_equal(tied$components$zero.pairs[["x"]], 1)
  expect_true(is.finite(tied$components$T.HT))
  expect_error(
    feng_spatial_rank_proportionality_test(x[1:3, ], y),
    "at least 4"
  )
  identical.x <- matrix(1, 5, 2)
  identical.y <- matrix(2, 6, 2)
  expect_error(
    feng_spatial_rank_proportionality_test(identical.x, identical.y),
    "variance"
  )
  invalid <- NULL
  expect_warning(
    invalid <- feng_spatial_rank_proportionality_test(
      identical.x, identical.y, strict = FALSE
    ),
    "No variance floor"
  )
  expect_false(invalid$diagnostics$calibrated)
  expect_true(is.na(invalid$p.value))
  expect_true(is.na(unname(invalid$statistic)))
})


ch3pp_axis_fixture <- function() {
  rbind(c(1, 0), c(-1, 0), c(0, 1), c(0, -1))
}


ch3pp_lp_vertex_reference <- function(a, target, lambda) {
  boundary.matrix <- rbind(a[1, ], a[1, ], a[2, ], a[2, ],
                           c(1, 0), c(0, 1))
  boundary.value <- c(
    target[1] - lambda, target[1] + lambda,
    target[2] - lambda, target[2] + lambda, 0, 0
  )
  candidates <- list()
  for (i in seq_len(5)) for (j in seq.int(i + 1L, 6L)) {
    system <- boundary.matrix[c(i, j), , drop = FALSE]
    if (abs(det(system)) < 1e-12) next
    value <- solve(system, boundary.value[c(i, j)])
    if (max(abs(a %*% value - target)) <= lambda + 1e-10) {
      candidates[[length(candidates) + 1L]] <- value
    }
  }
  objective <- vapply(candidates, function(value) sum(abs(value)), numeric(1))
  candidates[[which.min(objective)]]
}


test_that("SCLIME internal solver recovers the analytic identity solution", {
  lambda <- 0.2
  fit <- cpp_ch3pp_sclime(diag(3), lambda, 1e-8, 200000L)
  expect_true(fit$all_converged)
  expect_equal(fit$raw_solution, diag(1 - lambda, 3), tolerance = 2e-8)
  expect_equal(fit$solution, diag(1 - lambda, 3), tolerance = 2e-8)
  expect_true(all(fit$primal_violation <= 1e-8))
  expect_true(all(fit$dual_violation <= 1e-8))
  expect_true(all(fit$stationarity_residual <= 1e-8))
  expect_true(all(fit$relative_gap <= 1e-8))
  expect_equal(fit$primal_objective, rep(1 - lambda, 3), tolerance = 2e-8)
  expect_equal(fit$dual_objective, fit$primal_objective, tolerance = 2e-8)
})


test_that("SCLIME agrees with an independent two-dimensional LP vertex search", {
  a <- matrix(c(1, 0.3, 0.3, 0.8), 2, 2)
  lambda <- 0.18
  fit <- cpp_ch3pp_sclime(a, lambda, 1e-8, 300000L)
  references <- cbind(
    ch3pp_lp_vertex_reference(a, c(1, 0), lambda),
    ch3pp_lp_vertex_reference(a, c(0, 1), lambda)
  )
  expect_true(fit$all_converged)
  expect_equal(fit$raw_solution, references, tolerance = 3e-8)
  expect_equal(colSums(abs(fit$raw_solution)),
               colSums(abs(references)), tolerance = 3e-8)
  expect_true(all(apply(
    a %*% fit$raw_solution - diag(2), 2L,
    function(value) max(abs(value)) <= lambda + 1e-8
  )))
  expect_true(all(fit$stationarity_residual <= 1e-8))
  expect_true(all(fit$relative_gap <= 1e-8))
})


test_that("SGLASSO penalizes the diagonal as the primary objective requires", {
  lambda <- 0.2
  fit <- cpp_ch3pp_sglasso(diag(3), lambda, 1e-9, 10000L, 1, 100L)
  expected <- diag(1 / (1 + lambda), 3)
  expect_true(fit$converged)
  expect_false(fit$backtracking_failed)
  expect_equal(fit$solution, expected, tolerance = 2e-9)
  expect_equal(fit$inverse, diag(1 + lambda, 3), tolerance = 3e-9)
  expect_lte(fit$kkt_residual, 1e-9)
  expect_gt(fit$minimum_eigenvalue, 0)
  expect_false(isTRUE(all.equal(fit$solution, diag(3))))
})


test_that("public precision APIs reproduce analytic SSCM identity fits", {
  x <- ch3pp_axis_fixture()
  sclime <- spatial_sign_precision(
    x, lambda = 0.2, method = "sclime", solver_tol = 1e-7
  )
  sglasso <- spatial_sign_precision(
    x, lambda = 0.2, method = "sglasso", solver_tol = 1e-8
  )

  expect_s3_class(sclime, "spatial_sign_precision_fit")
  expect_true(sclime$valid)
  expect_true(sglasso$valid)
  expect_equal(as.numeric(sclime$center), c(0, 0), tolerance = 1e-14)
  expect_equal(sclime$sscm, diag(0.5, 2), tolerance = 1e-14)
  expect_equal(sclime$coefficient, diag(2), tolerance = 1e-14)
  expect_equal(sclime$estimate, diag(0.8, 2), tolerance = 2e-7)
  expect_equal(sglasso$estimate, diag(1 / 1.2, 2), tolerance = 2e-8)
  expect_true(sclime$diagnostics$solver$all.columns.certified)
  expect_true(sglasso$diagnostics$solver$kkt.certified)
  expect_true(sglasso$diagnostics$diagonal.penalty)
  expect_match(sclime$diagnostics$target, "Lambda0")
})


test_that("SCLIME symmetrization is the smaller-absolute primary rule", {
  a <- matrix(c(1, 0.35, 0.1, 0.35, 1, 0.25, 0.1, 0.25, 1), 3, 3)
  fit <- cpp_ch3pp_sclime(a, 0.25, 2e-7, 300000L)
  expect_true(fit$all_converged)
  raw <- fit$raw_solution
  expected <- raw
  for (i in seq_len(3)) for (j in seq_len(3)) {
    expected[i, j] <- if (abs(raw[i, j]) <= abs(raw[j, i]))
      raw[i, j] else raw[j, i]
  }
  expect_equal(fit$solution, expected, tolerance = 1e-14)
  expect_equal(fit$solution, t(fit$solution), tolerance = 1e-14)
  expect_true(all(fit$primal_violation <= 2e-7))
  expect_true(all(fit$stationarity_residual <= 2e-7))
  expect_true(all(fit$relative_gap <= 2e-7))
  expect_true(is.finite(fit$final_feasibility_violation))
})


test_that("precision estimates have translation, global-scale, and permutation equivariance", {
  x <- ch3pp_axis_fixture()
  permutation <- c(2, 1)
  for (method in c("sclime", "sglasso")) {
    base <- spatial_sign_precision(
      x, 0.25, method = method, solver_tol = 2e-7
    )
    transformed <- spatial_sign_precision(
      11 * sweep(x, 2L, c(100, -70), "+"),
      0.25, method = method, solver_tol = 2e-7
    )
    permuted <- spatial_sign_precision(
      x[, permutation, drop = FALSE],
      0.25, method = method, solver_tol = 2e-7
    )
    expect_true(base$valid)
    expect_true(transformed$valid)
    expect_true(permuted$valid)
    expect_equal(transformed$sscm, base$sscm, tolerance = 3e-12)
    expect_equal(transformed$estimate, base$estimate, tolerance = 3e-7)
    expect_equal(permuted$estimate,
                 base$estimate[permutation, permutation, drop = FALSE],
                 tolerance = 3e-7)
  }
})


test_that("failed optimization never returns a precision-looking matrix", {
  x <- ch3pp_axis_fixture()
  expect_error(
    spatial_sign_precision(
      x, 0.2, method = "sclime", solver_tol = 1e-12,
      solver_max_iter = 1L
    ),
    "certificate"
  )
  failed.sclime <- NULL
  expect_warning(
    failed.sclime <- spatial_sign_precision(
      x, 0.2, method = "sclime", solver_tol = 1e-12,
      solver_max_iter = 1L, strict = FALSE
    ),
    "No precision estimate"
  )
  expect_false(failed.sclime$valid)
  expect_null(failed.sclime$estimate)
  expect_identical(failed.sclime$diagnostics$failure.stage,
                   "optimization certificate")
  expect_true(is.list(failed.sclime$diagnostics$last.iterate))

  failed.sglasso <- NULL
  expect_warning(
    failed.sglasso <- spatial_sign_precision(
      x, 0.2, method = "sglasso", solver_tol = 1e-14,
      solver_max_iter = 1L, strict = FALSE
    ),
    "No precision estimate"
  )
  expect_false(failed.sglasso$valid)
  expect_null(failed.sglasso$estimate)
  expect_false(failed.sglasso$diagnostics$solver$kkt.certified)
})


test_that("degenerate SSCM and median failures are explicit and unrepaired", {
  identical.data <- matrix(1, 6, 3)
  expect_error(
    spatial_sign_precision(identical.data, 0.3, method = "sclime"),
    "no directional information"
  )
  invalid <- NULL
  expect_warning(
    invalid <- spatial_sign_precision(
      identical.data, 0.3, method = "sglasso", strict = FALSE
    ),
    "no directional information"
  )
  expect_false(invalid$valid)
  expect_null(invalid$estimate)
  expect_identical(invalid$diagnostics$failure.stage, "degenerate SSCM")
  expect_match(invalid$diagnostics$no.repair, "No ridge")

  x <- ch3pp_fixture_x(9)[, 1:3]
  failed.median <- NULL
  expect_warning(
    failed.median <- spatial_sign_precision(
      x, 0.4, method = "sglasso", median_max_iter = 1L,
      strict = FALSE
    ),
    "spatial median did not converge"
  )
  expect_false(failed.median$valid)
  expect_null(failed.median$estimate)
  expect_identical(failed.median$diagnostics$failure.stage,
                   "spatial median")
})


test_that("threshold support uses >=, includes diagonals, and rejects invalid fits", {
  fit <- structure(
    list(
      valid = TRUE,
      estimate = matrix(c(0.5, -0.2, -0.2, 0.1), 2, 2),
      method = "sclime"
    ),
    class = "spatial_sign_precision_fit"
  )
  thresholded <- threshold_spatial_sign_precision(fit, 0.2)
  expect_s3_class(thresholded, "thresholded_spatial_sign_precision")
  expect_equal(thresholded$estimate,
               matrix(c(0.5, -0.2, -0.2, 0), 2, 2))
  expect_true(thresholded$support[1, 1])
  expect_true(thresholded$support[1, 2])
  expect_false(thresholded$support[2, 2])
  expect_identical(thresholded$comparison, ">=")
  expect_true(thresholded$diagonal.thresholded)

  invalid <- fit
  invalid$valid <- FALSE
  invalid$estimate <- NULL
  expect_error(threshold_spatial_sign_precision(invalid, 0.2), "invalid")
  expect_error(threshold_spatial_sign_precision(fit, -1), "non-negative")
})


test_that("public validators reject malformed inputs and tuning controls", {
  x <- ch3pp_axis_fixture()
  expect_error(cheng_sscm_equality_test(x, x[, 1, drop = FALSE]),
               "same number")
  expect_error(feng_spatial_rank_proportionality_test(x[1:3, ], x),
               "at least 4")
  expect_error(spatial_sign_precision(x, 0), "positive")
  expect_error(spatial_sign_precision(x, NA_real_), "positive")
  expect_error(spatial_sign_precision(x, 0.2, solver_tol = 0),
               "solver_tol")
  expect_error(spatial_sign_precision(x, 0.2, solver_max_iter = 1.5),
               "positive integer")
  expect_error(spatial_sign_precision(x, 0.2, initial_step = 0),
               "initial_step")
  expect_error(spatial_sign_precision(x, 0.2, strict = NA), "strict")
})


test_that("SCLIME rejects an infeasible symmetrized precision estimate", {
  coefficient <- structure(
    c(
      1, 0.648397841714642, -0.313576247774984,
      0.648397841714642, 1, -0.529847057204945,
      -0.313576247774984, -0.529847057204945, 1
    ),
    dim = c(3L, 3L)
  )
  eig <- eigen(coefficient, symmetric = TRUE)
  signs <- eig$vectors %*% diag(sqrt(eig$values), 3L) %*%
    t(eig$vectors)
  x <- rbind(signs, -signs)

  native <- cpp_ch3pp_sclime(
    coefficient, 0.05, 1e-5, 50000L
  )
  expect_true(native$all_converged)
  expect_lte(max(native$primal_violation), 1e-5)
  expect_gt(native$final_feasibility_violation, 1e-5)

  expect_error(
    spatial_sign_precision(
      x, lambda = 0.05, method = "sclime", solver_tol = 1e-5,
      solver_max_iter = 50000L
    ),
    "certificate"
  )
  expect_warning(
    invalid <- spatial_sign_precision(
      x, lambda = 0.05, method = "sclime", solver_tol = 1e-5,
      solver_max_iter = 50000L, strict = FALSE
    ),
    "certificate"
  )
  expect_false(invalid$valid)
  expect_null(invalid$estimate)
  expect_gt(
    invalid$diagnostics$solver$symmetrized.feasibility.violation,
    1e-5
  )
})
