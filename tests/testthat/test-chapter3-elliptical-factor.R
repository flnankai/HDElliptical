ch3ef_fixture <- function() {
  z <- rbind(
    c(3.0, 2.0, 0.5, -0.2),
    c(2.0, 1.3, -0.4, 0.8),
    c(1.2, 0.9, 1.1, -0.7),
    c(0.6, -0.3, 0.7, 1.4),
    c(1.6, 1.1, -1.2, -0.5),
    c(0.8, 0.2, 1.5, 0.4)
  )
  x <- rbind(z, -z)
  colnames(x) <- paste0("v", seq_len(ncol(x)))
  x
}


ch3ef_sign_reference <- function(x, center) {
  residual <- sweep(x, 2L, center, "-")
  radius <- sqrt(rowSums(residual^2))
  residual / radius
}


test_that("POET-SS pilot and reconstruction follow the literal formulas", {
  x <- ch3ef_fixture()
  fit <- spatial_sign_poet(
    x, factors = 1, threshold = 0, tol = 1e-10, max_iter = 2000L
  )
  u <- ch3ef_sign_reference(x, fit$center)
  pilot <- ncol(x) * crossprod(u) / nrow(x)

  expect_s3_class(fit, "elliptical_factor_fit")
  expect_true(fit$valid)
  expect_equal(fit$signs, u, tolerance = 2e-10, ignore_attr = TRUE)
  expect_equal(fit$pilot, pilot, tolerance = 2e-10)
  expect_equal(sum(diag(fit$pilot)), ncol(x), tolerance = 2e-10)
  expect_equal(fit$estimate, fit$pilot, tolerance = 2e-10)
  expect_equal(
    fit$components$low.rank + fit$components$idiosyncratic.raw,
    fit$pilot, tolerance = 2e-12
  )
  expect_equal(
    fit$components$low.rank +
      fit$components$idiosyncratic.thresholded,
    fit$estimate, tolerance = 2e-12
  )
  expect_equal(fit$diagnostics$positive.definiteness.repair, "none")
  expect_equal(colnames(fit$estimate), colnames(x))
  expect_equal(names(fit$center), colnames(x))
})


test_that("common threshold and generalized maps affect only residual off diagonal", {
  x <- ch3ef_fixture()
  raw <- spatial_sign_poet(x, factors = 1, threshold = 0)
  hard <- spatial_sign_poet(x, factors = 1, threshold = 1e6, rule = "hard")
  soft <- spatial_sign_poet(x, factors = 1, threshold = 0.05, rule = "soft")
  residual <- raw$components$idiosyncratic.raw

  expect_equal(
    diag(hard$components$idiosyncratic.thresholded), diag(residual),
    tolerance = 1e-13
  )
  expect_equal(
    hard$components$idiosyncratic.thresholded[row(hard$estimate) !=
                                                col(hard$estimate)],
    rep(0, ncol(x) * (ncol(x) - 1L))
  )
  expected.soft <- residual
  off <- row(expected.soft) != col(expected.soft)
  expected.soft[off] <- sign(residual[off]) *
    pmax(abs(residual[off]) - 0.05, 0)
  expect_equal(
    soft$components$idiosyncratic.thresholded, expected.soft,
    tolerance = 1e-13
  )

  rate.fit <- spatial_sign_poet(x, factors = 1, constant = 2)
  expected.rate <- sqrt(log(ncol(x)) / nrow(x)) +
    sqrt(log(nrow(x)) / nrow(x))
  expect_equal(unname(rate.fit$threshold[1, 2]), 2 * expected.rate)
  expect_equal(rate.fit$diagnostics$threshold.source,
               "constant * {sqrt(log(p)/n) + sqrt(log(n)/n)}")
  expect_equal(rate.fit$diagnostics$tuning.selection, "none")
})


test_that("ER and GR factor counts reproduce their primary criteria", {
  x <- ch3ef_fixture()
  er <- elliptical_factor_number(x, max_factors = 3, method = "er")
  expected.er <- er$eigenvalues[1:3] / er$eigenvalues[2:4]
  expect_equal(er$criterion, expected.er, tolerance = 1e-12)
  expect_equal(er$selected, which.max(expected.er))
  expect_equal(er$max.factors, 3L)

  gr <- elliptical_factor_number(x, max_factors = 2, method = "gr")
  r <- min(nrow(x), ncol(x))
  v <- vapply(0:2, function(j) sum(gr$eigenvalues[(j + 1):(r - 1)]),
              numeric(1))
  expected.gr <- vapply(1:2, function(j) {
    log1p(gr$eigenvalues[j] / v[j]) /
      log1p(gr$eigenvalues[j + 1] / v[j + 1])
  }, numeric(1))
  expect_equal(gr$tail.sums, v, tolerance = 1e-13)
  expect_equal(gr$criterion, expected.gr, tolerance = 1e-12)
  expect_equal(gr$selected, which.max(expected.gr))
  expect_error(
    elliptical_factor_number(x[, 1:2], 1, method = "gr"),
    "min\\(n, p\\)"
  )
})


test_that("the one-step Tyler pilot matches an independent matrix formula", {
  x <- ch3ef_fixture()
  precision <- matrix(c(
    1.5, 0.1, 0, 0,
    0.1, 1.2, 0.05, 0,
    0, 0.05, 0.9, 0.08,
    0, 0, 0.08, 1.1
  ), 4, 4, byrow = TRUE)
  fit <- poet_tme(
    x, factors = 1, threshold = 0,
    preliminary_precision = precision,
    tol = 1e-10, max_iter = 2000L
  )
  residual <- sweep(x, 2L, fit$center, "-")
  direct <- Reduce(`+`, lapply(seq_len(nrow(x)), function(i) {
    z <- residual[i, ]
    tcrossprod(z) / drop(crossprod(z, precision %*% z))
  })) * ncol(x) / nrow(x)

  expect_equal(fit$pilot, direct, tolerance = 2e-12, ignore_attr = TRUE)
  expect_equal(fit$estimate, fit$pilot, tolerance = 2e-12)
  expect_equal(fit$diagnostics$refinement.steps, 1L)
  expect_equal(fit$diagnostics$preliminary.precision.source,
               "supplied matrix")
  expect_gt(fit$diagnostics$tyler$minimum.standardized.quadratic, 0)
  expect_equal(fit$diagnostics$positive.definiteness.repair, "none")
})


test_that("sign and Tyler factor pilots obey their stated invariances", {
  x <- ch3ef_fixture()
  shift <- c(1e6, -2e6, 3e6, -4e6)
  sign0 <- spatial_sign_poet(x, factors = 1, threshold = 0)
  sign1 <- spatial_sign_poet(
    sweep(x, 2L, shift, "+") * 1e120,
    factors = 1, threshold = 0, tol = 1e-8, max_iter = 3000L
  )
  expect_equal(sign1$pilot, sign0$pilot, tolerance = 2e-8)

  tyler0 <- poet_tme(
    x, factors = 1, threshold = 0, preliminary_precision = diag(4)
  )
  tyler1 <- poet_tme(
    sweep(x, 2L, shift, "+") * 1e120,
    factors = 1, threshold = 0, preliminary_precision = diag(4),
    max_iter = 3000L
  )
  expect_equal(tyler1$pilot, tyler0$pilot, tolerance = 2e-8)

  permutation <- c(3, 1, 4, 2)
  permuted <- spatial_sign_poet(
    x[, permutation], factors = 1, threshold = 0
  )
  expect_equal(
    permuted$pilot,
    sign0$pilot[permutation, permutation], tolerance = 2e-10,
    ignore_attr = TRUE
  )
})


test_that("Woodbury factor precision uses raw residual and certified solver", {
  x <- ch3ef_fixture()
  factor.fit <- spatial_sign_poet(x, factors = 1, threshold = 0.15)
  precision.fit <- elliptical_factor_precision(
    factor.fit, lambda = 0.5, method = "sglasso",
    solver_tol = 1e-7, solver_max_iter = 20000L
  )
  expect_s3_class(precision.fit, "elliptical_factor_precision_fit")
  expect_true(precision.fit$valid)
  vu <- precision.fit$idiosyncratic.precision
  gamma <- factor.fit$components$leading.eigenvectors
  lambda <- factor.fit$components$leading.eigenvalues
  middle <- diag(1 / lambda, length(lambda)) +
    crossprod(gamma, vu %*% gamma)
  direct <- vu - vu %*% gamma %*% solve(middle) %*%
    t(gamma) %*% vu
  expect_equal(precision.fit$estimate, direct, tolerance = 2e-9)
  expect_equal(
    precision.fit$diagnostics$residual.source,
    "raw unthresholded idiosyncratic complement"
  )
  expect_true(precision.fit$diagnostics$solver$kkt.certified)
  expect_lte(precision.fit$diagnostics$solver$kkt.residual, 1e-7)
  expect_true(all(is.finite(precision.fit$estimate)))
  expect_equal(rownames(precision.fit$estimate), colnames(x))
})


test_that("failed precision certificates never return a plausible matrix", {
  x <- ch3ef_fixture()
  factor.fit <- spatial_sign_poet(x, factors = 1, threshold = 0.1)
  expect_warning(
    failed <- elliptical_factor_precision(
      factor.fit, lambda = 1e-8, method = "sclime",
      solver_tol = 1e-14, solver_max_iter = 1L, strict = FALSE
    ),
    "failed"
  )
  expect_false(failed$valid)
  expect_null(failed$estimate)
  expect_null(failed$idiosyncratic.precision)
  expect_match(failed$diagnostics$no.repair, "No ridge")
  expect_error(
    elliptical_factor_precision(
      factor.fit, lambda = 1e-8, method = "sclime",
      solver_tol = 1e-14, solver_max_iter = 1L
    ),
    "failed"
  )
})


test_that("factor APIs reject undefined inputs without numerical repair", {
  x <- ch3ef_fixture()
  expect_error(spatial_sign_poet(x, factors = 4, threshold = 0.1),
               "factors")
  expect_error(spatial_sign_poet(x, factors = 1, threshold = -1),
               "threshold")
  expect_error(spatial_sign_poet(x, factors = 1, rule = "scad", scad_a = 2),
               "scad_a")
  expect_error(
    poet_tme(
      x, factors = 1, threshold = 0.1,
      preliminary_precision = diag(c(1, 1, 1, -1))
    ),
    "positive definite"
  )
  expect_error(
    poet_tme(
      rbind(x, rep(0, ncol(x))), factors = 1, threshold = 0.1,
      preliminary_precision = diag(ncol(x)), zero_tol = 0
    ),
    "residual"
  )
  expect_error(elliptical_factor_number(x, 0), "max_factors")
  expect_error(elliptical_factor_number(x, 4), "max_factors")
  expect_error(elliptical_factor_precision(list(), 0.1), "valid")
})
