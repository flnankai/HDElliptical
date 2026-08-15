sscca_fixture <- function() {
  index <- seq_len(18)
  x <- cbind(
    x1 = sin(index / 2),
    x2 = cos(index / 3) + 0.08 * sin(index),
    x3 = sin(index / 5) - 0.12 * cos(index / 2)
  )
  y <- cbind(
    y1 = 0.85 * x[, 1] + 0.18 * cos(index / 4),
    y2 = -0.65 * x[, 2] + 0.22 * sin(index / 3)
  )
  rownames(x) <- rownames(y) <- paste0("r", index)
  list(x = x, y = y)
}


sscca_direct_blocks <- function(x, y, center) {
  signs <- spatial_sign(cbind(x, y), center = center)
  covariance <- crossprod(signs) / nrow(x)
  px <- ncol(x)
  py <- ncol(y)
  iy <- px + seq_len(py)
  list(
    signs = signs,
    xx = covariance[seq_len(px), seq_len(px), drop = FALSE],
    xy = covariance[seq_len(px), iy, drop = FALSE],
    yy = covariance[iy, iy, drop = FALSE]
  )
}


test_that("primary SSCCA uses the literal p-scaled spatial-sign blocks", {
  dat <- sscca_fixture()
  fit <- sscca(dat$x, dat$y, lambda_x = 0.01, lambda_y = 0.01)
  direct <- sscca_direct_blocks(dat$x, dat$y, fit$center$joint)
  dimension <- ncol(dat$x) + ncol(dat$y)

  expect_s3_class(fit, "sscca_fit")
  expect_true(fit$valid)
  expect_equal(unname(fit$operator$xx), unname(dimension * direct$xx),
               tolerance = 1e-13)
  expect_equal(unname(fit$operator$xy), unname(dimension * direct$xy),
               tolerance = 1e-13)
  expect_equal(unname(fit$operator$yy), unname(dimension * direct$yy),
               tolerance = 1e-13)
  expect_equal(
    drop(crossprod(fit$x.coefficients,
                   fit$operator$xx %*% fit$x.coefficients)),
    1, tolerance = 2e-11
  )
  expect_equal(
    drop(crossprod(fit$y.coefficients,
                   fit$operator$yy %*% fit$y.coefficients)),
    1, tolerance = 2e-11
  )
  expect_equal(
    unname(fit$canonical.correlations),
    drop(crossprod(fit$x.coefficients,
                   fit$operator$xy %*% fit$y.coefficients)),
    tolerance = 2e-12
  )
  expect_equal(fit$x.scores,
               drop(sweep(dat$x, 2, fit$center$x) %*% fit$x.coefficients),
               tolerance = 1e-13)
  expect_lte(fit$diagnostics$pair.kkt.x,
             fit$diagnostics$certificate.limit)
  expect_lte(fit$diagnostics$pair.kkt.y,
             fit$diagnostics$certificate.limit)
})


test_that("metric coordinate updates divide by non-unit diagonals", {
  ax <- diag(c(4, 1))
  ay <- matrix(1, 1, 1)
  cross <- matrix(c(0.5, 0.3), 2, 1)
  native <- cpp_ch6_sscca_primary(
    ax, ay, cross, 0.1, 0, 0L, 20L,
    1e-12, 100L, 1e-13, 1000L, 0
  )
  expected.raw <- c((0.5 - 0.1) / 4, (0.3 - 0.1) / 1)

  expect_true(native$valid)
  expect_true(native$converged)
  expect_equal(abs(as.numeric(native$raw_x)), expected.raw, tolerance = 1e-12)
  expect_false(isTRUE(all.equal(abs(as.numeric(native$raw_x)), c(0.4, 0.2))))
  expect_lt(native$pair_kkt_x, 1e-11)
  expect_equal(native$metric_norm_x, 1, tolerance = 1e-13)
})


test_that("both primary BIC formulas use the unnormalised lasso solution", {
  ax <- diag(c(4, 1))
  ay <- matrix(1, 1, 1)
  cross <- matrix(c(0.5, 0.3), 2, 1)
  lambda <- c(0.2, 0.1, 0)
  n <- 20L
  fit1 <- cpp_ch6_sscca_primary(
    ax, ay, cross, lambda, 0, 1L, n,
    1e-12, 100L, 1e-13, 1000L, 0
  )
  fit2 <- cpp_ch6_sscca_primary(
    ax, ay, cross, lambda, 0, 2L, n,
    1e-12, 100L, 1e-13, 1000L, 0
  )
  df1 <- vapply(lambda, function(value) {
    sum(abs(c((0.5 - value) / 4, max(0, 0.3 - value))) > 0)
  }, integer(1))

  expect_equal(
    fit1$bic_values_x,
    fit1$rss_values_x + df1 * log(n) / n,
    tolerance = 1e-12
  )
  expect_equal(
    fit2$bic_values_x,
    log(n * fit2$rss_values_x / (n - df1)) + df1 * log(n) / n,
    tolerance = 1e-12
  )
  expect_equal(fit1$selected_index_x,
               which.min(fit1$bic_values_x))
  expect_equal(fit2$selected_index_x,
               which.min(fit2$bic_values_x))
})


test_that("R BIC paths expose every supplied candidate and certificate", {
  dat <- sscca_fixture()
  fit <- sscca(
    dat$x, dat$y,
    lambda_x = c(0.04, 0.02, 0.01),
    lambda_y = c(0.03, 0.015, 0.005),
    selection = "bic1"
  )

  expect_equal(nrow(fit$tuning$candidate.x), 3L)
  expect_equal(nrow(fit$tuning$candidate.y), 3L)
  expect_equal(fit$tuning$candidate.x$lambda, c(0.04, 0.02, 0.01))
  expect_true(all(fit$tuning$candidate.x$converged))
  expect_true(fit$tuning$lambda.x %in% fit$tuning$candidate.x$lambda)
  expect_true(fit$tuning$lambda.y %in% fit$tuning$candidate.y$lambda)
})


test_that("primary SSCCA respects its exact paired-data invariances", {
  dat <- sscca_fixture()
  base <- sscca(dat$x, dat$y, 0.01, 0.01)
  shifted <- sscca(
    sweep(dat$x, 2, c(1e8, -2e8, 3e8), FUN = "+"),
    sweep(dat$y, 2, c(-4e8, 5e8), FUN = "+"),
    0.01, 0.01
  )
  scaled <- sscca(dat$x * 1e100, dat$y * 1e100, 0.01, 0.01)
  ordering <- c(18:10, 1:9)
  permuted <- sscca(dat$x[ordering, ], dat$y[ordering, ], 0.01, 0.01)
  signs.x <- diag(c(-1, 1, -1))
  signs.y <- diag(c(-1, 1))
  reflected <- sscca(dat$x %*% signs.x, dat$y %*% signs.y, 0.01, 0.01)

  expect_equal(shifted$canonical.correlations,
               base$canonical.correlations, tolerance = 2e-7)
  expect_equal(scaled$canonical.correlations,
               base$canonical.correlations, tolerance = 2e-10)
  expect_equal(permuted$canonical.correlations,
               base$canonical.correlations, tolerance = 2e-10)
  expect_equal(reflected$canonical.correlations,
               base$canonical.correlations, tolerance = 2e-10)
  expect_equal(unname(abs(drop(signs.x %*% reflected$x.coefficients))),
               unname(abs(base$x.coefficients)), tolerance = 2e-8)
  expect_equal(unname(abs(drop(signs.y %*% reflected$y.coefficients))),
               unname(abs(base$y.coefficients)), tolerance = 2e-8)
})


test_that("book PMD reduces to the leading SVD pair when l1 bounds are inactive", {
  cross <- matrix(c(0.8, 0.1, -0.2, 0.5, 0.3, 0.4), 3, 2)
  native <- cpp_ch6_sign_whitened_pmd(
    cross, sqrt(3), sqrt(2), 1e-11, 500L
  )
  reference <- svd(cross, nu = 1, nv = 1)

  expect_true(native$valid)
  expect_true(native$converged)
  expect_equal(native$association, reference$d[1], tolerance = 1e-10)
  expect_equal(abs(as.numeric(native$x)), abs(reference$u[, 1]),
               tolerance = 1e-9)
  expect_equal(abs(as.numeric(native$y)), abs(reference$v[, 1]),
               tolerance = 1e-9)
  expect_equal(native$l2_x, 1, tolerance = 1e-12)
  expect_equal(native$l2_y, 1, tolerance = 1e-12)
})


test_that("book PMD enforces active l1 bounds and handles a tied maximum", {
  cross <- matrix(c(1, 0.4, 0.2, 0.3, 0.9, -0.1), 3, 2)
  native <- cpp_ch6_sign_whitened_pmd(cross, 1.1, 1.05, 1e-10, 1000L)
  tie <- cpp_ch6_sign_whitened_pmd(matrix(1, 2, 2), 1.1, 1.1,
                                   1e-10, 1000L)

  expect_true(native$valid)
  expect_lte(native$l1_x, 1.1 + 2e-9)
  expect_lte(native$l1_y, 1.05 + 2e-9)
  expect_lte(native$l2_x, 1 + 2e-12)
  expect_lte(native$l2_y, 1 + 2e-12)
  expect_lte(native$constraint_violation_x, 2e-9)
  expect_true(tie$valid)
  expect_true(tie$top_tie_slack_x || tie$top_tie_slack_y)
  expect_lte(tie$l1_x, 1.1 + 2e-9)
})


test_that("book wrapper matches literal whitening and returns both geometries", {
  dat <- sscca_fixture()
  fit <- sign_whitened_sparse_cca(
    dat$x, dat$y, c_x = sqrt(3), c_y = sqrt(2)
  )
  direct <- sscca_direct_blocks(dat$x, dat$y, fit$center$joint)
  rx <- eigen(direct$xx, symmetric = TRUE)
  ry <- eigen(direct$yy, symmetric = TRUE)
  ix <- rx$vectors %*% diag(1 / sqrt(rx$values)) %*% t(rx$vectors)
  iy <- ry$vectors %*% diag(1 / sqrt(ry$values)) %*% t(ry$vectors)
  operator <- ix %*% direct$xy %*% iy

  expect_s3_class(fit, "sign_whitened_sparse_cca_fit")
  expect_equal(fit$operator, operator, tolerance = 2e-12)
  expect_equal(
    unname(fit$x.coefficients),
    unname(drop(ix %*% fit$x.whitened.directions)), tolerance = 2e-12
  )
  expect_equal(
    unname(fit$y.coefficients),
    unname(drop(iy %*% fit$y.whitened.directions)), tolerance = 2e-12
  )
  expect_equal(
    unname(fit$canonical.correlations),
    drop(crossprod(fit$x.whitened.directions,
                   fit$operator %*% fit$y.whitened.directions)),
    tolerance = 2e-12
  )
  expect_match(fit$diagnostics$variant, "not the penalized metric")
})


test_that("explicit ridge is reported and no hidden whitening repair occurs", {
  index <- seq_len(10)
  x <- cbind(a = index, b = index)
  y <- cbind(c = sin(index), d = cos(index))

  expect_error(
    sign_whitened_sparse_cca(x, y, sqrt(2), sqrt(2)),
    "not full-rank SPD"
  )
  fit <- sign_whitened_sparse_cca(
    x, y, sqrt(2), sqrt(2), ridge_x = 0.01
  )
  expect_true(fit$valid)
  expect_equal(fit$tuning$ridge.x, 0.01)
  expect_equal(
    fit$diagnostics$reciprocal.condition.x,
    min(eigen(fit$spatial.sign.covariance$xx + diag(0.01, 2),
              symmetric = TRUE, only.values = TRUE)$values) /
      max(eigen(fit$spatial.sign.covariance$xx + diag(0.01, 2),
                symmetric = TRUE, only.values = TRUE)$values),
    tolerance = 1e-12
  )
})


test_that("zero signs, solver failure, and invalid controls are explicit", {
  x <- matrix(0, 6, 2)
  y <- matrix(0, 6, 2)
  expect_error(sscca(x, y, 0.01, 0.01), "diagonal entry")
  expect_warning(
    invalid <- sscca(x, y, 0.01, 0.01, strict = FALSE),
    "diagonal entry"
  )
  expect_false(invalid$valid)
  expect_null(invalid$estimate)

  dat <- sscca_fixture()
  expect_error(sscca(dat$x, dat$y, c(0.1, 0.01), 0.01), "scalar")
  expect_error(sscca(dat$x, dat$y, -0.1, 0.01), "non-negative")
  expect_error(sscca(dat$x[-1, ], dat$y, 0.01, 0.01), "same number")
  expect_error(
    sign_whitened_sparse_cca(dat$x, dat$y, 0.9, sqrt(2)),
    "must lie"
  )
  expect_error(
    sign_whitened_sparse_cca(dat$x, dat$y, sqrt(3), sqrt(2),
                             ridge_x = -1),
    "non-negative"
  )
})


test_that("names, row pairing, supports, and score orientation are stable", {
  dat <- sscca_fixture()
  fit <- sscca(dat$x, dat$y, 0.02, 0.02)

  expect_identical(names(fit$x.coefficients), colnames(dat$x))
  expect_identical(names(fit$y.coefficients), colnames(dat$y))
  expect_identical(names(fit$x.scores), rownames(dat$x))
  expect_identical(fit$variable.names$x, colnames(dat$x))
  expect_gte(fit$x.coefficients[fit$diagnostics$sign.anchor.x], 0)
  expect_equal(fit$diagnostics$degrees.freedom.x,
               sum(abs(fit$diagnostics$raw.coefficients$x) > 0))

  changed <- dat$y
  rownames(changed) <- rev(rownames(changed))
  expect_error(sscca(dat$x, changed, 0.02, 0.02), "row names")
})
