ch6_cca_fixture <- function() {
  h1 <- sqrt(5 / 2) * c(1, -1, 0, 0, 0)
  h2 <- sqrt(5 / 6) * c(1, 1, -2, 0, 0)
  h3 <- sqrt(5 / 12) * c(1, 1, 1, -3, 0)
  h4 <- sqrt(5 / 20) * c(1, 1, 1, 1, -4)
  x <- cbind(x1 = h1, x2 = h2)
  y <- cbind(
    y1 = 0.8 * h1 + 0.6 * h3,
    y2 = 0.3 * h2 + sqrt(0.91) * h4
  )
  rownames(x) <- rownames(y) <- paste0("r", seq_len(5))
  list(x = x, y = y, h = cbind(h1, h2, h3, h4))
}


ch6_factor_fixture <- function() {
  x <- rbind(
    c(-3, -1, 0), c(-2, 1, 1), c(-1, -2, 0),
    c(1, 2, 0), c(2, -1, -1), c(3, 1, 0)
  )
  colnames(x) <- c("a", "b", "c")
  rownames(x) <- paste0("t", seq_len(nrow(x)))
  x
}


ch6_axis_fixture <- function() {
  x <- rbind(
    c(2, 0), c(-2, 0), c(1, 0), c(-1, 0),
    c(0, 1), c(0, -1)
  )
  colnames(x) <- c("first", "second")
  x
}


test_that("Chapter 6 public formals are explicit and stable", {
  expect_identical(
    names(formals(classical_pca)),
    c("x", "components", "center", "covariance_divisor", "eigen_tol")
  )
  expect_identical(
    names(formals(classical_cca)),
    c("x", "y", "components", "center", "covariance_divisor", "rank_tol")
  )
  expect_identical(names(formals(cca_bartlett_test)),
                   c("object", "null_rank"))
  expect_identical(
    names(formals(robust_factor_subspace)),
    c("x", "factors", "method", "center", "tol", "max_iter",
      "zero_tol", "zero_action", "eigen_tol")
  )
  expect_identical(
    names(formals(rts_factor)),
    c("x", "factors", "center", "zero_tol", "zero_pair_action",
      "eigen_tol")
  )
  expect_identical(
    names(formals(kendall_factor_number)),
    c("x", "kmax", "c", "method", "zero_tol", "zero_pair_action",
      "eigen_tol")
  )
})


test_that("classical PCA matches the book's divisor-n hand calculation", {
  x <- rbind(c(-2, 0), c(0, -1), c(0, 1), c(2, 0))
  colnames(x) <- c("wide", "narrow")
  rownames(x) <- paste0("o", seq_len(nrow(x)))
  fit <- classical_pca(
    x, components = 2, center = "mean", covariance_divisor = "n"
  )

  expect_s3_class(fit, "classical_pca_fit")
  expect_s3_class(fit, "hd_pca_fit")
  expect_true(isTRUE(fit$valid))
  expect_true(all(c(
    "method", "eigenvalues", "loadings", "scores", "center",
    "operator", "rank", "n", "p", "variable.names", "diagnostics", "call"
  ) %in% names(fit)))
  expect_equal(fit$operator, diag(c(2, 0.5)), tolerance = 1e-14,
               ignore_attr = TRUE)
  expect_equal(unname(fit$eigenvalues), c(2, 0.5), tolerance = 1e-14)
  expect_equal(unname(fit$loadings), diag(2), tolerance = 1e-14)
  expect_equal(unname(fit$scores), unname(x), tolerance = 1e-14)
  expect_equal(fit$reconstruction, x, tolerance = 1e-14)
  expect_identical(rownames(fit$loadings), colnames(x))
  expect_identical(colnames(fit$loadings), c("PC1", "PC2"))
  expect_identical(colnames(fit$scores), c("PC1", "PC2"))
  expect_identical(rownames(fit$scores), rownames(x))
  expect_identical(fit$diagnostics$covariance.divisor, "n")
  expect_identical(fit$diagnostics$sign.anchor.coordinates, c(1L, 2L))
})


test_that("PCA divisor changes eigenvalue scale but not directions or scores", {
  x <- ch6_factor_fixture()
  fit.n <- classical_pca(x, components = 2, covariance_divisor = "n")
  fit.n1 <- classical_pca(x, components = 2, covariance_divisor = "n-1")

  expect_equal(fit.n1$eigenvalues,
               fit.n$eigenvalues * nrow(x) / (nrow(x) - 1),
               tolerance = 1e-12)
  expect_equal(fit.n1$loadings, fit.n$loadings, tolerance = 1e-12)
  expect_equal(fit.n1$scores, fit.n$scores, tolerance = 1e-12)
  expect_equal(fit.n1$projector, fit.n$projector, tolerance = 1e-12)
  expect_equal(sum(fit.n$explained.variance),
               sum(fit.n1$explained.variance), tolerance = 1e-14)
  expect_identical(fit.n$diagnostics$covariance.divisor.value,
                   as.numeric(nrow(x)))
  expect_identical(fit.n1$diagnostics$covariance.divisor.value,
                   as.numeric(nrow(x) - 1))
})


test_that("PCA centering is explicit, named, and translation equivariant", {
  x <- ch6_factor_fixture()
  shift <- c(a = 10, b = -7, c = 3)
  shifted <- sweep(x, 2L, shift, "+")
  mean.fit <- classical_pca(x, components = 2, center = "mean")
  shifted.fit <- classical_pca(shifted, components = 2, center = "mean")

  expect_equal(shifted.fit$operator, mean.fit$operator, tolerance = 1e-12)
  expect_equal(shifted.fit$loadings, mean.fit$loadings, tolerance = 1e-12)
  expect_equal(shifted.fit$scores, mean.fit$scores, tolerance = 1e-12)
  expect_equal(shifted.fit$center, mean.fit$center + shift,
               tolerance = 1e-14)
  expect_equal(shifted.fit$reconstruction,
               sweep(mean.fit$reconstruction, 2L, shift, "+"),
               tolerance = 1e-12)

  supplied <- c(c = 30, a = 10, b = 20)
  named.fit <- classical_pca(x, components = 2, center = supplied)
  expect_identical(names(named.fit$center), colnames(x))
  expect_equal(named.fit$center, c(a = 10, b = 20, c = 30))
  expect_identical(named.fit$diagnostics$centering, "supplied")
  expect_error(classical_pca(x, center = c(wrong = 1, b = 2, c = 3)),
               "match")
})


test_that("PCA repeated eigenspaces use projectors and cannot be split", {
  x <- rbind(
    c(sqrt(2), 0), c(-sqrt(2), 0),
    c(0, sqrt(2)), c(0, -sqrt(2))
  )
  fit <- classical_pca(x, covariance_divisor = "n")

  expect_identical(fit$rank, 2L)
  expect_equal(unname(fit$eigenvalues), c(1, 1), tolerance = 1e-14)
  expect_equal(fit$projector, diag(2), tolerance = 1e-14,
               ignore_attr = TRUE)
  expect_length(fit$diagnostics$repeated.eigenvalue.blocks, 1L)
  expect_identical(fit$diagnostics$repeated.eigenvalue.blocks[[1]], 1:2)
  expect_error(
    classical_pca(x, components = 1, covariance_divisor = "n"),
    "repeated-eigenvalue"
  )
})


test_that("PCA fails strictly at invalid rank and data boundaries", {
  expect_error(classical_pca(matrix(1, 4, 2)), "rank zero")
  expect_error(classical_pca(matrix(1:8, 4, 2), components = 0),
               "components")
  expect_error(classical_pca(matrix(1:8, 4, 2), components = 3),
               "components")
  expect_error(classical_pca(matrix(c(1, NA, 3, 4), 2, 2)), "finite")
  bad.names <- matrix(1:8, 4, 2,
                      dimnames = list(NULL, c("same", "same")))
  expect_error(classical_pca(bad.names), "unique")
  expect_error(classical_pca(matrix(1:8, 4, 2), eigen_tol = 0),
               "strictly positive")
  expect_error(classical_pca(matrix(1:2, 1, 2)), "at least")
})


test_that("native PCA moment kernel equals direct centered cross-products", {
  x <- ch6_factor_fixture()
  center <- c(0.5, -0.25, 0.75)
  out <- HDElliptical:::cpp_ch6cf_center_moments(x, center, 6)
  centered <- sweep(x, 2L, center, "-")

  expect_equal(out$centered, centered, tolerance = 0, ignore_attr = TRUE)
  expect_equal(out$covariance, crossprod(centered) / 6,
               tolerance = 1e-14, ignore_attr = TRUE)
  expect_identical(out$divisor, 6)
  expect_error(
    HDElliptical:::cpp_ch6cf_center_moments(x, center[1:2], 6),
    "one value per column"
  )
  expect_error(HDElliptical:::cpp_ch6cf_center_moments(x, center, 0),
               "strictly positive")
})


test_that("classical CCA matches an exact orthogonal-basis calculation", {
  fixture <- ch6_cca_fixture()
  fit <- classical_cca(
    fixture$x, fixture$y, center = "mean", covariance_divisor = "n"
  )

  expect_s3_class(fit, "classical_cca_fit")
  expect_s3_class(fit, "hd_cca_fit")
  expect_true(isTRUE(fit$valid))
  expect_equal(unname(fit$all.canonical.correlations), c(0.8, 0.3),
               tolerance = 1e-12)
  expect_equal(fit$covariance$xx, diag(2), tolerance = 1e-12,
               ignore_attr = TRUE)
  expect_equal(fit$covariance$yy, diag(2), tolerance = 1e-12,
               ignore_attr = TRUE)
  expect_equal(fit$covariance$xy, diag(c(0.8, 0.3)), tolerance = 1e-12,
               ignore_attr = TRUE)
  expect_equal(crossprod(fit$x.scores) / 5, diag(2), tolerance = 1e-12,
               ignore_attr = TRUE)
  expect_equal(crossprod(fit$y.scores) / 5, diag(2), tolerance = 1e-12,
               ignore_attr = TRUE)
  expect_equal(crossprod(fit$x.scores, fit$y.scores) / 5,
               diag(c(0.8, 0.3)), tolerance = 1e-12,
               ignore_attr = TRUE)
  expect_identical(colnames(fit$x.coefficients), c("CC1", "CC2"))
  expect_identical(rownames(fit$x.coefficients), colnames(fixture$x))
  expect_identical(fit$diagnostics$joint.rank, 4L)
  expect_true(fit$diagnostics$joint.full.rank)
})


test_that("CCA common divisor and translations preserve canonical correlations", {
  fixture <- ch6_cca_fixture()
  fit.n <- classical_cca(fixture$x, fixture$y, covariance_divisor = "n")
  fit.n1 <- classical_cca(fixture$x, fixture$y,
                          covariance_divisor = "n-1")
  coefficient.scale <- sqrt((nrow(fixture$x) - 1) / nrow(fixture$x))

  expect_equal(fit.n1$canonical.correlations, fit.n$canonical.correlations,
               tolerance = 1e-12)
  expect_equal(fit.n1$x.coefficients,
               coefficient.scale * fit.n$x.coefficients,
               tolerance = 1e-12)
  expect_equal(fit.n1$y.coefficients,
               coefficient.scale * fit.n$y.coefficients,
               tolerance = 1e-12)
  expect_equal(fit.n1$projector$x.whitened, fit.n$projector$x.whitened,
               tolerance = 1e-12)

  shifted <- classical_cca(
    sweep(fixture$x, 2L, c(8, -3), "+"),
    sweep(fixture$y, 2L, c(-4, 12), "+"),
    covariance_divisor = "n"
  )
  expect_equal(shifted$canonical.correlations, fit.n$canonical.correlations,
               tolerance = 1e-12)
  expect_equal(shifted$x.scores, fit.n$x.scores, tolerance = 1e-12)
  expect_equal(shifted$y.scores, fit.n$y.scores, tolerance = 1e-12)
})


test_that("CCA is symmetric under view exchange and agrees with correlation", {
  fixture <- ch6_cca_fixture()
  xy <- classical_cca(fixture$x, fixture$y, covariance_divisor = "n")
  yx <- classical_cca(fixture$y, fixture$x, covariance_divisor = "n")

  expect_equal(yx$canonical.correlations, xy$canonical.correlations,
               tolerance = 1e-12)
  expect_equal(yx$projector$x.whitened, xy$projector$y.whitened,
               tolerance = 1e-12, ignore_attr = TRUE)
  expect_equal(yx$projector$y.whitened, xy$projector$x.whitened,
               tolerance = 1e-12, ignore_attr = TRUE)

  a <- c(-3, -1, 0, 2, 5)
  b <- c(2, 1, 0, -1, -4)
  one <- classical_cca(matrix(a), matrix(b), covariance_divisor = "n")
  expect_equal(unname(one$canonical.correlations), abs(cor(a, b)),
               tolerance = 1e-14)
})


test_that("CCA rejects rank defects, pairing errors, and split repeated roots", {
  fixture <- ch6_cca_fixture()
  duplicated <- cbind(fixture$x[, 1], fixture$x[, 1])
  expect_error(classical_cca(duplicated, fixture$y),
               "strictly positive definite")
  expect_error(classical_cca(fixture$x[-1, ], fixture$y),
               "same number")
  bad.rows <- fixture$y
  rownames(bad.rows) <- rev(rownames(bad.rows))
  expect_error(classical_cca(fixture$x, bad.rows), "row names")
  expect_error(classical_cca(fixture$x, fixture$y, rank_tol = 0),
               "strictly positive")

  h <- fixture$h
  tied.y <- cbind(
    0.5 * h[, 1] + sqrt(0.75) * h[, 3],
    0.5 * h[, 2] + sqrt(0.75) * h[, 4]
  )
  expect_error(
    classical_cca(fixture$x, tied.y, components = 1,
                  covariance_divisor = "n"),
    "repeated-eigenvalue"
  )
  tied <- classical_cca(fixture$x, tied.y, components = 2,
                        covariance_divisor = "n")
  expect_equal(tied$projector$x.whitened, diag(2), tolerance = 1e-12,
               ignore_attr = TRUE)
  expect_length(tied$diagnostics$repeated.correlation.blocks, 1L)
})


test_that("Bartlett tail tests match exact Wilks calculations", {
  fixture <- ch6_cca_fixture()
  fit <- classical_cca(fixture$x, fixture$y, covariance_divisor = "n")
  global <- cca_bartlett_test(fit, null_rank = 0)
  tail <- cca_bartlett_test(fit, null_rank = 1)

  expect_s3_class(global, "htest")
  expect_equal(global$wilks.lambda, 0.3276, tolerance = 1e-14)
  expect_equal(unname(global$statistic), 1.6739428901, tolerance = 1e-9)
  expect_identical(unname(global$parameter), 4L)
  expect_identical(global$sample.size.factor, 1.5)
  expect_equal(tail$wilks.lambda, 0.91, tolerance = 1e-14)
  expect_equal(unname(tail$statistic), 0.1414660192, tolerance = 1e-9)
  expect_identical(unname(tail$parameter), 1L)
  expect_equal(unname(tail$tested.canonical.correlations), 0.3,
               tolerance = 1e-12)
  expect_equal(global$p.value,
               pchisq(unname(global$statistic), 4, lower.tail = FALSE),
               tolerance = 0)
})


test_that("Bartlett contract is strict about centering, joint rank, and logs", {
  fixture <- ch6_cca_fixture()
  independent <- classical_cca(
    fixture$x, fixture$h[, 3:4], covariance_divisor = "n"
  )
  zero <- cca_bartlett_test(independent)
  expect_equal(unname(zero$statistic), 0, tolerance = 1e-14)
  expect_equal(zero$p.value, 1, tolerance = 1e-14)

  nonmean <- classical_cca(fixture$x, fixture$y, center = "none",
                           covariance_divisor = "n")
  expect_error(cca_bartlett_test(nonmean), "mean-centered")
  dependent <- classical_cca(fixture$x, fixture$x,
                             covariance_divisor = "n")
  expect_false(dependent$diagnostics$joint.full.rank)
  expect_error(cca_bartlett_test(dependent), "jointly")
  expect_error(cca_bartlett_test(list(valid = TRUE)), "valid")
  expect_error(cca_bartlett_test(independent, null_rank = 2), "null_rank")

  bad.log <- independent
  bad.log$all.canonical.correlations[1] <- 1
  expect_error(cca_bartlett_test(bad.log), "log\\(1 - rho\\^2\\)")
  bad.n <- independent
  bad.n$n <- 2
  expect_error(cca_bartlett_test(bad.n), "sample-size factor")
})


test_that("spatial-sign factor subspace equals its direct operator", {
  x <- ch6_axis_fixture()
  fit <- robust_factor_subspace(
    x, factors = 1, method = "spatial_sign", center = "none"
  )

  expect_s3_class(fit, "robust_factor_subspace_fit")
  expect_s3_class(fit, "hd_factor_fit")
  expect_true(isTRUE(fit$valid))
  expect_equal(fit$operator, diag(c(4 / 6, 2 / 6)), tolerance = 1e-14,
               ignore_attr = TRUE)
  expect_equal(unname(fit$principal_subspace), matrix(c(1, 0), 2, 1),
               tolerance = 1e-14)
  expect_equal(fit$projector, diag(c(1, 0)), tolerance = 1e-14,
               ignore_attr = TRUE)
  expect_equal(unname(fit$scores), matrix(x[, 1], ncol = 1),
               tolerance = 1e-14)
  expect_identical(fit$factor.scores, NULL)
  expect_match(fit$diagnostics$target, "no finite-sample equality")
  expect_identical(fit$diagnostics$zero.directions$normalization,
                   "sum of sign outer products divided by n")
})


test_that("robust factor operators obey their center and invariance contracts", {
  x <- ch6_factor_fixture()
  kendall <- robust_factor_subspace(
    x, factors = 1, method = "kendall", center = "mean"
  )
  shifted <- robust_factor_subspace(
    sweep(x, 2L, c(5, -9, 2), "+"), factors = 1,
    method = "kendall", center = "none"
  )
  scaled <- robust_factor_subspace(
    7 * x, factors = 1, method = "kendall", center = "mean"
  )
  permuted <- robust_factor_subspace(
    x[c(6, 2, 4, 1, 5, 3), ], factors = 1,
    method = "kendall", center = "mean"
  )

  expect_equal(kendall$operator, spatial_kendall(x), tolerance = 1e-14,
               ignore_attr = TRUE)
  expect_equal(shifted$operator, kendall$operator, tolerance = 1e-14)
  expect_equal(scaled$operator, kendall$operator, tolerance = 1e-14)
  expect_equal(permuted$operator, kendall$operator, tolerance = 1e-14)
  expect_equal(shifted$projector, kendall$projector, tolerance = 1e-12)
  expect_equal(scaled$projector, kendall$projector, tolerance = 1e-12)
  expect_true(kendall$diagnostics$zero.directions$translation.invariant.operator)

  spatial <- robust_factor_subspace(
    ch6_axis_fixture(), factors = 1, method = "spatial_sign",
    center = "spatial"
  )
  expect_equal(spatial$center, c(first = 0, second = 0), tolerance = 1e-10)
  expect_true(spatial$diagnostics$spatial.median$converged)
})


test_that("robust subspace zero directions and eigengaps are explicit", {
  with.zero <- rbind(ch6_axis_fixture(), c(0, 0))
  expect_error(
    robust_factor_subspace(with.zero, 1, method = "spatial_sign",
                           center = "none"),
    "zero spatial residual"
  )
  retained <- robust_factor_subspace(
    with.zero, 1, method = "spatial_sign", center = "none",
    zero_action = "zero"
  )
  expect_identical(retained$diagnostics$zero.directions$zero.count, 1L)
  expect_equal(retained$operator, diag(c(4 / 7, 2 / 7)), tolerance = 1e-14,
               ignore_attr = TRUE)

  with.tie <- rbind(ch6_factor_fixture(), ch6_factor_fixture()[1, ])
  expect_error(
    robust_factor_subspace(with.tie, 1, method = "kendall",
                           center = "mean"),
    "tied pairwise"
  )
  retained.tie <- robust_factor_subspace(
    with.tie, 1, method = "kendall", center = "mean",
    zero_action = "zero"
  )
  expect_identical(retained.tie$diagnostics$zero.directions$zero.count, 1L)

  isotropic <- rbind(c(1, 0), c(-1, 0), c(0, 1), c(0, -1))
  expect_error(
    robust_factor_subspace(isotropic, 1, method = "spatial_sign",
                           center = "none"),
    "repeated-eigenvalue"
  )
  expect_error(robust_factor_subspace(ch6_factor_fixture(), 3), "factors")
  expect_error(robust_factor_subspace(matrix(1:6, 6, 1), 1),
               "at least two variables")
})


test_that("RTS returns sqrt-p loadings and cross-sectional OLS scores", {
  x <- ch6_factor_fixture()
  fit <- rts_factor(x, factors = 1, center = "mean")
  centered <- sweep(x, 2L, colMeans(x), "-")

  expect_s3_class(fit, "rts_factor_fit")
  expect_s3_class(fit, "hd_factor_fit")
  expect_true(isTRUE(fit$valid))
  expect_identical(fit$factor.number, 1L)
  expect_equal(crossprod(fit$loadings) / ncol(x), diag(1),
               tolerance = 1e-14, ignore_attr = TRUE)
  expect_equal(fit$factor.scores, centered %*% fit$loadings / ncol(x),
               tolerance = 1e-14, ignore_attr = TRUE)
  expect_equal(fit$common.component,
               fit$factor.scores %*% t(fit$loadings),
               tolerance = 1e-14, ignore_attr = TRUE)
  expect_equal(fit$common.component, centered %*% fit$projector,
               tolerance = 1e-14, ignore_attr = TRUE)
  expect_equal(fit$residuals, centered - fit$common.component,
               tolerance = 1e-14, ignore_attr = TRUE)
  expect_equal(fit$residuals %*% fit$loadings, matrix(0, nrow(x), 1),
               tolerance = 1e-12, ignore_attr = TRUE)
  expect_equal(fit$fitted, fit$common.component +
                 matrix(fit$center, nrow(x), ncol(x), byrow = TRUE),
               tolerance = 1e-14, ignore_attr = TRUE)
  expect_match(fit$diagnostics$loading.normalization,
               "crossprod\\(loadings\\) / p = I")
  expect_match(fit$diagnostics$target, "not asserted equal to span\\(B\\)")
})


test_that("RTS separates translation-invariant subspace from explicit scores", {
  x <- ch6_factor_fixture()
  shift <- c(a = 4, b = -2, c = 7)
  base <- rts_factor(x, factors = 1, center = "mean")
  moved <- rts_factor(sweep(x, 2L, shift, "+"), factors = 1,
                      center = "mean")
  none <- rts_factor(sweep(x, 2L, shift, "+"), factors = 1,
                     center = "none")

  expect_equal(moved$operator, base$operator, tolerance = 1e-14)
  expect_equal(moved$principal_subspace, base$principal_subspace,
               tolerance = 1e-12)
  expect_equal(moved$factor.scores, base$factor.scores, tolerance = 1e-12)
  expect_equal(moved$common.component, base$common.component,
               tolerance = 1e-12)
  expect_equal(moved$center, base$center + shift, tolerance = 1e-14)
  expect_equal(none$factor.scores,
               sweep(x, 2L, shift, "+") %*% none$loadings / ncol(x),
               tolerance = 1e-14, ignore_attr = TRUE)
})


test_that("RTS rejects ties, unresolved factor boundaries, and invalid ranks", {
  x <- ch6_factor_fixture()
  with.tie <- rbind(x, x[1, ])
  expect_error(rts_factor(with.tie, 1), "tied pairwise")
  allowed <- rts_factor(with.tie, 1, zero_pair_action = "zero")
  expect_identical(allowed$diagnostics$zero.kendall.pairs, 1L)
  expect_error(rts_factor(matrix(1, 4, 3), 1,
                          zero_pair_action = "zero"),
               "All pairwise differences")
  expect_error(rts_factor(matrix(1:6, 6, 1), 1),
               "at least two variables")
  expect_error(rts_factor(x, 0), "factors")
  expect_error(rts_factor(x, 3), "factors")

  isotropic <- rbind(c(1, 0), c(-1, 0), c(0, 1), c(0, -1))
  expect_error(rts_factor(isotropic, 1), "repeated-eigenvalue")
})


test_that("MKER uses the stabilized adjacent-eigenvalue ratio", {
  x <- ch6_factor_fixture()
  fit <- kendall_factor_number(x, kmax = 2, c = 0.25, method = "mker")
  adjusted <- unname(fit$eigenvalues) + 0.25 / sqrt(min(nrow(x), ncol(x)))
  manual <- adjusted[1:2] / adjusted[2:3]

  expect_s3_class(fit, "kendall_factor_number")
  expect_true(isTRUE(fit$valid))
  expect_identical(fit$method, "MKER")
  expect_equal(unname(fit$adjusted.eigenvalues), adjusted,
               tolerance = 1e-14)
  expect_equal(unname(fit$criterion), manual, tolerance = 1e-14)
  expect_identical(fit$selected, as.integer(which(manual == max(manual))[1]))
  expect_identical(fit$factor.number, fit$selected)
  expect_identical(fit$tail.sums, NULL)
  expect_equal(fit$stabilization, 0.25 / sqrt(3), tolerance = 1e-14)
  expect_identical(fit$diagnostics$tie.break,
                   "smallest j attaining the exact maximum")
})


test_that("MKTCR uses the stated V-j tail sums and logarithmic ratio", {
  x <- ch6_factor_fixture()
  fit <- kendall_factor_number(x, kmax = 2, c = 0.4, method = "mktcr")
  adjusted <- unname(fit$adjusted.eigenvalues)
  m <- length(adjusted)
  tails <- vapply(0:(m - 1), function(j) sum(adjusted[(j + 1):m]),
                  numeric(1))
  manual <- vapply(1:2, function(j) {
    log1p(adjusted[j] / tails[j]) /
      log1p(adjusted[j + 1] / tails[j + 1])
  }, numeric(1))

  expect_identical(fit$method, "MKTCR")
  expect_equal(unname(fit$tail.sums), tails, tolerance = 1e-14)
  expect_equal(unname(fit$criterion), manual, tolerance = 1e-14)
  expect_identical(fit$selected, as.integer(which(manual == max(manual))[1]))
  expect_identical(names(fit$tail.sums), paste0("V", 0:(m - 1)))
  expect_identical(fit$m, min(nrow(x), ncol(x)))
  expect_identical(fit$diagnostics$N, ncol(x))
  expect_identical(fit$diagnostics$T, nrow(x))
})


test_that("Kendall factor-number selection is affine-scale and row invariant", {
  x <- ch6_factor_fixture()
  base <- kendall_factor_number(x, 2, 0.2, "mker")
  moved <- kendall_factor_number(
    sweep(x, 2L, c(5, -2, 8), "+"), 2, 0.2, "mker"
  )
  scaled <- kendall_factor_number(11 * x, 2, 0.2, "mker")
  permuted <- kendall_factor_number(x[c(4, 2, 6, 1, 5, 3), ],
                                    2, 0.2, "mker")

  expect_equal(moved$operator, base$operator, tolerance = 1e-14)
  expect_equal(scaled$operator, base$operator, tolerance = 1e-14)
  expect_equal(permuted$operator, base$operator, tolerance = 1e-14)
  expect_equal(moved$criterion, base$criterion, tolerance = 1e-14)
  expect_equal(scaled$criterion, base$criterion, tolerance = 1e-14)
  expect_equal(permuted$criterion, base$criterion, tolerance = 1e-14)
  expect_identical(moved$selected, base$selected)
  expect_identical(scaled$selected, base$selected)
  expect_identical(permuted$selected, base$selected)
})


test_that("factor-number selector tuning, ties, and zero pairs are strict", {
  x <- ch6_factor_fixture()
  expect_error(kendall_factor_number(x, 0, 0.1), "kmax")
  expect_error(kendall_factor_number(x, 3, 0.1), "kmax")
  expect_error(kendall_factor_number(x, 2.5, 0.1), "kmax")
  expect_error(kendall_factor_number(x, 2, 0), "strictly positive")
  expect_error(kendall_factor_number(x, 2, -1), "strictly positive")
  expect_error(kendall_factor_number(x, 2, NA_real_), "strictly positive")
  expect_error(
    kendall_factor_number(matrix(1, 5, 3), 2, 0.1,
                          zero_pair_action = "zero"),
    "All pairwise differences"
  )

  with.tie <- rbind(x, x[1, ])
  expect_error(kendall_factor_number(with.tie, 2, 0.1), "tied pairwise")
  allowed <- kendall_factor_number(
    with.tie, 2, 0.1, zero_pair_action = "zero"
  )
  expect_identical(allowed$diagnostics$zero.kendall.pairs, 1L)

  isotropic <- rbind(diag(3), -diag(3))
  tied <- kendall_factor_number(isotropic, 2, 0.2, method = "mker")
  expect_equal(unname(tied$criterion), c(1, 1), tolerance = 1e-14)
  expect_identical(tied$selected, 1L)
})


test_that("native CCA and RTS kernels equal direct matrix identities", {
  fixture <- ch6_cca_fixture()
  mx <- colMeans(fixture$x)
  my <- colMeans(fixture$y)
  moments <- HDElliptical:::cpp_ch6cf_cca_moments(
    fixture$x, fixture$y, mx, my, 5
  )
  xc <- sweep(fixture$x, 2L, mx, "-")
  yc <- sweep(fixture$y, 2L, my, "-")

  expect_equal(moments$centered_x, xc, tolerance = 0, ignore_attr = TRUE)
  expect_equal(moments$centered_y, yc, tolerance = 0, ignore_attr = TRUE)
  expect_equal(moments$covariance_xx, crossprod(xc) / 5,
               tolerance = 1e-14, ignore_attr = TRUE)
  expect_equal(moments$covariance_xy, crossprod(xc, yc) / 5,
               tolerance = 1e-14, ignore_attr = TRUE)
  expect_equal(moments$covariance_yy, crossprod(yc) / 5,
               tolerance = 1e-14, ignore_attr = TRUE)
  expect_error(
    HDElliptical:::cpp_ch6cf_cca_moments(fixture$x, fixture$y[-1, ],
                                        mx, my, 5),
    "same number"
  )

  x <- ch6_factor_fixture()
  loadings <- sqrt(3) * matrix(c(1, 0, 0), 3, 1)
  out <- HDElliptical:::cpp_ch6cf_rts_components(
    x, colMeans(x), loadings
  )
  centered <- sweep(x, 2L, colMeans(x), "-")
  scores <- centered %*% loadings / 3
  common <- scores %*% t(loadings)
  expect_equal(out$factor_scores, scores, tolerance = 1e-14,
               ignore_attr = TRUE)
  expect_equal(out$common_component, common, tolerance = 1e-14,
               ignore_attr = TRUE)
  expect_equal(out$residuals, centered - common, tolerance = 1e-14,
               ignore_attr = TRUE)
  expect_equal(out$loading_gram_scaled, diag(1), tolerance = 1e-14,
               ignore_attr = TRUE)
  expect_equal(out$normal_equation_residual, 0, tolerance = 1e-14)
  expect_error(
    HDElliptical:::cpp_ch6cf_rts_components(x, colMeans(x), matrix(1, 2, 1)),
    "p by K"
  )
})
