# Deterministic formula and contract tests for Chapter 6 sparse PCA/CCA.

c6spc_projector <- function(loadings) {
  tcrossprod(unname(loadings))
}


c6spc_axis_data <- function() {
  x <- rbind(
    c(4, 0, 0), c(-4, 0, 0), c(2, 0, 0), c(-2, 0, 0),
    c(0, 3, 0), c(0, -3, 0)
  )
  colnames(x) <- c("first", "second", "third")
  rownames(x) <- paste0("r", seq_len(nrow(x)))
  x
}


c6spc_orthogonal_scores <- function() {
  h1 <- sqrt(5 / 2) * c(1, -1, 0, 0, 0)
  h2 <- sqrt(5 / 6) * c(1, 1, -2, 0, 0)
  h3 <- sqrt(5 / 12) * c(1, 1, 1, -3, 0)
  h4 <- sqrt(5 / 20) * c(1, 1, 1, 1, -4)
  cbind(h1, h2, h3, h4)
}


test_that("sparse Chapter 6 public formals expose all tuning controls", {
  expect_identical(
    names(formals(truncated_power_pca)),
    c("operator", "sparsity", "components", "initial", "solver_tol",
      "solver_max_iter", "symmetry_tol", "strict", "keep_operator")
  )
  expect_identical(
    names(formals(sparse_spatial_sign_pca)),
    c("x", "sparsity", "components", "center", "divisor", "zero_action",
      "initial", "median_tol", "median_max_iter", "zero_tol",
      "solver_tol", "solver_max_iter", "symmetry_tol", "strict",
      "keep_operator")
  )
  expect_identical(
    names(formals(fantope_pca)),
    c("x", "rank", "tau", "center", "scale", "covariance_divisor",
      "rho", "solver_tol", "solver_max_iter", "symmetry_tol", "strict",
      "keep_operator")
  )
  expect_identical(
    names(formals(pmd_sparse_pca)),
    c("x", "l1_bound", "components", "center", "scale",
      "covariance_divisor", "initial", "solver_tol", "solver_max_iter",
      "symmetry_tol", "strict", "keep_operator")
  )
  expect_identical(
    names(formals(pmd_sparse_cca)),
    c("x", "y", "l1_x", "l1_y", "components", "center", "scale",
      "covariance_divisor", "initial_x", "initial_y", "solver_tol",
      "solver_max_iter", "strict", "keep_operator")
  )
})


test_that("truncated power solves a diagonal sparse eigenproblem exactly", {
  operator <- diag(c(7, 4, 1))
  dimnames(operator) <- list(c("a", "b", "c"), c("a", "b", "c"))
  fit <- truncated_power_pca(operator, sparsity = 1, components = 2)

  expect_s3_class(fit, "truncated_power_pca_fit")
  expect_s3_class(fit, "hd_pca_fit")
  expect_true(fit$valid)
  expect_equal(unname(fit$loadings), rbind(c(1, 0), c(0, 1), c(0, 0)),
               tolerance = 0)
  expect_equal(unname(fit$eigenvalues), c(7, 4), tolerance = 0)
  expect_equal(unname(fit$deflated.operator %*% fit$loadings),
               matrix(0, 3, 2), tolerance = 1e-14)
  expect_identical(fit$variable.names, c("a", "b", "c"))
  expect_null(fit$scores)
  expect_null(fit$center)
  expect_true(all(vapply(fit$diagnostics$solver,
                         function(z) isTRUE(z$certified), logical(1))))
  expect_match(fit$diagnostics$no.implicit.regularization, "no ridge")
})


test_that("truncated power resolves support and loading-sign ties deterministically", {
  fit <- truncated_power_pca(
    diag(2), sparsity = 1, initial = c(1, 1), solver_tol = 1e-12
  )
  fit.neg <- truncated_power_pca(
    diag(2), sparsity = 1, initial = c(-1, -1), solver_tol = 1e-12
  )

  expect_equal(unname(fit$loadings[, 1]), c(1, 0), tolerance = 0)
  expect_equal(fit.neg$loadings, fit$loadings, tolerance = 0)
  expect_identical(fit$diagnostics$solver[[1]]$support_size, 1L)
  expect_equal(fit$diagnostics$solver[[1]]$fixed_point_residual, 0,
               tolerance = 0)
  expect_gte(fit$loadings[fit$diagnostics$sign.anchor.coordinates[1], 1], 0)
})


test_that("untruncated power agrees with the leading eigenspace", {
  operator <- matrix(c(4, 1, 0, 1, 2, 0, 0, 0, 1), 3, 3)
  fit <- truncated_power_pca(operator, sparsity = 3, solver_tol = 1e-10)
  reference <- eigen(operator, symmetric = TRUE)$vectors[, 1, drop = FALSE]

  expect_equal(c6spc_projector(fit$loadings), c6spc_projector(reference),
               tolerance = 1e-8)
  expect_equal(unname(fit$eigenvalues),
               eigen(operator, symmetric = TRUE)$values[1], tolerance = 1e-9)
  expect_lte(fit$diagnostics$solver[[1]]$fixed_point_residual, 1e-10)
})


test_that("truncated power rejects non-PSD input and exposes solver failure", {
  expect_error(truncated_power_pca(diag(c(1, -1)), 1),
               "not positive semidefinite")
  expect_error(truncated_power_pca(matrix(c(1, 1e-3, 0, 1), 2), 1,
                                   symmetry_tol = 1e-8), "not symmetric")
  expect_error(truncated_power_pca(diag(2), 0), "between 1 and 2")
  expect_error(truncated_power_pca(diag(2), 1, initial = c(0, 0)), "nonzero")
  expect_error(truncated_power_pca(matrix(0, 2, 2), 1),
               "matrix-vector product is zero")
  expect_warning(
    invalid <- truncated_power_pca(matrix(0, 2, 2), 1, strict = FALSE),
    "matrix-vector product is zero"
  )
  expect_false(invalid$valid)
  expect_identical(invalid$rank, 0L)
  expect_identical(invalid$diagnostics$failure.stage,
                   "truncated-power solver")
})


test_that("sparse spatial-sign PCA equals the explicit SSCM plus TPM formula", {
  x <- c6spc_axis_data()
  fit <- sparse_spatial_sign_pca(
    x, sparsity = 1, center = "none", solver_tol = 1e-12
  )
  signs <- spatial_sign(x, center = rep(0, ncol(x)))
  manual <- crossprod(unclass(signs)) / nrow(x)

  expect_s3_class(fit, "sparse_spatial_sign_pca_fit")
  expect_true(fit$valid)
  expect_equal(unname(fit$operator), unname(manual), tolerance = 0)
  expect_equal(unname(diag(fit$operator)), c(2 / 3, 1 / 3, 0),
               tolerance = 2e-16)
  expect_equal(unname(fit$loadings[, 1]), c(1, 0, 0), tolerance = 0)
  expect_equal(fit$scores, x %*% fit$loadings, tolerance = 0)
  expect_identical(fit$diagnostics$divisor.value, 6L)
  expect_equal(fit$diagnostics$operator.trace, 1, tolerance = 2e-16)
})


test_that("sparse spatial-sign PCA shares zero and divisor semantics", {
  x <- rbind(c(0, 0), c(2, 0), c(-2, 0), c(0, 1))
  fit.n <- sparse_spatial_sign_pca(x, 1, center = "none")
  fit.nz <- sparse_spatial_sign_pca(
    x, 1, center = "none", divisor = "nonzero"
  )

  expect_identical(fit.n$diagnostics$n.zero.residuals, 1L)
  expect_equal(sum(diag(fit.n$operator)), 3 / 4, tolerance = 0)
  expect_equal(sum(diag(fit.nz$operator)), 1, tolerance = 0)
  expect_error(sparse_spatial_sign_pca(
    x, 1, center = "none", zero_action = "error"
  ), "zero residual")
  expect_error(sparse_spatial_sign_pca(
    matrix(0, 3, 2), 1, center = "none", divisor = "nonzero"
  ), "divisor is zero")
})


test_that("sparse spatial-sign subspaces are translation and scale invariant", {
  x <- rbind(
    c(3, 1, 0), c(-2, 1, 1), c(1, -3, 2),
    c(-1, 2, -2), c(2, 2, 1), c(-3, -1, -1)
  )
  base <- sparse_spatial_sign_pca(
    x, sparsity = 3, center = "mean", solver_tol = 1e-10
  )
  changed <- sparse_spatial_sign_pca(
    5 * sweep(x, 2, c(8, -3, 11), "+"),
    sparsity = 3, center = "mean", solver_tol = 1e-10
  )
  dense <- spatial_sign_pca(x, rank = 1, center = "mean")

  expect_equal(changed$operator, base$operator, tolerance = 3e-15)
  expect_equal(c6spc_projector(changed$loadings),
               c6spc_projector(base$loadings), tolerance = 1e-8)
  expect_equal(c6spc_projector(base$loadings),
               c6spc_projector(dense$loadings), tolerance = 1e-8)
})


test_that("Fantope tau zero returns the exact leading projector", {
  x <- c6spc_axis_data()
  fit <- fantope_pca(
    x, rank = 2, tau = 0, center = "none", solver_tol = 1e-12
  )
  target <- diag(c(1, 1, 0))

  expect_s3_class(fit, "fantope_pca_fit")
  expect_s3_class(fit, "hd_pca_fit")
  expect_true(fit$valid)
  expect_equal(unname(fit$relaxed.projector), target, tolerance = 1e-11)
  expect_equal(c6spc_projector(fit$loadings), target, tolerance = 1e-11)
  expect_equal(sum(diag(fit$relaxed.projector)), 2, tolerance = 1e-12)
  expect_gte(min(eigen(fit$relaxed.projector, symmetric = TRUE)$values),
             -1e-12)
  expect_lte(max(eigen(fit$relaxed.projector, symmetric = TRUE)$values),
             1 + 1e-12)
  expect_lte(fit$diagnostics$scaled.primal.residual, 1e-12)
  expect_lte(fit$diagnostics$scaled.dual.residual, 1e-12)
  expect_equal(fit$scores,
               sweep(x, 2, fit$center, "-") %*% fit$loadings,
               tolerance = 0)
})


test_that("Fantope retains the relaxed estimate at repeated values", {
  x <- rbind(diag(3), -diag(3))
  fit <- fantope_pca(
    x, rank = 1, tau = 0, center = "none", solver_tol = 1e-12
  )

  expect_equal(unname(fit$relaxed.projector), diag(1 / 3, 3),
               tolerance = 1e-12)
  expect_equal(sum(diag(fit$relaxed.projector)), 1, tolerance = 1e-13)
  expect_equal(sum(fit$relaxed.projector^2), 1 / 3, tolerance = 1e-12)
  expect_true(length(fit$diagnostics$repeated.relaxed.eigenvalue.groups) > 0L)
  expect_match(fit$diagnostics$relaxed.estimate.contract,
               "primary convex Fantope estimate")
})


test_that("Fantope checks both ADMM residuals and strict failure semantics", {
  x <- c6spc_axis_data()
  expect_warning(
    invalid <- fantope_pca(
      x, rank = 1, tau = 0, center = "none",
      solver_tol = 1e-14, solver_max_iter = 1, strict = FALSE
    ),
    "did not pass"
  )
  expect_false(invalid$valid)
  expect_gt(invalid$diagnostics$scaled.dual.residual, 1e-14)
  expect_identical(invalid$diagnostics$failure.stage, "Fantope ADMM")
  expect_error(fantope_pca(x, 0, 0), "rank")
  expect_error(fantope_pca(x, 1, -1), "tau")
  expect_error(fantope_pca(x, 1, 0, rho = 0), "rho")
})


test_that("PMD sparse PCA enforces the actual l1 bound and tie optimum", {
  z <- c(-2, -1, 1, 2)
  x <- outer(z, c(1, 1, 0))
  colnames(x) <- c("a", "b", "c")
  fit <- pmd_sparse_pca(
    x, l1_bound = 1, center = "none", solver_tol = 1e-12
  )
  certificate <- fit$diagnostics$solver[[1]]

  expect_s3_class(fit, "pmd_sparse_pca_fit")
  expect_true(fit$valid)
  expect_equal(unname(fit$loadings[, 1]), c(1, 0, 0), tolerance = 0)
  expect_equal(sum(abs(fit$loadings[, 1])), 1, tolerance = 0)
  expect_true(certificate$maximum_tie_branch)
  expect_true(certificate$active_l1_constraint)
  expect_lte(certificate$constraint_residual, 1e-12)
  expect_lte(certificate$bound_residual, 1e-12)
  expect_lte(certificate$scaled_v_kkt_residual, 1e-11)
})


test_that("unpenalized PMD sparse PCA agrees with the SVD projector", {
  x <- rbind(
    c(3, 1, 0), c(-2, 0, 1), c(1, -3, 1),
    c(-1, 2, -2), c(2, 1, 0), c(-3, -1, 0)
  )
  centered <- sweep(x, 2, colMeans(x), "-")
  fit <- pmd_sparse_pca(x, l1_bound = sqrt(ncol(x)), solver_tol = 1e-10)
  reference <- svd(centered, nu = 0, nv = 1)$v[, 1, drop = FALSE]

  expect_equal(c6spc_projector(fit$loadings), c6spc_projector(reference),
               tolerance = 1e-9)
  expect_equal(unname(fit$eigenvalues),
               svd(centered, nu = 0, nv = 0)$d[1]^2 / nrow(x),
               tolerance = 1e-9)
  expect_lte(fit$diagnostics$solver[[1]]$fixed_point_residual, 1e-10)
})


test_that("PMD sparse PCA residual deflation returns successive axis factors", {
  h <- c6spc_orthogonal_scores()
  x <- cbind(a = 3 * h[, 1], b = 2 * h[, 2], c = h[, 3])
  fit <- pmd_sparse_pca(
    x, l1_bound = 1, components = 2, center = "none",
    solver_tol = 1e-12
  )

  expect_equal(unname(fit$loadings), rbind(c(1, 0), c(0, 1), c(0, 0)),
               tolerance = 1e-13)
  expect_equal(unname(fit$singular.values), c(3, 2) * sqrt(5),
               tolerance = 1e-12)
  expect_equal(unname(fit$eigenvalues), c(9, 4), tolerance = 1e-12)
  expect_equal(unname(fit$residual.matrix %*% fit$loadings),
               matrix(0, 5, 2), tolerance = 1e-12)
  expect_identical(fit$diagnostics$deflation, "R <- R - d u v'")
})


test_that("PMD sparse PCA never promotes an uncertified iterate", {
  z <- c(-2, -1, 1, 2)
  x <- outer(z, c(1, 1, 0))
  expect_warning(
    invalid <- pmd_sparse_pca(
      x, l1_bound = 1, center = "none", solver_tol = 1e-14,
      solver_max_iter = 1, strict = FALSE
    ),
    "did not pass"
  )
  expect_false(invalid$valid)
  expect_identical(invalid$rank, 0L)
  expect_error(pmd_sparse_pca(x, l1_bound = 0.99), "\\[1, sqrt")
  expect_error(pmd_sparse_pca(x, l1_bound = 2), "\\[1, sqrt")
  expect_warning(
    zero <- pmd_sparse_pca(matrix(0, 4, 2), 1, center = "none",
                           strict = FALSE),
    "leading-SVD initialization failed"
  )
  expect_false(zero$valid)
})


test_that("PMD sparse CCA matches an explicit diagonal cross operator", {
  h <- c6spc_orthogonal_scores()
  x <- cbind(x1 = h[, 1], x2 = h[, 2])
  y <- cbind(y1 = 3 * h[, 1], y2 = 2 * h[, 2])
  fit <- pmd_sparse_cca(
    x, y, l1_x = 1, l1_y = 1, components = 2,
    center = FALSE, scale = FALSE, solver_tol = 1e-12
  )

  expect_s3_class(fit, "pmd_sparse_cca_fit")
  expect_s3_class(fit, "hd_pca_fit")
  expect_true(fit$valid)
  expect_equal(unname(fit$operator), diag(c(3, 2)), tolerance = 2e-15)
  expect_equal(unname(fit$loadings), diag(2), tolerance = 1e-13)
  expect_equal(unname(fit$y.loadings), diag(2), tolerance = 1e-13)
  expect_equal(unname(fit$singular.values), c(3, 2), tolerance = 1e-13)
  expect_equal(unname(fit$canonical.correlations), c(1, 1),
               tolerance = 2e-15)
  expect_false(isTRUE(all.equal(fit$singular.values,
                                fit$canonical.correlations)))
  expect_true(all(vapply(fit$diagnostics$solver,
                         function(z) isTRUE(z$certified), logical(1))))
  expect_identical(fit$diagnostics$cross.operator.evaluation,
                   "matrix-free X'(Yv) and Y'(Xu), plus explicit rank-one deflations")
})


test_that("PMD sparse CCA handles maximum ties without a norm hack", {
  z <- c(-2, -1, 1, 2)
  x <- cbind(a = z, b = z)
  y <- cbind(c = z, d = z)
  fit <- pmd_sparse_cca(
    x, y, l1_x = 1, l1_y = 1,
    center = FALSE, scale = FALSE, solver_tol = 1e-12
  )
  certificate <- fit$diagnostics$solver[[1]]

  expect_equal(unname(fit$loadings[, 1]), c(1, 0), tolerance = 0)
  expect_equal(unname(fit$y.loadings[, 1]), c(1, 0), tolerance = 0)
  expect_true(certificate$x_maximum_tie_branch)
  expect_true(certificate$y_maximum_tie_branch)
  expect_lte(certificate$x_bound_residual, 1e-12)
  expect_lte(certificate$y_bound_residual, 1e-12)
})


test_that("scaled PMD CCA is invariant to shifts and positive feature units", {
  h <- c6spc_orthogonal_scores()
  x <- cbind(x1 = h[, 1] + 0.3 * h[, 3], x2 = h[, 2])
  y <- cbind(y1 = 0.8 * h[, 1] + 0.6 * h[, 3],
             y2 = 0.3 * h[, 2] + sqrt(0.91) * h[, 4])
  base <- pmd_sparse_cca(
    x, y, l1_x = sqrt(2), l1_y = sqrt(2), scale = TRUE,
    solver_tol = 1e-10
  )
  changed <- pmd_sparse_cca(
    sweep(sweep(x, 2, c(4, 0.5), "*"), 2, c(10, -7), "+"),
    sweep(sweep(y, 2, c(2, 5), "*"), 2, c(-3, 20), "+"),
    l1_x = sqrt(2), l1_y = sqrt(2), scale = TRUE,
    solver_tol = 1e-10
  )

  expect_equal(c6spc_projector(changed$loadings),
               c6spc_projector(base$loadings), tolerance = 1e-8)
  expect_equal(c6spc_projector(changed$y.loadings),
               c6spc_projector(base$y.loadings), tolerance = 1e-8)
  expect_equal(changed$canonical.correlations,
               base$canonical.correlations, tolerance = 1e-10)
})


test_that("PMD sparse CCA validates paired data and solver degeneracy", {
  h <- c6spc_orthogonal_scores()
  x <- cbind(a = h[, 1])
  y <- cbind(b = h[, 2])
  expect_error(pmd_sparse_cca(x, y[-1, , drop = FALSE], 1, 1),
               "same number of rows")
  expect_error(pmd_sparse_cca(x, y, 1, 1, initial_x = 1),
               "supplied together")
  expect_error(pmd_sparse_cca(x, y, 0.5, 1), "\\[1, sqrt")
  expect_warning(
    invalid <- pmd_sparse_cca(
      x, y, 1, 1, center = FALSE, scale = FALSE, strict = FALSE
    ),
    "cross-covariance operator is zero"
  )
  expect_false(invalid$valid)
  expect_identical(invalid$diagnostics$failure.stage,
                   "PMD CCA alternating solver")
  expect_error(pmd_sparse_cca(
    cbind(x, constant = 1), cbind(y, other = 1), sqrt(2), sqrt(2),
    scale = TRUE
  ), "constant")
})


test_that("all successful fits implement the common hd_pca_fit contract", {
  x <- c6spc_axis_data()
  fits <- list(
    truncated_power_pca(diag(c(3, 2, 1)), 1),
    sparse_spatial_sign_pca(x, 1, center = "none"),
    fantope_pca(x, 1, 0, center = "none"),
    pmd_sparse_pca(x, 1, center = "none")
  )
  required <- c(
    "method", "eigenvalues", "loadings", "scores", "center", "operator",
    "rank", "n", "p", "variable.names", "diagnostics", "call"
  )

  for (fit in fits) {
    expect_s3_class(fit, "hd_pca_fit")
    expect_true(fit$valid)
    expect_true(all(required %in% names(fit)))
    for (j in seq_len(ncol(fit$loadings))) {
      anchor <- which(abs(fit$loadings[, j]) ==
                        max(abs(fit$loadings[, j])))[1L]
      expect_gte(fit$loadings[anchor, j], 0)
    }
  }
})
