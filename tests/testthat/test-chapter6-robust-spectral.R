# Deterministic tests for Chapter 6 robust spectral estimators.

c6rs_projector <- function(loadings) {
  unname(tcrossprod(loadings))
}


c6rs_plain_matrix <- function(x) {
  matrix(as.numeric(x), nrow(x), ncol(x))
}


c6rs_manual_kendall <- function(x, zero_tol = 0, nonzero_divisor = FALSE) {
  p <- ncol(x)
  total <- choose(nrow(x), 2)
  zero <- 0L
  answer <- matrix(0, p, p)
  for (i in seq_len(nrow(x) - 1L)) {
    for (j in seq.int(i + 1L, nrow(x))) {
      difference <- x[i, ] - x[j, ]
      radius <- sqrt(sum(difference^2))
      if (radius <= zero_tol) {
        zero <- zero + 1L
      } else {
        direction <- difference / radius
        answer <- answer + tcrossprod(direction)
      }
    }
  }
  denominator <- if (nonzero_divisor) total - zero else total
  list(matrix = answer / denominator, pairs = total, zero = zero)
}


c6rs_axis_data <- function(radii) {
  answer <- cbind(radii, rep(0, length(radii)))
  colnames(answer) <- c("x", "y")
  answer
}


test_that("spatial-sign PCA is exactly the Chapter 1 SSCM eigendecomposition", {
  x <- rbind(
    c(3, 0, 0), c(-3, 0, 0), c(0, 2, 0),
    c(0, -2, 0), c(0, 0, 1), c(0, 0, -1)
  )
  colnames(x) <- c("a", "b", "c")
  fit <- spatial_sign_pca(x, rank = 2, center = "none")
  manual <- sscm(x, center = "none")

  expect_s3_class(fit, "spatial_sign_pca_fit")
  expect_s3_class(fit, "hd_pca_fit")
  expect_equal(unname(fit$operator), c6rs_plain_matrix(manual),
               tolerance = 0)
  expect_equal(fit$eigenvalues, eigen(manual, symmetric = TRUE)$values,
               tolerance = 1e-14)
  expect_equal(fit$scores, x %*% fit$loadings, tolerance = 0)
  expect_identical(rownames(fit$loadings), colnames(x))
  expect_identical(colnames(fit$loadings), c("PC1", "PC2"))
  expect_identical(fit$diagnostics$divisor.value, 6L)
  expect_identical(fit$diagnostics$divisor.rule, "n")
  expect_equal(fit$diagnostics$operator.trace, 1, tolerance = 1e-15)
  for (j in seq_len(ncol(fit$loadings))) {
    anchor <- which(abs(fit$loadings[, j]) ==
                      max(abs(fit$loadings[, j])))[1L]
    expect_gte(fit$loadings[anchor, j], 0)
  }
})


test_that("spatial-sign zero contributions and both divisors are explicit", {
  x <- rbind(c(0, 0), c(2, 0), c(-2, 0), c(0, 1))
  default <- spatial_sign_pca(x, rank = 1, center = "none")
  nonzero <- spatial_sign_pca(
    x, rank = 1, center = "none", divisor = "nonzero"
  )

  expect_identical(default$diagnostics$n.zero.residuals, 1L)
  expect_equal(sum(diag(default$operator)), 3 / 4, tolerance = 0)
  expect_equal(sum(diag(nonzero$operator)), 1, tolerance = 0)
  expect_identical(nonzero$diagnostics$divisor.value, 3L)
  expect_error(
    spatial_sign_pca(x, rank = 1, center = "none", zero_action = "error"),
    "zero residual"
  )
  expect_error(
    spatial_sign_pca(
      matrix(0, 3, 2), rank = 1, center = "none", divisor = "nonzero"
    ),
    "divisor is zero"
  )
  boundary <- spatial_sign_pca(
    c6rs_axis_data(c(1, 2)), rank = 1, center = "none", zero_tol = 1
  )
  expect_identical(boundary$diagnostics$n.zero.residuals, 1L)
  expect_equal(unname(boundary$operator), diag(c(1 / 2, 0)),
               tolerance = 0)

})


test_that("spatial-sign PCA is translation, positive-scale, and row invariant", {
  x <- rbind(
    c(3, 1, 0), c(-2, 1, 1), c(1, -3, 2),
    c(-1, 2, -2), c(2, 2, 1), c(-3, -1, -1)
  )
  colnames(x) <- c("a", "b", "c")
  shift <- c(10, -7, 4)
  base <- spatial_sign_pca(x, rank = 2, center = "mean")
  moved.x <- 5 * sweep(x, 2, shift, "+")
  moved <- spatial_sign_pca(moved.x, rank = 2, center = "mean")
  permutation <- c(4, 1, 6, 2, 5, 3)
  permuted <- spatial_sign_pca(
    x[permutation, ], rank = 2, center = "mean"
  )

  expect_equal(moved$operator, base$operator, tolerance = 2e-14)
  expect_equal(c6rs_projector(moved$loadings),
               c6rs_projector(base$loadings), tolerance = 2e-13)
  expect_equal(moved$scores %*% t(moved$loadings),
               5 * (base$scores %*% t(base$loadings)), tolerance = 2e-12)
  expect_equal(permuted$operator, base$operator, tolerance = 2e-14)
})


test_that("repeated spatial-sign eigenvalues are tested through their projector", {
  x <- rbind(
    c(1, 0, 0), c(-1, 0, 0), c(2, 0, 0), c(-2, 0, 0),
    c(0, 1, 0), c(0, -1, 0), c(0, 2, 0), c(0, -2, 0),
    c(0, 0, 1), c(0, 0, -1)
  )
  fit <- spatial_sign_pca(x, rank = 2, center = "none")

  expect_equal(fit$eigenvalues, c(0.4, 0.4, 0.2), tolerance = 1e-15)
  expect_equal(c6rs_projector(fit$loadings), diag(c(1, 1, 0)),
               tolerance = 1e-14)
  expect_true(any(vapply(
    fit$diagnostics$repeated.eigenvalue.groups,
    function(index) identical(as.integer(index), 1:2), logical(1)
  )))
})


test_that("Kendall PCA matches direct enumeration of every unordered pair", {
  x <- rbind(c(2, 0), c(-1, 0), c(0, 2), c(0, -3))
  colnames(x) <- c("u", "v")
  manual <- c6rs_manual_kendall(x)
  fit <- kendall_pca(x, rank = 2, center = "mean")

  expect_equal(unname(fit$operator), manual$matrix, tolerance = 2e-16)
  expect_identical(fit$diagnostics$n.pairs, 6L)
  expect_identical(fit$diagnostics$n.zero.pairs, 0L)
  expect_equal(sum(diag(fit$operator)), 1, tolerance = 2e-16)
  expect_equal(fit$scores,
               sweep(x, 2, colMeans(x), "-") %*% fit$loadings,
               tolerance = 2e-16)
  expect_identical(
    fit$diagnostics$operator.center,
    "none: pairwise differences are translation invariant"
  )
})


test_that("Kendall ties have zero, error, and nonzero-divisor contracts", {
  x <- rbind(c(0, 0), c(0, 0), c(1, 0))
  default <- kendall_pca(x, rank = 1, center = "none")
  omitted <- kendall_pca(
    x, rank = 1, center = "none", divisor = "nonzero_pairs"
  )

  expect_equal(unname(default$operator), diag(c(2 / 3, 0)), tolerance = 0)
  expect_equal(unname(omitted$operator), diag(c(1, 0)), tolerance = 0)
  expect_identical(default$diagnostics$n.zero.pairs, 1L)
  expect_error(kendall_pca(x, rank = 1, ties = "error"), "tied")
  expect_error(
    kendall_pca(
      matrix(1, 3, 2), rank = 1, divisor = "nonzero_pairs"
    ),
    "divisor is zero"
  )
})


test_that("Kendall operator and subspace obey affine nuisance invariances", {
  x <- rbind(
    c(2, 0, 1), c(-1, 2, 0), c(0, -3, 2),
    c(4, 1, -1), c(-2, -1, 3)
  )
  colnames(x) <- c("a", "b", "c")
  base <- kendall_pca(x, rank = 2)
  order <- c(5, 2, 4, 1, 3)
  changed.x <- 7 * sweep(x[order, ], 2, c(20, -8, 11), "+")
  changed <- kendall_pca(changed.x, rank = 2)

  expect_equal(changed$operator, base$operator, tolerance = 3e-15)
  expect_equal(c6rs_projector(changed$loadings),
               c6rs_projector(base$loadings), tolerance = 2e-14)
  expect_equal(
    changed$scores %*% t(changed$loadings),
    7 * (base$scores[order, ] %*% t(base$loadings)),
    tolerance = 2e-13
  )
})


test_that("Kendall PCA is equivariant under signed feature permutations", {
  x <- rbind(
    c(3, 1, 0), c(-2, 2, 1), c(1, -4, 2),
    c(-1, 1, -3), c(2, 3, 1)
  )
  transform <- matrix(c(
    0, 0, -1,
    1, 0, 0,
    0, 1, 0
  ), 3, 3, byrow = TRUE)
  base <- kendall_pca(x, rank = 2, center = "none")
  changed <- kendall_pca(x %*% transform, rank = 2, center = "none")

  expect_equal(unname(changed$operator),
               unname(t(transform) %*% base$operator %*% transform),
               tolerance = 3e-15)
  expect_equal(c6rs_projector(changed$loadings),
               t(transform) %*% c6rs_projector(base$loadings) %*% transform,
               tolerance = 2e-14)
  expect_true(all(vapply(seq_len(2), function(j) {
    anchor <- which(abs(changed$loadings[, j]) ==
                      max(abs(changed$loadings[, j])))[1L]
    changed$loadings[anchor, j] >= 0
  }, logical(1))))
})


test_that("all GSSCM radial families obey their exact endpoint formulas", {
  x <- c6rs_axis_data(c(1, 2, 3, 4, 5))
  cutoffs <- c(Q1 = 1, Q2 = 2, Q3 = 3, Q3_star = 4)
  expected <- list(
    winsor = c(1, 2, 2, 2, 2),
    quadratic = c(1, 2, 4 / 3, 1, 4 / 5),
    ball = c(1, 2, 0, 0, 0),
    shell = c(1, 2, 3, 0, 0),
    linear_redescending = c(1, 2, 3 / 2, 0, 0)
  )

  for (weight in names(expected)) {
    fit <- generalized_sign_pca(
      x, rank = 1, weight = weight, center = "none",
      cutoff = "user", cutoffs = cutoffs
    )
    expect_equal(fit$diagnostics$transformed.norms, expected[[weight]],
                 tolerance = 2e-15, info = weight)
    expect_equal(unname(fit$operator[1, 1]),
                 mean(expected[[weight]]^2), tolerance = 2e-15,
                 info = weight)
  }
  shell <- generalized_sign_pca(
    x, rank = 1, weight = "shell", center = "none",
    cutoff = "user", cutoffs = cutoffs
  )
  expect_identical(
    shell$diagnostics$boundary.rule,
    "Ball includes Q2; Shell includes Q1 and Q3; linear redescending includes its inner endpoint"
  )
})


test_that("the default GSPCA cutoffs use ordinary median and raw MAD", {
  transformed.radii <- c(1, 2, 3, 4, 100)
  radii <- transformed.radii^(3 / 2)
  x <- c6rs_axis_data(radii)
  fit <- generalized_sign_pca(
    x, rank = 1, weight = "ball", center = "none",
    cutoff = "median_mad"
  )
  expected <- c(
    Q1 = 2^(3 / 2), Q2 = 3^(3 / 2), Q3 = 4^(3 / 2),
    Q3_star = (3 + 1.4826)^(3 / 2)
  )

  expect_equal(fit$diagnostics$cutoffs, expected, tolerance = 2e-12)
  expect_equal(fit$diagnostics$transformed.radius.center, 3,
               tolerance = 2e-15)
  expect_equal(fit$diagnostics$transformed.radius.mad, 1,
               tolerance = 2e-15)
  expect_identical(fit$diagnostics$cutoff.method, "median_mad")
  expect_match(fit$diagnostics$cutoff.source, "2024")
  even.transformed <- c(1, 2, 4, 8)
  even <- generalized_sign_pca(
    c6rs_axis_data(even.transformed^(3 / 2)), rank = 1,
    weight = "ball", center = "none"
  )
  expect_equal(unname(even$diagnostics$cutoffs["Q2"]),
               stats::median(even.transformed^(3 / 2)), tolerance = 2e-14)
  expect_equal(even$diagnostics$transformed.radius.center, 3, tolerance = 2e-15)

})


test_that("the original h-order cutoff is exact and fails when h exceeds n", {
  transformed.radii <- 1:6
  x <- c6rs_axis_data(transformed.radii^(3 / 2))
  fit <- generalized_sign_pca(
    x, rank = 1, weight = "winsor", center = "none",
    cutoff = "original_h_order"
  )
  expected <- c(
    Q1 = 2^(3 / 2), Q2 = 4^(3 / 2), Q3 = 6^(3 / 2),
    Q3_star = (4 + 1.4826 * 2)^(3 / 2)
  )

  expect_identical(fit$diagnostics$original.h, 4L)
  expect_equal(fit$diagnostics$cutoffs, expected, tolerance = 3e-14)
  expect_match(fit$diagnostics$cutoff.source, "2019")

  wide <- matrix(seq_len(40), 4, 10)
  expect_error(
    generalized_sign_pca(
      wide, rank = 1, weight = "ball", center = "none",
      cutoff = "original_h_order"
    ),
    "exceeds n"
  )
})


test_that("zero MAD uses an exact limiting rule or an explicit error", {
  x <- rbind(c(2, 0), c(-2, 0), c(0, 2), c(0, -2))
  lr <- generalized_sign_pca(
    x, rank = 2, weight = "linear_redescending", center = "none"
  )
  shell <- generalized_sign_pca(
    x, rank = 2, weight = "shell", center = "none"
  )

  expect_true(lr$diagnostics$zero.mad)
  expect_true(all(lr$diagnostics$cutoffs == 2))
  expect_equal(lr$diagnostics$transformed.norms, rep(2, 4), tolerance = 0)
  expect_equal(shell$diagnostics$transformed.norms, rep(2, 4), tolerance = 0)
  expect_match(lr$diagnostics$zero.mad.rule, "coincident")
  expect_error(
    generalized_sign_pca(
      x, rank = 1, weight = "ball", center = "none",
      zero_mad = "error"
    ),
    "zero MAD"
  )
})


test_that("Q1 lower-bound and user-cutoff choices are never implicit", {
  lowered <- HDElliptical:::.c6rs_q1(-1, "zero")
  expect_equal(lowered$value, 0, tolerance = 0)
  expect_true(lowered$adjusted)
  expect_error(HDElliptical:::.c6rs_q1(-1, "error"), "negative")

  x <- c6rs_axis_data(1:4)
  expect_error(
    generalized_sign_pca(
      x, rank = 1, weight = "ball", center = "none", cutoff = "user"
    ),
    "required"
  )
  expect_error(
    generalized_sign_pca(
      x, rank = 1, weight = "ball", center = "none", cutoff = "user",
      cutoffs = c(Q1 = 2, Q2 = 1, Q3 = 3, Q3_star = 4)
    ),
    "Q1 <= Q2"
  )
  expect_error(
    generalized_sign_pca(
      x, rank = 1, weight = "ball", center = "none", cutoff = "user",
      cutoffs = c(1, 2, 3)
    ),
    "four"
  )
  expect_error(
    generalized_sign_pca(
      x, rank = 1, weight = "ball", center = "none", cutoff = "user",
      cutoffs = c(1, 2, 3, 4)
    ),
    "must be named"
  )
  expect_error(
    generalized_sign_pca(x, rank = 1, weight = "identity", center = "none",
                         cutoff = "median_mad"), "not used"
  )

})


test_that("identity GSPCA is the n-divisor centered second moment", {
  x <- rbind(c(3, 1), c(-2, 2), c(1, -4), c(-1, 0), c(2, 3))
  colnames(x) <- c("a", "b")
  center <- colMeans(x)
  fit <- generalized_sign_pca(
    x, rank = 2, weight = "identity", center = center
  )
  centered <- sweep(x, 2, center, "-")
  expected <- crossprod(centered) / nrow(x)
  scaled <- generalized_sign_pca(
    6 * x, rank = 2, weight = "identity", center = 6 * center
  )

  expect_equal(unname(fit$operator), unname(expected), tolerance = 3e-15)
  expect_equal(fit$scores, centered %*% fit$loadings, tolerance = 3e-15)
  expect_equal(unname(scaled$operator), 36 * unname(expected), tolerance = 1e-13)
  expect_equal(c6rs_projector(scaled$loadings),
               c6rs_projector(fit$loadings), tolerance = 2e-14)
  expect_null(fit$diagnostics$cutoffs)
  expect_identical(fit$diagnostics$cutoff.method, "none")
})


test_that("spatial GSPCA reduces exactly to spatial-sign PCA", {
  x <- rbind(
    c(3, 1, 0), c(-2, 2, 1), c(1, -4, 2),
    c(-1, 1, -3), c(2, 3, 1)
  )
  center <- colMeans(x)
  sign <- spatial_sign_pca(x, rank = 2, center = center)
  generalized <- generalized_sign_pca(
    x, rank = 2, weight = "spatial", center = center
  )

  expect_equal(generalized$operator, sign$operator, tolerance = 1e-15)
  expect_equal(c6rs_projector(generalized$loadings),
               c6rs_projector(sign$loadings), tolerance = 2e-14)
  expect_equal(generalized$scores %*% t(generalized$loadings),
               sign$scores %*% t(sign$loadings), tolerance = 2e-14)
  expect_true(all(generalized$diagnostics$transformed.norms == 1))
  extreme <- generalized_sign_pca(
    cbind(c(1e20, -1e20, 1e-300, -1e-300), 0),
    rank = 1, weight = "spatial", center = "none"
  )
  expect_true(all(is.finite(extreme$operator)))
  expect_equal(unname(extreme$operator), diag(c(1, 0)), tolerance = 0)
  expect_true(all(is.finite(extreme$diagnostics$radial.multiplier)))
  expect_equal(extreme$diagnostics$radial.multiplier[3:4],
               rep(1e300, 2), tolerance = 5e-5)

})


test_that("k-step LTS center follows deterministic C-steps", {
  x <- rbind(
    c(0, 0), c(1, 0), c(0, 1), c(1, 1), c(2, 2), c(30, -10)
  )
  initial <- spatial_median(x, tol = 1e-10, warn = FALSE)
  current <- as.numeric(initial)
  h <- floor((nrow(x) + 1) / 2)
  for (step in 1:2) {
    signed <- spatial_sign(x, center = current)
    radii <- attr(signed, "norms")
    selected <- order(radii, seq_along(radii))[seq_len(h)]
    current <- colMeans(x[selected, , drop = FALSE])
  }
  fit <- generalized_sign_pca(
    x, rank = 1, weight = "identity", center = "kstep_lts",
    lts_steps = 2, tol = 1e-10
  )
  shifted <- generalized_sign_pca(
    sweep(x, 2, c(20, -7), "+"), rank = 1, weight = "identity",
    center = "kstep_lts", lts_steps = 2, tol = 1e-10
  )

  expect_equal(unname(fit$center), current, tolerance = 2e-14)
  expect_identical(fit$diagnostics$center$h, 3L)
  expect_length(fit$diagnostics$center$path, 2L)
  expect_equal(unname(shifted$center), current + c(20, -7),
               tolerance = 2e-13)
  default.fit <- generalized_sign_pca(
    x, rank = 1, weight = "identity", tol = 1e-10
  )
  expect_identical(default.fit$diagnostics$center$type, "k-step LTS")
  expect_identical(default.fit$diagnostics$center$steps, 2L)

})


test_that("GSPCA cutoffs, operators, and projectors are equivariant", {
  x <- rbind(
    c(3, 1, 0), c(-3, -1, 0), c(1, -4, 2),
    c(-1, 4, -2), c(2, 2, 1), c(-2, -2, -1)
  )
  transform <- matrix(c(
    0, -1, 0,
    0, 0, 1,
    1, 0, 0
  ), 3, 3, byrow = TRUE)
  base <- generalized_sign_pca(
    x, rank = 2, weight = "winsor", center = "none"
  )
  scaled <- generalized_sign_pca(
    7 * x, rank = 2, weight = "winsor", center = "none"
  )
  rotated <- generalized_sign_pca(
    x %*% transform, rank = 2, weight = "winsor", center = "none"
  )

  expect_equal(scaled$diagnostics$cutoffs,
               7 * base$diagnostics$cutoffs, tolerance = 2e-13)
  expect_equal(unname(scaled$operator), 49 * unname(base$operator),
               tolerance = 3e-12)
  expect_equal(c6rs_projector(scaled$loadings),
               c6rs_projector(base$loadings), tolerance = 2e-13)
  expect_equal(unname(rotated$operator),
               unname(t(transform) %*% base$operator %*% transform),
               tolerance = 2e-14)
  expect_equal(c6rs_projector(rotated$loadings),
               t(transform) %*% c6rs_projector(base$loadings) %*% transform,
               tolerance = 2e-13)
})


test_that("all three fits obey the shared minimal PCA object contract", {
  x <- rbind(
    c(3, 0, 1), c(-3, 0, -1), c(0, 2, 1),
    c(0, -2, -1), c(1, 1, 2), c(-1, -1, -2)
  )
  colnames(x) <- c("alpha", "beta", "gamma")
  rownames(x) <- paste0("r", seq_len(nrow(x)))
  fits <- list(
    spatial_sign_pca(x, rank = 2, center = "none", keep_operator = FALSE),
    kendall_pca(x, rank = 2, center = "none", keep_operator = FALSE),
    generalized_sign_pca(
      x, rank = 2, weight = "ball", center = "none",
      keep_operator = FALSE
    )
  )

  for (fit in fits) {
    expect_s3_class(fit, "hd_pca_fit")
    expect_identical(class(fit)[length(class(fit))], "list")
    expect_length(fit$eigenvalues, 3L)
    expect_identical(dim(fit$loadings), c(3L, 2L))
    expect_identical(dim(fit$scores), c(6L, 2L))
    expect_identical(rownames(fit$loadings), colnames(x))
    expect_identical(rownames(fit$scores), rownames(x))
    expect_identical(colnames(fit$loadings), c("PC1", "PC2"))
    expect_identical(colnames(fit$scores), c("PC1", "PC2"))
    expect_identical(fit$variable.names, colnames(x))
    expect_identical(fit$rank, 2L)
    expect_identical(fit$n, 6L)
    expect_identical(fit$p, 3L)
    expect_null(fit$operator)
    expect_match(fit$diagnostics$no.implicit.regularization,
                 "no ridge")
  }
})


test_that("robust spectral APIs reject invalid controls and data", {
  x <- matrix(seq_len(18), 6, 3)
  expect_error(spatial_sign_pca(x, rank = 0), "at least")
  expect_error(kendall_pca(x, rank = 4), "cannot exceed")
  expect_error(generalized_sign_pca(x, rank = 1, weight = "bad"),
               "arg")
  expect_error(spatial_sign_pca(x, rank = 1, zero_tol = -1),
               "non-negative")
  expect_error(kendall_pca(x, rank = 1, keep_operator = NA),
               "TRUE or FALSE")
  expect_error(generalized_sign_pca(
    x, rank = 1, center = "kstep_lts", lts_steps = 0
  ), "at least")
  bad <- x
  bad[1, 1] <- Inf
  expect_error(spatial_sign_pca(bad, rank = 1), "finite")
  names <- x
  colnames(names) <- c("a", "a", "b")
  expect_error(kendall_pca(names, rank = 1), "unique")
})


test_that("uncertified spatial medians are never promoted", {
  x <- rbind(
    c(0, 0), c(1, 0), c(0, 1), c(20, 3), c(-4, 10), c(2, -7)
  )
  expect_error(
    spatial_sign_pca(
      x, rank = 1, center = "spatial", tol = 1e-15, max_iter = 1
    ),
    "spatial median"
  )
  expect_error(
    generalized_sign_pca(
      x, rank = 1, center = "kstep_lts", tol = 1e-15, max_iter = 1
    ),
    "spatial median"
  )
})


test_that("a zero operator remains zero without a hidden eigenvalue floor", {
  x <- matrix(rep(c(1, 2), 4), 4, 2, byrow = TRUE)
  fit <- spatial_sign_pca(x, rank = 2)

  expect_equal(unname(fit$operator), matrix(0, 2, 2), tolerance = 0)
  expect_equal(fit$eigenvalues, c(0, 0), tolerance = 0)
  expect_equal(c6rs_projector(fit$loadings), diag(2), tolerance = 0)
  expect_equal(unname(fit$scores), matrix(0, 4, 2), tolerance = 0)
  expect_true(any(vapply(
    fit$diagnostics$repeated.eigenvalue.groups,
    function(index) identical(as.integer(index), 1:2), logical(1)
  )))
})


test_that("GSPCA reports the radial term in its population scope", {
  x <- rbind(c(2, 0), c(-2, 0), c(0, 1), c(0, -1))
  fit <- generalized_sign_pca(
    x, rank = 1, weight = "winsor", center = "none"
  )
  expect_match(fit$diagnostics$population.scope, "radial variable")
  expect_match(fit$diagnostics$population.scope, "not asserted")
})
