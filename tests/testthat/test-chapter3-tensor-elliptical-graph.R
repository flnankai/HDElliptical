teg_fixture_matrix <- function(dims = c(2L, 3L), half_n = 6L) {
  pstar <- prod(dims)
  base <- outer(seq_len(half_n), seq_len(pstar), function(i, j) {
    sin((i + 1.3) * (j + 0.7)) +
      cos((2 * i + j) / 3) + i * j / (11 * half_n * pstar)
  })
  rbind(base, -base)
}

teg_as_observation_array <- function(x, dims) {
  array(as.numeric(x), dim = c(nrow(x), dims))
}

teg_as_tensor_list <- function(x, dims) {
  lapply(seq_len(nrow(x)), function(i) array(x[i, ], dim = dims))
}

teg_unfold <- function(x, dims, mode) {
  order <- c(mode, setdiff(seq_along(dims), mode))
  matrix(aperm(array(x, dim = dims), order), nrow = dims[mode])
}

teg_mode_transform <- function(x, dims, matrices) {
  transform <- Reduce(kronecker, rev(matrices))
  x %*% t(transform)
}

teg_frob_normalize <- function(x) {
  x / sqrt(sum(x * x))
}


test_that("column-major matricization and pilot crossproducts match literals", {
  signs <- matrix(1:12, nrow = 2, byrow = TRUE)
  dims <- c(2L, 3L)
  got <- HDElliptical:::cpp_ch3teg_mode_crossproducts(signs, dims)

  first <- matrix(1:6, nrow = 2)
  second <- matrix(7:12, nrow = 2)
  mode1 <- 2 / 2 * (tcrossprod(first) + tcrossprod(second))
  mode2 <- 3 / 2 * (crossprod(first) + crossprod(second))
  expect_equal(got[[1]], mode1, tolerance = 1e-14)
  expect_equal(got[[2]], mode2, tolerance = 1e-14)
  expect_equal(teg_unfold(1:6, dims, 1), first)
  expect_equal(teg_unfold(1:6, dims, 2), t(first))
  expect_equal(as.vector(first), 1:6)
})


test_that("other-mode whitening matches a literal two-mode reference", {
  signs <- rbind(1:6, c(2, -1, 3, -2, 4, -3)) / 10
  dims <- c(2L, 3L)
  a1 <- matrix(c(1.2, 0.1, 0.1, 0.8), 2)
  a2 <- matrix(c(1.1, 0.1, 0, 0.1, 0.9, 0.05,
                 0, 0.05, 1.3), 3)
  got <- HDElliptical:::cpp_ch3teg_whitened_mode_scatter(
    signs, dims, list(a1, a2), TRUE
  )

  m <- lapply(seq_len(nrow(signs)), function(i) {
    matrix(signs[i, ], nrow = 2)
  })
  v1 <- lapply(m, function(value) value %*% t(a2))
  v2 <- lapply(m, function(value) t(a1 %*% value))
  s1 <- 2 / nrow(signs) * Reduce(`+`, lapply(v1, tcrossprod))
  s2 <- 3 / nrow(signs) * Reduce(`+`, lapply(v2, tcrossprod))

  expect_equal(got$scatter[[1]], s1, tolerance = 1e-14)
  expect_equal(got$scatter[[2]], s2, tolerance = 1e-14)
  expect_equal(got$vectors[[1]][1, ], as.vector(v1[[1]]), tolerance = 1e-14)
  expect_equal(
    teg_unfold(got$vectors[[2]][1, ], dims, 2), v2[[1]],
    tolerance = 1e-14
  )
  expect_equal(got$scatter[[1]], t(got$scatter[[1]]), tolerance = 0)
  expect_equal(got$scatter[[2]], t(got$scatter[[2]]), tolerance = 0)
})


test_that("three-mode whitening follows the reverse Kronecker vec identity", {
  dims <- c(2L, 2L, 2L)
  signs <- rbind(seq_len(8), c(2, -1, 4, -3, 6, -5, 8, -7)) / 20
  roots <- list(
    matrix(c(1.1, 0.1, 0.1, 0.9), 2),
    matrix(c(0.8, -0.05, -0.05, 1.2), 2),
    matrix(c(1.3, 0.2, 0.2, 0.7), 2)
  )
  got <- HDElliptical:::cpp_ch3teg_whitened_mode_scatter(
    signs, dims, roots, TRUE
  )
  for (target in seq_along(dims)) {
    multipliers <- roots
    multipliers[[target]] <- diag(dims[target])
    transform <- Reduce(kronecker, rev(multipliers))
    whitened <- signs %*% t(transform)
    reference <- dims[target] / nrow(signs) *
      Reduce(`+`, lapply(seq_len(nrow(signs)), function(i) {
        tcrossprod(teg_unfold(whitened[i, ], dims, target))
      }))
    expect_equal(got$vectors[[target]], whitened, tolerance = 2e-14)
    expect_equal(got$scatter[[target]], reference, tolerance = 2e-14)
  }
})


test_that("off-diagonal glasso matches the two-dimensional exact solution", {
  scatter <- matrix(c(1, 0.4, 0.4, 2), 2)
  rho <- 0.1
  fit <- HDElliptical:::cpp_ch3teg_offdiag_glasso(
    scatter, rho, 1e-8, 10000L, 1, 100L
  )
  dual <- matrix(c(1, 0.3, 0.3, 2), 2)
  reference <- solve(dual)
  objective <- sum(diag(scatter %*% reference)) -
    as.numeric(determinant(reference, logarithm = TRUE)$modulus) +
    rho * sum(abs(reference[row(reference) != col(reference)]))
  gradient <- scatter - solve(reference)

  expect_true(fit$converged)
  expect_true(fit$objective_descent)
  expect_false(fit$backtracking_failed)
  expect_equal(fit$solution, reference, tolerance = 3e-7)
  expect_equal(fit$objective, objective, tolerance = 1e-10)
  expect_lte(max(diff(fit$objective_history)), 1e-12)
  expect_gt(fit$minimum_eigenvalue, 0)
  expect_lte(fit$relative_update, 1e-8)
  expect_lte(fit$kkt_diagonal, 1e-8)
  expect_lte(fit$kkt_off_diagonal, 1e-8)
  expect_equal(diag(solve(fit$solution)), diag(scatter), tolerance = 1e-8)
  expect_equal(
    gradient[1, 2] + rho * sign(reference[1, 2]), 0,
    tolerance = 1e-12
  )
})


test_that("zero-penalty glasso is an exact SPD inverse without repair", {
  scatter <- matrix(c(2, 0.25, 0.25, 1), 2)
  fit <- HDElliptical:::cpp_ch3teg_offdiag_glasso(
    scatter, 0, 1e-10, 10L, 1, 10L
  )
  expect_true(fit$converged)
  expect_equal(fit$iterations, 0L)
  expect_equal(fit$solution, solve(scatter), tolerance = 1e-13)
  expect_equal(fit$inverse, scatter, tolerance = 1e-13)
  expect_lte(fit$kkt_maximum, 1e-10)
  expect_error(
    HDElliptical:::cpp_ch3teg_offdiag_glasso(
      matrix(c(1, 1, 1, 1), 2), 0, 1e-8, 10L, 1, 10L
    ),
    "strictly positive definite.*no pseudoinverse or ridge"
  )
  expect_error(
    HDElliptical:::cpp_ch3teg_offdiag_glasso(
      matrix(c(0, 0, 0, 1), 2), 0.1, 1e-8, 10L, 1, 10L
    ),
    "not coercive.*no ridge"
  )
})


test_that("K equals one reduces to the vector spatial-sign formula", {
  dims <- 3L
  x <- teg_fixture_matrix(dims, half_n = 6L)
  data <- teg_as_observation_array(x, dims)
  fit <- tensor_spatial_sign_precision(
    data, lambda = 0, center = rep(0, dims), solver_tol = 1e-9
  )
  signs <- t(apply(x, 1, function(value) value / sqrt(sum(value^2))))
  scatter <- dims / nrow(x) * crossprod(signs)
  raw <- solve(scatter)

  expect_true(fit$valid)
  expect_identical(fit$dims, 3L)
  expect_equal(fit$pstar, 3L)
  expect_identical(unname(fit$pilot.branch), "inverse")
  expect_equal(fit$sign.matrix, signs, tolerance = 1e-14)
  expect_equal(fit$pilot.scatter[[1]], scatter, tolerance = 1e-13)
  expect_equal(fit$mode.scatter[[1]], scatter, tolerance = 1e-13)
  expect_equal(fit$raw.precision[[1]], raw, tolerance = 1e-8)
  expect_equal(fit$estimate[[1]], teg_frob_normalize(raw), tolerance = 1e-8)
  expect_equal(sqrt(sum(fit$estimate[[1]]^2)), 1, tolerance = 1e-13)
  expect_equal(fit$effective.penalty, c(mode1 = 0))
  expect_true(fit$solver[[1]]$certified)
})


test_that("two-mode public fit reproduces all paper formula blocks", {
  dims <- c(2L, 3L)
  x <- teg_fixture_matrix(dims, half_n = 7L)
  data <- teg_as_observation_array(x, dims)
  lambda <- c(0.03, 0.02)
  fit <- tensor_spatial_sign_precision(
    data, lambda, center = array(0, dims), keep_whitened = TRUE,
    solver_tol = 2e-7
  )

  signs <- t(apply(x, 1, function(value) value / sqrt(sum(value^2))))
  pilot.source <- lapply(seq_along(dims), function(k) {
    dims[k] / nrow(x) * Reduce(`+`, lapply(seq_len(nrow(x)), function(i) {
      tcrossprod(teg_unfold(signs[i, ], dims, k))
    }))
  })
  pilot.raw <- lapply(pilot.source, solve)
  pilot <- lapply(pilot.raw, teg_frob_normalize)

  expect_true(fit$valid)
  expect_identical(fit$dims, dims)
  expect_equal(dim(fit$signs), c(nrow(x), dims))
  expect_equal(fit$sign.matrix, signs, tolerance = 1e-13)
  expect_equal(fit$pilot.scatter, pilot.source, tolerance = 1e-12,
               ignore_attr = TRUE)
  expect_equal(fit$pilot.raw, pilot.raw, tolerance = 1e-10,
               ignore_attr = TRUE)
  expect_equal(fit$pilot, pilot, tolerance = 1e-10, ignore_attr = TRUE)
  expect_equal(vapply(fit$pilot, function(z) sqrt(sum(z^2)), numeric(1)),
               c(mode1 = 1, mode2 = 1), tolerance = 1e-13)
  expect_equal(vapply(fit$estimate, function(z) sqrt(sum(z^2)), numeric(1)),
               c(mode1 = 1, mode2 = 1), tolerance = 1e-13)
  expect_equal(unname(fit$effective.penalty), dims * lambda)
  expect_true(all(vapply(fit$solver, `[[`, logical(1), "certified")))
  expect_true(all(vapply(fit$solver, function(z) z$minimum_eigenvalue > 0,
                         logical(1))))
  expect_true(all(vapply(fit$solver, function(z) z$kkt_maximum <= 2e-7,
                         logical(1))))
  expect_true(all(vapply(fit$solver, function(z) z$relative_update <= 2e-7,
                         logical(1))))
  expect_equal(
    vapply(fit$solver, `[[`, numeric(1), "paper.objective"),
    vapply(seq_along(dims), function(k) {
      omega <- fit$raw.precision[[k]]
      sum(diag(fit$mode.scatter[[k]] %*% omega)) / dims[k] -
        as.numeric(determinant(omega, logarithm = TRUE)$modulus) / dims[k] +
        lambda[k] * sum(abs(omega[row(omega) != col(omega)]))
    }, numeric(1)),
    tolerance = 1e-10, ignore_attr = TRUE
  )
})


test_that("three-mode fixtures are valid and preserve the separable contract", {
  dims <- c(2L, 3L, 2L)
  x <- teg_fixture_matrix(dims, half_n = 8L)
  fit <- tensor_spatial_sign_precision(
    teg_as_tensor_list(x, dims), lambda = c(0.5, 0.4, 0.3),
    center = array(0, dims)
  )
  expect_true(fit$valid)
  expect_equal(fit$n, 16L)
  expect_identical(fit$dims, dims)
  expect_equal(fit$pstar, 12L)
  expect_length(fit$pilot, 3L)
  expect_length(fit$mode.scatter, 3L)
  expect_length(fit$raw.precision, 3L)
  expect_length(fit$estimate, 3L)
  expect_named(fit$estimate, c("mode1", "mode2", "mode3"))
  expect_true(all(vapply(fit$solver, `[[`, logical(1), "certified")))
  expect_equal(vapply(fit$pilot, function(z) sqrt(sum(z^2)), numeric(1)),
               c(mode1 = 1, mode2 = 1, mode3 = 1), tolerance = 1e-13)
  expect_equal(vapply(fit$estimate, function(z) sqrt(sum(z^2)), numeric(1)),
               c(mode1 = 1, mode2 = 1, mode3 = 1), tolerance = 1e-13)
})


test_that("list and observation-first array inputs are identical", {
  dims <- c(2L, 3L)
  x <- teg_fixture_matrix(dims, half_n = 6L)
  array.fit <- tensor_spatial_sign_precision(
    teg_as_observation_array(x, dims), 0.3, center = array(0, dims)
  )
  list.fit <- tensor_spatial_sign_precision(
    teg_as_tensor_list(x, dims), 0.3, center = array(0, dims)
  )
  expect_true(array.fit$valid)
  expect_true(list.fit$valid)
  expect_equal(array.fit$sign.matrix, list.fit$sign.matrix, tolerance = 0)
  expect_equal(array.fit$pilot, list.fit$pilot, tolerance = 0)
  expect_equal(array.fit$mode.scatter, list.fit$mode.scatter, tolerance = 0)
  expect_equal(array.fit$raw.precision, list.fit$raw.precision, tolerance = 0)
  expect_equal(array.fit$estimate, list.fit$estimate, tolerance = 0)
  expect_match(array.fit$input.type, "observation-first")
  expect_match(list.fit$input.type, "list")
})


test_that("mode permutation permutes the fitted mode list", {
  dims <- c(2L, 3L, 2L)
  x <- teg_fixture_matrix(dims, half_n = 7L)
  data <- teg_as_observation_array(x, dims)
  lambda <- c(0.4, 0.35, 0.3)
  original <- tensor_spatial_sign_precision(data, lambda, center = array(0, dims))
  mode.order <- c(2L, 3L, 1L)
  permuted.data <- aperm(data, c(1L, mode.order + 1L))
  permuted <- tensor_spatial_sign_precision(
    permuted.data, lambda[mode.order],
    center = array(0, dims[mode.order])
  )

  expect_true(original$valid)
  expect_true(permuted$valid)
  expect_identical(permuted$dims, dims[mode.order])
  for (k in seq_along(dims)) {
    expect_equal(permuted$pilot[[k]], original$pilot[[mode.order[k]]],
                 tolerance = 2e-10)
    expect_equal(permuted$mode.scatter[[k]],
                 original$mode.scatter[[mode.order[k]]], tolerance = 2e-10)
    expect_equal(permuted$estimate[[k]], original$estimate[[mode.order[k]]],
                 tolerance = 2e-8)
  }
})


test_that("within-mode signed permutations transform the corresponding graph", {
  dims <- c(2L, 3L)
  x <- teg_fixture_matrix(dims, half_n = 8L)
  p1 <- matrix(c(0, -1, 1, 0), 2)
  p2 <- diag(3)
  transformed.x <- teg_mode_transform(x, dims, list(p1, p2))
  original <- tensor_spatial_sign_precision(
    teg_as_observation_array(x, dims), c(0.03, 0.03),
    center = array(0, dims), solver_tol = 2e-7
  )
  transformed <- tensor_spatial_sign_precision(
    teg_as_observation_array(transformed.x, dims), c(0.03, 0.03),
    center = array(0, dims), solver_tol = 2e-7
  )

  expect_true(original$valid)
  expect_true(transformed$valid)
  expect_equal(transformed$pilot[[1]], p1 %*% original$pilot[[1]] %*% t(p1),
               tolerance = 2e-10)
  expect_equal(transformed$mode.scatter[[1]],
               p1 %*% original$mode.scatter[[1]] %*% t(p1),
               tolerance = 2e-10)
  expect_equal(transformed$raw.precision[[1]],
               p1 %*% original$raw.precision[[1]] %*% t(p1),
               tolerance = 2e-6)
  expect_equal(transformed$estimate[[1]],
               p1 %*% original$estimate[[1]] %*% t(p1),
               tolerance = 2e-6)
  expect_equal(transformed$pilot[[2]], original$pilot[[2]], tolerance = 2e-10)
  expect_equal(transformed$mode.scatter[[2]], original$mode.scatter[[2]],
               tolerance = 2e-10)
  expect_equal(transformed$estimate[[2]], original$estimate[[2]],
               tolerance = 2e-6)
})


test_that("translation, global scale, sign, and row order invariances hold", {
  dims <- c(2L, 2L)
  x <- teg_fixture_matrix(dims, half_n = 7L)
  data <- teg_as_observation_array(x, dims)
  shift <- c(3, -5, 7, 11)
  translated <- teg_as_observation_array(
    sweep(x, 2, shift, `+`), dims
  )
  permutation <- c(8:14, 1:7)
  fits <- list(
    base = tensor_spatial_sign_precision(data, 0.4),
    translated = tensor_spatial_sign_precision(translated, 0.4),
    large = tensor_spatial_sign_precision(data * 1e150, 0.4),
    small = tensor_spatial_sign_precision(data * 1e-150, 0.4),
    negative = tensor_spatial_sign_precision(-data, 0.4),
    reordered = tensor_spatial_sign_precision(
      teg_as_observation_array(x[permutation, , drop = FALSE], dims), 0.4
    )
  )

  expect_true(all(vapply(fits, `[[`, logical(1), "valid")))
  expect_equal(as.numeric(fits$base$center), rep(0, 4), tolerance = 5e-8)
  expect_equal(as.numeric(fits$translated$center), shift, tolerance = 5e-8)
  expect_equal(fits$translated$estimate, fits$base$estimate, tolerance = 2e-9)
  expect_equal(fits$large$estimate, fits$base$estimate, tolerance = 2e-9)
  expect_equal(fits$small$estimate, fits$base$estimate, tolerance = 2e-9)
  expect_equal(fits$negative$estimate, fits$base$estimate, tolerance = 2e-9)
  expect_equal(fits$reordered$estimate, fits$base$estimate, tolerance = 2e-9)
  expect_equal(fits$large$sign.matrix, fits$base$sign.matrix, tolerance = 2e-13)
  expect_equal(fits$small$sign.matrix, fits$base$sign.matrix, tolerance = 2e-13)
  expect_equal(fits$negative$sign.matrix, -fits$base$sign.matrix,
               tolerance = 5e-8)
  expect_equal(fits$reordered$sign.matrix,
               fits$base$sign.matrix[permutation, , drop = FALSE],
               tolerance = 2e-13)
})


test_that("the exact paper sample-size branch including equality is locked", {
  dims <- c(4L, 3L)
  x <- matrix(
    c(sin(seq_len(24)) + cos(seq_len(24) / 3)),
    nrow = 2
  )
  fit <- tensor_spatial_sign_precision(
    teg_as_observation_array(x, dims), lambda = 2,
    center = array(0, dims)
  )
  expect_true(fit$valid)
  expect_equal(2 * prod(dims), dims[1]^2 * (dims[1] - 1) / 2)
  expect_identical(unname(fit$pilot.branch[1]), "identity")
  expect_identical(unname(fit$pilot.branch[2]), "inverse")
  expect_equal(fit$pilot.raw[[1]], diag(4), tolerance = 0)
  expect_equal(fit$pilot[[1]], diag(4) / 2, tolerance = 0)
  expect_equal(fit$pilot.sqrt[[1]], diag(4) / sqrt(2), tolerance = 1e-14)
  expect_equal(fit$pilot.raw[[2]], solve(fit$pilot.scatter[[2]]),
               tolerance = 1e-10)
})


test_that("exact zero residuals are mapped to zero and diagnosed", {
  dims <- c(2L, 2L)
  x <- rbind(teg_fixture_matrix(dims, half_n = 6L), rep(0, 4))
  fit <- tensor_spatial_sign_precision(
    teg_as_observation_array(x, dims), lambda = 0.5,
    center = array(0, dims)
  )
  expect_true(fit$valid)
  expect_equal(fit$diagnostics$zero.residuals$count, 1L)
  expect_equal(fit$diagnostics$zero.residuals$indices, nrow(x))
  expect_equal(fit$sign.matrix[nrow(x), ], rep(0, 4), tolerance = 0)
  expect_equal(sum(fit$sign.matrix[nrow(x), ]^2), 0)
  expect_equal(
    rowSums(fit$sign.matrix[-nrow(x), , drop = FALSE]^2),
    rep(1, nrow(x) - 1L), tolerance = 1e-14
  )

  identical.data <- array(3, dim = c(6, 2, 2))
  expect_warning(
    all.zero <- tensor_spatial_sign_precision(
      identical.data, lambda = 0.5, strict = FALSE
    ),
    "pilot failed"
  )
  expect_false(all.zero$valid)
  expect_null(all.zero$estimate)
  expect_equal(all.zero$diagnostics$zero.residuals$count, 6L)
  expect_equal(all.zero$sign.matrix, matrix(0, 6, 4), tolerance = 0)
})


test_that("singular pilots and failed solvers never produce pseudo-estimates", {
  collinear <- cbind(-3:3, 2 * (-3:3))
  expect_error(
    tensor_spatial_sign_precision(collinear, 0.2, center = c(0, 0)),
    "pilot failed.*No inverse or SPD repair"
  )
  expect_warning(
    invalid.pilot <- tensor_spatial_sign_precision(
      collinear, 0.2, center = c(0, 0), strict = FALSE
    ),
    "pilot failed"
  )
  expect_false(invalid.pilot$valid)
  expect_null(invalid.pilot$estimate)
  expect_match(invalid.pilot$diagnostics$no.repair, "No ridge")

  dims <- c(2L, 3L)
  x <- teg_fixture_matrix(dims, half_n = 7L)
  expect_warning(
    invalid.solver <- tensor_spatial_sign_precision(
      teg_as_observation_array(x, dims), lambda = 1e-4,
      center = array(0, dims), solver_max_iter = 1L, strict = FALSE
    ),
    "did not satisfy every SPD.*certificate"
  )
  expect_false(invalid.solver$valid)
  expect_null(invalid.solver$estimate)
  expect_match(invalid.solver$failure.stage, "solver certificate")
  expect_true(any(vapply(invalid.solver$solver, function(z) {
    !is.null(z) && !isTRUE(z$certified)
  }, logical(1))))
})


test_that("thresholding uses absolute greater-than-or-equal and keeps diagonals", {
  omega1 <- matrix(c(0.1, -0.5, -0.5, 0.2), 2)
  omega2 <- matrix(c(0.01, 0.49, 0.49, 0.02), 2)
  fit <- structure(
    list(valid = TRUE, estimate = list(mode1 = omega1, mode2 = omega2)),
    class = "tensor_spatial_sign_precision_fit"
  )
  got <- threshold_tensor_spatial_sign_precision(fit, c(0.5, 0.5))
  expect_true(got$valid)
  expect_equal(got$estimate[[1]][1, 2], -0.5)
  expect_equal(got$estimate[[1]][2, 1], -0.5)
  expect_equal(diag(got$estimate[[1]]), c(0.1, 0.2))
  expect_equal(got$estimate[[2]][1, 2], 0)
  expect_equal(got$estimate[[2]][2, 1], 0)
  expect_equal(diag(got$estimate[[2]]), c(0.01, 0.02))
  expect_equal(got$tau, c(mode1 = 0.5, mode2 = 0.5))
  expect_match(got$comparison, "absolute value >= tau")
})


test_that("input and control contracts reject malformed requests", {
  good <- teg_as_observation_array(teg_fixture_matrix(c(2L, 2L), 4L),
                                   c(2L, 2L))
  expect_error(tensor_spatial_sign_precision(1:10, 0.1),
               "list.*observation-first")
  expect_error(tensor_spatial_sign_precision(list(array(1:4, c(2, 2))), 0.1),
               "at least two")
  expect_error(tensor_spatial_sign_precision(
    list(array(1:4, c(2, 2)), array(1:6, c(2, 3))), 0.1
  ), "identical dimensions")
  expect_error(tensor_spatial_sign_precision(
    list(array(1:4, c(2, 2)), matrix(letters[1:4], 2)), 0.1
  ), "numeric.*identical")
  expect_error(tensor_spatial_sign_precision(array(1:4, c(1, 2, 2)), 0.1),
               "at least two observations")
  bad.na <- good
  bad.na[1] <- NA_real_
  expect_error(tensor_spatial_sign_precision(bad.na, 0.1), "must not contain")
  bad.inf <- good
  bad.inf[1] <- Inf
  expect_error(tensor_spatial_sign_precision(bad.inf, 0.1), "must not contain")
  expect_error(tensor_spatial_sign_precision(good, -0.1), "non-negative")
  expect_error(tensor_spatial_sign_precision(good, c(0.1, 0.2, 0.3)),
               "one-per-mode")
  expect_error(tensor_spatial_sign_precision(good, Inf), "non-negative")
  expect_error(tensor_spatial_sign_precision(good, 1e308),
               "rescaled penalties.*finite")
  expect_error(tensor_spatial_sign_precision(good, 0.1, center = 1:3),
               "one entry per tensor coordinate")
  expect_error(tensor_spatial_sign_precision(
    good, 0.1, center = array(0, c(4, 1))
  ), "mode dimensions")
  expect_error(tensor_spatial_sign_precision(
    good, 0.1, center = 1:3, strict = FALSE
  ), "one entry per tensor coordinate")
  expect_error(tensor_spatial_sign_precision(good, 0.1, solver_tol = 0),
               "strictly positive")
  expect_error(tensor_spatial_sign_precision(good, 0.1, solver_max_iter = 1.5),
               "positive integer")
  expect_error(tensor_spatial_sign_precision(good, 0.1, initial_step = NA),
               "strictly positive")
  expect_error(tensor_spatial_sign_precision(good, 0.1, max_backtracking = 0),
               "positive integer")
  expect_error(tensor_spatial_sign_precision(good, 0.1, keep_whitened = NA),
               "TRUE or FALSE")
  expect_error(tensor_spatial_sign_precision(good, 0.1, strict = 1),
               "TRUE or FALSE")
  expect_error(threshold_tensor_spatial_sign_precision(list(), 0.1),
               "valid tensor")
  fake <- structure(list(valid = TRUE, estimate = list(diag(2))),
                    class = "tensor_spatial_sign_precision_fit")
  expect_error(threshold_tensor_spatial_sign_precision(fake, -1),
               "non-negative")
  expect_error(threshold_tensor_spatial_sign_precision(fake, c(1, 2)),
               "one-per-mode")
})


test_that("native kernels reject dimension and SPD contract violations", {
  expect_error(
    HDElliptical:::cpp_ch3teg_mode_crossproducts(matrix(1:6, 1), c(2L, 2L)),
    "product of `dims`"
  )
  expect_error(
    HDElliptical:::cpp_ch3teg_whitened_mode_scatter(
      matrix(1:4, 1), c(2L, 2L), list(diag(2)), FALSE
    ),
    "one matrix per tensor mode"
  )
  expect_error(
    HDElliptical:::cpp_ch3teg_whitened_mode_scatter(
      matrix(1:4, 1), c(2L, 2L),
      list(matrix(c(1, 2, 0, 1), 2), diag(2)), FALSE
    ),
    "must be symmetric"
  )
  expect_error(
    HDElliptical:::cpp_ch3teg_whitened_mode_scatter(
      matrix(1:4, 1), c(2L, 2L),
      list(diag(c(1, 0)), diag(2)), FALSE
    ),
    "strictly positive definite.*no eigenvalue repair"
  )
  expect_error(
    HDElliptical:::cpp_ch3teg_offdiag_glasso(
      matrix(c(1, 0.2, 0.1, 1), 2), 0.1, 1e-8, 10L, 1, 10L
    ),
    "must be symmetric"
  )
  expect_error(
    HDElliptical:::cpp_ch3teg_offdiag_glasso(
      diag(2), -0.1, 1e-8, 10L, 1, 10L
    ),
    "Invalid.*controls"
  )
})
