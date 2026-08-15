# Deterministic formula, boundary, invariance, and certificate tests for the
# Chapter 7 spatial-clustering methods. These are method tests, not paper
# simulation reproductions.

ch7sc_sparse_fixture <- function() {
  x <- rbind(
    c(-4, 0.0,  2), c(-3, 1.0,  2), c(-2, 2.0,  2),
    c( 2, 0.2, -2), c( 3, 1.2, -2), c( 4, 2.2, -2)
  )
  colnames(x) <- c("signal1", "weak", "signal2")
  x
}


ch7sc_permutations <- function(n, p) {
  stopifnot(n == 6L, p == 3L)
  list(
    cbind(n:1L, c(3L:n, 1:2), c(2L:n, 1L)),
    cbind(c(2L:n, 1L), n:1L, c(4L:n, 1:3))
  )
}


ch7sc_manual_tau_objective <- function(x, fit, tol = 1e-6,
                                        max_iter = 2000L) {
  active <- fit$active_set
  overall <- spatial_median(
    x[, active, drop = FALSE], tol = tol, max_iter = max_iter, warn = FALSE
  )
  stopifnot(isTRUE(attr(overall, "converged")))
  sizes <- tabulate(fit$labels, nrow(fit$algorithm_centers))
  difference <- sweep(
    fit$algorithm_centers[, active, drop = FALSE], 2L, overall, "-"
  )
  sum(sizes * rowSums(difference^2))
}


test_that("Chapter 7 spatial public APIs expose all method decisions", {
  expect_identical(
    names(formals(k_spatial_median)),
    c("x", "K", "init", "first_index", "tol", "max_iter",
      "spatial_max_iter", "zero_tol", "empty_action", "cycle_action",
      "ties", "keep_distances")
  )
  expect_identical(
    names(formals(sm_sscm)),
    c("x", "K", "lambda", "init", "first_index", "tol", "max_iter",
      "spatial_max_iter", "zero_tol", "empty_action", "cycle_action",
      "ties", "keep_distances", "keep_signs")
  )
  expect_identical(
    names(formals(sparse_k_spatial_median)),
    c("x", "K", "tau", "init", "first_index", "tol", "max_iter",
      "spatial_max_iter", "zero_tol", "empty_action", "empty_active",
      "cycle_action", "ties", "reset_excluded", "keep_distances")
  )
  expect_true(all(
    c("tau_grid", "B", "permutations", "seed", "selection_ties") %in%
      names(formals(sparse_sm_select_tau))
  ))
  expect_true(all(
    c("k_grid", "tau_grid", "B", "permutations", "seed",
      "selection_ties") %in% names(formals(sparse_sm_select_k))
  ))
})


test_that("K-spatial-median matches the unsquared-distance formula", {
  x <- matrix(c(0, 1, 2, 8, 9, 10), ncol = 1)
  fit <- k_spatial_median(x, 2, init = c(1, 6), keep_distances = TRUE)

  expect_s3_class(fit, "k_spatial_median_fit")
  expect_identical(fit$labels, c(1L, 1L, 1L, 2L, 2L, 2L))
  expect_equal(unname(fit$centers[, 1]), c(1, 9), tolerance = 0)
  manual <- sum(abs(x[, 1] - fit$centers[fit$labels, 1]))
  expect_equal(fit$objective, manual, tolerance = 0)
  expect_equal(fit$objective, 4, tolerance = 0)
  expect_equal(
    unname(fit$distances),
    unname(abs(outer(x[, 1], fit$centers[, 1], "-"))),
    tolerance = 1e-14
  )
  expect_true(fit$diagnostics$converged)
  expect_true(fit$diagnostics$objective.monotone)
  expect_true(all(diff(fit$diagnostics$objective.history) <= 1e-12))
  expect_true(all(vapply(
    fit$diagnostics$center.certificates, `[[`, logical(1), "converged"
  )))

  one <- k_spatial_median(x, 1, init = 1)
  expect_identical(one$labels, rep.int(1L, nrow(x)))
  expect_equal(one$objective, sum(abs(x[, 1] - one$centers[1, 1])))
})


test_that("K-spatial-median assignment, max-min ties, and repair are explicit", {
  tie <- HDElliptical:::cpp_ch7sc_assign_euclidean(
    matrix(0, 1, 1), matrix(c(-1, 1), 2, 1), 0L, TRUE
  )
  expect_identical(tie$labels, 1L)
  expect_equal(as.numeric(tie$assigned_distance), 1, tolerance = 0)
  expect_identical(tie$tie_count, 1L)

  maxmin <- k_spatial_median(
    matrix(c(0, -1, 1), ncol = 1), 2,
    init = "maxmin", first_index = 1
  )
  expect_identical(maxmin$diagnostics$initialization.indices, c(1L, 2L))

  x <- matrix(0:2, ncol = 1)
  duplicated.centers <- matrix(0, 2, 1)
  expect_error(
    k_spatial_median(x, 2, init = duplicated.centers),
    "empty cluster.*no repair"
  )
  repaired <- k_spatial_median(
    x, 2, init = duplicated.centers, empty_action = "farthest"
  )
  expect_identical(repaired$labels, c(1L, 1L, 2L))
  expect_identical(repaired$diagnostics$repairs[[1]]$row, 3L)
  expect_identical(repaired$diagnostics$repairs[[1]]$from, 1L)
  expect_identical(repaired$diagnostics$repairs[[1]]$to, 2L)

  expect_error(
    k_spatial_median(matrix(c(0, 0, 1), ncol = 1), 3),
    "at least `K` distinct rows"
  )
})


test_that("K-spatial-median is translation, scale, and rotation equivariant", {
  x <- rbind(c(-4, 0), c(-3, 0), c(-2, 0),
             c( 2, 0), c( 3, 0), c( 4, 0))
  base <- k_spatial_median(x, 2, init = c(1, 6))
  angle <- 0.37
  rotation <- matrix(c(cos(angle), -sin(angle),
                       sin(angle),  cos(angle)), 2, 2, byrow = TRUE)
  shift <- c(7, -4)
  transformed.x <- sweep(x %*% rotation, 2L, shift, "+")
  transformed <- k_spatial_median(transformed.x, 2, init = c(1, 6))

  expect_identical(transformed$labels, base$labels)
  expect_equal(
    unname(transformed$centers),
    unname(sweep(base$centers %*% rotation, 2L, shift, "+")),
    tolerance = 1e-12
  )
  expect_equal(transformed$objective, base$objective, tolerance = 1e-12)

  scaled <- k_spatial_median(3 * x, 2, init = c(1, 6))
  expect_identical(scaled$labels, base$labels)
  expect_equal(unname(scaled$centers), 3 * unname(base$centers),
               tolerance = 1e-12)
  expect_equal(scaled$objective, 3 * base$objective, tolerance = 1e-12)
})


test_that("SSCM kernel implements U(0)=0, ridge trace, and exact inverse", {
  x <- rbind(c(0, 0), c(1, 0), c(2, 0), c(2, 0))
  centers <- rbind(c(0, 0), c(2, 0))
  metric <- HDElliptical:::cpp_ch7sc_sscm_metric(
    x, centers, c(1L, 1L, 2L, 2L), 0.1, 0, TRUE
  )

  expect_identical(metric$n_zero, 3L)
  expect_equal(metric$signs[1, ], c(0, 0), tolerance = 0)
  expect_equal(metric$signs[2, ], c(1, 0), tolerance = 0)
  expect_equal(metric$metric, diag(c(0.35, 0.1)), tolerance = 1e-15)
  expect_equal(metric$trace_expected, 3 / 4 - 2 / 4 + 0.2,
               tolerance = 1e-15)
  expect_equal(sum(diag(metric$metric)), metric$trace_expected,
               tolerance = 1e-15)
  expect_lte(metric$trace_error, 1e-15)
  expect_gte(metric$minimum_eigenvalue, 0.1 - 1e-15)
  expect_lte(metric$inverse_residual, 1e-14)
  expect_equal(metric$metric %*% metric$inverse, diag(2), tolerance = 1e-14)
  expect_match(metric$inverse_method, "no pseudoinverse")

  dropped <- HDElliptical:::cpp_ch7sc_sscm_metric(
    x, centers, c(1L, 1L, 2L, 2L), 0.1, 0, FALSE
  )
  expect_null(dropped$signs)
})


test_that("SSCM geometry is orthogonally equivariant and works for p over n", {
  x <- rbind(c(-1, 0), c(-2, 0), c(1, 2), c(2, 2))
  centers <- rbind(c(0, 0), c(1, 1))
  labels <- c(1L, 1L, 2L, 2L)
  angle <- 0.61
  q <- matrix(c(cos(angle), -sin(angle),
                sin(angle),  cos(angle)), 2, 2, byrow = TRUE)
  original <- HDElliptical:::cpp_ch7sc_sscm_metric(
    x, centers, labels, 0.2, 0, FALSE
  )
  rotated <- HDElliptical:::cpp_ch7sc_sscm_metric(
    x %*% q, centers %*% q, labels, 0.2, 0, FALSE
  )
  expect_equal(rotated$metric, t(q) %*% original$metric %*% q,
               tolerance = 1e-14)
  expect_equal(rotated$inverse, t(q) %*% original$inverse %*% q,
               tolerance = 1e-13)

  high.dimensional <- cbind(c(-2, -1, 1, 2), matrix(0, 4, 5))
  fit <- sm_sscm(
    high.dimensional, 2, lambda = 0.2, init = c(1, 4),
    keep_signs = TRUE, keep_distances = TRUE
  )
  expect_s3_class(fit, "sm_sscm_fit")
  expect_identical(fit$labels, c(1L, 1L, 2L, 2L))
  expect_gt(min(eigen(fit$metric, symmetric = TRUE, only.values = TRUE)$values),
            0)
  expect_lte(fit$diagnostics$inverse.residual, 1e-13)
  expect_match(fit$diagnostics$global.objective.monotonicity,
               "not claimed")
  expect_identical(dim(fit$signs), dim(high.dimensional))
  expect_identical(dim(fit$distances), c(4L, 2L))
})


test_that("Sparse--SM scores and active-boundary rules are exact", {
  centers <- rbind(c(-2, 1, 3), c(2, 1, -1))
  score <- HDElliptical:::cpp_ch7sc_feature_scores(centers)
  expect_equal(as.numeric(score$average_center), c(0, 1, 1), tolerance = 0)
  expect_equal(as.numeric(score$scores), c(4, 0, 4), tolerance = 0)

  boundary <- HDElliptical:::.c7sc_active(c(4, 0, 4), 4, "error")
  expect_identical(boundary$active, c(1L, 3L))
  expect_error(
    HDElliptical:::.c7sc_active(c(4, 0, 4), 5, "error"),
    "empty active set.*no fallback"
  )
  largest <- HDElliptical:::.c7sc_active(c(4, 0, 4), 5, "largest")
  expect_identical(largest$active, 1L)
  all.features <- HDElliptical:::.c7sc_active(c(4, 0, 4), 5, "all")
  expect_identical(all.features$active, 1:3)

  x <- ch7sc_sparse_fixture()
  fit <- sparse_k_spatial_median(x, 2, tau = 4, init = c(1, 6))
  expect_identical(fit$labels, c(1L, 1L, 1L, 2L, 2L, 2L))
  expect_equal(fit$scores, c(6, 0.2, 4), tolerance = 1e-14)
  expect_identical(fit$active_set, c(1L, 3L))
  expect_match(fit$diagnostics$active.boundary.rule, ">=")
  expect_match(fit$diagnostics$global.objective.monotonicity,
               "not claimed")

  expect_error(
    sparse_k_spatial_median(x, 2, tau = 7, init = c(1, 6)),
    "empty active set.*no fallback"
  )
  fallback <- sparse_k_spatial_median(
    x, 2, tau = 7, init = c(1, 6), empty_active = "largest"
  )
  expect_identical(fallback$active_set, 1L)
  expect_match(
    fallback$diagnostics$empty.active.fallback.history[1], "largest score"
  )
  keep.all <- sparse_k_spatial_median(
    x, 2, tau = 7, init = c(1, 6), empty_active = "all"
  )
  expect_identical(keep.all$active_set, 1:3)
})


test_that("Sparse--SM ignores inactive coordinates and resets output only", {
  x <- ch7sc_sparse_fixture()
  none <- sparse_k_spatial_median(
    x, 2, tau = 4, init = c(1, 6), reset_excluded = "none"
  )
  reset <- sparse_k_spatial_median(
    x, 2, tau = 4, init = c(1, 6),
    reset_excluded = "overall_spatial_median"
  )

  expect_identical(reset$labels, none$labels)
  expect_identical(reset$active_set, none$active_set)
  expect_equal(reset$algorithm_centers, none$algorithm_centers, tolerance = 0)
  expect_equal(
    reset$centers[, reset$active_set, drop = FALSE],
    none$centers[, none$active_set, drop = FALSE], tolerance = 0
  )
  expect_equal(
    unname(reset$centers[, 2]),
    rep(reset$diagnostics$reset.baseline[2], 2), tolerance = 0
  )
  expect_equal(reset$diagnostics$reset.baseline[2], 1.1, tolerance = 1e-7)
  expect_false(isTRUE(all.equal(reset$centers[, 2], none$centers[, 2])))
  expect_true(reset$diagnostics$reset.output.only)

  centers <- rbind(x[1, ], x[6, ])
  base <- HDElliptical:::cpp_ch7sc_assign_euclidean(
    x, centers, c(0L, 2L), TRUE
  )
  changed.x <- x
  changed.centers <- centers
  changed.x[, 2] <- c(-1e6, 2e6, -3e6, 4e6, -5e6, 6e6)
  changed.centers[, 2] <- c(9e8, -9e8)
  changed <- HDElliptical:::cpp_ch7sc_assign_euclidean(
    changed.x, changed.centers, c(0L, 2L), TRUE
  )
  expect_identical(changed$labels, base$labels)
  expect_equal(changed$distance_matrix, base$distance_matrix, tolerance = 0)
})


test_that("Sparse--SM final assignment honors explicit farthest repair", {
  x <- structure(
    c(3L, 1L, -1L, 1L, -1L, 2L, 3L, 0L,
      -2L, 2L, 2L, -2L, -1L, -2L, 3L, -1L),
    dim = c(8L, 2L)
  )
  fit <- sparse_k_spatial_median(
    x, 3L, tau = 2.1435062466189265, init = c(1L, 6L, 8L),
    max_iter = 1L, empty_action = "farthest", empty_active = "largest",
    cycle_action = "return"
  )
  expect_false(fit$diagnostics$converged)
  expect_identical(fit$diagnostics$termination, "max_iter")
  expect_true(length(fit$diagnostics$next.assignment.repairs) > 0L)
})


test_that("Sparse--SM cycle return exposes a self-consistent final state", {
  x <- structure(
    c(4L, -4L, -2L, 2L, 4L, -1L, 2L, -3L,
      -2L, 1L, 4L, 4L, 2L, -2L, -1L, 4L,
      0L, 1L, 3L, -3L, 3L, 2L, -2L, -4L),
    dim = c(8L, 3L)
  )
  tau <- 6.0433367919176817
  args <- list(
    x = x, K = 3L, tau = tau, init = c(1L, 5L, 6L),
    max_iter = 40L, empty_action = "farthest",
    empty_active = "largest"
  )
  fit <- do.call(
    sparse_k_spatial_median, c(args, list(cycle_action = "return"))
  )
  expect_false(fit$diagnostics$converged)
  expect_identical(fit$diagnostics$termination, "cycle")
  expect_identical(fit$diagnostics$cycle.length, 6L)
  expected.centers <- t(vapply(seq_len(3L), function(k) {
    as.numeric(spatial_median(
      x[fit$labels == k, , drop = FALSE], warn = FALSE
    ))
  }, numeric(ncol(x))))
  expect_equal(
    unname(fit$algorithm_centers), unname(expected.centers), tolerance = 1e-8
  )
  expected.scores <- as.numeric(
    HDElliptical:::cpp_ch7sc_feature_scores(expected.centers)$scores
  )
  expect_equal(fit$scores, expected.scores, tolerance = 1e-12)
  expect_identical(fit$active_set, which(expected.scores >= tau))
  expect_error(
    do.call(
      sparse_k_spatial_median, c(args, list(cycle_action = "error"))
    ),
    "cycle length 6"
  )
})


test_that("Selector geometry matches the stated retained-subspace formulas", {
  x <- matrix(c(0, 1, 3, 4), ncol = 1)
  centers <- matrix(c(1, 3), ncol = 1)
  geometry <- HDElliptical:::cpp_ch7sc_geometry(
    x, centers, c(1L, 1L, 2L, 2L), 2, 0L
  )
  expect_equal(as.integer(geometry$sizes), c(2L, 2L))
  expect_equal(geometry$within_sum, 2, tolerance = 0)
  expect_equal(geometry$within_mean, 0.5, tolerance = 0)
  expect_equal(geometry$between_dispersion, 4, tolerance = 0)
  expect_equal(geometry$average_between, 2, tolerance = 0)
  expect_identical(geometry$pair_count, 1L)
})


test_that("Sparse--SM tau selector computes the exact permutation Gap", {
  x <- ch7sc_sparse_fixture()
  permutations <- ch7sc_permutations(nrow(x), ncol(x))
  selection <- sparse_sm_select_tau(
    x, 2, tau_grid = c(1, 0, 1), permutations = permutations,
    init = c(1, 6), tol = 1e-6, spatial_max_iter = 2000,
    empty_active = "largest", keep_fits = TRUE,
    keep_permutations = TRUE
  )

  expect_s3_class(selection, "sparse_sm_tau_selection")
  expect_identical(selection$table$tau, c(0, 1))
  expected.gap <- log(selection$table$objective) -
    colMeans(log(selection$reference_objectives))
  expect_equal(selection$table$gap, expected.gap, tolerance = 1e-14)
  expect_identical(selection$selected_index, which.max(expected.gap))
  expect_identical(
    selection$selected_tau, selection$table$tau[selection$selected_index]
  )
  manual <- vapply(
    selection$fits, ch7sc_manual_tau_objective, numeric(1),
    x = x, tol = 1e-6, max_iter = 2000L
  )
  expect_equal(selection$table$objective, manual, tolerance = 1e-12)
  expect_identical(selection$permutations, permutations)
  expect_identical(selection$diagnostics$selection.tie.rule, "smallest tau")
  expect_identical(
    selection$diagnostics$coordinate.convention, "retained active subspace"
  )
})


test_that("Permutation generation is seeded, local, and validated", {
  set.seed(918)
  before <- .Random.seed
  first <- HDElliptical:::.c7sc_permutation_plan(6L, 3L, 2, NULL, 41)
  expect_identical(.Random.seed, before)
  second <- HDElliptical:::.c7sc_permutation_plan(6L, 3L, 2, NULL, 41)
  expect_identical(first, second)
  expect_identical(.Random.seed, before)
  expect_identical(first$seed, 41L)

  singleton.array <- array(
    c(seq_len(6L), 6:1), dim = c(6L, 1L, 2L)
  )
  singleton <- HDElliptical:::.c7sc_permutation_plan(
    6L, 1L, 2L, singleton.array, NULL
  )
  expect_identical(lapply(singleton$plan, dim), rep(list(c(6L, 1L)), 2L))
  expect_identical(singleton$plan[[1L]], matrix(seq_len(6L), ncol = 1L))
  expect_identical(singleton$plan[[2L]], matrix(6:1, ncol = 1L))

  valid <- ch7sc_permutations(6, 3)
  expect_error(
    HDElliptical:::.c7sc_permutation_plan(6L, 3L, 2, valid, 41),
    "either `permutations` or `seed`"
  )
  invalid <- valid
  invalid[[1]][1, 1] <- invalid[[1]][2, 1]
  expect_error(
    HDElliptical:::.c7sc_permutation_plan(6L, 3L, 2, invalid, NULL),
    "must permute 1:n"
  )
})


test_that("Sparse--SM K selector retunes tau and computes BWDM exactly", {
  x <- rbind(
    c(-7, -0.2), c(-6, 0.1), c(-5, 0.3),
    c(-1, -0.3), c( 0, 0.1), c( 1, 0.4),
    c( 5, -0.2), c( 6, 0.2), c( 7, 0.5)
  )
  n <- nrow(x)
  permutations <- list(cbind(c(2L:n, 1L), c(4L:n, 1:3)))
  selection <- sparse_sm_select_k(
    x, k_grid = c(3, 2, 3), tau_grid = 0,
    permutations = permutations, init = "maxmin", tol = 1e-6,
    spatial_max_iter = 2000, empty_action = "farthest",
    keep_selections = TRUE, keep_permutations = TRUE
  )

  expect_s3_class(selection, "sparse_sm_k_selection")
  expect_identical(selection$table$K, c(2L, 3L))
  expected <- with(
    selection$table,
    (ABDM / (K - 1)) / (AWDM / (n - K))
  )
  expect_equal(selection$table$BWDM, expected, tolerance = 1e-14)
  expect_identical(selection$selected_index, which.max(expected))
  expect_identical(
    selection$selected_k, selection$table$K[selection$selected_index]
  )
  expect_true(selection$diagnostics$tau.reselected.for.each.K)
  expect_identical(selection$diagnostics$selection.tie.rule, "smallest K")
  expect_identical(selection$permutations, permutations)
  expect_length(selection$selections, 2L)
  expect_true(all(vapply(
    selection$selections,
    function(value) inherits(value, "sparse_sm_tau_selection"), logical(1)
  )))
})


test_that("Selectors reject zero-log and zero-within degeneracies", {
  identical.x <- matrix(0, 4, 2)
  identity.permutation <- list(cbind(1:4, 1:4))
  expect_error(
    sparse_sm_select_tau(
      identical.x, 2, tau_grid = 0,
      permutations = identity.permutation,
      init = matrix(0, 2, 2), empty_action = "farthest"
    ),
    "strictly positive before taking its logarithm"
  )
  expect_error(
    sparse_sm_select_tau(
      ch7sc_sparse_fixture(), 1, tau_grid = 0,
      permutations = ch7sc_permutations(6, 3)
    ),
    "between 2 and"
  )

  zero.within <- rbind(
    c(-1, 0), c(-1, 0), c(0, 0),
    c(0, 0), c(1, 0), c(1, 0)
  )
  expect_error(
    sparse_sm_select_k(
      zero.within, 3, tau_grid = 0,
      permutations = list(cbind(1:6, 1:6)), init = "maxmin"
    ),
    "AWDM is not strictly positive"
  )
})


test_that("Spatial distance kernels preserve finite extreme geometry", {
  left <- matrix(c(1e100 + 2e84, 1e300), 1L, 2L)
  right <- matrix(c(1e100, 1e300), 1L, 2L)
  delta <- left[1L, 1L] - right[1L, 1L]
  expect_gt(abs(delta), 0)

  euclidean <- HDElliptical:::cpp_ch7sc_assign_euclidean(
    left, right, 0:1, TRUE
  )
  expect_equal(
    as.numeric(euclidean$assigned_distance) / abs(delta), 1,
    tolerance = 1e-14
  )
  metric <- HDElliptical:::cpp_ch7sc_assign_metric(
    left, right, diag(2), TRUE
  )
  expect_equal(
    as.numeric(metric$assigned_distance) / (delta * delta), 1,
    tolerance = 1e-14
  )
  sign.fit <- HDElliptical:::cpp_ch7sc_sscm_metric(
    left, right, 1L, 0.1, 0.8 * abs(delta), TRUE
  )
  expect_identical(as.integer(sign.fit$n_zero), 0L)
  expect_equal(drop(sign.fit$signs), c(1, 0), tolerance = 0)

  accumulated <- HDElliptical:::cpp_ch7sc_assign_metric(
    matrix(c(1e-200, 1e-200), 1L, 2L), matrix(0, 1L, 2L),
    diag(c(1e308, 1e308)), TRUE
  )
  expect_gt(as.numeric(accumulated$assigned_distance), 0)
  expect_equal(
    as.numeric(accumulated$assigned_distance) / 2e-92, 1,
    tolerance = 1e-14
  )

  maximum <- .Machine$double.xmax
  opposite.x <- matrix(c(maximum, 0), 1L, 2L)
  opposite.center <- matrix(c(-maximum, 0), 1L, 2L)
  expected.root <- maximum * (2 * sqrt(1e-309))
  expect_true(is.finite(expected.root^2))
  opposite.metric <- HDElliptical:::cpp_ch7sc_assign_metric(
    opposite.x, opposite.center, diag(c(1e-309, 1)), TRUE
  )
  expect_equal(
    as.numeric(opposite.metric$assigned_distance) / expected.root^2, 1,
    tolerance = 1e-13
  )
  opposite.sign <- HDElliptical:::cpp_ch7sc_sscm_metric(
    opposite.x, opposite.center, 1L, 0.1, 0, TRUE
  )
  expect_equal(drop(opposite.sign$signs), c(1, 0), tolerance = 0)

  general.metric <- matrix(c(2, 0.4, 0.4, 1), 2L, 2L)
  general.difference <- c(3, -2)
  general <- HDElliptical:::cpp_ch7sc_assign_metric(
    matrix(general.difference, 1L), matrix(0, 1L, 2L),
    general.metric, TRUE
  )
  expected <- drop(t(general.difference) %*%
                     general.metric %*% general.difference)
  expect_equal(as.numeric(general$assigned_distance), expected,
               tolerance = 1e-14)
})


test_that("Spatial clustering rejects invalid boundary contracts", {
  x <- ch7sc_sparse_fixture()
  expect_error(k_spatial_median(x, 0), "between 1 and")
  expect_error(k_spatial_median(x, 2, init = c(1, 1)), "distinct valid")
  expect_error(k_spatial_median(x, 2, ties = "random"), "ties")
  expect_error(k_spatial_median(x, 2, empty_action = "random"),
               "empty_action")
  expect_error(k_spatial_median(x, 2, cycle_action = "ignore"),
               "cycle_action")
  expect_error(k_spatial_median(x, 2, zero_tol = -1), "non-negative")
  expect_error(sm_sscm(x, 2, lambda = 0), "strictly positive")
  expect_error(sm_sscm(x, 2, lambda = Inf), "strictly positive")
  expect_error(sparse_k_spatial_median(x, 2, tau = -1), "non-negative")
  expect_error(sparse_sm_select_tau(
    x, 2, tau_grid = numeric(), permutations = ch7sc_permutations(6, 3)
  ), "tau_grid")
  expect_error(sparse_sm_select_k(
    x, c(1, 2), tau_grid = 0, permutations = ch7sc_permutations(6, 3)
  ), "2 <= K < n")
  expect_error(sparse_sm_select_k(
    x, nrow(x), tau_grid = 0, permutations = ch7sc_permutations(6, 3)
  ), "2 <= K < n")
  expect_error(k_spatial_median(cbind(x, NA_real_), 2), "finite")

  expect_error(
    HDElliptical:::cpp_ch7sc_assign_metric(
      matrix(c(0, 0), 1, 2), matrix(c(1, 1), 1, 2),
      matrix(c(1, 0.1, 0, 1), 2, 2), TRUE
    ),
    "exactly symmetric and SPD"
  )

  finite.distance <- HDElliptical:::cpp_ch7sc_assign_metric(
    matrix(1e200, 1, 1), matrix(0, 1, 1), matrix(1e-308, 1, 1), TRUE
  )
  expect_equal(
    as.numeric(finite.distance$assigned_distance) / 1e92, 1, tolerance = 1e-14
  )
})
