# Formula, invariance, degeneracy, RNG, and certificate tests for Chapter 7
# classical/sparse clustering.  These are small deterministic method tests,
# not reproductions of any paper simulation.

ch7cs_partition <- function(cluster) outer(cluster, cluster, `==`)


ch7cs_ordered_bcss <- function(x, cluster) {
  vapply(seq_len(ncol(x)), function(j) {
    total <- sum(outer(x[, j], x[, j], "-")^2) / nrow(x)
    within <- sum(vapply(sort(unique(cluster)), function(k) {
      values <- x[cluster == k, j]
      sum(outer(values, values, "-")^2) / length(values)
    }, numeric(1)))
    total - within
  }, numeric(1))
}


ch7cs_half_bcss <- function(x, cluster) {
  total <- colSums(sweep(x, 2L, colMeans(x), "-")^2)
  within <- numeric(ncol(x))
  for (k in sort(unique(cluster))) {
    group <- x[cluster == k, , drop = FALSE]
    within <- within + colSums(sweep(group, 2L, colMeans(group), "-")^2)
  }
  total - within
}


ch7cs_weight_reference <- function(scores, bound, tolerance = 1e-13) {
  if (bound == 1) {
    answer <- numeric(length(scores))
    answer[which(scores == max(scores))[1L]] <- 1
    return(answer)
  }
  if (sum(scores) / sqrt(sum(scores^2)) <= bound) {
    return(scores / sqrt(sum(scores^2)))
  }
  lower <- 0
  upper <- max(scores)
  for (iteration in seq_len(500L)) {
    threshold <- (lower + upper) / 2
    soft <- pmax(scores - threshold, 0)
    ratio <- if (sum(soft^2) == 0) 1 else sum(soft) / sqrt(sum(soft^2))
    if (abs(ratio - bound) <= tolerance) break
    if (ratio > bound) lower <- threshold else upper <- threshold
  }
  soft / sqrt(sum(soft^2))
}


ch7cs_kmedian_scores <- function(x, cluster) {
  global <- apply(x, 2L, stats::median)
  total <- colSums(abs(sweep(x, 2L, global, "-")))
  within <- numeric(ncol(x))
  for (k in sort(unique(cluster))) {
    group <- x[cluster == k, , drop = FALSE]
    center <- apply(group, 2L, stats::median)
    within <- within + colSums(abs(sweep(group, 2L, center, "-")))
  }
  total - within
}


test_that("Chapter 7 public formals expose initialization and repairs", {
  expect_identical(
    names(formals(lloyd_kmeans)),
    c("x", "clusters", "initial", "initialization", "first_index", "seed",
      "ties", "empty_action", "solver_tol", "solver_max_iter", "strict",
      "keep_distances")
  )
  expect_identical(
    names(formals(gaussian_mixture_em)),
    c("x", "means", "covariances", "proportions", "ridge", "solver_tol",
      "monotone_tol", "rank_tol", "symmetry_tol", "solver_max_iter",
      "strict", "keep_responsibilities")
  )
  expect_identical(
    names(formals(sparse_kmeans)),
    c("x", "clusters", "s", "initial", "initialization", "first_index",
      "seed", "ties", "empty_action", "weight_tol", "weight_max_iter",
      "solver_tol", "solver_max_iter", "strict", "keep_distances")
  )
  expect_identical(
    names(formals(sparse_kmedian)),
    c("x", "clusters", "s_w", "initial", "initialization", "first_index",
      "seed", "ties", "empty_action", "score_tol", "weight_tol",
      "weight_max_iter", "solver_tol", "solver_max_iter", "strict",
      "keep_distances")
  )
  expect_true(all(c("permutation_indices", "permutation_seed",
                    "clustering_seed", "selection") %in%
                  names(formals(sparse_kmeans_select_s))))
})


test_that("Lloyd K-means matches the explicit mean and WCSS formulas", {
  x <- rbind(c(-3, 1), c(-1, -1), c(2, 1), c(4, -1))
  colnames(x) <- c("a", "b")
  fit <- lloyd_kmeans(x, 2, initial = c(1, 3), keep_distances = TRUE)
  reference.centers <- rbind(colMeans(x[1:2, ]), colMeans(x[3:4, ]))
  reference.objective <- sum(vapply(seq_len(nrow(x)), function(i) {
    sum((x[i, ] - reference.centers[fit$cluster[i], ])^2)
  }, numeric(1)))

  expect_s3_class(fit, "lloyd_kmeans_fit")
  expect_identical(unname(fit$cluster), c(1L, 1L, 2L, 2L))
  expect_equal(unname(fit$centers), unname(reference.centers), tolerance = 0)
  expect_equal(fit$objective, reference.objective, tolerance = 0)
  expect_true(all(diff(fit$objective_trace) <= 1e-12))
  expect_true(fit$converged)
  expect_true(fit$valid)
  expect_equal(fit$diagnostics$assignment_residual, 0, tolerance = 0)
  expect_identical(dim(fit$distances), c(4L, 2L))
})


test_that("Lloyd exact assignment ties select the smallest cluster", {
  x <- matrix(c(-2, 0, 2), ncol = 1)
  fit <- lloyd_kmeans(x, 2, initial = matrix(c(-1, 1), ncol = 1))

  expect_identical(unname(fit$cluster), c(1L, 1L, 2L))
  expect_gte(fit$diagnostics$assignment_ties, 1L)
  expect_equal(unname(fit$centers), matrix(c(-1, 2), ncol = 1),
               tolerance = 0)
  expect_error(lloyd_kmeans(x, 2, initial = c(1, 3), ties = "random"),
               "smallest index")
})


test_that("max-min seeding has deterministic row-index tie rules", {
  x <- rbind(c(0, 0), c(1, 0), c(-1, 0), c(0, 0))
  fit <- lloyd_kmeans(x, 2, first_index = 1)

  expect_identical(fit$initialization$indices, c(1L, 2L))
  expect_match(fit$initialization$source, "max-min")
  expect_error(lloyd_kmeans(matrix(1, 4, 2), 2), "unique rows")
  expect_error(lloyd_kmeans(x, 5), "between 1 and 4")
})


test_that("seeded initialization is reproducible and isolates caller RNG", {
  x <- cbind(seq_len(8), c(0, 1, 0, 1, 5, 6, 5, 6))
  set.seed(807)
  old <- .Random.seed
  first <- lloyd_kmeans(
    x, 2, initialization = "random", seed = 31, empty_action = "farthest"
  )
  expect_identical(.Random.seed, old)
  second <- lloyd_kmeans(
    x, 2, initialization = "random", seed = 31, empty_action = "farthest"
  )

  expect_identical(first$initialization$indices,
                   second$initialization$indices)
  expect_equal(first$centers, second$centers, tolerance = 0)
  expect_error(lloyd_kmeans(x, 2, initialization = "random"),
               "explicit integer `seed`")
  expect_error(lloyd_kmeans(x, 2, seed = 1), "must be NULL")
})


test_that("empty clusters error by default and farthest repair is explicit", {
  x <- matrix(c(0, 1, 10), ncol = 1)
  duplicate.centers <- matrix(c(0, 0), ncol = 1)
  expect_error(lloyd_kmeans(x, 2, initial = duplicate.centers),
               "empty cluster")
  repaired <- lloyd_kmeans(
    x, 2, initial = duplicate.centers, empty_action = "farthest"
  )

  expect_identical(unname(repaired$cluster), c(1L, 1L, 2L))
  expect_identical(unname(repaired$size), c(2L, 1L))
  expect_equal(repaired$centers[2, 1], 10, tolerance = 0)
  expect_identical(repaired$diagnostics$empty_repairs, 1L)
  expect_true(repaired$valid)
  many <- lloyd_kmeans(
    matrix(c(0, 1, 10, 11), ncol = 1), 3,
    initial = matrix(0, 3, 1), empty_action = "farthest"
  )
  expect_true(all(many$size > 0L))
  expect_identical(many$diagnostics$empty_repairs, 2L)
  expect_true(many$valid)
})


test_that("Lloyd partition is rigid-motion and common-scale equivariant", {
  x <- rbind(c(-3, 0), c(-2, 1), c(2, -1), c(4, 0))
  rotation <- matrix(c(0, -1, 1, 0), 2, 2)
  base <- lloyd_kmeans(x, 2, initial = c(1, 3))
  moved <- lloyd_kmeans(
    x %*% rotation + matrix(c(7, -4), nrow(x), 2, byrow = TRUE),
    2, initial = c(1, 3)
  )
  scaled <- lloyd_kmeans(3 * x, 2, initial = c(1, 3))

  expect_identical(ch7cs_partition(moved$cluster),
                   ch7cs_partition(base$cluster))
  expect_equal(unname(moved$centers),
               unname(base$centers %*% rotation +
                        matrix(c(7, -4), 2, 2, byrow = TRUE)),
               tolerance = 1e-14)
  expect_identical(ch7cs_partition(scaled$cluster),
                   ch7cs_partition(base$cluster))
  expect_equal(scaled$objective, 9 * base$objective, tolerance = 1e-13)
})


test_that("strict Lloyd mode rejects an unfinished iterate", {
  x <- matrix(c(0, 1, 2, 10), ncol = 1)
  expect_error(lloyd_kmeans(x, 2, initial = c(1, 2), solver_max_iter = 1),
               "without its solver certificate")
  expect_warning(
    loose <- lloyd_kmeans(
      x, 2, initial = c(1, 2), solver_max_iter = 1, strict = FALSE
    ),
    "without its solver certificate"
  )
  expect_false(loose$valid)
  expect_false(loose$converged)
})


test_that("Gaussian E/M update matches a pure-R one-dimensional reference", {
  x <- matrix(c(-2, -1, 1, 3), ncol = 1)
  means <- matrix(c(-1.5, 2), ncol = 1)
  covariance <- array(c(1.2, 0.8), c(1, 1, 2))
  proportions <- c(0.4, 0.6)
  density <- cbind(
    proportions[1] * stats::dnorm(x[, 1], means[1, 1], sqrt(covariance[1, 1, 1])),
    proportions[2] * stats::dnorm(x[, 1], means[2, 1], sqrt(covariance[1, 1, 2]))
  )
  tau <- density / rowSums(density)
  mass <- colSums(tau)
  expected.means <- colSums(tau * x[, 1]) / mass
  expected.variance <- vapply(seq_len(2), function(k) {
    sum(tau[, k] * (x[, 1] - expected.means[k])^2) / mass[k]
  }, numeric(1))
  core <- cpp_ch7cs_gmm_core(
    x, means, covariance, proportions, 0, 1e-16, 1e-12, 1e-12, 1L, TRUE
  )

  expect_equal(as.numeric(core$means), expected.means, tolerance = 2e-14)
  expect_equal(as.numeric(core$covariances), expected.variance,
               tolerance = 2e-14)
  expect_equal(as.numeric(core$proportions), mass / nrow(x),
               tolerance = 2e-14)
  expect_true(all(diff(core$log_likelihood_trace) >= -1e-12))
})


test_that("Gaussian EM exposes likelihood and fixed-point certificates", {
  x <- matrix(c(-2.2, -2, -1.8, 1.8, 2, 2.2), ncol = 1)
  fit <- gaussian_mixture_em(
    x, matrix(c(-2, 2), ncol = 1), array(c(.2, .2), c(1, 1, 2)),
    c(.5, .5), solver_tol = 1e-7
  )

  expect_s3_class(fit, "gaussian_mixture_em_fit")
  expect_true(fit$valid)
  expect_true(fit$converged)
  expect_equal(unname(rowSums(fit$responsibilities)), rep(1, nrow(x)),
               tolerance = 2e-15)
  expect_equal(sum(fit$proportions), 1, tolerance = 2e-15)
  expect_equal(sum(fit$diagnostics$effective_masses), nrow(x),
               tolerance = 2e-15)
  expect_true(all(diff(fit$log_likelihood_trace) >= -1e-10))
  expect_lte(fit$diagnostics$fixed_point_residual, 1e-7)
  expect_match(fit$diagnostics$inversion, "no inverse or pseudo-inverse")
})


test_that("Gaussian responsibilities respect affine orthogonal re-expression", {
  x <- rbind(c(-3, 0), c(-2, 1), c(-2, -1),
             c(2, 0), c(3, 1), c(3, -1))
  means <- rbind(c(-2.5, 0), c(2.5, 0))
  covariances <- array(0, c(2, 2, 2))
  covariances[, , 1] <- diag(c(1, .8))
  covariances[, , 2] <- diag(c(.9, 1.1))
  base <- gaussian_mixture_em(
    x, means, covariances, c(.5, .5), solver_tol = 1e-7
  )
  q <- matrix(c(0, -1, 1, 0), 2, 2)
  shift <- c(9, -2)
  transformed.cov <- array(0, c(2, 2, 2))
  for (k in 1:2) transformed.cov[, , k] <-
    t(q) %*% covariances[, , k] %*% q
  moved <- gaussian_mixture_em(
    x %*% q + matrix(shift, nrow(x), 2, byrow = TRUE),
    means %*% q + matrix(shift, 2, 2, byrow = TRUE),
    transformed.cov, c(.5, .5), solver_tol = 1e-7
  )

  expect_equal(moved$responsibilities, base$responsibilities,
               tolerance = 2e-9)
  expect_identical(ch7cs_partition(moved$cluster),
                   ch7cs_partition(base$cluster))
  expect_equal(moved$log_likelihood, base$log_likelihood, tolerance = 2e-9)
})


test_that("Gaussian EM never hides singular covariance repair", {
  x <- rbind(c(0, 0, 0), c(1, 0, 0), c(0, 1, 0))
  expect_error(gaussian_mixture_em(
    x, means = matrix(colMeans(x), nrow = 1), covariances = diag(3),
    proportions = 1
  ), "singular")
  regularized <- gaussian_mixture_em(
    x, means = matrix(colMeans(x), nrow = 1), covariances = diag(3),
    proportions = 1, ridge = .05, solver_tol = 1e-9
  )

  expect_true(regularized$valid)
  expect_true(regularized$diagnostics$regularized)
  expect_false(regularized$diagnostics$monotone)
  expect_equal(regularized$diagnostics$ridge, .05, tolerance = 0)
  expect_identical(unname(regularized$diagnostics$raw_covariance_ranks), 2L)
  expect_gte(min(regularized$diagnostics$minimum_eigenvalues), .05 - 1e-12)
})


test_that("Gaussian EM validates the full-covariance input contract", {
  x <- matrix(c(-1, 0, 1, 2), ncol = 1)
  expect_error(gaussian_mixture_em(x, matrix(0, 1, 1), 1, .9),
               "sum exactly")
  expect_error(gaussian_mixture_em(x, matrix(0, 1, 1), matrix(0, 1, 1), 1),
               "positive definite")
  expect_error(gaussian_mixture_em(
    cbind(x, x^2), matrix(c(0, 0), 1, 2),
    matrix(c(1, .1, 0, 1), 2), 1, symmetry_tol = 0
  ), "not symmetric")
  expect_error(gaussian_mixture_em(x, matrix(0, 1, 1), 1, 1, ridge = -1),
               "non-negative")
})


test_that("sparse K-means returns the displayed ordered-pair BCSS", {
  x <- rbind(c(-3, -1), c(-2, 0), c(2, 1), c(4, 2))
  colnames(x) <- c("signal", "secondary")
  fit <- sparse_kmeans(x, 2, sqrt(2), initial = c(1, 3))
  ordered <- ch7cs_ordered_bcss(x, fit$cluster)
  half <- ch7cs_half_bcss(x, fit$cluster)

  expect_s3_class(fit, "sparse_kmeans_fit")
  expect_equal(unname(fit$feature_scores), ordered, tolerance = 2e-14)
  expect_equal(unname(fit$feature_scores), 2 * unname(half), tolerance = 2e-14)
  expect_equal(unname(fit$diagnostics$half_scaled_feature_scores), unname(half),
               tolerance = 2e-14)
  expect_equal(fit$objective, sum(fit$weights * ordered), tolerance = 2e-14)
  expect_match(fit$diagnostics$feature_score_identity, "2 \\* \\(TSS - WCSS\\)")
  expect_true(fit$valid)
})


test_that("sparse K-means active weight block matches pure-R soft thresholding", {
  x <- rbind(c(-3, -1), c(-2, -1), c(2, 1), c(3, 1))
  bound <- 1.1
  fit <- sparse_kmeans(x, 2, bound, initial = c(1, 3),
                       weight_tol = 1e-12)
  reference <- ch7cs_weight_reference(fit$feature_scores, bound)

  expect_equal(unname(fit$weights), unname(reference), tolerance = 2e-11)
  expect_equal(sum(fit$weights), bound, tolerance = 2e-10)
  expect_equal(sqrt(sum(fit$weights^2)), 1, tolerance = 2e-14)
  expect_true(fit$diagnostics$weight_certificate$constraint_active)
  expect_true(fit$diagnostics$weight_certificate$certified)
  expect_lte(fit$diagnostics$weight_certificate$stationarity_residual, 1e-9)
  expect_true(all(fit$weights >= 0))
})


test_that("sparse K-means s=1 ties use the first feature", {
  x <- rbind(c(-2, -2), c(-1, -1), c(1, 1), c(2, 2))
  fit <- sparse_kmeans(x, 2, 1, initial = c(1, 3))

  expect_equal(unname(fit$feature_scores),
               rep(unname(fit$feature_scores[1]), 2),
               tolerance = 0)
  expect_identical(unname(fit$weights), c(1, 0))
  expect_true(fit$diagnostics$weight_certificate$s_one_tie_branch)
  expect_equal(fit$diagnostics$weight_certificate$threshold,
               unname(fit$feature_scores[1]), tolerance = 0)
  expect_equal(fit$diagnostics$weight_certificate$gamma, 0, tolerance = 0)
  expect_equal(fit$diagnostics$weight_certificate$stationarity_residual, 0)
  expect_error(sparse_kmeans(x, 1, 1), "feature scores are zero")
  expect_error(sparse_kmeans(x, 2, .9), "must lie")
  expect_error(sparse_kmeans(x, 2, 2), "must lie")
})


test_that("sparse K-means is translation and feature-permutation equivariant", {
  x <- rbind(c(-4, -1), c(-2, 0), c(2, 1), c(5, 2))
  base <- sparse_kmeans(x, 2, 1.25, initial = c(1, 3))
  moved <- sparse_kmeans(
    sweep(x, 2L, c(10, -7), "+"), 2, 1.25, initial = c(1, 3)
  )
  permutation <- c(2, 1)
  permuted <- sparse_kmeans(
    x[, permutation], 2, 1.25, initial = c(1, 3)
  )

  expect_identical(ch7cs_partition(moved$cluster),
                   ch7cs_partition(base$cluster))
  expect_equal(moved$weights, base$weights, tolerance = 2e-14)
  expect_equal(moved$feature_scores, base$feature_scores, tolerance = 2e-13)
  expect_identical(ch7cs_partition(permuted$cluster),
                   ch7cs_partition(base$cluster))
  expect_equal(unname(permuted$weights), unname(base$weights[permutation]),
               tolerance = 2e-12)
  expect_equal(unname(permuted$feature_scores),
               unname(base$feature_scores[permutation]), tolerance = 2e-12)
})


test_that("sparse K-means exposes empty-cluster and unfinished-solver states", {
  x <- rbind(c(0, 0), c(1, 0), c(10, 1))
  centers <- rbind(c(0, 0), c(0, 0))
  expect_error(sparse_kmeans(x, 2, 1, initial = centers), "empty cluster")
  repaired <- sparse_kmeans(
    x, 2, 1, initial = centers, empty_action = "farthest"
  )
  expect_gte(repaired$diagnostics$empty_repairs, 1L)
  expect_true(repaired$valid)
  expect_error(sparse_kmeans(
    matrix(c(0, 1, 2, 10), ncol = 1), 2, 1,
    initial = c(1, 2), solver_max_iter = 1
  ), "without its solver certificate")
})


test_that("permutation selector matches explicit objective and Gap formulas", {
  x <- rbind(c(-3, 0), c(-2, 1), c(2, 0), c(3, 1))
  indices <- array(c(
    4, 3, 2, 1, 2, 1, 4, 3,
    2, 4, 1, 3, 3, 1, 4, 2
  ), dim = c(4, 2, 2))
  fit <- sparse_kmeans_select_s(
    x, 2, c(1, sqrt(2)), permutation_indices = indices,
    selection = "max", initial = c(1, 3), empty_action = "farthest",
    keep_permutations = TRUE
  )
  manual.gap <- log(fit$calibration$observed_objective) -
    colMeans(log(fit$permutation_objectives))

  expect_s3_class(fit, "sparse_kmeans_selection")
  expect_equal(unname(fit$calibration$gap), unname(manual.gap),
               tolerance = 2e-15)
  expect_equal(unname(fit$calibration$standard_error),
               unname(sqrt(1 + 1 / 2) *
                        apply(log(fit$permutation_objectives), 2, stats::sd)),
               tolerance = 2e-15)
  expect_equal(fit$permutation_indices, indices, tolerance = 0)
  expect_identical(fit$selected_index,
                   which(fit$calibration$gap == max(fit$calibration$gap))[1L])
  expect_equal(fit$selected_s, fit$calibration$s[fit$selected_index],
               tolerance = 0)
})


test_that("generated permutation calibration is reproducible and RNG-isolated", {
  x <- rbind(c(-3, 0), c(-2, 1), c(2, 0), c(3, 1))
  set.seed(1907)
  old <- .Random.seed
  generated <- sparse_kmeans_select_s(
    x, 2, c(1, sqrt(2)), n_permutations = 3, permutation_seed = 21,
    selection = "one_se", initial = c(1, 3), empty_action = "farthest",
    keep_permutations = TRUE
  )
  expect_identical(.Random.seed, old)
  supplied <- sparse_kmeans_select_s(
    x, 2, c(1, sqrt(2)),
    permutation_indices = generated$permutation_indices,
    selection = "one_se", initial = c(1, 3), empty_action = "farthest",
    keep_permutations = TRUE
  )

  expect_equal(generated$calibration, supplied$calibration, tolerance = 0)
  expect_identical(generated$selected_s, supplied$selected_s)
  expect_equal(generated$permutation_objectives,
               supplied$permutation_objectives, tolerance = 0)
  expect_true(generated$selected_s %in% generated$calibration$s)
})


test_that("selector tie and validation rules are explicit", {
  x <- cbind(c(-3, -2, 2, 3), 0)
  indices <- array(rep(c(4, 3, 2, 1, 2, 1, 4, 3), 2), c(4, 2, 2))
  fit <- sparse_kmeans_select_s(
    x, 2, c(1, sqrt(2)), permutation_indices = indices,
    selection = "one_se", initial = c(1, 3), empty_action = "farthest"
  )

  expect_identical(fit$selected_s, 1)
  expect_error(sparse_kmeans_select_s(
    x, 2, 1, n_permutations = 1, permutation_seed = 1,
    selection = "one_se", initial = c(1, 3)
  ), "at least two")
  bad <- indices
  bad[1, 1, 1] <- 1.5
  expect_error(sparse_kmeans_select_s(
    x, 2, 1, permutation_indices = bad, initial = c(1, 3)
  ), "finite integer")
  bad <- indices
  bad[, 1, 1] <- 1:4
  bad[4, 1, 1] <- 3
  expect_error(sparse_kmeans_select_s(
    x, 2, 1, permutation_indices = bad, initial = c(1, 3)
  ), "not a permutation")
})


test_that("sparse K-median matches midpoint-median improvement formulas", {
  x <- rbind(c(-4, 0), c(-2, 2), c(2, 0), c(8, 2))
  colnames(x) <- c("signal", "noise")
  fit <- sparse_kmedian(x, 2, 1, initial = c(1, 3))
  expected <- ch7cs_kmedian_scores(x, fit$cluster)

  expect_s3_class(fit, "sparse_kmedian_fit")
  expect_identical(unname(fit$cluster), c(1L, 1L, 2L, 2L))
  expect_equal(unname(fit$centers), rbind(c(-3, 1), c(5, 1)), tolerance = 0)
  expect_equal(unname(fit$global_median), c(0, 1), tolerance = 0)
  expect_equal(unname(fit$feature_scores), unname(expected), tolerance = 0)
  expect_identical(unname(fit$weights), c(1, 0))
  expect_equal(fit$objective, sum(fit$weights * expected), tolerance = 0)
  expect_match(fit$diagnostics$median_convention, "midpoint")
})


test_that("sparse K-median active weights match pure-R soft thresholding", {
  x <- rbind(c(-4, -2), c(-2, -1), c(2, 1), c(6, 3))
  bound <- 1.2
  fit <- sparse_kmedian(x, 2, bound, initial = c(1, 3),
                        weight_tol = 1e-12)
  reference <- ch7cs_weight_reference(fit$feature_scores, bound)

  expect_equal(unname(fit$weights), unname(reference), tolerance = 2e-11)
  expect_equal(sum(fit$weights), bound, tolerance = 2e-10)
  expect_equal(sqrt(sum(fit$weights^2)), 1, tolerance = 2e-14)
  expect_true(fit$diagnostics$weight_certificate$constraint_active)
  expect_true(fit$diagnostics$weight_certificate$certified)
  expect_equal(fit$diagnostics$assignment_residual, 0, tolerance = 0)
})


test_that("sparse K-median has its stated scale and coordinate geometry", {
  x <- rbind(c(-5, -1), c(-2, 0), c(2, 1), c(9, 2))
  base <- sparse_kmedian(x, 2, 1.25, initial = c(1, 3))
  scaled <- sparse_kmedian(3 * x, 2, 1.25, initial = c(1, 3))
  permutation <- c(2, 1)
  permuted <- sparse_kmedian(
    x[, permutation], 2, 1.25, initial = c(1, 3)
  )

  expect_identical(ch7cs_partition(scaled$cluster),
                   ch7cs_partition(base$cluster))
  expect_equal(scaled$weights, base$weights, tolerance = 2e-12)
  expect_equal(scaled$feature_scores, 3 * base$feature_scores,
               tolerance = 2e-12)
  expect_equal(scaled$objective, 3 * base$objective, tolerance = 2e-12)
  expect_identical(ch7cs_partition(permuted$cluster),
                   ch7cs_partition(base$cluster))
  expect_equal(unname(permuted$weights), unname(base$weights[permutation]),
               tolerance = 2e-12)
  expect_false(base$diagnostics$rotation_equivariant)
})


test_that("sparse K-median fails explicitly on degenerate objectives and empties", {
  x <- rbind(c(-2, 0), c(-1, 0), c(1, 0), c(2, 0))
  expect_error(sparse_kmedian(x, 1, 1), "feature scores are zero")
  duplicate.centers <- rbind(c(-2, 0), c(-2, 0))
  expect_error(sparse_kmedian(x, 2, 1, initial = duplicate.centers),
               "empty cluster")
  repaired <- sparse_kmedian(
    x, 2, 1, initial = duplicate.centers, empty_action = "farthest"
  )

  expect_gte(repaired$diagnostics$empty_repairs, 1L)
  expect_true(repaired$valid)
  expect_identical(repaired$diagnostics$method_scope,
                   "book-defined sparse K-median baseline")
  expect_error(sparse_kmedian(x, 2, 1, initial = c(1, 3), score_tol = 0),
               "strictly positive")
})


test_that("all Chapter 7 classical/sparse examples are deterministic calls", {
  x <- rbind(c(-2, 0), c(-1, 0), c(2, 0), c(3, 0))
  expect_no_error(lloyd_kmeans(x, 2, initial = c(1, 3)))
  expect_no_error(sparse_kmeans(x, 2, 1, initial = c(1, 3)))
  expect_no_error(sparse_kmedian(x, 2, 1, initial = c(1, 3)))
  expect_no_error(gaussian_mixture_em(
    matrix(c(-2.2, -2, -1.8, 1.8, 2, 2.2), ncol = 1),
    matrix(c(-2, 2), ncol = 1), array(c(.2, .2), c(1, 1, 2)),
    c(.5, .5), solver_tol = 1e-6
  ))
})
