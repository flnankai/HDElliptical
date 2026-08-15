# SEMC formula, certificate, fixed-fixture, RNG, and boundary tests.
#
# Provenance: Feng and Zhuang (2026), arXiv:2605.08995, and the authors'
# MIT-licensed GEMcluster implementation at commit
# 10fce04fe690fe274dd5d237cfcd3d5c6a4139f6. No simulation is run here.


.semc_test_fixture <- function() {
  x <- matrix(c(
    -3.0, -2.0,
    -2.0, -3.0,
    -2.0, -1.0,
    -1.0, -2.0,
    -2.5, -2.0,
    -1.5, -2.2,
     3.0,  2.0,
     2.0,  3.0,
     2.0,  1.0,
     1.0,  2.0,
     2.5,  2.0,
     1.5,  2.2
  ), ncol = 2L, byrow = TRUE,
  dimnames = list(NULL, c("feature1", "feature2")))
  list(x = x, labels = rep(1:2, each = 6L))
}


.semc_test_args <- function(shape = "tyler") {
  fixture <- .semc_test_fixture()
  list(
    shape = shape,
    initialization = "labels",
    initial_labels = fixture$labels,
    init_tau = 0,
    init_nstart = 1L,
    outer_nstart = 1L,
    poet_factor_selection = "supplied",
    poet_factors = 0L,
    poet_threshold = 0.1,
    poet_ridge = 0.05,
    tyler_ridge = 0.1,
    tyler_tol = 1e-5,
    tyler_max_iter = 200L,
    glasso_lambda = if (shape == "glasso") 0.2 else NULL,
    glasso_tol = 1e-5,
    glasso_max_iter = 20000L,
    max_iter = 10L,
    convergence_tol = 0.2,
    generator_grid_size = 60L,
    strict = TRUE,
    keep_path = TRUE
  )
}


test_that("SEMC quadratic radii and posterior softmax match literal formulas", {
  x <- matrix(c(1, 2, -1, 0.5, 2, -0.5), ncol = 2L, byrow = TRUE)
  centers <- matrix(c(0, 0, 1, 1), ncol = 2L, byrow = TRUE)
  precision <- matrix(c(2, 0.25, 0.25, 1.5), ncol = 2L)
  observed <- HDElliptical:::cpp_ch7_semc_delta(x, centers, precision)
  expected <- vapply(seq_len(nrow(centers)), function(k) {
    residual <- sweep(x, 2L, centers[k, ], "-")
    rowSums((residual %*% precision) * residual)
  }, numeric(nrow(x)))
  expect_equal(observed, expected, tolerance = 1e-13)

  log_scores <- matrix(c(1000, 999, -2, 3, 1, 1), ncol = 2L, byrow = TRUE)
  result <- HDElliptical:::cpp_ch7_semc_softmax(log_scores)
  maximum <- apply(log_scores, 1L, max)
  unnormalized <- exp(log_scores - maximum)
  expected_probability <- unnormalized / rowSums(unnormalized)
  expect_equal(result$probabilities, expected_probability,
               tolerance = 1e-15)
  expect_equal(as.numeric(result$log_normalizers),
               maximum + log(rowSums(unnormalized)), tolerance = 1e-15)
  expect_equal(rowSums(result$probabilities), rep(1, nrow(log_scores)),
               tolerance = 1e-15)
})


test_that("SEMC sign scatter and one Tyler map match their equations", {
  residuals <- matrix(c(
    1, 0,
    0, 2,
    0, 0
  ), ncol = 2L, byrow = TRUE)
  weights <- c(1, 2, 3)
  radial_floor <- 0.5
  radii_squared <- rowSums(residuals^2)
  coefficients <- weights / pmax(radii_squared, radial_floor)
  expected <- crossprod(residuals, residuals * coefficients) / sum(weights)
  sign_fit <- HDElliptical:::cpp_ch7_semc_weighted_sign_scatter(
    residuals, weights, radial_floor
  )
  expect_equal(sign_fit$estimate, expected, tolerance = 1e-15)
  expect_equal(sign_fit$radial_floor_uses, 1L)

  tyler_residuals <- matrix(c(
     1, 0,
     0, 2,
     1, 1,
    -1, 1
  ), ncol = 2L, byrow = TRUE)
  tyler_weights <- c(1, 2, 1, 3)
  initial_shape <- diag(c(1.5, 0.5))
  ridge <- 0.2
  inverse <- solve(initial_shape)
  quadratic <- rowSums((tyler_residuals %*% inverse) * tyler_residuals)
  coefficient <- 2 * tyler_weights / (pmax(quadratic, 0.1) *
                                        sum(tyler_weights))
  proposal <- crossprod(tyler_residuals, tyler_residuals * coefficient)
  proposal <- (1 - ridge) * proposal + ridge * diag(2L)
  proposal <- proposal * 2 / sum(diag(proposal))
  tyler <- HDElliptical:::cpp_ch7_semc_weighted_tyler(
    tyler_residuals, tyler_weights, initial_shape,
    ridge, 0.1, 1e6, 1L, 0
  )
  expect_true(tyler$converged)
  expect_equal(tyler$iterations, 1L)
  expect_equal(tyler$shape, proposal, tolerance = 1e-13)
  expect_equal(sum(diag(tyler$shape)), 2, tolerance = 1e-14)
  expect_gt(min(eigen(tyler$shape, symmetric = TRUE)$values), 0)
})


test_that("SEMC KDE and off-diagonal glasso expose exact native contracts", {
  grid <- c(-1, 0, 1)
  observations <- c(-0.5, 0.75)
  weights <- c(1, 3)
  bandwidth <- 0.4
  expected <- vapply(grid, function(point) {
    sum(weights * stats::dnorm((point - observations) / bandwidth)) /
      (bandwidth * sum(weights))
  }, numeric(1L))
  density <- HDElliptical:::cpp_ch7_semc_weighted_kde(
    grid, observations, weights, bandwidth
  )
  expect_equal(as.numeric(density), expected, tolerance = 1e-15)

  diagonal_scatter <- diag(c(2, 4))
  diagonal_fit <- HDElliptical:::cpp_ch7_semc_offdiag_glasso(
    diagonal_scatter, 0.3, 5e-6, 20000L, 1, 100L, 0, 1e-12
  )
  expect_true(diagonal_fit$converged)
  expect_false(diagonal_fit$backtracking_failed)
  expect_false(diagonal_fit$diagonal_penalty)
  expect_lte(diagonal_fit$kkt_residual, 5e-6)
  expect_equal(diagonal_fit$solution, diag(c(0.5, 0.25)),
               tolerance = 1e-6)
  expect_gt(diagonal_fit$minimum_eigenvalue, 0)

  uncertified <- HDElliptical:::cpp_ch7_semc_offdiag_glasso(
    matrix(c(2, 0.8, 0.8, 1), 2L), 0.1, 1e-14,
    1L, 1, 100L, 0, 1e-12
  )
  expect_false(uncertified$converged)
  expect_true(uncertified$backtracking_failed ||
                uncertified$kkt_residual > 1e-14)
})


test_that("all three SEMC shape pipelines pass strict fixed-fixture certificates", {
  fixture <- .semc_test_fixture()
  pipelines <- c(
    tyler = "weighted Tyler",
    poet = "weighted Tyler -> POET",
    glasso = "off-diagonal graphical lasso"
  )
  fits <- lapply(names(pipelines), function(shape) {
    do.call(semc_fit, c(
      list(x = fixture$x, K = 2L),
      .semc_test_args(shape)
    ))
  })
  names(fits) <- names(pipelines)

  for (shape in names(fits)) {
    fit <- fits[[shape]]
    expect_s3_class(fit, "semc_fit")
    expect_true(fit$valid)
    expect_true(fit$converged)
    expect_equal(rowSums(fit$posterior), rep(1, nrow(fixture$x)),
                 tolerance = 1e-12)
    expect_equal(sum(fit$mixing), 1, tolerance = 1e-14)
    expect_equal(sum(diag(fit$shape)), ncol(fixture$x),
                 tolerance = 1e-12)
    expect_gt(min(eigen(fit$shape, symmetric = TRUE)$values), 0)
    expect_gt(min(eigen(fit$precision, symmetric = TRUE)$values), 0)
    expect_equal(unname(fit$shape %*% fit$precision), diag(ncol(fixture$x)),
                 tolerance = 1e-10)
    expect_true(fit$shape.fit$certificate$valid)
    expect_match(fit$shape.fit$controls$nested.pipeline, pipelines[[shape]],
                 fixed = TRUE)
    expect_equal(predict(fit, type = "posterior"), fit$posterior)
    expect_identical(predict(fit), fit$cluster)
    expect_true(all(fit$cluster[1:6] == fit$cluster[1L]))
    expect_true(all(fit$cluster[7:12] == fit$cluster[7L]))
    expect_false(fit$cluster[1L] == fit$cluster[7L])

    generator <- fit$generator
    p <- ncol(fixture$x)
    integral <- pi^(p / 2) / gamma(p / 2) *
      HDElliptical:::.semc_trapezoid(
        generator$u.grid,
        generator$u.grid^(p / 2 - 1) * exp(generator$log.g.grid)
      )
    expect_equal(integral, 1, tolerance = 2e-12)
    expect_gte(min(generator$score.grid),
               fit$controls$generator$score.clip[1L])
    expect_lte(max(generator$score.grid),
               fit$controls$generator$score.clip[2L])
    expect_match(fit$diagnostics$no.hidden.repair, "no eigenvalue projection",
                 ignore.case = TRUE)
    expect_equal(fit$controls$shape$poet.max.factors, 0L)
    expect_equal(fit$controls$shape$poet.max.factors.requested, 8L)
  }

  glasso <- fits$glasso$shape.fit$glasso
  expect_false(glasso$diagonal.penalty)
  expect_false(glasso$fit$diagonal_penalty)
  expect_true(glasso$fit$converged)
  expect_lte(glasso$fit$kkt_residual,
             fits$glasso$controls$shape$glasso.tol)
  expect_equal(glasso$selected.lambda, 0.2)
  expect_identical(
    fits$tyler$provenance$official.software,
    paste0("flnankai/GEMcluster@",
           "10fce04fe690fe274dd5d237cfcd3d5c6a4139f6")
  )
})


test_that("paper and official-software Gap dispersions stay distinct", {
  x <- matrix(c(-2, 0, 0, 1, 3, -1), ncol = 2L, byrow = TRUE)
  fit <- list(
    centers = matrix(c(-1, 0, 2, 0), ncol = 2L, byrow = TRUE),
    precision = matrix(c(1.5, 0.2, 0.2, 1), 2L),
    cluster = c(1L, 1L, 2L),
    posterior = matrix(c(0.8, 0.2, 0.6, 0.4, 0.1, 0.9),
                       ncol = 2L, byrow = TRUE)
  )
  delta <- HDElliptical:::cpp_ch7_semc_delta(
    x, fit$centers, fit$precision
  )
  paper <- mean(log1p(delta[cbind(seq_len(nrow(x)), fit$cluster)]))
  software <- sum(fit$posterior * delta)
  expect_equal(
    HDElliptical:::.semc_gap_dispersion(x, fit, "paper_hard_log1p"),
    paper, tolerance = 1e-15
  )
  expect_equal(
    HDElliptical:::.semc_gap_dispersion(x, fit, "software_soft_delta"),
    software, tolerance = 1e-15
  )
  expect_false(isTRUE(all.equal(paper, software)))
})


test_that("Gap-LSE uses fixed permutations, paper formulas, and local RNG", {
  fixture <- .semc_test_fixture()
  n <- nrow(fixture$x)
  p <- ncol(fixture$x)
  B <- 2L
  plan <- array(NA_integer_, dim = c(n, p, B))
  for (b in seq_len(B)) {
    for (j in seq_len(p)) {
      plan[, j, b] <- seq_len(n)
    }
  }
  control <- .semc_test_args("tyler")
  set.seed(709L)
  before <- .Random.seed
  selection <- semc_select_k_gap(
    fixture$x, k_grid = 2L, B = B, control = control,
    dispersion = "paper_hard_log1p", rule = "lse",
    permutation_indices = plan, seed = 817L,
    keep_reference_fits = FALSE, keep_permutations = TRUE
  )
  expect_identical(.Random.seed, before)
  expect_s3_class(selection, "semc_gap_selection")
  expect_equal(selection$selected.k, 2L)
  expect_equal(selection$lse.k, 2L)
  expect_equal(selection$maximum.gap.k, 2L)
  expect_equal(selection$summary$mean.log.reference,
               rowMeans(selection$reference.log.dispersion),
               tolerance = 1e-15)
  expect_equal(selection$summary$gap,
               rowMeans(selection$reference.log.dispersion) -
                 log(selection$summary$observed.dispersion),
               tolerance = 1e-15)
  expect_equal(selection$summary$standard.error,
               sqrt(1 + 1 / B) *
                 apply(selection$reference.log.dispersion, 1L, stats::sd),
               tolerance = 1e-15)
  expect_equal(selection$summary$gap, 0, tolerance = 1e-11)
  expect_identical(selection$permutation.indices, plan)
  expect_true(selection$diagnostics$paper.software.conflict.exposed)
  expect_match(selection$diagnostics$formula,
               "mean_i log(1 + delta_i,C_i)", fixed = TRUE)
})


test_that("SEMC stochastic contracts isolate RNG and require explicit seeds", {
  fixture <- .semc_test_fixture()
  args <- .semc_test_args("tyler")
  args$initialization <- "software_random"
  args$initial_labels <- NULL
  args$init_nstart <- 2L
  args$init_empty_action <- "farthest"
  args$seed <- 91L

  missing_seed <- args
  missing_seed$seed <- NULL
  expect_error(
    do.call(semc_fit, c(list(x = fixture$x, K = 2L), missing_seed)),
    "explicit seed"
  )

  set.seed(1301L)
  before <- .Random.seed
  first <- do.call(semc_fit, c(list(x = fixture$x, K = 2L), args))
  expect_identical(.Random.seed, before)
  second <- do.call(semc_fit, c(list(x = fixture$x, K = 2L), args))
  expect_identical(.Random.seed, before)
  expect_identical(first$cluster, second$cluster)
  expect_equal(first$centers, second$centers, tolerance = 0)
  expect_true(first$controls$initialization$software.contract)
})


test_that("SEMC fails at explicit SPD and convergence boundaries without repair", {
  rank_one_residuals <- matrix(c(
    -2, 0,
    -1, 0,
     1, 0,
     2, 0
  ), ncol = 2L, byrow = TRUE)
  expect_error(
    HDElliptical:::cpp_ch7_semc_weighted_tyler(
      rank_one_residuals, rep(1, 4L), diag(2L),
      0, 1e-6, 1e6, 1L, 0
    ),
    "not strictly positive definite"
  )
  explicit_ridge <- HDElliptical:::cpp_ch7_semc_weighted_tyler(
    rank_one_residuals, rep(1, 4L), diag(2L),
    0.1, 1e-6, 1e6, 1L, 0
  )
  expect_true(explicit_ridge$converged)
  expect_equal(explicit_ridge$ridge, 0.1)
  expect_gt(explicit_ridge$minimum_eigenvalue, 0)

  fixture <- .semc_test_fixture()
  fit <- do.call(semc_fit, c(
    list(x = fixture$x, K = 2L), .semc_test_args("tyler")
  ))
  expect_error(predict(fit, matrix(0, 2L, 3L)), "feature dimension")
  renamed <- fixture$x[1:2, , drop = FALSE]
  colnames(renamed) <- rev(colnames(renamed))
  expect_error(predict(fit, renamed), "column names and order")
  expect_error(predict(fit, extra = TRUE), "No additional")

  expect_error(
    semc_fit(fixture$x, 2L, glasso_lambda = 0.1,
             glasso_lambda_grid = c(0.2, 0.1)),
    "at most one"
  )
  bad_plan <- array(1L, dim = c(nrow(fixture$x), 2L, 2L))
  expect_error(
    semc_select_k_gap(fixture$x, 2L, B = 2L,
                      permutation_indices = bad_plan, seed = 1L),
    "must permute"
  )
})


test_that("SEMC runs on a frozen subset of the pinned official fixture", {
  # Literal official rows 1, 2, 4, 7, 3, 9, 11, 12, 5, 6, 8, and 10
  # from the set.seed(2), n=36, p=6 Gaussian fixture in pinned test-gem.R.
  # The numbers are frozen here; no simulation is executed by this test.
  x <- matrix(c(
     2.03831275788377,  2.39462180375117,  0.832539591177201,
    -0.18016835935994, -0.856093149616075, 1.34022251299227,
     1.46270626172199, -0.942181108916423, -0.529365285184869,
     0.489209200923289, -0.477893411082912, -0.64782304106809,
     1.82760400942345,  2.20028825124934,  3.08845161457641,
     2.66443728539639, -1.19579798202605, -1.47638616020004,
    -0.178572652382754, -1.01744139504381, 1.36766109206865,
     0.774466828236058, 0.543433604952253, -0.507474615729431,
     1.4233120196421, 2.44519099795236, 1.68742393071414,
     0.0617325888625494, 0.70573830213935, 0.518328081148286,
     0.464307554999709, 0.128946054368327, 0.135438333176472,
    -1.70190147862167, -1.32809182629611, -1.80020264197458,
    -0.615684235217505, -0.664403867013466, 1.65755088068686,
    -0.252054571478647, 2.44563233073367, 0.481456507007696,
    -0.074189307135252, 0.10427660135395, 2.28107481275925,
     1.3312556119978, 0.975325236524682, 0.260476587640124,
    -2.04351445007947, 0.61796195014417, -1.05218440267279,
    -0.861864049782676, 0.434894671603151, -1.92423026347599,
    -0.27227643185354, -0.596184098950831, -0.0409355859849468,
    -1.56370787332449, 0.145434735415286, -0.155853477414398,
    -1.92897327289182, 1.00424217865491, -1.20628478707983,
    -1.89545206259417, 1.05169906257342, -1.26951259603692,
    -0.778507369496471, 0.183728033200902, 0.00425347908939661,
    -0.708431657371514, -0.0773912830230504, -1.89243157197245
  ), nrow = 12L, ncol = 6L)
  labels <- rep(1:3, each = 4L)
  expect_equal(sum(x), 8.8601750343878471, tolerance = 1e-14)
  expect_equal(sum(x^2), 117.93265446310033, tolerance = 1e-12)

  fit <- semc_fit(
    x, K = 3L, shape = "tyler",
    initialization = "labels", initial_labels = labels,
    init_tau = 0, init_empty_action = "farthest",
    poet_factor_selection = "supplied", poet_factors = 0L,
    poet_threshold = 0.1, poet_ridge = 0.1,
    tyler_ridge = 0.1, tyler_tol = 1e-5, tyler_max_iter = 500L,
    max_iter = 1L, convergence_tol = 10,
    generator_grid_size = 50L, strict = TRUE
  )
  expect_true(fit$valid)
  expect_equal(dim(fit$posterior), c(12L, 3L))
  expect_equal(sum(diag(fit$shape)), 6, tolerance = 1e-12)
  expect_gt(fit$shape.fit$certificate$minimum.shape.eigenvalue, 0)
  expect_identical(
    fit$provenance$official.software,
    paste0("flnankai/GEMcluster@",
           "10fce04fe690fe274dd5d237cfcd3d5c6a4139f6")
  )
})


test_that("strict FALSE returns a diagnostic-only nonconverged iterate", {
  fixture <- .semc_test_fixture()
  args <- .semc_test_args("tyler")
  args$max_iter <- 1L
  args$convergence_tol <- 1e-15
  args$strict <- FALSE
  args$keep_path <- TRUE

  fit <- do.call(semc_fit, c(
    list(x = fixture$x, K = 2L), args
  ))
  expect_s3_class(fit, "semc_fit")
  expect_false(fit$converged)
  expect_false(fit$valid)
  expect_equal(fit$iterations, 1L)
  expect_true(fit$shape.fit$certificate$valid)
  expect_equal(nrow(fit$path), 2L)
  final_change <- unlist(tail(
    fit$path[c("center.change", "precision.change", "mixing.change")], 1L
  ))
  expect_gt(max(final_change), args$convergence_tol)
  expect_error(predict(fit), "valid certified semc_fit")
})

