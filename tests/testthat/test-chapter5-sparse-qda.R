c5qda_fixture <- function() {
  x1 <- rbind(
    c(-3.0, -0.5), c(-2.2, 1.1), c(-1.4, -1.3),
    c(-2.7, 1.8), c(-0.9, 0.4), c(-1.8, -2.0)
  )
  x2 <- rbind(
    c(2.8, 0.2), c(1.9, 1.7), c(1.2, -1.8),
    c(2.5, 2.2), c(0.7, -0.2), c(1.6, -2.4)
  )
  colnames(x1) <- colnames(x2) <- c("first", "second")
  list(
    x = rbind(x1, x2),
    y = factor(rep(c("class1", "class2"), each = 6L)),
    x1 = x1, x2 = x2
  )
}


c5qda_mle_cov <- function(x) {
  centered <- sweep(x, 2L, colMeans(x), "-")
  crossprod(centered) / nrow(x)
}


c5qda_quadratic_score <- function(x, model) {
  rowSums((x %*% model$quadratic) * x) +
    drop(x %*% model$linear) + model$intercept
}


test_that("Li-Shao uses MLE covariances and all three exact inequalities", {
  fixture <- c5qda_fixture()
  sample1 <- c5qda_mle_cov(fixture$x1)
  sample2 <- c5qda_mle_cov(fixture$x2)
  mean_difference <- colMeans(fixture$x2) - colMeans(fixture$x1)
  pooled <- (nrow(fixture$x1) * sample1 + nrow(fixture$x2) * sample2) /
    (nrow(fixture$x1) + nrow(fixture$x2))

  mean_boundary <- HDElliptical:::cpp_c5qda_li_shao(
    fixture$x1, fixture$x2, abs(mean_difference[1L]), 0, 0, 0.1
  )
  expect_equal(unname(mean_boundary$sample_covariance1), unname(sample1),
               tolerance = 2e-14)
  expect_equal(unname(mean_boundary$sample_covariance2), unname(sample2),
               tolerance = 2e-14)
  expect_equal(mean_boundary$thresholded_difference[1L], 0, tolerance = 0)

  difference_boundary <- abs(sample1[1L, 2L] - sample2[1L, 2L])
  pooled_fit <- HDElliptical:::cpp_c5qda_li_shao(
    fixture$x1, fixture$x2, 0, difference_boundary,
    abs(pooled[1L, 2L]), 0.1
  )
  expect_equal(unname(pooled_fit$pooled_covariance), unname(pooled),
               tolerance = 2e-14)
  expect_equal(pooled_fit$covariance1[1L, 2L], 0, tolerance = 0)
  expect_equal(pooled_fit$covariance2[1L, 2L], 0, tolerance = 0)
  expect_gt(pooled_fit$covariance1[1L, 1L], 0)
  expect_gt(pooled_fit$covariance2[2L, 2L], 0)
})


test_that("Li-Shao public score is the primary centered formula", {
  fixture <- c5qda_fixture()
  fit <- li_shao_sparse_qda(
    fixture$x, fixture$y, threshold_mean = 0.25,
    threshold_difference = 0.2, threshold_covariance = 0.1,
    ridge = 0.05
  )
  expect_s3_class(fit, "hd_classifier_fit")
  expect_true(fit$valid)
  estimate <- fit$estimate
  z <- fixture$x[c(1L, 7L, 10L), , drop = FALSE]
  centered <- sweep(z, 2L, estimate$mean1, "-")
  dhat <- estimate$thresholded.mean.difference
  manual <- rowSums(
    (centered %*% (estimate$precision2 - estimate$precision1)) * centered
  ) - 2 * drop(centered %*% estimate$precision2 %*% dhat) +
    as.numeric(crossprod(dhat, estimate$precision2 %*% dhat)) -
    fit$diagnostics$log.determinant[1L] +
    fit$diagnostics$log.determinant[2L]
  expect_equal(predict(fit, z, type = "score"), manual, tolerance = 2e-11)
  expect_equal(c5qda_quadratic_score(z, fit$score.model), manual,
               tolerance = 2e-11)
  expect_identical(fit$score.model$tie, "class1")
  expect_match(fit$primary.orientation, "Q < 0", fixed = TRUE)
})


test_that("Li-Shao paper bisection is deterministic and preserves RNG", {
  fixture <- c5qda_fixture()
  set.seed(913)
  state <- .Random.seed
  first <- li_shao_sparse_qda(
    fixture$x, fixture$y, selection = "paper_bisection",
    bisection_tol = 0.5, ridge = 0.05
  )
  expect_identical(.Random.seed, state)
  second <- li_shao_sparse_qda(
    fixture$x, fixture$y, selection = "paper_bisection",
    bisection_tol = 0.5, ridge = 0.05
  )
  expect_equal(first$tuning$thresholds, second$tuning$thresholds, tolerance = 0)
  bounds <- HDElliptical:::.c5qda_li_bounds(
    fixture$x, as.integer(fixture$y)
  )
  expect_equal(first$tuning$bisection$initial.bounds, bounds, tolerance = 0)
  expect_true(length(first$tuning$bisection$history) >= 1L)
  expect_true(all(first$tuning$thresholds >= 0))
})


test_that("Li-Shao never invents a ridge or inverse", {
  x <- cbind(seq_len(8), 2 * seq_len(8), 3 * seq_len(8))
  y <- factor(rep(c("one", "two"), each = 4L))
  expect_error(
    li_shao_sparse_qda(
      x, y, threshold_mean = 0, threshold_difference = 0,
      threshold_covariance = 0
    ),
    "not invertible"
  )
  expect_warning(
    fit <- li_shao_sparse_qda(
      x, y, threshold_mean = 0, threshold_difference = 0,
      threshold_covariance = 0, strict = FALSE
    ),
    "not invertible"
  )
  expect_false(fit$valid)
  expect_identical(fit$diagnostics$failure.stage,
                   "thresholded_covariance_inversion")
})


test_that("Jiang matrix and vector losses match separable diagonal solutions", {
  s1 <- diag(c(1, 2))
  s2 <- diag(c(2, 1))
  lambda <- 0.2
  matrix_fit <- HDElliptical:::cpp_c5qda_jiang_matrix(
    s1, s2, lambda, 1, 1e-9, 10000L
  )
  target <- diag(s1 - s2)
  curvature <- diag(s1) * diag(s2)
  expected_d <- sign(target) * pmax(abs(target) - lambda, 0) / curvature
  expect_true(matrix_fit$converged)
  expect_equal(matrix_fit$solution, diag(expected_d), tolerance = 2e-8)
  expect_lte(matrix_fit$kkt_residual, 1e-9)

  h <- diag(c(3, 4))
  gamma <- c(1, -2)
  vector_fit <- HDElliptical:::cpp_c5qda_jiang_vector(
    h, gamma, lambda, 1e-10, 10000L
  )
  expected_beta <- sign(gamma) * pmax(abs(gamma) - lambda, 0) / diag(h)
  expect_true(vector_fit$converged)
  expect_equal(as.numeric(vector_fit$solution), expected_beta,
               tolerance = 2e-9)
  expect_lte(vector_fit$kkt_residual, 1e-10)
})


test_that("Jiang eta checks all score breakpoints with deterministic ties", {
  raw <- c(-2, -0.5, 0.25, 1.5, 3)
  group <- c(2L, 1L, 2L, 1L, 1L)
  fit <- HDElliptical:::.c5qda_eta(raw, group)
  brute_loss <- vapply(fit$candidates, function(eta) {
    mean(ifelse(raw + eta > 0, 1L, 2L) != group)
  }, numeric(1))
  expect_equal(fit$losses, brute_loss, tolerance = 0)
  expect_equal(fit$error, min(brute_loss), tolerance = 0)
  winners <- which(brute_loss == min(brute_loss))
  expected <- fit$candidates[winners][
    order(abs(fit$candidates[winners]), fit$candidates[winners])[1L]
  ]
  expect_equal(fit$eta, expected, tolerance = 0)
  expect_true(all(-raw %in% fit$candidates))
})


test_that("Jiang public estimator uses losses, symmetrization, and strict tie", {
  fixture <- c5qda_fixture()
  fit <- jiang_da_qda(
    fixture$x, fixture$y, lambda_interaction = 0.2,
    lambda_linear = 0.2, solver_tol = 1e-6
  )
  expect_true(fit$valid)
  expect_equal(fit$estimate$interaction, t(fit$estimate$interaction),
               tolerance = 0)
  expect_identical(fit$score.model$tie, "class2")
  expect_match(fit$diagnostics$estimator, "not Dantzig", fixed = TRUE)
  midpoint <- fit$estimate$midpoint
  centered <- sweep(fixture$x, 2L, midpoint, "-")
  manual <- rowSums(
    (centered %*% fit$estimate$interaction) * centered
  ) + drop(centered %*% fit$estimate$linear) + fit$estimate$eta
  expect_equal(predict(fit, fixture$x, type = "score"), manual,
               tolerance = 2e-10)
  gamma <- 4 * (fit$estimate$mean1 - fit$estimate$mean2) +
    drop((fit$estimate$covariance1 - fit$estimate$covariance2) %*%
           fit$estimate$interaction %*%
           (fit$estimate$mean1 - fit$estimate$mean2))
  expect_equal(fit$estimate$gamma, gamma, tolerance = 2e-13)
})


test_that("joint tuning is deterministic and no default penalty is invented", {
  fixture <- c5qda_fixture()
  grid <- data.frame(
    lambda_interaction = c(0.15, 0.4),
    lambda_linear = c(0.2, 0.5)
  )
  set.seed(101)
  state <- .Random.seed
  fit <- jiang_da_qda(
    fixture$x, fixture$y, parameter_grid = grid, folds = 3L,
    solver_tol = 2e-6
  )
  expect_identical(.Random.seed, state)
  expect_true(fit$valid)
  expect_identical(fit$tuning$selection, "joint_cv")
  expect_equal(nrow(fit$tuning$cv), 2L)
  expect_true(fit$tuning$winner %in% 1:2)
  expect_error(jiang_da_qda(fixture$x, fixture$y), "No default")
  expect_error(
    jiang_da_qda(fixture$x, fixture$y, lambda_interaction = 0.2),
    "supplied together"
  )
})


test_that("matrix Dantzig operator and adjoint obey the inner-product identity", {
  s1 <- matrix(c(2, 0.4, 0.4, 1), 2)
  s2 <- matrix(c(1.5, -0.2, -0.2, 3), 2)
  d <- matrix(c(1, 2, -3, 0.5), 2)
  dual <- matrix(c(-1, 0.2, 1.4, 2), 2)
  applied <- HDElliptical:::cpp_c5qda_operator(s1, s2, d, FALSE)
  adjoint <- HDElliptical:::cpp_c5qda_operator(s1, s2, dual, TRUE)
  expect_equal(sum(applied * dual), sum(d * adjoint), tolerance = 2e-13)
  reference <- 0.5 * (s1 %*% d %*% s2 + s2 %*% d %*% s1)
  expect_equal(applied, reference, tolerance = 2e-14)
})


test_that("matrix-free Dantzig solutions match diagonal linear programs", {
  s1 <- diag(c(1, 2))
  s2 <- diag(c(2, 1))
  lambda <- 0.2
  matrix_fit <- HDElliptical:::cpp_c5qda_dantzig_matrix(
    s1, s2, lambda, 1e-7, 20000L
  )
  target <- diag(s1 - s2)
  expected <- sign(target) * pmax(abs(target) - lambda, 0) /
    (diag(s1) * diag(s2))
  expect_true(matrix_fit$converged)
  expect_true(matrix_fit$matrix_free)
  expect_equal(matrix_fit$solution, diag(expected), tolerance = 2e-6)
  expect_equal(matrix_fit$solution, t(matrix_fit$solution), tolerance = 0)
  residual <- 0.5 * (
    s1 %*% matrix_fit$solution %*% s2 +
      s2 %*% matrix_fit$solution %*% s1
  ) - s1 + s2
  expect_lte(max(abs(residual)), lambda + 1e-7)
  expect_equal(matrix_fit$feasible_maxnorm, max(abs(residual)),
               tolerance = 2e-14)

  target_beta <- c(1, -1)
  vector_fit <- HDElliptical:::cpp_c5qda_dantzig_vector(
    s2, target_beta, lambda, 1e-7, 20000L
  )
  expected_beta <- sign(target_beta) *
    pmax(abs(target_beta) - lambda, 0) / diag(s2)
  expect_true(vector_fit$converged)
  expect_equal(as.numeric(vector_fit$solution), expected_beta,
               tolerance = 2e-6)
  expect_lte(max(abs(s2 %*% vector_fit$solution - target_beta)),
             lambda + 1e-7)
})


test_that("SDAR stores certified primary score and is translation invariant", {
  fixture <- c5qda_fixture()
  fit <- sdar_qda(
    fixture$x, fixture$y, lambda_D = 0.3, lambda_beta = 0.3,
    solver_tol = 1e-6, feasibility_tol = 2e-6
  )
  expect_true(fit$valid)
  expect_true(fit$diagnostics$matrix.free)
  expect_true(fit$diagnostics$feasibility.rechecked.after.symmetrization)
  expect_gt(fit$diagnostics$signed.determinant$sign, 0)
  z <- fixture$x[c(2L, 8L, 11L), , drop = FALSE]
  centered1 <- sweep(z, 2L, fit$estimate$mean1, "-")
  midpoint <- fit$estimate$midpoint
  manual <- rowSums(
    (centered1 %*% fit$estimate$interaction) * centered1
  ) - 2 * drop(sweep(z, 2L, midpoint, "-") %*% fit$estimate$beta) -
    fit$diagnostics$signed.determinant$log_determinant
  expect_equal(predict(fit, z, type = "score"), manual,
               tolerance = 3e-9)

  shift <- c(17, -9)
  translated <- sdar_qda(
    sweep(fixture$x, 2L, shift, "+"), fixture$y,
    lambda_D = 0.3, lambda_beta = 0.3,
    solver_tol = 1e-6, feasibility_tol = 2e-6
  )
  expect_true(translated$valid)
  expect_equal(
    predict(translated, sweep(z, 2L, shift, "+"), type = "score"),
    predict(fit, z, type = "score"), tolerance = 2e-7
  )
  expect_equal(translated$estimate$interaction, fit$estimate$interaction,
               tolerance = 3e-12)
  expect_equal(translated$estimate$beta, fit$estimate$beta,
               tolerance = 3e-10)
})


test_that("signed log determinant rejects a negative determinant without abs", {
  negative <- HDElliptical:::cpp_c5qda_signed_logdet(diag(c(-2, 1)))
  positive <- HDElliptical:::cpp_c5qda_signed_logdet(diag(c(-2, -3)))
  expect_false(negative$valid)
  expect_lt(negative$sign, 0)
  expect_identical(negative$failure_stage, "nonpositive_determinant")
  expect_true(positive$valid)
  expect_equal(positive$log_determinant, log(6), tolerance = 2e-14)
})


test_that("SSQDA O(np) trace equals the ordered triple-U sum term by term", {
  x <- rbind(c(0, 1), c(2, -1), c(3, 4), c(-2, 2), c(1, -3))
  center <- c(0.25, 0.4)
  core <- HDElliptical:::cpp_c5qda_ssqda_moments(x, center, 0)
  ordered <- 0
  n <- nrow(x)
  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      for (k in seq_len(n)) {
        if (length(unique(c(i, j, k))) == 3L) {
          ordered <- ordered + sum((x[i, ] - x[j, ]) *
                                     (x[k, ] - x[j, ]))
        }
      }
    }
  }
  ordered <- ordered / (n * (n - 1) * (n - 2))
  identity <- sum(sweep(x, 2L, colMeans(x), "-")^2) / (n - 1)
  expect_true(core$valid)
  expect_equal(core$trace_estimate, ordered, tolerance = 2e-14)
  expect_equal(core$trace_estimate, identity, tolerance = 2e-14)
  residual <- sweep(x, 2L, center, "-")
  signs <- residual / sqrt(rowSums(residual^2))
  expect_equal(core$sscm, crossprod(signs) / n, tolerance = 2e-14)
  expect_equal(core$covariance, core$trace_estimate * core$sscm,
               tolerance = 2e-14)
  expect_identical(core$trace_algorithm, "O(np) triple-U identity")
})


test_that("SSQDA moments transform correctly and reject zero residuals", {
  x <- rbind(c(0, 1), c(2, -1), c(3, 4), c(-2, 2), c(1, -3))
  center <- c(0.25, 0.4)
  transform <- matrix(c(0, -1, 1, 0), 2)
  base <- HDElliptical:::cpp_c5qda_ssqda_moments(x, center, 0)
  rotated <- HDElliptical:::cpp_c5qda_ssqda_moments(
    x %*% transform, drop(center %*% transform), 0
  )
  expect_equal(rotated$trace_estimate, base$trace_estimate, tolerance = 2e-14)
  expect_equal(rotated$sscm, t(transform) %*% base$sscm %*% transform,
               tolerance = 2e-14)
  expect_equal(rotated$covariance,
               t(transform) %*% base$covariance %*% transform,
               tolerance = 3e-14)

  zero <- HDElliptical:::cpp_c5qda_ssqda_moments(x, x[1L, ], 0)
  expect_false(zero$valid)
  expect_equal(zero$zero_residuals, 1)
  expect_identical(zero$failure_stage, "ssqda_zero_spatial_residual")
})


test_that("SSQDA public fit stores robust moments and primary score", {
  fixture <- c5qda_fixture()
  fit <- ssqda(
    fixture$x, fixture$y, lambda_D = 1, lambda_beta = 1,
    solver_tol = 2e-6, feasibility_tol = 3e-6
  )
  expect_true(fit$valid)
  expect_equal(fit$prior$probabilities, c(class1 = 0.5, class2 = 0.5),
               tolerance = 0)
  expect_identical(
    fit$diagnostics$robust.class1$trace.algorithm,
    "O(np) triple-U identity"
  )
  expect_true(fit$diagnostics$robust.class1$median$converged)
  expect_true(fit$diagnostics$robust.class2$median$converged)
  expect_equal(fit$estimate$covariance1,
               fit$diagnostics$robust.class1$trace * fit$estimate$sscm1,
               tolerance = 2e-13)
  expect_equal(fit$estimate$covariance2,
               fit$diagnostics$robust.class2$trace * fit$estimate$sscm2,
               tolerance = 2e-13)
  expect_gt(fit$diagnostics$signed.determinant$sign, 0)
  expect_true(fit$diagnostics$matrix.free)
})


test_that("SSQDA enforces class size and median contracts", {
  x <- rbind(c(0, 0), c(1, 0), c(2, 0), c(3, 0), c(4, 0))
  y <- factor(c("a", "a", "b", "b", "b"))
  expect_error(ssqda(x, y, lambda_D = 1, lambda_beta = 1),
               "at least 3")

  fixture <- c5qda_fixture()
  expect_warning(
    failed <- ssqda(
      fixture$x, fixture$y, lambda_D = 1, lambda_beta = 1,
      median_max_iter = 1L, median_tol = 1e-15, strict = FALSE
    ),
    "median"
  )
  expect_false(failed$valid)
  expect_identical(failed$diagnostics$failure.stage, "spatial_median")
})


test_that("certified sparse plug-in QDA equals direct log-density algebra", {
  fixture <- c5qda_fixture()
  mean1 <- c(-2, 0.2)
  mean2 <- c(1.7, -0.1)
  covariance1 <- matrix(c(2, 0.3, 0.3, 1), 2)
  covariance2 <- matrix(c(1.2, -0.2, -0.2, 2.5), 2)
  prior <- c(class1 = 0.7, class2 = 0.3)
  fit <- sparse_plugin_qda(
    fixture$x, fixture$y, mean1, mean2, covariance1, covariance2,
    prior = prior
  )
  expect_true(fit$valid)
  expect_true(fit$diagnostics$covariance.certificate$cholesky)
  z <- fixture$x[c(1L, 7L, 12L), , drop = FALSE]
  precision1 <- solve(covariance1)
  precision2 <- solve(covariance2)
  residual1 <- sweep(z, 2L, mean1, "-")
  residual2 <- sweep(z, 2L, mean2, "-")
  direct <- -determinant(covariance1, logarithm = TRUE)$modulus +
    determinant(covariance2, logarithm = TRUE)$modulus -
    rowSums((residual1 %*% precision1) * residual1) +
    rowSums((residual2 %*% precision2) * residual2) +
    2 * log(0.7 / 0.3)
  expect_equal(predict(fit, z, type = "score"), as.numeric(direct),
               tolerance = 2e-13)
  expect_equal(fit$estimate$precision1, precision1, tolerance = 2e-14)
  expect_equal(fit$estimate$precision2, precision2, tolerance = 2e-14)
})


test_that("plug-in QDA preserves feature matching and rejects uncertified SPD", {
  fixture <- c5qda_fixture()
  fit <- sparse_plugin_qda(
    fixture$x, fixture$y, c(-2, 0), c(2, 0), diag(2), diag(c(2, 1))
  )
  z <- fixture$x[1:3, , drop = FALSE]
  reversed <- z[, 2:1, drop = FALSE]
  expect_equal(predict(fit, reversed, type = "score"),
               predict(fit, z, type = "score"), tolerance = 0)

  bad <- matrix(c(1, 2, 2, 1), 2)
  expect_error(
    sparse_plugin_qda(
      fixture$x, fixture$y, c(-2, 0), c(2, 0), bad, diag(2)
    ),
    "positive definite"
  )
  expect_warning(
    invalid <- sparse_plugin_qda(
      fixture$x, fixture$y, c(-2, 0), c(2, 0), bad, diag(2),
      strict = FALSE
    ),
    "positive definite"
  )
  expect_false(invalid$valid)
  expect_identical(invalid$diagnostics$failure.stage,
                   "covariance_certificate")
})


test_that("solver nonconvergence obeys strict and non-strict contracts", {
  fixture <- c5qda_fixture()
  expect_error(
    sdar_qda(
      fixture$x, fixture$y, lambda_D = 0.01, lambda_beta = 0.01,
      solver_tol = 1e-14, solver_max_iter = 1L
    ),
    "did not certify"
  )
  expect_warning(
    invalid <- sdar_qda(
      fixture$x, fixture$y, lambda_D = 0.01, lambda_beta = 0.01,
      solver_tol = 1e-14, solver_max_iter = 1L, strict = FALSE
    ),
    "did not certify"
  )
  expect_false(invalid$valid)
  expect_identical(invalid$diagnostics$failure.stage,
                   "matrix_dantzig_solver")
})


test_that("sparse QDA APIs expose only explicit tuning routes", {
  fixture <- c5qda_fixture()
  expect_error(sdar_qda(fixture$x, fixture$y), "No default")
  expect_error(ssqda(fixture$x, fixture$y), "No default")
  expect_error(
    sdar_qda(fixture$x, fixture$y, lambda_D = 1), "supplied together"
  )
  expect_error(
    li_shao_sparse_qda(
      fixture$x, fixture$y, selection = "paper_bisection"
    ),
    "bisection_tol"
  )
  expect_true(all(c(
    "li_shao_sparse_qda", "jiang_da_qda", "sdar_qda", "ssqda",
    "sparse_plugin_qda"
  ) %in% getNamespaceExports("HDElliptical")))
})
