ch4ind_panel_fixture <- function() {
  structure(
    matrix(c(
      -2.0,  0.5,  1.4, -0.7,
      -1.1,  1.8, -0.4,  2.2,
       0.3, -1.2,  2.5,  0.8,
       1.7,  0.2, -1.6, -2.1,
       2.4,  1.1,  0.6,  1.5,
      -0.6,  2.7, -2.2,  0.4
    ), ncol = 4L, byrow = TRUE),
    dimnames = list(paste0("t", 1:6), paste0("u", 1:4))
  )
}

ch4ind_rank_fixture <- function() {
  list(
    x = structure(
      matrix(c(
        1, 7,
        4, 2,
        2, 6,
        7, 1,
        3, 5,
        6, 3,
        5, 4
      ), ncol = 2L, byrow = TRUE),
      dimnames = list(paste0("row", 1:7), c("x.a", "x.b"))
    ),
    y = structure(
      matrix(c(
        6, 3,
        1, 7,
        5, 2,
        2, 6,
        7, 4,
        3, 1,
        4, 5
      ), ncol = 2L, byrow = TRUE),
      dimnames = list(paste0("row", 1:7), c("y.a", "y.b"))
    )
  )
}

ch4ind_corr_reference <- function(residuals) {
  squared.norms <- colSums(residuals^2)
  crossprod(residuals) / sqrt(outer(squared.norms, squared.norms))
}

ch4ind_extreme_p_reference <- function(statistic, constant) {
  -expm1(-exp(-statistic / 2) * constant)
}

ch4ind_fjlx_reference <- function(panel, regressors = NULL) {
  time <- nrow(panel)
  units <- ncol(panel)
  scaled <- sweep(panel, 2L, apply(abs(panel), 2L, max), "/")
  rank <- 0L
  qs <- vector("list", units)
  if (is.null(regressors)) {
    residuals <- scaled
  } else if (is.list(regressors)) {
    residuals <- matrix(0, time, units)
    for (unit in seq_len(units)) {
      fit <- qr(regressors[[unit]], LAPACK = FALSE)
      rank <- fit$rank
      qs[[unit]] <- qr.Q(fit, complete = FALSE)[, seq_len(rank),
                                                  drop = FALSE]
      residuals[, unit] <- qr.resid(fit, scaled[, unit])
    }
  } else {
    fit <- qr(regressors, LAPACK = FALSE)
    rank <- fit$rank
    q <- qr.Q(fit, complete = FALSE)[, seq_len(rank), drop = FALSE]
    qs <- rep(list(q), units)
    residuals <- qr.resid(fit, scaled)
  }
  correlations <- ch4ind_corr_reference(residuals)
  upper <- correlations[upper.tri(correlations)]
  overlap <- 0
  for (first in seq_len(units - 1L)) for (second in (first + 1L):units) {
    overlap <- overlap + if (rank == 0L) {
      time
    } else {
      time - 2 * rank + sum(crossprod(qs[[first]], qs[[second]])^2)
    }
  }
  residual.df <- time - rank
  mu <- time * overlap / residual.df^2
  sum.raw <- time * sum(upper^2)
  sum.z <- (sum.raw - mu) / units
  maximum <- max(abs(upper))
  max.gumbel <- time * maximum^2 - 4 * log(units) + log(log(units))
  p.max <- ch4ind_extreme_p_reference(
    max.gumbel, 1 / sqrt(8 * pi)
  )
  p.sum <- pnorm(sum.z, lower.tail = FALSE)
  minimum.p <- min(p.max, p.sum)
  list(
    residuals = residuals,
    correlations = correlations,
    overlap = overlap,
    mu = mu,
    sum.raw = sum.raw,
    sum.z = sum.z,
    p.sum = p.sum,
    maximum = maximum,
    max.gumbel = max.gumbel,
    p.max = p.max,
    C.N = minimum.p,
    p.combined = 2 * minimum.p - minimum.p^2
  )
}

ch4ind_serial_reference <- function(panel, temporal) {
  residuals <- panel / max(abs(panel))
  directions <- sweep(
    residuals, 2L, sqrt(colSums(residuals^2)), "/"
  )
  correlations <- crossprod(directions)
  units <- ncol(panel)
  upper <- correlations[upper.tri(correlations)]
  sum.statistic <- sqrt(2 / (units * (units - 1))) * sum(upper)
  variance.sum <- 0
  for (first in seq_len(units - 1L)) for (second in (first + 1L):units) {
    bar <- rowMeans(directions[, -c(first, second), drop = FALSE])
    variance.sum <- variance.sum +
      sum(directions[, second] * (directions[, first] - bar)) *
      sum(directions[, first] * (directions[, second] - bar))
  }
  sum.variance <- 2 * variance.sum / (units * (units - 1))
  sum.z <- sum.statistic / sqrt(sum.variance)
  p.sum <- pnorm(sum.z, lower.tail = FALSE)
  maximum.square <- max(upper^2)
  effective.dimension <- sum(diag(temporal))^2 / sum(temporal^2)
  max.gumbel <- effective.dimension * maximum.square -
    4 * log(units) + log(log(units))
  p.max <- ch4ind_extreme_p_reference(
    max.gumbel, 1 / sqrt(8 * pi)
  )
  fisher <- -2 * (log(p.sum) + log(p.max))
  list(
    residuals = residuals,
    directions = directions,
    correlations = correlations,
    sum.statistic = sum.statistic,
    sum.variance = sum.variance,
    sum.z = sum.z,
    p.sum = p.sum,
    maximum.square = maximum.square,
    effective.dimension = effective.dimension,
    max.gumbel = max.gumbel,
    p.max = p.max,
    fisher = fisher,
    p.fisher = pchisq(fisher, 4, lower.tail = FALSE)
  )
}

ch4ind_temporal_reference <- function(panel, nu) {
  residuals <- panel / max(abs(panel))
  time <- nrow(residuals)
  units <- ncol(residuals)
  centered <- sweep(residuals, 1L, rowMeans(residuals), "-")
  sigma.hat <- tcrossprod(centered) / (units - 1)
  trace.hat <- sum(diag(sigma.hat))
  u.hat <- crossprod(residuals) / trace.hat
  p.hat <- (sum(u.hat^2) - sum(diag(u.hat))^2 / time) / units
  threshold <- nu * sqrt(p.hat * log(time) / units)
  sigma.tilde <- diag(diag(sigma.hat), time)
  for (first in seq_len(time - 1L)) for (second in (first + 1L):time) {
    theta <- sigma.hat[first, second] /
      sqrt(sigma.hat[first, first] * sigma.hat[second, second])
    score <- if (1 - theta^2 <= 0) Inf else abs(theta) / (1 - theta^2)
    if (score >= threshold) {
      sigma.tilde[first, second] <- sigma.hat[first, second]
      sigma.tilde[second, first] <- sigma.hat[first, second]
    }
  }
  list(
    sigma.hat = sigma.hat,
    u.hat = u.hat,
    p.hat = p.hat,
    threshold = threshold,
    sigma.tilde = sigma.tilde,
    effective.dimension = sum(diag(sigma.tilde))^2 / sum(sigma.tilde^2)
  )
}

ch4ind_rank_reference <- function(x, y, measure, permutations) {
  n <- nrow(x)
  comparisons <- ncol(x) * ncol(y)
  null.second <- if (measure == "spearman") {
    1 / (n - 1)
  } else {
    2 * (2 * n + 5) / (9 * n * (n - 1))
  }
  correlations <- cor(x, y, method = measure)
  statistic.from <- function(current.x) {
    sum(cor(current.x, y, method = measure)^2) -
      comparisons * null.second
  }
  permutation.statistics <- apply(
    permutations, 2L,
    function(index) statistic.from(x[index, , drop = FALSE])
  )
  sum.statistic <- sum(correlations^2) - comparisons * null.second
  permutation.variance <- var(permutation.statistics)
  sum.z <- sum.statistic / sqrt(permutation.variance)
  p.sum <- pnorm(sum.z, lower.tail = FALSE)
  maximum <- max(abs(correlations))
  max.gumbel <- maximum^2 / null.second -
    2 * log(comparisons) + log(log(comparisons))
  p.max <- ch4ind_extreme_p_reference(max.gumbel, 1 / sqrt(pi))
  fisher <- -2 * (log(p.sum) + log(p.max))
  list(
    correlations = correlations,
    null.second = null.second,
    maximum = maximum,
    sum.statistic = sum.statistic,
    permutation.statistics = permutation.statistics,
    permutation.mean = mean(permutation.statistics),
    permutation.variance = permutation.variance,
    sum.z = sum.z,
    p.sum = p.sum,
    max.gumbel = max.gumbel,
    p.max = p.max,
    fisher = fisher,
    p.fisher = pchisq(fisher, 4, lower.tail = FALSE)
  )
}


test_that("Gaussian Wilks matches determinant and canonical references", {
  x <- matrix(c(
    1, 4, 2, 1, 3, 7, 4, 2,
    5, 8, 6, 3, 7, 6, 8, 5
  ), ncol = 2L, byrow = TRUE)
  y <- matrix(c(
    2, 7, 8, 3, 5, 1, 1, 6,
    7, 4, 3, 8, 6, 2, 4, 5
  ), ncol = 2L, byrow = TRUE)
  colnames(x) <- c("left.one", "left.two")
  colnames(y) <- c("right.one", "right.two")
  observed <- gaussian_wilks_independence_test(x, y)

  x.scaled <- sweep(x, 2L, apply(abs(x), 2L, max), "/")
  y.scaled <- sweep(y, 2L, apply(abs(y), 2L, max), "/")
  covariance <- cov(cbind(x.scaled, y.scaled))
  p <- ncol(x)
  s.xx <- covariance[seq_len(p), seq_len(p), drop = FALSE]
  s.yy <- covariance[p + seq_len(ncol(y)), p + seq_len(ncol(y)),
                     drop = FALSE]
  lambda <- det(covariance) / (det(s.xx) * det(s.yy))
  multiplier <- nrow(x) - 1 - (ncol(x) + ncol(y) + 1) / 2
  statistic <- -multiplier * log(lambda)

  expect_s3_class(observed, "hd_independence_test")
  expect_equal(unname(observed$estimate), lambda, tolerance = 1e-12)
  expect_equal(unname(observed$statistic), statistic, tolerance = 1e-12)
  expect_equal(observed$p.value,
               pchisq(statistic, ncol(x) * ncol(y), lower.tail = FALSE),
               tolerance = 1e-12)
  expect_equal(observed$components$canonical_correlations,
               cancor(x.scaled, y.scaled)$cor, tolerance = 1e-12,
               ignore_attr = TRUE)
  expect_equal(observed$components$covariance_x, s.xx,
               tolerance = 1e-12, ignore_attr = TRUE)
  expect_equal(observed$components$covariance_y, s.yy,
               tolerance = 1e-12, ignore_attr = TRUE)
  expect_equal(observed$components$cross_covariance,
               covariance[seq_len(p), p + seq_len(ncol(y)), drop = FALSE],
               tolerance = 1e-12, ignore_attr = TRUE)
  expect_identical(dimnames(observed$components$covariance_x),
                   list(colnames(x), colnames(x)))
  expect_identical(unname(observed$parameter), 4L)
})

test_that("Gaussian Wilks is block-affine invariant and fails singularly", {
  x <- matrix(c(
    1, 4, 2, 1, 3, 7, 4, 2,
    5, 8, 6, 3, 7, 6, 8, 5
  ), ncol = 2L, byrow = TRUE)
  y <- matrix(c(
    2, 7, 8, 3, 5, 1, 1, 6,
    7, 4, 3, 8, 6, 2, 4, 5
  ), ncol = 2L, byrow = TRUE)
  transform.x <- matrix(c(2, 0.4, -0.3, 1.5), 2L, 2L)
  transform.y <- matrix(c(1.2, -0.2, 0.7, 1.8), 2L, 2L)
  baseline <- gaussian_wilks_independence_test(x, y)
  transformed <- gaussian_wilks_independence_test(
    sweep(x %*% transform.x, 2L, c(10, -4), "+"),
    sweep(y %*% transform.y, 2L, c(-3, 7), "+")
  )
  expect_equal(transformed$components$lambda,
               baseline$components$lambda, tolerance = 1e-11)
  expect_equal(transformed$statistic, baseline$statistic,
               tolerance = 1e-10)
  expect_error(
    gaussian_wilks_independence_test(cbind(1:6, 2 * (1:6)),
                                     matrix(c(2, 1, 4, 3, 6, 5), 6L, 1L)),
    "not positive definite", fixed = TRUE
  )
  expect_error(
    gaussian_wilks_independence_test(matrix(1:12, 3L, 4L), 1:3),
    "n > p + q", fixed = TRUE
  )
  expect_error(
    gaussian_wilks_independence_test(matrix(1:8, 4L, 2L), 1:3),
    "same number of rows", fixed = TRUE
  )
})

test_that("Pesaran CD matches the uncentered residual formula", {
  residuals <- ch4ind_panel_fixture()
  correlations <- ch4ind_corr_reference(residuals)
  upper <- correlations[upper.tri(correlations)]
  time <- nrow(residuals)
  units <- ncol(residuals)
  statistic <- sqrt(2 * time / (units * (units - 1))) * sum(upper)
  observed <- pesaran_cd_test(residuals, keep_correlations = TRUE)

  expect_equal(unname(observed$statistic), statistic, tolerance = 1e-13)
  expect_equal(observed$p.value, 2 * pnorm(abs(statistic), lower.tail = FALSE),
               tolerance = 1e-13)
  expect_equal(observed$components$correlations, correlations,
               tolerance = 1e-13, ignore_attr = TRUE)
  expect_identical(dimnames(observed$components$correlations),
                   list(colnames(residuals), colnames(residuals)))
  expect_equal(observed$components$sum, sum(upper), tolerance = 1e-13)
  expect_equal(unname(observed$estimate), mean(upper), tolerance = 1e-13)
  expect_false(observed$diagnostics$residuals.centered.by.function)
  expect_null(pesaran_cd_test(residuals)$components$correlations)

  positive.scaled <- sweep(residuals, 2L, c(2, 0.5, 3, 4), "*")
  permuted <- residuals[, c(4, 2, 1, 3)]
  expect_equal(pesaran_cd_test(positive.scaled)$statistic,
               observed$statistic, tolerance = 1e-13)
  expect_equal(pesaran_cd_test(permuted)$statistic,
               observed$statistic, tolerance = 1e-13)
})

test_that("Pesaran CD validates residual and retention contracts", {
  expect_error(pesaran_cd_test(matrix(0, 4L, 3L)),
               "positive finite scale", fixed = TRUE)
  expect_error(pesaran_cd_test(matrix(1:8, 4L, 2L), keep_correlations = NA),
               "must be TRUE or FALSE", fixed = TRUE)
  expect_error(pesaran_cd_test(matrix(c(1:7, NA), 4L, 2L)),
               "finite values", fixed = TRUE)
  expect_error(pesaran_cd_test(matrix(1:4, 4L, 1L)),
               "at least 2 rows and 2 columns", fixed = TRUE)
})

test_that("Feng-Jiang-Liu-Xiong matches every primary component", {
  panel <- ch4ind_panel_fixture()
  reference <- ch4ind_fjlx_reference(panel)
  observed <- feng_jiang_liu_xiong_panel_independence_test(
    panel, keep_correlations = TRUE
  )

  expect_equal(observed$raw.statistic[["S.N"]], reference$sum.raw,
               tolerance = 1e-12)
  expect_equal(observed$raw.statistic[["L.N"]], reference$maximum,
               tolerance = 1e-12)
  expect_equal(observed$components$projection.overlap.sum,
               reference$overlap, tolerance = 1e-12)
  expect_equal(observed$components$mu, reference$mu, tolerance = 1e-12)
  expect_equal(observed$components$sum.z, reference$sum.z,
               tolerance = 1e-12)
  expect_equal(observed$components$p.sum, reference$p.sum,
               tolerance = 1e-12)
  expect_equal(observed$components$max.gumbel, reference$max.gumbel,
               tolerance = 1e-12)
  expect_equal(observed$components$p.max, reference$p.max,
               tolerance = 1e-12)
  expect_equal(observed$components$C.N, reference$C.N,
               tolerance = 1e-12)
  expect_equal(observed$p.value, reference$p.combined, tolerance = 1e-12)
  expect_equal(observed$components$correlations, reference$correlations,
               tolerance = 1e-12, ignore_attr = TRUE)
  expect_identical(observed$diagnostics$serial.correlation.allowed, FALSE)

  maximum <- feng_jiang_liu_xiong_panel_independence_test(
    panel, component = "max"
  )
  sum.test <- feng_jiang_liu_xiong_panel_independence_test(
    panel, component = "sum"
  )
  expect_equal(unname(maximum$statistic), reference$max.gumbel,
               tolerance = 1e-12)
  expect_equal(maximum$p.value, reference$p.max, tolerance = 1e-12)
  expect_equal(unname(sum.test$statistic), reference$sum.z,
               tolerance = 1e-12)
  expect_equal(sum.test$p.value, reference$p.sum, tolerance = 1e-12)
})

test_that("Feng-Jiang-Liu-Xiong uses exact common and unit projections", {
  panel <- ch4ind_panel_fixture()
  time.index <- seq_len(nrow(panel)) - mean(seq_len(nrow(panel)))
  common <- cbind(intercept = 1, trend = time.index)
  common.reference <- ch4ind_fjlx_reference(panel, common)
  common.observed <- feng_jiang_liu_xiong_panel_independence_test(
    panel, common
  )
  expect_equal(common.observed$components$projection.overlap.sum,
               common.reference$overlap, tolerance = 1e-12)
  expect_equal(common.observed$components$mu, common.reference$mu,
               tolerance = 1e-12)
  expect_equal(common.observed$components$sum.raw, common.reference$sum.raw,
               tolerance = 1e-12)
  expect_equal(common.observed$components$sum.z, common.reference$sum.z,
               tolerance = 1e-12)

  designs <- list(
    cbind(1, c(-3, -2, -1, 1, 2, 3)),
    cbind(1, c(-2, -1, 2, 3, 1, -3)),
    cbind(1, c(1, -3, 2, -1, -2, 3)),
    cbind(1, c(-1, 3, -3, 2, 1, -2))
  )
  unit.reference <- ch4ind_fjlx_reference(panel, designs)
  unit.observed <- feng_jiang_liu_xiong_panel_independence_test(
    panel, designs
  )
  expect_equal(unit.observed$components$projection.overlap.sum,
               unit.reference$overlap, tolerance = 1e-12)
  expect_equal(unit.observed$components$mu, unit.reference$mu,
               tolerance = 1e-12)
  expect_equal(unit.observed$components$sum.raw, unit.reference$sum.raw,
               tolerance = 1e-12)
  expect_equal(unit.observed$components$sum.z, unit.reference$sum.z,
               tolerance = 1e-12)
  expect_match(unit.observed$diagnostics$design.route, "unit-specific",
               fixed = TRUE)

  signed.scaled <- sweep(panel, 2L, c(-2, 0.5, 3, -4), "*")
  transformed <- feng_jiang_liu_xiong_panel_independence_test(signed.scaled)
  baseline <- feng_jiang_liu_xiong_panel_independence_test(panel)
  expect_equal(transformed$raw.statistic, baseline$raw.statistic,
               tolerance = 1e-12)
  expect_equal(transformed$p.value, baseline$p.value, tolerance = 1e-12)
})

test_that("Feng-Jiang-Liu-Xiong rejects undefined design cases", {
  panel <- ch4ind_panel_fixture()
  expect_error(
    feng_jiang_liu_xiong_panel_independence_test(
      panel, regressors = list(diag(6), diag(6), diag(6), diag(6))
    ),
    "fewer than T columns", fixed = TRUE
  )
  expect_error(
    feng_jiang_liu_xiong_panel_independence_test(
      panel, regressors = cbind(rep(1, nrow(panel)), rep(1, nrow(panel)))
    ),
    "full column rank", fixed = TRUE
  )
  expect_error(
    feng_jiang_liu_xiong_panel_independence_test(
      panel, regressors = list(matrix(1, 6L, 1L),
                              cbind(1, seq_len(6)),
                              matrix(1, 6L, 1L),
                              matrix(1, 6L, 1L))
    ),
    "same column rank", fixed = TRUE
  )
  expect_error(
    feng_jiang_liu_xiong_panel_independence_test(
      panel, regressors = list(matrix(1, 6L, 1L))
    ),
    "one T by p matrix per panel unit", fixed = TRUE
  )
  expect_error(
    feng_jiang_liu_xiong_panel_independence_test(
      panel, component = "unknown"
    ),
    "arg", fixed = TRUE
  )
})

test_that("serial-panel test matches signed sum and supplied covariance", {
  panel <- ch4ind_panel_fixture()
  temporal <- toeplitz(0.45^(0:(nrow(panel) - 1L)))
  reference <- ch4ind_serial_reference(panel, temporal)
  observed <- wang_liu_feng_ma_serial_panel_test(
    panel, temporal_covariance = temporal, keep_matrices = TRUE
  )

  expect_s3_class(observed, "hd_independence_test")
  expect_equal(observed$components$correlations, reference$correlations,
               tolerance = 1e-12, ignore_attr = TRUE)
  expect_equal(observed$raw.statistic[["S.N"]], reference$sum.statistic,
               tolerance = 1e-12)
  expect_equal(observed$components$sum_variance, reference$sum.variance,
               tolerance = 1e-12)
  expect_equal(observed$components$sum.z, reference$sum.z,
               tolerance = 1e-12)
  expect_equal(observed$components$p.sum, reference$p.sum,
               tolerance = 1e-12)
  expect_equal(observed$raw.statistic[["L.N"]], reference$maximum.square,
               tolerance = 1e-12)
  expect_equal(observed$components$effective_dimension,
               reference$effective.dimension, tolerance = 1e-12)
  expect_equal(observed$components$max.gumbel, reference$max.gumbel,
               tolerance = 1e-12)
  expect_equal(observed$components$p.max, reference$p.max,
               tolerance = 1e-12)
  expect_equal(unname(observed$statistic), reference$fisher,
               tolerance = 1e-12)
  expect_equal(observed$p.value, reference$p.fisher, tolerance = 1e-12)
  expect_identical(unname(observed$parameter), 4)
  expect_equal(observed$components$sigma_tilde, temporal,
               tolerance = 1e-14, ignore_attr = TRUE)
  expect_null(observed$components$sigma_hat)
  expect_null(observed$components$u_hat)
  expect_true(is.na(observed$components$p_hat))
  expect_true(is.na(observed$components$threshold))
  expect_match(observed$diagnostics$temporal.route,
               "supplied positive-definite", fixed = TRUE)

  maximum <- wang_liu_feng_ma_serial_panel_test(
    panel, temporal_covariance = temporal, component = "max"
  )
  sum.test <- wang_liu_feng_ma_serial_panel_test(
    panel, temporal_covariance = temporal, component = "sum"
  )
  expect_equal(unname(maximum$statistic), reference$max.gumbel,
               tolerance = 1e-12)
  expect_equal(maximum$p.value, reference$p.max, tolerance = 1e-12)
  expect_equal(unname(sum.test$statistic), reference$sum.z,
               tolerance = 1e-12)
  expect_equal(sum.test$p.value, reference$p.sum, tolerance = 1e-12)
})

test_that("literal serial temporal estimator matches every threshold step", {
  panel <- ch4ind_panel_fixture()
  nu <- 1.51
  reference <- ch4ind_temporal_reference(panel, nu)
  observed <- wang_liu_feng_ma_serial_panel_test(
    panel, nu = nu, keep_matrices = TRUE
  )

  expect_equal(observed$components$sigma_hat, reference$sigma.hat,
               tolerance = 1e-12, ignore_attr = TRUE)
  expect_equal(observed$components$u_hat, reference$u.hat,
               tolerance = 1e-12, ignore_attr = TRUE)
  expect_equal(observed$components$p_hat, reference$p.hat,
               tolerance = 1e-12)
  expect_equal(observed$components$threshold, reference$threshold,
               tolerance = 1e-12)
  expect_equal(observed$components$sigma_tilde, reference$sigma.tilde,
               tolerance = 1e-12, ignore_attr = TRUE)
  expect_equal(observed$components$effective_dimension,
               reference$effective.dimension, tolerance = 1e-12)
  expect_identical(dimnames(observed$components$sigma_hat),
                   list(rownames(panel), rownames(panel)))
  expect_identical(dimnames(observed$components$u_hat),
                   list(colnames(panel), colnames(panel)))
  expect_match(observed$diagnostics$temporal.route,
               "literal primary", fixed = TRUE)
  expect_equal(observed$diagnostics$threshold.nu, nu)

  compact <- wang_liu_feng_ma_serial_panel_test(panel)
  expect_null(compact$components$correlations)
  expect_null(compact$components$sigma_hat)
  expect_null(compact$components$sigma_tilde)
  expect_null(compact$components$u_hat)
})

test_that("serial-panel aggregate is invariant to legitimate re-expression", {
  panel <- ch4ind_panel_fixture()
  temporal <- toeplitz(0.35^(0:(nrow(panel) - 1L)))
  baseline <- wang_liu_feng_ma_serial_panel_test(
    panel, temporal_covariance = temporal
  )
  rescaled <- wang_liu_feng_ma_serial_panel_test(
    -7.5 * panel, temporal_covariance = temporal
  )
  permuted <- wang_liu_feng_ma_serial_panel_test(
    panel[, c(4, 1, 3, 2)], temporal_covariance = temporal
  )
  expect_equal(rescaled$raw.statistic, baseline$raw.statistic,
               tolerance = 1e-12)
  expect_equal(rescaled$p.value, baseline$p.value, tolerance = 1e-12)
  expect_equal(permuted$raw.statistic, baseline$raw.statistic,
               tolerance = 1e-12)
  expect_equal(permuted$p.value, baseline$p.value, tolerance = 1e-12)

  common <- cbind(1, seq_len(nrow(panel)))
  projected <- wang_liu_feng_ma_serial_panel_test(
    panel, regressors = common, temporal_covariance = temporal
  )
  manual.residuals <- qr.resid(
    qr(common, LAPACK = FALSE), panel / max(abs(panel))
  )
  manual <- wang_liu_feng_ma_serial_panel_test(
    manual.residuals, temporal_covariance = temporal
  )
  expect_equal(projected$raw.statistic, manual$raw.statistic,
               tolerance = 1e-12)
  expect_equal(projected$p.value, manual$p.value, tolerance = 1e-12)
})

test_that("serial-panel fails rather than repairing undefined estimators", {
  panel <- ch4ind_panel_fixture()
  time <- nrow(panel)
  expect_error(
    wang_liu_feng_ma_serial_panel_test(panel[, 1:2],
                                      temporal_covariance = diag(time)),
    "at least three panel units", fixed = TRUE
  )
  expect_error(
    wang_liu_feng_ma_serial_panel_test(panel, nu = sqrt(2)),
    "strictly greater than sqrt(2)", fixed = TRUE
  )
  expect_error(
    wang_liu_feng_ma_serial_panel_test(
      panel, temporal_covariance = diag(time - 1L)
    ),
    "at least 6 rows and 6 columns", fixed = TRUE
  )
  asymmetric <- diag(time)
  asymmetric[1, 2] <- 0.2
  expect_error(
    wang_liu_feng_ma_serial_panel_test(
      panel, temporal_covariance = asymmetric
    ),
    "symmetric T by T", fixed = TRUE
  )
  semidefinite <- diag(c(rep(1, time - 1L), 0))
  expect_error(
    wang_liu_feng_ma_serial_panel_test(
      panel, temporal_covariance = semidefinite
    ),
    "positive definite", fixed = TRUE
  )
  repeated <- matrix(rep(c(1, 2, 4, 3, 5, 6), 4L), 6L, 4L)
  expect_error(
    wang_liu_feng_ma_serial_panel_test(
      repeated, temporal_covariance = diag(time)
    ),
    "sum-variance estimate is not strictly positive", fixed = TRUE
  )
  spherical <- rbind(
    c(1, 1, -1, -1),
    c(1, -1, 1, -1)
  )
  expect_error(
    wang_liu_feng_ma_serial_panel_test(spherical),
    "P-hat_N is not positive", fixed = TRUE
  )
})

test_that("Spearman vector test matches exact ranks and intrinsic permutations", {
  fixture <- ch4ind_rank_fixture()
  B <- 29L
  caller.seed <- 812L
  permutation.seed <- 271L
  set.seed(caller.seed)
  caller.state <- .Random.seed
  observed <- wang_liu_feng_vector_independence_test(
    fixture$x, fixture$y, measure = "spearman", B = B,
    seed = permutation.seed, keep_correlations = TRUE,
    keep_permutation = TRUE
  )
  expect_identical(.Random.seed, caller.state)

  set.seed(permutation.seed)
  permutations <- replicate(B, sample.int(nrow(fixture$x)))
  storage.mode(permutations) <- "integer"
  reference <- ch4ind_rank_reference(
    fixture$x, fixture$y, "spearman", permutations
  )
  expect_equal(observed$components$correlations, reference$correlations,
               tolerance = 1e-13, ignore_attr = TRUE)
  expect_equal(observed$components$null_second_moment,
               reference$null.second, tolerance = 1e-14)
  expect_equal(observed$raw.statistic[["maximum.absolute.rank.correlation"]],
               reference$maximum, tolerance = 1e-13)
  expect_equal(observed$raw.statistic[["centered.sum.squares"]],
               reference$sum.statistic, tolerance = 1e-13)
  expect_equal(observed$components$permutation_statistics,
               reference$permutation.statistics, tolerance = 1e-13,
               ignore_attr = TRUE)
  expect_equal(observed$components$permutation_mean,
               reference$permutation.mean, tolerance = 1e-13)
  expect_equal(observed$components$permutation_variance,
               reference$permutation.variance, tolerance = 1e-13)
  expect_equal(observed$components$sum.z, reference$sum.z,
               tolerance = 1e-13)
  expect_equal(observed$components$p.sum, reference$p.sum,
               tolerance = 1e-13)
  expect_equal(observed$components$max.gumbel, reference$max.gumbel,
               tolerance = 1e-13)
  expect_equal(observed$components$p.max, reference$p.max,
               tolerance = 1e-13)
  expect_equal(unname(observed$statistic), reference$fisher,
               tolerance = 1e-13)
  expect_equal(observed$p.value, reference$p.fisher, tolerance = 1e-13)
  expect_identical(unname(observed$components$permutation.indices),
                   unname(permutations))
  expect_identical(dimnames(observed$components$correlations),
                   list(colnames(fixture$x), colnames(fixture$y)))
  expect_identical(unname(observed$parameter), 4)
})

test_that("Kendall vector test matches pair-sign reference and routing", {
  fixture <- ch4ind_rank_fixture()
  B <- 31L
  seed <- 99L
  set.seed(seed)
  permutations <- replicate(B, sample.int(nrow(fixture$x)))
  storage.mode(permutations) <- "integer"
  reference <- ch4ind_rank_reference(
    fixture$x, fixture$y, "kendall", permutations
  )
  observed <- wang_liu_feng_vector_independence_test(
    fixture$x, fixture$y, measure = "kendall", B = B, seed = seed,
    keep_correlations = TRUE, keep_permutation = TRUE
  )
  expect_equal(observed$components$correlations, reference$correlations,
               tolerance = 1e-13, ignore_attr = TRUE)
  expect_equal(observed$components$null_second_moment,
               reference$null.second, tolerance = 1e-14)
  expect_equal(observed$components$permutation_statistics,
               reference$permutation.statistics, tolerance = 1e-13,
               ignore_attr = TRUE)
  expect_equal(observed$components$permutation_variance,
               reference$permutation.variance, tolerance = 1e-13)
  expect_equal(observed$components$sum.z, reference$sum.z,
               tolerance = 1e-13)
  expect_equal(observed$components$p.sum, reference$p.sum,
               tolerance = 1e-13)
  expect_equal(observed$components$max.gumbel, reference$max.gumbel,
               tolerance = 1e-13)
  expect_equal(observed$components$p.max, reference$p.max,
               tolerance = 1e-13)
  expect_equal(observed$p.value, reference$p.fisher, tolerance = 1e-13)

  maximum <- wang_liu_feng_vector_independence_test(
    fixture$x, fixture$y, measure = "kendall", component = "max",
    B = B, seed = seed
  )
  sum.test <- wang_liu_feng_vector_independence_test(
    fixture$x, fixture$y, measure = "kendall", component = "sum",
    B = B, seed = seed
  )
  expect_equal(unname(maximum$statistic), reference$max.gumbel,
               tolerance = 1e-13)
  expect_equal(maximum$p.value, reference$p.max, tolerance = 1e-13)
  expect_equal(unname(sum.test$statistic), reference$sum.z,
               tolerance = 1e-13)
  expect_equal(sum.test$p.value, reference$p.sum, tolerance = 1e-13)
  expect_length(observed$diagnostics$review.only, 2L)
})

test_that("rank vector test is monotone invariant and obeys RNG contracts", {
  fixture <- ch4ind_rank_fixture()
  baseline <- wang_liu_feng_vector_independence_test(
    fixture$x, fixture$y, B = 23L, seed = 44L
  )
  transformed <- wang_liu_feng_vector_independence_test(
    fixture$x^3, log(fixture$y), B = 23L, seed = 44L
  )
  permuted <- wang_liu_feng_vector_independence_test(
    fixture$x[, 2:1], fixture$y[, 2:1], B = 23L, seed = 44L
  )
  expect_equal(transformed$raw.statistic, baseline$raw.statistic,
               tolerance = 1e-13)
  expect_equal(transformed$p.value, baseline$p.value, tolerance = 1e-13)
  expect_equal(permuted$raw.statistic, baseline$raw.statistic,
               tolerance = 1e-13)
  expect_equal(permuted$p.value, baseline$p.value, tolerance = 1e-13)

  set.seed(733L)
  invisible(wang_liu_feng_vector_independence_test(
    fixture$x, fixture$y, B = 17L, seed = NULL
  ))
  observed.state <- .Random.seed
  set.seed(733L)
  invisible(replicate(17L, sample.int(nrow(fixture$x))))
  expected.state <- .Random.seed
  expect_identical(observed.state, expected.state)

  compact <- wang_liu_feng_vector_independence_test(
    fixture$x, fixture$y, B = 17L, seed = 8L
  )
  expect_null(compact$components$correlations)
  expect_null(compact$components$permutation_statistics)
  expect_null(compact$components$permutation.indices)
})

test_that("rank vector calibration rejects ties and malformed controls", {
  fixture <- ch4ind_rank_fixture()
  tied <- fixture$x
  tied[2, 1] <- tied[1, 1]
  expect_error(
    wang_liu_feng_vector_independence_test(
      tied, fixture$y, B = 11L, seed = 1L
    ),
    "exact tie", fixed = TRUE
  )
  expect_error(
    wang_liu_feng_vector_independence_test(
      fixture$x[, 1], fixture$y[, 1], B = 11L, seed = 1L
    ),
    "p * q >= 2", fixed = TRUE
  )
  expect_error(
    wang_liu_feng_vector_independence_test(
      fixture$x, fixture$y[-1, ], B = 11L, seed = 1L
    ),
    "same number of rows", fixed = TRUE
  )
  expect_error(
    wang_liu_feng_vector_independence_test(
      fixture$x, fixture$y, B = 1L, seed = 1L
    ),
    "between 2", fixed = TRUE
  )
  expect_error(
    wang_liu_feng_vector_independence_test(
      fixture$x, fixture$y, B = 11L, seed = -1L
    ),
    "between 0", fixed = TRUE
  )
  expect_error(
    wang_liu_feng_vector_independence_test(
      fixture$x, fixture$y, B = 11L, seed = 1L,
      keep_permutation = 1
    ),
    "must be TRUE or FALSE", fixed = TRUE
  )
  expect_error(
    wang_liu_feng_vector_independence_test(
      fixture$x, fixture$y, measure = "pearson", B = 11L, seed = 1L
    ),
    "arg", fixed = TRUE
  )
})
