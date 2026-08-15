afc_fixture <- function(scale = 1) {
  tt <- seq_len(20)
  u <- (tt - 10.5) / 10
  basis <- cbind(linear = u, quadratic = u^2)
  factors <- cbind(market = sin(tt / 4) + 0.03 * tt)
  error <- cbind(
    0.7 * sin(0.71 * tt) + 0.21 * cos(0.19 * tt),
    0.6 * cos(0.47 * tt) - 0.17 * sin(0.83 * tt),
    0.8 * sin(0.31 * tt + 0.4) + 0.13 * cos(0.97 * tt),
    0.5 * cos(0.29 * tt - 0.2) + 0.19 * sin(1.07 * tt)
  )
  alpha_curve <- cbind(
    0.08 + 0.04 * u,
    -0.03 + 0.02 * u^2,
    0.01 - 0.03 * u,
    0.06 + 0.015 * u^2
  )
  beta_curve <- cbind(
    0.5 + 0.1 * u,
    -0.3 + 0.05 * u,
    0.2 - 0.08 * u,
    0.4 + 0.03 * u^2
  )
  returns <- scale * (
    error + alpha_curve + as.numeric(factors[, 1L]) * beta_curve
  )
  colnames(returns) <- c("A", "B", "C", "D")
  rownames(returns) <- paste0("t", tt)
  list(returns = returns, factors = factors, basis = basis)
}

afc_fit_pair <- function(scale = 1) {
  x <- afc_fixture(scale)
  null.design <- conditional_alpha_sieve_design(x$basis, x$factors)
  trace.design <- conditional_alpha_sieve_design(
    x$basis, x$factors, center_alpha = FALSE
  )
  list(
    data = x,
    null.design = null.design,
    trace.design = trace.design,
    fit = conditional_alpha_sieve_fit(x$returns, null.design),
    trace.fit = conditional_alpha_sieve_fit(x$returns, trace.design)
  )
}

afc_project_r <- function(y, z) {
  coefficient <- solve(crossprod(z), crossprod(z, y))
  one.coefficient <- solve(crossprod(z), crossprod(z, rep(1, nrow(z))))
  list(
    residuals = y - z %*% coefficient,
    h = as.numeric(1 - z %*% one.coefficient),
    coefficient = coefficient
  )
}

afc_light_r <- function(fit, factor.count) {
  e <- fit$residuals.scaled
  h <- fit$h
  T <- nrow(e)
  N <- ncol(e)
  q <- fit$design.columns
  S <- sum(colSums(e)^2) / (N * T)
  mu <- sum(e^2 * h^2) / (N * T)
  centered <- sweep(e, 2L, colMeans(e), "-")
  sigma <- crossprod(centered) / T
  trace.raw <- sum(sigma^2)
  trace.one <- sum(diag(sigma))
  trace.hat <- T^2 / ((T + q - 1) * (T - q)) *
    (trace.raw - trace.one^2 / (T - q))
  h.off <- sum(outer(h^2, h^2)) - sum(h^4)
  variance <- 2 * trace.hat * h.off / (N^2 * T^2)
  marginal <- colSums(e^2) / (T - factor.count - 1)
  coordinate <- colSums(e)^2 / (T * marginal)
  list(
    S = S, mu = mu, trace.raw = trace.raw, trace.one = trace.one,
    trace.hat = trace.hat, h.off = h.off, variance = variance,
    z = (S - mu) / sqrt(variance), marginal = marginal,
    coordinate = coordinate, maximum = max(coordinate)
  )
}

afc_sign_rows <- function(x) {
  norms <- sqrt(rowSums(x^2))
  answer <- matrix(0, nrow(x), ncol(x))
  keep <- norms > 0
  answer[keep, ] <- x[keep, , drop = FALSE] / norms[keep]
  answer
}

afc_css_r <- function(fit, trace.fit) {
  u <- afc_sign_rows(fit$residuals.scaled)
  utilde <- afc_sign_rows(trace.fit$residuals.scaled)
  h <- fit$h
  h2 <- sum(h^2)
  quadratic <- sum(drop(crossprod(u, h))^2) / h2
  gram <- tcrossprod(utilde)
  weights <- outer(h^2, h^2)
  trace.numerator <- sum(weights[row(weights) != col(weights)] *
                           gram[row(gram) != col(gram)]^2)
  trace.denominator <- h2 * (h2 - 1)
  trace.hat <- trace.numerator / trace.denominator
  list(
    quadratic = quadratic,
    numerator = quadratic - 1,
    trace.numerator = trace.numerator,
    trace.denominator = trace.denominator,
    trace.hat = trace.hat,
    statistic = (quadratic - 1) / sqrt(trace.hat)
  )
}

afc_kendall_r <- function(x) {
  p <- ncol(x)
  answer <- matrix(0, p, p)
  count <- 0
  for (i in seq_len(nrow(x) - 1L)) {
    for (j in seq.int(i + 1L, nrow(x))) {
      difference <- x[i, ] - x[j, ]
      norm <- sqrt(sum(difference^2))
      direction <- if (norm == 0) numeric(p) else difference / norm
      answer <- answer + tcrossprod(direction)
      count <- count + 1L
    }
  }
  answer / count
}


test_that("conditional sieve design follows the literal block construction", {
  tt <- seq_len(8)
  basis <- cbind(tt / 8, (tt / 8)^2)
  factors <- cbind(market = (-1)^tt / (tt + 1))
  observed <- conditional_alpha_sieve_design(basis, factors)
  expected <- cbind(
    sweep(basis, 2L, colMeans(basis), "-"),
    basis * factors[, 1]
  )

  expect_equal(unclass(observed), unclass(expected), ignore_attr = TRUE)
  expect_equal(attr(observed, "basis.columns"), 2L)
  expect_equal(attr(observed, "factor.count"), 1L)
  expect_true(attr(observed, "center.alpha"))

  raw <- conditional_alpha_sieve_design(
    basis, factors, center_alpha = FALSE
  )
  expect_equal(unclass(raw)[, 1:2], basis, ignore_attr = TRUE)
  expect_false(attr(raw, "center.alpha"))
})


test_that("spline partition dependence is never silently repaired", {
  tt <- seq(0, 1, length.out = 10)
  partition <- cbind(1 - tt, tt)
  expect_error(
    conditional_alpha_sieve_design(partition),
    "full column rank"
  )
  contrast <- matrix(c(1, -1), 2L, 1L)
  observed <- conditional_alpha_sieve_design(
    partition, alpha_contrast = contrast
  )
  expected <- sweep(partition, 2L, colMeans(partition), "-") %*%
    contrast
  expect_equal(unclass(observed), unclass(expected), ignore_attr = TRUE)
  expect_equal(ncol(observed), 1L)
})


test_that("conditional sieve fit matches strict pure-R projection", {
  x <- afc_fixture()
  design <- conditional_alpha_sieve_design(x$basis, x$factors)
  observed <- conditional_alpha_sieve_fit(x$returns, design)
  reference <- afc_project_r(x$returns, unclass(design))

  expect_s3_class(observed, "conditional_alpha_sieve_fit")
  expect_equal(observed$residuals, reference$residuals, tolerance = 1e-12)
  expect_equal(observed$h, reference$h, tolerance = 1e-12,
               ignore_attr = TRUE)
  expect_equal(observed$h2, sum(reference$h^2), tolerance = 1e-12)
  expect_equal(observed$factor.count, 1L)
  expect_equal(observed$design.columns, ncol(design))
  expect_lt(observed$diagnostics$normal.equation.residual, 1e-12)
  expect_identical(observed$diagnostics$generalized.inverse, "none")
})


test_that("HDA feasible centering and trace match the literal formulas", {
  x <- afc_fit_pair()
  reference <- afc_light_r(x$fit, factor.count = 1L)
  observed <- ma_lan_su_tsai_conditional_alpha_sum_test(x$fit)

  expect_s3_class(observed, "htest")
  expect_equal(unname(observed$statistic), reference$z,
               tolerance = 2e-11)
  expect_equal(observed$components$S.NT.scaled, reference$S,
               tolerance = 2e-12)
  expect_equal(observed$components$mu.hat.scaled, reference$mu,
               tolerance = 2e-12)
  expect_equal(
    observed$components$trace.Sigma.squared.raw,
    reference$trace.raw, tolerance = 2e-11
  )
  expect_equal(
    observed$components$trace.Sigma.squared.estimate,
    reference$trace.hat, tolerance = 2e-11
  )
  expect_equal(observed$components$variance.hat.scaled,
               reference$variance, tolerance = 2e-11)
  expect_equal(observed$p.value,
               pnorm(reference$z, lower.tail = FALSE), tolerance = 1e-14)
  expect_true(observed$diagnostics$book.oracle.standardization.conflict)
})


test_that("conditional maximum and Fisher adaptive tests match primary formulas", {
  x <- afc_fit_pair()
  reference <- afc_light_r(x$fit, factor.count = 1L)
  maximum <- ma_feng_wang_bao_conditional_alpha_test(
    x$fit, component = "max"
  )
  centered <- reference$maximum - 2 * log(x$fit$N) +
    log(log(x$fit$N))
  intensity <- exp(-0.5 * log(pi) - centered / 2)
  p.max <- -expm1(-intensity)

  expect_equal(maximum$components$coordinate.t.squared,
               reference$coordinate, tolerance = 2e-12,
               ignore_attr = TRUE)
  expect_equal(maximum$components$marginal.variance.scaled,
               reference$marginal, tolerance = 2e-12,
               ignore_attr = TRUE)
  expect_equal(maximum$components$marginal.df, x$fit$T - 2)
  expect_equal(unname(maximum$statistic), centered, tolerance = 2e-12)
  expect_equal(maximum$p.value, p.max, tolerance = 2e-14)

  adaptive <- ma_feng_wang_bao_conditional_alpha_test(
    x$fit, component = "adaptive"
  )
  p.sum <- pnorm(reference$z, lower.tail = FALSE)
  fisher <- -2 * (log(p.max) + log(p.sum))
  expect_equal(unname(adaptive$statistic), fisher, tolerance = 2e-12)
  expect_equal(adaptive$p.value,
               pchisq(fisher, 4, lower.tail = FALSE), tolerance = 2e-14)
  expect_identical(
    adaptive$diagnostics$adaptive.calibration,
    "primary Fisher chi-square with 4 df"
  )
  expect_true(adaptive$diagnostics$book.cauchy.attribution.conflict)
})


test_that("CSS statistic and off-diagonal trace match the literal formula", {
  x <- afc_fit_pair()
  reference <- afc_css_r(x$fit, x$trace.fit)
  observed <- zhao_conditional_spatial_sign_sum_test(
    x$fit, x$trace.fit
  )

  expect_equal(observed$components$sign.quadratic,
               reference$quadratic, tolerance = 2e-12)
  expect_equal(observed$components$centered.numerator,
               reference$numerator, tolerance = 2e-12)
  expect_equal(observed$components$trace.numerator,
               reference$trace.numerator, tolerance = 2e-12)
  expect_equal(observed$components$trace.denominator,
               reference$trace.denominator, tolerance = 2e-12)
  expect_equal(observed$components$trace.Sigma.u.squared,
               reference$trace.hat, tolerance = 2e-12)
  expect_equal(unname(observed$statistic), reference$statistic,
               tolerance = 2e-12)
  expect_equal(observed$p.value,
               pnorm(reference$statistic, lower.tail = FALSE),
               tolerance = 2e-14)
  expect_equal(observed$diagnostics$primary.variance.factor, 1)
  expect_true(observed$diagnostics$book.sqrt.two.conflict)
})


test_that("CSM radial correction and truncated Cauchy match primary formulas", {
  x <- afc_fit_pair()
  maximum <- zhao_wang_conditional_spatial_sign_test(
    x$fit, component = "max", tol = 1e-6, max_iter = 2000L
  )
  r <- maximum$components$radii
  omega <- sum(x$fit$h^2) / x$fit$T
  inverse <- mean(1 / r)
  denominator <- 1 - 2 * (1 - omega) * inverse * mean(r) +
    (1 - omega) * mean(r^2) * inverse^2
  zeta <- x$fit$N * inverse^2 / denominator
  raw <- x$fit$T * zeta *
    max(maximum$components$standardized.theta^2)
  centered <- raw - 2 * log(x$fit$N) + log(log(x$fit$N))

  expect_equal(maximum$components$omega.ratio, omega,
               tolerance = 2e-13)
  expect_equal(maximum$components$zeta.denominator, denominator,
               tolerance = 2e-12)
  expect_equal(maximum$components$zeta, zeta, tolerance = 2e-12)
  expect_equal(maximum$components$CSM.raw, raw, tolerance = 2e-12)
  expect_equal(unname(maximum$statistic), centered, tolerance = 2e-12)

  combined <- zhao_wang_conditional_spatial_sign_test(
    x$fit, x$trace.fit, component = "combined",
    tol = 1e-6, max_iter = 2000L
  )
  p <- c(
    combined$components$CSS$p.value,
    combined$components$CSM.p.value
  )
  terms <- ifelse(p < 0.5, 0.5 / tan(pi * p), 0)
  score <- sum(terms)
  expected.p <- atan2(1, score) / pi
  expect_equal(unname(combined$statistic), score, tolerance = 2e-12)
  expect_equal(combined$p.value, expected.p, tolerance = 2e-14)
  expect_equal(combined$components$truncated.Cauchy$active, p < 0.5,
               ignore_attr = TRUE)
})


test_that("SS-BH uses restricted no-intercept residuals and full varsigma", {
  x <- afc_fixture()
  observed <- wang_zhao_feng_wang_mutual_fund_fdr(
    x$returns, x$factors, q = 0.2, tol = 1e-6, max_iter = 2000L
  )
  reference <- afc_project_r(x$returns, x$factors)
  scaled.reference <- reference$residuals / max(abs(x$returns))
  r <- observed$radii
  omega <- sum(reference$h^2) / nrow(x$returns)
  inverse <- mean(1 / r)
  denominator <- 1 - 2 * (1 - omega) * inverse * mean(r) +
    (1 - omega) * mean(r^2) * inverse^2
  varsigma <- ncol(x$returns) * inverse^2 / denominator
  statistic <- sqrt(nrow(x$returns) * varsigma) *
    observed$standardized.theta
  p <- pnorm(statistic, lower.tail = FALSE)
  ordering <- order(p, seq_along(p))
  passes <- p[ordering] <= 0.2 * seq_along(p) / length(p)
  k <- if (any(passes)) max(which(passes)) else 0L

  expect_s3_class(observed, "mutual_fund_fdr")
  expect_equal(observed$restricted.residuals.scaled, scaled.reference,
               tolerance = 2e-12, ignore_attr = TRUE)
  expect_equal(observed$omega.ratio, omega, tolerance = 2e-12)
  expect_equal(observed$varsigma.denominator, denominator,
               tolerance = 2e-12)
  expect_equal(observed$varsigma, varsigma, tolerance = 2e-12)
  expect_equal(observed$statistic, statistic, tolerance = 2e-12,
               ignore_attr = TRUE)
  expect_equal(observed$p.value, p, tolerance = 2e-14,
               ignore_attr = TRUE)
  expect_equal(observed$bh.k, k)
  expect_identical(
    observed$diagnostics$restricted.residual.source,
    "returns projected on observed factors without an intercept"
  )
  expect_false(observed$diagnostics$observable.projection.contains.intercept)
  expect_true(observed$diagnostics$book.varsigma.omission)
})


test_that("spatial-Kendall FSS factor extraction is literal and explicit", {
  x <- afc_fixture()
  observed <- wang_zhao_feng_wang_mutual_fund_fdr(
    x$returns, x$factors, q = 0.1,
    adjustment = "spatial_kendall", n_factors = 1L,
    tol = 1e-6, max_iter = 2000L
  )
  restricted <- afc_project_r(x$returns, x$factors)$residuals
  restricted <- restricted / max(abs(x$returns))
  expected.kendall <- afc_kendall_r(restricted)
  eig <- eigen(expected.kendall, symmetric = TRUE)
  loadings <- sqrt(ncol(restricted)) * eig$vectors[, 1L, drop = FALSE]
  scores <- restricted %*% loadings / ncol(restricted)
  adjusted <- restricted - scores %*% t(loadings)
  adjusted <- adjusted / max(abs(adjusted))

  expect_equal(observed$factor.adjustment$spatial.kendall,
               expected.kendall, tolerance = 2e-12)
  expect_equal(abs(observed$factor.adjustment$loadings), abs(loadings),
               tolerance = 2e-10)
  expect_equal(observed$factor.adjusted.residuals.scaled, adjusted,
               tolerance = 2e-10, ignore_attr = TRUE)
  expect_equal(observed$factor.adjustment$n.factors, 1L)
  expect_identical(observed$factor.adjustment$selection, "fixed")

  expect_error(
    wang_zhao_feng_wang_mutual_fund_fdr(
      x$returns, x$factors, adjustment = "spatial_kendall"
    ),
    "supply n_factors or an explicit k_max|Supply n_factors or an explicit k_max"
  )
})


test_that("factor nuisance terms and row order have the required invariances", {
  x <- afc_fixture()
  base <- wang_zhao_feng_wang_mutual_fund_fdr(
    x$returns, x$factors, tol = 1e-6, max_iter = 2000L
  )
  factor.shift <- x$factors %*% matrix(c(0.9, -0.7, 0.4, 1.1), 1L)
  shifted <- wang_zhao_feng_wang_mutual_fund_fdr(
    x$returns + factor.shift, x$factors,
    tol = 1e-6, max_iter = 2000L
  )
  expect_equal(shifted$statistic, base$statistic, tolerance = 2e-9)
  expect_equal(shifted$p.value, base$p.value, tolerance = 2e-10)

  order <- c(seq(2, 20, 2), seq(1, 19, 2))
  permuted <- wang_zhao_feng_wang_mutual_fund_fdr(
    x$returns[order, ], x$factors[order, , drop = FALSE],
    tol = 1e-6, max_iter = 2000L
  )
  expect_equal(permuted$statistic, base$statistic, tolerance = 2e-9)
  expect_equal(permuted$p.value, base$p.value, tolerance = 2e-10)
})


test_that("conditional tests are invariant to common scale and joint row order", {
  base <- afc_fit_pair(1)
  huge <- afc_fit_pair(1e150)
  sum.base <- ma_lan_su_tsai_conditional_alpha_sum_test(base$fit)
  sum.huge <- ma_lan_su_tsai_conditional_alpha_sum_test(huge$fit)
  max.base <- ma_feng_wang_bao_conditional_alpha_test(
    base$fit, component = "max"
  )
  max.huge <- ma_feng_wang_bao_conditional_alpha_test(
    huge$fit, component = "max"
  )
  css.base <- zhao_conditional_spatial_sign_sum_test(
    base$fit, base$trace.fit
  )
  css.huge <- zhao_conditional_spatial_sign_sum_test(
    huge$fit, huge$trace.fit
  )
  csm.base <- zhao_wang_conditional_spatial_sign_test(
    base$fit, tol = 1e-6, max_iter = 2000L
  )
  csm.huge <- zhao_wang_conditional_spatial_sign_test(
    huge$fit, tol = 1e-6, max_iter = 2000L
  )

  expect_equal(unname(sum.huge$statistic), unname(sum.base$statistic),
               tolerance = 2e-10)
  expect_equal(unname(max.huge$statistic), unname(max.base$statistic),
               tolerance = 2e-10)
  expect_equal(unname(css.huge$statistic), unname(css.base$statistic),
               tolerance = 2e-10)
  expect_equal(unname(csm.huge$statistic), unname(csm.base$statistic),
               tolerance = 2e-8)
  expect_true(all(is.finite(c(
    sum.huge$p.value, max.huge$p.value, css.huge$p.value, csm.huge$p.value
  ))))

  order <- rev(seq_len(base$fit$T))
  null.design <- unclass(base$null.design)[order, , drop = FALSE]
  trace.design <- unclass(base$trace.design)[order, , drop = FALSE]
  attr(null.design, "factor.count") <- 1L
  attr(null.design, "center.alpha") <- TRUE
  attr(trace.design, "factor.count") <- 1L
  attr(trace.design, "center.alpha") <- FALSE
  fit.p <- conditional_alpha_sieve_fit(
    base$data$returns[order, ], null.design
  )
  trace.p <- conditional_alpha_sieve_fit(
    base$data$returns[order, ], trace.design
  )
  sum.p <- ma_lan_su_tsai_conditional_alpha_sum_test(fit.p)
  css.p <- zhao_conditional_spatial_sign_sum_test(fit.p, trace.p)
  expect_equal(unname(sum.p$statistic), unname(sum.base$statistic),
               tolerance = 2e-10)
  expect_equal(unname(css.p$statistic), unname(css.base$statistic),
               tolerance = 2e-10)
})


test_that("CSM is equivariant to positive assetwise rescaling", {
  x <- afc_fixture()
  design <- conditional_alpha_sieve_design(x$basis, x$factors)
  fit1 <- conditional_alpha_sieve_fit(x$returns, design)
  fit2 <- conditional_alpha_sieve_fit(
    sweep(x$returns, 2L, c(1e-8, 3, 2e5, 0.4), "*"), design
  )
  one <- zhao_wang_conditional_spatial_sign_test(
    fit1, tol = 1e-6, max_iter = 2000L
  )
  two <- zhao_wang_conditional_spatial_sign_test(
    fit2, tol = 1e-6, max_iter = 2000L
  )
  expect_equal(unname(two$statistic), unname(one$statistic),
               tolerance = 2e-6)
  expect_equal(two$p.value, one$p.value, tolerance = 2e-7)
})


test_that("supplied FDR residual and latent-fit contracts are explicit", {
  x <- afc_fixture()
  restricted <- afc_project_r(x$returns, x$factors)
  direct <- wang_zhao_feng_wang_mutual_fund_fdr(
    restricted_residuals = restricted$residuals,
    omega_ratio = sum(restricted$h^2) / nrow(x$returns),
    tol = 1e-6, max_iter = 2000L
  )
  supplied <- wang_zhao_feng_wang_mutual_fund_fdr(
    restricted_residuals = restricted$residuals,
    omega_ratio = sum(restricted$h^2) / nrow(x$returns),
    adjustment = "supplied",
    latent_fit = list(
      adjusted.residuals = restricted$residuals,
      n.factors = 0L
    ), tol = 1e-6, max_iter = 2000L
  )
  expect_equal(supplied$statistic, direct$statistic, tolerance = 2e-9)
  expect_equal(supplied$p.value, direct$p.value, tolerance = 2e-10)
  expect_error(
    wang_zhao_feng_wang_mutual_fund_fdr(
      restricted_residuals = restricted$residuals
    ),
    "omega_ratio"
  )
  expect_error(
    wang_zhao_feng_wang_mutual_fund_fdr(
      restricted_residuals = restricted$residuals,
      omega_ratio = 0.9, adjustment = "supplied", latent_fit = list()
    ),
    "must contain"
  )
})


test_that("invalid, rank-deficient, and degenerate inputs fail without repair", {
  x <- afc_fixture()
  expect_error(conditional_alpha_sieve_design(x$basis, cbind(1, 1)),
               "same number of rows|at least")
  expect_error(
    conditional_alpha_sieve_design(cbind(x$basis, x$basis[, 1])),
    "full column rank"
  )
  expect_error(
    conditional_alpha_sieve_fit(
      x$returns, cbind(x$basis[, 1], x$basis[, 1])
    ),
    "full column rank"
  )
  bad <- x$returns
  bad[1, 1] <- NA_real_
  expect_error(
    conditional_alpha_sieve_fit(
      bad, conditional_alpha_sieve_design(x$basis, x$factors)
    ),
    "finite"
  )

  pair <- afc_fit_pair()
  expect_error(
    ma_feng_wang_bao_conditional_alpha_test(
      pair$fit, factor_count = 2L, component = "max"
    ),
    "conflicts"
  )
  expect_error(
    zhao_conditional_spatial_sign_sum_test(pair$fit, pair$fit),
    "uncentered"
  )
  expect_error(
    zhao_wang_conditional_spatial_sign_test(
      pair$fit, component = "combined"
    ),
    "trace_fit"
  )
  expect_error(
    zhao_wang_conditional_spatial_sign_test(pair$fit, tol = 0),
    "tol"
  )
  expect_error(
    wang_zhao_feng_wang_mutual_fund_fdr(
      x$returns, x$factors, q = 1
    ),
    "q"
  )

  zero.row <- rbind(0, diag(4), -diag(4))
  expect_error(
    wang_zhao_feng_wang_mutual_fund_fdr(
      restricted_residuals = zero.row, omega_ratio = 1,
      max_iter = 2000L
    ),
    "radius|inverse-radius|standardised residual"
  )
})

test_that("spatial-Kendall eigenvalue-ratio selection is literal", {
  x <- afc_fixture()
  observed <- wang_zhao_feng_wang_mutual_fund_fdr(
    x$returns, x$factors, adjustment = "spatial_kendall", k_max = 2L,
    tol = 1e-6, max_iter = 2000L
  )
  values <- observed$factor.adjustment$eigenvalues
  ratios <- values[1:2] / values[2:3]

  expect_equal(observed$factor.adjustment$eigenvalue.ratios, ratios,
               tolerance = 2e-12)
  expect_equal(observed$factor.adjustment$n.factors, which.max(ratios))
  expect_identical(
    observed$factor.adjustment$selection,
    "paper eigenvalue-ratio rule with explicit k_max"
  )
  expect_identical(observed$factor.adjustment$k.max, 2L)
})


test_that("FDR nuisance arguments cannot silently contradict each other", {
  x <- afc_fixture()
  restricted <- afc_project_r(x$returns, x$factors)
  residuals <- restricted$residuals
  dimnames(residuals) <- dimnames(x$returns)
  omega <- sum(restricted$h^2) / nrow(residuals)

  consistent <- wang_zhao_feng_wang_mutual_fund_fdr(
    restricted_residuals = residuals, factors = x$factors,
    omega_ratio = omega, tol = 1e-6, max_iter = 2000L
  )
  expect_equal(consistent$omega.ratio, omega, tolerance = 2e-12)
  expect_error(
    wang_zhao_feng_wang_mutual_fund_fdr(
      restricted_residuals = residuals, factors = x$factors,
      omega_ratio = 0.97 * omega
    ),
    "conflicts with the ratio implied by factors"
  )
  expect_error(
    wang_zhao_feng_wang_mutual_fund_fdr(
      x$returns, x$factors, n_factors = 1L
    ),
    "not used when"
  )
  expect_error(
    wang_zhao_feng_wang_mutual_fund_fdr(
      restricted_residuals = residuals, omega_ratio = omega,
      adjustment = "supplied",
      latent_fit = list(adjusted.residuals = residuals), k_max = 1L
    ),
    "used only with"
  )
})


test_that("paired and supplied nuisance components enforce row and asset order", {
  x <- afc_fixture()
  named.factors <- x$factors
  rownames(named.factors) <- rownames(x$returns)
  expect_error(
    wang_zhao_feng_wang_mutual_fund_fdr(
      x$returns, named.factors[rev(seq_len(nrow(named.factors))), ,
                               drop = FALSE]
    ),
    "row names"
  )

  restricted <- afc_project_r(x$returns, x$factors)
  residuals <- restricted$residuals
  dimnames(residuals) <- dimnames(x$returns)
  omega <- sum(restricted$h^2) / nrow(residuals)
  expect_error(
    wang_zhao_feng_wang_mutual_fund_fdr(
      restricted_residuals = residuals, omega_ratio = omega,
      adjustment = "supplied",
      latent_fit = list(
        adjusted.residuals = residuals[rev(seq_len(nrow(residuals))), ,
                                       drop = FALSE]
      )
    ),
    "row names"
  )

  pair <- afc_fit_pair()
  row.bad <- pair$trace.fit
  order <- rev(seq_len(row.bad$T))
  row.bad$residuals.scaled <- row.bad$residuals.scaled[order, , drop = FALSE]
  row.bad$h <- row.bad$h[order]
  row.bad$observation.names <- row.bad$observation.names[order]
  expect_error(
    zhao_conditional_spatial_sign_sum_test(pair$fit, row.bad),
    "different observation orders"
  )

  asset.bad <- pair$trace.fit
  order <- c(2L, 1L, 3L, 4L)
  asset.bad$residuals.scaled <-
    asset.bad$residuals.scaled[, order, drop = FALSE]
  asset.bad$asset.names <- asset.bad$asset.names[order]
  expect_error(
    zhao_conditional_spatial_sign_sum_test(pair$fit, asset.bad),
    "different asset orders"
  )
})


test_that("tampered conditional sieve fits are rejected", {
  fit <- afc_fit_pair()$fit

  bad <- fit
  bad$h2 <- 1.01 * bad$h2
  expect_error(
    ma_lan_su_tsai_conditional_alpha_sum_test(bad),
    "valid conditional_alpha_sieve_fit"
  )

  bad <- fit
  bad$asset.names <- rev(bad$asset.names)
  expect_error(
    ma_feng_wang_bao_conditional_alpha_test(bad, component = "max"),
    "valid conditional_alpha_sieve_fit"
  )

  bad <- fit
  bad$T <- bad$T + 1L
  expect_error(
    zhao_wang_conditional_spatial_sign_test(bad),
    "valid conditional_alpha_sieve_fit"
  )
})


test_that("all public conditional-alpha methods leave the RNG state unchanged", {
  x <- afc_fixture()
  set.seed(40404)
  before <- .Random.seed

  null.design <- conditional_alpha_sieve_design(x$basis, x$factors)
  trace.design <- conditional_alpha_sieve_design(
    x$basis, x$factors, center_alpha = FALSE
  )
  fit <- conditional_alpha_sieve_fit(x$returns, null.design)
  trace.fit <- conditional_alpha_sieve_fit(x$returns, trace.design)
  invisible(wang_zhao_feng_wang_mutual_fund_fdr(
    x$returns, x$factors, tol = 1e-6, max_iter = 2000L
  ))
  invisible(ma_lan_su_tsai_conditional_alpha_sum_test(fit))
  invisible(ma_feng_wang_bao_conditional_alpha_test(
    fit, component = "adaptive"
  ))
  invisible(zhao_conditional_spatial_sign_sum_test(fit, trace.fit))
  invisible(zhao_wang_conditional_spatial_sign_test(
    fit, trace.fit, component = "combined",
    tol = 1e-6, max_iter = 2000L
  ))

  expect_identical(.Random.seed, before)
})


test_that("SS-BH is equivariant to positive assetwise rescaling", {
  x <- afc_fixture()
  one <- wang_zhao_feng_wang_mutual_fund_fdr(
    x$returns, x$factors, q = 0.2, tol = 1e-6, max_iter = 2000L
  )
  scaled.returns <- sweep(
    x$returns, 2L, c(1e-7, 0.3, 8, 2e4), "*"
  )
  two <- wang_zhao_feng_wang_mutual_fund_fdr(
    scaled.returns, x$factors, q = 0.2,
    tol = 1e-6, max_iter = 2000L
  )

  expect_equal(two$statistic, one$statistic, tolerance = 3e-6)
  expect_equal(two$p.value, one$p.value, tolerance = 3e-7)
  expect_identical(two$rejected.names, one$rejected.names)
})


test_that("truncated Cauchy uses the strict one-half activation boundary", {
  boundary <- HDElliptical:::.ch4_afc_truncated_cauchy(
    c(0, 0.25, 0.5, 1)
  )
  expect_identical(boundary$active, c(TRUE, TRUE, FALSE, FALSE))
  expect_equal(boundary$terms[3:4], c(0, 0))
  expect_true(is.infinite(boundary$statistic))
  expect_equal(boundary$p.value, 0)

  finite <- HDElliptical:::.ch4_afc_truncated_cauchy(c(0.25, 0.5))
  expect_equal(finite$statistic, 0.5, tolerance = 1e-15)
  expect_equal(finite$p.value, atan2(1, 0.5) / pi, tolerance = 1e-15)
  expect_error(
    HDElliptical:::.ch4_afc_truncated_cauchy(c(0.2, 1.01)),
    "must lie in"
  )
})
