alpha_fixture <- function() {
  time <- seq_len(12)
  factors <- matrix(seq(-0.6, 1.4, length.out = 12), ncol = 1L)
  colnames(factors) <- "market"
  error <- cbind(
    sin(0.7 * time) + 0.2 * cos(1.3 * time),
    1.2 * cos(0.4 * time) + 0.3 * sin(1.1 * time),
    0.8 * sin(0.3 * time + 1) + 0.6 * cos(0.9 * time)
  )
  returns <- sweep(error, 2L, c(0.15, -0.08, 0.05), "+") +
    factors %*% matrix(c(0.6, -0.4, 0.2), nrow = 1L)
  colnames(returns) <- c("A", "B", "C")
  list(returns = returns, factors = factors)
}

alpha_literal_ols <- function(returns, factors) {
  design <- cbind(intercept = 1, factors)
  coefficients <- solve(crossprod(design), crossprod(design, returns))
  residuals <- returns - design %*% coefficients
  h <- rep(1, nrow(returns)) -
    factors %*% solve(crossprod(factors),
                      crossprod(factors, rep(1, nrow(returns))))
  h <- as.numeric(h)
  v <- nrow(returns) - ncol(factors) - 1
  t.squared <- coefficients[1, ]^2 * sum(h^2) /
    (colSums(residuals^2) / v)
  correlation <- stats::cov2cor(crossprod(residuals))
  list(
    coefficients = coefficients,
    residuals = residuals,
    h = h,
    h2 = sum(h^2),
    v = v,
    t.squared = t.squared,
    correlation = correlation
  )
}

alpha_literal_direction <- function(residual, diagonal) {
  standardized <- residual / sqrt(diagonal)
  standardized / sqrt(sum(standardized^2))
}

alpha_literal_lfm_trace <- function(returns, factors, h, diagonal) {
  total <- 0
  n <- nrow(returns)
  p <- ncol(returns)
  for (first in seq_len(n)) {
    for (second in seq_len(n)) {
      if (first == second) next
      remaining <- setdiff(seq_len(n), c(first, second))
      split <- floor(length(remaining) / 2)
      first.half <- remaining[seq_len(split)]
      second.half <- remaining[seq.int(split + 1L, length(remaining))]
      slope1 <- solve(
        crossprod(factors[first.half, , drop = FALSE]),
        crossprod(
          factors[first.half, , drop = FALSE],
          returns[first.half, , drop = FALSE]
        )
      )
      slope2 <- solve(
        crossprod(factors[second.half, , drop = FALSE]),
        crossprod(
          factors[second.half, , drop = FALSE],
          returns[second.half, , drop = FALSE]
        )
      )
      residual1 <- returns[first, ] -
        factors[first, , drop = FALSE] %*% slope1
      residual2 <- returns[second, ] -
        factors[second, , drop = FALSE] %*% slope2
      u1 <- alpha_literal_direction(as.numeric(residual1), diagonal)
      u2 <- alpha_literal_direction(as.numeric(residual2), diagonal)
      total <- total + h[first]^2 * h[second]^2 * sum(u1 * u2)^2
    }
  }
  p^2 * total / (sum(h^2) * (sum(h^2) - 1))
}


test_that("GRS matches joint-OLS primary formula and locks the divisor", {
  x <- alpha_fixture()
  literal <- alpha_literal_ols(x$returns, x$factors)
  observed <- grs_alpha_test(x$returns, x$factors)
  T <- nrow(x$returns)
  N <- ncol(x$returns)
  K <- ncol(x$factors)
  scatter <- crossprod(literal$residuals) / T
  quadratic <- drop(
    literal$coefficients[1, ] %*% solve(scatter) %*%
      literal$coefficients[1, ]
  )
  expected <- (T - N - K) / N * literal$h2 / T * quadratic

  expect_s3_class(observed, "htest")
  expect_equal(unname(observed$statistic), expected, tolerance = 1e-12)
  expect_equal(
    observed$p.value,
    unname(stats::pf(expected, N, T - N - K, lower.tail = FALSE)),
    tolerance = 1e-12
  )
  expect_equal(unname(observed$estimate),
               unname(literal$coefficients[1, ]), tolerance = 1e-12)
  expect_equal(observed$components$h, literal$h, tolerance = 1e-12)
  expect_identical(observed$diagnostics$covariance.divisor, "T")
  expect_true(observed$diagnostics$book.slope.conflict)
  expect_true(observed$diagnostics$book.covariance.divisor.conflict)

  book.mixed <- (T - N - K) / N * literal$h2 / T * drop(
    literal$coefficients[1, ] %*%
      solve(crossprod(literal$residuals) / literal$v) %*%
      literal$coefficients[1, ]
  )
  expect_equal(book.mixed, expected * literal$v / T, tolerance = 1e-12)
  expect_gt(abs(book.mixed - expected), 1e-6)
})


test_that("PY and FLLM max match literal finite-v and Gumbel formulas", {
  x <- alpha_fixture()
  literal <- alpha_literal_ols(x$returns, x$factors)
  N <- ncol(x$returns)
  v <- literal$v
  p0 <- 0.2
  delta <- 1.3
  threshold <- qnorm(1 - p0 / (2 * N^delta))
  rho <- literal$correlation[upper.tri(literal$correlation)]
  keep <- abs(sqrt(v) * rho) > threshold
  rho2 <- 2 * sum(rho[keep]^2) / (N * (N - 1))
  center <- v / (v - 2)
  scale <- center * sqrt(
    2 * (v - 1) / (v - 4) * (1 + (N - 1) * rho2)
  )
  expected.py <- sum(literal$t.squared - center) /
    sqrt(N) / scale

  py <- pesaran_yamagata_alpha_test(
    x$returns, x$factors, p0 = p0, delta = delta
  )
  expect_equal(unname(py$statistic), expected.py, tolerance = 1e-12)
  expect_equal(py$components$rho.tilde.squared, rho2,
               tolerance = 1e-14)
  expect_identical(py$components$retained.upper.triangle, keep)
  expect_equal(py$p.value, pnorm(expected.py, lower.tail = FALSE),
               tolerance = 1e-14)
  expect_false(py$diagnostics$book.generic.sum.is.primary.PY)

  maximum <- max(literal$t.squared)
  centered <- maximum - 2 * log(N) + log(log(N))
  intensity <- exp(-centered / 2) / sqrt(pi)
  max.test <- feng_lan_liu_ma_alpha_max_test(
    x$returns, x$factors
  )
  expect_equal(unname(max.test$raw.statistic), maximum,
               tolerance = 1e-12)
  expect_equal(unname(max.test$statistic), centered,
               tolerance = 1e-12)
  expect_equal(max.test$p.value, -expm1(-intensity),
               tolerance = 1e-14)
})


test_that("Gaussian COM is the primary Bonferroni minimum-p rule", {
  x <- alpha_fixture()
  observed <- gaussian_alpha_combination_test(
    x$returns, x$factors, p0 = 0.15, delta = 0.8
  )
  component <- observed$components$p.values
  expect_equal(observed$p.value, min(1, 2 * min(component)),
               tolerance = 0)
  expect_equal(unname(observed$statistic), min(component),
               tolerance = 0)
  expect_identical(observed$components$multiplicity.factor, 2)
  expect_false(observed$diagnostics$book.generic.Cauchy.implemented)
})


test_that("LFM Q, estimating equation, and ordered split trace are literal", {
  x <- alpha_fixture()
  y.scale <- apply(abs(x$returns), 2L, max)
  f.scale <- apply(abs(x$factors), 2L, max)
  y <- sweep(x$returns, 2L, y.scale, "/")
  f <- sweep(x$factors, 2L, f.scale, "/")
  observed <- liu_feng_ma_spatial_sign_alpha_test(
    x$returns, x$factors, bias = "none", tol = 1e-9,
    max_iter = 3000L
  )
  U <- observed$components$directions
  h <- observed$components$h
  h2 <- sum(h^2)
  Q <- ncol(y) / h2 * (
    sum(colSums(U * h)^2) -
      sum(h^2 * rowSums(U^2))
  )
  expected.trace <- alpha_literal_lfm_trace(
    y, f, h,
    observed$components$scale.diagonal.scaled.coordinates
  )

  expect_equal(observed$components$Q, Q, tolerance = 1e-11)
  expect_equal(
    unname(ncol(y) * colMeans(U^2)), rep(1, ncol(y)),
    tolerance = 2e-9
  )
  expect_equal(rowSums(U^2), rep(1, nrow(U)), tolerance = 1e-12)
  expect_equal(observed$components$trace.R.squared,
               expected.trace, tolerance = 2e-10)
  expect_equal(observed$components$trace.denominator,
               h2 * (h2 - 1), tolerance = 1e-12)
  expect_gt(
    abs(observed$components$trace.denominator - h2^2),
    1
  )
  expect_equal(
    unname(observed$statistic),
    Q / sqrt(2 * expected.trace),
    tolerance = 2e-10
  )
  expect_true(observed$diagnostics$book.trace.denominator.conflict)
  expect_true(observed$diagnostics$book.bias.omission)
  expect_true(observed$diagnostics$book.leaveout.scale.conflict)
  expect_match(observed$diagnostics$weighted.INST.status, "review-only")
  expect_match(observed$diagnostics$dependent.Lq.status, "review-only")
})


test_that("LFM bootstrap uses iid asset-time signs and isolates explicit seed", {
  x <- alpha_fixture()
  set.seed(812)
  state <- .Random.seed
  observed <- liu_feng_ma_spatial_sign_alpha_test(
    x$returns, x$factors, bias = "wild_bootstrap",
    bootstrap_reps = 3L, seed = 47, keep_bootstrap = TRUE,
    tol = 1e-7, max_iter = 1000L
  )
  expect_identical(.Random.seed, state)
  expect_length(observed$components$bootstrap.Q, 3L)
  expect_equal(
    observed$components$delta.Q,
    mean(observed$components$bootstrap.Q),
    tolerance = 0
  )
  repeated <- liu_feng_ma_spatial_sign_alpha_test(
    x$returns, x$factors, bias = "wild_bootstrap",
    bootstrap_reps = 3L, seed = 47, keep_bootstrap = TRUE,
    tol = 1e-7, max_iter = 1000L
  )
  expect_identical(.Random.seed, state)
  expect_equal(repeated$components$bootstrap.Q,
               observed$components$bootstrap.Q, tolerance = 0)
  expect_true(observed$diagnostics$rng.isolated)

  y.scale <- apply(abs(x$returns), 2L, max)
  f.scale <- apply(abs(x$factors), 2L, max)
  y <- sweep(x$returns, 2L, y.scale, "/")
  f <- sweep(x$factors, 2L, f.scale, "/")
  unrestricted <- lm.fit(cbind(1, f), y)
  restricted.slope <- solve(crossprod(f), crossprod(f, y))
  restricted.fitted <- f %*% restricted.slope
  set.seed(47)
  literal.q <- numeric(3)
  for (b in seq_len(3)) {
    multipliers <- matrix(
      sample(c(-1, 1), nrow(y) * ncol(y), replace = TRUE),
      nrow = nrow(y), ncol = ncol(y)
    )
    core <- HDElliptical:::cpp_ch4_lfm_spatial_sign_core(
      restricted.fitted + unrestricted$residuals * multipliers,
      f, 1e-7, 1000L, 0, FALSE
    )
    literal.q[b] <- core$Q
  }
  expect_equal(observed$components$bootstrap.Q, literal.q,
               tolerance = 1e-12)

  set.seed(99)
  before <- .Random.seed
  invisible(liu_feng_ma_spatial_sign_alpha_test(
    x$returns, x$factors, bias = "wild_bootstrap",
    bootstrap_reps = 1L, seed = NULL, tol = 1e-7,
    max_iter = 1000L
  ))
  expect_false(identical(.Random.seed, before))
})


test_that("NULL and vector factor interfaces have the stated contracts", {
  x <- alpha_fixture()
  vector.factor <- drop(x$factors)
  matrix.calls <- list(
    grs_alpha_test(x$returns, x$factors),
    pesaran_yamagata_alpha_test(x$returns, x$factors),
    feng_lan_liu_ma_alpha_max_test(x$returns, x$factors),
    gaussian_alpha_combination_test(x$returns, x$factors),
    liu_feng_ma_spatial_sign_alpha_test(
      x$returns, x$factors, bias = "none", max_iter = 2000L
    ),
    zhao_feng_wang_wang_robust_alpha_test(
      x$returns, x$factors, max_iter = 2000L
    )
  )
  vector.calls <- list(
    grs_alpha_test(x$returns, vector.factor),
    pesaran_yamagata_alpha_test(x$returns, vector.factor),
    feng_lan_liu_ma_alpha_max_test(x$returns, vector.factor),
    gaussian_alpha_combination_test(x$returns, vector.factor),
    liu_feng_ma_spatial_sign_alpha_test(
      x$returns, vector.factor, bias = "none", max_iter = 2000L
    ),
    zhao_feng_wang_wang_robust_alpha_test(
      x$returns, vector.factor, max_iter = 2000L
    )
  )
  for (i in seq_along(matrix.calls)) {
    expect_equal(
      unname(vector.calls[[i]]$statistic),
      unname(matrix.calls[[i]]$statistic),
      tolerance = 1e-10
    )
  }

  no.factor <- zhao_feng_wang_wang_robust_alpha_test(
    x$returns, factors = NULL, max_iter = 2000L
  )
  expect_equal(no.factor$components$eta, 0, tolerance = 0)
  expect_equal(no.factor$components$omega, 1, tolerance = 0)
  expect_equal(
    no.factor$components$zeta,
    unname(ncol(x$returns) *
      no.factor$components$radial.moments["mean.inverse.r"]^2),
    tolerance = 1e-12
  )
})


test_that("ZFWW zeta, robust max, and truncated Cauchy are primary formulas", {
  x <- alpha_fixture()
  max.test <- zhao_feng_wang_wang_robust_alpha_test(
    x$returns, x$factors, component = "max",
    tol = 1e-7, max_iter = 2000L
  )
  radial <- max.test$components$radial.moments
  eta <- max.test$components$eta
  denominator <- 1 -
    2 * eta * radial["mean.inverse.r"] * radial["mean.r"] +
    eta * radial["mean.inverse.r.squared"] * radial["mean.r.squared"]
  zeta <- ncol(x$returns) * radial["mean.inverse.r"]^2 / denominator
  raw <- nrow(x$returns) * zeta *
    max(max.test$components$standardized.theta^2)
  centered <- raw - 2 * log(ncol(x$returns)) +
    log(log(ncol(x$returns)))

  expect_equal(max.test$components$zeta.denominator,
               unname(denominator), tolerance = 1e-13)
  expect_equal(max.test$components$zeta, unname(zeta),
               tolerance = 1e-13)
  expect_equal(max.test$components$robust.maximum, unname(raw),
               tolerance = 1e-12)
  expect_equal(unname(max.test$statistic), unname(centered),
               tolerance = 1e-12)
  expect_true(max.test$diagnostics$book.zeta.inverse.moment.conflict)
  expect_true(max.test$diagnostics$book.zeta.eta.squared.conflict)

  book.zeta <- ncol(x$returns) / radial["mean.inverse.r"]^2 /
    (1 - 2 * eta * radial["mean.inverse.r"] * radial["mean.r"] +
       eta^2 * radial["mean.inverse.r.squared"] *
         radial["mean.r.squared"])
  expect_gt(abs(unname(book.zeta - zeta)), 1e-4)

  combined <- zhao_feng_wang_wang_robust_alpha_test(
    x$returns, x$factors, component = "combined", bias = "none",
    tol = 1e-7, max_iter = 2000L
  )
  p <- c(combined$components$p.SS, combined$components$p.SM)
  terms <- ifelse(p < 0.5, 0.5 / tan(pi * p), 0)
  score <- sum(terms)
  expect_equal(unname(combined$statistic), score, tolerance = 1e-14)
  expect_equal(combined$p.value, atan2(1, score) / pi,
               tolerance = 1e-14)
  expect_identical(
    combined$diagnostics$combined.calibration,
    "primary truncated-Cauchy combination"
  )
  expect_false(combined$diagnostics$book.generic.Cauchy.implemented)
})


test_that("all six procedures obey their scale and factor-span invariances", {
  x <- alpha_fixture()
  asset.scale <- c(-2.5, 0.4, 3.2)
  factor.scale <- -3.7
  shifted <- x$returns +
    x$factors %*% matrix(c(1.1, -0.7, 0.3), nrow = 1L)

  calls <- list(
    GRS = function(y, f) grs_alpha_test(y, f),
    PY = function(y, f) pesaran_yamagata_alpha_test(y, f),
    MAX = function(y, f) feng_lan_liu_ma_alpha_max_test(y, f),
    COM = function(y, f) gaussian_alpha_combination_test(y, f),
    SS = function(y, f) liu_feng_ma_spatial_sign_alpha_test(
      y, f, bias = "none", tol = 1e-7, max_iter = 2000L
    ),
    SM = function(y, f) zhao_feng_wang_wang_robust_alpha_test(
      y, f, component = "max", tol = 1e-7, max_iter = 2000L
    )
  )
  for (name in names(calls)) {
    baseline <- calls[[name]](x$returns, x$factors)
    rescaled <- calls[[name]](
      sweep(x$returns, 2L, asset.scale, "*"),
      x$factors * factor.scale
    )
    factor.shift <- calls[[name]](shifted, x$factors)
    expect_equal(
      unname(rescaled$statistic), unname(baseline$statistic),
      tolerance = 3e-7, info = name
    )
    expect_equal(
      unname(factor.shift$statistic), unname(baseline$statistic),
      tolerance = 3e-7, info = name
    )
  }
})


test_that("column pre-scaling handles extreme global magnitudes", {
  x <- alpha_fixture()
  functions <- list(
    function(y, f) grs_alpha_test(y, f),
    function(y, f) pesaran_yamagata_alpha_test(y, f),
    function(y, f) feng_lan_liu_ma_alpha_max_test(y, f),
    function(y, f) gaussian_alpha_combination_test(y, f),
    function(y, f) liu_feng_ma_spatial_sign_alpha_test(
      y, f, bias = "none", tol = 1e-7, max_iter = 2000L
    ),
    function(y, f) zhao_feng_wang_wang_robust_alpha_test(
      y, f, tol = 1e-7, max_iter = 2000L
    )
  )
  for (fun in functions) {
    baseline <- fun(x$returns, x$factors)
    huge <- fun(x$returns * 1e150, x$factors * 1e150)
    tiny <- fun(x$returns * 1e-150, x$factors * 1e-150)
    expect_true(is.finite(unname(huge$statistic)))
    expect_true(is.finite(unname(tiny$statistic)))
    expect_equal(unname(huge$statistic), unname(baseline$statistic),
                 tolerance = 3e-7)
    expect_equal(unname(tiny$statistic), unname(baseline$statistic),
                 tolerance = 3e-7)
  }
})


test_that("alpha APIs reject undefined and silently repaired cases", {
  x <- alpha_fixture()
  bad <- x$returns
  bad[1, 1] <- NA_real_
  expect_error(grs_alpha_test(bad, x$factors), "finite")
  expect_error(grs_alpha_test(x$returns, x$factors[-1, ]),
               "at least")
  expect_error(grs_alpha_test(x$returns, matrix(0, 12, 1)),
               "positive finite magnitude")
  expect_error(grs_alpha_test(
    x$returns, cbind(x$factors, x$factors)
  ), "rank deficient")
  expect_error(zhao_feng_wang_wang_robust_alpha_test(
    x$returns, matrix(1, 12, 1)
  ), "rank deficient")
  expect_error(
    grs_alpha_test(cbind(x$returns, x$returns[, 1:9 %% 3 + 1]),
                   x$factors),
    "T > N"
  )
  expect_error(grs_alpha_test(cbind(x$returns[, 1], x$returns[, 1]),
                              x$factors),
               "positive definite")
  expect_error(pesaran_yamagata_alpha_test(x$returns[, 1, drop = FALSE],
                                           x$factors),
               "at least")
  expect_error(
    pesaran_yamagata_alpha_test(x$returns[1:6, ], x$factors[1:6, ]),
    "v > 4"
  )
  expect_error(pesaran_yamagata_alpha_test(x$returns, x$factors, p0 = 1),
               "strictly between")
  expect_error(pesaran_yamagata_alpha_test(x$returns, x$factors, delta = 0),
               "positive")
  expect_error(feng_lan_liu_ma_alpha_max_test(
    x$returns[, 1, drop = FALSE], x$factors
  ), "at least")

  expect_error(liu_feng_ma_spatial_sign_alpha_test(
    x$returns, x$factors, bias = "supplied"
  ), "delta_q")
  expect_error(liu_feng_ma_spatial_sign_alpha_test(
    x$returns, x$factors, bias = "none", delta_q = 0
  ), "must be NULL")
  expect_error(liu_feng_ma_spatial_sign_alpha_test(
    x$returns, x$factors, bootstrap_reps = 0
  ), "positive integer")
  expect_error(liu_feng_ma_spatial_sign_alpha_test(
    x$returns, x$factors, seed = -1
  ), "seed")
  expect_error(liu_feng_ma_spatial_sign_alpha_test(
    x$returns, x$factors, keep_bootstrap = NA
  ), "TRUE or FALSE")
  expect_error(liu_feng_ma_spatial_sign_alpha_test(
    x$returns, x$factors, tol = 0
  ), "positive")
  expect_error(liu_feng_ma_spatial_sign_alpha_test(
    x$returns, x$factors, max_iter = 1L, bias = "none",
    tol = 1e-15
  ), "did not satisfy")
  expect_error(liu_feng_ma_spatial_sign_alpha_test(
    matrix(0, 8, 2), NULL, bias = "none"
  ), "positive finite magnitude")

  expect_error(zhao_feng_wang_wang_robust_alpha_test(
    matrix(1, 8, 2), NULL
  ), "marginal sample variation")
  expect_error(zhao_feng_wang_wang_robust_alpha_test(
    x$returns, x$factors, component = "bad"
  ), "arg")
  expect_false(exists("zhao_chen_zi_weighted_alpha_test",
                      mode = "function"))
  expect_false(exists("dependent_alpha_test", mode = "function"))
  expect_false(exists("lq_alpha_test", mode = "function"))
})
