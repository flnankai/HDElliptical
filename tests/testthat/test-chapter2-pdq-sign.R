pdq_test_fixture <- function() {
  v <- rbind(
    c(1, 2, 3),
    c(2, -3, 1),
    c(-4, 1, 2),
    c(3, 4, -2)
  )
  w <- v %*% matrix(c(
    1.2, 0.1, -0.2,
    0.3, 0.9, 0.25,
    -0.15, 0.2, 1.1
  ), 3, 3, byrow = TRUE)
  x <- sweep(rbind(v, -v), 2L, c(0.4, -0.3, 0.2), "+")
  y <- sweep(rbind(w, -w), 2L, c(-0.2, 0.5, -0.4), "+")
  colnames(x) <- colnames(y) <- c("u", "v", "w")
  rownames(x) <- paste0("A", seq_len(nrow(x)))
  rownames(y) <- paste0("B", seq_len(nrow(y)))
  list(x = x, y = y)
}

pdq_ref_sign <- function(x) {
  radius <- sqrt(rowSums(x^2))
  answer <- matrix(0, nrow(x), ncol(x))
  positive <- radius > 0
  answer[positive, ] <- x[positive, , drop = FALSE] / radius[positive]
  list(sign = answer, radius = radius)
}

pdq_ref_quantile <- function(x, probability) {
  pair.count <- choose(nrow(x), 2)
  rank <- which(seq_len(pair.count) / pair.count >= probability)[1L]
  value <- vapply(seq_len(ncol(x)), function(j) {
    sort(abs(outer(x[, j], x[, j], "-"))[lower.tri(x = matrix(
      0, nrow(x), nrow(x)
    ))])[rank]
  }, numeric(1))
  list(value = value, pair.count = pair.count, rank = rank)
}

pdq_ref_median <- function(x, tol = 1e-12, max_iter = 10000L) {
  location <- colMeans(x)
  for (iteration in seq_len(max_iter)) {
    residual <- sweep(x, 2L, location, "-")
    score <- pdq_ref_sign(residual)
    if (any(score$radius == 0)) {
      stop("pure-R reference fixture hit an observation")
    }
    equation <- sqrt(sum(colMeans(score$sign)^2))
    if (equation <= tol) break
    weights <- 1 / score$radius
    next.location <- colSums(x * weights) / sum(weights)
    location <- next.location
  }
  list(location = location, iterations = iteration, residual = equation)
}

pdq_ref_fit_group <- function(x, probability) {
  quantile <- pdq_ref_quantile(x, probability)
  standardized <- sweep(x, 2L, quantile$value, "/")
  median <- pdq_ref_median(standardized)
  location <- quantile$value * median$location
  residual <- sweep(x, 2L, location, "-")
  residual <- sweep(residual, 2L, quantile$value, "/")
  score <- pdq_ref_sign(residual)
  p <- ncol(x)
  Omega <- crossprod(score$sign) / nrow(x)
  G <- matrix(0, p, p)
  for (i in seq_len(nrow(x))) {
    G <- G + (diag(p) - tcrossprod(score$sign[i, ])) /
      score$radius[i]
  }
  G <- G / nrow(x)
  list(
    quantile = quantile,
    location = location,
    standardized = residual,
    signs = score$sign,
    radii = score$radius,
    Omega = Omega,
    G = G
  )
}

pdq_ref_formula <- function(x, y, probability) {
  fit1 <- pdq_ref_fit_group(x, probability)
  fit2 <- pdq_ref_fit_group(y, probability)
  A12 <- diag(fit2$quantile$value / fit1$quantile$value)
  A21 <- diag(fit1$quantile$value / fit2$quantile$value)
  M1 <- fit2$G %*% A21 %*% solve(fit1$G)
  M2 <- solve(fit2$G) %*% t(A12) %*% fit1$G
  K1 <- (M1 + t(M1)) / 2
  K2 <- (M2 + t(M2)) / 2
  C12 <- solve(fit2$G) %*% t(A12) %*% fit1$G %*%
    fit2$G %*% A21 %*% solve(fit1$G)
  K3 <- diag(ncol(x)) + t(C12)

  cross1 <- sweep(x, 2L, fit2$location, "-")
  cross1 <- sweep(cross1, 2L, fit1$quantile$value, "/")
  cross2 <- sweep(y, 2L, fit1$location, "-")
  cross2 <- sweep(cross2, 2L, fit2$quantile$value, "/")
  cross1 <- pdq_ref_sign(cross1)$sign
  cross2 <- pdq_ref_sign(cross2)$sign
  cross.inner <- cross1 %*% t(cross2)
  R.hat <- -mean(cross.inner)

  mean1 <- colMeans(fit1$signs)
  mean2 <- colMeans(fit2$signs)
  Q.hat <- drop(t(mean1) %*% K1 %*% mean1) +
    drop(t(mean2) %*% K2 %*% mean2) -
    drop(t(mean1) %*% K3 %*% mean2)
  H1 <- fit1$signs %*% K1 %*% t(fit1$signs) / nrow(x)^2
  H2 <- fit2$signs %*% K2 %*% t(fit2$signs) / nrow(y)^2
  H12 <- -fit1$signs %*% K3 %*% t(fit2$signs) /
    (2 * nrow(x) * nrow(y))
  bias <- sum(diag(H1)) + sum(diag(H2))
  H1.zero <- H1
  H2.zero <- H2
  diag(H1.zero) <- diag(H2.zero) <- 0
  tau.sq <- 2 * (sum(H1.zero^2) + sum(H2.zero^2) + 2 * sum(H12^2))
  list(
    fit1 = fit1, fit2 = fit2, A12 = A12, A21 = A21,
    K1 = K1, K2 = K2, K3 = K3,
    cross1 = cross1, cross2 = cross2, cross.inner = cross.inner,
    R.hat = R.hat, Q.hat = Q.hat, H1 = H1, H2 = H2, H12 = H12,
    bias = bias, statistic = R.hat - bias,
    fitted.diagonal.deleted = Q.hat - bias, tau.sq = tau.sq
  )
}

pdq_ref_u32_add <- function(a, b) (a + b) %% 2^32

pdq_ref_u32_multiply <- function(a, b) {
  a0 <- a %% 2^16
  a1 <- floor(a / 2^16)
  b0 <- b %% 2^16
  b1 <- floor(b / 2^16)
  low <- a0 * b0
  middle <- (a0 * b1 + a1 * b0) %% 2^16
  (low + middle * 2^16) %% 2^32
}

pdq_ref_u32_xor <- function(a, b) {
  powers <- 2^(0:31)
  bits.a <- floor(a / powers) %% 2
  bits.b <- floor(b / powers) %% 2
  sum(((bits.a + bits.b) %% 2) * powers)
}

pdq_ref_mix32 <- function(value) {
  value <- pdq_ref_u32_xor(value, floor(value / 2^16))
  value <- pdq_ref_u32_multiply(value, 2146121005)
  value <- pdq_ref_u32_xor(value, floor(value / 2^15))
  value <- pdq_ref_u32_multiply(value, 2221713035)
  pdq_ref_u32_xor(value, floor(value / 2^16))
}

pdq_ref_counter_word <- function(counter, seed) {
  low <- counter %% 2^32
  high <- floor(counter / 2^32)
  high.key <- pdq_ref_mix32(pdq_ref_u32_add(high, 2654435769))
  keyed <- pdq_ref_u32_add(
    pdq_ref_u32_xor(seed, high.key),
    pdq_ref_u32_multiply(2654435769, pdq_ref_u32_add(low, 1))
  )
  pdq_ref_mix32(keyed)
}

pdq_ref_counter_signs <- function(B, dimension, seed) {
  answer <- matrix(-1L, B, dimension)
  for (b in seq_len(B)) {
    for (j in seq_len(dimension)) {
      counter <- (b - 1) * dimension + j - 1
      answer[b, j] <- if (
        pdq_ref_counter_word(counter, seed) %% 2 == 1
      ) 1L else -1L
    }
  }
  answer
}

pdq_ref_bootstrap <- function(reference, multipliers) {
  n1 <- nrow(reference$fit1$signs)
  n2 <- nrow(reference$fit2$signs)
  vapply(seq_len(nrow(multipliers)), function(b) {
    mean1 <- colSums(
      reference$fit1$signs * multipliers[b, seq_len(n1)]
    ) / n1
    mean2 <- colSums(
      reference$fit2$signs * multipliers[b, n1 + seq_len(n2)]
    ) / n2
    drop(t(mean1) %*% reference$K1 %*% mean1) +
      drop(t(mean2) %*% reference$K2 %*% mean2) -
      drop(t(mean1) %*% reference$K3 %*% mean2) - reference$bias
  }, numeric(1))
}

test_that("PDQ observed statistic matches the literal small-sample formula", {
  data <- pdq_test_fixture()
  reference <- pdq_ref_formula(data$x, data$y, 0.25)
  result <- feng_wang_pdq_two_sample_test(
    data$x, data$y, B = 31L, seed = 2605, keep_bootstrap = TRUE
  )

  expect_s3_class(result, "htest")
  expect_s3_class(result, "hd_location_test")
  expect_identical(result$alternative, "two.sided")
  expect_equal(result$components$scales$group1$quantile.input,
               reference$fit1$quantile$value, tolerance = 2e-13,
               ignore_attr = TRUE)
  expect_equal(result$components$scales$group2$quantile.input,
               reference$fit2$quantile$value, tolerance = 2e-13,
               ignore_attr = TRUE)
  expect_equal(result$components$medians$group1$location,
               reference$fit1$location, tolerance = 2e-10,
               ignore_attr = TRUE)
  expect_equal(result$components$medians$group2$location,
               reference$fit2$location, tolerance = 2e-10,
               ignore_attr = TRUE)
  expect_equal(result$components$fitted.signs$group1,
               reference$fit1$signs, tolerance = 2e-10,
               ignore_attr = TRUE)
  expect_equal(result$components$fitted.signs$group2,
               reference$fit2$signs, tolerance = 2e-10,
               ignore_attr = TRUE)
  expect_equal(result$components$matrices$Omega1,
               reference$fit1$Omega, tolerance = 2e-10,
               ignore_attr = TRUE)
  expect_equal(result$components$matrices$Omega2,
               reference$fit2$Omega, tolerance = 2e-10,
               ignore_attr = TRUE)
  expect_equal(result$components$matrices$K1, reference$K1,
               tolerance = 5e-10, ignore_attr = TRUE)
  expect_equal(result$components$matrices$K2, reference$K2,
               tolerance = 5e-10, ignore_attr = TRUE)
  expect_equal(result$components$matrices$K3, reference$K3,
               tolerance = 5e-10, ignore_attr = TRUE)
  expect_equal(result$components$cross.inner.products,
               reference$cross.inner, tolerance = 2e-10,
               ignore_attr = TRUE)
  expect_equal(result$components$observed$R.PDQ, reference$R.hat,
               tolerance = 2e-10)
  expect_equal(result$components$observed$b.hat, reference$bias,
               tolerance = 5e-10)
  expect_equal(unname(result$statistic), reference$statistic,
               tolerance = 5e-10)
  expect_equal(result$components$observed$Q.fitted,
               reference$Q.hat, tolerance = 5e-10)
  expect_false(isTRUE(all.equal(
    result$components$observed$R.PDQ,
    result$components$observed$Q.fitted, tolerance = 1e-8
  )))
  expect_true(result$components$observed$fitted.quadratic.is.not.observed)
})

test_that("U-quantiles use the empirical inverse CDF without interpolation", {
  data <- pdq_test_fixture()
  probability <- 0.4
  q1 <- pdq_ref_quantile(data$x, probability)
  q2 <- pdq_ref_quantile(data$y, probability)
  result <- feng_wang_pdq_two_sample_test(
    data$x, data$y, quantile_prob = probability,
    B = 7L, seed = 11
  )

  expect_equal(q1$pair.count, 28)
  expect_equal(q1$rank, 12)
  expect_equal(q2$rank, 12)
  expect_equal(result$components$scales$group1$pair.count, 28)
  expect_equal(result$components$scales$group2$pair.count, 28)
  expect_equal(result$components$scales$group1$order.index, 12)
  expect_equal(result$components$scales$group2$order.index, 12)
  expect_equal(result$components$scales$group1$quantile.input,
               q1$value, tolerance = 2e-13, ignore_attr = TRUE)
  expect_equal(result$components$scales$group2$quantile.input,
               q2$value, tolerance = 2e-13, ignore_attr = TRUE)
  expect_identical(result$diagnostics$quantile.definition,
                   "unordered-pair empirical inverse CDF; no interpolation")
  expect_equal(result$diagnostics$quantile.order.index,
               c(group1 = 12, group2 = 12))
})

test_that("K matrices, diagonal deletion, and variance retain published factors", {
  data <- pdq_test_fixture()
  reference <- pdq_ref_formula(data$x, data$y, 0.25)
  result <- feng_wang_pdq_two_sample_test(
    data$x, data$y, B = 9L, seed = 19
  )
  deletion <- result$components$diagonal.deletion

  expect_equal(result$components$matrices$A12.diagonal.input,
               diag(reference$A12), tolerance = 2e-13,
               ignore_attr = TRUE)
  expect_equal(result$components$matrices$A21.diagonal.input,
               diag(reference$A21), tolerance = 2e-13,
               ignore_attr = TRUE)
  expect_equal(deletion$within.kernel1, reference$H1,
               tolerance = 5e-10, ignore_attr = TRUE)
  expect_equal(deletion$within.kernel2, reference$H2,
               tolerance = 5e-10, ignore_attr = TRUE)
  expect_equal(deletion$cross.kernel, reference$H12,
               tolerance = 5e-10, ignore_attr = TRUE)
  expect_equal(sum(deletion$contribution1), sum(diag(reference$H1)),
               tolerance = 5e-10)
  expect_equal(sum(deletion$contribution2), sum(diag(reference$H2)),
               tolerance = 5e-10)
  expect_equal(deletion$b.hat, reference$bias, tolerance = 5e-10)
  expect_equal(deletion$within.denominators,
               c(group1 = nrow(data$x)^2, group2 = nrow(data$y)^2))
  expect_equal(unname(result$variance), reference$tau.sq,
               tolerance = 1e-9)
  expect_equal(result$components$bootstrap$variance.diagonal.deleted,
               reference$tau.sq, tolerance = 1e-9)
})

test_that("Rademacher draws and both finite-B decisions are literal", {
  data <- pdq_test_fixture()
  reference <- pdq_ref_formula(data$x, data$y, 0.25)
  B <- 37L
  seed <- 314159
  level <- 0.2
  result <- feng_wang_pdq_two_sample_test(
    data$x, data$y, level = level, B = B, seed = seed,
    keep_bootstrap = TRUE
  )
  expected.signs <- pdq_ref_counter_signs(
    B, nrow(data$x) + nrow(data$y), seed
  )
  expected.bootstrap <- pdq_ref_bootstrap(reference, expected.signs)
  exceedances <- sum(expected.bootstrap >= reference$statistic)
  critical.rank <- which(seq_len(B) / B >= 1 - level)[1L]
  critical <- sort(expected.bootstrap)[critical.rank]

  expect_equal(unname(result$components$bootstrap$multipliers),
               expected.signs)
  expect_equal(unname(result$components$bootstrap$statistics),
               expected.bootstrap, tolerance = 2e-9)
  expect_equal(result$components$bootstrap$exceedances, exceedances)
  expect_equal(result$p.value, (1 + exceedances) / (B + 1),
               tolerance = 0)
  expect_equal(result$components$bootstrap$critical.order.index,
               critical.rank)
  expect_equal(result$components$bootstrap$critical.value, critical,
               tolerance = 2e-9)
  expect_identical(result$components$bootstrap$primary.rejected,
                   reference$statistic > critical)
  expect_identical(result$components$bootstrap$auxiliary.p.rejected,
                   result$p.value <= level)
  expect_identical(result$components$bootstrap$decision.disagreement,
                   xor(reference$statistic > critical,
                       result$p.value <= level))
  expect_true(result$diagnostics$plus.one.correction)
  expect_false(result$diagnostics$paper.finite.B.p.value.specified)
  expect_identical(result$diagnostics$p.value.tie.rule,
                   ">= observed statistic")
  expect_identical(result$diagnostics$primary.tie.rule,
                   "> critical value")
  expect_false(result$null.distribution$finite.sample.exact)
  expect_false(result$null.distribution$nuisance.refitted)
})

test_that("translation, signed coordinate scaling, and observed group swap hold", {
  data <- pdq_test_fixture()
  baseline <- feng_wang_pdq_two_sample_test(
    data$x, data$y, B = 31L, seed = 123, keep_bootstrap = TRUE
  )
  shift <- c(1000, -2000, 3000)
  translated <- feng_wang_pdq_two_sample_test(
    sweep(data$x, 2L, shift, "+"),
    sweep(data$y, 2L, shift, "+"),
    B = 31L, seed = 123, keep_bootstrap = TRUE
  )
  scale <- c(-3.5, 0.2, -7)
  scaled <- feng_wang_pdq_two_sample_test(
    sweep(data$x, 2L, scale, "*"),
    sweep(data$y, 2L, scale, "*"),
    B = 31L, seed = 123, keep_bootstrap = TRUE
  )
  swapped <- feng_wang_pdq_two_sample_test(
    data$y, data$x, B = 31L, seed = 123, keep_bootstrap = TRUE
  )

  expect_equal(unname(translated$statistic), unname(baseline$statistic),
               tolerance = 2e-10)
  expect_equal(translated$p.value, baseline$p.value, tolerance = 0)
  expect_equal(translated$components$bootstrap$statistics,
               baseline$components$bootstrap$statistics,
               tolerance = 2e-10, ignore_attr = TRUE)
  expect_equal(unname(scaled$statistic), unname(baseline$statistic),
               tolerance = 2e-10)
  expect_equal(scaled$p.value, baseline$p.value, tolerance = 0)
  expect_equal(scaled$components$bootstrap$statistics,
               baseline$components$bootstrap$statistics,
               tolerance = 2e-10, ignore_attr = TRUE)
  expect_equal(unname(swapped$statistic), unname(baseline$statistic),
               tolerance = 2e-10)
  expect_equal(swapped$components$observed$R.PDQ,
               baseline$components$observed$R.PDQ, tolerance = 2e-10)
  expect_equal(swapped$components$observed$b.hat,
               baseline$components$observed$b.hat, tolerance = 2e-10)
  expect_equal(swapped$components$matrices$K1,
               baseline$components$matrices$K2,
               tolerance = 2e-9, ignore_attr = TRUE)
  expect_equal(swapped$components$matrices$K2,
               baseline$components$matrices$K1,
               tolerance = 2e-9, ignore_attr = TRUE)
  expect_equal(swapped$components$matrices$K3,
               t(baseline$components$matrices$K3),
               tolerance = 2e-9, ignore_attr = TRUE)
  expect_match(swapped$diagnostics$finite.counter.stream.group.swap,
               "conditional bootstrap law")
})

test_that("explicit seeds isolate R RNG and NULL seeds remain reproducible", {
  data <- pdq_test_fixture()
  set.seed(801)
  state <- .Random.seed
  first <- feng_wang_pdq_two_sample_test(
    data$x, data$y, B = 13L, seed = 4294967295,
    keep_bootstrap = TRUE
  )
  expect_identical(.Random.seed, state)
  second <- feng_wang_pdq_two_sample_test(
    data$x, data$y, B = 13L, seed = 4294967295,
    keep_bootstrap = TRUE
  )
  expect_identical(first$components$bootstrap$multipliers,
                   second$components$bootstrap$multipliers)
  expect_identical(first$components$bootstrap$statistics,
                   second$components$bootstrap$statistics)
  expect_identical(first$p.value, second$p.value)

  set.seed(802)
  null.first <- feng_wang_pdq_two_sample_test(
    data$x, data$y, B = 11L, seed = NULL, keep_bootstrap = TRUE
  )
  set.seed(802)
  null.second <- feng_wang_pdq_two_sample_test(
    data$x, data$y, B = 11L, seed = NULL, keep_bootstrap = TRUE
  )
  expect_identical(null.first$diagnostics$seed.used,
                   null.second$diagnostics$seed.used)
  expect_identical(null.first$components$bootstrap$multipliers,
                   null.second$components$bootstrap$multipliers)
  compact <- feng_wang_pdq_two_sample_test(
    data$x, data$y, B = 5L, seed = 1, keep_bootstrap = FALSE
  )
  expect_null(compact$components$bootstrap$statistics)
  expect_null(compact$components$bootstrap$multipliers)
})

test_that("extreme common units retain the statistic and disclose representation", {
  data <- pdq_test_fixture()
  baseline <- feng_wang_pdq_two_sample_test(
    data$x, data$y, B = 23L, seed = 91, keep_bootstrap = TRUE
  )
  extreme.scale <- c(1e155, -1e-200, 1e80)
  extreme <- feng_wang_pdq_two_sample_test(
    sweep(data$x, 2L, extreme.scale, "*"),
    sweep(data$y, 2L, extreme.scale, "*"),
    B = 23L, seed = 91, keep_bootstrap = TRUE
  )

  expect_true(is.finite(unname(extreme$statistic)))
  expect_equal(unname(extreme$statistic), unname(baseline$statistic),
               tolerance = 5e-10)
  expect_equal(extreme$p.value, baseline$p.value, tolerance = 0)
  expect_equal(extreme$components$bootstrap$statistics,
               baseline$components$bootstrap$statistics,
               tolerance = 5e-10, ignore_attr = TRUE)
  expect_gt(sum(extreme$diagnostics$nonrepresentable.input.diagonals), 0)
  expect_true(all(is.finite(
    extreme$components$scales$group1$log.quantile.input
  )))
  expect_true(all(is.finite(
    extreme$components$scales$group1$D.diagonal.input.canonical
  )))
  expect_identical(extreme$diagnostics$internal.preconditioning,
                   "common coordinatewise midrange/range affine transformation")
  expect_match(extreme$diagnostics$no.repair.policy, "no ridge")
})

test_that("U(0)=0 is limited to cross signs and fitted zeros fail", {
  directions.x <- rbind(
    c(1, 1, 1), c(0.4, -0.7, 0.2),
    c(-0.8, 0.3, 0.6), c(0.2, 0.9, -0.5)
  )
  center.x <- c(1, 1, 1)
  x.zero <- rbind(
    sweep(directions.x, 2L, center.x, "+"),
    sweep(-directions.x, 2L, center.x, "+")
  )
  directions.y <- rbind(
    c(1, 2, 3), c(2, -3, 1), c(-4, 1, 2), c(3, 4, -2)
  )
  y.zero <- 10 * rbind(directions.y, -directions.y)
  cross.zero <- feng_wang_pdq_two_sample_test(
    x.zero, y.zero, B = 7L, seed = 5
  )
  expect_equal(cross.zero$diagnostics$cross.zero.sign.count[["group1"]], 1)
  expect_equal(unname(cross.zero$components$cross.signs$group1[5L, ]),
               c(0, 0, 0), tolerance = 0)
  expect_identical(cross.zero$diagnostics$U.zero.convention,
                   "U(0) = 0 only where no inverse radius is needed")

  x.coincident <- rbind(c(0, 0), c(1, 0.1), c(2, 0.3))
  y.regular <- rbind(
    c(-2, -1), c(-1, 2), c(0.5, -2),
    c(2, 1), c(1, 3), c(-0.5, -3)
  )
  expect_error(
    feng_wang_pdq_two_sample_test(
      x.coincident, y.regular, B = 3L, seed = 1
    ),
    "inverse-radius definition"
  )
})

test_that("nonconvergence and no-repair boundaries are explicit", {
  x <- matrix(c(
    -2.1, 0.2, 1.3,
    -1.0, 2.4, -0.7,
    0.1, -1.7, 2.2,
    1.4, 0.8, -2.5,
    2.8, -0.4, 0.5,
    4.2, 1.6, 3.1
  ), 6, 3, byrow = TRUE)
  y <- matrix(c(
    -3.0, 1.1, -1.2,
    -1.7, -2.2, 0.4,
    -0.2, 2.7, 1.8,
    1.2, -0.8, -2.1,
    2.5, 1.9, 0.9,
    3.8, -1.3, 2.6
  ), 6, 3, byrow = TRUE)
  expect_error(
    feng_wang_pdq_two_sample_test(
      x, y, B = 3L, seed = 1, tol = 1e-14,
      max_iter = 1L, strict = TRUE
    ),
    "did not meet"
  )
  expect_warning(
    loose <- feng_wang_pdq_two_sample_test(
      x, y, B = 3L, seed = 1, tol = 1e-14,
      max_iter = 1L, strict = FALSE
    ),
    "did not meet"
  )
  expect_false(all(loose$diagnostics$median.converged))
  expect_false(loose$diagnostics$median.strict)
  expect_match(loose$diagnostics$no.repair.policy, "generalized inverse")

  data <- pdq_test_fixture()
  constant <- data$x
  constant[, 1L] <- 1
  expect_error(
    feng_wang_pdq_two_sample_test(constant, data$y, B = 3L, seed = 1),
    "U-quantile"
  )
  pooled.constant.x <- data$x
  pooled.constant.y <- data$y
  pooled.constant.x[, 2L] <- pooled.constant.y[, 2L] <- 4
  expect_error(
    feng_wang_pdq_two_sample_test(
      pooled.constant.x, pooled.constant.y, B = 3L, seed = 1
    ),
    "pooled variation"
  )
})

test_that("invalid dimensions and controls fail before calibration", {
  data <- pdq_test_fixture()
  expect_error(
    feng_wang_pdq_two_sample_test(data$x[1:2, ], data$y),
    "at least 3 row"
  )
  expect_error(
    feng_wang_pdq_two_sample_test(data$x, data$y[1:2, ]),
    "at least 3 row"
  )
  expect_error(
    feng_wang_pdq_two_sample_test(data$x[, 1L, drop = FALSE],
                                  data$y[, 1L, drop = FALSE]),
    "at least two variables"
  )
  expect_error(
    feng_wang_pdq_two_sample_test(data$x, data$y[, 1:2]),
    "same number of columns"
  )
  bad <- data$x
  bad[1, 1] <- Inf
  expect_error(feng_wang_pdq_two_sample_test(bad, data$y), "finite")
  expect_error(
    feng_wang_pdq_two_sample_test(data$x, data$y, quantile_prob = 0),
    "quantile_prob"
  )
  expect_error(
    feng_wang_pdq_two_sample_test(data$x, data$y, quantile_prob = 1),
    "quantile_prob"
  )
  expect_error(
    feng_wang_pdq_two_sample_test(data$x, data$y, level = 0),
    "level"
  )
  expect_error(feng_wang_pdq_two_sample_test(data$x, data$y, B = 0),
               "positive integer")
  expect_error(feng_wang_pdq_two_sample_test(data$x, data$y, seed = -1),
               "2\\^32")
  expect_error(feng_wang_pdq_two_sample_test(data$x, data$y, seed = 2^32),
               "2\\^32")
  expect_error(feng_wang_pdq_two_sample_test(data$x, data$y, seed = 1.5),
               "2\\^32")
  expect_error(feng_wang_pdq_two_sample_test(data$x, data$y, tol = 0),
               "finite positive")
  expect_error(feng_wang_pdq_two_sample_test(data$x, data$y, max_iter = 0),
               "positive integer")
  expect_error(
    feng_wang_pdq_two_sample_test(data$x, data$y, strict = NA),
    "TRUE.*FALSE"
  )
  expect_error(
    feng_wang_pdq_two_sample_test(data$x, data$y, keep_bootstrap = 1),
    "TRUE.*FALSE"
  )
})
