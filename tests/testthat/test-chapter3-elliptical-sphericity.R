sph_fixture <- function() {
  structure(
    matrix(c(
      -2.1,  0.4,  1.0,
      -1.3,  1.7, -0.6,
      -0.7, -1.8,  0.2,
       0.1,  0.8,  2.0,
       0.9, -0.5, -1.5,
       1.6,  1.2,  0.4,
       2.4, -1.1,  1.3,
       3.0,  0.2, -0.8
    ), ncol = 3L, byrow = TRUE),
    dimnames = list(paste0("r", 1:8), c("a", "b", "c"))
  )
}

sph_unit <- function(value) {
  maximum <- max(abs(value))
  if (maximum == 0) return(numeric(length(value)))
  scaled <- value / maximum
  scaled / sqrt(sum(scaled^2))
}

sph_sign_reference <- function(x, center) {
  t(vapply(
    seq_len(nrow(x)),
    function(i) sph_unit(x[i, ] - center),
    numeric(ncol(x))
  ))
}

sph_delta_reference <- function(radii) {
  n <- length(radii)
  m1 <- mean(radii^(-1))
  a2 <- mean(radii^(-2)) / m1^2
  a3 <- mean(radii^(-3)) / m1^3
  c(
    delta = (2 - 2 * a2 + a2^2) / n^2 +
      (8 * a2 - 6 * a2^2 + 2 * a2 * a3 - 2 * a3) / n^3,
    ratio2 = a2,
    ratio3 = a3
  )
}

sph_rank_reference <- function(x, method) {
  n <- nrow(x)
  p <- ncol(x)
  directions <- array(0, dim = c(n, n, p))
  for (i in seq_len(n)) {
    for (j in seq_len(n)) {
      directions[i, j, ] <- sph_unit(x[i, ] - x[j, ])
    }
  }
  total <- 0
  for (i in seq_len(n)) for (j in seq_len(n)) {
    for (k in seq_len(n)) for (l in seq_len(n)) {
      if (length(unique(c(i, j, k, l))) != 4L) next
      first <- sum(directions[i, j, ] * directions[k, l, ])
      if (method == "spearman") {
        second <- sum(directions[k, j, ] * directions[i, l, ])
        total <- total + first * second
      } else {
        total <- total + first^2
      }
    }
  }
  denominator <- n * (n - 1) * (n - 2) * (n - 3)
  trace.hat <- if (method == "spearman") {
    total / (2 * denominator)
  } else {
    total / denominator
  }
  q <- if (method == "spearman") {
    4 * p * trace.hat - 1
  } else {
    p * trace.hat - 1
  }
  list(total = total, denominator = denominator,
       trace = trace.hat, q = q)
}

sph_max_reference <- function(signs) {
  n <- nrow(signs)
  p <- ncol(signs)
  psi <- crossprod(signs) / n
  diagonal <- n * p * (p + 2) * (diag(psi) - 1 / p)^2 /
    (2 * (1 - 1 / p))
  off <- matrix(NA_real_, p, p)
  off[upper.tri(off)] <- n * p * (p + 2) *
    psi[upper.tri(psi)]^2
  maximum <- max(diagonal, off, na.rm = TRUE)
  comparisons <- p * (p + 1) / 2
  statistic <- maximum - 2 * log(comparisons) + log(log(comparisons))
  p.value <- -expm1(-exp(-statistic / 2) / sqrt(pi))
  list(psi = psi, diagonal = diagonal, off = off,
       maximum = maximum, statistic = statistic, p.value = p.value)
}


test_that("Zou-Peng-Feng-Wang residual calibration matches a literal formula", {
  x <- sph_fixture()
  result <- zou_peng_feng_wang_sphericity_test(x, tol = 1e-11)
  center <- as.numeric(result$components$fitted.location)
  signs <- sph_sign_reference(unname(x), center)
  radii <- sqrt(rowSums((sweep(unname(x), 2, center))^2))
  delta <- sph_delta_reference(radii)
  n <- nrow(x)
  p <- ncol(x)
  q <- p * (sum((signs %*% t(signs))^2) -
              sum(rowSums(signs^2)^2)) / (n * (n - 1)) - 1
  sigma2 <- 4 * (p - 1) / (n * (n - 1) * (p + 2))
  z <- (q - p * delta[["delta"]]) / sqrt(sigma2)

  expect_s3_class(result, "htest")
  expect_s3_class(result, "hd_sphericity_test")
  expect_equal(unname(result$components$fitted.signs), unname(signs),
               tolerance = 2e-13)
  expect_equal(result$components$Q.tilde, q, tolerance = 2e-13)
  expect_equal(result$components$delta.hat, delta[["delta"]],
               tolerance = 2e-13)
  expect_equal(result$components$inverse.moment.ratios$ratio2,
               delta[["ratio2"]], tolerance = 2e-13)
  expect_equal(result$components$inverse.moment.ratios$ratio3,
               delta[["ratio3"]], tolerance = 2e-13)
  expect_equal(result$components$sigma0.squared, sigma2, tolerance = 1e-15)
  expect_equal(unname(result$statistic), z, tolerance = 2e-13)
  expect_equal(result$p.value, pnorm(z, lower.tail = FALSE),
               tolerance = 2e-15)
  expect_identical(result$diagnostics$bias.estimator, "residual")
  expect_true(result$diagnostics$translation.invariant.calibration)
})


test_that("normal-limit and literal second-order bias paths are locked", {
  x <- sph_fixture() / 7
  normal <- zou_peng_feng_wang_sphericity_test(
    x, bias_estimator = "normal_limit", tol = 1e-11
  )
  n <- nrow(x)
  expect_equal(normal$components$delta.hat, n^(-2) + 2 * n^(-3),
               tolerance = 0)
  expect_null(normal$components$radii.used.scaled)

  second <- zou_peng_feng_wang_sphericity_test(
    x, bias_estimator = "second_order", tol = 1e-11
  )
  center <- as.numeric(second$components$fitted.location)
  signs <- unname(second$components$fitted.signs)
  residual.scale <- second$diagnostics$centering.scale
  extra.scale <- second$diagnostics$second.order.extra.scale
  radii <- sqrt(rowSums((sweep(unname(x), 2, center))^2)) /
    residual.scale
  center.scaled <- center / residual.scale
  corrected <- (radii + as.numeric(signs %*% center.scaled) -
    sum(center.scaled^2) / (2 * radii)) / extra.scale
  delta <- sph_delta_reference(corrected)

  expect_equal(
    as.numeric(second$components$second.order.corrected.radii.scaled),
    as.numeric(corrected), tolerance = 2e-13
  )
  expect_equal(second$components$delta.hat, delta[["delta"]],
               tolerance = 2e-13)
  expect_equal(second$components$inverse.moment.ratios$ratio2,
               delta[["ratio2"]], tolerance = 2e-13)
  expect_equal(second$components$inverse.moment.ratios$ratio3,
               delta[["ratio3"]], tolerance = 2e-13)
  expect_false(second$diagnostics$translation.invariant.calibration)
  expect_match(second$diagnostics$bias.source, "literal 2014")
})


test_that("Feng-Liu ordered sums and denominators match literal quadruple loops", {
  x <- sph_fixture()[1:6, , drop = FALSE]
  for (method in c("spearman", "kendall")) {
    reference <- sph_rank_reference(unname(x), method)
    result <- feng_liu_rank_sphericity_test(
      x, method = method, keep_pair_signs = TRUE
    )
    expect_equal(result$components$ordered.sum, reference$total,
                 tolerance = 2e-12, info = method)
    expect_equal(result$components$ordered.denominator,
                 reference$denominator, tolerance = 0, info = method)
    expect_equal(result$components$trace.estimate, reference$trace,
                 tolerance = 2e-14, info = method)
    expect_equal(result$components$Q.tilde, reference$q,
                 tolerance = 2e-13, info = method)
    expect_equal(nrow(result$components$pair.signs), choose(nrow(x), 2),
                 info = method)
    expect_equal(nrow(result$components$pair.endpoints), choose(nrow(x), 2),
                 info = method)
  }
  spearman <- feng_liu_rank_sphericity_test(x, "spearman")
  kendall <- feng_liu_rank_sphericity_test(x, "kendall")
  expect_equal(spearman$components$trace.denominator,
               2 * spearman$components$ordered.denominator, tolerance = 0)
  expect_equal(kendall$components$trace.denominator,
               kendall$components$ordered.denominator, tolerance = 0)
  expect_null(spearman$components$pair.signs)
})


test_that("spatial-sign max statistic and Gumbel tail match a literal SSCM", {
  x <- sph_fixture()
  result <- zhao_yang_zhang_feng_wang_sign_max_test(
    x, tol = 1e-11, keep_sscm = TRUE
  )
  reference <- sph_max_reference(unname(result$components$fitted.signs))

  expect_s3_class(result, "htest")
  expect_equal(unname(result$components$sscm), unname(reference$psi),
               tolerance = 2e-14)
  expect_equal(result$components$maximum.standardized.square,
               reference$maximum, tolerance = 2e-13)
  expect_equal(result$components$T.SM, reference$statistic,
               tolerance = 2e-13)
  expect_equal(result$p.value, reference$p.value, tolerance = 2e-15)
  expect_equal(result$components$p.SM, reference$p.value, tolerance = 2e-15)
  expect_equal(result$components$comparison.count,
               ncol(x) * (ncol(x) + 1) / 2, tolerance = 0)
  expect_equal(result$diagnostics$sscm.trace,
               sum(diag(reference$psi)), tolerance = 2e-15)
  expect_identical(
    result$components$reject,
    result$components$T.SM >= result$components$critical.value
  )
})


test_that("adaptive output separates Cauchy score from the final probability", {
  x <- sph_fixture()
  result <- zhao_yang_zhang_feng_wang_adaptive_sphericity_test(
    x, tol = 1e-11, keep_sscm = TRUE
  )
  p.ss <- result$components$sum$p.SS
  p.sm <- result$components$max$p.SM
  terms <- ifelse(c(p.ss, p.sm) < 0.5,
                  0.5 / tan(pi * c(p.ss, p.sm)), 0)
  score <- sum(terms)
  combined <- atan2(1, score) / pi

  expect_equal(unname(result$statistic), score, tolerance = 2e-14)
  expect_equal(result$components$Cauchy.score, score, tolerance = 2e-14)
  expect_equal(result$components$combined.p.value, combined,
               tolerance = 2e-15)
  expect_equal(result$p.value, combined, tolerance = 2e-15)
  expect_equal(unname(result$components$Cauchy.terms), unname(terms),
               tolerance = 2e-14)
  expect_equal(result$components$sum$Q.tilde,
               zou_peng_feng_wang_sphericity_test(x, tol = 1e-11)$components$Q.tilde,
               tolerance = 2e-13)
  expect_equal(result$components$max$T.SM,
               zhao_yang_zhang_feng_wang_sign_max_test(x, tol = 1e-11)$components$T.SM,
               tolerance = 2e-13)
  expect_match(result$diagnostics$cauchy.output.contract, "p.value is 1 minus")

  empty <- .sph_truncated_cauchy(c(0.5, 0.9))
  expect_equal(empty$score, 0, tolerance = 0)
  expect_equal(empty$p.value, 0.5, tolerance = 0)
  expect_true(empty$empty.active.set)
  expect_identical(empty$active, c(FALSE, FALSE))

  boundary <- .sph_truncated_cauchy(c(0, 0.8))
  expect_true(is.infinite(boundary$score))
  expect_equal(boundary$p.value, 0, tolerance = 0)
})


test_that("elliptical procedures have the stated similarity invariances", {
  x <- sph_fixture()
  shift <- c(10, -7, 3)
  rotation <- qr.Q(qr(matrix(c(
    2, -1, 0,
    1,  2, 1,
    0, -1, 2
  ), 3, 3, byrow = TRUE)))
  transformed <- sweep(unname(x) %*% rotation * -3.5, 2, shift, "+")
  signed.permuted <- sweep(
    unname(x)[, c(3, 1, 2)] * matrix(c(-3.5, 3.5, -3.5),
                                     nrow(x), 3, byrow = TRUE),
    2, shift, "+"
  )

  sign.base <- zou_peng_feng_wang_sphericity_test(x, tol = 1e-11)
  sign.changed <- zou_peng_feng_wang_sphericity_test(
    transformed, tol = 1e-11
  )
  max.base <- zhao_yang_zhang_feng_wang_sign_max_test(x, tol = 1e-11)
  max.changed <- zhao_yang_zhang_feng_wang_sign_max_test(
    signed.permuted, tol = 1e-11
  )
  adaptive.base <- zhao_yang_zhang_feng_wang_adaptive_sphericity_test(
    x, tol = 1e-11
  )
  adaptive.changed <- zhao_yang_zhang_feng_wang_adaptive_sphericity_test(
    signed.permuted, tol = 1e-11
  )

  expect_equal(sign.changed$statistic, sign.base$statistic, tolerance = 3e-8)
  expect_equal(max.changed$statistic, max.base$statistic, tolerance = 3e-8)
  expect_equal(adaptive.changed$p.value, adaptive.base$p.value,
               tolerance = 3e-9)
  rotated.max <- zhao_yang_zhang_feng_wang_sign_max_test(
    unname(x) %*% rotation, tol = 1e-11
  )
  expect_gt(abs(unname(rotated.max$statistic - max.base$statistic)), 1e-3)
  for (method in c("spearman", "kendall")) {
    base <- feng_liu_rank_sphericity_test(x, method)
    changed <- feng_liu_rank_sphericity_test(transformed, method)
    permuted <- feng_liu_rank_sphericity_test(x[c(8:1), ], method)
    expect_equal(changed$statistic, base$statistic, tolerance = 2e-12,
                 info = method)
    expect_equal(permuted$statistic, base$statistic, tolerance = 2e-12,
                 info = method)
  }

  row.permuted <- zou_peng_feng_wang_sphericity_test(
    x[c(3, 8, 1, 7, 2, 6, 4, 5), ], tol = 1e-11
  )
  expect_equal(row.permuted$statistic, sign.base$statistic, tolerance = 3e-8)
})


test_that("coordinatewise scaling is correctly not treated as a symmetry", {
  x <- sph_fixture()
  changed <- sweep(unname(x), 2, c(0.2, 3, 8), "*")
  sign.base <- zou_peng_feng_wang_sphericity_test(x, tol = 1e-11)
  sign.changed <- zou_peng_feng_wang_sphericity_test(changed, tol = 1e-11)
  rank.base <- feng_liu_rank_sphericity_test(x, "kendall")
  rank.changed <- feng_liu_rank_sphericity_test(changed, "kendall")

  expect_gt(abs(unname(sign.base$statistic - sign.changed$statistic)), 1e-3)
  expect_gt(abs(unname(rank.base$statistic - rank.changed$statistic)), 1e-3)
})


test_that("ties, zero residuals, and no-repair contracts are explicit", {
  tied <- rbind(
    c(0, 0), c(0, 0), c(1, 0), c(0, 1), c(-1, 0), c(0, -2)
  )
  rank <- feng_liu_rank_sphericity_test(tied, "kendall")
  reference <- sph_rank_reference(tied, "kendall")
  expect_equal(rank$diagnostics$zero.pair.directions, 1, tolerance = 0)
  expect_equal(rank$components$ordered.denominator,
               nrow(tied) * 5 * 4 * 3, tolerance = 0)
  expect_equal(rank$components$ordered.sum, reference$total,
               tolerance = 2e-12)

  coincident <- matrix(rep(c(2, -1), 6), ncol = 2, byrow = TRUE)
  expect_error(
    zou_peng_feng_wang_sphericity_test(coincident),
    "Inverse radial moments are undefined"
  )
  expect_error(
    zhao_yang_zhang_feng_wang_adaptive_sphericity_test(coincident),
    "Inverse radial moments are undefined"
  )
  normal <- zou_peng_feng_wang_sphericity_test(
    coincident, bias_estimator = "normal_limit"
  )
  maximum <- zhao_yang_zhang_feng_wang_sign_max_test(coincident)
  expect_true(is.finite(unname(normal$statistic)))
  expect_equal(normal$diagnostics$zero.residuals, nrow(coincident), tolerance = 0)
  expect_equal(maximum$diagnostics$zero.residuals, nrow(coincident), tolerance = 0)
  expect_equal(maximum$diagnostics$sscm.trace, 0, tolerance = 0)
  expect_identical(maximum$diagnostics$ridge, "none")
  expect_identical(maximum$diagnostics$numerical.floor, "none")
})


test_that("extreme common scales remain finite and preserve results", {
  x <- sph_fixture()
  sign.base <- zou_peng_feng_wang_sphericity_test(x, tol = 1e-11)
  rank.base <- feng_liu_rank_sphericity_test(x, "spearman")
  max.base <- zhao_yang_zhang_feng_wang_sign_max_test(x, tol = 1e-11)

  for (scale in c(1e250, 1e-250)) {
    sign.scaled <- zou_peng_feng_wang_sphericity_test(x * scale, tol = 1e-11)
    rank.scaled <- feng_liu_rank_sphericity_test(x * scale, "spearman")
    max.scaled <- zhao_yang_zhang_feng_wang_sign_max_test(
      x * scale, tol = 1e-11
    )
    expect_true(is.finite(unname(sign.scaled$statistic)))
    expect_true(is.finite(unname(rank.scaled$statistic)))
    expect_true(is.finite(unname(max.scaled$statistic)))
    expect_equal(sign.scaled$statistic, sign.base$statistic, tolerance = 3e-8)
    expect_equal(rank.scaled$statistic, rank.base$statistic, tolerance = 2e-12)
    expect_equal(max.scaled$statistic, max.base$statistic, tolerance = 3e-8)
  }
})


test_that("input, convergence, and second-order failures are diagnosed", {
  x <- sph_fixture()
  expect_error(zou_peng_feng_wang_sphericity_test(x[, 1, drop = FALSE]),
               "at least two columns")
  expect_error(feng_liu_rank_sphericity_test(x[1:3, ]), "at least 4")
  expect_error(feng_liu_rank_sphericity_test(x, "other"), "arg")
  expect_error(zhao_yang_zhang_feng_wang_sign_max_test(x, alpha = 0),
               "strictly between")
  expect_error(zhao_yang_zhang_feng_wang_sign_max_test(x, keep_sscm = NA),
               "TRUE or FALSE")
  expect_error(zou_peng_feng_wang_sphericity_test(x, tol = 0),
               "finite positive")
  expect_error(zou_peng_feng_wang_sphericity_test(x, max_iter = 0),
               "positive integer")
  bad <- x
  bad[1, 1] <- Inf
  expect_error(feng_liu_rank_sphericity_test(bad), "finite")

  expect_error(
    zhao_yang_zhang_feng_wang_sign_max_test(
      x, tol = 1e-16, max_iter = 1L, strict = TRUE
    ),
    "did not converge"
  )
  expect_warning(
    loose <- zhao_yang_zhang_feng_wang_sign_max_test(
      x, tol = 1e-16, max_iter = 1L, strict = FALSE
    ),
    "last finite iterate"
  )
  expect_false(loose$diagnostics$spatial.median$converged)

  far.from.origin <- sweep(unname(x), 2, c(100, -70, 30), "+")
  expect_error(
    zou_peng_feng_wang_sphericity_test(
      far.from.origin, bias_estimator = "second_order", tol = 1e-11
    ),
    "second-order corrected radii"
  )
})
