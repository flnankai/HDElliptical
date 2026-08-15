composite_greedy_reference <- function(covariance, block_size) {
  correlation <- stats::cov2cor(covariance)
  remaining <- seq_len(ncol(covariance))
  blocks <- list()
  while (length(remaining) > 0L) {
    if (length(remaining) <= block_size) {
      blocks[[length(blocks) + 1L]] <- remaining
      break
    }
    if (block_size == 1L) {
      selected <- remaining[[1L]]
    } else {
      pairs <- utils::combn(remaining, 2L)
      pair.score <- abs(correlation[cbind(pairs[1L, ], pairs[2L, ])])
      selected <- pairs[, which.max(pair.score)]
      while (length(selected) < block_size) {
        candidates <- setdiff(remaining, selected)
        addition.score <- vapply(candidates, function(candidate) {
          sum(abs(correlation[candidate, selected]))
        }, numeric(1L))
        selected <- c(selected, candidates[[which.max(addition.score)]])
      }
    }
    selected <- sort(selected)
    blocks[[length(blocks) + 1L]] <- selected
    remaining <- setdiff(remaining, selected)
  }
  blocks
}

composite_block_inverse_reference <- function(covariance, block_size) {
  blocks <- composite_greedy_reference(covariance, block_size)
  inverse <- matrix(0, nrow(covariance), ncol(covariance))
  for (block in blocks) {
    inverse[block, block] <- solve(
      covariance[block, block, drop = FALSE]
    )
  }
  list(inverse = inverse, blocks = blocks)
}

composite_scatter_reference <- function(z) {
  centered <- sweep(z, 2L, colMeans(z), "-")
  crossprod(centered)
}

composite_pooled_reference <- function(x, y) {
  (composite_scatter_reference(x) + composite_scatter_reference(y)) /
    (nrow(x) + nrow(y) - 2)
}

composite_q_reference <- function(x, y, block_size) {
  n1 <- nrow(x)
  n2 <- nrow(y)
  ordered.numerator <- 0
  for (i1 in seq_len(n1)) {
    for (i2 in setdiff(seq_len(n1), i1)) {
      keep1 <- setdiff(seq_len(n1), c(i1, i2))
      for (j1 in seq_len(n2)) {
        for (j2 in setdiff(seq_len(n2), j1)) {
          keep2 <- setdiff(seq_len(n2), c(j1, j2))
          covariance <- composite_pooled_reference(
            x[keep1, , drop = FALSE], y[keep2, , drop = FALSE]
          )
          inverse <- composite_block_inverse_reference(
            covariance, block_size
          )$inverse
          ordered.numerator <- ordered.numerator + as.numeric(
            crossprod(x[i1, ] - y[j1, ],
                      inverse %*% (x[i2, ] - y[j2, ]))
          )
        }
      }
    }
  }
  ordered.denominator <- n1 * (n1 - 1) * n2 * (n2 - 1)
  list(
    estimate = ordered.numerator / ordered.denominator,
    ordered.numerator = ordered.numerator,
    ordered.denominator = ordered.denominator
  )
}

composite_trace_reference <- function(x, block_size) {
  n <- nrow(x)
  ordered.numerator <- 0
  for (i1 in seq_len(n)) {
    for (i2 in setdiff(seq_len(n), i1)) {
      for (i3 in setdiff(seq_len(n), c(i1, i2))) {
        for (i4 in setdiff(seq_len(n), c(i1, i2, i3))) {
          keep <- setdiff(seq_len(n), c(i1, i2, i3, i4))
          covariance <- stats::cov(x[keep, , drop = FALSE])
          inverse <- composite_block_inverse_reference(
            covariance, block_size
          )$inverse
          first <- as.numeric(crossprod(
            x[i1, ] - x[i2, ], inverse %*% (x[i3, ] - x[i4, ])
          ))
          second <- as.numeric(crossprod(
            x[i1, ] - x[i4, ], inverse %*% (x[i3, ] - x[i2, ])
          ))
          ordered.numerator <- ordered.numerator + first * second
        }
      }
    }
  }
  P4 <- n * (n - 1) * (n - 2) * (n - 3)
  list(
    estimate = ordered.numerator / (2 * P4),
    ordered.numerator = ordered.numerator,
    ordered.denominator = 2 * P4,
    P4 = P4
  )
}


test_that("Composite T2 matches the literal published two-sample formulas", {
  set.seed(2701)
  x <- matrix(stats::rnorm(24), 8, 3)
  y <- matrix(stats::rnorm(18, 0.25, 1.2), 6, 3)
  q.reference <- composite_q_reference(x, y, 2L)
  trace.reference <- composite_trace_reference(x, 2L)
  result <- composite_t2_two_sample_test(x, y, block_size = 2L)
  coefficient <- 2 * (1 / nrow(x) + 1 / nrow(y))^2
  variance <- coefficient * trace.reference$estimate
  z <- q.reference$estimate / sqrt(variance)

  expect_equal(result$components$Q.n, q.reference$estimate,
               tolerance = 2e-10)
  expect_equal(result$components$Q.ordered.numerator,
               q.reference$ordered.numerator, tolerance = 3e-9)
  expect_equal(result$components$Q.ordered.denominator,
               q.reference$ordered.denominator)
  expect_equal(result$components$trace.Lambda.K.squared,
               trace.reference$estimate, tolerance = 3e-9)
  expect_equal(result$components$trace.ordered.numerator,
               trace.reference$ordered.numerator, tolerance = 5e-8)
  expect_equal(result$components$trace.ordered.denominator,
               trace.reference$ordered.denominator)
  expect_equal(result$components$P4.n1, trace.reference$P4)
  expect_equal(result$components$variance.coefficient, coefficient,
               tolerance = 1e-15)
  expect_equal(result$variance[[1L]], variance, tolerance = 3e-9)
  expect_equal(unname(result$statistic), z, tolerance = 3e-10)
  expect_equal(result$p.value, stats::pnorm(z, lower.tail = FALSE),
               tolerance = 1e-14)
})


test_that("Composite T2 exposes exact ordered denominators and counts", {
  set.seed(2702)
  x <- matrix(stats::rnorm(21), 7, 3)
  y <- matrix(stats::rnorm(15), 5, 3)
  result <- composite_t2_two_sample_test(x, y)

  expect_equal(result$components$Q.ordered.denominator,
               7 * 6 * 5 * 4)
  expect_equal(result$components$Q.unordered.pair.combinations,
               choose(7, 2) * choose(5, 2))
  expect_equal(result$components$trace.ordered.denominator,
               2 * 7 * 6 * 5 * 4)
  expect_equal(result$components$trace.unordered.quadruples,
               choose(7, 4))
  expect_equal(result$components$P4.n1, 7 * 6 * 5 * 4)
  expect_equal(
    result$components$Q.n,
    result$components$Q.ordered.numerator /
      result$components$Q.ordered.denominator,
    tolerance = 1e-14
  )
  expect_equal(
    result$components$trace.Lambda.K.squared,
    result$components$trace.ordered.numerator /
      result$components$trace.ordered.denominator,
    tolerance = 1e-14
  )
  expect_equal(
    result$components$denominator.variance,
    2 * (1 / 7 + 1 / 5)^2 *
      result$components$trace.Lambda.K.squared,
    tolerance = 1e-14
  )
})


test_that("Composite T2 reports the paper greedy partition and remainder", {
  set.seed(2703)
  x <- matrix(stats::rnorm(32), 8, 4)
  y <- matrix(stats::rnorm(24), 6, 4)
  x[, 2] <- 0.85 * x[, 1] + 0.35 * x[, 2]
  y[, 2] <- 0.85 * y[, 1] + 0.35 * y[, 2]
  result <- composite_t2_two_sample_test(x, y, block_size = 3L)
  pooled <- composite_pooled_reference(x, y)
  expected <- composite_greedy_reference(pooled, 3L)

  expect_equal(
    unname(lapply(result$components$full.sample.blocks, unname)),
    expected
  )
  expect_equal(unname(lengths(result$components$full.sample.blocks)),
               c(3L, 1L))
  expect_identical(result$components$selection, "paper_greedy")
  expect_match(result$diagnostics$block.selection, "lexicographic")
  expect_true(result$diagnostics$dynamic.block.selection)
  expect_gte(result$diagnostics$q.unique.partitions, 1)
  expect_gte(result$diagnostics$trace.unique.partitions, 1)
})


test_that("Composite T2 keeps the paper's group-1 trace calibration", {
  set.seed(2704)
  x <- matrix(stats::rnorm(24, sd = 0.8), 8, 3)
  y <- matrix(stats::rnorm(27, 0.2, 1.5), 9, 3)
  forward <- composite_t2_two_sample_test(x, y)
  reverse <- composite_t2_two_sample_test(y, x)

  expect_equal(forward$components$Q.n, reverse$components$Q.n,
               tolerance = 2e-10)
  expect_equal(reverse$estimate, -forward$estimate, tolerance = 1e-14)
  expect_equal(forward$components$variance.coefficient,
               reverse$components$variance.coefficient,
               tolerance = 1e-15)
  expect_gt(abs(forward$components$trace.Lambda.K.squared -
                  reverse$components$trace.Lambda.K.squared), 1e-8)
  expect_identical(forward$diagnostics$trace.calibration.sample,
                   "group1 only")
  expect_false(forward$diagnostics$finite.sample.group.label.symmetric)
})


test_that("Composite T2 is translation and diagonal-scale invariant", {
  set.seed(2705)
  x <- matrix(stats::rnorm(21), 7, 3)
  y <- matrix(stats::rnorm(15, 0.1, 1.3), 5, 3)
  baseline <- composite_t2_two_sample_test(x, y)
  shift <- c(30, -17, 4.5)
  translated <- composite_t2_two_sample_test(
    sweep(x, 2L, shift, "+"), sweep(y, 2L, shift, "+")
  )
  units <- c(1e200, -1e-200, 3e150)
  rescaled <- composite_t2_two_sample_test(
    sweep(x, 2L, units, "*"), sweep(y, 2L, units, "*")
  )
  permuted <- composite_t2_two_sample_test(
    x[c(6, 2, 7, 1, 4, 3, 5), , drop = FALSE],
    y[c(4, 1, 5, 2, 3), , drop = FALSE]
  )
  combined.maximum <- apply(rbind(x, y), 2L, function(z) max(abs(z)))
  bounded.x <- sweep(x, 2L, 2 * combined.maximum, "/")
  bounded.y <- sweep(y, 2L, 2 * combined.maximum, "/")
  bounded <- composite_t2_two_sample_test(bounded.x, bounded.y)
  extreme.units <- c(1.7e308, -1e-308, 8e250)
  extreme <- composite_t2_two_sample_test(
    sweep(bounded.x, 2L, extreme.units, "*"),
    sweep(bounded.y, 2L, extreme.units, "*")
  )

  expect_lte(
    abs(unname(translated$statistic - baseline$statistic)),
    3e-10
  )
  expect_equal(translated$components$Q.n, baseline$components$Q.n,
               tolerance = 2e-10)
  expect_lte(
    abs(unname(rescaled$statistic - baseline$statistic)),
    3e-10
  )
  expect_equal(rescaled$components$trace.Lambda.K.squared,
               baseline$components$trace.Lambda.K.squared,
               tolerance = 3e-9)
  expect_equal(rescaled$components$full.sample.blocks,
               baseline$components$full.sample.blocks)
  expect_lte(
    abs(unname(permuted$statistic - baseline$statistic)),
    4e-10
  )
  expect_equal(permuted$p.value, baseline$p.value, tolerance = 1e-12)
  expect_equal(extreme$statistic, bounded$statistic, tolerance = 4e-10)
  expect_equal(extreme$components$Q.n, bounded$components$Q.n,
               tolerance = 3e-10)
})


test_that("Composite T2 supports K one and the minimum scalar sizes", {
  x <- matrix(c(-3, -1, 0.5, 2, 4, 7), ncol = 1)
  y <- matrix(c(-2, 1, 5), ncol = 1)
  result <- composite_t2_two_sample_test(x, y, block_size = 1L)
  q.reference <- composite_q_reference(x, y, 1L)
  trace.reference <- composite_trace_reference(x, 1L)

  expect_s3_class(result, "hd_location_test")
  expect_s3_class(result, "htest")
  expect_equal(result$components$Q.n, q.reference$estimate,
               tolerance = 2e-12)
  expect_equal(result$components$trace.Lambda.K.squared,
               trace.reference$estimate, tolerance = 2e-11)
  expect_equal(result$n, c(group1 = 6, group2 = 3))
  expect_equal(result$p, 1L)
  expect_equal(unname(lengths(result$components$full.sample.blocks)), 1L)
})


test_that("Composite T2 returns a documented htest contract", {
  set.seed(2706)
  colnames.x <- c("alpha", "beta", "gamma")
  x <- matrix(stats::rnorm(21), 7, 3, dimnames = list(NULL, colnames.x))
  y <- matrix(stats::rnorm(18), 6, 3, dimnames = list(NULL, colnames.x))
  result <- composite_t2_two_sample_test(x, y)

  expect_named(result$statistic, "Z")
  expect_named(result$raw.statistic, "Q.n")
  expect_equal(names(result$estimate), colnames.x)
  expect_identical(result$null.distribution$family, "normal")
  expect_identical(result$null.distribution$tail, "upper")
  expect_identical(result$diagnostics$covariance.model, "common")
  expect_false(result$diagnostics$behrens.fisher)
  expect_false(result$diagnostics$component.test.combination)
  expect_identical(result$diagnostics$component.correlation.calibration,
                   "not applicable")
  expect_false(grepl("Behrens", result$method, fixed = TRUE))
  expect_true(all(is.finite(
    result$diagnostics$minimum.reciprocal.condition
  )))
  expect_true(all(result$diagnostics$minimum.reciprocal.condition > 0))
  expect_identical(result$diagnostics$regularization, "none")
  expect_identical(result$diagnostics$variance.repair, "none")
})


test_that("Composite T2 validates inputs and paper size conditions", {
  set.seed(2707)
  x <- matrix(stats::rnorm(21), 7, 3)
  y <- matrix(stats::rnorm(15), 5, 3)

  expect_error(composite_t2_two_sample_test(letters[1:7], y),
               "numeric matrix")
  expect_error(composite_t2_two_sample_test(x, y[, 1:2]),
               "same number of columns")
  x.named <- x
  y.named <- y
  colnames(x.named) <- c("a", "b", "c")
  colnames(y.named) <- c("a", "c", "b")
  expect_error(composite_t2_two_sample_test(x.named, y.named),
               "same names")
  x.bad <- x
  x.bad[1, 1] <- Inf
  expect_error(composite_t2_two_sample_test(x.bad, y), "finite")
  expect_error(composite_t2_two_sample_test(x, y, block_size = 0),
               "positive integer")
  expect_error(composite_t2_two_sample_test(x, y, block_size = 1.5),
               "positive integer")
  expect_error(composite_t2_two_sample_test(x, y, block_size = 4),
               "must not exceed")
  expect_error(composite_t2_two_sample_test(x, y, selection = "exact"),
               "paper_greedy")
  expect_error(composite_t2_two_sample_test(x[1:6, ], y),
               "nrow\\(x\\) >= block_size \\+ 5")
  expect_error(composite_t2_two_sample_test(x, y[1:2, ]),
               "at least 3 row")
})


test_that("Composite T2 rejects degenerate covariance blocks without repair", {
  varying <- seq_len(8)
  x.zero <- cbind(varying[1:7], constant = 1)
  y.zero <- cbind(varying[1:5] + 0.2, constant = 1)
  expect_error(
    composite_t2_two_sample_test(x.zero, y.zero, block_size = 1L),
    "marginal variance.*strictly positive"
  )

  x.base <- c(-3, -2, -1, 0, 1, 2, 4)
  y.base <- c(-2.5, -1.5, 0.5, 1.5, 3)
  x.singular <- cbind(x.base, 2 * x.base)
  y.singular <- cbind(y.base, 2 * y.base)
  colnames(x.singular) <- colnames(y.singular) <- c("a", "b")
  expect_error(
    composite_t2_two_sample_test(x.singular, y.singular, block_size = 2L),
    "positive definite|numerically invertible"
  )
})
