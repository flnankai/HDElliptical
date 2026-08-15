fzl_fixture <- function() {
  list(
    x = rbind(
      c(0.7, -1.1, 0.2),
      c(-0.4, 0.8, 1.2),
      c(1.1, 0.3, -0.9),
      c(-1.2, -0.5, 0.7),
      c(0.2, 1.4, -0.4),
      c(1.5, -0.2, 0.5)
    ),
    y = rbind(
      c(-0.3, 0.7, -0.1),
      c(1.2, -0.6, 0.8),
      c(-1.1, 0.2, 1.3),
      c(0.5, 1.1, -0.7),
      c(0.9, -1.3, 0.4),
      c(-0.6, -0.4, -1.2)
    )
  )
}


fzl_ref_direction <- function(difference, diagonal) {
  standardized <- difference / sqrt(diagonal)
  radius <- sqrt(sum(standardized^2))
  if (radius == 0) {
    return(rep(0, length(difference)))
  }
  standardized / radius
}


fzl_ref_prepare <- function(x, y, identification) {
  p <- ncol(x)
  anchor <- x[1L, ]
  global <- max(abs(x), abs(y))
  x.out <- y.out <- NULL
  x.out <- matrix(0, nrow(x), p)
  y.out <- matrix(0, nrow(y), p)
  log.scale <- coordinate.range <- numeric(p)
  for (j in seq_len(p)) {
    operand <- if (identification == "geometric") {
      max(abs(x[, j]), abs(y[, j]))
    } else {
      global
    }
    x.out[, j] <- x[, j] / operand - anchor[j] / operand
    y.out[, j] <- y[, j] / operand - anchor[j] / operand
    coordinate.range[j] <- max(abs(x.out[, j]), abs(y.out[, j]))
    if (identification == "geometric") {
      x.out[, j] <- x.out[, j] / coordinate.range[j]
      y.out[, j] <- y.out[, j] / coordinate.range[j]
      log.scale[j] <- log(operand) + log(coordinate.range[j])
    } else {
      log.scale[j] <- log(global)
    }
  }
  list(
    x = x.out,
    y = y.out,
    anchor = anchor,
    log.scale = log.scale,
    coordinate.range = coordinate.range
  )
}


fzl_ref_identify <- function(diagonal, identification) {
  if (identification == "geometric") {
    diagonal / exp(mean(log(diagonal)))
  } else {
    ncol <- length(diagonal)
    ncol * diagonal / sum(diagonal)
  }
}


fzl_ref_fit <- function(data, exclude = integer(), identification,
                        tol = 1e-7, max_iter = 500L) {
  keep <- setdiff(seq_len(nrow(data)), exclude)
  p <- ncol(data)
  m <- length(keep)
  diagonal <- apply(data[keep, , drop = FALSE], 2L, stats::var)
  diagonal <- fzl_ref_identify(diagonal, identification)

  for (iteration in 0:max_iter) {
    ranks <- matrix(0, m, p)
    zero.pairs <- 0L
    if (m >= 2L) {
      for (a in seq_len(m - 1L)) {
        for (b in (a + 1L):m) {
          direction <- fzl_ref_direction(
            data[keep[a], ] - data[keep[b], ], diagonal
          )
          if (all(direction == 0)) zero.pairs <- zero.pairs + 1L
          ranks[a, ] <- ranks[a, ] + direction / m
          ranks[b, ] <- ranks[b, ] - direction / m
        }
      }
    }
    score <- colMeans(ranks^2)
    residual <- max(abs(p * score / sum(score) - 1))
    if (residual <= tol) {
      return(list(
        diagonal = diagonal,
        log.diagonal = log(diagonal),
        iterations = iteration,
        residual = residual,
        zero.pairs = zero.pairs
      ))
    }
    if (iteration == max_iter) stop("reference did not converge")
    diagonal <- fzl_ref_identify(diagonal * score, identification)
  }
  stop("unreachable")
}


fzl_ref_permutations4 <- function() {
  output <- vector("list", 24L)
  index <- 1L
  for (a in 1:4) {
    for (b in setdiff(1:4, a)) {
      for (cc in setdiff(1:4, c(a, b))) {
        d <- setdiff(1:4, c(a, b, cc))
        output[[index]] <- c(a, b, cc, d)
        index <- index + 1L
      }
    }
  }
  do.call(rbind, output)
}


fzl_ref_within_trace <- function(data, quads, fits) {
  p <- ncol(data)
  permutations <- fzl_ref_permutations4()
  contributions <- matrix(0, nrow(quads), 24L)
  for (q in seq_len(nrow(quads))) {
    diagonal <- fits[[q]]$diagonal
    for (r in seq_len(24L)) {
      index <- quads[q, permutations[r, ]]
      u12 <- fzl_ref_direction(
        data[index[1L], ] - data[index[2L], ], diagonal
      )
      u34 <- fzl_ref_direction(
        data[index[3L], ] - data[index[4L], ], diagonal
      )
      u32 <- fzl_ref_direction(
        data[index[3L], ] - data[index[2L], ], diagonal
      )
      u14 <- fzl_ref_direction(
        data[index[1L], ] - data[index[4L], ], diagonal
      )
      contributions[q, r] <- sum(u12 * u34) * sum(u32 * u14)
    }
  }
  ordered.sum <- sum(contributions)
  n <- nrow(data)
  list(
    estimate = 2 * p^2 * ordered.sum /
      (n * (n - 1) * (n - 2) * (n - 3)),
    ordered.sum = ordered.sum,
    contributions = contributions
  )
}


fzl_reference <- function(x, y, identification = "geometric",
                          tol = 1e-7, max_iter = 500L) {
  prepared <- fzl_ref_prepare(x, y, identification)
  x.data <- prepared$x
  y.data <- prepared$y
  n1 <- nrow(x)
  n2 <- nrow(y)
  p <- ncol(x)
  pairs1 <- t(utils::combn(n1, 2L))
  pairs2 <- t(utils::combn(n2, 2L))
  quads1 <- t(utils::combn(n1, 4L))
  quads2 <- t(utils::combn(n2, 4L))
  full1 <- fzl_ref_fit(x.data, identification = identification,
                       tol = tol, max_iter = max_iter)
  full2 <- fzl_ref_fit(y.data, identification = identification,
                       tol = tol, max_iter = max_iter)
  pair.fits1 <- lapply(seq_len(nrow(pairs1)), function(index) {
    fzl_ref_fit(x.data, pairs1[index, ], identification, tol, max_iter)
  })
  pair.fits2 <- lapply(seq_len(nrow(pairs2)), function(index) {
    fzl_ref_fit(y.data, pairs2[index, ], identification, tol, max_iter)
  })
  quad.fits1 <- lapply(seq_len(nrow(quads1)), function(index) {
    fzl_ref_fit(x.data, quads1[index, ], identification, tol, max_iter)
  })
  quad.fits2 <- lapply(seq_len(nrow(quads2)), function(index) {
    fzl_ref_fit(y.data, quads2[index, ], identification, tol, max_iter)
  })

  block.main <- block.trace3 <- matrix(0, nrow(pairs1), nrow(pairs2))
  weight1 <- n1 / (n1 + n2)
  weight2 <- n2 / (n1 + n2)
  for (a in seq_len(nrow(pairs1))) {
    i <- pairs1[a, 1L]
    j <- pairs1[a, 2L]
    for (b in seq_len(nrow(pairs2))) {
      s <- pairs2[b, 1L]
      tt <- pairs2[b, 2L]
      pooled <- weight1 * pair.fits1[[a]]$diagonal +
        weight2 * pair.fits2[[b]]$diagonal
      uis <- fzl_ref_direction(x.data[i, ] - y.data[s, ], pooled)
      ujt <- fzl_ref_direction(x.data[j, ] - y.data[tt, ], pooled)
      uit <- fzl_ref_direction(x.data[i, ] - y.data[tt, ], pooled)
      ujs <- fzl_ref_direction(x.data[j, ] - y.data[s, ], pooled)
      block.main[a, b] <- 2 * (sum(uis * ujt) + sum(uit * ujs))
      ux <- fzl_ref_direction(x.data[i, ] - x.data[j, ], pooled)
      uy <- fzl_ref_direction(y.data[s, ] - y.data[tt, ], pooled)
      block.trace3[a, b] <- 4 * sum(ux * uy)^2
    }
  }
  main.sum <- sum(block.main)
  T.n <- main.sum / (n1 * (n1 - 1) * n2 * (n2 - 1))
  trace1 <- fzl_ref_within_trace(x.data, quads1, quad.fits1)
  trace2 <- fzl_ref_within_trace(y.data, quads2, quad.fits2)
  trace3.sum <- sum(block.trace3)
  trace3 <- p^2 * trace3.sum / (n1^2 * n2^2)
  variance.term1 <- trace1$estimate / (2 * n1 * (n1 - 1) * p^2)
  variance.term2 <- trace2$estimate / (2 * n2 * (n2 - 1) * p^2)
  variance.term3 <- trace3 / (n1 * n2 * p^2)
  sigma2 <- variance.term1 + variance.term2 + variance.term3

  list(
    prepared = prepared,
    full1 = full1,
    full2 = full2,
    pair.fits1 = pair.fits1,
    pair.fits2 = pair.fits2,
    quad.fits1 = quad.fits1,
    quad.fits2 = quad.fits2,
    pairs1 = pairs1,
    pairs2 = pairs2,
    quads1 = quads1,
    quads2 = quads2,
    main.blocks = block.main,
    trace3.blocks = block.trace3,
    main.sum = main.sum,
    T.n = T.n,
    trace1 = trace1,
    trace2 = trace2,
    trace3.sum = trace3.sum,
    trace3 = trace3,
    variance.term1 = variance.term1,
    variance.term2 = variance.term2,
    variance.term3 = variance.term3,
    sigma2 = sigma2,
    statistic = T.n / sqrt(sigma2)
  )
}


test_that("Feng-Zhang-Liu matches both literal pure-R identifications", {
  data <- fzl_fixture()
  for (identification in c("geometric", "paper_trace")) {
    reference <- fzl_reference(
      data$x, data$y, identification = identification, tol = 1e-7
    )
    result <- feng_zhang_liu_spatial_rank_test(
      data$x, data$y, tol = 1e-7,
      scale_identification = identification
    )
    components <- result$components

    expect_equal(unname(result$raw.statistic), reference$T.n,
                 tolerance = 3e-10)
    expect_equal(unname(result$statistic), reference$statistic,
                 tolerance = 3e-9)
    expect_equal(components$main.ordered.sum, reference$main.sum,
                 tolerance = 3e-9)
    expect_equal(unname(components$main.pair.block.contributions),
                 unname(reference$main.blocks), tolerance = 3e-9)
    expect_equal(components$trace.R1.squared.hat,
                 reference$trace1$estimate, tolerance = 3e-9)
    expect_equal(components$trace.R2.squared.hat,
                 reference$trace2$estimate, tolerance = 3e-9)
    expect_equal(components$trace.R1.R2.hat, reference$trace3,
                 tolerance = 3e-9)
    expect_equal(components$trace1.ordered.sum,
                 reference$trace1$ordered.sum, tolerance = 3e-9)
    expect_equal(components$trace2.ordered.sum,
                 reference$trace2$ordered.sum, tolerance = 3e-9)
    expect_equal(components$trace3.ordered.sum, reference$trace3.sum,
                 tolerance = 3e-9)
    expect_equal(unname(components$trace1.permutation.contributions),
                 unname(reference$trace1$contributions), tolerance = 3e-9)
    expect_equal(unname(components$trace2.permutation.contributions),
                 unname(reference$trace2$contributions), tolerance = 3e-9)
    expect_equal(unname(components$trace3.pair.block.contributions),
                 unname(reference$trace3.blocks), tolerance = 3e-9)
    expect_equal(components$variance.term1, reference$variance.term1,
                 tolerance = 3e-11)
    expect_equal(components$variance.term2, reference$variance.term2,
                 tolerance = 3e-11)
    expect_equal(components$variance.term3, reference$variance.term3,
                 tolerance = 3e-11)
    expect_equal(components$sigma2.hat, reference$sigma2,
                 tolerance = 3e-11)
    expect_equal(
      unname(components$full.fit$group1$scale.diagonal.standardized),
      reference$full1$diagonal, tolerance = 3e-8
    )
    expect_equal(
      unname(components$full.fit$group2$scale.diagonal.standardized),
      reference$full2$diagonal, tolerance = 3e-8
    )
    expect_equal(
      unname(
        components$leave.two.out.fit$group1$scale.diagonal.standardized
      ),
      do.call(rbind, lapply(reference$pair.fits1, `[[`, "diagonal")),
      tolerance = 3e-8
    )
    expect_equal(
      unname(
        components$leave.four.out.fit$group2$scale.diagonal.standardized
      ),
      do.call(rbind, lapply(reference$quad.fits2, `[[`, "diagonal")),
      tolerance = 3e-8
    )
    expect_identical(result$diagnostics$scale.identification,
                     identification)
  }
})


test_that("ordered factors and the corrected feasible variance are locked", {
  data <- fzl_fixture()
  result <- feng_zhang_liu_spatial_rank_test(data$x, data$y)
  components <- result$components
  n1 <- nrow(data$x)
  n2 <- nrow(data$y)
  p <- ncol(data$x)

  expect_equal(components$main.ordered.denominator,
               n1 * (n1 - 1) * n2 * (n2 - 1))
  expect_equal(components$trace1.ordered.denominator,
               n1 * (n1 - 1) * (n1 - 2) * (n1 - 3))
  expect_equal(components$trace2.ordered.denominator,
               n2 * (n2 - 1) * (n2 - 2) * (n2 - 3))
  expect_equal(components$trace3.ordered.denominator, n1^2 * n2^2)
  expect_equal(components$T.n,
               sum(components$main.pair.block.contributions) /
                 components$main.ordered.denominator,
               tolerance = 1e-14)
  expect_equal(components$trace.R1.squared.hat,
               2 * p^2 * components$trace1.ordered.sum /
                 components$trace1.ordered.denominator,
               tolerance = 1e-13)
  expect_equal(components$trace.R2.squared.hat,
               2 * p^2 * components$trace2.ordered.sum /
                 components$trace2.ordered.denominator,
               tolerance = 1e-13)
  expect_equal(components$trace.R1.R2.hat,
               p^2 * components$trace3.ordered.sum /
                 components$trace3.ordered.denominator,
               tolerance = 1e-13)
  expect_equal(components$variance.term2,
               components$trace.R2.squared.hat /
                 (2 * n2 * (n2 - 1) * p^2), tolerance = 1e-15)
  expect_false(isTRUE(all.equal(
    components$variance.term2,
    components$trace.R2.squared.hat / (2 * n2 * (n2 - 1)),
    tolerance = 1e-12
  )))
  expect_true(result$diagnostics$second.variance.term.p.squared.restored)
})


test_that("geometric identification has the claimed exact invariances", {
  data <- fzl_fixture()
  base <- feng_zhang_liu_spatial_rank_test(data$x, data$y)
  multiplier <- c(-7, 0.2, 13)
  shift <- c(10, -4, 0.7)
  scaled <- feng_zhang_liu_spatial_rank_test(
    sweep(sweep(data$x, 2L, multiplier, "*"), 2L, shift, "+"),
    sweep(sweep(data$y, 2L, multiplier, "*"), 2L, shift, "+")
  )
  permuted <- feng_zhang_liu_spatial_rank_test(
    data$x[c(4, 1, 6, 2, 5, 3), , drop = FALSE],
    data$y[c(2, 6, 3, 1, 5, 4), , drop = FALSE]
  )
  swapped <- feng_zhang_liu_spatial_rank_test(data$y, data$x)

  expect_equal(unname(scaled$raw.statistic), unname(base$raw.statistic),
               tolerance = 2e-10)
  expect_equal(unname(scaled$statistic), unname(base$statistic),
               tolerance = 2e-9)
  expect_equal(unname(scaled$variance), unname(base$variance),
               tolerance = 2e-10)
  expect_equal(unname(permuted$raw.statistic), unname(base$raw.statistic),
               tolerance = 2e-10)
  expect_equal(unname(permuted$variance), unname(base$variance),
               tolerance = 2e-10)
  expect_equal(unname(swapped$raw.statistic), unname(base$raw.statistic),
               tolerance = 2e-10)
  expect_equal(unname(swapped$variance), unname(base$variance),
               tolerance = 2e-10)
  expect_match(base$diagnostics$scale.identification.note,
               "coordinatewise-scale")
})


test_that("literal paper trace normalization has a finite-sample counterexample", {
  data <- fzl_fixture()
  multiplier <- c(-7, 0.2, 13)
  base <- feng_zhang_liu_spatial_rank_test(
    data$x, data$y, scale_identification = "paper_trace"
  )
  scaled <- feng_zhang_liu_spatial_rank_test(
    sweep(data$x, 2L, multiplier, "*"),
    sweep(data$y, 2L, multiplier, "*"),
    scale_identification = "paper_trace"
  )

  expect_gt(abs(unname(base$raw.statistic - scaled$raw.statistic)), 1e-6)
  expect_gt(abs(unname(base$statistic - scaled$statistic)), 1e-6)
  expect_identical(base$diagnostics$scale.identification, "paper_trace")
  expect_match(base$diagnostics$scale.identification.note,
               "violate.*coordinatewise-scale")
})


test_that("extreme finite coordinate units are handled in geometric mode", {
  data <- fzl_fixture()
  base <- feng_zhang_liu_spatial_rank_test(data$x, data$y)
  multiplier <- c(1e300, -1e-300, 5e150)
  extreme <- feng_zhang_liu_spatial_rank_test(
    sweep(data$x, 2L, multiplier, "*"),
    sweep(data$y, 2L, multiplier, "*")
  )

  expect_true(is.finite(unname(extreme$statistic)))
  expect_true(is.finite(unname(extreme$variance)))
  expect_equal(unname(extreme$raw.statistic), unname(base$raw.statistic),
               tolerance = 3e-10)
  expect_equal(unname(extreme$statistic), unname(base$statistic),
               tolerance = 3e-9)
  expect_identical(extreme$diagnostics$variance.repair, "none")
})


test_that("U(0)=0 is explicit and cross-sample ties remain callable", {
  data <- fzl_fixture()
  data$y[1L, ] <- data$x[1L, ]
  result <- feng_zhang_liu_spatial_rank_test(data$x, data$y)

  expect_true(is.finite(unname(result$statistic)))
  expect_gt(result$diagnostics$zero.directions$main.evaluations, 0)
  expect_match(result$diagnostics$sign.at.zero, "U\\(0\\) = 0")
  expect_identical(result$diagnostics$regularization, "none")
})


test_that("the htest and applicability contracts are explicit", {
  data <- fzl_fixture()
  result <- feng_zhang_liu_spatial_rank_test(data$x, data$y)

  expect_s3_class(result, "hd_location_test")
  expect_s3_class(result, "htest")
  expect_identical(result$alternative, "two.sided")
  expect_identical(result$null.distribution$tail, "upper")
  expect_equal(result$p.value,
               stats::pnorm(unname(result$statistic), lower.tail = FALSE))
  expect_true(result$diagnostics$applicability$formal.common.scatter.required)
  expect_false(
    result$diagnostics$applicability$unequal.scatter.theory.available.in.paper
  )
  expect_match(result$diagnostics$calibration, "local-alternative.*not used")
  expect_identical(result$diagnostics$variance.repair, "none")
})


test_that("invalid and degenerate inputs fail without silent repair", {
  data <- fzl_fixture()
  expect_error(
    feng_zhang_liu_spatial_rank_test(data$x[-1L, ], data$y),
    "at least 6 row"
  )
  expect_error(
    feng_zhang_liu_spatial_rank_test(data$x, data$y[, -1L]),
    "same number of columns"
  )
  expect_error(
    feng_zhang_liu_spatial_rank_test(
      data$x, data$y, scale_identification = "unknown"
    ),
    "arg"
  )
  x.constant <- data$x
  x.constant[, 1L] <- 1
  expect_error(
    feng_zhang_liu_spatial_rank_test(x.constant, data$y),
    "marginal variance is non-positive"
  )
  expect_error(
    feng_zhang_liu_spatial_rank_test(
      data$x, data$y, tol = 1e-14, max_iter = 1L
    ),
    "failed to converge"
  )
  expect_error(
    feng_zhang_liu_spatial_rank_test(data$x, data$y, tol = 0),
    "finite positive"
  )
})
