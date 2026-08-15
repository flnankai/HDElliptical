tinst_rows <- function(x, omitted = NULL) {
  if (is.null(omitted)) x else x[-omitted, , drop = FALSE]
}


tinst_joint_reference <- function(x, omitted = NULL, tol = 1e-8,
                                  max_iter = 500L, zero_tol = 0) {
  z <- tinst_rows(x, omitted)
  p <- ncol(z)
  theta <- colMeans(z)
  diagonal <- apply(z, 2L, stats::var)
  if (any(!is.finite(diagonal)) || any(diagonal <= 0)) {
    stop("non-positive reference variance")
  }
  initial.diagonal <- diagonal
  converged <- FALSE
  iterations <- 0L
  location.change <- Inf
  diagonal.change <- Inf

  for (update in 0:(max_iter - 1L)) {
    epsilon <- sweep(sweep(z, 2L, theta, "-"), 2L, sqrt(diagonal), "/")
    radius <- sqrt(rowSums(epsilon^2))
    if (any(radius <= zero_tol)) stop("zero reference residual")
    direction <- epsilon / radius
    location.residual <- sqrt(sum(colMeans(direction)^2))
    diagonal.residual <- max(abs(p * colMeans(direction^2) - 1))
    minimum.radius <- min(radius)
    standardized.step <- colSums(direction) * minimum.radius /
      sum(minimum.radius / radius)
    next.theta <- theta + sqrt(diagonal) * standardized.step
    next.diagonal <- p * diagonal * colMeans(direction^2)
    location.change <- max(abs(
      (next.theta - theta) / sqrt(initial.diagonal)
    ))
    diagonal.change <- max(abs(log(next.diagonal / diagonal)))
    theta <- next.theta
    diagonal <- next.diagonal
    iterations <- update + 1L
    if (max(location.change, diagonal.change) <= tol) {
      converged <- TRUE
      break
    }
  }
  list(
    theta = theta,
    diagonal = diagonal,
    iterations = iterations,
    converged = converged,
    location.change = location.change,
    diagonal.change = diagonal.change,
    location.residual = location.residual,
    diagonal.residual = diagonal.residual
  )
}


tinst_weighted_reference <- function(x, diagonal, omitted = NULL,
                                     tol = 1e-8, max_iter = 500L,
                                     zero_tol = 0) {
  z <- tinst_rows(x, omitted)
  theta <- colMeans(z)
  converged <- FALSE
  iterations <- 0L
  location.change <- Inf

  for (update in 0:(max_iter - 1L)) {
    epsilon <- sweep(sweep(z, 2L, theta, "-"),
                     2L, sqrt(diagonal), "/")
    radius <- sqrt(rowSums(epsilon^2))
    if (any(radius <= zero_tol)) stop("zero reference residual")
    direction <- epsilon / radius
    minimum.radius <- min(radius)
    inverse.radius.weights <- minimum.radius / radius
    equation.numerator <- colSums(direction * inverse.radius.weights)
    location.residual <- sqrt(sum(
      (equation.numerator / sum(inverse.radius.weights))^2
    ))
    standardized.step <- minimum.radius * equation.numerator /
      sum(inverse.radius.weights^2)
    next.theta <- theta + sqrt(diagonal) * standardized.step
    location.change <- max(abs((next.theta - theta) / sqrt(diagonal)))
    theta <- next.theta
    iterations <- update + 1L
    if (location.change <= tol) {
      converged <- TRUE
      break
    }
  }
  list(
    theta = theta,
    iterations = iterations,
    converged = converged,
    location.change = location.change,
    location.residual = location.residual
  )
}


tinst_reference <- function(x, y, tol = 1e-8, max_iter = 500L,
                            zero_tol = 0) {
  n1 <- nrow(x)
  n2 <- nrow(y)
  p <- ncol(x)
  full1 <- tinst_joint_reference(
    x, tol = tol, max_iter = max_iter, zero_tol = zero_tol
  )
  full2 <- tinst_joint_reference(
    y, tol = tol, max_iter = max_iter, zero_tol = zero_tol
  )

  joint1 <- weighted1 <- vector("list", n1)
  joint2 <- weighted2 <- vector("list", n2)
  direction1 <- matrix(NA_real_, n1, p)
  direction2 <- matrix(NA_real_, n2, p)
  radius1 <- numeric(n1)
  radius2 <- numeric(n2)
  for (i in seq_len(n1)) {
    joint1[[i]] <- tinst_joint_reference(
      x, i, tol, max_iter, zero_tol
    )
    weighted1[[i]] <- tinst_weighted_reference(
      x, joint1[[i]]$diagonal, i, tol, max_iter, zero_tol
    )
    epsilon <- (x[i, ] - weighted1[[i]]$theta) /
      sqrt(joint1[[i]]$diagonal)
    radius1[i] <- sqrt(sum(epsilon^2))
    direction1[i, ] <- epsilon / radius1[i]
  }
  for (i in seq_len(n2)) {
    joint2[[i]] <- tinst_joint_reference(
      y, i, tol, max_iter, zero_tol
    )
    weighted2[[i]] <- tinst_weighted_reference(
      y, joint2[[i]]$diagonal, i, tol, max_iter, zero_tol
    )
    epsilon <- (y[i, ] - weighted2[[i]]$theta) /
      sqrt(joint2[[i]]$diagonal)
    radius2[i] <- sqrt(sum(epsilon^2))
    direction2[i, ] <- epsilon / radius2[i]
  }

  cross.sum <- 0
  for (i in seq_len(n1)) {
    for (j in seq_len(n2)) {
      epsilon1 <- (x[i, ] - weighted2[[j]]$theta) /
        sqrt(joint1[[i]]$diagonal)
      epsilon2 <- (y[j, ] - weighted1[[i]]$theta) /
        sqrt(joint2[[j]]$diagonal)
      radius.cross1 <- sqrt(sum(epsilon1^2))
      radius.cross2 <- sqrt(sum(epsilon2^2))
      value <- if (radius.cross1 == 0 || radius.cross2 == 0) {
        0
      } else {
        sum((epsilon1 / radius.cross1) * (epsilon2 / radius.cross2)) /
          (radius.cross1 * radius.cross2)
      }
      cross.sum <- cross.sum + value
    }
  }
  T <- -cross.sum / (n1 * n2)

  ratio1 <- sqrt(full1$diagonal / full2$diagonal)
  ratio2 <- sqrt(full2$diagonal / full1$diagonal)
  a1.sum <- 0
  a2.sum <- 0
  a12.sum <- 0
  for (i in seq_len(n1)) {
    for (j in seq_len(n1)) {
      if (i != j) {
        a1.sum <- a1.sum + sum(direction1[j, ] * ratio1 *
                                 direction1[i, ])^2
      }
    }
  }
  for (i in seq_len(n2)) {
    for (j in seq_len(n2)) {
      if (i != j) {
        a2.sum <- a2.sum + sum(direction2[j, ] * ratio2 *
                                 direction2[i, ])^2
      }
    }
  }
  for (i in seq_len(n1)) {
    for (j in seq_len(n2)) {
      a12.sum <- a12.sum + sum(direction1[i, ] * direction2[j, ])^2
    }
  }
  A1 <- p^2 * a1.sum / (n1 * (n1 - 1))
  A2 <- p^2 * a2.sum / (n2 * (n2 - 1))
  A12 <- p^2 * a12.sum / (n1 * n2)
  nu12 <- mean(radius1^-2)
  nu22 <- mean(radius2^-2)
  c1 <- nu12
  c2 <- nu22
  variance.term1 <- 2 * nu12^2 * c2^2 / c1^2 * A1 /
    (n1 * (n1 - 1) * p^2)
  variance.term2 <- 2 * nu22^2 * c1^2 / c2^2 * A2 /
    (n2 * (n2 - 1) * p^2)
  variance.term12 <- 4 * nu12 * nu22 * A12 / (n1 * n2 * p^2)
  variance <- variance.term1 + variance.term2 + variance.term12

  list(
    T = T,
    A1 = A1,
    A2 = A2,
    A12 = A12,
    nu12 = nu12,
    nu22 = nu22,
    c1 = c1,
    c2 = c2,
    variance.term1 = variance.term1,
    variance.term2 = variance.term2,
    variance.term12 = variance.term12,
    variance = variance,
    z = T / sqrt(variance),
    full1 = full1,
    full2 = full2,
    joint1 = joint1,
    joint2 = joint2,
    weighted1 = weighted1,
    weighted2 = weighted2,
    direction1 = direction1,
    direction2 = direction2,
    radius1 = radius1,
    radius2 = radius2,
    ratio1 = ratio1,
    ratio2 = ratio2
  )
}


test_that("tINST matches a literal n1 != n2 reference", {
  set.seed(25101)
  x <- matrix(stats::rnorm(8 * 3), 8, 3)
  y <- matrix(stats::rnorm(11 * 3, 0.2), 11, 3)
  x <- sweep(x, 2L, c(0.5, 2, 4), "*")
  y <- sweep(y, 2L, c(0.5, 2, 4), "*")
  tolerance <- 1e-7
  reference <- tinst_reference(x, y, tol = tolerance, max_iter = 2000L)
  result <- tinst_two_sample_test(
    x, y, tol = tolerance, max_iter = 2000L
  )

  expect_s3_class(result, "hd_location_test")
  expect_s3_class(result, "htest")
  expect_equal(unname(result$statistic), reference$z, tolerance = 3e-7)
  expect_equal(result$components$T.tINST, reference$T, tolerance = 2e-7)
  expect_equal(result$components$A1, reference$A1, tolerance = 2e-7)
  expect_equal(result$components$A2, reference$A2, tolerance = 2e-7)
  expect_equal(result$components$A12, reference$A12, tolerance = 2e-7)
  expect_equal(result$components$nu12, reference$nu12, tolerance = 2e-7)
  expect_equal(result$components$nu22, reference$nu22, tolerance = 2e-7)
  expect_identical(result$components$nu12, result$components$c1)
  expect_identical(result$components$nu22, result$components$c2)
  expect_equal(result$components$variance.term1,
               reference$variance.term1, tolerance = 4e-7)
  expect_equal(result$components$variance.term2,
               reference$variance.term2, tolerance = 4e-7)
  expect_equal(result$components$variance.term12,
               reference$variance.term12, tolerance = 4e-7)
  expect_equal(result$variance[[1L]], reference$variance, tolerance = 4e-7)
  expect_equal(
    result$variance[[1L]],
    result$components$variance.term1 +
      result$components$variance.term2 +
      result$components$variance.term12,
    tolerance = 2e-15
  )
  expect_equal(
    unname(result$statistic),
    result$components$T.tINST / sqrt(result$variance[[1L]]),
    tolerance = 2e-15
  )
  expect_equal(result$p.value,
               stats::pnorm(reference$z, lower.tail = FALSE),
               tolerance = 3e-8)
  expect_true(result$diagnostics$iteration.stable)
  expect_equal(result$n, c(x = 8, y = 11))
})


test_that("tINST trace factors use ordered pairs and the corrected n2 limit", {
  set.seed(25102)
  x <- matrix(stats::rnorm(7 * 3), 7, 3)
  y <- matrix(stats::rnorm(10 * 3, -0.15), 10, 3)
  reference <- tinst_reference(x, y, tol = 1e-7, max_iter = 2000L)
  result <- tinst_two_sample_test(
    x, y, tol = 1e-7, max_iter = 2000L
  )

  unordered1 <- 0
  for (i in seq_len(nrow(x) - 1L)) {
    for (j in (i + 1L):nrow(x)) {
      unordered1 <- unordered1 + sum(
        reference$direction1[j, ] * reference$ratio1 *
          reference$direction1[i, ]
      )^2
    }
  }
  expect_equal(
    result$components$A1,
    2 * ncol(x)^2 * unordered1 / (nrow(x) * (nrow(x) - 1)),
    tolerance = 2e-7
  )

  wrong.a2.sum <- 0
  for (i in seq_len(nrow(x))) {
    for (j in seq_len(nrow(y))) {
      if (i != j) {
        wrong.a2.sum <- wrong.a2.sum + sum(
          reference$direction2[j, ] * reference$ratio2 *
            reference$direction2[i, ]
        )^2
      }
    }
  }
  wrong.A2 <- ncol(x)^2 * wrong.a2.sum /
    (nrow(y) * (nrow(y) - 1))
  expect_equal(result$components$A2, reference$A2, tolerance = 2e-7)
  expect_gt(abs(result$components$A2 - wrong.A2), 1e-4)
  expect_identical(result$diagnostics$trace.pair.convention,
                   "ordered i != j pairs")
  expect_match(result$diagnostics$corrected.source.indices[[3L]], "n2")
})


test_that("one-update diagnostics regress the n_k and D_k source typos", {
  set.seed(25103)
  x <- matrix(stats::rnorm(8 * 3), 8, 3)
  y <- matrix(stats::rnorm(12 * 3, 0.1), 12, 3)
  expect_warning(
    result <- tinst_two_sample_test(
      x, y, tol = 1e-14, max_iter = 1L, strict = FALSE
    ),
    "iteration-stable"
  )
  anchor <- result$diagnostics$internal.anchor
  scale <- result$diagnostics$internal.column.scale
  xs <- sweep(sweep(x, 2L, anchor, "-"), 2L, scale, "/")
  ys <- sweep(sweep(y, 2L, anchor, "-"), 2L, scale, "/")
  first1 <- tinst_joint_reference(
    xs, omitted = 1L, tol = 1e-14, max_iter = 1L
  )
  first2 <- tinst_joint_reference(
    ys, omitted = 1L, tol = 1e-14, max_iter = 1L
  )

  expect_equal(result$components$joint.location1.scaled[1, ],
               first1$theta, tolerance = 2e-13)
  expect_equal(result$components$leaveout.diagonal1.scaled[1, ],
               first1$diagonal, tolerance = 2e-13)
  expect_equal(result$components$joint.location2.scaled[1, ],
               first2$theta, tolerance = 2e-13)
  expect_equal(result$components$leaveout.diagonal2.scaled[1, ],
               first2$diagonal, tolerance = 2e-13)
  expect_equal(length(result$diagnostics$leaveout.joint$group1$iterations), 8)
  expect_equal(length(result$diagnostics$leaveout.joint$group2$iterations), 12)
})


test_that("tINST is invariant to group order, shifts, and diagonal units", {
  set.seed(25104)
  x <- matrix(stats::rnorm(10 * 4), 10, 4)
  y <- matrix(stats::rnorm(13 * 4, 0.15), 13, 4)
  base <- tinst_two_sample_test(x, y, tol = 1e-7, max_iter = 2000L)
  swapped <- tinst_two_sample_test(y, x, tol = 1e-7, max_iter = 2000L)

  expect_equal(unname(base$statistic), unname(swapped$statistic),
               tolerance = 2e-10)
  expect_equal(base$components$T.tINST, swapped$components$T.tINST,
               tolerance = 2e-10)
  expect_equal(base$variance, swapped$variance, tolerance = 2e-10)
  expect_equal(base$components$A1, swapped$components$A2,
               tolerance = 2e-10)
  expect_equal(base$components$A2, swapped$components$A1,
               tolerance = 2e-10)
  expect_equal(base$components$nu12, swapped$components$nu22,
               tolerance = 2e-10)
  expect_equal(base$components$nu22, swapped$components$nu12,
               tolerance = 2e-10)

  multipliers <- c(-3, 0.2, 8, -0.75)
  translation <- c(1e8, -2e8, 3e8, -4e8)
  scaled <- tinst_two_sample_test(
    sweep(x, 2L, multipliers, "*"),
    sweep(y, 2L, multipliers, "*"),
    tol = 1e-7, max_iter = 2000L
  )
  shifted <- tinst_two_sample_test(
    sweep(x, 2L, translation, "+"),
    sweep(y, 2L, translation, "+"),
    tol = 1e-7, max_iter = 2000L
  )
  expect_equal(unname(base$statistic), unname(scaled$statistic),
               tolerance = 2e-10)
  expect_equal(base$variance, scaled$variance, tolerance = 2e-10)
  expect_equal(unname(base$statistic), unname(shifted$statistic),
               tolerance = 2e-6)
  expect_equal(base$variance, shifted$variance, tolerance = 2e-6)

  transformed.x <- sweep(
    sweep(x, 2L, multipliers, "*"), 2L, translation, "+"
  )
  transformed.y <- sweep(
    sweep(y, 2L, multipliers, "*"), 2L, translation, "+"
  )
  transformed <- tinst_two_sample_test(
    transformed.x, transformed.y, tol = 1e-7, max_iter = 2000L
  )
  expect_equal(unname(base$statistic), unname(transformed$statistic),
               tolerance = 2e-6)
  expect_equal(base$components$T.tINST,
               transformed$components$T.tINST, tolerance = 2e-6)
  expect_equal(base$variance, transformed$variance, tolerance = 2e-6)
  expect_equal(base$p.value, transformed$p.value, tolerance = 2e-6)
})


test_that("tINST returns its upper-tail rejection and iteration diagnostics", {
  set.seed(25105)
  x <- matrix(stats::rnorm(9 * 3), 9, 3)
  y <- matrix(stats::rnorm(12 * 3, 0.4), 12, 3)
  result <- tinst_two_sample_test(
    x, y, alpha = 0.1, tol = 1e-7, max_iter = 2000L
  )

  expect_equal(result$diagnostics$rejection$critical.value,
               stats::qnorm(0.9), tolerance = 1e-15)
  expect_identical(
    result$diagnostics$rejection$reject,
    isTRUE(unname(result$statistic) > stats::qnorm(0.9))
  )
  expect_identical(result$null.distribution$tail, "upper")
  expect_true(all(
    result$diagnostics$leaveout.joint$group1$iteration.stable
  ))
  expect_true(all(
    result$diagnostics$leaveout.joint$group2$iteration.stable
  ))
  expect_true(all(
    result$diagnostics$leaveout.weighted.location$group1$iteration.stable
  ))
  expect_true(all(
    result$diagnostics$leaveout.weighted.location$group2$iteration.stable
  ))
  expect_true(all(!is.na(
    result$diagnostics$leaveout.weighted.location$group1$
      score.residual
  )))
  expect_true(all(
    result$diagnostics$leaveout.weighted.location$group1$relative.update <=
      1e-7
  ))
  expect_gt(
    result$diagnostics$leaveout.weighted.location$group1$
      worst.score.residual,
    0.1
  )
  expect_identical(
    result$diagnostics$convergence.basis,
    "relative update as specified by paper"
  )
  expect_equal(
    result$diagnostics$leaveout.joint$group1$worst.relative.update,
    max(result$diagnostics$leaveout.joint$group1$relative.update)
  )
  expect_equal(
    result$diagnostics$leaveout.weighted.location$group2$
      smallest.residual.distance,
    min(result$diagnostics$leaveout.weighted.location$group2$
          minimum.residual.distance)
  )
  expect_equal(
    result$diagnostics$worst.leaveout$relative.update,
    max(
      result$diagnostics$leaveout.joint$group1$worst.relative.update,
      result$diagnostics$leaveout.joint$group2$worst.relative.update,
      result$diagnostics$leaveout.weighted.location$group1$
        worst.relative.update,
      result$diagnostics$leaveout.weighted.location$group2$
        worst.relative.update
    )
  )
  expect_equal(
    result$diagnostics$minimum.residual.distance,
    min(
      result$diagnostics$full.sample$group1$minimum.residual.distance,
      result$diagnostics$full.sample$group2$minimum.residual.distance,
      result$diagnostics$worst.leaveout$minimum.residual.distance
    )
  )
  expect_identical(result$diagnostics$regularization, "none")
  expect_identical(result$diagnostics$variance.repair, "none")
  expect_identical(result$diagnostics$zero.residual.perturbation, "none")
})


test_that("tINST has explicit degeneracy and nonconvergence contracts", {
  set.seed(25106)
  x <- matrix(stats::rnorm(8 * 3), 8, 3)
  y <- matrix(stats::rnorm(10 * 3), 10, 3)

  expect_error(tinst_two_sample_test(x[1:2, ], y), "at least 3")
  expect_error(tinst_two_sample_test(x, y[, 1:2]), "same number")
  expect_error(tinst_two_sample_test(x, y, alpha = 0), "strictly between")
  expect_error(tinst_two_sample_test(x, y, tol = 0), "positive")
  expect_error(tinst_two_sample_test(x, y, max_iter = 1.5),
               "positive integer")
  expect_error(tinst_two_sample_test(x, y, zero_tol = -1), "non-negative")
  expect_error(tinst_two_sample_test(x, y, strict = NA), "TRUE.*FALSE")

  constant.x <- cbind(stats::rnorm(8), 1, stats::rnorm(8))
  constant.y <- cbind(stats::rnorm(10), 1, stats::rnorm(10))
  expect_error(
    tinst_two_sample_test(constant.x, constant.y),
    "strictly positive; no ridge"
  )

  coincident.mean <- matrix(c(-1, 0, 1), ncol = 1)
  ordinary.y <- matrix(c(-2, -0.4, 0.7, 2), ncol = 1)
  expect_error(
    tinst_two_sample_test(coincident.mean, ordinary.y),
    "zero standardised residual.*no perturbation"
  )

  expect_error(
    tinst_two_sample_test(x, y, tol = 1e-14, max_iter = 1L),
    "iteration-stable"
  )
  expect_warning(
    loose <- tinst_two_sample_test(
      x, y, tol = 1e-14, max_iter = 1L, strict = FALSE
    ),
    "iteration-stable"
  )
  expect_false(loose$diagnostics$iteration.stable)
  expect_gt(loose$diagnostics$stability.failures, 0)
  expect_true(all(
    loose$diagnostics$leaveout.joint$group1$iterations <= 1L
  ))
  expect_true(all(
    loose$diagnostics$leaveout.weighted.location$group2$iterations <= 1L
  ))
})
