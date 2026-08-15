ch3g_highdim_fixture <- function() {
  matrix(c(
    -2.0,  0.3,  1.0,
    -1.4,  1.6, -0.8,
    -0.8, -1.3,  0.2,
    -0.1,  0.9,  1.8,
     0.5, -0.7, -1.5,
     1.0,  1.2,  0.6,
     1.7, -1.0,  1.4,
     2.3,  0.1, -0.9,
     2.8,  1.5,  0.3,
     3.2, -1.6, -0.2
  ), ncol = 3L, byrow = TRUE)
}


ch3g_trace_u_reference <- function(x) {
  n <- nrow(x)
  gram <- tcrossprod(x)
  y1 <- mean(diag(gram))
  y2.total <- y3.total <- y4.total <- y5.total <- 0
  for (i in seq_len(n)) for (j in seq_len(n)) {
    if (i == j) next
    y2.total <- y2.total + gram[i, j]^2
    y3.total <- y3.total + gram[i, j]
    for (k in seq_len(n)) {
      if (k == i || k == j) next
      y4.total <- y4.total + gram[i, j] * gram[j, k]
      for (l in seq_len(n)) {
        if (l == i || l == j || l == k) next
        y5.total <- y5.total + gram[i, j] * gram[k, l]
      }
    }
  }
  p2 <- n * (n - 1)
  p3 <- p2 * (n - 2)
  p4 <- p3 * (n - 3)
  values <- c(
    Y1 = y1,
    Y2 = y2.total / p2,
    Y3 = y3.total / p2,
    Y4 = y4.total / p3,
    Y5 = y5.total / p4
  )
  c(values, T1 = values[["Y1"]] - values[["Y3"]],
    T2 = values[["Y2"]] - 2 * values[["Y4"]] + values[["Y5"]])
}


ch3g_cross_u_reference <- function(x, y) {
  n1 <- nrow(x)
  n2 <- nrow(y)
  gram <- x %*% t(y)
  first <- sum(gram^2) / (n1 * n2)
  second <- third <- fourth <- 0
  for (i in seq_len(n1)) for (k in seq_len(n1)) {
    if (i == k) next
    for (j in seq_len(n2)) {
      second <- second + gram[i, j] * gram[k, j]
    }
  }
  for (i in seq_len(n1)) for (j in seq_len(n2)) {
    for (l in seq_len(n2)) {
      if (j == l) next
      third <- third + gram[i, j] * gram[i, l]
    }
  }
  for (i in seq_len(n1)) for (k in seq_len(n1)) {
    if (i == k) next
    for (j in seq_len(n2)) for (l in seq_len(n2)) {
      if (j == l) next
      fourth <- fourth + gram[i, j] * gram[k, l]
    }
  }
  first - second / (n1 * n2 * (n1 - 1)) -
    third / (n1 * n2 * (n2 - 1)) +
    fourth / (n1 * n2 * (n1 - 1) * (n2 - 1))
}


test_that("Wang-Yao corrected LRT matches both mean conventions", {
  x <- ch3g_highdim_fixture()
  beta <- 0.45
  for (center in c(FALSE, TRUE)) {
    result <- wang_yao_corrected_lrt(x, center = center, beta = beta)
    n <- nrow(x)
    p <- ncol(x)
    m <- n - as.integer(center)
    residual <- if (center) sweep(x, 2L, colMeans(x), "-") else x
    S <- crossprod(residual) / m
    values <- eigen(S, symmetric = TRUE, only.values = TRUE)$values
    L <- -sum(log(values)) + p * log(mean(values))
    y <- p / m
    corrected <- L + (p - m) * log1p(-y) - p
    null.mean <- -log1p(-y) / 2 + beta * y / 2
    null.variance <- -2 * log1p(-y) - 2 * y
    z <- (corrected - null.mean) / sqrt(null.variance)

    expect_equal(result$components$L, L, tolerance = 3e-13)
    expect_equal(result$components$corrected.L, corrected,
                 tolerance = 3e-13)
    expect_equal(result$components$null.mean, null.mean,
                 tolerance = 2e-15)
    expect_equal(result$components$null.variance, null.variance,
                 tolerance = 2e-15)
    expect_equal(unname(result$statistic), z, tolerance = 3e-13)
    expect_equal(result$p.value, pnorm(z, lower.tail = FALSE),
                 tolerance = 2e-15)
    expect_equal(result$diagnostics$effective.covariance.df, m)
  }

  base <- wang_yao_corrected_lrt(x, beta = 0)
  large <- wang_yao_corrected_lrt(1e150 * x, beta = 0)
  Q <- qr.Q(qr(matrix(c(1, 2, 0, -1, 1, 3, 2, 0, 1), 3L)))
  rotated <- wang_yao_corrected_lrt(x %*% Q, beta = 0)
  expect_equal(unname(large$statistic), unname(base$statistic),
               tolerance = 5e-12)
  expect_equal(unname(rotated$statistic), unname(base$statistic),
               tolerance = 5e-12)
})


test_that("Wang-Yao corrected John statistic locks the n versus n-1 shift", {
  x <- ch3g_highdim_fixture()
  beta <- -0.2
  for (center in c(FALSE, TRUE)) {
    result <- wang_yao_corrected_john_test(x, center = center, beta = beta)
    n <- nrow(x)
    p <- ncol(x)
    m <- n - as.integer(center)
    residual <- if (center) sweep(x, 2L, colMeans(x), "-") else x
    S <- crossprod(residual) / m
    U <- p * sum(S^2) / sum(diag(S))^2 - 1
    spectral.center <- if (center) n * p / (n - 1) else p
    centered <- n * U - spectral.center
    z <- (centered - (1 + beta)) / 2
    expect_equal(result$components$U, U, tolerance = 3e-14)
    expect_equal(result$components$spectral.center, spectral.center,
                 tolerance = 0)
    expect_equal(result$components$centered.nU, centered,
                 tolerance = 3e-13)
    expect_equal(unname(result$statistic), z, tolerance = 3e-13)
  }

  estimated <- wang_yao_corrected_john_test(x, beta = NULL)
  residual <- x / max(abs(x))
  beta.hat <- mean(residual^4) / mean(residual^2)^2 - 3
  expect_equal(estimated$components$beta, beta.hat, tolerance = 3e-15)
  expect_identical(estimated$diagnostics$beta.source,
                   "scale-standardized empirical fourth moment")
  expect_true(is.finite(wang_yao_corrected_john_test(
    matrix(seq_len(60), 5L, 12L), beta = 0
  )$statistic))
})


test_that("Chen-Zhang-Zhong five U-statistics match ordered index sums", {
  x <- ch3g_highdim_fixture()[1:6, , drop = FALSE]
  scale <- max(abs(x))
  residual <- sweep(x / scale, 2L, colMeans(x / scale), "-")
  reference <- ch3g_trace_u_reference(residual)
  sphericity <- chen_zhang_zhong_covariance_test(x, "sphericity")
  for (name in paste0("Y", 1:5)) {
    expect_equal(sphericity$components[[paste0(name, ".scaled")]],
                 reference[[name]], tolerance = 3e-14)
  }
  expect_equal(sphericity$components$T1.scaled, reference[["T1"]],
               tolerance = 3e-14)
  expect_equal(sphericity$components$T2.scaled, reference[["T2"]],
               tolerance = 4e-14)
  U <- ncol(x) * reference[["T2"]] / reference[["T1"]]^2 - 1
  expect_equal(sphericity$raw.statistic[["U.n"]], U, tolerance = 4e-14)
  expect_equal(unname(sphericity$statistic), nrow(x) * U / 2,
               tolerance = 4e-13)

  identity <- chen_zhang_zhong_covariance_test(x, "identity")
  T1 <- reference[["T1"]] * scale^2
  T2 <- reference[["T2"]] * scale^4
  V <- T2 / ncol(x) - 2 * T1 / ncol(x) + 1
  expect_equal(identity$raw.statistic[["V.n"]], V, tolerance = 5e-14)
  expect_equal(unname(identity$statistic), nrow(x) * V / 2,
               tolerance = 5e-13)

  shifted <- chen_zhang_zhong_covariance_test(
    sweep(x, 2L, c(10, -7, 4), "+"), "sphericity"
  )
  large <- chen_zhang_zhong_covariance_test(1e150 * x, "sphericity")
  expect_equal(unname(shifted$statistic), unname(sphericity$statistic),
               tolerance = 2e-11)
  expect_equal(unname(large$statistic), unname(sphericity$statistic),
               tolerance = 2e-11)
})


test_that("Fisher-Sun-Gallagher trace correction matches its exact coefficients", {
  x <- ch3g_highdim_fixture()[1:8, , drop = FALSE]
  N <- nrow(x)
  n <- N - 1
  p <- ncol(x)
  z <- sweep(x, 2L, colMeans(x), "-")
  values <- svd(z, nu = 0L, nv = 0L)$d^2 / n
  traces <- c(tr1 = sum(values), tr2 = sum(values^2),
              tr3 = sum(values^3), tr4 = sum(values^4))
  b <- -4 / n
  c.star <- -(2 * n^2 + 3 * n - 6) / (n * (n^2 + n + 2))
  d <- 2 * (5 * n + 6) / (n * (n^2 + n + 2))
  e <- -(5 * n + 6) / (n^2 * (n^2 + n + 2))
  tau <- n^5 * (n^2 + n + 2) /
    ((n + 1) * (n + 2) * (n + 4) * (n + 6) *
       (n - 1) * (n - 2) * (n - 3))
  a2 <- n^2 / ((n - 1) * (n + 2) * p) *
    (traces[["tr2"]] - traces[["tr1"]]^2 / n)
  a4 <- tau / p * (traces[["tr4"]] +
    b * traces[["tr3"]] * traces[["tr1"]] +
    c.star * traces[["tr2"]]^2 +
    d * traces[["tr2"]] * traces[["tr1"]]^2 +
    e * traces[["tr1"]]^4)
  psi <- a4 / a2^2
  multiplier <- sqrt(n * p / (8 * (8 + 12 * p / n + (p / n)^2)))
  fit <- fisher_sun_gallagher_sphericity_test(x)

  expect_equal(fit$components$sample.trace.powers.scaled *
                 max(abs(x))^(2 * seq_len(4)),
               traces, tolerance = 5e-12)
  expect_equal(fit$components$a2.scaled * max(abs(x))^4,
               a2, tolerance = 5e-13)
  expect_equal(fit$components$a4.scaled * max(abs(x))^8,
               a4, tolerance = 2e-11)
  expect_equal(fit$components$coefficients,
               c(b = b, c.star = c.star, d = d, e = e, tau = tau),
               tolerance = 3e-15)
  expect_equal(fit$components$psi.2, psi, tolerance = 2e-12)
  expect_equal(unname(fit$statistic), multiplier * (psi - 1),
               tolerance = 2e-12)

  shifted <- fisher_sun_gallagher_sphericity_test(
    sweep(x, 2L, c(4, -8, 2), "+")
  )
  tiny <- fisher_sun_gallagher_sphericity_test(1e-120 * x)
  expect_equal(unname(shifted$statistic), unname(fit$statistic),
               tolerance = 2e-11)
  expect_equal(unname(tiny$statistic), unname(fit$statistic),
               tolerance = 2e-11)
})


test_that("Li-Chen leave-four-out and cross traces match literal sums", {
  all <- ch3g_highdim_fixture()
  x <- all[1:4, , drop = FALSE]
  y <- rbind(all[5:10, , drop = FALSE], c(3.8, 0.7, -1.1))
  A1 <- ch3g_trace_u_reference(x)[["T2"]]
  A2 <- ch3g_trace_u_reference(y)[["T2"]]
  cross <- ch3g_cross_u_reference(x, y)
  T.LC <- A1 + A2 - 2 * cross
  se <- 2 * A1 / nrow(y) + 2 * A2 / nrow(x)
  fit <- li_chen_covariance_test(x, y)

  expect_equal(fit$components$A1, A1, tolerance = 2e-12)
  expect_equal(fit$components$A2, A2, tolerance = 2e-12)
  expect_equal(fit$components$C, cross, tolerance = 2e-12)
  expect_equal(fit$components$T.LC, T.LC, tolerance = 3e-12)
  expect_equal(fit$components$estimated.null.standard.error, se,
               tolerance = 3e-12)
  expect_equal(unname(fit$statistic), T.LC / se, tolerance = 3e-12)
  expect_equal(fit$p.value, pnorm(T.LC / se, lower.tail = FALSE),
               tolerance = 2e-15)

  swapped <- li_chen_covariance_test(y, x)
  expect_equal(unname(swapped$statistic), unname(fit$statistic),
               tolerance = 3e-12)
  shifted <- li_chen_covariance_test(
    sweep(x, 2L, c(20, -7, 3), "+"),
    sweep(y, 2L, c(-11, 9, 5), "+")
  )
  expect_equal(unname(shifted$statistic), unname(fit$statistic),
               tolerance = 2e-10)
  large <- li_chen_covariance_test(1e70 * x, 1e70 * y)
  expect_equal(unname(large$statistic), unname(fit$statistic),
               tolerance = 3e-11)
})


test_that("high-dimensional covariance tests reject degenerate contracts", {
  x <- ch3g_highdim_fixture()
  expect_error(wang_yao_corrected_lrt(matrix(1:24, 4L, 6L), beta = 0),
               "p / effective.df")
  expect_error(wang_yao_corrected_lrt(x, beta = Inf), "finite")
  expect_error(wang_yao_corrected_john_test(matrix(0, 5L, 8L), beta = 0),
               "nonconstant")
  expect_error(chen_zhang_zhong_covariance_test(x[1:3, ]),
               "at least 4")
  expect_error(chen_zhang_zhong_covariance_test(matrix(1, 5L, 2L)),
               "degenerate")
  expect_error(fisher_sun_gallagher_sphericity_test(x[1:4, ]),
               "at least 5")
  expect_error(fisher_sun_gallagher_sphericity_test(matrix(1:6, ncol = 1L)),
               "two variables")
  expect_error(li_chen_covariance_test(x[1:3, ], x[4:8, ]),
               "at least 4")
  expect_error(li_chen_covariance_test(x[1:5, ], x[6:10, 1:2]),
               "same number")
  expect_error(li_chen_covariance_test(matrix(0, 4L, 2L),
                                       matrix(0, 4L, 2L)), "nonzero")
})
