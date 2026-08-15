c7_chime_fixture <- function() {
  x <- rbind(
    c(-1.5, -0.7), c(-1.1, -0.9), c(-1.3, -0.4), c(-0.8, -0.6),
    c( 0.9,  0.5), c( 1.4,  0.8), c( 1.1,  0.9), c( 0.7,  0.4)
  )
  colnames(x) <- c("a", "b")
  rownames(x) <- paste0("s", seq_len(nrow(x)))
  list(
    x = x,
    initial = list(
      omega = 0.45,
      mu1 = c(a = -0.9, b = -0.5),
      mu2 = c(a = 0.8, b = 0.5),
      covariance = structure(
        matrix(c(1.2, 0.15, 0.15, 0.8), 2, 2),
        dimnames = list(c("a", "b"), c("a", "b"))
      )
    ),
    lambda = c(0.35, 0.28, 0.22)
  )
}


c7_if_fixture <- function() {
  z <- seq(-2.5, 2.5, length.out = 12)
  x <- cbind(
    z,
    z^2,
    z^3,
    sin(z),
    cos(z),
    c(-2, -1, 0, 1, 2, 3, -3, -2, -1, 0, 1, 2)
  )
  colnames(x) <- paste0("f", seq_len(ncol(x)))
  rownames(x) <- paste0("r", seq_len(nrow(x)))
  x
}


c7_ref_soft <- function(value, penalty) {
  sign(value) * pmax(abs(value) - penalty, 0)
}


c7_ref_kkt <- function(metric, target, coefficient, lambda) {
  gradient <- drop(metric %*% coefficient - target)
  component <- ifelse(
    coefficient != 0,
    abs(gradient + lambda * sign(coefficient)),
    pmax(abs(gradient) - lambda, 0)
  )
  max(component)
}


c7_ref_lasso <- function(metric, target, coefficient, lambda,
                         tolerance, max_iterations) {
  residual <- drop(target - metric %*% coefficient)
  kkt.scale <- 1 + max(abs(target)) + lambda
  converged <- FALSE
  relative <- Inf
  kkt <- Inf
  iterations <- 0L
  for (iteration in seq_len(max_iterations)) {
    old <- coefficient
    for (j in seq_along(coefficient)) {
      partial <- residual[[j]] + metric[j, j] * coefficient[[j]]
      updated <- c7_ref_soft(partial, lambda) / metric[j, j]
      change <- updated - coefficient[[j]]
      if (change != 0) {
        coefficient[[j]] <- updated
        residual <- residual - metric[, j] * change
      }
    }
    iterations <- iteration
    relative <- max(abs(coefficient - old)) /
      max(1, abs(old), abs(coefficient))
    kkt <- c7_ref_kkt(metric, target, coefficient, lambda)
    if (relative <= tolerance && kkt <= tolerance * kkt.scale) {
      converged <- TRUE
      break
    }
  }
  list(
    coefficient = coefficient,
    iterations = iterations,
    converged = converged,
    relative = relative,
    kkt = kkt,
    kkt.scale = kkt.scale,
    objective = drop(
      crossprod(coefficient, metric %*% coefficient) / 2 -
        crossprod(coefficient, target) + lambda * sum(abs(coefficient))
    )
  )
}


c7_ref_responsibility <- function(x, omega, mu1, mu2, beta) {
  eta <- drop((x - rep((mu1 + mu2) / 2, each = nrow(x))) %*% beta)
  stats::plogis(stats::qlogis(omega) - eta)
}


c7_ref_chime <- function(x, initial, lambda, tolerance = 1e-10,
                         max_sweeps = 20000L) {
  stages <- length(lambda)
  omega <- numeric(stages)
  mu1 <- matrix(0, stages, ncol(x))
  mu2 <- matrix(0, stages, ncol(x))
  beta <- matrix(0, stages, ncol(x))
  covariance <- vector("list", stages)
  gamma <- matrix(0, nrow(x), stages - 1L)
  fits <- vector("list", stages)
  omega[[1L]] <- initial$omega
  mu1[1L, ] <- initial$mu1
  mu2[1L, ] <- initial$mu2
  covariance[[1L]] <- initial$covariance
  fits[[1L]] <- c7_ref_lasso(
    covariance[[1L]], mu1[1L, ] - mu2[1L, ], rep(0, ncol(x)),
    lambda[[1L]], tolerance, max_sweeps
  )
  beta[1L, ] <- fits[[1L]]$coefficient
  for (stage in 2:stages) {
    gamma[, stage - 1L] <- c7_ref_responsibility(
      x, omega[[stage - 1L]], mu1[stage - 1L, ],
      mu2[stage - 1L, ], beta[stage - 1L, ]
    )
    second <- sum(gamma[, stage - 1L])
    first <- nrow(x) - second
    omega[[stage]] <- second / nrow(x)
    mu1[stage, ] <- colSums(x * (1 - gamma[, stage - 1L])) / first
    mu2[stage, ] <- colSums(x * gamma[, stage - 1L]) / second
    first.residual <- x - rep(mu1[stage, ], each = nrow(x))
    second.residual <- x - rep(mu2[stage, ], each = nrow(x))
    covariance[[stage]] <-
      crossprod(first.residual * sqrt(1 - gamma[, stage - 1L])) / nrow(x) +
      crossprod(second.residual * sqrt(gamma[, stage - 1L])) / nrow(x)
    fits[[stage]] <- c7_ref_lasso(
      covariance[[stage]], mu1[stage, ] - mu2[stage, ],
      beta[stage - 1L, ], lambda[[stage]], tolerance, max_sweeps
    )
    beta[stage, ] <- fits[[stage]]$coefficient
  }
  list(
    omega = omega, mu1 = mu1, mu2 = mu2, beta = beta,
    covariance = covariance, gamma = gamma, fits = fits
  )
}


c7_ref_ks <- function(value) {
  value <- sort(value)
  n <- length(value)
  distribution <- stats::pnorm(value)
  d.plus <- max(seq_len(n) / n - distribution)
  d.minus <- max(distribution - (seq_len(n) - 1L) / n)
  c(score = sqrt(n) * max(d.plus, d.minus),
    d.plus = d.plus, d.minus = d.minus)
}


test_that("CHIME native stages agree with a line-by-line R reference", {
  fixture <- c7_chime_fixture()
  native <- cpp_ch7_chime_fit(
    fixture$x, fixture$initial$omega, fixture$initial$mu1,
    fixture$initial$mu2, fixture$initial$covariance, fixture$lambda,
    1e-10, 20000L
  )
  reference <- c7_ref_chime(
    fixture$x, fixture$initial, fixture$lambda, 1e-10, 20000L
  )
  expect_equal(drop(native$omega), reference$omega, tolerance = 2e-9)
  expect_equal(native$mu1, reference$mu1, tolerance = 2e-9)
  expect_equal(native$mu2, reference$mu2, tolerance = 2e-9)
  expect_equal(native$beta, reference$beta, tolerance = 2e-8)
  expect_equal(
    native$responsibility_history, reference$gamma, tolerance = 2e-9
  )
  for (stage in seq_along(reference$covariance)) {
    expect_equal(
      unname(native$covariance_history[[stage]]),
      unname(reference$covariance[[stage]]),
      tolerance = 2e-9
    )
  }
  expect_true(all(native$lasso_converged == 1L))
})


test_that("CHIME reports the printed objective and exact-zero KKT equations", {
  fixture <- c7_chime_fixture()
  fit <- chime_clustering(
    fixture$x, fixture$lambda, fixture$initial,
    max_iter = 2L, lasso_tol = 1e-10
  )
  expect_s3_class(fit, "chime_fit")
  expect_true(fit$valid)
  expect_identical(names(fit$cluster), rownames(fixture$x))
  expect_identical(names(fit$score), rownames(fixture$x))
  expect_identical(names(fit$responsibility), rownames(fixture$x))
  expect_equal(length(fit$history$covariance), 3L)
  for (stage in seq_len(3L)) {
    metric <- fit$history$covariance[[stage]]
    target <- fit$history$mu1[stage, ] - fit$history$mu2[stage, ]
    coefficient <- fit$history$beta[stage, ]
    expected.kkt <- c7_ref_kkt(
      metric, target, coefficient, fixture$lambda[[stage]]
    )
    expected.objective <- drop(
      crossprod(coefficient, metric %*% coefficient) / 2 -
        crossprod(coefficient, target) +
        fixture$lambda[[stage]] * sum(abs(coefficient))
    )
    expect_equal(fit$history$lasso.kkt[[stage]], expected.kkt,
                 tolerance = 2e-10)
    expect_equal(fit$history$lasso.objective[[stage]], expected.objective,
                 tolerance = 2e-10)
    expect_lte(
      fit$history$lasso.kkt[[stage]],
      1e-10 * fit$history$lasso.kkt.scale[[stage]]
    )
  }
  expect_match(fit$diagnostics$no.repair, "No ridge")
  expect_equal(fit$diagnostics$support.size, sum(fit$support))
})


test_that("CHIME initial diagonal metric has the analytic soft-threshold fit", {
  x <- matrix(c(-3, -1, -2, 0, 2, 0, 3, 1), ncol = 2, byrow = TRUE)
  covariance <- diag(c(2, 1))
  mu1 <- c(-1, -0.5)
  mu2 <- c(1, 0.5)
  lambda <- 0.2
  core <- cpp_ch7_chime_fit(
    x, 0.5, mu1, mu2, covariance, lambda, 1e-12, 100L
  )
  expected <- c7_ref_soft(mu1 - mu2, lambda) / diag(covariance)
  expect_equal(drop(core$beta), expected, tolerance = 1e-14)
  expect_equal(drop(core$lasso_kkt), 0, tolerance = 1e-14)
  expect_equal(core$cluster, ifelse(core$score >= core$threshold, 1L, 2L))
})


test_that("CHIME stable logistic remains finite for extreme scores", {
  x <- matrix(c(-1001, -999, -998, 998, 999, 1001), ncol = 1)
  core <- cpp_ch7_chime_fit(
    x, 0.5, -1000, 1000, matrix(1, 1, 1), 0, 1e-12, 100L
  )
  expect_true(all(is.finite(core$responsibility)))
  expect_true(all(core$responsibility >= 0 & core$responsibility <= 1))
  expect_identical(drop(core$cluster), c(1L, 1L, 1L, 2L, 2L, 2L))
  expect_true(all(is.finite(core$score)))
})


test_that("CHIME is translation invariant", {
  fixture <- c7_chime_fixture()
  base <- chime_clustering(
    fixture$x, fixture$lambda, fixture$initial, max_iter = 2L,
    lasso_tol = 1e-10
  )
  shift <- c(a = 4.5, b = -2.25)
  moved.initial <- fixture$initial
  moved.initial$mu1 <- moved.initial$mu1 + shift
  moved.initial$mu2 <- moved.initial$mu2 + shift
  moved <- chime_clustering(
    sweep(fixture$x, 2L, shift, "+"), fixture$lambda, moved.initial,
    max_iter = 2L, lasso_tol = 1e-10
  )
  expect_equal(moved$history$omega, base$history$omega, tolerance = 2e-9)
  expect_equal(moved$history$beta, base$history$beta, tolerance = 2e-8)
  expect_equal(moved$responsibility, base$responsibility, tolerance = 2e-9)
  expect_identical(moved$cluster, base$cluster)
  expect_equal(moved$estimate$mu1 - shift, base$estimate$mu1,
               tolerance = 2e-9)
  expect_equal(moved$estimate$mu2 - shift, base$estimate$mu2,
               tolerance = 2e-9)
})


test_that("CHIME is equivariant to feature permutation", {
  fixture <- c7_chime_fixture()
  base <- chime_clustering(
    fixture$x, fixture$lambda, fixture$initial, max_iter = 2L,
    lasso_tol = 1e-10
  )
  permutation <- c(2L, 1L)
  changed.initial <- fixture$initial
  changed.initial$mu1 <- changed.initial$mu1[permutation]
  changed.initial$mu2 <- changed.initial$mu2[permutation]
  changed.initial$covariance <-
    changed.initial$covariance[permutation, permutation]
  changed <- chime_clustering(
    fixture$x[, permutation], fixture$lambda, changed.initial,
    max_iter = 2L, lasso_tol = 1e-10
  )
  expect_identical(changed$cluster, base$cluster)
  expect_equal(changed$responsibility, base$responsibility, tolerance = 2e-9)
  expect_equal(
    changed$estimate$beta[c("a", "b")], base$estimate$beta,
    tolerance = 2e-8
  )
  expect_equal(
    changed$estimate$covariance[c("a", "b"), c("a", "b")],
    base$estimate$covariance,
    tolerance = 2e-9
  )
})


test_that("CHIME component labels are canonicalized at omega one half", {
  fixture <- c7_chime_fixture()
  base <- chime_clustering(
    fixture$x, fixture$lambda, fixture$initial, max_iter = 2L,
    lasso_tol = 1e-10
  )
  swapped.initial <- fixture$initial
  swapped.initial$omega <- 1 - fixture$initial$omega
  swapped.initial$mu1 <- fixture$initial$mu2
  swapped.initial$mu2 <- fixture$initial$mu1
  swapped <- chime_clustering(
    fixture$x, fixture$lambda, swapped.initial, max_iter = 2L,
    lasso_tol = 1e-10
  )
  expect_true(all(base$history$omega <= 0.5))
  expect_true(all(swapped$history$omega <= 0.5))
  expect_false(base$history$label.swapped[[1L]])
  expect_true(swapped$history$label.swapped[[1L]])
  expect_equal(swapped$history$omega, base$history$omega,
               tolerance = 2e-9)
  expect_equal(swapped$history$mu1, base$history$mu1, tolerance = 2e-9)
  expect_equal(swapped$history$mu2, base$history$mu2, tolerance = 2e-9)
  expect_equal(swapped$history$beta, base$history$beta, tolerance = 2e-8)
  expect_equal(swapped$responsibility, base$responsibility,
               tolerance = 2e-9)
  expect_equal(swapped$score, base$score, tolerance = 2e-8)
  expect_equal(swapped$threshold, base$threshold, tolerance = 2e-9)
  expect_identical(swapped$cluster, base$cluster)
})


test_that("CHIME canonicalizes updates whose raw mass exceeds one half", {
  x <- matrix(c(-3, -2.8, 2, 2.1, 2.2, 2.3, 2.4, 2.5), ncol = 1)
  initial <- list(
    omega = 0.49, mu1 = -3, mu2 = 2,
    covariance = diag(1)
  )
  lambda <- c(0.1, 0.08, 0.06)
  reference <- c7_ref_chime(x, initial, lambda, tolerance = 1e-10)
  fit <- chime_clustering(x, lambda, initial, max_iter = 2L,
                          lasso_tol = 1e-10)
  raw.gamma <- reference$gamma[, 1L]
  expect_identical(fit$history$label.swapped, c(FALSE, TRUE, FALSE))
  expect_gt(mean(raw.gamma), 0.5)
  expect_equal(fit$history$responsibility[, 1L], 1 - raw.gamma,
               tolerance = 2e-9)
  expect_equal(unname(fit$history$mu1[2L, ]),
               unname(reference$mu2[2L, ]), tolerance = 2e-9)
  expect_equal(unname(fit$history$mu2[2L, ]),
               unname(reference$mu1[2L, ]), tolerance = 2e-9)
  expect_equal(
    unname(fit$history$covariance[[2L]]),
    unname(reference$covariance[[2L]]),
    tolerance = 2e-9
  )
  expect_equal(
    unname(fit$history$beta[2L, ]),
    unname(-reference$fits[[2L]]$coefficient),
    tolerance = 2e-8
  )
  expect_true(all(fit$history$omega <= 0.5))
})


test_that("CHIME requires the full path and supports only K equals two", {
  fixture <- c7_chime_fixture()
  expect_error(
    chime_clustering(
      fixture$x, fixture$lambda[-1L], fixture$initial, max_iter = 2L
    ),
    "explicit"
  )
  expect_error(
    chime_clustering(
      fixture$x, fixture$lambda, fixture$initial, K = 3L, max_iter = 2L
    ),
    "only two"
  )
  expect_error(
    chime_clustering(
      fixture$x, c(0.3, NA, 0.1), fixture$initial, max_iter = 2L
    ),
    "finite"
  )
  expect_error(
    chime_clustering(fixture$x, fixture$lambda, list(), max_iter = 2L),
    "missing"
  )
})


test_that("CHIME fails rather than repairing invalid metrics or solvers", {
  fixture <- c7_chime_fixture()
  indefinite <- fixture$initial
  indefinite$covariance <- matrix(c(1, 2, 2, 1), 2, 2)
  expect_error(
    chime_clustering(
      fixture$x, fixture$lambda, indefinite, max_iter = 2L
    ),
    "indefinite"
  )
  zero.diagonal <- fixture$initial
  zero.diagonal$covariance <- diag(c(1, 0))
  expect_error(
    chime_clustering(
      fixture$x, fixture$lambda, zero.diagonal, max_iter = 2L
    ),
    "diagonal"
  )
  expect_error(
    chime_clustering(
      fixture$x, fixture$lambda, fixture$initial, max_iter = 2L,
      lasso_tol = 1e-16, lasso_max_iter = 1L
    ),
    "did not meet"
  )
})


test_that("IF-PCA paper scores agree with a line-by-line R reference", {
  x <- c7_if_fixture()
  core <- cpp_ch7_if_scores(x, FALSE)
  expected.center <- colMeans(x)
  expected.scale <- apply(x, 2L, stats::sd)
  expected.w <- sweep(sweep(x, 2L, expected.center), 2L, expected.scale, "/")
  expected <- vapply(seq_len(ncol(x)), function(j) {
    c7_ref_ks(expected.w[, j])
  }, numeric(3))
  expect_equal(unname(drop(core$center)), unname(expected.center),
               tolerance = 1e-14)
  expect_equal(unname(drop(core$scale)), unname(expected.scale),
               tolerance = 1e-14)
  expect_equal(unname(core$standardized), unname(expected.w),
               tolerance = 1e-14)
  expect_equal(drop(core$score), expected["score", ], tolerance = 2e-14)
  expect_equal(drop(core$d_plus), expected["d.plus", ], tolerance = 2e-14)
  expect_equal(drop(core$d_minus), expected["d.minus", ], tolerance = 2e-14)
})


test_that("IF-PCA software KS convention is explicit and leaves PCA W alone", {
  x <- c7_if_fixture()
  paper <- cpp_ch7_if_scores(x, FALSE)
  software <- cpp_ch7_if_scores(x, TRUE)
  factor <- sqrt(1 - 1 / nrow(x))
  expected <- vapply(seq_len(ncol(x)), function(j) {
    c7_ref_ks(paper$standardized[, j] / factor)[["score"]]
  }, numeric(1))
  expect_equal(software$standardized, paper$standardized, tolerance = 0)
  expect_equal(drop(software$score), expected, tolerance = 2e-14)
  expect_equal(software$software_ks_scale, factor, tolerance = 1e-15)
  expect_false(isTRUE(all.equal(
    drop(software$score), drop(paper$score), tolerance = 1e-15
  )))
  fit <- if_pca(
    x, 2L, threshold = 0, empirical_null = "none",
    ks_convention = "software"
  )
  expect_equal(fit$diagnostics$standard.deviation.divisor, nrow(x) - 1L)
  expect_match(fit$diagnostics$software.ks.rescaling, "divide W")
})


test_that("IF-PCA empirical-null and fixed threshold are exact", {
  x <- c7_if_fixture()
  raw <- drop(cpp_ch7_if_scores(x, FALSE)$score)
  expected <- (raw - mean(raw)) / stats::sd(raw)
  fit <- if_pca(
    x, 2L, selection = "fixed", threshold = 0,
    empirical_null = "mean_sd"
  )
  expect_equal(unname(drop(fit$raw.score)), unname(raw), tolerance = 2e-14)
  expect_equal(unname(drop(fit$adjusted.score)), unname(expected),
               tolerance = 2e-14)
  expect_identical(fit$selected, which(expected >= 0))
  expect_equal(fit$empirical.null$location, mean(raw), tolerance = 1e-15)
  expect_equal(fit$empirical.null$scale, stats::sd(raw), tolerance = 1e-15)
  expect_identical(
    fit$ranking, order(-expected, seq_along(expected))
  )
  expect_match(fit$diagnostics$selection.rule, ">=")
})


test_that("IF-PCA supplied-null HCT matches a hand calculation", {
  x <- c7_if_fixture()
  null <- seq(-1, 2, length.out = 301)
  fit <- if_pca(
    x, 2L, selection = "hct", empirical_null = "none",
    null_scores = null, hct_convention = "paper"
  )
  score <- drop(fit$adjusted.score)
  p.value <- vapply(score, function(value) mean(null > value), numeric(1))
  ordering <- order(p.value, seq_along(p.value))
  sorted <- p.value[ordering]
  rank.probability <- seq_along(score) / length(score)
  contrast <- rank.probability - sorted
  hc <- sqrt(length(score)) * contrast /
    sqrt(rank.probability + pmax(sqrt(nrow(x)) * contrast, 0))
  admissible <- seq_along(score) <= floor(length(score) / 2) &
    sorted > log(length(score)) / length(score)
  expected.rank <- which(admissible)[
    which.max(hc[admissible])
  ]
  expected.order <- order(-score, seq_along(score))
  expect_equal(unname(fit$hct$p.value), unname(p.value), tolerance = 0)
  expect_equal(fit$hct$rank.probability, rank.probability, tolerance = 0)
  expect_equal(unname(fit$hct$hc), unname(hc), tolerance = 2e-15)
  expect_identical(fit$hct$selected.rank, as.integer(expected.rank))
  expect_identical(fit$selected, expected.order[seq_len(expected.rank)])
  expect_equal(length(fit$selected), fit$hct$selected.rank)
  expect_match(fit$hct$selection.rule, "top j-hat")
})


test_that("IF-PCA HCT paper and software search and maximum ties are explicit", {
  score <- seq(20, 1)
  p.paper <- seq_len(20) / 20
  p.software <- seq_len(20) / 21
  paper <- .c7_if_hct(score, p.paper, 25L, "paper")
  software <- .c7_if_hct(score, p.software, 25L, "software")
  expect_identical(paper$selected.rank, 3L)
  expect_identical(software$selected.rank, 10L)
  expect_equal(paper$rank.probability, seq_len(20) / 20)
  expect_equal(software$rank.probability, seq_len(20) / 21)
  expect_match(paper$maximum.tie.rule, "first")
  expect_match(software$maximum.tie.rule, "last")
  expect_identical(paper$feature.tie.rule,
                   software$feature.tie.rule)
})


test_that("IF-PCA HCT demands an identified null and admissible rank", {
  x <- c7_if_fixture()
  expect_error(
    if_pca(x, 2L, selection = "hct", empirical_null = "none"),
    "exactly one"
  )
  expect_error(
    if_pca(
      x, 2L, selection = "hct", empirical_null = "none",
      null_scores = 1:10, null_cdf = stats::pnorm
    ),
    "exactly one"
  )
  expect_error(
    if_pca(
      x, 2L, selection = "hct", empirical_null = "none",
      null_cdf = function(value) 1
    ),
    "no admissible"
  )
  expect_error(
    if_pca(
      x, 2L, selection = "hct", empirical_null = "none",
      null_reps = 20L
    ),
    "seed"
  )
})


test_that("IF-PCA explicit null simulation is reproducible and restores RNG", {
  x <- c7_if_fixture()
  set.seed(707)
  before <- .Random.seed
  first <- if_pca(
    x, 2L, selection = "hct", empirical_null = "none",
    null_reps = 400L, seed = 901L
  )
  expect_identical(.Random.seed, before)
  second <- if_pca(
    x, 2L, selection = "hct", empirical_null = "none",
    null_reps = 400L, seed = 901L
  )
  expect_identical(first$null.calibration$raw.scores,
                   second$null.calibration$raw.scores)
  expect_identical(first$hct$p.value, second$hct$p.value)
  expect_identical(first$selected, second$selected)
  expect_identical(first$cluster, second$cluster)
  expect_identical(first$null.calibration$repetitions, 400L)
  expect_match(first$null.calibration$software.reference, "100\\*p")
  expect_match(first$null.calibration$paper.numeric.reference, "2000\\*p")
})


test_that("IF-PCA is invariant to feature translation and positive scaling", {
  x <- c7_if_fixture()
  base <- if_pca(x, 2L, threshold = 0, empirical_null = "none")
  moved <- sweep(x, 2L, seq_len(ncol(x)) * 10, "+")
  moved <- sweep(moved, 2L, seq_len(ncol(x)) + 0.5, "*")
  changed <- if_pca(moved, 2L, threshold = 0, empirical_null = "none")
  expect_equal(changed$raw.score, base$raw.score, tolerance = 2e-14)
  expect_equal(changed$standardized, base$standardized, tolerance = 2e-14)
  expect_identical(changed$selected, base$selected)
  expect_identical(
    outer(changed$cluster, changed$cluster, "=="),
    outer(base$cluster, base$cluster, "==")
  )
})


test_that("IF-PCA is equivariant to feature and row permutations", {
  x <- c7_if_fixture()
  base <- if_pca(x, 2L, threshold = 0, empirical_null = "none")
  columns <- c(6L, 2L, 5L, 1L, 4L, 3L)
  changed <- if_pca(
    x[, columns], 2L, threshold = 0, empirical_null = "none"
  )
  expect_equal(
    changed$raw.score[colnames(x)], base$raw.score, tolerance = 2e-14
  )
  expect_setequal(changed$selected.names, base$selected.names)
  expect_identical(
    outer(changed$cluster, changed$cluster, "=="),
    outer(base$cluster, base$cluster, "==")
  )
  rows <- c(12L, 1L, 11L, 2L, 10L, 3L, 9L, 4L, 8L, 5L, 7L, 6L)
  permuted <- if_pca(
    x[rows, ], 2L, threshold = 0, empirical_null = "none"
  )
  restored <- permuted$cluster[match(rownames(x), rownames(x)[rows])]
  expect_identical(
    outer(restored, restored, "=="),
    outer(base$cluster, base$cluster, "==")
  )
})


test_that("IF-PCA deterministic k-means has a fixed distance tie rule", {
  embedding <- matrix(c(-1, 1, 0), ncol = 1)
  fit <- cpp_ch7_if_deterministic_kmeans(
    embedding, 2L, 1e-12, 100L
  )
  expect_true(fit$valid)
  expect_true(fit$converged)
  expect_identical(drop(fit$cluster), c(1L, 2L, 1L))
  expect_match(fit$tie_rule, "smallest")
  again <- cpp_ch7_if_deterministic_kmeans(
    embedding, 2L, 1e-12, 100L
  )
  expect_identical(fit$cluster, again$cluster)
  expect_identical(fit$centers, again$centers)
})


test_that("IF-PCA fixed mode consumes no RNG and theoretical truncation is bounded", {
  x <- c7_if_fixture()
  set.seed(808)
  before <- .Random.seed
  fit <- if_pca(
    x, 2L, threshold = 0, empirical_null = "none", truncate = TRUE
  )
  expect_identical(.Random.seed, before)
  bound <- log(ncol(x)) / sqrt(nrow(x))
  expect_lte(max(abs(fit$embedding)), bound)
  expect_equal(fit$diagnostics$truncation.threshold, bound, tolerance = 0)
  expect_true(fit$diagnostics$truncation)
})


test_that("IF-PCA fails deterministically on zero variance and rank defects", {
  x <- c7_if_fixture()
  zero <- cbind(x, constant = 1)
  expect_error(
    if_pca(zero, 2L, threshold = 0, empirical_null = "none"),
    "zero or non-finite"
  )
  identical.features <- cbind(f1 = x[, 1L], f2 = x[, 1L])
  expect_error(
    if_pca(
      identical.features, 3L, threshold = 0, empirical_null = "none"
    ),
    "certified rank"
  )
  repeated <- x[, rep(1L, 3L), drop = FALSE]
  colnames(repeated) <- paste0("g", seq_len(ncol(repeated)))
  expect_error(
    if_pca(repeated, 2L, threshold = 0, empirical_null = "mean_sd"),
    "empirical-null scale"
  )
  expect_error(
    if_pca(x, 2L, threshold = 100, empirical_null = "none"),
    "selected no features"
  )
})


test_that("IF-PCA validates null arguments and names without silent fallback", {
  x <- c7_if_fixture()
  expect_error(
    if_pca(
      x, 2L, threshold = 0, empirical_null = "none",
      null_scores = 1:5
    ),
    "only used"
  )
  expect_error(
    if_pca(
      x, 2L, selection = "hct", empirical_null = "none",
      null_cdf = function(value) c(0.2, 0.3)
    ),
    "one finite"
  )
  duplicated <- x
  colnames(duplicated)[2L] <- colnames(duplicated)[1L]
  expect_error(
    if_pca(duplicated, 2L, threshold = 0, empirical_null = "none"),
    "unique"
  )
  expect_error(
    if_pca(x, nrow(x), threshold = 0, empirical_null = "none"),
    "smaller"
  )
  fit <- if_pca(x, 2L, threshold = 0, empirical_null = "none")
  expect_s3_class(fit, "if_pca_fit")
  expect_identical(fit$selected.names, colnames(x)[fit$selected])
  expect_match(fit$diagnostics$no.repair, "No zero-variance")
})
