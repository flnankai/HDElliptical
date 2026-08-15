zfrd_fixture <- function() {
  matrix(c(
     1.8,  0.2, -0.7,
    -1.1,  1.5,  0.4,
     0.6, -1.7,  1.3,
     2.2,  0.9,  0.8,
    -0.5, -0.8, -2.0,
     1.1,  2.1, -0.4,
    -2.3,  0.4,  1.0,
     0.3, -2.4, -1.1,
     1.5, -0.6,  2.3,
    -1.7, -1.3,  0.6
  ), ncol = 3, byrow = TRUE,
  dimnames = list(paste0("row", 1:10), c("a", "b", "c")))
}

zfrd_reference <- function(x, location, shape) {
  decomposition <- eigen(shape, symmetric = TRUE)
  inverse.root <- decomposition$vectors %*%
    diag(1 / sqrt(decomposition$values), nrow(shape)) %*%
    t(decomposition$vectors)
  standardized <- sweep(x, 2, location, "-") %*% inverse.root
  radii <- sqrt(rowSums(standardized^2))
  directions <- standardized / radii
  correlations <- vapply(
    seq_len(ncol(x)),
    function(j) cor(log(radii), directions[, j]),
    numeric(1)
  )
  n <- nrow(x)
  p <- ncol(x)
  list(
    standardized = standardized,
    radii = radii,
    log.radii = log(radii),
    directions = directions,
    correlations = correlations,
    T.sum = n * sum(correlations^2),
    T.max.raw = n * max(correlations^2),
    T.max = n * max(correlations^2) - 2 * log(p) + log(log(p))
  )
}

zfrd_bootstrap_reference <- function(log.radii, B, seed) {
  n <- length(log.radii)
  p <- 3L
  set.seed(seed)
  sum.values <- max.values <- numeric(B)
  for (b in seq_len(B)) {
    indices <- integer(n)
    directions <- matrix(NA_real_, n, p)
    for (i in seq_len(n)) {
      indices[i] <- floor(runif(1) * n) + 1L
      row <- rnorm(p)
      directions[i, ] <- row / sqrt(sum(row^2))
    }
    correlations <- vapply(
      seq_len(p),
      function(j) cor(log.radii[indices], directions[, j]),
      numeric(1)
    )
    sum.values[b] <- n * sum(correlations^2)
    max.values[b] <- n * max(correlations^2) -
      2 * log(p) + log(log(p))
  }
  list(sum = sum.values, max = max.values)
}

test_that("radial-directional statistics match a literal full-matrix reference", {
  x <- zfrd_fixture()
  location <- c(0.15, -0.2, 0.1)
  shape <- matrix(c(1.4, 0.25, -0.1,
                    0.25, 0.9, 0.18,
                    -0.1, 0.18, 1.2), 3, 3, byrow = TRUE)
  reference <- zfrd_reference(x, location, shape)
  result <- zhang_feng_radial_directional_test(
    x, location = location, shape = shape
  )

  expect_s3_class(result, "radial_directional_test")
  expect_s3_class(result, "htest")
  expect_equal(unname(result$estimate), reference$correlations,
               tolerance = 2e-13)
  expect_equal(result$components$T.sum, reference$T.sum,
               tolerance = 2e-13)
  expect_equal(result$components$T.max.raw, reference$T.max.raw,
               tolerance = 2e-13)
  expect_equal(result$components$T.max, reference$T.max,
               tolerance = 2e-13)
  expect_equal(result$components$log.radii, reference$log.radii,
               tolerance = 2e-13, ignore_attr = TRUE)
  expect_equal(unname(result$components$directions),
               unname(reference$directions), tolerance = 2e-13)
  z <- (reference$T.sum - 3) / sqrt(6)
  expect_equal(result$components$Z.sum, z, tolerance = 1e-14)
  expect_equal(result$components$p.values[["sum"]],
               pnorm(z, lower.tail = FALSE), tolerance = 1e-15)
  expected.max <- -expm1(
    -exp(-reference$T.max / 2) / sqrt(pi)
  )
  expect_equal(result$components$p.values[["max"]], expected.max,
               tolerance = 1e-15)
  expected.cauchy <- 0.5 - atan(
    0.5 * tanpi(0.5 - pnorm(z, lower.tail = FALSE)) +
      0.5 * tanpi(0.5 - expected.max)
  ) / pi
  expect_equal(result$p.value, expected.cauchy, tolerance = 2e-15)
  expect_identical(result$components$largest.coordinate,
                   names(which.max(abs(result$estimate))))
})

test_that("sum and max components can be selected without changing the fit", {
  x <- zfrd_fixture()
  common <- list(x = x, location = c(0.15, -0.2, 0.1), shape = diag(3))
  combined <- do.call(zhang_feng_radial_directional_test, common)
  sum.test <- do.call(
    zhang_feng_radial_directional_test,
    c(common, list(component = "sum"))
  )
  max.test <- do.call(
    zhang_feng_radial_directional_test,
    c(common, list(component = "max"))
  )
  expect_equal(sum.test$p.value, combined$components$p.values[["sum"]])
  expect_equal(max.test$p.value, combined$components$p.values[["max"]])
  expect_equal(sum.test$estimate, combined$estimate)
  expect_equal(max.test$estimate, combined$estimate)
  expect_named(sum.test$statistic, "Normal.Z")
  expect_named(max.test$statistic, "Gumbel.centered")
  expect_named(combined$statistic, "Cauchy.angle")
})

test_that("exact and basis-sensitive transformations match the documented boundary", {
  x <- zfrd_fixture()
  location <- c(0.15, -0.2, 0.1)
  shape <- matrix(c(1.4, 0.25, -0.1,
                    0.25, 0.9, 0.18,
                    -0.1, 0.18, 1.2), 3, 3, byrow = TRUE)
  baseline <- zhang_feng_radial_directional_test(
    x, location = location, shape = shape
  )
  huge <- zhang_feng_radial_directional_test(
    x * 1e150, location = location * 1e150,
    shape = shape * 1e300
  )
  expect_equal(huge$estimate, baseline$estimate, tolerance = 2e-13)
  expect_equal(huge$components$T.sum, baseline$components$T.sum,
               tolerance = 2e-13)
  expect_equal(huge$components$T.max, baseline$components$T.max,
               tolerance = 2e-13)

  transform <- matrix(c(1.3, 0.2, -0.1,
                        0.4, 0.9, 0.3,
                        -0.2, 0.1, 1.1), 3, 3, byrow = TRUE)
  shift <- c(4, -3, 2)
  transformed.x <- sweep(x %*% t(transform), 2, shift, "+")
  transformed.location <- as.vector(transform %*% location) + shift
  transformed.shape <- transform %*% shape %*% t(transform)
  transformed <- zhang_feng_radial_directional_test(
    transformed.x, location = transformed.location,
    shape = transformed.shape, component = "sum"
  )
  # Coordinatewise empirical variance normalisation prevents exact finite-
  # sample rotation invariance, even though the null theory is asymptotically
  # rotation-compatible after a consistently transformed shape.
  expect_gt(abs(transformed$components$T.sum -
                  baseline$components$T.sum), 1e-4)

  signed.permutation <- matrix(c(0, -1, 0,
                                 1,  0, 0,
                                 0,  0, 1), 3, 3, byrow = TRUE)
  permuted <- zhang_feng_radial_directional_test(
    x %*% t(signed.permutation),
    location = as.vector(signed.permutation %*% location),
    shape = signed.permutation %*% shape %*% t(signed.permutation)
  )
  expect_equal(permuted$components$T.max, baseline$components$T.max,
               tolerance = 2e-13)
  expect_equal(permuted$components$T.sum, baseline$components$T.sum,
               tolerance = 2e-13)
})

test_that("validated and internally fitted HR routes agree with explicit fits", {
  x <- zfrd_fixture()
  fitted <- high_dimensional_hr(
    x, diag(3), bandwidth = 2L, scale_estimator = "none",
    tol = 1e-7, max_iter = 1000L
  )
  via.fit <- zhang_feng_radial_directional_test(x, fit = fitted)
  explicit <- zhang_feng_radial_directional_test(
    x, location = fitted$location, shape = fitted$shape
  )
  internal <- zhang_feng_radial_directional_test(
    x, pilot_precision = diag(3), bandwidth = 2L,
    tol = 1e-7, max_iter = 1000L
  )
  expect_equal(via.fit$estimate, explicit$estimate, tolerance = 1e-14)
  expect_equal(via.fit$components$T.sum, explicit$components$T.sum,
               tolerance = 1e-14)
  expect_equal(internal$estimate, via.fit$estimate, tolerance = 1e-14)
  expect_match(via.fit$diagnostics$standardisation.source,
               "high_dimensional_hr_fit")
  expect_match(internal$diagnostics$standardisation.source,
               "internally fitted")
})

test_that("intrinsic bootstrap matches the literal radial-directional resampling", {
  x <- zfrd_fixture()
  B <- 24L
  seed <- 281L
  set.seed(77)
  state <- .Random.seed
  result <- zhang_feng_radial_directional_test(
    x, location = c(0.15, -0.2, 0.1), shape = diag(3),
    calibration = "radial_bootstrap", B = B, seed = seed,
    keep_bootstrap = TRUE
  )
  expect_identical(.Random.seed, state)
  reference <- zfrd_bootstrap_reference(
    unname(result$components$log.radii), B, seed
  )
  expect_equal(result$components$bootstrap$values$T.sum,
               reference$sum, tolerance = 2e-13)
  expect_equal(result$components$bootstrap$values$T.max,
               reference$max, tolerance = 2e-13)
  expect_equal(result$components$bootstrap$mean.sum,
               mean(reference$sum), tolerance = 2e-14)
  expect_equal(result$components$bootstrap$sd.sum,
               sd(reference$sum), tolerance = 2e-14)
  expect_equal(result$components$bootstrap$mean.max,
               mean(reference$max), tolerance = 2e-14)
  expect_equal(result$components$bootstrap$sd.max,
               sd(reference$max), tolerance = 2e-14)

  z <- (result$components$T.sum - mean(reference$sum)) /
    sd(reference$sum)
  expect_equal(result$components$Z.sum, z, tolerance = 2e-14)
  mu.g <- 2 * 0.5772156649015328606 - log(pi)
  sigma.g <- sqrt(2 * pi^2 / 3)
  g <- mu.g + sigma.g / sd(reference$max) *
    (result$components$T.max - mean(reference$max))
  expect_equal(result$components$Gumbel.argument, g, tolerance = 2e-14)
  expect_true(result$diagnostics$intrinsic.calibration)
  expect_false(result$diagnostics$simulation.replication)
})

test_that("NULL seed uses the caller stream normally", {
  x <- zfrd_fixture()
  set.seed(812)
  before <- .Random.seed
  invisible(zhang_feng_radial_directional_test(
    x, location = c(0.15, -0.2, 0.1), shape = diag(3),
    calibration = "radial_bootstrap", B = 4L
  ))
  expect_false(identical(.Random.seed, before))
})

test_that("standardisation and calibration contracts reject ambiguous inputs", {
  x <- zfrd_fixture()
  expect_error(
    zhang_feng_radial_directional_test(x),
    "Supply a valid"
  )
  expect_error(
    zhang_feng_radial_directional_test(
      x, location = numeric(3), shape = diag(3),
      pilot_precision = diag(3)
    ),
    "cannot be combined"
  )
  expect_error(
    zhang_feng_radial_directional_test(
      x, location = numeric(3), shape = diag(3), bandwidth = 1L
    ),
    "used only"
  )
  expect_error(
    zhang_feng_radial_directional_test(x, location = numeric(3)),
    "supplied together"
  )
  expect_error(
    zhang_feng_radial_directional_test(
      x, location = numeric(3), shape = diag(c(1, 1, -1))
    ),
    "positive definite"
  )
  nonsymmetric <- diag(3)
  nonsymmetric[1, 2] <- 0.2
  expect_error(
    zhang_feng_radial_directional_test(
      x, location = numeric(3), shape = nonsymmetric
    ),
    "symmetric"
  )
  expect_error(
    zhang_feng_radial_directional_test(
      x[, 1, drop = FALSE], location = 0, shape = matrix(1)
    ),
    "at least two variables"
  )
  expect_error(
    zhang_feng_radial_directional_test(
      x, location = numeric(3), shape = diag(3),
      calibration = "radial_bootstrap", B = 1
    ),
    "at least 2"
  )
  expect_error(
    zhang_feng_radial_directional_test(
      x, location = numeric(3), shape = diag(3), seed = 1
    ),
    "used only"
  )
  expect_error(
    zhang_feng_radial_directional_test(
      x, location = numeric(3), shape = diag(3), keep_bootstrap = TRUE
    ),
    "requires radial-bootstrap"
  )
})

test_that("undefined radial and directional variances are never repaired", {
  circle <- rbind(c(1, 0), c(0, 1), c(-1, 0), c(0, -1),
                  c(0.6, 0.8), c(-0.6, -0.8))
  expect_error(
    zhang_feng_radial_directional_test(
      circle, location = c(0, 0), shape = diag(2)
    ),
    "log radii have zero"
  )
  ray <- cbind(1:8, 0)
  expect_error(
    zhang_feng_radial_directional_test(
      ray, location = c(0, 0), shape = diag(2)
    ),
    "direction coordinate has zero"
  )
  x <- zfrd_fixture()
  x[1, ] <- 0
  expect_error(
    zhang_feng_radial_directional_test(
      x, location = c(0, 0, 0), shape = diag(3)
    ),
    "residual is exactly zero"
  )
  tiny <- zfrd_fixture()
  tiny[1, ] <- tiny[1, ] * 1e-200
  expect_error(
    zhang_feng_radial_directional_test(
      tiny, location = numeric(3), shape = diag(3),
      radius_zero_tol = 1e-150
    ),
    "radius_zero_tol"
  )
})

test_that("diagnostics disclose target, invariance, and no-repair boundary", {
  result <- zhang_feng_radial_directional_test(
    zfrd_fixture(), location = c(0.15, -0.2, 0.1), shape = diag(3)
  )
  expect_match(result$diagnostics$method.target, "necessary elliptical")
  expect_match(result$diagnostics$invariance, "basis")
  expect_match(result$diagnostics$primary.numerical.boundary,
               "not used")
  expect_match(result$diagnostics$book.completeness, "omits")
  expect_match(result$diagnostics$no.repair, "No ridge")
  expect_equal(result$components$p.values[["combined"]], result$p.value)
  expect_true(is.finite(result$components$largest.absolute.correlation))
})
