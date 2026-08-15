.hdhr_test_base_x <- function() {
  rbind(
    c(2, 0), c(-2, 0), c(0, 1), c(0, -1),
    c(1, 1), c(-1, -1), c(1, -1), c(-1, 1)
  )
}

.hdhr_test_one_step_x <- function() {
  rbind(
    c(2, 0), c(-2, 0), c(1, 0), c(-1, 0),
    c(0, 1), c(0, -1), c(1, 1), c(-1, -1)
  )
}

.hdhr_test_invsqrt <- function(a) {
  eig <- eigen(a, symmetric = TRUE)
  eig$vectors %*% diag(eig$values^(-0.5), nrow(a)) %*%
    t(eig$vectors)
}

.hdhr_test_sqrt <- function(a) {
  eig <- eigen(a, symmetric = TRUE)
  eig$vectors %*% diag(sqrt(eig$values), nrow(a)) %*%
    t(eig$vectors)
}

.hdhr_test_trace_p <- function(a) {
  nrow(a) * a / sum(diag(a))
}


test_that("HDHR implements the literal one-step Algorithm 2 shape map", {
  x <- .hdhr_test_one_step_x()
  u <- x / sqrt(rowSums(x^2))
  s <- crossprod(u) / nrow(x)
  raw.ref <- .hdhr_test_trace_p(s)
  band.ref <- .hdhr_test_trace_p(diag(diag(s), 2L))

  expect_warning(
    fit <- high_dimensional_hr(
      x, diag(2), bandwidth = 0, max_iter = 1L,
      tol = 1e-14, strict = FALSE, scale_estimator = "none"
    ),
    "did not stabilize"
  )
  expect_false(fit$valid)
  expect_null(fit$estimate)
  expect_equal(fit$diagnostics$last.iterate$shape, band.ref,
               tolerance = 2e-13)
  expect_equal(fit$diagnostics$last.iterate$raw.shape, raw.ref,
               tolerance = 2e-13)
  expect_equal(fit$diagnostics$map$score.residual.l2, 0,
               tolerance = 1e-15)
  expect_equal(fit$diagnostics$map$shape.relative.update,
               norm(band.ref - diag(2), "F") / norm(diag(2), "F"),
               tolerance = 2e-13)
})


test_that("final SSCM, raw map, banding, and equation diagnostics are literal", {
  x <- .hdhr_test_base_x()
  fit <- high_dimensional_hr(x, diag(2), bandwidth = 0, tol = 1e-11)
  invsqrt <- .hdhr_test_invsqrt(fit$shape)
  sqrt.shape <- .hdhr_test_sqrt(fit$shape)
  residual <- sweep(x, 2, fit$location, "-") %*% invsqrt
  u <- residual / sqrt(rowSums(residual^2))
  s.ref <- crossprod(u) / nrow(x)
  b.ref <- diag(diag(s.ref), 2L)
  raw.ref <- .hdhr_test_trace_p(sqrt.shape %*% s.ref %*% sqrt.shape)
  shape.map <- .hdhr_test_trace_p(
    sqrt.shape %*% b.ref %*% sqrt.shape
  )

  expect_true(fit$valid)
  expect_true(fit$diagnostics$iteration.stable)
  expect_equal(fit$sscm, s.ref, tolerance = 3e-12)
  expect_equal(fit$banded.sscm, b.ref, tolerance = 3e-12)
  expect_equal(fit$raw.shape, raw.ref, tolerance = 3e-12)
  expect_equal(fit$diagnostics$shape.equation.residual,
               norm(shape.map - fit$shape, "F") /
                 max(1, norm(fit$shape, "F")),
               tolerance = 3e-12)
  expect_equal(fit$diagnostics$score.residual$l2,
               sqrt(sum(colMeans(u)^2)), tolerance = 3e-12)
  expect_equal(fit$diagnostics$score.residual$infinity,
               max(abs(colMeans(u))), tolerance = 3e-12)
})


test_that("shape trace, inverse, and scale identification are self-consistent", {
  x <- .hdhr_test_base_x()
  fit <- high_dimensional_hr(x, 7 * diag(2), bandwidth = 0)
  centered <- sweep(x, 2, colMeans(x), "-")
  trace.ref <- sum(centered^2) / (nrow(x) - 1)

  expect_true(fit$valid)
  expect_equal(sum(diag(fit$shape)), ncol(x), tolerance = 2e-13)
  expect_equal(fit$precision %*% fit$shape, diag(2), tolerance = 3e-11)
  expect_equal(fit$scatter.scale, trace.ref / ncol(x),
               tolerance = 2e-13)
  expect_equal(fit$scatter.trace, trace.ref, tolerance = 2e-13)
  expect_equal(fit$scatter, fit$scatter.scale * fit$shape,
               tolerance = 2e-13)
  expect_identical(fit$diagnostics$scale$robust, FALSE)
  expect_match(fit$diagnostics$scale$formula, "n-1", fixed = TRUE)
  expect_match(fit$diagnostics$scale$identification, "trace-p shape")
})


test_that("the default bandwidth is dimension-resolved and endpoints are explicit", {
  x2 <- .hdhr_test_base_x()
  fit2 <- high_dimensional_hr(x2, diag(2))
  x1 <- matrix(c(-3, -1, 1, 3), ncol = 1)
  fit1 <- high_dimensional_hr(x1, matrix(1, 1, 1))

  expect_identical(fit2$bandwidth$requested, 3L)
  expect_identical(fit2$bandwidth$resolved, 1L)
  expect_true(fit2$bandwidth$default.resolved.for.dimension)
  expect_identical(fit1$bandwidth$resolved, 0L)
  expect_equal(fit1$shape, matrix(1, 1, 1), tolerance = 1e-12)
  expect_equal(fit1$precision, matrix(1, 1, 1), tolerance = 1e-12)
  expect_error(high_dimensional_hr(x2, diag(2), bandwidth = 2L),
               "zero through p - 1")
  expect_error(high_dimensional_hr(x2, diag(2), bandwidth = -1L),
               "zero through p - 1")
})


test_that("translation and common scaling obey the identified transformations", {
  x <- .hdhr_test_base_x()
  shift <- c(1e8, -3e8)
  base <- high_dimensional_hr(x, diag(2), bandwidth = 0, tol = 1e-10)
  translated <- high_dimensional_hr(
    sweep(x, 2, shift, "+"), diag(2), bandwidth = 0, tol = 1e-10
  )
  big <- high_dimensional_hr(
    x * 1e150, diag(2), bandwidth = 0, tol = 1e-10
  )
  small <- high_dimensional_hr(
    x * 1e-150, diag(2), bandwidth = 0, tol = 1e-10
  )

  expect_equal(translated$location, base$location + shift,
               tolerance = 2e-8)
  expect_equal(translated$shape, base$shape, tolerance = 2e-9)
  expect_equal(translated$precision, base$precision, tolerance = 2e-9)
  expect_equal(translated$scatter.scale, base$scatter.scale,
               tolerance = 2e-9)
  expect_equal(big$location / 1e150, base$location, tolerance = 2e-9)
  expect_equal(big$shape, base$shape, tolerance = 2e-9)
  expect_equal(big$precision, base$precision, tolerance = 2e-9)
  expect_equal(big$scatter.scale / 1e300, base$scatter.scale,
               tolerance = 2e-9)
  expect_equal(small$location / 1e-150, base$location, tolerance = 2e-9)
  expect_equal(small$shape, base$shape, tolerance = 2e-9)
  expect_equal(small$scatter.scale / 1e-300, base$scatter.scale,
               tolerance = 2e-9)
})


test_that("signed permutations preserve the bandwidth-zero estimator", {
  x <- cbind(
    c(-3, -2, -1, 0, 1, 2, 3, 4),
    c(2, -1, 3, -2, 4, -3, 1, -4),
    c(-4, 1, -2, 3, -1, 4, -3, 2)
  )
  pilot <- diag(c(1.4, 0.8, 1.9))
  transform <- diag(c(-1, 1, -1))[c(3, 1, 2), ]
  transformed.x <- x %*% t(transform)
  transformed.pilot <- transform %*% pilot %*% t(transform)
  fit <- high_dimensional_hr(x, pilot, bandwidth = 0, tol = 2e-9)
  transformed.fit <- high_dimensional_hr(
    transformed.x, transformed.pilot, bandwidth = 0, tol = 2e-9
  )

  expect_true(fit$valid)
  expect_true(transformed.fit$valid)
  expect_equal(transformed.fit$location,
               as.numeric(transform %*% fit$location), tolerance = 2e-7)
  expect_equal(transformed.fit$shape,
               transform %*% fit$shape %*% t(transform), tolerance = 2e-7)
  expect_equal(transformed.fit$precision,
               transform %*% fit$precision %*% t(transform),
               tolerance = 3e-7)
  expect_equal(transformed.fit$scatter.scale, fit$scatter.scale,
               tolerance = 2e-12)
})


test_that("pilot common scale is immaterial but its certificate is retained", {
  x <- .hdhr_test_base_x()
  base <- high_dimensional_hr(x, diag(2), bandwidth = 0)
  rescaled <- high_dimensional_hr(x, 1e100 * diag(2), bandwidth = 0)
  tiny.pilot <- high_dimensional_hr(x, 1e-300 * diag(2), bandwidth = 0)
  certified <- structure(
    list(
      valid = TRUE,
      estimate = diag(2),
      method = "sglasso",
      diagnostics = list(solver = list(kkt.certified = TRUE))
    ),
    class = "spatial_sign_precision_fit"
  )
  from.fit <- high_dimensional_hr(x, certified, bandwidth = 0)

  expect_equal(rescaled$location, base$location, tolerance = 2e-10)
  expect_equal(rescaled$shape, base$shape, tolerance = 2e-10)
  expect_equal(rescaled$precision, base$precision, tolerance = 2e-10)
  expect_equal(tiny.pilot$shape, base$shape, tolerance = 2e-10)
  expect_equal(tiny.pilot$location, base$location, tolerance = 2e-10)
  expect_equal(from.fit$shape, base$shape, tolerance = 2e-10)
  expect_identical(from.fit$diagnostics$pilot$source,
                   "spatial_sign_precision_fit:sglasso")
  expect_true(from.fit$diagnostics$kkt$pilot$kkt.certified)
  expect_false(from.fit$diagnostics$kkt$applicable.to.hdhr)
  expect_true(from.fit$diagnostics$pilot$common.scale.normalized)
})


test_that("invalid precision fits and singular pilots never yield estimates", {
  x <- .hdhr_test_base_x()
  invalid.fit <- structure(
    list(valid = FALSE, estimate = NULL, diagnostics = list(solver = NULL)),
    class = "spatial_sign_precision_fit"
  )
  expect_warning(
    bad.fit <- high_dimensional_hr(
      x, invalid.fit, bandwidth = 0, strict = FALSE
    ),
    "invalid or uncertified"
  )
  expect_false(bad.fit$valid)
  expect_null(bad.fit$estimate)

  singular <- diag(c(1, 0))
  expect_error(high_dimensional_hr(x, singular, bandwidth = 0),
               "not positive definite")
  expect_warning(
    bad <- high_dimensional_hr(
      x, singular, bandwidth = 0, strict = FALSE
    ),
    "pilot precision"
  )
  expect_false(bad$valid)
  expect_null(bad$estimate)
  expect_identical(bad$diagnostics$failure.stage, "pilot precision")
  expect_match(bad$diagnostics$no.repair, "pseudoinverse")
})


test_that("nonsymmetric and mismatched pilots are not silently repaired", {
  x <- .hdhr_test_base_x()
  nonsymmetric <- matrix(c(1, 0.2, 0, 1), 2, 2)
  expect_warning(
    bad <- high_dimensional_hr(
      x, nonsymmetric, bandwidth = 0, strict = FALSE
    ),
    "not symmetric"
  )
  expect_false(bad$valid)
  expect_null(bad$estimate)
  tiny.nonsymmetric <- 1e-300 * matrix(c(1, 0.2, 0, 1), 2, 2)
  expect_warning(
    tiny.bad <- high_dimensional_hr(
      x, tiny.nonsymmetric, bandwidth = 0, strict = FALSE
    ),
    "not symmetric"
  )
  expect_false(tiny.bad$valid)
  expect_null(tiny.bad$estimate)
  expect_error(high_dimensional_hr(x, diag(3), bandwidth = 0),
               "p by p")
  expect_error(high_dimensional_hr(x, matrix(NA_real_, 2, 2), bandwidth = 0),
               "finite")
  expect_error(high_dimensional_hr(x, bandwidth = 0),
               "pilot_precision")
})


test_that("banding-induced non-PD and unbanded p greater than n fail", {
  indefinite.x <- rbind(
    c(1, 1, 1),
    c(-1, -1, -1),
    c(2, 2, 2),
    c(-2, -2, -2)
  )
  expect_warning(
    bad.band <- high_dimensional_hr(
      indefinite.x, diag(3), bandwidth = 1,
      strict = FALSE, scale_estimator = "none"
    ),
    "banded standardized-sign SSCM"
  )
  expect_false(bad.band$valid)
  expect_null(bad.band$estimate)
  expect_identical(bad.band$diagnostics$failure.stage,
                   "banded sign SSCM")
  expect_match(bad.band$diagnostics$failure, "No ridge")

  wide.x <- rbind(
    c(1, 2, 3, 4, 5), c(-1, -2, -2, -4, -3),
    c(2, -1, 4, -3, 1), c(-2, 1, -4, 3, -1)
  )
  expect_warning(
    wide <- high_dimensional_hr(
      wide.x, diag(5), bandwidth = 4,
      strict = FALSE, scale_estimator = "none"
    ),
    "structurally singular"
  )
  expect_false(wide$valid)
  expect_null(wide$estimate)
  expect_identical(wide$bandwidth$resolved, 4L)
})


test_that("p at least n remains executable with diagonal banding", {
  wide.x <- rbind(
    c(1, 2, 3, 4, 5, 6),
    c(-2, 1, -4, 3, -6, 5),
    c(3, -4, 1, -6, 2, -5),
    c(-4, -3, 6, 1, 5, 2)
  )
  fit <- high_dimensional_hr(
    wide.x, diag(6), bandwidth = 0,
    tol = 1e-7, max_iter = 3000L, scale_estimator = "none"
  )
  expect_true(fit$valid)
  expect_equal(dim(fit$shape), c(6L, 6L))
  expect_equal(sum(diag(fit$shape)), 6, tolerance = 2e-11)
  expect_true(all(eigen(fit$shape, symmetric = TRUE,
                        only.values = TRUE)$values > 0))
  expect_null(fit$scatter)
  expect_true(is.na(fit$scatter.scale))
  expect_true(fit$diagnostics$iteration.stable)
})


test_that("exact zero residuals are reported without perturbation", {
  x <- rbind(c(0, 0), c(1, 0), c(-1, 0), c(0, 1), c(0, -1))
  expect_error(high_dimensional_hr(x, diag(2), bandwidth = 0),
               "zero or non-finite")
  expect_warning(
    bad <- high_dimensional_hr(
      x, diag(2), bandwidth = 0, strict = FALSE
    ),
    "inverse-radius update is undefined"
  )
  expect_false(bad$valid)
  expect_null(bad$estimate)
  expect_identical(bad$diagnostics$failure.stage,
                   "standardized residual")
  expect_match(bad$diagnostics$no.repair, "zero-residual")
})


test_that("joint and spatial-median nonconvergence have explicit contracts", {
  x <- .hdhr_test_one_step_x()
  expect_warning(
    joint <- high_dimensional_hr(
      x, diag(2), bandwidth = 0, tol = 1e-15, max_iter = 1L,
      strict = FALSE, scale_estimator = "none"
    ),
    "did not stabilize"
  )
  expect_false(joint$valid)
  expect_null(joint$estimate)
  expect_identical(joint$diagnostics$failure.stage,
                   "joint fixed-point iteration")
  expect_identical(joint$diagnostics$iterations, 1L)
  expect_false(joint$diagnostics$iteration.stable)

  asymmetric <- rbind(c(0, 0), c(10, 0), c(0, 2), c(7, 9), c(-2, 5))
  expect_warning(
    median <- high_dimensional_hr(
      asymmetric, diag(2), bandwidth = 0,
      median_tol = 1e-16, median_max_iter = 1L,
      strict = FALSE, scale_estimator = "none"
    ),
    "spatial median"
  )
  expect_false(median$valid)
  expect_null(median$estimate)
  expect_identical(median$diagnostics$failure.stage,
                   "initial spatial median")
  expect_false(median$diagnostics$median$stable)
})


test_that("extreme unrepresentable scale is explicit and shape-only is usable", {
  x <- .hdhr_test_base_x() * 1e200
  expect_warning(
    scaled <- high_dimensional_hr(
      x, diag(2), bandwidth = 0, strict = FALSE
    ),
    "outside the finite positive double range"
  )
  expect_false(scaled$valid)
  expect_null(scaled$estimate)
  expect_identical(scaled$diagnostics$failure.stage, "QDA scatter scale")

  shape.only <- high_dimensional_hr(
    x, diag(2), bandwidth = 0, scale_estimator = "none"
  )
  reference <- high_dimensional_hr(
    .hdhr_test_base_x(), diag(2), bandwidth = 0,
    scale_estimator = "none"
  )
  expect_true(shape.only$valid)
  expect_equal(shape.only$location / 1e200, reference$location,
               tolerance = 3e-9)
  expect_equal(shape.only$shape, reference$shape, tolerance = 3e-9)
  expect_null(shape.only$scatter)
  expect_true(is.na(shape.only$scatter.scale))
})


test_that("underflowing QDA scale is likewise never floored", {
  x <- .hdhr_test_base_x() * 1e-200
  expect_warning(
    scaled <- high_dimensional_hr(
      x, diag(2), bandwidth = 0, strict = FALSE
    ),
    "outside the finite positive double range"
  )
  expect_false(scaled$valid)
  expect_null(scaled$estimate)
  expect_identical(scaled$diagnostics$failure.stage, "QDA scatter scale")
  expect_match(scaled$diagnostics$no.repair, "eigenvalue floor")
})


test_that("unrepresentable within-matrix dynamic range fails during scaling", {
  x <- rbind(c(0, 0), c(1e308, 0), c(1e-308, 1))
  expect_warning(
    bad <- high_dimensional_hr(
      x, diag(2), bandwidth = 0, strict = FALSE,
      scale_estimator = "none"
    ),
    "underflows double precision"
  )
  expect_false(bad$valid)
  expect_null(bad$estimate)
  expect_identical(bad$diagnostics$failure.stage, "data scaling")
  expect_match(bad$diagnostics$no.repair, "No ridge")
})


test_that("constant data and malformed data are rejected", {
  constant <- matrix(4, 6, 3)
  expect_error(high_dimensional_hr(constant, diag(3), bandwidth = 0),
               "All observations coincide")
  expect_warning(
    bad <- high_dimensional_hr(
      constant, diag(3), bandwidth = 0, strict = FALSE
    ),
    "All observations coincide"
  )
  expect_false(bad$valid)
  expect_null(bad$estimate)
  expect_error(high_dimensional_hr(matrix(1, 1, 2), diag(2)),
               "at least two rows")
  expect_error(high_dimensional_hr(matrix(c(1, NA, 2, 3), 2), diag(2)),
               "finite")
  expect_error(high_dimensional_hr(data.frame(a = 1:3, b = letters[1:3]),
                                   diag(2)), "Every column")
})


test_that("control arguments are validated exactly", {
  x <- .hdhr_test_base_x()
  expect_error(high_dimensional_hr(x, diag(2), bandwidth = 0, tol = 0),
               "tol")
  expect_error(high_dimensional_hr(x, diag(2), bandwidth = 0, max_iter = 1.5),
               "max_iter")
  expect_error(high_dimensional_hr(x, diag(2), bandwidth = 0,
                                   median_tol = Inf), "median_tol")
  expect_error(high_dimensional_hr(x, diag(2), bandwidth = 0,
                                   median_max_iter = 0), "median_max_iter")
  expect_error(high_dimensional_hr(x, diag(2), bandwidth = 0,
                                   zero_tol = -1), "zero_tol")
  expect_error(high_dimensional_hr(x, diag(2), bandwidth = 0,
                                   scale_estimator = "invented"),
               "arg")
  expect_error(high_dimensional_hr(x, diag(2), bandwidth = 0,
                                   strict = NA), "strict")
})


test_that("names, diagnostics, and errata remain explicit", {
  x <- .hdhr_test_base_x()
  colnames(x) <- c("alpha", "beta")
  rownames(x) <- paste0("obs", seq_len(nrow(x)))
  fit <- high_dimensional_hr(x, diag(2), bandwidth = 0)

  expect_identical(names(fit$location), colnames(x))
  expect_identical(dimnames(fit$shape), list(colnames(x), colnames(x)))
  expect_identical(dimnames(fit$precision), list(colnames(x), colnames(x)))
  expect_identical(colnames(fit$signs), colnames(x))
  expect_identical(rownames(fit$signs), rownames(x))
  expect_null(fit$diagnostics$failure)
  expect_match(fit$diagnostics$convergence.basis, "software stopping rule")
  expect_match(fit$diagnostics$primary.algorithm,
               "band the standardized-sign SSCM")
  expect_match(paste(fit$diagnostics$book.errata, collapse = " "),
               "does not band the final inverse")
  expect_match(paste(fit$diagnostics$book.errata, collapse = " "),
               "not a finite-sample guarantee")
  expect_match(fit$diagnostics$invariance, "not arbitrary")
  expect_match(fit$diagnostics$no.repair, "jitter")
  expect_true(is.finite(fit$diagnostics$positive.definiteness$shape.rcond))
  expect_gt(fit$diagnostics$positive.definiteness$shape.rcond, 0)
})
