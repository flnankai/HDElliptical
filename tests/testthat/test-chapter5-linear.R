# Deterministic tests for Chapter 5 sparse linear classifiers.

c5lin_fixture <- function() {
  first <- rbind(
    c(2.0, 1.0, 0.0),
    c(1.0, 2.0, 1.0),
    c(3.0, 1.0, -1.0),
    c(2.0, 3.0, 0.5),
    c(1.5, 1.0, -0.5),
    c(2.5, 2.0, 1.5)
  )
  second <- -first
  x <- rbind(first, second)
  colnames(x) <- c("a", "b", "c")
  y <- factor(
    rep(c("class-a", "class-b"), each = nrow(first)),
    levels = c("class-a", "class-b")
  )
  list(x = x, y = y, first = first, second = second)
}


c5lin_manual_signs <- function(x, center) {
  residual <- sweep(x, 2L, center, "-")
  radii <- sqrt(rowSums(residual^2))
  answer <- matrix(0, nrow(x), ncol(x))
  positive <- radii > 0
  answer[positive, ] <- residual[positive, , drop = FALSE] /
    radii[positive]
  answer
}


c5lin_kkt_residual <- function(x, response, coefficient, intercept,
                               lambda) {
  residual <- response - drop(intercept + x %*% coefficient)
  gradient <- -2 * drop(crossprod(x, residual)) / nrow(x)
  components <- numeric(length(coefficient))
  positive <- coefficient > 0
  negative <- coefficient < 0
  zero <- !(positive | negative)
  components[positive] <- abs(gradient[positive] + lambda)
  components[negative] <- abs(gradient[negative] - lambda)
  components[zero] <- pmax(0, abs(gradient[zero]) - lambda)
  max(abs(-2 * mean(residual)), components)
}


test_that("Dantzig Rcpp kernel matches a diagonal problem exactly", {
  fit <- HDElliptical:::cpp_c5lin_dantzig(
    diag(c(2, 4)), c(1, -3), lambda = 0.2,
    tolerance = 1e-9, maximum_iterations = 100000L
  )

  expect_true(fit$converged)
  expect_equal(drop(fit$solution), c(0.4, -0.7), tolerance = 2e-8)
  expect_lte(fit$primal_violation, 1e-9)
  expect_lte(fit$dual_violation, 1e-9)
  expect_lte(fit$stationarity_residual, 1e-9)
  expect_lte(fit$relative_gap, 1e-9)
  expect_equal(
    max(abs(diag(c(2, 4)) %*% fit$solution - c(1, -3))),
    0.2, tolerance = 2e-8
  )
})


test_that("Dantzig zero and zero-operator boundaries are explicit", {
  zero.optimum <- HDElliptical:::cpp_c5lin_dantzig(
    diag(2), c(0.1, -0.2), lambda = 0.2,
    tolerance = 1e-9, maximum_iterations = 10L
  )
  expect_true(zero.optimum$converged)
  expect_equal(drop(zero.optimum$solution), c(0, 0), tolerance = 0)
  expect_identical(zero.optimum$iterations, 0L)
  expect_equal(zero.optimum$primal_objective, 0, tolerance = 0)

  feasible.zero.operator <- HDElliptical:::cpp_c5lin_dantzig(
    matrix(0, 2, 2), c(0.1, -0.2), lambda = 0.2,
    tolerance = 1e-9, maximum_iterations = 10L
  )
  expect_true(feasible.zero.operator$converged)
  expect_true(feasible.zero.operator$operator_zero)
  expect_false(feasible.zero.operator$infeasible_zero_operator)

  infeasible <- HDElliptical:::cpp_c5lin_dantzig(
    matrix(0, 2, 2), c(0.1, -0.3), lambda = 0.2,
    tolerance = 1e-9, maximum_iterations = 10L
  )
  expect_false(infeasible$converged)
  expect_true(infeasible$operator_zero)
  expect_true(infeasible$infeasible_zero_operator)
  expect_gt(infeasible$primal_violation, 0)
})


test_that("LPD uses the unbiased pooled covariance and certified constraint", {
  z <- c5lin_fixture()
  fit <- lpd_classifier(
    z$x, z$y, lambda = 0.2, prior = c(0.7, 0.3),
    ridge = 0.15, solver_tol = 1e-8
  )
  mean1 <- colMeans(z$first)
  mean2 <- colMeans(z$second)
  residual1 <- sweep(z$first, 2L, mean1, "-")
  residual2 <- sweep(z$second, 2L, mean2, "-")
  pooled <- (crossprod(residual1) + crossprod(residual2)) /
    (nrow(z$x) - 2)
  operator <- pooled
  diag(operator) <- diag(operator) + 0.15
  direction <- fit$estimate$direction

  expect_s3_class(fit, "hd_classifier_fit")
  expect_true(fit$valid)
  expect_equal(unname(fit$estimate$pooled.covariance), pooled,
               tolerance = 1e-14)
  expect_equal(unname(fit$estimate$constraint.operator), operator,
               tolerance = 1e-14)
  expect_equal(
    max(abs(drop(operator %*% direction) - (mean1 - mean2))),
    0.2, tolerance = 2e-7
  )
  expect_true(fit$diagnostics$solver$certified)
  expect_equal(
    fit$score.model$intercept,
    -sum((mean1 + mean2) / 2 * direction) + log(0.7 / 0.3),
    tolerance = 1e-13
  )
  expect_identical(fit$diagnostics$pooled.covariance.divisor, 10)
  expect_true(fit$tuning$ridge.changes.operator)
})


test_that("LPD score and prior threshold use class1-positive orientation", {
  z <- c5lin_fixture()
  fit <- lpd_classifier(
    z$x, z$y, lambda = 0.2, prior = c(0.8, 0.2),
    ridge = 0.1
  )
  manual <- drop(z$x %*% fit$estimate$direction) +
    fit$score.model$intercept
  score <- predict(fit, z$x, type = "score")
  predicted <- predict(fit, z$x, type = "class")

  expect_equal(score, manual, tolerance = 1e-13)
  expect_identical(levels(predicted), levels(z$y))
  expect_true(all(predicted[seq_len(6)] == "class-a"))
  expect_true(all(predicted[7:12] == "class-b"))
  expect_identical(fit$primary.orientation,
                   "positive score selects class1")
  expect_match(fit$score.scale, "log\\(pi1/pi2\\)")
})


test_that("LPD translation and feature-name contracts are exact", {
  z <- c5lin_fixture()
  shift <- c(a = 100, b = -50, c = 7)
  moved <- sweep(z$x, 2L, shift, "+")
  base <- lpd_classifier(z$x, z$y, 0.2, ridge = 0.1)
  translated <- lpd_classifier(moved, z$y, 0.2, ridge = 0.1)

  expect_equal(translated$estimate$direction,
               base$estimate$direction, tolerance = 2e-8)
  expect_equal(
    predict(translated, moved, type = "score"),
    predict(base, z$x, type = "score"),
    tolerance = 2e-7
  )
  reordered <- moved[, c("b", "a", "c")]
  expect_equal(
    predict(translated, reordered, type = "score"),
    predict(translated, moved, type = "score"),
    tolerance = 0
  )
})


test_that("LPD never promotes an uncertified last iterate", {
  z <- c5lin_fixture()
  expect_error(
    lpd_classifier(
      z$x, z$y, lambda = 0.01, ridge = 0.1,
      solver_tol = 1e-14, solver_max_iter = 1L
    ),
    "failed at Dantzig optimization"
  )
  invalid <- NULL
  expect_warning(
    invalid <- lpd_classifier(
      z$x, z$y, lambda = 0.01, ridge = 0.1,
      solver_tol = 1e-14, solver_max_iter = 1L,
      strict = FALSE
    ),
    "failed at Dantzig optimization"
  )
  expect_s3_class(invalid, "hd_classifier_fit")
  expect_false(invalid$valid)
  expect_length(invalid$estimate, 0L)
  expect_false(is.null(invalid$diagnostics$solver$last.iterate))
  expect_error(predict(invalid, z$x), "invalid|failed",
               ignore.case = TRUE)
})


test_that("DSDA uses the primary response coding and n-inverse objective", {
  z <- c5lin_fixture()
  lambda <- 0.15
  fit <- dsda_classifier(
    z$x, z$y, lambda = lambda, solver_tol = 1e-9
  )
  response <- c(rep(-2, 6), rep(2, 6))
  beta <- fit$estimate$regression.direction
  intercept <- fit$estimate$regression.intercept
  residual <- response - drop(intercept + z$x %*% beta)
  objective <- mean(residual^2) + lambda * sum(abs(beta))

  expect_true(fit$valid)
  expect_equal(fit$estimate$response, response, tolerance = 0)
  expect_equal(fit$diagnostics$response.coding,
               c(class1 = -2, class2 = 2), tolerance = 0)
  expect_equal(fit$diagnostics$solver$primal.objective,
               objective, tolerance = 2e-12)
  expect_lt(abs(
    c5lin_kkt_residual(z$x, response, beta, intercept, lambda) -
      fit$diagnostics$solver$kkt.residual
  ), 1e-12)
  expect_lte(fit$diagnostics$solver$kkt.residual, 1e-9)
  expect_lte(fit$diagnostics$solver$relative.gap, 1e-9)
  expect_match(fit$diagnostics$objective, "n\\^\\{-1\\}")
})


test_that("DSDA prediction reverses only the class2-oriented direction", {
  z <- c5lin_fixture()
  fit <- dsda_classifier(z$x, z$y, lambda = 0.15)
  beta <- fit$estimate$regression.direction
  midpoint <- (
    fit$estimate$locations$class1 + fit$estimate$locations$class2
  ) / 2
  expected <- -drop(sweep(z$x, 2L, midpoint, "-") %*% beta)

  expect_equal(fit$estimate$direction, -beta, tolerance = 0)
  expect_equal(predict(fit, z$x, type = "score"),
               expected, tolerance = 2e-13)
  expect_true(all(predict(fit, z$x)[1:6] == "class-a"))
  expect_true(all(predict(fit, z$x)[7:12] == "class-b"))
  expect_equal(
    fit$estimate$regression.intercept,
    -sum(colMeans(z$x) * beta), tolerance = 2e-13
  )
  expect_equal(
    fit$estimate$classification.intercept,
    -sum(midpoint * fit$estimate$direction), tolerance = 2e-13
  )
  expect_match(fit$diagnostics$classification.intercept,
               "not used for prediction")
})


test_that("DSDA equal-prior and zero-direction boundaries are explicit", {
  z <- c5lin_fixture()
  expect_error(
    dsda_classifier(z$x, z$y, lambda = 0.2,
                    prior = c(0.7, 0.3)),
    "equal", ignore.case = TRUE
  )
  zero <- dsda_classifier(z$x, z$y, lambda = 100)
  expect_equal(unname(zero$estimate$direction), rep(0, 3), tolerance = 0)
  expect_equal(predict(zero, z$x, type = "score"),
               rep(0, 12), tolerance = 0)
  expect_true(all(predict(zero, z$x) == "class-a"))

  zero.class2 <- dsda_classifier(
    z$x, z$y, lambda = 100, tie = "class2"
  )
  expect_true(all(predict(zero.class2, z$x) == "class-b"))
})


test_that("DSDA refuses a non-certified coordinate-descent result", {
  z <- c5lin_fixture()
  invalid <- NULL
  expect_warning(
    invalid <- dsda_classifier(
      z$x, z$y, lambda = 0.01, solver_tol = 1e-14,
      solver_max_iter = 1L, strict = FALSE
    ),
    "coordinate-descent certificate"
  )
  expect_false(invalid$valid)
  expect_length(invalid$estimate, 0L)
  expect_false(is.null(invalid$diagnostics$solver$last.iterate))
})


test_that("SSLDA classwise SSCMs and explicit ridge operator cross-lock", {
  z <- c5lin_fixture()
  fit <- sslda(
    z$x, z$y, lambda = 0.25, ridge = 0.2,
    median_tol = 1e-10, solver_tol = 1e-8
  )
  location1 <- fit$estimate$locations$class1
  location2 <- fit$estimate$locations$class2
  signs1 <- c5lin_manual_signs(z$first, location1)
  signs2 <- c5lin_manual_signs(z$second, location2)
  sscm1 <- crossprod(signs1) / 6
  sscm2 <- crossprod(signs2) / 6
  pooled <- (6 * sscm1 + 6 * sscm2) / 12
  operator <- 3 * pooled
  diag(operator) <- diag(operator) + 0.2
  direction <- fit$estimate$direction

  expect_true(fit$valid)
  expect_equal(unname(fit$estimate$classwise.sscm$class1),
               sscm1, tolerance = 2e-12)
  expect_equal(unname(fit$estimate$classwise.sscm$class2),
               sscm2, tolerance = 2e-12)
  expect_equal(unname(fit$estimate$pooled.sscm), pooled,
               tolerance = 2e-12)
  expect_equal(unname(fit$estimate$constraint.operator), operator,
               tolerance = 2e-12)
  expect_equal(
    max(abs(drop(operator %*% direction) -
              (location1 - location2))),
    0.25, tolerance = 2e-7
  )
  expect_true(fit$diagnostics$solver$certified)
  expect_true(fit$tuning$ridge.changes.operator)
  expect_identical(fit$diagnostics$sscm.pooling,
                   "(n1*S1+n2*S2)/(n1+n2)")
})


test_that("SSLDA has its translation and positive-scale decision invariance", {
  z <- c5lin_fixture()
  shift <- c(11, -7, 5)
  base <- sslda(z$x, z$y, lambda = 0.25, ridge = 0.2)
  moved.x <- 4 * sweep(z$x, 2L, shift, "+")
  moved <- sslda(
    moved.x, z$y, lambda = 4 * 0.25, ridge = 0.2
  )

  expect_identical(
    as.character(predict(moved, moved.x)),
    as.character(predict(base, z$x))
  )
  expect_equal(
    sign(predict(moved, moved.x, type = "score")),
    sign(predict(base, z$x, type = "score")),
    tolerance = 0
  )
  expect_equal(
    moved$estimate$locations$class1,
    4 * (base$estimate$locations$class1 + shift),
    tolerance = 2e-7
  )
})


test_that("SSLDA enforces equal priors and spatial-median certificates", {
  z <- c5lin_fixture()
  expect_error(
    sslda(z$x, z$y, lambda = 0.25,
          prior = c(0.6, 0.4)),
    "equal", ignore.case = TRUE
  )

  hard <- z$x
  hard[1, ] <- c(20, -10, 7)
  invalid <- NULL
  expect_warning(
    invalid <- sslda(
      hard, z$y, lambda = 0.25, median_tol = 1e-15,
      median_max_iter = 1L, strict = FALSE
    ),
    "spatial medians"
  )
  expect_false(invalid$valid)
  expect_identical(invalid$diagnostics$failure.stage,
                   "classwise spatial medians")
})


test_that("sparse plug-in LDA matches supplied precision and covariance", {
  z <- c5lin_fixture()
  precision <- diag(c(2, 4, 8))
  dimnames(precision) <- list(colnames(z$x), colnames(z$x))
  from.precision <- sparse_plugin_lda(
    z$x, z$y, precision = precision, prior = c(0.75, 0.25)
  )
  covariance <- diag(c(0.5, 0.25, 0.125))
  dimnames(covariance) <- list(colnames(z$x), colnames(z$x))
  from.covariance <- sparse_plugin_lda(
    z$x, z$y, covariance = covariance,
    prior = c(0.75, 0.25)
  )
  delta <- colMeans(z$first) - colMeans(z$second)
  expected.direction <- drop(precision %*% delta)

  expect_equal(from.precision$estimate$direction,
               expected.direction, tolerance = 0)
  expect_equal(from.covariance$estimate$direction,
               expected.direction, tolerance = 2e-14)
  expect_equal(from.covariance$estimate$precision,
               precision, tolerance = 2e-14)
  expect_equal(
    from.precision$score.model$intercept,
    -sum((colMeans(z$first) + colMeans(z$second)) / 2 *
           expected.direction) + log(3),
    tolerance = 1e-14
  )
  expect_identical(from.precision$diagnostics$inversion, "none")
  expect_identical(from.covariance$diagnostics$inversion,
                   "ordinary Cholesky inverse")
})


test_that("sparse plug-in LDA never repairs invalid matrix inputs", {
  z <- c5lin_fixture()
  expect_error(sparse_plugin_lda(z$x, z$y), "exactly one")
  expect_error(
    sparse_plugin_lda(
      z$x, z$y, precision = diag(3), covariance = diag(3)
    ),
    "exactly one"
  )
  nonsymmetric <- diag(3)
  nonsymmetric[1, 2] <- 0.2
  expect_error(
    sparse_plugin_lda(z$x, z$y, precision = nonsymmetric),
    "symmetric"
  )
  singular <- diag(c(1, 1, 0))
  expect_error(
    sparse_plugin_lda(z$x, z$y, covariance = singular),
    "positive definite"
  )
  indefinite <- diag(c(1, 1, -1))
  expect_error(
    sparse_plugin_lda(z$x, z$y, precision = indefinite),
    "positive definite"
  )
  wrong.names <- diag(3)
  dimnames(wrong.names) <- list(c("b", "a", "c"),
                                c("b", "a", "c"))
  expect_error(
    sparse_plugin_lda(z$x, z$y, precision = wrong.names),
    "match the training features"
  )
})


test_that("sparse plug-in LDA accepts explicit locations without alteration", {
  z <- c5lin_fixture()
  locations <- list(
    c(a = 2, b = 1, c = 0.5),
    c(a = -1, b = -2, c = -0.5)
  )
  fit <- sparse_plugin_lda(
    z$x, z$y, precision = diag(3), locations = locations
  )
  expect_equal(fit$estimate$locations$class1,
               locations[[1]], tolerance = 0)
  expect_equal(fit$estimate$locations$class2,
               locations[[2]], tolerance = 0)
  expect_equal(fit$estimate$direction,
               locations[[1]] - locations[[2]], tolerance = 0)
  expect_identical(fit$estimate$locations$source,
                   "supplied classwise locations")
})


c5lin_certified_precision <- function(p = 3L,
                                      method = "sglasso",
                                      certified = TRUE) {
  solver <- if (identical(method, "sclime")) {
    list(all.columns.certified = certified)
  } else {
    list(kkt.certified = certified)
  }
  structure(
    list(
      estimate = diag(p),
      valid = TRUE,
      method = method,
      lambda = 0.3,
      diagnostics = list(solver = solver)
    ),
    class = "spatial_sign_precision_fit"
  )
}


test_that("precision plug-in LDA consumes only a certified precision fit", {
  z <- c5lin_fixture()
  locations <- list(colMeans(z$first), colMeans(z$second))
  fit <- spatial_sign_precision_lda(
    z$x, z$y, c5lin_certified_precision(),
    locations = locations
  )
  delta <- locations[[1]] - locations[[2]]

  expect_true(fit$valid)
  expect_equal(unname(fit$estimate$direction), delta, tolerance = 0)
  expect_identical(fit$diagnostics$precision.certificate,
                   "SGLASSO KKT certified")
  expect_identical(fit$estimate$precision.method, "sglasso")
  expect_match(fit$score.scale, "scale-free")

  sclime <- spatial_sign_precision_lda(
    z$x, z$y, c5lin_certified_precision(method = "sclime"),
    locations = locations
  )
  expect_identical(sclime$diagnostics$precision.certificate,
                   "all SCLIME columns certified")
})


test_that("precision plug-in LDA rejects plain, invalid, and uncertified fits", {
  z <- c5lin_fixture()
  locations <- list(colMeans(z$first), colMeans(z$second))
  expect_error(
    spatial_sign_precision_lda(
      z$x, z$y, diag(3), locations = locations
    ),
    "valid SCLIME or SGLASSO"
  )
  invalid <- c5lin_certified_precision()
  invalid$valid <- FALSE
  expect_error(
    spatial_sign_precision_lda(
      z$x, z$y, invalid, locations = locations
    ),
    "valid SCLIME or SGLASSO"
  )
  uncertified <- c5lin_certified_precision(certified = FALSE)
  expect_error(
    spatial_sign_precision_lda(
      z$x, z$y, uncertified, locations = locations
    ),
    "certificate"
  )
  wrong.method <- c5lin_certified_precision(method = "other")
  expect_error(
    spatial_sign_precision_lda(
      z$x, z$y, wrong.method, locations = locations
    ),
    "valid SCLIME or SGLASSO"
  )
  expect_error(
    spatial_sign_precision_lda(
      z$x, z$y, c5lin_certified_precision(),
      locations = locations, prior = c(0.6, 0.4)
    ),
    "equal", ignore.case = TRUE
  )
})


test_that("precision plug-in LDA computes classwise spatial medians by default", {
  z <- c5lin_fixture()
  fit <- spatial_sign_precision_lda(
    z$x, z$y, c5lin_certified_precision(),
    median_tol = 1e-10
  )
  expected1 <- spatial_median(z$first, tol = 1e-10)
  expected2 <- spatial_median(z$second, tol = 1e-10)

  expect_equal(unname(fit$estimate$locations$class1),
               as.numeric(expected1), tolerance = 2e-8)
  expect_equal(unname(fit$estimate$locations$class2),
               as.numeric(expected2), tolerance = 2e-8)
  expect_identical(fit$diagnostics$location.source,
                   "classwise sample spatial medians")
  expect_true(all(vapply(
    fit$diagnostics$spatial.medians,
    function(item) isTRUE(item$converged),
    logical(1)
  )))
})


test_that("linear classifiers validate tuning and response inputs", {
  z <- c5lin_fixture()
  expect_error(lpd_classifier(z$x, z$y, 0), "greater than")
  expect_error(lpd_classifier(z$x, z$y, NA_real_), "finite")
  expect_error(lpd_classifier(z$x, z$y, 0.2, ridge = -1),
               "at least")
  expect_error(lpd_classifier(z$x, z$y, 0.2, solver_tol = 0),
               "greater than")
  expect_error(lpd_classifier(z$x, z$y, 0.2,
                              solver_max_iter = 1.5),
               "positive integer")
  expect_error(lpd_classifier(z$x, z$y, 0.2, strict = NA),
               "TRUE or FALSE")
  expect_error(lpd_classifier(z$x, z$y, 0.2, tie = "bad"),
               "arg")
  expect_error(
    lpd_classifier(z$x, factor(rep("one", 12)), 0.2),
    "two|binary", ignore.case = TRUE
  )
  bad <- z$x
  bad[1, 1] <- NA_real_
  expect_error(lpd_classifier(bad, z$y, 0.2),
               "finite|missing", ignore.case = TRUE)
})


test_that("all public linear fits expose one shared score contract", {
  z <- c5lin_fixture()
  locations <- list(colMeans(z$first), colMeans(z$second))
  fits <- list(
    lpd_classifier(z$x, z$y, 0.2, ridge = 0.1),
    dsda_classifier(z$x, z$y, 0.15),
    sslda(z$x, z$y, 0.25, ridge = 0.2),
    sparse_plugin_lda(z$x, z$y, precision = diag(3)),
    spatial_sign_precision_lda(
      z$x, z$y, c5lin_certified_precision(),
      locations = locations
    )
  )

  for (fit in fits) {
    expect_s3_class(fit, "hd_classifier_fit")
    expect_true(fit$valid)
    expect_identical(fit$score.model$type, "linear")
    expect_length(fit$score.model$coefficients, 3L)
    expect_true(is.finite(fit$score.model$intercept))
    expect_identical(fit$levels, levels(z$y))
    expect_identical(fit$feature.names, colnames(z$x))
    expect_equal(
      predict(fit, z$x, type = "score"),
      drop(z$x %*% fit$score.model$coefficients) +
        fit$score.model$intercept,
      tolerance = 2e-12
    )
    expect_identical(fit$primary.orientation,
                     "positive score selects class1")
    expect_true(is.character(fit$score.scale))
    expect_true(nzchar(fit$score.scale))
  }
})
