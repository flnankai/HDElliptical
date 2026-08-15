c5_gqda_fixture <- function() {
  x1 <- rbind(
    c(-2.4, -0.7), c(-2.0, 0.4), c(-1.7, 1.1), c(-1.2, -1.0),
    c(-2.8, 0.8), c(-1.4, 0.1), c(-2.2, -1.3), c(-0.9, 0.9)
  )
  x2 <- rbind(
    c(2.6, -0.4), c(1.9, 0.7), c(1.2, 1.5), c(2.3, -1.2),
    c(3.1, 1.0), c(1.5, -0.2), c(2.7, -1.5), c(0.8, 0.8)
  )
  x <- rbind(x1, x2)
  colnames(x) <- c("first", "second")
  rownames(x) <- paste0("g", seq_len(nrow(x)))
  newdata <- rbind(c(-1.3, 0.2), c(0, 0), c(2.1, -0.1))
  colnames(newdata) <- c("first", "second")
  list(
    x = x,
    y = factor(rep(c("one", "two"), each = 8L),
               levels = c("one", "two")),
    newdata = newdata
  )
}


c5_gqda_class_fits <- function(x, y) {
  lapply(levels(y), function(level) {
    z <- x[y == level, , drop = FALSE]
    list(
      valid = TRUE,
      location = colMeans(z),
      scatter = stats::cov(z)
    )
  })
}


c5_gqda_manual_score <- function(newdata, class.fits, c) {
  scatter <- lapply(class.fits, `[[`, "scatter")
  precision <- lapply(scatter, solve)
  logdet <- vapply(scatter, function(value) {
    determinant(value, logarithm = TRUE)$modulus[[1L]]
  }, numeric(1))
  squared <- lapply(seq_len(2L), function(group) {
    residual <- sweep(newdata, 2L, class.fits[[group]]$location, "-")
    rowSums((residual %*% precision[[group]]) * residual)
  })
  squared[[2L]] - squared[[1L]] - c * (logdet[[1L]] - logdet[[2L]])
}


test_that("GQDA fixed scores equal the literal primary inequality", {
  fixture <- c5_gqda_fixture()
  class.fits <- c5_gqda_class_fits(fixture$x, fixture$y)
  for (constant in c(0, 0.37, 1)) {
    fit <- gqda_classifier(fixture$x, fixture$y, c = constant)
    observed <- predict(fit, fixture$newdata, type = "score")
    expected <- c5_gqda_manual_score(fixture$newdata, class.fits, constant)
    expect_equal(unname(observed), expected, tolerance = 2e-12)
    expect_equal(fit$c, constant)
    expect_identical(fit$tuning$selection, "fixed")
  }
})


test_that("c=1 is twice the equal-prior canonical Gaussian log likelihood ratio", {
  fixture <- c5_gqda_fixture()
  fit <- gqda_classifier(fixture$x, fixture$y, c = 1)
  class.fits <- c5_gqda_class_fits(fixture$x, fixture$y)
  score <- predict(fit, fixture$newdata, type = "score")
  canonical <- c5_gqda_manual_score(fixture$newdata, class.fits, 1) / 2
  expect_equal(unname(score), 2 * canonical, tolerance = 2e-12)
  expect_identical(fit$score.scale, "method_threshold")
  expect_match(fit$diagnostics$primary.rule, "class1 iff")
})


test_that("breakpoint enumeration handles signed and zero log determinants", {
  difference <- c(0.15, 0.72, -0.18, -0.83, 0.05, -0.02)
  logdet <- c(1, 1, -1, -1, 0, 0)
  truth <- c(1L, 2L, 1L, 2L, 1L, 2L)
  selected <- HDElliptical:::.c5gqda_select_c(difference, logdet, truth)
  brute.errors <- vapply(selected$candidates, function(constant) {
    sum(ifelse(difference >= constant * logdet, 1L, 2L) != truth)
  }, numeric(1))
  expect_equal(selected$errors, brute.errors)
  expect_equal(selected$c, selected$candidates[which.min(brute.errors)])
  expect_true(any(selected$breakpoints > 0 & selected$breakpoints < 1))

  fixture <- c5_gqda_fixture()
  equal.det <- list(
    list(valid = TRUE, location = c(-1, 0), scatter = diag(c(2, 0.5))),
    list(valid = TRUE, location = c(1, 0), scatter = diag(c(2, 0.5)))
  )
  fit <- robust_gqda(fixture$x, fixture$y, equal.det)
  expect_equal(fit$log.determinant.contrast, 0, tolerance = 1e-14)
  expect_false(fit$tuning$c.identified)
  expect_equal(fit$c, 0)
})


test_that("the primary weak inequality sends exact boundary ties to class 1", {
  fixture <- c5_gqda_fixture()
  same <- list(
    valid = TRUE,
    location = c(first = 0, second = 0),
    scatter = diag(2)
  )
  fit <- robust_gqda(fixture$x, fixture$y, list(same, same), c = 0.8)
  boundary <- matrix(c(0, 0, 2, -3), 2, 2, byrow = TRUE,
                     dimnames = list(NULL, c("first", "second")))
  expect_equal(predict(fit, boundary, type = "score"), c(0, 0))
  expect_identical(as.character(predict(fit, boundary)), c("one", "one"))
})


test_that("certified robust fits reproduce classical GQDA exactly", {
  fixture <- c5_gqda_fixture()
  class.fits <- c5_gqda_class_fits(fixture$x, fixture$y)
  classical <- gqda_classifier(fixture$x, fixture$y, c = 0.41)
  robust <- robust_gqda(fixture$x, fixture$y, class.fits, c = 0.41)
  expect_equal(
    predict(robust, fixture$newdata, type = "score"),
    predict(classical, fixture$newdata, type = "score"),
    tolerance = 2e-12
  )
  expect_match(robust$method, "certified robust")
})


test_that("GQDA is equivariant under class exchange and coordinate transforms", {
  fixture <- c5_gqda_fixture()
  constant <- 0.63
  original <- gqda_classifier(fixture$x, fixture$y, c = constant)
  score <- predict(original, fixture$newdata, type = "score")

  reversed.y <- factor(fixture$y, levels = rev(levels(fixture$y)))
  reversed <- gqda_classifier(fixture$x, reversed.y, c = constant)
  expect_equal(
    predict(reversed, fixture$newdata, type = "score"),
    -score, tolerance = 3e-12
  )

  transform <- c(first = -3.5, second = 0.4)
  shift <- c(first = 17, second = -9)
  transformed.x <- sweep(sweep(fixture$x, 2L, transform, "*"),
                         2L, shift, "+")
  transformed.new <- sweep(sweep(fixture$newdata, 2L, transform, "*"),
                           2L, shift, "+")
  transformed <- gqda_classifier(transformed.x, fixture$y, c = constant)
  expect_equal(
    unname(predict(transformed, transformed.new, type = "score")),
    unname(score), tolerance = 2e-11
  )
  expect_equal(
    predict(original, fixture$newdata[, 2:1, drop = FALSE], type = "score"),
    score, tolerance = 2e-12
  )
})


test_that("finite extreme common scaling preserves the audited score", {
  fixture <- c5_gqda_fixture()
  baseline <- gqda_classifier(fixture$x, fixture$y, c = 0.2)
  baseline.score <- predict(baseline, fixture$newdata, type = "score")
  for (scale in c(1e-100, 1e100)) {
    fit <- gqda_classifier(fixture$x * scale, fixture$y, c = 0.2)
    value <- predict(fit, fixture$newdata * scale, type = "score")
    expect_equal(unname(value), unname(baseline.score), tolerance = 2e-11)
  }
})


test_that("classical cross-validation is stratified, leakage-free, and RNG local", {
  fixture <- c5_gqda_fixture()
  set.seed(815)
  before <- .Random.seed
  fit <- gqda_classifier(
    fixture$x, fixture$y, selection = "cross_validation",
    folds = 4, seed = 291
  )
  expect_identical(.Random.seed, before)
  assignment <- fit$tuning$folds$assignment
  expect_equal(tabulate(assignment[fixture$y == "one"], 4), rep(2L, 4))
  expect_equal(tabulate(assignment[fixture$y == "two"], 4), rep(2L, 4))
  for (fold in seq_len(4L)) {
    diagnostic <- fit$diagnostics$cross.validation[[fold]]
    expect_length(intersect(diagnostic$validation, diagnostic$fitting), 0L)
  }
  repeat.fit <- gqda_classifier(
    fixture$x, fixture$y, selection = "cross_validation",
    folds = 4, seed = 291
  )
  expect_equal(fit$c, repeat.fit$c)
  expect_identical(fit$tuning$path$errors,
                   repeat.fit$tuning$path$errors)
})


test_that("invalid certificates fail without numerical repair", {
  fixture <- c5_gqda_fixture()
  good <- c5_gqda_class_fits(fixture$x, fixture$y)
  singular <- good
  singular[[1L]]$scatter <- diag(c(1, 0))
  expect_error(
    robust_gqda(fixture$x, fixture$y, singular, c = 0),
    "positive definite"
  )
  expect_warning(
    invalid <- robust_gqda(
      fixture$x, fixture$y, singular, c = 0, strict = FALSE
    ),
    "positive definite"
  )
  expect_false(invalid$valid)
  expect_match(invalid$diagnostics$no.repair, "No ridge")
  expect_error(predict(invalid, fixture$newdata), "invalid")

  inconsistent <- good
  inconsistent[[1L]]$precision <- diag(2)
  expect_error(
    robust_gqda(fixture$x, fixture$y, inconsistent, c = 0),
    "inconsistent scatter and precision"
  )
  expect_error(
    gqda_classifier(fixture$x, fixture$y, c = 0,
                    selection = "resubstitution"),
    "incompatible"
  )
})


test_that("HR-GQDA uses certified paper-QDA scales and alias semantics", {
  angle <- 2 * pi * (0:7) / 8
  cloud <- cbind(cos(angle), sin(angle))
  x <- rbind(
    sweep(cloud, 2L, c(-2, 0), "+"),
    sweep(sweep(cloud, 2L, c(1.4, 0.7), "*"), 2L, c(2, 0), "+")
  )
  colnames(x) <- c("first", "second")
  y <- factor(rep(c("left", "right"), each = 8L),
              levels = c("left", "right"))
  fit <- hr_gqda(
    x, y, pilot_precision = diag(2), bandwidth = 1L, c = 0,
    tol = 1e-7, max_iter = 2000L,
    median_tol = 1e-8, median_max_iter = 2000L
  )
  expect_true(fit$valid)
  expect_true(all(vapply(fit$hr.fits, inherits, logical(1),
                         what = "high_dimensional_hr_fit")))
  expect_true(all(vapply(fit$hr.fits, function(value) {
    identical(value$scale.estimator, "paper_qda")
  }, logical(1))))
  expect_true(all(is.finite(predict(fit, x, type = "score"))))

  alias <- hr_qda(
    x, y, pilot_precision = diag(2), bandwidth = 1L, c = 0,
    tol = 1e-7, max_iter = 2000L,
    median_tol = 1e-8, median_max_iter = 2000L
  )
  expect_equal(predict(alias, x, type = "score"),
               predict(fit, x, type = "score"), tolerance = 2e-11)
  expect_match(alias$diagnostics$alias, "generalized QDA")
})

