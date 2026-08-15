ch5_classical_fixture <- function() {
  x <- rbind(
    c(2.0, 0.2, 1.0),
    c(1.0, 1.4, 0.2),
    c(3.0, -0.5, 1.7),
    c(2.5, 0.8, -0.3),
    c(0.5, -1.0, 0.8),
    c(-2.0, 0.0, -0.7),
    c(-1.0, -1.2, 0.4),
    c(-3.0, 0.7, -1.5),
    c(-2.5, -0.8, 0.6),
    c(-0.5, 1.1, -0.1)
  )
  colnames(x) <- c("signal", "aux", "third")
  rownames(x) <- paste0("r", seq_len(nrow(x)))
  list(x = x, y = factor(rep(c("A", "B"), each = 5),
                         levels = c("A", "B")))
}


ch5_manual_moments <- function(x, y) {
  levels <- levels(y)
  x1 <- x[y == levels[1], , drop = FALSE]
  x2 <- x[y == levels[2], , drop = FALSE]
  mean1 <- colMeans(x1)
  mean2 <- colMeans(x2)
  centered1 <- sweep(x1, 2L, mean1, "-")
  centered2 <- sweep(x2, 2L, mean2, "-")
  list(
    mean1 = mean1, mean2 = mean2,
    sscp1 = crossprod(centered1), sscp2 = crossprod(centered2),
    n1 = nrow(x1), n2 = nrow(x2)
  )
}


test_that("shared classifier helper formals and fields are stable", {
  expect_identical(
    names(formals(HDElliptical:::.clf_prepare_xy)),
    c("x", "y", "prior", "equal_prior", "min_class")
  )
  expect_identical(
    names(formals(HDElliptical:::.clf_validate_newdata)),
    c("newdata", "object")
  )
  expect_identical(
    names(formals(HDElliptical:::.clf_new_fit)),
    c(
      "method", "training", "levels", "feature_names", "n_features",
      "score_model", "estimate", "tuning", "diagnostics", "call",
      "primary_orientation", "score_scale", "tie", "data_name", "valid"
    )
  )

  fixture <- ch5_classical_fixture()
  training <- HDElliptical:::.clf_prepare_xy(
    fixture$x, fixture$y, prior = "empirical"
  )
  expect_identical(training$class1, 1:5)
  expect_identical(training$class2, 6:10)
  expect_identical(training$levels, c("A", "B"))
  expect_identical(training$feature.names, colnames(fixture$x))
  expect_identical(training$observation.names, rownames(fixture$x))
  expect_equal(training$prior$probabilities, c(A = 0.5, B = 0.5))

  fit <- classical_lda_classifier(fixture$x, fixture$y)
  expect_identical(
    names(fit),
    c(
      "valid", "method", "levels", "feature.names", "p", "prior",
      "score.model", "estimate", "tuning", "diagnostics", "call",
      "primary.orientation", "score.scale", "data.name"
    )
  )
  expect_identical(fit$score.model$type, "linear")
  expect_identical(fit$score.model$tie, "class1")
})


test_that("training labels and priors have deterministic contracts", {
  x <- matrix(seq_len(16), 8, 2)
  y <- c("second", "second", "second", "second",
         "first", "first", "first", "first")
  prepared <- HDElliptical:::.clf_prepare_xy(
    x, y, prior = c(first = 1, second = 3)
  )
  expect_identical(prepared$levels, c("second", "first"))
  expect_equal(prepared$prior$probabilities,
               c(second = 0.75, first = 0.25))
  expect_equal(unname(prepared$prior$log.ratio), log(3))

  y.factor <- factor(y, levels = c("first", "second", "unused"))
  prepared.factor <- HDElliptical:::.clf_prepare_xy(x, y.factor)
  expect_identical(prepared.factor$levels, c("first", "second"))
  expect_error(
    HDElliptical:::.clf_prepare_xy(x, y, prior = c(first = 1, wrong = 1)),
    "match"
  )
  expect_error(
    HDElliptical:::.clf_prepare_xy(x, y, prior = c(0.7, 0.3),
                                  equal_prior = TRUE),
    "equal-prior"
  )
  expect_error(
    HDElliptical:::.clf_prepare_xy(x, rep("one", 8)),
    "Exactly two"
  )
})


test_that("low-level fit and prediction support custom, linear, and quadratic scores", {
  custom <- hd_classifier_fit(
    function(z) z[, 1] - z[, 2], 2,
    levels = c("yes", "no"), feature_names = c("a", "b")
  )
  z <- matrix(c(0, 0, 2, 1, 0, 1), 3, 2, byrow = TRUE,
              dimnames = list(c("tie", "positive", "negative"), c("a", "b")))
  expect_equal(predict(custom, z, type = "score"), c(tie = 0, positive = 1,
                                                       negative = -1))
  expect_identical(as.character(predict(custom, z)), c("yes", "yes", "no"))
  expect_error(predict(custom, z[, 1, drop = FALSE]), "exactly 2")

  reversed <- z[, c("b", "a"), drop = FALSE]
  expect_equal(predict(custom, reversed, type = "score"),
               predict(custom, z, type = "score"))
  expect_error(
    predict(custom, setNames(as.data.frame(z), c("a", "wrong"))),
    "match"
  )

  quadratic <- HDElliptical:::.clf_new_fit(
    method = "quadratic check", levels = c("one", "two"), n_features = 2L,
    score_model = list(
      type = "quadratic", quadratic = diag(c(1, -1)),
      linear = c(2, 3), intercept = -4
    ),
    score_scale = "method_threshold"
  )
  expected <- rowSums((z %*% diag(c(1, -1))) * z) +
    as.numeric(z %*% c(2, 3)) - 4
  expect_equal(
    unname(predict(quadratic, unname(z), type = "score")),
    unname(expected)
  )

  tie.second <- hd_classifier_fit(
    function(z) numeric(nrow(z)), 2, tie = "class2"
  )
  expect_identical(as.character(predict(tie.second, matrix(0, 1, 2))),
                   "class2")
})


test_that("Gaussian oracle LDA equals its defining log likelihood ratio", {
  mu1 <- c(x = 1, y = -0.5)
  mu2 <- c(x = -0.5, y = 0.25)
  sigma <- matrix(c(2, 0.4, 0.4, 1.2), 2)
  prior <- c(class1 = 0.7, class2 = 0.3)
  z <- rbind(a = c(0, 0), b = c(2, -1), c = c(-1, 1))
  fit <- gaussian_lda_oracle(
    mu1, mu2, covariance = sigma, prior = prior,
    feature_names = c("x", "y")
  )
  omega <- solve(sigma)
  expected <- as.numeric((z - rep((mu1 + mu2) / 2, each = nrow(z))) %*%
                           omega %*% (mu1 - mu2)) + log(0.7 / 0.3)
  expect_equal(unname(predict(fit, z, type = "score")), expected,
               tolerance = 1e-12)
  fit.precision <- gaussian_lda_oracle(mu1, mu2, precision = omega,
                                       prior = prior)
  expect_equal(unname(predict(fit.precision, unname(z), type = "score")),
               expected, tolerance = 1e-12)
  expect_identical(fit$score.scale, "canonical_log_lr")
  expect_error(
    gaussian_lda_oracle(mu1, mu2, covariance = sigma, precision = omega),
    "exactly one"
  )
  expect_error(
    gaussian_lda_oracle(mu1, mu2, covariance = matrix(c(1, 2, 2, 1), 2)),
    "positive definite"
  )
})


test_that("Gaussian oracle QDA is canonical rather than the twice-scaled display", {
  mu1 <- c(0.5, -0.2)
  mu2 <- c(-0.4, 0.6)
  sigma1 <- matrix(c(1.2, 0.2, 0.2, 0.8), 2)
  sigma2 <- matrix(c(0.7, -0.1, -0.1, 1.5), 2)
  z <- rbind(c(0, 0), c(1, -1), c(-1, 2))
  fit <- gaussian_qda_oracle(
    mu1, mu2, covariance1 = sigma1, covariance2 = sigma2,
    prior = c(0.4, 0.6)
  )
  omega1 <- solve(sigma1)
  omega2 <- solve(sigma2)
  d1 <- rowSums(((z - rep(mu1, each = nrow(z))) %*% omega1) *
                  (z - rep(mu1, each = nrow(z))))
  d2 <- rowSums(((z - rep(mu2, each = nrow(z))) %*% omega2) *
                  (z - rep(mu2, each = nrow(z))))
  expected <- -determinant(sigma1, logarithm = TRUE)$modulus / 2 +
    determinant(sigma2, logarithm = TRUE)$modulus / 2 - d1 / 2 + d2 / 2 +
    log(0.4 / 0.6)
  expect_equal(unname(predict(fit, z, type = "score")), as.numeric(expected),
               tolerance = 1e-12)
  expect_match(fit$diagnostics$book.errata, "twice")
  expect_identical(fit$score.scale, "canonical_log_lr")
})


test_that("elliptical oracle evaluates the supplied generators exactly", {
  mu1 <- c(1, 0)
  mu2 <- c(-1, 0)
  shape <- diag(2)
  gaussian <- function(d) -d / 2
  z <- rbind(c(0.2, 0), c(0.2, 3), c(-0.5, 1))
  elliptical <- elliptical_oracle_classifier(
    mu1, mu2, shape, shape, gaussian, prior = c(0.7, 0.3)
  )
  qda <- gaussian_qda_oracle(
    mu1, mu2, covariance1 = shape, covariance2 = shape,
    prior = c(0.7, 0.3)
  )
  expect_equal(predict(elliptical, z, type = "score"),
               predict(qda, z, type = "score"), tolerance = 1e-12)

  nu <- 4
  student <- function(d) -(nu + 2) * log1p(d / nu) / 2
  t.fit <- elliptical_oracle_classifier(
    mu1, mu2, shape, shape, student, prior = c(0.7, 0.3)
  )
  t.score <- predict(t.fit, z[1:2, , drop = FALSE], type = "score")
  expect_false(isTRUE(all.equal(unname(t.score[1]), unname(t.score[2]))))
  expect_match(t.fit$diagnostics$book.errata, "equal priors")

  bad <- elliptical_oracle_classifier(
    mu1, mu2, shape, shape, function(d) rep(NA_real_, length(d))
  )
  expect_error(predict(bad, z), "log-generator")
})


test_that("classical LDA uses pooled divisor n minus two and is translation equivariant", {
  fixture <- ch5_classical_fixture()
  moments <- ch5_manual_moments(fixture$x, fixture$y)
  covariance <- (moments$sscp1 + moments$sscp2) /
    (nrow(fixture$x) - 2)
  direction <- solve(covariance, moments$mean1 - moments$mean2)
  midpoint <- (moments$mean1 + moments$mean2) / 2
  expected <- as.numeric(fixture$x %*% direction - sum(midpoint * direction))
  fit <- classical_lda_classifier(fixture$x, fixture$y)
  expect_equal(predict(fit, fixture$x, type = "score"), expected,
               tolerance = 1e-11, ignore_attr = TRUE)
  expect_equal(fit$estimate$covariance, covariance, tolerance = 1e-12)
  expect_identical(fit$diagnostics$covariance.divisor,
                   nrow(fixture$x) - 2L)

  shift <- c(10, -3, 2)
  shifted <- sweep(fixture$x, 2L, shift, "+")
  shifted.fit <- classical_lda_classifier(shifted, fixture$y)
  expect_equal(
    predict(shifted.fit, shifted, type = "score"),
    predict(fit, fixture$x, type = "score"), tolerance = 1e-10
  )
})


test_that("classical QDA uses class divisors n_k minus one", {
  fixture <- ch5_classical_fixture()
  moments <- ch5_manual_moments(fixture$x, fixture$y)
  sigma1 <- moments$sscp1 / (moments$n1 - 1)
  sigma2 <- moments$sscp2 / (moments$n2 - 1)
  fit <- classical_qda_classifier(fixture$x, fixture$y)
  expect_equal(fit$estimate$covariance1, sigma1, tolerance = 1e-12)
  expect_equal(fit$estimate$covariance2, sigma2, tolerance = 1e-12)
  expect_equal(
    unname(fit$diagnostics$covariance.divisor),
    c(moments$n1 - 1, moments$n2 - 1)
  )
  z <- fixture$x[c(2, 8), , drop = FALSE]
  manual <- function(row) {
    d1 <- row - moments$mean1
    d2 <- row - moments$mean2
    -as.numeric(determinant(sigma1, logarithm = TRUE)$modulus) / 2 +
      as.numeric(determinant(sigma2, logarithm = TRUE)$modulus) / 2 -
      sum(d1 * solve(sigma1, d1)) / 2 +
      sum(d2 * solve(sigma2, d2)) / 2
  }
  expect_equal(
    unname(predict(fit, z, type = "score")),
    unname(apply(z, 1L, manual)), tolerance = 1e-10
  )
})


test_that("strict covariance failures never trigger hidden repair", {
  x <- cbind(rep(c(1, 2, -1, -2), each = 2),
             rep(c(1, 2, -1, -2), each = 2))
  y <- rep(c("A", "B"), each = 4)
  expect_error(classical_lda_classifier(x, y), "positive-definite")
  expect_warning(
    invalid <- classical_lda_classifier(x, y, strict = FALSE),
    "No ridge"
  )
  expect_false(invalid$valid)
  expect_null(invalid$estimate)
  expect_match(invalid$diagnostics$no.repair, "No ridge")
  expect_error(predict(invalid, matrix(0, 1, 2)), "invalid")
  expect_error(classical_qda_classifier(x, y), "positive definite")
})


test_that("independence rule exposes its marginal variance convention", {
  fixture <- ch5_classical_fixture()
  moments <- ch5_manual_moments(fixture$x, fixture$y)
  delta <- moments$mean1 - moments$mean2
  variance.unbiased <- diag(moments$sscp1 + moments$sscp2) /
    (nrow(fixture$x) - 2)
  fit <- independence_classifier(fixture$x, fixture$y)
  expect_equal(fit$estimate$marginal.variance, variance.unbiased)
  expect_equal(fit$estimate$direction, delta / variance.unbiased)
  expect_identical(fit$diagnostics$covariance.divisor,
                   nrow(fixture$x) - 2L)

  mle <- independence_classifier(
    fixture$x, fixture$y, variance_divisor = "mle"
  )
  expect_equal(mle$estimate$marginal.variance,
               diag(moments$sscp1 + moments$sscp2) / nrow(fixture$x))

  constant <- cbind(fixture$x, constant = 1)
  expect_error(independence_classifier(constant, fixture$y),
               "marginal variance")
  dropped <- independence_classifier(
    constant, fixture$y, zero_variance = "drop"
  )
  expect_identical(unname(dropped$diagnostics$dropped.features), 4L)
  expect_identical(unname(dropped$estimate$direction[4]), 0)
})


test_that("FAIR uses Welch ranking and the average class variance", {
  fixture <- ch5_classical_fixture()
  moments <- ch5_manual_moments(fixture$x, fixture$y)
  variance1 <- diag(moments$sscp1) / (moments$n1 - 1)
  variance2 <- diag(moments$sscp2) / (moments$n2 - 1)
  delta <- moments$mean1 - moments$mean2
  welch <- delta / sqrt(variance1 / moments$n1 + variance2 / moments$n2)
  ordering <- order(-abs(welch), seq_along(welch))
  fit <- fair_classifier(fixture$x, fixture$y, selection = "m", m = 2)
  expect_equal(fit$estimate$t.statistics, welch, tolerance = 1e-12)
  expect_equal(fit$estimate$marginal.variance,
               (variance1 + variance2) / 2, tolerance = 1e-12)
  expect_identical(unname(fit$diagnostics$feature.order), ordering)
  expect_identical(
    unname(fit$diagnostics$selected.features), ordering[1:2]
  )
  expect_true(all(fit$estimate$direction[-ordering[1:2]] == 0))

  at.boundary <- fair_classifier(
    fixture$x, fixture$y, selection = "threshold",
    t_threshold = abs(welch[ordering[1]])
  )
  expect_false(ordering[1] %in% at.boundary$diagnostics$selected.features)
  expect_match(at.boundary$diagnostics$threshold.comparison, ">", fixed = TRUE)
})


test_that("FAIR paper feature-count criterion matches equation 4.3", {
  fixture <- ch5_classical_fixture()
  fit <- fair_classifier(fixture$x, fixture$y, selection = "paper")
  moments <- ch5_manual_moments(fixture$x, fixture$y)
  variance1 <- diag(moments$sscp1) / (moments$n1 - 1)
  variance2 <- diag(moments$sscp2) / (moments$n2 - 1)
  t.value <- (moments$mean1 - moments$mean2) /
    sqrt(variance1 / moments$n1 + variance2 / moments$n2)
  ordering <- order(-abs(t.value), seq_along(t.value))
  sscp <- moments$sscp1 + moments$sscp2
  correlation <- sscp / sqrt(outer(diag(sscp), diag(sscp)))
  correlation <- correlation[ordering, ordering, drop = FALSE]
  n <- moments$n1 + moments$n2
  criterion <- numeric(length(ordering))
  lambda <- numeric(length(ordering))
  for (m in seq_along(ordering)) {
    lambda[m] <- max(eigen(correlation[seq_len(m), seq_len(m), drop = FALSE],
                           symmetric = TRUE, only.values = TRUE)$values)
    sum.t <- sum(t.value[ordering[seq_len(m)]]^2)
    criterion[m] <- n * (sum.t + m * (moments$n1 - moments$n2) / n)^2 /
      (lambda[m] * (m * moments$n1 * moments$n2 +
                      moments$n1 * moments$n2 * sum.t))
  }
  expect_equal(fit$tuning$criterion, criterion, tolerance = 1e-11)
  expect_equal(fit$tuning$truncated.lambda.max, lambda, tolerance = 1e-11)
  expect_identical(fit$tuning$m, which.max(criterion))
})


test_that("FAIR has explicit equal-prior and zero-variance contracts", {
  fixture <- ch5_classical_fixture()
  unbalanced <- fixture$x[-1, , drop = FALSE]
  labels <- droplevels(fixture$y[-1])
  expect_error(
    fair_classifier(unbalanced, labels, prior = "empirical"),
    "equal-prior"
  )
  constant <- cbind(fixture$x, constant = 1)
  expect_error(fair_classifier(constant, fixture$y), "marginal variance")
  dropped <- fair_classifier(
    constant, fixture$y, zero_variance = "drop", selection = "m", m = 1
  )
  expect_identical(unname(dropped$diagnostics$dropped.features), 4L)
})


test_that("Shao threshold LDA uses divisor n and strict hard thresholds", {
  fixture <- ch5_classical_fixture()
  fit <- shao_threshold_lda(
    fixture$x, fixture$y, M_cov = 0, M_mean = 0, alpha = 0.25
  )
  moments <- ch5_manual_moments(fixture$x, fixture$y)
  covariance <- (moments$sscp1 + moments$sscp2) / nrow(fixture$x)
  expect_equal(fit$estimate$covariance, covariance, tolerance = 1e-12)
  expect_identical(fit$diagnostics$covariance.divisor, nrow(fixture$x))
  expect_equal(fit$estimate$thresholded.covariance, covariance,
               tolerance = 1e-12)
  expect_equal(fit$estimate$thresholded.mean.difference,
               moments$mean1 - moments$mean2)

  matrix0 <- matrix(c(2, 0.5, 0.5, 1), 2)
  thresholded <- HDElliptical:::cpp_ch5_classical_hard_threshold_covariance(
    matrix0, 0.5
  )
  expect_equal(thresholded, diag(c(2, 1)))
  expect_match(fit$diagnostics$covariance.threshold.comparison,
               "absolute value > threshold", fixed = TRUE)
})


test_that("Shao explicit-grid LOO is deterministic and leakage-free", {
  fixture <- ch5_classical_fixture()
  grid <- data.frame(M_cov = c(0, 5), M_mean = c(0, 0.5))
  first <- shao_threshold_lda(
    fixture$x, fixture$y, alpha = 0.25, parameter_grid = grid
  )
  second <- shao_threshold_lda(
    fixture$x, fixture$y, alpha = 0.25, parameter_grid = grid
  )
  expect_identical(first$tuning$correct, second$tuning$correct)
  expect_identical(first$tuning$eligible, second$tuning$eligible)
  expect_identical(first$tuning$selected.row, second$tuning$selected.row)
  expect_identical(first$tuning$source,
                   "explicit grid with leave-one-out classification")
  expect_match(first$diagnostics$cv, "all nuisance estimates refit")
  expect_error(
    shao_threshold_lda(
      fixture$x, fixture$y, M_cov = 0, M_mean = 0, alpha = 0.25,
      parameter_grid = grid
    ),
    "either"
  )
  expect_error(
    shao_threshold_lda(fixture$x, fixture$y, alpha = 0.25),
    "Supply both"
  )
})


test_that("Shao contracts reject invented priors and positive-definite repair", {
  fixture <- ch5_classical_fixture()
  expect_error(
    shao_threshold_lda(
      fixture$x[-1, ], droplevels(fixture$y[-1]),
      M_cov = 0, M_mean = 0, alpha = 0.25, prior = "empirical"
    ),
    "equal-prior"
  )
  singular <- cbind(fixture$x[, 1], fixture$x[, 1])
  expect_error(
    shao_threshold_lda(
      singular, fixture$y, M_cov = 0, M_mean = 0, alpha = 0.25
    ),
    "No ridge"
  )
  expect_warning(
    invalid <- shao_threshold_lda(
      singular, fixture$y, M_cov = 0, M_mean = 0, alpha = 0.25,
      strict = FALSE
    ),
    "No ridge"
  )
  expect_false(invalid$valid)
  expect_error(
    shao_threshold_lda(
      fixture$x, fixture$y, M_cov = 0, M_mean = 0, alpha = 0.5
    ),
    "strictly between"
  )
})


test_that("native Chapter 5 classical kernels match direct R calculations", {
  fixture <- ch5_classical_fixture()
  coefficients <- c(1, -2, 0.5)
  intercept <- -0.25
  expect_equal(
    as.numeric(HDElliptical:::cpp_ch5_classical_linear_scores(
      fixture$x, coefficients, intercept
    )),
    as.numeric(fixture$x %*% coefficients + intercept),
    tolerance = 1e-13
  )
  quadratic <- matrix(c(1, 0.2, 0, 0.2, -0.5, 0.1, 0, 0.1, 0.3), 3)
  expected <- rowSums((fixture$x %*% quadratic) * fixture$x) +
    as.numeric(fixture$x %*% coefficients) + intercept
  expect_equal(
    as.numeric(HDElliptical:::cpp_ch5_classical_quadratic_scores(
      fixture$x, quadratic, coefficients, intercept
    )),
    unname(expected), tolerance = 1e-12
  )
  moments <- HDElliptical:::cpp_ch5_classical_two_class_moments(
    fixture$x, as.integer(fixture$y)
  )
  reference <- ch5_manual_moments(fixture$x, fixture$y)
  expect_equal(as.numeric(moments$mean1), as.numeric(reference$mean1))
  expect_equal(as.numeric(moments$mean2), as.numeric(reference$mean2))
  expect_equal(unname(moments$sscp1), unname(reference$sscp1),
               tolerance = 1e-13)
  expect_equal(unname(moments$sscp2), unname(reference$sscp2),
               tolerance = 1e-13)
})


test_that("argument validation is strict and feature-safe", {
  fixture <- ch5_classical_fixture()
  duplicated <- fixture$x
  colnames(duplicated) <- c("x", "x", "z")
  expect_error(classical_lda_classifier(duplicated, fixture$y), "unique")
  expect_error(classical_lda_classifier(fixture$x, fixture$y, strict = NA),
               "TRUE or FALSE")
  expect_error(independence_classifier(fixture$x, fixture$y,
                                       variance_divisor = "other"))
  expect_error(fair_classifier(fixture$x, fixture$y, selection = "m"),
               "must be supplied")
  expect_error(
    fair_classifier(fixture$x, fixture$y, selection = "paper", m = 1),
    "used only"
  )
  expect_error(
    shao_threshold_lda(
      fixture$x, fixture$y, alpha = 0.25,
      parameter_grid = data.frame(M_cov = c(0, 0), M_mean = c(0, 0))
    ),
    "duplicate"
  )
})
