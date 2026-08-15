# Chapter 5 direct and robust sparse linear classifiers.

.c5lin_scalar <- function(value, name, lower = 0, strict_lower = FALSE) {
  value <- as.numeric(value)
  bad.lower <- if (strict_lower) value <= lower else value < lower
  if (length(value) != 1L || is.na(value) || !is.finite(value) ||
      bad.lower) {
    relation <- if (strict_lower) "greater than" else "at least"
    stop(sprintf("'%s' must be one finite number %s %s.",
                 name, relation, format(lower)), call. = FALSE)
  }
  value
}


.c5lin_positive_integer <- function(value, name) {
  if (!is.numeric(value) || length(value) != 1L || is.na(value) ||
      !is.finite(value) || value < 1 || value != floor(value) ||
      value > .Machine$integer.max) {
    stop(sprintf("'%s' must be a positive integer.", name), call. = FALSE)
  }
  as.integer(value)
}


.c5lin_controls <- function(lambda, ridge, solver_tol, solver_max_iter,
                            strict, tie) {
  if (!is.logical(strict) || length(strict) != 1L || is.na(strict)) {
    stop("'strict' must be TRUE or FALSE.", call. = FALSE)
  }
  list(
    lambda = .c5lin_scalar(lambda, "lambda", 0, TRUE),
    ridge = .c5lin_scalar(ridge, "ridge", 0, FALSE),
    solver.tol = .c5lin_scalar(solver_tol, "solver_tol", 0, TRUE),
    solver.max.iter = .c5lin_positive_integer(
      solver_max_iter, "solver_max_iter"
    ),
    strict = strict,
    tie = match.arg(tie, c("class1", "class2"))
  )
}


.c5lin_named <- function(value, feature.names) {
  value <- as.numeric(value)
  if (!is.null(feature.names)) names(value) <- feature.names
  value
}


.c5lin_matrix <- function(value, p, name, feature.names,
                          positive.definite = TRUE) {
  if (is.data.frame(value)) value <- data.matrix(value)
  if (!is.matrix(value) || !is.numeric(value) ||
      !identical(dim(value), c(p, p))) {
    stop(sprintf("'%s' must be a numeric %d by %d matrix.",
                 name, p, p), call. = FALSE)
  }
  storage.mode(value) <- "double"
  if (anyNA(value) || any(!is.finite(value))) {
    stop(sprintf("'%s' must contain only finite values.", name),
         call. = FALSE)
  }
  symmetry.bound <- 64 * .Machine$double.eps *
    max(1, max(abs(value)))
  symmetry.error <- max(abs(value - t(value)))
  if (!is.finite(symmetry.error) || symmetry.error > symmetry.bound) {
    stop(sprintf("'%s' must be symmetric; it is not symmetrized.",
                 name), call. = FALSE)
  }
  if (!is.null(feature.names)) {
    if (!is.null(rownames(value)) &&
        !identical(rownames(value), feature.names)) {
      stop(sprintf("Row names of '%s' must match the training features.",
                   name), call. = FALSE)
    }
    if (!is.null(colnames(value)) &&
        !identical(colnames(value), feature.names)) {
      stop(sprintf("Column names of '%s' must match the training features.",
                   name), call. = FALSE)
    }
    dimnames(value) <- list(feature.names, feature.names)
  }
  if (positive.definite) {
    factor <- tryCatch(chol(value), error = identity)
    if (inherits(factor, "condition")) {
      stop(sprintf(
        "'%s' must be strictly positive definite; no ridge, eigenvalue floor, or pseudoinverse is applied.",
        name
      ), call. = FALSE)
    }
  }
  attr(value, "symmetry.error") <- symmetry.error
  attr(value, "symmetry.bound") <- symmetry.bound
  value
}


.c5lin_class_parts <- function(training) {
  x1 <- training$x[training$class1, , drop = FALSE]
  x2 <- training$x[training$class2, , drop = FALSE]
  list(x1 = x1, x2 = x2)
}


.c5lin_sample_locations <- function(training) {
  parts <- .c5lin_class_parts(training)
  list(
    class1 = .c5lin_named(colMeans(parts$x1), training$feature.names),
    class2 = .c5lin_named(colMeans(parts$x2), training$feature.names),
    source = "classwise sample means"
  )
}


.c5lin_validate_locations <- function(locations, training) {
  if (!is.list(locations) || length(locations) != 2L) {
    stop("'locations' must be NULL or a list of two location vectors.",
         call. = FALSE)
  }
  first <- .as_location(locations[[1L]], training$p, "locations[[1]]")
  second <- .as_location(locations[[2L]], training$p, "locations[[2]]")
  for (location in list(first, second)) {
    if (!is.null(names(location)) && !is.null(training$feature.names) &&
        !identical(names(location), training$feature.names)) {
      stop("Named supplied locations must match the training features.",
           call. = FALSE)
    }
  }
  list(
    class1 = .c5lin_named(first, training$feature.names),
    class2 = .c5lin_named(second, training$feature.names),
    source = "supplied classwise locations"
  )
}


.c5lin_spatial_locations <- function(training, median_tol,
                                     median_max_iter, zero_tol) {
  controls <- .validate_iteration_controls(
    median_tol, median_max_iter, zero_tol
  )
  parts <- .c5lin_class_parts(training)
  first <- .ch3pp_spatial_fit(
    parts$x1, controls$tol, controls$max_iter, controls$zero_tol
  )
  second <- .ch3pp_spatial_fit(
    parts$x2, controls$tol, controls$max_iter, controls$zero_tol
  )
  list(
    valid = isTRUE(first$diagnostics$converged) &&
      isTRUE(second$diagnostics$converged),
    class1 = .c5lin_named(first$center, training$feature.names),
    class2 = .c5lin_named(second$center, training$feature.names),
    signs1 = first$signs,
    signs2 = second$signs,
    fits = list(class1 = first$diagnostics, class2 = second$diagnostics),
    zero.radii = c(class1 = first$n.zero, class2 = second$n.zero),
    radial.scale = c(class1 = first$radial.scale,
                     class2 = second$radial.scale),
    controls = controls,
    source = "classwise sample spatial medians"
  )
}


.c5lin_dantzig <- function(operator, target, controls) {
  answer <- tryCatch(
    cpp_c5lin_dantzig(
      operator, target, controls$lambda, controls$solver.tol,
      controls$solver.max.iter
    ),
    error = identity
  )
  if (inherits(answer, "condition")) {
    return(list(
      valid = FALSE,
      failure = conditionMessage(answer),
      certificate = list(error = conditionMessage(answer)),
      solution = NULL
    ))
  }
  gap.scale <- max(
    1, abs(answer$primal_objective), abs(answer$dual_objective)
  )
  valid <- isTRUE(answer$converged) &&
    !isTRUE(answer$infeasible_zero_operator) &&
    is.finite(answer$primal_violation) &&
    answer$primal_violation <= controls$solver.tol &&
    is.finite(answer$dual_violation) &&
    answer$dual_violation <= controls$solver.tol &&
    is.finite(answer$stationarity_residual) &&
    answer$stationarity_residual <= controls$solver.tol &&
    is.finite(answer$relative_gap) &&
    answer$relative_gap <= controls$solver.tol &&
    is.finite(answer$duality_gap) &&
    answer$duality_gap >= -controls$solver.tol * gap.scale
  certificate <- list(
    certified = valid,
    converged = isTRUE(answer$converged),
    iterations = as.integer(answer$iterations),
    relative.update = as.numeric(answer$relative_update),
    primal.violation = as.numeric(answer$primal_violation),
    stationarity.residual = as.numeric(answer$stationarity_residual),
    dual.violation = as.numeric(answer$dual_violation),
    primal.objective = as.numeric(answer$primal_objective),
    dual.objective = as.numeric(answer$dual_objective),
    duality.gap = as.numeric(answer$duality_gap),
    relative.gap = as.numeric(answer$relative_gap),
    operator.norm = as.numeric(answer$operator_norm),
    primal.dual.step = as.numeric(answer$primal_dual_step),
    zero.operator = isTRUE(answer$operator_zero),
    infeasible.zero.operator =
      isTRUE(answer$infeasible_zero_operator),
    dual = as.numeric(answer$dual),
    last.iterate = if (valid) NULL else as.numeric(answer$solution)
  )
  list(
    valid = valid,
    failure = if (valid) NULL else paste(
      "the Dantzig optimizer did not pass primal feasibility, dual",
      "feasibility, l1 stationarity, and relative-gap checks"
    ),
    certificate = certificate,
    solution = if (valid) as.numeric(answer$solution) else NULL
  )
}


.c5lin_failure <- function(method, training, call, data.name, tie,
                           score.scale, stage, message, diagnostics,
                           tuning, strict) {
  full.message <- paste0(method, " failed at ", stage, ": ", message,
                         ". No fitted classifier was returned and no ",
                         "numerical repair was applied.")
  if (isTRUE(strict)) stop(full.message, call. = FALSE)
  warning(full.message, call. = FALSE)
  .clf_new_fit(
    method = method,
    training = training,
    score_model = list(
      type = "linear",
      coefficients = rep.int(NA_real_, training$p),
      intercept = NA_real_
    ),
    estimate = list(),
    tuning = tuning,
    diagnostics = c(
      list(
        failure.stage = stage,
        failure = full.message,
        no.repair = paste(
          "No ridge beyond the explicit ridge argument, eigenvalue floor,",
          "pseudoinverse, constraint relaxation, or post-hoc repair"
        )
      ),
      diagnostics
    ),
    call = call,
    primary_orientation = "positive score selects class1",
    score_scale = score.scale,
    tie = tie,
    data_name = data.name,
    valid = FALSE
  )
}


.c5lin_success <- function(method, training, score.model, estimate,
                           tuning, diagnostics, call, data.name,
                           score.scale, tie) {
  .clf_new_fit(
    method = method,
    training = training,
    score_model = score.model,
    estimate = estimate,
    tuning = tuning,
    diagnostics = c(
      diagnostics,
      list(
        no.implicit.repair = paste(
          "Only explicitly requested regularization is used; no",
          "pseudoinverse, eigenvalue clipping, or constraint relaxation"
        )
      )
    ),
    call = call,
    primary_orientation = "positive score selects class1",
    score_scale = score.scale,
    tie = tie,
    data_name = data.name,
    valid = TRUE
  )
}


#' Linear programming discriminant classifier
#'
#' Fits the Cai--Liu linear programming discriminant (LPD) direction
#' \deqn{\widehat\gamma\in\arg\min_\gamma\|\gamma\|_1
#' :\|A\gamma-(\bar x_1-\bar x_2)\|_\infty\leq\lambda,}
#' where \eqn{A} is the unbiased pooled within-class covariance plus the
#' explicitly supplied \eqn{ridge I}. The package score is
#' \deqn{(z-(\bar x_1+\bar x_2)/2)^\top\widehat\gamma+
#' \log(\pi_1/\pi_2),}
#' so a positive score selects class 1. The Dantzig solution is returned only
#' after primal feasibility, dual feasibility, l1 stationarity, and a
#' primal--dual gap are all certified.
#'
#' @param x Numeric training matrix with observations in rows.
#' @param y Binary response. Its first factor level is class 1.
#' @param lambda Positive Dantzig constraint radius; it is never selected
#'   implicitly.
#' @param prior Two positive class probabilities summing to one. NULL uses
#'   the shared classifier default.
#' @param ridge Explicit non-negative ridge added to the pooled covariance.
#'   Zero implements the unmodified paper operator.
#' @param solver_tol Positive tolerance required by every solver certificate.
#' @param solver_max_iter Positive maximum number of primal--dual iterations.
#' @param strict If TRUE, a failed certificate is an error. If FALSE, an
#'   explicitly invalid, non-predictable fit is returned with a warning.
#' @param tie Which class receives a score exactly equal to zero.
#'
#' @return An object of class hd_classifier_fit.
#'
#' @references
#' Cai, T. and Liu, W. (2011). A direct estimation approach to sparse linear
#' discriminant analysis. Journal of the American Statistical Association,
#' 106, 1566--1577. \doi{10.1198/jasa.2011.tm11199}.
#'
#' @examples
#' x <- rbind(
#'   c(2, 1), c(1, 2), c(2, 2), c(3, 1),
#'   c(-2, -1), c(-1, -2), c(-2, -2), c(-3, -1)
#' )
#' y <- factor(rep(c("first", "second"), each = 4))
#' fit <- lpd_classifier(x, y, lambda = 0.2, ridge = 0.1)
#' predict(fit, x, type = "score")
#'
#' @export
lpd_classifier <- function(
    x, y, lambda, prior = NULL, ridge = 0,
    solver_tol = 1e-7, solver_max_iter = 100000L,
    strict = TRUE, tie = c("class1", "class2")) {
  call <- match.call()
  data.name <- deparse1(substitute(x))
  training <- .clf_prepare_xy(
    x, y, prior = prior, equal_prior = FALSE, min_class = 2L
  )
  controls <- .c5lin_controls(
    lambda, ridge, solver_tol, solver_max_iter, strict, tie
  )
  parts <- .c5lin_class_parts(training)
  locations <- .c5lin_sample_locations(training)
  residual1 <- sweep(parts$x1, 2L, locations$class1, "-")
  residual2 <- sweep(parts$x2, 2L, locations$class2, "-")
  pooled <- (crossprod(residual1) + crossprod(residual2)) /
    (training$n - 2)
  if (any(!is.finite(pooled))) {
    return(.c5lin_failure(
      "Cai-Liu LPD classifier", training, call, data.name, controls$tie,
      "LPD covariance-inverse linear score", "pooled covariance",
      "the unbiased pooled covariance is non-finite in original units",
      list(), controls, controls$strict
    ))
  }
  operator <- pooled
  diag(operator) <- diag(operator) + controls$ridge
  if (any(!is.finite(operator))) {
    return(.c5lin_failure(
      "Cai-Liu LPD classifier", training, call, data.name, controls$tie,
      "LPD covariance-inverse linear score", "explicit ridge operator",
      "the covariance-plus-ridge operator is non-finite",
      list(), controls, controls$strict
    ))
  }
  if (!is.null(training$feature.names)) {
    dimnames(pooled) <- dimnames(operator) <-
      list(training$feature.names, training$feature.names)
  }
  delta <- locations$class1 - locations$class2
  solver <- .c5lin_dantzig(operator, delta, controls)
  if (!solver$valid) {
    return(.c5lin_failure(
      "Cai-Liu LPD classifier", training, call, data.name, controls$tie,
      "LPD covariance-inverse linear score", "Dantzig optimization",
      solver$failure, list(solver = solver$certificate),
      list(lambda = controls$lambda, ridge = controls$ridge,
           solver.tolerance = controls$solver.tol,
           solver.maximum.iterations = controls$solver.max.iter),
      controls$strict
    ))
  }
  direction <- .c5lin_named(solver$solution, training$feature.names)
  midpoint <- (locations$class1 + locations$class2) / 2
  intercept <- -sum(midpoint * direction) + training$prior$log.ratio
  .c5lin_success(
    "Cai-Liu LPD classifier", training,
    list(type = "linear", coefficients = direction,
         intercept = intercept),
    estimate = list(
      direction = direction,
      locations = locations[c("class1", "class2", "source")],
      midpoint = midpoint,
      mean.difference = delta,
      pooled.covariance = pooled,
      constraint.operator = operator
    ),
    tuning = list(
      lambda = controls$lambda,
      ridge = controls$ridge,
      ridge.changes.operator = controls$ridge > 0,
      solver.tolerance = controls$solver.tol,
      solver.maximum.iterations = controls$solver.max.iter
    ),
    diagnostics = list(
      solver = solver$certificate,
      pooled.covariance.divisor = training$n - 2,
      prior.threshold = training$prior$log.ratio,
      primary.problem = paste(
        "min ||gamma||_1 subject to",
        "||operator gamma - (mean1-mean2)||_infinity <= lambda"
      )
    ),
    call = call, data.name = data.name,
    score.scale = paste(
      "LPD covariance-inverse direction with additive",
      "log(pi1/pi2) threshold"
    ),
    tie = controls$tie
  )
}


#' Direct sparse discriminant analysis
#'
#' Fits the Mai--Zou--Yuan DSDA lasso with the exact response coding
#' \eqn{-n/n_1} for class 1 and \eqn{+n/n_2} for class 2 and objective
#' \deqn{n^{-1}\sum_i(y_i-a-x_i^\top\beta)^2+
#' \lambda\|\beta\|_1.}
#' The regression coefficient points toward class 2. The package reverses it
#' only when constructing its documented score, so positive always means
#' class 1. Following the primary equal-prior rule, classification uses the
#' midpoint intercept rather than the penalized regression intercept. The
#' coordinate-descent solution must pass both KKT and Fenchel dual-gap checks.
#'
#' @inheritParams lpd_classifier
#' @param prior Equal class probabilities. Non-equal probabilities are
#'   rejected in this first implementation because the primary DSDA
#'   classification intercept has a scale-dependent unequal-prior correction.
#'
#' @return An object of class hd_classifier_fit.
#'
#' @references
#' Mai, Q., Zou, H. and Yuan, M. (2012). A direct approach to sparse
#' discriminant analysis in ultra-high dimensions. Biometrika, 99, 29--42.
#' \doi{10.1093/biomet/asr066}.
#'
#' @examples
#' x <- rbind(
#'   c(2, 1), c(1, 2), c(2, 2), c(3, 1),
#'   c(-2, -1), c(-1, -2), c(-2, -2), c(-3, -1)
#' )
#' y <- factor(rep(c("first", "second"), each = 4))
#' dsda_classifier(x, y, lambda = 0.2)
#'
#' @export
dsda_classifier <- function(
    x, y, lambda, prior = c(0.5, 0.5),
    solver_tol = 1e-7, solver_max_iter = 100000L,
    strict = TRUE, tie = c("class1", "class2")) {
  call <- match.call()
  data.name <- deparse1(substitute(x))
  training <- .clf_prepare_xy(
    x, y, prior = prior, equal_prior = TRUE, min_class = 2L
  )
  controls <- .c5lin_controls(
    lambda, 0, solver_tol, solver_max_iter, strict, tie
  )
  response <- numeric(training$n)
  response[training$class1] <- -training$n / training$n1
  response[training$class2] <- training$n / training$n2
  solver <- tryCatch(
    cpp_c5lin_dsda(
      training$x, response, controls$lambda, controls$solver.tol,
      controls$solver.max.iter
    ),
    error = identity
  )
  if (inherits(solver, "condition")) {
    return(.c5lin_failure(
      "Mai-Zou-Yuan DSDA classifier", training, call, data.name,
      controls$tie, "negative DSDA regression direction",
      "coordinate-descent optimization", conditionMessage(solver),
      list(solver = list(error = conditionMessage(solver))),
      list(lambda = controls$lambda,
           solver.tolerance = controls$solver.tol,
           solver.maximum.iterations = controls$solver.max.iter),
      controls$strict
    ))
  }
  gap.scale <- max(
    1, abs(solver$primal_objective), abs(solver$dual_objective)
  )
  certified <- isTRUE(solver$converged) &&
    is.finite(solver$kkt_residual) &&
    solver$kkt_residual <= controls$solver.tol &&
    is.finite(solver$dual_feasibility) &&
    solver$dual_feasibility <= controls$solver.tol &&
    is.finite(solver$relative_gap) &&
    solver$relative_gap <= controls$solver.tol &&
    is.finite(solver$duality_gap) &&
    solver$duality_gap >= -controls$solver.tol * gap.scale
  certificate <- list(
    certified = certified,
    converged = isTRUE(solver$converged),
    iterations = as.integer(solver$iterations),
    relative.update = as.numeric(solver$relative_update),
    kkt.residual = as.numeric(solver$kkt_residual),
    primal.objective = as.numeric(solver$primal_objective),
    dual.objective = as.numeric(solver$dual_objective),
    duality.gap = as.numeric(solver$duality_gap),
    relative.gap = as.numeric(solver$relative_gap),
    dual.feasibility = as.numeric(solver$dual_feasibility),
    dual.scale = as.numeric(solver$dual_scale),
    dual = as.numeric(solver$dual),
    last.iterate = if (certified) NULL else
      as.numeric(solver$coefficient)
  )
  if (!certified) {
    return(.c5lin_failure(
      "Mai-Zou-Yuan DSDA classifier", training, call, data.name,
      controls$tie, "negative DSDA regression direction",
      "coordinate-descent certificate",
      paste("the optimizer did not pass the KKT and Fenchel",
            "primal-dual gap checks"),
      list(solver = certificate),
      list(lambda = controls$lambda,
           solver.tolerance = controls$solver.tol,
           solver.maximum.iterations = controls$solver.max.iter),
      controls$strict
    ))
  }
  regression.direction <- .c5lin_named(
    solver$coefficient, training$feature.names
  )
  direction <- -regression.direction
  locations <- .c5lin_sample_locations(training)
  midpoint <- (locations$class1 + locations$class2) / 2
  classification.intercept <- -sum(midpoint * direction)
  .c5lin_success(
    "Mai-Zou-Yuan DSDA classifier", training,
    list(type = "linear", coefficients = direction,
         intercept = classification.intercept),
    estimate = list(
      direction = direction,
      regression.direction = regression.direction,
      regression.intercept = as.numeric(solver$regression_intercept),
      classification.intercept = classification.intercept,
      locations = locations[c("class1", "class2", "source")],
      midpoint = midpoint,
      response = response,
      fitted.regression = as.numeric(solver$fitted),
      regression.residual = as.numeric(solver$residual)
    ),
    tuning = list(
      lambda = controls$lambda,
      solver.tolerance = controls$solver.tol,
      solver.maximum.iterations = controls$solver.max.iter
    ),
    diagnostics = list(
      solver = certificate,
      response.coding = c(class1 = -training$n / training$n1,
                          class2 = training$n / training$n2),
      objective = paste(
        "n^{-1} sum (y-a-x'beta)^2 + lambda ||beta||_1"
      ),
      classification.intercept = paste(
        "equal-prior midpoint intercept; the regression intercept",
        "is retained separately and is not used for prediction"
      ),
      package.orientation = "negative of the class2-oriented DSDA direction"
    ),
    call = call, data.name = data.name,
    score.scale = "negative DSDA regression-direction scale",
    tie = controls$tie
  )
}


#' Spatial-sign direct sparse linear discriminant analysis
#'
#' Computes classwise spatial medians and SSCMs, pools the SSCMs with their
#' class sample sizes, and solves
#' \deqn{\min_\gamma\|\gamma\|_1:
#' \|(p\widetilde S+ridge I)\gamma-
#' (\widetilde\mu_1-\widetilde\mu_2)\|_\infty\leq\lambda.}
#' The explicit ridge changes the constraint operator and is therefore stored
#' as part of the method rather than treated as a numerical repair. This
#' implementation follows the equal-prior SSLDA rule.
#'
#' @inheritParams lpd_classifier
#' @param prior Equal class probabilities; non-equal probabilities are
#'   rejected because the SSLDA theory and score in the primary paper assume
#'   equal priors.
#' @param median_tol,median_max_iter,zero_tol Spatial-median and zero-sign
#'   controls.
#'
#' @return An object of class hd_classifier_fit.
#'
#' @references
#' Zhuang, D. and Feng, L. (2025). Spatial sign based direct sparse linear
#' discriminant analysis for high dimensional data.
#' \doi{10.48550/arXiv.2504.11117}.
#'
#' @examples
#' x <- rbind(
#'   c(2, 1), c(1, 2), c(2, 2), c(3, 1),
#'   c(-2, -1), c(-1, -2), c(-2, -2), c(-3, -1)
#' )
#' y <- factor(rep(c("first", "second"), each = 4))
#' sslda(x, y, lambda = 0.25, ridge = 0.1)
#'
#' @export
sslda <- function(
    x, y, lambda, prior = c(0.5, 0.5), ridge = 0,
    median_tol = 1e-8, median_max_iter = 500L, zero_tol = 0,
    solver_tol = 1e-7, solver_max_iter = 100000L,
    strict = TRUE, tie = c("class1", "class2")) {
  call <- match.call()
  data.name <- deparse1(substitute(x))
  training <- .clf_prepare_xy(
    x, y, prior = prior, equal_prior = TRUE, min_class = 2L
  )
  controls <- .c5lin_controls(
    lambda, ridge, solver_tol, solver_max_iter, strict, tie
  )
  locations <- .c5lin_spatial_locations(
    training, median_tol, median_max_iter, zero_tol
  )
  if (!locations$valid) {
    return(.c5lin_failure(
      "Zhuang-Feng SSLDA classifier", training, call, data.name,
      controls$tie, "scale-free spatial-sign linear score",
      "classwise spatial medians",
      "at least one classwise spatial median did not converge",
      list(spatial.medians = locations$fits),
      list(lambda = controls$lambda, ridge = controls$ridge),
      controls$strict
    ))
  }
  sscm1 <- crossprod(locations$signs1) / training$n1
  sscm2 <- crossprod(locations$signs2) / training$n2
  pooled.sscm <- (
    training$n1 * sscm1 + training$n2 * sscm2
  ) / training$n
  operator <- training$p * pooled.sscm
  diag(operator) <- diag(operator) + controls$ridge
  if (any(!is.finite(operator))) {
    return(.c5lin_failure(
      "Zhuang-Feng SSLDA classifier", training, call, data.name,
      controls$tie, "scale-free spatial-sign linear score",
      "SSCM constraint operator", "the operator is non-finite",
      list(spatial.medians = locations$fits),
      list(lambda = controls$lambda, ridge = controls$ridge),
      controls$strict
    ))
  }
  if (!is.null(training$feature.names)) {
    dims <- list(training$feature.names, training$feature.names)
    dimnames(sscm1) <- dimnames(sscm2) <-
      dimnames(pooled.sscm) <- dimnames(operator) <- dims
  }
  delta <- locations$class1 - locations$class2
  solver <- .c5lin_dantzig(operator, delta, controls)
  if (!solver$valid) {
    return(.c5lin_failure(
      "Zhuang-Feng SSLDA classifier", training, call, data.name,
      controls$tie, "scale-free spatial-sign linear score",
      "Dantzig optimization", solver$failure,
      list(solver = solver$certificate,
           spatial.medians = locations$fits),
      list(lambda = controls$lambda, ridge = controls$ridge,
           solver.tolerance = controls$solver.tol,
           solver.maximum.iterations = controls$solver.max.iter),
      controls$strict
    ))
  }
  direction <- .c5lin_named(solver$solution, training$feature.names)
  midpoint <- (locations$class1 + locations$class2) / 2
  intercept <- -sum(midpoint * direction)
  .c5lin_success(
    "Zhuang-Feng SSLDA classifier", training,
    list(type = "linear", coefficients = direction,
         intercept = intercept),
    estimate = list(
      direction = direction,
      locations = locations[c("class1", "class2", "source")],
      midpoint = midpoint,
      mean.difference = delta,
      classwise.sscm = list(class1 = sscm1, class2 = sscm2),
      pooled.sscm = pooled.sscm,
      constraint.operator = operator
    ),
    tuning = list(
      lambda = controls$lambda,
      ridge = controls$ridge,
      ridge.changes.operator = controls$ridge > 0,
      median.tolerance = locations$controls$tol,
      median.maximum.iterations = locations$controls$max_iter,
      zero.tolerance = locations$controls$zero_tol,
      solver.tolerance = controls$solver.tol,
      solver.maximum.iterations = controls$solver.max.iter
    ),
    diagnostics = list(
      solver = solver$certificate,
      spatial.medians = locations$fits,
      zero.radii = locations$zero.radii,
      radial.scale = locations$radial.scale,
      sscm.pooling = "(n1*S1+n2*S2)/(n1+n2)",
      target.shape.normalization = "trace(Lambda)=p",
      prior.restriction = "equal priors"
    ),
    call = call, data.name = data.name,
    score.scale = "scale-free spatial-sign inverse-shape direction",
    tie = controls$tie
  )
}


#' Sparse matrix plug-in LDA classifier
#'
#' Applies the Gaussian plug-in LDA score using exactly one supplied sparse
#' precision or covariance matrix and classwise sample means (or two explicit
#' locations). A supplied covariance is inverted only after an ordinary
#' Cholesky factorization proves it is strictly positive definite. A supplied
#' precision must itself be strictly positive definite. The function never
#' uses a pseudoinverse, ridge, symmetrization, or eigenvalue clipping.
#'
#' @param x Numeric training matrix with observations in rows.
#' @param y Binary response. Its first factor level is class 1.
#' @param precision Optional supplied sparse precision matrix.
#' @param covariance Optional supplied sparse covariance matrix. Exactly one
#'   of precision and covariance must be supplied.
#' @param locations NULL for classwise sample means, or a list of two finite
#'   classwise location vectors.
#' @param prior Two positive class probabilities summing to one. NULL uses
#'   the shared classifier default.
#' @param tie Which class receives a score exactly equal to zero.
#'
#' @return An object of class hd_classifier_fit.
#'
#' @references
#' Anderson, T. W. (2003). *An Introduction to Multivariate Statistical
#' Analysis*, 3rd ed. Wiley.
#'
#' Feng, L. (2026). *High-Dimensional Data Analysis for Elliptical Symmetric
#' Distributions*, Chapter 5 (book manuscript). The supplied sparse-matrix
#' adapter is a package construction, not a separate primary method.
#'
#' @examples
#' x <- rbind(
#'   c(2, 1), c(1, 2), c(2, 2), c(3, 1),
#'   c(-2, -1), c(-1, -2), c(-2, -2), c(-3, -1)
#' )
#' y <- factor(rep(c("first", "second"), each = 4))
#' sparse_plugin_lda(x, y, precision = diag(2))
#'
#' @export
sparse_plugin_lda <- function(
    x, y, precision = NULL, covariance = NULL, locations = NULL,
    prior = NULL, tie = c("class1", "class2")) {
  call <- match.call()
  data.name <- deparse1(substitute(x))
  training <- .clf_prepare_xy(
    x, y, prior = prior, equal_prior = FALSE, min_class = 2L
  )
  tie <- match.arg(tie)
  if (is.null(precision) == is.null(covariance)) {
    stop("Supply exactly one of 'precision' and 'covariance'.",
         call. = FALSE)
  }
  if (is.null(locations)) {
    locations <- .c5lin_sample_locations(training)
  } else {
    locations <- .c5lin_validate_locations(locations, training)
  }
  if (!is.null(precision)) {
    omega <- .c5lin_matrix(
      precision, training$p, "precision", training$feature.names,
      positive.definite = TRUE
    )
    supplied <- omega
    source <- "supplied precision"
  } else {
    sigma <- .c5lin_matrix(
      covariance, training$p, "covariance", training$feature.names,
      positive.definite = TRUE
    )
    omega <- chol2inv(chol(sigma))
    if (!is.null(training$feature.names)) {
      dimnames(omega) <-
        list(training$feature.names, training$feature.names)
    }
    supplied <- sigma
    source <- "inverse of supplied positive-definite covariance"
  }
  delta <- locations$class1 - locations$class2
  direction <- .c5lin_named(
    as.numeric(omega %*% delta), training$feature.names
  )
  midpoint <- (locations$class1 + locations$class2) / 2
  intercept <- -sum(midpoint * direction) + training$prior$log.ratio
  .c5lin_success(
    "Sparse matrix plug-in LDA classifier", training,
    list(type = "linear", coefficients = direction,
         intercept = intercept),
    estimate = list(
      direction = direction,
      locations = locations[c("class1", "class2", "source")],
      midpoint = midpoint,
      mean.difference = delta,
      precision = omega,
      supplied.matrix = supplied,
      matrix.source = source
    ),
    tuning = list(),
    diagnostics = list(
      matrix.source = source,
      prior.threshold = training$prior$log.ratio,
      matrix.requirement = "symmetric strictly positive definite",
      inversion = if (is.null(precision))
        "ordinary Cholesky inverse" else "none"
    ),
    call = call, data.name = data.name,
    score.scale = "canonical Gaussian log-likelihood-ratio scale",
    tie = tie
  )
}


.c5lin_precision_certificate <- function(fit, training) {
  if (!inherits(fit, "spatial_sign_precision_fit") ||
      !isTRUE(fit$valid) || is.null(fit$estimate) ||
      !fit$method %in% c("sclime", "sglasso")) {
    stop(
      "'precision_fit' must be a valid SCLIME or SGLASSO fit returned by spatial_sign_precision().",
      call. = FALSE
    )
  }
  solver <- fit$diagnostics$solver
  certified <- if (identical(fit$method, "sclime")) {
    isTRUE(solver$all.columns.certified)
  } else {
    isTRUE(solver$kkt.certified)
  }
  if (!certified) {
    stop("'precision_fit' does not retain a valid solver certificate.",
         call. = FALSE)
  }
  estimate <- .c5lin_matrix(
    fit$estimate, training$p, "precision_fit$estimate",
    training$feature.names, positive.definite = FALSE
  )
  list(
    estimate = estimate,
    method = fit$method,
    lambda = fit$lambda,
    solver = solver,
    certificate = if (identical(fit$method, "sclime"))
      "all SCLIME columns certified" else "SGLASSO KKT certified"
  )
}


#' Spatial-sign precision plug-in LDA
#'
#' Combines classwise robust locations with a certified SCLIME or SGLASSO
#' inverse-shape fit from spatial_sign_precision(). Plain matrices,
#' invalid fits, thresholded fits, and fits that have lost their original
#' feasibility/KKT certificate are rejected. The equal-prior score is
#' \deqn{(z-(\widetilde\mu_1+\widetilde\mu_2)/2)^\top
#' \widehat\Omega(\widetilde\mu_1-\widetilde\mu_2).}
#' Because the decision is equal-prior, the inverse-shape scale is irrelevant.
#'
#' @param x Numeric training matrix with observations in rows.
#' @param y Binary response. Its first factor level is class 1.
#' @param precision_fit A valid, certified object returned by
#'   spatial_sign_precision() with method sclime or sglasso.
#' @param locations NULL to compute classwise spatial medians, or a list of
#'   two supplied finite robust location vectors.
#' @param prior Equal class probabilities; non-equal probabilities are
#'   rejected for this scale-free elliptical rule.
#' @param median_tol,median_max_iter,zero_tol Spatial-median controls used
#'   only when locations is NULL.
#' @param strict If TRUE, spatial-median nonconvergence is an error. If FALSE,
#'   it produces an explicitly invalid fit.
#' @param tie Which class receives a score exactly equal to zero.
#'
#' @return An object of class hd_classifier_fit.
#'
#' @references
#' Lu, Z. and Feng, L. (2025). Robust sparse precision matrix estimation and
#' its applications. \doi{10.48550/arXiv.2503.03575}.
#'
#' @examples
#' x <- rbind(
#'   c(2, 1), c(1, 2), c(2, 2), c(3, 1),
#'   c(-2, -1), c(-1, -2), c(-2, -2), c(-3, -1)
#' )
#' y <- factor(rep(c("first", "second"), each = 4))
#' residual <- rbind(
#'   sweep(x[1:4, ], 2, colMeans(x[1:4, ]), "-"),
#'   sweep(x[5:8, ], 2, colMeans(x[5:8, ]), "-")
#' )
#' precision_fit <- spatial_sign_precision(
#'   residual, lambda = 0.5, method = "sglasso"
#' )
#' spatial_sign_precision_lda(x, y, precision_fit)
#'
#' @export
spatial_sign_precision_lda <- function(
    x, y, precision_fit, locations = NULL,
    prior = c(0.5, 0.5), median_tol = 1e-8,
    median_max_iter = 500L, zero_tol = 0, strict = TRUE,
    tie = c("class1", "class2")) {
  call <- match.call()
  data.name <- deparse1(substitute(x))
  training <- .clf_prepare_xy(
    x, y, prior = prior, equal_prior = TRUE, min_class = 2L
  )
  tie <- match.arg(tie)
  if (!is.logical(strict) || length(strict) != 1L || is.na(strict)) {
    stop("'strict' must be TRUE or FALSE.", call. = FALSE)
  }
  precision <- .c5lin_precision_certificate(precision_fit, training)
  spatial <- NULL
  if (is.null(locations)) {
    spatial <- .c5lin_spatial_locations(
      training, median_tol, median_max_iter, zero_tol
    )
    if (!spatial$valid) {
      return(.c5lin_failure(
        "Spatial-sign precision plug-in LDA classifier",
        training, call, data.name, tie,
        "scale-free certified inverse-shape linear score",
        "classwise spatial medians",
        "at least one classwise spatial median did not converge",
        list(spatial.medians = spatial$fits,
             precision.certificate = precision),
        list(), strict
      ))
    }
    locations <- spatial[c("class1", "class2", "source")]
  } else {
    locations <- .c5lin_validate_locations(locations, training)
  }
  delta <- locations$class1 - locations$class2
  direction <- .c5lin_named(
    as.numeric(precision$estimate %*% delta),
    training$feature.names
  )
  midpoint <- (locations$class1 + locations$class2) / 2
  intercept <- -sum(midpoint * direction)
  .c5lin_success(
    "Spatial-sign precision plug-in LDA classifier", training,
    list(type = "linear", coefficients = direction,
         intercept = intercept),
    estimate = list(
      direction = direction,
      locations = locations,
      midpoint = midpoint,
      mean.difference = delta,
      precision = precision$estimate,
      precision.method = precision$method
    ),
    tuning = list(
      precision.lambda = precision$lambda,
      median.tolerance = if (is.null(spatial)) NULL else
        spatial$controls$tol,
      median.maximum.iterations = if (is.null(spatial)) NULL else
        spatial$controls$max_iter,
      zero.tolerance = if (is.null(spatial)) NULL else
        spatial$controls$zero_tol
    ),
    diagnostics = list(
      precision.certificate = precision$certificate,
      precision.solver = precision$solver,
      location.source = locations$source,
      spatial.medians = if (is.null(spatial)) NULL else spatial$fits,
      zero.radii = if (is.null(spatial)) NULL else spatial$zero.radii,
      prior.restriction = "equal priors",
      scale.identification = paste(
        "the certified precision estimates inverse trace-normalized shape;",
        "its scalar is irrelevant to this equal-prior sign rule"
      )
    ),
    call = call, data.name = data.name,
    score.scale = "scale-free certified inverse-shape linear score",
    tie = tie
  )
}
