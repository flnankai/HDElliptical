# Chapter 5: shared classifier infrastructure and classical benchmarks.

.clf_validate_strict <- function(strict) {
  if (!is.logical(strict) || length(strict) != 1L || is.na(strict)) {
    stop("`strict` must be TRUE or FALSE.", call. = FALSE)
  }
  strict
}


.clf_positive_integer <- function(x, name, minimum = 1L) {
  value <- as.numeric(x)
  if (length(value) != 1L || is.na(value) || !is.finite(value) ||
      value < minimum || value != floor(value) ||
      value > .Machine$integer.max) {
    stop(sprintf("`%s` must be one integer at least %d.", name, minimum),
         call. = FALSE)
  }
  as.integer(value)
}


.clf_numeric_matrix <- function(x, name, min_rows = 1L) {
  if (is.data.frame(x)) {
    if (!all(vapply(x, is.numeric, logical(1)))) {
      stop(sprintf("Every column of `%s` must be numeric.", name),
           call. = FALSE)
    }
    x <- as.matrix(x)
  }
  if (is.numeric(x) && is.vector(x)) {
    value.names <- names(x)
    x <- matrix(as.numeric(x), nrow = 1L)
    if (!is.null(value.names)) colnames(x) <- value.names
  }
  if (!is.matrix(x) || !is.numeric(x)) {
    stop(sprintf("`%s` must be a numeric matrix or all-numeric data frame.",
                 name), call. = FALSE)
  }
  storage.mode(x) <- "double"
  if (nrow(x) < min_rows || ncol(x) < 1L) {
    stop(sprintf("`%s` must have at least %d row(s) and one column.",
                 name, min_rows), call. = FALSE)
  }
  if (anyNA(x) || any(!is.finite(x))) {
    stop(sprintf("`%s` must contain only finite values.", name),
         call. = FALSE)
  }
  variable.names <- colnames(x)
  if (!is.null(variable.names) &&
      (anyNA(variable.names) || any(!nzchar(variable.names)) ||
       anyDuplicated(variable.names))) {
    stop(sprintf("Column names of `%s` must be non-empty and unique.", name),
         call. = FALSE)
  }
  x
}


.clf_prior <- function(prior, counts, levels, equal_prior) {
  if (is.null(prior)) prior <- "equal"
  if (is.character(prior)) {
    if (length(prior) != 1L || is.na(prior) ||
        !prior %in% c("equal", "empirical")) {
      stop(
        paste(
          "`prior` must be \"equal\", \"empirical\", or a positive",
          "numeric vector of length two."
        ),
        call. = FALSE
      )
    }
    probabilities <- if (prior == "equal") c(0.5, 0.5) else
      as.numeric(counts) / sum(counts)
    source <- prior
  } else {
    prior.names <- names(prior)
    probabilities <- as.numeric(prior)
    if (length(probabilities) != 2L || anyNA(probabilities) ||
        any(!is.finite(probabilities)) || any(probabilities <= 0)) {
      stop("A numeric `prior` must contain two finite positive values.",
           call. = FALSE)
    }
    if (!is.null(prior.names)) {
      if (length(prior.names) != 2L || anyNA(prior.names) ||
          any(!nzchar(prior.names)) || anyDuplicated(prior.names) ||
          !setequal(prior.names, levels)) {
        stop("Named `prior` entries must match the two class levels exactly.",
             call. = FALSE)
      }
      probabilities <- probabilities[match(levels, prior.names)]
    }
    probabilities <- probabilities / sum(probabilities)
    source <- "specified"
  }
  names(probabilities) <- levels
  if (isTRUE(equal_prior) &&
      max(abs(probabilities - c(0.5, 0.5))) >
      100 * .Machine$double.eps) {
    stop("This method has a strict equal-prior contract.", call. = FALSE)
  }
  list(
    probabilities = probabilities,
    source = source,
    log.ratio = log(probabilities[1L] / probabilities[2L])
  )
}


# Stable shared contract used by all Chapter 5 classifier modules.
.clf_prepare_xy <- function(
    x, y, prior = NULL, equal_prior = FALSE, min_class = 2L) {
  x <- .clf_numeric_matrix(x, "x", min_rows = 2L)
  if (length(y) != nrow(x) || anyNA(y)) {
    stop("`y` must contain one non-missing class label per row of `x`.",
         call. = FALSE)
  }
  if (is.factor(y)) {
    y <- droplevels(y)
  } else {
    if (is.list(y) || is.matrix(y) || is.data.frame(y)) {
      stop("`y` must be an atomic vector or factor.", call. = FALSE)
    }
    labels <- unique(as.character(y))
    y <- factor(as.character(y), levels = labels)
  }
  if (nlevels(y) != 2L) {
    stop("Exactly two observed classes are required.", call. = FALSE)
  }
  min_class <- .clf_positive_integer(min_class, "min_class")
  class.integer <- as.integer(y)
  class1 <- which(class.integer == 1L)
  class2 <- which(class.integer == 2L)
  if (length(class1) < min_class || length(class2) < min_class) {
    stop(sprintf("Each class must contain at least %d observations.", min_class),
         call. = FALSE)
  }
  levels <- levels(y)
  prior.fit <- .clf_prior(
    prior, c(length(class1), length(class2)), levels,
    equal_prior = equal_prior
  )
  list(
    x = x,
    y = y,
    class = class.integer,
    levels = levels,
    class1 = class1,
    class2 = class2,
    n = nrow(x),
    n1 = length(class1),
    n2 = length(class2),
    p = ncol(x),
    feature.names = colnames(x),
    observation.names = rownames(x),
    prior = prior.fit
  )
}


.clf_validate_newdata <- function(newdata, object) {
  if (!inherits(object, "hd_classifier_fit")) {
    stop("`object` must inherit from `hd_classifier_fit`.", call. = FALSE)
  }
  x <- .clf_numeric_matrix(newdata, "newdata", min_rows = 1L)
  if (ncol(x) != object$p) {
    stop(sprintf("`newdata` must have exactly %d columns.", object$p),
         call. = FALSE)
  }
  expected <- object$feature.names
  observed <- colnames(x)
  if (!is.null(expected)) {
    if (!is.null(observed)) {
      if (!setequal(observed, expected)) {
        stop(
          "Column names of `newdata` must match the training features exactly.",
          call. = FALSE
        )
      }
      x <- x[, match(expected, observed), drop = FALSE]
    } else {
      colnames(x) <- expected
    }
  }
  x
}


.clf_validate_spd <- function(x, name, p = NULL) {
  if (!is.matrix(x) || !is.numeric(x)) x <- as.matrix(x)
  storage.mode(x) <- "double"
  if (length(dim(x)) != 2L || nrow(x) != ncol(x) ||
      (!is.null(p) && nrow(x) != p) || anyNA(x) || any(!is.finite(x))) {
    stop(sprintf("`%s` must be a finite square matrix of the required size.",
                 name), call. = FALSE)
  }
  symmetry.tolerance <- 100 * .Machine$double.eps * max(1, max(abs(x)))
  if (max(abs(x - t(x))) > symmetry.tolerance) {
    stop(sprintf("`%s` must be symmetric.", name), call. = FALSE)
  }
  x <- (x + t(x)) / 2
  factor <- tryCatch(chol(x), error = function(e) NULL)
  if (is.null(factor)) {
    stop(sprintf("`%s` must be strictly positive definite.", name),
         call. = FALSE)
  }
  reciprocal.condition <- rcond(x)
  if (!is.finite(reciprocal.condition) || reciprocal.condition <= 0) {
    stop(sprintf("`%s` has no finite positive reciprocal condition number.",
                 name), call. = FALSE)
  }
  list(
    matrix = x,
    chol = factor,
    log.determinant = 2 * sum(log(diag(factor))),
    reciprocal.condition = reciprocal.condition
  )
}


.clf_solve_chol <- function(factor, rhs) {
  backsolve(factor, forwardsolve(t(factor), rhs))
}


.clf_validate_score_model <- function(score_model, p, valid) {
  if (!isTRUE(valid)) return(score_model)
  if (!is.list(score_model) || length(score_model$type) != 1L ||
      !is.character(score_model$type) || is.na(score_model$type)) {
    stop("`score_model` must be a list with one character `type`.",
         call. = FALSE)
  }
  type <- score_model$type
  if (type == "linear") {
    coefficients <- as.numeric(score_model$coefficients)
    intercept <- as.numeric(score_model$intercept)
    if (length(coefficients) != p || anyNA(coefficients) ||
        any(!is.finite(coefficients))) {
      stop("A linear score needs `p` finite coefficients.", call. = FALSE)
    }
    if (length(intercept) != 1L || is.na(intercept) ||
        !is.finite(intercept)) {
      stop("A linear score needs one finite intercept.", call. = FALSE)
    }
    score_model$coefficients <- coefficients
    score_model$intercept <- intercept
  } else if (type == "quadratic") {
    quadratic <- score_model$quadratic
    if (!is.matrix(quadratic) || !is.numeric(quadratic) ||
        !identical(dim(quadratic), c(p, p)) || anyNA(quadratic) ||
        any(!is.finite(quadratic))) {
      stop("A quadratic score needs a finite `p` by `p` matrix.",
           call. = FALSE)
    }
    tolerance <- 100 * .Machine$double.eps *
      max(1, max(abs(quadratic)))
    if (max(abs(quadratic - t(quadratic))) > tolerance) {
      stop("The quadratic score matrix must be symmetric.", call. = FALSE)
    }
    linear <- as.numeric(score_model$linear)
    intercept <- as.numeric(score_model$intercept)
    if (length(linear) != p || anyNA(linear) || any(!is.finite(linear))) {
      stop("A quadratic score needs `p` finite linear coefficients.",
           call. = FALSE)
    }
    if (length(intercept) != 1L || is.na(intercept) ||
        !is.finite(intercept)) {
      stop("A quadratic score needs one finite intercept.", call. = FALSE)
    }
    score_model$quadratic <- (quadratic + t(quadratic)) / 2
    score_model$linear <- linear
    score_model$intercept <- intercept
  } else if (type == "elliptical") {
    required <- c(
      "location1", "location2", "precision1", "precision2",
      "log.determinant1", "log.determinant2", "log.generator1",
      "log.generator2", "log.prior.ratio"
    )
    if (!all(required %in% names(score_model))) {
      stop("The elliptical score model is incomplete.", call. = FALSE)
    }
    if (!is.function(score_model$log.generator1) ||
        !is.function(score_model$log.generator2)) {
      stop("Elliptical log-generators must be functions.", call. = FALSE)
    }
  } else if (type == "custom") {
    if (!is.function(score_model$fun)) {
      stop("A custom score model needs a function `fun`.", call. = FALSE)
    }
    if (length(score_model$description) != 1L ||
        !is.character(score_model$description) ||
        is.na(score_model$description) || !nzchar(score_model$description)) {
      stop("A custom score model needs one non-empty `description`.",
           call. = FALSE)
    }
  } else {
    stop(sprintf("Unsupported classifier score model type `%s`.", type),
         call. = FALSE)
  }
  score_model
}


# Stable shared constructor used by all Chapter 5 classifier modules.
.clf_new_fit <- function(
    method, training = NULL, levels = NULL, feature_names = NULL,
    n_features = NULL, score_model, estimate = list(), tuning = list(),
    diagnostics = list(), call = NULL,
    primary_orientation = "positive score selects class1", score_scale,
    tie = c("class1", "class2"), data_name = NULL, valid = TRUE) {
  if (length(method) != 1L || !is.character(method) || is.na(method) ||
      !nzchar(method)) {
    stop("`method` must be one non-empty character string.", call. = FALSE)
  }
  if (!is.logical(valid) || length(valid) != 1L || is.na(valid)) {
    stop("`valid` must be TRUE or FALSE.", call. = FALSE)
  }
  tie <- match.arg(tie)
  if (!is.null(training)) {
    levels <- training$levels
    feature_names <- training$feature.names
    n_features <- training$p
    prior <- training$prior
  } else {
    if (length(levels) != 2L || anyNA(levels) || anyDuplicated(levels)) {
      stop("`levels` must contain two distinct non-missing labels.",
           call. = FALSE)
    }
    levels <- as.character(levels)
    n_features <- .clf_positive_integer(n_features, "n_features")
    if (!is.null(feature_names)) {
      feature_names <- as.character(feature_names)
      if (length(feature_names) != n_features || anyNA(feature_names) ||
          any(!nzchar(feature_names)) || anyDuplicated(feature_names)) {
        stop("`feature_names` must contain one unique non-empty name per feature.",
             call. = FALSE)
      }
    }
    prior <- score_model$prior
    if (is.null(prior)) {
      prior <- .clf_prior("equal", c(1L, 1L), levels, FALSE)
    }
  }
  if (length(score_scale) != 1L || !is.character(score_scale) ||
      is.na(score_scale) || !nzchar(score_scale)) {
    stop("`score_scale` must be one non-empty character string.",
         call. = FALSE)
  }
  score_model <- .clf_validate_score_model(score_model, n_features, valid)
  score_model$tie <- tie
  structure(
    list(
      valid = valid,
      method = method,
      levels = levels,
      feature.names = feature_names,
      p = n_features,
      prior = prior,
      score.model = score_model,
      estimate = estimate,
      tuning = tuning,
      diagnostics = diagnostics,
      call = call,
      primary.orientation = primary_orientation,
      score.scale = score_scale,
      data.name = data_name
    ),
    class = "hd_classifier_fit"
  )
}


.clf_failure <- function(message, strict, method, training, call,
                         score_scale, stage, tuning = list(),
                         diagnostics = list(), data_name = NULL) {
  strict <- .clf_validate_strict(strict)
  if (strict) stop(message, call. = FALSE)
  warning(message, call. = FALSE)
  diagnostics$failure.stage <- stage
  diagnostics$failure <- message
  diagnostics$no.repair <- paste(
    "No ridge, jitter, pseudoinverse, eigenvalue floor, determinant",
    "absolute value, or post-hoc positive-definite repair was applied."
  )
  .clf_new_fit(
    method = method, training = training,
    score_model = list(type = "invalid"), estimate = NULL,
    tuning = tuning, diagnostics = diagnostics, call = call,
    score_scale = score_scale, data_name = data_name, valid = FALSE
  )
}


.clf_oracle_training <- function(levels, p, feature.names, prior) {
  if (length(levels) != 2L || anyNA(levels) || anyDuplicated(levels)) {
    stop("`levels` must contain two distinct non-missing labels.",
         call. = FALSE)
  }
  levels <- as.character(levels)
  if (!is.null(feature.names)) {
    feature.names <- as.character(feature.names)
    if (length(feature.names) != p || anyNA(feature.names) ||
        any(!nzchar(feature.names)) || anyDuplicated(feature.names)) {
      stop("`feature_names` must contain one unique non-empty name per feature.",
           call. = FALSE)
    }
  }
  list(
    levels = levels, feature.names = feature.names, p = p,
    prior = .clf_prior(prior, c(1L, 1L), levels, FALSE)
  )
}


.clf_score <- function(object, newdata) {
  model <- object$score.model
  if (model$type == "linear") {
    score <- cpp_ch5_classical_linear_scores(
      newdata, model$coefficients, model$intercept
    )
  } else if (model$type == "quadratic") {
    score <- cpp_ch5_classical_quadratic_scores(
      newdata, model$quadratic, model$linear, model$intercept
    )
  } else if (model$type == "elliptical") {
    distances <- cpp_ch5_classical_mahalanobis_pairs(
      newdata, model$location1, model$precision1,
      model$location2, model$precision2
    )
    log1 <- model$log.generator1(distances$distance1)
    log2 <- model$log.generator2(distances$distance2)
    log1 <- as.numeric(log1)
    log2 <- as.numeric(log2)
    n <- nrow(newdata)
    if (length(log1) != n || length(log2) != n || anyNA(log1) ||
        anyNA(log2) || any(is.nan(log1)) || any(is.nan(log2)) ||
        any(log1 == Inf) || any(log2 == Inf)) {
      stop(
        paste(
          "Each log-generator must return one finite or negative-infinite",
          "numeric value per squared Mahalanobis distance."
        ),
        call. = FALSE
      )
    }
    score <- model$log.prior.ratio - model$log.determinant1 / 2 +
      model$log.determinant2 / 2 + log1 - log2
    if (anyNA(score) || any(is.nan(score))) {
      stop("The elliptical log-density contrast is undefined.", call. = FALSE)
    }
  } else if (model$type == "custom") {
    score <- as.numeric(model$fun(newdata))
    if (length(score) != nrow(newdata) || anyNA(score) ||
        any(!is.finite(score))) {
      stop("The custom scorer must return one finite score per row.",
           call. = FALSE)
    }
  } else {
    stop("This fit has no usable score model.", call. = FALSE)
  }
  score <- as.numeric(score)
  if (!is.null(rownames(newdata))) names(score) <- rownames(newdata)
  score
}


#' Construct a classifier from a scoring function
#'
#' This low-level constructor creates a two-class `hd_classifier_fit` from a
#' user-supplied scoring function. The function must accept a finite numeric
#' matrix with observations in rows and return one finite score per row.
#' Positive scores select the first class. Package methods use transparent
#' structured linear or quadratic score models instead of this custom route.
#'
#' @param score Function mapping an observation matrix to numeric scores.
#' @param n_features Positive number of input features.
#' @param levels Two distinct class labels; the first is selected by a positive
#'   score.
#' @param feature_names Optional unique feature names.
#' @param method Descriptive method name.
#' @param score_scale Description of the numerical score scale.
#' @param tie Whether score zero selects the first or second class.
#'
#' @return An object of class `hd_classifier_fit`.
#' @references
#' Feng, L. (2026). *High-Dimensional Data Analysis for Elliptically
#' Symmetric Distributions*, Chapter 5 (book manuscript). This low-level
#' constructor contract is package infrastructure, not a distinct method.
#' @examples
#' fit <- hd_classifier_fit(function(z) z[, 1] - z[, 2], 2)
#' predict(fit, rbind(c(2, 1), c(0, 1)))
#' @export
hd_classifier_fit <- function(
    score, n_features, levels = c("class1", "class2"),
    feature_names = NULL, method = "Custom two-class classifier",
    score_scale = "method_threshold", tie = c("class1", "class2")) {
  call <- match.call()
  if (!is.function(score)) stop("`score` must be a function.", call. = FALSE)
  p <- .clf_positive_integer(n_features, "n_features")
  .clf_new_fit(
    method = method, levels = levels, feature_names = feature_names,
    n_features = p,
    score_model = list(type = "custom", fun = score, description = method),
    diagnostics = list(
      custom.scorer = TRUE,
      no.repair = "The package does not modify custom scores."
    ),
    call = call, score_scale = score_scale, tie = tie
  )
}


#' Predict from a Chapter 5 classifier
#'
#' @param object A valid object inheriting from `hd_classifier_fit`.
#' @param newdata Numeric matrix or all-numeric data frame with observations
#'   in rows.
#' @param type Return class labels or numerical scores.
#' @param ... Reserved for method compatibility; currently unused.
#'
#' @return For `type = "score"`, a numeric vector. For `type = "class"`, a
#'   factor with the training class levels.
#' @export
predict.hd_classifier_fit <- function(
    object, newdata, type = c("class", "score"), ...) {
  dots <- list(...)
  if (length(dots)) stop("Unused arguments were supplied in `...`.", call. = FALSE)
  if (!inherits(object, "hd_classifier_fit")) {
    stop("`object` must inherit from `hd_classifier_fit`.", call. = FALSE)
  }
  if (!isTRUE(object$valid)) {
    failure <- object$diagnostics$failure
    if (is.null(failure)) failure <- "unknown fitting failure"
    stop(paste0("Cannot predict from an invalid classifier fit: ", failure),
         call. = FALSE)
  }
  type <- match.arg(type)
  newdata <- .clf_validate_newdata(newdata, object)
  score <- .clf_score(object, newdata)
  if (type == "score") return(score)
  tie <- object$score.model$tie
  first <- if (identical(tie, "class2")) score > 0 else score >= 0
  answer <- factor(ifelse(first, object$levels[1L], object$levels[2L]),
                   levels = object$levels)
  if (!is.null(rownames(newdata))) names(answer) <- rownames(newdata)
  answer
}


.clf_locations <- function(location1, location2, feature_names = NULL) {
  names1 <- names(location1)
  names2 <- names(location2)
  location1 <- as.numeric(location1)
  location2 <- as.numeric(location2)
  if (!length(location1) || length(location2) != length(location1) ||
      anyNA(location1) || anyNA(location2) ||
      any(!is.finite(location1)) || any(!is.finite(location2))) {
    stop("Class locations must be finite numeric vectors of equal positive length.",
         call. = FALSE)
  }
  if (is.null(feature_names)) {
    if (!is.null(names1)) feature_names <- names1 else if (!is.null(names2))
      feature_names <- names2
  }
  p <- length(location1)
  if (!is.null(feature_names)) {
    feature_names <- as.character(feature_names)
    if (length(feature_names) != p || anyNA(feature_names) ||
        any(!nzchar(feature_names)) || anyDuplicated(feature_names)) {
      stop("`feature_names` must contain one unique non-empty name per feature.",
           call. = FALSE)
    }
    names(location1) <- names(location2) <- feature_names
  }
  list(location1 = location1, location2 = location2,
       feature.names = feature_names, p = p)
}


.clf_matrix_source <- function(covariance, precision, covariance_name,
                               precision_name, p) {
  if (is.null(covariance) == is.null(precision)) {
    stop(sprintf("Supply exactly one of `%s` and `%s`.",
                 covariance_name, precision_name), call. = FALSE)
  }
  if (!is.null(covariance)) {
    covariance.fit <- .clf_validate_spd(covariance, covariance_name, p)
    precision.matrix <- chol2inv(covariance.fit$chol)
    list(
      covariance = covariance.fit$matrix,
      precision = precision.matrix,
      log.covariance.determinant = covariance.fit$log.determinant,
      reciprocal.condition = covariance.fit$reciprocal.condition,
      source = "covariance"
    )
  } else {
    precision.fit <- .clf_validate_spd(precision, precision_name, p)
    list(
      covariance = NULL,
      precision = precision.fit$matrix,
      log.covariance.determinant = -precision.fit$log.determinant,
      reciprocal.condition = precision.fit$reciprocal.condition,
      source = "precision"
    )
  }
}


#' Oracle Gaussian linear discriminant classifier
#'
#' Constructs the exact Gaussian log-likelihood-ratio score
#' \deqn{(z-(\mu_1+\mu_2)/2)'\Omega(\mu_1-\mu_2)
#'       +\log(\pi_1/\pi_2).}
#' Exactly one of `covariance` and `precision` must be supplied. No ridge,
#' generalized inverse, or positive-definite repair is used.
#'
#' @param location1,location2 Finite class-location vectors.
#' @param covariance Optional common positive-definite covariance matrix.
#' @param precision Optional common positive-definite precision matrix.
#' @param prior `"equal"` or a positive numeric vector of length two.
#' @param levels Two class labels.
#' @param feature_names Optional unique feature names.
#'
#' @return A valid `hd_classifier_fit`.
#' @references
#' Anderson, T. W. (2003). *An Introduction to Multivariate Statistical
#' Analysis*, 3rd ed. Wiley.
#' @examples
#' fit <- gaussian_lda_oracle(c(1, 0), c(-1, 0), covariance = diag(2))
#' predict(fit, rbind(c(2, 0), c(-2, 0)))
#' @export
gaussian_lda_oracle <- function(
    location1, location2, covariance = NULL, precision = NULL,
    prior = "equal", levels = c("class1", "class2"),
    feature_names = NULL) {
  call <- match.call()
  locations <- .clf_locations(location1, location2, feature_names)
  matrix.fit <- .clf_matrix_source(
    covariance, precision, "covariance", "precision", locations$p
  )
  training <- .clf_oracle_training(
    levels, locations$p, locations$feature.names, prior
  )
  delta <- locations$location1 - locations$location2
  direction <- as.numeric(matrix.fit$precision %*% delta)
  midpoint <- (locations$location1 + locations$location2) / 2
  intercept <- -sum(midpoint * direction) + training$prior$log.ratio
  if (any(!is.finite(direction)) || !is.finite(intercept)) {
    stop("The oracle LDA score is outside the finite double range.",
         call. = FALSE)
  }
  if (!is.null(locations$feature.names)) names(direction) <-
    locations$feature.names
  .clf_new_fit(
    method = "Oracle Gaussian linear discriminant analysis",
    training = training,
    score_model = list(
      type = "linear", coefficients = direction, intercept = intercept
    ),
    estimate = list(
      location1 = locations$location1,
      location2 = locations$location2,
      covariance = matrix.fit$covariance,
      precision = if (matrix.fit$source == "precision")
        matrix.fit$precision else NULL,
      direction = direction
    ),
    diagnostics = list(
      matrix.source = matrix.fit$source,
      reciprocal.condition = matrix.fit$reciprocal.condition,
      prior.term = training$prior$log.ratio,
      regularization = "none",
      no.repair = paste(
        "Strict Cholesky validation; no ridge, pseudoinverse, jitter,",
        "or eigenvalue floor."
      )
    ),
    call = call, score_scale = "canonical_log_lr", tie = "class1"
  )
}


#' Oracle Gaussian quadratic discriminant classifier
#'
#' Uses the canonical Gaussian log-likelihood ratio
#' \deqn{-\tfrac12\log|\Sigma_1|+\tfrac12\log|\Sigma_2|
#' -\tfrac12d_1^2+\tfrac12d_2^2+\log(\pi_1/\pi_2).}
#' The direct quadratic formula displayed later in the book is twice this
#' score. It has the same boundary but is not the same numerical score.
#'
#' @inheritParams gaussian_lda_oracle
#' @param covariance1,covariance2 Optional class covariance matrices.
#' @param precision1,precision2 Optional class precision matrices. For each
#'   class, supply exactly one covariance or precision matrix.
#'
#' @return A valid `hd_classifier_fit` with a canonical-log-ratio score.
#' @references
#' Anderson, T. W. (2003). *An Introduction to Multivariate Statistical
#' Analysis*, 3rd ed. Wiley.
#' @examples
#' fit <- gaussian_qda_oracle(
#'   c(1, 0), c(-1, 0), covariance1 = diag(2),
#'   covariance2 = diag(c(2, 1))
#' )
#' predict(fit, rbind(c(1, 0), c(-1, 0)), type = "score")
#' @export
gaussian_qda_oracle <- function(
    location1, location2, covariance1 = NULL, covariance2 = NULL,
    precision1 = NULL, precision2 = NULL, prior = "equal",
    levels = c("class1", "class2"), feature_names = NULL) {
  call <- match.call()
  locations <- .clf_locations(location1, location2, feature_names)
  class1.fit <- .clf_matrix_source(
    covariance1, precision1, "covariance1", "precision1", locations$p
  )
  class2.fit <- .clf_matrix_source(
    covariance2, precision2, "covariance2", "precision2", locations$p
  )
  training <- .clf_oracle_training(
    levels, locations$p, locations$feature.names, prior
  )
  omega1 <- class1.fit$precision
  omega2 <- class2.fit$precision
  quadratic <- (omega2 - omega1) / 2
  linear <- as.numeric(omega1 %*% locations$location1 -
                         omega2 %*% locations$location2)
  intercept <- -sum(locations$location1 *
                      as.numeric(omega1 %*% locations$location1)) / 2 +
    sum(locations$location2 *
          as.numeric(omega2 %*% locations$location2)) / 2 -
    class1.fit$log.covariance.determinant / 2 +
    class2.fit$log.covariance.determinant / 2 +
    training$prior$log.ratio
  if (any(!is.finite(quadratic)) || any(!is.finite(linear)) ||
      !is.finite(intercept)) {
    stop("The oracle QDA score is outside the finite double range.",
         call. = FALSE)
  }
  .clf_new_fit(
    method = "Oracle Gaussian quadratic discriminant analysis",
    training = training,
    score_model = list(
      type = "quadratic", quadratic = quadratic,
      linear = linear, intercept = intercept
    ),
    estimate = list(
      location1 = locations$location1,
      location2 = locations$location2,
      covariance1 = class1.fit$covariance,
      covariance2 = class2.fit$covariance,
      precision1 = omega1, precision2 = omega2
    ),
    diagnostics = list(
      matrix.source = c(class1.fit$source, class2.fit$source),
      log.covariance.determinant = c(
        class1.fit$log.covariance.determinant,
        class2.fit$log.covariance.determinant
      ),
      reciprocal.condition = c(
        class1.fit$reciprocal.condition,
        class2.fit$reciprocal.condition
      ),
      book.errata = paste(
        "The book's direct D/beta QDA display is twice the canonical",
        "Gaussian log-likelihood-ratio score used here."
      ),
      regularization = "none",
      no.repair = paste(
        "Strict Cholesky validation; no ridge, pseudoinverse, jitter,",
        "eigenvalue floor, or absolute determinant."
      )
    ),
    call = call, score_scale = "canonical_log_lr", tie = "class1"
  )
}


#' Exact oracle classifier for two elliptical populations
#'
#' Computes the exact generator-aware score
#' \deqn{\log\pi_1-\tfrac12\log|\Lambda_1|+\log g_1(d_1)
#' -\log\pi_2+\tfrac12\log|\Lambda_2|-\log g_2(d_2),}
#' where \eqn{d_k=(z-\mu_k)'\Lambda_k^{-1}(z-\mu_k)}. The
#' log-generator functions must be vectorized. They may return `-Inf` for a
#' zero density, but not `NA`, `NaN`, or `Inf`.
#'
#' A common generator and common shape do not generally reduce this rule to a
#' midpoint linear rule under unequal priors. That shortcut is exact for equal
#' priors, and for the exponential radial generator underlying Gaussian LDA.
#'
#' @param location1,location2 Finite class locations.
#' @param shape1,shape2 Positive-definite class shape matrices, on the scales
#'   used by the corresponding generators.
#' @param log_generator1,log_generator2 Vectorized functions evaluating
#'   \eqn{\log g_k(d)}. By default the second generator equals the first.
#' @inheritParams gaussian_lda_oracle
#'
#' @return A generator-aware `hd_classifier_fit`.
#' @references
#' Fang, K.-T. and Anderson, T. W. (1990). *Statistical Inference in
#' Elliptically Contoured and Related Distributions*. Allerton Press.
#' Wakaki, H. (1994). Discriminant analysis under elliptical populations.
#' *Hiroshima Mathematical Journal*, 24, 257--298.
#' @examples
#' normal_log_generator <- function(d) -d / 2
#' fit <- elliptical_oracle_classifier(
#'   c(1, 0), c(-1, 0), diag(2), diag(2), normal_log_generator
#' )
#' predict(fit, rbind(c(1, 0), c(-1, 0)))
#' @export
elliptical_oracle_classifier <- function(
    location1, location2, shape1, shape2, log_generator1,
    log_generator2 = log_generator1, prior = "equal",
    levels = c("class1", "class2"), feature_names = NULL) {
  call <- match.call()
  locations <- .clf_locations(location1, location2, feature_names)
  if (!is.function(log_generator1) || !is.function(log_generator2)) {
    stop("`log_generator1` and `log_generator2` must be functions.",
         call. = FALSE)
  }
  shape1.fit <- .clf_validate_spd(shape1, "shape1", locations$p)
  shape2.fit <- .clf_validate_spd(shape2, "shape2", locations$p)
  training <- .clf_oracle_training(
    levels, locations$p, locations$feature.names, prior
  )
  score.model <- list(
    type = "elliptical",
    location1 = locations$location1,
    location2 = locations$location2,
    precision1 = chol2inv(shape1.fit$chol),
    precision2 = chol2inv(shape2.fit$chol),
    log.determinant1 = shape1.fit$log.determinant,
    log.determinant2 = shape2.fit$log.determinant,
    log.generator1 = log_generator1,
    log.generator2 = log_generator2,
    log.prior.ratio = training$prior$log.ratio
  )
  .clf_new_fit(
    method = "Exact oracle elliptical discriminant classifier",
    training = training, score_model = score.model,
    estimate = list(
      location1 = locations$location1,
      location2 = locations$location2,
      shape1 = shape1.fit$matrix, shape2 = shape2.fit$matrix
    ),
    diagnostics = list(
      generator.argument = "squared Mahalanobis distance",
      reciprocal.condition = c(
        shape1.fit$reciprocal.condition,
        shape2.fit$reciprocal.condition
      ),
      unequal.prior.contract = paste(
        "The exact generator score is retained; no midpoint-linear",
        "shortcut is made."
      ),
      book.errata = paste(
        "Common generator and common shape imply a midpoint-linear rule",
        "only for equal priors or an exponential radial generator."
      ),
      regularization = "none",
      no.repair = paste(
        "Strict Cholesky validation; no ridge, pseudoinverse, jitter,",
        "eigenvalue floor, or absolute determinant."
      )
    ),
    call = call, score_scale = "canonical_log_lr", tie = "class1"
  )
}

.clf_two_class_moments <- function(training) {
  answer <- tryCatch(
    cpp_ch5_classical_two_class_moments(training$x, training$class),
    error = identity
  )
  if (inherits(answer, "condition")) return(answer)
  variable.names <- training$feature.names
  for (name in c("mean1", "mean2")) {
    value <- as.numeric(answer[[name]])
    if (!is.null(variable.names)) names(value) <- variable.names
    answer[[name]] <- value
  }
  for (name in c("sscp1", "sscp2")) {
    value <- as.matrix(answer[[name]])
    if (!is.null(variable.names)) {
      dimnames(value) <- list(variable.names, variable.names)
    }
    answer[[name]] <- value
  }
  answer
}


.clf_classical_failure <- function(core, method, training, strict, call,
                                    score_scale, data.name) {
  if (inherits(core, "condition")) {
    message <- paste0(
      method, " moment computation failed: ", conditionMessage(core),
      " No numerical repair was applied."
    )
    return(.clf_failure(
      message, strict, method, training, call, score_scale,
      "two-class moment computation", data_name = data.name
    ))
  }
  NULL
}


.clf_name_estimates <- function(estimate, training) {
  variables <- training$feature.names
  if (is.null(variables)) return(estimate)
  vector.names <- c("location1", "location2", "direction",
                    "mean.difference", "marginal.variance",
                    "t.statistics")
  matrix.names <- c("covariance", "covariance1", "covariance2",
                    "thresholded.covariance", "correlation")
  for (name in intersect(names(estimate), vector.names)) {
    value <- estimate[[name]]
    if (length(value) == training$p) names(value) <- variables
    estimate[[name]] <- value
  }
  for (name in intersect(names(estimate), matrix.names)) {
    value <- estimate[[name]]
    if (is.matrix(value) && identical(dim(value), c(training$p, training$p))) {
      dimnames(value) <- list(variables, variables)
      estimate[[name]] <- value
    }
  }
  estimate
}


#' Classical plug-in linear discriminant classifier
#'
#' Estimates the common covariance with the pooled unbiased estimator and
#' inserts it into the Gaussian LDA log-likelihood ratio. The covariance must
#' be strictly positive definite. No generalized inverse or ridge is used.
#'
#' @param x Numeric matrix or all-numeric data frame, observations in rows.
#' @param y Two-class label vector or factor.
#' @param prior `"equal"`, `"empirical"`, or a positive numeric vector of
#'   length two. Named entries must match the class levels.
#' @param strict If `TRUE`, a method failure is an error. If `FALSE`, it is a
#'   warning followed by an invalid `hd_classifier_fit`.
#'
#' @return An `hd_classifier_fit`.
#' @references
#' Anderson, T. W. (2003). *An Introduction to Multivariate Statistical
#' Analysis*, 3rd ed. Wiley.
#' @examples
#' x <- rbind(c(2, 0), c(1, 1), c(1, -1),
#'            c(-2, 0), c(-1, 1), c(-1, -1))
#' classical_lda_classifier(x, rep(c("A", "B"), each = 3))
#' @export
classical_lda_classifier <- function(
    x, y, prior = "equal", strict = TRUE) {
  call <- match.call()
  data.name <- deparse1(substitute(x))
  strict <- .clf_validate_strict(strict)
  training <- .clf_prepare_xy(x, y, prior = prior, min_class = 2L)
  method <- "Classical plug-in Gaussian linear discriminant analysis"
  core <- .clf_two_class_moments(training)
  failed <- .clf_classical_failure(
    core, method, training, strict, call, "canonical_log_lr", data.name
  )
  if (!is.null(failed)) return(failed)
  covariance <- (core$sscp1 + core$sscp2) / (training$n - 2L)
  covariance.fit <- tryCatch(
    .clf_validate_spd(covariance, "pooled covariance", training$p),
    error = identity
  )
  if (inherits(covariance.fit, "condition")) {
    return(.clf_failure(
      paste0(
        "Classical LDA requires a strictly positive-definite pooled ",
        "covariance: ", conditionMessage(covariance.fit),
        " No ridge or generalized inverse was used."
      ),
      strict, method, training, call, "canonical_log_lr",
      "pooled covariance factorization", data_name = data.name
    ))
  }
  delta <- core$mean1 - core$mean2
  direction <- as.numeric(.clf_solve_chol(covariance.fit$chol, delta))
  midpoint <- (core$mean1 + core$mean2) / 2
  intercept <- -sum(midpoint * direction) + training$prior$log.ratio
  if (any(!is.finite(direction)) || !is.finite(intercept)) {
    return(.clf_failure(
      paste(
        "The classical LDA direction or intercept is outside the finite",
        "double range; no clipping was applied."
      ),
      strict, method, training, call, "canonical_log_lr",
      "score construction", data_name = data.name
    ))
  }
  estimate <- .clf_name_estimates(list(
    location1 = core$mean1,
    location2 = core$mean2,
    covariance = covariance,
    direction = direction,
    mean.difference = delta
  ), training)
  .clf_new_fit(
    method = method, training = training,
    score_model = list(
      type = "linear", coefficients = direction, intercept = intercept
    ),
    estimate = estimate,
    diagnostics = list(
      covariance.divisor = training$n - 2L,
      residual.degrees.of.freedom = training$n - 2L,
      reciprocal.condition = covariance.fit$reciprocal.condition,
      prior.term = training$prior$log.ratio,
      book.errata = paste(
        "For fixed-p score asymptotics, the pooled-covariance variance",
        "contribution is (a'Omega delta)^2 +",
        "(a'Omega a)(delta'Omega delta), not the single term in the draft."
      ),
      regularization = "none",
      no.repair = paste(
        "Strict Cholesky solve; no ridge, pseudoinverse, jitter, or",
        "eigenvalue floor."
      )
    ),
    call = call, data_name = data.name,
    score_scale = "canonical_log_lr", tie = "class1"
  )
}


#' Classical plug-in quadratic discriminant classifier
#'
#' Uses classwise unbiased covariance matrices with divisors `n1 - 1` and
#' `n2 - 1`, then forms the canonical Gaussian QDA log-likelihood ratio.
#' Both covariance matrices must be strictly positive definite. Ordinary
#' positive determinants are used; no pseudo-determinant, absolute value,
#' ridge, or generalized inverse is substituted.
#'
#' @inheritParams classical_lda_classifier
#' @return An `hd_classifier_fit`.
#' @references
#' Anderson, T. W. (2003). *An Introduction to Multivariate Statistical
#' Analysis*, 3rd ed. Wiley.
#' @examples
#' x <- rbind(c(2, 0), c(1, 1), c(1, -1), c(2, 2),
#'            c(-2, 0), c(-1, 2), c(-1, -2), c(-2, -1))
#' classical_qda_classifier(x, rep(c("A", "B"), each = 4))
#' @export
classical_qda_classifier <- function(
    x, y, prior = "equal", strict = TRUE) {
  call <- match.call()
  data.name <- deparse1(substitute(x))
  strict <- .clf_validate_strict(strict)
  training <- .clf_prepare_xy(x, y, prior = prior, min_class = 2L)
  method <- "Classical plug-in Gaussian quadratic discriminant analysis"
  core <- .clf_two_class_moments(training)
  failed <- .clf_classical_failure(
    core, method, training, strict, call, "canonical_log_lr", data.name
  )
  if (!is.null(failed)) return(failed)
  covariance1 <- core$sscp1 / (training$n1 - 1L)
  covariance2 <- core$sscp2 / (training$n2 - 1L)
  fit1 <- tryCatch(
    .clf_validate_spd(covariance1, "class-1 covariance", training$p),
    error = identity
  )
  fit2 <- tryCatch(
    .clf_validate_spd(covariance2, "class-2 covariance", training$p),
    error = identity
  )
  if (inherits(fit1, "condition") || inherits(fit2, "condition")) {
    detail <- paste(
      if (inherits(fit1, "condition")) conditionMessage(fit1) else NULL,
      if (inherits(fit2, "condition")) conditionMessage(fit2) else NULL,
      collapse = "; "
    )
    return(.clf_failure(
      paste0(
        "Classical QDA requires both class covariance matrices to be ",
        "strictly positive definite: ", detail,
        " No ridge or generalized inverse was used."
      ),
      strict, method, training, call, "canonical_log_lr",
      "class covariance factorization", data_name = data.name
    ))
  }
  omega1 <- chol2inv(fit1$chol)
  omega2 <- chol2inv(fit2$chol)
  quadratic <- (omega2 - omega1) / 2
  linear <- as.numeric(omega1 %*% core$mean1 - omega2 %*% core$mean2)
  intercept <- -sum(core$mean1 * as.numeric(omega1 %*% core$mean1)) / 2 +
    sum(core$mean2 * as.numeric(omega2 %*% core$mean2)) / 2 -
    fit1$log.determinant / 2 + fit2$log.determinant / 2 +
    training$prior$log.ratio
  if (any(!is.finite(quadratic)) || any(!is.finite(linear)) ||
      !is.finite(intercept)) {
    return(.clf_failure(
      paste(
        "The classical QDA score is outside the finite double range;",
        "no clipping was applied."
      ),
      strict, method, training, call, "canonical_log_lr",
      "score construction", data_name = data.name
    ))
  }
  estimate <- .clf_name_estimates(list(
    location1 = core$mean1, location2 = core$mean2,
    covariance1 = covariance1, covariance2 = covariance2,
    precision1 = omega1, precision2 = omega2
  ), training)
  .clf_new_fit(
    method = method, training = training,
    score_model = list(
      type = "quadratic", quadratic = quadratic,
      linear = linear, intercept = intercept
    ),
    estimate = estimate,
    diagnostics = list(
      covariance.divisor = c(class1 = training$n1 - 1L,
                             class2 = training$n2 - 1L),
      log.covariance.determinant = c(
        class1 = fit1$log.determinant, class2 = fit2$log.determinant
      ),
      reciprocal.condition = c(
        class1 = fit1$reciprocal.condition,
        class2 = fit2$reciprocal.condition
      ),
      book.errata = c(
        paste(
          "The direct D/beta QDA display is twice the canonical score",
          "returned by this function."
        ),
        paste(
          "In diverging dimension, operator-norm covariance consistency",
          "alone does not control log determinants; eigenvalue and",
          "log-determinant conditions are also required."
        )
      ),
      regularization = "none",
      determinant = "ordinary signed positive-definite determinant",
      no.repair = paste(
        "Strict Cholesky solves; no ridge, pseudoinverse, jitter,",
        "eigenvalue floor, or absolute determinant."
      )
    ),
    call = call, data_name = data.name,
    score_scale = "canonical_log_lr", tie = "class1"
  )
}


#' Diagonal-covariance independence classifier
#'
#' Fits the Gaussian independence rule after replacing the common covariance
#' by its pooled diagonal. `variance_divisor = "unbiased"` uses `n - 2`;
#' `"mle"` uses `n`. Zero marginal variances are never floored. They either
#' cause a failure or are dropped explicitly.
#'
#' @inheritParams classical_lda_classifier
#' @param variance_divisor Pooled marginal-variance convention.
#' @param zero_variance Whether an exactly non-positive pooled marginal
#'   variance is an error or is explicitly dropped.
#'
#' @return An `hd_classifier_fit`.
#' @references
#' Bickel, P. J. and Levina, E. (2004). Some theory for Fisher's linear
#' discriminant function, naive Bayes, and some alternatives when there are
#' many more variables than observations. *Bernoulli*, 10, 989--1010.
#' @examples
#' x <- rbind(c(2, 1), c(1, 2), c(3, 0),
#'            c(-2, -1), c(-1, -2), c(-3, 0))
#' independence_classifier(x, rep(c("A", "B"), each = 3))
#' @export
independence_classifier <- function(
    x, y, prior = "equal",
    variance_divisor = c("unbiased", "mle"),
    zero_variance = c("error", "drop"), strict = TRUE) {
  call <- match.call()
  data.name <- deparse1(substitute(x))
  strict <- .clf_validate_strict(strict)
  variance.divisor <- match.arg(variance_divisor)
  zero.variance <- match.arg(zero_variance)
  training <- .clf_prepare_xy(x, y, prior = prior, min_class = 2L)
  method <- "Diagonal-covariance Gaussian independence classifier"
  core <- .clf_two_class_moments(training)
  failed <- .clf_classical_failure(
    core, method, training, strict, call, "canonical_log_lr", data.name
  )
  if (!is.null(failed)) return(failed)
  divisor <- if (variance.divisor == "unbiased") training$n - 2L else
    training$n
  variances <- diag(core$sscp1 + core$sscp2) / divisor
  bad <- !is.finite(variances) | variances <= 0
  if (any(bad) && zero.variance == "error") {
    return(.clf_failure(
      paste(
        "The independence rule encountered a non-positive or non-finite",
        "pooled marginal variance; no variance floor was applied."
      ),
      strict, method, training, call, "canonical_log_lr",
      "marginal variance", data_name = data.name
    ))
  }
  keep <- which(!bad)
  if (!length(keep)) {
    return(.clf_failure(
      "No positive finite marginal variance remains after explicit dropping.",
      strict, method, training, call, "canonical_log_lr",
      "feature selection", data_name = data.name
    ))
  }
  delta <- core$mean1 - core$mean2
  direction <- numeric(training$p)
  direction[keep] <- delta[keep] / variances[keep]
  midpoint <- (core$mean1 + core$mean2) / 2
  intercept <- -sum(midpoint * direction) + training$prior$log.ratio
  if (any(!is.finite(direction)) || !is.finite(intercept)) {
    return(.clf_failure(
      paste(
        "The independence direction is outside the finite double range;",
        "no clipping was applied."
      ),
      strict, method, training, call, "canonical_log_lr",
      "score construction", data_name = data.name
    ))
  }
  estimate <- .clf_name_estimates(list(
    location1 = core$mean1, location2 = core$mean2,
    marginal.variance = variances, mean.difference = delta,
    direction = direction
  ), training)
  .clf_new_fit(
    method = method, training = training,
    score_model = list(
      type = "linear", coefficients = direction, intercept = intercept
    ),
    estimate = estimate,
    tuning = list(
      variance.divisor = variance.divisor,
      zero.variance = zero.variance
    ),
    diagnostics = list(
      covariance.divisor = divisor,
      selected.features = keep,
      dropped.features = which(bad),
      prior.term = training$prior$log.ratio,
      covariance.model = "diagonal common covariance",
      regularization = "none",
      no.repair = paste(
        "No variance floor, ridge, pseudoinverse, or imputation was used;",
        "dropping occurs only when explicitly requested."
      )
    ),
    call = call, data_name = data.name,
    score_scale = "canonical_log_lr", tie = "class1"
  )
}

.clf_fair_selection <- function(core, training, usable, t.statistics,
                                selection, m, t_threshold) {
  ordering <- usable[order(-abs(t.statistics[usable]), usable)]
  if (selection == "m") {
    if (is.null(m)) {
      stop("`m` must be supplied when `selection = \"m\"`.", call. = FALSE)
    }
    m <- .clf_positive_integer(m, "m")
    if (m > length(ordering)) {
      stop("`m` cannot exceed the number of usable features.", call. = FALSE)
    }
    selected <- ordering[seq_len(m)]
    return(list(
      selected = selected, ordering = ordering, selected.m = m,
      criterion = NULL, lambda.max = NULL, threshold = NULL
    ))
  }
  if (selection == "threshold") {
    threshold <- as.numeric(t_threshold)
    if (length(threshold) != 1L || is.na(threshold) ||
        !is.finite(threshold) || threshold < 0) {
      stop("`t_threshold` must be one finite non-negative number.",
           call. = FALSE)
    }
    selected <- ordering[abs(t.statistics[ordering]) > threshold]
    return(list(
      selected = selected, ordering = ordering,
      selected.m = length(selected), criterion = NULL,
      lambda.max = NULL, threshold = threshold
    ))
  }

  pooled.sscp <- core$sscp1 + core$sscp2
  marginal.sscp <- diag(pooled.sscp)[ordering]
  denominator <- sqrt(outer(marginal.sscp, marginal.sscp))
  correlation <- pooled.sscp[ordering, ordering, drop = FALSE] / denominator
  correlation <- (correlation + t(correlation)) / 2
  diag(correlation) <- 1
  paper <- cpp_ch5_classical_fair_criterion(
    t.statistics[ordering]^2, correlation, training$n1, training$n2
  )
  criterion <- as.numeric(paper$criterion)
  lambda.max <- as.numeric(paper$lambda_max)
  if (length(criterion) != length(ordering) || anyNA(criterion) ||
      any(!is.finite(criterion))) {
    stop("The primary FAIR feature-count criterion is not finite.",
         call. = FALSE)
  }
  selected.m <- which.max(criterion)
  list(
    selected = ordering[seq_len(selected.m)], ordering = ordering,
    selected.m = selected.m, criterion = criterion,
    lambda.max = lambda.max, threshold = NULL,
    correlation = correlation
  )
}


#' Features annealed independence rule (FAIR)
#'
#' Implements Fan and Fan's FAIR classifier. Features are ranked by the
#' absolute Welch two-sample statistic
#' \deqn{T_j=(\bar X_j-\bar Y_j)/
#' \sqrt{S_{1j}^2/n_1+S_{2j}^2/n_2},}
#' while classification uses the marginal variance
#' \eqn{(S_{1j}^2+S_{2j}^2)/2}. Feature-count selection by `"paper"`
#' maximizes equation (4.3), including the largest eigenvalue of every
#' truncated pooled within-class correlation matrix. Ties in absolute
#' t-statistics are broken by the original feature index, and a tie in the
#' criterion selects the smallest feature count.
#'
#' FAIR is implemented under its equal-prior common-covariance contract. An
#' exactly zero marginal variance is never floored: it is either an error or
#' is dropped explicitly.
#'
#' @inheritParams classical_lda_classifier
#' @param selection Select an explicit feature count, an explicit absolute
#'   t-statistic threshold, or the primary paper criterion.
#' @param m Positive feature count used only by `selection = "m"`.
#' @param t_threshold Non-negative threshold used only by
#'   `selection = "threshold"`; the comparison is strict.
#' @param zero_variance Whether non-positive marginal variances are errors or
#'   are explicitly dropped.
#'
#' @return An `hd_classifier_fit`.
#' @references Fan, J. and Fan, Y. (2008). High-dimensional classification
#'   using features annealed independence rules. *Annals of Statistics*, 36,
#'   2605--2637. \doi{10.1214/07-AOS504}.
#' @examples
#' x <- rbind(c(3, 1, 0), c(2, 2, 1), c(4, 0, -1),
#'            c(-3, -1, 0), c(-2, -2, -1), c(-4, 0, 1))
#' fair_classifier(x, rep(c("A", "B"), each = 3), selection = "m", m = 1)
#' @export
fair_classifier <- function(
    x, y, selection = c("paper", "m", "threshold"),
    m = NULL, t_threshold = NULL, prior = "equal",
    zero_variance = c("error", "drop"), strict = TRUE) {
  call <- match.call()
  data.name <- deparse1(substitute(x))
  strict <- .clf_validate_strict(strict)
  selection <- match.arg(selection)
  zero.variance <- match.arg(zero_variance)
  if (selection != "m" && !is.null(m)) {
    stop("`m` is used only when `selection = \"m\"`.", call. = FALSE)
  }
  if (selection != "threshold" && !is.null(t_threshold)) {
    stop(
      "`t_threshold` is used only when `selection = \"threshold\"`.",
      call. = FALSE
    )
  }
  training <- .clf_prepare_xy(
    x, y, prior = prior, equal_prior = TRUE, min_class = 2L
  )
  method <- "Fan--Fan features annealed independence rule (FAIR)"
  core <- .clf_two_class_moments(training)
  failed <- .clf_classical_failure(
    core, method, training, strict, call, "method_threshold", data.name
  )
  if (!is.null(failed)) return(failed)
  variance1 <- diag(core$sscp1) / (training$n1 - 1L)
  variance2 <- diag(core$sscp2) / (training$n2 - 1L)
  classification.variance <- (variance1 + variance2) / 2
  bad <- !is.finite(classification.variance) | classification.variance <= 0
  if (any(bad) && zero.variance == "error") {
    return(.clf_failure(
      paste(
        "FAIR encountered a non-positive or non-finite marginal variance;",
        "no variance floor was applied."
      ),
      strict, method, training, call, "method_threshold",
      "marginal variance", data_name = data.name
    ))
  }
  usable <- which(!bad)
  if (!length(usable)) {
    return(.clf_failure(
      "FAIR has no usable feature after explicit zero-variance dropping.",
      strict, method, training, call, "method_threshold",
      "feature screening", data_name = data.name
    ))
  }
  delta <- core$mean1 - core$mean2
  welch.denominator <- sqrt(
    variance1 / training$n1 + variance2 / training$n2
  )
  t.statistics <- rep.int(NA_real_, training$p)
  t.statistics[usable] <- delta[usable] / welch.denominator[usable]
  if (any(!is.finite(t.statistics[usable]))) {
    return(.clf_failure(
      "The FAIR Welch statistics are not finite; no clipping was applied.",
      strict, method, training, call, "method_threshold",
      "Welch statistics", data_name = data.name
    ))
  }
  selected.fit <- tryCatch(
    .clf_fair_selection(
      core, training, usable, t.statistics, selection, m, t_threshold
    ),
    error = identity
  )
  if (inherits(selected.fit, "condition")) {
    if (selection %in% c("m", "threshold")) stop(selected.fit)
    return(.clf_failure(
      paste0(
        "The primary FAIR selection criterion failed: ",
        conditionMessage(selected.fit), " No spectral repair was applied."
      ),
      strict, method, training, call, "method_threshold",
      "paper feature-count criterion", data_name = data.name
    ))
  }
  selected <- selected.fit$selected
  direction <- numeric(training$p)
  if (length(selected)) {
    direction[selected] <- delta[selected] /
      classification.variance[selected]
  }
  midpoint <- (core$mean1 + core$mean2) / 2
  intercept <- -sum(midpoint * direction)
  if (any(!is.finite(direction)) || !is.finite(intercept)) {
    return(.clf_failure(
      "The FAIR score is outside the finite double range; no clipping was applied.",
      strict, method, training, call, "method_threshold",
      "score construction", data_name = data.name
    ))
  }
  estimate <- .clf_name_estimates(list(
    location1 = core$mean1, location2 = core$mean2,
    marginal.variance = classification.variance,
    mean.difference = delta, t.statistics = t.statistics,
    direction = direction
  ), training)
  .clf_new_fit(
    method = method, training = training,
    score_model = list(
      type = "linear", coefficients = direction, intercept = intercept
    ),
    estimate = estimate,
    tuning = list(
      selection = selection, m = selected.fit$selected.m,
      t.threshold = selected.fit$threshold,
      criterion = selected.fit$criterion,
      truncated.lambda.max = selected.fit$lambda.max
    ),
    diagnostics = list(
      variance.formula = "(unbiased class-1 variance + unbiased class-2 variance) / 2",
      welch.formula = "difference / sqrt(S1^2/n1 + S2^2/n2)",
      selected.features = selected,
      dropped.features = which(bad),
      feature.order = selected.fit$ordering,
      feature.tie = "original feature index",
      criterion.tie = "smallest m",
      threshold.comparison = "absolute t statistic > t_threshold",
      prior.contract = "equal",
      regularization = "none",
      no.repair = paste(
        "No variance floor, ridge, pseudoinverse, eigenvalue floor, or",
        "imputation was used; dropping occurs only when requested."
      )
    ),
    call = call, data_name = data.name,
    score_scale = "method_threshold", tie = "class1"
  )
}


.clf_nonnegative_scalar <- function(x, name) {
  x <- as.numeric(x)
  if (length(x) != 1L || is.na(x) || !is.finite(x) || x < 0) {
    stop(sprintf("`%s` must be one finite non-negative number.", name),
         call. = FALSE)
  }
  x
}


.clf_shao_core <- function(x, class, M_cov, M_mean, alpha) {
  n <- nrow(x)
  p <- ncol(x)
  moments <- tryCatch(
    cpp_ch5_classical_two_class_moments(x, class), error = identity
  )
  if (inherits(moments, "condition")) {
    return(list(valid = FALSE, failure = conditionMessage(moments),
                stage = "two-class moments"))
  }
  covariance <- (moments$sscp1 + moments$sscp2) / n
  covariance.threshold <- M_cov * sqrt(log(p) / n)
  thresholded <- tryCatch(
    cpp_ch5_classical_hard_threshold_covariance(
      covariance, covariance.threshold
    ),
    error = identity
  )
  if (inherits(thresholded, "condition")) {
    return(list(valid = FALSE, failure = conditionMessage(thresholded),
                stage = "covariance thresholding"))
  }
  covariance.fit <- tryCatch(
    .clf_validate_spd(
      as.matrix(thresholded), "thresholded covariance", p
    ),
    error = identity
  )
  if (inherits(covariance.fit, "condition")) {
    return(list(
      valid = FALSE, failure = conditionMessage(covariance.fit),
      stage = "thresholded covariance factorization",
      covariance = covariance, thresholded.covariance = thresholded,
      covariance.threshold = covariance.threshold
    ))
  }
  delta <- as.numeric(moments$mean1 - moments$mean2)
  mean.threshold <- M_mean * (log(p) / n)^alpha
  thresholded.delta <- delta * (abs(delta) > mean.threshold)
  direction <- as.numeric(
    .clf_solve_chol(covariance.fit$chol, thresholded.delta)
  )
  midpoint <- as.numeric(moments$mean1 + moments$mean2) / 2
  intercept <- -sum(midpoint * direction)
  if (any(!is.finite(direction)) || !is.finite(intercept)) {
    return(list(
      valid = FALSE,
      failure = "The thresholded LDA score is outside the finite double range.",
      stage = "score construction"
    ))
  }
  list(
    valid = TRUE, mean1 = as.numeric(moments$mean1),
    mean2 = as.numeric(moments$mean2), covariance = covariance,
    thresholded.covariance = as.matrix(covariance.fit$matrix),
    covariance.threshold = covariance.threshold,
    delta = delta, thresholded.delta = thresholded.delta,
    mean.threshold = mean.threshold, direction = direction,
    intercept = intercept,
    reciprocal.condition = covariance.fit$reciprocal.condition
  )
}


.clf_shao_grid <- function(parameter_grid) {
  if (is.matrix(parameter_grid)) {
    if (!is.numeric(parameter_grid) || ncol(parameter_grid) != 2L) {
      stop("A matrix `parameter_grid` must have two numeric columns.",
           call. = FALSE)
    }
    parameter_grid <- as.data.frame(parameter_grid)
    names(parameter_grid) <- c("M_cov", "M_mean")
  }
  if (!is.data.frame(parameter_grid) ||
      !all(c("M_cov", "M_mean") %in% names(parameter_grid)) ||
      nrow(parameter_grid) < 1L) {
    stop(
      "`parameter_grid` must have at least one row and columns `M_cov` and `M_mean`.",
      call. = FALSE
    )
  }
  grid <- data.frame(
    M_cov = as.numeric(parameter_grid$M_cov),
    M_mean = as.numeric(parameter_grid$M_mean)
  )
  if (anyNA(grid) || any(!is.finite(as.matrix(grid))) ||
      any(as.matrix(grid) < 0)) {
    stop("Both columns of `parameter_grid` must be finite and non-negative.",
         call. = FALSE)
  }
  if (anyDuplicated(grid)) {
    stop("`parameter_grid` must not contain duplicate parameter pairs.",
         call. = FALSE)
  }
  grid
}


.clf_shao_loo <- function(training, grid, alpha) {
  correct <- integer(nrow(grid))
  valid <- rep.int(TRUE, nrow(grid))
  fold.failure <- vector("list", nrow(grid))
  for (g in seq_len(nrow(grid))) {
    failure <- character()
    for (i in seq_len(training$n)) {
      keep <- seq_len(training$n) != i
      core <- .clf_shao_core(
        training$x[keep, , drop = FALSE], training$class[keep],
        grid$M_cov[g], grid$M_mean[g], alpha
      )
      if (!isTRUE(core$valid)) {
        valid[g] <- FALSE
        failure <- c(failure, sprintf("row %d: %s", i, core$failure))
        break
      }
      score <- cpp_ch5_classical_linear_scores(
        training$x[i, , drop = FALSE], core$direction, core$intercept
      )
      predicted <- if (score[1L] >= 0) 1L else 2L
      correct[g] <- correct[g] + as.integer(predicted == training$class[i])
    }
    fold.failure[[g]] <- failure
  }
  list(correct = correct, valid = valid, fold.failure = fold.failure)
}


#' Shao--Wang--Deng--Wang thresholded sparse LDA
#'
#' Implements the primary hard-threshold LDA benchmark. The pooled covariance
#' uses divisor `n`; only its off-diagonal entries are retained when
#' \eqn{|s_{ij}|>M_{cov}\sqrt{\log(p)/n}}, while the diagonal is always
#' retained. Mean differences are retained when
#' \eqn{|\bar X_j-\bar Y_j|>M_{mean}\{\log(p)/n\}^{\alpha}}, with
#' \eqn{0<\alpha<1/2}. Equality is thresholded out in both cases.
#'
#' Supply either both constants or an explicit two-column grid. A grid is
#' selected by exact leave-one-out correct classification count, recomputing
#' all moments and thresholds in every fold. Candidates whose thresholded
#' covariance is not positive definite in any fold are ineligible. A score
#' tie selects class 1; a tuning tie selects the first supplied grid row.
#'
#' @inheritParams classical_lda_classifier
#' @param M_cov,M_mean Explicit non-negative threshold constants.
#' @param alpha Exponent strictly between zero and one-half.
#' @param parameter_grid Optional data frame or numeric matrix with columns
#'   `M_cov` and `M_mean`. No default grid is invented.
#' @param cv The primary leave-one-out scheme; currently only `"loo"`.
#'
#' @return An `hd_classifier_fit`.
#' @references Shao, J., Wang, Y., Deng, X., and Wang, S. (2011). Sparse
#'   linear discriminant analysis by thresholding for high dimensional data.
#'   *Annals of Statistics*, 39, 1241--1265.
#' @examples
#' x <- rbind(c(3, 1), c(2, 2), c(4, -0.5),
#'            c(-3, -1), c(-2, -2), c(-4, 0.5))
#' shao_threshold_lda(
#'   x, rep(c("A", "B"), each = 3), M_cov = 0, M_mean = 0,
#'   alpha = 0.25
#' )
#' @export
shao_threshold_lda <- function(
    x, y, M_cov = NULL, M_mean = NULL, alpha,
    parameter_grid = NULL, cv = "loo", prior = "equal",
    strict = TRUE) {
  call <- match.call()
  data.name <- deparse1(substitute(x))
  strict <- .clf_validate_strict(strict)
  alpha <- as.numeric(alpha)
  if (length(alpha) != 1L || is.na(alpha) || !is.finite(alpha) ||
      alpha <= 0 || alpha >= 0.5) {
    stop("`alpha` must be one finite number strictly between 0 and 0.5.",
         call. = FALSE)
  }
  if (length(cv) != 1L || !is.character(cv) || is.na(cv) || cv != "loo") {
    stop("`cv` must be \"loo\".", call. = FALSE)
  }
  using.grid <- !is.null(parameter_grid)
  if (using.grid && (!is.null(M_cov) || !is.null(M_mean))) {
    stop("Supply either `parameter_grid` or `M_cov` and `M_mean`, not both.",
         call. = FALSE)
  }
  if (!using.grid && (is.null(M_cov) || is.null(M_mean))) {
    stop(
      "Supply both `M_cov` and `M_mean`, or an explicit `parameter_grid`.",
      call. = FALSE
    )
  }
  training <- .clf_prepare_xy(
    x, y, prior = prior, equal_prior = TRUE,
    min_class = if (using.grid) 3L else 2L
  )
  method <- "Shao--Wang--Deng--Wang thresholded sparse LDA"
  if (using.grid) {
    grid <- .clf_shao_grid(parameter_grid)
    cv.fit <- .clf_shao_loo(training, grid, alpha)
    eligible <- which(cv.fit$valid)
    if (!length(eligible)) {
      return(.clf_failure(
        paste(
          "Every supplied threshold pair failed strict positive-definite",
          "validation in leave-one-out fitting; no ridge was used."
        ),
        strict, method, training, call, "method_threshold",
        "leave-one-out tuning", tuning = list(
          grid = grid, correct = cv.fit$correct,
          eligible = cv.fit$valid, fold.failure = cv.fit$fold.failure
        ), data_name = data.name
      ))
    }
    best.correct <- max(cv.fit$correct[eligible])
    chosen <- eligible[which(cv.fit$correct[eligible] == best.correct)[1L]]
    M.cov <- grid$M_cov[chosen]
    M.mean <- grid$M_mean[chosen]
    tuning <- list(
      source = "explicit grid with leave-one-out classification",
      M.cov = M.cov, M.mean = M.mean, alpha = alpha,
      grid = grid, correct = cv.fit$correct, eligible = cv.fit$valid,
      selected.row = chosen, tie.rule = "first supplied grid row",
      fold.failure = cv.fit$fold.failure
    )
  } else {
    M.cov <- .clf_nonnegative_scalar(M_cov, "M_cov")
    M.mean <- .clf_nonnegative_scalar(M_mean, "M_mean")
    tuning <- list(
      source = "specified", M.cov = M.cov, M.mean = M.mean,
      alpha = alpha
    )
  }
  core <- .clf_shao_core(
    training$x, training$class, M.cov, M.mean, alpha
  )
  if (!isTRUE(core$valid)) {
    return(.clf_failure(
      paste0(
        "Thresholded sparse LDA failed at ", core$stage, ": ",
        core$failure, " No ridge or generalized inverse was used."
      ),
      strict, method, training, call, "method_threshold", core$stage,
      tuning = tuning, data_name = data.name
    ))
  }
  if (!is.null(training$feature.names)) {
    for (name in c("mean1", "mean2", "delta", "thresholded.delta",
                   "direction")) {
      names(core[[name]]) <- training$feature.names
    }
    for (name in c("covariance", "thresholded.covariance")) {
      dimnames(core[[name]]) <- list(
        training$feature.names, training$feature.names
      )
    }
  }
  tuning$covariance.threshold <- core$covariance.threshold
  tuning$mean.threshold <- core$mean.threshold
  .clf_new_fit(
    method = method, training = training,
    score_model = list(
      type = "linear", coefficients = core$direction,
      intercept = core$intercept
    ),
    estimate = list(
      location1 = core$mean1, location2 = core$mean2,
      covariance = core$covariance,
      thresholded.covariance = core$thresholded.covariance,
      mean.difference = core$delta,
      thresholded.mean.difference = core$thresholded.delta,
      direction = core$direction
    ),
    tuning = tuning,
    diagnostics = list(
      covariance.divisor = training$n,
      covariance.threshold.comparison = "off-diagonal absolute value > threshold; diagonal retained",
      mean.threshold.comparison = "absolute value > threshold",
      reciprocal.condition = core$reciprocal.condition,
      prior.contract = "equal",
      cv = if (using.grid) "leave-one-out with all nuisance estimates refit" else
        "not used; constants specified",
      regularization = "none",
      no.repair = paste(
        "No ridge, pseudoinverse, jitter, eigenvalue floor, or",
        "post-threshold positive-definite repair was used."
      )
    ),
    call = call, data_name = data.name,
    score_scale = "method_threshold", tie = "class1"
  )
}

