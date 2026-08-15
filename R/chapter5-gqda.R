.c5gqda_flag <- function(value, name) {
  if (!is.logical(value) || length(value) != 1L || is.na(value)) {
    stop(sprintf("`%s` must be TRUE or FALSE.", name), call. = FALSE)
  }
  value
}

.c5gqda_positive_integer <- function(value, name, minimum = 1L) {
  value <- as.numeric(value)
  if (length(value) != 1L || is.na(value) || !is.finite(value) ||
      value != floor(value) || value < minimum ||
      value > .Machine$integer.max) {
    stop(sprintf("`%s` must be one integer not smaller than %d.",
                 name, minimum), call. = FALSE)
  }
  as.integer(value)
}

.c5gqda_unit_interval <- function(value, name = "c") {
  value <- as.numeric(value)
  if (length(value) != 1L || is.na(value) || !is.finite(value) ||
      value < 0 || value > 1) {
    stop(sprintf("`%s` must be one finite number in [0, 1].", name),
         call. = FALSE)
  }
  value
}

.c5gqda_align_location <- function(value, training, name) {
  original.names <- names(value)
  value <- as.numeric(value)
  if (length(value) != training$p || anyNA(value) || any(!is.finite(value))) {
    stop(sprintf("`%s` must be a finite vector of length %d.",
                 name, training$p), call. = FALSE)
  }
  if (!is.null(training$feature.names) && !is.null(original.names)) {
    if (!setequal(original.names, training$feature.names)) {
      stop(sprintf("Names on `%s` do not match the training variables.", name),
           call. = FALSE)
    }
    value <- value[match(training$feature.names, original.names)]
  }
  names(value) <- training$feature.names
  value
}

.c5gqda_align_square <- function(value, training, name) {
  if (!is.matrix(value) || !is.numeric(value) ||
      !identical(dim(value), base::c(training$p, training$p))) {
    stop(sprintf("`%s` must be a numeric %d by %d matrix.",
                 name, training$p, training$p), call. = FALSE)
  }
  storage.mode(value) <- "double"
  if (anyNA(value) || any(!is.finite(value))) {
    stop(sprintf("`%s` must contain only finite values.", name),
         call. = FALSE)
  }
  variable.names <- training$feature.names
  if (!is.null(variable.names)) {
    rn <- rownames(value)
    cn <- colnames(value)
    if (xor(is.null(rn), is.null(cn))) {
      stop(sprintf("`%s` must have both row and column names, or neither.",
                   name), call. = FALSE)
    }
    if (!is.null(rn)) {
      if (!setequal(rn, variable.names) || !setequal(cn, variable.names)) {
        stop(sprintf("Dimnames on `%s` do not match the training variables.",
                     name), call. = FALSE)
      }
      value <- value[match(variable.names, rn), match(variable.names, cn),
                     drop = FALSE]
    }
  }
  asymmetry <- max(abs(value - t(value)))
  matrix.scale <- max(1, max(abs(value)))
  symmetry.tolerance <- 512 * .Machine$double.eps * training$p * matrix.scale
  if (asymmetry > symmetry.tolerance) {
    stop(sprintf("`%s` is not numerically symmetric.", name), call. = FALSE)
  }
  value <- value / 2 + t(value) / 2
  dimnames(value) <- list(variable.names, variable.names)
  value
}

.c5gqda_spd <- function(value, training, name) {
  value <- .c5gqda_align_square(value, training, name)
  factor <- tryCatch(chol(value), error = identity)
  if (inherits(factor, "condition")) {
    stop(sprintf("`%s` must be strictly positive definite; no repair was used.",
                 name), call. = FALSE)
  }
  reciprocal.condition <- rcond(value)
  if (!is.finite(reciprocal.condition) || reciprocal.condition <= 0) {
    stop(sprintf("`%s` has no certified positive reciprocal condition.", name),
         call. = FALSE)
  }
  list(
    matrix = value,
    factor = factor,
    log.determinant = 2 * sum(log(diag(factor))),
    rcond = reciprocal.condition
  )
}

.c5gqda_extract_class_fit <- function(fit, training, label) {
  if (!is.list(fit)) {
    stop(sprintf("%s must be a certified list-like class fit.", label),
         call. = FALSE)
  }
  if (!is.null(fit$valid) && !isTRUE(fit$valid)) {
    stop(sprintf("%s is marked invalid.", label), call. = FALSE)
  }
  location <- fit$location
  if (is.null(location) && is.list(fit$estimate)) {
    location <- fit$estimate$location
  }
  scatter <- fit$scatter
  if (is.null(scatter) && is.list(fit$estimate)) {
    scatter <- fit$estimate$scatter
  }
  if (is.null(location) || is.null(scatter)) {
    stop(sprintf("%s must provide certified `location` and `scatter` fields.",
                 label), call. = FALSE)
  }
  location <- .c5gqda_align_location(location, training,
                                      paste0(label, "$location"))
  scatter.info <- .c5gqda_spd(scatter, training,
                               paste0(label, "$scatter"))

  precision <- NULL
  precision.source <- "Cholesky inverse of certified scatter"
  if (inherits(fit, "high_dimensional_hr_fit")) {
    if (!is.numeric(fit$scatter.scale) || length(fit$scatter.scale) != 1L ||
        !is.finite(fit$scatter.scale) || fit$scatter.scale <= 0 ||
        is.null(fit$precision)) {
      stop(sprintf("%s lacks the primary QDA scatter scale.", label),
           call. = FALSE)
    }
    precision <- as.matrix(fit$precision) / fit$scatter.scale
    precision.source <- "HDHR trace-scaled shape precision"
  } else if (!is.null(fit$scatter.precision)) {
    precision <- fit$scatter.precision
    precision.source <- "supplied certified scatter precision"
  } else if (!is.null(fit$precision)) {
    precision <- fit$precision
    precision.source <- "supplied precision checked against scatter"
  }
  if (is.null(precision)) {
    precision <- chol2inv(scatter.info$factor)
  }
  precision.info <- .c5gqda_spd(precision, training,
                                 paste0(label, "$precision"))
  inverse.residual <- max(abs(
    scatter.info$matrix %*% precision.info$matrix - diag(training$p)
  ))
  inverse.tolerance <- min(
    1e-5,
    4096 * .Machine$double.eps * training$p /
      min(scatter.info$rcond, precision.info$rcond)
  )
  if (!is.finite(inverse.residual) || inverse.residual > inverse.tolerance) {
    stop(sprintf(
      "%s has inconsistent scatter and precision (residual %.3g).",
      label, inverse.residual
    ), call. = FALSE)
  }
  list(
    location = location,
    scatter = scatter.info$matrix,
    precision = precision.info$matrix,
    log.determinant = scatter.info$log.determinant,
    scatter.rcond = scatter.info$rcond,
    precision.rcond = precision.info$rcond,
    inverse.residual = inverse.residual,
    inverse.tolerance = inverse.tolerance,
    precision.source = precision.source,
    original.fit = fit
  )
}

.c5gqda_classical_fit <- function(x, training, label) {
  if (nrow(x) < 2L) {
    stop(sprintf("%s needs at least two observations.", label),
         call. = FALSE)
  }
  scatter <- stats::cov(x)
  location <- colMeans(x)
  list(
    valid = TRUE,
    location = location,
    scatter = scatter,
    diagnostics = list(
      source = "classical unbiased sample covariance",
      divisor = nrow(x) - 1L
    )
  )
}

.c5gqda_components <- function(newdata, fits) {
  answer <- cpp_ch5_gqda_components(
    newdata,
    fits[[1L]]$location,
    fits[[2L]]$location,
    fits[[1L]]$precision,
    fits[[2L]]$precision
  )
  if (length(answer$difference) != nrow(newdata) ||
      anyNA(answer$difference) || any(!is.finite(answer$difference))) {
    stop(paste(
      "The GQDA Mahalanobis-distance difference is not finitely",
      "representable; no clipping or rescaling was applied."
    ), call. = FALSE)
  }
  answer
}

.c5gqda_select_c <- function(difference, logdet, truth) {
  difference <- as.numeric(difference)
  logdet <- rep_len(as.numeric(logdet), length(difference))
  truth <- as.integer(truth)
  if (!length(difference) || length(truth) != length(difference) ||
      anyNA(difference) || any(!is.finite(difference)) ||
      anyNA(logdet) || any(!is.finite(logdet)) ||
      any(!truth %in% 1:2)) {
    stop("The GQDA tuning path requires finite scores and binary truths.",
         call. = FALSE)
  }
  tolerance <- 512 * .Machine$double.eps * pmax(1, abs(logdet))
  active <- abs(logdet) > tolerance
  ratios <- difference[active] / logdet[active]
  breakpoints <- ratios[is.finite(ratios) & ratios >= 0 & ratios <= 1]
  knots <- sort(unique(base::c(0, breakpoints, 1)))
  midpoints <- if (length(knots) > 1L) {
    knots[-length(knots)] / 2 + knots[-1L] / 2
  } else {
    numeric()
  }
  candidates <- sort(unique(base::c(knots, midpoints)))
  errors <- vapply(candidates, function(candidate) {
    prediction <- ifelse(difference >= candidate * logdet, 1L, 2L)
    sum(prediction != truth)
  }, numeric(1))
  selected <- which(errors == min(errors))[[1L]]
  list(
    c = candidates[[selected]],
    identified = any(active),
    candidates = candidates,
    errors = errors,
    error.rate = errors / length(truth),
    breakpoints = sort(unique(breakpoints)),
    tie.rule = "smallest c among equal empirical error counts",
    inequality = "difference >= c * signed log-determinant contrast"
  )
}

.c5gqda_with_seed <- function(seed, expression) {
  if (is.null(seed)) return(force(expression))
  seed <- .c5gqda_positive_integer(seed, "seed", minimum = 0L)
  had.seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had.seed) old.seed <- get(".Random.seed", envir = .GlobalEnv)
  on.exit({
    if (had.seed) {
      assign(".Random.seed", old.seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv,
                      inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  set.seed(seed)
  force(expression)
}

.c5gqda_folds <- function(training, folds, seed) {
  folds <- .c5gqda_positive_integer(folds, "folds", minimum = 2L)
  if (folds > min(training$n1, training$n2)) {
    stop("`folds` cannot exceed the smaller class size.", call. = FALSE)
  }
  assign.class <- function(index) {
    order <- if (is.null(seed)) index else sample(index, length(index))
    answer <- integer(length(index))
    answer[match(order, index)] <- rep(seq_len(folds), length.out = length(index))
    answer
  }
  assignment <- .c5gqda_with_seed(seed, {
    answer <- integer(training$n)
    answer[training$class1] <- assign.class(training$class1)
    answer[training$class2] <- assign.class(training$class2)
    answer
  })
  list(
    assignment = assignment,
    folds = folds,
    seed = seed,
    rule = if (is.null(seed))
      "deterministic stratified round-robin in input order" else
      "seeded stratified permutation with RNG state restored"
  )
}

.c5gqda_cross_validated_components <- function(training, fold.info,
                                                fit.function) {
  difference <- rep(NA_real_, training$n)
  logdet <- rep(NA_real_, training$n)
  fold.diagnostics <- vector("list", fold.info$folds)
  for (fold in seq_len(fold.info$folds)) {
    validation <- which(fold.info$assignment == fold)
    fitting <- which(fold.info$assignment != fold)
    first <- fit.function(
      training$x[intersect(fitting, training$class1), , drop = FALSE], 1L
    )
    second <- fit.function(
      training$x[intersect(fitting, training$class2), , drop = FALSE], 2L
    )
    fits <- list(
      .c5gqda_extract_class_fit(first, training, "fold class-1 fit"),
      .c5gqda_extract_class_fit(second, training, "fold class-2 fit")
    )
    components <- .c5gqda_components(
      training$x[validation, , drop = FALSE], fits
    )
    difference[validation] <- components$difference
    logdet[validation] <-
      fits[[1L]]$log.determinant - fits[[2L]]$log.determinant
    fold.diagnostics[[fold]] <- list(
      validation = validation,
      fitting = fitting,
      log.determinant.contrast = unique(logdet[validation]),
      class.fit = lapply(fits, function(value) {
        list(
          scatter.rcond = value$scatter.rcond,
          precision.rcond = value$precision.rcond,
          inverse.residual = value$inverse.residual
        )
      })
    )
  }
  if (anyNA(difference) || anyNA(logdet)) {
    stop("Cross-validation did not score every training observation.",
         call. = FALSE)
  }
  list(difference = difference, logdet = logdet,
       diagnostics = fold.diagnostics)
}

.c5gqda_invalid <- function(message, strict, training, method, call,
                            data.name, failure.stage) {
  if (strict) stop(message, call. = FALSE)
  warning(message, call. = FALSE)
  .clf_new_fit(
    method = method,
    training = training,
    score_model = list(
      type = "custom",
      fun = function(newdata) rep(NA_real_, nrow(newdata)),
      description = "invalid GQDA fit"
    ),
    estimate = list(),
    tuning = list(),
    diagnostics = list(
      failure.stage = failure.stage,
      failure = message,
      no.repair = paste(
        "No ridge, eigenvalue floor, pseudoinverse, jitter, determinant",
        "absolute value, or covariance projection was applied."
      )
    ),
    call = call,
    primary_orientation =
      "positive score selects class1; invalid fit has no usable score",
    score_scale = "method_threshold",
    tie = "class1",
    data_name = data.name,
    valid = FALSE
  )
}

.c5gqda_construct <- function(training, fits, selection, requested.c,
                              tuning.components = NULL, fold.info = NULL,
                              fold.diagnostics = NULL,
                              method, call, data.name, strict) {
  tryCatch({
    logdet <- fits[[1L]]$log.determinant - fits[[2L]]$log.determinant
    logdet.tolerance <- 512 * .Machine$double.eps *
      max(1, abs(fits[[1L]]$log.determinant),
          abs(fits[[2L]]$log.determinant))
    logdet.identified <- abs(logdet) > logdet.tolerance

    if (selection == "fixed") {
      selected <- list(
        c = .c5gqda_unit_interval(requested.c),
        identified = logdet.identified,
        candidates = .c5gqda_unit_interval(requested.c),
        errors = NA_real_, error.rate = NA_real_, breakpoints = numeric(),
        tie.rule = "not applicable: c supplied by user",
        inequality = "difference >= c * signed log-determinant contrast"
      )
    } else {
      if (is.null(tuning.components)) {
        tuning.components <- .c5gqda_components(training$x, fits)
        tuning.logdet <- logdet
      } else {
        tuning.logdet <- tuning.components$logdet
      }
      selected <- .c5gqda_select_c(
        tuning.components$difference,
        tuning.logdet,
        ifelse(seq_len(training$n) %in% training$class1, 1L, 2L)
      )
    }
    requested.c.original <- selected$c
    if (!logdet.identified) {
      selected$c <- 0
      selected$identified <- FALSE
    }
    resolved.c <- selected$c
    location1 <- fits[[1L]]$location
    location2 <- fits[[2L]]$location
    precision1 <- fits[[1L]]$precision
    precision2 <- fits[[2L]]$precision
    scorer <- local({
      first.location <- location1
      second.location <- location2
      first.precision <- precision1
      second.precision <- precision2
      constant <- resolved.c
      determinant.contrast <- logdet
      function(newdata) {
        components <- cpp_ch5_gqda_components(
          newdata, first.location, second.location,
          first.precision, second.precision
        )
        score <- as.numeric(components$difference) -
          constant * determinant.contrast
        if (anyNA(score) || any(!is.finite(score))) {
          stop("GQDA prediction produced a non-finite score.", call. = FALSE)
        }
        score
      }
    })
    score.model <- list(
      type = "custom",
      fun = scorer,
      description = paste(
        "squared Mahalanobis distance to class2 minus class1, minus",
        "c times log(det(scatter1)/det(scatter2))"
      )
    )
    expanded <- list(
      quadratic = precision2 - precision1,
      linear = 2 * (precision1 %*% location1 -
                      precision2 %*% location2),
      intercept = as.numeric(
        crossprod(location2, precision2 %*% location2) -
          crossprod(location1, precision1 %*% location1) -
          resolved.c * logdet
      )
    )
    fit <- .clf_new_fit(
      method = method,
      training = training,
      score_model = score.model,
      estimate = list(
        location = list(class1 = location1, class2 = location2),
        scatter = list(class1 = fits[[1L]]$scatter,
                       class2 = fits[[2L]]$scatter),
        precision = list(class1 = precision1,
                         class2 = precision2),
        expanded.quadratic.score = expanded
      ),
      tuning = list(
        selection = selection,
        c = resolved.c,
        c.before.zero.logdet.identification = requested.c.original,
        c.identified = logdet.identified,
        path = selected,
        folds = fold.info
      ),
      diagnostics = list(
        failure.stage = NULL,
        failure = NULL,
        log.determinant = base::c(
          class1 = fits[[1L]]$log.determinant,
          class2 = fits[[2L]]$log.determinant
        ),
        signed.log.determinant.contrast = logdet,
        log.determinant.zero.tolerance = logdet.tolerance,
        class.fit = lapply(fits, function(value) {
          list(
            scatter.rcond = value$scatter.rcond,
            precision.rcond = value$precision.rcond,
            inverse.residual = value$inverse.residual,
            inverse.tolerance = value$inverse.tolerance,
            precision.source = value$precision.source
          )
        }),
        cross.validation = fold.diagnostics,
        primary.rule = paste(
          "class1 iff Delta_d^2 >= c log{det(scatter1)/det(scatter2)};",
          "c=0 is MMD and c=1 is Gaussian QDA under equal priors"
        ),
        primary.errata = base::c(
          paste(
            "The book's Gaussian log-likelihood HR-QDA display is not the",
            "Yan--Feng--Zhang primary classifier; the primary uses GQDA."
          ),
          paste(
            "The published sorted-ratio selector reverses direction when",
            "the signed log-determinant contrast is negative. This",
            "implementation evaluates the original inequality at all",
            "breakpoints, adjacent midpoints, and endpoints."
          ),
          paste(
            "The published separated-interval upper endpoint r1(n1) is",
            "logically r1(1); direct empirical-loss enumeration avoids it."
          )
        ),
        no.repair = paste(
          "All scatters and precisions are certified SPD and mutually",
          "inverse. No ridge, floor, pseudoinverse, jitter, nearest-PD",
          "projection, or absolute determinant was used."
        )
      ),
      call = call,
      primary_orientation = paste(
        "positive Delta_d^2 - c log{det(scatter1)/det(scatter2)}",
        "selects class1"
      ),
      score_scale = "method_threshold",
      tie = "class1",
      data_name = data.name,
      valid = TRUE
    )
    fit$c <- resolved.c
    fit$location <- list(class1 = location1, class2 = location2)
    fit$scatter <- list(class1 = fits[[1L]]$scatter,
                        class2 = fits[[2L]]$scatter)
    fit$precision <- list(class1 = precision1,
                          class2 = precision2)
    fit$log.determinant.contrast <- logdet
    fit$class.fits <- lapply(fits, `[[`, "original.fit")
    fit
  }, error = function(error) {
    .c5gqda_invalid(
      paste0("GQDA fitting failed: ", conditionMessage(error)),
      strict, training, method, call, data.name, "GQDA construction"
    )
  })
}

#' Generalized quadratic discriminant analysis for elliptical populations
#'
#' Fits the Bose--Pal--Saha Ray--Nayak generalized QDA rule using classical
#' class moments or two explicitly supplied certified location--scatter fits.
#' The rule assigns an observation to class 1 when
#' \deqn{\Delta_d^2(x) \ge c\log\{|S_1|/|S_2|\},\qquad 0\le c\le1.}
#' A supplied `c` is used directly. Otherwise the resubstitution selector (or
#' an explicitly requested leakage-free classical cross-validation selector)
#' enumerates both endpoints, every legal decision breakpoint, and adjacent
#' midpoints. This direct loss calculation remains correct for positive,
#' negative, and zero signed log-determinant contrasts.
#'
#' @param x Numeric training matrix with observations in rows.
#' @param y Two-class response. Factor levels define class 1 then class 2;
#'   otherwise first appearance defines the order.
#' @param class_fits Optional length-two list of certified fits, each providing
#'   finite `location` and strictly positive-definite `scatter`. A compatible
#'   precision may also be supplied and is checked against the scatter.
#' @param c Optional fixed constant in `[0,1]`.
#' @param selection One of `"resubstitution"`, `"fixed"`, or
#'   `"cross_validation"`. Supplying `c` selects `"fixed"`. Cross-validation
#'   is available only for internally refitted classical moments.
#' @param folds Number of stratified folds for `"cross_validation"`.
#' @param seed Optional seed for stratified fold shuffling. The caller's RNG
#'   state is restored. With `NULL`, input-order round-robin folds are used.
#' @param strict If `TRUE`, numerical/certificate failure is an error;
#'   otherwise an invalid, non-predictable fit is returned with a warning.
#'
#' @return An `hd_classifier_fit`. Package scores are positive for class 1.
#'
#' @details This is an equal-prior classifier. The current book's HR-QDA
#'   subsection writes an ordinary Gaussian plug-in QDA score, but the cited
#'   Yan--Feng--Zhang primary applies the generalized threshold above. No
#'   simulation tuning grid, ridge, determinant absolute value, or covariance
#'   repair is used.
#'
#' @references Bose, S., Pal, A., SahaRay, R., and Nayak, J. (2015).
#'   Generalized quadratic discriminant analysis. *Pattern Recognition*,
#'   48(8), 2676--2684. \doi{10.1016/j.patcog.2015.02.016}.
#'
#' @examples
#' x <- rbind(c(-2, 0), c(-1, 1), c(-1, -1), c(-2, 1),
#'            c(2, 0), c(1, 1), c(1, -1), c(2, -1))
#' y <- factor(rep(c("left", "right"), each = 4),
#'             levels = c("left", "right"))
#' fit <- gqda_classifier(x, y, c = 0)
#' predict(fit, x)
#'
#' @export
gqda_classifier <- function(
    x, y, class_fits = NULL, c = NULL,
    selection = base::c("resubstitution", "fixed", "cross_validation"),
    folds = 5L, seed = NULL, strict = TRUE) {
  call <- match.call()
  data.name <- paste(deparse1(substitute(x)), deparse1(substitute(y)),
                     sep = ", ")
  selection.was.missing <- missing(selection)
  selection <- match.arg(selection)
  strict <- .c5gqda_flag(strict, "strict")
  if (!is.null(c)) {
    if (!selection.was.missing && selection != "fixed") {
      stop("Supplying `c` is incompatible with a non-fixed `selection`.",
           call. = FALSE)
    }
    selection <- "fixed"
  } else if (selection == "fixed") {
    stop("`selection = \"fixed\"` requires `c`.", call. = FALSE)
  }
  training <- .clf_prepare_xy(
    x, y, prior = base::c(0.5, 0.5), equal_prior = TRUE, min_class = 2L
  )
  tryCatch({
    supplied <- !is.null(class_fits)
    if (supplied) {
      if (!is.list(class_fits) || length(class_fits) != 2L) {
        stop("`class_fits` must be a length-two list.", call. = FALSE)
      }
      if (selection == "cross_validation") {
        stop(paste(
          "Cross-validation cannot refit arbitrary supplied class fits;",
          "use fixed or resubstitution selection."
        ), call. = FALSE)
      }
      raw.fits <- class_fits
    } else {
      raw.fits <- list(
        .c5gqda_classical_fit(training$x[training$class1, , drop = FALSE],
                              training, "class 1"),
        .c5gqda_classical_fit(training$x[training$class2, , drop = FALSE],
                              training, "class 2")
      )
    }
    fits <- list(
      .c5gqda_extract_class_fit(raw.fits[[1L]], training, "class_fits[[1]]"),
      .c5gqda_extract_class_fit(raw.fits[[2L]], training, "class_fits[[2]]")
    )
    fold.info <- NULL
    fold.diagnostics <- NULL
    tuning.components <- NULL
    if (selection == "cross_validation") {
      fold.info <- .c5gqda_folds(training, folds, seed)
      cv <- .c5gqda_cross_validated_components(
        training, fold.info,
        function(group.x, group) {
          .c5gqda_classical_fit(group.x, training,
                                paste("fold class", group))
        }
      )
      tuning.components <- list(
        difference = cv$difference,
        logdet = cv$logdet
      )
      fold.diagnostics <- cv$diagnostics
    }
    .c5gqda_construct(
      training, fits, selection, c, tuning.components,
      fold.info, fold.diagnostics,
      "Generalized quadratic discriminant analysis",
      call, data.name, strict
    )
  }, error = function(error) {
    .c5gqda_invalid(
      paste0("GQDA fitting failed: ", conditionMessage(error)),
      strict, training, "Generalized quadratic discriminant analysis",
      call, data.name, "class fits or tuning"
    )
  })
}

#' Generalized QDA from certified robust class fits
#'
#' This adapter applies the audited GQDA decision and tuning rule to two
#' certified robust location--scatter fits. It deliberately does not invent a
#' new M-, MVE-, MCD-, S-, or SD-estimation algorithm.
#'
#' @inheritParams gqda_classifier
#'
#' @references Bose, S., Pal, A., SahaRay, R., and Nayak, J. (2015).
#'   Generalized quadratic discriminant analysis. *Pattern Recognition*,
#'   48(8), 2676--2684. \doi{10.1016/j.patcog.2015.02.016}.
#'
#' @return An `hd_classifier_fit`.
#'
#' @examples
#' x <- rbind(c(-2, 0), c(-1, 1), c(-1, -1), c(-2, 1),
#'            c(2, 0), c(1, 1), c(1, -1), c(2, -1))
#' y <- factor(rep(c("left", "right"), each = 4))
#' fits <- lapply(split(seq_len(nrow(x)), y), function(ii) {
#'   list(valid = TRUE, location = colMeans(x[ii, , drop = FALSE]),
#'        scatter = stats::cov(x[ii, , drop = FALSE]))
#' })
#' robust_gqda(x, y, fits, c = 0)
#'
#' @export
robust_gqda <- function(
    x, y, class_fits, c = NULL,
    selection = base::c("resubstitution", "fixed"), strict = TRUE) {
  call <- match.call()
  selection.was.missing <- missing(selection)
  selection <- match.arg(selection)
  if (!is.null(c) && selection.was.missing) selection <- "fixed"
  fit <- gqda_classifier(
    x = x, y = y, class_fits = class_fits, c = c,
    selection = selection, strict = strict
  )
  fit$method <- "Generalized QDA from certified robust class fits"
  fit$call <- call
  fit
}

.c5gqda_hr_pilots <- function(pilot_precision) {
  certified.class <- base::c(
    "spatial_sign_precision_fit", "tensor_spatial_sign_precision_fit"
  )
  if (is.matrix(pilot_precision) ||
      any(vapply(certified.class, function(class.name) {
        inherits(pilot_precision, class.name)
      }, logical(1)))) {
    return(list(pilot_precision, pilot_precision))
  }
  if (!is.list(pilot_precision) || length(pilot_precision) != 2L) {
    stop(paste(
      "`pilot_precision` must be one matrix/certified fit or a length-two",
      "list of class-specific pilots."
    ), call. = FALSE)
  }
  pilot_precision
}

.c5gqda_hr_bandwidth <- function(bandwidth) {
  bandwidth <- as.numeric(bandwidth)
  if (!length(bandwidth) || length(bandwidth) > 2L || anyNA(bandwidth) ||
      any(!is.finite(bandwidth)) || any(bandwidth < 0) ||
      any(bandwidth != floor(bandwidth))) {
    stop("`bandwidth` must contain one or two non-negative integers.",
         call. = FALSE)
  }
  rep_len(as.integer(bandwidth), 2L)
}

#' High-dimensional HR generalized quadratic discriminant analysis
#'
#' Fits the primary Yan--Feng--Zhang classwise high-dimensional HR estimator,
#' including its separate QDA second-moment scale, and plugs the two certified
#' fits into generalized QDA. This is not the ordinary Gaussian HR-QDA score
#' currently displayed in the book.
#'
#' @inheritParams gqda_classifier
#' @param pilot_precision One pilot precision matrix/certified fit shared by
#'   both classes, or a length-two list of class-specific pilots.
#' @param bandwidth One or two hard-banding widths passed to
#'   [high_dimensional_hr()].
#' @param tol,max_iter,median_tol,median_max_iter,zero_tol Iteration controls
#'   passed to every full-sample and cross-validation HR refit.
#'
#' @return An `hd_classifier_fit` containing both certified HR fits.
#'
#' @references Yan, G., Feng, L., and Zhang, X. (2025).
#'   *High-Dimensional Hettmansperger--Randles Estimator and its Applications*.
#'   arXiv:2505.01669. \url{https://arxiv.org/abs/2505.01669}
#'
#' @examples
#' x <- rbind(c(-2, 0), c(-1, 1), c(-1, -1), c(-2, 1),
#'            c(2, 0), c(1, 1), c(1, -1), c(2, -1),
#'            c(-1.5, .4), c(1.5, -.4), c(-1.7, -.3), c(1.7, .3))
#' y <- factor(rep(c("left", "right"), each = 6))
#' hr_gqda(x, y, pilot_precision = diag(2), bandwidth = 0, c = 0)
#'
#' @export
hr_gqda <- function(
    x, y, pilot_precision, bandwidth = 3L, c = NULL,
    selection = base::c("resubstitution", "fixed", "cross_validation"),
    folds = 5L, seed = NULL,
    tol = 1e-8, max_iter = 1000L,
    median_tol = 1e-8, median_max_iter = 1000L,
    zero_tol = 0, strict = TRUE) {
  call <- match.call()
  data.name <- paste(deparse1(substitute(x)), deparse1(substitute(y)),
                     sep = ", ")
  selection.was.missing <- missing(selection)
  selection <- match.arg(selection)
  strict <- .c5gqda_flag(strict, "strict")
  if (!is.null(c)) {
    if (!selection.was.missing && selection != "fixed") {
      stop("Supplying `c` is incompatible with a non-fixed `selection`.",
           call. = FALSE)
    }
    selection <- "fixed"
  } else if (selection == "fixed") {
    stop("`selection = \"fixed\"` requires `c`.", call. = FALSE)
  }
  training <- .clf_prepare_xy(
    x, y, prior = base::c(0.5, 0.5), equal_prior = TRUE, min_class = 2L
  )
  tryCatch({
    pilots <- .c5gqda_hr_pilots(pilot_precision)
    bandwidth <- .c5gqda_hr_bandwidth(bandwidth)
    fit.hr <- function(group.x, group) {
      high_dimensional_hr(
        group.x,
        pilot_precision = pilots[[group]],
        bandwidth = bandwidth[[group]],
        tol = tol, max_iter = max_iter,
        median_tol = median_tol, median_max_iter = median_max_iter,
        zero_tol = zero_tol,
        scale_estimator = "paper_qda",
        strict = TRUE
      )
    }
    raw.fits <- list(
      fit.hr(training$x[training$class1, , drop = FALSE], 1L),
      fit.hr(training$x[training$class2, , drop = FALSE], 2L)
    )
    fits <- list(
      .c5gqda_extract_class_fit(raw.fits[[1L]], training, "class-1 HR fit"),
      .c5gqda_extract_class_fit(raw.fits[[2L]], training, "class-2 HR fit")
    )
    fold.info <- NULL
    fold.diagnostics <- NULL
    tuning.components <- NULL
    if (selection == "cross_validation") {
      fold.info <- .c5gqda_folds(training, folds, seed)
      cv <- .c5gqda_cross_validated_components(
        training, fold.info, fit.hr
      )
      tuning.components <- list(
        difference = cv$difference,
        logdet = cv$logdet
      )
      fold.diagnostics <- cv$diagnostics
    }
    answer <- .c5gqda_construct(
      training, fits, selection, c, tuning.components,
      fold.info, fold.diagnostics,
      "Yan--Feng--Zhang high-dimensional HR GQDA",
      call, data.name, strict
    )
    answer$hr.fits <- raw.fits
    answer$diagnostics$hr <- lapply(raw.fits, function(value) {
      list(
        valid = value$valid,
        bandwidth = value$bandwidth,
        scale.estimator = value$scale.estimator,
        iteration.stable = value$diagnostics$iteration.stable,
        score.residual = value$diagnostics$score.residual,
        shape.equation.residual =
          value$diagnostics$shape.equation.residual
      )
    })
    answer
  }, error = function(error) {
    .c5gqda_invalid(
      paste0("HR-GQDA fitting failed: ", conditionMessage(error)),
      strict, training, "Yan--Feng--Zhang high-dimensional HR GQDA",
      call, data.name, "HDHR class fit or tuning"
    )
  })
}

#' @rdname hr_gqda
#' @param ... Arguments passed unchanged to [hr_gqda()].
#' @export
hr_qda <- function(...) {
  call <- match.call()
  answer <- hr_gqda(...)
  answer$method <- paste0(
    "Yan--Feng--Zhang high-dimensional HR GQDA ",
    "(`hr_qda` compatibility alias)"
  )
  answer$call <- call
  answer$diagnostics$alias <- paste(
    "The cited primary method is generalized QDA with an estimated c;",
    "this alias does not implement the book's ordinary Gaussian QDA display."
  )
  answer
}
