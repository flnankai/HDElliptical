# Chapter 5: sparse quadratic discriminant classifiers.

.c5qda_nonnegative <- function(x, name, allow_null = FALSE) {
  if (is.null(x) && allow_null) return(NULL)
  x <- as.numeric(x)
  if (length(x) != 1L || is.na(x) || !is.finite(x) || x < 0) {
    stop(sprintf("`%s` must be one finite non-negative number.", name),
         call. = FALSE)
  }
  x
}


.c5qda_positive <- function(x, name) {
  x <- .c5qda_nonnegative(x, name)
  if (x == 0) {
    stop(sprintf("`%s` must be strictly positive.", name), call. = FALSE)
  }
  x
}


.c5qda_count <- function(x, name, minimum = 1L) {
  .clf_positive_integer(x, name, minimum = minimum)
}


.c5qda_mle_moments <- function(x, class) {
  x1 <- x[class == 1L, , drop = FALSE]
  x2 <- x[class == 2L, , drop = FALSE]
  mean1 <- colMeans(x1)
  mean2 <- colMeans(x2)
  centered1 <- sweep(x1, 2L, mean1, "-")
  centered2 <- sweep(x2, 2L, mean2, "-")
  list(
    x1 = x1,
    x2 = x2,
    mean1 = mean1,
    mean2 = mean2,
    covariance1 = crossprod(centered1) / nrow(x1),
    covariance2 = crossprod(centered2) / nrow(x2)
  )
}


.c5qda_score <- function(model, x) {
  if (model$type != "quadratic") {
    stop("Internal sparse-QDA score model is not quadratic.", call. = FALSE)
  }
  rowSums((x %*% model$quadratic) * x) +
    drop(x %*% model$linear) + model$intercept
}


.c5qda_model_from_center <- function(quadratic, linear_centered, center,
                                     intercept_centered) {
  quadratic <- (quadratic + t(quadratic)) / 2
  linear_centered <- as.numeric(linear_centered)
  center <- as.numeric(center)
  list(
    type = "quadratic",
    quadratic = quadratic,
    linear = linear_centered - 2 * drop(quadratic %*% center),
    intercept = as.numeric(
      crossprod(center, quadratic %*% center) -
        crossprod(linear_centered, center) + intercept_centered
    )
  )
}


.c5qda_parameter_grid <- function(parameter_grid, first, second) {
  if (is.null(parameter_grid)) return(NULL)
  parameter_grid <- as.data.frame(parameter_grid)
  required <- c(first, second)
  if (!all(required %in% names(parameter_grid))) {
    stop(sprintf(
      "`parameter_grid` must contain columns `%s` and `%s`.", first, second
    ), call. = FALSE)
  }
  answer <- parameter_grid[, required, drop = FALSE]
  answer[[first]] <- as.numeric(answer[[first]])
  answer[[second]] <- as.numeric(answer[[second]])
  if (!nrow(answer) || anyNA(answer) || any(!is.finite(as.matrix(answer))) ||
      any(as.matrix(answer) < 0)) {
    stop("Every tuning-grid entry must be finite and non-negative.",
         call. = FALSE)
  }
  unique(answer)
}


.c5qda_tuning_choice <- function(first_value, second_value, parameter_grid,
                                 first_name, second_name, training, folds,
                                 fit_core) {
  first_value <- .c5qda_nonnegative(first_value, first_name, allow_null = TRUE)
  second_value <- .c5qda_nonnegative(
    second_value, second_name, allow_null = TRUE
  )
  grid <- .c5qda_parameter_grid(parameter_grid, first_name, second_name)
  supplied_pair <- !is.null(first_value) || !is.null(second_value)
  if (supplied_pair && (is.null(first_value) || is.null(second_value))) {
    stop(sprintf("`%s` and `%s` must be supplied together.",
                 first_name, second_name), call. = FALSE)
  }
  if (supplied_pair && !is.null(grid)) {
    stop("Supply either one explicit tuning pair or `parameter_grid`, not both.",
         call. = FALSE)
  }
  if (!supplied_pair && is.null(grid)) {
    stop(
      paste0(
        "No default regularization is invented: supply both `", first_name,
        "` and `", second_name, "`, or provide `parameter_grid`."
      ),
      call. = FALSE
    )
  }
  if (supplied_pair) {
    return(list(
      first = first_value, second = second_value, selection = "specified",
      folds = NULL, cv = NULL
    ))
  }

  folds <- .c5qda_count(folds, "folds", minimum = 2L)
  if (folds > min(training$n1, training$n2)) {
    stop("`folds` cannot exceed the smaller class size.", call. = FALSE)
  }
  fold_id <- integer(training$n)
  fold_id[training$class1] <-
    (seq_along(training$class1) - 1L) %% folds + 1L
  fold_id[training$class2] <-
    (seq_along(training$class2) - 1L) %% folds + 1L
  cv_error <- rep.int(Inf, nrow(grid))
  failures <- integer(nrow(grid))
  for (g in seq_len(nrow(grid))) {
    mistakes <- 0L
    observations <- 0L
    failed <- FALSE
    for (fold in seq_len(folds)) {
      validation <- fold_id == fold
      core <- tryCatch(
        fit_core(
          training$x[!validation, , drop = FALSE],
          training$class[!validation],
          grid[[first_name]][g], grid[[second_name]][g]
        ),
        error = function(error) list(
          valid = FALSE, failure.stage = "cross_validation_fold",
          failure = conditionMessage(error)
        )
      )
      if (!isTRUE(core$valid)) {
        failed <- TRUE
        failures[g] <- failures[g] + 1L
        break
      }
      score <- .c5qda_score(
        core$score.model, training$x[validation, , drop = FALSE]
      )
      predicted <- if (identical(core$tie, "class2")) {
        ifelse(score > 0, 1L, 2L)
      } else {
        ifelse(score >= 0, 1L, 2L)
      }
      mistakes <- mistakes + sum(predicted != training$class[validation])
      observations <- observations + sum(validation)
    }
    if (!failed) cv_error[g] <- mistakes / observations
  }
  ordering <- order(cv_error, grid[[first_name]], grid[[second_name]],
                    seq_len(nrow(grid)))
  winner <- ordering[1L]
  if (!is.finite(cv_error[winner])) {
    return(list(
      first = NA_real_, second = NA_real_, selection = "joint_cv",
      folds = folds,
      cv = cbind(grid, error = cv_error, failed.folds = failures),
      failure = "Every tuning pair failed in deterministic stratified CV."
    ))
  }
  list(
    first = grid[[first_name]][winner],
    second = grid[[second_name]][winner],
    selection = "joint_cv",
    folds = folds,
    cv = cbind(grid, error = cv_error, failed.folds = failures),
    winner = winner
  )
}


.c5qda_eta <- function(raw_score, class) {
  breakpoints <- sort(unique(-as.numeric(raw_score)))
  span <- max(1, abs(breakpoints), diff(range(breakpoints)))
  middle <- if (length(breakpoints) > 1L) {
    (breakpoints[-1L] + breakpoints[-length(breakpoints)]) / 2
  } else {
    numeric()
  }
  candidates <- sort(unique(c(
    breakpoints[1L] - span, breakpoints, middle,
    breakpoints[length(breakpoints)] + span
  )))
  loss <- vapply(candidates, function(eta) {
    mean(ifelse(raw_score + eta > 0, 1L, 2L) != class)
  }, numeric(1))
  winner <- order(loss, abs(candidates), candidates)[1L]
  list(
    eta = candidates[winner],
    error = loss[winner],
    candidates = candidates,
    losses = loss,
    tie_rule = "minimum error, then smallest absolute eta, then smaller eta"
  )
}


.c5qda_failure <- function(message, strict, method, training, call,
                           score_scale, stage, tuning = list(),
                           diagnostics = list()) {
  .clf_failure(
    message = message, strict = strict, method = method,
    training = training, call = call, score_scale = score_scale,
    stage = stage, tuning = tuning, diagnostics = diagnostics,
    data_name = deparse(call$x)
  )
}


.c5qda_li_model <- function(core) {
  quadratic <- core$precision2 - core$precision1
  dhat <- as.numeric(core$thresholded_difference)
  mean1 <- as.numeric(core$mean1)
  centered_linear <- -2 * drop(core$precision2 %*% dhat)
  centered_intercept <- as.numeric(
    crossprod(dhat, core$precision2 %*% dhat) -
      core$logdet1 + core$logdet2
  )
  .c5qda_model_from_center(
    quadratic, centered_linear, mean1, centered_intercept
  )
}


.c5qda_li_fit_core <- function(x, class, thresholds, ridge) {
  core <- cpp_c5qda_li_shao(
    x[class == 1L, , drop = FALSE], x[class == 2L, , drop = FALSE],
    thresholds[1L], thresholds[2L], thresholds[3L], ridge
  )
  if (!isTRUE(core$valid)) {
    return(list(valid = FALSE, failure.stage = core$failure_stage,
                components = core))
  }
  list(
    valid = TRUE,
    score.model = .c5qda_li_model(core),
    components = core,
    tie = "class1"
  )
}


.c5qda_li_bounds <- function(x, class) {
  moments <- .c5qda_mle_moments(x, class)
  off_diagonal <- row(moments$covariance1) != col(moments$covariance1)
  bound3 <- if (any(off_diagonal)) {
    max(abs(c(
      moments$covariance1[off_diagonal],
      moments$covariance2[off_diagonal]
    )))
  } else 0
  c(
    mean = max(abs(moments$mean2 - moments$mean1)),
    difference = max(abs(moments$covariance1 - moments$covariance2)),
    covariance = bound3
  )
}


.c5qda_li_bisection <- function(training, tolerance, ridge) {
  lower <- c(mean = 0, difference = 0, covariance = 0)
  upper <- .c5qda_li_bounds(training$x, training$class)
  max_iterations <- if (max(upper) == 0) 1L else
    max(1L, ceiling(log2(max(upper) / tolerance)) + 1L)
  history <- vector("list", max_iterations)
  winner <- lower
  for (iteration in seq_len(max_iterations)) {
    candidates <- unique(expand.grid(
      mean = c(lower[1L], upper[1L]),
      difference = c(lower[2L], upper[2L]),
      covariance = c(lower[3L], upper[3L]),
      KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
    ))
    error <- rep.int(Inf, nrow(candidates))
    invalid <- integer(nrow(candidates))
    for (candidate in seq_len(nrow(candidates))) {
      mistakes <- 0L
      for (i in seq_len(training$n)) {
        fit <- .c5qda_li_fit_core(
          training$x[-i, , drop = FALSE], training$class[-i],
          as.numeric(candidates[candidate, ]), ridge
        )
        if (!isTRUE(fit$valid)) {
          invalid[candidate] <- invalid[candidate] + 1L
          next
        }
        score <- .c5qda_score(
          fit$score.model, training$x[i, , drop = FALSE]
        )
        mistakes <- mistakes + (ifelse(score >= 0, 1L, 2L) !=
                                   training$class[i])
      }
      if (!invalid[candidate]) error[candidate] <- mistakes / training$n
    }
    choice <- order(error, candidates$mean, candidates$difference,
                    candidates$covariance, seq_len(nrow(candidates)))[1L]
    winner <- as.numeric(candidates[choice, ])
    names(winner) <- names(lower)
    history[[iteration]] <- cbind(
      iteration = iteration, candidates, error = error,
      invalid.leave.one.out = invalid
    )
    if (!is.finite(error[choice])) {
      return(list(
        valid = FALSE,
        failure = "Every Li-Shao endpoint failed during leave-one-out tuning.",
        bounds = upper, history = history[seq_len(iteration)]
      ))
    }
    midpoint <- (lower + upper) / 2
    choose_lower <- winner == lower
    upper[choose_lower] <- midpoint[choose_lower]
    lower[!choose_lower] <- midpoint[!choose_lower]
    if (max(upper - lower) <= tolerance) {
      history <- history[seq_len(iteration)]
      break
    }
  }
  list(
    valid = TRUE, thresholds = winner, lower = lower, upper = upper,
    initial.bounds = .c5qda_li_bounds(training$x, training$class),
    tolerance = tolerance, history = history
  )
}


#' Li-Shao thresholded sparse QDA
#'
#' Fits the three-threshold rule of Li and Shao. Class covariances use the
#' maximum-likelihood divisor \eqn{n_k}. Mean differences are retained only
#' when strictly above the first threshold. Entries whose class difference is
#' at most the second threshold are pooled; off-diagonal entries are then kept
#' only when strictly above the third threshold. Diagonals are not removed.
#'
#' Paper bisection uses deterministic leave-one-out endpoint searches. Its
#' tolerance must be supplied. A ridge is used only when explicitly positive.
#'
#' @param x Numeric training matrix with observations in rows.
#' @param y Two-class response. Its first observed level is class 1.
#' @param threshold_mean,threshold_difference,threshold_covariance Three finite
#'   non-negative thresholds, required with specified selection.
#' @param selection Either \code{"specified"} or \code{"paper_bisection"}.
#' @param bisection_tol Explicit positive bisection stopping tolerance.
#' @param ridge Explicit non-negative diagonal ridge after thresholding.
#' @param strict Whether numerical failure is an error; otherwise an invalid
#'   \code{hd_classifier_fit} is returned with a warning.
#'
#' @return An \code{hd_classifier_fit}; non-negative scores select class 1.
#' @export
#'
#' @references
#' Li, J. and Shao, J. (2015). Sparse quadratic discriminant analysis for
#' high dimensional data. *Statistica Sinica*, 25, 457-473.
#'
#' @examples
#' x <- rbind(
#'   c(-2, 0), c(-1, 1), c(-1, -1), c(-2, 1),
#'   c(2, 0), c(1, 2), c(1, -2), c(2, 1)
#' )
#' y <- factor(rep(c("left", "right"), each = 4))
#' fit <- li_shao_sparse_qda(
#'   x, y, threshold_mean = 0, threshold_difference = 0,
#'   threshold_covariance = 0, ridge = 0.1
#' )
li_shao_sparse_qda <- function(
    x, y, threshold_mean = NULL, threshold_difference = NULL,
    threshold_covariance = NULL,
    selection = c("specified", "paper_bisection"),
    bisection_tol = NULL, ridge = 0, strict = TRUE) {
  call <- match.call()
  selection <- match.arg(selection)
  strict <- .clf_validate_strict(strict)
  ridge <- .c5qda_nonnegative(ridge, "ridge")
  min_class <- if (selection == "paper_bisection") 3L else 2L
  training <- .clf_prepare_xy(
    x, y, prior = "equal", equal_prior = TRUE, min_class = min_class
  )
  tuning_details <- list(selection = selection, ridge = ridge)
  if (selection == "specified") {
    thresholds <- c(
      mean = .c5qda_nonnegative(threshold_mean, "threshold_mean"),
      difference = .c5qda_nonnegative(
        threshold_difference, "threshold_difference"
      ),
      covariance = .c5qda_nonnegative(
        threshold_covariance, "threshold_covariance"
      )
    )
    if (!is.null(bisection_tol)) {
      stop("'bisection_tol' is used only for paper bisection.", call. = FALSE)
    }
  } else {
    if (!all(vapply(
      list(threshold_mean, threshold_difference, threshold_covariance),
      is.null, logical(1)
    ))) {
      stop("Do not supply fixed thresholds with paper bisection.",
           call. = FALSE)
    }
    bisection_tol <- .c5qda_positive(bisection_tol, "bisection_tol")
    tuned <- .c5qda_li_bisection(training, bisection_tol, ridge)
    tuning_details$bisection <- tuned
    if (!isTRUE(tuned$valid)) {
      return(.c5qda_failure(
        tuned$failure, strict, "Li-Shao sparse QDA", training, call,
        "Li-Shao primary Q score", "leave_one_out_bisection",
        tuning = tuning_details
      ))
    }
    thresholds <- tuned$thresholds
  }
  tuning_details$thresholds <- thresholds
  core <- .c5qda_li_fit_core(
    training$x, training$class, thresholds, ridge
  )
  if (!isTRUE(core$valid)) {
    return(.c5qda_failure(
      paste(
        "The thresholded Li-Shao covariance matrix is not invertible and",
        "positive definite under the explicitly requested ridge."
      ),
      strict, "Li-Shao sparse QDA", training, call,
      "Li-Shao primary Q score", core$failure.stage,
      tuning = tuning_details,
      diagnostics = list(no.implicit.ridge = TRUE)
    ))
  }
  components <- core$components
  .clf_new_fit(
    method = "Li-Shao sparse QDA", training = training,
    score_model = core$score.model,
    estimate = list(
      mean1 = as.numeric(components$mean1),
      mean2 = as.numeric(components$mean2),
      thresholded.mean.difference =
        as.numeric(components$thresholded_difference),
      sample.covariance1 = components$sample_covariance1,
      sample.covariance2 = components$sample_covariance2,
      pooled.covariance = components$pooled_covariance,
      covariance1 = components$covariance1,
      covariance2 = components$covariance2,
      precision1 = components$precision1,
      precision2 = components$precision2
    ),
    tuning = tuning_details,
    diagnostics = list(
      covariance.divisor = c(class1 = training$n1, class2 = training$n2),
      log.determinant = c(components$logdet1, components$logdet2),
      no.implicit.ridge = ridge == 0,
      threshold.inequalities = c(
        mean = "strict greater than",
        pooling = "less than or equal",
        covariance = "strict greater than"
      )
    ),
    call = call,
    primary_orientation = "class 2 iff Q < 0; non-negative Q selects class1",
    score_scale = "Li-Shao primary Q score", tie = "class1",
    data_name = deparse(call$x)
  )
}


.c5qda_jiang_core <- function(x, class, lambda_interaction, lambda_linear,
                              rho, solver_tol, solver_max_iter) {
  moments <- .c5qda_mle_moments(x, class)
  matrix_solver <- cpp_c5qda_jiang_matrix(
    moments$covariance1, moments$covariance2, lambda_interaction,
    rho, solver_tol, solver_max_iter
  )
  if (!isTRUE(matrix_solver$converged)) {
    return(list(
      valid = FALSE, failure.stage = matrix_solver$failure_stage,
      matrix.solver = matrix_solver
    ))
  }
  interaction <- (matrix_solver$solution + t(matrix_solver$solution)) / 2
  mean_difference <- moments$mean1 - moments$mean2
  gamma <- 4 * mean_difference +
    drop((moments$covariance1 - moments$covariance2) %*%
           interaction %*% mean_difference)
  vector_solver <- cpp_c5qda_jiang_vector(
    moments$covariance1 + moments$covariance2, gamma, lambda_linear,
    solver_tol, solver_max_iter
  )
  if (!isTRUE(vector_solver$converged)) {
    return(list(
      valid = FALSE, failure.stage = vector_solver$failure_stage,
      matrix.solver = matrix_solver, vector.solver = vector_solver
    ))
  }
  beta <- as.numeric(vector_solver$solution)
  midpoint <- (moments$mean1 + moments$mean2) / 2
  centered <- sweep(x, 2L, midpoint, "-")
  raw_score <- rowSums((centered %*% interaction) * centered) +
    drop(centered %*% beta)
  eta_fit <- .c5qda_eta(raw_score, class)
  list(
    valid = TRUE,
    score.model = .c5qda_model_from_center(
      interaction, beta, midpoint, eta_fit$eta
    ),
    tie = "class2",
    mean1 = moments$mean1, mean2 = moments$mean2,
    covariance1 = moments$covariance1,
    covariance2 = moments$covariance2,
    interaction = interaction, beta = beta, gamma = gamma,
    midpoint = midpoint, eta = eta_fit,
    matrix.solver = matrix_solver, vector.solver = vector_solver
  )
}


#' Jiang-Wang-Leng direct sparse QDA
#'
#' Implements the direct quadratic-discriminant estimator of Jiang, Wang and
#' Leng. The interaction estimate minimizes a penalized quadratic trace loss,
#' is symmetrized after optimization, and the linear coefficient minimizes its
#' lasso quadratic loss. This is deliberately not the Dantzig program sometimes
#' misattributed to this paper. The intercept exhaustively checks sorted raw
#' training-score breakpoints and intervening intervals for minimum 0-1 loss.
#'
#' @inheritParams li_shao_sparse_qda
#' @param lambda_interaction,lambda_linear Explicit non-negative penalties.
#'   Supply both, or instead provide \code{parameter_grid}.
#' @param parameter_grid Paired \code{lambda_interaction} and
#'   \code{lambda_linear} columns for deterministic stratified joint CV.
#' @param folds Number of folds used only with a parameter grid.
#' @param rho Positive ADMM penalty for the interaction loss.
#' @param solver_tol Positive optimization and certificate tolerance.
#' @param solver_max_iter Positive optimization iteration limit.
#'
#' @return An \code{hd_classifier_fit}. Strictly positive scores select class 1;
#'   exact zero selects class 2, following the primary decision rule.
#' @export
#'
#' @references
#' Jiang, B., Wang, X. and Leng, C. (2018). A direct approach for sparse
#' quadratic discriminant analysis. *Journal of Machine Learning Research*,
#' 19(31), 1-37.
#'
#' @examples
#' x <- rbind(
#'   c(-2, 0), c(-1, 1), c(-1, -1), c(-2, 1),
#'   c(2, 0), c(1, 2), c(1, -2), c(2, 1)
#' )
#' y <- factor(rep(c("left", "right"), each = 4))
#' fit <- jiang_da_qda(x, y, 0.2, 0.2, solver_tol = 1e-6)
jiang_da_qda <- function(
    x, y, lambda_interaction = NULL, lambda_linear = NULL,
    parameter_grid = NULL, folds = 5L, rho = 1,
    solver_tol = 1e-7, solver_max_iter = 10000L, strict = TRUE) {
  call <- match.call()
  strict <- .clf_validate_strict(strict)
  rho <- .c5qda_positive(rho, "rho")
  solver_tol <- .c5qda_positive(solver_tol, "solver_tol")
  solver_max_iter <- .c5qda_count(
    solver_max_iter, "solver_max_iter", minimum = 1L
  )
  grid_mode <- is.null(lambda_interaction) && is.null(lambda_linear) &&
    !is.null(parameter_grid)
  training <- .clf_prepare_xy(
    x, y, prior = "equal", equal_prior = TRUE,
    min_class = if (grid_mode) 3L else 2L
  )
  fit_core <- function(z, group, first, second) {
    .c5qda_jiang_core(
      z, group, first, second, rho, solver_tol, solver_max_iter
    )
  }
  chosen <- .c5qda_tuning_choice(
    lambda_interaction, lambda_linear, parameter_grid,
    "lambda_interaction", "lambda_linear", training, folds, fit_core
  )
  if (!is.null(chosen$failure)) {
    return(.c5qda_failure(
      chosen$failure, strict, "Jiang-Wang-Leng direct sparse QDA",
      training, call, "Jiang primary quadratic score",
      "joint_cross_validation", tuning = chosen
    ))
  }
  core <- fit_core(training$x, training$class, chosen$first, chosen$second)
  if (!isTRUE(core$valid)) {
    return(.c5qda_failure(
      "The penalized Jiang direct-QDA optimization did not certify convergence.",
      strict, "Jiang-Wang-Leng direct sparse QDA", training, call,
      "Jiang primary quadratic score", core$failure.stage,
      tuning = chosen,
      diagnostics = list(
        matrix.solver = core$matrix.solver,
        vector.solver = core$vector.solver
      )
    ))
  }
  .clf_new_fit(
    method = "Jiang-Wang-Leng direct sparse QDA", training = training,
    score_model = core$score.model,
    estimate = list(
      mean1 = core$mean1, mean2 = core$mean2,
      covariance1 = core$covariance1, covariance2 = core$covariance2,
      interaction = core$interaction, linear = core$beta,
      gamma = core$gamma, midpoint = core$midpoint, eta = core$eta$eta
    ),
    tuning = c(chosen, list(
      lambda.interaction = chosen$first,
      lambda.linear = chosen$second, rho = rho
    )),
    diagnostics = list(
      matrix.solver = core$matrix.solver,
      vector.solver = core$vector.solver,
      eta = core$eta,
      interaction.symmetry.residual =
        max(abs(core$interaction - t(core$interaction))),
      estimator = "penalized quadratic and lasso losses; not Dantzig"
    ),
    call = call,
    primary_orientation = "strictly positive score selects class1",
    score_scale = "Jiang primary quadratic score", tie = "class2",
    data_name = deparse(call$x)
  )
}


.c5qda_ssqda_class_moments <- function(
    x, median_tol, median_max_iter, zero_tol) {
  if (nrow(x) < 3L) {
    return(list(
      valid = FALSE, failure.stage = "ssqda_class_size",
      failure = "SSQDA requires at least three observations in each fitted class."
    ))
  }
  location <- spatial_median(
    x, tol = median_tol, max_iter = median_max_iter,
    zero_tol = zero_tol, warn = FALSE
  )
  median_diagnostics <- lapply(
    c("objective", "iterations", "converged",
      "relative_change", "equation_residual"),
    function(name) attr(location, name)
  )
  names(median_diagnostics) <- c(
    "objective", "iterations", "converged",
    "relative.change", "equation.residual"
  )
  if (!isTRUE(median_diagnostics$converged)) {
    return(list(
      valid = FALSE, failure.stage = "spatial_median",
      failure = "The SSQDA spatial median did not certify convergence.",
      median = median_diagnostics
    ))
  }
  robust <- cpp_c5qda_ssqda_moments(x, location, zero_tol)
  if (!isTRUE(robust$valid)) {
    return(list(
      valid = FALSE, failure.stage = robust$failure_stage,
      failure = paste(
        "SSQDA has a zero spatial residual; its SSCM is undefined under",
        "the requested zero tolerance."
      ),
      median = median_diagnostics, robust = robust
    ))
  }
  robust$signs <- NULL
  list(
    valid = TRUE, mean = as.numeric(location),
    covariance = robust$covariance,
    sscm = robust$sscm, trace = robust$trace_estimate,
    zero.residuals = robust$zero_residuals,
    trace.algorithm = robust$trace_algorithm,
    median = median_diagnostics
  )
}


.c5qda_dantzig_core <- function(
    x, class, lambda_D, lambda_beta, solver_tol, feasibility_tol,
    solver_max_iter, robust = FALSE, median_tol = NULL,
    median_max_iter = NULL, zero_tol = NULL) {
  if (robust) {
    class1 <- .c5qda_ssqda_class_moments(
      x[class == 1L, , drop = FALSE],
      median_tol, median_max_iter, zero_tol
    )
    if (!isTRUE(class1$valid)) {
      return(list(
        valid = FALSE, failure.stage = class1$failure.stage,
        failure = class1$failure, robust.class1 = class1
      ))
    }
    class2 <- .c5qda_ssqda_class_moments(
      x[class == 2L, , drop = FALSE],
      median_tol, median_max_iter, zero_tol
    )
    if (!isTRUE(class2$valid)) {
      return(list(
        valid = FALSE, failure.stage = class2$failure.stage,
        failure = class2$failure, robust.class1 = class1,
        robust.class2 = class2
      ))
    }
    mean1 <- class1$mean
    mean2 <- class2$mean
    covariance1 <- class1$covariance
    covariance2 <- class2$covariance
  } else {
    moments <- .c5qda_mle_moments(x, class)
    mean1 <- moments$mean1
    mean2 <- moments$mean2
    covariance1 <- moments$covariance1
    covariance2 <- moments$covariance2
    class1 <- class2 <- NULL
  }
  matrix_solver <- cpp_c5qda_dantzig_matrix(
    covariance1, covariance2, lambda_D, solver_tol, solver_max_iter
  )
  matrix_certified <- isTRUE(matrix_solver$converged) &&
    is.finite(matrix_solver$feasible_maxnorm) &&
    matrix_solver$feasible_maxnorm <= lambda_D + feasibility_tol
  if (!matrix_certified) {
    return(list(
      valid = FALSE, failure.stage = "matrix_dantzig_solver",
      failure = paste(
        "The matrix-free Dantzig interaction solver did not certify",
        "convergence and post-symmetrization feasibility."
      ),
      matrix.solver = matrix_solver, robust.class1 = class1,
      robust.class2 = class2
    ))
  }
  difference <- mean2 - mean1
  vector_solver <- cpp_c5qda_dantzig_vector(
    covariance2, difference, lambda_beta,
    solver_tol, solver_max_iter
  )
  vector_certified <- isTRUE(vector_solver$converged) &&
    is.finite(vector_solver$feasible_maxnorm) &&
    vector_solver$feasible_maxnorm <= lambda_beta + feasibility_tol
  if (!vector_certified) {
    return(list(
      valid = FALSE, failure.stage = "vector_dantzig_solver",
      failure = "The vector Dantzig solver did not certify feasibility.",
      matrix.solver = matrix_solver, vector.solver = vector_solver,
      robust.class1 = class1, robust.class2 = class2
    ))
  }
  interaction <- matrix_solver$solution
  beta <- as.numeric(vector_solver$solution)
  determinant <- cpp_c5qda_signed_logdet(
    diag(ncol(x)) + interaction %*% covariance1
  )
  if (!isTRUE(determinant$valid)) {
    return(list(
      valid = FALSE, failure.stage = determinant$failure_stage,
      failure = paste(
        "The signed determinant of I + D S1 is not strictly positive;",
        "no absolute value or numerical repair is applied."
      ),
      matrix.solver = matrix_solver, vector.solver = vector_solver,
      determinant = determinant, robust.class1 = class1,
      robust.class2 = class2
    ))
  }
  midpoint <- (mean1 + mean2) / 2
  centered_linear <- -2 * beta
  centered_intercept <- as.numeric(
    2 * crossprod(beta, midpoint - mean1) -
      determinant$log_determinant
  )
  list(
    valid = TRUE,
    score.model = .c5qda_model_from_center(
      interaction, centered_linear, mean1, centered_intercept
    ),
    tie = "class1",
    mean1 = mean1, mean2 = mean2,
    covariance1 = covariance1, covariance2 = covariance2,
    interaction = interaction, beta = beta, midpoint = midpoint,
    determinant = determinant,
    matrix.solver = matrix_solver, vector.solver = vector_solver,
    robust.class1 = class1, robust.class2 = class2
  )
}


.c5qda_sparse_dantzig_fit <- function(
    training, call, method, score_scale, lambda_D, lambda_beta,
    parameter_grid, folds, solver_tol, feasibility_tol, solver_max_iter,
    strict, robust = FALSE, median_tol = NULL, median_max_iter = NULL,
    zero_tol = NULL) {
  fit_core <- function(z, group, first, second) {
    .c5qda_dantzig_core(
      z, group, first, second, solver_tol, feasibility_tol,
      solver_max_iter, robust, median_tol, median_max_iter, zero_tol
    )
  }
  chosen <- .c5qda_tuning_choice(
    lambda_D, lambda_beta, parameter_grid,
    "lambda_D", "lambda_beta", training, folds, fit_core
  )
  if (!is.null(chosen$failure)) {
    return(.c5qda_failure(
      chosen$failure, strict, method, training, call, score_scale,
      "joint_cross_validation", tuning = chosen
    ))
  }
  core <- fit_core(training$x, training$class, chosen$first, chosen$second)
  if (!isTRUE(core$valid)) {
    return(.c5qda_failure(
      core$failure, strict, method, training, call, score_scale,
      core$failure.stage, tuning = chosen,
      diagnostics = list(
        matrix.solver = core$matrix.solver,
        vector.solver = core$vector.solver,
        signed.determinant = core$determinant,
        robust.class1 = core$robust.class1,
        robust.class2 = core$robust.class2
      )
    ))
  }
  .clf_new_fit(
    method = method, training = training, score_model = core$score.model,
    estimate = list(
      mean1 = core$mean1, mean2 = core$mean2,
      covariance1 = core$covariance1, covariance2 = core$covariance2,
      interaction = core$interaction, beta = core$beta,
      midpoint = core$midpoint,
      sscm1 = if (robust) core$robust.class1$sscm else NULL,
      sscm2 = if (robust) core$robust.class2$sscm else NULL
    ),
    tuning = c(chosen, list(
      lambda.D = chosen$first, lambda.beta = chosen$second,
      solver.tolerance = solver_tol,
      feasibility.tolerance = feasibility_tol
    )),
    diagnostics = list(
      matrix.solver = core$matrix.solver,
      vector.solver = core$vector.solver,
      signed.determinant = core$determinant,
      determinant.matrix = "I + D S1",
      matrix.free = isTRUE(core$matrix.solver$matrix_free),
      feasibility.rechecked.after.symmetrization = TRUE,
      equal.prior = TRUE,
      robust.class1 = core$robust.class1,
      robust.class2 = core$robust.class2,
      no.numerical.repair = TRUE
    ),
    call = call,
    primary_orientation = "non-negative twice-canonical score selects class1",
    score_scale = score_scale, tie = "class1",
    data_name = deparse(call$x)
  )
}


#' Sparse discriminant analysis with quadratic interactions
#'
#' Fits the SDAR rule through two Dantzig programs. The interaction operator is
#' applied directly as one half of \eqn{S_1 D S_2 + S_2 D S_1}; no Kronecker
#' matrix is formed. The interaction is symmetrized and its feasibility is then
#' checked again. The method has an equal-prior contract. The signed
#' determinant of \eqn{I + D S_1} must be positive.
#'
#' @inheritParams jiang_da_qda
#' @param lambda_D,lambda_beta Explicit non-negative Dantzig bounds.
#' @param parameter_grid Paired \code{lambda_D} and \code{lambda_beta} columns
#'   for deterministic stratified joint cross-validation.
#' @param feasibility_tol Positive tolerance for final constraint feasibility.
#'
#' @return An \code{hd_classifier_fit}; non-negative scores select class 1.
#' @export
#'
#' @references
#' Cai, T. T. and Zhang, L. (2021). A convex optimization approach to
#' high-dimensional sparse quadratic discriminant analysis. *Annals of
#' Statistics*, 49, 1537-1568.
#'
#' @examples
#' x <- rbind(
#'   c(-2, 0), c(-1, 1), c(-1, -1), c(-2, 1),
#'   c(2, 0), c(1, 2), c(1, -2), c(2, 1)
#' )
#' y <- factor(rep(c("left", "right"), each = 4))
#' fit <- sdar_qda(x, y, lambda_D = 0.5, lambda_beta = 0.5,
#'                 solver_tol = 1e-6)
sdar_qda <- function(
    x, y, lambda_D = NULL, lambda_beta = NULL,
    parameter_grid = NULL, folds = 5L, solver_tol = 1e-7,
    feasibility_tol = 1e-7, solver_max_iter = 10000L, strict = TRUE) {
  call <- match.call()
  strict <- .clf_validate_strict(strict)
  solver_tol <- .c5qda_positive(solver_tol, "solver_tol")
  feasibility_tol <- .c5qda_positive(feasibility_tol, "feasibility_tol")
  solver_max_iter <- .c5qda_count(
    solver_max_iter, "solver_max_iter", minimum = 1L
  )
  grid_mode <- is.null(lambda_D) && is.null(lambda_beta) &&
    !is.null(parameter_grid)
  training <- .clf_prepare_xy(
    x, y, prior = "equal", equal_prior = TRUE,
    min_class = if (grid_mode) 3L else 2L
  )
  .c5qda_sparse_dantzig_fit(
    training, call, "SDAR sparse quadratic discriminant analysis",
    "SDAR twice-canonical Q score", lambda_D, lambda_beta,
    parameter_grid, folds, solver_tol, feasibility_tol,
    solver_max_iter, strict
  )
}


#' Spatial-sign sparse quadratic discriminant analysis
#'
#' Fits SSQDA with classwise spatial medians and covariance surrogates equal to
#' a trace estimate times the spatial-sign covariance matrix. The ordered
#' triple-U trace is evaluated through its exact \eqn{O(np)} identity,
#' \eqn{\sum_i \lVert X_i-\bar X\rVert^2/(n-1)}. Each class must contain at
#' least three observations. The Dantzig programs, determinant certificate,
#' equal-prior contract, and score orientation are the same as in SDAR.
#'
#' @inheritParams sdar_qda
#' @param median_tol Positive spatial-median equation tolerance.
#' @param median_max_iter Positive spatial-median iteration limit.
#' @param zero_tol Non-negative threshold for a zero spatial residual. No
#'   perturbation is made when such a residual occurs.
#'
#' @return An \code{hd_classifier_fit}; non-negative scores select class 1.
#' @export
#'
#' @references
#' Feng, L. (2025). Spatial sign based sparse quadratic discriminant analysis
#' for high-dimensional elliptical distributions. arXiv:2504.11187.
#'
#' @examples
#' x <- rbind(
#'   c(-3, -0.5), c(-2.2, 1.1), c(-1.4, -1.3),
#'   c(-2.7, 1.8), c(-0.9, 0.4), c(-1.8, -2),
#'   c(2.8, 0.2), c(1.9, 1.7), c(1.2, -1.8),
#'   c(2.5, 2.2), c(0.7, -0.2), c(1.6, -2.4)
#' )
#' y <- factor(rep(c("left", "right"), each = 6))
#' fit <- ssqda(x, y, lambda_D = 1, lambda_beta = 1,
#'              solver_tol = 1e-6)
ssqda <- function(
    x, y, lambda_D = NULL, lambda_beta = NULL,
    parameter_grid = NULL, folds = 10L, median_tol = 1e-8,
    median_max_iter = 1000L, zero_tol = 0, solver_tol = 1e-7,
    feasibility_tol = 1e-7, solver_max_iter = 10000L, strict = TRUE) {
  call <- match.call()
  strict <- .clf_validate_strict(strict)
  median_tol <- .c5qda_positive(median_tol, "median_tol")
  median_max_iter <- .c5qda_count(
    median_max_iter, "median_max_iter", minimum = 1L
  )
  zero_tol <- .c5qda_nonnegative(zero_tol, "zero_tol")
  solver_tol <- .c5qda_positive(solver_tol, "solver_tol")
  feasibility_tol <- .c5qda_positive(feasibility_tol, "feasibility_tol")
  solver_max_iter <- .c5qda_count(
    solver_max_iter, "solver_max_iter", minimum = 1L
  )
  grid_mode <- is.null(lambda_D) && is.null(lambda_beta) &&
    !is.null(parameter_grid)
  training <- .clf_prepare_xy(
    x, y, prior = "equal", equal_prior = TRUE,
    min_class = if (grid_mode) 4L else 3L
  )
  .c5qda_sparse_dantzig_fit(
    training, call, "Spatial-sign sparse QDA",
    "SSQDA twice-canonical Q score", lambda_D, lambda_beta,
    parameter_grid, folds, solver_tol, feasibility_tol,
    solver_max_iter, strict, robust = TRUE,
    median_tol = median_tol, median_max_iter = median_max_iter,
    zero_tol = zero_tol
  )
}


.c5qda_location <- function(x, p, name) {
  x <- as.numeric(x)
  if (length(x) != p || anyNA(x) || any(!is.finite(x))) {
    stop(sprintf("'%s' must contain one finite value per feature.", name),
         call. = FALSE)
  }
  x
}


#' Certified sparse plug-in QDA
#'
#' Constructs a quadratic discriminant rule from explicitly supplied class
#' means and covariance matrices. Both covariance matrices must pass strict
#' symmetry, Cholesky, and reciprocal-condition certificates. The function
#' never estimates or repairs a covariance, making it suitable for plugging in
#' externally obtained sparse positive-definite estimators.
#'
#' @param x Numeric training matrix used to define features and class levels.
#' @param y Two-class response.
#' @param mean1,mean2 Explicit finite class mean vectors.
#' @param covariance1,covariance2 Explicit symmetric positive-definite class
#'   covariance matrices.
#' @param prior Equal, empirical, or explicitly supplied positive class prior.
#' @param strict Whether covariance-certificate failure is an error; otherwise
#'   return an invalid fit with a warning.
#'
#' @return An \code{hd_classifier_fit} whose score is twice the log posterior
#'   density ratio of class 1 to class 2.
#'
#' @references
#' Anderson, T. W. (2003). *An Introduction to Multivariate Statistical
#' Analysis*, 3rd ed. Wiley.
#'
#' Feng, L. (2026). *High-Dimensional Data Analysis for Elliptical Symmetric
#' Distributions*, Chapter 5 (book manuscript). The supplied sparse-matrix
#' adapter is a package construction, not a separate primary method.
#' @export
#'
#' @examples
#' x <- rbind(c(-1, 0), c(-2, 1), c(1, 0), c(2, -1))
#' y <- factor(rep(c("left", "right"), each = 2))
#' fit <- sparse_plugin_qda(
#'   x, y, mean1 = c(-1.5, 0.5), mean2 = c(1.5, -0.5),
#'   covariance1 = diag(2), covariance2 = diag(c(1, 2))
#' )
sparse_plugin_qda <- function(
    x, y, mean1, mean2, covariance1, covariance2,
    prior = "equal", strict = TRUE) {
  call <- match.call()
  strict <- .clf_validate_strict(strict)
  training <- .clf_prepare_xy(
    x, y, prior = prior, equal_prior = FALSE, min_class = 2L
  )
  mean1 <- .c5qda_location(mean1, training$p, "mean1")
  mean2 <- .c5qda_location(mean2, training$p, "mean2")
  covariance_fit <- tryCatch(
    list(
      first = .clf_validate_spd(
        covariance1, "covariance1", p = training$p
      ),
      second = .clf_validate_spd(
        covariance2, "covariance2", p = training$p
      )
    ),
    error = function(error) error
  )
  if (inherits(covariance_fit, "error")) {
    return(.c5qda_failure(
      conditionMessage(covariance_fit), strict, "Certified sparse plug-in QDA",
      training, call, "twice log posterior density ratio",
      "covariance_certificate",
      diagnostics = list(explicit.covariance.input = TRUE)
    ))
  }
  precision1 <- .clf_solve_chol(
    covariance_fit$first$chol, diag(training$p)
  )
  precision2 <- .clf_solve_chol(
    covariance_fit$second$chol, diag(training$p)
  )
  quadratic <- precision2 - precision1
  linear <- 2 * drop(precision1 %*% mean1) -
    2 * drop(precision2 %*% mean2)
  intercept <- as.numeric(
    -crossprod(mean1, precision1 %*% mean1) +
      crossprod(mean2, precision2 %*% mean2) -
      covariance_fit$first$log.determinant +
      covariance_fit$second$log.determinant +
      2 * training$prior$log.ratio
  )
  .clf_new_fit(
    method = "Certified sparse plug-in QDA", training = training,
    score_model = list(
      type = "quadratic", quadratic = quadratic,
      linear = linear, intercept = intercept
    ),
    estimate = list(
      mean1 = mean1, mean2 = mean2,
      covariance1 = covariance_fit$first$matrix,
      covariance2 = covariance_fit$second$matrix,
      precision1 = precision1, precision2 = precision2
    ),
    diagnostics = list(
      covariance.certificate = list(
        reciprocal.condition = c(
          covariance_fit$first$reciprocal.condition,
          covariance_fit$second$reciprocal.condition
        ),
        cholesky = TRUE, explicitly.supplied = TRUE
      ),
      no.numerical.repair = TRUE
    ),
    call = call,
    primary_orientation = "positive posterior-density contrast selects class1",
    score_scale = "twice log posterior density ratio", tie = "class1",
    data_name = deparse(call$x)
  )
}
