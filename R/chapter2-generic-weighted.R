.ch2gw_matrix <- function(x, name = "x", min_rows = 2L) {
  if (is.data.frame(x)) {
    if (!all(vapply(x, is.numeric, logical(1)))) {
      stop(sprintf("`%s` must contain only numeric columns.", name),
           call. = FALSE)
    }
    x <- as.matrix(x)
  }
  if (!is.matrix(x) || !is.numeric(x) || length(dim(x)) != 2L) {
    stop(sprintf("`%s` must be a numeric matrix or numeric data frame.",
                 name), call. = FALSE)
  }
  if (nrow(x) < min_rows || ncol(x) < 1L) {
    stop(sprintf(
      "`%s` must have at least %d rows and one column.", name, min_rows
    ), call. = FALSE)
  }
  if (anyNA(x) || any(!is.finite(x))) {
    stop(sprintf("`%s` must contain only finite values.", name),
         call. = FALSE)
  }
  storage.mode(x) <- "double"
  x
}


.ch2gw_vector <- function(value, length.out, name, positive = FALSE) {
  if (!is.numeric(value) || !is.null(dim(value)) ||
      length(value) != length.out || anyNA(value) ||
      any(!is.finite(value)) || (positive && any(value <= 0))) {
    qualifier <- if (positive) "strictly positive finite numeric" else
      "finite numeric"
    stop(sprintf("`%s` must be a length-%d %s vector.",
                 name, length.out, qualifier), call. = FALSE)
  }
  as.numeric(value)
}


.ch2gw_scalar <- function(value, name, positive = FALSE,
                          nonnegative = FALSE, integer = FALSE) {
  if (!is.numeric(value) || length(value) != 1L || is.na(value) ||
      !is.finite(value) || (positive && value <= 0) ||
      (nonnegative && value < 0) ||
      (integer && (value != floor(value) ||
                   value > .Machine$integer.max))) {
    description <- if (integer) {
      "one finite positive integer"
    } else if (positive) {
      "one strictly positive finite number"
    } else if (nonnegative) {
      "one non-negative finite number"
    } else {
      "one finite number"
    }
    stop(sprintf("`%s` must be %s.", name, description), call. = FALSE)
  }
  if (integer) as.integer(value) else as.numeric(value)
}


.ch2gw_flag <- function(value, name) {
  if (!is.logical(value) || length(value) != 1L || is.na(value)) {
    stop(sprintf("`%s` must be TRUE or FALSE.", name), call. = FALSE)
  }
  value
}


.ch2gw_weight_spec <- function(K, power, power.supplied) {
  if (is.function(K)) {
    if (power.supplied) {
      stop("`power` is only used when `K = \"power\"`; it is not silently ignored for an R callback.",
           call. = FALSE)
    }
    return(list(kind = "callback", label = "user-supplied scalar R callback",
                callback = K, power = NULL))
  }
  if (!is.character(K) || length(K) != 1L || is.na(K)) {
    stop("`K` must be one of \"constant\", \"inverse_norm\", \"power\", or an R function.",
         call. = FALSE)
  }
  kind <- match.arg(K, c("constant", "inverse_norm", "power"))
  if (kind != "power" && power.supplied) {
    stop("`power` is only used when `K = \"power\"`; it is not silently ignored.",
         call. = FALSE)
  }
  if (kind == "power") {
    power <- .ch2gw_scalar(power, "power")
    return(list(kind = kind, label = sprintf("K(r) = r^(%.17g)", power),
                callback = NULL, power = power))
  }
  list(
    kind = kind,
    label = if (kind == "constant") "K(r) = 1" else "K(r) = 1/r",
    callback = NULL,
    power = if (kind == "constant") 0 else -1
  )
}


.ch2gw_evaluate_weights <- function(radii, specification) {
  radii <- as.numeric(radii)
  if (!length(radii) || anyNA(radii) || any(!is.finite(radii)) ||
      any(radii <= 0)) {
    stop("Internal radii must be strictly positive and finite.",
         call. = FALSE)
  }
  weights <- switch(
    specification$kind,
    constant = rep(1, length(radii)),
    inverse_norm = 1 / radii,
    power = {
      if (specification$power == 0) {
        rep(1, length(radii))
      } else if (specification$power == -1) {
        1 / radii
      } else {
        radii^specification$power
      }
    },
    callback = vapply(seq_along(radii), function(index) {
      value <- tryCatch(
        specification$callback(radii[index]),
        error = function(error) {
          stop(sprintf(
            "The radial-weight callback failed at radius %d: %s",
            index, conditionMessage(error)
          ), call. = FALSE)
        }
      )
      if (!is.numeric(value) || !is.null(dim(value)) ||
          length(value) != 1L || is.na(value) || !is.finite(value)) {
        stop(sprintf(
          "The radial-weight callback must return one finite numeric scalar at radius %d.",
          index
        ), call. = FALSE)
      }
      as.numeric(value)
    }, numeric(1))
  )
  if (anyNA(weights) || any(!is.finite(weights))) {
    stop("The radial weights are not all finite; no cap, truncation, or replacement is applied.",
         call. = FALSE)
  }
  as.numeric(weights)
}


.ch2gw_variable_names <- function(x) {
  if (is.null(colnames(x))) paste0("variable", seq_len(ncol(x))) else
    colnames(x)
}


.ch2gw_observation_names <- function(x) {
  if (is.null(rownames(x))) paste0("observation", seq_len(nrow(x))) else
    rownames(x)
}


#' Generic weighted Hettmansperger--Randles location estimator
#'
#' Solves the formula-complete weighted location equation in Chapter 2 while
#' retaining the unweighted diagonal Hettmansperger--Randles (HR) scale
#' equation. At iteration \eqn{m}, let
#' \deqn{e_i^{(m)}=(D^{(m)})^{-1/2}(X_i-\theta^{(m)}),\quad
#' r_i^{(m)}=\|e_i^{(m)}\|,\quad U_i^{(m)}=e_i^{(m)}/r_i^{(m)}.}
#' The updates are
#' \deqn{\theta^{(m+1)}=\theta^{(m)}+(D^{(m)})^{1/2}
#' \frac{\sum_iK(r_i^{(m)})U_i^{(m)}}
#' {\sum_iK(r_i^{(m)})/r_i^{(m)}}}
#' and
#' \deqn{D^{(m+1)}=pD^{(m)}\operatorname{diag}
#' \{n^{-1}\sum_iU_i^{(m)}U_i^{(m)T}\}.}
#' Thus weights enter the location equation only; the scale update is always
#' the unweighted HR update printed in the manuscript.
#'
#' `K` may be an R function or one of the auditable built-ins `"constant"`,
#' `"inverse_norm"`, and `"power"`. A user callback is called separately on
#' each scalar radius and must return exactly one finite numeric scalar. The
#' recursion rejects zero radii, a non-positive weighted denominator,
#' non-positive scale iterates, and non-convergence. It never floors a radius,
#' caps a weight, flips a denominator sign, adds a ridge, or returns an
#' unconverged last iterate.
#'
#' @param x Numeric matrix or data frame with observations in rows.
#' @param K Radial-weight name or scalar R callback.
#' @param power Finite exponent used only when `K = "power"`.
#' @param initial_location Optional finite initial location. The default is the
#'   vector of sample column means.
#' @param initial_diagonal Optional strictly positive initial diagonal scale.
#'   The default is the vector of unbiased marginal sample variances.
#' @param tol Strictly positive relative iterate tolerance.
#' @param max_iter Positive maximum number of simultaneous location/scale
#'   updates.
#' @param zero_tol Non-negative threshold for a singular standardized radius.
#' @param keep_history Whether to retain all location and diagonal iterates.
#' @return A list of class `generic_weighted_hr_location` containing the fitted
#'   location and diagonal, final directions/radii/weights, equation residuals,
#'   and strict convergence diagnostics.
#' @references Feng, L., Liu, B. and Ma, Y. (2021). An inverse norm sign test
#' for location parameters in high-dimensional data. Journal of Business &
#' Economic Statistics 39, 807--815.
#' @examples
#' x <- matrix(c(-2, 1, 0, 3, -1, 2, 1, -3, 2, 0, 4, -2), ncol = 2)
#' generic_weighted_hr_location(x, K = "constant", tol = 1e-6)
#' @export
generic_weighted_hr_location <- function(
    x, K = "constant", power = 0, initial_location = NULL,
    initial_diagonal = NULL, tol = 1e-8, max_iter = 500L,
    zero_tol = 0, keep_history = FALSE) {
  call <- match.call()
  data.name <- deparse1(substitute(x))
  power.supplied <- !missing(power)
  x <- .ch2gw_matrix(x, "x", min_rows = 2L)
  specification <- .ch2gw_weight_spec(K, power, power.supplied)
  tol <- .ch2gw_scalar(tol, "tol", positive = TRUE)
  max_iter <- .ch2gw_scalar(
    max_iter, "max_iter", positive = TRUE, integer = TRUE
  )
  zero_tol <- .ch2gw_scalar(
    zero_tol, "zero_tol", nonnegative = TRUE
  )
  keep_history <- .ch2gw_flag(keep_history, "keep_history")
  initial <- cpp_ch2_generic_weighted_initial(x)
  if (!is.null(initial_location)) {
    initial$location <- .ch2gw_vector(
      initial_location, ncol(x), "initial_location"
    )
  }
  if (!is.null(initial_diagonal)) {
    initial$diagonal <- .ch2gw_vector(
      initial_diagonal, ncol(x), "initial_diagonal", positive = TRUE
    )
  }
  location <- as.numeric(initial$location)
  diagonal <- as.numeric(initial$diagonal)
  initial.location <- location
  initial.diagonal <- diagonal
  location.history <- list(location)
  diagonal.history <- list(diagonal)
  converged <- FALSE
  relative.update <- Inf
  location.update <- Inf
  diagonal.update <- Inf

  for (iteration in seq_len(max_iter)) {
    geometry <- cpp_ch2_generic_weighted_geometry(
      x, location, diagonal, zero_tol
    )
    weights <- .ch2gw_evaluate_weights(
      geometry$radii, specification
    )
    step <- cpp_ch2_generic_weighted_step(
      location, diagonal, geometry$directions,
      geometry$radii, weights
    )
    next.location <- as.numeric(step$next_location)
    next.diagonal <- as.numeric(geometry$next_diagonal)
    location.update <- max(
      abs(next.location - location) / sqrt(initial.diagonal)
    )
    diagonal.update <- max(abs(log(next.diagonal / diagonal)))
    relative.update <- max(location.update, diagonal.update)
    if (!is.finite(relative.update)) {
      stop("The weighted HR relative update is not finite; no repair is applied.",
           call. = FALSE)
    }
    location <- next.location
    diagonal <- next.diagonal
    if (keep_history) {
      location.history[[length(location.history) + 1L]] <- location
      diagonal.history[[length(diagonal.history) + 1L]] <- diagonal
    }
    if (relative.update <= tol) {
      converged <- TRUE
      break
    }
  }
  if (!converged) {
    stop(sprintf(
      paste0(
        "The generic weighted HR recursion did not converge within %d ",
        "updates (last relative update %.6g); no unconverged iterate is ",
        "returned."
      ),
      max_iter, relative.update
    ), call. = FALSE)
  }

  final.geometry <- cpp_ch2_generic_weighted_geometry(
    x, location, diagonal, zero_tol
  )
  final.weights <- .ch2gw_evaluate_weights(
    final.geometry$radii, specification
  )
  final.step <- cpp_ch2_generic_weighted_step(
    location, diagonal, final.geometry$directions,
    final.geometry$radii, final.weights
  )
  diagonal.residual <- max(abs(final.geometry$diagonal_equation - 1))
  variable.names <- .ch2gw_variable_names(x)
  observation.names <- .ch2gw_observation_names(x)
  names(location) <- variable.names
  names(diagonal) <- variable.names
  names(initial.location) <- variable.names
  names(initial.diagonal) <- variable.names
  names(final.step$numerator) <- variable.names
  names(final.geometry$diagonal_equation) <- variable.names
  colnames(final.geometry$directions) <- variable.names
  rownames(final.geometry$directions) <- observation.names
  names(final.geometry$radii) <- observation.names
  names(final.weights) <- observation.names

  result <- list(
    location = location,
    scale.diagonal = diagonal,
    log.scale.diagonal = log(diagonal),
    directions = final.geometry$directions,
    radii = final.geometry$radii,
    weights = final.weights,
    weight = specification$label,
    power = specification$power,
    initial.location = initial.location,
    initial.scale.diagonal = initial.diagonal,
    weighted.location.numerator = final.step$numerator,
    weighted.location.denominator = as.numeric(final.step$denominator),
    normalized.location.equation.residual = as.numeric(
      final.step$normalized_location_residual
    ),
    diagonal.equation = final.geometry$diagonal_equation,
    diagonal.equation.residual = diagonal.residual,
    diagnostics = list(
      converged = TRUE,
      iterations = iteration,
      relative.update = relative.update,
      location.relative.update = location.update,
      log.diagonal.relative.update = diagonal.update,
      convergence.basis = paste(
        "maximum of initial-scale-standardized location update and",
        "log-diagonal update"
      ),
      location.update = "weighted Chapter 2 equation",
      diagonal.update = "unweighted HR equation",
      diagonal.common.scale.identification = paste(
        "inherited from the initial marginal variance scale;",
        "no trace renormalization"
      ),
      callback.evaluation = if (specification$kind == "callback") {
        "one scalar radius per R call"
      } else {
        "built-in vector evaluation"
      },
      minimum.radius = as.numeric(final.geometry$minimum_radius),
      zero.radius.policy = "error",
      numerical.repair = "none"
    ),
    n = nrow(x),
    p = ncol(x),
    data.name = data.name,
    call = call
  )
  if (keep_history) {
    result$location.history <- do.call(rbind, location.history)
    result$scale.diagonal.history <- do.call(rbind, diagonal.history)
    colnames(result$location.history) <- variable.names
    colnames(result$scale.diagonal.history) <- variable.names
    rownames(result$location.history) <- paste0(
      "iterate", seq_len(nrow(result$location.history)) - 1L
    )
    rownames(result$scale.diagonal.history) <- rownames(
      result$location.history
    )
  }
  class(result) <- c("generic_weighted_hr_location", "list")
  result
}


#' Oracle generic weighted-sign sum statistic
#'
#' Evaluates the Chapter 2 oracle statistic from a supplied reference location
#' and diagonal shape. With
#' \deqn{r_i=\|D^{-1/2}(X_i-\theta)\|,\quad
#' V_i(K)=K(r_i)U\{D^{-1/2}(X_i-\theta)\},}
#' the raw score is
#' \deqn{T_n(K)=\frac{2}{n(n-1)}\sum_{i<j}V_i(K)^TV_j(K).}
#'
#' This function deliberately says `oracle`. It does not estimate the
#' reference location, diagonal, radial moment, or correlation trace. A valid
#' upper-tail normal p-value is returned only when the caller supplies either
#' `null_sd`, or both `nu2` and `trace_R2`, in which case
#' \deqn{\sigma_{n,K}=\left\{
#' \frac{2\nu_{2,K}^2\operatorname{tr}(R^2)}
#' {n(n-1)p^2}\right\}^{1/2}.}
#' With no calibration inputs the raw score remains available and
#' `calibrated` is `FALSE`; no feasible leave-out or plug-in estimator is
#' guessed.
#'
#' @param x Numeric matrix or data frame with observations in rows.
#' @param theta Required finite reference location vector.
#' @param diagonal Required strictly positive reference diagonal shape.
#' @param K Radial-weight name or scalar R callback, as in
#'   `generic_weighted_hr_location()`.
#' @param power Finite exponent used only when `K = "power"`.
#' @param null_sd Optional strictly positive supplied null standard deviation.
#' @param nu2 Optional strictly positive supplied population moment
#'   \eqn{\nu_{2,K}=E\{K^2(r_i)\}}.
#' @param trace_R2 Optional strictly positive supplied
#'   \eqn{\operatorname{tr}(R^2)}. It must be supplied together with `nu2`.
#' @param zero_tol Non-negative threshold for a singular standardized radius.
#' @param keep_scores Whether to retain directions, radii, weights, and
#'   weighted score means.
#' @return A list of class `oracle_weighted_sign_sum_test`. The `p.value` field
#'   is absent unless explicit calibration is complete.
#' @references Feng, L., Liu, B. and Ma, Y. (2021). An inverse norm sign test
#' for location parameters in high-dimensional data. Journal of Business &
#' Economic Statistics 39, 807--815.
#' @examples
#' x <- matrix(c(-2, 1, 0, 3, -1, 2, 1, -3, 2, 0, 4, -2), ncol = 2)
#' oracle_weighted_sign_sum_test(
#'   x, theta = c(0, 0), diagonal = c(1, 1), K = "constant"
#' )
#' @export
oracle_weighted_sign_sum_test <- function(
    x, theta, diagonal, K = "constant", power = 0,
    null_sd = NULL, nu2 = NULL, trace_R2 = NULL,
    zero_tol = 0, keep_scores = FALSE) {
  call <- match.call()
  data.name <- deparse1(substitute(x))
  power.supplied <- !missing(power)
  x <- .ch2gw_matrix(x, "x", min_rows = 2L)
  theta <- .ch2gw_vector(theta, ncol(x), "theta")
  diagonal <- .ch2gw_vector(
    diagonal, ncol(x), "diagonal", positive = TRUE
  )
  specification <- .ch2gw_weight_spec(K, power, power.supplied)
  zero_tol <- .ch2gw_scalar(
    zero_tol, "zero_tol", nonnegative = TRUE
  )
  keep_scores <- .ch2gw_flag(keep_scores, "keep_scores")

  has.null.sd <- !is.null(null_sd)
  has.nu2 <- !is.null(nu2)
  has.trace <- !is.null(trace_R2)
  if (has.null.sd && (has.nu2 || has.trace)) {
    stop("Supply either `null_sd` or the pair `nu2` and `trace_R2`, not conflicting calibration routes.",
         call. = FALSE)
  }
  if (xor(has.nu2, has.trace)) {
    stop("`nu2` and `trace_R2` must be supplied together; no missing calibration component is estimated.",
         call. = FALSE)
  }
  if (has.null.sd) {
    null_sd <- .ch2gw_scalar(null_sd, "null_sd", positive = TRUE)
  }
  if (has.nu2) {
    nu2 <- .ch2gw_scalar(nu2, "nu2", positive = TRUE)
    trace_R2 <- .ch2gw_scalar(
      trace_R2, "trace_R2", positive = TRUE
    )
  }

  geometry <- cpp_ch2_generic_weighted_geometry(
    x, theta, diagonal, zero_tol
  )
  weights <- .ch2gw_evaluate_weights(
    geometry$radii, specification
  )
  core <- cpp_ch2_generic_weighted_quadratic(
    geometry$directions, weights
  )
  calibrated <- has.null.sd || has.nu2
  calibration.route <- "none"
  if (has.nu2) {
    log.null.sd <- 0.5 * log(2) + log(nu2) +
      0.5 * log(trace_R2) -
      0.5 * (log(nrow(x)) + log(nrow(x) - 1) + 2 * log(ncol(x)))
    if (!is.finite(log.null.sd) ||
        log.null.sd > log(.Machine$double.xmax) ||
        log.null.sd < log(.Machine$double.xmin)) {
      stop("The supplied-moment null standard deviation is not representable; no floor or cap is applied.",
           call. = FALSE)
    }
    null_sd <- exp(log.null.sd)
    calibration.route <- "supplied nu2 and trace_R2"
  } else if (has.null.sd) {
    calibration.route <- "supplied null_sd"
  }
  score <- as.numeric(core$score)
  statistic <- stats::setNames(score, "T.n.K")
  if (calibrated) {
    standardized <- score / null_sd
    if (!is.finite(standardized)) {
      stop("The standardized oracle weighted-sign statistic is not finite; no repair is applied.",
           call. = FALSE)
    }
    statistic <- stats::setNames(standardized, "oracle weighted-sign Z")
  }
  variable.names <- .ch2gw_variable_names(x)
  observation.names <- .ch2gw_observation_names(x)
  names(theta) <- variable.names
  names(diagonal) <- variable.names

  result <- list(
    statistic = statistic,
    score = stats::setNames(score, "T.n.K"),
    raw.statistic = stats::setNames(score, "T.n.K"),
    calibrated = calibrated,
    null.sd = if (calibrated) null_sd else NULL,
    alternative = "greater",
    reference.location = theta,
    reference.diagonal = diagonal,
    weight = specification$label,
    power = specification$power,
    empirical.nu2 = as.numeric(core$empirical_nu2),
    supplied.nu2 = if (has.nu2) nu2 else NULL,
    supplied.trace.R2 = if (has.trace) trace_R2 else NULL,
    diagnostics = list(
      scope = "oracle supplied location and diagonal only",
      calibration = calibration.route,
      p.value.available = calibrated,
      feasible.leaveout.or.plugin = FALSE,
      callback.evaluation = if (specification$kind == "callback") {
        "one scalar radius per R call"
      } else {
        "built-in vector evaluation"
      },
      minimum.radius = as.numeric(geometry$minimum_radius),
      zero.radius.policy = "error",
      numerical.repair = "none"
    ),
    n = nrow(x),
    p = ncol(x),
    data.name = data.name,
    call = call
  )
  if (calibrated) {
    result$p.value <- stats::pnorm(
      unname(statistic), lower.tail = FALSE
    )
  }
  if (keep_scores) {
    colnames(geometry$directions) <- variable.names
    rownames(geometry$directions) <- observation.names
    names(geometry$radii) <- observation.names
    names(weights) <- observation.names
    names(core$weighted_score_mean) <- variable.names
    result$directions <- geometry$directions
    result$radii <- geometry$radii
    result$weights <- weights
    result$weighted.score.mean <- core$weighted_score_mean
  }
  result <- result[!vapply(result, is.null, logical(1))]
  class(result) <- c("oracle_weighted_sign_sum_test", "list")
  result
}




