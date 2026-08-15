.ch3ef_strict <- function(strict) {
  if (!is.logical(strict) || length(strict) != 1L || is.na(strict)) {
    stop("`strict` must be TRUE or FALSE.", call. = FALSE)
  }
  strict
}


.ch3ef_warn_or_stop <- function(message, strict) {
  if (isTRUE(strict)) {
    stop(message, call. = FALSE)
  }
  warning(message, call. = FALSE)
  invisible(FALSE)
}


.ch3ef_integer <- function(value, name, minimum = 0L, maximum = Inf) {
  if (!is.numeric(value) || length(value) != 1L || is.na(value) ||
      !is.finite(value) || value != floor(value) || value < minimum ||
      value > maximum || value > .Machine$integer.max) {
    stop(sprintf(
      "`%s` must be one integer between %s and %s.",
      name, format(minimum), format(maximum)
    ), call. = FALSE)
  }
  as.integer(value)
}


.ch3ef_sign_fit <- function(x, tol, max_iter, zero_tol, strict) {
  controls <- .validate_iteration_controls(tol, max_iter, zero_tol)
  fit <- .ch3pp_spatial_fit(
    x, controls$tol, controls$max_iter, controls$zero_tol
  )
  if (!isTRUE(fit$diagnostics$converged)) {
    .ch3ef_warn_or_stop(
      paste(
        "The spatial median did not satisfy its equation-residual",
        "convergence criterion; no factor estimator is certified."
      ), strict
    )
  }
  fit$valid <- isTRUE(fit$diagnostics$converged)
  fit$controls <- controls
  fit
}


.ch3ef_eigen <- function(matrix, factors) {
  matrix <- (matrix + t(matrix)) / 2
  decomposition <- eigen(matrix, symmetric = TRUE)
  values <- as.numeric(decomposition$values)
  scale <- max(1, max(abs(values)))
  material.negative <- values < -512 * .Machine$double.eps * scale
  if (any(material.negative)) {
    stop(
      "The factor pilot has a materially negative eigenvalue; no positive-" %+%
        "definiteness repair is applied.",
      call. = FALSE
    )
  }
  rounded <- values < 0
  values[rounded] <- 0
  vectors <- decomposition$vectors
  for (j in seq_len(ncol(vectors))) {
    first <- which(abs(vectors[, j]) > 0)[1L]
    if (!is.na(first) && vectors[first, j] < 0) {
      vectors[, j] <- -vectors[, j]
    }
  }
  p <- nrow(matrix)
  if (factors > 0L) {
    leading.values <- values[seq_len(factors)]
    leading.vectors <- vectors[, seq_len(factors), drop = FALSE]
    low.rank <- tcrossprod(
      sweep(leading.vectors, 2L, sqrt(leading.values), "*")
    )
  } else {
    leading.values <- numeric()
    leading.vectors <- matrix(numeric(), p, 0L)
    low.rank <- matrix(0, p, p)
  }
  list(
    values = values,
    vectors = vectors,
    leading.values = leading.values,
    leading.vectors = leading.vectors,
    low.rank = low.rank,
    idiosyncratic = (matrix - low.rank + t(matrix - low.rank)) / 2,
    negative.roundoff.clipped = sum(rounded)
  )
}


.ch3ef_factor_limit <- function(n, p) max(0L, min(n, p) - 1L)


.ch3ef_threshold <- function(threshold, constant, n, p) {
  if (is.null(threshold)) {
    constant <- as.numeric(constant)
    if (length(constant) != 1L || is.na(constant) ||
        !is.finite(constant) || constant < 0) {
      stop("`constant` must be one finite non-negative number.",
           call. = FALSE)
    }
    rate <- sqrt(log(p) / n) + sqrt(log(n) / n)
    list(
      value = constant * rate,
      source = "constant * {sqrt(log(p)/n) + sqrt(log(n)/n)}",
      constant = constant,
      rate = rate
    )
  } else {
    threshold <- as.numeric(threshold)
    if (length(threshold) != 1L || is.na(threshold) ||
        !is.finite(threshold) || threshold < 0) {
      stop("`threshold` must be one finite non-negative number.",
           call. = FALSE)
    }
    list(
      value = threshold, source = "supplied", constant = NULL,
      rate = sqrt(log(p) / n) + sqrt(log(n) / n)
    )
  }
}


.ch3ef_poet_core <- function(
    pilot, center, signs, n, factors, threshold.fit, rule.fit,
    pilot.type, median.diagnostics, zero.radii, data.name, call,
    extra.diagnostics = list()) {
  p <- ncol(pilot)
  maximum <- .ch3ef_factor_limit(n, p)
  factors <- .ch3ef_integer(factors, "factors", 0L, maximum)
  eig <- .ch3ef_eigen(pilot, factors)
  threshold.matrix <- matrix(threshold.fit$value, p, p)
  thresholded <- cpp_ch3_gaussian_threshold_matrix(
    eig$idiosyncratic, threshold.matrix, rule.fit$rule,
    rule.fit$shape, FALSE
  )
  estimate <- (eig$low.rank + thresholded)
  estimate <- (estimate + t(estimate)) / 2
  variable.names <- colnames(pilot)
  if (!is.null(variable.names)) {
    matrices <- list(pilot, estimate, eig$low.rank, eig$idiosyncratic,
                     thresholded, threshold.matrix)
    matrices <- lapply(matrices, function(value) {
      dimnames(value) <- list(variable.names, variable.names)
      value
    })
    pilot <- matrices[[1L]]
    estimate <- matrices[[2L]]
    eig$low.rank <- matrices[[3L]]
    eig$idiosyncratic <- matrices[[4L]]
    thresholded <- matrices[[5L]]
    threshold.matrix <- matrices[[6L]]
    rownames(eig$leading.vectors) <- variable.names
    names(center) <- variable.names
  }
  minimum.eigenvalue <- min(eigen(
    (estimate + t(estimate)) / 2, symmetric = TRUE, only.values = TRUE
  )$values)
  center.output <- as.numeric(center)
  names(center.output) <- variable.names
  structure(
    list(
      estimate = estimate,
      valid = TRUE,
      center = center.output,
      signs = signs,
      pilot = pilot,
      factors = factors,
      threshold = threshold.matrix,
      rule = rule.fit$rule,
      data.name = data.name,
      components = list(
        eigenvalues = eig$values,
        leading.eigenvalues = eig$leading.values,
        leading.eigenvectors = eig$leading.vectors,
        low.rank = eig$low.rank,
        idiosyncratic.raw = eig$idiosyncratic,
        idiosyncratic.thresholded = thresholded
      ),
      diagnostics = c(list(
        n = n, p = p, pilot.type = pilot.type,
        factor.count.source = "supplied",
        threshold.source = threshold.fit$source,
        threshold.constant = threshold.fit$constant,
        theoretical.rate = threshold.fit$rate,
        threshold.diagonal = FALSE,
        minimum.estimate.eigenvalue = minimum.eigenvalue,
        positive.definiteness.repair = "none",
        negative.eigenvalue.roundoff.count =
          eig$negative.roundoff.clipped,
        spatial.median = median.diagnostics,
        zero.fitted.radii = zero.radii,
        tuning.selection = "none",
        call = call
      ), extra.diagnostics)
    ),
    class = c("elliptical_factor_fit", "hd_covariance_estimator", "list")
  )
}


.ch3ef_spd_inverse <- function(matrix, name) {
  matrix <- as.matrix(matrix)
  if (!is.numeric(matrix) || nrow(matrix) != ncol(matrix) ||
      anyNA(matrix) || any(!is.finite(matrix))) {
    stop(sprintf("`%s` must be a finite square matrix.", name),
         call. = FALSE)
  }
  scale <- max(1, max(abs(matrix)))
  error <- max(abs(matrix - t(matrix)))
  if (error > 512 * .Machine$double.eps * scale) {
    stop(sprintf("`%s` must be numerically symmetric.", name),
         call. = FALSE)
  }
  matrix <- (matrix + t(matrix)) / 2
  factor <- tryCatch(chol(matrix), error = identity)
  if (inherits(factor, "condition")) {
    stop(sprintf(
      "`%s` must be positive definite; no ridge or eigenvalue floor is applied.",
      name
    ), call. = FALSE)
  }
  list(matrix = matrix, inverse = chol2inv(factor), rcond = rcond(matrix))
}


#' Select the number of spikes in an elliptical factor model
#'
#' Computes the trace-\eqn{p} spatial-sign pilot centered at the sample spatial
#' median and applies the eigenvalue-ratio (ER) or growth-ratio (GR) selector
#' of Xu, Ma, Wang and Feng.  For ordered eigenvalues \eqn{\lambda_j},
#' \deqn{\widehat m_{ER}=\arg\max_{1\le j\le M}
#' \lambda_j/\lambda_{j+1}.}
#' The GR criterion uses
#' \eqn{V_j=\sum_{l=j+1}^{\min(n,p)-1}\lambda_l} exactly as in the
#' primary paper. `max_factors` is deliberately supplied: the paper treats
#' \eqn{M} as a predetermined upper bound and its numerical choice is not a
#' universal tuning rule.
#'
#' @param x Numeric matrix or data frame with observations in rows.
#' @param max_factors Predetermined positive integer upper bound \eqn{M}.
#' @param method Either eigenvalue ratio (`"er"`) or growth ratio (`"gr"`).
#' @param tol,max_iter,zero_tol Spatial-median controls.
#' @param strict If `TRUE`, fail when the spatial median is not certified;
#'   otherwise return an explicitly invalid diagnostic object.
#'
#' @return An `elliptical_factor_number` object containing the selected count,
#'   all pilot eigenvalues, criterion values, and spatial-median diagnostics.
#' @references Xu, X., Ma, H., Wang, H. and Feng, L. (2025). High dimensional
#'   matrix estimation through elliptical factor models. arXiv:2512.19325.
#'   \url{https://arxiv.org/abs/2512.19325}.
#' @examples
#' x <- rbind(c(-3, -2, 0), c(-2, -1, 1), c(-1, 0, -1),
#'            c(1, 0, 1), c(2, 1, -1), c(3, 2, 0))
#' elliptical_factor_number(x, max_factors = 1)
#' @export
elliptical_factor_number <- function(
    x, max_factors, method = c("er", "gr"),
    tol = 1e-8, max_iter = 500L, zero_tol = 0, strict = TRUE) {
  call <- match.call()
  data.name <- deparse1(substitute(x))
  x <- .as_data_matrix(x, "x", min_rows = 2L)
  strict <- .ch3ef_strict(strict)
  method <- match.arg(method)
  n <- nrow(x)
  p <- ncol(x)
  rank.bound <- min(n, p)
  maximum <- if (method == "er") rank.bound - 1L else rank.bound - 2L
  if (maximum < 1L) {
    stop(sprintf(
      "The `%s` selector needs min(n, p) of at least %d.",
      method, if (method == "er") 2L else 3L
    ), call. = FALSE)
  }
  max_factors <- .ch3ef_integer(
    max_factors, "max_factors", 1L, maximum
  )
  fit <- .ch3ef_sign_fit(x, tol, max_iter, zero_tol, strict)
  if (!fit$valid) {
    return(structure(list(
      selected = NA_integer_, valid = FALSE, method = method,
      max.factors = max_factors, data.name = data.name,
      diagnostics = list(spatial.median = fit$diagnostics, call = call)
    ), class = "elliptical_factor_number"))
  }
  pilot <- p * crossprod(fit$signs) / n
  eig <- .ch3ef_eigen(pilot, 0L)
  values <- eig$values
  if (method == "er") {
    denominators <- values[seq.int(2L, max_factors + 1L)]
    criterion <- values[seq_len(max_factors)] / denominators
    criterion[denominators <= 0] <- Inf
    tails <- NULL
  } else {
    r <- rank.bound
    tails <- vapply(0:max_factors, function(j) {
      lower <- j + 1L
      upper <- r - 1L
      if (lower > upper) 0 else sum(values[lower:upper])
    }, numeric(1))
    criterion <- vapply(seq_len(max_factors), function(j) {
      numerator.tail <- tails[j]
      denominator.tail <- tails[j + 1L]
      if (!(numerator.tail > 0) || !(denominator.tail > 0)) return(NA_real_)
      numerator <- log1p(values[j] / numerator.tail)
      denominator <- log1p(values[j + 1L] / denominator.tail)
      if (!(denominator > 0)) NA_real_ else numerator / denominator
    }, numeric(1))
    if (all(!is.finite(criterion))) {
      stop(
        "The growth-ratio criterion is undefined because its tail sums " %+%
          "are non-positive.", call. = FALSE
      )
    }
  }
  criterion.for.selection <- criterion
  criterion.for.selection[is.na(criterion.for.selection) |
                            is.nan(criterion.for.selection)] <- -Inf
  selected <- which.max(criterion.for.selection)
  structure(list(
    selected = as.integer(selected), valid = TRUE, method = method,
    max.factors = max_factors, criterion = criterion,
    eigenvalues = values, tail.sums = tails, pilot = pilot,
    center = fit$center, signs = fit$signs, data.name = data.name,
    diagnostics = list(
      n = n, p = p, rank.bound = rank.bound,
      factor.upper.bound.source = "supplied",
      spatial.median = fit$diagnostics,
      zero.fitted.radii = fit$n.zero,
      tie.break = "smallest index returned by which.max",
      call = call
    )
  ), class = "elliptical_factor_number")
}


#' Spatial-sign POET estimator for an elliptical factor model
#'
#' Implements POET-SS from Xu et al. The pilot is
#' \deqn{\widehat\Sigma_0=\frac{p}{n}\sum_i
#' U(X_i-\widehat\mu)U(X_i-\widehat\mu)^T,}
#' its supplied leading `factors` eigenpairs form the low-rank component, and
#' a generalized thresholding rule is applied only to off-diagonal entries of
#' the raw idiosyncratic complement.  When `threshold = NULL`, the common
#' threshold is `constant * (sqrt(log(p)/n) + sqrt(log(n)/n))`.  The unknown
#' sufficiently-large theoretical constant is not estimated from the same
#' data or borrowed from a paper simulation.
#'
#' The result may be indefinite in finite samples because thresholding is
#' used literally. No eigenvalue floor, ridge, nearest-PD projection, or
#' pseudoinverse is applied.
#'
#' @inheritParams elliptical_factor_number
#' @param factors Supplied non-negative integer number of factors.
#' @param threshold Optional common threshold in trace-normalized scatter
#'   units. `NULL` uses the primary rate with `constant`.
#' @param constant Non-negative multiplier for the primary rate threshold.
#' @param rule One of `"hard"`, `"soft"`, `"scad"`, or
#'   `"adaptive_lasso"`.
#' @param scad_a,adaptive_eta Generalized-threshold rule parameters.
#'
#' @return An `elliptical_factor_fit` with the pilot, factor eigensystem, raw
#'   and thresholded idiosyncratic components, and reconstructed scatter.
#' @references Xu, X., Ma, H., Wang, H. and Feng, L. (2025). arXiv:2512.19325.
#' @examples
#' x <- rbind(c(-3, -2, 0), c(-2, -1, 1), c(-1, 0, -1),
#'            c(1, 0, 1), c(2, 1, -1), c(3, 2, 0))
#' spatial_sign_poet(x, factors = 1, threshold = 0.1)
#' @export
spatial_sign_poet <- function(
    x, factors, threshold = NULL, constant = 1,
    rule = c("hard", "soft", "scad", "adaptive_lasso"),
    scad_a = 3.7, adaptive_eta = 1,
    tol = 1e-8, max_iter = 500L, zero_tol = 0, strict = TRUE) {
  call <- match.call()
  data.name <- deparse1(substitute(x))
  x <- .as_data_matrix(x, "x", min_rows = 2L)
  strict <- .ch3ef_strict(strict)
  fit <- .ch3ef_sign_fit(x, tol, max_iter, zero_tol, strict)
  if (!fit$valid) {
    return(structure(list(
      estimate = NULL, valid = FALSE, data.name = data.name,
      diagnostics = list(spatial.median = fit$diagnostics, call = call)
    ), class = c("elliptical_factor_fit", "list")))
  }
  n <- nrow(x)
  p <- ncol(x)
  threshold.fit <- .ch3ef_threshold(threshold, constant, n, p)
  rule.fit <- .ch3g_threshold_rule(rule, scad_a, adaptive_eta)
  pilot <- p * crossprod(fit$signs) / n
  dimnames(pilot) <- list(colnames(x), colnames(x))
  .ch3ef_poet_core(
    pilot, fit$center, fit$signs, n, factors, threshold.fit, rule.fit,
    "spatial-sign covariance", fit$diagnostics, fit$n.zero,
    data.name, call,
    extra.diagnostics = list(
      radial.scale = fit$radial.scale,
      subtraction.overflow.fallback = fit$overflow.fallback
    )
  )
}


#' One-step Tyler POET estimator for an elliptical factor model
#'
#' Implements the primary POET-TME construction.  With a positive-definite
#' preliminary inverse shape \eqn{\widehat V_S}, it first computes
#' \deqn{\widehat\Sigma_T=\frac{p}{n}\sum_i
#' \frac{(X_i-\widehat\mu)(X_i-\widehat\mu)^T}
#' {(X_i-\widehat\mu)^T\widehat V_S(X_i-\widehat\mu)}}
#' and then applies the same supplied-factor POET decomposition once.  This is
#' the one-step estimator recommended in the primary paper, not a new
#' convergence iteration.
#'
#' If `preliminary_precision = NULL`, a POET-SS fit with the same factor and
#' threshold controls is constructed and inverted only when it is strictly
#' positive definite. A supplied matrix, or the certified `estimate` from an
#' `elliptical_factor_precision_fit`, may be used instead. No generalized
#' inverse, ridge, or positive-definiteness repair is performed.
#'
#' @inheritParams spatial_sign_poet
#' @param preliminary_precision Optional positive-definite preliminary inverse
#'   shape matrix or a valid `elliptical_factor_precision_fit`.
#'
#' @return An `elliptical_factor_fit` whose pilot is the one-step Tyler matrix.
#' @references Xu, X., Ma, H., Wang, H. and Feng, L. (2025). arXiv:2512.19325.
#' @examples
#' x <- rbind(c(-3, -2), c(-2, -1), c(-1, 1),
#'            c(1, -1), c(2, 1), c(3, 2))
#' poet_tme(x, factors = 0, threshold = 0.2,
#'          preliminary_precision = diag(2))
#' @export
poet_tme <- function(
    x, factors, threshold = NULL, constant = 1,
    rule = c("hard", "soft", "scad", "adaptive_lasso"),
    scad_a = 3.7, adaptive_eta = 1,
    preliminary_precision = NULL,
    tol = 1e-8, max_iter = 500L, zero_tol = 0, strict = TRUE) {
  call <- match.call()
  data.name <- deparse1(substitute(x))
  x <- .as_data_matrix(x, "x", min_rows = 2L)
  strict <- .ch3ef_strict(strict)
  fit <- .ch3ef_sign_fit(x, tol, max_iter, zero_tol, strict)
  if (!fit$valid) {
    return(structure(list(
      estimate = NULL, valid = FALSE, data.name = data.name,
      diagnostics = list(spatial.median = fit$diagnostics, call = call)
    ), class = c("elliptical_factor_fit", "list")))
  }
  n <- nrow(x)
  p <- ncol(x)
  threshold.fit <- .ch3ef_threshold(threshold, constant, n, p)
  rule.fit <- .ch3g_threshold_rule(rule, scad_a, adaptive_eta)
  preliminary.fit <- NULL
  source <- "supplied matrix"
  if (is.null(preliminary_precision)) {
    sign.pilot <- p * crossprod(fit$signs) / n
    dimnames(sign.pilot) <- list(colnames(x), colnames(x))
    preliminary.fit <- .ch3ef_poet_core(
      sign.pilot, fit$center, fit$signs, n, factors, threshold.fit, rule.fit,
      "spatial-sign covariance", fit$diagnostics, fit$n.zero,
      data.name, call
    )
    inverse <- .ch3ef_spd_inverse(
      preliminary.fit$estimate, "the preliminary POET-SS estimate"
    )
    precision <- inverse$inverse
    source <- "inverse of internally fitted POET-SS scatter"
    preliminary.rcond <- inverse$rcond
  } else {
    if (inherits(preliminary_precision, "elliptical_factor_precision_fit")) {
      if (!isTRUE(preliminary_precision$valid) ||
          is.null(preliminary_precision$estimate)) {
        stop("`preliminary_precision` is not a valid certified fit.",
             call. = FALSE)
      }
      preliminary_precision <- preliminary_precision$estimate
      source <- "certified elliptical_factor_precision_fit"
    }
    if (!is.matrix(preliminary_precision) ||
        !identical(dim(preliminary_precision), c(p, p))) {
      stop("`preliminary_precision` must be a p by p matrix.",
           call. = FALSE)
    }
    inverse <- .ch3ef_spd_inverse(
      preliminary_precision, "preliminary_precision"
    )
    precision <- inverse$matrix
    preliminary.rcond <- inverse$rcond
  }
  tyler <- cpp_ch3ef_tyler_one_step(
    x, as.numeric(fit$center), precision, fit$controls$zero_tol
  )
  tyler.pilot <- as.matrix(tyler$estimate)
  dimnames(tyler.pilot) <- list(colnames(x), colnames(x))
  .ch3ef_poet_core(
    tyler.pilot, fit$center, fit$signs, n, factors, threshold.fit, rule.fit,
    "one-step self-normalized Tyler", fit$diagnostics, fit$n.zero,
    data.name, call,
    extra.diagnostics = list(
      preliminary.precision.source = source,
      preliminary.precision.rcond = preliminary.rcond,
      preliminary.spatial.sign.fit = preliminary.fit,
      tyler = list(
        minimum.standardized.quadratic =
          tyler$minimum_standardized_quadratic,
        maximum.standardized.quadratic =
          tyler$maximum_standardized_quadratic,
        minimum.residual.radius = tyler$minimum_residual_radius,
        trace = tyler$trace,
        subtraction.overflow.fallback =
          tyler$subtraction_overflow_fallback
      ),
      refinement.steps = 1L,
      no.repair = paste(
        "No ridge, eigenvalue floor, pseudoinverse, or nearest-PD",
        "projection is used"
      )
    )
  )
}


#' Sparse precision estimation for an elliptical factor fit
#'
#' Applies the primary CLIME or GLASSO problem to the *raw* idiosyncratic
#' complement of a `spatial_sign_poet()` or `poet_tme()` fit, and reconstructs
#' the full inverse scatter through the Woodbury identity,
#' \deqn{V=V_u-V_u\Gamma(\Lambda^{-1}+\Gamma^T V_u\Gamma)^{-1}
#' \Gamma^T V_u.}
#' CLIME must pass per-column primal, dual, stationarity, gap, and final
#' symmetrized-feasibility certificates. GLASSO must remain positive definite
#' and satisfy its full elementwise-\eqn{\ell_1} KKT equation. The tuning
#' parameter is explicit because the theoretical constant multiplying the
#' rate is not an observable default.
#'
#' @param fit A valid `elliptical_factor_fit`.
#' @param lambda Positive CLIME constraint or GLASSO penalty.
#' @param method Either `"sclime"` or `"sglasso"`.
#' @param solver_tol,solver_max_iter,initial_step,max_backtracking Solver
#'   controls; see `spatial_sign_precision()`.
#' @param strict If `TRUE`, fail on any solver or Woodbury certificate failure.
#'   If `FALSE`, return `estimate = NULL` and explicit diagnostics.
#'
#' @return An `elliptical_factor_precision_fit` containing the idiosyncratic
#'   and full precision estimates and complete solver certificates.
#' @references Xu, X., Ma, H., Wang, H. and Feng, L. (2025). arXiv:2512.19325.
#' @examples
#' x <- rbind(c(-3, -2), c(-2, -1), c(-1, 1),
#'            c(1, -1), c(2, 1), c(3, 2))
#' f <- spatial_sign_poet(x, factors = 0, threshold = 0.2)
#' elliptical_factor_precision(f, lambda = 0.5, method = "sglasso")
#' @export
elliptical_factor_precision <- function(
    fit, lambda, method = c("sclime", "sglasso"),
    solver_tol = 1e-7, solver_max_iter = 100000L,
    initial_step = 1, max_backtracking = 100L, strict = TRUE) {
  call <- match.call()
  if (!inherits(fit, "elliptical_factor_fit") || !isTRUE(fit$valid) ||
      is.null(fit$estimate)) {
    stop("`fit` must be a valid elliptical factor fit.", call. = FALSE)
  }
  method <- match.arg(method)
  strict <- .ch3ef_strict(strict)
  controls <- .ch3pp_validate_precision_controls(
    lambda, solver_tol, solver_max_iter, initial_step, max_backtracking
  )
  residual <- fit$components$idiosyncratic.raw
  p <- nrow(residual)
  failed <- function(stage, message, solver = NULL) {
    .ch3ef_warn_or_stop(message, strict)
    structure(list(
      estimate = NULL, idiosyncratic.precision = NULL,
      valid = FALSE, method = method, lambda = controls$lambda,
      source = fit,
      diagnostics = list(
        failure.stage = stage, failure = message, solver = solver,
        no.repair = paste(
          "No ridge, eigenvalue floor, pseudoinverse, constraint",
          "relaxation, or post-hoc KKT repair"
        ), call = call
      )
    ), class = "elliptical_factor_precision_fit")
  }

  if (method == "sclime") {
    solver <- tryCatch(cpp_ch3pp_sclime(
      residual, controls$lambda, controls$solver_tol,
      controls$solver_max_iter
    ), error = identity)
    if (inherits(solver, "condition")) {
      return(failed(
        "optimization error",
        paste("Elliptical-factor SCLIME failed:", conditionMessage(solver)),
        list(error = conditionMessage(solver))
      ))
    }
    valid <- isTRUE(solver$all_converged) &&
      all(solver$primal_violation <= controls$solver_tol) &&
      all(solver$dual_violation <= controls$solver_tol) &&
      all(solver$stationarity_residual <= controls$solver_tol) &&
      all(solver$relative_gap <= controls$solver_tol) &&
      is.finite(solver$final_feasibility_violation) &&
      solver$final_feasibility_violation <= controls$solver_tol
    certificate <- list(
      all.columns.certified = valid,
      column.converged = as.logical(solver$converged),
      iterations = as.integer(solver$iterations),
      primal.violation = as.numeric(solver$primal_violation),
      dual.violation = as.numeric(solver$dual_violation),
      stationarity.residual = as.numeric(solver$stationarity_residual),
      relative.gap = as.numeric(solver$relative_gap),
      raw.maximum.feasibility.violation = max(solver$primal_violation),
      symmetrized.feasibility.violation =
        as.numeric(solver$final_feasibility_violation),
      minimum.eigenvalue = as.numeric(solver$final_minimum_eigenvalue)
    )
    vu <- as.matrix(solver$solution)
    last <- list(raw = as.matrix(solver$raw_solution),
                 symmetrized = vu,
                 dual = as.matrix(solver$dual_solution))
  } else {
    solver <- tryCatch(cpp_ch3pp_sglasso(
      residual, controls$lambda, controls$solver_tol,
      controls$solver_max_iter, controls$initial_step,
      controls$max_backtracking
    ), error = identity)
    if (inherits(solver, "condition")) {
      return(failed(
        "optimization error",
        paste("Elliptical-factor SGLASSO failed:", conditionMessage(solver)),
        list(error = conditionMessage(solver))
      ))
    }
    valid <- isTRUE(solver$converged) &&
      !isTRUE(solver$backtracking_failed) &&
      is.finite(solver$minimum_eigenvalue) &&
      solver$minimum_eigenvalue > 0 &&
      is.finite(solver$kkt_residual) &&
      solver$kkt_residual <= controls$solver_tol
    certificate <- list(
      kkt.certified = valid, converged = isTRUE(solver$converged),
      backtracking.failed = isTRUE(solver$backtracking_failed),
      iterations = as.integer(solver$iterations),
      relative.update = as.numeric(solver$relative_update),
      kkt.residual = as.numeric(solver$kkt_residual),
      objective = as.numeric(solver$objective),
      minimum.eigenvalue = as.numeric(solver$minimum_eigenvalue),
      accepted.step = as.numeric(solver$accepted_step)
    )
    vu <- as.matrix(solver$solution)
    last <- list(symmetrized = vu, inverse = as.matrix(solver$inverse))
  }
  if (!valid) {
    return(failed(
      "optimization certificate",
      paste0(
        "Elliptical-factor ", toupper(method),
        " failed its feasibility/KKT certificate; no estimate was returned."
      ), certificate
    ))
  }

  gamma <- fit$components$leading.eigenvectors
  eigenvalues <- fit$components$leading.eigenvalues
  factors <- length(eigenvalues)
  if (factors > 0L) {
    if (any(!is.finite(eigenvalues)) || any(eigenvalues <= 0)) {
      return(failed(
        "Woodbury reconstruction",
        "Leading factor eigenvalues must be finite and positive."
      ))
    }
    vu.gamma <- vu %*% gamma
    middle <- diag(1 / eigenvalues, factors) + crossprod(gamma, vu.gamma)
    middle.factor <- tryCatch(chol((middle + t(middle)) / 2),
                              error = identity)
    if (inherits(middle.factor, "condition")) {
      return(failed(
        "Woodbury reconstruction",
        paste(
          "The Woodbury middle matrix is not positive definite; no",
          "pseudoinverse or ridge is applied."
        )
      ))
    }
    middle.inverse <- chol2inv(middle.factor)
    estimate <- vu - vu.gamma %*% middle.inverse %*% t(vu.gamma)
    middle.rcond <- rcond(middle)
  } else {
    estimate <- vu
    middle <- matrix(numeric(), 0L, 0L)
    middle.inverse <- middle
    middle.rcond <- NA_real_
  }
  estimate <- (estimate + t(estimate)) / 2
  if (any(!is.finite(estimate))) {
    return(failed(
      "Woodbury reconstruction",
      "The reconstructed precision is non-finite; no clipping is applied."
    ))
  }
  variable.names <- rownames(residual)
  if (!is.null(variable.names)) {
    dimnames(vu) <- dimnames(estimate) <-
      list(variable.names, variable.names)
  }
  minimum <- min(eigen(estimate, symmetric = TRUE, only.values = TRUE)$values)
  structure(list(
    estimate = estimate, idiosyncratic.precision = vu,
    valid = TRUE, method = method, lambda = controls$lambda,
    source = fit,
    components = list(
      factor.loadings = gamma,
      factor.eigenvalues = eigenvalues,
      woodbury.middle = middle,
      woodbury.middle.inverse = middle.inverse
    ),
    diagnostics = list(
      solver = certificate, last.iterate = last,
      residual.source = "raw unthresholded idiosyncratic complement",
      woodbury.middle.rcond = middle.rcond,
      full.minimum.eigenvalue = minimum,
      positive.definiteness.required = method == "sglasso",
      no.repair = paste(
        "No ridge, eigenvalue floor, pseudoinverse, constraint relaxation,",
        "or post-hoc KKT repair"
      ), call = call
    )
  ), class = "elliptical_factor_precision_fit")
}
