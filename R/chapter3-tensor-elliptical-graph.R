.ch3teg_mode_names <- function(dims) {
  paste0("mode", seq_along(dims))
}

.ch3teg_parse_data <- function(data) {
  if (is.list(data) && !is.data.frame(data)) {
    if (length(data) < 2L) {
      stop("A tensor list must contain at least two observations.",
           call. = FALSE)
    }
    first <- data[[1L]]
    if (!is.numeric(first)) {
      stop("Every tensor-list element must be numeric.", call. = FALSE)
    }
    dims <- dim(first)
    if (is.null(dims)) dims <- length(first)
    dims <- as.integer(dims)
    if (!length(dims) || anyNA(dims) || any(dims < 1L)) {
      stop("Every tensor must have positive mode dimensions.", call. = FALSE)
    }
    expected.dim <- dim(first)
    values <- lapply(seq_along(data), function(i) {
      value <- data[[i]]
      if (!is.numeric(value) || !identical(dim(value), expected.dim) ||
          length(value) != prod(dims)) {
        stop("All tensor-list elements must be numeric and have identical " %+%
               "dimensions.", call. = FALSE)
      }
      value <- as.numeric(value)
      if (anyNA(value) || any(!is.finite(value))) {
        stop("Tensor data must not contain NA, NaN, or infinite values.",
             call. = FALSE)
      }
      value
    })
    x <- do.call(rbind, values)
    input.type <- "list of column-major tensor arrays"
  } else {
    if (!is.numeric(data) || is.null(dim(data)) || length(dim(data)) < 2L) {
      stop("`data` must be a list of numeric tensors or an observation-first " %+%
             "numeric array.", call. = FALSE)
    }
    array.dim <- dim(data)
    if (anyNA(array.dim) || any(array.dim < 1L) || array.dim[1L] < 2L) {
      stop("An observation-first array must contain at least two " %+%
             "observations and positive mode dimensions.", call. = FALSE)
    }
    if (anyNA(data) || any(!is.finite(data))) {
      stop("Tensor data must not contain NA, NaN, or infinite values.",
           call. = FALSE)
    }
    dims <- as.integer(array.dim[-1L])
    x <- matrix(as.numeric(data), nrow = array.dim[1L])
    input.type <- "observation-first column-major numeric array"
  }

  pstar <- prod(as.double(dims))
  if (!is.finite(pstar) || pstar > .Machine$integer.max ||
      pstar != ncol(x)) {
    stop("The product of tensor mode dimensions is invalid or too large.",
         call. = FALSE)
  }
  storage.mode(x) <- "double"
  list(x = x, n = nrow(x), dims = dims, pstar = as.integer(pstar),
       input.type = input.type)
}

.ch3teg_validate_logical <- function(x, name) {
  if (!is.logical(x) || length(x) != 1L || is.na(x)) {
    stop(sprintf("`%s` must be TRUE or FALSE.", name), call. = FALSE)
  }
  x
}

.ch3teg_positive_integer <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || !is.finite(x) ||
      x < 1 || x != floor(x) || x > .Machine$integer.max) {
    stop(sprintf("`%s` must be one positive integer.", name), call. = FALSE)
  }
  as.integer(x)
}

.ch3teg_positive_scalar <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || !is.finite(x) ||
      x <= 0) {
    stop(sprintf("`%s` must be one finite, strictly positive number.", name),
         call. = FALSE)
  }
  as.numeric(x)
}

.ch3teg_expand_nonnegative <- function(x, length.out, name) {
  if (!is.numeric(x) || !(length(x) %in% c(1L, length.out)) || anyNA(x) ||
      any(!is.finite(x)) || any(x < 0)) {
    stop(sprintf("`%s` must contain either one or one-per-mode finite " %+%
                   "non-negative value.", name), call. = FALSE)
  }
  rep(as.numeric(x), length.out = length.out)
}

.ch3teg_normalize_frobenius <- function(x, context) {
  scale <- max(abs(x))
  if (!is.finite(scale) || scale <= 0) {
    stop(context %+% " has zero or non-finite Frobenius norm; no repair " %+%
           "was applied.", call. = FALSE)
  }
  scaled <- x / scale
  unit.norm <- sqrt(sum(scaled * scaled))
  if (!is.finite(unit.norm) || unit.norm <= 0) {
    stop(context %+% " has an invalid Frobenius norm; no repair was " %+%
           "applied.", call. = FALSE)
  }
  list(matrix = scaled / unit.norm, norm = scale * unit.norm)
}

.ch3teg_spd_sqrt <- function(x, context) {
  decomposition <- tryCatch(
    eigen(x, symmetric = TRUE),
    error = identity
  )
  if (inherits(decomposition, "condition") ||
      any(!is.finite(decomposition$values)) ||
      any(decomposition$values <= 0)) {
    stop(context %+% " is not strictly positive definite; no ridge, " %+%
           "eigenvalue floor, or pseudoinverse was used.", call. = FALSE)
  }
  root <- decomposition$vectors %*%
    (sqrt(decomposition$values) * t(decomposition$vectors))
  root <- (root + t(root)) / 2
  if (any(!is.finite(root)) || inherits(try(chol(root), silent = TRUE),
                                         "try-error")) {
    stop(context %+% " has no numerically valid SPD square root; no " %+%
           "repair was used.", call. = FALSE)
  }
  root
}

.ch3teg_failure <- function(stage, message, strict, partial, call) {
  if (strict) stop(message, call. = FALSE)
  warning(message, call. = FALSE)
  partial$estimate <- NULL
  partial$valid <- FALSE
  partial$failure.stage <- stage
  partial$failure <- message
  partial$call <- call
  partial$diagnostics <- c(
    partial$diagnostics,
    list(
      failure.stage = stage,
      failure = message,
      no.repair = paste(
        "No ridge, jitter, eigenvalue floor, pseudoinverse, objective",
        "relaxation, or post-hoc KKT repair was used."
      )
    )
  )
  structure(partial, class = "tensor_spatial_sign_precision_fit")
}

.ch3teg_center_and_sign <- function(parsed, center, median.controls) {
  x <- parsed$x
  if (is.null(center)) {
    center.vector <- spatial_median(
      x, tol = median.controls$tol,
      max_iter = median.controls$max_iter,
      zero_tol = median.controls$zero_tol,
      warn = FALSE
    )
    median.diagnostics <- list(
      estimated = TRUE,
      converged = isTRUE(attr(center.vector, "converged")),
      iterations = as.integer(attr(center.vector, "iterations")),
      objective = as.numeric(attr(center.vector, "objective")),
      relative.update = as.numeric(attr(center.vector, "relative_change")),
      equation.residual = as.numeric(attr(center.vector,
                                          "equation_residual"))
    )
  } else {
    if (!is.numeric(center) || length(center) != parsed$pstar ||
        anyNA(center) || any(!is.finite(center))) {
      stop("`center` must be NULL or a finite numeric tensor/vector with " %+%
             "one entry per tensor coordinate.", call. = FALSE)
    }
    if (!is.null(dim(center)) &&
        !identical(as.integer(dim(center)), parsed$dims)) {
      stop("An array-valued `center` must have the tensor mode dimensions.",
           call. = FALSE)
    }
    center.vector <- as.numeric(center)
    median.diagnostics <- list(
      estimated = FALSE, converged = NA, iterations = 0L,
      objective = NA_real_, relative.update = NA_real_,
      equation.residual = NA_real_
    )
  }
  if (any(!is.finite(center.vector))) {
    stop("The spatial median is non-finite; no tensor estimate was formed.",
         call. = FALSE)
  }
  signs.with.attributes <- spatial_sign(
    x, center = center.vector, zero_tol = median.controls$zero_tol
  )
  norms <- as.numeric(attr(signs.with.attributes, "norms"))
  n.zero <- as.integer(attr(signs.with.attributes, "n_zero"))
  signs <- matrix(as.numeric(signs.with.attributes), nrow = parsed$n,
                  ncol = parsed$pstar)
  if (any(!is.finite(signs))) {
    stop("Tensor spatial signs are non-finite; no estimate was formed.",
         call. = FALSE)
  }
  list(
    center.vector = as.numeric(center.vector),
    center.tensor = array(as.numeric(center.vector), dim = parsed$dims),
    signs = signs,
    tensor.signs = array(as.numeric(signs), dim = c(parsed$n, parsed$dims)),
    norms = norms,
    n.zero = n.zero,
    zero.indices = which(rowSums(abs(signs)) == 0),
    median = median.diagnostics
  )
}


#' Sparse tensor-elliptical precision matrices from tensor spatial signs
#'
#' Implements the Spatial-Sign Separate tensor lasso of Liu, Lu, Zhou, Feng,
#' and Wang. Observations may be a list of identically dimensioned numeric
#' arrays or one numeric array whose first dimension indexes observations.
#' R's column-major vectorization is used: the first tensor mode varies
#' fastest, and mode-\eqn{k} matricization follows Kolda--Bader ordering.
#'
#' The ordinary spatial median of the vectorized observations is used unless
#' `center` is supplied. Write \eqn{U_{i,(l)}} for the mode-\eqn{l}
#' matricization of tensor sign \eqn{U_i}. The mode-\eqn{l} pilot is
#' \deqn{\widetilde\Omega_l^{raw}=
#' \left\{\frac{p_l}{n}\sum_i U_{i,(l)}U_{i,(l)}^T\right\}^{-1}}
#' exactly when
#' \eqn{n p_* > p_l^2(p_l-1)/2}; otherwise it is \eqn{I_{p_l}}. Every pilot
#' is Frobenius-normalized. Other modes are whitened with the symmetric
#' positive-definite square roots of these normalized pilots, yielding
#' \deqn{\widehat S_k=\frac{p_k}{n}\sum_iV_i^{(k)}V_i^{(k)T}.}
#'
#' Each raw precision matrix minimizes the primary-paper objective
#' \deqn{\frac{1}{p_k}\{\operatorname{tr}(\widehat S_k\Omega)
#' -\log\det(\Omega)\}+\lambda_k\sum_{a\ne b}|\Omega_{ab}|.}
#' Thus the equivalent ordinary graphical-lasso penalty is
#' \eqn{\rho_k=p_k\lambda_k}; diagonal entries are not penalized. A
#' self-contained positive-definite proximal-gradient solver is used. A fit is
#' valid only if its objective is non-increasing, every raw solution is SPD,
#' the diagonal and off-diagonal KKT residuals do not exceed `solver_tol`, and
#' the relative update meets `solver_tol`. The reported normalized estimates
#' are raw solutions divided by their Frobenius norms; KKT certificates refer
#' to the raw solutions.
#'
#' No ridge, jitter, eigenvalue floor, pseudoinverse, or post-hoc repair is
#' applied. A singular pilot on the paper's inverse branch is therefore a
#' failure. With `strict = FALSE`, any such failure or any uncertified solver
#' run returns `estimate = NULL` and `valid = FALSE`, together with available
#' diagnostics; it never exposes a last iterate as a valid estimate.
#'
#' @param data A list of same-dimension numeric tensors, or a numeric array
#'   with observations in its first dimension. A numeric matrix represents
#'   \eqn{n} observations of an order-one tensor.
#' @param lambda One non-negative tuning value or one value per mode. The
#'   value is \eqn{\lambda_k} in the paper, not the rescaled glasso penalty.
#' @param center `NULL` for the ordinary sample spatial median, or a finite
#'   numeric vector/array with the tensor mode dimensions.
#' @param median_tol,median_max_iter,zero_tol Controls for the ordinary spatial
#'   median and tensor spatial signs. Exact zero residuals map to zero signs.
#' @param solver_tol Positive tolerance required simultaneously for KKT and
#'   relative-update certificates.
#' @param solver_max_iter Positive graphical-lasso iteration limit.
#' @param initial_step Positive initial proximal-gradient step.
#' @param max_backtracking Positive line-search reduction limit per iteration.
#' @param keep_whitened If `TRUE`, retain the vectorized tensors after
#'   other-mode pilot whitening for each target mode. This can be large.
#' @param strict If `TRUE`, stop on any median, pilot, or solver failure. If
#'   `FALSE`, return an explicitly invalid fit with no estimate.
#'
#' @return An object of class `tensor_spatial_sign_precision_fit`. A valid fit
#'   contains `estimate`, `raw.precision`, the center and tensor signs, mode
#'   dimensions, pilot source scatters and branches, pilot square roots,
#'   whitened mode scatters, paper and effective penalties, and complete SPD,
#'   descent, KKT, and relative-update certificates.
#'
#' @references
#' Liu, J., Lu, Z., Zhou, L., Feng, L., and Wang, Z. (2025). Tensor
#' Elliptical Graphic Model. arXiv:2508.00333.
#' \url{https://arxiv.org/abs/2508.00333}.
#'
#' @examples
#' x <- array(c(
#'   -2, -1, 1, 2, -1, 1, 2, -2,
#'   1, -2, 2, -1, 2, 1, -1, -2,
#'   2, 1, -2, -1, -2, 2, 1, -1,
#'   -1, 2, -2, 1, 1, -2, 2, -1
#' ), dim = c(8, 2, 2))
#' fit <- tensor_spatial_sign_precision(
#'   x, lambda = 1, center = array(0, c(2, 2))
#' )
#' fit$estimate
#'
#' @export
tensor_spatial_sign_precision <- function(
    data, lambda, center = NULL,
    median_tol = 1e-8, median_max_iter = 500L, zero_tol = 0,
    solver_tol = 1e-7, solver_max_iter = 10000L,
    initial_step = 1, max_backtracking = 100L,
    keep_whitened = FALSE, strict = TRUE) {
  call <- match.call()
  data.name <- deparse1(substitute(data))
  parsed <- .ch3teg_parse_data(data)
  median.controls <- .validate_iteration_controls(
    median_tol, median_max_iter, zero_tol
  )
  solver_tol <- .ch3teg_positive_scalar(solver_tol, "solver_tol")
  solver_max_iter <- .ch3teg_positive_integer(
    solver_max_iter, "solver_max_iter"
  )
  initial_step <- .ch3teg_positive_scalar(initial_step, "initial_step")
  max_backtracking <- .ch3teg_positive_integer(
    max_backtracking, "max_backtracking"
  )
  keep_whitened <- .ch3teg_validate_logical(keep_whitened, "keep_whitened")
  strict <- .ch3teg_validate_logical(strict, "strict")
  lambda <- .ch3teg_expand_nonnegative(lambda, length(parsed$dims), "lambda")
  mode.names <- .ch3teg_mode_names(parsed$dims)
  names(lambda) <- mode.names
  if (!is.null(center)) {
    if (!is.numeric(center) || length(center) != parsed$pstar ||
        anyNA(center) || any(!is.finite(center))) {
      stop("`center` must be NULL or a finite numeric tensor/vector with " %+%
             "one entry per tensor coordinate.", call. = FALSE)
    }
    if (!is.null(dim(center)) &&
        !identical(as.integer(dim(center)), parsed$dims)) {
      stop("An array-valued `center` must have the tensor mode dimensions.",
           call. = FALSE)
    }
  }
  effective.penalty <- lambda * parsed$dims
  if (any(!is.finite(effective.penalty))) {
    stop("The rescaled penalties `dims * lambda` must remain finite.",
         call. = FALSE)
  }
  names(effective.penalty) <- mode.names

  centered <- tryCatch(
    .ch3teg_center_and_sign(parsed, center, median.controls),
    error = identity
  )
  base <- list(
    estimate = NULL, valid = FALSE, data.name = data.name,
    n = parsed$n, dims = parsed$dims, pstar = parsed$pstar,
    lambda = lambda, input.type = parsed$input.type,
    vectorization = paste(
      "R column-major vec; mode 1 varies fastest; observations occupy",
      "the first dimension of returned tensor signs"
    ),
    call = call, diagnostics = list()
  )
  if (inherits(centered, "condition")) {
    return(.ch3teg_failure(
      "center and signs",
      paste0("Tensor spatial-sign preprocessing failed: ",
             conditionMessage(centered)),
      strict, base, call
    ))
  }
  base$center <- centered$center.tensor
  base$center.vector <- centered$center.vector
  base$signs <- centered$tensor.signs
  base$sign.matrix <- centered$signs
  base$diagnostics$spatial.median <- centered$median
  base$diagnostics$zero.residuals <- list(
    count = centered$n.zero, indices = centered$zero.indices,
    norms = centered$norms, zero.tolerance = median.controls$zero_tol
  )
  if (isTRUE(centered$median$estimated) &&
      !isTRUE(centered$median$converged)) {
    return(.ch3teg_failure(
      "spatial median",
      "The ordinary sample spatial median did not converge; no tensor " %+%
        "precision estimate was formed.",
      strict, base, call
    ))
  }

  pilot.scatter <- tryCatch(
    cpp_ch3teg_mode_crossproducts(centered$signs, parsed$dims),
    error = identity
  )
  if (inherits(pilot.scatter, "condition")) {
    return(.ch3teg_failure(
      "pilot scatter",
      paste0("Tensor pilot scatter construction failed: ",
             conditionMessage(pilot.scatter)),
      strict, base, call
    ))
  }
  names(pilot.scatter) <- mode.names
  pilot.raw <- pilot <- pilot.sqrt <- vector("list", length(parsed$dims))
  branch <- character(length(parsed$dims))
  pilot.raw.norm <- numeric(length(parsed$dims))
  names(pilot.raw) <- names(pilot) <- names(pilot.sqrt) <-
    names(branch) <- names(pilot.raw.norm) <- mode.names

  for (k in seq_along(parsed$dims)) {
    pk <- parsed$dims[k]
    inverse.branch <- as.double(parsed$n) * as.double(parsed$pstar) >
      as.double(pk)^2 * as.double(pk - 1L) / 2
    branch[k] <- if (inverse.branch) "inverse" else "identity"
    component <- tryCatch({
      if (inverse.branch) {
        factor <- chol(pilot.scatter[[k]])
        raw <- chol2inv(factor)
        if (any(!is.finite(raw))) {
          stop("the inverse contains non-finite entries")
        }
      } else {
        raw <- diag(pk)
      }
      normalized <- .ch3teg_normalize_frobenius(
        raw, paste0("Mode-", k, " pilot")
      )
      root <- .ch3teg_spd_sqrt(
        normalized$matrix, paste0("Mode-", k, " normalized pilot")
      )
      list(raw = raw, normalized = normalized$matrix,
           norm = normalized$norm, root = root)
    }, error = identity)
    if (inherits(component, "condition")) {
      base$pilot.scatter <- pilot.scatter
      base$pilot.branch <- branch
      return(.ch3teg_failure(
        paste0("pilot mode ", k),
        paste0("Mode-", k, " pilot failed on the paper's ", branch[k],
               " branch: ", conditionMessage(component),
               ". No inverse or SPD repair was substituted."),
        strict, base, call
      ))
    }
    pilot.raw[[k]] <- component$raw
    pilot[[k]] <- component$normalized
    pilot.sqrt[[k]] <- component$root
    pilot.raw.norm[k] <- component$norm
  }

  whitened <- tryCatch(
    cpp_ch3teg_whitened_mode_scatter(
      centered$signs, parsed$dims, pilot.sqrt, keep_whitened
    ),
    error = identity
  )
  base$pilot.scatter <- pilot.scatter
  base$pilot.raw <- pilot.raw
  base$pilot <- pilot
  base$pilot.sqrt <- pilot.sqrt
  base$pilot.branch <- branch
  base$pilot.raw.frobenius <- pilot.raw.norm
  if (inherits(whitened, "condition")) {
    return(.ch3teg_failure(
      "other-mode whitening",
      paste0("Tensor pilot whitening failed: ", conditionMessage(whitened)),
      strict, base, call
    ))
  }
  mode.scatter <- whitened$scatter
  names(mode.scatter) <- mode.names
  base$mode.scatter <- mode.scatter
  if (keep_whitened) {
    names(whitened$vectors) <- mode.names
    base$whitened.sign.matrix <- whitened$vectors
  }

  raw.precision <- estimate <- solver <- vector("list", length(parsed$dims))
  names(raw.precision) <- names(estimate) <- names(solver) <- mode.names
  for (k in seq_along(parsed$dims)) {
    fit <- tryCatch(
      cpp_ch3teg_offdiag_glasso(
        mode.scatter[[k]], effective.penalty[k], solver_tol,
        solver_max_iter, initial_step, max_backtracking
      ),
      error = identity
    )
    if (inherits(fit, "condition")) {
      base$raw.precision <- raw.precision
      base$solver <- solver
      base$effective.penalty <- effective.penalty
      return(.ch3teg_failure(
        paste0("solver mode ", k),
        paste0("Mode-", k, " off-diagonal graphical lasso failed: ",
               conditionMessage(fit), " No numerical repair was applied."),
        strict, base, call
      ))
    }
    fit$paper.objective <- fit$objective / parsed$dims[k]
    fit$paper.objective.history <- fit$objective_history / parsed$dims[k]
    fit$paper.lambda <- lambda[k]
    fit$effective.penalty <- effective.penalty[k]
    certified <- isTRUE(fit$converged) &&
      !isTRUE(fit$backtracking_failed) &&
      isTRUE(fit$objective_descent) &&
      is.finite(fit$minimum_eigenvalue) && fit$minimum_eigenvalue > 0 &&
      is.finite(fit$kkt_diagonal) && fit$kkt_diagonal <= solver_tol &&
      is.finite(fit$kkt_off_diagonal) &&
      fit$kkt_off_diagonal <= solver_tol &&
      is.finite(fit$relative_update) && fit$relative_update <= solver_tol
    fit$certified <- certified
    solver[[k]] <- fit
    if (!certified) {
      base$raw.precision <- raw.precision
      base$solver <- solver
      base$effective.penalty <- effective.penalty
      return(.ch3teg_failure(
        paste0("solver certificate mode ", k),
        paste0("Mode-", k, " graphical lasso did not satisfy every SPD, ",
               "objective-descent, KKT, and relative-update certificate ",
               "within the requested limits. No last iterate was returned ",
               "as an estimate."),
        strict, base, call
      ))
    }
    raw.precision[[k]] <- fit$solution
    estimate[[k]] <- .ch3teg_normalize_frobenius(
      fit$solution, paste0("Mode-", k, " raw precision")
    )$matrix
  }

  base$estimate <- estimate
  base$raw.precision <- raw.precision
  base$valid <- TRUE
  base$effective.penalty <- effective.penalty
  base$solver <- solver
  base$diagnostics$solver.tolerance <- solver_tol
  base$diagnostics$solver.max.iterations <- solver_max_iter
  base$diagnostics$all.modes.certified <- TRUE
  base$diagnostics$no.repair <- paste(
    "No ridge, jitter, eigenvalue floor, pseudoinverse, objective",
    "relaxation, or post-hoc KKT repair was used."
  )
  structure(base, class = "tensor_spatial_sign_precision_fit")
}


#' Threshold fitted tensor spatial-sign precision matrices
#'
#' Applies the book's corrected threshold rule
#' \eqn{\widehat\Omega_{ab} I\{|\widehat\Omega_{ab}|\ge\tau_k\}} to every
#' off-diagonal entry of a valid tensor spatial-sign precision fit. Equality is
#' retained and diagonal entries are always preserved. The arXiv display omits
#' the absolute-value bars, but its proof explicitly treats both positive and
#' negative edges; the absolute-value rule is therefore used here.
#'
#' @param fit A valid object returned by
#'   \code{\link{tensor_spatial_sign_precision}}.
#' @param tau One non-negative threshold or one threshold per mode.
#'
#' @return A list with thresholded modewise precision estimates, thresholds,
#'   and the original unthresholded estimates.
#'
#' @references
#' Liu, J., Lu, Z., Zhou, L., Feng, L., and Wang, Z. (2025). Tensor
#' Elliptical Graphic Model. arXiv:2508.00333.
#' \url{https://arxiv.org/abs/2508.00333}.
#'
#' @examples
#' x <- array(
#'   sin(seq_len(48)) + cos(seq_len(48) / 3), dim = c(12, 2, 2)
#' )
#' fit <- tensor_spatial_sign_precision(x, lambda = 1, center = c(0, 0, 0, 0))
#' threshold_tensor_spatial_sign_precision(fit, tau = 0.05)
#'
#' @export
threshold_tensor_spatial_sign_precision <- function(fit, tau) {
  call <- match.call()
  if (!inherits(fit, "tensor_spatial_sign_precision_fit") ||
      !isTRUE(fit$valid) || is.null(fit$estimate)) {
    stop("`fit` must be a valid tensor spatial-sign precision fit.",
         call. = FALSE)
  }
  tau <- .ch3teg_expand_nonnegative(tau, length(fit$estimate), "tau")
  names(tau) <- names(fit$estimate)
  thresholded <- lapply(seq_along(fit$estimate), function(k) {
    value <- fit$estimate[[k]]
    keep <- abs(value) >= tau[k]
    diag(keep) <- TRUE
    value[!keep] <- 0
    value
  })
  names(thresholded) <- names(fit$estimate)
  structure(
    list(
      estimate = thresholded, tau = tau, original = fit$estimate,
      valid = TRUE, comparison = "absolute value >= tau; diagonal retained",
      call = call
    ),
    class = "thresholded_tensor_spatial_sign_precision_fit"
  )
}
