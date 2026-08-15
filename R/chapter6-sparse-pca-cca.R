# Chapter 6: sparse PCA and sparse CCA.
#
# Every numerical tuning parameter in this file is user visible.  Iterative
# methods return a directly checkable fixed-point, ADMM, or block-KKT
# certificate.  An uncertified iterate is never promoted to a fitted model.

.c6spc_flag <- function(x, name) {
  if (!is.logical(x) || length(x) != 1L || is.na(x)) {
    stop(sprintf("`%s` must be TRUE or FALSE.", name), call. = FALSE)
  }
  x
}


.c6spc_count <- function(x, name, minimum = 1L, maximum = Inf) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || !is.finite(x) ||
      x != floor(x) || x < minimum || x > maximum ||
      x > .Machine$integer.max) {
    stop(sprintf(
      "`%s` must be one integer between %s and %s.",
      name, format(minimum), format(maximum)
    ), call. = FALSE)
  }
  as.integer(x)
}


.c6spc_scalar <- function(x, name, allow_zero = FALSE) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || !is.finite(x) ||
      if (allow_zero) x < 0 else x <= 0) {
    qualifier <- if (allow_zero) "non-negative" else "strictly positive"
    stop(sprintf("`%s` must be one finite %s number.", name, qualifier),
         call. = FALSE)
  }
  as.numeric(x)
}


.c6spc_data <- function(x, name = "x", min_rows = 2L) {
  x <- .as_data_matrix(x, name = name, min_rows = min_rows)
  variable.names <- colnames(x)
  if (is.null(variable.names)) {
    variable.names <- paste0("V", seq_len(ncol(x)))
    colnames(x) <- variable.names
  }
  if (anyNA(variable.names) || any(!nzchar(variable.names)) ||
      anyDuplicated(variable.names)) {
    stop(sprintf(
      "Column names of `%s`, when supplied, must be non-empty and unique.",
      name
    ), call. = FALSE)
  }
  x
}


.c6spc_named_vector <- function(x, p, variable.names, name,
                                positive = FALSE) {
  supplied.names <- names(x)
  x <- as.numeric(x)
  if (length(x) != p || anyNA(x) || any(!is.finite(x)) ||
      (positive && any(x <= 0))) {
    requirement <- if (positive) "finite strictly positive" else "finite"
    stop(sprintf("`%s` must contain %d %s values.", name, p, requirement),
         call. = FALSE)
  }
  if (!is.null(supplied.names)) {
    if (anyNA(supplied.names) || any(!nzchar(supplied.names)) ||
        anyDuplicated(supplied.names) ||
        !setequal(supplied.names, variable.names)) {
      stop(sprintf(
        "Names of `%s` must match the data column names exactly.", name
      ), call. = FALSE)
    }
    x <- x[match(variable.names, supplied.names)]
  }
  setNames(x, variable.names)
}


.c6spc_center_scale <- function(x, center, scale, name = "x") {
  n <- nrow(x)
  p <- ncol(x)
  variable.names <- colnames(x)
  if (is.character(center)) {
    center <- match.arg(center, c("mean", "none"))
    center.value <- if (center == "mean") colMeans(x) else numeric(p)
    center.source <- center
  } else {
    center.value <- .c6spc_named_vector(
      center, p, variable.names, paste0("center for ", name)
    )
    center.source <- "supplied"
  }
  center.value <- setNames(as.numeric(center.value), variable.names)
  centered <- sweep(x, 2L, center.value, "-", check.margin = FALSE)
  if (any(!is.finite(centered))) {
    stop(sprintf("Centered `%s` is outside the finite double range.", name),
         call. = FALSE)
  }

  if (is.logical(scale)) {
    scale <- .c6spc_flag(scale, paste0("scale for ", name))
    if (scale) {
      if (n < 2L) {
        stop("Sample standard deviations require at least two rows.",
             call. = FALSE)
      }
      scale.value <- sqrt(colSums(centered^2) / (n - 1))
      if (any(!is.finite(scale.value)) || any(scale.value <= 0)) {
        stop(sprintf(
          "`%s` contains a constant or non-finite column; no scale floor is applied.",
          name
        ), call. = FALSE)
      }
      scale.source <- "sample standard deviation (divisor n - 1)"
    } else {
      scale.value <- rep.int(1, p)
      scale.source <- "none"
    }
  } else {
    scale.value <- .c6spc_named_vector(
      scale, p, variable.names, paste0("scale for ", name), positive = TRUE
    )
    scale.source <- "supplied"
  }
  scale.value <- setNames(as.numeric(scale.value), variable.names)
  processed <- sweep(centered, 2L, scale.value, "/", check.margin = FALSE)
  if (any(!is.finite(processed))) {
    stop(sprintf("Scaled `%s` is outside the finite double range.", name),
         call. = FALSE)
  }
  dimnames(processed) <- dimnames(x)
  list(
    x = processed, center = center.value, scale = scale.value,
    center.source = center.source, scale.source = scale.source
  )
}


.c6spc_divisor <- function(covariance_divisor, n) {
  covariance_divisor <- match.arg(covariance_divisor, c("n", "n-1"))
  value <- if (covariance_divisor == "n") n else n - 1L
  if (value <= 0L) {
    stop("The selected covariance divisor is not positive.", call. = FALSE)
  }
  list(label = covariance_divisor, value = as.numeric(value))
}


.c6spc_operator <- function(operator, symmetry_tol, require_psd = TRUE) {
  symmetry_tol <- .c6spc_scalar(symmetry_tol, "symmetry_tol")
  operator <- as.matrix(operator)
  storage.mode(operator) <- "double"
  if (!is.numeric(operator) || nrow(operator) < 1L ||
      nrow(operator) != ncol(operator) || anyNA(operator) ||
      any(!is.finite(operator))) {
    stop("`operator` must be a finite non-empty square numeric matrix.",
         call. = FALSE)
  }
  p <- nrow(operator)
  rn <- rownames(operator)
  cn <- colnames(operator)
  if (!is.null(rn) && !is.null(cn) && !identical(rn, cn)) {
    stop("Row and column names of `operator` must agree exactly.",
         call. = FALSE)
  }
  variable.names <- if (!is.null(cn)) cn else if (!is.null(rn)) rn else
    paste0("V", seq_len(p))
  if (anyNA(variable.names) || any(!nzchar(variable.names)) ||
      anyDuplicated(variable.names)) {
    stop("Names of `operator` must be non-empty and unique.", call. = FALSE)
  }
  matrix.scale <- max(1, max(abs(operator)))
  symmetry.error <- max(abs(operator - t(operator)))
  if (!is.finite(symmetry.error) ||
      symmetry.error > symmetry_tol * matrix.scale) {
    stop("`operator` is not symmetric within `symmetry_tol`.",
         call. = FALSE)
  }
  evaluated <- (operator + t(operator)) / 2
  dimnames(evaluated) <- list(variable.names, variable.names)
  spectrum <- eigen(evaluated, symmetric = TRUE, only.values = TRUE)$values
  spectral.scale <- max(1, max(abs(spectrum)))
  if (require_psd && min(spectrum) < -symmetry_tol * spectral.scale) {
    stop("`operator` is not positive semidefinite within `symmetry_tol`; no eigenvalue floor is applied.",
         call. = FALSE)
  }
  list(
    matrix = evaluated, values = as.numeric(spectrum),
    minimum.eigenvalue = min(spectrum), symmetry.error = symmetry.error,
    tolerance = symmetry_tol, variable.names = variable.names
  )
}


.c6spc_anchor <- function(loadings, paired = NULL) {
  loadings <- as.matrix(loadings)
  anchors <- integer(ncol(loadings))
  for (j in seq_len(ncol(loadings))) {
    magnitude <- abs(loadings[, j])
    anchor <- which(magnitude == max(magnitude))[1L]
    anchors[j] <- anchor
    if (loadings[anchor, j] < 0) {
      loadings[, j] <- -loadings[, j]
      if (!is.null(paired)) paired[, j] <- -paired[, j]
    }
  }
  list(loadings = loadings, paired = paired, anchors = anchors)
}


.c6spc_repeat_groups <- function(values, tolerance) {
  if (!length(values)) return(list())
  tolerance <- .c6spc_scalar(tolerance, "symmetry_tol")
  threshold <- tolerance * max(1, max(abs(values)))
  groups <- split(seq_along(values), cumsum(c(TRUE, abs(diff(values)) > threshold)))
  Filter(function(index) length(index) > 1L, groups)
}


.c6spc_fit <- function(method, method.class, eigenvalues, loadings, scores,
                       center, operator, rank, n, p, variable.names,
                       diagnostics, call, valid = TRUE, extra = list()) {
  structure(c(list(
    valid = isTRUE(valid),
    method = method,
    eigenvalues = eigenvalues,
    loadings = loadings,
    scores = scores,
    center = center,
    operator = operator,
    rank = as.integer(rank),
    n = as.integer(n),
    p = as.integer(p),
    variable.names = variable.names,
    diagnostics = c(diagnostics, list(
      no.implicit.regularization = paste(
        "no ridge, pseudoinverse, eigenvalue floor, denominator floor,",
        "or hidden tuning parameter is used"
      ),
      loading.sign.rule = paste(
        "the first coordinate attaining the largest absolute loading",
        "is positive"
      ),
      repeated.value.contract = paste(
        "individual directions are not identified at repeated spectral",
        "values; compare the corresponding orthogonal projector"
      )
    )),
    call = call
  ), extra), class = c(method.class, "hd_pca_fit", "list"))
}


.c6spc_invalid <- function(method, method.class, n, p, variable.names,
                           center, operator, requested.rank, stage, message,
                           call, diagnostics = list(), extra = list()) {
  loadings <- matrix(numeric(), p, 0L,
                     dimnames = list(variable.names, character()))
  scores <- if (is.na(n)) NULL else
    matrix(numeric(), n, 0L, dimnames = list(NULL, character()))
  .c6spc_fit(
    method = method, method.class = method.class,
    eigenvalues = numeric(), loadings = loadings, scores = scores,
    center = center, operator = operator, rank = 0L, n = n, p = p,
    variable.names = variable.names,
    diagnostics = c(diagnostics, list(
      requested.rank = as.integer(requested.rank),
      failure.stage = stage, failure.message = message,
      certified = FALSE
    )),
    call = call, valid = FALSE, extra = extra
  )
}


.c6spc_fail <- function(strict, message, ...) {
  if (strict) stop(message, call. = FALSE)
  warning(message, call. = FALSE)
  .c6spc_invalid(message = message, ...)
}


.c6spc_initial <- function(initial, p, components, name = "initial") {
  if (is.null(initial)) return(matrix(numeric(), 0L, 0L))
  initial <- as.matrix(initial)
  storage.mode(initial) <- "double"
  if (!identical(dim(initial), c(p, components)) || anyNA(initial) ||
      any(!is.finite(initial))) {
    stop(sprintf("`%s` must be a finite %d by %d matrix.",
                 name, p, components), call. = FALSE)
  }
  if (any(colSums(initial^2) == 0)) {
    stop(sprintf("Every column of `%s` must be nonzero.", name),
         call. = FALSE)
  }
  initial
}


.c6spc_sparsity <- function(sparsity, components, p, name = "sparsity") {
  if (!is.numeric(sparsity) || !length(sparsity) || anyNA(sparsity) ||
      any(!is.finite(sparsity)) || any(sparsity != floor(sparsity)) ||
      any(sparsity < 1) || any(sparsity > p)) {
    stop(sprintf("`%s` must contain integers between 1 and %d.", name, p),
         call. = FALSE)
  }
  if (length(sparsity) == 1L) sparsity <- rep.int(sparsity, components)
  if (length(sparsity) != components) {
    stop(sprintf("`%s` must have length one or `components`.", name),
         call. = FALSE)
  }
  as.integer(sparsity)
}


.c6spc_l1_bounds <- function(bound, components, p, name) {
  if (!is.numeric(bound) || !length(bound) || anyNA(bound) ||
      any(!is.finite(bound))) {
    stop(sprintf("`%s` must contain finite numeric values.", name),
         call. = FALSE)
  }
  if (length(bound) == 1L) bound <- rep.int(bound, components)
  if (length(bound) != components) {
    stop(sprintf("`%s` must have length one or `components`.", name),
         call. = FALSE)
  }
  upper <- sqrt(p)
  if (any(bound < 1) || any(bound > upper)) {
    stop(sprintf("`%s` must lie in [1, sqrt(%d)].", name, p),
         call. = FALSE)
  }
  as.numeric(bound)
}


.c6spc_solver_ok <- function(result, components) {
  isTRUE(result$converged) &&
    identical(as.integer(result$completed), as.integer(components)) &&
    length(result$certificates) == components &&
    all(vapply(result$certificates, function(x) {
      isTRUE(x$converged) && isTRUE(x$certified)
    }, logical(1)))
}


#' Truncated-power sparse principal components
#'
#' Applies the truncated-power iteration to a supplied positive-semidefinite
#' operator.  At every iteration only the `sparsity` largest absolute entries
#' of the matrix-vector product are retained; exact ties are resolved by the
#' smallest coordinate index.  Multiple components use the Mackey deflation
#' \eqn{(I-vv')M(I-vv')}.
#'
#' The returned certificate reports unit-norm and support violations, the
#' sign-invariant truncated fixed-point residual, objective monotonicity, and
#' iteration count.  Because an operator rather than observations is supplied,
#' `scores`, `center`, and `n` are `NULL`, `NULL`, and `NA`.
#'
#' @param operator Finite symmetric positive-semidefinite matrix.
#' @param sparsity Required support size, scalar or one value per component.
#' @param components Number of components.
#' @param initial Optional `p` by `components` deterministic starting matrix.
#' @param solver_tol Positive fixed-point tolerance.
#' @param solver_max_iter Positive maximum number of iterations per component.
#' @param symmetry_tol Positive relative symmetry/PSD certification tolerance.
#' @param strict If `TRUE`, an uncertified solver stops; otherwise a warning and
#'   an invalid `hd_pca_fit` are returned.
#' @param keep_operator Retain the supplied operator in the fit.
#'
#' @return A `truncated_power_pca_fit` inheriting from `hd_pca_fit`.
#' @export
#'
#' @references
#' Yuan, X.-T. and Zhang, T. (2013). Truncated power method for sparse
#' eigenvalue problems. *Journal of Machine Learning Research*, 14, 899--925.
#'
#' @examples
#' truncated_power_pca(diag(c(4, 2, 1)), sparsity = 1)
truncated_power_pca <- function(
    operator, sparsity, components = 1L, initial = NULL,
    solver_tol = 1e-8, solver_max_iter = 1000L,
    symmetry_tol = sqrt(.Machine$double.eps), strict = TRUE,
    keep_operator = TRUE) {
  call <- match.call()
  strict <- .c6spc_flag(strict, "strict")
  keep_operator <- .c6spc_flag(keep_operator, "keep_operator")
  solver_tol <- .c6spc_scalar(solver_tol, "solver_tol")
  solver_max_iter <- .c6spc_count(solver_max_iter, "solver_max_iter")
  checked <- .c6spc_operator(operator, symmetry_tol, require_psd = TRUE)
  p <- nrow(checked$matrix)
  components <- .c6spc_count(components, "components", maximum = p)
  sparsity <- .c6spc_sparsity(sparsity, components, p)
  initial <- .c6spc_initial(initial, p, components)
  result <- cpp_ch6spc_tpm(
    checked$matrix, sparsity, initial, solver_tol, solver_max_iter
  )
  retained.operator <- if (keep_operator) checked$matrix else NULL
  base.diagnostics <- list(
    solver = result$certificates,
    solver.tolerance = solver_tol,
    solver.maximum.iterations = solver_max_iter,
    sparsity = sparsity,
    deflation = "(I - vv') M (I - vv')",
    deterministic.support.ties = "decreasing absolute value, then coordinate index",
    operator.symmetry.error = checked$symmetry.error,
    operator.minimum.eigenvalue = checked$minimum.eigenvalue,
    operator.symmetry.tolerance = checked$tolerance
  )
  if (!.c6spc_solver_ok(result, components)) {
    message <- sprintf(
      "Truncated-power component %d did not pass its fixed-point certificate: %s",
      as.integer(result$failure_component), as.character(result$failure_message)
    )
    return(.c6spc_fail(
      strict, message,
      method = "truncated-power sparse PCA",
      method.class = "truncated_power_pca_fit",
      n = NA_integer_, p = p, variable.names = checked$variable.names,
      center = NULL, operator = retained.operator,
      requested.rank = components, stage = "truncated-power solver",
      call = call, diagnostics = base.diagnostics
    ))
  }
  loadings <- as.matrix(result$loadings)
  anchored <- .c6spc_anchor(loadings)
  loadings <- anchored$loadings
  component.names <- paste0("SPC", seq_len(components))
  dimnames(loadings) <- list(checked$variable.names, component.names)
  eigenvalues <- vapply(seq_len(components), function(j) {
    as.numeric(crossprod(loadings[, j], checked$matrix %*% loadings[, j]))
  }, numeric(1))
  names(eigenvalues) <- component.names
  .c6spc_fit(
    method = "truncated-power sparse PCA",
    method.class = "truncated_power_pca_fit",
    eigenvalues = eigenvalues, loadings = loadings, scores = NULL,
    center = NULL, operator = retained.operator, rank = components,
    n = NA_integer_, p = p, variable.names = checked$variable.names,
    diagnostics = c(base.diagnostics, list(
      certified = TRUE,
      sign.anchor.coordinates = setNames(anchored$anchors, component.names),
      eigenvalue.definition = "Rayleigh value in the original supplied operator",
      completed.components = as.integer(result$completed)
    )),
    call = call,
    extra = list(deflated.operator = as.matrix(result$deflated_operator))
  )
}


#' Sparse spatial-sign principal component analysis
#'
#' Forms the Chapter 1 sample spatial-sign covariance matrix and applies the
#' certified truncated-power solver.  Zero residuals have the same explicit
#' semantics as [spatial_sign_pca()]; scores project the centered original
#' observations rather than their signs.
#'
#' @param x Numeric matrix or data frame with observations in rows.
#' @param sparsity Required support size, scalar or one value per component.
#' @param components Number of sparse components.
#' @param center Numeric center or one of `"spatial"`, `"mean"`, and `"none"`.
#' @param divisor `"n"` or the number of nonzero residuals.
#' @param zero_action Either retain exact zero sign contributions or reject.
#' @param initial Optional deterministic loading starts.
#' @param median_tol,median_max_iter Spatial-median controls.
#' @param zero_tol Non-negative zero-residual tolerance.
#' @param solver_tol,solver_max_iter Truncated-power controls.
#' @param symmetry_tol Positive relative operator certification tolerance.
#' @param strict If `TRUE`, an uncertified location or sparse solver stops;
#'   otherwise an invalid fit is returned with a warning.
#' @param keep_operator Retain the SSCM.
#'
#' @return A `sparse_spatial_sign_pca_fit` inheriting from `hd_pca_fit`.
#' @export
#'
#' @references
#' Zhao, Y., Wang, L., and Feng, L. (2024). Spatial-sign based high-dimensional
#' principal component analysis. arXiv:2409.13267.
#'
#' @examples
#' x <- rbind(c(4, 0, 0), c(-4, 0, 0), c(0, 2, 0), c(0, -2, 0))
#' sparse_spatial_sign_pca(x, sparsity = 1, center = "none")
sparse_spatial_sign_pca <- function(
    x, sparsity, components = 1L,
    center = c("spatial", "mean", "none"),
    divisor = c("n", "nonzero"), zero_action = c("zero", "error"),
    initial = NULL, median_tol = 1e-8, median_max_iter = 500L,
    zero_tol = 0, solver_tol = 1e-8, solver_max_iter = 1000L,
    symmetry_tol = sqrt(.Machine$double.eps), strict = TRUE,
    keep_operator = TRUE) {
  call <- match.call()
  if (missing(center)) center <- "spatial"
  strict <- .c6spc_flag(strict, "strict")
  keep_operator <- .c6spc_flag(keep_operator, "keep_operator")
  x <- .c6spc_data(x)
  n <- nrow(x)
  p <- ncol(x)
  components <- .c6spc_count(
    components, "components", maximum = min(n, p)
  )
  sparsity <- .c6spc_sparsity(sparsity, components, p)
  initial <- .c6spc_initial(initial, p, components)
  divisor <- match.arg(divisor)
  zero_action <- match.arg(zero_action)
  median_tol <- .c6spc_scalar(median_tol, "median_tol")
  median_max_iter <- .c6spc_count(median_max_iter, "median_max_iter")
  zero_tol <- .validate_zero_tol(zero_tol)
  solver_tol <- .c6spc_scalar(solver_tol, "solver_tol")
  solver_max_iter <- .c6spc_count(solver_max_iter, "solver_max_iter")

  center.fit <- tryCatch(
    .c6rs_center(
      x, center, median_tol, median_max_iter, zero_tol,
      allow_lts = FALSE
    ),
    error = identity
  )
  if (inherits(center.fit, "error")) {
    message <- conditionMessage(center.fit)
    return(.c6spc_fail(
      strict, message,
      method = "sparse spatial-sign PCA",
      method.class = "sparse_spatial_sign_pca_fit",
      n = n, p = p, variable.names = colnames(x), center = NULL,
      operator = NULL, requested.rank = components,
      stage = "spatial-median center", call = call,
      diagnostics = list(
        sparsity = sparsity, median.tolerance = median_tol,
        median.maximum.iterations = median_max_iter
      )
    ))
  }
  signs <- spatial_sign(x, center = center.fit$value, zero_tol = zero_tol)
  n.zero <- as.integer(attr(signs, "n_zero"))
  if (n.zero > 0L && zero_action == "error") {
    stop(sprintf("The centered sample contains %d zero residual(s).", n.zero),
         call. = FALSE)
  }
  denominator <- if (divisor == "n") n else n - n.zero
  if (denominator <= 0L) {
    stop("The nonzero-residual divisor is zero.", call. = FALSE)
  }
  operator <- crossprod(unclass(signs)) / denominator
  dimnames(operator) <- list(colnames(x), colnames(x))
  checked <- .c6spc_operator(operator, symmetry_tol, require_psd = TRUE)
  result <- cpp_ch6spc_tpm(
    checked$matrix, sparsity, initial, solver_tol, solver_max_iter
  )
  retained.operator <- if (keep_operator) checked$matrix else NULL
  base.diagnostics <- list(
    center.source = center.fit$source,
    center = center.fit$diagnostics,
    divisor.rule = divisor, divisor.value = as.integer(denominator),
    zero.action = zero_action, zero.tolerance = zero_tol,
    n.zero.residuals = n.zero, operator.trace = sum(diag(operator)),
    sparsity = sparsity, solver = result$certificates,
    solver.tolerance = solver_tol,
    solver.maximum.iterations = solver_max_iter,
    deflation = "(I - vv') M (I - vv')",
    deterministic.support.ties = "decreasing absolute value, then coordinate index",
    score.definition = "centered original observations times sparse loadings"
  )
  if (!.c6spc_solver_ok(result, components)) {
    message <- sprintf(
      "Sparse spatial-sign component %d did not pass its fixed-point certificate: %s",
      as.integer(result$failure_component), as.character(result$failure_message)
    )
    return(.c6spc_fail(
      strict, message,
      method = "sparse spatial-sign PCA",
      method.class = "sparse_spatial_sign_pca_fit",
      n = n, p = p, variable.names = colnames(x),
      center = center.fit$value, operator = retained.operator,
      requested.rank = components, stage = "truncated-power solver",
      call = call, diagnostics = base.diagnostics
    ))
  }
  anchored <- .c6spc_anchor(as.matrix(result$loadings))
  loadings <- anchored$loadings
  component.names <- paste0("SPC", seq_len(components))
  dimnames(loadings) <- list(colnames(x), component.names)
  scores <- .c6rs_scores(x, center.fit$value, loadings)
  eigenvalues <- vapply(seq_len(components), function(j) {
    as.numeric(crossprod(loadings[, j], operator %*% loadings[, j]))
  }, numeric(1))
  names(eigenvalues) <- component.names
  .c6spc_fit(
    method = "sparse spatial-sign PCA",
    method.class = "sparse_spatial_sign_pca_fit",
    eigenvalues = eigenvalues, loadings = loadings, scores = scores,
    center = center.fit$value, operator = retained.operator,
    rank = components, n = n, p = p, variable.names = colnames(x),
    diagnostics = c(base.diagnostics, list(
      certified = TRUE,
      sign.anchor.coordinates = setNames(anchored$anchors, component.names),
      operator.symmetry.error = checked$symmetry.error,
      operator.minimum.eigenvalue = checked$minimum.eigenvalue
    )),
    call = call,
    extra = list(deflated.operator = as.matrix(result$deflated_operator))
  )
}


#' Fantope projection and selection sparse PCA
#'
#' Solves
#' \deqn{\max_H \langle S,H\rangle-\tau\lVert H\rVert_{1,1},\quad
#' 0\preceq H\preceq I,\quad \mathrm{tr}(H)=r}
#' by a certified two-block ADMM.  The `H` update is the Euclidean Fantope
#' projection and the split update is entrywise soft thresholding.  The
#' returned `relaxed.projector` is the primary convex estimate; `loadings` are
#' its leading `rank` eigenvectors and are not claimed to be individually
#' identified across repeated eigenvalues.
#'
#' @param x Numeric data with observations in rows.
#' @param rank Required Fantope trace/rank.
#' @param tau Required non-negative entrywise lasso penalty.
#' @param center Numeric center or `"mean"`/`"none"`.
#' @param scale Logical or supplied positive scale vector.
#' @param covariance_divisor `"n"` (book convention) or `"n-1"`.
#' @param rho Positive ADMM penalty.
#' @param solver_tol Positive scaled primal/dual residual tolerance.
#' @param solver_max_iter Positive ADMM iteration limit.
#' @param symmetry_tol Positive operator certification tolerance.
#' @param strict If `TRUE`, failure is an error; otherwise return an invalid fit.
#' @param keep_operator Retain the covariance operator.
#'
#' @return A `fantope_pca_fit` inheriting from `hd_pca_fit`.
#' @export
#'
#' @references
#' Vu, V. Q., Cho, J., Lei, J., and Rohe, K. (2013). Fantope projection and
#' selection: a near-optimal convex relaxation of sparse PCA. *NeurIPS 26*.
#'
#' @examples
#' x <- rbind(c(3, 0), c(-3, 0), c(0, 1), c(0, -1))
#' fantope_pca(x, rank = 1, tau = 0.1)
fantope_pca <- function(
    x, rank, tau, center = c("mean", "none"), scale = FALSE,
    covariance_divisor = c("n", "n-1"), rho = 1,
    solver_tol = 1e-7, solver_max_iter = 5000L,
    symmetry_tol = sqrt(.Machine$double.eps), strict = TRUE,
    keep_operator = TRUE) {
  call <- match.call()
  if (missing(center)) center <- "mean"
  strict <- .c6spc_flag(strict, "strict")
  keep_operator <- .c6spc_flag(keep_operator, "keep_operator")
  x <- .c6spc_data(x)
  n <- nrow(x)
  p <- ncol(x)
  rank <- .c6spc_count(rank, "rank", maximum = p)
  tau <- .c6spc_scalar(tau, "tau", allow_zero = TRUE)
  rho <- .c6spc_scalar(rho, "rho")
  solver_tol <- .c6spc_scalar(solver_tol, "solver_tol")
  solver_max_iter <- .c6spc_count(solver_max_iter, "solver_max_iter")
  processed <- .c6spc_center_scale(x, center, scale)
  divisor <- .c6spc_divisor(covariance_divisor, n)
  operator <- crossprod(processed$x) / divisor$value
  dimnames(operator) <- list(colnames(x), colnames(x))
  checked <- .c6spc_operator(operator, symmetry_tol, require_psd = TRUE)
  result <- cpp_ch6spc_fantope(
    checked$matrix, rank, tau, rho, solver_tol, solver_max_iter
  )
  retained.operator <- if (keep_operator) checked$matrix else NULL
  base.diagnostics <- list(
    tuning.source = "explicit rank and tau",
    tau = tau, rho = rho,
    solver.tolerance = solver_tol,
    solver.maximum.iterations = solver_max_iter,
    iterations = as.integer(result$iterations),
    primal.residual = as.numeric(result$primal_residual),
    scaled.primal.residual = as.numeric(result$scaled_primal_residual),
    dual.residual = as.numeric(result$dual_residual),
    scaled.dual.residual = as.numeric(result$scaled_dual_residual),
    trace.residual = as.numeric(result$trace_residual),
    minimum.relaxed.eigenvalue = as.numeric(result$minimum_eigenvalue),
    maximum.relaxed.eigenvalue = as.numeric(result$maximum_eigenvalue),
    fantope.lower.violation = as.numeric(result$lower_violation),
    fantope.upper.violation = as.numeric(result$upper_violation),
    objective = as.numeric(result$objective),
    centering = processed$center.source,
    scaling = processed$scale.source,
    scale = processed$scale,
    covariance.divisor = divisor$label,
    covariance.divisor.value = divisor$value,
    operator.symmetry.error = checked$symmetry.error,
    operator.minimum.eigenvalue = checked$minimum.eigenvalue
  )
  certified <- isTRUE(result$converged) && isTRUE(result$certified)
  if (!certified) {
    message <- paste(
      "Fantope ADMM did not pass its primal/dual and Fantope certificate",
      sprintf("after %d iterations.", as.integer(result$iterations))
    )
    return(.c6spc_fail(
      strict, message,
      method = "Fantope projection and selection PCA",
      method.class = "fantope_pca_fit",
      n = n, p = p, variable.names = colnames(x),
      center = processed$center, operator = retained.operator,
      requested.rank = rank, stage = "Fantope ADMM",
      call = call, diagnostics = base.diagnostics,
      extra = list(scale = processed$scale, relaxed.projector = NULL)
    ))
  }
  relaxed <- as.matrix(result$H)
  dimnames(relaxed) <- list(colnames(x), colnames(x))
  decomposition <- eigen(relaxed, symmetric = TRUE)
  selected <- seq_len(rank)
  anchored <- .c6spc_anchor(decomposition$vectors[, selected, drop = FALSE])
  loadings <- anchored$loadings
  component.names <- paste0("FPC", selected)
  dimnames(loadings) <- list(colnames(x), component.names)
  scores <- processed$x %*% loadings
  dimnames(scores) <- list(rownames(x), component.names)
  eigenvalues <- vapply(selected, function(j) {
    as.numeric(crossprod(loadings[, j], operator %*% loadings[, j]))
  }, numeric(1))
  names(eigenvalues) <- component.names
  relaxed.values <- as.numeric(decomposition$values)
  .c6spc_fit(
    method = "Fantope projection and selection PCA",
    method.class = "fantope_pca_fit",
    eigenvalues = eigenvalues, loadings = loadings, scores = scores,
    center = processed$center, operator = retained.operator,
    rank = rank, n = n, p = p, variable.names = colnames(x),
    diagnostics = c(base.diagnostics, list(
      certified = TRUE,
      relaxed.eigenvalues = relaxed.values,
      repeated.relaxed.eigenvalue.groups =
        .c6spc_repeat_groups(relaxed.values, symmetry_tol),
      sign.anchor.coordinates = setNames(anchored$anchors, component.names),
      eigenvalue.definition = "Rayleigh variance in the sample operator",
      relaxed.estimate.contract = paste(
        "relaxed.projector is the primary convex Fantope estimate;",
        "leading vectors are a derived representation"
      )
    )),
    call = call,
    extra = list(scale = processed$scale, relaxed.projector = relaxed,
                 split.matrix = as.matrix(result$Z))
  )
}


#' Penalized-matrix-decomposition sparse PCA
#'
#' For each residual data matrix `R`, solves the PMD rank-one problem
#' \deqn{\max_{u,v} u'Rv,\quad \lVert u\rVert_2\leq1,\quad
#' \lVert v\rVert_2\leq1,\quad \lVert v\rVert_1\leq c,}
#' then applies `R <- R - d u v'`.  The supplied `l1_bound` is the actual
#' \eqn{c\in[1,\sqrt p]}, not a rescaled UI fraction.  The alternating solver
#' reports both block-KKT and fixed-point residuals.
#'
#' @param x Numeric data with observations in rows.
#' @param l1_bound Required PMD loading bound, scalar or one per component.
#' @param components Number of components.
#' @param center Numeric center or `"mean"`/`"none"`.
#' @param scale Logical or supplied positive scale vector.
#' @param covariance_divisor Divisor used to report component variances.
#' @param initial Optional deterministic loading starts.
#' @param solver_tol,solver_max_iter Alternating-solver controls.
#' @param symmetry_tol Positive covariance certification tolerance.
#' @param strict If `TRUE`, failure is an error; otherwise return an invalid fit.
#' @param keep_operator Retain the processed sample covariance.
#'
#' @return A `pmd_sparse_pca_fit` inheriting from `hd_pca_fit`.
#' @export
#'
#' @references
#' Witten, D. M., Tibshirani, R., and Hastie, T. (2009). A penalized matrix
#' decomposition, with applications to sparse principal components and
#' canonical correlation analysis. *Biostatistics*, 10, 515--534.
#'
#' @examples
#' x <- rbind(c(3, 0), c(-3, 0), c(0, 1), c(0, -1))
#' pmd_sparse_pca(x, l1_bound = 1)
pmd_sparse_pca <- function(
    x, l1_bound, components = 1L, center = c("mean", "none"),
    scale = FALSE, covariance_divisor = c("n", "n-1"),
    initial = NULL, solver_tol = 1e-8, solver_max_iter = 1000L,
    symmetry_tol = sqrt(.Machine$double.eps), strict = TRUE,
    keep_operator = TRUE) {
  call <- match.call()
  if (missing(center)) center <- "mean"
  strict <- .c6spc_flag(strict, "strict")
  keep_operator <- .c6spc_flag(keep_operator, "keep_operator")
  x <- .c6spc_data(x)
  n <- nrow(x)
  p <- ncol(x)
  components <- .c6spc_count(
    components, "components", maximum = min(n, p)
  )
  l1.bound <- .c6spc_l1_bounds(l1_bound, components, p, "l1_bound")
  initial <- .c6spc_initial(initial, p, components)
  solver_tol <- .c6spc_scalar(solver_tol, "solver_tol")
  solver_max_iter <- .c6spc_count(solver_max_iter, "solver_max_iter")
  processed <- .c6spc_center_scale(x, center, scale)
  divisor <- .c6spc_divisor(covariance_divisor, n)
  operator <- crossprod(processed$x) / divisor$value
  dimnames(operator) <- list(colnames(x), colnames(x))
  checked <- .c6spc_operator(operator, symmetry_tol, require_psd = TRUE)
  result <- cpp_ch6spc_pmd_pca(
    processed$x, l1.bound, initial, solver_tol, solver_max_iter
  )
  retained.operator <- if (keep_operator) checked$matrix else NULL
  base.diagnostics <- list(
    tuning.source = "explicit l1_bound",
    l1.bound = l1.bound,
    solver = result$certificates,
    solver.tolerance = solver_tol,
    solver.maximum.iterations = solver_max_iter,
    deflation = "R <- R - d u v'",
    centering = processed$center.source,
    scaling = processed$scale.source,
    scale = processed$scale,
    covariance.divisor = divisor$label,
    covariance.divisor.value = divisor$value,
    operator.symmetry.error = checked$symmetry.error,
    operator.minimum.eigenvalue = checked$minimum.eigenvalue
  )
  if (!.c6spc_solver_ok(result, components)) {
    message <- sprintf(
      "PMD sparse-PCA component %d did not pass its block certificate: %s",
      as.integer(result$failure_component), as.character(result$failure_message)
    )
    return(.c6spc_fail(
      strict, message,
      method = "PMD sparse PCA", method.class = "pmd_sparse_pca_fit",
      n = n, p = p, variable.names = colnames(x),
      center = processed$center, operator = retained.operator,
      requested.rank = components, stage = "PMD alternating solver",
      call = call, diagnostics = base.diagnostics,
      extra = list(scale = processed$scale)
    ))
  }
  left <- as.matrix(result$left_vectors)
  anchored <- .c6spc_anchor(as.matrix(result$loadings), paired = left)
  loadings <- anchored$loadings
  left <- anchored$paired
  component.names <- paste0("SPC", seq_len(components))
  dimnames(loadings) <- list(colnames(x), component.names)
  dimnames(left) <- list(rownames(x), component.names)
  scores <- processed$x %*% loadings
  dimnames(scores) <- list(rownames(x), component.names)
  eigenvalues <- as.numeric(result$singular_values)^2 / divisor$value
  names(eigenvalues) <- component.names
  .c6spc_fit(
    method = "PMD sparse PCA", method.class = "pmd_sparse_pca_fit",
    eigenvalues = eigenvalues, loadings = loadings, scores = scores,
    center = processed$center, operator = retained.operator,
    rank = components, n = n, p = p, variable.names = colnames(x),
    diagnostics = c(base.diagnostics, list(
      certified = TRUE,
      sign.anchor.coordinates = setNames(anchored$anchors, component.names),
      eigenvalue.definition = "squared PMD residual singular value divided by the covariance divisor"
    )),
    call = call,
    extra = list(
      scale = processed$scale,
      left.vectors = left,
      singular.values = setNames(as.numeric(result$singular_values),
                                 component.names),
      residual.matrix = as.matrix(result$residual_matrix)
    )
  )
}


#' Penalized-matrix-decomposition sparse CCA
#'
#' Solves the PMD sparse CCA problem
#' \deqn{\max_{u,v}u'Mv,\quad \lVert u\rVert_2,\lVert v\rVert_2\leq1,
#' \quad \lVert u\rVert_1\leq c_x,\quad \lVert v\rVert_1\leq c_y}
#' with the cross-covariance operator applied matrix-free.  Later pairs use
#' `M <- M - d u v'`.  Both actual l1 bounds must be supplied; no permutation
#' tuning or default grid is generated.
#'
#' For the shared `hd_pca_fit` contract, `loadings`, `scores`, `center`, `p`,
#' and `variable.names` refer to the x view.  The y-view counterparts are
#' returned separately.
#'
#' @param x,y Numeric paired data matrices with observations in rows.
#' @param l1_x,l1_y Required actual l1 bounds, scalar or one per component.
#' @param components Number of canonical pairs.
#' @param center Logical or supplied centers. A length-two logical vector may
#'   control x and y separately.
#' @param scale Logical or supplied list `list(x=..., y=...)` of positive
#'   scales. Logical `TRUE` uses sample standard deviations.
#' @param covariance_divisor Common cross-covariance divisor.
#' @param initial_x,initial_y Optional deterministic coefficient starts; both
#'   must be supplied together.
#' @param solver_tol,solver_max_iter Alternating-solver controls.
#' @param strict If `TRUE`, failure is an error; otherwise return an invalid fit.
#' @param keep_operator Retain the explicit cross-covariance matrix. The solver
#'   itself remains matrix-free.
#'
#' @return A `pmd_sparse_cca_fit` inheriting from `hd_pca_fit`.
#' @export
#'
#' @references
#' Witten, D. M., Tibshirani, R., and Hastie, T. (2009). A penalized matrix
#' decomposition, with applications to sparse principal components and
#' canonical correlation analysis. *Biostatistics*, 10, 515--534.
#'
#' @examples
#' x <- cbind(a = c(-2, -1, 1, 2), b = c(1, -1, -1, 1))
#' y <- cbind(c = x[, 1], d = c(-1, 1, 1, -1))
#' pmd_sparse_cca(x, y, l1_x = 1, l1_y = 1, scale = FALSE)
pmd_sparse_cca <- function(
    x, y, l1_x, l1_y, components = 1L, center = TRUE, scale = TRUE,
    covariance_divisor = c("n", "n-1"), initial_x = NULL,
    initial_y = NULL, solver_tol = 1e-8, solver_max_iter = 1000L,
    strict = TRUE, keep_operator = TRUE) {
  call <- match.call()
  strict <- .c6spc_flag(strict, "strict")
  keep_operator <- .c6spc_flag(keep_operator, "keep_operator")
  x <- .c6spc_data(x, "x")
  y <- .c6spc_data(y, "y")
  if (nrow(x) != nrow(y)) {
    stop("`x` and `y` must have the same number of rows.", call. = FALSE)
  }
  if (!is.null(rownames(x)) && !is.null(rownames(y)) &&
      !identical(rownames(x), rownames(y))) {
    stop("Row names of paired `x` and `y` must agree exactly.",
         call. = FALSE)
  }
  n <- nrow(x)
  px <- ncol(x)
  py <- ncol(y)
  components <- .c6spc_count(
    components, "components", maximum = min(n, px, py)
  )
  l1.x <- .c6spc_l1_bounds(l1_x, components, px, "l1_x")
  l1.y <- .c6spc_l1_bounds(l1_y, components, py, "l1_y")
  if (xor(is.null(initial_x), is.null(initial_y))) {
    stop("`initial_x` and `initial_y` must be supplied together.",
         call. = FALSE)
  }
  initial.x <- .c6spc_initial(initial_x, px, components, "initial_x")
  initial.y <- .c6spc_initial(initial_y, py, components, "initial_y")
  solver_tol <- .c6spc_scalar(solver_tol, "solver_tol")
  solver_max_iter <- .c6spc_count(solver_max_iter, "solver_max_iter")

  if (is.logical(center)) {
    if (!length(center) %in% c(1L, 2L) || anyNA(center)) {
      stop("Logical `center` must have length one or two.", call. = FALSE)
    }
    if (length(center) == 1L) center <- rep.int(center, 2L)
    center.x <- if (center[1L]) "mean" else "none"
    center.y <- if (center[2L]) "mean" else "none"
  } else if (is.list(center) && setequal(names(center), c("x", "y"))) {
    center.x <- center$x
    center.y <- center$y
  } else {
    stop("`center` must be logical or a named list with x and y entries.",
         call. = FALSE)
  }
  if (is.logical(scale)) {
    if (!length(scale) %in% c(1L, 2L) || anyNA(scale)) {
      stop("Logical `scale` must have length one or two.", call. = FALSE)
    }
    if (length(scale) == 1L) scale <- rep.int(scale, 2L)
    scale.x <- scale[1L]
    scale.y <- scale[2L]
  } else if (is.list(scale) && setequal(names(scale), c("x", "y"))) {
    scale.x <- scale$x
    scale.y <- scale$y
  } else {
    stop("`scale` must be logical or a named list with x and y entries.",
         call. = FALSE)
  }
  processed.x <- .c6spc_center_scale(x, center.x, scale.x, "x")
  processed.y <- .c6spc_center_scale(y, center.y, scale.y, "y")
  divisor <- .c6spc_divisor(covariance_divisor, n)
  root.divisor <- sqrt(divisor$value)
  result <- cpp_ch6spc_pmd_cca(
    processed.x$x / root.divisor, processed.y$x / root.divisor,
    l1.x, l1.y, initial.x, initial.y, solver_tol, solver_max_iter
  )
  operator <- if (keep_operator) {
    crossprod(processed.x$x, processed.y$x) / divisor$value
  } else NULL
  if (!is.null(operator)) {
    dimnames(operator) <- list(colnames(x), colnames(y))
  }
  base.diagnostics <- list(
    tuning.source = "explicit l1_x and l1_y",
    l1.x = l1.x, l1.y = l1.y,
    solver = result$certificates,
    solver.tolerance = solver_tol,
    solver.maximum.iterations = solver_max_iter,
    cross.operator.evaluation = paste(
      "matrix-free X'(Yv) and Y'(Xu), plus explicit rank-one deflations"
    ),
    deflation = "M <- M - d u v'",
    x.centering = processed.x$center.source,
    y.centering = processed.y$center.source,
    x.scaling = processed.x$scale.source,
    y.scaling = processed.y$scale.source,
    covariance.divisor = divisor$label,
    covariance.divisor.value = divisor$value
  )
  if (!.c6spc_solver_ok(result, components)) {
    message <- sprintf(
      "PMD sparse-CCA pair %d did not pass its block certificate: %s",
      as.integer(result$failure_component), as.character(result$failure_message)
    )
    return(.c6spc_fail(
      strict, message,
      method = "PMD sparse CCA", method.class = "pmd_sparse_cca_fit",
      n = n, p = px, variable.names = colnames(x),
      center = processed.x$center, operator = operator,
      requested.rank = components, stage = "PMD CCA alternating solver",
      call = call, diagnostics = base.diagnostics,
      extra = list(
        scale = processed.x$scale,
        y.center = processed.y$center, y.scale = processed.y$scale,
        y.p = py, y.variable.names = colnames(y),
        y.loadings = matrix(numeric(), py, 0L),
        y.scores = matrix(numeric(), n, 0L)
      )
    ))
  }
  x.loadings <- as.matrix(result$x_loadings)
  y.loadings <- as.matrix(result$y_loadings)
  anchored <- .c6spc_anchor(x.loadings, paired = y.loadings)
  x.loadings <- anchored$loadings
  y.loadings <- anchored$paired
  component.names <- paste0("CC", seq_len(components))
  dimnames(x.loadings) <- list(colnames(x), component.names)
  dimnames(y.loadings) <- list(colnames(y), component.names)
  x.scores <- processed.x$x %*% x.loadings
  y.scores <- processed.y$x %*% y.loadings
  dimnames(x.scores) <- list(rownames(x), component.names)
  dimnames(y.scores) <- list(rownames(y), component.names)
  strengths <- setNames(as.numeric(result$singular_values), component.names)
  correlations <- vapply(seq_len(components), function(j) {
    sx <- x.scores[, j]
    sy <- y.scores[, j]
    if (sum(sx^2) == 0 || sum(sy^2) == 0) NA_real_ else
      stats::cor(sx, sy)
  }, numeric(1))
  names(correlations) <- component.names
  .c6spc_fit(
    method = "PMD sparse CCA", method.class = "pmd_sparse_cca_fit",
    eigenvalues = strengths, loadings = x.loadings, scores = x.scores,
    center = processed.x$center, operator = operator,
    rank = components, n = n, p = px, variable.names = colnames(x),
    diagnostics = c(base.diagnostics, list(
      certified = TRUE,
      sign.anchor.coordinates = setNames(anchored$anchors, component.names),
      eigenvalue.definition = "PMD deflated cross-covariance strength",
      correlation.definition = "ordinary correlation of the paired sparse scores"
    )),
    call = call,
    extra = list(
      scale = processed.x$scale,
      y.center = processed.y$center, y.scale = processed.y$scale,
      y.p = py, y.variable.names = colnames(y),
      y.loadings = y.loadings, y.scores = y.scores,
      canonical.correlations = correlations,
      singular.values = strengths
    )
  )
}
