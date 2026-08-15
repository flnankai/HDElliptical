# Chapter 6 spatial-sign sparse canonical correlation analysis.

.sscca_positive <- function(value, name, allow_zero = FALSE) {
  value <- as.numeric(value)
  bad <- length(value) != 1L || is.na(value) || !is.finite(value) ||
    if (allow_zero) value < 0 else value <= 0
  if (bad) {
    qualifier <- if (allow_zero) "non-negative" else "strictly positive"
    stop(sprintf("`%s` must be one finite %s number.", name, qualifier),
         call. = FALSE)
  }
  value
}


.sscca_count <- function(value, name) {
  if (!is.numeric(value) || length(value) != 1L || is.na(value) ||
      !is.finite(value) || value < 1 || value != floor(value) ||
      value > .Machine$integer.max) {
    stop(sprintf("`%s` must be a positive integer.", name), call. = FALSE)
  }
  as.integer(value)
}


.sscca_strict <- function(strict) {
  if (!is.logical(strict) || length(strict) != 1L || is.na(strict)) {
    stop("`strict` must be TRUE or FALSE.", call. = FALSE)
  }
  strict
}


.sscca_pair <- function(x, y, min_rows = 3L) {
  x <- .as_data_matrix(x, "x", min_rows = min_rows)
  y <- .as_data_matrix(y, "y", min_rows = min_rows)
  if (nrow(x) != nrow(y)) {
    stop("`x` and `y` must have the same number of paired rows.",
         call. = FALSE)
  }
  if (!is.null(rownames(x)) && !is.null(rownames(y)) &&
      !identical(rownames(x), rownames(y))) {
    stop("Non-null row names of paired `x` and `y` must be identical.",
         call. = FALSE)
  }
  for (entry in list(x = x, y = y)) {
    feature.names <- colnames(entry)
    if (!is.null(feature.names) &&
        (anyNA(feature.names) || any(!nzchar(feature.names)) ||
         anyDuplicated(feature.names))) {
      stop("Column names in each data block must be non-empty and unique.",
           call. = FALSE)
    }
  }
  list(x = x, y = y)
}


.sscca_lambda <- function(value, name, selection) {
  value <- as.numeric(value)
  if (!length(value) || anyNA(value) || any(!is.finite(value)) ||
      any(value < 0)) {
    stop(sprintf("`%s` must contain finite non-negative values.", name),
         call. = FALSE)
  }
  if (selection == "fixed" && length(value) != 1L) {
    stop(sprintf("`%s` must be scalar when `selection = \"fixed\"`.", name),
         call. = FALSE)
  }
  sort(unique(value), decreasing = TRUE)
}


.sscca_metric_diagnostics <- function(metric, tolerance, name) {
  metric <- (metric + t(metric)) / 2
  values <- eigen(metric, symmetric = TRUE, only.values = TRUE)$values
  scale <- max(1, max(abs(values)))
  if (any(values < -tolerance * scale)) {
    stop(sprintf(
      "The `%s` spatial-sign metric is materially indefinite; no PSD repair is applied.",
      name
    ), call. = FALSE)
  }
  if (any(diag(metric) <= 0)) {
    stop(sprintf(
      "Every diagonal entry of the `%s` spatial-sign metric must be positive.",
      name
    ), call. = FALSE)
  }
  positive <- values > tolerance * scale
  list(
    eigenvalues = values,
    rank = sum(positive),
    reciprocal.condition = if (all(positive)) min(values) / max(values) else 0,
    minimum.eigenvalue = min(values),
    tolerance = tolerance,
    no.repair = TRUE
  )
}


.sscca_inverse_sqrt <- function(metric, ridge, tolerance, name) {
  ridge <- .sscca_positive(ridge, paste0("ridge_", name), allow_zero = TRUE)
  regularized <- (metric + t(metric)) / 2 + diag(ridge, nrow(metric))
  decomposition <- eigen(regularized, symmetric = TRUE)
  values <- as.numeric(decomposition$values)
  scale <- max(1, max(abs(values)))
  if (any(values <= tolerance * scale)) {
    stop(sprintf(
      "The `%s` sign block plus its explicit ridge is not full-rank SPD at `rank_tol`; no pseudoinverse or eigenvalue floor is used.",
      name
    ), call. = FALSE)
  }
  root <- decomposition$vectors %*%
    (diag(1 / sqrt(values), length(values))) %*%
    t(decomposition$vectors)
  list(
    inverse.sqrt = (root + t(root)) / 2,
    regularized = regularized,
    eigenvalues = values,
    reciprocal.condition = min(values) / max(values),
    ridge = ridge
  )
}


.sscca_invalid <- function(method, class_name, message, strict, call,
                            data.name, diagnostics = list()) {
  if (strict) stop(message, call. = FALSE)
  warning(message, call. = FALSE)
  structure(list(
    valid = FALSE,
    method = method,
    estimate = NULL,
    canonical.correlations = NA_real_,
    diagnostics = c(list(
      failure = message,
      failure.stage = "method fitting",
      no.repair = TRUE
    ), diagnostics),
    call = call,
    data.name = data.name
  ), class = c(class_name, "hd_cca_fit", "list"))
}


.sscca_joint_sign_blocks <- function(x, y, median_tol, median_max_iter,
                                      zero_tol, zero_action) {
  joint <- cbind(x, y)
  center <- spatial_median(
    joint, tol = median_tol, max_iter = median_max_iter,
    zero_tol = zero_tol, warn = FALSE
  )
  median.diagnostics <- list(
    converged = isTRUE(attr(center, "converged")),
    iterations = attr(center, "iterations"),
    relative.change = attr(center, "relative_change"),
    equation.residual = attr(center, "equation_residual"),
    objective = attr(center, "objective")
  )
  if (!median.diagnostics$converged) {
    return(list(valid = FALSE,
                message = "The joint spatial median did not converge.",
                median.diagnostics = median.diagnostics))
  }
  signs <- spatial_sign(joint, center = center, zero_tol = zero_tol)
  n.zero <- as.integer(attr(signs, "n_zero"))
  if (zero_action == "error" && n.zero > 0L) {
    return(list(
      valid = FALSE,
      message = paste(
        "At least one observation equals the fitted joint spatial median;",
        "choose `zero_action = \"zero\"` to use the book's U(0)=0 convention."
      ),
      median.diagnostics = median.diagnostics
    ))
  }
  px <- ncol(x)
  py <- ncol(y)
  n <- nrow(x)
  covariance <- crossprod(signs) / n
  xx <- covariance[seq_len(px), seq_len(px), drop = FALSE]
  yy.index <- px + seq_len(py)
  yy <- covariance[yy.index, yy.index, drop = FALSE]
  xy <- covariance[seq_len(px), yy.index, drop = FALSE]
  dimnames(xx) <- list(colnames(x), colnames(x))
  dimnames(yy) <- list(colnames(y), colnames(y))
  dimnames(xy) <- list(colnames(x), colnames(y))
  center.x <- center[seq_len(px)]
  center.y <- center[yy.index]
  names(center.x) <- colnames(x)
  names(center.y) <- colnames(y)
  list(
    valid = TRUE,
    signs = signs,
    covariance = list(xx = xx, xy = xy, yx = t(xy), yy = yy),
    center = list(x = center.x, y = center.y, joint = center),
    n.zero = n.zero,
    norms = attr(signs, "norms"),
    median.diagnostics = median.diagnostics
  )
}


.sscca_candidate_table <- function(lambda, bic, rss, kkt, iterations,
                                    converged) {
  data.frame(
    lambda = as.numeric(lambda),
    bic = as.numeric(bic),
    residual.sum.of.squares = as.numeric(rss),
    kkt.residual = as.numeric(kkt),
    iterations = as.integer(iterations),
    converged = as.logical(converged),
    row.names = NULL,
    check.names = FALSE
  )
}


#' Primary spatial-sign sparse canonical correlation analysis
#'
#' Fits the sparse CCA estimator proposed by Qian, Liu, and Feng.  For paired
#' rows, a joint spatial median and joint sample spatial-sign covariance matrix
#' are computed, then the first sparse pair solves
#' \deqn{\max_{w_x,w_y}\;w_x^T(pS_{xy})w_y-
#' \lambda_x\|w_x\|_1-\lambda_y\|w_y\|_1}
#' subject to \eqn{w_x^T(pS_{xx})w_x\leq1} and
#' \eqn{w_y^T(pS_{yy})w_y\leq1}, where
#' \eqn{p=p_x+p_y}.
#'
#' Every conditional metric-lasso block is solved with its actual diagonal.
#' This corrects the referenced `mixedCCA` coordinate update, which omits
#' division by the metric diagonal and is exact only when that diagonal is one.
#' With `selection = "bic1"` or `"bic2"`, each alternating block chooses from
#' the explicitly supplied lambda sequence using the two criteria printed in
#' the primary paper.  No lambda grid, ridge, pseudoinverse, or matrix repair is
#' invented by this function.
#'
#' @param x,y Paired finite numeric matrices, observations in rows.
#' @param lambda_x,lambda_y Non-negative penalties.  Scalars are required for
#'   fixed fitting; explicit sequences are allowed for BIC selection.
#' @param selection One of `"fixed"`, `"bic1"`, or `"bic2"`.
#' @param tol Positive relative tolerance for the alternating directions.
#' @param max_iter Maximum alternating iterations.
#' @param inner_tol Positive KKT/update tolerance for each metric-lasso block.
#' @param inner_max_iter Maximum coordinate-descent sweeps per candidate.
#' @param median_tol,median_max_iter Controls passed to the joint spatial
#'   median fit.
#' @param zero_tol Non-negative exact/near-zero spatial-sign tolerance.
#' @param zero_action Use the package convention `"zero"` for \eqn{U(0)=0}, or
#'   fail with `"error"` when a fitted residual is zero.
#' @param support_tol Non-negative absolute threshold used only to count BIC
#'   degrees of freedom and report support.  It does not alter coefficients.
#' @param metric_tol Positive relative tolerance for PSD/rank diagnostics; no
#'   eigenvalue is floored or projected.
#' @param strict If `TRUE`, method failures are errors.  Otherwise an invalid
#'   `hd_cca_fit` is returned and no last iterate is presented as an estimate.
#'
#' @return An `sscca_fit` and `hd_cca_fit` object containing the two canonical
#'   coefficient vectors, raw and sign scores, spatial-sign blocks, selected
#'   penalties, BIC paths, convergence and KKT certificates.
#' @export
#'
#' @references
#' Qian, J., Liu, W., and Feng, L. (2025). High dimensional sparse canonical
#' correlation analysis for elliptical symmetric distributions.
#' arXiv:2504.13018.
#'
#' @examples
#' t <- seq_len(12)
#' x <- cbind(x1 = sin(t), x2 = cos(t / 2))
#' y <- cbind(y1 = sin(t) + 0.2 * cos(t), y2 = cos(t / 2) - 0.1 * sin(t))
#' fit <- sscca(x, y, lambda_x = 0.01, lambda_y = 0.01)
#' fit$canonical.correlations
sscca <- function(
    x, y, lambda_x, lambda_y,
    selection = c("fixed", "bic1", "bic2"),
    tol = 1e-7, max_iter = 500L,
    inner_tol = 1e-9, inner_max_iter = 5000L,
    median_tol = 1e-8, median_max_iter = 500L,
    zero_tol = 0, zero_action = c("zero", "error"),
    support_tol = 0, metric_tol = sqrt(.Machine$double.eps),
    strict = TRUE) {
  call <- match.call()
  data.name <- paste(deparse1(substitute(x)), "and", deparse1(substitute(y)))
  pair <- .sscca_pair(x, y)
  x <- pair$x
  y <- pair$y
  selection <- match.arg(selection)
  zero_action <- match.arg(zero_action)
  strict <- .sscca_strict(strict)
  tol <- .sscca_positive(tol, "tol")
  max_iter <- .sscca_count(max_iter, "max_iter")
  inner_tol <- .sscca_positive(inner_tol, "inner_tol")
  inner_max_iter <- .sscca_count(inner_max_iter, "inner_max_iter")
  median_tol <- .sscca_positive(median_tol, "median_tol")
  median_max_iter <- .sscca_count(median_max_iter, "median_max_iter")
  zero_tol <- .sscca_positive(zero_tol, "zero_tol", allow_zero = TRUE)
  support_tol <- .sscca_positive(
    support_tol, "support_tol", allow_zero = TRUE
  )
  metric_tol <- .sscca_positive(metric_tol, "metric_tol")
  lambda_x <- .sscca_lambda(lambda_x, "lambda_x", selection)
  lambda_y <- .sscca_lambda(lambda_y, "lambda_y", selection)

  blocks <- .sscca_joint_sign_blocks(
    x, y, median_tol, median_max_iter, zero_tol, zero_action
  )
  if (!isTRUE(blocks$valid)) {
    return(.sscca_invalid(
      "Primary spatial-sign sparse CCA", "sscca_fit", blocks$message,
      strict, call, data.name,
      diagnostics = list(spatial.median = blocks$median.diagnostics)
    ))
  }
  dimension <- ncol(x) + ncol(y)
  metric.x <- dimension * blocks$covariance$xx
  metric.y <- dimension * blocks$covariance$yy
  cross <- dimension * blocks$covariance$xy
  metric.x.diagnostics <- tryCatch(
    .sscca_metric_diagnostics(metric.x, metric_tol, "x"),
    error = identity
  )
  metric.y.diagnostics <- tryCatch(
    .sscca_metric_diagnostics(metric.y, metric_tol, "y"),
    error = identity
  )
  if (inherits(metric.x.diagnostics, "error") ||
      inherits(metric.y.diagnostics, "error")) {
    message <- if (inherits(metric.x.diagnostics, "error")) {
      conditionMessage(metric.x.diagnostics)
    } else {
      conditionMessage(metric.y.diagnostics)
    }
    return(.sscca_invalid(
      "Primary spatial-sign sparse CCA", "sscca_fit", message,
      strict, call, data.name
    ))
  }

  native <- tryCatch(
    cpp_ch6_sscca_primary(
      metric.x, metric.y, cross, lambda_x, lambda_y,
      match(selection, c("fixed", "bic1", "bic2")) - 1L,
      nrow(x), tol, max_iter, inner_tol, inner_max_iter, support_tol
    ),
    error = identity
  )
  if (inherits(native, "error") || !isTRUE(native$valid)) {
    message <- if (inherits(native, "error")) conditionMessage(native) else
      native$message
    return(.sscca_invalid(
      "Primary spatial-sign sparse CCA", "sscca_fit", message,
      strict, call, data.name,
      diagnostics = list(native = if (inherits(native, "error")) NULL else native)
    ))
  }
  certificate.limit <- 20 * (tol + inner_tol) *
    (1 + max(abs(metric.x), abs(metric.y), abs(cross)))
  certified <- isTRUE(native$converged) &&
    isTRUE(native$inner_converged_x) && isTRUE(native$inner_converged_y) &&
    native$pair_kkt_x <= certificate.limit &&
    native$pair_kkt_y <= certificate.limit
  if (!certified) {
    return(.sscca_invalid(
      "Primary spatial-sign sparse CCA", "sscca_fit",
      paste(
        "The alternating SSCCA fit did not meet its update/KKT certificate;",
        "no uncertified last iterate is returned as an estimate."
      ), strict, call, data.name,
      diagnostics = list(native = native,
                         certificate.limit = certificate.limit)
    ))
  }

  x.coefficients <- as.numeric(native$x)
  y.coefficients <- as.numeric(native$y)
  names(x.coefficients) <- colnames(x)
  names(y.coefficients) <- colnames(y)
  centered.x <- sweep(x, 2L, blocks$center$x, check.margin = FALSE)
  centered.y <- sweep(y, 2L, blocks$center$y, check.margin = FALSE)
  x.scores <- drop(centered.x %*% x.coefficients)
  y.scores <- drop(centered.y %*% y.coefficients)
  names(x.scores) <- rownames(x)
  names(y.scores) <- rownames(y)
  signs.x <- blocks$signs[, seq_len(ncol(x)), drop = FALSE]
  signs.y <- blocks$signs[, ncol(x) + seq_len(ncol(y)), drop = FALSE]
  x.sign.scores <- drop(signs.x %*% x.coefficients)
  y.sign.scores <- drop(signs.y %*% y.coefficients)
  names(x.sign.scores) <- rownames(x)
  names(y.sign.scores) <- rownames(y)
  association <- as.numeric(native$association)
  canonical <- stats::setNames(association, "CC1")
  selected.x <- as.integer(native$selected_index_x)
  selected.y <- as.integer(native$selected_index_y)
  candidate.x <- .sscca_candidate_table(
    lambda_x, native$bic_values_x, native$rss_values_x,
    native$candidate_kkt_x, native$candidate_iterations_x,
    native$candidate_converged_x
  )
  candidate.y <- .sscca_candidate_table(
    lambda_y, native$bic_values_y, native$rss_values_y,
    native$candidate_kkt_y, native$candidate_iterations_y,
    native$candidate_converged_y
  )

  structure(list(
    valid = TRUE,
    method = "Primary spatial-sign sparse canonical correlation analysis",
    canonical.correlations = canonical,
    x.coefficients = x.coefficients,
    y.coefficients = y.coefficients,
    x.scores = x.scores,
    y.scores = y.scores,
    x.sign.scores = x.sign.scores,
    y.sign.scores = y.sign.scores,
    center = blocks$center,
    spatial.sign.covariance = blocks$covariance,
    operator = list(xx = metric.x, xy = cross, yx = t(cross), yy = metric.y),
    rank = 1L,
    n = nrow(x),
    p.x = ncol(x),
    p.y = ncol(y),
    variable.names = list(x = colnames(x), y = colnames(y)),
    observation.names = if (!is.null(rownames(x))) rownames(x) else rownames(y),
    tuning = list(
      selection = selection,
      lambda.x = as.numeric(native$lambda_x),
      lambda.y = as.numeric(native$lambda_y),
      candidate.x = candidate.x,
      candidate.y = candidate.y,
      selected.index.x = selected.x,
      selected.index.y = selected.y,
      support.tolerance = support_tol
    ),
    diagnostics = list(
      spatial.median = blocks$median.diagnostics,
      zero.action = zero_action,
      zero.signs = blocks$n.zero,
      outside.continuous.theory = blocks$n.zero > 0L,
      metric.x = metric.x.diagnostics,
      metric.y = metric.y.diagnostics,
      iterations = as.integer(native$iterations),
      converged = isTRUE(native$converged),
      relative.update = native$relative_update,
      objective = native$objective,
      objective.trace = native$objective_trace,
      update.trace = native$update_trace,
      lambda.x.trace = native$lambda_x_trace,
      lambda.y.trace = native$lambda_y_trace,
      metric.norm.x = native$metric_norm_x,
      metric.norm.y = native$metric_norm_y,
      inner.kkt.x = native$inner_kkt_x,
      inner.kkt.y = native$inner_kkt_y,
      pair.kkt.x = native$pair_kkt_x,
      pair.kkt.y = native$pair_kkt_y,
      certificate.limit = certificate.limit,
      degrees.freedom.x = as.integer(native$degrees_freedom_x),
      degrees.freedom.y = as.integer(native$degrees_freedom_y),
      raw.coefficients = list(x = native$raw_x, y = native$raw_y),
      raw.metric.scale = c(x = native$raw_metric_scale_x,
                           y = native$raw_metric_scale_y),
      sign.anchor.x = as.integer(native$sign_anchor_x),
      corrected.coordinate.update = paste(
        "each coordinate is divided by its actual spatial-sign metric diagonal;",
        "the referenced mixedCCA update omits this divisor"
      ),
      numerical.repair = "none; no ridge, pseudoinverse, floor, or projection"
    ),
    call = call,
    data.name = data.name
  ), class = c("sscca_fit", "hd_cca_fit", "list"))
}


#' Book sign-whitened sparse CCA variant
#'
#' Implements the constrained program printed in Chapter 6 of the book.  It
#' first whitens the two diagonal blocks of the joint sample spatial-sign
#' covariance and then solves
#' \deqn{\max_{u,v} u^T K_Sv,\quad
#' \|u\|_2\leq1,\ \|v\|_2\leq1,\
#' \|u\|_1\leq c_x,\ \|v\|_1\leq c_y.}
#'
#' This is intentionally a separately named book variant: it is not the
#' penalized metric program proposed in Qian, Liu, and Feng (2025), which is
#' implemented by [sscca()].  Any ridge used to make a whitening block SPD is
#' explicit and changes the reported operator; the defaults apply no ridge.
#'
#' @param x,y Paired finite numeric matrices, observations in rows.
#' @param c_x,c_y Explicit l1 bounds in `[1, sqrt(ncol(block))]`.
#' @param ridge_x,ridge_y Explicit non-negative diagonal additions used only
#'   in the two whitening metrics.
#' @param tol,max_iter Alternating PMD update controls.
#' @param median_tol,median_max_iter Joint spatial-median controls.
#' @param zero_tol,zero_action Spatial-sign zero convention; see [sscca()].
#' @param rank_tol Positive relative SPD tolerance for whitening.  No
#'   pseudoinverse or eigenvalue floor is used.
#' @param strict If `TRUE`, method failures are errors; otherwise an invalid
#'   fit is returned without exposing an uncertified estimate.
#'
#' @return A `sign_whitened_sparse_cca_fit` and `hd_cca_fit` object containing
#'   whitened and original-coordinate directions, scores, constraints and
#'   convergence diagnostics.
#'
#' @references
#' Feng, L. (2026). *High-Dimensional Data Analysis for Elliptical Symmetric
#' Distributions*, Chapter 6 (book manuscript).
#'
#' For the distinct primary metric-penalized method, see Qian, J., Liu, W., and
#' Feng, L. (2025). High dimensional sparse canonical correlation analysis for
#' elliptical symmetric distributions. \url{https://arxiv.org/abs/2504.13018}.
#' @export
#'
#' @examples
#' t <- seq_len(12)
#' x <- cbind(x1 = sin(t), x2 = cos(t / 2))
#' y <- cbind(y1 = sin(t) + 0.2 * cos(t), y2 = cos(t / 2) - 0.1 * sin(t))
#' fit <- sign_whitened_sparse_cca(x, y, c_x = sqrt(2), c_y = sqrt(2))
#' fit$canonical.correlations
sign_whitened_sparse_cca <- function(
    x, y, c_x, c_y, ridge_x = 0, ridge_y = 0,
    tol = 1e-7, max_iter = 500L,
    median_tol = 1e-8, median_max_iter = 500L,
    zero_tol = 0, zero_action = c("zero", "error"),
    rank_tol = sqrt(.Machine$double.eps), strict = TRUE) {
  call <- match.call()
  data.name <- paste(deparse1(substitute(x)), "and", deparse1(substitute(y)))
  pair <- .sscca_pair(x, y)
  x <- pair$x
  y <- pair$y
  zero_action <- match.arg(zero_action)
  strict <- .sscca_strict(strict)
  c_x <- .sscca_positive(c_x, "c_x")
  c_y <- .sscca_positive(c_y, "c_y")
  if (c_x < 1 || c_x > sqrt(ncol(x)) ||
      c_y < 1 || c_y > sqrt(ncol(y))) {
    stop("`c_x` and `c_y` must lie in [1, sqrt(number of block variables)].",
         call. = FALSE)
  }
  ridge_x <- .sscca_positive(ridge_x, "ridge_x", allow_zero = TRUE)
  ridge_y <- .sscca_positive(ridge_y, "ridge_y", allow_zero = TRUE)
  tol <- .sscca_positive(tol, "tol")
  max_iter <- .sscca_count(max_iter, "max_iter")
  median_tol <- .sscca_positive(median_tol, "median_tol")
  median_max_iter <- .sscca_count(median_max_iter, "median_max_iter")
  zero_tol <- .sscca_positive(zero_tol, "zero_tol", allow_zero = TRUE)
  rank_tol <- .sscca_positive(rank_tol, "rank_tol")

  blocks <- .sscca_joint_sign_blocks(
    x, y, median_tol, median_max_iter, zero_tol, zero_action
  )
  if (!isTRUE(blocks$valid)) {
    return(.sscca_invalid(
      "Book sign-whitened sparse CCA variant",
      "sign_whitened_sparse_cca_fit", blocks$message, strict, call,
      data.name, diagnostics = list(spatial.median = blocks$median.diagnostics)
    ))
  }
  root.x <- tryCatch(
    .sscca_inverse_sqrt(blocks$covariance$xx, ridge_x, rank_tol, "x"),
    error = identity
  )
  root.y <- tryCatch(
    .sscca_inverse_sqrt(blocks$covariance$yy, ridge_y, rank_tol, "y"),
    error = identity
  )
  if (inherits(root.x, "error") || inherits(root.y, "error")) {
    message <- if (inherits(root.x, "error")) conditionMessage(root.x) else
      conditionMessage(root.y)
    return(.sscca_invalid(
      "Book sign-whitened sparse CCA variant",
      "sign_whitened_sparse_cca_fit", message, strict, call, data.name
    ))
  }
  operator <- root.x$inverse.sqrt %*% blocks$covariance$xy %*%
    root.y$inverse.sqrt
  native <- tryCatch(
    cpp_ch6_sign_whitened_pmd(operator, c_x, c_y, tol, max_iter),
    error = identity
  )
  if (inherits(native, "error") || !isTRUE(native$valid)) {
    message <- if (inherits(native, "error")) conditionMessage(native) else
      native$message
    return(.sscca_invalid(
      "Book sign-whitened sparse CCA variant",
      "sign_whitened_sparse_cca_fit", message, strict, call, data.name,
      diagnostics = list(native = if (inherits(native, "error")) NULL else native)
    ))
  }
  certified <- isTRUE(native$converged) &&
    native$constraint_violation_x <= 20 * tol &&
    native$constraint_violation_y <= 20 * tol
  if (!certified) {
    return(.sscca_invalid(
      "Book sign-whitened sparse CCA variant",
      "sign_whitened_sparse_cca_fit",
      paste(
        "The sign-whitened PMD fit did not meet its update/constraint",
        "certificate; no uncertified last iterate is returned."
      ), strict, call, data.name,
      diagnostics = list(native = native)
    ))
  }

  whitened.x <- as.numeric(native$x)
  whitened.y <- as.numeric(native$y)
  coefficients.x <- drop(root.x$inverse.sqrt %*% whitened.x)
  coefficients.y <- drop(root.y$inverse.sqrt %*% whitened.y)
  names(whitened.x) <- names(coefficients.x) <- colnames(x)
  names(whitened.y) <- names(coefficients.y) <- colnames(y)
  centered.x <- sweep(x, 2L, blocks$center$x, check.margin = FALSE)
  centered.y <- sweep(y, 2L, blocks$center$y, check.margin = FALSE)
  x.scores <- drop(centered.x %*% coefficients.x)
  y.scores <- drop(centered.y %*% coefficients.y)
  names(x.scores) <- rownames(x)
  names(y.scores) <- rownames(y)
  signs.x <- blocks$signs[, seq_len(ncol(x)), drop = FALSE]
  signs.y <- blocks$signs[, ncol(x) + seq_len(ncol(y)), drop = FALSE]
  association <- as.numeric(native$association)

  structure(list(
    valid = TRUE,
    method = "Book sign-whitened spatial-sign sparse CCA variant",
    canonical.correlations = stats::setNames(association, "CC1"),
    x.coefficients = coefficients.x,
    y.coefficients = coefficients.y,
    x.whitened.directions = whitened.x,
    y.whitened.directions = whitened.y,
    x.scores = x.scores,
    y.scores = y.scores,
    x.sign.scores = drop(signs.x %*% coefficients.x),
    y.sign.scores = drop(signs.y %*% coefficients.y),
    center = blocks$center,
    spatial.sign.covariance = blocks$covariance,
    operator = operator,
    rank = 1L,
    n = nrow(x),
    p.x = ncol(x),
    p.y = ncol(y),
    variable.names = list(x = colnames(x), y = colnames(y)),
    observation.names = if (!is.null(rownames(x))) rownames(x) else rownames(y),
    tuning = list(
      c.x = c_x, c.y = c_y, ridge.x = ridge_x, ridge.y = ridge_y
    ),
    diagnostics = list(
      variant = paste(
        "book whitened l1-constrained program; not the penalized metric",
        "SSCCA estimator in Qian, Liu, and Feng (2025)"
      ),
      spatial.median = blocks$median.diagnostics,
      zero.action = zero_action,
      zero.signs = blocks$n.zero,
      outside.continuous.theory = blocks$n.zero > 0L,
      reciprocal.condition.x = root.x$reciprocal.condition,
      reciprocal.condition.y = root.y$reciprocal.condition,
      iterations = as.integer(native$iterations),
      converged = isTRUE(native$converged),
      relative.update = native$relative_update,
      objective.trace = native$objective_trace,
      l1.x = native$l1_x,
      l1.y = native$l1_y,
      l2.x = native$l2_x,
      l2.y = native$l2_y,
      threshold.x = native$threshold_x,
      threshold.y = native$threshold_y,
      top.tie.slack.x = native$top_tie_slack_x,
      top.tie.slack.y = native$top_tie_slack_y,
      constraint.violation.x = native$constraint_violation_x,
      constraint.violation.y = native$constraint_violation_y,
      sign.anchor.x = as.integer(native$sign_anchor_x),
      numerical.repair = paste(
        "none; only the user-supplied ridge is included in the reported",
        "whitening operator"
      )
    ),
    call = call,
    data.name = data.name
  ), class = c("sign_whitened_sparse_cca_fit", "hd_cca_fit", "list"))
}
