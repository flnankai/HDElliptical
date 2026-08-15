# Chapter 6: classical PCA/CCA and robust factor analysis.

.ch6cf_data <- function(x, name = "x", min_rows = 2L) {
  x <- .as_data_matrix(x, name = name, min_rows = min_rows)
  feature.names <- colnames(x)
  if (!is.null(feature.names) &&
      (anyNA(feature.names) || any(!nzchar(feature.names)) ||
       anyDuplicated(feature.names))) {
    stop(sprintf("Column names of `%s` must be non-empty and unique.", name),
         call. = FALSE)
  }
  x
}


.ch6cf_count <- function(x, name, minimum = 1L, maximum = Inf) {
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


.ch6cf_positive <- function(x, name, allow_zero = FALSE) {
  x <- as.numeric(x)
  bad <- length(x) != 1L || is.na(x) || !is.finite(x) ||
    if (allow_zero) x < 0 else x <= 0
  if (bad) {
    qualifier <- if (allow_zero) "non-negative" else "strictly positive"
    stop(sprintf("`%s` must be one finite %s number.", name, qualifier),
         call. = FALSE)
  }
  x
}


.ch6cf_center <- function(x, center, choices = c("mean", "none"),
                          name = "center") {
  p <- ncol(x)
  feature.names <- colnames(x)
  if (is.character(center)) {
    source <- match.arg(center, choices)
    value <- switch(source, mean = colMeans(x), none = numeric(p))
  } else {
    supplied.names <- names(center)
    center <- as.numeric(center)
    if (length(center) != p || anyNA(center) || any(!is.finite(center))) {
      stop(sprintf("`%s` must contain one finite value per column.", name),
           call. = FALSE)
    }
    if (!is.null(supplied.names)) {
      if (anyNA(supplied.names) || any(!nzchar(supplied.names)) ||
          anyDuplicated(supplied.names)) {
        stop(sprintf("Names of `%s` must be non-empty and unique.", name),
             call. = FALSE)
      }
      if (is.null(feature.names) || !setequal(supplied.names, feature.names)) {
        stop(sprintf(
          "Named `%s` must match the data column names exactly.", name
        ), call. = FALSE)
      }
      center <- center[match(feature.names, supplied.names)]
    }
    value <- center
    source <- "supplied"
  }
  value <- as.numeric(value)
  names(value) <- feature.names
  list(value = value, source = source)
}


.ch6cf_divisor <- function(covariance_divisor, n) {
  covariance_divisor <- match.arg(covariance_divisor, c("n-1", "n"))
  divisor <- if (covariance_divisor == "n-1") n - 1 else n
  if (divisor <= 0) {
    stop("The selected covariance divisor is not positive.", call. = FALSE)
  }
  list(label = covariance_divisor, value = as.numeric(divisor))
}


.ch6cf_anchor_columns <- function(vectors, paired = NULL) {
  vectors <- as.matrix(vectors)
  if (!ncol(vectors)) return(list(primary = vectors, paired = paired,
                                   anchors = integer()))
  anchors <- integer(ncol(vectors))
  for (j in seq_len(ncol(vectors))) {
    magnitudes <- abs(vectors[, j])
    anchor <- which(magnitudes == max(magnitudes))[1L]
    anchors[j] <- anchor
    if (vectors[anchor, j] < 0) {
      vectors[, j] <- -vectors[, j]
      if (!is.null(paired)) paired[, j] <- -paired[, j]
    }
  }
  list(primary = vectors, paired = paired, anchors = anchors)
}


.ch6cf_tie_blocks <- function(values, tolerance) {
  if (!length(values)) return(list())
  scale <- max(1, max(abs(values)))
  breaks <- c(TRUE, abs(diff(values)) > tolerance * scale)
  split(seq_along(values), cumsum(breaks))
}


.ch6cf_psd_eigen <- function(operator, name, tolerance) {
  tolerance <- .ch6cf_positive(tolerance, "eigen_tol")
  operator <- as.matrix(operator)
  if (!is.numeric(operator) || nrow(operator) != ncol(operator) ||
      nrow(operator) < 1L || anyNA(operator) || any(!is.finite(operator))) {
    stop(sprintf("`%s` must be a finite non-empty square matrix.", name),
         call. = FALSE)
  }
  scale <- max(1, max(abs(operator)))
  symmetry.error <- max(abs(operator - t(operator)))
  if (symmetry.error > tolerance * scale) {
    stop(sprintf("`%s` is not symmetric within `eigen_tol`.", name),
         call. = FALSE)
  }
  operator <- (operator + t(operator)) / 2
  decomposition <- eigen(operator, symmetric = TRUE)
  values <- as.numeric(decomposition$values)
  eigen.scale <- max(1, max(abs(values)))
  if (any(values < -tolerance * eigen.scale)) {
    stop(sprintf(
      "`%s` has a materially negative eigenvalue; no PSD repair is applied.",
      name
    ), call. = FALSE)
  }
  rounded <- values < 0
  values[rounded] <- 0
  anchored <- .ch6cf_anchor_columns(decomposition$vectors)
  rank <- sum(values > tolerance * max(1, values[1L]))
  list(
    operator = operator,
    values = values,
    vectors = anchored$primary,
    anchors = anchored$anchors,
    rank = as.integer(rank),
    tie.blocks = .ch6cf_tie_blocks(values, tolerance),
    symmetry.error = symmetry.error,
    negative.roundoff.clipped = sum(rounded),
    tolerance = tolerance
  )
}


.ch6cf_spd_inverse_sqrt <- function(matrix, name, tolerance) {
  spectral <- .ch6cf_psd_eigen(matrix, name, tolerance)
  threshold <- tolerance * max(1, spectral$values[1L])
  if (spectral$rank != nrow(matrix) || any(spectral$values <= threshold)) {
    stop(sprintf(
      "`%s` must be strictly positive definite and full rank at `rank_tol`; no ridge or pseudoinverse is used.",
      name
    ), call. = FALSE)
  }
  inverse.sqrt <- spectral$vectors %*%
    (diag(1 / sqrt(spectral$values), nrow(matrix))) %*%
    t(spectral$vectors)
  list(
    inverse.sqrt = (inverse.sqrt + t(inverse.sqrt)) / 2,
    eigenvalues = spectral$values,
    reciprocal.condition = min(spectral$values) / max(spectral$values),
    tolerance = tolerance
  )
}


.ch6cf_component_count <- function(components, maximum, default, name) {
  if (is.null(components)) return(as.integer(default))
  .ch6cf_count(components, name, minimum = 1L, maximum = maximum)
}


.ch6cf_boundary_gap <- function(values, rank, tolerance, name,
                                require_gap = FALSE) {
  if (rank >= length(values)) {
    return(list(gap = Inf, scale = max(1, max(abs(values))),
                separated = TRUE))
  }
  scale <- max(1, max(abs(values)))
  gap <- values[rank] - values[rank + 1L]
  separated <- is.finite(gap) && gap > tolerance * scale
  if (require_gap && !separated) {
    stop(sprintf(
      "The requested `%s` boundary splits an unresolved repeated-eigenvalue block; use a separated boundary and compare repeated directions through their projector.",
      name
    ), call. = FALSE)
  }
  list(gap = gap, scale = scale, separated = separated)
}


.ch6cf_dimnames_square <- function(matrix, feature.names) {
  if (!is.null(feature.names)) dimnames(matrix) <- list(feature.names,
                                                        feature.names)
  matrix
}


#' Classical principal component analysis with explicit covariance convention
#'
#' Computes PCA from the centered second-moment matrix with a user-visible
#' divisor.  The book's displayed estimator uses `covariance_divisor = "n"`;
#' `"n-1"` gives the usual unbiased sample covariance.  No variance floor or
#' rank repair is applied.  `components = NULL` returns all numerically usable
#' positive-variance components.
#'
#' Eigenvector signs are anchored by making the first coordinate attaining the
#' largest absolute loading positive.  Individual eigenvectors inside a
#' repeated-eigenvalue block are not identified; the returned projector is the
#' invariant object for such a block.
#'
#' @param x Numeric matrix or data frame with observations in rows.
#' @param components Number of leading components. `NULL` uses the numerical
#'   rank certified at `eigen_tol`.
#' @param center Either `"mean"`, `"none"`, or a supplied finite vector with
#'   one entry per variable. Named supplied centers are reordered to the data.
#' @param covariance_divisor Either `"n-1"` or `"n"`.
#' @param eigen_tol Positive relative tolerance for PSD, numerical-rank, and
#'   repeated-eigenvalue diagnostics.
#'
#' @return A `classical_pca_fit` and `hd_pca_fit` object containing loadings,
#'   scores, covariance operator, eigenvalues, projector, reconstruction, and
#'   complete centering/divisor diagnostics.
#' @export
#'
#' @references
#' Anderson, T. W. (2003). *An Introduction to Multivariate Statistical
#' Analysis*, 3rd ed. Wiley.
#'
#' @examples
#' x <- rbind(c(2, 0), c(1, 1), c(-1, -1), c(-2, 0))
#' fit <- classical_pca(x, components = 1, covariance_divisor = "n")
#' fit$loadings
classical_pca <- function(
    x, components = NULL, center = c("mean", "none"),
    covariance_divisor = c("n-1", "n"),
    eigen_tol = sqrt(.Machine$double.eps)) {
  call <- match.call()
  data.name <- deparse1(substitute(x))
  x <- .ch6cf_data(x, "x", min_rows = 2L)
  n <- nrow(x)
  p <- ncol(x)
  center.fit <- .ch6cf_center(x, center, c("mean", "none"))
  divisor.fit <- .ch6cf_divisor(covariance_divisor, n)
  moments <- cpp_ch6cf_center_moments(
    x, center.fit$value, divisor.fit$value
  )
  covariance <- .ch6cf_dimnames_square(
    as.matrix(moments$covariance), colnames(x)
  )
  spectral <- .ch6cf_psd_eigen(
    covariance, "the PCA covariance operator", eigen_tol
  )
  if (spectral$rank < 1L) {
    stop("The PCA operator has numerical rank zero at `eigen_tol`.",
         call. = FALSE)
  }
  rank <- .ch6cf_component_count(
    components, maximum = spectral$rank, default = spectral$rank,
    name = "components"
  )
  boundary <- .ch6cf_boundary_gap(
    spectral$values, rank, spectral$tolerance, "components",
    require_gap = TRUE
  )
  loadings <- spectral$vectors[, seq_len(rank), drop = FALSE]
  component.names <- paste0("PC", seq_len(rank))
  dimnames(loadings) <- list(colnames(x), component.names)
  centered <- as.matrix(moments$centered)
  dimnames(centered) <- dimnames(x)
  scores <- centered %*% loadings
  dimnames(scores) <- list(rownames(x), component.names)
  reconstruction.centered <- scores %*% t(loadings)
  reconstruction <- sweep(
    reconstruction.centered, 2L, center.fit$value, "+"
  )
  dimnames(reconstruction) <- dimnames(x)
  projector <- tcrossprod(loadings)
  projector <- .ch6cf_dimnames_square(projector, colnames(x))
  total <- sum(spectral$values)
  explained <- spectral$values[seq_len(rank)] / total
  names(explained) <- component.names
  eigenvalues <- spectral$values
  names(eigenvalues) <- paste0("PC", seq_along(eigenvalues))
  structure(list(
    valid = TRUE,
    method = "Classical principal component analysis",
    eigenvalues = eigenvalues,
    loadings = loadings,
    scores = scores,
    center = center.fit$value,
    operator = covariance,
    rank = rank,
    n = n,
    p = p,
    variable.names = colnames(x),
    observation.names = rownames(x),
    projector = projector,
    reconstruction = reconstruction,
    explained.variance = explained,
    cumulative.explained.variance = cumsum(explained),
    diagnostics = list(
      centering = center.fit$source,
      covariance.divisor = divisor.fit$label,
      covariance.divisor.value = divisor.fit$value,
      numerical.rank = spectral$rank,
      eigen.tolerance = spectral$tolerance,
      boundary.gap = boundary$gap,
      boundary.separated = boundary$separated,
      repeated.eigenvalue.blocks = spectral$tie.blocks,
      sign.anchor.coordinates = spectral$anchors[seq_len(rank)],
      sign.anchor.rule = paste(
        "first coordinate attaining maximum absolute loading is positive"
      ),
      negative.roundoff.clipped = spectral$negative.roundoff.clipped,
      repeated.direction.contract = paste(
        "individual eigenvectors are not identified inside a repeated block;",
        "compare the returned projector"
      ),
      numerical.repair = "none; only certified negative roundoff is set to zero"
    ),
    call = call,
    data.name = data.name
  ), class = c("classical_pca_fit", "hd_pca_fit", "list"))
}


#' Classical canonical correlation analysis
#'
#' Forms all three covariance blocks with one explicit divisor, constructs the
#' symmetric inverse square roots of the two within-block covariance matrices,
#' and applies an SVD to the whitened cross-covariance matrix. Both within-block
#' matrices must be strictly positive definite and full rank at `rank_tol`.
#' No ridge, pseudoinverse, or eigenvalue floor is used.
#'
#' @param x,y Paired numeric data matrices with equal row counts.
#' @param components Number of leading canonical pairs. `NULL` returns all
#'   `min(ncol(x), ncol(y))` pairs after the within-block rank certificates.
#' @param center Either `"mean"` or `"none"`, applied to both blocks.
#' @param covariance_divisor Either `"n-1"` or `"n"` for every covariance
#'   block; canonical correlations are invariant to this common choice.
#' @param rank_tol Positive relative tolerance for strict within-block rank and
#'   repeated-canonical-correlation diagnostics.
#'
#' @return A `classical_cca_fit` object with canonical coefficients, scores,
#'   structure loadings, correlations, covariance blocks, and whitened
#'   singular-vector projectors.
#' @export
#'
#' @references
#' Hotelling, H. (1936). Relations between two sets of variates.
#' *Biometrika*, 28, 321--377.
#'
#' @examples
#' basis <- stats::poly(seq_len(8), degree = 4)
#' x <- basis[, 1:2, drop = FALSE]
#' y <- cbind(0.8 * basis[, 1] + 0.6 * basis[, 3],
#'            0.3 * basis[, 2] + sqrt(0.91) * basis[, 4])
#' colnames(x) <- c("x1", "x2")
#' colnames(y) <- c("y1", "y2")
#' fit <- classical_cca(x, y, components = 1)
#' fit$canonical.correlations
classical_cca <- function(
    x, y, components = NULL, center = c("mean", "none"),
    covariance_divisor = c("n-1", "n"),
    rank_tol = sqrt(.Machine$double.eps)) {
  call <- match.call()
  data.name <- paste(deparse1(substitute(x)), "and", deparse1(substitute(y)))
  x <- .ch6cf_data(x, "x", min_rows = 2L)
  y <- .ch6cf_data(y, "y", min_rows = 2L)
  if (nrow(x) != nrow(y)) {
    stop("`x` and `y` must have the same number of paired rows.",
         call. = FALSE)
  }
  if (!is.null(rownames(x)) && !is.null(rownames(y)) &&
      !identical(rownames(x), rownames(y))) {
    stop("Non-null row names of paired `x` and `y` must be identical.",
         call. = FALSE)
  }
  center <- match.arg(center)
  n <- nrow(x)
  px <- ncol(x)
  py <- ncol(y)
  center.x <- .ch6cf_center(x, center, c("mean", "none"), "center")
  center.y <- .ch6cf_center(y, center, c("mean", "none"), "center")
  divisor.fit <- .ch6cf_divisor(covariance_divisor, n)
  rank_tol <- .ch6cf_positive(rank_tol, "rank_tol")
  moments <- cpp_ch6cf_cca_moments(
    x, y, center.x$value, center.y$value, divisor.fit$value
  )
  sxx <- .ch6cf_dimnames_square(as.matrix(moments$covariance_xx),
                                 colnames(x))
  syy <- .ch6cf_dimnames_square(as.matrix(moments$covariance_yy),
                                 colnames(y))
  sxy <- as.matrix(moments$covariance_xy)
  dimnames(sxy) <- list(colnames(x), colnames(y))
  fit.x <- .ch6cf_spd_inverse_sqrt(sxx, "the x covariance block", rank_tol)
  fit.y <- .ch6cf_spd_inverse_sqrt(syy, "the y covariance block", rank_tol)
  whitened <- fit.x$inverse.sqrt %*% sxy %*% fit.y$inverse.sqrt
  if (anyNA(whitened) || any(!is.finite(whitened))) {
    stop("The whitened cross-covariance is outside the finite double range.",
         call. = FALSE)
  }
  maximum <- min(px, py)
  decomposition <- svd(whitened, nu = maximum, nv = maximum)
  correlations <- as.numeric(decomposition$d[seq_len(maximum)])
  if (any(correlations > 1 + rank_tol)) {
    stop(paste(
      "A canonical correlation exceeds one beyond `rank_tol`;",
      "the joint covariance blocks failed their numerical consistency check."
    ), call. = FALSE)
  }
  rounded <- correlations > 1
  correlations[rounded] <- 1
  rank <- .ch6cf_component_count(
    components, maximum = maximum, default = maximum, name = "components"
  )
  boundary <- .ch6cf_boundary_gap(
    correlations, rank, rank_tol, "components", require_gap = TRUE
  )
  u <- decomposition$u[, seq_len(rank), drop = FALSE]
  v <- decomposition$v[, seq_len(rank), drop = FALSE]
  x.coefficients <- fit.x$inverse.sqrt %*% u
  y.coefficients <- fit.y$inverse.sqrt %*% v
  anchored <- .ch6cf_anchor_columns(x.coefficients, y.coefficients)
  x.coefficients <- anchored$primary
  y.coefficients <- anchored$paired
  u <- fit.x$inverse.sqrt %*% x.coefficients
  u <- sxx %*% u
  v <- fit.y$inverse.sqrt %*% y.coefficients
  v <- syy %*% v
  pair.names <- paste0("CC", seq_len(rank))
  dimnames(x.coefficients) <- list(colnames(x), pair.names)
  dimnames(y.coefficients) <- list(colnames(y), pair.names)
  centered.x <- as.matrix(moments$centered_x)
  centered.y <- as.matrix(moments$centered_y)
  dimnames(centered.x) <- dimnames(x)
  dimnames(centered.y) <- dimnames(y)
  joint.singular.values <- svd(
    cbind(centered.x, centered.y), nu = 0L, nv = 0L
  )$d
  joint.rank <- sum(
    joint.singular.values > rank_tol * max(1, joint.singular.values[1L])
  )
  x.scores <- centered.x %*% x.coefficients
  y.scores <- centered.y %*% y.coefficients
  dimnames(x.scores) <- list(rownames(x), pair.names)
  dimnames(y.scores) <- list(rownames(y), pair.names)
  x.loadings <- sweep(sxx %*% x.coefficients, 1L, sqrt(diag(sxx)), "/")
  y.loadings <- sweep(syy %*% y.coefficients, 1L, sqrt(diag(syy)), "/")
  dimnames(x.loadings) <- list(colnames(x), pair.names)
  dimnames(y.loadings) <- list(colnames(y), pair.names)
  x.whitened <- fit.x$inverse.sqrt %*% x.coefficients
  x.whitened <- sxx %*% x.whitened
  y.whitened <- fit.y$inverse.sqrt %*% y.coefficients
  y.whitened <- syy %*% y.whitened
  x.projector <- tcrossprod(x.whitened)
  y.projector <- tcrossprod(y.whitened)
  dimnames(x.projector) <- list(colnames(x), colnames(x))
  dimnames(y.projector) <- list(colnames(y), colnames(y))
  all.correlations <- correlations
  names(all.correlations) <- paste0("CC", seq_along(all.correlations))
  structure(list(
    valid = TRUE,
    method = "Classical canonical correlation analysis",
    canonical.correlations = all.correlations[seq_len(rank)],
    all.canonical.correlations = all.correlations,
    x.coefficients = x.coefficients,
    y.coefficients = y.coefficients,
    x.scores = x.scores,
    y.scores = y.scores,
    x.loadings = x.loadings,
    y.loadings = y.loadings,
    center = list(x = center.x$value, y = center.y$value),
    covariance = list(xx = sxx, xy = sxy, yx = t(sxy), yy = syy),
    operator = whitened,
    rank = rank,
    n = n,
    p.x = px,
    p.y = py,
    variable.names = list(x = colnames(x), y = colnames(y)),
    observation.names = if (!is.null(rownames(x))) rownames(x) else rownames(y),
    projector = list(x.whitened = x.projector,
                     y.whitened = y.projector),
    diagnostics = list(
      centering = center,
      covariance.divisor = divisor.fit$label,
      covariance.divisor.value = divisor.fit$value,
      rank.tolerance = rank_tol,
      reciprocal.condition.x = fit.x$reciprocal.condition,
      reciprocal.condition.y = fit.y$reciprocal.condition,
      within.rank.x = px,
      within.rank.y = py,
      joint.rank = as.integer(joint.rank),
      joint.full.rank = joint.rank == px + py,
      boundary.gap = boundary$gap,
      boundary.separated = boundary$separated,
      repeated.correlation.blocks = .ch6cf_tie_blocks(
        correlations, rank_tol
      ),
      sign.anchor.coordinates = anchored$anchors,
      sign.anchor.rule = paste(
        "the x coefficient's first maximum-absolute coordinate is positive;",
        "the paired y coefficient is flipped jointly"
      ),
      correlations.above.one.roundoff.clipped = sum(rounded),
      repeated.direction.contract = paste(
        "canonical vectors inside a repeated singular-value block are not",
        "identified; compare the whitened projectors"
      ),
      numerical.repair = "none; no ridge, pseudoinverse, or eigenvalue floor"
    ),
    call = call,
    data.name = data.name
  ), class = c("classical_cca_fit", "hd_cca_fit", "list"))
}


#' Bartlett likelihood-ratio approximation for canonical correlations
#'
#' Tests the tail null that canonical correlations after `null_rank` are zero.
#' For `null_rank = k`, the statistic uses
#' \deqn{\Lambda_k=\prod_{j=k+1}^{m}(1-\widehat\rho_j^2),\qquad
#' -\{n-1-(p_x+p_y+1)/2\}\log\Lambda_k,}
#' with chi-squared degrees of freedom `(p_x-k)(p_y-k)`.  The special case
#' `k = 0` is the full independence test printed in the book.
#'
#' @param object A valid `classical_cca_fit` produced with mean centering.
#' @param null_rank Non-negative integer rank under the tail null.
#'
#' @return An object of class `htest` with Wilks' Lambda and the Bartlett
#'   approximation diagnostics.
#' @export
#'
#' @references
#' Anderson, T. W. (2003). *An Introduction to Multivariate Statistical
#' Analysis*, 3rd ed. Wiley.
#'
#' @examples
#' basis <- stats::poly(seq_len(12), degree = 4)
#' x <- basis[, 1:2, drop = FALSE]
#' y <- cbind(0.8 * basis[, 1] + 0.6 * basis[, 3],
#'            0.3 * basis[, 2] + sqrt(0.91) * basis[, 4])
#' colnames(x) <- c("x1", "x2")
#' colnames(y) <- c("y1", "y2")
#' fit <- classical_cca(x, y)
#' cca_bartlett_test(fit, null_rank = 0)
cca_bartlett_test <- function(object, null_rank = 0L) {
  call <- match.call()
  if (!inherits(object, "classical_cca_fit") || !isTRUE(object$valid)) {
    stop("`object` must be a valid `classical_cca_fit`.", call. = FALSE)
  }
  if (!identical(object$diagnostics$centering, "mean")) {
    stop(paste(
      "Bartlett's stated sample-size correction requires a mean-centered",
      "classical CCA fit."
    ), call. = FALSE)
  }
  correlations <- as.numeric(object$all.canonical.correlations)
  m <- length(correlations)
  null_rank <- .ch6cf_count(
    null_rank, "null_rank", minimum = 0L, maximum = m - 1L
  )
  if (!isTRUE(object$diagnostics$joint.full.rank)) {
    stop(paste(
      "Bartlett's likelihood-ratio approximation requires the jointly",
      "centered [X, Y] matrix to have full column rank; no generalized",
      "determinant or pseudoinverse substitute is used."
    ), call. = FALSE)
  }
  tail <- correlations[seq.int(null_rank + 1L, m)]
  terms <- 1 - tail^2
  if (anyNA(terms) || any(!is.finite(terms)) || any(terms <= 0)) {
    stop(paste(
      "`log(1 - rho^2)` is not finite and defined for every tested",
      "canonical correlation; Bartlett's approximation is not reported."
    ), call. = FALSE)
  }
  sample.factor <- object$n - 1 - (object$p.x + object$p.y + 1) / 2
  if (!is.finite(sample.factor) || sample.factor <= 0) {
    stop(paste(
      "The Bartlett sample-size factor must be strictly positive;",
      "the fixed-dimensional approximation is unavailable for this fit."
    ), call. = FALSE)
  }
  log.wilks <- sum(log1p(-tail^2))
  statistic <- -sample.factor * log.wilks
  degrees <- (object$p.x - null_rank) * (object$p.y - null_rank)
  p.value <- stats::pchisq(statistic, degrees, lower.tail = FALSE)
  structure(list(
    statistic = c("X-squared" = statistic),
    parameter = c(df = degrees),
    p.value = p.value,
    method = sprintf(
      "Bartlett approximation to Wilks' Lambda for CCA rank <= %d",
      null_rank
    ),
    data.name = object$data.name,
    null.value = c(rank = null_rank),
    alternative = sprintf("canonical rank is greater than %d", null_rank),
    wilks.lambda = exp(log.wilks),
    log.wilks.lambda = log.wilks,
    sample.size.factor = sample.factor,
    tested.canonical.correlations = tail,
    call = call
  ), class = "htest")
}

.ch6cf_robust_center <- function(x, center, tol, max_iter, zero_tol) {
  if (is.character(center)) {
    center <- match.arg(center, c("spatial", "mean", "none"))
  }
  if (is.character(center) && identical(center, "spatial")) {
    controls <- .validate_iteration_controls(tol, max_iter, zero_tol)
    location <- spatial_median(
      x, tol = controls$tol, max_iter = controls$max_iter,
      zero_tol = controls$zero_tol, warn = FALSE
    )
    diagnostics <- list(
      objective = as.numeric(attr(location, "objective")),
      iterations = as.integer(attr(location, "iterations")),
      converged = isTRUE(attr(location, "converged")),
      relative.change = as.numeric(attr(location, "relative_change")),
      equation.residual = as.numeric(attr(location, "equation_residual"))
    )
    if (!diagnostics$converged) {
      stop(paste(
        "The spatial median did not satisfy its equation-residual",
        "certificate; no robust principal subspace is returned."
      ), call. = FALSE)
    }
    value <- as.numeric(location)
    names(value) <- colnames(x)
    return(list(value = value, source = "spatial", diagnostics = diagnostics,
                controls = controls))
  }
  fit <- .ch6cf_center(x, center, c("mean", "none"))
  list(value = fit$value, source = fit$source, diagnostics = NULL,
       controls = .validate_iteration_controls(tol, max_iter, zero_tol))
}


.ch6cf_robust_operator <- function(x, method, center.fit, zero_action) {
  zero_action <- match.arg(zero_action, c("error", "zero"))
  if (method == "spatial_sign") {
    operator <- sscm(
      x, center = center.fit$value,
      zero_tol = center.fit$controls$zero_tol
    )
    n.zero <- as.integer(attr(operator, "n_zero"))
    if (n.zero > 0L && zero_action == "error") {
      stop(sprintf(
        paste(
          "The selected center produces %d zero spatial residual(s);",
          "choose `zero_action = \"zero\"` to retain their explicit zero",
          "contribution with denominator n."
        ), n.zero
      ), call. = FALSE)
    }
    diagnostics <- list(
      zero.count = n.zero,
      total.count = nrow(x),
      normalization = "sum of sign outer products divided by n",
      zero.contribution = zero_action,
      translation.invariant.operator = FALSE
    )
  } else {
    operator <- spatial_kendall(
      x, zero_tol = center.fit$controls$zero_tol
    )
    n.zero <- as.integer(attr(operator, "n_zero_pairs"))
    n.pairs <- as.integer(attr(operator, "n_pairs"))
    if (n.zero > 0L && zero_action == "error") {
      stop(sprintf(
        paste(
          "The data contain %d tied pairwise difference(s); choose",
          "`zero_action = \"zero\"` to retain their explicit zero",
          "contribution with the all-pairs denominator."
        ), n.zero
      ), call. = FALSE)
    }
    diagnostics <- list(
      zero.count = n.zero,
      total.count = n.pairs,
      normalization = "2 / {n(n-1)} over all unordered pairs",
      zero.contribution = zero_action,
      translation.invariant.operator = TRUE
    )
  }
  list(operator = operator, diagnostics = diagnostics)
}


#' Robust principal subspace for an elliptical factor model
#'
#' Estimates a supplied-dimensional principal subspace from either the sample
#' spatial-sign covariance matrix or the exact multivariate Kendall matrix.
#' The result is deliberately called a `principal_subspace`: a finite-sample
#' equality with the span of an unknown factor-loading matrix is not asserted.
#'
#' The requested boundary must have a certified eigengap. Signs are anchored
#' deterministically, while directions within a repeated eigenvalue block are
#' represented by their projector. Zero spatial residuals or tied Kendall
#' pairs are errors unless their zero contribution is explicitly requested.
#'
#' @param x Numeric matrix or data frame with observations in rows.
#' @param factors Supplied positive principal-subspace dimension.
#' @param method Either `"spatial_sign"` or `"kendall"`.
#' @param center `"spatial"`, `"mean"`, `"none"`, or a supplied finite center.
#'   Kendall's operator is translation invariant, but this center still defines
#'   the returned component scores.
#' @param tol,max_iter Spatial-median iteration controls.
#' @param zero_tol Non-negative zero-residual or tied-pair tolerance.
#' @param zero_action Either error on zero directions or retain their explicit
#'   zero contribution under the primary denominator.
#' @param eigen_tol Positive relative eigengap and PSD tolerance.
#'
#' @return A `robust_factor_subspace_fit` and `hd_factor_fit` object containing
#'   the orthonormal principal-subspace basis, projector, robust operator, and
#'   centered component scores.
#' @export
#'
#' @references
#' Han, F. and Liu, H. (2018). ECA: High-dimensional elliptical component
#' analysis in non-Gaussian distributions. *JASA*, 113, 252--268.
#'
#' @examples
#' x <- rbind(c(-3, -1), c(-2, 1), c(-1, -2),
#'            c(1, 2), c(2, -1), c(3, 1))
#' fit <- robust_factor_subspace(x, factors = 1, method = "kendall",
#'                               center = "mean")
#' fit$principal_subspace
robust_factor_subspace <- function(
    x, factors, method = c("spatial_sign", "kendall"),
    center = c("spatial", "mean", "none"),
    tol = 1e-8, max_iter = 500L, zero_tol = 0,
    zero_action = c("error", "zero"),
    eigen_tol = sqrt(.Machine$double.eps)) {
  call <- match.call()
  data.name <- deparse1(substitute(x))
  x <- .ch6cf_data(x, "x", min_rows = 2L)
  method <- match.arg(method)
  zero_action <- match.arg(zero_action)
  p <- ncol(x)
  n <- nrow(x)
  if (p < 2L) {
    stop("A factor principal subspace requires at least two variables.",
         call. = FALSE)
  }
  factors <- .ch6cf_count(factors, "factors", 1L, p - 1L)
  center.fit <- .ch6cf_robust_center(
    x, center, tol, max_iter, zero_tol
  )
  robust <- .ch6cf_robust_operator(
    x, method, center.fit, zero_action
  )
  operator <- .ch6cf_dimnames_square(
    as.matrix(robust$operator), colnames(x)
  )
  spectral <- .ch6cf_psd_eigen(
    operator, "the robust factor operator", eigen_tol
  )
  if (factors > spectral$rank) {
    stop(sprintf(
      "`factors` exceeds the robust operator's numerical rank (%d).",
      spectral$rank
    ), call. = FALSE)
  }
  boundary <- .ch6cf_boundary_gap(
    spectral$values, factors, spectral$tolerance, "factors",
    require_gap = TRUE
  )
  basis <- spectral$vectors[, seq_len(factors), drop = FALSE]
  factor.names <- paste0("RF", seq_len(factors))
  dimnames(basis) <- list(colnames(x), factor.names)
  projector <- tcrossprod(basis)
  projector <- .ch6cf_dimnames_square(projector, colnames(x))
  centered <- sweep(x, 2L, center.fit$value, "-")
  scores <- centered %*% basis
  dimnames(scores) <- list(rownames(x), factor.names)
  eigenvalues <- spectral$values
  names(eigenvalues) <- paste0("EV", seq_along(eigenvalues))
  structure(list(
    valid = TRUE,
    method = if (method == "spatial_sign")
      "Spatial-sign robust factor principal subspace" else
      "Kendall robust factor principal subspace",
    principal_subspace = basis,
    loadings = basis,
    factor.scores = NULL,
    common.component = NULL,
    scores = scores,
    projector = projector,
    eigenvalues = eigenvalues,
    factor.number = factors,
    center = center.fit$value,
    operator = operator,
    rank = factors,
    n = n,
    p = p,
    variable.names = colnames(x),
    observation.names = rownames(x),
    diagnostics = list(
      operator = method,
      centering = center.fit$source,
      spatial.median = center.fit$diagnostics,
      zero.directions = robust$diagnostics,
      numerical.rank = spectral$rank,
      eigen.tolerance = spectral$tolerance,
      boundary.gap = boundary$gap,
      repeated.eigenvalue.blocks = spectral$tie.blocks,
      sign.anchor.coordinates = spectral$anchors[seq_len(factors)],
      sign.anchor.rule = paste(
        "first coordinate attaining maximum absolute loading is positive"
      ),
      target = paste(
        "sample principal subspace of the selected robust operator; no",
        "finite-sample equality with span(B) is claimed"
      ),
      repeated.direction.contract = paste(
        "directions inside a repeated block are represented by the projector"
      ),
      numerical.repair = "none"
    ),
    call = call,
    data.name = data.name
  ), class = c("robust_factor_subspace_fit", "hd_factor_fit", "list"))
}


.ch6cf_kendall_spectral <- function(
    x, zero_tol, zero_pair_action, eigen_tol, name) {
  zero_tol <- .validate_zero_tol(zero_tol)
  zero_pair_action <- match.arg(zero_pair_action, c("error", "zero"))
  operator <- spatial_kendall(x, zero_tol = zero_tol)
  n.pairs <- as.integer(attr(operator, "n_pairs"))
  n.zero <- as.integer(attr(operator, "n_zero_pairs"))
  if (n.zero > 0L && zero_pair_action == "error") {
    stop(sprintf(
      paste(
        "The data contain %d tied pairwise difference(s); choose",
        "`zero_pair_action = \"zero\"` to retain their explicit zero",
        "contribution under the all-pairs denominator."
      ), n.zero
    ), call. = FALSE)
  }
  if (n.zero == n.pairs) {
    stop("All pairwise differences are zero; the Kendall operator is undefined as directional information.",
         call. = FALSE)
  }
  operator <- .ch6cf_dimnames_square(as.matrix(operator), colnames(x))
  spectral <- .ch6cf_psd_eigen(operator, name, eigen_tol)
  list(
    operator = operator, spectral = spectral, n.pairs = n.pairs,
    n.zero = n.zero, zero.action = zero_pair_action, zero.tol = zero_tol
  )
}


#' Robust two-step factor estimator based on multivariate Kendall's tau
#'
#' Implements the primary RTS estimator.  If `Gamma` contains the supplied top
#' Kendall eigenvectors, loadings are `sqrt(p) * Gamma`, factor scores are the
#' cross-sectional OLS estimates `X_centered %*% loadings / p`, and the common
#' component is `factor_scores %*% t(loadings)`.  Centering affects the latter
#' two quantities and is therefore explicit even though the Kendall operator is
#' translation invariant.
#'
#' @param x Numeric matrix or data frame with time/observations in rows and
#'   cross-sectional variables in columns.
#' @param factors Supplied positive factor number.
#' @param center `"mean"`, `"none"`, or a supplied finite center vector.
#' @param zero_tol Non-negative tolerance for tied pairwise differences.
#' @param zero_pair_action Either error on ties or retain their explicit zero
#'   contribution under the all-pairs Kendall denominator.
#' @param eigen_tol Positive relative PSD, numerical-rank, and eigengap
#'   tolerance.
#'
#' @return An `rts_factor_fit` and `hd_factor_fit` object with loading matrix,
#'   OLS factor scores, common component, residuals, fitted data, principal
#'   subspace, and complete Kendall diagnostics.
#' @export
#'
#' @references
#' He, Y., Kong, X., Yu, L. and Zhang, X. (2022). Large-dimensional factor
#' analysis without moment constraints. *Journal of Business & Economic
#' Statistics*, 40, 302--312.
#'
#' @examples
#' x <- rbind(c(-3, -1, 0), c(-2, 1, 1), c(-1, -2, 0),
#'            c(1, 2, 0), c(2, -1, -1), c(3, 1, 0))
#' fit <- rts_factor(x, factors = 1)
#' crossprod(fit$loadings) / ncol(x)
rts_factor <- function(
    x, factors, center = c("mean", "none"), zero_tol = 0,
    zero_pair_action = c("error", "zero"),
    eigen_tol = sqrt(.Machine$double.eps)) {
  call <- match.call()
  data.name <- deparse1(substitute(x))
  x <- .ch6cf_data(x, "x", min_rows = 2L)
  n <- nrow(x)
  p <- ncol(x)
  if (p < 2L) {
    stop("RTS factor analysis requires at least two variables.",
         call. = FALSE)
  }
  factors <- .ch6cf_count(
    factors, "factors", minimum = 1L, maximum = p - 1L
  )
  center.fit <- .ch6cf_center(x, center, c("mean", "none"))
  kendall <- .ch6cf_kendall_spectral(
    x, zero_tol, zero_pair_action, eigen_tol,
    "the RTS Kendall operator"
  )
  if (factors > kendall$spectral$rank) {
    stop(sprintf(
      "`factors` exceeds the Kendall operator's numerical rank (%d).",
      kendall$spectral$rank
    ), call. = FALSE)
  }
  boundary <- .ch6cf_boundary_gap(
    kendall$spectral$values, factors, kendall$spectral$tolerance,
    "factors", require_gap = TRUE
  )
  principal <- kendall$spectral$vectors[, seq_len(factors), drop = FALSE]
  factor.names <- paste0("F", seq_len(factors))
  dimnames(principal) <- list(colnames(x), factor.names)
  loadings <- sqrt(p) * principal
  dimnames(loadings) <- list(colnames(x), factor.names)
  components <- cpp_ch6cf_rts_components(x, center.fit$value, loadings)
  factor.scores <- as.matrix(components$factor_scores)
  common.component <- as.matrix(components$common_component)
  residuals <- as.matrix(components$residuals)
  fitted <- sweep(common.component, 2L, center.fit$value, "+")
  dimnames(factor.scores) <- list(rownames(x), factor.names)
  dimnames(common.component) <- dimnames(x)
  dimnames(residuals) <- dimnames(x)
  dimnames(fitted) <- dimnames(x)
  projector <- tcrossprod(principal)
  projector <- .ch6cf_dimnames_square(projector, colnames(x))
  eigenvalues <- kendall$spectral$values
  names(eigenvalues) <- paste0("EV", seq_along(eigenvalues))
  structure(list(
    valid = TRUE,
    method = "Robust two-step Kendall factor estimator",
    principal_subspace = principal,
    loadings = loadings,
    factor.scores = factor.scores,
    common.component = common.component,
    residuals = residuals,
    fitted = fitted,
    projector = projector,
    eigenvalues = eigenvalues,
    factor.number = factors,
    center = center.fit$value,
    operator = kendall$operator,
    rank = factors,
    n = n,
    p = p,
    variable.names = colnames(x),
    observation.names = rownames(x),
    diagnostics = list(
      centering = center.fit$source,
      kendall.normalization = "2 / {n(n-1)} over all unordered pairs",
      kendall.pairs = kendall$n.pairs,
      zero.kendall.pairs = kendall$n.zero,
      zero.pair.action = kendall$zero.action,
      zero.tolerance = kendall$zero.tol,
      numerical.rank = kendall$spectral$rank,
      eigen.tolerance = kendall$spectral$tolerance,
      boundary.gap = boundary$gap,
      repeated.eigenvalue.blocks = kendall$spectral$tie.blocks,
      sign.anchor.coordinates =
        kendall$spectral$anchors[seq_len(factors)],
      loading.normalization = "crossprod(loadings) / p = I",
      factor.score.estimator = "cross-sectional ordinary least squares",
      normal.equation.residual =
        as.numeric(components$normal_equation_residual),
      target = paste(
        "Kendall principal subspace; identified up to orthogonal rotation",
        "and sign, not asserted equal to span(B) in finite samples"
      ),
      numerical.repair = "none"
    ),
    call = call,
    data.name = data.name
  ), class = c("rts_factor_fit", "hd_factor_fit", "list"))
}


#' Robust factor-number selection from multivariate Kendall eigenvalues
#'
#' Implements the modified Kendall eigenvalue-ratio (MKER) and transformed
#' contribution-ratio (MKTCR) selectors.  With `m = min(N, T)` for `N = p`
#' variables and `T = n` observations, the stabilized eigenvalues are
#' \deqn{\widetilde\lambda_j=\lambda_j+c/\sqrt{m}.}
#' MKER maximizes `tilde_lambda[j] / tilde_lambda[j + 1]`.  For MKTCR,
#' \deqn{V_j=\sum_{i=j+1}^{m}\widetilde\lambda_i,\qquad
#' \frac{\log(1+\widetilde\lambda_j/V_{j-1})}
#' {\log(1+\widetilde\lambda_{j+1}/V_j)}}
#' is maximized. `kmax` and strictly positive `c` are always supplied, and an
#' exact tie selects the smallest index.
#'
#' @param x Numeric matrix or data frame with observations in rows.
#' @param kmax Supplied positive upper bound, at most `min(n, p) - 1`.
#' @param c Supplied strictly positive stabilization constant.
#' @param method Either `"mker"` or `"mktcr"`.
#' @param zero_tol Non-negative tolerance for tied pairwise differences.
#' @param zero_pair_action Either error on ties or retain their explicit zero
#'   contribution under the all-pairs denominator.
#' @param eigen_tol Positive relative PSD tolerance.
#'
#' @return A `kendall_factor_number` object with the selected factor number,
#'   raw and stabilized eigenvalues, criterion path, and MKTCR tail sums.
#' @export
#'
#' @references
#' Yu, L., He, Y. and Zhang, X. (2019). Robust factor number specification for
#' large-dimensional elliptical factor models. arXiv:1808.09107.
#'
#' @examples
#' x <- rbind(c(-3, -1, 0), c(-2, 1, 1), c(-1, -2, 0),
#'            c(1, 2, 0), c(2, -1, -1), c(3, 1, 0))
#' kendall_factor_number(x, kmax = 1, c = 0.1, method = "mker")
kendall_factor_number <- function(
    x, kmax, c, method = c("mker", "mktcr"), zero_tol = 0,
    zero_pair_action = c("error", "zero"),
    eigen_tol = sqrt(.Machine$double.eps)) {
  call <- match.call()
  data.name <- deparse1(substitute(x))
  x <- .ch6cf_data(x, "x", min_rows = 2L)
  method <- match.arg(method)
  n <- nrow(x)
  p <- ncol(x)
  m <- min(n, p)
  if (m < 2L) {
    stop("Factor-number selection requires `min(n, p) >= 2`.",
         call. = FALSE)
  }
  kmax <- .ch6cf_count(kmax, "kmax", minimum = 1L, maximum = m - 1L)
  constant <- .ch6cf_positive(c, "c")
  kendall <- .ch6cf_kendall_spectral(
    x, zero_tol, zero_pair_action, eigen_tol,
    "the factor-number Kendall operator"
  )
  raw <- kendall$spectral$values[seq_len(m)]
  stabilization <- constant / sqrt(m)
  adjusted <- raw + stabilization
  if (anyNA(adjusted) || any(!is.finite(adjusted)) || any(adjusted <= 0)) {
    stop("The stabilized Kendall eigenvalues are not finite and positive.",
         call. = FALSE)
  }
  if (method == "mker") {
    criterion <- adjusted[seq_len(kmax)] /
      adjusted[seq.int(2L, kmax + 1L)]
    tail.sums <- NULL
  } else {
    tail.sums <- vapply(0:(m - 1L), function(j) {
      sum(adjusted[seq.int(j + 1L, m)])
    }, numeric(1))
    names(tail.sums) <- paste0("V", 0:(m - 1L))
    criterion <- vapply(seq_len(kmax), function(j) {
      numerator <- log1p(adjusted[j] / tail.sums[j])
      denominator <- log1p(adjusted[j + 1L] / tail.sums[j + 1L])
      if (!is.finite(numerator) || !is.finite(denominator) ||
          denominator <= 0) {
        stop("An MKTCR logarithmic ratio is not finite and positive.",
             call. = FALSE)
      }
      numerator / denominator
    }, numeric(1))
  }
  names(raw) <- paste0("lambda", seq_len(m))
  names(adjusted) <- paste0("tilde.lambda", seq_len(m))
  names(criterion) <- paste0("j", seq_len(kmax))
  maximum <- max(criterion)
  selected <- which(criterion == maximum)[1L]
  structure(list(
    valid = TRUE,
    selected = as.integer(selected),
    factor.number = as.integer(selected),
    method = toupper(method),
    kmax = kmax,
    c = constant,
    stabilization = stabilization,
    eigenvalues = raw,
    adjusted.eigenvalues = adjusted,
    criterion = criterion,
    tail.sums = tail.sums,
    operator = kendall$operator,
    n = n,
    p = p,
    m = m,
    variable.names = colnames(x),
    diagnostics = list(
      stabilization.formula = "c / sqrt(min(N, T))",
      N = p,
      T = n,
      tie.break = "smallest j attaining the exact maximum",
      kendall.normalization = "2 / {T(T-1)} over all unordered pairs",
      kendall.pairs = kendall$n.pairs,
      zero.kendall.pairs = kendall$n.zero,
      zero.pair.action = kendall$zero.action,
      zero.tolerance = kendall$zero.tol,
      eigen.tolerance = kendall$spectral$tolerance,
      negative.roundoff.clipped =
        kendall$spectral$negative.roundoff.clipped,
      tuning.source = "explicit kmax and c",
      numerical.repair = "none"
    ),
    call = call,
    data.name = data.name
  ), class = c("kendall_factor_number", "list"))
}

