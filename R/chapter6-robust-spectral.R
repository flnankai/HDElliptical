.c6rs_flag <- function(x, name) {
  if (!is.logical(x) || length(x) != 1L || is.na(x)) {
    stop(sprintf("`%s` must be TRUE or FALSE.", name), call. = FALSE)
  }
  x
}


.c6rs_positive_integer <- function(x, name, minimum = 1L) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) ||
      !is.finite(x) || x != floor(x) || x < minimum ||
      x > .Machine$integer.max) {
    stop(sprintf("`%s` must be an integer at least %d.", name, minimum),
         call. = FALSE)
  }
  as.integer(x)
}


.c6rs_rank <- function(rank, n, p) {
  if (is.null(rank)) {
    return(as.integer(min(n, p)))
  }
  rank <- .c6rs_positive_integer(rank, "rank")
  if (rank > p) {
    stop("`rank` cannot exceed the number of variables.", call. = FALSE)
  }
  rank
}


.c6rs_names <- function(x) {
  variable.names <- colnames(x)
  if (is.null(variable.names)) {
    variable.names <- paste0("V", seq_len(ncol(x)))
  }
  if (anyNA(variable.names) || any(!nzchar(variable.names)) ||
      anyDuplicated(variable.names)) {
    stop("Column names, when supplied, must be non-empty and unique.",
         call. = FALSE)
  }
  variable.names
}


.c6rs_spatial_center <- function(x, tol, max_iter, zero_tol) {
  controls <- .validate_iteration_controls(tol, max_iter, zero_tol)
  estimate <- spatial_median(
    x, tol = controls$tol, max_iter = controls$max_iter,
    zero_tol = controls$zero_tol, warn = FALSE
  )
  diagnostics <- list(
    converged = isTRUE(attr(estimate, "converged")),
    iterations = as.integer(attr(estimate, "iterations")),
    objective = as.numeric(attr(estimate, "objective")),
    relative.change = as.numeric(attr(estimate, "relative_change")),
    equation.residual = as.numeric(attr(estimate, "equation_residual")),
    tolerance = controls$tol,
    maximum.iterations = controls$max_iter,
    zero.tolerance = controls$zero_tol
  )
  if (!diagnostics$converged) {
    stop(
      "The spatial median did not pass its convergence certificate; " %+%
        "increase `max_iter` or change the center explicitly.",
      call. = FALSE
    )
  }
  list(value = as.numeric(estimate), diagnostics = diagnostics)
}


.c6rs_center <- function(x, center, tol, max_iter, zero_tol,
                         allow_lts = FALSE, lts_steps = 2L) {
  p <- ncol(x)
  variable.names <- .c6rs_names(x)
  if (is.numeric(center)) {
    value <- .as_location(center, p)
    names(value) <- variable.names
    return(list(
      value = value, source = "supplied numeric center",
      diagnostics = list(type = "supplied", converged = TRUE)
    ))
  }
  choices <- c("spatial", "spatial_median", "mean", "none")
  if (allow_lts) {
    choices <- c("kstep_lts", choices)
  }
  center <- match.arg(center, choices)
  if (identical(center, "mean")) {
    value <- colMeans(x)
    names(value) <- variable.names
    return(list(
      value = value, source = "sample mean",
      diagnostics = list(type = "mean", converged = TRUE)
    ))
  }
  if (identical(center, "none")) {
    value <- setNames(rep.int(0, p), variable.names)
    return(list(
      value = value, source = "origin (no centering)",
      diagnostics = list(type = "none", converged = TRUE)
    ))
  }

  spatial <- .c6rs_spatial_center(x, tol, max_iter, zero_tol)
  names(spatial$value) <- variable.names
  if (center %in% c("spatial", "spatial_median")) {
    return(list(
      value = spatial$value, source = "sample spatial median",
      diagnostics = c(list(type = "spatial median"), spatial$diagnostics)
    ))
  }

  lts_steps <- .c6rs_positive_integer(lts_steps, "lts_steps")
  h <- as.integer(floor((nrow(x) + 1) / 2))
  current <- spatial$value
  path <- vector("list", lts_steps)
  for (step in seq_len(lts_steps)) {
    signed <- spatial_sign(x, center = current, zero_tol = zero_tol)
    radii <- as.numeric(attr(signed, "norms"))
    selected <- order(radii, seq_along(radii))[seq_len(h)]
    next.center <- colMeans(x[selected, , drop = FALSE])
    names(next.center) <- variable.names
    path[[step]] <- list(
      step = step, selected = as.integer(selected),
      objective = sum(radii[selected]^2), center = next.center
    )
    current <- next.center
  }
  list(
    value = current,
    source = sprintf("%d-step LTS initialized at the spatial median", lts_steps),
    diagnostics = list(
      type = "k-step LTS", converged = TRUE, h = h,
      steps = lts_steps, path = path,
      spatial.median = spatial$diagnostics
    )
  )
}


.c6rs_scores <- function(x, center, loadings) {
  scaled <- .center_and_scale(x, center, zero_tol = 0)
  scores <- (scaled$x %*% loadings) * scaled$scale
  if (any(!is.finite(scores))) {
    stop("The requested principal-component scores are non-finite.",
         call. = FALSE)
  }
  dimnames(scores) <- list(rownames(x), colnames(loadings))
  scores
}


.c6rs_anchor <- function(loadings) {
  anchors <- integer(ncol(loadings))
  for (j in seq_len(ncol(loadings))) {
    magnitudes <- abs(loadings[, j])
    anchor <- which(magnitudes == max(magnitudes))[1L]
    anchors[j] <- anchor
    if (loadings[anchor, j] < 0) {
      loadings[, j] <- -loadings[, j]
    }
  }
  list(loadings = loadings, anchors = anchors)
}


.c6rs_repeated_groups <- function(values) {
  if (!length(values)) {
    return(list())
  }
  tolerance <- sqrt(.Machine$double.eps) * max(1, max(abs(values)))
  group <- cumsum(c(TRUE, abs(diff(values)) > tolerance))
  groups <- split(seq_along(values), group)
  Filter(function(index) length(index) > 1L, groups)
}


.c6rs_decompose <- function(operator, rank, variable.names) {
  if (!is.matrix(operator) || nrow(operator) != ncol(operator) ||
      any(!is.finite(operator))) {
    stop("The spectral operator must be a finite square matrix.",
         call. = FALSE)
  }
  symmetry.error <- max(abs(operator - t(operator)))
  if (!is.finite(symmetry.error) || symmetry.error != 0) {
    stop("The spectral operator is not exactly symmetric; no implicit " %+%
           "symmetrization is performed.", call. = FALSE)
  }
  decomposition <- eigen(operator, symmetric = TRUE)
  if (any(!is.finite(decomposition$values)) ||
      any(!is.finite(decomposition$vectors))) {
    stop("The symmetric eigendecomposition returned non-finite values.",
         call. = FALSE)
  }
  selected <- seq_len(rank)
  anchored <- .c6rs_anchor(decomposition$vectors[, selected, drop = FALSE])
  component.names <- paste0("PC", selected)
  dimnames(anchored$loadings) <- list(variable.names, component.names)
  list(
    eigenvalues = as.numeric(decomposition$values),
    loadings = anchored$loadings,
    anchors = setNames(anchored$anchors, component.names),
    repeated.groups = .c6rs_repeated_groups(decomposition$values),
    symmetry.error = symmetry.error,
    minimum.eigenvalue = min(decomposition$values)
  )
}


.c6rs_fit <- function(method, method.class, x, rank, center, operator,
                      decomposition, diagnostics, call, keep_operator) {
  scores <- .c6rs_scores(x, center, decomposition$loadings)
  structure(
    list(
      method = method,
      eigenvalues = decomposition$eigenvalues,
      loadings = decomposition$loadings,
      scores = scores,
      center = center,
      operator = if (keep_operator) operator else NULL,
      rank = rank,
      n = nrow(x),
      p = ncol(x),
      variable.names = .c6rs_names(x),
      diagnostics = c(
        diagnostics,
        list(
          sign.anchor.indices = decomposition$anchors,
          sign.anchor.rule = paste(
            "first coordinate attaining the largest absolute loading is positive"
          ),
          repeated.eigenvalue.groups = decomposition$repeated.groups,
          repeated.eigenvalue.contract = paste(
            "individual loading vectors are not identified within a repeated",
            "eigenspace; compare its orthogonal projector"
          ),
          operator.symmetry.error = decomposition$symmetry.error,
          operator.minimum.eigenvalue = decomposition$minimum.eigenvalue,
          no.implicit.regularization = paste(
            "no ridge, pseudoinverse, eigenvalue floor, or hidden",
            "symmetrization is used"
          )
        )
      ),
      call = call
    ),
    class = c(method.class, "hd_pca_fit", "list")
  )
}


.c6rs_validate_cutoffs <- function(values) {
  required <- c("Q1", "Q2", "Q3", "Q3_star")
  original.names <- names(values)
  values <- as.numeric(values)
  if (length(values) != 4L || anyNA(values) || any(!is.finite(values)) ||
      any(values < 0)) {
    stop("`cutoffs` must contain four finite non-negative values.",
         call. = FALSE)
  }
  if (is.null(original.names)) {
    stop("`cutoffs` must be named Q1, Q2, Q3, and Q3_star.",
         call. = FALSE)
  }
  if (!setequal(original.names, required) || anyDuplicated(original.names)) {
    stop("Named `cutoffs` must be exactly Q1, Q2, Q3, and Q3_star.",
         call. = FALSE)
  }
  values <- setNames(values, original.names)[required]
  if (values["Q1"] > values["Q2"] ||
      values["Q2"] > values["Q3"] ||
      values["Q3"] > values["Q3_star"]) {
    stop("`cutoffs` must satisfy Q1 <= Q2 <= Q3 <= Q3_star.",
         call. = FALSE)
  }
  values
}


.c6rs_q1 <- function(base, policy) {
  adjusted <- FALSE
  if (base < 0) {
    if (identical(policy, "error")) {
      stop("The Q1 Wilson--Hilferty base is negative.", call. = FALSE)
    }
    base <- 0
    adjusted <- TRUE
  }
  list(value = base^(3 / 2), base = base, adjusted = adjusted)
}


.c6rs_cutoffs <- function(radii, n, p, method, user,
                          q1_lower_bound, zero_mad) {
  if (identical(method, "user")) {
    if (is.null(user)) {
      stop("`cutoffs` is required when `cutoff = \"user\"`.",
           call. = FALSE)
    }
    return(list(
      values = .c6rs_validate_cutoffs(user),
      source = "user supplied", transformed.center = NA_real_,
      transformed.mad = NA_real_, h = NA_integer_,
      zero.mad = NA, q1.adjusted = FALSE, q1.base = NA_real_
    ))
  }
  transformed <- radii^(2 / 3)
  if (identical(method, "median_mad")) {
    transformed.center <- stats::median(transformed)
    transformed.mad <- stats::median(abs(transformed - transformed.center))
    q2 <- stats::median(radii)
    h <- NA_integer_
    source <- "ordinary median and raw MAD (2024 GSPCA convention)"
  } else {
    h <- as.integer(floor((n + p + 1) / 2))
    if (h > n) {
      stop(
        "The original h-order cutoff is undefined because " %+%
          sprintf("h = %d exceeds n = %d.", h, n),
        call. = FALSE
      )
    }
    transformed.center <- sort(transformed, partial = h)[h]
    deviations <- abs(transformed - transformed.center)
    transformed.mad <- sort(deviations, partial = h)[h]
    q2 <- transformed.center^(3 / 2)
    source <- "original h-order median and h-order MAD (2019 convention)"
  }
  is.zero.mad <- identical(as.numeric(transformed.mad), 0)
  if (is.zero.mad && identical(zero_mad, "error")) {
    stop("The transformed radii have zero MAD.", call. = FALSE)
  }
  q1 <- .c6rs_q1(transformed.center - transformed.mad,
                  q1_lower_bound)
  values <- c(
    Q1 = q1$value,
    Q2 = q2,
    Q3 = (transformed.center + transformed.mad)^(3 / 2),
    Q3_star = (transformed.center + 1.4826 * transformed.mad)^(3 / 2)
  )
  values <- .c6rs_validate_cutoffs(values)
  list(
    values = values, source = source,
    transformed.center = transformed.center,
    transformed.mad = transformed.mad, h = h,
    zero.mad = is.zero.mad,
    zero.mad.rule = if (is.zero.mad) {
      paste(
        "coincident cutoffs use the closed-inner/zero-outer limiting rule;",
        "no denominator is perturbed"
      )
    } else {
      "not applicable"
    },
    q1.adjusted = q1$adjusted,
    q1.base = q1$base
  )
}


#' Spatial-sign principal component analysis
#'
#' Decomposes the sample spatial-sign covariance matrix (SSCM). When the
#' spatial-median center is certified, the default operator conventions match
#' [sscm()]: zero residuals map to zero and the divisor is `n`. Scores are
#' projections of the centered original observations, not projections of
#' their spatial signs.
#'
#' @param x Numeric matrix or data frame with observations in rows.
#' @param rank Number of loading vectors. `NULL` uses `min(n, p)`.
#' @param center A numeric center or one of `"spatial"`, `"mean"`, and
#'   `"none"`.
#' @param divisor Either `"n"`, the book and Chapter 1 definition, or
#'   `"nonzero"`, which divides by the number of nonzero residual signs.
#' @param zero_action `"zero"` maps residuals no larger than `zero_tol` to
#'   zero; `"error"` rejects such a sample.
#' @param tol,max_iter Spatial-median convergence controls.
#' @param zero_tol Non-negative tolerance for a zero residual. The default
#'   detects exact zeros and preserves scale equivariance.
#' @param keep_operator If `TRUE`, retain the decomposed SSCM.
#'
#' @return An object inheriting from `hd_pca_fit`. `eigenvalues` contains the
#'   full spectrum and `loadings` contains the requested leading vectors.
#'
#' @references
#' Taskinen, S., Kankainen, A., and Oja, H. (2012). Sign covariance matrix
#' estimate with an application to principal components. *Statistics &
#' Probability Letters*, 82, 1153--1161.
#'
#' @examples
#' x <- rbind(c(3, 0), c(-3, 0), c(0, 1), c(0, -1))
#' spatial_sign_pca(x, rank = 1, center = "none")
#'
#' @importFrom stats setNames
#' @export
spatial_sign_pca <- function(
    x, rank = NULL, center = c("spatial", "mean", "none"),
    divisor = c("n", "nonzero"), zero_action = c("zero", "error"),
    tol = 1e-8, max_iter = 500L, zero_tol = 0,
    keep_operator = TRUE) {
  call <- match.call()
  if (missing(center)) center <- "spatial"
  x <- .as_data_matrix(x, min_rows = 2L)
  rank <- .c6rs_rank(rank, nrow(x), ncol(x))
  divisor <- match.arg(divisor)
  zero_action <- match.arg(zero_action)
  keep_operator <- .c6rs_flag(keep_operator, "keep_operator")
  center.fit <- .c6rs_center(
    x, center, tol, max_iter, zero_tol, allow_lts = FALSE
  )
  signs <- spatial_sign(
    x, center = center.fit$value, zero_tol = zero_tol
  )
  n.zero <- as.integer(attr(signs, "n_zero"))
  if (n.zero > 0L && identical(zero_action, "error")) {
    stop(sprintf("The centered sample contains %d zero residual(s).", n.zero),
         call. = FALSE)
  }
  denominator <- if (identical(divisor, "n")) {
    nrow(x)
  } else {
    nrow(x) - n.zero
  }
  if (denominator <= 0) {
    stop("The nonzero-residual divisor is zero.", call. = FALSE)
  }
  operator <- crossprod(unclass(signs)) / denominator
  variable.names <- .c6rs_names(x)
  dimnames(operator) <- list(variable.names, variable.names)
  decomposition <- .c6rs_decompose(operator, rank, variable.names)
  .c6rs_fit(
    method = "spatial-sign PCA", method.class = "spatial_sign_pca_fit",
    x = x, rank = rank, center = center.fit$value,
    operator = operator, decomposition = decomposition,
    diagnostics = list(
      center.source = center.fit$source,
      center = center.fit$diagnostics,
      divisor.rule = divisor,
      divisor.value = denominator,
      zero.action = zero_action,
      zero.tolerance = .validate_zero_tol(zero_tol),
      n.zero.residuals = n.zero,
      operator.trace = sum(diag(operator)),
      score.definition = "centered original observations times loadings"
    ),
    call = call, keep_operator = keep_operator
  )
}


#' Multivariate-Kendall principal component analysis
#'
#' Decomposes the exact multivariate Kendall U-statistic over all unordered
#' pairs. Its default tie convention matches [spatial_kendall()]: a tied
#' difference contributes the zero matrix while the divisor remains
#' `choose(n, 2)`. The operator is translation invariant; `center` controls
#' only the reported scores.
#'
#' @inheritParams spatial_sign_pca
#' @param center A numeric score center or one of `"mean"`, `"spatial"`, and
#'   `"none"`. It does not enter the Kendall operator.
#' @param divisor `"all_pairs"` uses every unordered pair;
#'   `"nonzero_pairs"` excludes tied pairs from the divisor.
#' @param ties `"zero"` gives tied differences zero contribution;
#'   `"error"` rejects any tied pair.
#' @param zero_tol Non-negative tolerance for a tied pairwise difference.
#'
#' @return An object inheriting from `hd_pca_fit`. `eigenvalues` contains the
#'   full spectrum and `loadings` contains the requested leading vectors.
#'
#' @references
#' Han, F. and Liu, H. (2018). ECA: High-dimensional elliptical component
#' analysis in non-Gaussian distributions. *Journal of the American
#' Statistical Association*, 113, 252--268.
#'
#' @examples
#' x <- rbind(c(2, 0), c(-2, 0), c(0, 1), c(0, -1))
#' kendall_pca(x, rank = 1)
#'
#' @export
kendall_pca <- function(
    x, rank = NULL, center = c("mean", "spatial", "none"),
    divisor = c("all_pairs", "nonzero_pairs"),
    ties = c("zero", "error"), tol = 1e-8, max_iter = 500L,
    zero_tol = 0, keep_operator = TRUE) {
  call <- match.call()
  if (missing(center)) center <- "mean"
  x <- .as_data_matrix(x, min_rows = 2L)
  rank <- .c6rs_rank(rank, nrow(x), ncol(x))
  divisor <- match.arg(divisor)
  ties <- match.arg(ties)
  keep_operator <- .c6rs_flag(keep_operator, "keep_operator")
  zero_tol <- .validate_zero_tol(zero_tol)
  kendall <- spatial_kendall(x, zero_tol = zero_tol)
  n.pairs <- as.integer(attr(kendall, "n_pairs"))
  n.zero <- as.integer(attr(kendall, "n_zero_pairs"))
  if (n.zero > 0L && identical(ties, "error")) {
    stop(sprintf("The sample contains %d tied pairwise difference(s).", n.zero),
         call. = FALSE)
  }
  denominator <- if (identical(divisor, "all_pairs")) {
    n.pairs
  } else {
    n.pairs - n.zero
  }
  if (denominator <= 0L) {
    stop("The nonzero-pair divisor is zero.", call. = FALSE)
  }
  operator <- matrix(as.numeric(kendall), ncol(x), ncol(x))
  if (identical(divisor, "nonzero_pairs")) {
    operator <- operator * (n.pairs / denominator)
  }
  variable.names <- .c6rs_names(x)
  dimnames(operator) <- list(variable.names, variable.names)
  center.fit <- .c6rs_center(
    x, center, tol, max_iter, zero_tol, allow_lts = FALSE
  )
  decomposition <- .c6rs_decompose(operator, rank, variable.names)
  .c6rs_fit(
    method = "multivariate-Kendall ECA", method.class = "kendall_pca_fit",
    x = x, rank = rank, center = center.fit$value,
    operator = operator, decomposition = decomposition,
    diagnostics = list(
      operator.center = "none: pairwise differences are translation invariant",
      score.center.source = center.fit$source,
      score.center = center.fit$diagnostics,
      divisor.rule = divisor,
      divisor.value = denominator,
      n.pairs = n.pairs,
      n.zero.pairs = n.zero,
      ties.policy = ties,
      zero.tolerance = zero_tol,
      operator.trace = sum(diag(operator)),
      score.definition = "centered original observations times loadings"
    ),
    call = call, keep_operator = keep_operator
  )
}


#' Generalized spatial-sign principal component analysis
#'
#' Forms the generalized spatial-sign covariance matrix
#' \deqn{n^{-1}\sum_i g(x_i-T)g(x_i-T)',\qquad
#' g(t)=\xi(\lVert t\rVert_2)t,}
#' and decomposes it. Available radial multipliers are `"winsor"`,
#' `"quadratic"`, `"ball"`, `"shell"`, `"linear_redescending"`,
#' `"spatial"` (ordinary spatial signs), and `"identity"` (the centered
#' second moment).
#'
#' The default cutoff convention is the ordinary median/raw-MAD construction
#' stated in Chapter 6 of the book. Because the cited 2024 manuscript was not
#' available for formula verification, `"median_mad"` is documented here as a
#' book-defined convention and is not claimed to reproduce that manuscript.
#' `"original_h_order"` implements the 2019 order statistic with
#' \eqn{h=\lfloor(n+p+1)/2\rfloor} and fails when `h > n`.
#' `"user"` requires all four named cutoffs. Ball includes `Q2`; Shell includes
#' both `Q1` and `Q3`. When the MAD is zero, the default exact limiting rule
#' uses coincident cutoffs and never perturbs a denominator.
#' The `"spatial"` and `"identity"` families use no cutoffs and reject an
#' explicitly supplied `cutoff` or `cutoffs` argument.
#'
#' For a general elliptical representation
#' \eqn{X-\mu=RAU}, \eqn{U=Z/\lVert Z\rVert}, and
#' \eqn{S=\sum_l\lambda_l Z_l^2}, the population generalized eigenvalue
#' contains the radial variable:
#' \deqn{E\{K^2(|R|\sqrt{S}/\lVert Z\rVert)
#' \lambda_j Z_j^2/S\}.}
#' The simpler expression without `R` is therefore not claimed outside the
#' Gaussian representation. Eigenvector and ordering results require the
#' conditions stated in Chapter 6 for this book-defined construction, not
#' merely a measurable weight.
#'
#' @inheritParams spatial_sign_pca
#' @param weight Radial-multiplier family.
#' @param center Numeric center or one of `"kstep_lts"`, `"spatial"`,
#'   `"mean"`, and `"none"`. K-step LTS starts from the certified Chapter 1
#'   spatial median.
#' @param cutoff One of the Chapter 6 book-defined `"median_mad"` rule, the
#'   2019 `"original_h_order"` rule, or `"user"`.
#' @param cutoffs For `cutoff = "user"`, four non-negative values named
#'   exactly `Q1`, `Q2`, `Q3`, and `Q3_star`; input order is arbitrary.
#' @param q1_lower_bound `"zero"` explicitly applies the software lower bound
#'   \eqn{\max(0,m-s)} before the power `3/2`; `"error"` rejects a negative
#'   Wilson--Hilferty base.
#' @param zero_mad `"limit"` uses the exact coincident-cutoff limit;
#'   `"error"` rejects zero transformed MAD.
#' @param lts_steps Number of deterministic LTS C-steps when
#'   `center = "kstep_lts"`.
#' @param divisor `"n"` uses the Chapter 6 book definition; `"nonzero"`
#'   excludes only zero residuals from the divisor. Observations rejected by a
#'   radial weight remain part of either divisor.
#'
#' @return An object inheriting from `hd_pca_fit`. `eigenvalues` contains the
#'   full spectrum and `loadings` contains the requested leading vectors.
#'
#' @references
#' Raymaekers, J. and Rousseeuw, P. J. (2019). A generalized spatial sign
#' covariance matrix. *Journal of Multivariate Analysis*, 171, 94--111.
#'
#' Chapter 6 cites Leyder, S., Raymaekers, J., and Verdonck, T. (2024),
#' Generalized spherical principal component analysis, *Statistics and
#' Computing*, 34, 104. That manuscript was unavailable for formula
#' verification; the `"median_mad"` construction above is therefore attributed
#' to the book rather than asserted as a reproduction of the cited article.
#'
#' @examples
#' x <- rbind(c(3, 0), c(-3, 0), c(0, 1), c(0, -1))
#' generalized_sign_pca(
#'   x, rank = 1, weight = "winsor", center = "none"
#' )
#'
#' @export
generalized_sign_pca <- function(
    x, rank = NULL,
    weight = c("winsor", "quadratic", "ball", "shell",
               "linear_redescending", "spatial", "identity"),
    center = c("kstep_lts", "spatial", "mean", "none"),
    cutoff = c("median_mad", "original_h_order", "user"),
    cutoffs = NULL, q1_lower_bound = c("zero", "error"),
    zero_mad = c("limit", "error"), lts_steps = 2L,
    divisor = c("n", "nonzero"), zero_action = c("zero", "error"),
    tol = 1e-8, max_iter = 500L, zero_tol = 0,
    keep_operator = TRUE) {
  call <- match.call()
  cutoff.supplied <- !missing(cutoff)
  if (missing(center)) center <- "kstep_lts"
  x <- .as_data_matrix(x, min_rows = 2L)
  rank <- .c6rs_rank(rank, nrow(x), ncol(x))
  weight <- match.arg(weight)
  cutoff <- match.arg(cutoff)
  uses.cutoff <- weight %in% c(
    "winsor", "quadratic", "ball", "shell", "linear_redescending"
  )
  if (!uses.cutoff && (cutoff.supplied || !is.null(cutoffs))) {
    stop("`cutoff` and `cutoffs` are not used by this radial family.",
         call. = FALSE)
  }
  q1_lower_bound <- match.arg(q1_lower_bound)
  zero_mad <- match.arg(zero_mad)
  divisor <- match.arg(divisor)
  zero_action <- match.arg(zero_action)
  keep_operator <- .c6rs_flag(keep_operator, "keep_operator")
  center.fit <- .c6rs_center(
    x, center, tol, max_iter, zero_tol,
    allow_lts = TRUE, lts_steps = lts_steps
  )
  scaled <- .center_and_scale(x, center.fit$value, zero_tol)
  sign.kernel <- cpp_spatial_sign(scaled$x, scaled$zero_tol)
  scaled.radii <- as.numeric(sign.kernel$norms)
  n.zero <- as.integer(sign.kernel$n_zero)
  if (n.zero > 0L && identical(zero_action, "error")) {
    stop(sprintf("The centered sample contains %d zero residual(s).", n.zero),
         call. = FALSE)
  }
  denominator <- if (identical(divisor, "n")) {
    nrow(x)
  } else {
    nrow(x) - n.zero
  }
  if (denominator <= 0L) {
    stop("The nonzero-residual divisor is zero.", call. = FALSE)
  }

  if (uses.cutoff) {
    scaled.user <- if (is.null(cutoffs)) NULL else
      .c6rs_validate_cutoffs(cutoffs) / scaled$scale
    cutoff.fit <- .c6rs_cutoffs(
      scaled.radii, nrow(x), ncol(x), cutoff, scaled.user,
      q1_lower_bound, zero_mad
    )
    scaled.cutoffs <- cutoff.fit$values
  } else {
    scaled.cutoffs <- setNames(rep(0, 4L),
                               c("Q1", "Q2", "Q3", "Q3_star"))
    cutoff.fit <- list(
      values = scaled.cutoffs, source = "not used by this radial family",
      transformed.center = NA_real_, transformed.mad = NA_real_,
      h = NA_integer_, zero.mad = NA, zero.mad.rule = "not applicable",
      q1.adjusted = FALSE, q1.base = NA_real_
    )
  }
  kernel <- cpp_ch6rs_radial_transform(
    scaled$x, weight, as.numeric(scaled.cutoffs), scaled$zero_tol
  )
  if (as.integer(kernel$n_zero) != n.zero) {
    stop("Internal zero-residual counts disagree.", call. = FALSE)
  }
  original.radii <- scaled.radii * scaled$scale
  if (identical(weight, "spatial")) {
    operator <- kernel$matrix_sum / denominator
    transformed.norms <- as.numeric(kernel$transformed_norms)
    radial.multiplier <- numeric(length(original.radii))
    retained <- transformed.norms > 0
    radial.multiplier[retained] <- 1 / original.radii[retained]
  } else {
    transformed <- kernel$transformed * scaled$scale
    if (any(!is.finite(transformed))) {
      stop("The original-scale radial transform is non-finite.",
           call. = FALSE)
    }
    operator <- crossprod(transformed) / denominator
    radial.multiplier <- as.numeric(kernel$radial_multiplier)
    transformed.norms <- as.numeric(kernel$transformed_norms) * scaled$scale
  }
  if (any(!is.finite(operator))) {
    stop("The generalized spatial-sign operator is non-finite.",
         call. = FALSE)
  }
  original.cutoffs <- scaled.cutoffs * scaled$scale
  transformed.scale <- scaled$scale^(2 / 3)
  original.transformed.center <-
    cutoff.fit$transformed.center * transformed.scale
  original.transformed.mad <- cutoff.fit$transformed.mad * transformed.scale
  original.q1.base <- cutoff.fit$q1.base * transformed.scale
  variable.names <- .c6rs_names(x)
  dimnames(operator) <- list(variable.names, variable.names)
  decomposition <- .c6rs_decompose(operator, rank, variable.names)
  .c6rs_fit(
    method = sprintf("generalized spatial-sign PCA (%s)", weight),
    method.class = "generalized_sign_pca_fit",
    x = x, rank = rank, center = center.fit$value,
    operator = operator, decomposition = decomposition,
    diagnostics = list(
      center.source = center.fit$source,
      center = center.fit$diagnostics,
      radial.family = weight,
      radial.definition = "g(t) = xi(||t||) t",
      cutoffs = if (uses.cutoff) original.cutoffs else NULL,
      cutoff.source = cutoff.fit$source,
      cutoff.method = if (uses.cutoff) cutoff else "none",
      transformed.radius.center = original.transformed.center,
      transformed.radius.mad = original.transformed.mad,
      original.h = cutoff.fit$h,
      q1.lower.bound.policy = q1_lower_bound,
      q1.lower.bound.adjusted = cutoff.fit$q1.adjusted,
      q1.transformed.base = original.q1.base,
      zero.mad = cutoff.fit$zero.mad,
      zero.mad.policy = zero_mad,
      zero.mad.rule = cutoff.fit$zero.mad.rule,
      boundary.rule = paste(
        "Ball includes Q2; Shell includes Q1 and Q3;",
        "linear redescending includes its inner endpoint"
      ),
      divisor.rule = divisor,
      divisor.value = denominator,
      zero.action = zero_action,
      zero.tolerance = .validate_zero_tol(zero_tol),
      n.zero.residuals = n.zero,
      n.nonzero.transformed = as.integer(kernel$n_nonzero_transform),
      radii = original.radii,
      radial.multiplier = radial.multiplier,
      transformed.norms = transformed.norms,
      operator.trace = sum(diag(operator)),
      score.definition = "centered original observations times loadings",
      population.scope = paste(
        "general elliptical population formulas retain the independent radial",
        "variable; the Gaussian-only radial simplification is not asserted"
      )
    ),
    call = call, keep_operator = keep_operator
  )
}
