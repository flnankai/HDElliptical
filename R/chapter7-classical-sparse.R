# Chapter 7: classical and sparse prototype clustering.
#
# All initialization, random-number, tie, empty-cluster, regularization, and
# stopping choices are public arguments.  The strict defaults never repair an
# empty cluster and never regularize a Gaussian covariance.

.ch7cs_flag <- function(x, name) {
  if (!is.logical(x) || length(x) != 1L || is.na(x)) {
    stop(sprintf("`%s` must be TRUE or FALSE.", name), call. = FALSE)
  }
  x
}


.ch7cs_count <- function(x, name, minimum = 1L, maximum = Inf) {
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


.ch7cs_scalar <- function(x, name, allow_zero = FALSE) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) || !is.finite(x) ||
      if (allow_zero) x < 0 else x <= 0) {
    qualifier <- if (allow_zero) "non-negative" else "strictly positive"
    stop(sprintf("`%s` must be one finite %s number.", name, qualifier),
         call. = FALSE)
  }
  as.numeric(x)
}


.ch7cs_data <- function(x, name = "x") {
  if (is.data.frame(x)) x <- data.matrix(x)
  x <- as.matrix(x)
  storage.mode(x) <- "double"
  if (!is.numeric(x) || nrow(x) < 2L || ncol(x) < 1L || anyNA(x) ||
      any(!is.finite(x))) {
    stop(sprintf(
      "`%s` must be a finite numeric matrix with at least two rows and one column.",
      name
    ), call. = FALSE)
  }
  variable.names <- colnames(x)
  if (is.null(variable.names)) variable.names <- paste0("V", seq_len(ncol(x)))
  if (anyNA(variable.names) || any(!nzchar(variable.names)) ||
      anyDuplicated(variable.names)) {
    stop(sprintf("Column names of `%s` must be non-empty and unique.", name),
         call. = FALSE)
  }
  colnames(x) <- variable.names
  if (is.null(rownames(x))) rownames(x) <- paste0("row", seq_len(nrow(x)))
  x
}


.ch7cs_seed <- function(seed, name = "seed") {
  .ch7cs_count(seed, name, minimum = 0L, maximum = .Machine$integer.max)
}


.ch7cs_with_seed <- function(seed, expression) {
  seed <- .ch7cs_seed(seed)
  had.seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had.seed) old.seed <- get(".Random.seed", envir = .GlobalEnv,
                                inherits = FALSE)
  on.exit({
    if (had.seed) {
      assign(".Random.seed", old.seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  set.seed(seed)
  force(expression)
}


.ch7cs_ties <- function(ties) {
  if (!is.character(ties) || length(ties) != 1L || is.na(ties) ||
      ties != "first") {
    stop("`ties` must be \"first\"; exact ties use the smallest index.",
         call. = FALSE)
  }
  ties
}


.ch7cs_empty_action <- function(empty_action) {
  match.arg(empty_action, c("error", "farthest"))
}


.ch7cs_initial <- function(x, clusters, initial, initialization, first_index,
                           seed, metric = c("squared", "manhattan")) {
  metric <- match.arg(metric)
  n <- nrow(x)
  p <- ncol(x)
  clusters <- .ch7cs_count(clusters, "clusters", 1L, n)
  initialization <- match.arg(initialization, c("maxmin", "random"))

  if (!is.null(initial)) {
    if (!is.null(seed)) {
      stop("`seed` must be NULL when `initial` is supplied.", call. = FALSE)
    }
    if (is.matrix(initial) || is.data.frame(initial)) {
      centers <- as.matrix(initial)
      storage.mode(centers) <- "double"
      if (!identical(dim(centers), c(clusters, p)) || anyNA(centers) ||
          any(!is.finite(centers))) {
        stop(sprintf(
          "A matrix `initial` must be a finite %d by %d center matrix.",
          clusters, p
        ), call. = FALSE)
      }
      colnames(centers) <- colnames(x)
      return(list(
        centers = centers, indices = NULL, source = "supplied centers",
        initialization = "supplied"
      ))
    }
    indices <- as.numeric(initial)
    if (length(indices) != clusters || anyNA(indices) ||
        any(!is.finite(indices)) || any(indices != floor(indices)) ||
        any(indices < 1) || any(indices > n) || anyDuplicated(indices)) {
      stop(sprintf(
        "A vector `initial` must contain %d distinct row indices in 1, ..., %d.",
        clusters, n
      ), call. = FALSE)
    }
    indices <- as.integer(indices)
    return(list(
      centers = x[indices, , drop = FALSE], indices = indices,
      source = "supplied row indices", initialization = "supplied"
    ))
  }

  if (initialization == "random") {
    if (is.null(seed)) {
      stop("Random initialization requires an explicit integer `seed`.",
           call. = FALSE)
    }
    indices <- .ch7cs_with_seed(seed, sample.int(n, clusters,
                                                 replace = FALSE))
    return(list(
      centers = x[indices, , drop = FALSE], indices = as.integer(indices),
      source = "seeded random row indices", initialization = "random"
    ))
  }
  if (!is.null(seed)) {
    stop("`seed` must be NULL for deterministic max-min initialization.",
         call. = FALSE)
  }
  first_index <- .ch7cs_count(first_index, "first_index", 1L, n)
  if (sum(!duplicated(x)) < clusters) {
    stop("Max-min initialization requires at least `clusters` unique rows.",
         call. = FALSE)
  }
  indices <- integer(clusters)
  indices[1L] <- first_index
  minimum.distance <- rep(Inf, n)
  for (k in seq_len(clusters)) {
    if (k > 1L) {
      maximum <- max(minimum.distance)
      candidates <- which(minimum.distance == maximum)
      indices[k] <- candidates[1L]
    }
    difference <- sweep(x, 2L, x[indices[k], ], "-", check.margin = FALSE)
    distance <- if (metric == "squared") rowSums(difference^2) else
      rowSums(abs(difference))
    minimum.distance <- pmin(minimum.distance, distance)
    minimum.distance[indices[seq_len(k)]] <- -Inf
  }
  list(
    centers = x[indices, , drop = FALSE], indices = indices,
    source = paste0("deterministic max-min (", metric, ")"),
    initialization = "maxmin"
  )
}


.ch7cs_bound <- function(s, p, name) {
  s <- .ch7cs_scalar(s, name)
  if (s < 1 || s > sqrt(p)) {
    stop(sprintf("`%s` must lie in [1, sqrt(p)].", name), call. = FALSE)
  }
  s
}


.ch7cs_strict <- function(core, strict, method) {
  strict <- .ch7cs_flag(strict, "strict")
  if (!isTRUE(core$certified)) {
    message <- sprintf("%s stopped without its solver certificate.", method)
    if (strict) stop(message, call. = FALSE)
    warning(message, call. = FALSE)
  }
  strict
}


.ch7cs_names <- function(core, x) {
  rownames(core$centers) <- paste0("cluster", seq_len(nrow(core$centers)))
  colnames(core$centers) <- colnames(x)
  names(core$cluster) <- rownames(x)
  names(core$size) <- rownames(core$centers)
  if (!is.null(core$distances)) {
    dimnames(core$distances) <- list(rownames(x), rownames(core$centers))
  }
  core
}


#' Lloyd's deterministic K-means algorithm
#'
#' Minimize the Chapter 7 within-cluster sum of squared Euclidean distances by
#' alternating exact sample-mean and nearest-center updates. Exact assignment
#' ties use the smallest cluster index. The default `empty_action = "error"`
#' applies no repair. The optional `"farthest"` rule moves the smallest-row
#' maximizer of current within-cluster squared distance from a donor of size at
#' least two, then recomputes centers.
#'
#' `initial` may be a `clusters` by `p` center matrix or `clusters` distinct
#' row indices. Otherwise initialization is either deterministic max-min from
#' the explicit `first_index`, or seeded sampling without replacement. Seeded
#' initialization restores the caller's RNG state.
#'
#' @param x Numeric observation matrix, with observations in rows.
#' @param clusters Number of clusters.
#' @param initial Optional center matrix or vector of distinct row indices.
#' @param initialization Max-min or seeded random initialization when `initial`
#'   is `NULL`.
#' @param first_index First max-min seed row.
#' @param seed Explicit seed required for random initialization; otherwise
#'   `NULL`.
#' @param ties Exact tie rule; only deterministic `"first"` is implemented.
#' @param empty_action Error, or the explicit farthest-observation repair.
#' @param solver_tol Tolerance used only to certify objective monotonicity.
#' @param solver_max_iter Maximum Lloyd iterations.
#' @param strict Error rather than return an uncertified iterate.
#' @param keep_distances Retain the final squared-distance matrix.
#'
#' @return A `lloyd_kmeans_fit` object with labels, centers, sizes, objective
#'   trace, initialization record, and convergence/repair certificates.
#'
#' @references
#' Feng, L. (2026). *High-Dimensional Data Analysis for Elliptical Symmetric
#' Distributions*, Chapter 7 (book manuscript). This implementation follows
#' the manuscript's explicit deterministic Lloyd-update contract.
#'
#' For a foundational k-means formulation, see MacQueen, J. B. (1967). Some
#' methods for classification and analysis of multivariate observations.
#' *Proceedings of the Fifth Berkeley Symposium*, 1, 281--297.
#' @export
#' @examples
#' x <- rbind(c(-2, 0), c(-1, 0), c(2, 0), c(3, 0))
#' lloyd_kmeans(x, clusters = 2, initial = c(1, 3))
lloyd_kmeans <- function(
    x, clusters, initial = NULL,
    initialization = c("maxmin", "random"), first_index = 1L, seed = NULL,
    ties = "first", empty_action = c("error", "farthest"),
    solver_tol = 1e-10, solver_max_iter = 100L, strict = TRUE,
    keep_distances = FALSE) {
  call <- match.call()
  x <- .ch7cs_data(x)
  clusters <- .ch7cs_count(clusters, "clusters", 1L, nrow(x))
  ties <- .ch7cs_ties(ties)
  empty_action <- .ch7cs_empty_action(empty_action)
  solver_tol <- .ch7cs_scalar(solver_tol, "solver_tol")
  solver_max_iter <- .ch7cs_count(solver_max_iter, "solver_max_iter")
  keep_distances <- .ch7cs_flag(keep_distances, "keep_distances")
  initialization <- match.arg(initialization)
  start <- .ch7cs_initial(
    x, clusters, initial, initialization, first_index, seed, "squared"
  )
  core <- cpp_ch7cs_lloyd_core(
    x, start$centers, match(empty_action, c("error", "farthest")) - 1L,
    solver_tol, solver_max_iter, keep_distances
  )
  core$cluster <- as.integer(core$cluster)
  core$size <- as.integer(core$size)
  core <- .ch7cs_names(core, x)
  .ch7cs_strict(core, strict, "Lloyd K-means")
  answer <- list(
    call = call, method = "Lloyd K-means", cluster = core$cluster,
    centers = core$centers, size = core$size, objective = core$objective,
    objective_trace = core$objective_trace, iterations = core$iterations,
    converged = isTRUE(core$converged), valid = isTRUE(core$certified),
    initialization = start, distances = core$distances,
    diagnostics = core[c(
      "cycle_detected", "assignment_ties", "empty_repairs", "repair_ties",
      "maximum_objective_increase", "assignment_residual", "certified"
    )]
  )
  answer$initialization$centers <- start$centers
  class(answer) <- c("lloyd_kmeans_fit", "hd_clustering_fit")
  answer
}


.ch7cs_covariances <- function(covariances, p, clusters, symmetry_tol) {
  symmetry_tol <- .ch7cs_scalar(symmetry_tol, "symmetry_tol", allow_zero = TRUE)
  if (is.matrix(covariances) && clusters == 1L) {
    covariances <- array(covariances, dim = c(p, p, 1L))
  }
  if (!is.array(covariances) || length(dim(covariances)) != 3L ||
      !identical(dim(covariances), c(p, p, clusters))) {
    stop(sprintf("`covariances` must have dimensions %d by %d by %d.",
                 p, p, clusters), call. = FALSE)
  }
  storage.mode(covariances) <- "double"
  if (anyNA(covariances) || any(!is.finite(covariances))) {
    stop("`covariances` must be finite.", call. = FALSE)
  }
  errors <- numeric(clusters)
  for (k in seq_len(clusters)) {
    covariance <- covariances[, , k]
    scale <- max(1, max(abs(covariance)))
    errors[k] <- max(abs(covariance - t(covariance)))
    if (errors[k] > symmetry_tol * scale) {
      stop(sprintf(
        "Covariance %d is not symmetric within `symmetry_tol`.", k
      ), call. = FALSE)
    }
    # Arithmetic symmetrization is controlled by the explicit tolerance and
    # recorded; it is not an eigenvalue repair.
    covariances[, , k] <- (covariance + t(covariance)) / 2
  }
  list(value = covariances, errors = errors, tolerance = symmetry_tol)
}


#' Full-covariance Gaussian-mixture EM
#'
#' Implement the Chapter 7 component-specific full-covariance Gaussian EM
#' updates. Means, covariance arrays, and mixing proportions are all required:
#' this function has no hidden start. The E-step uses Cholesky solves and
#' log-sum-exp; no inverse or pseudo-inverse is formed. The M-step covariance
#' has maximum-likelihood divisor `sum(tau[, k])`.
#'
#' The strict default `ridge = 0` implements ordinary EM. A positive explicit
#' ridge adds `ridge * I` after every raw covariance M-step and is labelled
#' regularized EM; unpenalized likelihood monotonicity is then not claimed.
#' Singular initial or updated covariance matrices error instead of receiving
#' an eigenvalue floor. Gaussian-mixture likelihood is unbounded without
#' further restrictions, so a numerical EM fixed point is not a global-MLE
#' certificate.
#'
#' @param x Numeric observation matrix.
#' @param means Required `K` by `p` initial mean matrix.
#' @param covariances Required `p` by `p` by `K` initial covariance array (a
#'   matrix is accepted for `K = 1`).
#' @param proportions Required strictly positive length-`K` vector summing
#'   exactly to one.
#' @param ridge Explicit non-negative covariance ridge; zero means none.
#' @param solver_tol Relative parameter and log-likelihood fixed-point
#'   tolerance.
#' @param monotone_tol Numerical monotonicity tolerance for ordinary EM.
#' @param rank_tol Relative eigenvalue cutoff used only for raw-rank
#'   diagnostics, never for inversion or flooring.
#' @param symmetry_tol Explicit tolerance for arithmetic symmetrization of
#'   supplied covariance matrices.
#' @param solver_max_iter Maximum EM updates.
#' @param strict Error rather than return an uncertified iterate.
#' @param keep_responsibilities Retain the final posterior matrix.
#'
#' @return A `gaussian_mixture_em_fit` object containing fitted parameters,
#'   hard labels (posterior ties use the smallest component), likelihood trace,
#'   effective masses, raw covariance ranks, and solver certificates.
#'
#' @references
#' Feng, L. (2026). *High-Dimensional Data Analysis for Elliptical Symmetric
#' Distributions*, Chapter 7 (book manuscript). The implementation is the
#' explicitly documented standard Gaussian-mixture E/M construction.
#' @export
#' @examples
#' x <- matrix(c(-2.2, -2, -1.8, 1.8, 2, 2.2), ncol = 1)
#' gaussian_mixture_em(
#'   x, means = matrix(c(-2, 2), ncol = 1),
#'   covariances = array(c(0.2, 0.2), c(1, 1, 2)),
#'   proportions = c(0.5, 0.5), solver_tol = 1e-6
#' )
gaussian_mixture_em <- function(
    x, means, covariances, proportions, ridge = 0,
    solver_tol = 1e-8, monotone_tol = 1e-10, rank_tol = 1e-10,
    symmetry_tol = 1e-12, solver_max_iter = 200L, strict = TRUE,
    keep_responsibilities = TRUE) {
  call <- match.call()
  x <- .ch7cs_data(x)
  means <- as.matrix(means)
  storage.mode(means) <- "double"
  if (!is.numeric(means) || nrow(means) < 1L || nrow(means) > nrow(x) ||
      ncol(means) != ncol(x) || anyNA(means) || any(!is.finite(means))) {
    stop("`means` must be a finite K by p matrix with 1 <= K <= n.",
         call. = FALSE)
  }
  clusters <- nrow(means)
  colnames(means) <- colnames(x)
  rownames(means) <- paste0("cluster", seq_len(clusters))
  proportions <- as.numeric(proportions)
  if (length(proportions) != clusters || anyNA(proportions) ||
      any(!is.finite(proportions)) || any(proportions <= 0)) {
    stop("`proportions` must contain K finite strictly positive values.",
         call. = FALSE)
  }
  if (sum(proportions) != 1) {
    stop("`proportions` must sum exactly to one; normalize them explicitly.",
         call. = FALSE)
  }
  ridge <- .ch7cs_scalar(ridge, "ridge", allow_zero = TRUE)
  solver_tol <- .ch7cs_scalar(solver_tol, "solver_tol")
  monotone_tol <- .ch7cs_scalar(monotone_tol, "monotone_tol",
                                allow_zero = TRUE)
  rank_tol <- .ch7cs_scalar(rank_tol, "rank_tol")
  solver_max_iter <- .ch7cs_count(solver_max_iter, "solver_max_iter")
  keep_responsibilities <- .ch7cs_flag(
    keep_responsibilities, "keep_responsibilities"
  )
  covariance.input <- .ch7cs_covariances(
    covariances, ncol(x), clusters, symmetry_tol
  )
  core <- cpp_ch7cs_gmm_core(
    x, means, covariance.input$value, proportions, ridge, solver_tol,
    monotone_tol, rank_tol, solver_max_iter, keep_responsibilities
  )
  core$cluster <- as.integer(core$cluster)
  core$proportions <- as.numeric(core$proportions)
  core$effective_masses <- as.numeric(core$effective_masses)
  core$raw_covariance_ranks <- as.integer(core$raw_covariance_ranks)
  core$minimum_eigenvalues <- as.numeric(core$minimum_eigenvalues)
  names(core$cluster) <- rownames(x)
  dimnames(core$means) <- list(rownames(means), colnames(x))
  dimnames(core$covariances) <- list(
    colnames(x), colnames(x), rownames(means)
  )
  names(core$proportions) <- rownames(means)
  names(core$effective_masses) <- rownames(means)
  names(core$raw_covariance_ranks) <- rownames(means)
  names(core$minimum_eigenvalues) <- rownames(means)
  if (!is.null(core$responsibilities)) {
    dimnames(core$responsibilities) <- list(rownames(x), rownames(means))
  }
  .ch7cs_strict(core, strict, "Gaussian-mixture EM")
  answer <- list(
    call = call,
    method = if (ridge == 0) "full-covariance Gaussian-mixture EM" else
      "ridge-regularized full-covariance Gaussian-mixture EM",
    cluster = core$cluster, means = core$means,
    covariances = core$covariances, proportions = core$proportions,
    responsibilities = core$responsibilities,
    log_likelihood = core$log_likelihood,
    log_likelihood_trace = core$log_likelihood_trace,
    iterations = core$iterations, converged = isTRUE(core$converged),
    valid = isTRUE(core$certified),
    diagnostics = c(core[c(
      "effective_masses", "raw_covariance_ranks", "minimum_eigenvalues",
      "responsibility_row_sum_residual", "fixed_point_residual",
      "maximum_log_likelihood_decrease", "monotone", "label_ties",
      "certified"
    )], list(
      ridge = ridge, regularized = ridge > 0,
      supplied_covariance_symmetry_errors = covariance.input$errors,
      supplied_covariance_symmetry_tolerance = covariance.input$tolerance,
      covariance_divisor = "effective mass sum(tau[, k])",
      inversion = "Cholesky triangular solves; no inverse or pseudo-inverse",
      likelihood_warning = "finite Gaussian-mixture likelihood is unbounded"
    ))
  )
  class(answer) <- c("gaussian_mixture_em_fit", "hd_clustering_fit")
  answer
}


.ch7cs_sparse_fit <- function(
    x, clusters, bound, initial, initialization, first_index, seed, ties,
    empty_action, weight_tol, weight_max_iter, solver_tol, solver_max_iter,
    strict, keep_distances, median = FALSE, score_tol = NULL) {
  x <- .ch7cs_data(x)
  clusters <- .ch7cs_count(clusters, "clusters", 1L, nrow(x))
  bound.name <- if (median) "s_w" else "s"
  bound <- .ch7cs_bound(bound, ncol(x), bound.name)
  ties <- .ch7cs_ties(ties)
  empty_action <- .ch7cs_empty_action(empty_action)
  weight_tol <- .ch7cs_scalar(weight_tol, "weight_tol")
  weight_max_iter <- .ch7cs_count(weight_max_iter, "weight_max_iter")
  solver_tol <- .ch7cs_scalar(solver_tol, "solver_tol")
  solver_max_iter <- .ch7cs_count(solver_max_iter, "solver_max_iter")
  keep_distances <- .ch7cs_flag(keep_distances, "keep_distances")
  initialization <- match.arg(initialization, c("maxmin", "random"))
  start <- .ch7cs_initial(
    x, clusters, initial, initialization, first_index, seed,
    if (median) "manhattan" else "squared"
  )
  empty.code <- match(empty_action, c("error", "farthest")) - 1L
  if (median) {
    score_tol <- .ch7cs_scalar(score_tol, "score_tol")
    core <- cpp_ch7cs_sparse_kmedian_core(
      x, start$centers, bound, empty.code, score_tol, weight_tol,
      weight_max_iter, solver_tol, solver_max_iter, keep_distances
    )
  } else {
    core <- cpp_ch7cs_sparse_kmeans_core(
      x, start$centers, bound, empty.code, weight_tol, weight_max_iter,
      solver_tol, solver_max_iter, keep_distances
    )
  }
  core$cluster <- as.integer(core$cluster)
  core$size <- as.integer(core$size)
  core <- .ch7cs_names(core, x)
  core$weights <- as.numeric(core$weights)
  core$feature_scores <- as.numeric(core$feature_scores)
  names(core$weights) <- colnames(x)
  names(core$feature_scores) <- colnames(x)
  if (!is.null(core$half_scaled_feature_scores)) {
    core$half_scaled_feature_scores <-
      as.numeric(core$half_scaled_feature_scores)
    names(core$half_scaled_feature_scores) <- colnames(x)
  }
  if (!is.null(core$global_median)) {
    core$global_median <- as.numeric(core$global_median)
    names(core$global_median) <- colnames(x)
  }
  method <- if (median) "sparse K-median" else "sparse K-means"
  .ch7cs_strict(core, strict, method)
  list(x = x, core = core, start = start, bound = bound,
       method = method, empty_action = empty_action)
}


#' Sparse K-means with ordered-pair BCSS weights
#'
#' Optimize the Witten--Tibshirani sparse K-means criterion from Chapter 7.
#' For the returned partition, `feature_scores[j]` is exactly the displayed
#' ordered-pair score
#' `sum_ii' (x_ij-x_i'j)^2/n - sum_k sum_(i,i' in Ck)
#' (x_ij-x_i'j)^2/n_k`, hence it equals twice `TSS[j] - WCSS[j]`.
#' The half-scaled values are returned separately in diagnostics.
#'
#' The weight block is the normalized positive soft-threshold of these scores
#' under `norm(w, 2) <= 1`, `sum(w) <= s`, and `w >= 0`. At `s = 1`, a tied
#' maximum score selects the smallest feature index. All-zero scores error
#' rather than create uniform weights. Weighted assignment uses squared
#' distance after multiplying each column by `sqrt(w)`.
#'
#' @inheritParams lloyd_kmeans
#' @param s Required L1 weight bound in `[1, sqrt(p)]`.
#' @param weight_tol Weight-threshold root and KKT tolerance.
#' @param weight_max_iter Maximum bisection iterations for a weight update.
#'
#' @return A `sparse_kmeans_fit` object with feature scores, weights, labels,
#'   centers, weighted objective, and solver certificates.
#'
#' @references
#' Witten, D. M. and Tibshirani, R. (2010). A framework for feature selection
#' in clustering. *Journal of the American Statistical Association*, 105,
#' 713--726. \doi{10.1198/jasa.2010.tm09415}.
#' @export
#' @examples
#' x <- rbind(c(-3, 0), c(-2, 0), c(2, 0), c(3, 0))
#' sparse_kmeans(x, clusters = 2, s = 1, initial = c(1, 3))
sparse_kmeans <- function(
    x, clusters, s, initial = NULL,
    initialization = c("maxmin", "random"), first_index = 1L, seed = NULL,
    ties = "first", empty_action = c("error", "farthest"),
    weight_tol = 1e-10, weight_max_iter = 200L,
    solver_tol = 1e-10, solver_max_iter = 100L, strict = TRUE,
    keep_distances = FALSE) {
  call <- match.call()
  fit <- .ch7cs_sparse_fit(
    x, clusters, s, initial, initialization, first_index, seed, ties,
    empty_action, weight_tol, weight_max_iter, solver_tol, solver_max_iter,
    strict, keep_distances, median = FALSE
  )
  core <- fit$core
  answer <- list(
    call = call, method = fit$method, cluster = core$cluster,
    centers = core$centers, size = core$size, weights = core$weights,
    feature_scores = core$feature_scores, objective = core$objective,
    objective_trace = core$objective_trace, iterations = core$iterations,
    converged = isTRUE(core$converged), valid = isTRUE(core$certified),
    s = fit$bound, initialization = fit$start, distances = core$distances,
    diagnostics = c(core[c(
      "cycle_detected", "assignment_ties", "empty_repairs", "repair_ties",
      "maximum_objective_decrease", "assignment_residual",
      "weight_certificate", "weighted_within_loss", "certified"
    )], list(
      half_scaled_feature_scores = core$half_scaled_feature_scores,
      feature_score_identity =
        "book ordered-pair BCSS = 2 * (TSS - WCSS)",
      assignment_geometry = "squared Euclidean after column scaling by sqrt(w)",
      no_implicit_repair = fit$empty_action == "error"
    ))
  )
  class(answer) <- c("sparse_kmeans_fit", "hd_clustering_fit")
  answer
}


.ch7cs_permutation_array <- function(x, permutation_indices,
                                     n_permutations, seed, selection) {
  n <- nrow(x)
  p <- ncol(x)
  if (!is.null(permutation_indices)) {
    if (!is.null(n_permutations) || !is.null(seed)) {
      stop(paste0(
        "`n_permutations` and `seed` must be NULL when ",
        "`permutation_indices` is supplied."
      ), call. = FALSE)
    }
    indices <- permutation_indices
    if (!is.array(indices) || length(dim(indices)) != 3L ||
        dim(indices)[1L] != n || dim(indices)[2L] != p ||
        dim(indices)[3L] < 1L) {
      stop("`permutation_indices` must be an n by p by B integer array.",
           call. = FALSE)
    }
    if (!is.numeric(indices) || anyNA(indices) || any(!is.finite(indices)) ||
        any(indices != floor(indices))) {
      stop("`permutation_indices` must contain finite integer values.",
           call. = FALSE)
    }
    storage.mode(indices) <- "integer"
    B <- dim(indices)[3L]
  } else {
    if (is.null(n_permutations) || is.null(seed)) {
      stop(paste0(
        "Supply either `permutation_indices`, or both `n_permutations` ",
        "and an explicit `seed`."
      ), call. = FALSE)
    }
    B <- .ch7cs_count(n_permutations, "n_permutations")
    indices <- .ch7cs_with_seed(seed, {
      answer <- array(NA_integer_, dim = c(n, p, B))
      for (b in seq_len(B)) {
        for (j in seq_len(p)) answer[, j, b] <- sample.int(n)
      }
      answer
    })
  }
  if (selection == "one_se" && B < 2L) {
    stop("`selection = \"one_se\"` requires at least two permutations.",
         call. = FALSE)
  }
  target <- seq_len(n)
  for (b in seq_len(B)) {
    for (j in seq_len(p)) {
      index <- indices[, j, b]
      if (anyNA(index) || !identical(sort(as.integer(index)), target)) {
        stop(sprintf(
          "Permutation slice %d, feature %d is not a permutation of 1, ..., n.",
          b, j
        ), call. = FALSE)
      }
    }
  }
  indices
}


#' Select the sparse K-means L1 bound by permutation Gap calibration
#'
#' For every supplied `s_grid` value this function fits sparse K-means and
#' records its optimized weighted BCSS `O(s)`. Each reference data set
#' independently permutes the rows of every feature, preserving marginal
#' scales and tails while breaking joint clustering. The reported criterion is
#' `Gap(s) = log(O_obs(s)) - mean_b(log(O_perm,b(s)))`; every objective must be
#' finite and strictly positive.
#'
#' `selection = "max"` chooses the smallest `s` attaining the largest Gap.
#' `"one_se"` uses `se(s) = sqrt(1 + 1/B) * sd_b(log(O_perm,b(s)))` and chooses
#' the smallest `s` whose Gap is at least `max(Gap) - se(s_max)`, where `s_max`
#' is the smallest Gap maximizer. This fully states the implemented one-SE
#' convention; no paper simulation is reproduced.
#'
#' @inheritParams sparse_kmeans
#' @param s_grid Finite distinct candidate bounds in `[1, sqrt(p)]`.
#' @param permutation_indices Optional explicit `n` by `p` by `B` array; each
#'   feature column in each slice must be a permutation of `1:n`.
#' @param n_permutations Number of generated references when indices are not
#'   supplied.
#' @param permutation_seed Explicit seed for generated permutation indices.
#' @param selection Maximum-Gap or the stated one-SE rule.
#' @param clustering_seed Seed passed only to random clustering initialization,
#'   separate from `permutation_seed`.
#' @param keep_permutations Retain permutation objectives and index array.
#'
#' @return A `sparse_kmeans_selection` object with all observed fits, Gap/SE
#'   table, selected bound and fit, and optional reference details.
#'
#' @references
#' Witten, D. M. and Tibshirani, R. (2010). A framework for feature selection
#' in clustering. *Journal of the American Statistical Association*, 105,
#' 713--726. \doi{10.1198/jasa.2010.tm09415}.
#' @export
#' @examples
#' x <- rbind(c(-3, 0), c(-2, 1), c(2, 0), c(3, 1))
#' sparse_kmeans_select_s(
#'   x, clusters = 2, s_grid = c(1, sqrt(2)), n_permutations = 2,
#'   permutation_seed = 9, initial = c(1, 3), selection = "max"
#' )
sparse_kmeans_select_s <- function(
    x, clusters, s_grid, permutation_indices = NULL,
    n_permutations = NULL, permutation_seed = NULL,
    selection = c("max", "one_se"), initial = NULL,
    initialization = c("maxmin", "random"), first_index = 1L,
    clustering_seed = NULL, ties = "first",
    empty_action = c("error", "farthest"), weight_tol = 1e-10,
    weight_max_iter = 200L, solver_tol = 1e-10,
    solver_max_iter = 100L, strict = TRUE, keep_permutations = FALSE) {
  call <- match.call()
  x <- .ch7cs_data(x)
  clusters <- .ch7cs_count(clusters, "clusters", 1L, nrow(x))
  selection <- match.arg(selection)
  initialization <- match.arg(initialization)
  empty_action <- .ch7cs_empty_action(empty_action)
  ties <- .ch7cs_ties(ties)
  keep_permutations <- .ch7cs_flag(keep_permutations, "keep_permutations")
  s_grid <- as.numeric(s_grid)
  if (length(s_grid) < 1L || anyNA(s_grid) || any(!is.finite(s_grid)) ||
      anyDuplicated(s_grid)) {
    stop("`s_grid` must contain finite distinct candidate values.",
         call. = FALSE)
  }
  if (any(s_grid < 1 | s_grid > sqrt(ncol(x)))) {
    stop("Every `s_grid` value must lie in [1, sqrt(p)].", call. = FALSE)
  }
  s_grid <- sort(s_grid)
  indices <- .ch7cs_permutation_array(
    x, permutation_indices, n_permutations, permutation_seed, selection
  )
  B <- dim(indices)[3L]

  fit_one <- function(data, bound) sparse_kmeans(
    data, clusters = clusters, s = bound, initial = initial,
    initialization = initialization, first_index = first_index,
    seed = clustering_seed, ties = ties, empty_action = empty_action,
    weight_tol = weight_tol, weight_max_iter = weight_max_iter,
    solver_tol = solver_tol, solver_max_iter = solver_max_iter,
    strict = strict, keep_distances = FALSE
  )
  observed.fits <- lapply(s_grid, function(bound) fit_one(x, bound))
  observed.objective <- vapply(observed.fits, `[[`, numeric(1), "objective")
  permutation.objective <- matrix(
    NA_real_, nrow = B, ncol = length(s_grid),
    dimnames = list(paste0("permutation", seq_len(B)),
                    paste0("s=", format(s_grid)))
  )
  for (b in seq_len(B)) {
    reference <- x
    for (j in seq_len(ncol(x))) {
      reference[, j] <- x[indices[, j, b], j]
    }
    for (g in seq_along(s_grid)) {
      permutation.objective[b, g] <- fit_one(reference, s_grid[g])$objective
    }
  }
  all.objectives <- c(observed.objective, permutation.objective)
  if (any(!is.finite(all.objectives)) || any(all.objectives <= 0)) {
    stop("Every observed and permutation objective O(s) must be finite and positive.",
         call. = FALSE)
  }
  log.reference <- log(permutation.objective)
  gap <- log(observed.objective) - colMeans(log.reference)
  reference.sd <- if (B >= 2L) apply(log.reference, 2L, stats::sd) else
    rep.int(NA_real_, length(s_grid))
  standard.error <- sqrt(1 + 1 / B) * reference.sd
  maximum.index <- unname(which(gap == max(gap))[1L])
  selected.index <- if (selection == "max") maximum.index else {
    eligible <- which(gap >= gap[maximum.index] -
                        standard.error[maximum.index])
    unname(eligible[1L])
  }
  table <- data.frame(
    s = s_grid, observed_objective = observed.objective,
    mean_log_reference_objective = colMeans(log.reference),
    sd_log_reference_objective = reference.sd,
    standard_error = standard.error, gap = gap,
    row.names = NULL, check.names = FALSE
  )
  answer <- list(
    call = call, method = paste0("sparse K-means permutation Gap (", selection,
                                 ")"),
    selected_s = s_grid[selected.index], selected_index = selected.index,
    selected_fit = observed.fits[[selected.index]], fits = observed.fits,
    calibration = table, selection = selection,
    permutation_objectives = if (keep_permutations)
      permutation.objective else NULL,
    permutation_indices = if (keep_permutations) indices else NULL,
    diagnostics = list(
      permutations = B,
      permutation_source = if (is.null(permutation_indices))
        "featurewise seeded permutations" else "supplied index array",
      permutation_seed = if (is.null(permutation_indices))
        as.integer(permutation_seed) else NULL,
      one_se_definition = paste0(
        "sqrt(1 + 1/B) * sd(log O_b); compare to SE at smallest Gap maximizer"
      ),
      objective_scale = "book ordered-pair BCSS; common factor 2 cancels in Gap"
    )
  )
  class(answer) <- "sparse_kmeans_selection"
  answer
}


#' Coordinatewise sparse K-median clustering
#'
#' Implement the Chapter 7 sparse K-median baseline. Cluster centers and the
#' global center are coordinatewise medians, using the midpoint for even
#' samples. Feature `j` receives improvement
#' `D_j = sum_i |x_ij-median_j(x)| - sum_k sum_(i in Ck)
#' |x_ij-median_j(Ck)|`. The non-negative weight block obeys the same L2/L1
#' constraints and soft-threshold rule as sparse K-means; assignment minimizes
#' weighted L1 distance.
#'
#' This is the book-defined axis-aligned baseline, not a claim of a unique
#' primary sparse-median algorithm and not a rotation-equivariant spatial-
#' median method. All-zero improvements error. A negative score within the
#' explicit `score_tol` roundoff bound is set to zero and counted; a larger
#' violation errors.
#'
#' @inheritParams sparse_kmeans
#' @param s_w Required L1 weight bound in `[1, sqrt(p)]`.
#' @param score_tol Explicit relative roundoff tolerance for the theoretically
#'   non-negative median-improvement scores.
#'
#' @return A `sparse_kmedian_fit` object with midpoint medians, feature
#'   improvements, weights, labels, weighted-L1 objective and certificates.
#'
#' @references
#' Feng, L. (2026). *High-Dimensional Data Analysis for Elliptical Symmetric
#' Distributions*, Chapter 7 (book manuscript). This is the book-defined
#' axis-aligned baseline, not a separately attributed primary method.
#' @export
#' @examples
#' x <- rbind(c(-4, 0), c(-2, 0), c(2, 0), c(8, 0))
#' sparse_kmedian(x, clusters = 2, s_w = 1, initial = c(1, 3))
sparse_kmedian <- function(
    x, clusters, s_w, initial = NULL,
    initialization = c("maxmin", "random"), first_index = 1L, seed = NULL,
    ties = "first", empty_action = c("error", "farthest"),
    score_tol = 1e-12, weight_tol = 1e-10, weight_max_iter = 200L,
    solver_tol = 1e-10, solver_max_iter = 100L, strict = TRUE,
    keep_distances = FALSE) {
  call <- match.call()
  fit <- .ch7cs_sparse_fit(
    x, clusters, s_w, initial, initialization, first_index, seed, ties,
    empty_action, weight_tol, weight_max_iter, solver_tol, solver_max_iter,
    strict, keep_distances, median = TRUE, score_tol = score_tol
  )
  core <- fit$core
  answer <- list(
    call = call, method = fit$method, cluster = core$cluster,
    centers = core$centers, global_median = core$global_median,
    size = core$size, weights = core$weights,
    feature_scores = core$feature_scores, objective = core$objective,
    objective_trace = core$objective_trace, iterations = core$iterations,
    converged = isTRUE(core$converged), valid = isTRUE(core$certified),
    s_w = fit$bound, initialization = fit$start, distances = core$distances,
    diagnostics = c(core[c(
      "cycle_detected", "assignment_ties", "empty_repairs", "repair_ties",
      "score_roundoff_adjustments", "maximum_objective_decrease",
      "assignment_residual", "weight_certificate", "weighted_within_loss",
      "certified"
    )], list(
      median_convention = "midpoint of the two middle order statistics",
      assignment_geometry = "coordinatewise weighted L1",
      rotation_equivariant = FALSE,
      method_scope = "book-defined sparse K-median baseline",
      no_implicit_repair = fit$empty_action == "error"
    ))
  )
  class(answer) <- c("sparse_kmedian_fit", "hd_clustering_fit")
  answer
}
