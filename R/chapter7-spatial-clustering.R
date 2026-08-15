.c7sc_flag <- function(x, name) {
  if (!is.logical(x) || length(x) != 1L || is.na(x)) {
    stop(sprintf("`%s` must be TRUE or FALSE.", name), call. = FALSE)
  }
  x
}


.c7sc_integer <- function(x, name, minimum = 1L, maximum = Inf) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) ||
      !is.finite(x) || x != floor(x) || x < minimum || x > maximum ||
      x > .Machine$integer.max) {
    stop(sprintf("`%s` must be an integer between %s and %s.",
                 name, format(minimum), format(maximum)), call. = FALSE)
  }
  as.integer(x)
}


.c7sc_choice <- function(x, choices, name) {
  if (!is.character(x) || length(x) < 1L || anyNA(x)) {
    stop(sprintf("`%s` must be one of: %s.",
                 name, paste(choices, collapse = ", ")), call. = FALSE)
  }
  if (length(x) > 1L) {
    if (!identical(x, choices)) {
      stop(sprintf("`%s` must be one of: %s.",
                   name, paste(choices, collapse = ", ")), call. = FALSE)
    }
    x <- x[1L]
  }
  if (!x %in% choices) {
    stop(sprintf("`%s` must be one of: %s.",
                 name, paste(choices, collapse = ", ")), call. = FALSE)
  }
  x
}


.c7sc_validate_common <- function(x, K, tol, max_iter,
                                  spatial_max_iter, zero_tol,
                                  empty_action, cycle_action, ties,
                                  keep_distances) {
  x <- .as_data_matrix(x)
  K <- .c7sc_integer(K, "K", 1L, nrow(x))
  max_iter <- .c7sc_integer(max_iter, "max_iter")
  controls <- .validate_iteration_controls(tol, spatial_max_iter, zero_tol)
  empty_action <- .c7sc_choice(
    empty_action, c("error", "farthest"), "empty_action"
  )
  cycle_action <- .c7sc_choice(
    cycle_action, c("error", "return"), "cycle_action"
  )
  ties <- .c7sc_choice(ties, "first", "ties")
  keep_distances <- .c7sc_flag(keep_distances, "keep_distances")
  list(
    x = x, K = K, tol = controls$tol, max_iter = max_iter,
    spatial_max_iter = controls$max_iter, zero_tol = controls$zero_tol,
    empty_action = empty_action, cycle_action = cycle_action,
    ties = ties, keep_distances = keep_distances
  )
}


.c7sc_center_matrix <- function(centers, K, p, variable.names,
                                name = "init") {
  if (is.data.frame(centers)) {
    centers <- data.matrix(centers)
  }
  if (!is.matrix(centers) || !is.numeric(centers) ||
      !identical(dim(centers), c(K, p)) || anyNA(centers) ||
      any(!is.finite(centers))) {
    stop(sprintf("`%s` centers must be a finite %d by %d numeric matrix.",
                 name, K, p), call. = FALSE)
  }
  storage.mode(centers) <- "double"
  dimnames(centers) <- list(paste0("C", seq_len(K)), variable.names)
  centers
}


.c7sc_assignment_euclidean <- function(x, centers, active,
                                        empty_action, K) {
  raw <- cpp_ch7sc_assign_euclidean(
    x, centers, as.integer(active - 1L), TRUE
  )
  repaired <- .c7sc_repair_empty(
    as.integer(raw$labels), raw$distance_matrix, empty_action, K
  )
  rows <- seq_len(nrow(x))
  list(
    labels = repaired$labels,
    assigned_distance = as.numeric(
      raw$distance_matrix[cbind(rows, repaired$labels)]
    ),
    distance_matrix = raw$distance_matrix,
    tie_count = as.integer(raw$tie_count),
    repairs = repaired$repairs
  )
}


.c7sc_assignment_metric <- function(x, centers, inverse_metric,
                                     empty_action, K) {
  raw <- cpp_ch7sc_assign_metric(x, centers, inverse_metric, TRUE)
  repaired <- .c7sc_repair_empty(
    as.integer(raw$labels), raw$distance_matrix, empty_action, K
  )
  rows <- seq_len(nrow(x))
  list(
    labels = repaired$labels,
    assigned_distance = as.numeric(
      raw$distance_matrix[cbind(rows, repaired$labels)]
    ),
    distance_matrix = raw$distance_matrix,
    tie_count = as.integer(raw$tie_count),
    repairs = repaired$repairs
  )
}


.c7sc_repair_empty <- function(labels, distance_matrix, empty_action, K) {
  sizes <- tabulate(labels, nbins = K)
  empty <- which(sizes == 0L)
  if (!length(empty)) {
    return(list(labels = labels, repairs = list()))
  }
  if (identical(empty_action, "error")) {
    stop(sprintf(
      "The assignment produced empty cluster(s) %s; no repair was applied.",
      paste(empty, collapse = ", ")
    ), call. = FALSE)
  }
  repairs <- vector("list", length(empty))
  for (index in seq_along(empty)) {
    target <- empty[index]
    candidates <- which(sizes[labels] >= 2L)
    if (!length(candidates)) {
      stop("No donor cluster of size at least two is available for the " %+%
             "explicit farthest-point repair.", call. = FALSE)
    }
    donor_distance <- distance_matrix[cbind(candidates, labels[candidates])]
    selected <- candidates[which(donor_distance == max(donor_distance))[1L]]
    source <- labels[selected]
    repairs[[index]] <- list(
      row = as.integer(selected), from = as.integer(source),
      to = as.integer(target), distance = as.numeric(
        distance_matrix[selected, source]
      )
    )
    labels[selected] <- target
    sizes[source] <- sizes[source] - 1L
    sizes[target] <- sizes[target] + 1L
  }
  if (any(tabulate(labels, nbins = K) == 0L)) {
    stop("The explicit farthest-point repair did not remove every empty " %+%
           "cluster.", call. = FALSE)
  }
  list(labels = as.integer(labels), repairs = repairs)
}


.c7sc_initial_centers <- function(x, K, init, first_index) {
  n <- nrow(x)
  p <- ncol(x)
  variable.names <- colnames(x)
  if (is.null(variable.names)) {
    variable.names <- paste0("V", seq_len(p))
  }
  first_index <- .c7sc_integer(first_index, "first_index", 1L, n)
  if (is.character(init)) {
    if (length(init) != 1L || is.na(init) || !identical(init, "maxmin")) {
      stop("Character `init` must be exactly \"maxmin\".", call. = FALSE)
    }
    indices <- first_index
    if (K > 1L) {
      for (cluster in 2:K) {
        current <- x[indices, , drop = FALSE]
        assignment <- cpp_ch7sc_assign_euclidean(
          x, current, as.integer(seq_len(p) - 1L), FALSE
        )
        minimum <- as.numeric(assignment$assigned_distance)
        maximum <- max(minimum)
        if (!is.finite(maximum) || maximum <= 0) {
          stop("Deterministic max-min initialization needs at least `K` " %+%
                 "distinct rows.", call. = FALSE)
        }
        candidates <- which(minimum == maximum)
        indices <- c(indices, candidates[1L])
      }
    }
    centers <- x[indices, , drop = FALSE]
    dimnames(centers) <- list(paste0("C", seq_len(K)), variable.names)
    return(list(
      centers = centers, indices = as.integer(indices),
      source = "deterministic max-min; first maximum index resolves ties"
    ))
  }
  if (is.matrix(init) || is.data.frame(init)) {
    return(list(
      centers = .c7sc_center_matrix(init, K, p, variable.names),
      indices = NULL, source = "supplied center matrix"
    ))
  }
  indices <- as.numeric(init)
  if (length(indices) != K || anyNA(indices) || any(!is.finite(indices)) ||
      any(indices != floor(indices)) || any(indices < 1L) ||
      any(indices > n) || anyDuplicated(indices)) {
    stop("Numeric `init` must be a length-`K` vector of distinct valid row " %+%
           "indices, or a `K` by `p` center matrix.", call. = FALSE)
  }
  indices <- as.integer(indices)
  centers <- x[indices, , drop = FALSE]
  dimnames(centers) <- list(paste0("C", seq_len(K)), variable.names)
  list(centers = centers, indices = indices, source = "supplied row indices")
}


.c7sc_update_centers <- function(x, labels, K, tol,
                                  spatial_max_iter, zero_tol) {
  variable.names <- colnames(x)
  if (is.null(variable.names)) {
    variable.names <- paste0("V", seq_len(ncol(x)))
  }
  centers <- matrix(NA_real_, K, ncol(x),
                    dimnames = list(paste0("C", seq_len(K)), variable.names))
  certificates <- vector("list", K)
  for (cluster in seq_len(K)) {
    rows <- which(labels == cluster)
    if (!length(rows)) {
      stop(sprintf("Cluster %d is empty before its spatial-median update; " %+%
                     "no repair was applied.", cluster), call. = FALSE)
    }
    estimate <- spatial_median(
      x[rows, , drop = FALSE], tol = tol,
      max_iter = spatial_max_iter, zero_tol = zero_tol, warn = FALSE
    )
    certificate <- list(
      cluster = as.integer(cluster), size = as.integer(length(rows)),
      converged = isTRUE(attr(estimate, "converged")),
      iterations = as.integer(attr(estimate, "iterations")),
      objective = as.numeric(attr(estimate, "objective")),
      equation_residual = as.numeric(attr(estimate, "equation_residual"))
    )
    if (!certificate$converged) {
      stop(sprintf(
        "The spatial median for cluster %d failed its convergence certificate.",
        cluster
      ), call. = FALSE)
    }
    centers[cluster, ] <- as.numeric(estimate)
    certificates[[cluster]] <- certificate
  }
  list(centers = centers, certificates = certificates)
}


.c7sc_state_key <- function(labels, active = integer()) {
  paste0(paste(labels, collapse = ","), "|", paste(active, collapse = ","))
}


.c7sc_active <- function(scores, tau, empty_active) {
  active <- which(scores >= tau)
  fallback <- "none"
  if (!length(active)) {
    if (identical(empty_active, "error")) {
      stop("The threshold produced an empty active set; no fallback was " %+%
             "applied.", call. = FALSE)
    }
    if (identical(empty_active, "largest")) {
      active <- which(scores == max(scores))[1L]
      fallback <- "first feature attaining the largest score"
    } else {
      active <- seq_along(scores)
      fallback <- "all features"
    }
  }
  list(active = as.integer(active), fallback = fallback)
}


.c7sc_geometry <- function(x, centers, labels, active, overall = NULL) {
  if (is.null(overall)) {
    overall <- rep.int(0, ncol(x))
  }
  cpp_ch7sc_geometry(
    x, centers, as.integer(labels), as.numeric(overall),
    as.integer(active - 1L)
  )
}


.c7sc_certified_median <- function(x, tol, spatial_max_iter, zero_tol) {
  estimate <- spatial_median(
    x, tol = tol, max_iter = spatial_max_iter,
    zero_tol = zero_tol, warn = FALSE
  )
  if (!isTRUE(attr(estimate, "converged"))) {
    stop("A selector spatial median failed its convergence certificate.",
         call. = FALSE)
  }
  as.numeric(estimate)
}


.c7sc_common_diagnostics <- function(controls, initialization, iterations,
                                      converged, termination, tie_count,
                                      repairs, cycle_length = 0L) {
  list(
    converged = converged, iterations = as.integer(iterations),
    termination = termination, cycle.length = as.integer(cycle_length),
    initialization = initialization$source,
    initialization.indices = initialization$indices,
    first.index = as.integer(controls$first_index),
    tie.rule = "first minimum cluster, row, or feature index",
    tie.count = as.integer(tie_count),
    empty.action = controls$empty_action,
    repairs = repairs, cycle.action = controls$cycle_action,
    spatial.median.tolerance = controls$tol,
    spatial.median.maximum.iterations = controls$spatial_max_iter,
    zero.tolerance = controls$zero_tol,
    no.hidden.repair = paste(
      "empty clusters, empty active sets, cycles, and failed spatial-median",
      "certificates are not silently repaired"
    )
  )
}


#' K-spatial-median clustering
#'
#' Alternates exact sample spatial-median center updates with deterministic
#' nearest-center assignments. The fitted objective is the sum of unsquared
#' Euclidean distances. Empty-cluster repair is disabled by default.
#'
#' @param x Numeric observation-by-variable matrix or data frame.
#' @param K Number of clusters.
#' @param init Either `"maxmin"`, `K` distinct row indices, or a finite
#'   `K` by `p` center matrix.
#' @param first_index First max-min seed when `init = "maxmin"`.
#' @param tol Positive spatial-median equation tolerance.
#' @param max_iter Positive outer iteration limit.
#' @param spatial_max_iter Positive modified-Weiszfeld iteration limit.
#' @param zero_tol Non-negative zero-residual tolerance.
#' @param empty_action `"error"` or the explicit deterministic `"farthest"`
#'   repair. The latter only donates from a cluster of size at least two.
#' @param cycle_action Whether a detected repeated state is an `"error"` or is
#'   returned with a failed convergence certificate.
#' @param ties Deterministic tie rule; currently only `"first"` is supported.
#' @param keep_distances Whether to retain the final full distance matrix.
#'
#' @return A `k_spatial_median_fit` object with labels, centers, objective, and
#'   explicit convergence, tie, and repair diagnostics.
#' @export
#'
#' @references
#' Zhao, P., Zhuang, D., and Feng, L. (2026). Sparse K-spatial-median
#' clustering for high-dimensional data. arXiv:2605.00598.
#'
#' @examples
#' x <- matrix(c(0, 2, 8, 10), ncol = 1)
#' k_spatial_median(x, K = 2, init = c(1, 4))
k_spatial_median <- function(
    x, K, init = "maxmin", first_index = 1L, tol = 1e-8,
    max_iter = 100L, spatial_max_iter = 500L, zero_tol = 0,
    empty_action = c("error", "farthest"),
    cycle_action = c("error", "return"), ties = "first",
    keep_distances = FALSE) {
  call <- match.call()
  controls <- .c7sc_validate_common(
    x, K, tol, max_iter, spatial_max_iter, zero_tol,
    empty_action, cycle_action, ties, keep_distances
  )
  controls$first_index <- .c7sc_integer(
    first_index, "first_index", 1L, nrow(controls$x)
  )
  x <- controls$x
  initialization <- .c7sc_initial_centers(
    x, controls$K, init, controls$first_index
  )
  assignment <- .c7sc_assignment_euclidean(
    x, initialization$centers, seq_len(ncol(x)),
    controls$empty_action, controls$K
  )
  labels <- assignment$labels
  objective.history <- sum(assignment$assigned_distance)
  tie.count <- assignment$tie_count
  repairs <- assignment$repairs
  seen <- .c7sc_state_key(labels)
  converged <- FALSE
  termination <- "max_iter"
  cycle.length <- 0L
  center.certificates <- NULL
  centers <- initialization$centers
  iterations <- 0L
  for (iteration in seq_len(controls$max_iter)) {
    iterations <- iteration
    updated <- .c7sc_update_centers(
      x, labels, controls$K, controls$tol,
      controls$spatial_max_iter, controls$zero_tol
    )
    centers <- updated$centers
    center.certificates <- updated$certificates
    next.assignment <- .c7sc_assignment_euclidean(
      x, centers, seq_len(ncol(x)), controls$empty_action, controls$K
    )
    objective.history <- c(
      objective.history, sum(next.assignment$assigned_distance)
    )
    tie.count <- tie.count + next.assignment$tie_count
    repairs <- c(repairs, next.assignment$repairs)
    next.labels <- next.assignment$labels
    assignment <- next.assignment
    if (identical(next.labels, labels)) {
      converged <- TRUE
      termination <- "fixed_point"
      break
    }
    key <- .c7sc_state_key(next.labels)
    if (key %in% seen) {
      cycle.length <- length(seen) - match(key, seen) + 1L
      if (identical(controls$cycle_action, "error")) {
        stop(sprintf("A repeated clustering state (cycle length %d) was " %+%
                       "detected; no repair was applied.", cycle.length),
             call. = FALSE)
      }
      labels <- next.labels
      termination <- "cycle"
      break
    }
    seen <- c(seen, key)
    labels <- next.labels
  }
  if (!converged) {
    updated <- .c7sc_update_centers(
      x, labels, controls$K, controls$tol,
      controls$spatial_max_iter, controls$zero_tol
    )
    centers <- updated$centers
    center.certificates <- updated$certificates
    assignment <- .c7sc_assignment_euclidean(
      x, centers, seq_len(ncol(x)), controls$empty_action, controls$K
    )
  }
  geometry <- .c7sc_geometry(
    x, centers, labels, seq_len(ncol(x))
  )
  scale <- max(1, abs(objective.history))
  monotone <- !length(repairs) &&
    all(diff(objective.history) <= sqrt(.Machine$double.eps) * scale)
  diagnostics <- .c7sc_common_diagnostics(
    controls, initialization, iterations, converged, termination,
    tie.count, repairs, cycle.length
  )
  diagnostics$objective.history <- objective.history
  diagnostics$objective.monotone <- monotone
  diagnostics$center.certificates <- center.certificates
  diagnostics$next.labels <- assignment$labels
  diagnostics$next.assignment.repairs <- assignment$repairs
  structure(
    list(
      method = "K-spatial-median clustering", labels = as.integer(labels),
      centers = centers, sizes = as.integer(tabulate(labels, controls$K)),
      objective = as.numeric(geometry$within_sum),
      distances = if (controls$keep_distances) {
        assignment$distance_matrix
      } else NULL,
      diagnostics = diagnostics, call = call
    ),
    class = c("k_spatial_median_fit", "chapter7_clustering_fit", "list")
  )
}


#' Spatial-median clustering with a common SSCM metric
#'
#' Uses full-dimensional spatial medians and the explicitly regularized common
#' spatial-sign covariance matrix `crossprod(U) / n + lambda * I`. The inverse
#' is an exact SPD inverse; no pseudoinverse or eigenvalue repair is used.
#'
#' @inheritParams k_spatial_median
#' @param lambda Strictly positive ridge appearing in the stated SSCM method.
#' @param keep_signs Whether to retain final residual spatial signs.
#'
#' @return An `sm_sscm_fit` object. Its diagnostics explicitly do not claim
#'   monotonicity of a global objective.
#'
#' @references
#' Zhao, P., Zhuang, D., and Feng, L. (2026). Sparse K-spatial-median
#' clustering for high-dimensional data. arXiv:2605.00598.
#' \url{https://arxiv.org/abs/2605.00598}.
#' @export
#'
#' @examples
#' x <- rbind(c(-2, 0), c(-1, 0), c(1, 0), c(2, 0))
#' sm_sscm(x, K = 2, lambda = 0.1, init = c(1, 4))
sm_sscm <- function(
    x, K, lambda, init = "maxmin", first_index = 1L, tol = 1e-8,
    max_iter = 100L, spatial_max_iter = 500L, zero_tol = 0,
    empty_action = c("error", "farthest"),
    cycle_action = c("error", "return"), ties = "first",
    keep_distances = FALSE, keep_signs = FALSE) {
  call <- match.call()
  controls <- .c7sc_validate_common(
    x, K, tol, max_iter, spatial_max_iter, zero_tol,
    empty_action, cycle_action, ties, keep_distances
  )
  controls$first_index <- .c7sc_integer(
    first_index, "first_index", 1L, nrow(controls$x)
  )
  keep_signs <- .c7sc_flag(keep_signs, "keep_signs")
  lambda <- as.numeric(lambda)
  if (length(lambda) != 1L || is.na(lambda) || !is.finite(lambda) ||
      lambda <= 0) {
    stop("`lambda` must be a finite strictly positive number.", call. = FALSE)
  }
  x <- controls$x
  initialization <- .c7sc_initial_centers(
    x, controls$K, init, controls$first_index
  )
  initial.assignment <- .c7sc_assignment_euclidean(
    x, initialization$centers, seq_len(ncol(x)),
    controls$empty_action, controls$K
  )
  labels <- initial.assignment$labels
  tie.count <- initial.assignment$tie_count
  repairs <- initial.assignment$repairs
  seen <- .c7sc_state_key(labels)
  converged <- FALSE
  termination <- "max_iter"
  cycle.length <- 0L
  criterion.history <- numeric()
  centers <- initialization$centers
  metric <- NULL
  assignment <- NULL
  center.certificates <- NULL
  iterations <- 0L
  for (iteration in seq_len(controls$max_iter)) {
    iterations <- iteration
    updated <- .c7sc_update_centers(
      x, labels, controls$K, controls$tol,
      controls$spatial_max_iter, controls$zero_tol
    )
    centers <- updated$centers
    center.certificates <- updated$certificates
    metric <- cpp_ch7sc_sscm_metric(
      x, centers, as.integer(labels), lambda, controls$zero_tol, keep_signs
    )
    next.assignment <- .c7sc_assignment_metric(
      x, centers, metric$inverse, controls$empty_action, controls$K
    )
    criterion.history <- c(
      criterion.history, sum(next.assignment$assigned_distance)
    )
    tie.count <- tie.count + next.assignment$tie_count
    repairs <- c(repairs, next.assignment$repairs)
    next.labels <- next.assignment$labels
    assignment <- next.assignment
    if (identical(next.labels, labels)) {
      converged <- TRUE
      termination <- "fixed_point"
      break
    }
    key <- .c7sc_state_key(next.labels)
    if (key %in% seen) {
      cycle.length <- length(seen) - match(key, seen) + 1L
      if (identical(controls$cycle_action, "error")) {
        stop(sprintf("A repeated SM--SSCM state (cycle length %d) was " %+%
                       "detected; no repair was applied.", cycle.length),
             call. = FALSE)
      }
      labels <- next.labels
      termination <- "cycle"
      break
    }
    seen <- c(seen, key)
    labels <- next.labels
  }
  if (!converged) {
    updated <- .c7sc_update_centers(
      x, labels, controls$K, controls$tol,
      controls$spatial_max_iter, controls$zero_tol
    )
    centers <- updated$centers
    center.certificates <- updated$certificates
    metric <- cpp_ch7sc_sscm_metric(
      x, centers, as.integer(labels), lambda, controls$zero_tol, keep_signs
    )
    assignment <- .c7sc_assignment_metric(
      x, centers, metric$inverse, controls$empty_action, controls$K
    )
  }
  diagnostics <- .c7sc_common_diagnostics(
    controls, initialization, iterations, converged, termination,
    tie.count, repairs, cycle.length
  )
  diagnostics$center.certificates <- center.certificates
  diagnostics$assignment.criterion.history <- criterion.history
  diagnostics$global.objective.monotonicity <- paste(
    "not claimed: center and changing SSCM metric updates are not block",
    "minimizers of one stated objective"
  )
  diagnostics$n.zero.residuals <- as.integer(metric$n_zero)
  diagnostics$metric.trace.expected <- as.numeric(metric$trace_expected)
  diagnostics$metric.trace.error <- as.numeric(metric$trace_error)
  diagnostics$metric.minimum.eigenvalue <- as.numeric(
    metric$minimum_eigenvalue
  )
  diagnostics$inverse.residual <- as.numeric(metric$inverse_residual)
  diagnostics$next.labels <- assignment$labels
  diagnostics$next.assignment.repairs <- assignment$repairs
  structure(
    list(
      method = "SM--SSCM clustering", labels = as.integer(labels),
      centers = centers, sizes = as.integer(tabulate(labels, controls$K)),
      metric = metric$metric, inverse_metric = metric$inverse,
      signs = metric$signs,
      distances = if (controls$keep_distances) {
        assignment$distance_matrix
      } else NULL,
      lambda = lambda, diagnostics = diagnostics, call = call
    ),
    class = c("sm_sscm_fit", "chapter7_clustering_fit", "list")
  )
}


#' Sparse K-spatial-median clustering
#'
#' Updates full-dimensional spatial medians, hard-screens coordinates by
#' across-center separation, and assigns observations in the retained
#' subspace. This procedure is a deterministic block-coordinate heuristic and
#' is not represented as optimizing a single global objective.
#'
#' @inheritParams k_spatial_median
#' @param tau Finite non-negative hard-screening threshold. Equality is
#'   retained: a coordinate is active when its score is at least `tau`.
#' @param empty_active Action when no score reaches `tau`: `"error"`, the first
#'   `"largest"` score, or `"all"` coordinates.
#' @param reset_excluded Output-only treatment of excluded center coordinates.
#'   `"overall_spatial_median"` never feeds reset values back into iteration.
#'
#' @return A `sparse_k_spatial_median_fit` object with both reported centers
#'   and untouched `algorithm_centers`, the active set, scores, and diagnostics.
#'
#' @references
#' Zhao, P., Zhuang, D., and Feng, L. (2026). Sparse K-spatial-median
#' clustering for high-dimensional data. arXiv:2605.00598.
#' \url{https://arxiv.org/abs/2605.00598}.
#' @export
#'
#' @examples
#' x <- rbind(c(-2, 0), c(-1, 0), c(1, 0), c(2, 0))
#' sparse_k_spatial_median(x, K = 2, tau = 0.5, init = c(1, 4))
sparse_k_spatial_median <- function(
    x, K, tau, init = "maxmin", first_index = 1L, tol = 1e-8,
    max_iter = 100L, spatial_max_iter = 500L, zero_tol = 0,
    empty_action = c("error", "farthest"),
    empty_active = c("error", "largest", "all"),
    cycle_action = c("error", "return"), ties = "first",
    reset_excluded = c("none", "overall_spatial_median"),
    keep_distances = FALSE) {
  call <- match.call()
  controls <- .c7sc_validate_common(
    x, K, tol, max_iter, spatial_max_iter, zero_tol,
    empty_action, cycle_action, ties, keep_distances
  )
  controls$first_index <- .c7sc_integer(
    first_index, "first_index", 1L, nrow(controls$x)
  )
  empty_active <- .c7sc_choice(
    empty_active, c("error", "largest", "all"), "empty_active"
  )
  reset_excluded <- .c7sc_choice(
    reset_excluded, c("none", "overall_spatial_median"),
    "reset_excluded"
  )
  tau <- as.numeric(tau)
  if (length(tau) != 1L || is.na(tau) || !is.finite(tau) || tau < 0) {
    stop("`tau` must be a finite non-negative number.", call. = FALSE)
  }
  x <- controls$x
  initialization <- .c7sc_initial_centers(
    x, controls$K, init, controls$first_index
  )
  initial.assignment <- .c7sc_assignment_euclidean(
    x, initialization$centers, seq_len(ncol(x)),
    controls$empty_action, controls$K
  )
  labels <- initial.assignment$labels
  tie.count <- initial.assignment$tie_count
  repairs <- initial.assignment$repairs
  seen <- .c7sc_state_key(labels)
  converged <- FALSE
  termination <- "max_iter"
  cycle.length <- 0L
  criterion.history <- numeric()
  centers <- initialization$centers
  assignment <- NULL
  active <- seq_len(ncol(x))
  scores <- rep.int(0, ncol(x))
  active.history <- list()
  score.history <- list()
  fallback.history <- character()
  center.certificates <- NULL
  iterations <- 0L
  for (iteration in seq_len(controls$max_iter)) {
    iterations <- iteration
    updated <- .c7sc_update_centers(
      x, labels, controls$K, controls$tol,
      controls$spatial_max_iter, controls$zero_tol
    )
    centers <- updated$centers
    center.certificates <- updated$certificates
    score.result <- cpp_ch7sc_feature_scores(centers)
    scores <- as.numeric(score.result$scores)
    active.result <- .c7sc_active(scores, tau, empty_active)
    active <- active.result$active
    active.history[[iteration]] <- active
    score.history[[iteration]] <- scores
    fallback.history[iteration] <- active.result$fallback
    next.assignment <- .c7sc_assignment_euclidean(
      x, centers, active, controls$empty_action, controls$K
    )
    criterion.history <- c(
      criterion.history, sum(next.assignment$assigned_distance)
    )
    tie.count <- tie.count + next.assignment$tie_count
    repairs <- c(repairs, next.assignment$repairs)
    next.labels <- next.assignment$labels
    assignment <- next.assignment
    if (identical(next.labels, labels)) {
      converged <- TRUE
      termination <- "fixed_point"
      break
    }
    key <- .c7sc_state_key(next.labels, active)
    if (key %in% seen) {
      cycle.length <- length(seen) - match(key, seen) + 1L
      if (identical(controls$cycle_action, "error")) {
        stop(sprintf("A repeated Sparse--SM state (cycle length %d) was " %+%
                       "detected; no repair was applied.", cycle.length),
             call. = FALSE)
      }
      labels <- next.labels
      termination <- "cycle"
      break
    }
    seen <- c(seen, key)
    labels <- next.labels
  }
  if (!converged) {
    updated <- .c7sc_update_centers(
      x, labels, controls$K, controls$tol,
      controls$spatial_max_iter, controls$zero_tol
    )
    centers <- updated$centers
    center.certificates <- updated$certificates
    score.result <- cpp_ch7sc_feature_scores(centers)
    scores <- as.numeric(score.result$scores)
    active.result <- .c7sc_active(scores, tau, empty_active)
    active <- active.result$active
    assignment <- .c7sc_assignment_euclidean(
      x, centers, active, controls$empty_action, controls$K
    )
  }
  algorithm.centers <- centers
  reported.centers <- centers
  reset.baseline <- NULL
  if (identical(reset_excluded, "overall_spatial_median") &&
      length(active) < ncol(x)) {
    reset.baseline <- .c7sc_certified_median(
      x, controls$tol, controls$spatial_max_iter, controls$zero_tol
    )
    inactive <- setdiff(seq_len(ncol(x)), active)
    reported.centers[, inactive] <- matrix(
      reset.baseline[inactive], nrow = controls$K,
      ncol = length(inactive), byrow = TRUE
    )
  }
  diagnostics <- .c7sc_common_diagnostics(
    controls, initialization, iterations, converged, termination,
    tie.count, repairs, cycle.length
  )
  diagnostics$center.certificates <- center.certificates
  diagnostics$assignment.criterion.history <- criterion.history
  diagnostics$active.history <- active.history
  diagnostics$score.history <- score.history
  diagnostics$empty.active.action <- empty_active
  diagnostics$empty.active.fallback.history <- fallback.history
  diagnostics$active.boundary.rule <- "score >= tau is retained"
  diagnostics$global.objective.monotonicity <- paste(
    "not claimed: Sparse--SM alternates clustering and feature screening",
    "without one stated global objective"
  )
  diagnostics$reset.excluded <- reset_excluded
  diagnostics$reset.output.only <- TRUE
  diagnostics$reset.baseline <- reset.baseline
  diagnostics$next.labels <- assignment$labels
  diagnostics$next.assignment.repairs <- assignment$repairs
  structure(
    list(
      method = "Sparse K-spatial-median clustering",
      labels = as.integer(labels), centers = reported.centers,
      algorithm_centers = algorithm.centers,
      sizes = as.integer(tabulate(labels, controls$K)),
      active_set = as.integer(active), scores = scores, tau = tau,
      distances = if (controls$keep_distances) {
        assignment$distance_matrix
      } else NULL,
      diagnostics = diagnostics, call = call
    ),
    class = c(
      "sparse_k_spatial_median_fit", "chapter7_clustering_fit", "list"
    )
  )
}


.c7sc_validate_tau_grid <- function(tau_grid) {
  tau_grid <- as.numeric(tau_grid)
  if (!length(tau_grid) || anyNA(tau_grid) || any(!is.finite(tau_grid)) ||
      any(tau_grid < 0)) {
    stop("`tau_grid` must contain finite non-negative values.",
         call. = FALSE)
  }
  sort(unique(tau_grid))
}


.c7sc_validate_permutation <- function(permutation, n, p, index) {
  if (!is.matrix(permutation) || !is.numeric(permutation) ||
      !identical(dim(permutation), c(n, p)) || anyNA(permutation) ||
      any(!is.finite(permutation)) || any(permutation != floor(permutation))) {
    stop(sprintf("Permutation %d must be an integer-valued %d by %d matrix.",
                 index, n, p), call. = FALSE)
  }
  permutation <- matrix(as.integer(permutation), n, p)
  target <- seq_len(n)
  valid <- vapply(seq_len(p), function(j) {
    identical(sort(permutation[, j]), target)
  }, logical(1))
  if (!all(valid)) {
    stop(sprintf("Every column of permutation %d must permute 1:n.", index),
         call. = FALSE)
  }
  permutation
}


.c7sc_permutation_plan <- function(n, p, B, permutations, seed) {
  if (!is.null(permutations)) {
    if (!is.null(seed)) {
      stop("Supply either `permutations` or `seed`, not both.",
           call. = FALSE)
    }
    if (length(dim(permutations)) == 3L) {
      dimensions <- dim(permutations)
      if (!identical(dimensions[1:2], c(n, p))) {
        stop("A permutation array must have dimensions n by p by B.",
             call. = FALSE)
      }
      permutations <- lapply(seq_len(dimensions[3]), function(b) {
        matrix(permutations[, , b, drop = FALSE], nrow = n, ncol = p)
      })
    }
    if (!is.list(permutations) || !length(permutations)) {
      stop("`permutations` must be a non-empty list or n by p by B array.",
           call. = FALSE)
    }
    if (!is.null(B)) {
      B <- .c7sc_integer(B, "B")
      if (B != length(permutations)) {
        stop("`B` must equal the number of supplied permutations.",
             call. = FALSE)
      }
    }
    plan <- lapply(seq_along(permutations), function(b) {
      .c7sc_validate_permutation(permutations[[b]], n, p, b)
    })
    return(list(plan = plan, source = "supplied", seed = NULL))
  }
  if (is.null(B) || is.null(seed)) {
    stop("When `permutations` is absent, both `B` and `seed` are required.",
         call. = FALSE)
  }
  B <- .c7sc_integer(B, "B")
  seed <- .c7sc_integer(seed, "seed", 0L, .Machine$integer.max)
  had.seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had.seed) {
    old.seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  }
  on.exit({
    if (had.seed) {
      assign(".Random.seed", old.seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  set.seed(seed)
  plan <- lapply(seq_len(B), function(b) {
    vapply(seq_len(p), function(j) sample.int(n), integer(n))
  })
  list(plan = plan, source = "generated from explicit seed", seed = seed)
}


.c7sc_permute_columns <- function(x, permutation) {
  answer <- x
  for (j in seq_len(ncol(x))) {
    answer[, j] <- x[permutation[, j], j]
  }
  answer
}


.c7sc_tau_fit <- function(
    x, K, tau, init, first_index, tol, max_iter, spatial_max_iter,
    zero_tol, empty_action, empty_active, cycle_action, ties) {
  sparse_k_spatial_median(
    x = x, K = K, tau = tau, init = init,
    first_index = first_index, tol = tol, max_iter = max_iter,
    spatial_max_iter = spatial_max_iter, zero_tol = zero_tol,
    empty_action = empty_action, empty_active = empty_active,
    cycle_action = cycle_action, ties = ties,
    reset_excluded = "none", keep_distances = FALSE
  )
}


.c7sc_tau_objective <- function(x, fit, tol, spatial_max_iter, zero_tol) {
  active <- fit$active_set
  overall.active <- .c7sc_certified_median(
    x[, active, drop = FALSE], tol, spatial_max_iter, zero_tol
  )
  overall <- numeric(ncol(x))
  overall[active] <- overall.active
  geometry <- .c7sc_geometry(
    x, fit$algorithm_centers, fit$labels, active, overall
  )
  objective <- as.numeric(geometry$between_dispersion)
  if (!is.finite(objective) || objective <= 0) {
    stop("The active-subspace between-center objective must be finite and " %+%
           "strictly positive before taking its logarithm.",
         call. = FALSE)
  }
  list(value = objective, overall = overall.active, geometry = geometry)
}


#' Select the Sparse--SM threshold by the permutation Gap criterion
#'
#' Fits Sparse--SM over a finite threshold grid and maximizes
#' `log(O(tau)) - mean(log(O_perm(tau)))`, where every separation criterion is
#' evaluated in that fit's retained subspace. This is method-internal
#' calibration, not reproduction of a paper simulation.
#'
#' @inheritParams sparse_k_spatial_median
#' @param tau_grid Finite non-negative threshold grid.
#' @param B Number of reference permutations when `permutations` is absent.
#' @param permutations Optional list of `n` by `p` column-index permutation
#'   matrices, or an `n` by `p` by `B` array.
#' @param seed Required explicit RNG seed when permutations are generated.
#' @param selection_ties Tie rule for equal Gap values; only `"smallest"` is
#'   supported.
#' @param keep_fits Whether to retain every observed-data fit.
#' @param keep_permutations Whether to retain the complete permutation plan.
#'
#' @return A `sparse_sm_tau_selection` object containing the selected threshold,
#'   selected fit, and all observed/reference Gap ingredients.
#'
#' @references
#' Zhao, P., Zhuang, D., and Feng, L. (2026). Sparse K-spatial-median
#' clustering for high-dimensional data. arXiv:2605.00598.
#' \url{https://arxiv.org/abs/2605.00598}.
#'
#' @examples
#' \donttest{
#' x <- rbind(
#'   c(-4, 0, 2), c(-3, 1, 2), c(-2, 2, 2),
#'   c(2, 0.2, -2), c(3, 1.2, -2), c(4, 2.2, -2)
#' )
#' permutations <- list(
#'   cbind(6:1, c(3:6, 1:2), c(2:6, 1)),
#'   cbind(c(2:6, 1), 6:1, c(4:6, 1:3))
#' )
#' sparse_sm_select_tau(
#'   x, K = 2, tau_grid = c(0, 1), permutations = permutations,
#'   init = c(1, 6), tol = 1e-6, spatial_max_iter = 2000,
#'   empty_active = "largest"
#' )
#' }
#'
#' @export
sparse_sm_select_tau <- function(
    x, K, tau_grid, B = NULL, permutations = NULL, seed = NULL,
    init = "maxmin", first_index = 1L, tol = 1e-8, max_iter = 100L,
    spatial_max_iter = 500L, zero_tol = 0,
    empty_action = c("error", "farthest"),
    empty_active = c("error", "largest", "all"),
    cycle_action = c("error", "return"), ties = "first",
    selection_ties = "smallest", keep_fits = FALSE,
    keep_permutations = FALSE) {
  call <- match.call()
  x <- .as_data_matrix(x)
  K <- .c7sc_integer(K, "K", 2L, nrow(x))
  tau_grid <- .c7sc_validate_tau_grid(tau_grid)
  first_index <- .c7sc_integer(first_index, "first_index", 1L, nrow(x))
  max_iter <- .c7sc_integer(max_iter, "max_iter")
  controls <- .validate_iteration_controls(tol, spatial_max_iter, zero_tol)
  empty_action <- .c7sc_choice(
    empty_action, c("error", "farthest"), "empty_action"
  )
  empty_active <- .c7sc_choice(
    empty_active, c("error", "largest", "all"), "empty_active"
  )
  cycle_action <- .c7sc_choice(
    cycle_action, c("error", "return"), "cycle_action"
  )
  ties <- .c7sc_choice(ties, "first", "ties")
  selection_ties <- .c7sc_choice(
    selection_ties, "smallest", "selection_ties"
  )
  keep_fits <- .c7sc_flag(keep_fits, "keep_fits")
  keep_permutations <- .c7sc_flag(
    keep_permutations, "keep_permutations"
  )
  plan <- .c7sc_permutation_plan(
    nrow(x), ncol(x), B, permutations, seed
  )
  observed.fits <- vector("list", length(tau_grid))
  observed <- numeric(length(tau_grid))
  reference <- matrix(NA_real_, length(plan$plan), length(tau_grid))
  for (index in seq_along(tau_grid)) {
    tau <- tau_grid[index]
    observed.fits[[index]] <- .c7sc_tau_fit(
      x, K, tau, init, first_index, controls$tol, max_iter,
      controls$max_iter, controls$zero_tol, empty_action, empty_active,
      cycle_action, ties
    )
    if (!isTRUE(observed.fits[[index]]$diagnostics$converged)) {
      stop(sprintf("The observed Sparse--SM fit at tau=%s did not " %+%
                     "converge.", format(tau)), call. = FALSE)
    }
    observed[index] <- .c7sc_tau_objective(
      x, observed.fits[[index]], controls$tol,
      controls$max_iter, controls$zero_tol
    )$value
    for (b in seq_along(plan$plan)) {
      reference.x <- .c7sc_permute_columns(x, plan$plan[[b]])
      reference.fit <- .c7sc_tau_fit(
        reference.x, K, tau, init, first_index, controls$tol, max_iter,
        controls$max_iter, controls$zero_tol, empty_action, empty_active,
        cycle_action, ties
      )
      if (!isTRUE(reference.fit$diagnostics$converged)) {
        stop(sprintf("Reference fit b=%d at tau=%s did not converge.",
                     b, format(tau)), call. = FALSE)
      }
      reference[b, index] <- .c7sc_tau_objective(
        reference.x, reference.fit, controls$tol,
        controls$max_iter, controls$zero_tol
      )$value
    }
  }
  if (any(!is.finite(reference)) || any(reference <= 0)) {
    stop("Every reference objective must be finite and strictly positive.",
         call. = FALSE)
  }
  gap <- log(observed) - colMeans(log(reference))
  selected.index <- which(gap == max(gap))[1L]
  table <- data.frame(
    tau = tau_grid, objective = observed,
    mean_log_reference = colMeans(log(reference)), gap = gap,
    active_count = vapply(
      observed.fits, function(fit) length(fit$active_set), integer(1)
    ),
    stringsAsFactors = FALSE
  )
  structure(
    list(
      method = "Sparse--SM permutation Gap threshold selection",
      selected_tau = tau_grid[selected.index],
      selected_index = as.integer(selected.index),
      fit = observed.fits[[selected.index]], table = table,
      reference_objectives = reference,
      fits = if (keep_fits) observed.fits else NULL,
      permutations = if (keep_permutations) plan$plan else NULL,
      diagnostics = list(
        B = as.integer(length(plan$plan)),
        permutation.source = plan$source, seed = plan$seed,
        selection.tie.rule = "smallest tau",
        coordinate.convention = "retained active subspace",
        logarithm.contract = "all observed and reference O(tau) are > 0"
      ),
      call = call
    ),
    class = c("sparse_sm_tau_selection", "list")
  )
}


.c7sc_validate_k_grid <- function(k_grid, n) {
  values <- as.numeric(k_grid)
  if (!length(values) || anyNA(values) || any(!is.finite(values)) ||
      any(values != floor(values)) || any(values < 2L) || any(values >= n)) {
    stop("`k_grid` must contain integers satisfying 2 <= K < n.",
         call. = FALSE)
  }
  as.integer(sort(unique(values)))
}


.c7sc_retained_medians <- function(x, labels, active, K, tol,
                                    spatial_max_iter, zero_tol) {
  centers <- matrix(NA_real_, K, length(active))
  for (cluster in seq_len(K)) {
    rows <- which(labels == cluster)
    if (!length(rows)) {
      stop("A selected fit contains an empty cluster.", call. = FALSE)
    }
    centers[cluster, ] <- .c7sc_certified_median(
      x[rows, active, drop = FALSE], tol, spatial_max_iter, zero_tol
    )
  }
  centers
}


#' Select K for Sparse--SM by the BWDM rule
#'
#' Retunes `tau` for every candidate `K`, then evaluates average between-median
#' and within-median distances in the retained subspace. The selected `K`
#' maximizes the degree-of-freedom adjusted BWDM index.
#'
#' @inheritParams sparse_sm_select_tau
#' @param k_grid Candidate integers satisfying `2 <= K < n`.
#' @param keep_selections Whether to retain every candidate's threshold
#'   selection object.
#'
#' @return A `sparse_sm_k_selection` object with the selected `K`, selected
#'   fit, and exact ABDM, AWDM, and BWDM ingredients.
#'
#' @references
#' Zhao, P., Zhuang, D., and Feng, L. (2026). Sparse K-spatial-median
#' clustering for high-dimensional data. arXiv:2605.00598.
#' \url{https://arxiv.org/abs/2605.00598}.
#'
#' @examples
#' \donttest{
#' x <- rbind(
#'   c(-1, -0.3), c(0, 0.1), c(1, 0.4),
#'   c(5, -0.2), c(6, 0.2), c(7, 0.5)
#' )
#' permutations <- list(cbind(
#'   c(2:6, 1), c(4:6, 1:3)
#' ))
#' sparse_sm_select_k(
#'   x, k_grid = c(2, 3), tau_grid = 0, permutations = permutations,
#'   init = "maxmin", tol = 1e-6, spatial_max_iter = 2000,
#'   empty_action = "farthest"
#' )
#' }
#'
#' @export
sparse_sm_select_k <- function(
    x, k_grid, tau_grid, B = NULL, permutations = NULL, seed = NULL,
    init = "maxmin", first_index = 1L, tol = 1e-8, max_iter = 100L,
    spatial_max_iter = 500L, zero_tol = 0,
    empty_action = c("error", "farthest"),
    empty_active = c("error", "largest", "all"),
    cycle_action = c("error", "return"), ties = "first",
    selection_ties = "smallest", keep_selections = FALSE,
    keep_permutations = FALSE) {
  call <- match.call()
  x <- .as_data_matrix(x)
  k_grid <- .c7sc_validate_k_grid(k_grid, nrow(x))
  tau_grid <- .c7sc_validate_tau_grid(tau_grid)
  keep_selections <- .c7sc_flag(keep_selections, "keep_selections")
  keep_permutations <- .c7sc_flag(
    keep_permutations, "keep_permutations"
  )
  plan <- .c7sc_permutation_plan(
    nrow(x), ncol(x), B, permutations, seed
  )
  selections <- vector("list", length(k_grid))
  abdm <- awdm <- bwdm <- selected.tau <- numeric(length(k_grid))
  active.count <- integer(length(k_grid))
  for (index in seq_along(k_grid)) {
    K <- k_grid[index]
    selection <- sparse_sm_select_tau(
      x = x, K = K, tau_grid = tau_grid,
      B = length(plan$plan), permutations = plan$plan, seed = NULL,
      init = init, first_index = first_index, tol = tol,
      max_iter = max_iter, spatial_max_iter = spatial_max_iter,
      zero_tol = zero_tol, empty_action = empty_action,
      empty_active = empty_active, cycle_action = cycle_action,
      ties = ties, selection_ties = selection_ties,
      keep_fits = FALSE, keep_permutations = FALSE
    )
    fit <- selection$fit
    active <- fit$active_set
    retained.centers <- .c7sc_retained_medians(
      x, fit$labels, active, K, tol, spatial_max_iter, zero_tol
    )
    active.x <- x[, active, drop = FALSE]
    geometry <- cpp_ch7sc_geometry(
      active.x, retained.centers, as.integer(fit$labels),
      numeric(length(active)), as.integer(seq_along(active) - 1L)
    )
    abdm[index] <- as.numeric(geometry$average_between)
    awdm[index] <- as.numeric(geometry$within_mean)
    if (!is.finite(awdm[index]) || awdm[index] <= 0) {
      stop(sprintf("AWDM is not strictly positive for K=%d.", K),
           call. = FALSE)
    }
    bwdm[index] <- (abdm[index] / (K - 1)) /
      (awdm[index] / (nrow(x) - K))
    if (!is.finite(bwdm[index])) {
      stop(sprintf("BWDM is non-finite for K=%d.", K), call. = FALSE)
    }
    selected.tau[index] <- selection$selected_tau
    active.count[index] <- length(active)
    selections[[index]] <- selection
  }
  selected.index <- which(bwdm == max(bwdm))[1L]
  table <- data.frame(
    K = k_grid, tau = selected.tau, active_count = active.count,
    ABDM = abdm, AWDM = awdm, BWDM = bwdm,
    stringsAsFactors = FALSE
  )
  structure(
    list(
      method = "Sparse--SM BWDM cluster-count selection",
      selected_k = k_grid[selected.index],
      selected_index = as.integer(selected.index),
      fit = selections[[selected.index]]$fit, table = table,
      selections = if (keep_selections) selections else NULL,
      permutations = if (keep_permutations) plan$plan else NULL,
      diagnostics = list(
        B = as.integer(length(plan$plan)),
        permutation.source = plan$source, seed = plan$seed,
        selection.tie.rule = "smallest K",
        coordinate.convention = "retained active subspace",
        tau.reselected.for.each.K = TRUE,
        boundary.contract = "2 <= K < n and AWDM > 0"
      ),
      call = call
    ),
    class = c("sparse_sm_k_selection", "list")
  )
}
