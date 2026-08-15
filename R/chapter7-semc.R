# SEMC implementation for HDElliptical.
#
# Algorithm provenance: Feng and Zhuang (2026), arXiv:2605.08995, and the
# authors' MIT-licensed GEMcluster implementation at commit
# 10fce04fe690fe274dd5d237cfcd3d5c6a4139f6. This is an independent rewrite
# with every software-only tuning choice and numerical safeguard exposed.
# Upstream reference implementation copyright (c) 2026 Dan Zhuang and
# Long Feng. SPDX-License-Identifier: MIT


.semc_data <- function(x, name = "x", min_rows = 2L) {
  if (is.data.frame(x)) {
    x <- data.matrix(x)
  }
  if (!is.matrix(x) || !is.numeric(x) || length(dim(x)) != 2L) {
    stop(sprintf("`%s` must be a numeric matrix or data frame.", name),
         call. = FALSE)
  }
  storage.mode(x) <- "double"
  if (nrow(x) < min_rows || ncol(x) < 2L) {
    stop(sprintf("`%s` must have at least %d rows and two columns.",
                 name, min_rows), call. = FALSE)
  }
  if (anyNA(x) || any(!is.finite(x))) {
    stop(sprintf("`%s` must contain only finite values.", name),
         call. = FALSE)
  }
  x
}


.semc_scalar <- function(value, name, lower = -Inf, upper = Inf,
                         lower_open = FALSE, upper_open = FALSE) {
  value <- as.numeric(value)
  if (length(value) != 1L || is.na(value) || !is.finite(value) ||
      (if (lower_open) value <= lower else value < lower) ||
      (if (upper_open) value >= upper else value > upper)) {
    left <- if (lower_open) "(" else "["
    right <- if (upper_open) ")" else "]"
    stop(sprintf("`%s` must be one finite number in %s%s, %s%s.",
                 name, left, format(lower), format(upper), right),
         call. = FALSE)
  }
  value
}


.semc_integer <- function(value, name, lower = 1L, upper = .Machine$integer.max) {
  value <- as.numeric(value)
  if (length(value) != 1L || is.na(value) || !is.finite(value) ||
      value != floor(value) || value < lower || value > upper) {
    stop(sprintf("`%s` must be one integer in [%d, %d].",
                 name, as.integer(lower), as.integer(upper)), call. = FALSE)
  }
  as.integer(value)
}


.semc_flag <- function(value, name) {
  if (!is.logical(value) || length(value) != 1L || is.na(value)) {
    stop(sprintf("`%s` must be TRUE or FALSE.", name), call. = FALSE)
  }
  value
}


.semc_with_seed <- function(seed, code) {
  if (is.null(seed)) {
    return(code())
  }
  seed <- .semc_integer(seed, "seed", 0L)
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
  code()
}


.semc_spd <- function(matrix, name, spd_tol, symmetry_tol) {
  matrix <- as.matrix(matrix)
  if (!is.numeric(matrix) || nrow(matrix) == 0L ||
      nrow(matrix) != ncol(matrix) || anyNA(matrix) ||
      any(!is.finite(matrix))) {
    stop(sprintf("`%s` must be a finite square matrix.", name),
         call. = FALSE)
  }
  scale <- max(1, max(abs(matrix)))
  symmetry.error <- max(abs(matrix - t(matrix))) / scale
  if (!is.finite(symmetry.error) || symmetry.error > symmetry_tol) {
    stop(sprintf("`%s` is not symmetric within `symmetry_tol`.", name),
         call. = FALSE)
  }
  matrix <- (matrix + t(matrix)) / 2
  eigenvalues <- eigen(matrix, symmetric = TRUE, only.values = TRUE)$values
  minimum <- min(eigenvalues)
  factor <- tryCatch(chol(matrix), error = identity)
  if (inherits(factor, "condition") || !is.finite(minimum) ||
      minimum <= spd_tol) {
    stop(sprintf("`%s` is not strictly positive definite above `spd_tol`.",
                 name), call. = FALSE)
  }
  list(
    matrix = matrix, chol = factor, minimum.eigenvalue = minimum,
    symmetry.error = symmetry.error, condition.reciprocal = rcond(matrix)
  )
}


.semc_normalize_shape <- function(shape, name, spd_tol, symmetry_tol) {
  certificate <- .semc_spd(shape, name, spd_tol, symmetry_tol)
  trace <- sum(diag(certificate$matrix))
  if (!is.finite(trace) || trace <= 0) {
    stop(sprintf("`%s` has a nonpositive trace.", name), call. = FALSE)
  }
  p <- ncol(certificate$matrix)
  normalized <- certificate$matrix * (p / trace)
  normalized.certificate <- .semc_spd(
    normalized, paste0(name, " after trace normalization"),
    spd_tol, symmetry_tol
  )
  precision <- chol2inv(normalized.certificate$chol)
  precision <- (precision + t(precision)) / 2
  precision.certificate <- .semc_spd(
    precision, paste0(name, " inverse"), spd_tol, symmetry_tol
  )
  list(
    shape = normalized.certificate$matrix,
    precision = precision.certificate$matrix,
    shape.certificate = normalized.certificate,
    precision.certificate = precision.certificate,
    trace.before = trace,
    trace.after = sum(diag(normalized.certificate$matrix))
  )
}


.semc_l1_distances <- function(x, centers, active = seq_len(ncol(x))) {
  vapply(seq_len(nrow(centers)), function(k) {
    rowSums(abs(sweep(
      x[, active, drop = FALSE], 2L, centers[k, active], "-"
    )))
  }, numeric(nrow(x)))
}


.semc_coordinate_centers <- function(x, labels, K) {
  centers <- matrix(NA_real_, K, ncol(x))
  for (k in seq_len(K)) {
    rows <- which(labels == k)
    if (!length(rows)) {
      next
    }
    centers[k, ] <- apply(x[rows, , drop = FALSE], 2L, stats::median)
  }
  dimnames(centers) <- list(paste0("cluster", seq_len(K)), colnames(x))
  centers
}


.semc_maxmin_centers <- function(x, K, first_index = 1L) {
  n <- nrow(x)
  first_index <- .semc_integer(first_index, "first_index", 1L, n)
  selected <- integer(K)
  selected[1L] <- first_index
  if (K > 1L) {
    distance <- rowSums(abs(sweep(x, 2L, x[first_index, ], "-")))
    for (k in 2:K) {
      distance[selected[seq_len(k - 1L)]] <- -Inf
      selected[k] <- which.max(distance)
      candidate <- rowSums(abs(sweep(x, 2L, x[selected[k], ], "-")))
      distance <- pmin(distance, candidate)
    }
  }
  list(centers = x[selected, , drop = FALSE], indices = selected)
}


.semc_repair_empty <- function(labels, distances, action) {
  K <- ncol(distances)
  repairs <- integer()
  for (k in seq_len(K)) {
    if (any(labels == k)) {
      next
    }
    if (action == "error") {
      stop(sprintf("Initialization produced empty cluster %d.", k),
           call. = FALSE)
    }
    assigned.distance <- distances[cbind(seq_along(labels), labels)]
    donor.size <- tabulate(labels, nbins = K)[labels]
    candidates <- which(donor.size > 1L)
    if (!length(candidates)) {
      stop("No non-singleton cluster can supply an explicit empty repair.",
           call. = FALSE)
    }
    moved <- candidates[which.max(assigned.distance[candidates])]
    labels[moved] <- k
    repairs <- c(repairs, moved)
  }
  list(labels = labels, moved = repairs)
}


.semc_sparse_kmedian_once <- function(
    x, K, tau, centers, max_iter, empty_action, empty_active) {
  previous.labels <- NULL
  previous.active <- NULL
  objective.trace <- numeric()
  empty.repairs <- integer()
  converged <- FALSE
  active <- seq_len(ncol(x))
  labels <- rep.int(NA_integer_, nrow(x))
  signatures <- character()
  cycle <- FALSE

  for (iteration in seq_len(max_iter)) {
    center.mean <- colMeans(centers)
    scores <- colSums(abs(sweep(centers, 2L, center.mean, "-")))
    active <- which(scores >= tau)
    if (!length(active)) {
      if (empty_active == "error") {
        stop("The supplied initialization threshold excludes every feature.",
             call. = FALSE)
      }
      active <- seq_len(ncol(x))
    }
    distances <- .semc_l1_distances(x, centers, active)
    labels <- max.col(-distances, ties.method = "first")
    repair <- .semc_repair_empty(labels, distances, empty_action)
    labels <- repair$labels
    empty.repairs <- c(empty.repairs, repair$moved)
    centers.new <- .semc_coordinate_centers(x, labels, K)
    if (anyNA(centers.new)) {
      stop("An initialization center is undefined after empty handling.",
           call. = FALSE)
    }
    final.distances <- .semc_l1_distances(x, centers.new, active)
    objective <- sum(final.distances[cbind(seq_len(nrow(x)), labels)])
    objective.trace <- c(objective.trace, objective)
    signature <- paste(paste(labels, collapse = ","),
                       paste(active, collapse = ","), sep = "|")
    if (!is.null(previous.labels) && identical(labels, previous.labels) &&
        identical(active, previous.active)) {
      converged <- TRUE
      centers <- centers.new
      break
    }
    if (signature %in% signatures) {
      cycle <- TRUE
      centers <- centers.new
      break
    }
    signatures <- c(signatures, signature)
    previous.labels <- labels
    previous.active <- active
    centers <- centers.new
  }

  list(
    labels = as.integer(labels), centers = centers, active = active,
    objective = objective.trace[length(objective.trace)],
    objective.trace = objective.trace,
    feature.scores = scores, iterations = length(objective.trace),
    converged = converged, cycle.detected = cycle,
    empty.repairs = empty.repairs,
    certified = converged && !cycle && all(tabulate(labels, K) > 0L)
  )
}


.semc_initial_centers <- function(
    x, K, initialization, first_index, initial_labels, initial_centers) {
  if (initialization == "labels") {
    labels <- as.integer(initial_labels)
    if (length(labels) != nrow(x) || anyNA(labels) ||
        any(labels < 1L | labels > K) || any(tabulate(labels, K) == 0L)) {
      stop("`initial_labels` must use every integer label from 1 through K.",
           call. = FALSE)
    }
    return(list(
      centers = .semc_coordinate_centers(x, labels, K),
      labels = labels, source = "supplied labels", indices = NULL
    ))
  }
  if (initialization == "centers") {
    centers <- as.matrix(initial_centers)
    if (!is.numeric(centers) || !identical(dim(centers), c(K, ncol(x))) ||
        anyNA(centers) || any(!is.finite(centers))) {
      stop("`initial_centers` must be a finite K by p matrix.",
           call. = FALSE)
    }
    storage.mode(centers) <- "double"
    return(list(
      centers = centers, labels = NULL,
      source = "supplied centers", indices = NULL
    ))
  }
  if (initialization == "maxmin") {
    start <- .semc_maxmin_centers(x, K, first_index)
    return(list(
      centers = start$centers, labels = NULL,
      source = "deterministic max-min L1 centers", indices = start$indices
    ))
  }
  indices <- sample.int(nrow(x), K, replace = FALSE)
  list(
    centers = x[indices, , drop = FALSE], labels = NULL,
    source = "official-software random sparse-K-median start",
    indices = indices
  )
}


.semc_sparse_kmedian_multistart <- function(
    x, K, tau, nstart, max_iter, initialization, first_index,
    initial_labels, initial_centers, empty_action, empty_active) {
  fits <- vector("list", nstart)
  summaries <- vector("list", nstart)
  best <- NULL
  for (start.index in seq_len(nstart)) {
    start.method <- initialization
    if (start.index > 1L && initialization %in% c("labels", "centers")) {
      stop("Supplied labels or centers require `init_nstart = 1`.",
           call. = FALSE)
    }
    start <- .semc_initial_centers(
      x, K, start.method,
      ((first_index + start.index - 2L) %% nrow(x)) + 1L,
      initial_labels, initial_centers
    )
    fit <- tryCatch(
      .semc_sparse_kmedian_once(
        x, K, tau, start$centers, max_iter,
        empty_action, empty_active
      ),
      error = identity
    )
    if (inherits(fit, "condition")) {
      summaries[[start.index]] <- data.frame(
        start = start.index, objective = NA_real_, converged = FALSE,
        certified = FALSE, failure = conditionMessage(fit)
      )
      next
    }
    fit$start <- start
    fits[[start.index]] <- fit
    summaries[[start.index]] <- data.frame(
      start = start.index, objective = fit$objective,
      converged = fit$converged, certified = fit$certified,
      failure = NA_character_
    )
    if (isTRUE(fit$certified) &&
        (is.null(best) || fit$objective < best$objective)) {
      best <- fit
      best$best.start <- start.index
    }
  }
  if (is.null(best)) {
    failures <- vapply(summaries, function(item) item$failure, character(1L))
    stop(paste(
      "No sparse-K-median initialization start was certified:",
      paste(stats::na.omit(failures), collapse = "; ")
    ), call. = FALSE)
  }
  best$start.summaries <- do.call(rbind, summaries)
  best
}


.semc_permute_columns <- function(x) {
  out <- x
  for (j in seq_len(ncol(x))) {
    out[, j] <- x[sample.int(nrow(x)), j]
  }
  out
}


.semc_between_l1 <- function(x, labels, active) {
  overall <- apply(x[, active, drop = FALSE], 2L, stats::median)
  centers <- .semc_coordinate_centers(x, labels, max(labels))[, active,
                                                                    drop = FALSE]
  sum(tabulate(labels, max(labels)) * rowSums(abs(sweep(
    centers, 2L, overall, "-"
  ))))
}


.semc_select_initial_tau <- function(
    x, K, tau_grid, B, nstart, max_iter, initialization, first_index,
    empty_action, empty_active, dispersion_floor) {
  pilot <- .semc_sparse_kmedian_multistart(
    x, K, 0, max(1L, min(nstart, 2L)), max_iter,
    initialization, first_index, NULL, NULL, empty_action, empty_active
  )
  if (is.null(tau_grid)) {
    center.mean <- colMeans(pilot$centers)
    scores <- colSums(abs(sweep(pilot$centers, 2L, center.mean, "-")))
    tau_grid <- unique(as.numeric(stats::quantile(
      scores, probs = seq(0.2, 0.8, length.out = 6L), names = FALSE
    )))
    tau_grid <- tau_grid[is.finite(tau_grid) & tau_grid >= 0]
  }
  if (!is.numeric(tau_grid) || !length(tau_grid) || anyNA(tau_grid) ||
      any(!is.finite(tau_grid)) || any(tau_grid < 0)) {
    stop("`init_tau_grid` must contain finite nonnegative thresholds.",
         call. = FALSE)
  }
  tau_grid <- sort(unique(as.numeric(tau_grid)))
  gap <- numeric(length(tau_grid))
  observed <- numeric(length(tau_grid))
  reference <- matrix(NA_real_, length(tau_grid), B)
  fits <- vector("list", length(tau_grid))
  permutations <- lapply(seq_len(B), function(unused) .semc_permute_columns(x))
  for (index in seq_along(tau_grid)) {
    tau <- tau_grid[index]
    fit <- .semc_sparse_kmedian_multistart(
      x, K, tau, nstart, max_iter, initialization, first_index,
      NULL, NULL, empty_action, empty_active
    )
    fits[[index]] <- fit
    observed[index] <- .semc_between_l1(x, fit$labels, fit$active)
    for (b in seq_len(B)) {
      fit.reference <- .semc_sparse_kmedian_multistart(
        permutations[[b]], K, tau, max(1L, ceiling(nstart / 2)),
        max_iter, initialization, first_index,
        NULL, NULL, empty_action, empty_active
      )
      reference[index, b] <- .semc_between_l1(
        permutations[[b]], fit.reference$labels, fit.reference$active
      )
    }
    if (observed[index] <= dispersion_floor ||
        any(reference[index, ] <= dispersion_floor)) {
      stop(paste0(
        "Initialization Gap dispersion is not above the explicit ",
        "`init_dispersion_floor`."
      ), call. = FALSE)
    }
    gap[index] <- log(observed[index]) - mean(log(reference[index, ]))
  }
  selected <- which.max(gap)
  list(
    tau = tau_grid[selected], selected.index = selected,
    tau.grid = tau_grid, gap = gap, observed = observed,
    reference = reference, selected.fit = fits[[selected]],
    contract = paste(
      "official-software sparse-K-median quantile grid and",
      "column-permutation initialization Gap"
    )
  )
}


# Generator and shape helpers -------------------------------------------------

.semc_trapezoid <- function(x, y) {
  if (length(x) != length(y) || length(x) < 2L) {
    stop("Trapezoid inputs must have a common length of at least two.",
         call. = FALSE)
  }
  index <- 2:length(x)
  sum((x[index] - x[index - 1L]) *
        (y[index] + y[index - 1L]) / 2)
}


.semc_weighted_variance <- function(x, weights) {
  total <- sum(weights)
  mean <- sum(weights * x) / total
  sum(weights * (x - mean)^2) / total
}


.semc_build_generator <- function(
    delta, posterior, p, bandwidth, bandwidth_min, grid_size,
    radius_floor, density_floor, score_clip, spline_spar,
    extrapolation) {
  if (!identical(dim(delta), dim(posterior))) {
    stop("delta and posterior must have identical dimensions.",
         call. = FALSE)
  }
  if (any(!is.finite(delta)) || any(delta < 0) ||
      any(!is.finite(posterior)) || any(posterior < 0)) {
    stop("Generator inputs must be finite and nonnegative.", call. = FALSE)
  }
  y.observed <- log1p(as.vector(delta))
  weights <- as.vector(posterior)
  weight.sum <- sum(weights)
  if (!is.finite(weight.sum) || weight.sum <= 0) {
    stop("Generator weights must have a strictly positive sum.",
         call. = FALSE)
  }
  weights <- weights / weight.sum
  effective.n <- 1 / sum(weights^2)
  bandwidth.source <- "supplied"
  if (is.null(bandwidth)) {
    sd.y <- sqrt(.semc_weighted_variance(y.observed, weights))
    if (!is.finite(sd.y)) {
      stop("The weighted transformed-radius scale is not finite.",
           call. = FALSE)
    }
    bandwidth <- 1.06 * sd.y * effective.n^(-1 / 5)
    bandwidth.source <- "official-software weighted Silverman rule"
  }
  bandwidth.raw <- bandwidth
  bandwidth <- max(bandwidth, bandwidth_min)
  y.minimum <- max(min(y.observed) - 3 * bandwidth,
                   log1p(radius_floor))
  y.maximum <- max(y.observed) + 3 * bandwidth
  if (!is.finite(y.minimum) || !is.finite(y.maximum) ||
      y.maximum <= y.minimum) {
    stop("The transformed-radius generator grid is degenerate.",
         call. = FALSE)
  }
  y.grid <- seq(y.minimum, y.maximum, length.out = grid_size)
  u.grid <- expm1(y.grid)
  density <- cpp_ch7_semc_weighted_kde(
    y.grid, y.observed, weights, bandwidth
  )
  density.floor.uses <- sum(density < density_floor)
  log.q <- log(pmax(density, density_floor)) - log1p(u.grid)
  log.g.raw <- (1 - p / 2) * log(u.grid) + log.q
  spline <- tryCatch(
    stats::smooth.spline(y.grid, log.g.raw, spar = spline_spar),
    error = identity
  )
  if (inherits(spline, "condition")) {
    stop(paste(
      "The explicitly requested generator smoothing spline failed:",
      conditionMessage(spline)
    ), call. = FALSE)
  }
  log.g.grid <- stats::predict(spline, x = y.grid, deriv = 0)$y
  derivative <- stats::predict(spline, x = y.grid, deriv = 1)$y

  # Enforce the paper normalization after smoothing. This is a scalar
  # identification, not an eigenvalue or density repair.
  log.integrand <- (p / 2 - 1) * log(u.grid) + log.g.grid
  log.scale <- max(log.integrand)
  integral.scaled <- .semc_trapezoid(
    u.grid, exp(log.integrand - log.scale)
  )
  if (!is.finite(integral.scaled) || integral.scaled <= 0) {
    stop("The smoothed generator has a nonpositive radial integral.",
         call. = FALSE)
  }
  log.surface <- (p / 2) * log(pi) - lgamma(p / 2)
  log.normalizer <- log.surface + log.scale + log(integral.scaled)
  log.g.grid <- log.g.grid - log.normalizer
  score.unclipped <- -derivative / (1 + u.grid)
  score.grid <- pmin(pmax(score.unclipped, score_clip[1L]), score_clip[2L])
  score.clip.uses <- sum(score.unclipped < score_clip[1L] |
                           score.unclipped > score_clip[2L])
  list(
    y.grid = y.grid, u.grid = u.grid, log.g.grid = log.g.grid,
    score.grid = score.grid, bandwidth = bandwidth,
    bandwidth.raw = bandwidth.raw, bandwidth.source = bandwidth.source,
    effective.sample.size = effective.n,
    density.floor.uses = density.floor.uses,
    score.clip.uses = score.clip.uses,
    controls = list(
      bandwidth.min = bandwidth_min, grid.size = grid_size,
      radius.floor = radius_floor, density.floor = density_floor,
      score.clip = score_clip, spline.spar = spline_spar,
      extrapolation = extrapolation,
      normalization = paste(
        "pi^(p/2)/Gamma(p/2) times integral",
        "u^(p/2-1) g(u) du equals one"
      )
    )
  )
}


.semc_generator_evaluate <- function(generator, delta,
                                     type = c("log_g", "score")) {
  type <- match.arg(type)
  if (any(!is.finite(delta)) || any(delta < 0)) {
    stop("Generator evaluation radii must be finite and nonnegative.",
         call. = FALSE)
  }
  y <- log1p(pmax(delta, generator$controls$radius.floor))
  outside <- y < generator$y.grid[1L] |
    y > generator$y.grid[length(generator$y.grid)]
  if (any(outside) && generator$controls$extrapolation == "error") {
    stop("A transformed radius lies outside the fitted generator grid.",
         call. = FALSE)
  }
  values <- if (type == "log_g") generator$log.g.grid else
    generator$score.grid
  stats::approx(
    generator$y.grid, values, xout = y,
    rule = if (generator$controls$extrapolation == "constant") 2L else 1L,
    ties = "ordered"
  )$y
}


.semc_select_factor_count <- function(scatter, maximum, ratio_floor) {
  values <- sort(eigen(
    (scatter + t(scatter)) / 2, symmetric = TRUE, only.values = TRUE
  )$values, decreasing = TRUE)
  values <- pmax(values, ratio_floor)
  p <- length(values)
  if (p < 4L) {
    return(list(factors = 0L, ratios = numeric(), values = values))
  }
  maximum <- min(maximum, p - 2L)
  if (maximum < 1L) {
    return(list(factors = 0L, ratios = numeric(), values = values))
  }
  q1 <- p - 1L
  residual.sum <- vapply(0:(q1 - 1L), function(j) {
    if (j == 0L) sum(values[seq_len(q1)]) else
      sum(values[(j + 1L):q1])
  }, numeric(1L))
  ratios <- vapply(seq_len(maximum), function(j) {
    numerator <- log1p(values[j] / max(residual.sum[j], ratio_floor))
    denominator <- log1p(
      values[j + 1L] / max(residual.sum[j + 1L], ratio_floor)
    )
    numerator / max(denominator, ratio_floor)
  }, numeric(1L))
  list(factors = as.integer(which.max(ratios)), ratios = ratios,
       values = values)
}


.semc_poet <- function(
    scatter, factors, threshold, ridge, spd_tol, symmetry_tol, stage) {
  scatter <- (scatter + t(scatter)) / 2
  p <- ncol(scatter)
  eig <- eigen(scatter, symmetric = TRUE)
  order <- order(eig$values, decreasing = TRUE)
  values <- eig$values[order]
  vectors <- eig$vectors[, order, drop = FALSE]
  if (factors > 0L && any(values[seq_len(factors)] < 0)) {
    stop(sprintf("%s selects a negative leading eigenvalue.", stage),
         call. = FALSE)
  }
  if (factors > 0L) {
    leading.vectors <- vectors[, seq_len(factors), drop = FALSE]
    leading.values <- values[seq_len(factors)]
    low.rank <- tcrossprod(sweep(
      leading.vectors, 2L, sqrt(leading.values), "*"
    ))
  } else {
    leading.vectors <- matrix(numeric(), p, 0L)
    leading.values <- numeric()
    low.rank <- matrix(0, p, p)
  }
  residual <- scatter - low.rank
  thresholded <- residual
  off <- row(thresholded) != col(thresholded)
  thresholded[off] <- sign(thresholded[off]) *
    pmax(abs(thresholded[off]) - threshold, 0)
  diag(thresholded) <- diag(residual)
  estimate <- low.rank + thresholded + diag(ridge, p)
  normalized <- .semc_normalize_shape(
    estimate, stage, spd_tol, symmetry_tol
  )
  list(
    shape = normalized$shape, precision = normalized$precision,
    factors = factors, threshold = threshold, explicit.ridge = ridge,
    leading.eigenvalues = leading.values,
    leading.eigenvectors = leading.vectors,
    low.rank = low.rank, idiosyncratic.raw = residual,
    idiosyncratic.thresholded = thresholded,
    certificate = list(
      minimum.eigenvalue =
        normalized$shape.certificate$minimum.eigenvalue,
      trace = normalized$trace.after,
      positive.definiteness.repair = "none",
      ridge = ridge
    )
  )
}



.semc_glasso_grid <- function(
    scatter, lambda_grid, effective_n, ebic_gamma,
    tolerance, max_iterations, initial_step, max_backtracking,
    spd_tol, majorization_tolerance) {
  p <- ncol(scatter)
  fits <- vector("list", length(lambda_grid))
  table <- data.frame(
    lambda = lambda_grid, converged = FALSE,
    kkt = NA_real_, minimum_eigenvalue = NA_real_,
    edges = NA_integer_, ebic = Inf, failure = NA_character_
  )
  for (index in seq_along(lambda_grid)) {
    fit <- tryCatch(
      cpp_ch7_semc_offdiag_glasso(
        scatter, lambda_grid[index], tolerance, max_iterations,
        initial_step, max_backtracking, spd_tol,
        majorization_tolerance
      ),
      error = identity
    )
    if (inherits(fit, "condition")) {
      table$failure[index] <- conditionMessage(fit)
      next
    }
    certified <- isTRUE(fit$converged) &&
      !isTRUE(fit$backtracking_failed) &&
      is.finite(fit$minimum_eigenvalue) &&
      fit$minimum_eigenvalue > spd_tol &&
      is.finite(fit$kkt_residual) && fit$kkt_residual <= tolerance
    table$converged[index] <- certified
    table$kkt[index] <- fit$kkt_residual
    table$minimum_eigenvalue[index] <- fit$minimum_eigenvalue
    if (!certified) {
      table$failure[index] <- "SPD/KKT/convergence certificate failed"
      fits[[index]] <- fit
      next
    }
    precision <- as.matrix(fit$solution)
    edges <- sum(upper.tri(precision) & precision != 0)
    degrees <- p + edges
    loglik <- as.numeric(determinant(
      precision, logarithm = TRUE
    )$modulus) - sum(scatter * precision)
    ebic <- -effective_n * loglik + log(effective_n) * degrees +
      4 * ebic_gamma * log(p) * edges
    table$edges[index] <- edges
    table$ebic[index] <- ebic
    fits[[index]] <- fit
  }
  candidates <- which(table$converged & is.finite(table$ebic))
  if (!length(candidates)) {
    stop(paste0(
      "No graphical-lasso grid point passed every SPD, convergence, ",
      "backtracking, and KKT certificate."
    ), call. = FALSE)
  }
  selected <- candidates[which.min(table$ebic[candidates])]
  list(
    fit = fits[[selected]], selected.index = selected,
    selected.lambda = lambda_grid[selected], table = table,
    selection = if (length(lambda_grid) == 1L) "supplied singleton" else
      "official-software EBIC grid, independently certified",
    diagonal.penalty = FALSE
  )
}


.semc_shape_update <- function(
    x, centers, posterior, current_precision, shape_method,
    eta_precision, radial_floor, poet_factor_selection,
    poet_factors, poet_max_factors, factor_ratio_floor,
    poet_threshold, poet_threshold_scale, poet_ridge,
    tyler_ridge, tyler_tol, tyler_max_iter,
    glasso_lambda, glasso_lambda_grid, glasso_lambda_scale,
    glasso_ebic_gamma, glasso_tol, glasso_max_iter,
    glasso_initial_step, glasso_max_backtracking,
    glasso_majorization_tol, spd_tol, symmetry_tol) {
  K <- nrow(centers)
  residuals <- do.call(rbind, lapply(seq_len(K), function(k) {
    sweep(x, 2L, centers[k, ], "-")
  }))
  weights <- as.vector(posterior)
  effective.n <- sum(weights)^2 / sum(weights^2)
  sign.fit <- cpp_ch7_semc_weighted_sign_scatter(
    residuals, weights, radial_floor
  )
  sign.scatter <- as.matrix(sign.fit$estimate)
  factor.fit <- if (poet_factor_selection == "software_gr") {
    .semc_select_factor_count(
      sign.scatter, poet_max_factors, factor_ratio_floor
    )
  } else {
    list(
      factors = poet_factors, ratios = NULL,
      values = eigen(sign.scatter, symmetric = TRUE,
                     only.values = TRUE)$values
    )
  }
  factors <- factor.fit$factors
  threshold <- if (is.null(poet_threshold)) {
    poet_threshold_scale * sqrt(log(ncol(x)) / effective.n)
  } else {
    poet_threshold
  }
  first.poet <- .semc_poet(
    sign.scatter, factors, threshold, poet_ridge,
    spd_tol, symmetry_tol, "first SEMC POET shape"
  )
  tyler <- cpp_ch7_semc_weighted_tyler(
    residuals, weights, first.poet$shape,
    tyler_ridge, radial_floor, tyler_tol, tyler_max_iter, spd_tol
  )
  if (!isTRUE(tyler$converged) ||
      tyler$relative_change > tyler_tol ||
      tyler$minimum_eigenvalue <= spd_tol) {
    stop("The weighted Tyler iteration failed its convergence/SPD certificate.",
         call. = FALSE)
  }
  tyler.normalized <- .semc_normalize_shape(
    as.matrix(tyler$shape), "weighted Tyler shape",
    spd_tol, symmetry_tol
  )
  second.poet <- NULL
  glasso <- NULL
  if (shape_method == "tyler") {
    proposal.shape <- tyler.normalized$shape
  } else {
    second.poet <- .semc_poet(
      tyler.normalized$shape, factors, threshold, poet_ridge,
      spd_tol, symmetry_tol, "second SEMC POET shape"
    )
    proposal.shape <- second.poet$shape
    if (shape_method == "glasso") {
      if (is.null(glasso_lambda_grid)) {
        base <- if (is.null(glasso_lambda)) {
          glasso_lambda_scale * sqrt(log(ncol(x)) / effective.n)
        } else {
          glasso_lambda
        }
        glasso_lambda_grid <- if (is.null(glasso_lambda)) {
          base * exp(seq(log(1.8), log(0.45), length.out = 8L))
        } else {
          base
        }
      }
      glasso_lambda_grid <- sort(
        unique(as.numeric(glasso_lambda_grid)), decreasing = TRUE
      )
      glasso <- .semc_glasso_grid(
        second.poet$shape, glasso_lambda_grid, effective.n,
        glasso_ebic_gamma, glasso_tol, glasso_max_iter,
        glasso_initial_step, glasso_max_backtracking,
        spd_tol, glasso_majorization_tol
      )
      proposal.precision <- as.matrix(glasso$fit$solution)
      proposal.certificate <- .semc_spd(
        proposal.precision, "SEMC graphical-lasso precision",
        spd_tol, symmetry_tol
      )
      proposal.shape <- chol2inv(proposal.certificate$chol)
    }
  }
  proposal.normalized <- .semc_normalize_shape(
    proposal.shape, paste0("SEMC ", shape_method, " proposal shape"),
    spd_tol, symmetry_tol
  )
  proposal.precision <- proposal.normalized$precision
  if (is.null(current_precision)) {
    damped.precision <- proposal.precision
  } else {
    damped.precision <-
      (1 - eta_precision) * current_precision +
      eta_precision * proposal.precision
  }
  damped.certificate <- .semc_spd(
    damped.precision, "damped SEMC precision", spd_tol, symmetry_tol
  )
  damped.shape <- chol2inv(damped.certificate$chol)
  final <- .semc_normalize_shape(
    damped.shape, "damped SEMC shape", spd_tol, symmetry_tol
  )
  list(
    shape = final$shape, precision = final$precision,
    effective.sample.size = effective.n,
    sign = sign.fit, first.poet = first.poet,
    factor.selection = c(factor.fit, list(
      source = poet_factor_selection,
      ratio.floor = if (poet_factor_selection == "software_gr")
        factor_ratio_floor else NULL
    )),
    tyler = tyler, second.poet = second.poet, glasso = glasso,
    controls = list(
      method = shape_method, nested.pipeline = switch(
        shape_method,
        tyler = "sign pilot -> POET initialization -> weighted Tyler",
        poet = "sign pilot -> POET -> weighted Tyler -> POET",
        glasso = paste(
          "sign pilot -> POET -> weighted Tyler -> POET ->",
          "off-diagonal graphical lasso"
        )
      ),
      poet.threshold = threshold,
      poet.threshold.source = if (is.null(poet_threshold))
        "official-software rate" else "supplied",
      poet.ridge = poet_ridge, tyler.ridge = tyler_ridge,
      eta.precision = eta_precision,
      no.hidden.repair = paste(
        "No eigenvalue projection, pseudoinverse, unreported ridge,",
        "or graphical-lasso fallback"
      )
    ),
    certificate = list(
      minimum.shape.eigenvalue =
        final$shape.certificate$minimum.eigenvalue,
      minimum.precision.eigenvalue =
        final$precision.certificate$minimum.eigenvalue,
      trace.shape = sum(diag(final$shape)),
      symmetry.error = final$shape.certificate$symmetry.error,
      valid = TRUE
    )
  )
}



.semc_delta <- function(x, centers, precision) {
  delta <- cpp_ch7_semc_delta(x, centers, precision)
  if (any(!is.finite(delta)) || any(delta < 0)) {
    stop(paste0(
      "SEMC quadratic radii are not finite and nonnegative; no absolute-value ",
      "or zero-floor repair is applied."
    ), call. = FALSE)
  }
  delta
}


.semc_posterior <- function(delta, mixing, generator) {
  if (any(!is.finite(mixing)) || any(mixing <= 0) ||
      abs(sum(mixing) - 1) > 1e-10) {
    stop("SEMC mixing weights must be strictly positive and sum to one.",
         call. = FALSE)
  }
  log.g <- matrix(
    .semc_generator_evaluate(generator, as.vector(delta), "log_g"),
    nrow = nrow(delta), ncol = ncol(delta)
  )
  log.scores <- sweep(log.g, 2L, log(mixing), "+")
  softmax <- cpp_ch7_semc_softmax(log.scores)
  list(
    posterior = as.matrix(softmax$probabilities),
    log.normalizers = as.numeric(softmax$log_normalizers),
    log.scores = log.scores
  )
}


.semc_fit_single <- function(
    x, K, initialization_fit, shape_method,
    eta_mu, eta_precision, mixing_floor, center_weight_floor,
    max_iter, convergence_tol, bandwidth, bandwidth_min,
    generator_grid_size, generator_radius_floor, generator_density_floor,
    generator_score_clip, generator_spline_spar,
    generator_extrapolation, radial_floor,
    poet_factor_selection, poet_factors, poet_max_factors,
    factor_ratio_floor, poet_threshold, poet_threshold_scale, poet_ridge,
    tyler_ridge, tyler_tol, tyler_max_iter,
    glasso_lambda, glasso_lambda_grid, glasso_lambda_scale,
    glasso_ebic_gamma, glasso_tol, glasso_max_iter,
    glasso_initial_step, glasso_max_backtracking,
    glasso_majorization_tol, spd_tol, symmetry_tol, keep_path) {
  n <- nrow(x)
  centers <- initialization_fit$centers
  labels <- initialization_fit$labels
  posterior <- matrix(0, n, K)
  posterior[cbind(seq_len(n), labels)] <- 1
  mixing <- colMeans(posterior)
  if (any(mixing <= 0)) {
    stop("Initialization has an empty component.", call. = FALSE)
  }

  shape.fit <- .semc_shape_update(
    x, centers, posterior, NULL, shape_method,
    1, radial_floor, poet_factor_selection,
    poet_factors, poet_max_factors, factor_ratio_floor,
    poet_threshold, poet_threshold_scale, poet_ridge,
    tyler_ridge, tyler_tol, tyler_max_iter,
    glasso_lambda, glasso_lambda_grid, glasso_lambda_scale,
    glasso_ebic_gamma, glasso_tol, glasso_max_iter,
    glasso_initial_step, glasso_max_backtracking,
    glasso_majorization_tol, spd_tol, symmetry_tol
  )
  shape <- shape.fit$shape
  precision <- shape.fit$precision
  delta <- .semc_delta(x, centers, precision)
  generator <- .semc_build_generator(
    delta, posterior, ncol(x), bandwidth, bandwidth_min,
    generator_grid_size, generator_radius_floor,
    generator_density_floor, generator_score_clip,
    generator_spline_spar, generator_extrapolation
  )
  initial.posterior <- .semc_posterior(delta, mixing, generator)
  objective <- sum(initial.posterior$log.normalizers)
  path <- data.frame(
    iteration = 0L, objective = objective,
    center.change = NA_real_, precision.change = NA_real_,
    mixing.change = NA_real_, mixing.floor.uses = 0L,
    shape.minimum.eigenvalue =
      shape.fit$certificate$minimum.shape.eigenvalue,
    precision.minimum.eigenvalue =
      shape.fit$certificate$minimum.precision.eigenvalue,
    tyler.iterations = shape.fit$tyler$iterations,
    glasso.kkt = if (is.null(shape.fit$glasso)) NA_real_ else
      shape.fit$glasso$fit$kkt_residual
  )
  converged <- FALSE
  last.shape.fit <- shape.fit

  for (iteration in seq_len(max_iter)) {
    old.centers <- centers
    old.precision <- precision
    old.mixing <- mixing

    delta <- .semc_delta(x, centers, precision)
    e.step <- .semc_posterior(delta, mixing, generator)
    posterior.new <- e.step$posterior
    mixing.raw <- colMeans(posterior.new)
    mixing.floor.uses <- sum(mixing.raw < mixing_floor)
    mixing.new <- pmax(mixing.raw, mixing_floor)
    mixing.new <- mixing.new / sum(mixing.new)

    generator.new <- .semc_build_generator(
      delta, posterior.new, ncol(x), bandwidth, bandwidth_min,
      generator_grid_size, generator_radius_floor,
      generator_density_floor, generator_score_clip,
      generator_spline_spar, generator_extrapolation
    )
    score <- matrix(
      .semc_generator_evaluate(
        generator.new, as.vector(delta), type = "score"
      ),
      nrow = n, ncol = K
    )
    proposed.centers <- centers
    for (k in seq_len(K)) {
      center.weights <- posterior.new[, k] * score[, k]
      denominator <- sum(center.weights)
      if (!is.finite(denominator) || denominator <= center_weight_floor) {
        stop(sprintf(
          paste0(
            "Component %d center weight is not above the explicit ",
            "center_weight_floor."
          ), k
        ), call. = FALSE)
      }
      proposed.centers[k, ] <-
        colSums(x * center.weights) / denominator
    }
    centers.new <- (1 - eta_mu) * centers + eta_mu * proposed.centers

    shape.new <- .semc_shape_update(
      x, centers.new, posterior.new, precision, shape_method,
      eta_precision, radial_floor, poet_factor_selection,
      poet_factors, poet_max_factors, factor_ratio_floor,
      poet_threshold, poet_threshold_scale, poet_ridge,
      tyler_ridge, tyler_tol, tyler_max_iter,
      glasso_lambda, glasso_lambda_grid, glasso_lambda_scale,
      glasso_ebic_gamma, glasso_tol, glasso_max_iter,
      glasso_initial_step, glasso_max_backtracking,
      glasso_majorization_tol, spd_tol, symmetry_tol
    )
    centers <- centers.new
    shape <- shape.new$shape
    precision <- shape.new$precision
    posterior <- posterior.new
    mixing <- mixing.new
    generator <- generator.new
    last.shape.fit <- shape.new

    center.change <- max(abs(centers - old.centers))
    precision.change <- norm(
      precision - old.precision, type = "F"
    ) / max(1, norm(old.precision, type = "F"))
    mixing.change <- max(abs(mixing - old.mixing))
    delta.updated <- .semc_delta(x, centers, precision)
    objective.updated <- sum(.semc_posterior(
      delta.updated, mixing, generator
    )$log.normalizers)
    path <- rbind(path, data.frame(
      iteration = iteration, objective = objective.updated,
      center.change = center.change,
      precision.change = precision.change,
      mixing.change = mixing.change,
      mixing.floor.uses = mixing.floor.uses,
      shape.minimum.eigenvalue =
        shape.new$certificate$minimum.shape.eigenvalue,
      precision.minimum.eigenvalue =
        shape.new$certificate$minimum.precision.eigenvalue,
      tyler.iterations = shape.new$tyler$iterations,
      glasso.kkt = if (is.null(shape.new$glasso)) NA_real_ else
        shape.new$glasso$fit$kkt_residual
    ))
    objective <- objective.updated
    if (max(center.change, precision.change, mixing.change) <=
        convergence_tol) {
      converged <- TRUE
      break
    }
  }

  delta.final <- .semc_delta(x, centers, precision)
  generator.final <- .semc_build_generator(
    delta.final, posterior, ncol(x), bandwidth, bandwidth_min,
    generator_grid_size, generator_radius_floor,
    generator_density_floor, generator_score_clip,
    generator_spline_spar, generator_extrapolation
  )
  final.e.step <- .semc_posterior(delta.final, mixing, generator.final)
  posterior.final <- final.e.step$posterior
  cluster <- max.col(posterior.final, ties.method = "first")
  objective.final <- sum(final.e.step$log.normalizers)
  list(
    cluster = as.integer(cluster), posterior = posterior.final,
    mixing = mixing, centers = centers, shape = shape,
    precision = precision, generator = generator.final,
    delta = delta.final, objective = objective.final,
    iterations = iteration, converged = converged,
    valid = converged && isTRUE(last.shape.fit$certificate$valid),
    initialization = initialization_fit,
    shape.fit = last.shape.fit,
    path = if (keep_path) path else NULL,
    diagnostics = list(
      final.posterior.row.error = max(abs(rowSums(posterior.final) - 1)),
      final.mixing.residual = max(abs(colMeans(posterior.final) - mixing)),
      generator.final.reestimated = TRUE,
      no.hidden.repair = paste(
        "All radial/density/score/ridge/floor controls are explicit;",
        "no eigenvalue projection, pseudoinverse, or glasso fallback"
      )
    )
  )
}



#' Semiparametric elliptical mixture clustering
#'
#' Fits the Feng--Zhuang semiparametric generalized-EM model with common radial
#' generator and common trace-normalized shape. The nested shape pipelines are
#' "tyler" (sign pilot, POET initialization, weighted Tyler), "poet" (the
#' preceding pipeline plus a second POET step), and "glasso" (the full paper
#' pipeline plus off-diagonal graphical lasso).
#'
#' The paper does not uniquely determine initialization, bandwidth constants,
#' POET factor selection, penalty constants, generator clipping, or all
#' stopping rules. Defaults identified below as software contracts reproduce
#' the named choices at the pinned official implementation, while every
#' numerical floor, ridge, damping constant, and fallback policy is an explicit
#' argument. There is no positive-definite projection, pseudoinverse, hidden
#' ridge, or graphical-lasso fallback.
#'
#' @param x Numeric matrix with observations in rows.
#' @param K Integer number of mixture components, at least two and less than
#'   the number of observations.
#' @param shape One of "glasso", "poet", or "tyler".
#' @param initialization Initial sparse-K-median start: deterministic "maxmin",
#'   official-software-style "software_random", or supplied "labels" or
#'   "centers".
#' @param initial_labels,initial_centers Supplied initialization used by the
#'   corresponding initialization choice.
#' @param first_index First row for deterministic max-min initialization.
#' @param init_tau Nonnegative sparse-K-median feature threshold. NULL
#'   activates the official-software quantile-grid/permutation selector.
#' @param init_tau_grid Optional explicit grid when init_tau is NULL.
#' @param init_B Number of initialization reference permutations.
#' @param init_nstart,init_max_iter Sparse-K-median start and iteration limits.
#' @param init_empty_action Explicit empty-cluster policy.
#' @param init_empty_active Explicit all-features fallback policy when a
#'   threshold excludes every feature.
#' @param init_dispersion_floor Positive floor used only to reject degenerate
#'   initialization Gap logarithms; no objective is silently altered.
#' @param outer_nstart Number of full GEM starts.
#' @param eta_mu,eta_precision Damping constants in (0,1].
#' @param mixing_floor Explicit positive mixing-weight floor.
#' @param center_weight_floor Explicit positive denominator boundary for center
#'   updates.
#' @param max_iter,convergence_tol Outer GEM controls.
#' @param bandwidth Optional transformed-radius KDE bandwidth.
#' @param bandwidth_min Explicit lower bound for an automatic or supplied
#'   bandwidth.
#' @param generator_grid_size Number of transformed-radius grid points.
#' @param generator_radius_floor,generator_density_floor Explicit positive
#'   generator reconstruction floors.
#' @param generator_score_clip Two positive radial-score clipping endpoints.
#' @param generator_spline_spar Fixed smoothing-spline parameter.
#' @param generator_extrapolation Either constant endpoint extrapolation
#'   (official software contract) or strict error.
#' @param radial_floor Positive denominator floor in sign and Tyler maps.
#' @param poet_factor_selection Supplied factor count or the official-software
#'   GR selector.
#' @param poet_factors,poet_max_factors Factor controls.
#' @param factor_ratio_floor Explicit GR ratio floor.
#' @param poet_threshold Optional POET threshold. NULL uses
#'   poet_threshold_scale * sqrt(log(p) / n_eff).
#' @param poet_threshold_scale Explicit software-rate constant.
#' @param poet_ridge Explicit nonnegative POET diagonal ridge; zero means none.
#' @param tyler_ridge Explicit ridge in the paper's ridge-stabilized Tyler map.
#' @param tyler_tol,tyler_max_iter Tyler certificate controls.
#' @param glasso_lambda Optional singleton off-diagonal graphical-lasso penalty.
#' @param glasso_lambda_grid Optional explicit EBIC grid.
#' @param glasso_lambda_scale Software-rate constant when no penalty is supplied.
#' @param glasso_ebic_gamma EBIC gamma for a multi-penalty grid.
#' @param glasso_tol,glasso_max_iter,glasso_initial_step Graphical-lasso
#'   convergence and initial-step controls.
#' @param glasso_max_backtracking Maximum number of graphical-lasso
#'   SPD backtracking steps per iteration.
#' @param glasso_majorization_tol Nonnegative tolerance for the graphical-lasso
#'   majorization/KKT certificate.
#' @param spd_tol,symmetry_tol Strict shape/precision certificate tolerances.
#' @param seed Optional local seed. It is required for stochastic starts or
#'   automatic initialization-threshold selection and is restored on exit.
#' @param strict If TRUE, fail when the selected outer fit does not converge.
#' @param keep_path Retain the deterministic outer-iteration audit path.
#'
#' @return A semc_fit object containing labels, posterior probabilities,
#'   mixing weights, centers, trace-p shape, precision, generator grid,
#'   initialization and shape certificates, explicit controls, and provenance.
#'
#' @references
#' Feng, L. and Zhuang, D. (2026). Semiparametric Elliptical Mixture
#' Clustering for High-Dimensional Data. arXiv:2605.08995.
#' \url{https://arxiv.org/abs/2605.08995}.
#'
#' @examples
#' semc_x <- matrix(c(
#'   -3, -2, -2, -3, -2, -1, -1, -2, -2.5, -2, -1.5, -2.2,
#'    3,  2,  2,  3,  2,  1,  1,  2,  2.5,  2,  1.5,  2.2
#' ), ncol = 2L, byrow = TRUE)
#' semc_labels <- rep(1:2, each = 6L)
#' semc_model <- semc_fit(
#'   semc_x, K = 2L, shape = "tyler",
#'   initialization = "labels", initial_labels = semc_labels,
#'   first_index = 1L, init_tau = 0, init_tau_grid = NULL, init_B = 2L,
#'   init_nstart = 1L, init_max_iter = 20L,
#'   init_empty_action = "error", init_empty_active = "error",
#'   init_dispersion_floor = 1e-8, outer_nstart = 1L,
#'   eta_mu = 0.7, eta_precision = 0.7,
#'   mixing_floor = 1e-10, center_weight_floor = 1e-10,
#'   max_iter = 10L, convergence_tol = 0.2,
#'   bandwidth = 0.2, bandwidth_min = 0.05,
#'   generator_grid_size = 40L, generator_radius_floor = 1e-4,
#'   generator_density_floor = 1e-10,
#'   generator_score_clip = c(1e-3, 50),
#'   generator_spline_spar = 0.55, generator_extrapolation = "constant",
#'   radial_floor = 1e-6, poet_factor_selection = "supplied",
#'   poet_factors = 0L, poet_max_factors = 0L,
#'   factor_ratio_floor = 1e-8, poet_threshold = 0.1,
#'   poet_threshold_scale = 0.55, poet_ridge = 0.05,
#'   tyler_ridge = 0.1, tyler_tol = 1e-5, tyler_max_iter = 200L,
#'   spd_tol = 0, symmetry_tol = 1e-10, seed = 11L,
#'   strict = TRUE, keep_path = FALSE
#' )
#' predict(semc_model)
#'
#' @export
semc_fit <- function(
    x, K, shape = c("glasso", "poet", "tyler"),
    initialization = c("maxmin", "software_random", "labels", "centers"),
    initial_labels = NULL, initial_centers = NULL, first_index = 1L,
    init_tau = 0, init_tau_grid = NULL, init_B = 5L,
    init_nstart = 1L, init_max_iter = 50L,
    init_empty_action = c("error", "farthest"),
    init_empty_active = c("error", "all"),
    init_dispersion_floor = 1e-8, outer_nstart = 1L,
    eta_mu = 0.7, eta_precision = 0.7,
    mixing_floor = 1e-10, center_weight_floor = 1e-10,
    max_iter = 25L, convergence_tol = 2e-3,
    bandwidth = NULL, bandwidth_min = 0.05,
    generator_grid_size = 250L, generator_radius_floor = 1e-4,
    generator_density_floor = 1e-10,
    generator_score_clip = c(1e-3, 50),
    generator_spline_spar = 0.55,
    generator_extrapolation = c("constant", "error"),
    radial_floor = 1e-6,
    poet_factor_selection = c("supplied", "software_gr"),
    poet_factors = 0L, poet_max_factors = 8L,
    factor_ratio_floor = 1e-8,
    poet_threshold = NULL, poet_threshold_scale = 0.55,
    poet_ridge = 0, tyler_ridge = 0.02,
    tyler_tol = 1e-5, tyler_max_iter = 200L,
    glasso_lambda = NULL, glasso_lambda_grid = NULL,
    glasso_lambda_scale = 0.28, glasso_ebic_gamma = 0.5,
    glasso_tol = 1e-6, glasso_max_iter = 10000L,
    glasso_initial_step = 1, glasso_max_backtracking = 100L,
    glasso_majorization_tol = 1e-12,
    spd_tol = 0, symmetry_tol = 1e-10,
    seed = NULL, strict = TRUE, keep_path = TRUE) {
  call <- match.call()
  x.name <- deparse1(substitute(x))
  x <- .semc_data(x)
  n <- nrow(x)
  p <- ncol(x)
  K <- .semc_integer(K, "K", 2L, n - 1L)
  shape <- match.arg(shape)
  initialization <- match.arg(initialization)
  init_empty_action <- match.arg(init_empty_action)
  init_empty_active <- match.arg(init_empty_active)
  generator_extrapolation <- match.arg(generator_extrapolation)
  poet_factor_selection <- match.arg(poet_factor_selection)
  strict <- .semc_flag(strict, "strict")
  keep_path <- .semc_flag(keep_path, "keep_path")
  first_index <- .semc_integer(first_index, "first_index", 1L, n)
  init_B <- .semc_integer(init_B, "init_B", 1L)
  init_nstart <- .semc_integer(init_nstart, "init_nstart", 1L)
  init_max_iter <- .semc_integer(init_max_iter, "init_max_iter", 1L)
  outer_nstart <- .semc_integer(outer_nstart, "outer_nstart", 1L)
  max_iter <- .semc_integer(max_iter, "max_iter", 1L)
  generator_grid_size <- .semc_integer(
    generator_grid_size, "generator_grid_size", 20L
  )
  poet_factors <- .semc_integer(poet_factors, "poet_factors", 0L, p - 1L)
  poet_max_factors.requested <- .semc_integer(
    poet_max_factors, "poet_max_factors", 0L
  )
  poet_max_factors <- min(poet_max_factors.requested, max(0L, p - 2L))
  tyler_max_iter <- .semc_integer(
    tyler_max_iter, "tyler_max_iter", 1L
  )
  glasso_max_iter <- .semc_integer(
    glasso_max_iter, "glasso_max_iter", 1L
  )
  glasso_max_backtracking <- .semc_integer(
    glasso_max_backtracking, "glasso_max_backtracking", 1L
  )

  eta_mu <- .semc_scalar(eta_mu, "eta_mu", 0, 1, TRUE)
  eta_precision <- .semc_scalar(
    eta_precision, "eta_precision", 0, 1, TRUE
  )
  mixing_floor <- .semc_scalar(
    mixing_floor, "mixing_floor", 0, 1 / K, TRUE, TRUE
  )
  center_weight_floor <- .semc_scalar(
    center_weight_floor, "center_weight_floor", 0, Inf, TRUE
  )
  convergence_tol <- .semc_scalar(
    convergence_tol, "convergence_tol", 0, Inf, TRUE
  )
  bandwidth_min <- .semc_scalar(
    bandwidth_min, "bandwidth_min", 0, Inf, TRUE
  )
  generator_radius_floor <- .semc_scalar(
    generator_radius_floor, "generator_radius_floor", 0, Inf, TRUE
  )
  generator_density_floor <- .semc_scalar(
    generator_density_floor, "generator_density_floor", 0, Inf, TRUE
  )
  generator_spline_spar <- .semc_scalar(
    generator_spline_spar, "generator_spline_spar", 0, 1
  )
  radial_floor <- .semc_scalar(
    radial_floor, "radial_floor", 0, Inf, TRUE
  )
  factor_ratio_floor <- .semc_scalar(
    factor_ratio_floor, "factor_ratio_floor", 0, Inf, TRUE
  )
  poet_threshold_scale <- .semc_scalar(
    poet_threshold_scale, "poet_threshold_scale", 0
  )
  poet_ridge <- .semc_scalar(poet_ridge, "poet_ridge", 0)
  tyler_ridge <- .semc_scalar(
    tyler_ridge, "tyler_ridge", 0, 1, FALSE, TRUE
  )
  tyler_tol <- .semc_scalar(tyler_tol, "tyler_tol", 0, Inf, TRUE)
  glasso_lambda_scale <- .semc_scalar(
    glasso_lambda_scale, "glasso_lambda_scale", 0, Inf, TRUE
  )
  glasso_ebic_gamma <- .semc_scalar(
    glasso_ebic_gamma, "glasso_ebic_gamma", 0
  )
  glasso_tol <- .semc_scalar(
    glasso_tol, "glasso_tol", 0, Inf, TRUE
  )
  glasso_initial_step <- .semc_scalar(
    glasso_initial_step, "glasso_initial_step", 0, Inf, TRUE
  )
  glasso_majorization_tol <- .semc_scalar(
    glasso_majorization_tol, "glasso_majorization_tol", 0
  )
  spd_tol <- .semc_scalar(spd_tol, "spd_tol", 0)
  symmetry_tol <- .semc_scalar(
    symmetry_tol, "symmetry_tol", 0, Inf, TRUE
  )
  init_dispersion_floor <- .semc_scalar(
    init_dispersion_floor, "init_dispersion_floor", 0, Inf, TRUE
  )
  if (!is.null(bandwidth)) {
    bandwidth <- .semc_scalar(bandwidth, "bandwidth", 0, Inf, TRUE)
  }
  if (!is.null(init_tau)) {
    init_tau <- .semc_scalar(init_tau, "init_tau", 0)
  }
  if (!is.null(poet_threshold)) {
    poet_threshold <- .semc_scalar(poet_threshold, "poet_threshold", 0)
  }
  if (!is.null(glasso_lambda)) {
    glasso_lambda <- .semc_scalar(
      glasso_lambda, "glasso_lambda", 0, Inf, TRUE
    )
  }
  if (!is.null(glasso_lambda) && !is.null(glasso_lambda_grid)) {
    stop("Supply at most one of glasso_lambda and glasso_lambda_grid.",
         call. = FALSE)
  }
  if (!is.null(glasso_lambda_grid)) {
    glasso_lambda_grid <- as.numeric(glasso_lambda_grid)
    if (!length(glasso_lambda_grid) || anyNA(glasso_lambda_grid) ||
        any(!is.finite(glasso_lambda_grid)) ||
        any(glasso_lambda_grid <= 0)) {
      stop("glasso_lambda_grid must contain finite positive values.",
           call. = FALSE)
    }
  }
  generator_score_clip <- as.numeric(generator_score_clip)
  if (length(generator_score_clip) != 2L ||
      anyNA(generator_score_clip) ||
      any(!is.finite(generator_score_clip)) ||
      generator_score_clip[1L] <= 0 ||
      generator_score_clip[2L] <= generator_score_clip[1L]) {
    stop("generator_score_clip must contain two increasing positive values.",
         call. = FALSE)
  }
  if (initialization %in% c("labels", "centers") &&
      is.null(init_tau)) {
    stop("Supplied labels or centers require an explicit init_tau.",
         call. = FALSE)
  }
  if (initialization %in% c("labels", "centers") &&
      (init_nstart != 1L || outer_nstart != 1L)) {
    stop(paste(
      "Supplied labels or centers require init_nstart = 1 and",
      "outer_nstart = 1."
    ), call. = FALSE)
  }
  stochastic <- initialization == "software_random" || is.null(init_tau)
  if (stochastic && is.null(seed)) {
    stop(paste(
      "An explicit seed is required for software_random initialization",
      "or automatic initialization-threshold selection."
    ), call. = FALSE)
  }
  if (!is.null(seed)) {
    seed <- .semc_integer(seed, "seed", 0L)
  }

  result <- .semc_with_seed(seed, function() {
    tau.selection <- NULL
    selected.tau <- init_tau
    if (is.null(selected.tau)) {
      tau.selection <- .semc_select_initial_tau(
        x, K, init_tau_grid, init_B, init_nstart, init_max_iter,
        initialization, first_index, init_empty_action,
        init_empty_active, init_dispersion_floor
      )
      selected.tau <- tau.selection$tau
    }
    starts <- vector("list", outer_nstart)
    summaries <- vector("list", outer_nstart)
    best <- NULL
    for (outer in seq_len(outer_nstart)) {
      init.fit <- tryCatch(
        .semc_sparse_kmedian_multistart(
          x, K, selected.tau, init_nstart, init_max_iter,
          initialization,
          ((first_index + outer - 2L) %% n) + 1L,
          initial_labels, initial_centers,
          init_empty_action, init_empty_active
        ),
        error = identity
      )
      if (inherits(init.fit, "condition")) {
        summaries[[outer]] <- data.frame(
          start = outer, objective = NA_real_, iterations = NA_integer_,
          converged = FALSE, valid = FALSE,
          failure = conditionMessage(init.fit)
        )
        next
      }
      fit <- tryCatch(
        .semc_fit_single(
          x, K, init.fit, shape, eta_mu, eta_precision,
          mixing_floor, center_weight_floor, max_iter,
          convergence_tol, bandwidth, bandwidth_min,
          generator_grid_size, generator_radius_floor,
          generator_density_floor, generator_score_clip,
          generator_spline_spar, generator_extrapolation,
          radial_floor, poet_factor_selection, poet_factors,
          poet_max_factors, factor_ratio_floor, poet_threshold,
          poet_threshold_scale, poet_ridge, tyler_ridge,
          tyler_tol, tyler_max_iter, glasso_lambda,
          glasso_lambda_grid, glasso_lambda_scale,
          glasso_ebic_gamma, glasso_tol, glasso_max_iter,
          glasso_initial_step, glasso_max_backtracking,
          glasso_majorization_tol, spd_tol, symmetry_tol,
          keep_path
        ),
        error = identity
      )
      if (inherits(fit, "condition")) {
        summaries[[outer]] <- data.frame(
          start = outer, objective = NA_real_, iterations = NA_integer_,
          converged = FALSE, valid = FALSE,
          failure = conditionMessage(fit)
        )
        next
      }
      starts[[outer]] <- fit
      summaries[[outer]] <- data.frame(
        start = outer, objective = fit$objective,
        iterations = fit$iterations, converged = fit$converged,
        valid = fit$valid, failure = NA_character_
      )
      if ((isTRUE(fit$valid) || !strict) &&
          (is.null(best) || fit$objective > best$objective)) {
        best <- fit
        best$best.start <- outer
      }
    }
    if (is.null(best)) {
      summary.table <- do.call(rbind, summaries)
      stop(paste(
        "No outer SEMC start passed all convergence and shape certificates.",
        paste(stats::na.omit(summary.table$failure), collapse = "; ")
      ), call. = FALSE)
    }
    best$outer.starts <- starts
    best$start.summaries <- do.call(rbind, summaries)
    best$initialization$tau <- selected.tau
    best$initialization$tau.selection <- tau.selection
    best
  })

  if (!isTRUE(result$converged) && strict) {
    stop("The selected SEMC fit did not satisfy convergence_tol.",
         call. = FALSE)
  }
  variable.names <- colnames(x)
  if (!is.null(variable.names)) {
    colnames(result$centers) <- variable.names
    dimnames(result$shape) <- dimnames(result$precision) <-
      list(variable.names, variable.names)
  }
  rownames(result$centers) <- paste0("cluster", seq_len(K))
  colnames(result$posterior) <- paste0("cluster", seq_len(K))
  names(result$mixing) <- paste0("cluster", seq_len(K))
  result$call <- call
  result$data.name <- x.name
  result$method <- "Semiparametric elliptical mixture clustering"
  result$K <- K
  result$shape.method <- shape
  result$valid <- isTRUE(result$converged) &&
    isTRUE(result$shape.fit$certificate$valid)
  result$controls <- list(
    initialization = list(
      method = initialization, first.index = first_index,
      tau = init_tau, tau.grid = init_tau_grid, B = init_B,
      nstart = init_nstart, max.iter = init_max_iter,
      empty.action = init_empty_action,
      empty.active = init_empty_active,
      dispersion.floor = init_dispersion_floor,
      software.contract = initialization == "software_random" ||
        is.null(init_tau)
    ),
    outer = list(
      nstart = outer_nstart, eta.mu = eta_mu,
      eta.precision = eta_precision, max.iter = max_iter,
      convergence.tol = convergence_tol,
      mixing.floor = mixing_floor,
      center.weight.floor = center_weight_floor
    ),
    generator = list(
      bandwidth = bandwidth, bandwidth.min = bandwidth_min,
      grid.size = generator_grid_size,
      radius.floor = generator_radius_floor,
      density.floor = generator_density_floor,
      score.clip = generator_score_clip,
      spline.spar = generator_spline_spar,
      extrapolation = generator_extrapolation
    ),
    shape = list(
      method = shape, radial.floor = radial_floor,
      poet.factor.selection = poet_factor_selection,
      poet.factors = poet_factors,
      poet.max.factors = poet_max_factors,
      poet.max.factors.requested = poet_max_factors.requested,
      factor.ratio.floor = factor_ratio_floor,
      poet.threshold = poet_threshold,
      poet.threshold.scale = poet_threshold_scale,
      poet.ridge = poet_ridge, tyler.ridge = tyler_ridge,
      tyler.tol = tyler_tol, tyler.max.iter = tyler_max_iter,
      glasso.lambda = glasso_lambda,
      glasso.lambda.grid = glasso_lambda_grid,
      glasso.lambda.scale = glasso_lambda_scale,
      glasso.ebic.gamma = glasso_ebic_gamma,
      glasso.tol = glasso_tol,
      glasso.max.iter = glasso_max_iter,
      glasso.initial.step = glasso_initial_step,
      glasso.max.backtracking = glasso_max_backtracking,
      glasso.majorization.tol = glasso_majorization_tol,
      spd.tol = spd_tol, symmetry.tol = symmetry_tol
    ),
    seed = seed, strict = strict, keep.path = keep_path
  )
  result$provenance <- list(
    paper = "Feng and Zhuang (2026), arXiv:2605.08995",
    official.software = paste0(
      "flnankai/GEMcluster@",
      "10fce04fe690fe274dd5d237cfcd3d5c6a4139f6"
    ),
    license = "upstream MIT; HDElliptical implementation independently rewritten",
    software.contracts = c(
      "sparse-K-median initialization and optional permutation tau selection",
      "weighted Silverman bandwidth with explicit lower bound",
      "fixed smoothing spline and explicit radial-score clipping",
      "optional GR factor selector and EBIC glasso grid"
    ),
    deviations = c(
      "paper generator normalization is reimposed after smoothing",
      "off-diagonal graphical lasso is solved without GPL huge",
      "no project-PD or failure ridge fallback is applied"
    )
  )
  class(result) <- c("semc_fit", "hd_clustering_fit")
  result
}



#' Predict from a semiparametric elliptical mixture fit
#'
#' @param object A valid object returned by semc_fit().
#' @param newdata Optional numeric matrix. NULL returns fitted predictions.
#' @param type Either "class" or "posterior".
#' @param ... Unused.
#'
#' @return Integer class labels or a posterior-probability matrix.
#' @export
predict.semc_fit <- function(
    object, newdata = NULL, type = c("class", "posterior"), ...) {
  dots <- list(...)
  if (length(dots)) {
    stop("No additional prediction arguments are supported.", call. = FALSE)
  }
  if (!inherits(object, "semc_fit") || !isTRUE(object$valid)) {
    stop("object must be a valid certified semc_fit.", call. = FALSE)
  }
  type <- match.arg(type)
  if (is.null(newdata)) {
    posterior <- object$posterior
  } else {
    newdata <- .semc_data(newdata, "newdata", min_rows = 1L)
    if (ncol(newdata) != ncol(object$centers)) {
      stop("newdata does not match the fitted feature dimension.",
           call. = FALSE)
    }
    fitted.names <- colnames(object$centers)
    if (!is.null(fitted.names) &&
        !identical(colnames(newdata), fitted.names)) {
      stop("newdata column names and order must match the fitted data.",
           call. = FALSE)
    }
    delta <- .semc_delta(newdata, object$centers, object$precision)
    posterior <- .semc_posterior(
      delta, object$mixing, object$generator
    )$posterior
    colnames(posterior) <- names(object$mixing)
  }
  if (type == "posterior") {
    return(posterior)
  }
  as.integer(max.col(posterior, ties.method = "first"))
}


.semc_gap_dispersion <- function(x, fit, contract) {
  delta <- .semc_delta(x, fit$centers, fit$precision)
  if (contract == "paper_hard_log1p") {
    return(mean(log1p(delta[cbind(seq_len(nrow(x)), fit$cluster)])))
  }
  sum(fit$posterior * delta)
}


.semc_permutation_plan <- function(x, B, permutation_indices) {
  n <- nrow(x)
  p <- ncol(x)
  if (!is.null(permutation_indices)) {
    if (!is.array(permutation_indices) ||
        !identical(dim(permutation_indices), c(n, p, B)) ||
        !is.numeric(permutation_indices) || anyNA(permutation_indices) ||
        any(!is.finite(permutation_indices)) ||
        any(permutation_indices != floor(permutation_indices))) {
      stop("permutation_indices must be an integer n by p by B array.",
           call. = FALSE)
    }
    storage.mode(permutation_indices) <- "integer"
    for (b in seq_len(B)) {
      for (j in seq_len(p)) {
        if (!identical(sort(permutation_indices[, j, b]), seq_len(n))) {
          stop("Every permutation_indices column must permute 1:n.",
               call. = FALSE)
        }
      }
    }
    return(permutation_indices)
  }
  plan <- array(NA_integer_, dim = c(n, p, B))
  for (b in seq_len(B)) {
    for (j in seq_len(p)) {
      plan[, j, b] <- sample.int(n)
    }
  }
  plan
}


#' Select the SEMC cluster count by Gap-LSE
#'
#' The default follows equations (2.18)--(2.21) of Feng and Zhuang: hard fitted
#' labels, radial dispersion mean(log(1 + delta)), independently
#' column-permuted reference samples, and the lower-complexity one-standard-
#' error rule. The pinned official software instead defaults to posterior-
#' weighted sum(tau * delta); that distinct contract is available only through
#' the explicit dispersion = "software_soft_delta" choice.
#'
#' @param x Numeric matrix with observations in rows.
#' @param k_grid Sorted candidate component counts.
#' @param B At least two reference permutations.
#' @param control Named arguments passed to semc_fit(); x, K, and seed are
#'   controlled by this selector.
#' @param dispersion Paper hard-log1p or official-software soft-delta contract.
#' @param rule Lower-complexity one-standard-error or maximum-Gap selection.
#' @param permutation_indices Optional integer n by p by B permutation array.
#' @param seed Explicit local seed for permutation and fit seeds.
#' @param dispersion_floor Positive rejection boundary for logarithms.
#' @param keep_reference_fits Retain every reference fit.
#' @param keep_permutations Retain the permutation array.
#'
#' @return A semc_gap_selection object with both LSE and maximum-Gap choices,
#'   observed fits, all formula ingredients, local seeds, and provenance.
#'
#' @references
#' Feng, L. and Zhuang, D. (2026). Semiparametric Elliptical Mixture
#' Clustering for High-Dimensional Data. arXiv:2605.08995.
#' \url{https://arxiv.org/abs/2605.08995}.
#'
#' @examples
#' \donttest{
#' semc_x <- matrix(c(
#'   -3, -2, -2, -3, -2, -1, -1, -2, -2.5, -2, -1.5, -2.2,
#'    3,  2,  2,  3,  2,  1,  1,  2,  2.5,  2,  1.5,  2.2
#' ), ncol = 2L, byrow = TRUE)
#' semc_labels <- rep(1:2, each = 6L)
#' semc_B <- 2L
#' semc_plan <- array(
#'   rep(seq_len(nrow(semc_x)), ncol(semc_x) * semc_B),
#'   dim = c(nrow(semc_x), ncol(semc_x), semc_B)
#' )
#' semc_control <- list(
#'   shape = "tyler", initialization = "labels",
#'   initial_labels = semc_labels, first_index = 1L,
#'   init_tau = 0, init_tau_grid = NULL, init_B = 2L,
#'   init_nstart = 1L, init_max_iter = 20L,
#'   init_empty_action = "error", init_empty_active = "error",
#'   init_dispersion_floor = 1e-8, outer_nstart = 1L,
#'   eta_mu = 0.7, eta_precision = 0.7,
#'   mixing_floor = 1e-10, center_weight_floor = 1e-10,
#'   max_iter = 10L, convergence_tol = 0.2,
#'   bandwidth = 0.2, bandwidth_min = 0.05,
#'   generator_grid_size = 40L, generator_radius_floor = 1e-4,
#'   generator_density_floor = 1e-10,
#'   generator_score_clip = c(1e-3, 50),
#'   generator_spline_spar = 0.55, generator_extrapolation = "constant",
#'   radial_floor = 1e-6, poet_factor_selection = "supplied",
#'   poet_factors = 0L, poet_max_factors = 0L,
#'   factor_ratio_floor = 1e-8, poet_threshold = 0.1,
#'   poet_threshold_scale = 0.55, poet_ridge = 0.05,
#'   tyler_ridge = 0.1, tyler_tol = 1e-5, tyler_max_iter = 200L,
#'   spd_tol = 0, symmetry_tol = 1e-10,
#'   strict = TRUE, keep_path = FALSE
#' )
#' semc_gap <- semc_select_k_gap(
#'   semc_x, k_grid = 2L, B = semc_B, control = semc_control,
#'   dispersion = "paper_hard_log1p", rule = "lse",
#'   permutation_indices = semc_plan, seed = 19L,
#'   dispersion_floor = 1e-12, keep_reference_fits = FALSE,
#'   keep_permutations = FALSE
#' )
#' semc_gap$selected.k
#' }
#'
#' @export
semc_select_k_gap <- function(
    x, k_grid = 2:5, B = 20L, control = list(),
    dispersion = c("paper_hard_log1p", "software_soft_delta"),
    rule = c("lse", "max"), permutation_indices = NULL,
    seed, dispersion_floor = 1e-12,
    keep_reference_fits = FALSE, keep_permutations = FALSE) {
  call <- match.call()
  x <- .semc_data(x)
  n <- nrow(x)
  k_grid <- sort(unique(as.numeric(k_grid)))
  if (!length(k_grid) || anyNA(k_grid) || any(!is.finite(k_grid)) ||
      any(k_grid != floor(k_grid)) || any(k_grid < 2L) ||
      any(k_grid >= n)) {
    stop("k_grid must contain integers from 2 through n - 1.",
         call. = FALSE)
  }
  k_grid <- as.integer(k_grid)
  B <- .semc_integer(B, "B", 2L)
  seed <- .semc_integer(seed, "seed", 0L)
  dispersion_floor <- .semc_scalar(
    dispersion_floor, "dispersion_floor", 0, Inf, TRUE
  )
  keep_reference_fits <- .semc_flag(
    keep_reference_fits, "keep_reference_fits"
  )
  keep_permutations <- .semc_flag(
    keep_permutations, "keep_permutations"
  )
  dispersion <- match.arg(dispersion)
  rule <- match.arg(rule)
  if (!is.list(control) ||
      (length(control) > 0L &&
       (is.null(names(control)) || any(!nzchar(names(control)))))) {
    stop("control must be a named list.", call. = FALSE)
  }
  forbidden <- intersect(names(control), c("x", "K", "seed"))
  if (length(forbidden)) {
    stop("control must not contain x, K, or seed.", call. = FALSE)
  }

  result <- .semc_with_seed(seed, function() {
    plan <- .semc_permutation_plan(x, B, permutation_indices)
    fit.seeds <- matrix(
      sample.int(.Machine$integer.max, length(k_grid) * (B + 1L)),
      nrow = length(k_grid), ncol = B + 1L
    )
    observed.fits <- vector("list", length(k_grid))
    observed.dispersion <- numeric(length(k_grid))
    reference.log.dispersion <- matrix(
      NA_real_, nrow = length(k_grid), ncol = B
    )
    reference.fits <- if (keep_reference_fits) {
      matrix(vector("list", length(k_grid) * B),
             nrow = length(k_grid), ncol = B)
    } else {
      NULL
    }

    for (index in seq_along(k_grid)) {
      fit <- do.call(semc_fit, c(
        list(x = x, K = k_grid[index], seed = fit.seeds[index, 1L]),
        control
      ))
      if (!isTRUE(fit$valid)) {
        stop(sprintf("Observed SEMC fit for K=%d is not valid.",
                     k_grid[index]), call. = FALSE)
      }
      observed.fits[[index]] <- fit
      observed.dispersion[index] <- .semc_gap_dispersion(
        x, fit, dispersion
      )
      if (!is.finite(observed.dispersion[index]) ||
          observed.dispersion[index] <= dispersion_floor) {
        stop(paste0(
          "Observed dispersion is not above the explicit ",
          "dispersion_floor for K=", k_grid[index], "."
        ), call. = FALSE)
      }
    }
    for (b in seq_len(B)) {
      reference.x <- x
      for (j in seq_len(ncol(x))) {
        reference.x[, j] <- x[plan[, j, b], j]
      }
      for (index in seq_along(k_grid)) {
        fit <- do.call(semc_fit, c(
          list(
            x = reference.x, K = k_grid[index],
            seed = fit.seeds[index, b + 1L]
          ),
          control
        ))
        if (!isTRUE(fit$valid)) {
          stop(sprintf(
            "Reference SEMC fit for K=%d, permutation=%d is not valid.",
            k_grid[index], b
          ), call. = FALSE)
        }
        value <- .semc_gap_dispersion(reference.x, fit, dispersion)
        if (!is.finite(value) || value <= dispersion_floor) {
          stop(paste0(
            "Reference dispersion is not above the explicit ",
            "dispersion_floor for K=", k_grid[index],
            ", permutation=", b, "."
          ), call. = FALSE)
        }
        reference.log.dispersion[index, b] <- log(value)
        if (keep_reference_fits) {
          reference.fits[[index, b]] <- fit
        }
      }
    }

    log.observed <- log(observed.dispersion)
    gap <- rowMeans(reference.log.dispersion) - log.observed
    reference.sd <- apply(reference.log.dispersion, 1L, stats::sd)
    standard.error <- sqrt(1 + 1 / B) * reference.sd
    maximum.index <- which.max(gap)
    lse.index <- length(k_grid)
    if (length(k_grid) > 1L) {
      for (index in seq_len(length(k_grid) - 1L)) {
        if (gap[index] >= gap[index + 1L] -
            standard.error[index + 1L]) {
          lse.index <- index
          break
        }
      }
    }
    selected.index <- if (rule == "lse") lse.index else maximum.index
    summary <- data.frame(
      K = k_grid, observed.dispersion = observed.dispersion,
      log.observed.dispersion = log.observed,
      mean.log.reference = rowMeans(reference.log.dispersion),
      sd.log.reference = reference.sd,
      standard.error = standard.error, gap = gap,
      lse.margin = c(
        if (length(k_grid) > 1L) {
          gap[-length(gap)] -
            (gap[-1L] - standard.error[-1L])
        } else {
          numeric()
        },
        NA_real_
      )
    )
    list(
      selected.k = k_grid[selected.index],
      selected.index = selected.index,
      selected.fit = observed.fits[[selected.index]],
      lse.k = k_grid[lse.index],
      maximum.gap.k = k_grid[maximum.index],
      summary = summary,
      observed.fits = stats::setNames(
        observed.fits, paste0("K", k_grid)
      ),
      reference.log.dispersion = reference.log.dispersion,
      reference.fits = reference.fits,
      permutation.indices = if (keep_permutations) plan else NULL,
      fit.seeds = fit.seeds
    )
  })

  result$call <- call
  result$B <- B
  result$rule <- rule
  result$dispersion <- dispersion
  result$control <- control
  result$diagnostics <- list(
    formula = if (dispersion == "paper_hard_log1p") {
      "W_K = mean_i log(1 + delta_i,C_i)"
    } else {
      "official software: sum_i,k posterior_i,k * delta_i,k"
    },
    gap = "mean_b log(W_K^0,b) - log(W_K)",
    standard.error = "sqrt(1 + 1/B) * sd_b(log(W_K^0,b))",
    selection = paste(
      "smallest K_j with Gap(K_j) >= Gap(K_j+1) - s(K_j+1);",
      "largest candidate if none"
    ),
    rng.isolation = TRUE,
    dispersion.floor = dispersion_floor,
    paper.software.conflict.exposed = TRUE,
    official.software.commit =
      "10fce04fe690fe274dd5d237cfcd3d5c6a4139f6"
  )
  class(result) <- "semc_gap_selection"
  result
}









