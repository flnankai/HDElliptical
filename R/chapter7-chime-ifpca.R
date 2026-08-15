# Chapter 7 Gaussian-mixture and influential-feature clustering.

.c7_scalar <- function(value, name, lower = -Inf, strict_lower = FALSE) {
  value <- as.numeric(value)
  bad <- length(value) != 1L || is.na(value) || !is.finite(value) ||
    if (strict_lower) value <= lower else value < lower
  if (bad) {
    relation <- if (strict_lower) "larger than" else "not smaller than"
    stop(sprintf("'%s' must be one finite number %s %g.",
                 name, relation, lower), call. = FALSE)
  }
  value
}


.c7_count <- function(value, name, minimum = 1L) {
  value <- as.numeric(value)
  if (length(value) != 1L || is.na(value) || !is.finite(value) ||
      value != floor(value) || value < minimum ||
      value > .Machine$integer.max) {
    stop(sprintf("'%s' must be one integer not smaller than %d.",
                 name, minimum), call. = FALSE)
  }
  as.integer(value)
}


.c7_flag <- function(value, name) {
  if (!is.logical(value) || length(value) != 1L || is.na(value)) {
    stop(sprintf("'%s' must be TRUE or FALSE.", name), call. = FALSE)
  }
  value
}


.c7_data <- function(x, min_rows = 3L, min_columns = 1L) {
  x <- .as_data_matrix(x, "x", min_rows = min_rows)
  if (ncol(x) < min_columns) {
    stop(sprintf("'x' must have at least %d columns.", min_columns),
         call. = FALSE)
  }
  feature.names <- colnames(x)
  if (!is.null(feature.names) &&
      (anyNA(feature.names) || any(!nzchar(feature.names)) ||
       anyDuplicated(feature.names))) {
    stop("Column names of 'x' must be non-empty and unique.", call. = FALSE)
  }
  x
}


.c7_align_vector <- function(value, x, name) {
  original.names <- names(value)
  value <- as.numeric(value)
  p <- ncol(x)
  if (length(value) != p || anyNA(value) || any(!is.finite(value))) {
    stop(sprintf("'initial$%s' must be a finite vector of length %d.",
                 name, p), call. = FALSE)
  }
  feature.names <- colnames(x)
  if (!is.null(original.names) && !is.null(feature.names)) {
    if (!setequal(original.names, feature.names)) {
      stop(sprintf(
        "Names on 'initial$%s' do not match 'colnames(x)'.", name
      ), call. = FALSE)
    }
    value <- value[match(feature.names, original.names)]
  }
  names(value) <- feature.names
  value
}


.c7_align_covariance <- function(value, x, psd_tol) {
  p <- ncol(x)
  if (!is.matrix(value) || !is.numeric(value) ||
      !identical(dim(value), c(p, p))) {
    stop(sprintf(
      "'initial$covariance' must be a numeric %d by %d matrix.", p, p
    ), call. = FALSE)
  }
  storage.mode(value) <- "double"
  if (anyNA(value) || any(!is.finite(value))) {
    stop("'initial$covariance' must contain only finite values.",
         call. = FALSE)
  }
  feature.names <- colnames(x)
  rn <- rownames(value)
  cn <- colnames(value)
  if (xor(is.null(rn), is.null(cn))) {
    stop("'initial$covariance' must have both dimname vectors or neither.",
         call. = FALSE)
  }
  if (!is.null(rn) && !is.null(feature.names)) {
    if (!setequal(rn, feature.names) || !setequal(cn, feature.names)) {
      stop(paste(
        "Dimnames on 'initial$covariance' do not match 'colnames(x)'."
      ), call. = FALSE)
    }
    value <- value[match(feature.names, rn), match(feature.names, cn),
                   drop = FALSE]
  }
  scale <- max(1, max(abs(value)))
  symmetry.tolerance <- 512 * .Machine$double.eps * p * scale
  asymmetry <- max(abs(value - t(value)))
  if (asymmetry > symmetry.tolerance) {
    stop("'initial$covariance' must be numerically symmetric.",
         call. = FALSE)
  }
  value <- value / 2 + t(value) / 2
  eigenvalues <- eigen(value, symmetric = TRUE, only.values = TRUE)$values
  spectral.scale <- max(1, max(abs(eigenvalues)))
  if (min(eigenvalues) < -psd_tol * spectral.scale) {
    stop(paste(
      "'initial$covariance' is materially indefinite at 'psd_tol';",
      "no PSD projection or ridge was applied."
    ), call. = FALSE)
  }
  if (any(diag(value) <= 0)) {
    stop(paste(
      "Every diagonal of 'initial$covariance' must be positive for the",
      "quadratic-lasso solver; no diagonal floor was applied."
    ), call. = FALSE)
  }
  dimnames(value) <- list(feature.names, feature.names)
  list(
    matrix = value,
    eigenvalues = eigenvalues,
    minimum.eigenvalue = min(eigenvalues),
    rank = sum(eigenvalues > psd_tol * spectral.scale),
    psd.tolerance = psd_tol * spectral.scale,
    asymmetry = asymmetry,
    symmetry.tolerance = symmetry.tolerance
  )
}


.c7_chime_initial <- function(initial, x, psd_tol) {
  if (!is.list(initial)) {
    stop(paste(
      "'initial' must explicitly supply 'omega', 'mu1', 'mu2', and",
      "'covariance'; CHIME does not infer an initializer."
    ), call. = FALSE)
  }
  required <- c("omega", "mu1", "mu2", "covariance")
  missing <- setdiff(required, names(initial))
  if (length(missing)) {
    stop(sprintf("'initial' is missing: %s.", paste(missing, collapse = ", ")),
         call. = FALSE)
  }
  omega <- as.numeric(initial$omega)
  if (length(omega) != 1L || is.na(omega) || !is.finite(omega) ||
      omega <= 0 || omega >= 1) {
    stop(paste(
      "'initial$omega' must be one finite number strictly between 0 and 1."
    ), call. = FALSE)
  }
  list(
    omega = omega,
    mu1 = .c7_align_vector(initial$mu1, x, "mu1"),
    mu2 = .c7_align_vector(initial$mu2, x, "mu2"),
    covariance = .c7_align_covariance(initial$covariance, x, psd_tol)
  )
}


#' CHIME clustering for a two-component Gaussian mixture
#'
#' Fits the CHIME iteration of Cai, Ma, and Zhang for a two-component
#' Gaussian mixture with a common covariance. At stage \eqn{t}, the sparse
#' discriminant vector is the certified solution of
#' \deqn{\min_\beta\;\tfrac12\beta^T\widehat\Sigma_t\beta-
#' \beta^T(\widehat\mu_{1t}-\widehat\mu_{2t})+
#' \lambda_t\|\beta\|_1.}
#' The covariance is never inverted. Coordinate descent uses the actual
#' diagonal of the quadratic metric and must satisfy the exact lasso KKT
#' equations before an iterate is accepted.
#'
#' This function deliberately differs from the authors' numerical script in
#' three ways required for a reusable unsupervised method: it applies no hidden
#' covariance ridge, it never uses true labels to choose a penalty, and it
#' requires the complete penalty path (including the initial beta fit) from the
#' caller. The method is only defined here for K = 2, matching the stated
#' model and the available primary implementation.
#'
#' @param x Finite numeric matrix with observations in rows.
#' @param lambda Explicit non-negative vector of length max_iter + 1. Its first
#'   entry fits the initial discriminant vector and each remaining entry is
#'   used by one EM update. No path is generated or reordered.
#' @param initial Explicit list with omega, mu1, mu2, and covariance. Whenever
#'   omega exceeds one half, component labels are exchanged so the identifiable
#'   representation has omega at most one half.
#' @param K Number of components. Only 2 is implemented.
#' @param max_iter Number of EM updates. All are run; there is no unreported
#'   early stopping or reuse of a shorter lambda path.
#' @param lasso_tol Positive relative-update and scaled-KKT tolerance.
#' @param lasso_max_iter Maximum coordinate-descent sweeps at every stage.
#' @param support_tol Non-negative reporting threshold for coefficient support.
#'   It does not alter coefficients or KKT equations.
#' @param psd_tol Non-negative relative tolerance used only to diagnose the
#'   initial covariance as positive semidefinite. Eigenvalues are not changed.
#'
#' @return A chime_fit object containing labels, responsibilities, all
#'   parameter stages, the explicit lambda path, covariance matrices, and
#'   per-stage objective, convergence, and KKT certificates. Score ties use
#'   the printed weak inequality and are assigned to component 1.
#'
#' @references
#' Cai, T. T., Ma, J., and Zhang, L. (2019). CHIME: Clustering of
#' high-dimensional Gaussian mixtures with EM algorithm and its optimality.
#' \emph{Annals of Statistics}, \strong{47}, 1234--1267.
#' \doi{10.1214/18-AOS1711}.
#'
#' @examples
#' x <- rbind(
#'   c(-2.2, -1.0), c(-1.8, -1.2), c(-2.0, -0.8), c(-1.7, -1.1),
#'   c( 1.8,  1.0), c( 2.2,  1.1), c( 2.0,  0.9), c( 1.7,  1.2)
#' )
#' initial <- list(
#'   omega = 0.5, mu1 = c(-1.5, -0.8), mu2 = c(1.5, 0.8),
#'   covariance = diag(2)
#' )
#' chime_clustering(x, lambda = c(0.2, 0.15, 0.1), initial = initial,
#'                  max_iter = 2)
#'
#' @export
chime_clustering <- function(
    x, lambda, initial, K = 2L, max_iter = length(lambda) - 1L,
    lasso_tol = 1e-8, lasso_max_iter = 10000L, support_tol = 0,
    psd_tol = sqrt(.Machine$double.eps)) {
  call <- match.call()
  data.name <- deparse(substitute(x))
  x <- .c7_data(x, min_rows = 3L)
  K <- .c7_count(K, "K", minimum = 2L)
  if (K != 2L) {
    stop("'chime_clustering()' implements only two components.",
         call. = FALSE)
  }
  max_iter <- .c7_count(max_iter, "max_iter", minimum = 1L)
  lambda <- as.numeric(lambda)
  if (length(lambda) != max_iter + 1L || anyNA(lambda) ||
      any(!is.finite(lambda)) || any(lambda < 0)) {
    stop(paste(
      "'lambda' must be an explicit finite non-negative vector of length",
      "'max_iter + 1'; no implicit or recycled path is used."
    ), call. = FALSE)
  }
  lasso_tol <- .c7_scalar(lasso_tol, "lasso_tol", 0, strict_lower = TRUE)
  lasso_max_iter <- .c7_count(lasso_max_iter, "lasso_max_iter")
  support_tol <- .c7_scalar(support_tol, "support_tol", 0)
  psd_tol <- .c7_scalar(psd_tol, "psd_tol", 0)
  initial <- .c7_chime_initial(initial, x, psd_tol)

  core <- cpp_ch7_chime_fit(
    x, initial$omega, initial$mu1, initial$mu2,
    initial$covariance$matrix, lambda, lasso_tol, lasso_max_iter
  )
  stages <- 0:max_iter
  feature.names <- colnames(x)
  row.names <- rownames(x)
  colnames(core$mu1) <- feature.names
  colnames(core$mu2) <- feature.names
  colnames(core$beta) <- feature.names
  rownames(core$mu1) <- stages
  rownames(core$mu2) <- stages
  rownames(core$beta) <- stages
  for (index in seq_along(core$covariance_history)) {
    dimnames(core$covariance_history[[index]]) <-
      list(feature.names, feature.names)
  }
  dimnames(core$covariance) <- list(feature.names, feature.names)
  names(core$responsibility) <- row.names
  names(core$score) <- row.names
  names(core$cluster) <- row.names
  if (ncol(core$responsibility_history)) {
    dimnames(core$responsibility_history) <-
      list(row.names, seq_len(max_iter))
  }

  final <- max_iter + 1L
  estimate <- list(
    omega = core$omega[[final]],
    mu1 = setNames(as.numeric(core$mu1[final, ]), feature.names),
    mu2 = setNames(as.numeric(core$mu2[final, ]), feature.names),
    covariance = core$covariance,
    beta = setNames(as.numeric(core$beta[final, ]), feature.names)
  )
  support <- abs(estimate$beta) > support_tol
  structure(list(
    valid = TRUE,
    method = paste(
      "CHIME two-component Gaussian-mixture clustering",
      "(quadratic lasso with certified KKT equations)"
    ),
    cluster = setNames(as.integer(core$cluster), row.names),
    score = setNames(as.numeric(core$score), row.names),
    threshold = as.numeric(core$threshold),
    responsibility = setNames(as.numeric(core$responsibility), row.names),
    estimate = estimate,
    support = support,
    lambda = lambda,
    history = list(
      omega = core$omega,
      mu1 = core$mu1,
      mu2 = core$mu2,
      beta = core$beta,
      covariance = core$covariance_history,
      responsibility = core$responsibility_history,
      label.swapped = as.logical(core$label_swapped),
      parameter.change = core$parameter_change,
      lasso.iterations = core$lasso_iterations,
      lasso.converged = as.logical(core$lasso_converged),
      lasso.relative.update = core$lasso_relative_update,
      lasso.kkt = core$lasso_kkt,
      lasso.kkt.scale = core$lasso_kkt_scale,
      lasso.objective = core$lasso_objective
    ),
    diagnostics = list(
      components = 2L,
      iterations = max_iter,
      lasso.tolerance = lasso_tol,
      lasso.max.iterations = lasso_max_iter,
      maximum.scaled.kkt = max(core$lasso_kkt / core$lasso_kkt_scale),
      support.tolerance = support_tol,
      support.size = sum(support),
      initial.covariance = initial$covariance[-1L],
      stable.logistic = TRUE,
      label.canonicalization = paste(
        "At every stage omega > 1/2 is mapped to 1 - omega and component",
        "labels, responsibilities, means, and beta orientation are exchanged."
      ),
      tie.rule = "score >= log(omega / (1 - omega)) selects component 1",
      no.repair = paste(
        "No ridge, inverse, pseudoinverse, eigenvalue floor, probability",
        "clipping, component reset, or true-label tuning was used."
      )
    ),
    call = call,
    data.name = data.name
  ), class = c("chime_fit", "hd_clustering_fit", "list"))
}


.c7_if_adjust <- function(score, method, label) {
  score <- as.numeric(score)
  if (!length(score) || anyNA(score) || any(!is.finite(score))) {
    stop(sprintf("%s scores must be finite and non-empty.", label),
         call. = FALSE)
  }
  if (method == "none") {
    return(list(
      score = score, location = 0, scale = 1,
      rule = "no empirical-null translation or rescaling"
    ))
  }
  if (method == "mean_sd") {
    location <- mean(score)
    scale <- stats::sd(score)
    rule <- "(psi - mean(psi)) / sd(psi), with the R n-1 SD"
  } else {
    location <- stats::median(score)
    scale <- stats::mad(score, center = location, constant = 1.4826)
    rule <- paste(
      "(psi - median(psi)) / MAD(psi),",
      "with normal-consistency constant 1.4826"
    )
  }
  if (!is.finite(scale) || scale <= 0) {
    stop(sprintf(
      "%s scores have zero or non-finite empirical-null scale.", label
    ), call. = FALSE)
  }
  list(
    score = (score - location) / scale,
    location = location,
    scale = scale,
    rule = rule
  )
}


.c7_with_seed <- function(seed, expression) {
  seed <- .c7_count(seed, "seed", minimum = 0L)
  had.state <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had.state) old.state <- get(".Random.seed", envir = .GlobalEnv)
  on.exit({
    if (had.state) {
      assign(".Random.seed", old.state, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv,
                      inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  set.seed(seed)
  force(expression)
}


.c7_if_pvalues <- function(score, null_scores, null_cdf) {
  if (!is.null(null_cdf)) {
    if (!is.function(null_cdf)) {
      stop("'null_cdf' must be a function.", call. = FALSE)
    }
    distribution <- vapply(score, function(value) {
      result <- as.numeric(null_cdf(value))
      if (length(result) != 1L) NA_real_ else result
    }, numeric(1))
    if (anyNA(distribution) || any(!is.finite(distribution)) ||
        any(distribution < 0 | distribution > 1)) {
      stop("'null_cdf' must return one finite value in [0, 1] per score.",
           call. = FALSE)
    }
    return(list(
      p.value = as.numeric(1 - distribution),
      source = "caller-supplied null CDF"
    ))
  }
  null_scores <- as.numeric(null_scores)
  if (!length(null_scores) || anyNA(null_scores) ||
      any(!is.finite(null_scores))) {
    stop("'null_scores' must contain finite values.", call. = FALSE)
  }
  p.value <- vapply(score, function(value) {
    mean(null_scores > value)
  }, numeric(1))
  list(
    p.value = as.numeric(p.value),
    source = paste(
      "empirical 1-F0(t) = mean(null_scores > t),",
      "without an unreported pseudocount"
    )
  )
}


.c7_if_hct <- function(score, p.value, sample_size, convention) {
  score <- as.numeric(score)
  p.value <- as.numeric(p.value)
  p <- length(score)
  feature <- seq_len(p)
  p.order <- order(p.value, feature)
  sorted.p <- p.value[p.order]
  rank <- seq_len(p)
  rank.probability <- if (convention == "paper") {
    rank / p
  } else {
    rank / (p + 1)
  }
  contrast <- rank.probability - sorted.p
  denominator <- sqrt(
    rank.probability + pmax(sqrt(sample_size) * contrast, 0)
  )
  hc <- sqrt(p) * contrast / denominator
  p.value.lower <- log(p) / p
  maximum.rank <- floor(p / 2)
  admissible <- rank <= maximum.rank &
    sorted.p > p.value.lower & is.finite(hc)
  if (!any(admissible)) {
    stop(paste(
      "The IF-PCA HCT search has no admissible rank under",
      "p-value > log(p)/p and rank <= floor(p/2)."
    ), call. = FALSE)
  }
  candidates <- which(admissible)
  maximum <- max(hc[candidates])
  tied <- candidates[hc[candidates] == maximum]
  selected.rank <- as.integer(if (convention == "paper") {
    tied[[1L]]
  } else {
    tied[[length(tied)]]
  })
  score.order <- order(-score, feature)
  selected <- score.order[seq_len(selected.rank)]
  list(
    selected = selected,
    selected.rank = selected.rank,
    threshold = score[score.order[[selected.rank]]],
    score.order = score.order,
    p.value.order = p.order,
    sorted.p.value = sorted.p,
    rank.probability = rank.probability,
    hc = hc,
    admissible = admissible,
    p.value.lower = p.value.lower,
    maximum.rank = maximum.rank,
    convention = convention,
    maximum.tie.rule = if (convention == "paper") {
      "first (smallest) rank attaining the HC maximum"
    } else {
      "last (largest) rank attaining the HC maximum, as in the software"
    },
    feature.tie.rule =
      "decreasing adjusted score, then increasing feature index",
    selection.rule = paste(
      "exactly the top j-hat features; this resolves the paper's",
      "strict-threshold versus j-hat count ambiguity"
    )
  )
}


#' Influential features PCA clustering
#'
#' Implements influential features PCA (IF-PCA): standardize every feature
#' with its n - 1 sample standard deviation, rank features by
#' \deqn{\psi_{n,j}=\sqrt n\sup_t|\widehat F_{n,j}(t)-\Phi(t)|,}
#' optionally translate and rescale the scores by an empirical null, select
#' features, compute the leading K - 1 left singular vectors, and cluster
#' their rows.
#'
#' ks_convention = "paper" evaluates the displayed statistic on the n - 1
#' standardized data. "software" reproduces the distinct official MATLAB
#' score convention by dividing those standardized values once more by
#' sqrt(1 - 1/n), which is equivalent to using divisor n inside the KS
#' calculation. This distinction does not change the PCA matrix. The book's
#' printed n-divisor standardization is not used because the primary method
#' and software both form the PCA matrix with the n - 1 SD.
#'
#' With fixed selection, threshold is applied to the adjusted scores using a
#' greater-than-or-equal rule. HCT requires either null_scores already on the
#' same final scale, a supplied null_cdf, or explicit null_reps plus
#' seed. A standard KS CDF
#' is not silently substituted because centering and scale are estimated.
#' The original software used 100*p null replicates and the paper's numerical
#' study used 2000*p; neither computational choice is hidden as a default.
#'
#' @param x Finite numeric matrix with observations in rows.
#' @param K Known number of clusters, at least two.
#' @param selection Either "fixed" or higher-criticism thresholding ("hct").
#' @param threshold Required non-negative adjusted-score threshold for fixed
#'   selection.
#' @param empirical_null One of "mean_sd" (the primary empirical-null step),
#'   "median_mad" (the paper's robust variant), or "none".
#' @param ks_convention Either the primary-paper or official-software KS scale.
#' @param hct_convention For HCT, "paper" uses rank probability j/p and the
#'   first HC maximum; "software" uses j/(p+1) and the last maximum.
#' @param null_scores Optional finite reference scores already on the same final
#'   scale as the adjusted observed scores.
#' @param null_cdf Optional CDF function on that final scale.
#' @param null_reps Optional explicit number of null Monte Carlo replicates.
#'   There is deliberately no default.
#' @param seed Required when null_reps is supplied. The caller's RNG state is
#'   restored after calibration.
#' @param truncate If TRUE, apply the paper's theoretical entrywise bound
#'   log(p)/sqrt(n) to the left singular vectors. Numerical work in the
#'   primary paper and official software did not use it.
#' @param rank_tol Explicit non-negative relative singular-value tolerance used
#'   only to certify that selected data have rank at least K - 1.
#' @param kmeans_tol Positive tolerance for deterministic Lloyd updates.
#' @param kmeans_max_iter Maximum deterministic Lloyd updates.
#'
#' @return An if_pca_fit object with labels, selected indices and names,
#'   standardized data, raw and adjusted KS scores, the embedding,
#'   deterministic k-means diagnostics, and full HCT search details when used.
#'
#' @references
#' Jin, J. and Wang, W. (2016). Influential features PCA for high dimensional
#' clustering. \emph{Annals of Statistics}, \strong{44}, 2323--2359.
#' \doi{10.1214/15-AOS1423}.
#'
#' @examples
#' x <- cbind(
#'   c(-3, -2.5, -2, -1.5, 1.5, 2, 2.5, 3),
#'   c(-1, 0, 1, 0, -1, 0, 1, 0),
#'   c(0, 1, 0, -1, 0, 1, 0, -1)
#' )
#' if_pca(x, K = 2, selection = "fixed", threshold = 0,
#'        empirical_null = "none")
#'
#' @export
if_pca <- function(
    x, K, selection = c("fixed", "hct"), threshold = NULL,
    empirical_null = c("mean_sd", "median_mad", "none"),
    ks_convention = c("paper", "software"),
    hct_convention = c("paper", "software"), null_scores = NULL,
    null_cdf = NULL, null_reps = NULL, seed = NULL, truncate = FALSE,
    rank_tol = sqrt(.Machine$double.eps), kmeans_tol = 1e-10,
    kmeans_max_iter = 100L) {
  call <- match.call()
  data.name <- deparse(substitute(x))
  x <- .c7_data(x, min_rows = 3L, min_columns = 2L)
  n <- nrow(x)
  p <- ncol(x)
  K <- .c7_count(K, "K", minimum = 2L)
  if (K >= n) {
    stop("'K' must be smaller than 'nrow(x)'.", call. = FALSE)
  }
  selection <- match.arg(selection)
  empirical_null <- match.arg(empirical_null)
  ks_convention <- match.arg(ks_convention)
  hct_convention <- match.arg(hct_convention)
  truncate <- .c7_flag(truncate, "truncate")
  rank_tol <- .c7_scalar(rank_tol, "rank_tol", 0)
  kmeans_tol <- .c7_scalar(
    kmeans_tol, "kmeans_tol", 0, strict_lower = TRUE
  )
  kmeans_max_iter <- .c7_count(kmeans_max_iter, "kmeans_max_iter")

  score.core <- cpp_ch7_if_scores(x, ks_convention == "software")
  adjustment <- .c7_if_adjust(
    score.core$score, empirical_null, "Observed IF-PCA"
  )
  adjusted <- adjustment$score
  feature <- seq_len(p)
  ranking <- order(-adjusted, feature)
  hct <- NULL
  null.details <- NULL

  if (selection == "fixed") {
    if (!is.null(null_scores) || !is.null(null_cdf) ||
        !is.null(null_reps) || !is.null(seed)) {
      stop(paste(
        "Null calibration arguments are only used with selection = 'hct';",
        "none are silently ignored."
      ), call. = FALSE)
    }
    if (is.null(threshold)) {
      stop("'threshold' is required for fixed IF-PCA selection.",
           call. = FALSE)
    }
    threshold <- .c7_scalar(threshold, "threshold", 0)
    selected <- which(adjusted >= threshold)
    selection.rule <- "adjusted score >= caller-supplied threshold"
  } else {
    if (!is.null(threshold)) {
      stop("'threshold' must be NULL when selection = 'hct'.",
           call. = FALSE)
    }
    sources <- c(
      !is.null(null_scores), !is.null(null_cdf), !is.null(null_reps)
    )
    if (sum(sources) != 1L) {
      stop(paste(
        "HCT needs exactly one of 'null_scores', 'null_cdf', or",
        "explicit 'null_reps' plus 'seed'."
      ), call. = FALSE)
    }
    if (!is.null(null_reps)) {
      minimum <- if (empirical_null == "none") 1L else 2L
      null_reps <- .c7_count(null_reps, "null_reps", minimum = minimum)
      if (is.null(seed)) {
        stop("'seed' is required when 'null_reps' is supplied.",
             call. = FALSE)
      }
      raw.null <- .c7_with_seed(seed, cpp_ch7_if_null_scores(
        n, null_reps, ks_convention == "software"
      ))
      null.adjustment <- .c7_if_adjust(
        raw.null, empirical_null, "Simulated IF-PCA null"
      )
      null_scores <- null.adjustment$score
      null.details <- list(
        source = "explicit standard-normal Monte Carlo",
        repetitions = null_reps,
        seed = as.integer(seed),
        raw.scores = raw.null,
        adjustment = null.adjustment[-1L],
        rng.state.restored = TRUE,
        software.reference = "100*p repetitions (not used implicitly)",
        paper.numeric.reference = "2000*p repetitions (not used implicitly)"
      )
    } else {
      if (!is.null(seed)) {
        stop("'seed' is only meaningful with 'null_reps'.",
             call. = FALSE)
      }
      if (!is.null(null_scores)) null_scores <- as.numeric(null_scores)
    }
    p.values <- .c7_if_pvalues(adjusted, null_scores, null_cdf)
    hct <- .c7_if_hct(adjusted, p.values$p.value, n, hct_convention)
    hct$p.value <- p.values$p.value
    hct$null.source <- p.values$source
    selected <- hct$selected
    threshold <- hct$threshold
    ranking <- hct$score.order
    selection.rule <- hct$selection.rule
  }

  if (!length(selected)) {
    stop(paste(
      "IF-PCA selected no features; no automatic threshold relaxation or",
      "fallback to ordinary PCA was applied."
    ), call. = FALSE)
  }
  selected <- as.integer(selected)
  selected.matrix <- score.core$standardized[, selected, drop = FALSE]
  decomposition <- svd(selected.matrix, nu = K - 1L, nv = 0L)
  singular <- as.numeric(decomposition$d)
  singular.scale <- max(1, if (length(singular)) singular[[1L]] else 0)
  certified.rank <- sum(singular > rank_tol * singular.scale)
  if (certified.rank < K - 1L) {
    stop(sprintf(
      paste(
        "Selected IF-PCA data have certified rank %d, below K - 1 = %d",
        "at 'rank_tol'; no jitter or extra feature was added."
      ), certified.rank, K - 1L
    ), call. = FALSE)
  }
  embedding <- decomposition$u[, seq_len(K - 1L), drop = FALSE]
  for (component in seq_len(ncol(embedding))) {
    anchor <- which.max(abs(embedding[, component]))
    if (embedding[anchor, component] < 0) {
      embedding[, component] <- -embedding[, component]
    }
  }
  truncation.threshold <- Inf
  if (truncate) {
    truncation.threshold <- log(p) / sqrt(n)
    embedding[] <- pmax(
      -truncation.threshold, pmin(truncation.threshold, embedding)
    )
  }
  clustering <- cpp_ch7_if_deterministic_kmeans(
    embedding, K, kmeans_tol, kmeans_max_iter
  )
  if (!isTRUE(clustering$valid)) {
    stop(clustering$message, call. = FALSE)
  }
  if (!isTRUE(clustering$converged)) {
    stop(paste(
      "Deterministic IF-PCA k-means did not converge; no last iterate is",
      "reported as a fit."
    ), call. = FALSE)
  }

  feature.names <- colnames(x)
  if (is.null(feature.names)) feature.names <- paste0("V", feature)
  row.names <- rownames(x)
  dimnames(score.core$standardized) <- list(row.names, feature.names)
  names(score.core$center) <- feature.names
  names(score.core$scale) <- feature.names
  score.core$score <- setNames(as.numeric(score.core$score), feature.names)
  names(adjusted) <- feature.names
  rownames(embedding) <- row.names
  colnames(embedding) <- paste0("IFPC", seq_len(K - 1L))
  labels <- as.integer(clustering$cluster)
  names(labels) <- row.names

  structure(list(
    valid = TRUE,
    method = "Influential features PCA clustering",
    cluster = labels,
    selected = selected,
    selected.names = feature.names[selected],
    ranking = ranking,
    raw.score = score.core$score,
    adjusted.score = adjusted,
    standardized = score.core$standardized,
    embedding = embedding,
    singular.values = singular,
    threshold = threshold,
    selection = selection,
    empirical.null = list(
      method = empirical_null,
      location = adjustment$location,
      scale = adjustment$scale,
      rule = adjustment$rule
    ),
    hct = hct,
    null.calibration = null.details,
    kmeans = clustering,
    diagnostics = list(
      sample.size = n,
      features = p,
      clusters = K,
      standard.deviation.divisor = n - 1L,
      ks.convention = ks_convention,
      software.ks.rescaling = if (ks_convention == "software") {
        sprintf(
          "divide W by sqrt(1 - 1/n) = %.17g for KS only",
          score.core$software_ks_scale
        )
      } else {
        "none"
      },
      selection.rule = selection.rule,
      feature.tie.rule =
        "decreasing adjusted score, then increasing feature index",
      selected.rank = length(selected),
      certified.rank = certified.rank,
      rank.tolerance = rank_tol * singular.scale,
      truncation = truncate,
      truncation.threshold = truncation.threshold,
      kmeans.initialization = clustering$initialization,
      kmeans.tie.rule = clustering$tie_rule,
      no.repair = paste(
        "No zero-variance feature deletion, rank jitter, threshold",
        "relaxation, random center reset, or implicit null law was used."
      )
    ),
    call = call,
    data.name = data.name
  ), class = c("if_pca_fit", "hd_clustering_fit", "list"))
}
