# Chapter 3: high-dimensional Hettmansperger--Randles estimation.

.hdhr_data_matrix <- function(x) {
  if (is.data.frame(x)) {
    if (!all(vapply(x, is.numeric, logical(1)))) {
      stop("Every column of `x` must be numeric.", call. = FALSE)
    }
    x <- as.matrix(x)
  }
  if (!is.matrix(x) || !is.numeric(x)) {
    stop("`x` must be a numeric matrix or all-numeric data frame.",
         call. = FALSE)
  }
  storage.mode(x) <- "double"
  if (nrow(x) < 2L || ncol(x) < 1L) {
    stop("`x` must have at least two rows and one column.", call. = FALSE)
  }
  if (anyNA(x) || any(!is.finite(x))) {
    stop("`x` must contain only finite values.", call. = FALSE)
  }
  x
}

.hdhr_positive_scalar <- function(x, name, allow_zero = FALSE) {
  x <- as.numeric(x)
  lower_ok <- if (allow_zero) x >= 0 else x > 0
  if (length(x) != 1L || is.na(x) || !is.finite(x) || !lower_ok) {
    qualifier <- if (allow_zero) "non-negative" else "positive"
    stop(sprintf("`%s` must be one finite %s number.", name, qualifier),
         call. = FALSE)
  }
  x
}

.hdhr_positive_integer <- function(x, name) {
  value <- as.numeric(x)
  if (length(value) != 1L || is.na(value) || !is.finite(value) ||
      value < 1 || value != floor(value) || value > .Machine$integer.max) {
    stop(sprintf("`%s` must be one positive integer.", name),
         call. = FALSE)
  }
  as.integer(value)
}

.hdhr_strict <- function(strict) {
  if (!is.logical(strict) || length(strict) != 1L || is.na(strict)) {
    stop("`strict` must be TRUE or FALSE.", call. = FALSE)
  }
  strict
}

.hdhr_restore_from_log <- function(log_value) {
  smallest_log <- log(.Machine$double.xmin) -
    (.Machine$double.digits - 1) * log(2)
  largest_log <- log(.Machine$double.xmax)
  if (!is.finite(log_value) || log_value < smallest_log ||
      log_value > largest_log) {
    return(NA_real_)
  }
  value <- exp(log_value)
  if (!is.finite(value) || value <= 0) NA_real_ else value
}

.hdhr_invalid <- function(message, strict, call, data.name,
                          bandwidth, scale.estimator, pilot.source,
                          stage, core = NULL, pilot.certificate = NULL) {
  if (strict) {
    stop(message, call. = FALSE)
  }
  warning(message, call. = FALSE)
  last.iterate <- NULL
  if (is.list(core)) {
    last.iterate <- list(
      location.scaled = core$location_scaled,
      shape = core$shape_last,
      raw.shape = core$raw_shape_last,
      sscm = core$sscm_last,
      banded.sscm = core$banded_sscm_last
    )
  }
  structure(
    list(
      estimate = NULL,
      valid = FALSE,
      data.name = data.name,
      bandwidth = bandwidth,
      scale.estimator = scale.estimator,
      diagnostics = list(
        failure.stage = stage,
        failure = message,
        iteration.stable = FALSE,
        iterations = if (is.list(core) && !is.null(core$iterations))
          as.integer(core$iterations) else 0L,
        median = if (is.list(core)) core$median else NULL,
        map = if (is.list(core)) core$map else NULL,
        last.iterate = last.iterate,
        pilot.source = pilot.source,
        pilot.certificate = pilot.certificate,
        kkt = list(
          applicable.to.hdhr = FALSE,
          pilot = pilot.certificate,
          note = paste(
            "Algorithm 2 is a fixed-point iteration, not a KKT solver;",
            "a supplied precision fit retains its own certificate here."
          )
        ),
        no.repair = paste(
          "No ridge, eigenvalue floor, pseudoinverse, jitter, zero-residual",
          "perturbation, or post-hoc positive-definite repair was applied."
        ),
        call = call
      )
    ),
    class = "high_dimensional_hr_fit"
  )
}

.hdhr_extract_pilot <- function(pilot_precision, p) {
  if (inherits(pilot_precision, "spatial_sign_precision_fit")) {
    if (!isTRUE(pilot_precision$valid) || is.null(pilot_precision$estimate)) {
      return(list(
        valid = FALSE,
        message = paste(
          "`pilot_precision` is an invalid or uncertified",
          "spatial_sign_precision() fit."
        ),
        source = "spatial_sign_precision_fit",
        certificate = pilot_precision$diagnostics$solver
      ))
    }
    matrix <- pilot_precision$estimate
    source <- paste0("spatial_sign_precision_fit:", pilot_precision$method)
    certificate <- pilot_precision$diagnostics$solver
  } else {
    matrix <- pilot_precision
    source <- "user-supplied matrix"
    certificate <- NULL
  }
  if (!is.matrix(matrix) || !is.numeric(matrix)) {
    stop(
      paste(
        "`pilot_precision` must be a numeric p by p matrix or a valid",
        "spatial_sign_precision() fit."
      ),
      call. = FALSE
    )
  }
  storage.mode(matrix) <- "double"
  if (!identical(dim(matrix), c(p, p))) {
    stop("`pilot_precision` must be a p by p matrix matching `x`.",
         call. = FALSE)
  }
  if (anyNA(matrix) || any(!is.finite(matrix))) {
    stop("`pilot_precision` must contain only finite values.",
         call. = FALSE)
  }
  list(
    valid = TRUE,
    matrix = matrix,
    source = source,
    certificate = certificate
  )
}


#' High-dimensional Hettmansperger--Randles estimation
#'
#' Implements Algorithm 2 of Yan, Feng, and Zhang (2025).  Starting from the
#' sample spatial median and a supplied positive-definite pilot precision, the
#' method jointly updates location and trace-normalized shape.  At iteration
#' \eqn{k}, with
#' \eqn{e_i^{(k)}=\Sigma_k^{-1/2}(X_i-\mu_k)}, it applies
#' \deqn{\mu_{k+1}=\mu_k+\Sigma_k^{1/2}
#'   \frac{n^{-1}\sum_i U(e_i^{(k)})}
#'        {n^{-1}\sum_i\|e_i^{(k)}\|^{-1}}}
#' and
#' \deqn{\Sigma_{k+1}\ \mathrel{\propto}\
#'   \Sigma_k^{1/2}\mathcal B_h\left\{
#'   n^{-1}\sum_i U(e_i^{(k)})U(e_i^{(k)})^T\right\}
#'   \Sigma_k^{1/2},\qquad \operatorname{tr}(\Sigma_{k+1})=p.}
#'
#' The paper sets \eqn{h=3}; the formula permits \eqn{0\le h<p}.  When the
#' default is not explicitly supplied and \eqn{p\le3}, it resolves to
#' `min(3, p - 1)` and records that value.  The pilot's common scalar is not
#' identified by spatial signs, so its inverse is normalized to trace \eqn{p}
#' before iteration.
#'
#' The primary estimator is a shape estimator, not a covariance-scale
#' estimator.  With `scale_estimator = "paper_qda"`, this function additionally
#' implements the paper's QDA scale add-on
#' \deqn{\widehat{\operatorname{tr}(\Xi)}=
#'   \{\sum_i\|X_i\|^2-n\|\bar X\|^2\}/(n-1),}
#' returning `scatter.scale = trace.hat / p` and
#' `scatter = scatter.scale * shape`.  This second-moment scale is not part of
#' the robust sign fixed point and requires its finite double representation.
#'
#' The paper says to repeat until convergence but does not prescribe a norm.
#' Here, convergence means that the maximum of the relative location and
#' Frobenius shape-map updates is at most `tol`.  The spatial-sign score and
#' shape-map residual are returned separately; they are diagnostics rather
#' than an independently claimed estimating-equation certificate.
#'
#' A supplied `spatial_sign_precision()` fit is accepted only when its
#' feasibility/KKT certificate made it a valid fit.  A raw matrix must be
#' symmetric positive definite.  Hard banding itself can destroy positive
#' definiteness.  Such a failure, a singular pilot, an exact zero standardized
#' residual, overflow, or nonconvergence is reported without a ridge,
#' eigenvalue floor, pseudoinverse, jitter, or perturbation.  With
#' `strict = FALSE`, every such method failure returns `estimate = NULL` and
#' `valid = FALSE`; the last iterate appears only under diagnostics.
#'
#' The book's fixed-pilot score, post-hoc banded raw shape, second banding of
#' the inverse, and HR-centered divisor-\eqn{n} scale are not Algorithm 2.
#' Likewise, the book's displayed \eqn{r_n+h^{-\alpha}} theorem is a review
#' synthesis rather than a finite-sample guarantee stated in this primary.
#'
#' @param x Numeric \eqn{n\times p} data matrix.
#' @param pilot_precision A symmetric positive-definite \eqn{p\times p} pilot
#'   precision matrix, or a valid object returned by
#'   [spatial_sign_precision()].  No automatic tuning constant is invented.
#' @param bandwidth Non-negative integer less than \eqn{p}.  The paper uses 3.
#'   If omitted, the resolved default is `min(3, p - 1)`.
#' @param tol Positive relative tolerance for the joint location/shape fixed
#'   point.
#' @param max_iter Positive maximum number of joint updates.
#' @param median_tol Positive relative tolerance for the initial sample
#'   spatial median.
#' @param median_max_iter Positive maximum number of spatial-median updates.
#' @param zero_tol Non-negative radius threshold treated as coincidence.  The
#'   default zero applies no positive floor.
#' @param scale_estimator Either `"paper_qda"` for the primary paper's QDA
#'   covariance-trace add-on, or `"none"` for shape-only output.
#' @param strict If `TRUE`, method failures are errors.  If `FALSE`, they are
#'   warnings followed by an explicitly invalid fit with `estimate = NULL`.
#'
#' @return An object of class `high_dimensional_hr_fit`.  A valid object returns
#'   robust `location`, the final trace-\eqn{p} `shape`, the unbanded
#'   `raw.shape` fixed-point map, `precision = solve(shape)`, the final SSCMs,
#'   optional covariance `scatter.scale` and `scatter`, and detailed median,
#'   iteration, score, equation, positive-definiteness, reciprocal-condition,
#'   scale-identification, and pilot-certificate diagnostics.
#'
#' @references Yan, G., Feng, L., and Zhang, X. (2025).
#'   *High-Dimensional Hettmansperger-Randles Estimator and its Applications*.
#'   arXiv:2505.01669. \url{https://arxiv.org/abs/2505.01669}
#'
#' @examples
#' x <- rbind(
#'   c(2, 0), c(-2, 0), c(0, 1), c(0, -1),
#'   c(1, 1), c(-1, -1), c(1, -1), c(-1, 1)
#' )
#' high_dimensional_hr(x, diag(2), bandwidth = 0)
#'
#' @export
high_dimensional_hr <- function(
    x, pilot_precision, bandwidth = 3L,
    tol = 1e-8, max_iter = 1000L,
    median_tol = 1e-8, median_max_iter = 1000L,
    zero_tol = 0,
    scale_estimator = c("paper_qda", "none"),
    strict = TRUE) {
  call <- match.call()
  x.name <- deparse1(substitute(x))
  bandwidth.was.missing <- missing(bandwidth)
  x <- .hdhr_data_matrix(x)
  p <- ncol(x)
  strict <- .hdhr_strict(strict)
  scale.estimator <- match.arg(scale_estimator)
  tol <- .hdhr_positive_scalar(tol, "tol")
  max_iter <- .hdhr_positive_integer(max_iter, "max_iter")
  median_tol <- .hdhr_positive_scalar(median_tol, "median_tol")
  median_max_iter <- .hdhr_positive_integer(
    median_max_iter, "median_max_iter"
  )
  zero_tol <- .hdhr_positive_scalar(zero_tol, "zero_tol", allow_zero = TRUE)

  requested.bandwidth <- if (bandwidth.was.missing) 3L else bandwidth
  if (bandwidth.was.missing) {
    resolved.bandwidth <- min(3L, p - 1L)
  } else {
    value <- as.numeric(bandwidth)
    if (length(value) != 1L || is.na(value) || !is.finite(value) ||
        value < 0 || value != floor(value) || value >= p) {
      stop("`bandwidth` must be one integer from zero through p - 1.",
           call. = FALSE)
    }
    resolved.bandwidth <- as.integer(value)
  }
  bandwidth.info <- list(
    requested = as.integer(requested.bandwidth),
    resolved = resolved.bandwidth,
    default.resolved.for.dimension = bandwidth.was.missing && p <= 3L,
    paper.value = 3L
  )

  pilot <- .hdhr_extract_pilot(pilot_precision, p)
  if (!pilot$valid) {
    return(.hdhr_invalid(
      pilot$message, strict, call, x.name, bandwidth.info,
      scale.estimator, pilot$source, "pilot precision fit",
      pilot.certificate = pilot$certificate
    ))
  }

  core <- tryCatch(
    cpp_ch3_hdhr_fit(
      x, pilot$matrix, resolved.bandwidth,
      median_tol, median_max_iter, tol, max_iter, zero_tol
    ),
    error = identity
  )
  if (inherits(core, "condition")) {
    message <- paste0(
      "HDHR computation failed: ", conditionMessage(core),
      " No numerical repair was applied."
    )
    return(.hdhr_invalid(
      message, strict, call, x.name, bandwidth.info,
      scale.estimator, pilot$source, "C++ computation",
      pilot.certificate = pilot$certificate
    ))
  }
  if (!isTRUE(core$valid)) {
    message <- paste0(
      "HDHR estimation failed at ", core$failure_stage, ": ", core$failure
    )
    return(.hdhr_invalid(
      message, strict, call, x.name, bandwidth.info,
      scale.estimator, pilot$source, core$failure_stage, core,
      pilot$certificate
    ))
  }

  log.scatter.scale <- 2 * as.numeric(core$log_data_scale) +
    log(as.numeric(core$qda_trace_scaled)) - log(p)
  if (scale.estimator == "paper_qda") {
    scatter.scale <- .hdhr_restore_from_log(log.scatter.scale)
    if (is.na(scatter.scale)) {
      message <- paste(
        "HDHR shape estimation succeeded, but the primary paper's QDA",
        "scatter scale is outside the finite positive double range; no",
        "clipping or rescaling was applied."
      )
      return(.hdhr_invalid(
        message, strict, call, x.name, bandwidth.info,
        scale.estimator, pilot$source, "QDA scatter scale", core,
        pilot$certificate
      ))
    }
    scatter <- as.matrix(core$shape) * scatter.scale
    if (any(!is.finite(scatter)) || any(diag(scatter) <= 0)) {
      message <- paste(
        "The covariance-scale HDHR scatter cannot be represented as a",
        "finite positive-diagonal double matrix; no repair was applied."
      )
      return(.hdhr_invalid(
        message, strict, call, x.name, bandwidth.info,
        scale.estimator, pilot$source, "scatter restoration", core,
        pilot$certificate
      ))
    }
    scatter.trace <- p * scatter.scale
  } else {
    scatter.scale <- NA_real_
    scatter.trace <- NA_real_
    scatter <- NULL
  }

  location <- as.numeric(core$location)
  shape <- as.matrix(core$shape)
  raw.shape <- as.matrix(core$raw_shape)
  precision <- as.matrix(core$precision)
  sscm <- as.matrix(core$sscm)
  banded.sscm <- as.matrix(core$banded_sscm)
  signs <- as.matrix(core$signs)
  variable.names <- colnames(x)
  observation.names <- rownames(x)
  if (!is.null(variable.names)) {
    names(location) <- variable.names
    for (name in c("shape", "raw.shape", "precision", "sscm",
                   "banded.sscm")) {
      value <- get(name)
      dimnames(value) <- list(variable.names, variable.names)
      assign(name, value)
    }
    colnames(signs) <- variable.names
  }
  if (!is.null(observation.names)) {
    rownames(signs) <- observation.names
  }
  if (!is.null(scatter) && !is.null(variable.names)) {
    dimnames(scatter) <- list(variable.names, variable.names)
  }

  minimum.radius.log <- as.numeric(core$log_data_scale) +
    log(as.numeric(core$map$minimum.residual.distance.scaled))
  minimum.radius <- .hdhr_restore_from_log(minimum.radius.log)
  map <- core$map
  shape.equation.residual <- as.numeric(map$shape.relative.update)
  location.equation.residual <- list(
    l2 = as.numeric(map$score.residual.l2),
    infinity = as.numeric(map$score.residual.infinity)
  )
  estimate <- list(
    location = location,
    raw.shape = raw.shape,
    shape = shape,
    precision = precision,
    scatter.scale = scatter.scale,
    scatter = scatter
  )

  structure(
    list(
      estimate = estimate,
      valid = TRUE,
      location = location,
      raw.shape = raw.shape,
      shape = shape,
      precision = precision,
      scatter.scale = scatter.scale,
      scatter.trace = scatter.trace,
      scatter = scatter,
      sscm = sscm,
      banded.sscm = banded.sscm,
      signs = signs,
      bandwidth = bandwidth.info,
      scale.estimator = scale.estimator,
      data.name = x.name,
      diagnostics = list(
        failure.stage = NULL,
        failure = NULL,
        iteration.stable = isTRUE(core$iteration_stable),
        iterations = as.integer(core$iterations),
        convergence.basis = paste(
          "maximum relative location/shape fixed-point update",
          "as a documented software stopping rule"
        ),
        relative.update = as.numeric(map$relative.update),
        location.relative.update =
          as.numeric(map$location.relative.update),
        shape.relative.update = as.numeric(map$shape.relative.update),
        score.residual = location.equation.residual,
        shape.equation.residual = shape.equation.residual,
        score.certificate = paste(
          "recorded diagnostic; stability is certified by the relative",
          "fixed-point update, because the paper gives no separate tolerance"
        ),
        minimum.residual.distance = minimum.radius,
        minimum.residual.distance.log = minimum.radius.log,
        minimum.residual.distance.scaled =
          as.numeric(map$minimum.residual.distance.scaled),
        spatial.median = core$median,
        positive.definiteness = list(
          pilot.minimum.eigenvalue =
            as.numeric(core$pilot_minimum_eigenvalue),
          pilot.rcond = as.numeric(core$pilot_rcond),
          banded.sscm.minimum.eigenvalue =
            as.numeric(map$banded.sscm.minimum.eigenvalue),
          banded.sscm.rcond = as.numeric(map$banded.sscm.rcond),
          shape.minimum.eigenvalue =
            as.numeric(core$shape_minimum_eigenvalue),
          shape.rcond = as.numeric(core$shape_rcond)
        ),
        kkt = list(
          applicable.to.hdhr = FALSE,
          pilot = pilot$certificate,
          note = paste(
            "Algorithm 2 is a fixed-point iteration, not a KKT solver;",
            "the supplied precision fit retains its own certificate."
          )
        ),
        pilot = list(
          source = pilot$source,
          certificate = pilot$certificate,
          common.scale.normalized = TRUE,
          resulting.initial.shape.trace = p
        ),
        scale = list(
          identification = paste(
            "the sign equations identify a trace-p shape; covariance scale",
            "is a separate QDA second-moment add-on"
          ),
          estimator = scale.estimator,
          formula = if (scale.estimator == "paper_qda")
            "{sum ||Xi||^2 - n ||xbar||^2}/[(n-1)p]" else NULL,
          log.scatter.scale = log.scatter.scale,
          qda.trace.scaled = as.numeric(core$qda_trace_scaled),
          robust = FALSE
        ),
        primary.algorithm = paste(
          "Algorithm 2: joint location/shape fixed point; band the",
          "standardized-sign SSCM inside every shape update; normalize trace p"
        ),
        book.errata = c(
          "The fixed-pilot score/raw-SSCM template is not Algorithm 2.",
          paste(
            "Primary banding is inside every SSCM update; it does not band",
            "the final inverse a second time."
          ),
          paste(
            "Primary Algorithm 2 estimates trace-p shape; its QDA scale uses",
            "the sample mean and divisor n-1, not the HR center and divisor n."
          ),
          paste(
            "The book's r_n + h^{-alpha} theorem is not stated by the",
            "primary and is not a finite-sample guarantee."
          )
        ),
        invariance = paste(
          "translation and global scale; hard coordinate banding preserves",
          "signed-coordinate and order-reversal structure, not arbitrary",
          "affine transformations unless bandwidth = p - 1"
        ),
        no.repair = paste(
          "No ridge, eigenvalue floor, pseudoinverse, jitter, zero-residual",
          "perturbation, or post-hoc positive-definite repair was applied."
        ),
        call = call
      )
    ),
    class = "high_dimensional_hr_fit"
  )
}
