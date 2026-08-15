# Chapter 4: goodness-of-fit testing through radial--directional dependence.

.zfrd_logical <- function(value, name) {
  if (!is.logical(value) || length(value) != 1L || is.na(value)) {
    stop(sprintf("`%s` must be TRUE or FALSE.", name), call. = FALSE)
  }
  value
}

.zfrd_alpha <- function(alpha) {
  alpha <- as.numeric(alpha)
  if (length(alpha) != 1L || is.na(alpha) || !is.finite(alpha) ||
      alpha <= 0 || alpha >= 1) {
    stop("`alpha` must be one finite number strictly between zero and one.",
         call. = FALSE)
  }
  alpha
}

.zfrd_integer <- function(value, name, minimum = 1L) {
  value <- as.numeric(value)
  if (length(value) != 1L || is.na(value) || !is.finite(value) ||
      value < minimum || value != floor(value) ||
      value > .Machine$integer.max) {
    stop(sprintf("`%s` must be one integer at least %d.", name, minimum),
         call. = FALSE)
  }
  as.integer(value)
}

.zfrd_seed <- function(seed) {
  if (is.null(seed)) return(NULL)
  .zfrd_integer(seed, "seed", minimum = 0L)
}

.zfrd_relative_zero_tol <- function(value) {
  value <- as.numeric(value)
  if (length(value) != 1L || is.na(value) || !is.finite(value) ||
      value < 0 || value >= 1) {
    stop("`radius_zero_tol` must be one finite number in [0, 1).",
         call. = FALSE)
  }
  value
}

.zfrd_with_local_seed <- function(seed, expression) {
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

.zfrd_shape <- function(shape, p) {
  if (!is.matrix(shape) || !is.numeric(shape) ||
      !identical(dim(shape), c(p, p))) {
    stop("`shape` must be a numeric p by p matrix matching `x`.",
         call. = FALSE)
  }
  storage.mode(shape) <- "double"
  if (anyNA(shape) || any(!is.finite(shape))) {
    stop("`shape` must contain only finite values.", call. = FALSE)
  }
  scale <- max(abs(shape))
  if (!is.finite(scale) || scale <= 0) {
    stop("`shape` must have a finite positive numerical scale.",
         call. = FALSE)
  }
  relative.asymmetry <- max(abs(shape / scale - t(shape) / scale))
  if (!is.finite(relative.asymmetry) || relative.asymmetry > 1e-10) {
    stop("`shape` must be symmetric to relative tolerance 1e-10.",
         call. = FALSE)
  }
  # Average only a permitted roundoff-level asymmetry. Division before
  # addition avoids an unnecessary overflow for large finite entries.
  symmetrized <- (shape / scale / 2 + t(shape) / scale / 2) * scale
  list(matrix = symmetrized, relative.asymmetry = relative.asymmetry)
}

.zfrd_standardisation <- function(
    x, fit, location, shape, pilot_precision, bandwidth,
    tol, max_iter, median_tol, median_max_iter, zero_tol, strict) {
  p <- ncol(x)
  has.fit <- !is.null(fit)
  has.location <- !is.null(location)
  has.shape <- !is.null(shape)
  has.pilot <- !is.null(pilot_precision)
  if (has.fit && (has.location || has.shape || has.pilot)) {
    stop(paste(
      "Supply exactly one standardisation route: `fit`, both `location`",
      "and `shape`, or `pilot_precision`."
    ), call. = FALSE)
  }
  if (!is.null(bandwidth) && (has.fit || has.location || has.shape)) {
    stop(paste(
      "`bandwidth` is used only when standardisation is requested through",
      "`pilot_precision`."
    ), call. = FALSE)
  }
  if (xor(has.location, has.shape)) {
    stop("`location` and `shape` must be supplied together.",
         call. = FALSE)
  }
  if ((has.location || has.shape) && has.pilot) {
    stop(paste(
      "A supplied `location`/`shape` pair cannot be combined with",
      "`pilot_precision`."
    ), call. = FALSE)
  }

  if (has.fit) {
    if (!inherits(fit, "high_dimensional_hr_fit") ||
        !isTRUE(fit$valid) || is.null(fit$location) ||
        is.null(fit$shape)) {
      stop("`fit` must be a valid object returned by high_dimensional_hr().",
           call. = FALSE)
    }
    location <- .as_location(fit$location, p, "fit$location")
    checked.shape <- .zfrd_shape(fit$shape, p)
    return(list(
      location = location,
      shape = checked.shape$matrix,
      source = "validated high_dimensional_hr_fit",
      fit = fit,
      relative.asymmetry = checked.shape$relative.asymmetry
    ))
  }

  if (has.location && has.shape) {
    location <- .as_location(location, p, "location")
    checked.shape <- .zfrd_shape(shape, p)
    return(list(
      location = location,
      shape = checked.shape$matrix,
      source = "user-supplied location and shape",
      fit = NULL,
      relative.asymmetry = checked.shape$relative.asymmetry
    ))
  }

  if (!has.pilot) {
    stop(paste(
      "Supply a valid `fit`, both `location` and `shape`, or an explicit",
      "`pilot_precision` for high_dimensional_hr()."
    ), call. = FALSE)
  }
  arguments <- list(
    x = x,
    pilot_precision = pilot_precision,
    tol = tol,
    max_iter = max_iter,
    median_tol = median_tol,
    median_max_iter = median_max_iter,
    zero_tol = zero_tol,
    scale_estimator = "none",
    strict = strict
  )
  if (!is.null(bandwidth)) arguments$bandwidth <- bandwidth
  fitted <- do.call(high_dimensional_hr, arguments)
  if (!isTRUE(fitted$valid) || is.null(fitted$location) ||
      is.null(fitted$shape)) {
    stop(paste(
      "The internally requested high_dimensional_hr() standardisation is",
      "invalid; no test statistic is returned."
    ), call. = FALSE)
  }
  checked.shape <- .zfrd_shape(fitted$shape, p)
  list(
    location = .as_location(fitted$location, p, "fit$location"),
    shape = checked.shape$matrix,
    source = "internally fitted high_dimensional_hr() with supplied pilot",
    fit = fitted,
    relative.asymmetry = checked.shape$relative.asymmetry
  )
}

.zfrd_normal_tail <- function(z) {
  list(
    p.value = stats::pnorm(z, lower.tail = FALSE),
    log.p.value = stats::pnorm(z, lower.tail = FALSE, log.p = TRUE),
    log.one.minus.p.value = stats::pnorm(
      z, lower.tail = TRUE, log.p = TRUE
    )
  )
}

.zfrd_calibration <- function(T.sum, T.max, p, kind,
                              bootstrap = NULL) {
  if (kind == "analytic") {
    sum.argument <- (T.sum - p) / sqrt(2 * p)
    max.argument <- T.max
    details <- list(
      kind = "analytic",
      sum.center = p,
      sum.scale = sqrt(2 * p),
      max.location = NA_real_,
      max.scale = NA_real_
    )
  } else {
    sd.sum <- as.numeric(bootstrap$sd_sum)
    sd.max <- as.numeric(bootstrap$sd_max)
    if (!is.finite(sd.sum) || sd.sum <= 0 ||
        !is.finite(sd.max) || sd.max <= 0) {
      stop(paste(
        "The radial--directional bootstrap produced a non-positive or",
        "non-finite statistic standard deviation; no floor is applied."
      ), call. = FALSE)
    }
    euler.gamma <- 0.5772156649015328606
    mu.g <- 2 * euler.gamma - log(pi)
    sigma.g <- sqrt(2 * pi^2 / 3)
    sum.argument <- (T.sum - bootstrap$mean_sum) / sd.sum
    max.argument <- mu.g + sigma.g / sd.max *
      (T.max - bootstrap$mean_max)
    details <- list(
      kind = "radial_bootstrap",
      sum.center = as.numeric(bootstrap$mean_sum),
      sum.scale = sd.sum,
      max.center = as.numeric(bootstrap$mean_max),
      max.scale = sd.max,
      limiting.gumbel.mean = mu.g,
      limiting.gumbel.sd = sigma.g
    )
  }
  sum.tail <- .zfrd_normal_tail(sum.argument)
  max.tail <- .ssmax_gumbel_tail(max.argument)
  combination <- .ssmax_cauchy_combine_logtails(
    sum.tail$log.p.value, sum.tail$log.one.minus.p.value,
    max.tail$log.p.value, max.tail$log.one.minus.p.value
  )
  list(
    sum.argument = sum.argument,
    max.argument = max.argument,
    sum = sum.tail,
    max = max.tail,
    combined = combination,
    details = details
  )
}

#' Zhang--Feng radial--directional test of an elliptical model
#'
#' Tests the defining radial--directional implication of an elliptical model
#' after a supplied or explicitly fitted affine standardisation.  If
#' \eqn{Y_i=\widehat\Sigma^{-1/2}(X_i-\widehat\mu)},
#' \eqn{L_i=\log\|Y_i\|}, and \eqn{U_i=Y_i/\|Y_i\|}, the coordinate scores are
#' the ordinary empirical correlations \eqn{\widehat\gamma_j} between
#' \eqn{L_i} and \eqn{U_{ij}}.  The paper's statistics are
#' \deqn{T_{\rm sum}=n\sum_j\widehat\gamma_j^2,\qquad
#' T_{\max}=n\max_j\widehat\gamma_j^2-2\log p+\log\log p.}
#' Their analytic p-values use
#' \eqn{(T_{\rm sum}-p)/\sqrt{2p}\Rightarrow N(0,1)} and
#' \eqn{F_G(t)=\exp\{-\pi^{-1/2}\exp(-t/2)\}}.  The adaptive p-value is the
#' equal-weight Cauchy combination of the two marginal p-values.
#'
#' There are three deliberately exclusive standardisation routes.  Supply a
#' valid `high_dimensional_hr()` object through `fit`; supply both `location`
#' and a symmetric positive-definite `shape`; or supply `pilot_precision` and
#' let this function call `high_dimensional_hr()` with no hidden tuning.  The
#' latter route delegates `bandwidth`, convergence controls, and the strict
#' no-repair contract to that estimator.  The primary paper's numerical
#' section adds a ridge and an incompletely specified positive-definite
#' projection.  Those simulation choices are not package defaults and are not
#' reproduced here.
#'
#' With `calibration = "radial_bootstrap"`, the function implements Section
#' 2.4 of the primary paper: fitted radii are resampled with replacement and
#' paired independently with uniform directions on the sphere.  Bootstrap
#' means and standard deviations correct the normal and Gumbel arguments
#' exactly as displayed in the paper.  This is an intrinsic calibration of the
#' test, not a reproduction of the paper's simulation study.  A degenerate
#' resample or non-positive bootstrap standard deviation is reported without
#' deletion, redrawing, or flooring.
#'
#' The procedure tests the coordinatewise radial--directional correlation
#' restrictions targeted by the paper.  It is not an omnibus finite-sample
#' test against every non-elliptical distribution.  With a consistently
#' transformed supplied shape, translation, global scale, and signed
#' coordinate permutations are exact finite-sample invariances.  A general
#' affine transform induces an orthogonal rotation after standardisation, but
#' the paper's coordinatewise empirical variance normalisation means that
#' neither component is exactly rotation invariant in a finite sample.  The
#' null calibration is asymptotically rotation-compatible under uniform
#' directions.  The max diagnostic and an internally banded HR fit are also
#' explicitly basis/order-sensitive.
#'
#' @param x Numeric \eqn{n\times p} matrix with observations in rows.  At
#'   least three observations and two variables are required.
#' @param fit Optional valid object returned by [high_dimensional_hr()].
#' @param location,shape Optional jointly supplied location vector and
#'   symmetric positive-definite shape matrix.  Shape scale is immaterial.
#' @param pilot_precision Optional explicit pilot precision passed to
#'   [high_dimensional_hr()] when neither `fit` nor `location`/`shape` is
#'   supplied.
#' @param bandwidth Optional hard-banding width for an internally fitted HR
#'   estimator.  `NULL` delegates the dimension-aware paper default to
#'   [high_dimensional_hr()].
#' @param component Which result is the main `htest`: adaptive Cauchy
#'   `"combined"`, dense `"sum"`, or sparse `"max"`.  All three are always
#'   returned under `components`.
#' @param calibration Either the asymptotic `"analytic"` calibration or the
#'   paper's intrinsic `"radial_bootstrap"` mean--variance correction.
#' @param B Number of intrinsic bootstrap replicates; at least two when used.
#' @param seed Optional non-negative integer.  An explicit seed is localized
#'   and does not change the caller's R RNG state.  `NULL` uses the current
#'   stream normally.
#' @param keep_bootstrap Whether to retain the two bootstrap-statistic vectors.
#' @param alpha Test level strictly between zero and one.
#' @param tol,max_iter,median_tol,median_max_iter,zero_tol,strict Controls
#'   passed to an internally requested [high_dimensional_hr()] fit.
#' @param radius_zero_tol Non-negative relative fitted-radius threshold below
#'   one.  A positive value rejects a minimum radius no larger than this
#'   fraction of the maximum radius.  The relative definition preserves the
#'   arbitrary common scale of `shape`; zero rejects only exact coincidences.
#'   No ridge, eigenvalue floor, projection, pseudoinverse, or zero-radius
#'   perturbation is used.
#'
#' @return An object of class `c("radial_directional_test", "htest")` with the
#' selected main p-value, all three analytic or bootstrap-calibrated component
#' results, coordinate correlations and their largest coordinate, fitted log
#' radii and directions, the complete standardisation source, and numerical
#' diagnostics.
#'
#' @references Zhang, H. and Feng, L. (2026). *High-Dimensional Tests for
#' Elliptical Models via Radial--Directional Dependence*. arXiv:2605.03592.
#' \url{https://arxiv.org/abs/2605.03592}
#'
#' @examples
#' x <- rbind(
#'   c(2, 0), c(-2, 0), c(0, 1), c(0, -1),
#'   c(1, 2), c(-1, -2), c(2, -1), c(-2, 1)
#' )
#' zhang_feng_radial_directional_test(
#'   x, location = c(0, 0), shape = diag(2)
#' )
#'
#' @export
zhang_feng_radial_directional_test <- function(
    x, fit = NULL, location = NULL, shape = NULL,
    pilot_precision = NULL, bandwidth = NULL,
    component = c("combined", "sum", "max"),
    calibration = c("analytic", "radial_bootstrap"),
    B = 999L, seed = NULL, keep_bootstrap = FALSE,
    alpha = 0.05,
    tol = 1e-8, max_iter = 1000L,
    median_tol = 1e-8, median_max_iter = 1000L,
    zero_tol = 0, radius_zero_tol = 0, strict = TRUE) {
  call <- match.call()
  data.name <- deparse1(substitute(x))
  x <- .as_data_matrix(x, "x", min_rows = 3L)
  n <- nrow(x)
  p <- ncol(x)
  if (p < 2L) {
    stop("The radial--directional max statistic requires at least two variables.",
         call. = FALSE)
  }
  component <- match.arg(component)
  calibration <- match.arg(calibration)
  alpha <- .zfrd_alpha(alpha)
  strict <- .zfrd_logical(strict, "strict")
  keep_bootstrap <- .zfrd_logical(keep_bootstrap, "keep_bootstrap")
  controls <- .validate_iteration_controls(tol, max_iter, zero_tol)
  median.controls <- .validate_iteration_controls(
    median_tol, median_max_iter, zero_tol
  )
  radius.zero.tol <- .zfrd_relative_zero_tol(radius_zero_tol)
  seed <- .zfrd_seed(seed)
  if (calibration == "analytic") {
    if (!is.null(seed)) {
      stop("`seed` is used only with `calibration = \"radial_bootstrap\"`.",
           call. = FALSE)
    }
    if (keep_bootstrap) {
      stop("`keep_bootstrap` requires radial-bootstrap calibration.",
           call. = FALSE)
    }
    bootstrap.replicates <- 0L
  } else {
    bootstrap.replicates <- .zfrd_integer(B, "B", minimum = 2L)
  }

  standardisation <- .zfrd_standardisation(
    x, fit, location, shape, pilot_precision, bandwidth,
    controls$tol, controls$max_iter,
    median.controls$tol, median.controls$max_iter,
    controls$zero_tol, strict
  )
  run.core <- function() {
    cpp_zhang_feng_radial_directional(
      x, standardisation$location, standardisation$shape,
      radius.zero.tol, bootstrap.replicates
    )
  }
  core <- if (calibration == "radial_bootstrap" && !is.null(seed)) {
    .zfrd_with_local_seed(seed, run.core())
  } else {
    run.core()
  }

  T.sum <- as.numeric(core$T_sum)
  T.max.raw <- as.numeric(core$T_max_raw)
  T.max <- as.numeric(core$T_max)
  analytic <- .zfrd_calibration(T.sum, T.max, p, "analytic")
  selected <- if (calibration == "analytic") {
    analytic
  } else {
    .zfrd_calibration(T.sum, T.max, p, "radial_bootstrap", core$bootstrap)
  }
  correlations <- as.numeric(core$correlations)
  variable.names <- if (is.null(colnames(x))) {
    paste0("variable", seq_len(p))
  } else {
    colnames(x)
  }
  observation.names <- if (is.null(rownames(x))) {
    paste0("observation", seq_len(n))
  } else {
    rownames(x)
  }
  names(correlations) <- variable.names
  directions <- as.matrix(core$directions)
  dimnames(directions) <- list(observation.names, variable.names)
  log.radii <- stats::setNames(
    as.numeric(core$log_radii), observation.names
  )
  largest.index <- which.max(abs(correlations))
  largest.coordinate <- variable.names[[largest.index]]

  p.values <- c(
    sum = selected$sum$p.value,
    max = selected$max$p.value,
    combined = selected$combined$p.value
  )
  log.p.values <- c(
    sum = selected$sum$log.p.value,
    max = selected$max$log.p.value,
    combined = selected$combined$log.p.value
  )
  main.p <- unname(p.values[[component]])
  statistic <- switch(
    component,
    sum = c(Normal.Z = selected$sum.argument),
    max = c(Gumbel.centered = selected$max.argument),
    combined = c(Cauchy.angle = selected$combined$angle)
  )
  raw.statistic <- switch(
    component,
    sum = c(T.sum = T.sum),
    max = c(T.max = T.max),
    combined = c(Cauchy = selected$combined$statistic)
  )
  method.component <- switch(
    component,
    sum = "dense sum component",
    max = "sparse max component",
    combined = "adaptive Cauchy combination"
  )
  bootstrap.values <- if (calibration == "radial_bootstrap" &&
                          keep_bootstrap) {
    list(
      T.sum = as.numeric(core$bootstrap$T_sum),
      T.max = as.numeric(core$bootstrap$T_max)
    )
  } else {
    NULL
  }
  bootstrap.summary <- if (calibration == "radial_bootstrap") {
    list(
      replicates = bootstrap.replicates,
      mean.sum = as.numeric(core$bootstrap$mean_sum),
      sd.sum = as.numeric(core$bootstrap$sd_sum),
      mean.max = as.numeric(core$bootstrap$mean_max),
      sd.max = as.numeric(core$bootstrap$sd_max),
      seed = seed,
      explicit.seed.preserves.caller.RNG.state = !is.null(seed),
      retained = keep_bootstrap,
      values = bootstrap.values
    )
  } else {
    NULL
  }

  result <- list(
    statistic = statistic,
    parameter = c(n = n, p = p),
    p.value = main.p,
    alternative = "greater",
    method = paste(
      "Zhang--Feng radial--directional", method.component,
      sprintf("(%s calibration)", gsub("_", " ", calibration))
    ),
    data.name = data.name,
    raw.statistic = raw.statistic,
    null.value = c(`radial-directional correlations` = 0),
    estimate = correlations,
    components = list(
      correlations = correlations,
      largest.coordinate.index = largest.index,
      largest.coordinate = largest.coordinate,
      largest.absolute.correlation = abs(correlations[[largest.index]]),
      T.sum = T.sum,
      Z.sum = selected$sum.argument,
      T.max.raw = T.max.raw,
      T.max = T.max,
      Gumbel.argument = selected$max.argument,
      Cauchy.statistic = selected$combined$statistic,
      Cauchy.angle = selected$combined$angle,
      p.values = p.values,
      log.p.values = log.p.values,
      log.one.minus.p.values = c(
        sum = selected$sum$log.one.minus.p.value,
        max = selected$max$log.one.minus.p.value,
        combined = selected$combined$log.one.minus.p.value
      ),
      analytic = list(
        Z.sum = analytic$sum.argument,
        Gumbel.argument = analytic$max.argument,
        p.values = c(
          sum = analytic$sum$p.value,
          max = analytic$max$p.value,
          combined = analytic$combined$p.value
        )
      ),
      calibration = selected$details,
      bootstrap = bootstrap.summary,
      log.radii = log.radii,
      directions = directions,
      fitted.location = stats::setNames(
        as.numeric(standardisation$location), variable.names
      ),
      fitted.shape = standardisation$shape,
      n = n,
      p = p,
      alpha = alpha,
      rejected = main.p <= alpha
    ),
    diagnostics = list(
      standardisation.source = standardisation$source,
      standardisation.fit = standardisation$fit,
      supplied.shape.relative.asymmetry =
        standardisation$relative.asymmetry,
      shape.eigenvalues.scaled = as.numeric(
        core$shape_eigenvalues_scaled
      ),
      shape.scale = as.numeric(core$shape_scale),
      data.scale = as.numeric(core$data_scale),
      minimum.log.radius = as.numeric(core$minimum_log_radius),
      maximum.log.radius = as.numeric(core$maximum_log_radius),
      radius.zero.tolerance.relative = radius.zero.tol,
      log.radius.variance.sum =
        as.numeric(core$log_radius_variance_sum),
      direction.variance.sums = stats::setNames(
        as.numeric(core$direction_variance_sums), variable.names
      ),
      subtraction.overflow.fallbacks =
        as.integer(core$subtraction_fallbacks),
      correlation.roundoff.clips =
        as.integer(core$correlation_roundoff_clips),
      selected.component = component,
      selected.calibration = calibration,
      intrinsic.calibration = calibration == "radial_bootstrap",
      simulation.replication = FALSE,
      method.target = paste(
        "coordinatewise correlation between fitted log radius and fitted",
        "direction; a necessary elliptical implication, not a finite-sample",
        "omnibus characterization against every alternative"
      ),
      invariance = paste(
        "exact finite-sample: translation, global scale, and signed",
        "coordinate permutations; coordinatewise empirical variance",
        "normalisation makes both components basis-sensitive under a general",
        "rotation, while the null theory is asymptotically rotation-compatible;",
        "banded HR is additionally coordinate-order-sensitive"
      ),
      primary.numerical.boundary = paste(
        "the paper's numerical ridge and positive-definite projection are",
        "not used because the projection floor is not uniquely specified"
      ),
      book.completeness = paste(
        "the book gives the analytic statistics but omits the primary",
        "Section 2.4 radial--directional bootstrap calibration"
      ),
      no.repair = paste(
        "No ridge, eigenvalue floor, positive-definite projection,",
        "pseudoinverse, zero-radius perturbation, degenerate-resample",
        "deletion, redraw, or probability clipping was applied."
      ),
      call = call
    ),
    call = call
  )
  class(result) <- c("radial_directional_test", "htest")
  result
}
