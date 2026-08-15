.ch4_afc_matrix <- function(x, name, min_rows = 2L, min_cols = 1L) {
  if (is.data.frame(x)) {
    if (!all(vapply(x, is.numeric, logical(1)))) {
      stop(sprintf("%s must contain only numeric columns.", name),
           call. = FALSE)
    }
    x <- as.matrix(x)
  }
  if (!is.matrix(x) || !is.numeric(x) || length(dim(x)) != 2L) {
    stop(sprintf("%s must be a numeric matrix or numeric data frame.", name),
         call. = FALSE)
  }
  storage.mode(x) <- "double"
  if (nrow(x) < min_rows || ncol(x) < min_cols) {
    stop(sprintf(
      "%s must have at least %d rows and %d columns.",
      name, min_rows, min_cols
    ), call. = FALSE)
  }
  if (anyNA(x) || any(!is.finite(x))) {
    stop(sprintf("%s must contain only finite values.", name),
         call. = FALSE)
  }
  x
}

.ch4_afc_factors <- function(factors, n) {
  if (is.null(factors)) {
    return(matrix(numeric(), nrow = n, ncol = 0L))
  }
  if (is.numeric(factors) && is.vector(factors) && is.null(dim(factors))) {
    factors <- matrix(factors, ncol = 1L)
  }
  factors <- .ch4_afc_matrix(factors, "factors", n, 1L)
  if (nrow(factors) != n) {
    stop("returns and factors must have the same number of rows.",
         call. = FALSE)
  }
  factors
}

.ch4_afc_probability <- function(x, name) {
  x <- as.numeric(x)
  if (length(x) != 1L || is.na(x) || !is.finite(x) || x <= 0 || x >= 1) {
    stop(sprintf("%s must be one finite number strictly between zero and one.",
                 name), call. = FALSE)
  }
  x
}

.ch4_afc_integer <- function(x, name, minimum = 0L) {
  raw <- x
  x <- as.numeric(x)
  if (length(x) != 1L || is.na(x) || !is.finite(x) ||
      x != floor(x) || x < minimum || x > .Machine$integer.max) {
    stop(sprintf("%s must be one integer at least %d.", name, minimum),
         call. = FALSE)
  }
  as.integer(raw)
}

.ch4_afc_controls <- function(tol, max_iter, zero_tol) {
  tol <- as.numeric(tol)
  zero_tol <- as.numeric(zero_tol)
  max_iter <- .ch4_afc_integer(max_iter, "max_iter", 1L)
  if (length(tol) != 1L || is.na(tol) || !is.finite(tol) || tol <= 0) {
    stop("tol must be one finite positive number.", call. = FALSE)
  }
  if (length(zero_tol) != 1L || is.na(zero_tol) ||
      !is.finite(zero_tol) || zero_tol < 0) {
    stop("zero_tol must be one finite non-negative number.", call. = FALSE)
  }
  list(tol = tol, max_iter = max_iter, zero_tol = zero_tol)
}

.ch4_afc_common_scale <- function(x, name) {
  magnitude <- max(abs(x))
  if (!is.finite(magnitude) || magnitude <= 0) {
    stop(sprintf("%s must have positive finite magnitude.", name),
         call. = FALSE)
  }
  list(data = x / magnitude, scale = magnitude)
}

.ch4_afc_column_scale <- function(x, name) {
  if (!ncol(x)) {
    return(list(data = x, scale = numeric()))
  }
  magnitude <- apply(abs(x), 2L, max)
  if (any(!is.finite(magnitude)) || any(magnitude <= 0)) {
    stop(sprintf("Every column of %s must have positive finite magnitude.",
                 name), call. = FALSE)
  }
  list(data = sweep(x, 2L, magnitude, "/"), scale = as.numeric(magnitude))
}

.ch4_afc_flag <- function(x, name) {
  if (!is.logical(x) || length(x) != 1L || is.na(x)) {
    stop(sprintf("%s must be TRUE or FALSE.", name), call. = FALSE)
  }
  x
}

.ch4_afc_asset_names <- function(x) {
  if (is.null(colnames(x))) paste0("asset", seq_len(ncol(x))) else colnames(x)
}

.ch4_afc_observation_names <- function(x) {
  if (is.null(rownames(x))) {
    paste0("observation", seq_len(nrow(x)))
  } else {
    rownames(x)
  }
}

.ch4_afc_check_dimnames <- function(x, row_names = NULL, col_names = NULL,
                                    name) {
  if (!is.null(row_names) && !is.null(rownames(x)) &&
      !identical(as.character(rownames(x)), as.character(row_names))) {
    stop(sprintf("The row names of %s do not match the observation order.",
                 name), call. = FALSE)
  }
  if (!is.null(col_names) && !is.null(colnames(x)) &&
      !identical(as.character(colnames(x)), as.character(col_names))) {
    stop(sprintf("The column names of %s do not match the asset order.",
                 name), call. = FALSE)
  }
  invisible(x)
}

.ch4_afc_gumbel_tail <- function(centered) {
  centered <- as.numeric(centered)
  if (length(centered) != 1L || is.na(centered) || !is.finite(centered)) {
    stop("The centered Gumbel statistic must be finite.", call. = FALSE)
  }
  log_intensity <- -0.5 * log(pi) - 0.5 * centered
  if (log_intensity > log(.Machine$double.xmax)) {
    return(list(p.value = 1, log.p.value = 0,
                log.intensity = log_intensity))
  }
  if (log_intensity < log(.Machine$double.xmin)) {
    return(list(p.value = exp(log_intensity),
                log.p.value = log_intensity,
                log.intensity = log_intensity))
  }
  intensity <- exp(log_intensity)
  log_p <- if (intensity < 0.5) {
    log(-expm1(-intensity))
  } else {
    log1p(-exp(-intensity))
  }
  list(p.value = exp(log_p), log.p.value = log_p,
       log.intensity = log_intensity)
}

.ch4_afc_truncated_cauchy <- function(p.values) {
  p.values <- as.numeric(p.values)
  if (!length(p.values) || anyNA(p.values) || any(!is.finite(p.values)) ||
      any(p.values < 0) || any(p.values > 1)) {
    stop("Truncated-Cauchy component p-values must lie in [0,1].",
         call. = FALSE)
  }
  active <- p.values < 0.5
  terms <- numeric(length(p.values))
  terms[active] <- 0.5 / tan(pi * p.values[active])
  score <- sum(terms)
  if (is.nan(score) || score < 0) {
    stop("The truncated-Cauchy score is invalid.", call. = FALSE)
  }
  list(statistic = score, p.value = atan2(1, score) / pi,
       terms = terms, active = active)
}

.ch4_afc_htest <- function(statistic, p.value, method, data.name,
                            estimate, raw.statistic = NULL,
                            components = list(), diagnostics = list(),
                            call = NULL) {
  answer <- list(
    statistic = statistic,
    p.value = as.numeric(p.value),
    method = method,
    data.name = data.name,
    alternative = "two.sided",
    estimate = estimate,
    null.value = stats::setNames(0, "average conditional alpha"),
    raw.statistic = raw.statistic,
    components = components,
    diagnostics = diagnostics,
    call = call
  )
  answer <- answer[!vapply(answer, is.null, logical(1))]
  class(answer) <- c("conditional_alpha_test", "hd_alpha_test", "htest")
  answer
}

.ch4_afc_validate_design <- function(design, n) {
  design <- .ch4_afc_matrix(design, "design", n, 1L)
  if (nrow(design) != n) {
    stop("returns and design must have the same number of rows.",
         call. = FALSE)
  }
  if (ncol(design) >= n) {
    stop("design must have fewer columns than observations.",
         call. = FALSE)
  }
  scaled <- .ch4_afc_column_scale(design, "design")
  root <- tryCatch(chol(crossprod(scaled$data)), error = function(e) NULL)
  if (is.null(root)) {
    stop(paste0(
      "design must have full column rank; no generalized inverse, ridge, ",
      "or automatic spline-column deletion is used."
    ), call. = FALSE)
  }
  list(data = design, scaled = scaled$data, scale = scaled$scale)
}


#' Construct a conditional-alpha sieve design from a supplied basis
#'
#' This helper performs only the algebra specified by the conditional-alpha
#' papers. Rows are time observations. If `center_alpha = TRUE`, its first
#' block is the column-centered supplied basis; otherwise it is the original
#' basis. The remaining blocks are, in factor order, the rowwise products
#' \eqn{f_{jt}B(t/T)}. It does not choose spline order, knots, basis dimension,
#' or a BIC rule.
#'
#' A normalized B-spline basis usually contains the constant function, so all
#' of its centered columns are linearly dependent. The papers write an
#' ordinary inverse despite this identity. This function never silently drops
#' a column: pass an explicit full-column-rank `alpha_contrast`, or supply a
#' reduced basis. The contrast changes only the coordinates, not the spanned
#' centered-alpha space, when it has the intended range.
#'
#' @param basis Finite observation-by-basis numeric matrix.
#' @param factors `NULL`, a finite numeric vector, or an
#'   observation-by-factor numeric matrix.
#' @param center_alpha Whether to center the alpha-basis block. Use `TRUE` for
#'   the restricted null fit and `FALSE` for the unrestricted residual fit
#'   required by the Zhao CSS trace estimator.
#' @param alpha_contrast Optional finite matrix with `ncol(basis)` rows. The
#'   centered or uncentered alpha block is post-multiplied by this matrix.
#' @return A full-column-rank design matrix carrying its basis/factor metadata.
#' @references
#' Ma, S., Lan, W., Su, L. and Tsai, C.-L. (2020). Testing alphas in
#' conditional time-varying factor models with high-dimensional assets.
#' *Journal of Business & Economic Statistics*, 38, 214--227.
#' \doi{10.1080/07350015.2018.1482758}.
#' @examples
#' tt <- seq(0, 1, length.out = 12)
#' basis <- cbind(tt, tt^2)
#' factors <- cbind(market = sin(seq_len(12)))
#' conditional_alpha_sieve_design(basis, factors)
#' @export
conditional_alpha_sieve_design <- function(
    basis, factors = NULL, center_alpha = TRUE, alpha_contrast = NULL) {
  basis <- .ch4_afc_matrix(basis, "basis", 3L, 1L)
  factors <- .ch4_afc_factors(factors, nrow(basis))
  .ch4_afc_check_dimnames(factors, row_names = rownames(basis),
                          name = "factors")
  center_alpha <- .ch4_afc_flag(center_alpha, "center_alpha")
  alpha_block <- if (center_alpha) {
    sweep(basis, 2L, colMeans(basis), "-")
  } else {
    basis
  }
  if (!is.null(alpha_contrast)) {
    alpha_contrast <- .ch4_afc_matrix(
      alpha_contrast, "alpha_contrast", ncol(basis), 1L
    )
    if (nrow(alpha_contrast) != ncol(basis)) {
      stop("alpha_contrast must have one row per basis column.",
           call. = FALSE)
    }
    contrast_root <- tryCatch(
      chol(crossprod(alpha_contrast)), error = function(e) NULL
    )
    if (is.null(contrast_root)) {
      stop("alpha_contrast must have full column rank.", call. = FALSE)
    }
    alpha_block <- alpha_block %*% alpha_contrast
  }
  factor_blocks <- lapply(seq_len(ncol(factors)), function(j) {
    basis * factors[, j]
  })
  design <- do.call(cbind, c(list(alpha_block), factor_blocks))
  colnames(design) <- c(
    paste0(if (center_alpha) "alpha.centered." else "alpha.raw.",
           seq_len(ncol(alpha_block))),
    unlist(lapply(seq_len(ncol(factors)), function(j) {
      factor_name <- if (is.null(colnames(factors))) {
        paste0("factor", j)
      } else {
        colnames(factors)[j]
      }
      paste0(factor_name, ".basis.", seq_len(ncol(basis)))
    }), use.names = FALSE)
  )
  rownames(design) <- rownames(basis)
  .ch4_afc_validate_design(design, nrow(basis))
  attr(design, "basis.columns") <- ncol(basis)
  attr(design, "alpha.columns") <- ncol(alpha_block)
  attr(design, "factor.count") <- ncol(factors)
  attr(design, "center.alpha") <- center_alpha
  attr(design, "construction") <- "supplied basis; no knot/BIC selection"
  class(design) <- c("conditional_alpha_sieve_design", "matrix", "array")
  design
}


#' Fit the restricted conditional-alpha sieve regression
#'
#' Fits, independently for every asset, the no-intercept null regression on a
#' supplied full-rank sieve design. The residualized intercept
#' \eqn{h=M_Z1_T} is retained because it enters every feasible calibration.
#' The function uses one common response scale and independent design-column
#' scales for numerical stability; both leave the fitted residual subspace and
#' every reported test unchanged. It never uses a generalized inverse or
#' ridge.
#'
#' @param returns Finite observation-by-asset numeric matrix or data frame,
#'   with at least two assets.
#' @param design A finite full-column-rank observation-by-regressor matrix,
#'   normally produced by [conditional_alpha_sieve_design()].
#' @return A reusable `conditional_alpha_sieve_fit` object containing the
#'   restricted residuals, `h`, design metadata, and projection diagnostics.
#' @references
#' Ma, S., Lan, W., Su, L. and Tsai, C.-L. (2020). Testing alphas in
#' conditional time-varying factor models with high-dimensional assets.
#' *Journal of Business & Economic Statistics*, 38, 214--227.
#' \doi{10.1080/07350015.2018.1482758}.
#' @examples
#' tt <- seq(0, 1, length.out = 14)
#' z <- conditional_alpha_sieve_design(cbind(tt, tt^2))
#' y <- cbind(sin(1:14), cos(1:14), sin(1:14 / 2))
#' conditional_alpha_sieve_fit(y, z)
#' @export
conditional_alpha_sieve_fit <- function(returns, design) {
  call <- match.call()
  data.name <- deparse1(substitute(returns))
  returns <- .ch4_afc_matrix(returns, "returns", 3L, 2L)
  .ch4_afc_check_dimnames(design, row_names = rownames(returns),
                          name = "design")
  design_meta <- list(
    basis.columns = attr(design, "basis.columns", exact = TRUE),
    alpha.columns = attr(design, "alpha.columns", exact = TRUE),
    factor.count = attr(design, "factor.count", exact = TRUE),
    center.alpha = attr(design, "center.alpha", exact = TRUE),
    construction = attr(design, "construction", exact = TRUE)
  )
  design_checked <- .ch4_afc_validate_design(design, nrow(returns))
  response <- .ch4_afc_common_scale(returns, "returns")
  projected <- cpp_ch4_afc_project(response$data, design_checked$scaled)
  residuals_scaled <- as.matrix(projected$residuals)
  dimnames(residuals_scaled) <- list(
    .ch4_afc_observation_names(returns), .ch4_afc_asset_names(returns)
  )
  residuals <- residuals_scaled * response$scale
  residuals_representable <- all(is.finite(residuals))
  if (!residuals_representable) residuals <- NULL
  h <- as.numeric(projected$h)
  names(h) <- .ch4_afc_observation_names(returns)
  answer <- list(
    residuals = residuals,
    residuals.scaled = residuals_scaled,
    response.scale = response$scale,
    h = h,
    h2 = unname(projected$h2),
    coefficients.scaled = as.matrix(projected$coefficients),
    design = design_checked$data,
    design.column.scales = design_checked$scale,
    design.columns = ncol(design_checked$data),
    factor.count = design_meta$factor.count,
    basis.columns = design_meta$basis.columns,
    alpha.columns = design_meta$alpha.columns,
    center.alpha = design_meta$center.alpha,
    T = nrow(returns),
    N = ncol(returns),
    asset.names = .ch4_afc_asset_names(returns),
    observation.names = .ch4_afc_observation_names(returns),
    diagnostics = list(
      reciprocal.condition = unname(projected$reciprocal_condition),
      normal.equation.residual =
        unname(projected$normal_equation_residual),
      residuals.input.scale.representable = residuals_representable,
      construction = design_meta$construction,
      generalized.inverse = "none",
      ridge = "none",
      automatic.column.deletion = "none"
    ),
    data.name = data.name,
    call = call
  )
  class(answer) <- c("conditional_alpha_sieve_fit", "list")
  answer
}

.ch4_afc_fit <- function(fit, name = "fit") {
  if (!inherits(fit, "conditional_alpha_sieve_fit") ||
      !is.matrix(fit$residuals.scaled) || !is.numeric(fit$h) ||
      nrow(fit$residuals.scaled) != length(fit$h) ||
      ncol(fit$residuals.scaled) < 2L || anyNA(fit$residuals.scaled) ||
      any(!is.finite(fit$residuals.scaled)) || anyNA(fit$h) ||
      any(!is.finite(fit$h))) {
    stop(sprintf("%s must be a valid conditional_alpha_sieve_fit object.",
                 name), call. = FALSE)
  }
  T <- nrow(fit$residuals.scaled)
  N <- ncol(fit$residuals.scaled)
  h2 <- sum(fit$h^2)
  metadata_valid <-
    is.numeric(fit$T) && length(fit$T) == 1L && is.finite(fit$T) &&
    fit$T == T && is.numeric(fit$N) && length(fit$N) == 1L &&
    is.finite(fit$N) && fit$N == N && is.numeric(fit$h2) &&
    length(fit$h2) == 1L && is.finite(fit$h2) && fit$h2 > 0 &&
    is.finite(h2) && h2 > 0 &&
    isTRUE(all.equal(as.numeric(fit$h2), h2, tolerance = 1e-12)) &&
    is.numeric(fit$response.scale) && length(fit$response.scale) == 1L &&
    is.finite(fit$response.scale) && fit$response.scale > 0 &&
    is.numeric(fit$design.columns) && length(fit$design.columns) == 1L &&
    is.finite(fit$design.columns) &&
    fit$design.columns == floor(fit$design.columns) &&
    fit$design.columns >= 1 && fit$design.columns < T &&
    is.character(fit$asset.names) && length(fit$asset.names) == N &&
    !anyNA(fit$asset.names) &&
    is.character(fit$observation.names) &&
    length(fit$observation.names) == T && !anyNA(fit$observation.names) &&
    identical(as.character(colnames(fit$residuals.scaled)),
              as.character(fit$asset.names)) &&
    identical(as.character(rownames(fit$residuals.scaled)),
              as.character(fit$observation.names)) &&
    identical(as.character(names(fit$h)),
              as.character(fit$observation.names))
  if (!metadata_valid) {
    stop(sprintf("%s must be a valid conditional_alpha_sieve_fit object.",
                 name), call. = FALSE)
  }
  fit
}

.ch4_afc_omega_ratio <- function(x) {
  x <- as.numeric(x)
  if (length(x) != 1L || is.na(x) || !is.finite(x) || x <= 0 || x > 1) {
    stop("omega_ratio must be one finite number in (0, 1].",
         call. = FALSE)
  }
  x
}

.ch4_afc_fdr_base <- function(returns, factors, restricted_residuals,
                               omega_ratio) {
  if (is.null(restricted_residuals)) {
    if (is.null(returns)) {
      stop("Supply returns or restricted_residuals.", call. = FALSE)
    }
    returns <- .ch4_afc_matrix(returns, "returns", 3L, 2L)
    factors <- .ch4_afc_factors(factors, nrow(returns))
    .ch4_afc_check_dimnames(factors, row_names = rownames(returns),
                            name = "factors")
    response <- .ch4_afc_common_scale(returns, "returns")
    factor_scaled <- .ch4_afc_column_scale(factors, "factors")
    projected <- cpp_ch4_afc_project(response$data, factor_scaled$data)
    computed_omega <- unname(projected$h2) / nrow(returns)
    if (!is.null(omega_ratio) &&
        !isTRUE(all.equal(.ch4_afc_omega_ratio(omega_ratio),
                         computed_omega, tolerance = 1e-10))) {
      stop("omega_ratio conflicts with the ratio implied by factors.",
           call. = FALSE)
    }
    return(list(
      residuals.scaled = as.matrix(projected$residuals),
      unit.scale = response$scale,
      omega.ratio = computed_omega,
      factors.scaled = factor_scaled$data,
      source = "returns projected on observed factors without an intercept",
      asset.names = .ch4_afc_asset_names(returns),
      observation.names = .ch4_afc_observation_names(returns),
      projection = projected,
      returns.supplied = TRUE
    ))
  }
  if (!is.null(returns)) {
    stop("Supply either returns or restricted_residuals, not both.",
         call. = FALSE)
  }
  residuals <- .ch4_afc_matrix(
    restricted_residuals, "restricted_residuals", 3L, 2L
  )
  residual_scale <- .ch4_afc_common_scale(
    residuals, "restricted_residuals"
  )
  factors <- .ch4_afc_factors(factors, nrow(residuals))
  .ch4_afc_check_dimnames(factors, row_names = rownames(residuals),
                          name = "factors")
  factor_scaled <- .ch4_afc_column_scale(factors, "factors")
  if (is.null(omega_ratio)) {
    if (!ncol(factors)) {
      stop(paste0(
        "When restricted_residuals are supplied without factors, ",
        "omega_ratio must also be supplied; it cannot be inferred."
      ), call. = FALSE)
    }
    dummy <- matrix(0, nrow(residuals), 1L)
    projected <- cpp_ch4_afc_project(dummy, factor_scaled$data)
    omega_ratio <- unname(projected$h2) / nrow(residuals)
  } else {
    omega_ratio <- .ch4_afc_omega_ratio(omega_ratio)
    if (ncol(factors)) {
      dummy <- matrix(0, nrow(residuals), 1L)
      projected <- cpp_ch4_afc_project(dummy, factor_scaled$data)
      computed_omega <- unname(projected$h2) / nrow(residuals)
      if (!isTRUE(all.equal(omega_ratio, computed_omega,
                            tolerance = 1e-10))) {
        stop("omega_ratio conflicts with the ratio implied by factors.",
             call. = FALSE)
      }
    } else {
      projected <- NULL
    }
  }
  list(
    residuals.scaled = residual_scale$data,
    unit.scale = residual_scale$scale,
    omega.ratio = omega_ratio,
    factors.scaled = factor_scaled$data,
    source = "user-supplied restricted residuals",
    asset.names = .ch4_afc_asset_names(residuals),
    observation.names = .ch4_afc_observation_names(residuals),
    projection = projected,
    returns.supplied = FALSE
  )
}

.ch4_afc_latent_adjust <- function(base, adjustment, latent_fit,
                                    n_factors, k_max) {
  z <- base$residuals.scaled
  T <- nrow(z)
  N <- ncol(z)
  diagnostics <- list(method = adjustment, n.factors = 0L)
  if (adjustment == "none") {
    if (!is.null(latent_fit) || !is.null(n_factors) || !is.null(k_max)) {
      stop(paste0(
        "latent_fit, n_factors, and k_max are not used when ",
        "adjustment = 'none'."
      ), call. = FALSE)
    }
    return(list(data = z, unit.scale = base$unit.scale,
                diagnostics = diagnostics))
  }

  if (adjustment == "supplied") {
    if (!is.null(n_factors) || !is.null(k_max)) {
      stop(paste0(
        "n_factors and k_max are used only with ",
        "adjustment = 'spatial_kendall'."
      ), call. = FALSE)
    }
    if (!is.list(latent_fit)) {
      stop("latent_fit must be a list for adjustment = 'supplied'.",
           call. = FALSE)
    }
    if (!is.null(latent_fit$adjusted.residuals)) {
      adjusted <- .ch4_afc_matrix(
        latent_fit$adjusted.residuals,
        "latent_fit$adjusted.residuals", T, N
      )
      if (!identical(dim(adjusted), c(T, N))) {
        stop("latent_fit$adjusted.residuals has incompatible dimensions.",
             call. = FALSE)
      }
      .ch4_afc_check_dimnames(
        adjusted, base$observation.names, base$asset.names,
        "latent_fit$adjusted.residuals"
      )
      scaled <- .ch4_afc_common_scale(
        adjusted, "latent_fit$adjusted.residuals"
      )
      diagnostics$source <- "supplied final factor-adjusted residuals"
      diagnostics$n.factors <- if (is.null(latent_fit$n.factors)) {
        NA_integer_
      } else {
        .ch4_afc_integer(latent_fit$n.factors,
                         "latent_fit$n.factors", 0L)
      }
      return(list(data = scaled$data, unit.scale = scaled$scale,
                  diagnostics = diagnostics))
    }
    common <- latent_fit$common.component
    if (is.null(common) && !is.null(latent_fit$loadings) &&
        !is.null(latent_fit$scores)) {
      loadings <- .ch4_afc_matrix(
        latent_fit$loadings, "latent_fit$loadings", N, 1L
      )
      scores <- .ch4_afc_matrix(
        latent_fit$scores, "latent_fit$scores", T, ncol(loadings)
      )
      if (nrow(loadings) != N || nrow(scores) != T ||
          ncol(scores) != ncol(loadings)) {
        stop("The supplied latent loadings and scores have incompatible dimensions.",
             call. = FALSE)
      }
      .ch4_afc_check_dimnames(
        loadings, row_names = base$asset.names,
        name = "latent_fit$loadings"
      )
      .ch4_afc_check_dimnames(
        scores, row_names = base$observation.names,
        name = "latent_fit$scores"
      )
      if (!is.null(colnames(loadings)) && !is.null(colnames(scores)) &&
          !identical(colnames(loadings), colnames(scores))) {
        stop("The supplied latent loadings and scores use different factor orders.",
             call. = FALSE)
      }
      if (!is.null(latent_fit$n.factors) &&
          .ch4_afc_integer(latent_fit$n.factors,
                           "latent_fit$n.factors", 0L) != ncol(loadings)) {
        stop("latent_fit$n.factors conflicts with the supplied loadings.",
             call. = FALSE)
      }
      common <- scores %*% t(loadings)
      diagnostics$n.factors <- ncol(loadings)
      diagnostics$source <- "supplied loadings and scores"
    } else if (!is.null(common)) {
      common <- .ch4_afc_matrix(
        common, "latent_fit$common.component", T, N
      )
      if (!identical(dim(common), c(T, N))) {
        stop("latent_fit$common.component has incompatible dimensions.",
             call. = FALSE)
      }
      .ch4_afc_check_dimnames(
        common, base$observation.names, base$asset.names,
        "latent_fit$common.component"
      )
      diagnostics$n.factors <- if (is.null(latent_fit$n.factors)) {
        NA_integer_
      } else {
        .ch4_afc_integer(latent_fit$n.factors,
                         "latent_fit$n.factors", 0L)
      }
      diagnostics$source <- "supplied latent common component"
    } else {
      stop(paste0(
        "latent_fit must contain adjusted.residuals, common.component, ",
        "or compatible loadings and scores."
      ), call. = FALSE)
    }
    common.scaled <- common / base$unit.scale
    if (base$returns.supplied && ncol(base$factors.scaled)) {
      common.projected <- cpp_ch4_afc_project(
        common.scaled, base$factors.scaled
      )$residuals
    } else {
      common.projected <- common.scaled
      diagnostics$common.component.contract <- paste0(
        "With supplied restricted residuals, the supplied common component ",
        "must already be in the same restricted residual space."
      )
    }
    adjusted <- z - common.projected
    scaled <- .ch4_afc_common_scale(adjusted, "factor-adjusted residuals")
    return(list(data = scaled$data,
                unit.scale = base$unit.scale * scaled$scale,
                diagnostics = diagnostics))
  }

  if (!is.null(latent_fit)) {
    stop("latent_fit is used only with adjustment = 'supplied'.",
         call. = FALSE)
  }
  kendall <- cpp_ch4_afc_spatial_kendall(z)
  eig <- eigen(as.matrix(kendall$matrix), symmetric = TRUE)
  values <- as.numeric(eig$values)
  maximum_factors <- min(N - 1L, length(values) - 1L)
  selection <- "fixed"
  ratios <- NULL
  if (is.null(n_factors)) {
    if (is.null(k_max)) {
      stop(paste0(
        "For adjustment = 'spatial_kendall', supply n_factors or an ",
        "explicit k_max for the paper's eigenvalue-ratio rule."
      ), call. = FALSE)
    }
    k_max <- .ch4_afc_integer(k_max, "k_max", 1L)
    if (k_max > maximum_factors) {
      stop("k_max is too large for the spatial Kendall spectrum.",
           call. = FALSE)
    }
    denominators <- values[seq.int(2L, k_max + 1L)]
    if (any(!is.finite(denominators)) || any(denominators <= 0)) {
      stop(paste0(
        "The requested eigenvalue-ratio range has a non-positive ",
        "denominator; the factor count is not uniquely feasible."
      ), call. = FALSE)
    }
    ratios <- values[seq_len(k_max)] / denominators
    n_factors <- which.max(ratios)
    selection <- "paper eigenvalue-ratio rule with explicit k_max"
  } else {
    n_factors <- .ch4_afc_integer(n_factors, "n_factors", 1L)
    if (!is.null(k_max)) {
      stop("Supply n_factors or k_max, not both.", call. = FALSE)
    }
    if (n_factors > maximum_factors) {
      stop("n_factors is too large for the spatial Kendall spectrum.",
           call. = FALSE)
    }
  }
  loadings <- sqrt(N) * eig$vectors[, seq_len(n_factors), drop = FALSE]
  scores <- z %*% loadings / N
  common <- scores %*% t(loadings)
  adjusted <- z - common
  scaled <- .ch4_afc_common_scale(adjusted, "factor-adjusted residuals")
  diagnostics <- list(
    method = "spatial_kendall",
    source = "primary spatial Kendall elliptical PCA and least-squares scores",
    n.factors = n_factors,
    selection = selection,
    k.max = if (selection == "fixed") NULL else k_max,
    eigenvalues = values,
    eigenvalue.ratios = ratios,
    zero.pair.differences = as.integer(kendall$zero_pairs),
    spatial.kendall.trace = unname(kendall$trace),
    spatial.kendall = as.matrix(kendall$matrix),
    loadings = loadings,
    scores = scores,
    common.component.scaled = common
  )
  list(data = scaled$data,
       unit.scale = base$unit.scale * scaled$scale,
       diagnostics = diagnostics)
}

.ch4_afc_ssbh <- function(x, unit_scale, omega_ratio, q, controls,
                           asset_names, observation_names) {
  q <- .ch4_afc_probability(q, "q")
  omega_ratio <- .ch4_afc_omega_ratio(omega_ratio)
  fit <- cpp_scaled_spatial_median(
    x, numeric(ncol(x)), TRUE, controls$tol, controls$max_iter,
    controls$zero_tol
  )
  if (!isTRUE(fit$iteration_stable)) {
    stop(sprintf(
      paste0("The SS-BH scaled spatial median did not stabilize within %d ",
             "updates (score residual %.6g, relative update %.6g)."),
      as.integer(fit$max_iterations), as.numeric(fit$score_residual),
      as.numeric(fit$relative_update)
    ), call. = FALSE)
  }
  radii <- as.numeric(fit$radii)
  if (anyNA(radii) || any(!is.finite(radii)) || any(radii <= 0)) {
    stop(paste0(
      "Every fitted standardized radius must be strictly positive and ",
      "finite; no deletion, cap, or floor is used."
    ), call. = FALSE)
  }
  radial <- c(
    mean.inverse.r = mean(1 / radii),
    mean.r = mean(radii),
    mean.r.squared = mean(radii^2)
  )
  if (any(!is.finite(radial))) {
    stop("A required SS-BH radial moment is not finite.", call. = FALSE)
  }
  projection_loss <- 1 - omega_ratio
  denominator <- unname(
    1 - 2 * projection_loss * radial["mean.inverse.r"] *
      radial["mean.r"] + projection_loss * radial["mean.r.squared"] *
      radial["mean.inverse.r"]^2
  )
  if (!is.finite(denominator) || denominator <= 0) {
    stop(paste0(
      "The primary SS-BH varsigma denominator is not strictly positive ",
      "and finite; no absolute value or floor is used."
    ), call. = FALSE)
  }
  varsigma <- ncol(x) * unname(radial["mean.inverse.r"])^2 /
    denominator
  coordinates <- as.numeric(fit$standardized_location_coordinate)
  statistics <- sqrt(nrow(x) * varsigma) * coordinates
  if (any(!is.finite(statistics))) {
    stop("An SS-BH coordinate statistic is not finite.", call. = FALSE)
  }
  p.values <- stats::pnorm(statistics, lower.tail = FALSE)
  ordering <- order(p.values, seq_along(p.values))
  passes <- p.values[ordering] <= q * seq_along(p.values) / length(p.values)
  k <- if (any(passes)) max(which(passes)) else 0L
  threshold <- if (k) p.values[ordering[k]] else NA_real_
  rejected <- if (k) which(p.values <= threshold) else integer()
  adjusted <- stats::p.adjust(p.values, method = "BH")
  theta <- as.numeric(fit$location_relative_origin) * unit_scale
  log_diagonal <- as.numeric(fit$log_diagonal_input_canonical) +
    2 * log(unit_scale)
  log_diagonal <- log_diagonal - max(log_diagonal)
  diagonal <- exp(log_diagonal)
  directions <- as.matrix(fit$directions)
  dimnames(directions) <- list(observation_names, asset_names)
  names(statistics) <- names(p.values) <- names(adjusted) <-
    names(theta) <- names(coordinates) <- names(diagonal) <-
    names(log_diagonal) <- asset_names
  names(rejected) <- asset_names[rejected]
  list(
    statistic = statistics,
    p.value = p.values,
    p.adjusted = adjusted,
    rejected = rejected,
    rejected.names = asset_names[rejected],
    q = q,
    bh.k = k,
    bh.threshold = threshold,
    theta = theta,
    standardized.theta = coordinates,
    scale.diagonal = diagonal,
    log.scale.diagonal = log_diagonal,
    directions = directions,
    radii = radii,
    radial.moments = radial,
    varsigma = varsigma,
    varsigma.denominator = denominator,
    omega.ratio = omega_ratio,
    fit = fit
  )
}


#' Robust mutual-fund selection with SS-BH or FSS-BH
#'
#' Implements the one-sided procedures of Wang, Zhao, Feng and Wang for
#' \eqn{H_{0i}:\alpha_i\leq 0} against positive fund alpha. With observable
#' factors, the primary observations are the *restricted* residualized returns
#' \eqn{Z=M_FY}; the projection contains no intercept. A simultaneous scaled
#' spatial median gives \eqn{\widehat\theta} and
#' \eqn{\widehat D^{1/2}=\operatorname{diag}(\widehat d_i)}. If
#' \eqn{\widehat r_t=\|\widehat D^{-1/2}
#' (Z_t-\widehat\theta)\|}, then
#' \deqn{\widehat\varsigma=
#' \frac{N\overline{r^{-1}}^2}
#' {1-2(1-\omega_T/T)\overline{r^{-1}}\bar r+
#' (1-\omega_T/T)\overline{r^2}\,\overline{r^{-1}}^2}}
#' and \eqn{T_i^s=\sqrt{T\widehat\varsigma}
#' \widehat\theta_i/\widehat d_i}. One-sided normal p-values are passed to
#' the ordinary BH step-up rule.
#'
#' The book draft instead starts from unrestricted, intercept-containing OLS
#' residuals and replaces \eqn{\sqrt{\widehat\varsigma}} by the inverse-radius
#' mean. Both changes contradict the primary definition and erase the alpha
#' signal. This implementation follows the article and returns the full radial
#' correction. The published corrigendum changes only an affiliation.
#'
#' For FSS-BH, `adjustment = "spatial_kendall"` implements the primary spatial
#' Kendall matrix, \eqn{\widehat\Gamma=\sqrt N(\widehat\xi_1,\ldots,
#' \widehat\xi_r)}, least-squares scores, and factor removal. Supply a fixed
#' `n_factors`, or explicitly opt into the paper's eigenvalue-ratio selector by
#' supplying `k_max`; there is deliberately no tuning default. Alternatively,
#' `adjustment = "supplied"` accepts a validated nuisance fit. No simulation
#' tuning is embedded.
#'
#' @param returns `NULL` or a finite observation-by-fund matrix. Used when
#'   `restricted_residuals` is `NULL`.
#' @param factors Optional observed factor vector or matrix. The projection is
#'   through the origin, as required by the primary restricted fit.
#' @param restricted_residuals Optional already restricted observation-by-fund
#'   matrix. If factors are absent, `omega_ratio` is then mandatory.
#' @param omega_ratio Optional \eqn{\omega_T/T\in(0,1]}. It is inferred from
#'   factors for raw returns.
#' @param q Target BH FDR level in `(0,1)`.
#' @param adjustment One of `"none"`, `"supplied"`, or
#'   `"spatial_kendall"`.
#' @param latent_fit For supplied adjustment, a list containing either
#'   `adjusted.residuals`, `common.component`, or compatible `loadings` and
#'   `scores`, all in the input return units.
#' @param n_factors Explicit positive latent factor count.
#' @param k_max Explicit upper bound for the paper's eigenvalue-ratio rule;
#'   used only when `n_factors` is `NULL`.
#'   Both arguments are used only with `adjustment = "spatial_kendall"`;
#'   irrelevant tuning arguments are rejected rather than silently ignored.
#' @param tol Positive simultaneous scaled-median equation tolerance.
#' @param max_iter Positive maximum update count.
#' @param zero_tol Non-negative exact-zero radius tolerance. The default zero
#'   preserves every nonzero observation.
#' @return A `mutual_fund_fdr` object with fundwise statistics and p-values,
#'   BH decisions, the complete scaled-median fit, factor diagnostics, and all
#'   radial/projection quantities.
#' @references
#' Wang, H., Zhao, P., Feng, L. and Wang, Z. (2025). Robust mutual fund
#' selection with false discovery rate control. *Journal of Econometrics*,
#' 252, 106121. \doi{10.1016/j.jeconom.2025.106121}.
#' Preprint: \url{https://arxiv.org/abs/2411.14016}.
#'
#' Wang, H., Zhao, P., Feng, L. and Wang, Z. (2026). Corrigendum.
#' \doi{10.1016/j.jeconom.2025.106162}.
#' @examples
#' f <- cbind(market = seq(-1, 1, length.out = 14))
#' y <- cbind(f[, 1] + sin(1:14),
#'            0.15 - 0.3 * f[, 1] + cos(1:14),
#'            -0.1 + 0.2 * f[, 1] + sin(1:14 / 2))
#' wang_zhao_feng_wang_mutual_fund_fdr(y, f, q = 0.1)
#' @export
wang_zhao_feng_wang_mutual_fund_fdr <- function(
    returns = NULL, factors = NULL, restricted_residuals = NULL,
    omega_ratio = NULL, q = 0.05,
    adjustment = c("none", "supplied", "spatial_kendall"),
    latent_fit = NULL, n_factors = NULL, k_max = NULL,
    tol = 1e-7, max_iter = 500L, zero_tol = 0) {
  call <- match.call()
  data.name <- if (is.null(restricted_residuals)) {
    deparse1(substitute(returns))
  } else {
    deparse1(substitute(restricted_residuals))
  }
  adjustment <- match.arg(adjustment)
  controls <- .ch4_afc_controls(tol, max_iter, zero_tol)
  base <- .ch4_afc_fdr_base(
    returns, factors, restricted_residuals, omega_ratio
  )
  latent <- .ch4_afc_latent_adjust(
    base, adjustment, latent_fit, n_factors, k_max
  )
  core <- .ch4_afc_ssbh(
    latent$data, latent$unit.scale, base$omega.ratio, q, controls,
    base$asset.names, base$observation.names
  )
  answer <- c(core, list(
    method = if (adjustment == "none") {
      "Wang-Zhao-Feng-Wang spatial-sign BH mutual-fund selection (SS-BH)"
    } else {
      paste0("Wang-Zhao-Feng-Wang factor-adjusted spatial-sign BH ",
             "mutual-fund selection (FSS-BH)")
    },
    alternative = "greater",
    data.name = data.name,
    factor.adjustment = latent$diagnostics,
    restricted.residuals.scaled = base$residuals.scaled,
    factor.adjusted.residuals.scaled = latent$data,
    residual.unit.scale = latent$unit.scale,
    projection.h = if (is.null(base$projection)) NULL else
      as.numeric(base$projection$h),
    diagnostics = list(
      restricted.residual.source = base$source,
      observable.projection.contains.intercept = FALSE,
      calibration = "one-sided standard normal followed by BH step-up",
      book.unrestricted.residual.conflict = TRUE,
      book.varsigma.omission = TRUE,
      factor.count.tuning.default = "none",
      radial.floor = "none",
      ridge = "none",
      converged = isTRUE(core$fit$iteration_stable),
      iterations = as.integer(core$fit$iterations),
      score.residual = unname(core$fit$score_residual),
      relative.update = unname(core$fit$relative_update)
    ),
    call = call
  ))
  class(answer) <- c("mutual_fund_fdr", "list")
  answer
}

.ch4_afc_average_alpha <- function(fit) {
  estimate <- colSums(fit$residuals.scaled) / fit$h2 * fit$response.scale
  names(estimate) <- fit$asset.names
  estimate
}

.ch4_afc_factor_count <- function(fit, factor_count) {
  stored <- fit$factor.count
  if (is.null(factor_count)) {
    if (is.null(stored)) {
      stop(paste0(
        "factor_count is required because the supplied fit does not carry ",
        "factor metadata."
      ), call. = FALSE)
    }
    return(.ch4_afc_integer(stored, "fit$factor.count", 0L))
  }
  factor_count <- .ch4_afc_integer(factor_count, "factor_count", 0L)
  if (!is.null(stored) && as.integer(stored) != factor_count) {
    stop("factor_count conflicts with the metadata stored in fit.",
         call. = FALSE)
  }
  factor_count
}

.ch4_afc_light_components <- function(fit, factor_count = -1L,
                                       compute_max = FALSE,
                                       require_sum = TRUE) {
  fit <- .ch4_afc_fit(fit)
  core <- cpp_ch4_afc_light_components(
    fit$residuals.scaled, fit$h, as.integer(factor_count),
    as.integer(fit$design.columns), compute_max
  )
  trace_estimate <- unname(core$trace_sigma_squared_estimate)
  variance <- unname(core$variance_estimate)
  if (require_sum && (!is.finite(trace_estimate) || trace_estimate <= 0)) {
    stop(paste0(
      "The primary bias-corrected trace estimate is not strictly positive ",
      "and finite; no absolute value or floor is used."
    ), call. = FALSE)
  }
  if (require_sum && (!is.finite(variance) || variance <= 0)) {
    stop(paste0(
      "The primary conditional sum variance estimate is not strictly ",
      "positive and finite; no absolute value or floor is used."
    ), call. = FALSE)
  }
  if (require_sum) {
    sum_z <- (unname(core$sum_statistic) - unname(core$mean_estimate)) /
      sqrt(variance)
    if (!is.finite(sum_z)) {
      stop("The standardized conditional sum statistic is not finite.",
           call. = FALSE)
    }
    core$sum_z <- sum_z
    core$sum_p <- stats::pnorm(sum_z, lower.tail = FALSE)
    core$sum_log_p <- stats::pnorm(sum_z, lower.tail = FALSE, log.p = TRUE)
  } else {
    core$sum_z <- core$sum_p <- core$sum_log_p <- NA_real_
  }
  if (compute_max) {
    maximum <- unname(core$maximum_statistic)
    centered <- maximum - 2 * log(fit$N) + log(log(fit$N))
    tail <- .ch4_afc_gumbel_tail(centered)
    core$maximum_centered <- centered
    core$maximum_p <- tail$p.value
    core$maximum_log_p <- tail$log.p.value
    core$maximum_log_intensity <- tail$log.intensity
    names(core$coordinate_t_squared) <- fit$asset.names
    names(core$marginal_variance) <- fit$asset.names
  }
  core
}


#' Ma--Lan--Su--Tsai conditional high-dimensional alpha sum test
#'
#' Implements the complete feasible HDA statistic for a reusable restricted
#' sieve fit. With \eqn{h=M_Z1_T} and restricted residuals
#' \eqn{\widehat e_{it}},
#' \deqn{S_{NT}=\frac1{NT}\sum_i
#' (\widehat e_{i\cdot}'1_T)^2,\qquad
#' \widehat\mu_{NT}=\frac1{NT}\sum_{i,t}
#' \widehat e_{it}^2h_t^2.}
#' If \eqn{q=\operatorname{ncol}(Z)} and
#' \eqn{\widehat\Sigma=T^{-1}\sum_t
#' (\widehat e_t-\bar e)(\widehat e_t-\bar e)'}, the implemented primary
#' trace correction is
#' \deqn{\widehat{\operatorname{tr}(\Sigma^2)}=
#' \frac{T^2}{(T+q-1)(T-q)}
#' \left\{\operatorname{tr}(\widehat\Sigma^2)-
#' \frac{\operatorname{tr}^2(\widehat\Sigma)}{T-q}\right\}.}
#' The variance estimate is
#' \deqn{\widehat\sigma_{NT}^2=
#' \frac{2\widehat{\operatorname{tr}(\Sigma^2)}}{N^2T^2}
#' \sum_{t\ne s}h_t^2h_s^2.}
#'
#' The book draft's expression using only oracle
#' \eqn{N^{-1}\operatorname{tr}(\Omega)} and
#' \eqn{2N^{-2}\operatorname{tr}(\Omega^2)} omits `h` and the feasible
#' finite-sample trace correction. This function uses the primary estimator
#' and rejects non-positive estimates; it never takes an absolute value or
#' applies a floor.
#'
#' @param fit A valid object from [conditional_alpha_sieve_fit()].
#' @return A `conditional_alpha_test`/`htest` object with the upper-tail
#'   normal calibration and every feasible centering/variance component.
#' @references
#' Ma, S., Lan, W., Su, L. and Tsai, C.-L. (2020). Testing alphas in
#' conditional time-varying factor models with high-dimensional assets.
#' *Journal of Business & Economic Statistics*, 38, 214--227.
#' \doi{10.1080/07350015.2018.1482758}.
#' @examples
#' tt <- seq(0, 1, length.out = 18)
#' z <- conditional_alpha_sieve_design(cbind(tt, tt^2))
#' y <- cbind(sin(1:18), cos(1:18), sin(1:18 / 2))
#' fit <- conditional_alpha_sieve_fit(y, z)
#' ma_lan_su_tsai_conditional_alpha_sum_test(fit)
#' @export
ma_lan_su_tsai_conditional_alpha_sum_test <- function(fit) {
  call <- match.call()
  fit <- .ch4_afc_fit(fit)
  core <- .ch4_afc_light_components(fit, compute_max = FALSE)
  .ch4_afc_htest(
    statistic = stats::setNames(core$sum_z, "HDA Z"),
    p.value = core$sum_p,
    method = "Ma-Lan-Su-Tsai conditional high-dimensional alpha sum test",
    data.name = fit$data.name,
    estimate = .ch4_afc_average_alpha(fit),
    raw.statistic = c(
      S.NT.scaled = unname(core$sum_statistic),
      mu.hat.scaled = unname(core$mean_estimate),
      variance.hat.scaled = unname(core$variance_estimate)
    ),
    components = list(
      S.NT.scaled = unname(core$sum_statistic),
      mu.hat.scaled = unname(core$mean_estimate),
      variance.hat.scaled = unname(core$variance_estimate),
      trace.Sigma = unname(core$trace_sigma),
      trace.Sigma.squared.raw = unname(core$trace_sigma_squared_raw),
      trace.Sigma.squared.estimate =
        unname(core$trace_sigma_squared_estimate),
      trace.inner = unname(core$trace_inner),
      trace.correction.denominator =
        unname(core$trace_correction_denominator),
      sum.h.t.squared.h.s.squared = unname(core$off_diagonal_h4),
      h = fit$h,
      h2 = fit$h2,
      response.scale = fit$response.scale,
      design.columns = fit$design.columns,
      T = fit$T,
      N = fit$N
    ),
    diagnostics = list(
      calibration = "asymptotic upper-tail standard normal",
      trace.formula = "Ma-Lan-Su-Tsai feasible df-corrected estimator",
      book.oracle.standardization.conflict = TRUE,
      absolute.value.repair = "none",
      variance.floor = "none",
      projection = fit$diagnostics
    ),
    call = call
  )
}


#' Ma--Feng--Wang--Bao conditional maximum and adaptive alpha test
#'
#' For a supplied restricted sieve fit, the primary marginal statistics are
#' \deqn{t_i^2=T^{-1}\widehat\sigma_{ii}^{-1}
#' (\widehat e_{i\cdot}'1_T)^2,\qquad
#' \widehat\sigma_{ij}=\widehat e_{i\cdot}'
#' \widehat e_{j\cdot}/(T-d-1),}
#' where `d` is the observed factor count, not the number of sieve columns.
#' The maximum is centered by \eqn{2\log N-\log\log N} and calibrated with
#' cdf \eqn{F(x)=\exp\{-\pi^{-1/2}\exp(-x/2)\}}.
#'
#' The adaptive test combines this maximum p-value and the complete
#' [ma_lan_su_tsai_conditional_alpha_sum_test()] p-value by the primary Fisher
#' statistic \eqn{-2(\log p_M+\log p_S)} with a chi-squared distribution on
#' four degrees of freedom. The book draft, and a later review paragraph,
#' incorrectly attribute a Cauchy combination to this paper.
#'
#' @param fit A valid object from [conditional_alpha_sieve_fit()].
#' @param factor_count Optional non-negative observed factor count. It is
#'   inferred from a design constructed by
#'   [conditional_alpha_sieve_design()]; otherwise it is required.
#' @param component One of `"max"`, `"sum"`, or `"adaptive"`.
#' @return A `conditional_alpha_test`/`htest` object retaining both component
#'   statistics, log p-values, the primary marginal divisor, and calibration
#'   diagnostics.
#' @references
#' Ma, H., Feng, L., Wang, Z. and Bao, J. (2024). Adaptive testing for alphas
#' in conditional factor models with high dimensional assets. *Journal of
#' Business & Economic Statistics*, 42, 1356--1366.
#' \doi{10.1080/07350015.2024.2313543}.
#' Preprint: \url{https://arxiv.org/abs/2307.09397}.
#' @examples
#' tt <- seq(0, 1, length.out = 18)
#' f <- cbind(market = sin(1:18 / 3))
#' z <- conditional_alpha_sieve_design(cbind(tt, tt^2), f)
#' y <- cbind(sin(1:18), cos(1:18), sin(1:18 / 2))
#' fit <- conditional_alpha_sieve_fit(y, z)
#' ma_feng_wang_bao_conditional_alpha_test(fit, component = "max")
#' @export
ma_feng_wang_bao_conditional_alpha_test <- function(
    fit, factor_count = NULL, component = c("max", "sum", "adaptive")) {
  call <- match.call()
  component <- match.arg(component)
  fit <- .ch4_afc_fit(fit)
  if (component == "sum") {
    answer <- ma_lan_su_tsai_conditional_alpha_sum_test(fit)
    answer$call <- call
    return(answer)
  }
  factor_count <- .ch4_afc_factor_count(fit, factor_count)
  core <- .ch4_afc_light_components(
    fit, factor_count = factor_count, compute_max = TRUE,
    require_sum = component == "adaptive"
  )
  estimate <- .ch4_afc_average_alpha(fit)
  final.statistic <- stats::setNames(
    core$maximum_centered, "Gumbel centered"
  )
  final.p <- core$maximum_p
  final.method <- "Ma-Feng-Wang-Bao conditional maximum alpha test"
  fisher <- NULL
  if (component == "adaptive") {
    fisher.statistic <- -2 * (core$maximum_log_p + core$sum_log_p)
    fisher.p <- if (is.infinite(fisher.statistic)) {
      0
    } else {
      stats::pchisq(fisher.statistic, df = 4, lower.tail = FALSE)
    }
    fisher <- list(
      statistic = fisher.statistic,
      p.value = fisher.p,
      component.p.values = c(max = core$maximum_p, sum = core$sum_p),
      component.log.p.values = c(
        max = core$maximum_log_p, sum = core$sum_log_p
      )
    )
    final.statistic <- stats::setNames(fisher.statistic, "Fisher chi-square")
    final.p <- fisher.p
    final.method <- paste(
      "Ma-Feng-Wang-Bao adaptive conditional alpha test",
      "(primary Fisher max/sum combination)"
    )
  }
  .ch4_afc_htest(
    statistic = final.statistic,
    p.value = final.p,
    method = final.method,
    data.name = fit$data.name,
    estimate = estimate,
    raw.statistic = c(
      M.NT = unname(core$maximum_statistic),
      Gumbel.centered = core$maximum_centered,
      HDA.Z = core$sum_z
    ),
    components = list(
      coordinate.t.squared = core$coordinate_t_squared,
      marginal.variance.scaled = core$marginal_variance,
      marginal.df = unname(core$marginal_df),
      maximum = unname(core$maximum_statistic),
      maximum.centered = core$maximum_centered,
      maximum.p.value = core$maximum_p,
      maximum.log.p.value = core$maximum_log_p,
      sum.Z = core$sum_z,
      sum.p.value = core$sum_p,
      sum.log.p.value = core$sum_log_p,
      S.NT.scaled = unname(core$sum_statistic),
      mu.hat.scaled = unname(core$mean_estimate),
      variance.hat.scaled = unname(core$variance_estimate),
      trace.Sigma.squared.estimate =
        unname(core$trace_sigma_squared_estimate),
      fisher = fisher,
      h = fit$h,
      factor.count = factor_count,
      design.columns = fit$design.columns,
      response.scale = fit$response.scale,
      T = fit$T,
      N = fit$N
    ),
    diagnostics = list(
      maximum.calibration = "upper-tail type-I extreme value",
      sum.calibration = "upper-tail standard normal",
      adaptive.calibration = if (component == "adaptive") {
        "primary Fisher chi-square with 4 df"
      } else {
        NULL
      },
      primary.marginal.divisor = "T - d - 1",
      book.cauchy.attribution.conflict = TRUE,
      ridge = "none",
      variance.floor = "none",
      projection = fit$diagnostics
    ),
    call = call
  )
}

.ch4_afc_css_core <- function(fit, trace_fit) {
  fit <- .ch4_afc_fit(fit)
  trace_fit <- .ch4_afc_fit(trace_fit, "trace_fit")
  if (fit$T != trace_fit$T || fit$N != trace_fit$N) {
    stop("fit and trace_fit must describe the same observations and assets.",
         call. = FALSE)
  }
  if (!identical(fit$observation.names, trace_fit$observation.names)) {
    stop("fit and trace_fit use different observation orders.", call. = FALSE)
  }
  if (!identical(fit$asset.names, trace_fit$asset.names)) {
    stop("fit and trace_fit use different asset orders.", call. = FALSE)
  }
  if (!identical(trace_fit$center.alpha, FALSE)) {
    stop(paste0(
      "trace_fit must use the primary uncentered alpha-basis design ",
      "(center_alpha = FALSE)."
    ), call. = FALSE)
  }
  core <- cpp_ch4_afc_css_components(
    fit$residuals.scaled, trace_fit$residuals.scaled, fit$h
  )
  trace_estimate <- unname(core$trace_estimate)
  if (!is.finite(trace_estimate) || trace_estimate <= 0) {
    stop(paste0(
      "The primary CSS trace estimate is not strictly positive and finite; ",
      "no absolute value or floor is used."
    ), call. = FALSE)
  }
  statistic <- unname(core$numerator) / sqrt(trace_estimate)
  if (!is.finite(statistic)) {
    stop("The standardized CSS statistic is not finite.", call. = FALSE)
  }
  core$statistic <- statistic
  core$p.value <- stats::pnorm(statistic, lower.tail = FALSE)
  core$log.p.value <- stats::pnorm(
    statistic, lower.tail = FALSE, log.p = TRUE
  )
  core
}


#' Zhao robust spatial-sign sum test for conditional alpha
#'
#' Implements the spatial-sign sum test using a restricted null sieve `fit`
#' and the separate unrestricted, uncentered-alpha-basis `trace_fit` required
#' by the primary feasible trace. With signs \eqn{U_t} from the restricted
#' residuals and \eqn{h=M_Z1_T}, its centered numerator is
#' \deqn{(h'h)^{-1}h'UU'h-1.}
#' If \eqn{\widetilde U_t} denotes the signs from `trace_fit`, then
#' \deqn{\widehat{\operatorname{tr}(\Sigma_u^2)}=
#' \frac{\sum_{t\ne s}h_t^2h_s^2
#' (\widetilde U_t'\widetilde U_s)^2}
#' {(h'h)(h'h-1)}}
#' and the statistic divides the centered numerator by the square root of
#' this estimate.
#'
#' The book calls this a leave-two-out estimator but does not give the
#' required unrestricted residual construction, and it inserts an additional
#' factor of two under the square root. The primary Zhao formula, reproduced
#' explicitly in Zhao and Wang (2026), has no such factor. Exact zero rows use
#' the paper's definition \eqn{U(0)=0}; a non-positive final trace still fails.
#'
#' @param fit Restricted centered-basis fit from
#'   [conditional_alpha_sieve_fit()].
#' @param trace_fit Fit to the same returns using the corresponding design
#'   from [conditional_alpha_sieve_design()] with `center_alpha = FALSE`.
#' @return A `conditional_alpha_test`/`htest` object with the normal
#'   upper-tail calibration, both sign matrices, and the exact trace numerator
#'   and denominator.
#' @references
#' Zhao, P. (2023). Robust high-dimensional alpha test for conditional
#' time-varying factor models. *Statistics*, 57, 444--457.
#' \doi{10.1080/02331888.2023.2180003}.
#'
#' Zhao, P. and Wang, H. (2026). Robust spatial-sign-based testing of
#' high-dimensional alpha in conditional factor models.
#' \url{https://arxiv.org/abs/2604.12252}.
#' @examples
#' tt <- seq(0, 1, length.out = 18)
#' b <- cbind(tt, tt^2)
#' y <- cbind(sin(1:18), cos(1:18), sin(1:18 / 2))
#' null_fit <- conditional_alpha_sieve_fit(
#'   y, conditional_alpha_sieve_design(b)
#' )
#' trace_fit <- conditional_alpha_sieve_fit(
#'   y, conditional_alpha_sieve_design(b, center_alpha = FALSE)
#' )
#' zhao_conditional_spatial_sign_sum_test(null_fit, trace_fit)
#' @export
zhao_conditional_spatial_sign_sum_test <- function(fit, trace_fit) {
  call <- match.call()
  fit <- .ch4_afc_fit(fit)
  trace_fit <- .ch4_afc_fit(trace_fit, "trace_fit")
  core <- .ch4_afc_css_core(fit, trace_fit)
  numerator_signs <- as.matrix(core$numerator_signs)
  trace_signs <- as.matrix(core$trace_signs)
  dimnames(numerator_signs) <- list(
    fit$observation.names, fit$asset.names
  )
  dimnames(trace_signs) <- list(
    fit$observation.names, fit$asset.names
  )
  .ch4_afc_htest(
    statistic = stats::setNames(core$statistic, "CSS Z"),
    p.value = core$p.value,
    method = "Zhao robust spatial-sign sum test for conditional alpha",
    data.name = fit$data.name,
    estimate = .ch4_afc_average_alpha(fit),
    raw.statistic = c(
      sign.quadratic = unname(core$quadratic),
      centered.numerator = unname(core$numerator)
    ),
    components = list(
      sign.quadratic = unname(core$quadratic),
      centered.numerator = unname(core$numerator),
      trace.Sigma.u.squared = unname(core$trace_estimate),
      trace.numerator = unname(core$trace_numerator),
      trace.denominator = unname(core$trace_denominator),
      numerator.signs = numerator_signs,
      trace.signs = trace_signs,
      h = fit$h,
      h2 = fit$h2,
      numerator.zero.rows = as.integer(core$numerator_zero_rows),
      trace.zero.rows = as.integer(core$trace_zero_rows),
      T = fit$T,
      N = fit$N
    ),
    diagnostics = list(
      calibration = "asymptotic upper-tail standard normal",
      trace.residual.design =
        "unrestricted uncentered alpha-basis sieve fit",
      primary.trace.denominator = "(h'h)(h'h-1)",
      primary.variance.factor = 1,
      book.sqrt.two.conflict = TRUE,
      book.leave.two.out.description.conflict = TRUE,
      zero.sign.definition = "U(0)=0",
      trace.floor = "none",
      restricted.projection = fit$diagnostics,
      trace.projection = trace_fit$diagnostics
    ),
    call = call
  )
}

.ch4_afc_csm_core <- function(fit, controls) {
  fit <- .ch4_afc_fit(fit)
  spatial <- cpp_scaled_spatial_median(
    fit$residuals.scaled, numeric(fit$N), TRUE,
    controls$tol, controls$max_iter, controls$zero_tol
  )
  if (!isTRUE(spatial$iteration_stable)) {
    stop(sprintf(
      paste0("The conditional CSM scaled spatial median did not stabilize ",
             "within %d updates (score residual %.6g, relative update %.6g)."),
      as.integer(spatial$max_iterations),
      as.numeric(spatial$score_residual),
      as.numeric(spatial$relative_update)
    ), call. = FALSE)
  }
  radii <- as.numeric(spatial$radii)
  if (anyNA(radii) || any(!is.finite(radii)) || any(radii <= 0)) {
    stop(paste0(
      "Every conditional CSM standardized radius must be strictly positive ",
      "and finite; no deletion, cap, or floor is used."
    ), call. = FALSE)
  }
  radial <- c(
    mean.r.squared = mean(radii^2),
    mean.r = mean(radii),
    mean.inverse.r = mean(1 / radii)
  )
  if (any(!is.finite(radial))) {
    stop("A required conditional CSM radial moment is not finite.",
         call. = FALSE)
  }
  omega_ratio <- fit$h2 / fit$T
  projection_loss <- 1 - omega_ratio
  denominator <- unname(
    1 - 2 * projection_loss * radial["mean.inverse.r"] *
      radial["mean.r"] + projection_loss * radial["mean.r.squared"] *
      radial["mean.inverse.r"]^2
  )
  if (!is.finite(denominator) || denominator <= 0) {
    stop(paste0(
      "The primary conditional CSM zeta denominator is not strictly ",
      "positive and finite; no absolute value or floor is used."
    ), call. = FALSE)
  }
  zeta <- fit$N * unname(radial["mean.inverse.r"])^2 / denominator
  coordinate <- as.numeric(spatial$standardized_location_coordinate)
  raw <- fit$T * zeta * max(coordinate^2)
  centered <- raw - 2 * log(fit$N) + log(log(fit$N))
  if (!is.finite(zeta) || zeta <= 0 || !is.finite(raw) ||
      !is.finite(centered)) {
    stop("The conditional CSM statistic is not finite.", call. = FALSE)
  }
  tail <- .ch4_afc_gumbel_tail(centered)
  theta <- as.numeric(spatial$location_relative_origin) * fit$response.scale
  log_diagonal <- as.numeric(spatial$log_diagonal_input_canonical) +
    2 * log(fit$response.scale)
  log_diagonal <- log_diagonal - max(log_diagonal)
  diagonal <- exp(log_diagonal)
  directions <- as.matrix(spatial$directions)
  dimnames(directions) <- list(fit$observation.names, fit$asset.names)
  names(theta) <- names(coordinate) <- names(log_diagonal) <-
    names(diagonal) <- fit$asset.names
  list(
    raw = raw,
    centered = centered,
    p.value = tail$p.value,
    log.p.value = tail$log.p.value,
    log.intensity = tail$log.intensity,
    theta = theta,
    standardized.theta = coordinate,
    scale.diagonal = diagonal,
    log.scale.diagonal = log_diagonal,
    directions = directions,
    radii = radii,
    radial.moments = radial,
    omega.ratio = omega_ratio,
    zeta = zeta,
    zeta.denominator = denominator,
    fit = spatial
  )
}


#' Zhao--Wang robust conditional maximum and adaptive alpha test
#'
#' Implements the CSM maximum of Zhao and Wang (2026) on a restricted
#' conditional sieve fit. A simultaneous scaled spatial median gives
#' \eqn{\widehat\theta}, \eqn{\widehat D}, and standardized radii. With
#' \eqn{\omega_T=h'h}, the primary nuisance estimate is
#' \deqn{\widehat\zeta=
#' \frac{N\overline{r^{-1}}^2}
#' {1-2(1-\omega_T/T)\overline{r^{-1}}\bar r+
#' (1-\omega_T/T)\overline{r^2}\,\overline{r^{-1}}^2},}
#' and the centered statistic is
#' \deqn{T\widehat\zeta
#' \|\widehat D^{-1/2}\widehat\theta\|_\infty^2
#' -2\log N+\log\log N.}
#'
#' For `component = "combined"`, the CSM p-value and the exact Zhao CSS
#' component from [zhao_conditional_spatial_sign_sum_test()] are combined by
#' the paper's *truncated* Cauchy rule: a component contributes only when its
#' p-value is below one half. This current primary source resolves the book's
#' placeholder `ZhaoWang2024ConditionalMaxAlpha` citation. The book's CSM
#' radial formula is otherwise consistent, while its CSS square-root factor
#' is not.
#'
#' @param fit Restricted centered-basis conditional sieve fit.
#' @param trace_fit For `"sum"` or `"combined"`, the corresponding
#'   uncentered-alpha-basis fit required by CSS.
#' @param component One of `"max"`, `"sum"`, or `"combined"`.
#' @param tol Positive simultaneous scaled-median equation tolerance.
#' @param max_iter Positive maximum update count.
#' @param zero_tol Non-negative exact-zero radius tolerance.
#' @return A `conditional_alpha_test`/`htest` object containing the full CSM
#'   fit, radial correction, CSS component, and truncated-Cauchy diagnostics.
#' @references
#' Zhao, P. and Wang, H. (2026). Robust spatial-sign-based testing of
#' high-dimensional alpha in conditional factor models.
#' \url{https://arxiv.org/abs/2604.12252}.
#' @examples
#' tt <- seq(0, 1, length.out = 18)
#' b <- cbind(tt, tt^2)
#' y <- cbind(sin(1:18), cos(1:18), sin(1:18 / 2))
#' fit <- conditional_alpha_sieve_fit(
#'   y, conditional_alpha_sieve_design(b)
#' )
#' zhao_wang_conditional_spatial_sign_test(fit)
#' @export
zhao_wang_conditional_spatial_sign_test <- function(
    fit, trace_fit = NULL, component = c("max", "sum", "combined"),
    tol = 1e-7, max_iter = 500L, zero_tol = 0) {
  call <- match.call()
  component <- match.arg(component)
  fit <- .ch4_afc_fit(fit)
  if (component == "sum") {
    if (is.null(trace_fit)) {
      stop("trace_fit is required for component = 'sum'.", call. = FALSE)
    }
    answer <- zhao_conditional_spatial_sign_sum_test(fit, trace_fit)
    answer$call <- call
    return(answer)
  }
  controls <- .ch4_afc_controls(tol, max_iter, zero_tol)
  csm <- .ch4_afc_csm_core(fit, controls)
  css <- NULL
  combination <- NULL
  final.statistic <- stats::setNames(csm$centered, "Gumbel centered")
  final.p <- csm$p.value
  final.method <- "Zhao-Wang robust conditional spatial-sign maximum alpha test"
  if (component == "combined") {
    if (is.null(trace_fit)) {
      stop("trace_fit is required for component = 'combined'.",
           call. = FALSE)
    }
    trace_fit <- .ch4_afc_fit(trace_fit, "trace_fit")
    css <- .ch4_afc_css_core(fit, trace_fit)
    combination <- .ch4_afc_truncated_cauchy(c(
      CSS = css$p.value, CSM = csm$p.value
    ))
    names(combination$terms) <- names(combination$active) <- c("CSS", "CSM")
    final.statistic <- stats::setNames(
      combination$statistic, "truncated Cauchy"
    )
    final.p <- combination$p.value
    final.method <- paste(
      "Zhao-Wang robust adaptive conditional alpha test",
      "(primary truncated-Cauchy CSS/CSM combination)"
    )
  }
  .ch4_afc_htest(
    statistic = final.statistic,
    p.value = final.p,
    method = final.method,
    data.name = fit$data.name,
    estimate = csm$theta,
    raw.statistic = c(
      CSM.raw = csm$raw,
      Gumbel.centered = csm$centered
    ),
    components = list(
      CSM.raw = csm$raw,
      CSM.centered = csm$centered,
      CSM.p.value = csm$p.value,
      CSM.log.p.value = csm$log.p.value,
      theta = csm$theta,
      standardized.theta = csm$standardized.theta,
      scale.diagonal = csm$scale.diagonal,
      log.scale.diagonal = csm$log.scale.diagonal,
      directions = csm$directions,
      radii = csm$radii,
      radial.moments = csm$radial.moments,
      omega.ratio = csm$omega.ratio,
      zeta = csm$zeta,
      zeta.denominator = csm$zeta.denominator,
      CSS = if (is.null(css)) NULL else list(
        statistic = css$statistic,
        p.value = css$p.value,
        log.p.value = css$log.p.value,
        sign.quadratic = unname(css$quadratic),
        centered.numerator = unname(css$numerator),
        trace.Sigma.u.squared = unname(css$trace_estimate),
        trace.numerator = unname(css$trace_numerator),
        trace.denominator = unname(css$trace_denominator)
      ),
      truncated.Cauchy = combination,
      h = fit$h,
      h2 = fit$h2,
      T = fit$T,
      N = fit$N
    ),
    diagnostics = list(
      CSM.calibration = "upper-tail type-I extreme value",
      CSS.calibration = if (component == "combined") {
        "upper-tail standard normal"
      } else {
        NULL
      },
      combination = if (component == "combined") {
        "primary truncated Cauchy; components with p >= 1/2 contribute zero"
      } else {
        NULL
      },
      current.primary = "Zhao-Wang arXiv:2604.12252 (2026)",
      book.placeholder.citation.resolved = TRUE,
      CSS.book.sqrt.two.conflict = component == "combined",
      converged = isTRUE(csm$fit$iteration_stable),
      iterations = as.integer(csm$fit$iterations),
      score.residual = unname(csm$fit$score_residual),
      relative.update = unname(csm$fit$relative_update),
      radial.floor = "none",
      ridge = "none",
      restricted.projection = fit$diagnostics,
      trace.projection = if (is.null(trace_fit)) NULL else trace_fit$diagnostics
    ),
    call = call
  )
}
