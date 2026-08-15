.ch3g_covariance_estimator <- function(estimate, method, call,
                                       sample.covariance, threshold,
                                       components = list(),
                                       diagnostics = list()) {
  if (!is.matrix(estimate) || nrow(estimate) != ncol(estimate) ||
      anyNA(estimate) || any(!is.finite(estimate))) {
    stop("The covariance estimate is not a finite square matrix.", call. = FALSE)
  }
  result <- list(
    estimate = estimate,
    sample.covariance = sample.covariance,
    threshold = threshold,
    method = method,
    components = components,
    diagnostics = diagnostics,
    call = call
  )
  class(result) <- c("hd_covariance_estimator", "list")
  result
}


.ch3g_divisor <- function(divisor, n) {
  divisor <- match.arg(divisor, c("n", "n-1"))
  value <- if (divisor == "n") n else n - 1
  list(label = divisor, value = value)
}


.ch3g_threshold_rule <- function(rule, scad_a, adaptive_eta) {
  rule <- match.arg(rule, c("hard", "soft", "scad", "adaptive_lasso"))
  if (rule == "scad") {
    if (!is.numeric(scad_a) || length(scad_a) != 1L ||
        !is.finite(scad_a) || scad_a <= 2) {
      stop("`scad_a` must be one finite number greater than 2.", call. = FALSE)
    }
    shape <- scad_a
  } else if (rule == "adaptive_lasso") {
    if (!is.numeric(adaptive_eta) || length(adaptive_eta) != 1L ||
        !is.finite(adaptive_eta) || adaptive_eta < 0) {
      stop("`adaptive_eta` must be one finite non-negative number.",
           call. = FALSE)
    }
    shape <- adaptive_eta
  } else {
    shape <- 1
  }
  list(rule = rule, shape = shape)
}


.ch3g_threshold_inputs <- function(x, center, divisor) {
  x <- .as_data_matrix(x, "x", min_rows = 2L)
  center <- .ch3g_validate_center(center)
  n <- nrow(x)
  scale <- max(abs(x))
  if (!is.finite(scale) || scale <= 0) {
    stop("`x` must contain nonzero finite data.", call. = FALSE)
  }
  xs <- x / scale
  residual <- if (center) sweep(xs, 2L, colMeans(xs), "-") else xs
  divisor.fit <- .ch3g_divisor(divisor, n)
  covariance.scaled <- crossprod(residual) / divisor.fit$value
  list(
    x = x, residual = residual, covariance.scaled = covariance.scaled,
    scale = scale, n = n, p = ncol(x), center = center,
    divisor = divisor.fit
  )
}


.ch3g_restore_covariance_scale <- function(matrix, scale, what) {
  result <- matrix * scale * scale
  if (any(!is.finite(result))) {
    stop(sprintf("The %s overflows in the original measurement units.", what),
         call. = FALSE)
  }
  result
}


.ch3g_restore_fourth_order_scale <- function(value, scale, what) {
  result <- value * scale
  result <- result * scale
  result <- result * scale
  result <- result * scale
  invisible(what)
  result
}


#' Bickel--Levina hard-thresholded covariance estimator
#'
#' Forms the centered sample covariance with divisor \eqn{n}, as in Bickel
#' and Levina (2008), and applies
#' \deqn{s_{ij}\,1\{|s_{ij}|>\lambda\}.}
#' When `threshold = NULL`, \eqn{\lambda=C\sqrt{\log(p)/n}}. The book draft
#' first defines a divisor-\eqn{n-1} covariance and then cites the primary
#' divisor-\eqn{n} formula. `divisor` makes this finite-sample distinction
#' explicit; the default is the primary-paper convention. Thresholding need
#' not preserve positive definiteness and no eigenvalue repair is performed.
#'
#' @param x Numeric matrix with observations in rows.
#' @param threshold A finite non-negative threshold in covariance units.
#'   `NULL` uses the rate threshold controlled by `constant`.
#' @param constant A finite non-negative multiplier for the rate threshold.
#' @param center Whether to remove column means.
#' @param divisor Either `"n"` (primary-paper default) or `"n-1"`.
#' @param threshold_diagonal Whether to apply the rule to diagonal entries.
#'
#' @return An `hd_covariance_estimator` list whose `estimate` field is the
#'   thresholded covariance matrix.
#' @references Bickel, P. J. and Levina, E. (2008). *Annals of Statistics*,
#'   36, 2577--2604.
#' @examples
#' set.seed(31)
#' x <- matrix(rnorm(40), nrow = 10, ncol = 4)
#' bickel_levina_covariance_threshold(x, threshold = 0.2)
#' @export
bickel_levina_covariance_threshold <- function(
    x, threshold = NULL, constant = 1, center = TRUE,
    divisor = c("n", "n-1"), threshold_diagonal = TRUE) {
  call <- match.call()
  input <- .ch3g_threshold_inputs(x, center, divisor)
  if (is.null(threshold)) {
    if (!is.numeric(constant) || length(constant) != 1L ||
        !is.finite(constant) || constant < 0) {
      stop("`constant` must be one finite non-negative number.", call. = FALSE)
    }
    lambda <- constant * sqrt(log(input$p) / input$n)
    lambda.scaled <- lambda / input$scale / input$scale
    source <- "constant * sqrt(log(p) / n)"
  } else {
    if (!is.numeric(threshold) || length(threshold) != 1L ||
        !is.finite(threshold) || threshold < 0) {
      stop("`threshold` must be one finite non-negative number.", call. = FALSE)
    }
    lambda <- threshold
    lambda.scaled <- threshold / input$scale / input$scale
    source <- "supplied"
  }
  threshold.matrix <- matrix(lambda.scaled, input$p, input$p)
  estimate.scaled <- cpp_ch3_gaussian_threshold_matrix(
    input$covariance.scaled, threshold.matrix, "hard", 1,
    threshold_diagonal
  )
  estimate <- .ch3g_restore_covariance_scale(
    estimate.scaled, input$scale, "thresholded covariance estimate"
  )
  sample.covariance <- .ch3g_restore_covariance_scale(
    input$covariance.scaled, input$scale, "sample covariance"
  )
  .ch3g_covariance_estimator(
    estimate, "Bickel-Levina hard-thresholded covariance estimator", call,
    sample.covariance, matrix(lambda, input$p, input$p),
    components = list(rule = "hard"),
    diagnostics = list(
      n = input$n, p = input$p, center = input$center,
      covariance.divisor = input$divisor$value,
      divisor.convention = input$divisor$label,
      threshold.source = source,
      threshold.diagonal = threshold_diagonal,
      internal.data.scale = input$scale,
      fourth.order.components = paste(
        "theta is reported in original units when representable;",
        "the internally scaled value remains the computational source"
      ),
      positive.definiteness.repair = "none"
    )
  )
}


#' Rothman--Levina--Zhu generalized covariance thresholding
#'
#' Applies a generalized thresholding map entrywise to the divisor-\eqn{n}
#' centered sample covariance. The four rules in Rothman, Levina and Zhu
#' (2009) are hard, soft, SCAD, and adaptive lasso. Their adaptive-lasso rule
#' is
#' \deqn{\operatorname{sign}(z)
#' (|z|-\lambda^{\eta+1}|z|^{-\eta})_+.}
#' The accompanying book calls the fourth rule MCP; that attribution is
#' incorrect, so this API does not silently rename adaptive lasso as MCP.
#'
#' @inheritParams bickel_levina_covariance_threshold
#' @param rule One of `"hard"`, `"soft"`, `"scad"`, or
#'   `"adaptive_lasso"`.
#' @param scad_a SCAD shape parameter, greater than two.
#' @param adaptive_eta Non-negative adaptive-lasso exponent.
#'
#' @return An `hd_covariance_estimator` object.
#' @references Rothman, A. J., Levina, E. and Zhu, J. (2009). *Journal of
#'   the American Statistical Association*, 104, 177--186.
#'   \doi{10.1198/jasa.2009.0101}.
#' @examples
#' set.seed(32)
#' x <- matrix(rnorm(40), nrow = 10, ncol = 4)
#' rothman_levina_zhu_covariance_threshold(
#'   x, threshold = 0.2, rule = "soft"
#' )
#' @export
rothman_levina_zhu_covariance_threshold <- function(
    x, threshold = NULL, constant = 1,
    rule = c("hard", "soft", "scad", "adaptive_lasso"),
    scad_a = 3.7, adaptive_eta = 1, center = TRUE,
    divisor = c("n", "n-1"), threshold_diagonal = TRUE) {
  call <- match.call()
  input <- .ch3g_threshold_inputs(x, center, divisor)
  rule.fit <- .ch3g_threshold_rule(rule, scad_a, adaptive_eta)
  if (is.null(threshold)) {
    if (!is.numeric(constant) || length(constant) != 1L ||
        !is.finite(constant) || constant < 0) {
      stop("`constant` must be one finite non-negative number.", call. = FALSE)
    }
    lambda <- constant * sqrt(log(input$p) / input$n)
    lambda.scaled <- lambda / input$scale / input$scale
    source <- "constant * sqrt(log(p) / n)"
  } else {
    if (!is.numeric(threshold) || length(threshold) != 1L ||
        !is.finite(threshold) || threshold < 0) {
      stop("`threshold` must be one finite non-negative number.", call. = FALSE)
    }
    lambda <- threshold
    lambda.scaled <- threshold / input$scale / input$scale
    source <- "supplied"
  }
  threshold.matrix <- matrix(lambda.scaled, input$p, input$p)
  estimate.scaled <- cpp_ch3_gaussian_threshold_matrix(
    input$covariance.scaled, threshold.matrix, rule.fit$rule,
    rule.fit$shape, threshold_diagonal
  )
  estimate <- .ch3g_restore_covariance_scale(
    estimate.scaled, input$scale, "thresholded covariance estimate"
  )
  sample.covariance <- .ch3g_restore_covariance_scale(
    input$covariance.scaled, input$scale, "sample covariance"
  )
  .ch3g_covariance_estimator(
    estimate, "Rothman-Levina-Zhu generalized covariance thresholding", call,
    sample.covariance, matrix(lambda, input$p, input$p),
    components = list(
      rule = rule.fit$rule,
      shape.parameter = if (rule.fit$rule %in% c("scad", "adaptive_lasso"))
        rule.fit$shape else NULL
    ),
    diagnostics = list(
      n = input$n, p = input$p, center = input$center,
      covariance.divisor = input$divisor$value,
      divisor.convention = input$divisor$label,
      threshold.source = source,
      threshold.diagonal = threshold_diagonal,
      book.MCP.attribution = "not used; primary method is adaptive lasso",
      internal.data.scale = input$scale,
      positive.definiteness.repair = "none"
    )
  )
}


#' Cai--Liu adaptive covariance thresholding
#'
#' The primary-paper pilot covariance and variability estimate are
#' \deqn{\hat\sigma_{ij}=n^{-1}\sum_k Z_{ki}Z_{kj},\qquad
#' \hat\theta_{ij}=n^{-1}\sum_k(Z_{ki}Z_{kj}-\hat\sigma_{ij})^2,}
#' with entrywise threshold
#' \deqn{\lambda_{ij}=\delta
#' \sqrt{\hat\theta_{ij}\log(p)/n}.}
#' Both quantities therefore use the same divisor-\eqn{n} covariance by
#' default. Selecting `divisor = "n-1"` is explicit and is reported in the
#' diagnostics; it is not the primary finite-sample formula.
#'
#' @inheritParams rothman_levina_zhu_covariance_threshold
#' @param delta Finite non-negative adaptive-threshold multiplier.
#'
#' @return An `hd_covariance_estimator` object retaining
#'   \eqn{\hat\theta} and every entrywise threshold.
#' @references Cai, T. and Liu, W. (2011). *Journal of the American
#'   Statistical Association*, 106, 672--684. arXiv:1102.2237.
#' @examples
#' set.seed(33)
#' x <- matrix(rnorm(40), nrow = 10, ncol = 4)
#' cai_liu_adaptive_covariance_threshold(x, delta = 2)
#' @export
cai_liu_adaptive_covariance_threshold <- function(
    x, delta = 2,
    rule = c("hard", "soft", "scad", "adaptive_lasso"),
    scad_a = 3.7, adaptive_eta = 1, center = TRUE,
    divisor = c("n", "n-1"), threshold_diagonal = TRUE) {
  call <- match.call()
  input <- .ch3g_threshold_inputs(x, center, divisor)
  rule.fit <- .ch3g_threshold_rule(rule, scad_a, adaptive_eta)
  if (!is.numeric(delta) || length(delta) != 1L ||
      !is.finite(delta) || delta < 0) {
    stop("`delta` must be one finite non-negative number.", call. = FALSE)
  }
  theta.scaled <- cpp_ch3_gaussian_product_variability(
    input$residual, input$covariance.scaled
  )
  if (any(!is.finite(theta.scaled)) || any(theta.scaled < 0)) {
    stop("The Cai--Liu variability matrix is not finite and non-negative.",
         call. = FALSE)
  }
  threshold.scaled <- delta * sqrt(theta.scaled * log(input$p) / input$n)
  estimate.scaled <- cpp_ch3_gaussian_threshold_matrix(
    input$covariance.scaled, threshold.scaled, rule.fit$rule,
    rule.fit$shape, threshold_diagonal
  )
  estimate <- .ch3g_restore_covariance_scale(
    estimate.scaled, input$scale, "adaptive covariance estimate"
  )
  sample.covariance <- .ch3g_restore_covariance_scale(
    input$covariance.scaled, input$scale, "sample covariance"
  )
  theta <- .ch3g_restore_fourth_order_scale(
    theta.scaled, input$scale, "adaptive variability matrix"
  )
  threshold <- .ch3g_restore_covariance_scale(
    threshold.scaled, input$scale, "adaptive threshold matrix"
  )
  .ch3g_covariance_estimator(
    estimate, "Cai-Liu adaptive covariance thresholding", call,
    sample.covariance, threshold,
    components = list(
      rule = rule.fit$rule,
      delta = delta,
      theta = theta,
      shape.parameter = if (rule.fit$rule %in% c("scad", "adaptive_lasso"))
        rule.fit$shape else NULL
    ),
    diagnostics = list(
      n = input$n, p = input$p, center = input$center,
      covariance.divisor = input$divisor$value,
      divisor.convention = input$divisor$label,
      theta.outer.divisor = input$n,
      primary.finite.sample.formula = identical(input$divisor$label, "n"),
      threshold.diagonal = threshold_diagonal,
      internal.data.scale = input$scale,
      positive.definiteness.repair = "none"
    )
  )
}


#' POET covariance estimator with a supplied factor count
#'
#' Implements Principal Orthogonal complEment Thresholding (POET). The first
#' `factors` sample principal components form the low-rank part; a generalized
#' threshold is applied only to the off-diagonal entries of the orthogonal
#' residual covariance, and the two parts are added back. The factor count is
#' deliberately required: the cited display in the book does not specify a
#' unique data-driven selector.
#'
#' `threshold_type = "correlation"` uses primary equation (2.6),
#' \eqn{\lambda_{ij}=\tau\sqrt{r_{ii}r_{jj}}}. `"variance"` uses equation
#' (3.2), \eqn{C\sqrt{\hat\theta_{ij}}
#' \{p^{-1/2}+\sqrt{\log(p)/n}\}}. Residual diagonal entries are always
#' retained. No positive-definiteness repair, ridge, or automatic tuning is
#' applied.
#'
#' @param x Numeric matrix with observations in rows.
#' @param factors A supplied non-negative integer number of factors.
#' @param threshold Finite non-negative \eqn{\tau} or \eqn{C}, according to
#'   `threshold_type`.
#' @param threshold_type Either `"correlation"` or `"variance"`.
#' @param rule,scad_a,adaptive_eta Generalized thresholding rule controls.
#' @param center Whether to remove column means.
#'
#' @return An `hd_covariance_estimator` object containing the low-rank part,
#'   raw and thresholded orthogonal complements, eigenvalues/eigenvectors,
#'   and thresholds.
#' @references Fan, J., Liao, Y. and Mincheva, M. (2013). *Journal of the
#'   Royal Statistical Society: Series B*, 75, 603--680. arXiv:1201.0175.
#' @examples
#' set.seed(34)
#' x <- matrix(rnorm(60), nrow = 15, ncol = 4)
#' poet_covariance(x, factors = 1, threshold = 0.2)
#' @export
poet_covariance <- function(
    x, factors, threshold,
    threshold_type = c("correlation", "variance"),
    rule = c("hard", "soft", "scad", "adaptive_lasso"),
    scad_a = 3.7, adaptive_eta = 1, center = TRUE) {
  call <- match.call()
  x <- .as_data_matrix(x, "x", min_rows = 2L)
  center <- .ch3g_validate_center(center)
  n <- nrow(x)
  p <- ncol(x)
  if (!is.numeric(factors) || length(factors) != 1L ||
      !is.finite(factors) || factors != as.integer(factors) || factors < 0) {
    stop("`factors` must be one non-negative integer.", call. = FALSE)
  }
  factors <- as.integer(factors)
  maximum.factors <- min(n - as.integer(center), p)
  if (factors > maximum.factors) {
    stop("`factors` exceeds the residual sample rank bound.", call. = FALSE)
  }
  if (!is.numeric(threshold) || length(threshold) != 1L ||
      !is.finite(threshold) || threshold < 0) {
    stop("`threshold` must be one finite non-negative number.", call. = FALSE)
  }
  threshold_type <- match.arg(threshold_type)
  rule.fit <- .ch3g_threshold_rule(rule, scad_a, adaptive_eta)

  scale <- max(abs(x))
  if (!is.finite(scale) || scale <= 0) {
    stop("`x` must contain nonzero finite data.", call. = FALSE)
  }
  xs <- x / scale
  residual <- if (center) sweep(xs, 2L, colMeans(xs), "-") else xs
  decomposition <- svd(residual, nu = 0L, nv = maximum.factors)
  eigenvalues <- decomposition$d^2 / n
  vectors <- decomposition$v
  sample.scaled <- crossprod(residual) / n

  if (factors > 0L) {
    leading.vectors <- vectors[, seq_len(factors), drop = FALSE]
    leading.values <- eigenvalues[seq_len(factors)]
    lowrank.scaled <- tcrossprod(
      sweep(leading.vectors, 2L, sqrt(leading.values), "*")
    )
    idiosyncratic.data <- residual -
      residual %*% leading.vectors %*% t(leading.vectors)
  } else {
    leading.vectors <- matrix(numeric(), p, 0L)
    leading.values <- numeric()
    lowrank.scaled <- matrix(0, p, p)
    idiosyncratic.data <- residual
  }
  orthogonal.scaled <- crossprod(idiosyncratic.data) / n
  diagonal <- diag(orthogonal.scaled)
  if (any(!is.finite(diagonal)) || any(diagonal < 0)) {
    stop("The principal orthogonal complement has an invalid diagonal.",
         call. = FALSE)
  }
  if (threshold_type == "correlation") {
    threshold.scaled <- threshold * sqrt(outer(diagonal, diagonal))
    theta.scaled <- NULL
    omega <- NULL
  } else {
    theta.scaled <- cpp_ch3_gaussian_product_variability(
      idiosyncratic.data, orthogonal.scaled
    )
    omega <- 1 / sqrt(p) + sqrt(log(p) / n)
    threshold.scaled <- threshold * sqrt(theta.scaled) * omega
  }
  thresholded.scaled <- cpp_ch3_gaussian_threshold_matrix(
    orthogonal.scaled, threshold.scaled, rule.fit$rule,
    rule.fit$shape, FALSE
  )
  estimate.scaled <- lowrank.scaled + thresholded.scaled

  restore <- function(value, what) {
    .ch3g_restore_covariance_scale(value, scale, what)
  }
  estimate <- restore(estimate.scaled, "POET covariance estimate")
  sample.covariance <- restore(sample.scaled, "sample covariance")
  threshold.original <- restore(threshold.scaled, "POET threshold matrix")
  .ch3g_covariance_estimator(
    estimate, "POET covariance estimator", call,
    sample.covariance, threshold.original,
    components = list(
      factors = factors,
      rule = rule.fit$rule,
      leading.eigenvalues = .ch3g_restore_covariance_scale(
        matrix(leading.values, nrow = 1L), scale,
        "POET leading eigenvalues"
      )[1L, ],
      leading.eigenvectors = leading.vectors,
      low.rank = restore(lowrank.scaled, "POET low-rank component"),
      orthogonal.complement = restore(
        orthogonal.scaled, "POET orthogonal complement"
      ),
      thresholded.orthogonal.complement = restore(
        thresholded.scaled, "POET thresholded orthogonal complement"
      ),
      theta = if (is.null(theta.scaled)) NULL else
        .ch3g_restore_fourth_order_scale(
          theta.scaled, scale, "POET variability matrix"
        ),
      omega = omega
    ),
    diagnostics = list(
      n = n, p = p, center = center,
      covariance.divisor = n,
      factor.count.source = "supplied",
      threshold.type = threshold_type,
      threshold.diagonal = FALSE,
      internal.data.scale = scale,
      fourth.order.components = paste(
        "theta may overflow in original units; the internally scaled",
        "value remains the computational source"
      ),
      positive.definiteness.repair = "none",
      tuning.selection = "none"
    )
  )
}
