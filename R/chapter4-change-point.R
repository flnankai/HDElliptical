# Chapter 4: high-dimensional change-point methods.

.cp_data_matrix <- function(x, minimum_rows = 2L, minimum_columns = 1L) {
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
  if (nrow(x) < minimum_rows || ncol(x) < minimum_columns) {
    stop(sprintf(
      "`x` must have at least %d rows and %d columns.",
      minimum_rows, minimum_columns
    ), call. = FALSE)
  }
  if (anyNA(x) || any(!is.finite(x))) {
    stop("`x` must contain only finite values.", call. = FALSE)
  }
  x
}

.cp_scalar <- function(x, name, lower = -Inf, upper = Inf,
                       lower_open = FALSE, upper_open = FALSE) {
  x <- as.numeric(x)
  ok_lower <- if (lower_open) x > lower else x >= lower
  ok_upper <- if (upper_open) x < upper else x <= upper
  if (length(x) != 1L || is.na(x) || !is.finite(x) ||
      !ok_lower || !ok_upper) {
    stop(sprintf("`%s` is outside its permitted finite range.", name),
         call. = FALSE)
  }
  x
}

.cp_integer <- function(x, name, minimum = 1L) {
  value <- as.numeric(x)
  if (length(value) != 1L || is.na(value) || !is.finite(value) ||
      value != floor(value) || value < minimum ||
      value > .Machine$integer.max) {
    stop(sprintf("`%s` must be one integer not smaller than %d.",
                 name, minimum), call. = FALSE)
  }
  as.integer(value)
}

.cp_flag <- function(x, name) {
  if (!is.logical(x) || length(x) != 1L || is.na(x)) {
    stop(sprintf("`%s` must be TRUE or FALSE.", name), call. = FALSE)
  }
  x
}

.cp_iteration_controls <- function(tol, max_iter, zero_tol) {
  list(
    tol = .cp_scalar(tol, "tol", lower = 0, lower_open = TRUE),
    max_iter = .cp_integer(max_iter, "max_iter"),
    zero_tol = .cp_scalar(zero_tol, "zero_tol", lower = 0)
  )
}

.cp_row_norms <- function(x) {
  answer <- numeric(nrow(x))
  for (ii in seq_len(nrow(x))) {
    largest <- max(abs(x[ii, ]))
    if (!is.finite(largest)) {
      answer[[ii]] <- Inf
    } else if (largest == 0) {
      answer[[ii]] <- 0
    } else {
      answer[[ii]] <- largest * sqrt(sum((x[ii, ] / largest)^2))
    }
  }
  answer
}

.cp_with_seed <- function(seed, expression) {
  if (is.null(seed)) return(force(expression))
  seed <- .cp_integer(seed, "seed", minimum = 0L)
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv)
  on.exit({
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv,
                      inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  set.seed(seed)
  force(expression)
}

.cp_invalid <- function(message, strict, method, call, data.name,
                        diagnostics = list()) {
  if (strict) stop(message, call. = FALSE)
  warning(message, call. = FALSE)
  structure(
    list(
      statistic = c(invalid = NA_real_),
      parameter = numeric(),
      p.value = NA_real_,
      estimate = NULL,
      null.value = c(number_of_changes = 0),
      alternative = "at least one location change",
      method = method,
      data.name = data.name,
      valid = FALSE,
      diagnostics = c(
        list(
          failure = message,
          no.repair = paste(
            "No ridge beyond a method-specified ERHT ridge, variance floor,",
            "absolute-value repair, pseudoinverse, jitter, or residual",
            "perturbation was applied."
          ),
          call = call
        ),
        diagnostics
      )
    ),
    class = c("chapter4_change_point_test", "htest")
  )
}

.cp_gumbel_survival <- function(x, factor = 1) {
  y <- factor * exp(-x)
  out <- -expm1(-y)
  out[x > 745] <- 0
  out[x < -log(.Machine$double.xmax)] <- 1
  pmin(1, pmax(0, out))
}

.cp_gumbel_log_survival <- function(x, factor = 1) {
  if (x > 35) return(log(factor) - x)
  if (x < -log(.Machine$double.xmax)) return(0)
  log(-expm1(-factor * exp(-x)))
}

.cp_A <- function(x) {
  if (length(x) != 1L || !is.finite(x) || x <= 1) {
    stop("The argument of A(x)=sqrt(2 log x) must exceed one.",
         call. = FALSE)
  }
  sqrt(2 * log(x))
}

.cp_D <- function(x) {
  if (length(x) != 1L || !is.finite(x) || x <= 1) {
    stop("The argument of D(x) must exceed one.", call. = FALSE)
  }
  2 * log(x) + 0.5 * log(log(x)) - 0.5 * log(pi)
}

.cp_right_tail_mc <- function(observed, reference) {
  if (!length(reference) || anyNA(reference) || any(!is.finite(reference))) {
    stop("The calibration reference must be finite and non-empty.",
         call. = FALSE)
  }
  (1 + sum(reference >= observed)) / (length(reference) + 1)
}

.cp_psd_root <- function(covariance) {
  covariance <- (covariance + t(covariance)) / 2
  eig <- eigen(covariance, symmetric = TRUE)
  scale <- max(1, max(abs(eig$values)))
  tolerance <- 1000 * .Machine$double.eps * nrow(covariance) * scale
  if (min(eig$values) < -tolerance) {
    stop("A theoretical Gaussian-process covariance is not positive semidefinite.",
         call. = FALSE)
  }
  positive <- eig$values > tolerance
  root <- if (any(positive)) {
    eig$vectors[, positive, drop = FALSE] %*%
      diag(sqrt(eig$values[positive]), sum(positive))
  } else {
    matrix(0, nrow(covariance), 0L)
  }
  list(
    root = root,
    rank = sum(positive),
    minimum.eigenvalue = min(eig$values),
    zero.eigenvalue.tolerance = tolerance,
    note = paste(
      "Algebraic zero eigenvalues of the theoretical covariance are omitted",
      "from its rank factorization; no data-matrix eigenvalue repair is used."
    )
  )
}

.cp_gaussian_max_reference <- function(covariance, draws, absolute = FALSE) {
  draws <- .cp_integer(draws, "calibration_draws")
  factor <- .cp_psd_root(covariance)
  z <- matrix(
    stats::rnorm(draws * factor$rank), nrow = draws, ncol = factor$rank
  )
  paths <- if (factor$rank) z %*% t(factor$root) else
    matrix(0, draws, nrow(covariance))
  statistic <- if (absolute) apply(abs(paths), 1L, max) else
    apply(paths, 1L, max)
  list(statistic = statistic, factor = factor)
}

.cp_cauchy_component <- function(p.value, weight) {
  if (p.value == 0.5) return(list(sign = 0, log.abs = -Inf))
  if (p.value > 1e-8 && p.value < 1 - 1e-8) {
    transformed <- tan(pi * (0.5 - p.value))
    return(list(
      sign = sign(transformed),
      log.abs = log(weight) + log(abs(transformed))
    ))
  }
  if (p.value < 0.5) {
    return(list(
      sign = 1,
      log.abs = log(weight) - log(pi) - log(p.value)
    ))
  }
  list(
    sign = -1,
    log.abs = log(weight) - log(pi) - log1p(-p.value)
  )
}

.cp_signed_log_add <- function(first, second) {
  if (first$sign == 0) return(second)
  if (second$sign == 0) return(first)
  if (first$sign == second$sign) {
    largest <- max(first$log.abs, second$log.abs)
    smallest <- min(first$log.abs, second$log.abs)
    return(list(
      sign = first$sign,
      log.abs = if (is.infinite(largest)) Inf else
        largest + log1p(exp(smallest - largest))
    ))
  }
  if (is.infinite(first$log.abs) && is.infinite(second$log.abs)) {
    stop("The Cauchy combination has indeterminate opposite infinite transforms.",
         call. = FALSE)
  }
  if (identical(first$log.abs, second$log.abs)) {
    return(list(sign = 0, log.abs = -Inf))
  }
  largest <- if (first$log.abs > second$log.abs) first else second
  smallest <- if (first$log.abs > second$log.abs) second else first
  ratio <- exp(smallest$log.abs - largest$log.abs)
  if (ratio == 1) return(list(sign = 0, log.abs = -Inf))
  list(
    sign = largest$sign,
    log.abs = largest$log.abs + log1p(-ratio)
  )
}

.cp_cauchy_combine <- function(p.values, weights = NULL) {
  p.values <- as.numeric(p.values)
  if (!length(p.values) || anyNA(p.values) ||
      any(!is.finite(p.values)) || any(p.values < 0 | p.values > 1)) {
    stop("Cauchy combination requires finite p-values in [0,1].",
         call. = FALSE)
  }
  if (is.null(weights)) weights <- rep(1 / length(p.values), length(p.values))
  weights <- as.numeric(weights)
  if (length(weights) != length(p.values) || anyNA(weights) ||
      any(!is.finite(weights)) || any(weights <= 0) || sum(weights) <= 0) {
    stop("Cauchy weights must be positive, finite, and conformable.",
         call. = FALSE)
  }
  weights <- weights / sum(weights)
  lower <- any(p.values == 0)
  upper <- any(p.values == 1)
  if (lower && upper) {
    stop("The Cauchy combination is indeterminate with both 0 and 1 p-values.",
         call. = FALSE)
  }
  if (lower) return(list(statistic = Inf, p.value = 0))
  if (upper) return(list(statistic = -Inf, p.value = 1))
  combined <- list(sign = 0, log.abs = -Inf)
  for (ii in seq_along(p.values)) {
    combined <- .cp_signed_log_add(
      combined, .cp_cauchy_component(p.values[[ii]], weights[[ii]])
    )
  }
  if (combined$sign == 0) {
    return(list(
      statistic = 0, p.value = 0.5,
      statistic.sign = 0, log.absolute.statistic = -Inf,
      evaluation = "signed-log cotangents with atan2 tail"
    ))
  }
  statistic <- if (combined$log.abs > log(.Machine$double.xmax)) {
    combined$sign * Inf
  } else {
    combined$sign * exp(combined$log.abs)
  }
  inverse <- if (-combined$log.abs > log(.Machine$double.xmax)) Inf else
    exp(-combined$log.abs)
  complement.angle <- atan(inverse)
  p.value <- if (combined$sign > 0) {
    complement.angle / pi
  } else {
    1 - complement.angle / pi
  }
  list(
    statistic = statistic,
    p.value = p.value,
    statistic.sign = combined$sign,
    log.absolute.statistic = combined$log.abs,
    evaluation = "signed-log cotangents with atan2-equivalent tail"
  )
}

.cp_kolmogorov_survival <- function(x, tolerance = 1e-14) {
  if (!is.finite(x) || x < 0) return(NA_real_)
  if (x == 0) return(1)
  if (x < 1) {
    series <- 0
    for (k in seq_len(100000L)) {
      term <- exp(-((2 * k - 1)^2 * pi^2) / (8 * x^2))
      series <- series + term
      if (term < tolerance) break
    }
    cdf <- sqrt(2 * pi) / x * series
    return(pmin(1, pmax(0, 1 - cdf)))
  }
  answer <- 0
  for (k in seq_len(100000L)) {
    term <- 2 * (-1)^(k - 1L) * exp(-2 * k^2 * x^2)
    answer <- answer + term
    if (abs(term) < tolerance) break
  }
  pmin(1, pmax(0, answer))
}

.cp_kolmogorov_quantile <- function(probability) {
  probability <- .cp_scalar(
    probability, "probability", lower = 0, upper = 1,
    lower_open = TRUE, upper_open = TRUE
  )
  stats::uniroot(
    function(q) .cp_kolmogorov_survival(q) - (1 - probability),
    c(1e-8, 10), tol = 1e-11
  )$root
}


#' Classical CUSUM change-point test
#'
#' Computes the univariate CUSUM
#' \deqn{C_n(k)=n^{-1/2}\{S_k-kS_n/n\}}
#' and rejects for a large trimmed maximum of
#' \eqn{\{(k/n)(1-k/n)\}^{-\gamma}|C_n(k)|/\widehat\sigma}.
#' For `gamma = 0`, `calibration = "asymptotic"` uses the classical
#' two-sided Brownian-bridge (Kolmogorov) limit.  `calibration =
#' "grid-gaussian"` simulates the method's finite scan-grid Gaussian bridge;
#' this is null calibration, not replication of a paper simulation study.
#'
#' If `sigma` is not supplied, the scale is either the sample standard
#' deviation or Rice's adjacent-difference estimate.  This software choice is
#' returned explicitly because the book's classical display assumes a
#' population \eqn{\sigma} and does not prescribe its estimator.
#'
#' @param x Numeric vector or one-column matrix.
#' @param gamma Boundary weight in `[0, 1/2]`.
#' @param trim Number of candidate observations removed from each boundary.
#' @param sigma Optional positive known/null standard deviation.
#' @param variance Scale estimator used when `sigma = NULL`.
#' @param calibration Brownian-bridge calibration method.
#' @param calibration_draws Number of intrinsic Gaussian calibration draws.
#' @param alpha Test level.
#' @param seed Optional calibration seed; the previous RNG state is restored.
#' @param strict If `TRUE`, method failures error; otherwise they return an
#'   invalid object with `estimate = NULL`.
#' @return An `htest`-compatible object with the complete CUSUM path.
#' @references Page (1954); Csorgo and Horvath (1997).
#' @export
#' @examples
#' x <- c(-1, 0, 1, -1, 0, 3, 4, 3)
#' classical_cusum_test(x)
classical_cusum_test <- function(
    x,
    gamma = 0,
    trim = 1L,
    sigma = NULL,
    variance = c("difference", "sample"),
    calibration = c("asymptotic", "grid-gaussian"),
    calibration_draws = 4999L,
    alpha = 0.05,
    seed = NULL,
    strict = TRUE) {
  call <- match.call()
  data.name <- deparse(substitute(x))
  if (is.vector(x) && is.numeric(x)) x <- matrix(x, ncol = 1L)
  x <- .cp_data_matrix(x, minimum_rows = 3L, minimum_columns = 1L)
  if (ncol(x) != 1L) stop("The classical API requires univariate `x`.", call. = FALSE)
  gamma <- .cp_scalar(gamma, "gamma", lower = 0, upper = 0.5)
  trim <- .cp_integer(trim, "trim")
  alpha <- .cp_scalar(alpha, "alpha", lower = 0, upper = 1,
                      lower_open = TRUE, upper_open = TRUE)
  strict <- .cp_flag(strict, "strict")
  variance <- match.arg(variance)
  calibration <- match.arg(calibration)
  n <- nrow(x)
  k.grid <- seq.int(trim, n - trim)
  if (!length(k.grid) || min(k.grid) < 1L || max(k.grid) > n - 1L) {
    stop("`trim` leaves no candidate change location.", call. = FALSE)
  }
  if (calibration == "asymptotic" && gamma != 0) {
    stop("Asymptotic Kolmogorov calibration is available only for gamma = 0; use grid-gaussian.",
         call. = FALSE)
  }

  tryCatch({
    core <- ch4_cp_cusum_cpp(x)
    scale <- if (!is.null(sigma)) {
      .cp_scalar(sigma, "sigma", lower = 0, lower_open = TRUE)
    } else if (variance == "difference") {
      sqrt(core$difference_variance[[1L]])
    } else {
      stats::sd(x[, 1L])
    }
    if (!is.finite(scale) || scale <= 0) {
      stop("The selected CUSUM scale estimate is non-positive or non-finite.")
    }
    u <- k.grid / n
    path <- core$cusum0[k.grid, 1L] /
      (scale * (u * (1 - u))^gamma)
    selected <- which.max(abs(path))
    statistic <- abs(path[[selected]])
    k.hat <- k.grid[[selected]]

    calibration.details <- list(type = calibration)
    if (calibration == "asymptotic") {
      p.value <- .cp_kolmogorov_survival(statistic)
      critical <- .cp_kolmogorov_quantile(1 - alpha)
      calibration.details$limit <- "supremum absolute standard Brownian bridge"
    } else {
      covariance <- outer(u, u, function(a, b) pmin(a, b) - a * b)
      weights <- (u * (1 - u))^(-gamma)
      covariance <- covariance * tcrossprod(weights)
      reference <- .cp_with_seed(
        seed,
        .cp_gaussian_max_reference(
          covariance, calibration_draws, absolute = TRUE
        )
      )
      p.value <- .cp_right_tail_mc(statistic, reference$statistic)
      critical <- unname(stats::quantile(
        reference$statistic, 1 - alpha, type = 1
      ))
      calibration.details <- c(
        calibration.details,
        list(
          draws = length(reference$statistic),
          plus.one = TRUE,
          covariance.factor = reference$factor
        )
      )
    }

    structure(
      list(
        statistic = c(maximum.CUSUM = statistic),
        parameter = c(gamma = gamma, trim = trim, alpha = alpha),
        p.value = p.value,
        estimate = c(change.point = k.hat, change.fraction = k.hat / n),
        null.value = c(number_of_changes = 0),
        alternative = "at least one change in univariate location",
        method = "Classical univariate CUSUM change-point test",
        data.name = data.name,
        valid = TRUE,
        path = data.frame(k = k.grid, fraction = u, cusum = path),
        components = list(
          scale = scale,
          scale.source = if (!is.null(sigma)) "user supplied" else variance,
          critical.value = critical,
          reject = p.value <= alpha
        ),
        diagnostics = list(
          calibration = calibration.details,
          no.repair = "No scale floor or variance repair was applied.",
          call = call
        )
      ),
      class = c("chapter4_change_point_test", "htest")
    )
  }, error = function(error) {
    .cp_invalid(
      conditionMessage(error), strict,
      "Classical univariate CUSUM change-point test",
      call, data.name
    )
  })
}


#' Wang--Feng double-max-sum change-point test
#'
#' Implements the feasible DMS procedure of Wang and Feng (2023).  With Rice
#' difference variances \eqn{\widehat\sigma_j^2}, the sparse component is
#' either
#' \deqn{M=\max_{k,j}|C_{0,j}(k)|}
#' or the trimmed weighted version
#' \deqn{M^\dagger=\max_{\lambda\le k\le n-\lambda,j}
#' |C_{1/2,j}(k)|.}
#' The dense component is, literally,
#' \deqn{S=\sum_{k=1}^{n-1}\sum_{j=1}^p C_{1/2,j}(k)^2,}
#' not the book's displayed trimmed sum of `gamma = 0` CUSUMs.  It is centered
#' by `(n + 2) * p` and studentized by the primary paper's leave-four and
#' leave-three finite-difference estimators.
#'
#' The unweighted max pivot is `2 * M^2 - log(2 * p)`.  The weighted pivot is
#' `A(p * log(h)) * M - D(p * log(h))`, where
#' `h = ((lambda / n)^(-1) - 1)^2`.  These correct the nested-log and `h`
#' errors in the book.  `combination = "fisher"` is the primary DMS rule.
#' The optional Cauchy rule is clearly labelled a software extension.
#'
#' The leaveout diagonal is the paper's displayed
#' `A_m = {2,...,n} \\ {i_1,...,i_m}` estimator.  No bridging adjacent
#' difference is silently inserted.  A non-positive component variance or
#' feasible DMS variance is a method failure; it is never repaired by an
#' absolute value or floor.
#'
#' @param x Numeric `n` by `p` data matrix.
#' @param gamma Either `0` or `0.5` for the max component.
#' @param lambda Boundary removal integer required when `gamma = 0.5`.
#' @param combination Primary Fisher rule or a labelled Cauchy extension.
#' @param alpha Test level.
#' @param strict Failure contract; failed calculations return no estimate when
#'   `FALSE` and error when `TRUE`.
#' @return An `htest` object with max, sum, trace, variance, and localization
#'   components.
#' @references Wang, G. and Feng, L. (2023), JRSS B 85, 936--958;
#'   arXiv:2205.00709.
#' @export
#' @examples
#' x <- rbind(matrix(-0.5, 6, 3), matrix(0.5, 6, 3))
#' x <- x + matrix(seq_len(length(x)) %% 5, nrow(x), ncol(x)) / 10
#' wang_feng_dms_test(x, gamma = 0)
wang_feng_dms_test <- function(
    x,
    gamma = 0,
    lambda = NULL,
    combination = c("fisher", "cauchy"),
    alpha = 0.05,
    strict = TRUE) {
  call <- match.call()
  data.name <- deparse(substitute(x))
  x <- .cp_data_matrix(x, minimum_rows = 6L, minimum_columns = 1L)
  gamma <- .cp_scalar(gamma, "gamma", lower = 0, upper = 0.5)
  if (!gamma %in% c(0, 0.5)) stop("`gamma` must be exactly 0 or 0.5.", call. = FALSE)
  combination <- match.arg(combination)
  alpha <- .cp_scalar(alpha, "alpha", lower = 0, upper = 1,
                      lower_open = TRUE, upper_open = TRUE)
  strict <- .cp_flag(strict, "strict")
  n <- nrow(x)
  p <- ncol(x)
  if (gamma == 0.5) {
    if (is.null(lambda)) stop("`lambda` is required when gamma = 0.5.", call. = FALSE)
    lambda <- .cp_integer(lambda, "lambda")
    if (2L * lambda > n) stop("`lambda` leaves no weighted scan location.", call. = FALSE)
  }

  tryCatch({
    cusum <- ch4_cp_cusum_cpp(x)
    s2 <- as.numeric(cusum$difference_variance)
    if (any(!is.finite(s2)) || any(s2 <= 0)) {
      stop("At least one Rice difference variance is non-positive or non-finite.")
    }
    standardized.half <- sweep(cusum$cusum_half, 2L, sqrt(s2), "/")
    dense.by.k <- rowSums(standardized.half^2)
    S <- sum(dense.by.k)

    moments <- ch4_cp_dms_moments_cpp(x)
    if (!isTRUE(moments$valid)) stop(moments$failure)
    trace.hat <- moments$trace_hat
    fourth.hat <- moments$fourth_hat
    variance.hat <- ((2 * pi^2 - 18) / 3) * n^2 * trace.hat +
      ((15 - pi^2) / 3) * n * (fourth.hat - p^2)
    if (!is.finite(variance.hat) || variance.hat <= 0) {
      stop("The literal feasible DMS variance is non-positive or non-finite.")
    }
    sum.z <- (S - (n + 2) * p) / sqrt(variance.hat)
    log.p.sum <- stats::pnorm(sum.z, lower.tail = FALSE, log.p = TRUE)
    p.sum <- exp(log.p.sum)

    if (gamma == 0) {
      max.path <- sweep(cusum$cusum0, 2L, sqrt(s2), "/")
      scan.k <- seq_len(n - 1L)
      pivot <- 2 * max(abs(max.path))^2 - log(2 * p)
    } else {
      max.path <- standardized.half
      scan.k <- seq.int(lambda, n - lambda)
      M.temp <- max(abs(max.path[scan.k, , drop = FALSE]))
      h <- ((lambda / n)^(-1) - 1)^2
      argument <- p * log(h)
      pivot <- .cp_A(argument) * M.temp - .cp_D(argument)
    }
    scan.array <- abs(max.path[scan.k, , drop = FALSE])
    selected <- arrayInd(which.max(scan.array), dim(scan.array))
    max.k <- scan.k[selected[[1L]]]
    max.coordinate <- selected[[2L]]
    M <- scan.array[selected]
    log.p.max <- .cp_gumbel_log_survival(pivot)
    p.max <- exp(log.p.max)

    if (combination == "fisher") {
      statistic <- -2 * (log.p.max + log.p.sum)
      p.value <- stats::pchisq(statistic, df = 4, lower.tail = FALSE)
      combination.details <- list(
        primary = TRUE,
        distribution = "chi-squared with 4 degrees of freedom"
      )
    } else {
      combined <- .cp_cauchy_combine(c(p.max, p.sum))
      statistic <- combined$statistic
      p.value <- combined$p.value
      combination.details <- list(
        primary = FALSE,
        note = "Cauchy combination is a labelled software extension, not the paper's DMS rule."
      )
    }
    dense.k <- which.max(dense.by.k)

    structure(
      list(
        statistic = stats::setNames(statistic, paste0("DMS.", combination)),
        parameter = c(gamma = gamma, alpha = alpha),
        p.value = p.value,
        estimate = c(
          max.change.point = max.k,
          max.change.fraction = max.k / n,
          max.coordinate = max.coordinate,
          dense.change.point = dense.k,
          dense.change.fraction = dense.k / n
        ),
        null.value = c(number_of_changes = 0),
        alternative = "at least one high-dimensional mean change",
        method = "Wang--Feng double-max-sum change-point test",
        data.name = data.name,
        valid = TRUE,
        components = list(
          max = list(
            statistic = M,
            pivot = pivot,
            p.value = p.max,
            log.p.value = log.p.max,
            gamma = gamma,
            lambda = lambda
          ),
          sum = list(
            statistic = S,
            center = (n + 2) * p,
            z = sum.z,
            p.value = p.sum,
            log.p.value = log.p.sum,
            by.location = dense.by.k
          ),
          trace.hat = trace.hat,
          fourth.moment.hat = fourth.hat,
          variance.hat = variance.hat,
          variance.sd = sqrt(variance.hat),
          difference.variance = s2,
          combination = combination.details,
          reject = p.value <= alpha
        ),
        diagnostics = list(
          ordered.or.unordered = paste(
            "Primary DMS moment estimators are consecutive-difference",
            "leaveout averages, not pairwise U-statistics."
          ),
          leaveout.definition = moments$leaveout_definition,
          trace.terms = moments$trace_terms,
          fourth.terms = moments$fourth_terms,
          book.errata = c(
            "Dense S uses gamma=1/2 and k=1,...,n-1.",
            "Max(0) uses 2*M^2-log(2p).",
            "Weighted max uses A(p*log(h)) and h=((lambda/n)^(-1)-1)^2."
          ),
          no.repair = "No variance floor, absolute value, or ridge was applied.",
          call = call
        )
      ),
      class = c("chapter4_change_point_test", "htest")
    )
  }, error = function(error) {
    .cp_invalid(
      conditionMessage(error), strict,
      "Wang--Feng double-max-sum change-point test",
      call, data.name
    )
  })
}


.cp_scaled_hr_fit <- function(x, controls, label) {
  fit <- ch4_cp_scaled_hr_cpp(
    x, controls$tol, controls$max_iter, controls$zero_tol
  )
  if (!isTRUE(fit$valid) || !isTRUE(fit$stable)) {
    failure <- if (is.null(fit$failure)) "unknown scaled HR failure" else
      fit$failure
    stop(sprintf("%s failed: %s", label, failure), call. = FALSE)
  }
  fit$label <- label
  fit
}

.cp_standardized_signs <- function(x, location, diagonal,
                                   zero_tol, zero_is_failure = TRUE,
                                   sqrt_dimension = FALSE) {
  residual <- sweep(x, 2L, location, "-")
  residual <- sweep(residual, 2L, sqrt(diagonal), "/")
  radius <- .cp_row_norms(residual)
  bad <- !is.finite(radius) | radius <= zero_tol
  if (any(bad) && zero_is_failure) {
    stop("A standardized residual has zero or non-finite radius; no perturbation was applied.",
         call. = FALSE)
  }
  signs <- matrix(0, nrow(x), ncol(x))
  if (any(!bad)) {
    signs[!bad, ] <- sweep(residual[!bad, , drop = FALSE], 1L,
                           radius[!bad], "/")
  }
  if (sqrt_dimension) signs <- sqrt(ncol(x)) * signs
  list(signs = signs, radius = radius, zero = which(bad))
}

.cp_endpoint_trace <- function(x, endpoint.indices, endpoint.fit,
                               common.diagonal, zero_tol) {
  signs <- .cp_standardized_signs(
    x[endpoint.indices, , drop = FALSE],
    endpoint.fit$location,
    common.diagonal,
    zero_tol,
    zero_is_failure = TRUE
  )$signs
  ch4_cp_ordered_pair_square_sum_cpp(signs)
}

.cp_spatial_fit_summary <- function(fits) {
  iterations <- vapply(fits, function(fit) fit$iterations, integer(1))
  updates <- vapply(fits, function(fit) fit$relative_update, numeric(1))
  score.l2 <- vapply(fits, function(fit) fit$score_l2, numeric(1))
  score.inf <- vapply(fits, function(fit) fit$score_infinity, numeric(1))
  diagonal <- vapply(fits, function(fit) fit$diagonal_residual, numeric(1))
  minimum <- vapply(fits, function(fit) fit$minimum_distance, numeric(1))
  worst <- which.max(updates)
  list(
    all.stable = all(vapply(fits, function(fit) isTRUE(fit$stable), logical(1))),
    fit.count = length(fits),
    maximum.iterations = max(iterations),
    maximum.relative.update = max(updates),
    maximum.score.residual.l2 = max(score.l2),
    maximum.score.residual.infinity = max(score.inf),
    maximum.diagonal.residual = max(diagonal),
    minimum.residual.distance = min(minimum),
    worst.label = fits[[worst]]$label,
    convergence.basis = "relative location and log-diagonal update",
    scale.identification = "only relative diagonal scale is identified"
  )
}


#' Spatial-sign max-Linf, max-L2, and adaptive change-point test
#'
#' Implements the feasible procedure of Liu, Feng, Peng, and Wang (2025).
#' Endpoint blocks of size `floor(n * endpoint_fraction)` supply two joint
#' scaled-spatial-median/diagonal fits.  Their diagonals are normalized by
#' their first entries and averaged,
#' \deqn{\widehat D=\{\widehat D_1/\widehat d_{1,1}^2+
#' \widehat D_2/\widehat d_{2,1}^2\}/2,}
#' exactly as in the primary paper.  The same stable endpoint blocks estimate
#' \eqn{\zeta_1} and \eqn{\operatorname{tr}(R^2)}.  This differs from the
#' book's full-sample shortcut.
#'
#' For every candidate `k`, separate joint fits give the two segment location
#' estimates.  With the common endpoint diagonal, the max-Linf statistic is
#' based on
#' \deqn{C_\gamma(k)=\{u(1-u)\}^{1-\gamma}\sqrt n\,
#' \widehat D^{-1/2}(\widehat\theta_{1:k}-
#' \widehat\theta_{k+1:n}).}
#' The max-L2 statistic uses the full-sample sign partial sums.  For
#' `variant = "unweighted"`, the primary normalization is
#' `S / sqrt(2 * trace_hat)`, not the book's `p*S/(2*trace_hat)`.
#' Its Gaussian-process CDF is approximated on the actual trimmed scan grid.
#' For `variant = "weighted"`, the primary pivot uses
#' `abs(S_dagger) / sqrt(2 * trace_hat)`; the absolute value and square root
#' are both essential.
#'
#' Endpoint trace sums are ordered pairs: each block contributes
#' \deqn{\frac{p^2}{2m(m-1)}\sum_{i\ne j}(U_i^TU_j)^2.}
#' The adaptive p-value is the primary Fisher combination.  No paper size or
#' power simulation is reproduced; the Gaussian-process draws are intrinsic
#' null calibration only.
#'
#' @param x Numeric `n` by `p` matrix.
#' @param lambda Integer boundary removal parameter.
#' @param endpoint_fraction Stable endpoint proportion in `(0, 1/2)`.
#' @param variant Unweighted or weighted primary pair of max-Linf/max-L2 tests.
#' @param combination Fisher adaptive test, or `"none"` to report components
#'   with the max-Linf p-value as the `htest` p-value.
#' @param fv_draws Intrinsic Gaussian-process draws for the unweighted max-L2
#'   null CDF.
#' @param alpha Test level.
#' @param seed Optional calibration seed; previous RNG state is restored.
#' @param tol,max_iter,zero_tol Joint fixed-point controls.  The default
#'   `zero_tol = 0` follows the exact formula; positive values are explicit
#'   user-requested failure tolerances.
#' @param keep_reference Whether to retain Gaussian-process maxima.
#' @param strict Failure contract; no failed segment fit is silently dropped.
#' @return An `htest` object containing both primary component statistics,
#'   feasible radial/trace quantities, and worst-iteration diagnostics.
#' @references Liu, J., Feng, L., Peng, L. and Wang, Z. (2025),
#'   arXiv:2504.19306.
#' @export
#' @examples
#' x <- rbind(
#'   c(-2, 0), c(2, 0), c(0, -2), c(0, 2),
#'   c(-1, -1), c(1, 1), c(-1, 1), c(1, -1),
#'   c(-2, 1), c(2, -1), c(-1, 2), c(1, -2)
#' )
#' spatial_sign_change_point_test(
#'   x, lambda = 3, endpoint_fraction = 0.25,
#'   fv_draws = 99, seed = 1
#' )
spatial_sign_change_point_test <- function(
    x,
    lambda,
    endpoint_fraction = 0.2,
    variant = c("unweighted", "weighted"),
    combination = c("fisher", "none"),
    fv_draws = 4999L,
    alpha = 0.05,
    seed = NULL,
    tol = 1e-8,
    max_iter = 1000L,
    zero_tol = 0,
    keep_reference = FALSE,
    strict = TRUE) {
  call <- match.call()
  data.name <- deparse(substitute(x))
  x <- .cp_data_matrix(x, minimum_rows = 8L, minimum_columns = 2L)
  n <- nrow(x)
  p <- ncol(x)
  lambda <- .cp_integer(lambda, "lambda", minimum = 2L)
  if (2L * lambda > n) stop("`lambda` leaves no spatial-sign scan location.", call. = FALSE)
  endpoint_fraction <- .cp_scalar(
    endpoint_fraction, "endpoint_fraction", lower = 0, upper = 0.5,
    lower_open = TRUE, upper_open = TRUE
  )
  variant <- match.arg(variant)
  combination <- match.arg(combination)
  fv_draws <- .cp_integer(fv_draws, "fv_draws")
  alpha <- .cp_scalar(alpha, "alpha", lower = 0, upper = 1,
                      lower_open = TRUE, upper_open = TRUE)
  controls <- .cp_iteration_controls(tol, max_iter, zero_tol)
  keep_reference <- .cp_flag(keep_reference, "keep_reference")
  strict <- .cp_flag(strict, "strict")
  m <- floor(n * endpoint_fraction)
  if (m < 2L || 2L * m > n) {
    stop("`endpoint_fraction` must produce two disjoint blocks of at least two rows.",
         call. = FALSE)
  }
  k.grid <- seq.int(lambda, n - lambda)

  tryCatch({
    first.index <- seq_len(m)
    last.index <- seq.int(n - m + 1L, n)
    first.fit <- .cp_scaled_hr_fit(
      x[first.index, , drop = FALSE], controls, "first endpoint fit"
    )
    last.fit <- .cp_scaled_hr_fit(
      x[last.index, , drop = FALSE], controls, "last endpoint fit"
    )
    full.fit <- .cp_scaled_hr_fit(x, controls, "full-sample sign fit")

    first.normalized <- first.fit$diagonal / first.fit$diagonal[[1L]]
    last.normalized <- last.fit$diagonal / last.fit$diagonal[[1L]]
    common.diagonal <- (first.normalized + last.normalized) / 2
    if (any(!is.finite(common.diagonal)) || any(common.diagonal <= 0)) {
      stop("The endpoint-normalized common diagonal is invalid.")
    }

    segment.fits <- vector("list", 2L * length(k.grid))
    C.path <- matrix(NA_real_, length(k.grid), p)
    gamma <- if (variant == "unweighted") 0 else 0.5
    fit.id <- 1L
    for (ii in seq_along(k.grid)) {
      k <- k.grid[[ii]]
      left <- .cp_scaled_hr_fit(
        x[seq_len(k), , drop = FALSE], controls,
        sprintf("left segment k=%d", k)
      )
      right <- .cp_scaled_hr_fit(
        x[seq.int(k + 1L, n), , drop = FALSE], controls,
        sprintf("right segment k=%d", k)
      )
      segment.fits[[fit.id]] <- left
      segment.fits[[fit.id + 1L]] <- right
      fit.id <- fit.id + 2L
      u <- k / n
      C.path[ii, ] <- (u * (1 - u))^(1 - gamma) * sqrt(n) *
        (left$location - right$location) / sqrt(common.diagonal)
    }
    finite.correction <- 1 - n^(-0.5)
    row.max <- apply(abs(C.path), 1L, max)
    max.row <- which.max(row.max)
    max.coordinate <- which.max(abs(C.path[max.row, ]))
    M <- finite.correction * row.max[[max.row]]

    first.common <- .cp_standardized_signs(
      x[first.index, , drop = FALSE], first.fit$location,
      common.diagonal, controls$zero_tol, zero_is_failure = TRUE
    )
    last.common <- .cp_standardized_signs(
      x[last.index, , drop = FALSE], last.fit$location,
      common.diagonal, controls$zero_tol, zero_is_failure = TRUE
    )
    zeta.hat <- 0.5 * (mean(1 / first.common$radius) +
                        mean(1 / last.common$radius))
    if (!is.finite(zeta.hat) || zeta.hat <= 0) {
      stop("The endpoint inverse-radius estimate is invalid.")
    }

    trace.first.sum <- ch4_cp_ordered_pair_square_sum_cpp(first.common$signs)
    trace.last.sum <- ch4_cp_ordered_pair_square_sum_cpp(last.common$signs)
    trace.hat <- p^2 * (trace.first.sum + trace.last.sum) /
      (2 * m * (m - 1))
    if (!is.finite(trace.hat) || trace.hat <= 0) {
      stop("The ordered-pair endpoint trace estimate is non-positive or non-finite.")
    }

    full.common <- .cp_standardized_signs(
      x, full.fit$location, common.diagonal, controls$zero_tol,
      zero_is_failure = TRUE
    )
    signs <- full.common$signs
    partial <- apply(signs, 2L, cumsum)
    if (p == 1L) partial <- matrix(partial, ncol = 1L)
    total <- colSums(signs)
    l2.path <- numeric(length(k.grid))
    for (ii in seq_along(k.grid)) {
      k <- k.grid[[ii]]
      u <- k / n
      transformed <- sqrt(p / n) *
        (partial[k, ] - u * total) / (u * (1 - u))^gamma
      center <- if (variant == "unweighted") u * (1 - u) * p else p
      l2.path[[ii]] <- sum(transformed^2) - center
    }
    l2.row <- which.max(l2.path)
    S <- finite.correction * l2.path[[l2.row]]

    if (variant == "unweighted") {
      max.pivot <- 2 * p * zeta.hat^2 * M^2 - log(2 * p)
      p.max <- .cp_gumbel_survival(max.pivot)
      u.grid <- k.grid / n
      covariance <- outer(
        u.grid, u.grid,
        function(a, b) (pmin(a, b) - a * b)^2
      )
      reference <- .cp_with_seed(
        seed,
        .cp_gaussian_max_reference(covariance, fv_draws, absolute = FALSE)
      )
      l2.pivot <- S / sqrt(2 * trace.hat)
      p.l2 <- .cp_right_tail_mc(l2.pivot, reference$statistic)
      l2.critical <- unname(stats::quantile(
        reference$statistic, 1 - alpha, type = 1
      ))
      calibration <- list(
        distribution = "max of V(t) on the actual trimmed grid",
        covariance = "(min(s,t)-s*t)^2",
        draws = length(reference$statistic),
        plus.one = TRUE,
        factor = reference$factor,
        critical.value = l2.critical,
        reference = if (keep_reference) reference$statistic else NULL
      )
    } else {
      h <- ((lambda / n)^(-1) - 1)^2
      max.argument <- p * log(h)
      max.pivot <- sqrt(p) * zeta.hat * .cp_A(max.argument) * M -
        .cp_D(max.argument)
      p.max <- .cp_gumbel_survival(max.pivot)
      l2.argument <- log(n^2 / lambda^2)
      l2.pivot <- .cp_A(l2.argument) * abs(S) / sqrt(2 * trace.hat) -
        .cp_D(l2.argument)
      p.l2 <- .cp_gumbel_survival(l2.pivot, factor = 2)
      l2.critical <- (.cp_D(l2.argument) -
        log(-0.5 * log(1 - alpha))) / .cp_A(l2.argument)
      calibration <- list(
        distribution = "G2(x)=exp(-2 exp(-x))",
        argument = l2.argument,
        absolute.statistic = TRUE,
        critical.pivot = l2.critical
      )
    }

    if (combination == "fisher") {
      log.max <- if (p.max == 0) -Inf else log(p.max)
      log.l2 <- if (p.l2 == 0) -Inf else log(p.l2)
      adaptive.statistic <- -2 * (log.max + log.l2)
      p.value <- stats::pchisq(
        adaptive.statistic, df = 4, lower.tail = FALSE
      )
    } else {
      adaptive.statistic <- max.pivot
      p.value <- p.max
    }
    all.fits <- c(list(first.fit, last.fit, full.fit), segment.fits)
    fit.summary <- .cp_spatial_fit_summary(all.fits)

    structure(
      list(
        statistic = stats::setNames(
          adaptive.statistic,
          if (combination == "fisher") "adaptive.Fisher" else "max.Linf.pivot"
        ),
        parameter = c(
          lambda = lambda,
          endpoint_fraction = endpoint_fraction,
          alpha = alpha
        ),
        p.value = p.value,
        estimate = c(
          max.change.point = k.grid[[max.row]],
          max.change.fraction = k.grid[[max.row]] / n,
          max.coordinate = max.coordinate,
          l2.change.point = k.grid[[l2.row]],
          l2.change.fraction = k.grid[[l2.row]] / n
        ),
        null.value = c(number_of_changes = 0),
        alternative = "at least one elliptical location change",
        method = paste(
          "Liu--Feng--Peng--Wang spatial-sign",
          variant, "change-point test"
        ),
        data.name = data.name,
        valid = TRUE,
        components = list(
          max.Linf = list(
            statistic = M,
            pivot = max.pivot,
            p.value = p.max,
            path = C.path,
            gamma = gamma
          ),
          max.L2 = list(
            statistic = S,
            pivot = l2.pivot,
            p.value = p.l2,
            path = l2.path,
            gamma = gamma,
            calibration = calibration
          ),
          zeta1.hat = zeta.hat,
          trace.R2.hat = trace.hat,
          trace.ordered.sums = c(
            first = trace.first.sum,
            last = trace.last.sum
          ),
          common.diagonal = common.diagonal,
          common.diagonal.normalization =
            "average of endpoint diagonals after division by first entry",
          estimation = list(
            first.endpoint.location = first.fit$location,
            last.endpoint.location = last.fit$location,
            full.sample.location = full.fit$location,
            first.endpoint.diagonal.raw = first.fit$diagonal,
            last.endpoint.diagonal.raw = last.fit$diagonal,
            first.endpoint.signs = first.common$signs,
            last.endpoint.signs = last.common$signs,
            full.sample.signs = signs
          ),
          finite.sample.factor = finite.correction,
          reject = p.value <= alpha
        ),
        diagnostics = list(
          iteration = fit.summary,
          endpoint.block.size = m,
          segment.fit.count = length(segment.fits),
          scale.identification = paste(
            "The joint equations identify only relative diagonal scale;",
            "endpoint first-coordinate normalization fixes the reported common diagonal."
          ),
          trace.denominator = "2*m*(m-1) for each ordered-pair endpoint contribution",
          book.errata = c(
            "Feasible D and zeta1 use stable endpoint blocks, not the full sample.",
            "Unweighted max-L2 scale is sqrt(2*tr(R^2)), with no leading p.",
            "Weighted max-L2 uses abs(S_dagger)/sqrt(2*tr(R^2))."
          ),
          no.repair = paste(
            "No zero-radius perturbation, diagonal ridge, variance floor,",
            "absolute-value variance repair, or dropped failed split was used."
          ),
          call = call
        )
      ),
      class = c("chapter4_change_point_test", "htest")
    )
  }, error = function(error) {
    .cp_invalid(
      conditionMessage(error), strict,
      paste("Liu--Feng--Peng--Wang spatial-sign", variant,
            "change-point test"),
      call, data.name,
      diagnostics = list(
        endpoint.block.size = m,
        requested.splits = k.grid
      )
    )
  })
}


.cp_spatial_median_fit <- function(x, controls, label) {
  fit <- ch4_cp_spatial_median_cpp(
    x, controls$tol, controls$max_iter, controls$zero_tol
  )
  if (!isTRUE(fit$valid) || !isTRUE(fit$stable)) {
    failure <- if (is.null(fit$failure)) "unknown spatial-median failure" else
      fit$failure
    stop(sprintf("%s failed: %s", label, failure), call. = FALSE)
  }
  fit$label <- label
  fit
}

.cp_erht_signs <- function(x, location, zero_tol) {
  residual <- sweep(x, 2L, location, "-")
  radius <- .cp_row_norms(residual)
  bad <- !is.finite(radius)
  if (any(bad)) stop("An ERHT residual radius is non-finite.", call. = FALSE)
  signs <- matrix(0, nrow(x), ncol(x))
  nonzero <- radius > zero_tol
  if (any(nonzero)) {
    signs[nonzero, ] <- sqrt(ncol(x)) *
      sweep(residual[nonzero, , drop = FALSE], 1L, radius[nonzero], "/")
  }
  list(signs = signs, radius = radius, coincident = sum(!nonzero))
}

.cp_erht_inverse_distance <- function(x, location, zero_tol) {
  radius <- .cp_row_norms(sweep(x, 2L, location, "-"))
  if (any(!is.finite(radius))) {
    stop("An ERHT segment radius is non-finite.", call. = FALSE)
  }
  # The primary paper explicitly defines a coincident observation's weight as
  # zero.  This is a method convention, not a numerical floor.
  weight <- numeric(length(radius))
  nonzero <- radius > zero_tol
  weight[nonzero] <- sqrt(ncol(x)) / radius[nonzero]
  value <- mean(weight)
  if (!is.finite(value) || value <= 0) {
    stop("An ERHT inverse-distance average is non-positive or non-finite.",
         call. = FALSE)
  }
  list(value = value, radius = radius, coincident = sum(!nonzero))
}

.cp_erht_candidates <- function(n, epsilon, scan) {
  if (scan == "single") {
    k <- seq.int(ceiling(epsilon * n), floor((1 - epsilon) * n))
    k <- k[k >= 1L & k <= n - 1L]
    if (!length(k)) stop("`epsilon` leaves no single-change split.", call. = FALSE)
    return(data.frame(
      q1 = 0L, q2 = k, q3 = n,
      t1 = 0, t2 = k / n, t3 = 1
    ))
  }

  grid <- seq(0, 1, by = epsilon)
  if (utils::tail(grid, 1L) < 1 - 100 * .Machine$double.eps) grid <- c(grid, 1)
  grid <- unique(c(grid, 1))
  rows <- list()
  id <- 1L
  for (a in seq_len(length(grid) - 2L)) {
    for (b in seq.int(a + 1L, length(grid) - 1L)) {
      for (cc in seq.int(b + 1L, length(grid))) {
        if (grid[[b]] - grid[[a]] + 1e-12 < epsilon ||
            grid[[cc]] - grid[[b]] + 1e-12 < epsilon) next
        q <- floor(n * grid[base::c(a, b, cc)] + 1e-12)
        if (grid[[a]] == 0) q[[1L]] <- 0L
        if (grid[[cc]] == 1) q[[3L]] <- n
        if (q[[2L]] <= q[[1L]] || q[[3L]] <= q[[2L]]) next
        rows[[id]] <- data.frame(
          q1 = q[[1L]], q2 = q[[2L]], q3 = q[[3L]],
          t1 = q[[1L]] / n, t2 = q[[2L]] / n,
          t3 = q[[3L]] / n
        )
        id <- id + 1L
      }
    }
  }
  if (!length(rows)) stop("`epsilon` leaves no adjacent-triple candidate.", call. = FALSE)
  unique(do.call(rbind, rows))
}

.cp_interval_overlap <- function(a, b, c, d) {
  max(0, min(b, d) - max(a, c))
}

.cp_erht_temporal_inner <- function(left, right) {
  l1 <- left$t2 - left$t1
  l2 <- left$t3 - left$t2
  r1 <- right$t2 - right$t1
  r2 <- right$t3 - right$t2
  nl <- l1 * l2 / (l1 + l2)
  nr <- r1 * r2 / (r1 + r2)
  wl <- c(-sqrt(nl) / l1, sqrt(nl) / l2)
  wr <- c(-sqrt(nr) / r1, sqrt(nr) / r2)
  il <- matrix(c(left$t1, left$t2, left$t2, left$t3), 2L, byrow = TRUE)
  ir <- matrix(c(right$t1, right$t2, right$t2, right$t3), 2L, byrow = TRUE)
  value <- 0
  for (i in 1:2) for (j in 1:2) {
    value <- value + wl[[i]] * wr[[j]] *
      .cp_interval_overlap(il[i, 1L], il[i, 2L], ir[j, 1L], ir[j, 2L])
  }
  value
}

.cp_erht_temporal_covariance <- function(candidates) {
  count <- nrow(candidates)
  inner <- matrix(0, count, count)
  for (i in seq_len(count)) {
    for (j in i:count) {
      value <- .cp_erht_temporal_inner(candidates[i, ], candidates[j, ])
      inner[i, j] <- value
      inner[j, i] <- value
    }
  }
  diagonal <- diag(inner)
  if (any(!is.finite(diagonal)) || any(diagonal <= 0)) {
    stop("An ERHT temporal contrast has non-positive self inner product.",
         call. = FALSE)
  }
  covariance <- inner^2 / sqrt(outer(diagonal^2, diagonal^2))
  diag(covariance) <- 1
  covariance
}

.cp_erht_segment_cache <- function(x, candidates, controls) {
  keys <- unique(c(
    paste(candidates$q1, candidates$q2, sep = ":"),
    paste(candidates$q2, candidates$q3, sep = ":")
  ))
  cache <- new.env(parent = emptyenv())
  fit.list <- vector("list", length(keys))
  for (ii in seq_along(keys)) {
    limits <- as.integer(strsplit(keys[[ii]], ":", fixed = TRUE)[[1L]])
    index <- seq.int(limits[[1L]] + 1L, limits[[2L]])
    segment <- x[index, , drop = FALSE]
    fit <- .cp_spatial_median_fit(
      segment, controls, paste0("ERHT segment ", keys[[ii]])
    )
    inverse <- .cp_erht_inverse_distance(
      segment, fit$location, controls$zero_tol
    )
    value <- list(
      location = fit$location,
      e = inverse$value,
      n = length(index),
      fit = fit,
      coincident.weights = inverse$coincident
    )
    assign(keys[[ii]], value, envir = cache)
    fit.list[[ii]] <- fit
  }
  list(cache = cache, fits = fit.list)
}

.cp_erht_checked_inverse <- function(matrix, label) {
  chol.factor <- tryCatch(chol(matrix), error = identity)
  if (inherits(chol.factor, "error")) {
    stop(sprintf("%s is not numerically positive definite; no repair was applied.",
                 label), call. = FALSE)
  }
  inverse <- chol2inv(chol.factor)
  if (any(!is.finite(inverse))) {
    stop(sprintf("%s inverse is non-finite.", label), call. = FALSE)
  }
  inverse
}

.cp_erht_core <- function(x, ridge, candidates, controls) {
  n <- nrow(x)
  p <- ncol(x)
  pool.fit <- .cp_spatial_median_fit(x, controls, "ERHT common-pool median")
  pool.sign <- .cp_erht_signs(x, pool.fit$location, controls$zero_tol)
  Y <- pool.sign$signs
  R.hat <- crossprod(Y) / n
  segment <- .cp_erht_segment_cache(x, candidates, controls)
  count <- nrow(candidates)
  delta <- matrix(NA_real_, count, p)
  N.effective <- numeric(count)
  beta <- matrix(0, count, n)
  coincident <- integer(count)

  for (ii in seq_len(count)) {
    key1 <- paste(candidates$q1[[ii]], candidates$q2[[ii]], sep = ":")
    key2 <- paste(candidates$q2[[ii]], candidates$q3[[ii]], sep = ":")
    first <- get(key1, envir = segment$cache)
    second <- get(key2, envir = segment$cache)
    n1 <- first$n
    n2 <- second$n
    N.effective[[ii]] <- n1 * n2 / (n1 + n2)
    delta[ii, ] <- second$location - first$location
    root.N <- sqrt(N.effective[[ii]])
    first.index <- seq.int(candidates$q1[[ii]] + 1L,
                           candidates$q2[[ii]])
    second.index <- seq.int(candidates$q2[[ii]] + 1L,
                            candidates$q3[[ii]])
    beta[ii, first.index] <- -root.N / (n1 * first$e)
    beta[ii, second.index] <- root.N / (n2 * second$e)
    coincident[[ii]] <- first$coincident.weights + second$coincident.weights
  }

  z <- matrix(NA_real_, count, length(ridge))
  raw <- center <- variance <- matrix(NA_real_, count, length(ridge))
  kappa <- sigma2 <- matrix(NA_real_, count, length(ridge))
  companion.rcond <- numeric(length(ridge))
  Q.list <- vector("list", length(ridge))

  for (rr in seq_along(ridge)) {
    Q <- .cp_erht_checked_inverse(
      R.hat + diag(ridge[[rr]], p),
      sprintf("ERHT ridge matrix for rho=%g", ridge[[rr]])
    )
    Q.list[[rr]] <- Q
    companion <- (Y %*% Q %*% t(Y)) / n
    companion.rcond[[rr]] <- rcond(R.hat + diag(ridge[[rr]], p))
    for (ii in seq_len(count)) {
      raw[ii, rr] <- N.effective[[ii]] *
        drop(delta[ii, , drop = FALSE] %*% Q %*%
               matrix(delta[ii, ], ncol = 1L))
      moments <- ch4_cp_erht_moments_cpp(companion, beta[ii, ])
      if (!isTRUE(moments$valid)) {
        stop(sprintf(
          "ERHT candidate %d at ridge %g failed: %s",
          ii, ridge[[rr]], moments$failure
        ), call. = FALSE)
      }
      kappa[ii, rr] <- moments$kappa
      sigma2[ii, rr] <- moments$sigma2
      center[ii, rr] <- n * moments$kappa
      variance[ii, rr] <- n * moments$sigma2
      z[ii, rr] <- (raw[ii, rr] - center[ii, rr]) /
        sqrt(variance[ii, rr])
    }
  }
  if (any(!is.finite(z))) stop("At least one ERHT local statistic is non-finite.", call. = FALSE)

  fit.summary <- list(
    all.stable = TRUE,
    fit.count = 1L + length(segment$fits),
    maximum.iterations = max(c(pool.fit$iterations,
      vapply(segment$fits, function(fit) fit$iterations, integer(1)))),
    maximum.relative.update = max(c(pool.fit$relative_update,
      vapply(segment$fits, function(fit) fit$relative_update, numeric(1)))),
    maximum.score.residual = max(c(pool.fit$score_residual,
      vapply(segment$fits, function(fit) fit$score_residual, numeric(1)))),
    minimum.residual.distance = min(c(pool.fit$minimum_distance,
      vapply(segment$fits, function(fit) fit$minimum_distance, numeric(1)))),
    coincident.pool.signs = pool.sign$coincident,
    coincident.segment.weights = sum(coincident),
    convergence.basis = "modified Weiszfeld relative update with exact subgradient certificate"
  )

  list(
    z = z,
    raw = raw,
    center = center,
    variance = variance,
    kappa = kappa,
    sigma2 = sigma2,
    delta = delta,
    beta = beta,
    N.effective = N.effective,
    candidates = candidates,
    ridge = ridge,
    pool.location = pool.fit$location,
    pool.signs = Y,
    R.hat = R.hat,
    Q = Q.list,
    companion.rcond = companion.rcond,
    fit.summary = fit.summary,
    jacobian.convention =
      "e=mean(sqrt(p)/distance), with coincident weight zero; no (p-1)/p factor"
  )
}

.cp_erht_statistics <- function(core) {
  statistic <- apply(core$z, 2L, max)
  selected <- vapply(seq_along(statistic), function(j) {
    which.max(core$z[, j])
  }, integer(1))
  list(statistic = statistic, selected = selected)
}

.cp_erht_null_reference <- function(candidates, draws) {
  covariance <- .cp_erht_temporal_covariance(candidates)
  reference <- .cp_gaussian_max_reference(
    covariance, draws, absolute = FALSE
  )
  c(reference, list(covariance = covariance))
}


#' Elliptical regularized Hotelling change-point scan
#'
#' Implements the single-change path or the discretized adjacent-triple scan
#' of Song, Wen, and Feng (2026).  For adjacent segments with sizes `n1`, `n2`,
#' the raw statistic is
#' \deqn{V_\rho^{raw}=N\widehat\Delta^T
#' (\widehat R+\rho I)^{-1}\widehat\Delta,
#' \qquad N=n_1n_2/(n_1+n_2).}
#' The common pool is the full sample.  Its centered spatial signs are columns
#' of \eqn{Y}; with \eqn{A=Y^TQY/n}, the primary studentization is
#' \deqn{\widehat\kappa=\sum_i\beta_i^2A_{ii},\qquad
#' \widehat\sigma^2=2n\sum_{i\ne j}\beta_i^2\beta_j^2A_{ij}^2,}
#' \deqn{Z_\rho=(V_\rho^{raw}-n\widehat\kappa)/
#' \{n\widehat\sigma^2\}^{1/2}.}
#' The ordered off-diagonal variance is accumulated directly.  It is never
#' obtained by flooring a cancellation-prone difference.
#'
#' `ridge` contains the actual positive \eqn{\rho} values in the primary
#' formula.  The public ERHTCP reproduction code instead accepts empirical
#' ratios `rho / (p/n)`; this function does not silently make that conversion.
#' Users who want that convention can supply `ridge_ratio * p / n` explicitly,
#' and all resolved ridge values are returned.
#'
#' `calibration = "gaussian-supremum"` simulates the primary marginal
#' Gaussian-process limit on the exact candidate grid.  `"time-permutation"`
#' is the paper's practical exchangeability calibration and preserves each
#' multivariate row.  It is invalid under unaddressed serial dependence.
#' With multiple ridge values, the returned Cauchy transform combines the
#' marginal p-values, but `calibration.exact` is `FALSE`: the paper explicitly
#' distinguishes this analytic rule from exact calibration by a joint-limit
#' quantile involving the unknown cross-ridge correlation `r_E`.
#'
#' The primary inverse-distance average is exactly
#' `mean(sqrt(p) / distance)`, with a coincident observation assigned weight
#' zero.  The optional `(p - 1) / p` factor used by the public reproduction
#' repository is not part of the displayed primary statistic and is not used.
#'
#' @param x Numeric `n` by `p` data matrix.
#' @param ridge Positive actual ridge value or finite grid.
#' @param scan `"single"` or the primary discretized `"multiple"`
#'   adjacent-triple scan.
#' @param epsilon Trimming/minimum-segment fraction in `(0, 1/2)`; for
#'   `scan = "multiple"` it is also the primary grid spacing.
#' @param calibration Marginal Gaussian-supremum, time-permutation, or no
#'   p-value calibration.
#' @param calibration_draws Number of intrinsic null draws/permutations.
#' @param alpha Test level.
#' @param seed Optional calibration seed; previous RNG state is restored.
#' @param tol,max_iter,zero_tol Spatial-median controls.
#' @param cauchy_weights Optional positive ridge-combination weights.
#' @param keep_calibration Whether to retain simulated marginal null maxima.
#' @param strict Failure contract.  No invalid candidate is silently removed.
#' @return An `htest` object with every local raw/center/variance/Z component.
#' @references Song, Wen and Feng (2026), arXiv:2607.22162.
#' @export
#' @examples
#' x <- rbind(
#'   c(-2, 0), c(2, 0), c(0, -2), c(0, 2),
#'   c(-1, -1), c(1, 1), c(-1, 1), c(1, -1),
#'   c(-2, 1), c(2, -1), c(-1, 2), c(1, -2)
#' )
#' erht_change_point_test(
#'   x, ridge = 0.5, epsilon = 0.25,
#'   calibration_draws = 99, seed = 1
#' )
erht_change_point_test <- function(
    x,
    ridge,
    scan = c("single", "multiple"),
    epsilon = 0.1,
    calibration = c("gaussian-supremum", "time-permutation", "none"),
    calibration_draws = 4999L,
    alpha = 0.05,
    seed = NULL,
    tol = 1e-8,
    max_iter = 1000L,
    zero_tol = 0,
    cauchy_weights = NULL,
    keep_calibration = FALSE,
    strict = TRUE) {
  call <- match.call()
  data.name <- deparse(substitute(x))
  x <- .cp_data_matrix(x, minimum_rows = 6L, minimum_columns = 2L)
  ridge <- as.numeric(ridge)
  if (!length(ridge) || anyNA(ridge) || any(!is.finite(ridge)) ||
      any(ridge <= 0)) {
    stop("`ridge` must contain positive finite actual rho values.",
         call. = FALSE)
  }
  ridge <- unique(ridge)
  scan <- match.arg(scan)
  epsilon <- .cp_scalar(
    epsilon, "epsilon", lower = 0, upper = 0.5,
    lower_open = TRUE, upper_open = TRUE
  )
  calibration <- match.arg(calibration)
  calibration_draws <- .cp_integer(
    calibration_draws, "calibration_draws"
  )
  alpha <- .cp_scalar(alpha, "alpha", lower = 0, upper = 1,
                      lower_open = TRUE, upper_open = TRUE)
  controls <- .cp_iteration_controls(tol, max_iter, zero_tol)
  keep_calibration <- .cp_flag(keep_calibration, "keep_calibration")
  strict <- .cp_flag(strict, "strict")
  candidates <- .cp_erht_candidates(nrow(x), epsilon, scan)
  if (!is.null(cauchy_weights)) {
    cauchy_weights <- as.numeric(cauchy_weights)
    if (length(cauchy_weights) != length(ridge) ||
        anyNA(cauchy_weights) || any(!is.finite(cauchy_weights)) ||
        any(cauchy_weights <= 0)) {
      stop("`cauchy_weights` must be positive, finite, and match `ridge`.",
           call. = FALSE)
    }
    cauchy_weights <- cauchy_weights / sum(cauchy_weights)
  }

  tryCatch({
    observed <- .cp_erht_core(x, ridge, candidates, controls)
    observed.summary <- .cp_erht_statistics(observed)
    reference <- NULL
    p.ridge <- rep(NA_real_, length(ridge))
    critical <- rep(NA_real_, length(ridge))
    calibration.details <- list(type = calibration)

    if (calibration == "gaussian-supremum") {
      reference <- .cp_with_seed(
        seed,
        .cp_erht_null_reference(candidates, calibration_draws)
      )
      p.ridge <- vapply(
        observed.summary$statistic,
        .cp_right_tail_mc,
        numeric(1),
        reference = reference$statistic
      )
      critical[] <- unname(stats::quantile(
        reference$statistic, 1 - alpha, type = 1
      ))
      calibration.details <- c(
        calibration.details,
        list(
          draws = length(reference$statistic),
          plus.one = TRUE,
          temporal.covariance = reference$covariance,
          covariance.factor = reference$factor,
          exchangeability.required = FALSE
        )
      )
    } else if (calibration == "time-permutation") {
      permutation.statistic <- .cp_with_seed(seed, {
        answer <- matrix(NA_real_, calibration_draws, length(ridge))
        for (bb in seq_len(calibration_draws)) {
          permuted <- x[sample.int(nrow(x)), , drop = FALSE]
          permuted.core <- .cp_erht_core(
            permuted, ridge, candidates, controls
          )
          answer[bb, ] <- .cp_erht_statistics(permuted.core)$statistic
        }
        answer
      })
      p.ridge <- vapply(seq_along(ridge), function(j) {
        .cp_right_tail_mc(
          observed.summary$statistic[[j]], permutation.statistic[, j]
        )
      }, numeric(1))
      critical <- apply(
        permutation.statistic, 2L, stats::quantile,
        probs = 1 - alpha, type = 1, names = FALSE
      )
      reference <- list(statistic = permutation.statistic)
      calibration.details <- c(
        calibration.details,
        list(
          permutations = calibration_draws,
          plus.one = TRUE,
          exchangeability.required = TRUE,
          warning = paste(
            "Whole rows are permuted. This practical primary-paper calibration",
            "does not account for serial dependence."
          )
        )
      )
    }

    if (calibration == "none") {
      final.statistic <- max(observed.summary$statistic)
      p.value <- NA_real_
      combined <- NULL
    } else if (length(ridge) == 1L) {
      final.statistic <- observed.summary$statistic[[1L]]
      p.value <- p.ridge[[1L]]
      combined <- list(
        statistic = final.statistic,
        p.value = p.value,
        calibration.exact = TRUE,
        exact.scope = if (calibration == "gaussian-supremum")
          "marginal limiting scan law, up to finite Monte Carlo error" else
          "marginal randomization under row exchangeability"
      )
    } else {
      combined <- .cp_cauchy_combine(p.ridge, cauchy_weights)
      combined$calibration.exact <- FALSE
      combined$exact.scope <- paste(
        "Analytic Cauchy transformation of dependent ridge p-values;",
        "the primary exact joint-limit quantile requires cross-ridge r_E",
        "and is not claimed here."
      )
      final.statistic <- combined$statistic
      p.value <- combined$p.value
    }

    ridge.rows <- lapply(seq_along(ridge), function(j) {
      selected <- observed.summary$selected[[j]]
      data.frame(
        ridge = ridge[[j]],
        statistic = observed.summary$statistic[[j]],
        p.value = p.ridge[[j]],
        critical.value = critical[[j]],
        reject = if (is.na(p.ridge[[j]])) NA else p.ridge[[j]] <= alpha,
        candidate = selected,
        q1 = candidates$q1[[selected]],
        q2 = candidates$q2[[selected]],
        q3 = candidates$q3[[selected]],
        t1 = candidates$t1[[selected]],
        t2 = candidates$t2[[selected]],
        t3 = candidates$t3[[selected]]
      )
    })
    ridge.rows <- do.call(rbind, ridge.rows)
    best.ridge <- which.max(observed.summary$statistic)

    structure(
      list(
        statistic = stats::setNames(
          final.statistic,
          if (length(ridge) > 1L && calibration != "none")
            "analytic.Cauchy" else "ERHT.scan"
        ),
        parameter = c(
          ridge.count = length(ridge),
          epsilon = epsilon,
          alpha = alpha
        ),
        p.value = p.value,
        estimate = c(
          change.point = ridge.rows$q2[[best.ridge]],
          change.fraction = ridge.rows$t2[[best.ridge]],
          ridge = ridge[[best.ridge]]
        ),
        null.value = c(number_of_changes = 0),
        alternative = paste("at least one", scan, "elliptical location change"),
        method = paste("Song--Wen--Feng ERHT", scan, "change-point scan"),
        data.name = data.name,
        valid = TRUE,
        components = list(
          ridge.results = ridge.rows,
          combined = combined,
          local = list(
            candidates = candidates,
            Z = observed$z,
            raw = observed$raw,
            center = observed$center,
            variance = observed$variance,
            kappa = observed$kappa,
            sigma2 = observed$sigma2,
            delta = observed$delta,
            beta = observed$beta,
            N.effective = observed$N.effective
          ),
          pool.location = observed$pool.location,
          pool.signs = observed$pool.signs,
          spatial.sign.covariance = observed$R.hat,
          ridge.inverse = observed$Q,
          reject = if (is.na(p.value)) NA else p.value <= alpha
        ),
        diagnostics = list(
          iteration = observed$fit.summary,
          calibration = calibration.details,
          calibration.reference = if (keep_calibration && !is.null(reference))
            reference$statistic else NULL,
          calibration.exact = if (is.null(combined)) NA else
            combined$calibration.exact,
          ridge.parameterization = paste(
            "Actual rho values were used. No automatic multiplication by p/n",
            "or empirical rho/gamma convention was applied."
          ),
          resolved.ridge = ridge,
          companion.reciprocal.condition = observed$companion.rcond,
          jacobian = observed$jacobian.convention,
          book.errata = paste(
            "The analytic Cauchy transform is not the exact dependent-ridge",
            "joint-limit calibration; exactness requires a joint quantile with r_E."
          ),
          no.repair = paste(
            "Only the positive ridge explicitly requested by the method was",
            "used. No additional ridge, pseudoinverse, variance floor, jitter,",
            "or dropped invalid scan candidate was used."
          ),
          call = call
        )
      ),
      class = c("chapter4_change_point_test", "htest")
    )
  }, error = function(error) {
    .cp_invalid(
      conditionMessage(error), strict,
      paste("Song--Wen--Feng ERHT", scan, "change-point scan"),
      call, data.name,
      diagnostics = list(
        resolved.ridge = ridge,
        requested.candidates = nrow(candidates),
        calibration = calibration
      )
    )
  })
}


.cp_sequence_including_end <- function(from, to, by) {
  if (from > to) return(integer())
  value <- seq.int(from, to, by = by)
  unique(c(value, to))
}

.cp_wbs_candidate_triples <- function(length, epsilon, triple_step) {
  minimum <- ceiling(epsilon * length)
  rows <- list()
  id <- 1L
  first.boundary <- minimum
  last.boundary <- length - minimum
  if (first.boundary > last.boundary) {
    stop("A WBS interval has no admissible adjacent triple.", call. = FALSE)
  }
  for (q2 in seq.int(first.boundary, last.boundary)) {
    q1.max <- q2 - minimum
    q3.min <- q2 + minimum
    q1.grid <- .cp_sequence_including_end(0L, q1.max, triple_step)
    q3.grid <- .cp_sequence_including_end(q3.min, length, triple_step)
    for (q1 in q1.grid) for (q3 in q3.grid) {
      rows[[id]] <- data.frame(
        q1 = q1, q2 = q2, q3 = q3,
        t1 = q1 / length, t2 = q2 / length, t3 = q3 / length
      )
      id <- id + 1L
    }
  }
  do.call(rbind, rows)
}

.cp_wbs_random_intervals <- function(n, M) {
  answer <- matrix(NA_integer_, M, 2L)
  for (i in seq_len(M)) answer[i, ] <- sort(sample.int(n, 2L))
  unique(answer)
}

.cp_wbs_interval_score <- function(x, left, right, ridge, epsilon,
                                   triple_step, controls) {
  pool <- x[seq.int(left, right), , drop = FALSE]
  candidates <- .cp_wbs_candidate_triples(
    nrow(pool), epsilon, triple_step
  )
  core <- .cp_erht_core(pool, ridge, candidates, controls)
  selected <- which.max(core$z)
  index <- arrayInd(selected, dim(core$z))
  candidate <- index[[1L]]
  ridge.index <- index[[2L]]
  list(
    score = core$z[candidate, ridge.index],
    k = left - 1L + candidates$q2[[candidate]],
    a = left + candidates$q1[[candidate]],
    c = left - 1L + candidates$q3[[candidate]],
    ridge = ridge[[ridge.index]],
    ridge.index = ridge.index,
    candidate = candidate,
    candidate.count = nrow(candidates),
    core = core
  )
}

.cp_wbs_invalid <- function(message, strict, call, data.name,
                            diagnostics = list()) {
  if (strict) stop(message, call. = FALSE)
  warning(message, call. = FALSE)
  structure(
    list(
      estimate = NULL,
      valid = FALSE,
      method = "Song--Wen--Feng WBS--ERHT segmentation",
      data.name = data.name,
      diagnostics = c(
        list(
          failure = message,
          no.repair = paste(
            "No additional ridge, variance floor, jitter, pseudoinverse,",
            "failed-candidate deletion, or fabricated threshold was used."
          ),
          call = call
        ),
        diagnostics
      )
    ),
    class = "chapter4_change_point_fit"
  )
}


#' Wild binary segmentation with the ERHT local score
#'
#' Implements the primary WBS--ERHT algorithm of Song, Wen, and Feng (2026).
#' On an interval `I`, every (or, when `triple_step > 1`, an explicitly
#' thinned) adjacent triple compares two segments of length at least
#' `ceiling(epsilon * length(I))`; the common scatter pool is the whole
#' interval.  Scores are maximized over the supplied actual ridge grid.
#'
#' The recursion uses the paper's narrowest-over-threshold rule: among sampled
#' intervals contained in the current segment, plus the current segment
#' itself, retain those with score strictly greater than `threshold`; choose
#' the narrowest, break a length tie by the larger score, refine within the
#' requested radius, delete the requested neighborhood, and recurse.
#'
#' The primary theory does not choose the ridge grid, WBS threshold, number of
#' intervals, minimum length, refinement radius, or deletion radius.  Hence
#' this API requires those quantities (or explicit intervals) rather than
#' inventing automatic tuning constants.  Random WBS intervals are generated
#' independently of the observations from uniformly sampled unordered endpoint
#' pairs.  Scanning/WBS is intrinsic to the estimator and is not replication of
#' a paper simulation study.
#'
#' `triple_step = 1` is the literal exhaustive adjacent-triple collection.
#' Larger values are an explicit computational approximation and are recorded.
#' Any invalid spatial median or local variance invalidates the fit; failed
#' triples are never silently removed.
#'
#' @param x Numeric `n` by `p` data matrix.
#' @param ridge Positive actual ERHT ridge value or grid.
#' @param threshold Required strict WBS score threshold.
#' @param min_interval Required minimum recursive/WBS interval length.
#' @param refinement_radius Required non-negative local half-window radius.
#' @param deletion_radius Required non-negative recursion deletion radius.
#' @param intervals Optional two-column integer matrix of WBS intervals.
#' @param M Required number of uniformly sampled intervals when `intervals` is
#'   `NULL`.
#' @param epsilon Primary minimum adjacent-segment fraction.
#' @param triple_step Candidate endpoint thinning step; one is exact.
#' @param max_changes Optional explicit software stopping cap; `Inf` follows
#'   the primary recursion until no interval crosses the threshold.
#' @param seed Optional random-interval seed; previous RNG state is restored.
#' @param tol,max_iter,zero_tol Spatial-median controls.
#' @param strict Failure contract.  With `FALSE`, a method failure returns
#'   `estimate = NULL`, never a partially certified segmentation.
#' @return A change-point fit with sorted estimates and selection history.
#' @references Song, Wen and Feng (2026), arXiv:2607.22162, Algorithm 1.
#' @export
#' @examples
#' x <- rbind(
#'   c(-2, 0), c(2, 0), c(0, -2), c(0, 2),
#'   c(-1, -1), c(1, 1), c(-1, 1), c(1, -1),
#'   c(-2, 1), c(2, -1), c(-1, 2), c(1, -2)
#' )
#' erht_wbs(
#'   x, ridge = 0.5, threshold = 10,
#'   min_interval = 8, refinement_radius = 3, deletion_radius = 1,
#'   intervals = matrix(c(1, 12), ncol = 2), epsilon = 0.25
#' )
erht_wbs <- function(
    x,
    ridge,
    threshold,
    min_interval,
    refinement_radius,
    deletion_radius,
    intervals = NULL,
    M = NULL,
    epsilon = 0.1,
    triple_step = 1L,
    max_changes = Inf,
    seed = NULL,
    tol = 1e-8,
    max_iter = 1000L,
    zero_tol = 0,
    strict = TRUE) {
  call <- match.call()
  data.name <- deparse(substitute(x))
  x <- .cp_data_matrix(x, minimum_rows = 6L, minimum_columns = 2L)
  n <- nrow(x)
  ridge <- as.numeric(ridge)
  if (!length(ridge) || anyNA(ridge) || any(!is.finite(ridge)) ||
      any(ridge <= 0)) {
    stop("`ridge` must contain positive finite actual rho values.",
         call. = FALSE)
  }
  ridge <- unique(ridge)
  threshold <- .cp_scalar(threshold, "threshold")
  min_interval <- .cp_integer(min_interval, "min_interval", minimum = 2L)
  refinement_radius <- .cp_integer(
    refinement_radius, "refinement_radius", minimum = 0L
  )
  deletion_radius <- .cp_integer(
    deletion_radius, "deletion_radius", minimum = 0L
  )
  epsilon <- .cp_scalar(
    epsilon, "epsilon", lower = 0, upper = 0.5,
    lower_open = TRUE, upper_open = TRUE
  )
  triple_step <- .cp_integer(triple_step, "triple_step")
  if (is.finite(max_changes)) {
    max_changes <- .cp_integer(max_changes, "max_changes")
  } else if (!identical(as.numeric(max_changes), Inf)) {
    stop("`max_changes` must be a positive integer or Inf.", call. = FALSE)
  }
  controls <- .cp_iteration_controls(tol, max_iter, zero_tol)
  strict <- .cp_flag(strict, "strict")
  if (min_interval > n) stop("`min_interval` exceeds the sample size.", call. = FALSE)
  if (2L * ceiling(epsilon * min_interval) > min_interval) {
    stop("`min_interval` and `epsilon` cannot form two admissible sides.",
         call. = FALSE)
  }

  if (is.null(intervals)) {
    if (is.null(M)) {
      stop("Supply either explicit `intervals` or a positive `M`.",
           call. = FALSE)
    }
    M <- .cp_integer(M, "M")
    intervals <- .cp_with_seed(seed, .cp_wbs_random_intervals(n, M))
    interval.source <- "uniform unordered endpoint pairs"
  } else {
    if (!is.matrix(intervals) || ncol(intervals) != 2L ||
        !is.numeric(intervals) || anyNA(intervals) ||
        any(!is.finite(intervals)) || any(intervals != floor(intervals))) {
      stop("`intervals` must be a finite two-column integer matrix.",
           call. = FALSE)
    }
    storage.mode(intervals) <- "integer"
    intervals <- t(apply(intervals, 1L, sort))
    if (any(intervals < 1L) || any(intervals > n) ||
        any(intervals[, 1L] >= intervals[, 2L])) {
      stop("Every WBS interval must satisfy 1 <= left < right <= n.",
           call. = FALSE)
    }
    intervals <- unique(intervals)
    interval.source <- "user supplied"
    M <- nrow(intervals)
  }

  tryCatch({
    score.cache <- new.env(parent = emptyenv())
    score.interval <- function(left, right) {
      key <- paste(left, right, sep = ":")
      cached <- get0(key, envir = score.cache, inherits = FALSE)
      if (!is.null(cached)) return(cached)
      value <- .cp_wbs_interval_score(
        x, left, right, ridge, epsilon, triple_step, controls
      )
      assign(key, value, envir = score.cache)
      value
    }

    selected <- integer()
    history <- list()
    recurse <- function(left, right) {
      if (right - left + 1L < min_interval ||
          length(selected) >= max_changes) return(invisible(NULL))
      contained <- intervals[
        intervals[, 1L] >= left & intervals[, 2L] <= right &
          intervals[, 2L] - intervals[, 1L] + 1L >= min_interval,
        , drop = FALSE
      ]
      candidate.intervals <- unique(rbind(contained, c(left, right)))
      scores <- vector("list", nrow(candidate.intervals))
      for (ii in seq_len(nrow(candidate.intervals))) {
        scores[[ii]] <- score.interval(
          candidate.intervals[ii, 1L], candidate.intervals[ii, 2L]
        )
      }
      values <- vapply(scores, function(value) value$score, numeric(1))
      crossing <- which(values > threshold)
      if (!length(crossing)) return(invisible(NULL))
      lengths <- candidate.intervals[crossing, 2L] -
        candidate.intervals[crossing, 1L] + 1L
      order.index <- order(
        lengths,
        -values[crossing],
        candidate.intervals[crossing, 1L],
        candidate.intervals[crossing, 2L]
      )
      chosen.row <- crossing[order.index[[1L]]]
      chosen.interval <- candidate.intervals[chosen.row, ]
      rough <- scores[[chosen.row]]
      local.left <- max(left, rough$k - refinement_radius)
      local.right <- min(right, rough$k + refinement_radius)
      if (local.right - local.left + 1L >= min_interval) {
        refined <- score.interval(local.left, local.right)
      } else {
        refined <- rough
        local.left <- chosen.interval[[1L]]
        local.right <- chosen.interval[[2L]]
      }
      k.hat <- refined$k
      if (k.hat < left || k.hat >= right) {
        stop("A refined WBS split lies outside its recursive segment.")
      }
      selected <<- c(selected, k.hat)
      history[[length(history) + 1L]] <<- data.frame(
        order = length(history) + 1L,
        recursive.left = left,
        recursive.right = right,
        interval.left = chosen.interval[[1L]],
        interval.right = chosen.interval[[2L]],
        rough.change = rough$k,
        refined.left = local.left,
        refined.right = local.right,
        change.point = k.hat,
        score = refined$score,
        ridge = refined$ridge,
        triple.a = refined$a,
        triple.c = refined$c,
        candidate.count = refined$candidate.count
      )
      recurse(left, k.hat - deletion_radius)
      recurse(k.hat + deletion_radius + 1L, right)
      invisible(NULL)
    }

    recurse(1L, n)
    estimate <- sort(unique(selected))
    history.frame <- if (length(history)) do.call(rbind, history) else
      data.frame()
    cached.keys <- ls(score.cache, all.names = TRUE)
    cached.fits <- lapply(cached.keys, function(key) get(key, score.cache))
    maximum.update <- if (length(cached.fits)) max(vapply(
      cached.fits,
      function(value) value$core$fit.summary$maximum.relative.update,
      numeric(1)
    )) else NA_real_
    maximum.score.residual <- if (length(cached.fits)) max(vapply(
      cached.fits,
      function(value) value$core$fit.summary$maximum.score.residual,
      numeric(1)
    )) else NA_real_

    structure(
      list(
        estimate = estimate,
        change.fraction = estimate / n,
        valid = TRUE,
        method = "Song--Wen--Feng WBS--ERHT segmentation",
        data.name = data.name,
        history = history.frame,
        intervals = intervals,
        tuning = list(
          ridge = ridge,
          threshold = threshold,
          epsilon = epsilon,
          min.interval = min_interval,
          refinement.radius = refinement_radius,
          deletion.radius = deletion_radius,
          triple.step = triple_step,
          exhaustive.triples = identical(triple_step, 1L),
          max.changes = max_changes
        ),
        diagnostics = list(
          interval.source = interval.source,
          interval.count = nrow(intervals),
          scored.interval.count = length(cached.keys),
          threshold.rule = "strict score > threshold",
          selection.rule = paste(
            "narrowest over threshold; ties by larger score and then",
            "deterministic endpoint order"
          ),
          maximum.relative.update = maximum.update,
          maximum.score.residual = maximum.score.residual,
          jacobian = paste(
            "Primary e=mean(sqrt(p)/distance), coincident weight zero;",
            "no optional (p-1)/p factor."
          ),
          threshold.source = "user supplied; primary theory has no automatic finite-sample threshold",
          no.repair = paste(
            "Only the supplied ERHT ridge was used. No invalid candidate was",
            "dropped, and no variance floor, jitter, or pseudoinverse was used."
          ),
          call = call
        )
      ),
      class = "chapter4_change_point_fit"
    )
  }, error = function(error) {
    .cp_wbs_invalid(
      conditionMessage(error), strict, call, data.name,
      diagnostics = list(
        resolved.ridge = ridge,
        threshold = threshold,
        interval.source = interval.source,
        interval.count = nrow(intervals)
      )
    )
  })
}


#' @export
print.chapter4_change_point_fit <- function(x, ...) {
  cat("\n", x$method, "\n\n", sep = "")
  cat("data:  ", x$data.name, "\n", sep = "")
  if (!isTRUE(x$valid)) {
    cat("invalid fit: ", x$diagnostics$failure, "\n", sep = "")
    return(invisible(x))
  }
  if (length(x$estimate)) {
    cat("estimated changes: ", paste(x$estimate, collapse = ", "), "\n",
        sep = "")
  } else {
    cat("estimated changes: none\n")
  }
  invisible(x)
}
