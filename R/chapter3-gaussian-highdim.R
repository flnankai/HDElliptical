.ch3g_beta <- function(beta, residual) {
  if (is.null(beta)) {
    second <- mean(residual^2)
    if (!is.finite(second) || second <= 0) {
      stop("The empirical variance needed to estimate `beta` is non-positive.",
           call. = FALSE)
    }
    value <- mean(residual^4) / second^2 - 3
    source <- "scale-standardized empirical fourth moment"
  } else {
    if (!is.numeric(beta) || length(beta) != 1L || !is.finite(beta)) {
      stop("`beta` must be one finite number or NULL.", call. = FALSE)
    }
    value <- as.numeric(beta)
    source <- "supplied"
  }
  list(value = value, source = source)
}


.ch3g_wang_yao_inputs <- function(x, center, beta, require_below_one) {
  x <- .as_data_matrix(x, "x", min_rows = 3L)
  work <- .ch3g_centered_scaled(x, center)
  n <- nrow(x)
  p <- ncol(x)
  m <- work$effective.df
  if (require_below_one && p >= m) {
    stop("The corrected likelihood-ratio test requires p / effective.df < 1.",
         call. = FALSE)
  }
  beta.fit <- .ch3g_beta(beta, work$residual)
  S <- crossprod(work$residual) / m
  list(
    x = x, residual = work$residual, S = S, n = n, p = p, m = m,
    center = center, scale = work$scale, beta = beta.fit$value,
    beta.source = beta.fit$source
  )
}


#' Wang--Yao corrected likelihood-ratio test of sphericity
#'
#' For effective covariance degrees of freedom \eqn{m} and \eqn{y=p/m<1},
#' this function computes
#' \deqn{\mathcal L=-\log|S|+p\log\{\operatorname{tr}(S)/p\}}
#' and calibrates
#' \deqn{\mathcal L+(p-m)\log(1-p/m)-p}
#' by a normal distribution with mean
#' \eqn{-\log(1-y)/2+\beta y/2} and variance
#' \eqn{-2\log(1-y)-2y} for real data. `center = FALSE` is the known-zero-mean
#' formula in the book. `center = TRUE` uses Wang and Yao's unknown-mean
#' extension, replacing the spectral centering ratio by \eqn{p/(n-1)}.
#'
#' `beta = 0` gives the real Gaussian correction. A supplied value is the
#' standardized fourth cumulant \eqn{E Z^4-3}; `beta = NULL` uses the
#' scale-standardized empirical fourth moment under the spherical null.
#'
#' @param x Numeric matrix with observations in rows.
#' @param center Whether the mean is estimated and removed.
#' @param beta Finite fourth cumulant, or `NULL` for the empirical plug-in.
#'
#' @return An `hd_covariance_test` object with all centering terms retained.
#' @references Wang, Q. and Yao, J. (2013). *Electronic Journal of
#'   Statistics*, 7, 2164--2192. \doi{10.1214/13-EJS842}.
#' @examples
#' set.seed(40)
#' x <- matrix(rnorm(72), nrow = 24, ncol = 3)
#' wang_yao_corrected_lrt(x, center = TRUE)
#' @export
wang_yao_corrected_lrt <- function(x, center = FALSE, beta = 0) {
  call <- match.call()
  data.name <- deparse1(substitute(x))
  center <- .ch3g_validate_center(center)
  input <- .ch3g_wang_yao_inputs(x, center, beta, TRUE)
  values <- eigen(input$S, symmetric = TRUE, only.values = TRUE)$values
  if (any(!is.finite(values)) || any(values <= 0)) {
    stop("The sample covariance determinant is not strictly positive.",
         call. = FALSE)
  }
  log.ratio <- -sum(log(values)) + input$p * log(mean(values))
  y <- input$p / input$m
  corrected <- log.ratio + (input$p - input$m) * log1p(-y) - input$p
  null.mean <- -0.5 * log1p(-y) + 0.5 * input$beta * y
  null.variance <- -2 * log1p(-y) - 2 * y
  if (!is.finite(null.variance) || null.variance <= 0) {
    stop("The corrected LRT null variance is not positive.", call. = FALSE)
  }
  z <- (corrected - null.mean) / sqrt(null.variance)
  .ch3g_new_test(
    statistic = c(Z = z),
    p.value = stats::pnorm(z, lower.tail = FALSE),
    method = "Wang-Yao corrected high-dimensional sphericity LRT",
    data.name = data.name,
    raw.statistic = c(L.cal = corrected, L = log.ratio),
    null.distribution = list(
      family = "normal", parameters = c(mean = 0, sd = 1), exact = FALSE,
      tail = "upper",
      assumptions = paste(
        "real iid standardized entries with finite fourth moment;",
        "p / effective.df converges in (0, 1)"
      )
    ),
    components = list(
      L = log.ratio,
      corrected.L = corrected,
      aspect.ratio = y,
      beta = input$beta,
      null.mean = null.mean,
      null.variance = null.variance
    ),
    diagnostics = list(
      center.estimated = center,
      sample.size = input$n,
      effective.covariance.df = input$m,
      beta.source = input$beta.source,
      internal.data.scale = input$scale,
      attribution = paste(
        "sphericity formula is Wang-Yao (2013); Bai et al. (2009)",
        "develop the preceding given-covariance RMT LRT"
      ),
      regularization = "none"
    ),
    n = input$n, p = input$p, call = call
  )
}


#' Wang--Yao corrected John's test of sphericity
#'
#' For \eqn{U=p\operatorname{tr}(S^2)/\operatorname{tr}^2(S)-1}, the
#' known-zero-mean result is
#' \deqn{nU-p\ \Rightarrow\ N(1+\beta,4).}
#' With an estimated mean, the primary extension instead centers
#' \eqn{nU} by \eqn{np/(n-1)}. The test remains defined for \eqn{p\ge n}.
#'
#' @inheritParams wang_yao_corrected_lrt
#' @return An `hd_covariance_test` object.
#' @references
#' Wang, Q. and Yao, J. (2013). *Electronic Journal of Statistics*,
#' 7, 2164--2192. \doi{10.1214/13-EJS842}.
#' @examples
#' set.seed(41)
#' x <- matrix(rnorm(72), nrow = 24, ncol = 3)
#' wang_yao_corrected_john_test(x, center = TRUE)
#' @export
wang_yao_corrected_john_test <- function(x, center = FALSE, beta = 0) {
  call <- match.call()
  data.name <- deparse1(substitute(x))
  center <- .ch3g_validate_center(center)
  input <- .ch3g_wang_yao_inputs(x, center, beta, FALSE)
  trace.S <- sum(diag(input$S))
  if (!is.finite(trace.S) || trace.S <= 0) {
    stop("The sample covariance trace is non-positive.", call. = FALSE)
  }
  U <- input$p * sum(input$S * input$S) / trace.S^2 - 1
  spectral.center <- if (center) input$n * input$p / input$m else input$p
  centered <- input$n * U - spectral.center
  null.mean <- 1 + input$beta
  z <- (centered - null.mean) / 2
  .ch3g_new_test(
    statistic = c(Z = z),
    p.value = stats::pnorm(z, lower.tail = FALSE),
    method = "Wang-Yao corrected high-dimensional John's sphericity test",
    data.name = data.name,
    raw.statistic = c(nU.centered = centered, U = U),
    null.distribution = list(
      family = "normal", parameters = c(mean = 0, sd = 1), exact = FALSE,
      tail = "upper",
      assumptions = paste(
        "real iid standardized entries with finite fourth moment;",
        "p / effective.df has a finite positive limit"
      )
    ),
    components = list(
      U = U,
      nU = input$n * U,
      spectral.center = spectral.center,
      centered.nU = centered,
      beta = input$beta,
      null.mean = null.mean,
      null.variance = 4
    ),
    diagnostics = list(
      center.estimated = center,
      sample.size = input$n,
      effective.covariance.df = input$m,
      beta.source = input$beta.source,
      internal.data.scale = input$scale,
      regularization = "none"
    ),
    n = input$n, p = input$p, call = call
  )
}


#' Chen--Zhang--Zhong high-dimensional covariance test
#'
#' Constructs the location-invariant unbiased estimators
#' \deqn{T_{1n}=Y_{1n}-Y_{3n},\qquad
#' T_{2n}=Y_{2n}-2Y_{4n}+Y_{5n}}
#' of \eqn{\operatorname{tr}(\Sigma)} and
#' \eqn{\operatorname{tr}(\Sigma^2)}. For sphericity it uses
#' \eqn{U_n=pT_{2n}/T_{1n}^2-1}; for identity it uses
#' \eqn{V_n=T_{2n}/p-2T_{1n}/p+1}. In either case the primary null
#' calibration is \eqn{nU_n/2} or \eqn{nV_n/2}, respectively, against an
#' upper standard-normal tail.
#'
#' @param x Numeric matrix with observations in rows; at least four rows.
#' @param null Either `"sphericity"` or `"identity"`.
#'
#' @return An `hd_covariance_test` object retaining all five U-statistic
#'   building blocks.
#' @references Chen, S. X., Zhang, L.-X. and Zhong, P.-S. (2010).
#'   *Journal of the American Statistical Association*, 105, 810--819.
#'   \doi{10.1198/jasa.2010.tm09560}.
#' @examples
#' set.seed(42)
#' x <- matrix(rnorm(48), nrow = 12, ncol = 4)
#' chen_zhang_zhong_covariance_test(x, null = "sphericity")
#' @export
chen_zhang_zhong_covariance_test <- function(
    x, null = c("sphericity", "identity")) {
  call <- match.call()
  data.name <- deparse1(substitute(x))
  x <- .as_data_matrix(x, "x", min_rows = 4L)
  null <- match.arg(null)
  n <- nrow(x)
  p <- ncol(x)
  if (null == "sphericity" && p < 2L) {
    stop("Sphericity requires at least two variables.", call. = FALSE)
  }
  scale <- max(abs(x))
  if (!is.finite(scale) || scale <= 0) {
    stop("`x` must contain nonzero finite data.", call. = FALSE)
  }
  residual <- sweep(x / scale, 2L, colMeans(x / scale), "-")
  out <- cpp_ch3_gaussian_trace_u_statistics(residual)
  T1.scaled <- as.numeric(out$T1)
  T2.scaled <- as.numeric(out$T2)
  if (!is.finite(T1.scaled) || T1.scaled <= 0 || !is.finite(T2.scaled)) {
    stop("The unbiased trace estimates are degenerate or non-finite.",
         call. = FALSE)
  }
  if (null == "sphericity") {
    raw <- p * T2.scaled / T1.scaled^2 - 1
    z <- n * raw / 2
    method <- "Chen-Zhang-Zhong high-dimensional sphericity test"
    raw.name <- "U.n"
  } else {
    T1 <- T1.scaled * scale^2
    T2 <- T2.scaled * scale^4
    if (!is.finite(T1) || !is.finite(T2)) {
      stop("The identity-test trace estimates overflow in original units.",
           call. = FALSE)
    }
    raw <- T2 / p - 2 * T1 / p + 1
    z <- n * raw / 2
    method <- "Chen-Zhang-Zhong high-dimensional identity test"
    raw.name <- "V.n"
  }
  raw.statistic <- raw
  names(raw.statistic) <- raw.name
  original <- c(
    T1 = T1.scaled * scale^2,
    T2 = T2.scaled * scale^4
  )
  .ch3g_new_test(
    statistic = c(Z = z),
    p.value = stats::pnorm(z, lower.tail = FALSE),
    method = method,
    data.name = data.name,
    raw.statistic = raw.statistic,
    null.value = if (null == "identity") diag(p) else "positive multiple of I",
    null.distribution = list(
      family = "normal", parameters = c(mean = 0, sd = 1), exact = FALSE,
      tail = "upper",
      assumptions = paste(
        "Chen-Zhang-Zhong factor model with bounded eighth moments;",
        "trace(Sigma^4) / trace(Sigma^2)^2 tends to zero"
      )
    ),
    components = list(
      null = null,
      Y1.scaled = as.numeric(out$Y1),
      Y2.scaled = as.numeric(out$Y2),
      Y3.scaled = as.numeric(out$Y3),
      Y4.scaled = as.numeric(out$Y4),
      Y5.scaled = as.numeric(out$Y5),
      T1.scaled = T1.scaled,
      T2.scaled = T2.scaled,
      T1 = original[["T1"]],
      T2 = original[["T2"]]
    ),
    diagnostics = list(
      internal.data.scale = scale,
      centering = "algebraically neutral common translation removal",
      regularization = "none",
      trace.repair = "none"
    ),
    n = n, p = p, call = call
  )
}


#' Fisher--Sun--Gallagher fourth-to-second trace sphericity test
#'
#' For \eqn{N=n+1} observations, this implements the primary paper's unbiased
#' estimators \eqn{\hat a_2} and \eqn{\hat a_4} of
#' \eqn{p^{-1}\operatorname{tr}(\Sigma^2)} and
#' \eqn{p^{-1}\operatorname{tr}(\Sigma^4)}. With \eqn{c=p/n}, the statistic
#' \deqn{\sqrt{\frac{np}{8(8+12c+c^2)}}
#' \left(\frac{\hat a_4}{\hat a_2^2}-1\right)}
#' has an upper standard-normal null calibration. The calculation is performed
#' after a common scale normalization, which is algebraically neutral.
#'
#' @param x Numeric matrix with observations in rows; at least five rows.
#' @return An `hd_covariance_test` object retaining every trace and coefficient.
#' @references Fisher, T. J., Sun, X. and Gallagher, C. M. (2010). *Journal
#'   of Multivariate Analysis*, 101, 2554--2570.
#'   \doi{10.1016/j.jmva.2010.07.004}.
#' @examples
#' set.seed(43)
#' x <- matrix(rnorm(60), nrow = 15, ncol = 4)
#' fisher_sun_gallagher_sphericity_test(x)
#' @export
fisher_sun_gallagher_sphericity_test <- function(x) {
  call <- match.call()
  data.name <- deparse1(substitute(x))
  x <- .as_data_matrix(x, "x", min_rows = 5L)
  N <- nrow(x)
  n <- N - 1
  p <- ncol(x)
  if (p < 2L) stop("Sphericity requires at least two variables.", call. = FALSE)
  scale <- max(abs(x))
  if (!is.finite(scale) || scale <= 0) {
    stop("`x` must contain nonzero finite data.", call. = FALSE)
  }
  residual <- sweep(x / scale, 2L, colMeans(x / scale), "-")
  singular <- svd(residual, nu = 0L, nv = 0L)$d
  values <- singular^2 / n
  traces <- c(
    tr1 = sum(values),
    tr2 = sum(values^2),
    tr3 = sum(values^3),
    tr4 = sum(values^4)
  )
  b <- -4 / n
  c.star <- -(2 * n^2 + 3 * n - 6) / (n * (n^2 + n + 2))
  d <- 2 * (5 * n + 6) / (n * (n^2 + n + 2))
  e <- -(5 * n + 6) / (n^2 * (n^2 + n + 2))
  tau <- n^5 * (n^2 + n + 2) /
    ((n + 1) * (n + 2) * (n + 4) * (n + 6) *
       (n - 1) * (n - 2) * (n - 3))
  a2 <- n^2 / ((n - 1) * (n + 2) * p) *
    (traces[["tr2"]] - traces[["tr1"]]^2 / n)
  a4 <- tau / p * (
    traces[["tr4"]] + b * traces[["tr3"]] * traces[["tr1"]] +
      c.star * traces[["tr2"]]^2 +
      d * traces[["tr2"]] * traces[["tr1"]]^2 +
      e * traces[["tr1"]]^4
  )
  if (!is.finite(a2) || a2 <= 0 || !is.finite(a4)) {
    stop("The Fisher--Sun--Gallagher unbiased trace estimates are degenerate.",
         call. = FALSE)
  }
  psi <- a4 / a2^2
  concentration <- p / n
  multiplier <- sqrt(n * p / (8 * (8 + 12 * concentration + concentration^2)))
  z <- multiplier * (psi - 1)
  .ch3g_new_test(
    statistic = c(Z = z),
    p.value = stats::pnorm(z, lower.tail = FALSE),
    method = "Fisher-Sun-Gallagher high-dimensional sphericity test",
    data.name = data.name,
    raw.statistic = c(psi.2 = psi),
    null.distribution = list(
      family = "normal", parameters = c(mean = 0, sd = 1), exact = FALSE,
      tail = "upper",
      assumptions = paste(
        "Gaussian sampling; finite limiting spectral moments through order 16;",
        "p / (N - 1) has a finite positive limit"
      )
    ),
    components = list(
      sample.trace.powers.scaled = traces,
      a2.scaled = a2,
      a4.scaled = a4,
      psi.2 = psi,
      concentration = concentration,
      multiplier = multiplier,
      coefficients = c(b = b, c.star = c.star, d = d, e = e, tau = tau)
    ),
    diagnostics = list(
      observations = N,
      Wishart.df = n,
      covariance.divisor = n,
      internal.data.scale = scale,
      trace.repair = "none"
    ),
    n = N, p = p, call = call
  )
}


#' Li--Chen two-sample high-dimensional covariance equality test
#'
#' Computes the exact location-invariant leave-four-out statistics
#' \eqn{A_{n_1}}, \eqn{A_{n_2}}, and the two-by-two leave-out cross statistic
#' \eqn{C_{n_1n_2}} from Li and Chen (2012). Their combination
#' \deqn{T_{LC}=A_{n_1}+A_{n_2}-2C_{n_1n_2}}
#' is unbiased for \eqn{\operatorname{tr}(\Sigma_1-\Sigma_2)^2}. Under the
#' equality null, the primary feasible standard-error estimator is
#' \deqn{\hat\sigma_0=2A_{n_1}/n_2+2A_{n_2}/n_1.}
#' The published article labels this quantity once as \eqn{\hat\sigma_0^2};
#' its units, subsequent ratio-consistency theorem, and rejection rule show
#' that it is the standard error. A non-positive estimate is rejected without
#' an absolute value or floor.
#'
#' @param x,y Numeric matrices with observations in rows and matching columns;
#'   each group needs at least four observations.
#'
#' @return An `hd_covariance_test` object retaining original-unit and
#'   internally scaled trace estimates.
#' @references Li, J. and Chen, S. X. (2012). *Annals of Statistics*, 40,
#'   908--940. \doi{10.1214/12-AOS993}.
#' @examples
#' set.seed(44)
#' x <- matrix(rnorm(24), nrow = 8, ncol = 3)
#' y <- matrix(rnorm(27, mean = 0.1), nrow = 9, ncol = 3)
#' li_chen_covariance_test(x, y)
#' @export
li_chen_covariance_test <- function(x, y) {
  call <- match.call()
  x.name <- deparse1(substitute(x))
  y.name <- deparse1(substitute(y))
  x <- .as_data_matrix(x, "x", min_rows = 4L)
  y <- .as_data_matrix(y, "y", min_rows = 4L)
  .check_two_sample_variables(x, y)
  n1 <- nrow(x)
  n2 <- nrow(y)
  p <- ncol(x)
  scale <- max(abs(x), abs(y))
  if (!is.finite(scale) || scale <= 0) {
    stop("The two samples must contain nonzero finite data.", call. = FALSE)
  }
  xs <- sweep(x / scale, 2L, colMeans(x / scale), "-")
  ys <- sweep(y / scale, 2L, colMeans(y / scale), "-")
  out <- cpp_ch3_gaussian_li_chen(xs, ys)
  A1 <- as.numeric(out$A1)
  A2 <- as.numeric(out$A2)
  C <- as.numeric(out$C)
  T.LC <- as.numeric(out$T)
  standard.error <- 2 * A1 / n2 + 2 * A2 / n1
  if (!is.finite(standard.error) || standard.error <= 0) {
    stop("The Li--Chen estimated null standard error is not positive.",
         call. = FALSE)
  }
  z <- T.LC / standard.error
  unit.factor <- scale^4
  original <- c(
    A1 = A1 * unit.factor,
    A2 = A2 * unit.factor,
    C = C * unit.factor,
    T.LC = T.LC * unit.factor,
    standard.error = standard.error * unit.factor
  )
  .ch3g_new_test(
    statistic = c(Z = z),
    p.value = stats::pnorm(z, lower.tail = FALSE),
    method = "Li-Chen two-sample high-dimensional covariance equality test",
    data.name = paste(x.name, "and", y.name),
    raw.statistic = c(T.LC = original[["T.LC"]]),
    null.distribution = list(
      family = "normal", parameters = c(mean = 0, sd = 1), exact = FALSE,
      tail = "upper",
      assumptions = paste(
        "independent samples; Li-Chen factor-moment and trace regularity;",
        "positive feasible null standard error"
      )
    ),
    components = list(
      A1 = original[["A1"]],
      A2 = original[["A2"]],
      C = original[["C"]],
      T.LC = original[["T.LC"]],
      estimated.null.standard.error = original[["standard.error"]],
      A1.scaled = A1,
      A2.scaled = A2,
      C.scaled = C,
      T.LC.scaled = T.LC,
      estimated.null.standard.error.scaled = standard.error
    ),
    diagnostics = list(
      n1 = n1, n2 = n2,
      internal.common.scale = scale,
      centering = "algebraically neutral within-group centering",
      variance.label.erratum = paste(
        "primary display labels the feasible standard error as sigma-hat squared;",
        "dimensional analysis and the following theorem/rejection rule use sigma-hat"
      ),
      regularization = "none",
      standard.error.repair = "none"
    ),
    n = c(x = n1, y = n2), p = p, call = call
  )
}
