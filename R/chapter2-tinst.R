#' Huang--Liu--Zhou--Feng two-sample inverse norm sign test
#'
#' Tests equality of two multivariate location vectors with the two-sample
#' inverse norm sign test (tINST) of Huang, Liu, Zhou, and Feng (2023).  The
#' implementation uses the paper's feasible cross statistic, observation-wise
#' leave-one-out diagonal estimates, and the inverse norm weight
#' \eqn{\omega(r)=r^{-1}}.  The manuscript's three group-index typographical
#' errors are corrected: every leave-one-out sum uses the corresponding group
#' size \eqn{n_k}, every diagonal update is an update of \eqn{D_k}, and both
#' indices of the second within-group trace estimator run over group 2.
#'
#' Every full-sample and leave-one-out diagonal fit starts at the sample mean
#' and marginal sample variances.  The inverse-norm weighted location fit starts
#' at its leave-one-out sample mean.  The source article says to iterate until
#' convergence but does not prescribe a tolerance.  Here `iteration.stable`
#' means that the relative location update (and, for a joint fit, the
#' log-diagonal update) is at most `tol`.  The separately returned
#' `score.residual` is diagnostic only and is not presented as a certificate
#' that the estimating equation has been solved.  This distinction matters for
#' inverse-norm iteration because it can approach an included observation.
#' The diagonal fixed-point equations identify only relative coordinate scales;
#' this implementation retains the common scale carried by the sample-variance
#' initialization and applies no determinant or trace normalization.
#'
#' No zero marginal variance is replaced, no ridge is added, and a non-positive
#' or non-finite published variance estimate is an error.  With `strict = TRUE`
#' (the default), failure of any required iteration to become stable is also an
#' error.  Setting `strict = FALSE` returns the last iterates with a warning and
#' complete per-fit stability and score diagnostics.
#'
#' A common translation and common coordinatewise scaling are applied
#' internally before iteration.  This is algebraically neutral because tINST
#' is invariant to common shifts and nonsingular diagonal changes of units, and
#' it protects residual calculations from avoidable cancellation and overflow.
#'
#' @param x,y Numeric matrices or data frames with observations in rows and the
#'   same variables in columns.  Each group must have at least three rows.
#' @param alpha Significance level used for the returned upper-tail rejection
#'   decision.
#' @param tol Positive convergence tolerance for every joint and weighted
#'   location iteration.
#' @param max_iter Positive maximum number of updates for every individual fit.
#' @param zero_tol Non-negative tolerance for detecting a zero diagonally
#'   standardized residual.  The default detects exact zeros only.
#' @param strict If `TRUE`, error when any required iterative fit does not
#'   become stable.  If `FALSE`, warn and return the last iterates together with
#'   diagnostics.  A large score residual alone does not trigger this contract.
#'
#' @return An object of class `c("hd_location_test", "htest")`.  Its
#'   `components` field contains the feasible statistic, nuisance estimates,
#'   ordered-pair trace estimates, variance terms, leave-one-out radii and
#'   directions, and internally scaled fit estimates.  Its `diagnostics` field
#'   contains all iteration counts, update sizes, estimating-equation residuals,
#'   the rejection rule, and the explicit no-repair contract.
#'
#' @references
#' Huang, X., Liu, B., Zhou, Q., and Feng, L. (2023). A high-dimensional
#' inverse norm sign test for two-sample location problems. *Canadian Journal
#' of Statistics*, 51, 1004--1033. \doi{10.1002/cjs.11731}.
#'
#' @examples
#' set.seed(2023)
#' x <- matrix(stats::rt(24, df = 5), 8, 3)
#' y <- matrix(stats::rt(30, df = 5), 10, 3)
#' tinst_two_sample_test(x, y)
#'
#' @export
tinst_two_sample_test <- function(x, y, alpha = 0.05,
                                  tol = 1e-8, max_iter = 500L,
                                  zero_tol = 0, strict = TRUE) {
  call <- match.call()
  x.name <- deparse1(substitute(x))
  y.name <- deparse1(substitute(y))
  x <- .as_data_matrix(x, "x", min_rows = 3L)
  y <- .as_data_matrix(y, "y", min_rows = 3L)
  .check_two_sample_variables(x, y)

  alpha <- as.numeric(alpha)
  if (length(alpha) != 1L || is.na(alpha) || !is.finite(alpha) ||
      alpha <= 0 || alpha >= 1) {
    stop("`alpha` must be one finite number strictly between zero and one.",
         call. = FALSE)
  }
  controls <- .validate_iteration_controls(tol, max_iter, zero_tol)
  if (!is.logical(strict) || length(strict) != 1L || is.na(strict)) {
    stop("`strict` must be `TRUE` or `FALSE`.", call. = FALSE)
  }

  out <- cpp_tinst_two_sample(
    x, y, controls$tol, controls$max_iter, controls$zero_tol
  )
  if (!isTRUE(out$all_iteration_stable)) {
    message <- sprintf(
      paste0(
        "tINST did not become iteration-stable for %d required fit(s) in ",
        "`max_iter = %d` ",
        "updates. No ridge or iterate substitution was applied."
      ),
      as.integer(out$stability_failures), controls$max_iter
    )
    if (isTRUE(strict)) {
      stop(message, call. = FALSE)
    }
    warning(message, call. = FALSE)
  }

  n1 <- nrow(x)
  n2 <- nrow(y)
  p <- ncol(x)
  variable.names <- .hotelling_variable_names(x)
  z <- as.numeric(out$z)
  critical.value <- stats::qnorm(1 - alpha)
  reject <- isTRUE(z > critical.value)
  mean.x <- stats::setNames(as.numeric(out$mean_x), variable.names)
  mean.y <- stats::setNames(as.numeric(out$mean_y), variable.names)
  difference <- stats::setNames(as.numeric(out$difference), variable.names)
  full.diagonal1 <- stats::setNames(
    as.numeric(out$full_diagonal1_scaled), variable.names
  )
  full.diagonal2 <- stats::setNames(
    as.numeric(out$full_diagonal2_scaled), variable.names
  )
  worst.leaveout <- list(
    relative.update = max(
      out$joint_fits1$worst.relative.update,
      out$joint_fits2$worst.relative.update,
      out$weighted_fits1$worst.relative.update,
      out$weighted_fits2$worst.relative.update
    ),
    score.residual = max(
      out$joint_fits1$worst.score.residual,
      out$joint_fits2$worst.score.residual,
      out$weighted_fits1$worst.score.residual,
      out$weighted_fits2$worst.score.residual
    ),
    minimum.residual.distance = min(
      out$joint_fits1$smallest.residual.distance,
      out$joint_fits2$smallest.residual.distance,
      out$weighted_fits1$smallest.residual.distance,
      out$weighted_fits2$smallest.residual.distance
    ),
    maximum.iterations = max(
      out$joint_fits1$maximum.iterations,
      out$joint_fits2$maximum.iterations,
      out$weighted_fits1$maximum.iterations,
      out$weighted_fits2$maximum.iterations
    )
  )
  worst.overall.relative.update <- max(
    out$full_fit1$relative.update,
    out$full_fit2$relative.update,
    worst.leaveout$relative.update
  )
  worst.overall.score.residual <- max(
    out$full_fit1$score.residual,
    out$full_fit2$score.residual,
    worst.leaveout$score.residual
  )
  minimum.overall.residual.distance <- min(
    out$full_fit1$minimum.residual.distance,
    out$full_fit2$minimum.residual.distance,
    worst.leaveout$minimum.residual.distance
  )

  .new_hd_location_test(
    statistic = c(Z = z),
    p.value = stats::pnorm(z, lower.tail = FALSE),
    method = paste(
      "Huang-Liu-Zhou-Feng two-sample inverse norm sign test",
      "(tINST; asymptotic normal calibration)"
    ),
    data.name = paste(x.name, "and", y.name),
    alternative = "two.sided",
    raw.statistic = c(T.tINST = as.numeric(out$T)),
    estimate = difference,
    null.value = stats::setNames(numeric(p), variable.names),
    null.distribution = list(
      family = "normal",
      parameters = c(mean = 0, sd = 1),
      exact = FALSE,
      tail = "upper",
      assumptions = paste(
        "Huang-Liu-Zhou-Feng high-dimensional elliptical asymptotics;",
        "independent samples; finite inverse radial moments"
      )
    ),
    variance = c(denominator = as.numeric(out$variance)),
    components = list(
      mean.x = mean.x,
      mean.y = mean.y,
      difference = difference,
      n1 = as.numeric(out$n1),
      n2 = as.numeric(out$n2),
      p = as.numeric(out$p),
      T.tINST = as.numeric(out$T),
      nu12 = as.numeric(out$nu12),
      nu22 = as.numeric(out$nu22),
      c1 = as.numeric(out$c1),
      c2 = as.numeric(out$c2),
      A1 = as.numeric(out$A1),
      A2 = as.numeric(out$A2),
      A12 = as.numeric(out$A12),
      variance.term1 = as.numeric(out$variance_term1),
      variance.term2 = as.numeric(out$variance_term2),
      variance.term12 = as.numeric(out$variance_term12),
      denominator.variance = as.numeric(out$variance),
      full.diagonal1.scaled = full.diagonal1,
      full.diagonal2.scaled = full.diagonal2,
      leaveout.diagonal1.scaled = out$leaveout_diagonal1_scaled,
      leaveout.diagonal2.scaled = out$leaveout_diagonal2_scaled,
      joint.location1.scaled = out$joint_location1_scaled,
      joint.location2.scaled = out$joint_location2_scaled,
      weighted.location1.scaled = out$weighted_location1_scaled,
      weighted.location2.scaled = out$weighted_location2_scaled,
      own.direction1 = out$own_direction1,
      own.direction2 = out$own_direction2,
      own.radius1 = as.numeric(out$own_radius1),
      own.radius2 = as.numeric(out$own_radius2)
    ),
    diagnostics = list(
      iteration.stable = isTRUE(out$all_iteration_stable),
      stability.failures = as.integer(out$stability_failures),
      relative.update = worst.overall.relative.update,
      score.residual = worst.overall.score.residual,
      minimum.residual.distance = minimum.overall.residual.distance,
      worst.leaveout = worst.leaveout,
      convergence.basis = "relative update as specified by paper",
      score.residual.role = paste(
        "diagnostic only; it is not used as an estimating-equation root",
        "certificate"
      ),
      full.sample = list(
        group1 = out$full_fit1,
        group2 = out$full_fit2
      ),
      leaveout.joint = list(
        group1 = out$joint_fits1,
        group2 = out$joint_fits2
      ),
      leaveout.weighted.location = list(
        group1 = out$weighted_fits1,
        group2 = out$weighted_fits2
      ),
      inverse.norm.identity = list(
        group1 = isTRUE(all.equal(as.numeric(out$nu12),
                                  as.numeric(out$c1), tolerance = 0)),
        group2 = isTRUE(all.equal(as.numeric(out$nu22),
                                  as.numeric(out$c2), tolerance = 0))
      ),
      trace.pair.convention = "ordered i != j pairs",
      corrected.source.indices = c(
        "leaveout denominators use n_k",
        "diagonal iteration updates D_k",
        "A2 outer and inner indices use n2"
      ),
      rejection = list(
        alpha = alpha,
        critical.value = critical.value,
        rule = "reject when Z > qnorm(1 - alpha)",
        reject = reject
      ),
      cross.zero.vectors = as.numeric(out$cross_zero_vectors),
      internal.scaling = paste(
        "common per-coordinate translation and positive scaling across",
        "both groups"
      ),
      internal.anchor = stats::setNames(
        as.numeric(out$internal_anchor), variable.names
      ),
      internal.column.scale = stats::setNames(
        as.numeric(out$internal_column_scale), variable.names
      ),
      internal.overflow.fallback =
        isTRUE(out$internal_overflow_fallback),
      regularization = "none",
      zero.residual.perturbation = "none",
      variance.repair = "none"
    ),
    n = c(x = n1, y = n2),
    p = p,
    call = call
  )
}
