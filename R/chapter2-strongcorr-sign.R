.zfsc_validate_logical <- function(value, name) {
  if (!is.logical(value) || length(value) != 1L || is.na(value)) {
    stop(sprintf("`%s` must be TRUE or FALSE.", name), call. = FALSE)
  }
  value
}

.zfsc_validate_alpha <- function(alpha) {
  alpha <- as.numeric(alpha)
  if (length(alpha) != 1L || is.na(alpha) || !is.finite(alpha) ||
      alpha <= 0 || alpha >= 1) {
    stop("`alpha` must be one finite number strictly between zero and one.",
         call. = FALSE)
  }
  alpha
}

.zfsc_validate_B <- function(B) {
  if (!is.numeric(B) || length(B) != 1L || is.na(B) || !is.finite(B) ||
      B < 1 || B != floor(B) || B > .Machine$integer.max) {
    stop("`B` must be a positive integer no larger than R's integer limit.",
         call. = FALSE)
  }
  as.integer(B)
}

.zfsc_validate_seed <- function(seed) {
  if (is.null(seed)) return(NULL)
  if (!is.numeric(seed) || length(seed) != 1L || is.na(seed) ||
      !is.finite(seed) || seed < 0 || seed != floor(seed) ||
      seed > .Machine$integer.max) {
    stop(
      "`seed` must be NULL or one integer from zero through R's integer limit.",
      call. = FALSE
    )
  }
  as.integer(seed)
}

.zfsc_with_local_seed <- function(seed, expression) {
  had.state <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had.state) old.state <- get(".Random.seed", envir = .GlobalEnv)
  on.exit({
    if (had.state) {
      assign(".Random.seed", old.state, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  set.seed(seed)
  force(expression)
}

.zfsc_spatial_median_diagnostics <- function(location, tol, max_iter,
                                               strict) {
  list(
    converged = isTRUE(attr(location, "converged")),
    iterations = as.integer(attr(location, "iterations")),
    tolerance = tol,
    max.iterations.allowed = max_iter,
    objective = as.numeric(attr(location, "objective")),
    relative.change = as.numeric(attr(location, "relative_change")),
    equation.residual = as.numeric(attr(location, "equation_residual")),
    convergence.basis = paste(
      "modified Weiszfeld subgradient-equation convergence as returned by",
      "spatial_median()"
    ),
    strict = strict,
    zero.tolerance = 0,
    regularization = "none",
    numerical.floor = "none",
    perturbation = "none",
    ridge = "none"
  )
}

#' Zhao--Feng strong-correlation spatial-sign wild-bootstrap test
#'
#' Tests a one-sample location null with the wild-bootstrap calibration of
#' Zhao and Feng.  The observed raw statistic is
#' \deqn{S_n=\sum_{1\leq i<j\leq n}U(X_i-\mu_0)^\mathsf{T}
#' U(X_j-\mu_0),}
#' where \eqn{U(0)=0}.  Bootstrap signs are deliberately fitted differently:
#' if \eqn{\widehat\mu} is the ordinary Euclidean sample spatial median and
#' \eqn{\widehat U_i=U(X_i-\widehat\mu)}, then a replicate is
#' \deqn{S^*=\sum_{i<j}e_i e_j\widehat U_i^\mathsf{T}\widehat U_j.}
#' The multipliers are either independent Rademacher variables or independent
#' standard Gaussian variables.  The spatial median is fitted once and is not
#' refitted inside the bootstrap.
#'
#' The paper writes both statistics divided by
#' \eqn{\sqrt{\tau}\sqrt{\binom{n}{2}}}, where
#' \eqn{\tau=\operatorname{tr}(\Sigma_U^2)}, but explicitly notes that
#' \eqn{\tau} need not be estimated because this common factor cancels from
#' the bootstrap comparison.  Accordingly, the test uses raw pair sums and
#' also reports their \eqn{\sqrt{\binom{n}{2}}}-normalised, tau-free versions.
#' It does not invent a plug-in estimate of \eqn{\tau}.
#'
#' The paper prescribes the empirical \eqn{(1-\alpha)} quantile and rejects
#' when the observed statistic is strictly greater.  For reproducibility this
#' implementation defines that quantile as the type-1 inverse empirical cdf:
#' order statistic \eqn{\lceil(1-\alpha)B\rceil}.  The ordinary `p.value`
#' uses the separately labelled finite-Monte-Carlo plus-one convention
#' \deqn{\{1+\#(S_b^*\geq S_n)\}/(B+1),}
#' with ties counted in the upper tail.  Both rejection decisions are returned
#' because they can differ at finite `B`.
#'
#' The pair-sum identity used computationally is
#' \eqn{\{\|\sum_i U_i\|^2-\sum_i\|U_i\|^2\}/2}.  Thus an observation exactly
#' equal to a centre remains the literal zero sign; the implementation does
#' not silently replace the diagonal term by \eqn{n}, add jitter, or repair a
#' degenerate bootstrap distribution.  A degenerate distribution is returned
#' with an explicit diagnostic.
#'
#' The test is invariant to a common translation (when `mu` is translated),
#' orthogonal transformations, and a common nonzero scalar transformation.
#' It is not coordinatewise scale invariant.  `strict = TRUE` makes failure
#' of the ordinary spatial-median iteration an error.  With `strict = FALSE`,
#' the last finite iterate is used with a warning and its diagnostics are
#' retained.  No ridge, floor, perturbation, or pseudoinverse is applied.
#'
#' @param x A finite numeric matrix or data frame with observations in rows.
#'   At least two observations and one variable are required.
#' @param mu A finite null-location vector.  The default is the zero vector.
#' @param alpha A finite test level strictly between zero and one.
#' @param multiplier Either `"rademacher"` or `"gaussian"`.
#' @param B A positive integer number of wild-bootstrap replicates.
#' @param seed `NULL` or an integer from zero through R's integer limit.  An
#'   explicit seed makes the bootstrap reproducible without changing the
#'   caller's RNG state.  With `NULL`, one seed is first drawn from R's RNG and
#'   recorded, after which the replicate generation is isolated locally.
#' @param keep_bootstrap Whether to retain all raw and pair-normalised
#'   bootstrap statistics.
#' @param tol A finite positive tolerance passed to [spatial_median()].
#' @param max_iter A positive integer iteration limit for the spatial median.
#' @param strict Whether failure of the spatial-median iteration is an error
#'   (`TRUE`) or a warning followed by use of its last finite iterate
#'   (`FALSE`).
#'
#' @return An object of class `c("hd_location_test", "htest")`.  Its
#'   `p.value` is the auxiliary plus-one Monte Carlo upper-tail probability.
#'   `components` contains the raw and tau-free statistics, the paper
#'   type-1 critical rule, both decisions, null-centred and fitted signs, the
#'   fitted ordinary spatial median, and optional bootstrap draws.
#'   `diagnostics` records convergence, zero residuals, overflow-safe residual
#'   subtraction, Monte Carlo conventions, and the no-repair contract.
#'
#' @references
#' Zhao, P. and Feng, L. (2026). Note on High Dimensional Spatial-Sign Test
#' for One Sample Problem. arXiv:2601.08736.
#' \doi{10.48550/arXiv.2601.08736}.
#'
#' @examples
#' x <- rbind(c(-1.0, 0.4), c(0.2, -0.8), c(1.3, 0.6),
#'            c(-0.4, 1.1), c(0.8, -0.2))
#' zhao_feng_strongcorr_sign_test(x, B = 99, seed = 2601)
#'
#' @export
zhao_feng_strongcorr_sign_test <- function(
    x, mu = NULL, alpha = 0.05,
    multiplier = c("rademacher", "gaussian"), B = 9999L,
    seed = NULL, keep_bootstrap = FALSE, tol = 1e-8,
    max_iter = 1000L, strict = TRUE) {
  call <- match.call()
  data.name <- deparse1(substitute(x))
  x <- .as_data_matrix(x, "x", min_rows = 2L)
  n <- nrow(x)
  p <- ncol(x)
  variable.names <- .hotelling_variable_names(x)
  observation.names <- if (is.null(rownames(x))) {
    paste0("observation", seq_len(n))
  } else {
    rownames(x)
  }
  if (is.null(mu)) {
    mu <- numeric(p)
  } else {
    mu <- .as_location(mu, p, "mu")
  }
  names(mu) <- variable.names

  alpha <- .zfsc_validate_alpha(alpha)
  multiplier <- match.arg(multiplier)
  B <- .zfsc_validate_B(B)
  seed.requested <- .zfsc_validate_seed(seed)
  keep_bootstrap <- .zfsc_validate_logical(
    keep_bootstrap, "keep_bootstrap"
  )
  strict <- .zfsc_validate_logical(strict, "strict")
  controls <- .validate_iteration_controls(tol, max_iter, zero_tol = 0)

  fitted.location <- spatial_median(
    x, tol = controls$tol, max_iter = controls$max_iter,
    zero_tol = 0, warn = FALSE
  )
  median.diagnostics <- .zfsc_spatial_median_diagnostics(
    fitted.location, controls$tol, controls$max_iter, strict
  )
  if (!isTRUE(median.diagnostics$converged)) {
    message <- sprintf(
      paste0(
        "The ordinary spatial-median iteration did not converge within ",
        "%d iterations (equation residual %.6g, relative change %.6g)."
      ),
      controls$max_iter, median.diagnostics$equation.residual,
      median.diagnostics$relative.change
    )
    if (strict) stop(message, call. = FALSE)
    warning(paste0(message, " Returning the last finite iterate."),
            call. = FALSE)
  }
  fitted.location <- stats::setNames(
    as.numeric(fitted.location), variable.names
  )
  if (anyNA(fitted.location) || any(!is.finite(fitted.location))) {
    stop(
      paste0(
        "The ordinary spatial median is non-finite; no ridge, truncation, ",
        "or perturbation is applied."
      ),
      call. = FALSE
    )
  }

  seed.used <- if (is.null(seed.requested)) {
    as.integer(sample.int(.Machine$integer.max, 1L) - 1L)
  } else {
    seed.requested
  }
  fit <- .zfsc_with_local_seed(
    seed.used,
    cpp_zhao_feng_strongcorr_sign_bootstrap(
      x, as.numeric(mu), as.numeric(fitted.location), multiplier,
      B, alpha, keep_bootstrap
    )
  )

  observed.signs <- as.matrix(fit$observed_signs)
  fitted.signs <- as.matrix(fit$fitted_signs)
  dimnames(observed.signs) <- list(observation.names, variable.names)
  dimnames(fitted.signs) <- list(observation.names, variable.names)
  name.observations <- function(value) {
    stats::setNames(as.numeric(value), observation.names)
  }
  name.variables <- function(value) {
    stats::setNames(as.numeric(value), variable.names)
  }

  critical.reject <- isTRUE(fit$reject_by_critical)
  pvalue.reject <- isTRUE(fit$reject_by_p_value)
  randomization.degenerate <- isTRUE(fit$randomization_degenerate)
  bootstrap.raw <- if (keep_bootstrap) {
    as.numeric(fit$bootstrap_statistics_raw)
  } else {
    NULL
  }
  bootstrap.normalized <- if (keep_bootstrap) {
    as.numeric(fit$bootstrap_statistics_pair_normalized)
  } else {
    NULL
  }

  .new_hd_location_test(
    statistic = c(S.n = as.numeric(fit$observed_raw)),
    p.value = as.numeric(fit$p_value_plus_one),
    method = paste(
      "Zhao-Feng strong-correlation one-sample spatial-sign test",
      sprintf("(%s wild bootstrap)", multiplier)
    ),
    data.name = data.name,
    alternative = "two.sided",
    raw.statistic = c(S.n = as.numeric(fit$observed_raw)),
    estimate = fitted.location,
    null.value = mu,
    null.distribution = list(
      family = "wild-bootstrap empirical distribution",
      exact = FALSE,
      tail = "upper",
      multiplier = multiplier,
      replicates = B,
      paper.critical.quantile = paste0(
        "type-1 inverse ECDF order statistic ceiling((1-alpha)*B) = ",
        as.integer(fit$critical_index)
      ),
      paper.critical.value.raw = as.numeric(fit$critical_value_raw),
      paper.rejection = "observed raw pair sum strictly greater than critical",
      p.value.convention = paste(
        "auxiliary plus-one Monte Carlo upper tail; ties counted with >="
      ),
      theoretical.common.scale =
        "sqrt(tau) * sqrt(choose(n,2)); omitted from both sides",
      assumptions = paste(
        "Zhao-Feng elliptical-symmetry, radial-moment, spatial-median",
        "Bahadur, trace-regularity, and n*tau divergence conditions"
      )
    ),
    components = list(
      observed.raw = as.numeric(fit$observed_raw),
      observed.pair.normalized.tau.free =
        as.numeric(fit$observed_pair_normalized),
      pair.count = as.numeric(fit$pair_count),
      root.pair.count = as.numeric(fit$root_pair_count),
      tau.estimated = FALSE,
      tau.factor.used = FALSE,
      observed.signs.null.centered = observed.signs,
      observed.sign.sum = name.variables(fit$observed_sign_sum),
      observed.sign.diagonal = as.numeric(fit$observed_sign_diagonal),
      observed.residual.norms = name.observations(
        fit$observed_residual_norms
      ),
      observed.log.residual.norms = name.observations(
        fit$observed_log_residual_norms
      ),
      fitted.spatial.median = fitted.location,
      bootstrap.signs.median.centered = fitted.signs,
      fitted.sign.sum = name.variables(fit$fitted_sign_sum),
      fitted.sign.diagonal = as.numeric(fit$fitted_sign_diagonal),
      fitted.raw.pair.sum = as.numeric(fit$fitted_raw_pair_sum),
      fitted.residual.norms = name.observations(
        fit$fitted_residual_norms
      ),
      fitted.log.residual.norms = name.observations(
        fit$fitted_log_residual_norms
      ),
      bootstrap.statistics.raw = bootstrap.raw,
      bootstrap.statistics.pair.normalized.tau.free =
        bootstrap.normalized,
      bootstrap.summary.raw = c(
        mean = as.numeric(fit$bootstrap_mean_raw),
        variance.population = as.numeric(fit$bootstrap_variance_raw),
        minimum = as.numeric(fit$bootstrap_minimum_raw),
        maximum = as.numeric(fit$bootstrap_maximum_raw)
      ),
      critical.value.raw = as.numeric(fit$critical_value_raw),
      critical.value.pair.normalized.tau.free =
        as.numeric(fit$critical_value_pair_normalized),
      critical.order.index = as.integer(fit$critical_index),
      reject.paper.critical = critical.reject,
      p.value.plus.one = as.numeric(fit$p_value_plus_one),
      empirical.tail.without.plus.one = as.numeric(fit$empirical_tail),
      exceedances.including.ties = as.numeric(fit$exceedances),
      reject.plus.one.p.value = pvalue.reject,
      decisions.differ = isTRUE(fit$decisions_differ),
      monte.carlo.standard.error.unadjusted.tail =
        as.numeric(fit$mc_standard_error),
      multiplier = multiplier,
      B = B,
      alpha = alpha,
      seed.requested = seed.requested,
      seed.used = seed.used,
      n = n,
      p = p
    ),
    diagnostics = list(
      spatial.median = median.diagnostics,
      centering = list(
        observed = "U(X_i - mu0), centered at the null location",
        bootstrap = paste(
          "U(X_i - fitted ordinary spatial median); fitted once before",
          "bootstrap and not refitted per replicate"
        )
      ),
      scale = list(
        tau.estimated = FALSE,
        comparison = paste(
          "raw sums; the common sqrt(tau)*sqrt(choose(n,2)) factor",
          "cancels exactly from observed-versus-critical comparison"
        ),
        pair.normalized.reported = TRUE
      ),
      rejection = list(
        alpha = alpha,
        paper.wording = paste(
          "use their empirical (1-alpha)-quantiles as critical values;",
          "rejected whenever the observed test statistic exceeds the",
          "corresponding bootstrap critical value"
        ),
        paper.rule = "strict observed > type-1 empirical critical value",
        paper.critical.value.raw = as.numeric(fit$critical_value_raw),
        paper.reject = critical.reject,
        auxiliary.p.value.rule =
          "(1 + count(bootstrap >= observed)) / (B + 1)",
        auxiliary.p.value = as.numeric(fit$p_value_plus_one),
        auxiliary.reject = pvalue.reject,
        minimum.attainable.plus.one.p = 1 / (B + 1),
        finite.B.decisions.can.differ = TRUE,
        decisions.differ.here = isTRUE(fit$decisions_differ)
      ),
      randomization = list(
        multiplier = multiplier,
        B = B,
        seed.requested = seed.requested,
        seed.used = seed.used,
        explicit.seed.preserves.caller.RNG.state = !is.null(seed.requested),
        null.seed.draws.one.recorded.seed.from.caller.RNG =
          is.null(seed.requested),
        ties.counted.as.exceedances = TRUE,
        degenerate = randomization.degenerate
      ),
      observed.zero.residuals = as.numeric(
        fit$observed_zero_residuals
      ),
      fitted.zero.residuals = as.numeric(fit$fitted_zero_residuals),
      observed.overflow.fallback.rows = as.numeric(
        fit$observed_overflow_fallback_rows
      ),
      fitted.overflow.fallback.rows = as.numeric(
        fit$fitted_overflow_fallback_rows
      ),
      pair.sum.identity = paste(
        "0.5*(squared norm of summed signs - sum of squared sign norms);",
        "zero signs are retained literally"
      ),
      invariance = paste(
        "common translation with translated null, orthogonal transforms,",
        "and common nonzero scalar; not coordinatewise scale invariant"
      ),
      primary.source.notes = c(
        kappa4 = paste(
          "Eq. (2.7) defines inner-product fourth-moment kurtosis while",
          "Assumption 3.3 writes a trace-fourth ratio"
        ),
        assumption3.2 =
          "the matrix S in the displayed assumption is not defined there",
        zero.sign = paste(
          "the paper defines U(0)=0 although some proof identities assume",
          "all fitted signs have unit norm"
        )
      ),
      regularization = "none",
      numerical.floor = "none",
      absolute.value.repair = "none",
      perturbation = "none",
      pseudoinverse = "none"
    ),
    n = n,
    p = p,
    call = call
  )
}
