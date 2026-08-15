.ch3gp_validate_flag <- function(value, name) {
  if (!is.logical(value) || length(value) != 1L || is.na(value)) {
    stop(sprintf("`%s` must be TRUE or FALSE.", name), call. = FALSE)
  }
  value
}


.ch3gp_validate_number <- function(value, name, lower = 0,
                                   strict.lower = FALSE,
                                   upper = Inf) {
  value <- as.numeric(value)
  invalid.lower <- if (strict.lower) value <= lower else value < lower
  if (length(value) != 1L || is.na(value) || !is.finite(value) ||
      invalid.lower || value > upper) {
    interval <- if (strict.lower) "greater than" else "at least"
    stop(sprintf("`%s` must be one finite number %s %s%s.",
                 name, interval, format(lower),
                 if (is.finite(upper))
                   paste0(" and at most ", format(upper)) else ""),
         call. = FALSE)
  }
  value
}


.ch3gp_validate_count <- function(value, name) {
  if (!is.numeric(value) || length(value) != 1L || is.na(value) ||
      !is.finite(value) || value < 1 || value != floor(value) ||
      value > .Machine$integer.max) {
    stop(sprintf("`%s` must be one positive integer.", name),
         call. = FALSE)
  }
  as.integer(value)
}


.ch3gp_inputs <- function(x, center, divisor) {
  input <- .ch3g_threshold_inputs(x, center, divisor)
  sample.covariance <- .ch3g_restore_covariance_scale(
    input$covariance.scaled, input$scale, "sample covariance"
  )
  sample.covariance <- 0.5 *
    (sample.covariance + t(sample.covariance))
  if (any(diag(input$covariance.scaled) <= 0)) {
    stop(paste(
      "Every variable must have positive empirical variance; no zero-variance",
      "column was removed and no diagonal ridge was added."
    ), call. = FALSE)
  }
  variable.names <- colnames(input$x)
  if (is.null(variable.names)) {
    variable.names <- paste0("V", seq_len(input$p))
  }
  dimnames(sample.covariance) <- list(variable.names, variable.names)
  input$sample.covariance <- sample.covariance
  input$variable.names <- variable.names
  input
}


.ch3gp_stop_or_warn <- function(message, strict) {
  if (strict) {
    stop(message, call. = FALSE)
  }
  warning(message, call. = FALSE)
  invisible(NULL)
}


.ch3gp_precision_result <- function(estimate, valid, method, lambda,
                                    input, call, data.name, diagnostics) {
  if (!is.null(estimate)) {
    dimnames(estimate) <- list(input$variable.names, input$variable.names)
  }
  structure(
    list(
      estimate = estimate,
      valid = valid,
      method = method,
      lambda = lambda,
      sample.covariance = input$sample.covariance,
      center = input$center,
      data.name = data.name,
      diagnostics = diagnostics,
      call = call
    ),
    class = c("gaussian_precision_fit", "list")
  )
}


#' EC2 sparse covariance estimation with an eigenvalue constraint
#'
#' Implements the convex off-diagonal-lasso branch of Liu, Wang and Zhao's
#' EC2 estimator. The primary method first forms the empirical correlation
#' matrix \eqn{R}, solves
#' \deqn{\min_{C:\,\operatorname{diag}(C)=1,\,\lambda_{\min}(C)\geq\tau}
#' \frac12\|R-C\|_F^2+\lambda\sum_{i\ne j}|c_{ij}|,}
#' and restores marginal empirical standard deviations. This differs from the
#' book's shortened direct-covariance display. The default covariance divisor
#' is \eqn{n}; selecting `"n-1"` is explicit and changes only the marginal
#' covariance scale, not the sample correlation.
#'
#' The ISP/ADMM solver uses the paper's soft-threshold and spectral-projection
#' updates. A fit is valid only when the equality, fixed-point, off-diagonal
#' subgradient, spectral-dual, complementarity, diagonal, and eigenvalue
#' certificates all pass. The method-defining `tau` projection is not a
#' numerical repair. Adaptive EC2 and MC+ EC2 require additional weight or
#' shape choices and remain review-only; this function never guesses them.
#'
#' @param x Numeric observation-by-variable matrix or data frame.
#' @param lambda Finite non-negative off-diagonal lasso penalty on the
#'   correlation scale.
#' @param tau Finite minimum correlation-eigenvalue bound in \eqn{(0,1]}.
#' @param penalty Currently only `"l1"`. Adaptive and MC+ variants are
#'   deliberately review-only.
#' @param center Whether to subtract column means.
#' @param divisor Either `"n"` (formal default) or `"n-1"`.
#' @param rho Positive ISP/ADMM penalty multiplier.
#' @param solver_tol Positive tolerance for every feasibility and KKT
#'   certificate.
#' @param solver_max_iter Positive ISP/ADMM iteration limit.
#' @param strict If `TRUE`, a failed solver certificate is an error. If
#'   `FALSE`, the function warns and returns `estimate = NULL` with
#'   `valid = FALSE` and the uncertified last iterate in diagnostics.
#'
#' @return An `ec2_covariance_fit` list. `estimate` is the covariance matrix
#'   only for a certified fit; `correlation` is its EC2 correlation estimate.
#' @references Liu, H., Wang, L. and Zhao, T. (2014). Sparse covariance matrix
#'   estimation with eigenvalue constraints. *Journal of Computational and
#'   Graphical Statistics*, 23, 439--459.
#'   \doi{10.1080/10618600.2013.782818}.
#' @examples
#' x <- rbind(c(-2, -1), c(-1, 0), c(1, 0), c(2, 1))
#' ec2_covariance(x, lambda = 0.2, tau = 0.1)
#' @export
ec2_covariance <- function(
    x, lambda, tau, penalty = "l1", center = TRUE,
    divisor = c("n", "n-1"), rho = 1,
    solver_tol = 1e-7, solver_max_iter = 50000L,
    strict = TRUE) {
  call <- match.call()
  data.name <- deparse1(substitute(x))
  if (!is.character(penalty) || length(penalty) != 1L ||
      is.na(penalty) || penalty != "l1") {
    stop(paste(
      "Only `penalty = \"l1\"` is implemented as the uniquely specified",
      "convex EC2 branch. Adaptive and MC+ branches remain review-only."
    ), call. = FALSE)
  }
  lambda <- .ch3gp_validate_number(lambda, "lambda")
  tau <- .ch3gp_validate_number(
    tau, "tau", lower = 0, strict.lower = TRUE, upper = 1
  )
  rho <- .ch3gp_validate_number(
    rho, "rho", lower = 0, strict.lower = TRUE
  )
  solver_tol <- .ch3gp_validate_number(
    solver_tol, "solver_tol", lower = 0, strict.lower = TRUE
  )
  solver_max_iter <- .ch3gp_validate_count(
    solver_max_iter, "solver_max_iter"
  )
  strict <- .ch3gp_validate_flag(strict, "strict")
  input <- .ch3gp_inputs(x, center, divisor)

  marginal.sd.scaled <- sqrt(diag(input$covariance.scaled))
  correlation <- input$covariance.scaled /
    outer(marginal.sd.scaled, marginal.sd.scaled)
  correlation <- 0.5 * (correlation + t(correlation))
  diag(correlation) <- 1
  dimnames(correlation) <- list(input$variable.names, input$variable.names)

  core <- tryCatch(
    cpp_ch3gp_ec2_l1(
      correlation, lambda, tau, rho, solver_tol, solver_max_iter
    ),
    error = identity
  )
  if (inherits(core, "condition")) {
    message <- paste0(
      "EC2 optimization failed: ", conditionMessage(core),
      " No covariance estimate was returned and no repair was applied."
    )
    .ch3gp_stop_or_warn(message, strict)
    return(structure(
      list(
        estimate = NULL, correlation = NULL, valid = FALSE,
        method = "EC2 covariance (off-diagonal l1)",
        lambda = lambda, tau = tau,
        sample.covariance = input$sample.covariance,
        sample.correlation = correlation,
        center = input$center, data.name = data.name,
        diagnostics = list(
          failure.stage = "optimization error",
          error = conditionMessage(core),
          covariance.divisor = input$divisor$value,
          divisor.convention = input$divisor$label,
          no.repair = TRUE
        ),
        call = call
      ),
      class = c("ec2_covariance_fit", "list")
    ))
  }

  valid <- isTRUE(core$converged) && is.finite(core$kkt_maximum) &&
    core$kkt_maximum <= solver_tol &&
    is.finite(core$minimum_eigenvalue) &&
    core$minimum_eigenvalue > 0 &&
    core$eigenvalue_feasibility_violation <= solver_tol &&
    core$diagonal_residual <= solver_tol &&
    core$symmetry_residual <= solver_tol
  correlation.estimate <- as.matrix(core$solution)
  dimnames(correlation.estimate) <-
    list(input$variable.names, input$variable.names)
  estimate <- NULL
  if (valid) {
    estimate.scaled <- correlation.estimate *
      outer(marginal.sd.scaled, marginal.sd.scaled)
    estimate <- .ch3g_restore_covariance_scale(
      estimate.scaled, input$scale, "EC2 covariance estimate"
    )
    estimate <- 0.5 * (estimate + t(estimate))
    dimnames(estimate) <- list(input$variable.names, input$variable.names)
  } else {
    message <- paste0(
      "EC2 stopped without all feasibility and KKT certificates (maximum ",
      format(core$kkt_maximum, digits = 6), "). No covariance estimate was ",
      "returned and no iterate was repaired."
    )
    .ch3gp_stop_or_warn(message, strict)
  }

  structure(
    list(
      estimate = estimate,
      correlation = if (valid) correlation.estimate else NULL,
      valid = valid,
      method = "EC2 covariance (off-diagonal l1)",
      lambda = lambda,
      tau = tau,
      sample.covariance = input$sample.covariance,
      sample.correlation = correlation,
      center = input$center,
      data.name = data.name,
      diagnostics = list(
        failure.stage = if (valid) NULL else "optimization certificate",
        solver = list(
          iterations = as.integer(core$iterations),
          converged = isTRUE(core$converged),
          relative.update = as.numeric(core$relative_update),
          dual.residual = as.numeric(core$dual_residual),
          objective = as.numeric(core$objective),
          objective.history = as.numeric(core$objective_history),
          soft.threshold.feasible = isTRUE(core$sto_feasible),
          primal.equality.residual =
            as.numeric(core$primal_equality_residual),
          sparse.fixed.point.residual =
            as.numeric(core$sparse_fixed_point_residual),
          spectral.fixed.point.residual =
            as.numeric(core$spectral_fixed_point_residual),
          diagonal.residual = as.numeric(core$diagonal_residual),
          symmetry.residual = as.numeric(core$symmetry_residual),
          minimum.eigenvalue = as.numeric(core$minimum_eigenvalue),
          eigenvalue.feasibility.violation =
            as.numeric(core$eigenvalue_feasibility_violation),
          off.diagonal.stationarity.residual =
            as.numeric(core$off_diagonal_stationarity_residual),
          spectral.dual.violation =
            as.numeric(core$spectral_dual_violation),
          spectral.complementarity.residual =
            as.numeric(core$spectral_complementarity_residual),
          kkt.maximum = as.numeric(core$kkt_maximum)
        ),
        last.iterate = if (valid) NULL else list(
          sparse = as.matrix(core$solution),
          spectral = as.matrix(core$spectral_solution),
          dual = as.matrix(core$dual),
          spectral.multiplier = as.matrix(core$spectral_multiplier)
        ),
        covariance.divisor = input$divisor$value,
        divisor.convention = input$divisor$label,
        formulation = paste(
          "primary correlation EC2 with unit diagonal, followed by",
          "marginal-standard-deviation covariance reconstruction"
        ),
        book.direct.covariance.display = "not substituted for primary EC2",
        review.only.penalties = c("adaptive EC2", "MC+ EC2"),
        internal.data.scale = input$scale,
        positive.definiteness.repair = "none; tau is the model constraint",
        no.repair = paste(
          "No ridge, pseudoinverse, tolerance relaxation, or post-hoc",
          "eigenvalue/KKT repair"
        )
      ),
      call = call
    ),
    class = c("ec2_covariance_fit", "list")
  )
}


#' Gaussian graphical lasso with an off-diagonal penalty
#'
#' With the centered empirical covariance \eqn{S_n}, solves the Yuan--Lin
#' off-diagonal graphical-lasso program
#' \deqn{\min_{\Omega\succ0}
#' \operatorname{tr}(S_n\Omega)-\log\det(\Omega)
#' +\lambda\sum_{i\ne j}|\omega_{ij}|.}
#' Diagonal entries are never penalized. The optimizer uses an SPD-preserving
#' proximal-gradient step with explicit majorization backtracking and returns
#' the full diagonal/off-diagonal subgradient KKT residual. With `lambda = 0`,
#' the exact inverse is returned only when \eqn{S_n} is strictly positive
#' definite; no pseudoinverse or ridge is substituted.
#'
#' @inheritParams ec2_covariance
#' @param initial_step Positive initial proximal-gradient step.
#' @param max_backtracking Positive line-search reduction limit per iteration.
#'
#' @return A `gaussian_precision_fit` list. `estimate` is non-`NULL` only when
#'   SPD, objective-descent, relative-update, and KKT certificates all pass.
#' @references Yuan, M. and Lin, Y. (2007). Model selection and estimation in
#'   the Gaussian graphical model. *Biometrika*, 94, 19--35.
#'   \doi{10.1093/biomet/asm018}.
#' @examples
#' x <- rbind(c(-2, 0), c(-1, -1), c(1, 1), c(2, 0))
#' gaussian_graphical_lasso(x, lambda = 0.2)
#' @export
gaussian_graphical_lasso <- function(
    x, lambda, center = TRUE, divisor = c("n", "n-1"),
    solver_tol = 1e-7, solver_max_iter = 100000L,
    initial_step = 1, max_backtracking = 100L,
    strict = TRUE) {
  call <- match.call()
  data.name <- deparse1(substitute(x))
  lambda <- .ch3gp_validate_number(lambda, "lambda")
  solver_tol <- .ch3gp_validate_number(
    solver_tol, "solver_tol", lower = 0, strict.lower = TRUE
  )
  solver_max_iter <- .ch3gp_validate_count(
    solver_max_iter, "solver_max_iter"
  )
  initial_step <- .ch3gp_validate_number(
    initial_step, "initial_step", lower = 0, strict.lower = TRUE
  )
  max_backtracking <- .ch3gp_validate_count(
    max_backtracking, "max_backtracking"
  )
  strict <- .ch3gp_validate_flag(strict, "strict")
  input <- .ch3gp_inputs(x, center, divisor)

  core <- tryCatch(
    cpp_ch3gp_offdiag_glasso(
      input$sample.covariance, lambda, solver_tol,
      solver_max_iter, initial_step, max_backtracking
    ),
    error = identity
  )
  if (inherits(core, "condition")) {
    message <- paste0(
      "Gaussian graphical-lasso optimization failed: ",
      conditionMessage(core),
      " No precision estimate was returned and no repair was applied."
    )
    .ch3gp_stop_or_warn(message, strict)
    return(.ch3gp_precision_result(
      NULL, FALSE, "Gaussian graphical lasso (off-diagonal l1)",
      lambda, input, call, data.name,
      list(
        failure.stage = "optimization error",
        error = conditionMessage(core),
        covariance.divisor = input$divisor$value,
        divisor.convention = input$divisor$label,
        no.repair = TRUE
      )
    ))
  }

  valid <- isTRUE(core$converged) &&
    !isTRUE(core$backtracking_failed) &&
    isTRUE(core$objective_descent) &&
    is.finite(core$minimum_eigenvalue) &&
    core$minimum_eigenvalue > 0 &&
    is.finite(core$kkt_maximum) &&
    core$kkt_maximum <= solver_tol &&
    is.finite(core$relative_update) &&
    core$relative_update <= solver_tol
  estimate <- if (valid) as.matrix(core$solution) else NULL
  if (!valid) {
    message <- paste0(
      "Gaussian graphical lasso stopped without every SPD, descent, and KKT ",
      "certificate (maximum KKT residual ",
      format(core$kkt_maximum, digits = 6),
      "). No precision estimate was returned and no iterate was repaired."
    )
    .ch3gp_stop_or_warn(message, strict)
  }

  .ch3gp_precision_result(
    estimate, valid, "Gaussian graphical lasso (off-diagonal l1)",
    lambda, input, call, data.name,
    list(
      failure.stage = if (valid) NULL else "optimization certificate",
      solver = list(
        iterations = as.integer(core$iterations),
        converged = isTRUE(core$converged),
        backtracking.failed = isTRUE(core$backtracking_failed),
        objective.descent = isTRUE(core$objective_descent),
        relative.update = as.numeric(core$relative_update),
        objective = as.numeric(core$objective),
        objective.history = as.numeric(core$objective_history),
        diagonal.kkt.residual = as.numeric(core$kkt_diagonal),
        off.diagonal.kkt.residual =
          as.numeric(core$kkt_off_diagonal),
        kkt.maximum = as.numeric(core$kkt_maximum),
        minimum.eigenvalue = as.numeric(core$minimum_eigenvalue),
        accepted.step = as.numeric(core$accepted_step),
        total.backtracking = as.integer(core$total_backtracking)
      ),
      last.iterate = if (valid) NULL else list(
        precision = as.matrix(core$solution),
        inverse = as.matrix(core$inverse)
      ),
      covariance.divisor = input$divisor$value,
      divisor.convention = input$divisor$label,
      diagonal.penalty = FALSE,
      internal.data.scale = input$scale,
      positive.definiteness.repair = "none",
      no.repair = paste(
        "No ridge, pseudoinverse, eigenvalue floor, tolerance relaxation,",
        "or post-hoc KKT repair"
      )
    )
  )
}


#' Gaussian CLIME sparse precision estimator
#'
#' Solves the Cai--Liu--Luo CLIME program column by column,
#' \deqn{\min_b\|b\|_1\quad\text{subject to}\quad
#' \|S_n b-e_j\|_\infty\leq\lambda,}
#' using a primal--dual algorithm. Every raw column must pass primal
#' feasibility, dual feasibility, lasso stationarity, and relative duality-gap
#' certificates. The final symmetric estimator retains the entry of smaller
#' absolute value from each transposed pair; exact ties retain the row-column
#' entry before mirroring, making the convention deterministic.
#'
#' Primary CLIME does not guarantee that this post-LP symmetrization remains
#' feasible, nor that the final matrix is positive definite. Both facts are
#' reported and neither is repaired. Consequently `valid` certifies the raw
#' column programs, not an invented SPD condition. SCIO and scaled-lasso
#' precision estimation remain review-only because their additional programs
#' are not specified by the short book display.
#'
#' @inheritParams ec2_covariance
#'
#' @return A `gaussian_precision_fit` list. `estimate` is the primary
#'   smaller-absolute-value symmetrization of the certified raw columns.
#' @references Cai, T. T., Liu, W. and Luo, X. (2011). A constrained
#'   \eqn{\ell_1} minimization approach to sparse precision matrix estimation.
#'   *Journal of the American Statistical Association*, 106, 594--607.
#'   \doi{10.1198/jasa.2011.tm10155}.
#' @examples
#' x <- rbind(c(-2, 0), c(-1, -1), c(1, 1), c(2, 0))
#' clime_precision(x, lambda = 0)
#' @export
clime_precision <- function(
    x, lambda, center = TRUE, divisor = c("n", "n-1"),
    solver_tol = 1e-7, solver_max_iter = 100000L,
    strict = TRUE) {
  call <- match.call()
  data.name <- deparse1(substitute(x))
  lambda <- .ch3gp_validate_number(lambda, "lambda")
  solver_tol <- .ch3gp_validate_number(
    solver_tol, "solver_tol", lower = 0, strict.lower = TRUE
  )
  solver_max_iter <- .ch3gp_validate_count(
    solver_max_iter, "solver_max_iter"
  )
  strict <- .ch3gp_validate_flag(strict, "strict")
  input <- .ch3gp_inputs(x, center, divisor)

  core <- tryCatch(
    cpp_ch3gp_clime(
      input$sample.covariance, lambda, solver_tol, solver_max_iter
    ),
    error = identity
  )
  if (inherits(core, "condition")) {
    message <- paste0(
      "Gaussian CLIME optimization failed: ", conditionMessage(core),
      " No precision estimate was returned and no constraint was repaired."
    )
    .ch3gp_stop_or_warn(message, strict)
    return(.ch3gp_precision_result(
      NULL, FALSE, "Gaussian CLIME", lambda, input, call, data.name,
      list(
        failure.stage = "optimization error",
        error = conditionMessage(core),
        covariance.divisor = input$divisor$value,
        divisor.convention = input$divisor$label,
        no.repair = TRUE
      )
    ))
  }

  valid <- isTRUE(core$all_converged) &&
    all(is.finite(core$primal_violation)) &&
    all(core$primal_violation <= solver_tol) &&
    all(is.finite(core$dual_violation)) &&
    all(core$dual_violation <= solver_tol) &&
    all(is.finite(core$stationarity_residual)) &&
    all(core$stationarity_residual <= solver_tol) &&
    all(is.finite(core$relative_gap)) &&
    all(core$relative_gap <= solver_tol)
  estimate <- if (valid) as.matrix(core$solution) else NULL
  if (!valid) {
    message <- paste0(
      "Gaussian CLIME stopped without every raw-column primal-dual ",
      "certificate. No precision estimate was returned and no column was ",
      "repaired."
    )
    .ch3gp_stop_or_warn(message, strict)
  }

  raw <- as.matrix(core$raw_solution)
  dual <- as.matrix(core$dual_solution)
  dimnames(raw) <- dimnames(dual) <-
    list(input$variable.names, input$variable.names)
  .ch3gp_precision_result(
    estimate, valid, "Gaussian CLIME", lambda, input, call, data.name,
    list(
      failure.stage = if (valid) NULL else "optimization certificate",
      solver = list(
        all.columns.certified = valid,
        column.converged = as.logical(core$converged),
        iterations = as.integer(core$iterations),
        relative.update = as.numeric(core$relative_update),
        primal.violation = as.numeric(core$primal_violation),
        stationarity.residual = as.numeric(core$stationarity_residual),
        dual.violation = as.numeric(core$dual_violation),
        primal.objective = as.numeric(core$primal_objective),
        dual.objective = as.numeric(core$dual_objective),
        duality.gap = as.numeric(core$duality_gap),
        relative.gap = as.numeric(core$relative_gap),
        operator.norm = as.numeric(core$operator_norm),
        primal.dual.step = as.numeric(core$primal_dual_step)
      ),
      raw.solution = raw,
      dual.solution = dual,
      last.symmetrized.iterate = if (valid) NULL else
        as.matrix(core$solution),
      symmetrization = paste(
        "retain the smaller-absolute-value transposed entry;",
        "ties retain raw[i,j] before mirroring"
      ),
      symmetrized.feasibility.violation =
        as.numeric(core$symmetrized_feasibility_violation),
      symmetrized.constraint.feasible =
        is.finite(core$symmetrized_feasibility_violation) &&
        core$symmetrized_feasibility_violation <= solver_tol,
      symmetrized.minimum.eigenvalue =
        as.numeric(core$symmetrized_minimum_eigenvalue),
      symmetrized.positive.definite =
        is.finite(core$symmetrized_minimum_eigenvalue) &&
        core$symmetrized_minimum_eigenvalue > 0,
      covariance.divisor = input$divisor$value,
      divisor.convention = input$divisor$label,
      review.only.alternatives = c("SCIO", "scaled-lasso precision"),
      internal.data.scale = input$scale,
      positive.definiteness.repair = "none; not guaranteed by CLIME",
      no.repair = paste(
        "No ridge, pseudoinverse, feasibility relaxation, forced SPD,",
        "or post-hoc KKT repair"
      )
    )
  )
}
