#' Feng--Zhang--Liu high-dimensional two-sample spatial-rank test
#'
#' Tests equality of two high-dimensional elliptical location vectors using
#' the leave-two-out spatial-rank statistic of Feng, Zhang, and Liu (2020).
#' Observations are rows and variables are columns.  For distinct
#' \eqn{i,j} in group 1 and \eqn{s,t} in group 2, let
#' \eqn{\widehat D_{(i,j,s,t)}} be the weighted pool of the two diagonal
#' spatial-rank scale fits after deleting \eqn{i,j} and \eqn{s,t},
#' respectively.  The raw statistic is
#' \deqn{T_n=\{n_1(n_1-1)n_2(n_2-1)\}^{-1}
#' \sum_{i\ne j}\sum_{s\ne t}
#' U\{\widehat D_{(i,j,s,t)}^{-1/2}(X_{1i}-X_{2s})\}^{\mathsf T}
#' U\{\widehat D_{(i,j,s,t)}^{-1/2}(X_{1j}-X_{2t})\}.}
#' Here \eqn{U(v)=v/\lVert v\rVert} for nonzero \eqn{v}, and
#' \eqn{U(0)=0}.
#'
#' The p-value uses the feasible null calibration in the original article.
#' Its two within-group trace estimates use all ordered quadruples of
#' mutually distinct observations and their leave-four-out diagonal scale
#' fits.  The cross trace uses all ordered within-group pairs, pooled
#' leave-two-out scales, and the published denominator \eqn{n_1^2n_2^2}.
#' If these estimates are denoted by
#' \eqn{\widehat{\operatorname{tr}(R_1^2)}},
#' \eqn{\widehat{\operatorname{tr}(R_2^2)}}, and
#' \eqn{\widehat{\operatorname{tr}(R_1R_2)}}, then
#' \deqn{\widehat\sigma_n^2=
#' \frac{\widehat{\operatorname{tr}(R_1^2)}}
#' {2n_1(n_1-1)p^2}+
#' \frac{\widehat{\operatorname{tr}(R_2^2)}}
#' {2n_2(n_2-1)p^2}+
#' \frac{\widehat{\operatorname{tr}(R_1R_2)}}{n_1n_2p^2}.}
#' The \eqn{p^2} in the second term restores an evident typographical
#' omission in the author TeX: it is required by the symmetric oracle
#' variance and by consistency of the feasible estimator.  The reported
#' statistic is \eqn{Z=T_n/\widehat\sigma_n}; large positive values reject,
#' so the p-value is the upper standard-normal tail.  The local-alternative
#' oracle variance is not used.
#'
#' Each diagonal fit follows the published spatial-rank fixed-point update.
#' `scale_identification = "paper_trace"` applies the article's literal
#' \eqn{\operatorname{tr}(D)=p} normalization.  Because separately fitted
#' trace-normalized matrices acquire different scalar representatives after
#' coordinatewise rescaling, their weighted pool is not exactly
#' coordinatewise-scale equivariant in finite samples.  The default
#' `"geometric"` instead imposes unit geometric mean.  All leave-out fits
#' then acquire the same harmless scalar under a common coordinatewise
#' rescaling, making the complete statistic exactly coordinatewise-scale
#' invariant while leaving the published fixed-point equation unchanged.
#' Both choices and their distinction are recorded in `diagnostics`.
#'
#' At least six observations per group are necessary: a leave-four-out fit
#' must retain at least two observations.  More observations may be needed
#' for a particular data set if a retained subset has zero marginal variance
#' or zero spatial-rank energy.  Such degeneracy, nonconvergence, and a
#' non-positive feasible variance are errors.  No ridge, absolute value,
#' numerical floor, or random perturbation is applied.  Pairwise ties follow
#' \eqn{U(0)=0} and are counted.
#'
#' The formal asymptotic calibration in the paper assumes a common scatter
#' matrix and its stated high-dimensional trace conditions.  The article's
#' numerical study considered unequal scatter, but it did not establish the
#' null limit in that setting; this function therefore does not claim an
#' unequal-scatter guarantee and does not switch to a simulation or bootstrap
#' calibration.
#'
#' @param x,y Numeric matrices or data frames with observations in rows and
#'   the same variables in columns.  Each sample must have at least six rows.
#' @param tol Finite positive convergence tolerance for every full-sample,
#'   leave-two-out, and leave-four-out diagonal spatial-rank scale fit.
#' @param max_iter Positive integer maximum number of fixed-point updates for
#'   every scale fit.  Nonconvergence is an error.
#' @param scale_identification Scale representative used before pooling
#'   group fits.  The default `"geometric"` gives exact coordinatewise-scale
#'   invariance.  `"paper_trace"` reproduces the article's literal
#'   trace-normalized recursion for formula auditing.
#'
#' @return An object of class `c("hd_location_test", "htest")`.
#'   `components` contains the raw statistic, all three feasible trace and
#'   variance terms, exact ordered-sum denominators, block/permutation
#'   contributions, and every full/leave-two/leave-four scale fit.
#'   `diagnostics` records convergence, zero directions, preprocessing,
#'   identification, the corrected author-TeX factor, applicability, and the
#'   no-repair policy.
#'
#' @references
#' Feng, L., Zhang, X., and Liu, B. (2020). A high-dimensional spatial rank
#' test for the two-sample location problem. *Computational Statistics & Data
#' Analysis*, **144**, 106889. \doi{10.1016/j.csda.2019.106889}.
#'
#' @examples
#' x <- matrix(c(
#'   0.7, -1.1, 0.2, -0.4, 0.8, 1.2, 1.1, 0.3, -0.9,
#'   -1.2, -0.5, 0.7, 0.2, 1.4, -0.4, 1.5, -0.2, 0.5
#' ), ncol = 3, byrow = TRUE)
#' y <- matrix(c(
#'   -0.3, 0.7, -0.1, 1.2, -0.6, 0.8, -1.1, 0.2, 1.3,
#'   0.5, 1.1, -0.7, 0.9, -1.3, 0.4, -0.6, -0.4, -1.2
#' ), ncol = 3, byrow = TRUE)
#' feng_zhang_liu_spatial_rank_test(x, y, tol = 1e-6)
#'
#' @export
feng_zhang_liu_spatial_rank_test <- function(
    x, y, tol = 1e-7, max_iter = 500L,
    scale_identification = c("geometric", "paper_trace")) {
  call <- match.call()
  x.name <- deparse1(substitute(x))
  y.name <- deparse1(substitute(y))
  x <- .as_data_matrix(x, "x", min_rows = 6L)
  y <- .as_data_matrix(y, "y", min_rows = 6L)
  .check_two_sample_variables(x, y)
  controls <- .validate_iteration_controls(tol, max_iter, zero_tol = 0)
  scale_identification <- match.arg(scale_identification)

  fit <- cpp_feng_zhang_liu_spatial_rank(
    x, y, controls$tol, controls$max_iter, scale_identification
  )
  n1 <- nrow(x)
  n2 <- nrow(y)
  p <- ncol(x)
  variable.names <- .hotelling_variable_names(x)
  row.names1 <- rownames(x)
  row.names2 <- rownames(y)
  if (is.null(row.names1)) row.names1 <- paste0("x", seq_len(n1))
  if (is.null(row.names2)) row.names2 <- paste0("y", seq_len(n2))

  name.vector <- function(value) {
    stats::setNames(as.numeric(value), variable.names)
  }
  name.fit.matrix <- function(value, labels) {
    value <- as.matrix(value)
    dimnames(value) <- list(labels, variable.names)
    value
  }
  pair.labels <- function(indices, rows) {
    indices <- as.matrix(indices)
    paste(rows[indices[, 1L]], rows[indices[, 2L]], sep = "|")
  }
  quad.labels <- function(indices, rows) {
    indices <- as.matrix(indices)
    apply(indices, 1L, function(index) paste(rows[index], collapse = "|"))
  }
  pair.labels1 <- pair.labels(fit$pair_indices1, row.names1)
  pair.labels2 <- pair.labels(fit$pair_indices2, row.names2)
  quad.labels1 <- quad.labels(fit$quad_indices1, row.names1)
  quad.labels2 <- quad.labels(fit$quad_indices2, row.names2)
  rownames(fit$pair_indices1) <- pair.labels1
  rownames(fit$pair_indices2) <- pair.labels2
  rownames(fit$quad_indices1) <- quad.labels1
  rownames(fit$quad_indices2) <- quad.labels2
  colnames(fit$pair_indices1) <- colnames(fit$pair_indices2) <-
    c("first", "second")
  colnames(fit$quad_indices1) <- colnames(fit$quad_indices2) <-
    paste0("index", seq_len(4L))

  make.full.fit <- function(group) {
    suffix <- as.character(group)
    list(
      scale.diagonal.standardized = name.vector(
        fit[[paste0("full_D", suffix, "_standardised")]]
      ),
      log.scale.diagonal.standardized = name.vector(
        fit[[paste0("full_log_D", suffix, "_standardised")]]
      ),
      scale.diagonal.input.canonical = name.vector(
        fit[[paste0("full_D", suffix, "_input_canonical")]]
      ),
      log.scale.diagonal.input.canonical = name.vector(
        fit[[paste0("full_log_D", suffix, "_input_canonical")]]
      ),
      iterations = unname(fit[[paste0("full_iterations", suffix)]]),
      scale.equation.residual = unname(
        fit[[paste0("full_residual", suffix)]]
      ),
      zero.pair.directions = unname(
        fit[[paste0("full_zero_pair_directions", suffix)]]
      ),
      minimum.log.pair.radius = unname(
        fit[[paste0("full_minimum_log_pair_radius", suffix)]]
      )
    )
  }
  make.leave.fit <- function(order, group, labels) {
    suffix <- as.character(group)
    prefix <- paste0("leave", order, "_")
    list(
      indices = fit[[paste0(
        if (order == 2L) "pair_indices" else "quad_indices", suffix
      )]],
      scale.diagonal.standardized = name.fit.matrix(
        fit[[paste0(prefix, "D", suffix, "_standardised")]], labels
      ),
      log.scale.diagonal.standardized = name.fit.matrix(
        fit[[paste0(prefix, "log_D", suffix, "_standardised")]], labels
      ),
      scale.diagonal.standardized.canonical = name.fit.matrix(
        fit[[paste0(prefix, "D", suffix, "_standardised_canonical")]],
        labels
      ),
      scale.diagonal.input.canonical = name.fit.matrix(
        fit[[paste0(prefix, "D", suffix, "_input_canonical")]], labels
      ),
      log.scale.diagonal.input.canonical = name.fit.matrix(
        fit[[paste0(prefix, "log_D", suffix, "_input_canonical")]],
        labels
      ),
      iterations = stats::setNames(
        as.numeric(fit[[paste0(prefix, "iterations", suffix)]]), labels
      ),
      scale.equation.residual = stats::setNames(
        as.numeric(fit[[paste0(prefix, "residuals", suffix)]]), labels
      ),
      zero.pair.directions = stats::setNames(
        as.numeric(fit[[paste0(prefix, "zero_pair_directions", suffix)]]),
        labels
      ),
      minimum.log.pair.radius = stats::setNames(
        as.numeric(fit[[paste0(prefix,
          "minimum_log_pair_radii", suffix)]]), labels
      )
    )
  }

  full.fit <- list(group1 = make.full.fit(1L), group2 = make.full.fit(2L))
  leave2.fit <- list(
    group1 = make.leave.fit(2L, 1L, pair.labels1),
    group2 = make.leave.fit(2L, 2L, pair.labels2)
  )
  leave4.fit <- list(
    group1 = make.leave.fit(4L, 1L, quad.labels1),
    group2 = make.leave.fit(4L, 2L, quad.labels2)
  )

  mean.x <- name.vector(fit$sample_mean1)
  mean.y <- name.vector(fit$sample_mean2)
  difference <- stats::setNames(mean.x - mean.y, variable.names)
  z <- unname(fit$z)
  p.value <- stats::pnorm(z, lower.tail = FALSE)
  all.iterations <- c(
    full.fit$group1$iterations, full.fit$group2$iterations,
    leave2.fit$group1$iterations, leave2.fit$group2$iterations,
    leave4.fit$group1$iterations, leave4.fit$group2$iterations
  )
  all.residuals <- c(
    full.fit$group1$scale.equation.residual,
    full.fit$group2$scale.equation.residual,
    leave2.fit$group1$scale.equation.residual,
    leave2.fit$group2$scale.equation.residual,
    leave4.fit$group1$scale.equation.residual,
    leave4.fit$group2$scale.equation.residual
  )

  dimnames(fit$main_pair_block_contributions) <-
    list(pair.labels1, pair.labels2)
  dimnames(fit$trace3_pair_block_contributions) <-
    list(pair.labels1, pair.labels2)
  rownames(fit$trace1_permutation_contributions) <- quad.labels1
  rownames(fit$trace2_permutation_contributions) <- quad.labels2
  colnames(fit$trace1_permutation_contributions) <-
    colnames(fit$trace2_permutation_contributions) <-
      paste0("permutation", seq_len(24L))

  .new_hd_location_test(
    statistic = c(Z = z),
    p.value = p.value,
    method = paste(
      "Feng-Zhang-Liu high-dimensional two-sample spatial-rank test",
      "(leave-out statistic; feasible asymptotic normal calibration)"
    ),
    data.name = paste(x.name, "and", y.name),
    alternative = "two.sided",
    raw.statistic = c(T.n = unname(fit$T_n)),
    estimate = difference,
    null.value = stats::setNames(numeric(p), variable.names),
    null.distribution = list(
      family = "normal",
      parameters = c(mean = 0, sd = 1),
      exact = FALSE,
      tail = "upper",
      assumptions = paste(
        "Feng-Zhang-Liu high-dimensional trace conditions; independent",
        "samples with common scatter; converged diagonal rank-scale fits"
      )
    ),
    variance = c(estimated = unname(fit$sigma2_hat)),
    components = list(
      mean.x = mean.x,
      mean.y = mean.y,
      difference = difference,
      T.n = unname(fit$T_n),
      sigma2.hat = unname(fit$sigma2_hat),
      sigma.hat = unname(fit$sigma_hat),
      trace.R1.squared.hat = unname(fit$trace_R1_squared_hat),
      trace.R2.squared.hat = unname(fit$trace_R2_squared_hat),
      trace.R1.R2.hat = unname(fit$trace_R1_R2_hat),
      variance.term1 = unname(fit$variance_term1),
      variance.term2 = unname(fit$variance_term2),
      variance.term3 = unname(fit$variance_term3),
      main.ordered.sum = unname(fit$main_ordered_sum),
      main.ordered.denominator = unname(fit$main_ordered_denominator),
      main.pair.block.contributions = fit$main_pair_block_contributions,
      trace1.ordered.sum = unname(fit$trace1_ordered_sum),
      trace2.ordered.sum = unname(fit$trace2_ordered_sum),
      trace3.ordered.sum = unname(fit$trace3_ordered_sum),
      trace1.ordered.denominator = unname(
        fit$trace1_ordered_denominator
      ),
      trace2.ordered.denominator = unname(
        fit$trace2_ordered_denominator
      ),
      trace3.ordered.denominator = unname(
        fit$trace3_ordered_denominator
      ),
      trace1.permutation.contributions =
        fit$trace1_permutation_contributions,
      trace2.permutation.contributions =
        fit$trace2_permutation_contributions,
      trace3.pair.block.contributions =
        fit$trace3_pair_block_contributions,
      full.fit = full.fit,
      leave.two.out.fit = leave2.fit,
      leave.four.out.fit = leave4.fit,
      n1 = unname(fit$n1),
      n2 = unname(fit$n2),
      p = unname(fit$p),
      pair.count1 = unname(fit$pair_count1),
      pair.count2 = unname(fit$pair_count2),
      quadruple.count1 = unname(fit$quad_count1),
      quadruple.count2 = unname(fit$quad_count2)
    ),
    diagnostics = list(
      all.fits.converged = TRUE,
      tolerance = unname(fit$tolerance),
      max.iterations.allowed = unname(fit$max_iterations),
      iteration.summary = c(
        minimum = min(all.iterations),
        median = stats::median(all.iterations),
        mean = mean(all.iterations),
        maximum = max(all.iterations)
      ),
      maximum.scale.equation.residual = max(all.residuals),
      full.fit = full.fit,
      leave.two.out.fit = leave2.fit,
      leave.four.out.fit = leave4.fit,
      zero.directions = list(
        full.group1 = full.fit$group1$zero.pair.directions,
        full.group2 = full.fit$group2$zero.pair.directions,
        leave.two.group1 = leave2.fit$group1$zero.pair.directions,
        leave.two.group2 = leave2.fit$group2$zero.pair.directions,
        leave.four.group1 = leave4.fit$group1$zero.pair.directions,
        leave.four.group2 = leave4.fit$group2$zero.pair.directions,
        main.evaluations = unname(fit$main_zero_direction_evaluations),
        main.ordered.sign.uses = unname(fit$main_ordered_zero_sign_uses),
        trace1.evaluations = unname(
          fit$trace1_zero_direction_evaluations
        ),
        trace2.evaluations = unname(
          fit$trace2_zero_direction_evaluations
        ),
        trace3.evaluations = unname(
          fit$trace3_zero_direction_evaluations
        )
      ),
      sign.at.zero = "U(0) = 0 for pairwise and cross-sample directions",
      scale.identification = scale_identification,
      scale.identification.note = if (scale_identification == "geometric") {
        paste(
          "Unit-geometric-mean representatives are used for every fit;",
          "this coherent identification restores exact coordinatewise-scale",
          "invariance of the weighted pooled leave-out scale."
        )
      } else {
        paste(
          "Literal paper normalization tr(D)=p is used independently for",
          "each fit; weighted pooling can therefore violate the paper's",
          "claimed exact coordinatewise-scale invariance in finite samples."
        )
      },
      paper.trace.identification.available = TRUE,
      feasible.variance = paste(
        "two leave-four-out within-group ordered-quadruple traces plus",
        "one pooled leave-two-out cross trace"
      ),
      second.variance.term.p.squared.restored = TRUE,
      author.tex.erratum = paste(
        "The feasible variance's second term requires p^2 in its",
        "denominator; the author TeX omits it although the symmetric oracle",
        "formula and consistency require it."
      ),
      calibration = paste(
        "feasible null sigma.hat and upper normal tail; local-alternative",
        "oracle variance is not used"
      ),
      applicability = list(
        formal.common.scatter.required = TRUE,
        unequal.scatter.theory.available.in.paper = FALSE,
        note = paste(
          "The published null theorem assumes common scatter. Unequal",
          "scatter was examined numerically but has no stated null theorem;",
          "no simulation or bootstrap calibration is implemented here."
        )
      ),
      preprocessing = list(
        anchor = name.vector(fit$preprocessing_anchor),
        log.coordinate.scale = name.vector(
          fit$preprocessing_log_coordinate_scale
        ),
        scaled.coordinate.range = name.vector(
          fit$preprocessing_scaled_coordinate_range
        ),
        global.operand.scale = unname(
          fit$preprocessing_global_operand_scale
        ),
        standardized.x = fit$preprocessed_x,
        standardized.y = fit$preprocessed_y
      ),
      regularization = "none",
      variance.repair = "none"
    ),
    n = c(group1 = n1, group2 = n2),
    p = p,
    call = call
  )
}
