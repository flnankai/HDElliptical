#' Feng--Zou--Wang--Zhu Composite T-squared two-sample test
#'
#' Implements the two-sample Composite \eqn{T^2} test of Feng, Zou, Wang,
#' and Zhu (2017).  Observations are rows and variables are columns.  Despite
#' the wording at `ch2_location.tex:778--784` in the current book draft, this
#' is **not** a Behrens--Fisher test: the paper assumes that the two groups
#' have the same unknown covariance matrix.  Nor does the method combine
#' existing test p-values or estimate correlations among component tests.
#' Instead, it builds a block-diagonal approximation to the common pooled
#' covariance and combines the resulting block Hotelling quadratic forms.
#'
#' For a fixed leave-out covariance, the paper's practical block rule starts
#' with the pair having the largest absolute sample correlation and repeatedly
#' adds the variable with the largest sum of absolute correlations to the
#' current block.  Variables already assigned to a block are removed and the
#' procedure is repeated.  `selection = "paper_greedy"` implements that rule,
#' with lexicographic tie-breaking.  The final block contains all remaining
#' variables when fewer than `block_size` remain.  The paper's notation and
#' simulations use equal-size blocks; this deterministic remainder convention
#' extends the displayed rule to dimensions not divisible by the block size.
#'
#' Let \eqn{\widehat{\Sigma}_{O_K,i_1,i_2,j_1,j_2}^{-1}} be the
#' block-diagonal inverse obtained after deleting observations \eqn{i_1,i_2}
#' from group 1 and \eqn{j_1,j_2} from group 2, recomputing the unbiased pooled
#' covariance and its correlation matrix, and rebuilding the blocks.  The
#' published statistic is
#' \deqn{Q_n={1\over n_1n_2(n_1-1)(n_2-1)}
#' \sum_{i_1\ne i_2}\sum_{j_1\ne j_2}
#' (X_{1i_1}-X_{2j_1})^T
#' \widehat{\Sigma}_{O_K,i_1,i_2,j_1,j_2}^{-1}
#' (X_{1i_2}-X_{2j_2}).}
#'
#' The upper-tail normal calibration is
#' \deqn{Z={Q_n\over
#' {2(n_1^{-1}+n_2^{-1})^2
#' \widehat{\mathrm{tr}(\Lambda_K^2)}}^{1/2}}.}
#' Following Theorem 2 and the paragraph immediately after it, the feasible
#' trace estimate is the paper's one-sample leave-four-out estimator computed
#' from **group 1 only**:
#' \deqn{\widehat{\mathrm{tr}(\Lambda_K^2)}=
#' {1\over 2P_{n_1}^4}\sum^*
#' (X_{i_1}-X_{i_2})^T\widehat{\Sigma}_{O_K,\setminus4}^{-1}
#' (X_{i_3}-X_{i_4})
#' (X_{i_1}-X_{i_4})^T\widehat{\Sigma}_{O_K,\setminus4}^{-1}
#' (X_{i_3}-X_{i_2}).}
#' Thus \eqn{Q_n} and its population scaling are symmetric in the two samples,
#' but the paper's feasible finite-sample standardisation is label-asymmetric.
#' Exchanging the samples can change the reported \eqn{Z} and p-value.
#'
#' The first group must satisfy `nrow(x) >= block_size + 5`: after leaving out
#' four rows, a block of size `block_size` then has enough residual degrees of
#' freedom to be nonsingular.  The second group must contain at least three
#' rows for the pooled 2+2 leave-out covariance.  These count conditions are
#' necessary, not sufficient.  Every marginal variance used for selection and
#' every selected full or leave-out covariance block must be finite and
#' strictly positive definite.  Failure produces an error.  No ridge,
#' generalized inverse, absolute-value repair, or variance floor is applied.
#'
#' The C++ kernel uses sufficient statistics and sums unordered pairs and
#' quadruples as exact symmetry reductions of the displayed ordered formulas.
#' Both samples are internally translated by a common anchor and divided by a
#' common positive scale in each column.  This is algebraically neutral and
#' protects the scale-invariant method under very large or small units.
#'
#' @param x,y Numeric matrices or data frames with observations in rows and
#'   the same variables in columns.  The first sample supplies the feasible
#'   trace calibration.
#' @param block_size Positive integer block size \eqn{K}, no larger than the
#'   number of variables.  The paper recommends small fixed blocks and uses
#'   `2` in its main implementation.
#' @param selection Block construction rule.  The only supported value is
#'   `"paper_greedy"`, the practical algorithm in Remark 1 of the paper.
#'
#' @return An object of class `c("hd_location_test", "htest")`.  The
#'   statistic is the standardised \eqn{Z} and the p-value is its upper
#'   standard-normal tail probability.  `components` retains \eqn{Q_n}, the
#'   group-1 leave-four-out trace estimate, every calibration coefficient and
#'   ordered denominator, combination counts, and a full-sample diagnostic
#'   partition.  `diagnostics` records the common-covariance assumption,
#'   dynamic block selection, first-sample calibration, conditioning summaries,
#'   and the absence of numerical repair.
#'
#' @references
#' Feng, L., Zou, C., Wang, Z., and Zhu, L. (2017). Composite T-squared test
#' for high-dimensional data. *Statistica Sinica*, **27**, 1419--1436.
#' \doi{10.5705/ss.202015.0199}.
#'
#' @examples
#' set.seed(2717)
#' x <- matrix(rnorm(21), 7, 3)
#' y <- matrix(rnorm(15, 0.2), 5, 3)
#' composite_t2_two_sample_test(x, y)
#'
#' @export
composite_t2_two_sample_test <- function(
    x, y, block_size = 2L, selection = "paper_greedy") {
  call <- match.call()
  x.name <- deparse1(substitute(x))
  y.name <- deparse1(substitute(y))
  x <- .as_data_matrix(x, "x", min_rows = 1L)
  y <- .as_data_matrix(y, "y", min_rows = 3L)
  .check_two_sample_variables(x, y)

  if (!is.numeric(block_size) || length(block_size) != 1L ||
      is.na(block_size) || !is.finite(block_size) || block_size < 1 ||
      block_size != floor(block_size) ||
      block_size > .Machine$integer.max) {
    stop("`block_size` must be a positive integer.", call. = FALSE)
  }
  block_size <- as.integer(block_size)
  if (block_size > ncol(x)) {
    stop("`block_size` must not exceed the number of variables.",
         call. = FALSE)
  }
  selection <- match.arg(selection, "paper_greedy")
  if (nrow(x) < block_size + 5L) {
    stop(
      "Composite T-squared requires `nrow(x) >= block_size + 5` " %+%
        "for its group-1 leave-four-out trace calibration.",
      call. = FALSE
    )
  }

  n1 <- nrow(x)
  n2 <- nrow(y)
  p <- ncol(x)
  fit <- cpp_composite_t2_two_sample(
    x, y, block_size = block_size, selection = selection
  )
  z <- as.numeric(fit$z)
  p.value <- stats::pnorm(z, lower.tail = FALSE)
  variable.names <- .hotelling_variable_names(x)

  named.vector <- function(value) {
    stats::setNames(as.numeric(value), variable.names)
  }
  mean.x <- named.vector(fit$mean_x)
  mean.y <- named.vector(fit$mean_y)
  difference <- named.vector(fit$difference)
  column.scale <- named.vector(fit$column_scale)
  full.variance.scaled <- named.vector(fit$full_pooled_variance_scaled)

  full.blocks <- lapply(fit$full_blocks, function(index) {
    index <- as.integer(index)
    stats::setNames(index, variable.names[index])
  })
  names(full.blocks) <- paste0("block", seq_along(full.blocks))
  full.block.quadratic <- stats::setNames(
    as.numeric(fit$full_block_quadratic), names(full.blocks)
  )
  full.partition.signature <- paste(
    vapply(full.blocks, function(index) {
      paste(names(index), collapse = ",")
    }, character(1L)),
    collapse = " | "
  )

  .new_hd_location_test(
    statistic = c(Z = z),
    p.value = p.value,
    method = paste(
      "Feng-Zou-Wang-Zhu Composite T-squared two-sample test",
      "(common covariance; asymptotic normal calibration)"
    ),
    data.name = paste(x.name, "and", y.name),
    alternative = "two.sided",
    raw.statistic = c(Q.n = as.numeric(fit$Q_n)),
    estimate = difference,
    null.value = stats::setNames(numeric(p), variable.names),
    null.distribution = list(
      family = "normal",
      parameters = c(mean = 0, sd = 1),
      exact = FALSE,
      tail = "upper",
      assumptions = paste(
        "Feng-Zou-Wang-Zhu high-dimensional asymptotics; independent",
        "samples with a common covariance matrix; the paper's moment and",
        "dependence conditions; finite positive-definite selected leave-out",
        "blocks; strictly positive group-1 trace estimate"
      )
    ),
    variance = c(denominator = as.numeric(fit$variance)),
    components = list(
      mean.x = mean.x,
      mean.y = mean.y,
      difference = difference,
      Q.n = as.numeric(fit$Q_n),
      trace.Lambda.K.squared = as.numeric(fit$trace_lambda_squared),
      variance.coefficient = as.numeric(fit$variance_coefficient),
      denominator.variance = as.numeric(fit$variance),
      standard.error = as.numeric(fit$standard_error),
      n1 = as.numeric(fit$n1),
      n2 = as.numeric(fit$n2),
      p = as.numeric(fit$p),
      block.size = as.integer(fit$block_size),
      selection = as.character(fit$selection),
      Q.ordered.numerator = as.numeric(fit$q_ordered_numerator),
      Q.ordered.denominator = as.numeric(fit$q_ordered_denominator),
      Q.unordered.pair.combinations =
        as.numeric(fit$q_pair_combinations),
      trace.ordered.numerator =
        as.numeric(fit$trace_ordered_numerator),
      trace.ordered.denominator =
        as.numeric(fit$trace_ordered_denominator),
      trace.unordered.quadruples = as.numeric(fit$trace_quadruples),
      P4.n1 = as.numeric(fit$P4_n1),
      full.sample.blocks = full.blocks,
      full.sample.block.quadratic = full.block.quadratic,
      full.sample.partition.signature = full.partition.signature,
      full.pooled.variance.scaled = full.variance.scaled
    ),
    diagnostics = list(
      covariance.model = "common",
      behrens.fisher = FALSE,
      component.test.combination = FALSE,
      component.correlation.calibration = "not applicable",
      calibration = "upper-tail standard normal",
      trace.calibration.sample = "group1 only",
      finite.sample.group.label.symmetric = FALSE,
      Q.leaveout = paste(
        "delete two observations from each group; recompute unbiased",
        "pooled covariance, correlations, greedy blocks, and block inverses"
      ),
      trace.leaveout = paste(
        "delete four observations from group1; recompute its unbiased",
        "covariance, correlations, greedy blocks, and block inverses;",
        "denominator 2 * P_n1^4"
      ),
      block.selection = paste(
        "Remark 1 paper greedy rule on absolute correlations;",
        "lexicographic tie-breaking; final remainder block retained"
      ),
      dynamic.block.selection = TRUE,
      q.unique.partitions = as.numeric(fit$q_unique_partitions),
      trace.unique.partitions = as.numeric(fit$trace_unique_partitions),
      minimum.reciprocal.condition = c(
        Q.leave2.plus2 =
          as.numeric(fit$q_minimum_reciprocal_condition),
        trace.group1.leave4 =
          as.numeric(fit$trace_minimum_reciprocal_condition),
        full.pooled.diagnostic =
          as.numeric(fit$full_minimum_reciprocal_condition)
      ),
      minimum.cholesky.diagonal = c(
        Q.leave2.plus2 = as.numeric(fit$q_minimum_cholesky_diagonal),
        trace.group1.leave4 =
          as.numeric(fit$trace_minimum_cholesky_diagonal),
        full.pooled.diagnostic =
          as.numeric(fit$full_minimum_cholesky_diagonal)
      ),
      minimum.marginal.variance.scaled = c(
        Q.leave2.plus2 =
          as.numeric(fit$q_minimum_marginal_variance_scaled),
        trace.group1.leave4 =
          as.numeric(fit$trace_minimum_marginal_variance_scaled),
        full.pooled.diagnostic =
          as.numeric(fit$full_minimum_marginal_variance_scaled)
      ),
      minimum.sample.sizes = c(
        group1 = block_size + 5L,
        group2 = 3L
      ),
      internal.scaling = paste(
        "common per-column anchor and positive scale across both groups;",
        "long-double sufficient statistics"
      ),
      internal.scale.factors = column.scale,
      covariance.inverse = "strict Cholesky block solves",
      regularization = "none",
      generalized.inverse = "none",
      variance.repair = "none"
    ),
    n = c(group1 = n1, group2 = n2),
    p = p,
    call = call
  )
}
