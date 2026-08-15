// Feng--Zou--Wang--Zhu (2017) Composite T-squared two-sample test.
//
// The published two-sample method assumes a common covariance matrix.  For
// every ordered 2+2 leave-out tuple it recomputes the pooled covariance,
// rebuilds the correlation-driven blocks, and applies the corresponding
// block-diagonal inverse.  Its feasible trace calibration uses the first
// sample only and recomputes the blocks after leaving out four observations.
// The loops below use unordered pairs/quadruples only as exact symmetry
// reductions of those ordered formulas.

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp17)]]

#include <RcppArmadillo.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <numeric>
#include <set>
#include <sstream>
#include <string>
#include <vector>

namespace {

struct CompositeStandardizedData {
  arma::mat x;
  arma::mat y;
  std::vector<long double> origin;
  std::vector<long double> scale;
};

struct CompositeSufficientData {
  const arma::mat* values;
  arma::uword n;
  arma::uword p;
  std::vector<long double> sum;
  std::vector<long double> cross;
};

struct CompositeBlockFactor {
  std::vector<arma::uword> indices;
  arma::mat lower;
  double reciprocal_condition;
};

struct CompositePartition {
  std::vector<CompositeBlockFactor> blocks;
  std::string signature;
  double minimum_reciprocal_condition;
  double minimum_cholesky_diagonal;
  double minimum_marginal_variance;
};

struct CompositeAccumulator {
  long double value = 0.0L;
  long double correction = 0.0L;

  void add(long double increment) {
    const long double adjusted = increment - correction;
    const long double updated = value + adjusted;
    correction = (updated - value) - adjusted;
    value = updated;
  }
};

struct CompositeLoopDiagnostics {
  double minimum_reciprocal_condition =
    std::numeric_limits<double>::infinity();
  double minimum_cholesky_diagonal =
    std::numeric_limits<double>::infinity();
  double minimum_marginal_variance =
    std::numeric_limits<double>::infinity();
  std::set<std::string> partitions;
};

inline std::size_t composite_cross_index(arma::uword j,
                                         arma::uword k,
                                         arma::uword p) {
  return static_cast<std::size_t>(j) * static_cast<std::size_t>(p) +
    static_cast<std::size_t>(k);
}

void composite_require_finite_matrix(const arma::mat& x,
                                     const char* name) {
  if (!x.is_finite()) {
    Rcpp::stop("`%s` must contain only finite values.", name);
  }
}

double composite_checked_double(long double value,
                                const char* quantity) {
  const double result = static_cast<double>(value);
  if (!std::isfinite(result)) {
    Rcpp::stop(
      "Composite T-squared requires a finite, double-representable %s.",
      quantity
    );
  }
  return result;
}

double composite_checked_positive_double(long double value,
                                         const char* quantity) {
  const double result = static_cast<double>(value);
  if (!std::isfinite(result) || result <= 0.0) {
    Rcpp::stop(
      "Composite T-squared requires a finite, strictly positive and "
      "double-representable %s; no ridge, absolute-value repair, or "
      "numerical floor is applied.",
      quantity
    );
  }
  return result;
}

double composite_report_long_double(long double value) {
  const long double maximum = static_cast<long double>(
    std::numeric_limits<double>::max()
  );
  if (value > maximum) {
    return std::numeric_limits<double>::infinity();
  }
  if (value < -maximum) {
    return -std::numeric_limits<double>::infinity();
  }
  return static_cast<double>(value);
}

long double composite_choose2(arma::uword n) {
  const long double value = static_cast<long double>(n);
  return value * (value - 1.0L) / 2.0L;
}

long double composite_choose4(arma::uword n) {
  const long double value = static_cast<long double>(n);
  return value * (value - 1.0L) * (value - 2.0L) *
    (value - 3.0L) / 24.0L;
}

CompositeStandardizedData composite_standardize_columns(
    const arma::mat& x,
    const arma::mat& y) {
  CompositeStandardizedData out;
  const arma::uword n1 = x.n_rows;
  const arma::uword n2 = y.n_rows;
  const arma::uword p = x.n_cols;
  out.x.set_size(n1, p);
  out.y.set_size(n2, p);
  out.origin.resize(p);
  out.scale.resize(p);

  for (arma::uword j = 0; j < p; ++j) {
    const long double anchor = static_cast<long double>(x(0, j));
    long double maximum_difference = 0.0L;
    bool direct_finite = true;
    for (arma::uword i = 0; i < n1; ++i) {
      const long double difference =
        static_cast<long double>(x(i, j)) - anchor;
      direct_finite = direct_finite && std::isfinite(difference);
      maximum_difference = std::max(
        maximum_difference, std::abs(difference)
      );
    }
    for (arma::uword i = 0; i < n2; ++i) {
      const long double difference =
        static_cast<long double>(y(i, j)) - anchor;
      direct_finite = direct_finite && std::isfinite(difference);
      maximum_difference = std::max(
        maximum_difference, std::abs(difference)
      );
    }

    out.origin[j] = anchor;
    if (direct_finite && std::isfinite(maximum_difference)) {
      const long double scale = maximum_difference > 0.0L ?
        maximum_difference : 1.0L;
      out.scale[j] = scale;
      for (arma::uword i = 0; i < n1; ++i) {
        out.x(i, j) = static_cast<double>(
          (static_cast<long double>(x(i, j)) - anchor) / scale
        );
      }
      for (arma::uword i = 0; i < n2; ++i) {
        out.y(i, j) = static_cast<double>(
          (static_cast<long double>(y(i, j)) - anchor) / scale
        );
      }
      continue;
    }

    // On platforms whose long double has the double exponent range, direct
    // subtraction of opposite finite extremes can overflow.  Scaling the
    // operands first preserves the scale-invariant statistic.
    double operand_scale = std::abs(x(0, j));
    for (arma::uword i = 0; i < n1; ++i) {
      operand_scale = std::max(operand_scale, std::abs(x(i, j)));
    }
    for (arma::uword i = 0; i < n2; ++i) {
      operand_scale = std::max(operand_scale, std::abs(y(i, j)));
    }
    if (operand_scale == 0.0) {
      operand_scale = 1.0;
    }
    const long double operand_scale_ld =
      static_cast<long double>(operand_scale);
    const long double anchor_scaled = anchor / operand_scale_ld;
    maximum_difference = 0.0L;
    for (arma::uword i = 0; i < n1; ++i) {
      const long double difference =
        static_cast<long double>(x(i, j)) / operand_scale_ld -
        anchor_scaled;
      out.x(i, j) = static_cast<double>(difference);
      maximum_difference = std::max(
        maximum_difference, std::abs(difference)
      );
    }
    for (arma::uword i = 0; i < n2; ++i) {
      const long double difference =
        static_cast<long double>(y(i, j)) / operand_scale_ld -
        anchor_scaled;
      out.y(i, j) = static_cast<double>(difference);
      maximum_difference = std::max(
        maximum_difference, std::abs(difference)
      );
    }
    const long double second_scale = maximum_difference > 0.0L ?
      maximum_difference : 1.0L;
    for (arma::uword i = 0; i < n1; ++i) {
      out.x(i, j) = static_cast<double>(
        static_cast<long double>(out.x(i, j)) / second_scale
      );
    }
    for (arma::uword i = 0; i < n2; ++i) {
      out.y(i, j) = static_cast<double>(
        static_cast<long double>(out.y(i, j)) / second_scale
      );
    }
    out.scale[j] = operand_scale_ld * second_scale;
  }
  return out;
}

CompositeSufficientData composite_sufficient_data(const arma::mat& values) {
  CompositeSufficientData out;
  out.values = &values;
  out.n = values.n_rows;
  out.p = values.n_cols;
  out.sum.assign(out.p, 0.0L);
  out.cross.assign(
    static_cast<std::size_t>(out.p) * static_cast<std::size_t>(out.p),
    0.0L
  );

  for (arma::uword i = 0; i < out.n; ++i) {
    for (arma::uword j = 0; j < out.p; ++j) {
      const long double zj = static_cast<long double>(values(i, j));
      out.sum[j] += zj;
      for (arma::uword k = 0; k <= j; ++k) {
        out.cross[composite_cross_index(j, k, out.p)] +=
          zj * static_cast<long double>(values(i, k));
      }
    }
  }
  for (arma::uword j = 0; j < out.p; ++j) {
    for (arma::uword k = 0; k < j; ++k) {
      out.cross[composite_cross_index(k, j, out.p)] =
        out.cross[composite_cross_index(j, k, out.p)];
    }
  }
  return out;
}

template <std::size_t M>
arma::mat composite_scatter_excluding(
    const CompositeSufficientData& data,
    const std::array<arma::uword, M>& omitted) {
  const arma::uword remaining = data.n - static_cast<arma::uword>(M);
  if (remaining < 1) {
    Rcpp::stop("Composite T-squared leave-out covariance has no rows left.");
  }
  const long double remaining_ld = static_cast<long double>(remaining);
  std::vector<long double> sums(data.p);
  for (arma::uword j = 0; j < data.p; ++j) {
    long double value = data.sum[j];
    for (std::size_t r = 0; r < M; ++r) {
      value -= static_cast<long double>((*(data.values))(omitted[r], j));
    }
    sums[j] = value;
  }

  arma::mat scatter(data.p, data.p, arma::fill::zeros);
  for (arma::uword j = 0; j < data.p; ++j) {
    for (arma::uword k = 0; k <= j; ++k) {
      long double cross =
        data.cross[composite_cross_index(j, k, data.p)];
      for (std::size_t r = 0; r < M; ++r) {
        cross -= static_cast<long double>(
          (*(data.values))(omitted[r], j)
        ) * static_cast<long double>(
          (*(data.values))(omitted[r], k)
        );
      }
      const long double centered = cross -
        sums[j] * sums[k] / remaining_ld;
      const double value = static_cast<double>(centered);
      if (!std::isfinite(value)) {
        Rcpp::stop(
          "Composite T-squared produced a non-finite leave-out scatter."
        );
      }
      scatter(j, k) = value;
      scatter(k, j) = value;
    }
  }
  return scatter;
}

arma::mat composite_scatter_all(const CompositeSufficientData& data) {
  const long double n = static_cast<long double>(data.n);
  arma::mat scatter(data.p, data.p, arma::fill::zeros);
  for (arma::uword j = 0; j < data.p; ++j) {
    for (arma::uword k = 0; k <= j; ++k) {
      const long double centered =
        data.cross[composite_cross_index(j, k, data.p)] -
        data.sum[j] * data.sum[k] / n;
      const double value = static_cast<double>(centered);
      if (!std::isfinite(value)) {
        Rcpp::stop(
          "Composite T-squared produced a non-finite full-sample scatter."
        );
      }
      scatter(j, k) = value;
      scatter(k, j) = value;
    }
  }
  return scatter;
}

arma::mat composite_pooled_covariance_leave2(
    const CompositeSufficientData& group1,
    const CompositeSufficientData& group2,
    const std::array<arma::uword, 2>& omitted1,
    const std::array<arma::uword, 2>& omitted2) {
  const long double degrees = static_cast<long double>(
    group1.n + group2.n - 6
  );
  if (degrees <= 0.0L) {
    Rcpp::stop(
      "Composite T-squared pooled leave-2+2 covariance has no residual "
      "degrees of freedom."
    );
  }
  return (
    composite_scatter_excluding(group1, omitted1) +
    composite_scatter_excluding(group2, omitted2)
  ) / static_cast<double>(degrees);
}

arma::mat composite_group1_covariance_leave4(
    const CompositeSufficientData& group1,
    const std::array<arma::uword, 4>& omitted) {
  const long double degrees = static_cast<long double>(group1.n - 5);
  if (degrees <= 0.0L) {
    Rcpp::stop(
      "Composite T-squared group-1 leave-four covariance has no residual "
      "degrees of freedom."
    );
  }
  return composite_scatter_excluding(group1, omitted) /
    static_cast<double>(degrees);
}

arma::mat composite_full_pooled_covariance(
    const CompositeSufficientData& group1,
    const CompositeSufficientData& group2) {
  const long double degrees = static_cast<long double>(
    group1.n + group2.n - 2
  );
  return (composite_scatter_all(group1) +
          composite_scatter_all(group2)) /
    static_cast<double>(degrees);
}

std::string composite_partition_signature(
    const std::vector<std::vector<arma::uword>>& blocks) {
  std::vector<std::string> block_signatures;
  block_signatures.reserve(blocks.size());
  for (const auto& block : blocks) {
    std::ostringstream stream;
    for (std::size_t j = 0; j < block.size(); ++j) {
      if (j > 0) {
        stream << ',';
      }
      stream << block[j];
    }
    block_signatures.push_back(stream.str());
  }
  std::sort(block_signatures.begin(), block_signatures.end());
  std::ostringstream result;
  for (std::size_t b = 0; b < block_signatures.size(); ++b) {
    if (b > 0) {
      result << '|';
    }
    result << block_signatures[b];
  }
  return result.str();
}

std::vector<std::vector<arma::uword>> composite_greedy_blocks(
    const arma::mat& correlation,
    arma::uword block_size) {
  const arma::uword p = correlation.n_rows;
  std::vector<arma::uword> remaining(p);
  std::iota(remaining.begin(), remaining.end(), 0);
  std::vector<std::vector<arma::uword>> blocks;

  while (!remaining.empty()) {
    if (remaining.size() <= static_cast<std::size_t>(block_size)) {
      blocks.push_back(remaining);
      break;
    }

    std::vector<arma::uword> selected;
    if (block_size == 1) {
      selected.push_back(remaining.front());
    } else {
      arma::uword best_first = remaining[0];
      arma::uword best_second = remaining[1];
      double best_score = -1.0;
      for (std::size_t a = 0; a + 1 < remaining.size(); ++a) {
        for (std::size_t b = a + 1; b < remaining.size(); ++b) {
          const double score = std::abs(
            correlation(remaining[a], remaining[b])
          );
          if (score > best_score) {
            best_score = score;
            best_first = remaining[a];
            best_second = remaining[b];
          }
        }
      }
      selected.push_back(best_first);
      selected.push_back(best_second);

      while (selected.size() < static_cast<std::size_t>(block_size)) {
        arma::uword best_variable = 0;
        double best_addition = -1.0;
        bool initialized = false;
        for (arma::uword candidate : remaining) {
          if (std::find(selected.begin(), selected.end(), candidate) !=
              selected.end()) {
            continue;
          }
          double score = 0.0;
          for (arma::uword chosen : selected) {
            score += std::abs(correlation(candidate, chosen));
          }
          if (!initialized || score > best_addition) {
            initialized = true;
            best_addition = score;
            best_variable = candidate;
          }
        }
        selected.push_back(best_variable);
      }
    }

    std::sort(selected.begin(), selected.end());
    blocks.push_back(selected);
    std::vector<arma::uword> next;
    next.reserve(remaining.size() - selected.size());
    for (arma::uword candidate : remaining) {
      if (!std::binary_search(selected.begin(), selected.end(), candidate)) {
        next.push_back(candidate);
      }
    }
    remaining.swap(next);
  }
  return blocks;
}

CompositePartition composite_factor_partition(const arma::mat& covariance,
                                               arma::uword block_size,
                                               const char* stage) {
  if (covariance.n_rows < 1 || covariance.n_rows != covariance.n_cols) {
    Rcpp::stop("Composite T-squared received an invalid covariance matrix.");
  }
  const arma::uword p = covariance.n_rows;
  arma::mat correlation(p, p, arma::fill::eye);
  double minimum_variance = std::numeric_limits<double>::infinity();
  for (arma::uword j = 0; j < p; ++j) {
    const double variance = covariance(j, j);
    if (!std::isfinite(variance) || variance <= 0.0) {
      Rcpp::stop(
        "Composite T-squared requires every marginal variance in each %s "
        "covariance to be finite and strictly positive; no variance floor "
        "is applied.",
        stage
      );
    }
    minimum_variance = std::min(minimum_variance, variance);
  }
  for (arma::uword j = 0; j < p; ++j) {
    for (arma::uword k = 0; k < j; ++k) {
      const double denominator =
        std::sqrt(covariance(j, j)) * std::sqrt(covariance(k, k));
      const double value = covariance(j, k) / denominator;
      if (!std::isfinite(value)) {
        Rcpp::stop(
          "Composite T-squared produced a non-finite %s correlation.",
          stage
        );
      }
      correlation(j, k) = value;
      correlation(k, j) = value;
    }
  }

  const std::vector<std::vector<arma::uword>> block_indices =
    composite_greedy_blocks(correlation, block_size);
  CompositePartition result;
  result.signature = composite_partition_signature(block_indices);
  result.minimum_reciprocal_condition =
    std::numeric_limits<double>::infinity();
  result.minimum_cholesky_diagonal =
    std::numeric_limits<double>::infinity();
  result.minimum_marginal_variance = minimum_variance;

  for (const auto& indices : block_indices) {
    arma::uvec armadillo_indices(indices.size());
    for (std::size_t j = 0; j < indices.size(); ++j) {
      armadillo_indices(j) = indices[j];
    }
    const arma::mat block = covariance.submat(
      armadillo_indices, armadillo_indices
    );
    arma::mat lower;
    const bool success = arma::chol(lower, block, "lower");
    if (!success || !lower.is_finite()) {
      Rcpp::stop(
        "Composite T-squared requires every selected %s block to be "
        "strictly positive definite. Cholesky factorization failed; no "
        "ridge or generalized inverse is applied.",
        stage
      );
    }
    const double reciprocal_condition = arma::rcond(block);
    if (!std::isfinite(reciprocal_condition) ||
        reciprocal_condition <= 0.0) {
      Rcpp::stop(
        "Composite T-squared requires every selected %s block to be "
        "numerically invertible; no ridge or generalized inverse is "
        "applied.",
        stage
      );
    }
    const double minimum_cholesky = lower.diag().min();
    result.minimum_reciprocal_condition = std::min(
      result.minimum_reciprocal_condition, reciprocal_condition
    );
    result.minimum_cholesky_diagonal = std::min(
      result.minimum_cholesky_diagonal, minimum_cholesky
    );
    result.blocks.push_back({indices, lower, reciprocal_condition});
  }
  return result;
}

void composite_update_loop_diagnostics(
    CompositeLoopDiagnostics& diagnostics,
    const CompositePartition& partition) {
  diagnostics.minimum_reciprocal_condition = std::min(
    diagnostics.minimum_reciprocal_condition,
    partition.minimum_reciprocal_condition
  );
  diagnostics.minimum_cholesky_diagonal = std::min(
    diagnostics.minimum_cholesky_diagonal,
    partition.minimum_cholesky_diagonal
  );
  diagnostics.minimum_marginal_variance = std::min(
    diagnostics.minimum_marginal_variance,
    partition.minimum_marginal_variance
  );
  diagnostics.partitions.insert(partition.signature);
}

arma::mat composite_solve_lower_cholesky(
    const CompositeBlockFactor& block,
    const arma::mat& right_hand_side,
    const char* stage) {
  arma::mat intermediate;
  arma::mat solution;
  const bool first = arma::solve(
    intermediate,
    arma::trimatl(block.lower),
    right_hand_side,
    arma::solve_opts::fast
  );
  const bool second = first && arma::solve(
    solution,
    arma::trimatu(block.lower.t()),
    intermediate,
    arma::solve_opts::fast
  );
  if (!second || !solution.is_finite()) {
    Rcpp::stop(
      "Composite T-squared could not solve a selected %s block; no "
      "generalized inverse or numerical repair is applied.",
      stage
    );
  }
  return solution;
}

arma::rowvec composite_bilinear_columns(
    const CompositePartition& partition,
    const arma::mat& left,
    const arma::mat& right,
    const char* stage) {
  if (left.n_rows != right.n_rows || left.n_cols != right.n_cols) {
    Rcpp::stop("Composite T-squared received incompatible kernel vectors.");
  }
  arma::rowvec result(left.n_cols, arma::fill::zeros);
  for (const CompositeBlockFactor& block : partition.blocks) {
    arma::uvec indices(block.indices.size());
    for (std::size_t j = 0; j < block.indices.size(); ++j) {
      indices(j) = block.indices[j];
    }
    const arma::mat left_block = left.rows(indices);
    const arma::mat right_block = right.rows(indices);
    const arma::mat solution = composite_solve_lower_cholesky(
      block, left_block, stage
    );
    result += arma::sum(right_block % solution, 0);
  }
  if (!result.is_finite()) {
    Rcpp::stop("Composite T-squared produced a non-finite %s kernel.", stage);
  }
  return result;
}

arma::mat composite_local_gram(const CompositePartition& partition,
                               const arma::mat& rows,
                               const char* stage) {
  if (rows.n_rows != 4) {
    Rcpp::stop("Composite T-squared trace Gram matrix requires four rows.");
  }
  arma::mat centered = rows;
  centered.each_row() -= rows.row(0);
  arma::mat gram(4, 4, arma::fill::zeros);
  for (const CompositeBlockFactor& block : partition.blocks) {
    arma::uvec indices(block.indices.size());
    for (std::size_t j = 0; j < block.indices.size(); ++j) {
      indices(j) = block.indices[j];
    }
    const arma::mat row_block = centered.cols(indices);
    const arma::mat solution = composite_solve_lower_cholesky(
      block, row_block.t(), stage
    );
    gram += row_block * solution;
  }
  if (!gram.is_finite()) {
    Rcpp::stop(
      "Composite T-squared produced a non-finite leave-four-out Gram "
      "matrix."
    );
  }
  return gram;
}

double composite_difference_bilinear(const arma::mat& gram,
                                     int a,
                                     int b,
                                     int c,
                                     int d) {
  return gram(a, c) - gram(a, d) - gram(b, c) + gram(b, d);
}

Rcpp::List composite_blocks_to_list(const CompositePartition& partition) {
  Rcpp::List result(partition.blocks.size());
  for (std::size_t b = 0; b < partition.blocks.size(); ++b) {
    Rcpp::IntegerVector indices(partition.blocks[b].indices.size());
    for (std::size_t j = 0; j < partition.blocks[b].indices.size(); ++j) {
      indices[j] = static_cast<int>(partition.blocks[b].indices[j] + 1);
    }
    result[b] = indices;
  }
  return result;
}

arma::vec composite_standardized_means(const CompositeSufficientData& data) {
  arma::vec result(data.p);
  const long double n = static_cast<long double>(data.n);
  for (arma::uword j = 0; j < data.p; ++j) {
    result(j) = static_cast<double>(data.sum[j] / n);
  }
  return result;
}

arma::vec composite_full_block_quadratics(
    const CompositePartition& partition,
    const arma::vec& difference) {
  arma::vec result(partition.blocks.size());
  for (std::size_t b = 0; b < partition.blocks.size(); ++b) {
    const CompositeBlockFactor& block = partition.blocks[b];
    arma::uvec indices(block.indices.size());
    for (std::size_t j = 0; j < block.indices.size(); ++j) {
      indices(j) = block.indices[j];
    }
    arma::mat right = difference.elem(indices);
    const arma::mat solution = composite_solve_lower_cholesky(
      block, right, "full-sample diagnostic"
    );
    result(b) = arma::dot(right, solution);
  }
  return result;
}

}  // namespace


//' Feng--Zou--Wang--Zhu Composite T-squared two-sample kernel
//'
//' @param x,y Numeric observation-by-variable matrices.
//' @param block_size Positive block size no larger than the dimension.
//' @param selection Block construction rule. Currently `"paper_greedy"`.
//' @return Internal list of statistic, trace, and partition diagnostics.
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List cpp_composite_t2_two_sample(const arma::mat& x,
                                       const arma::mat& y,
                                       int block_size = 2,
                                       std::string selection =
                                         "paper_greedy") {
  composite_require_finite_matrix(x, "x");
  composite_require_finite_matrix(y, "y");
  const arma::uword n1 = x.n_rows;
  const arma::uword n2 = y.n_rows;
  const arma::uword p = x.n_cols;
  if (p < 1 || y.n_cols != p) {
    Rcpp::stop("`x` and `y` must have the same positive number of columns.");
  }
  if (block_size < 1 || static_cast<arma::uword>(block_size) > p) {
    Rcpp::stop("`block_size` must be a positive integer no larger than p.");
  }
  if (selection != "paper_greedy") {
    Rcpp::stop("`selection` must be \"paper_greedy\".");
  }
  if (n1 < static_cast<arma::uword>(block_size + 5)) {
    Rcpp::stop(
      "Composite T-squared requires n1 >= block_size + 5 so every "
      "group-1 leave-four-out block can have full rank."
    );
  }
  if (n2 < 3) {
    Rcpp::stop(
      "Composite T-squared requires at least three observations in the "
      "second group for the pooled leave-2+2 covariance."
    );
  }

  const CompositeStandardizedData standardized =
    composite_standardize_columns(x, y);
  const CompositeSufficientData group1 =
    composite_sufficient_data(standardized.x);
  const CompositeSufficientData group2 =
    composite_sufficient_data(standardized.y);

  CompositeAccumulator q_half_ordered_numerator;
  CompositeLoopDiagnostics q_diagnostics;
  std::size_t interrupt_counter = 0;
  for (arma::uword i1 = 0; i1 + 1 < n1; ++i1) {
    for (arma::uword i2 = i1 + 1; i2 < n1; ++i2) {
      const std::array<arma::uword, 2> omitted1 = {i1, i2};
      for (arma::uword j1 = 0; j1 + 1 < n2; ++j1) {
        for (arma::uword j2 = j1 + 1; j2 < n2; ++j2) {
          if ((++interrupt_counter & 1023U) == 0U) {
            Rcpp::checkUserInterrupt();
          }
          const std::array<arma::uword, 2> omitted2 = {j1, j2};
          const arma::mat covariance = composite_pooled_covariance_leave2(
            group1, group2, omitted1, omitted2
          );
          const CompositePartition partition = composite_factor_partition(
            covariance, static_cast<arma::uword>(block_size),
            "pooled leave-2+2"
          );
          composite_update_loop_diagnostics(q_diagnostics, partition);

          arma::mat left(p, 2);
          arma::mat right(p, 2);
          left.col(0) = standardized.x.row(i1).t() -
            standardized.y.row(j1).t();
          right.col(0) = standardized.x.row(i2).t() -
            standardized.y.row(j2).t();
          left.col(1) = standardized.x.row(i1).t() -
            standardized.y.row(j2).t();
          right.col(1) = standardized.x.row(i2).t() -
            standardized.y.row(j1).t();
          const arma::rowvec kernels = composite_bilinear_columns(
            partition, left, right, "Q_n leave-2+2"
          );
          q_half_ordered_numerator.add(
            static_cast<long double>(kernels(0)) +
            static_cast<long double>(kernels(1))
          );
        }
      }
    }
  }

  CompositeAccumulator trace_ordered_numerator;
  CompositeLoopDiagnostics trace_diagnostics;
  std::array<int, 4> permutation = {0, 1, 2, 3};
  for (arma::uword i1 = 0; i1 + 3 < n1; ++i1) {
    for (arma::uword i2 = i1 + 1; i2 + 2 < n1; ++i2) {
      for (arma::uword i3 = i2 + 1; i3 + 1 < n1; ++i3) {
        for (arma::uword i4 = i3 + 1; i4 < n1; ++i4) {
          if ((++interrupt_counter & 1023U) == 0U) {
            Rcpp::checkUserInterrupt();
          }
          const std::array<arma::uword, 4> omitted = {i1, i2, i3, i4};
          const arma::mat covariance =
            composite_group1_covariance_leave4(group1, omitted);
          const CompositePartition partition = composite_factor_partition(
            covariance, static_cast<arma::uword>(block_size),
            "group-1 leave-four"
          );
          composite_update_loop_diagnostics(trace_diagnostics, partition);
          arma::mat rows(4, p);
          rows.row(0) = standardized.x.row(i1);
          rows.row(1) = standardized.x.row(i2);
          rows.row(2) = standardized.x.row(i3);
          rows.row(3) = standardized.x.row(i4);
          const arma::mat gram = composite_local_gram(
            partition, rows, "trace leave-four"
          );

          permutation = {0, 1, 2, 3};
          do {
            const double first = composite_difference_bilinear(
              gram, permutation[0], permutation[1],
              permutation[2], permutation[3]
            );
            const double second = composite_difference_bilinear(
              gram, permutation[0], permutation[3],
              permutation[2], permutation[1]
            );
            trace_ordered_numerator.add(
              static_cast<long double>(first) *
              static_cast<long double>(second)
            );
          } while (std::next_permutation(
            permutation.begin(), permutation.end()
          ));
        }
      }
    }
  }

  const long double q_pair_combinations =
    composite_choose2(n1) * composite_choose2(n2);
  const long double q_ordered_denominator =
    4.0L * q_pair_combinations;
  const long double q_ordered_numerator =
    2.0L * q_half_ordered_numerator.value;
  const long double q_statistic = q_ordered_numerator /
    q_ordered_denominator;

  const long double trace_quadruples = composite_choose4(n1);
  const long double p4_n1 = 24.0L * trace_quadruples;
  const long double trace_denominator = 2.0L * p4_n1;
  const long double trace_estimate = trace_ordered_numerator.value /
    trace_denominator;
  if (!std::isfinite(trace_estimate) || trace_estimate <= 0.0L) {
    Rcpp::stop(
      "Composite T-squared requires the published group-1 leave-four-out "
      "trace estimate to be finite and strictly positive; no absolute-value "
      "repair or variance floor is applied."
    );
  }

  const long double n1_inverse =
    1.0L / static_cast<long double>(n1);
  const long double n2_inverse =
    1.0L / static_cast<long double>(n2);
  const long double variance_coefficient =
    2.0L * (n1_inverse + n2_inverse) *
    (n1_inverse + n2_inverse);
  const long double variance = variance_coefficient * trace_estimate;
  if (!std::isfinite(variance) || variance <= 0.0L) {
    Rcpp::stop(
      "Composite T-squared requires a finite, strictly positive normalising "
      "variance; no absolute-value repair or variance floor is applied."
    );
  }
  const long double standard_error = std::sqrt(variance);
  const long double z = q_statistic / standard_error;

  const arma::mat full_covariance = composite_full_pooled_covariance(
    group1, group2
  );
  const CompositePartition full_partition = composite_factor_partition(
    full_covariance, static_cast<arma::uword>(block_size),
    "full pooled diagnostic"
  );
  const arma::vec mean1_scaled = composite_standardized_means(group1);
  const arma::vec mean2_scaled = composite_standardized_means(group2);
  const arma::vec difference_scaled = mean1_scaled - mean2_scaled;
  const arma::vec full_block_quadratic =
    composite_full_block_quadratics(full_partition, difference_scaled);

  arma::vec mean_x(p);
  arma::vec mean_y(p);
  arma::vec difference(p);
  arma::vec column_scale(p);
  arma::vec full_variance_scaled = full_covariance.diag();
  for (arma::uword j = 0; j < p; ++j) {
    const long double scale = standardized.scale[j];
    mean_x(j) = composite_report_long_double(
      standardized.origin[j] +
      scale * static_cast<long double>(mean1_scaled(j))
    );
    mean_y(j) = composite_report_long_double(
      standardized.origin[j] +
      scale * static_cast<long double>(mean2_scaled(j))
    );
    difference(j) = composite_report_long_double(
      scale * static_cast<long double>(difference_scaled(j))
    );
    column_scale(j) = composite_report_long_double(scale);
  }

  return Rcpp::List::create(
    Rcpp::Named("z") = composite_checked_double(
      z, "standardised statistic"
    ),
    Rcpp::Named("Q_n") = composite_checked_double(
      q_statistic, "Q_n statistic"
    ),
    Rcpp::Named("trace_lambda_squared") =
      composite_checked_positive_double(
        trace_estimate, "group-1 leave-four-out trace estimate"
      ),
    Rcpp::Named("variance_coefficient") = composite_checked_positive_double(
      variance_coefficient, "sample-size variance coefficient"
    ),
    Rcpp::Named("variance") = composite_checked_positive_double(
      variance, "normalising variance"
    ),
    Rcpp::Named("standard_error") = composite_checked_positive_double(
      standard_error, "normalising standard error"
    ),
    Rcpp::Named("mean_x") = mean_x,
    Rcpp::Named("mean_y") = mean_y,
    Rcpp::Named("difference") = difference,
    Rcpp::Named("difference_scaled") = difference_scaled,
    Rcpp::Named("column_scale") = column_scale,
    Rcpp::Named("full_pooled_variance_scaled") = full_variance_scaled,
    Rcpp::Named("full_blocks") = composite_blocks_to_list(full_partition),
    Rcpp::Named("full_block_quadratic") = full_block_quadratic,
    Rcpp::Named("full_partition_signature") = full_partition.signature,
    Rcpp::Named("q_ordered_numerator") = composite_checked_double(
      q_ordered_numerator, "ordered Q_n numerator"
    ),
    Rcpp::Named("q_ordered_denominator") = composite_checked_positive_double(
      q_ordered_denominator, "ordered Q_n denominator"
    ),
    Rcpp::Named("q_pair_combinations") = composite_checked_positive_double(
      q_pair_combinations, "unordered 2+2 combination count"
    ),
    Rcpp::Named("trace_ordered_numerator") = composite_checked_double(
      trace_ordered_numerator.value, "ordered trace numerator"
    ),
    Rcpp::Named("trace_ordered_denominator") =
      composite_checked_positive_double(
        trace_denominator, "ordered trace denominator"
      ),
    Rcpp::Named("trace_quadruples") = composite_checked_positive_double(
      trace_quadruples, "unordered group-1 quadruple count"
    ),
    Rcpp::Named("P4_n1") = composite_checked_positive_double(
      p4_n1, "P_n1^4"
    ),
    Rcpp::Named("q_unique_partitions") = static_cast<double>(
      q_diagnostics.partitions.size()
    ),
    Rcpp::Named("trace_unique_partitions") = static_cast<double>(
      trace_diagnostics.partitions.size()
    ),
    Rcpp::Named("q_minimum_reciprocal_condition") =
      q_diagnostics.minimum_reciprocal_condition,
    Rcpp::Named("trace_minimum_reciprocal_condition") =
      trace_diagnostics.minimum_reciprocal_condition,
    Rcpp::Named("full_minimum_reciprocal_condition") =
      full_partition.minimum_reciprocal_condition,
    Rcpp::Named("q_minimum_cholesky_diagonal") =
      q_diagnostics.minimum_cholesky_diagonal,
    Rcpp::Named("trace_minimum_cholesky_diagonal") =
      trace_diagnostics.minimum_cholesky_diagonal,
    Rcpp::Named("full_minimum_cholesky_diagonal") =
      full_partition.minimum_cholesky_diagonal,
    Rcpp::Named("q_minimum_marginal_variance_scaled") =
      q_diagnostics.minimum_marginal_variance,
    Rcpp::Named("trace_minimum_marginal_variance_scaled") =
      trace_diagnostics.minimum_marginal_variance,
    Rcpp::Named("full_minimum_marginal_variance_scaled") =
      full_partition.minimum_marginal_variance,
    Rcpp::Named("n1") = static_cast<double>(n1),
    Rcpp::Named("n2") = static_cast<double>(n2),
    Rcpp::Named("p") = static_cast<double>(p),
    Rcpp::Named("block_size") = block_size,
    Rcpp::Named("selection") = selection
  );
}
