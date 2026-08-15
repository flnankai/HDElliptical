// Feng--Zhang--Liu high-dimensional two-sample spatial-rank test.
//
// This implements the leave-two-out statistic and the fully feasible trace
// estimators in Feng, Zhang and Liu (2020).  Ordered sums are evaluated by
// symmetry over unordered pairs/four-subsets, without changing their
// published denominators.

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp17)]]

#include <RcppArmadillo.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <string>
#include <utility>
#include <vector>

namespace {

void fzl_add(long double value, long double& total,
             long double& correction) {
  const long double updated = total + value;
  if (std::abs(total) >= std::abs(value)) {
    correction += (total - updated) + value;
  } else {
    correction += (value - updated) + total;
  }
  total = updated;
}

long double fzl_dot(const arma::vec& x, const arma::vec& y) {
  long double total = 0.0L;
  long double correction = 0.0L;
  for (arma::uword j = 0; j < x.n_elem; ++j) {
    fzl_add(static_cast<long double>(x(j)) *
              static_cast<long double>(y(j)),
            total, correction);
  }
  return total + correction;
}

arma::vec fzl_stable_mean(const arma::mat& x) {
  arma::vec output(x.n_cols, arma::fill::zeros);
  for (arma::uword j = 0; j < x.n_cols; ++j) {
    double scale = 0.0;
    for (arma::uword i = 0; i < x.n_rows; ++i) {
      scale = std::max(scale, std::abs(x(i, j)));
    }
    if (scale == 0.0) {
      continue;
    }
    long double total = 0.0L;
    long double correction = 0.0L;
    for (arma::uword i = 0; i < x.n_rows; ++i) {
      fzl_add(static_cast<long double>(x(i, j) / scale),
              total, correction);
    }
    output(j) = static_cast<double>(
      (total + correction) / static_cast<long double>(x.n_rows)
    ) * scale;
  }
  return output;
}

struct FzlPrepared {
  arma::mat x;
  arma::mat y;
  arma::vec anchor;
  arma::vec log_coordinate_scale;
  arma::vec scaled_coordinate_range;
  arma::vec mean_x;
  arma::vec mean_y;
  double global_operand_scale;
};

FzlPrepared fzl_prepare(const arma::mat& x, const arma::mat& y,
                        bool geometric) {
  const arma::uword n1 = x.n_rows;
  const arma::uword n2 = y.n_rows;
  const arma::uword p = x.n_cols;
  arma::mat xs(n1, p, arma::fill::zeros);
  arma::mat ys(n2, p, arma::fill::zeros);
  arma::vec anchor = x.row(0).t();
  arma::vec log_scale(p, arma::fill::zeros);
  arma::vec coordinate_range(p, arma::fill::zeros);

  double global_scale = 0.0;
  for (arma::uword j = 0; j < p; ++j) {
    for (arma::uword i = 0; i < n1; ++i) {
      global_scale = std::max(global_scale, std::abs(x(i, j)));
    }
    for (arma::uword i = 0; i < n2; ++i) {
      global_scale = std::max(global_scale, std::abs(y(i, j)));
    }
  }
  if (!(global_scale > 0.0) || !std::isfinite(global_scale)) {
    Rcpp::stop(
      "Feng-Zhang-Liu requires nonconstant finite pooled data."
    );
  }

  for (arma::uword j = 0; j < p; ++j) {
    double operand_scale = global_scale;
    if (geometric) {
      operand_scale = std::abs(anchor(j));
      for (arma::uword i = 0; i < n1; ++i) {
        operand_scale = std::max(operand_scale, std::abs(x(i, j)));
      }
      for (arma::uword i = 0; i < n2; ++i) {
        operand_scale = std::max(operand_scale, std::abs(y(i, j)));
      }
    }
    if (!(operand_scale > 0.0) || !std::isfinite(operand_scale)) {
      Rcpp::stop(
        "Feng-Zhang-Liu requires pooled positive variation in variable "
        "%llu.", static_cast<unsigned long long>(j + 1)
      );
    }

    const double anchor_scaled = anchor(j) / operand_scale;
    double max_difference = 0.0;
    for (arma::uword i = 0; i < n1; ++i) {
      xs(i, j) = x(i, j) / operand_scale - anchor_scaled;
      max_difference = std::max(max_difference, std::abs(xs(i, j)));
    }
    for (arma::uword i = 0; i < n2; ++i) {
      ys(i, j) = y(i, j) / operand_scale - anchor_scaled;
      max_difference = std::max(max_difference, std::abs(ys(i, j)));
    }
    if (!(max_difference > 0.0) || !std::isfinite(max_difference)) {
      Rcpp::stop(
        "Feng-Zhang-Liu requires pooled positive variation in variable "
        "%llu; the variation is zero or is not representable under the "
        "selected scale identification.",
        static_cast<unsigned long long>(j + 1)
      );
    }
    coordinate_range(j) = max_difference;
    if (geometric) {
      xs.col(j) /= max_difference;
      ys.col(j) /= max_difference;
      log_scale(j) = std::log(operand_scale) +
        std::log(max_difference);
    } else {
      log_scale(j) = std::log(global_scale);
    }
  }

  return FzlPrepared{
    xs, ys, anchor, log_scale, coordinate_range,
    fzl_stable_mean(x), fzl_stable_mean(y), global_scale
  };
}

double fzl_log_sum_exp(const arma::vec& values) {
  const double maximum = values.max();
  if (!std::isfinite(maximum)) {
    Rcpp::stop("Feng-Zhang-Liu scale normalisation became non-finite.");
  }
  long double total = 0.0L;
  long double correction = 0.0L;
  for (arma::uword j = 0; j < values.n_elem; ++j) {
    fzl_add(std::exp(static_cast<long double>(values(j) - maximum)),
            total, correction);
  }
  const long double sum = total + correction;
  if (!(sum > 0.0L) || !std::isfinite(sum)) {
    Rcpp::stop("Feng-Zhang-Liu scale normalisation is undefined.");
  }
  return maximum + std::log(static_cast<double>(sum));
}

void fzl_identify(arma::vec& log_diagonal, bool geometric) {
  double offset = 0.0;
  if (geometric) {
    long double total = 0.0L;
    long double correction = 0.0L;
    for (arma::uword j = 0; j < log_diagonal.n_elem; ++j) {
      fzl_add(static_cast<long double>(log_diagonal(j)),
              total, correction);
    }
    offset = static_cast<double>(
      (total + correction) /
        static_cast<long double>(log_diagonal.n_elem)
    );
  } else {
    offset = fzl_log_sum_exp(log_diagonal) -
      std::log(static_cast<double>(log_diagonal.n_elem));
  }
  log_diagonal -= offset;
  if (!log_diagonal.is_finite()) {
    Rcpp::stop(
      "Feng-Zhang-Liu diagonal scale became non-finite; no floor, ridge, "
      "or perturbation is applied."
    );
  }
}

bool fzl_direction(const arma::vec& difference,
                   const arma::vec& log_diagonal,
                   arma::vec& direction,
                   double& log_radius) {
  const arma::uword p = difference.n_elem;
  double maximum_log = -std::numeric_limits<double>::infinity();
  for (arma::uword j = 0; j < p; ++j) {
    if (difference(j) == 0.0) {
      direction(j) = 0.0;
      continue;
    }
    const double component_log = std::log(std::abs(difference(j))) -
      0.5 * log_diagonal(j);
    if (!std::isfinite(component_log)) {
      Rcpp::stop(
        "Feng-Zhang-Liu diagonal standardisation became non-finite; no "
        "numerical floor is applied."
      );
    }
    direction(j) = std::copysign(1.0, difference(j));
    maximum_log = std::max(maximum_log, component_log);
  }
  if (!std::isfinite(maximum_log)) {
    direction.zeros();
    log_radius = -std::numeric_limits<double>::infinity();
    return false;
  }

  long double square_total = 0.0L;
  long double square_correction = 0.0L;
  for (arma::uword j = 0; j < p; ++j) {
    if (difference(j) == 0.0) {
      continue;
    }
    const double component_log = std::log(std::abs(difference(j))) -
      0.5 * log_diagonal(j);
    direction(j) *= std::exp(component_log - maximum_log);
    const long double value = static_cast<long double>(direction(j));
    fzl_add(value * value, square_total, square_correction);
  }
  const long double squared_norm = square_total + square_correction;
  if (!(squared_norm > 0.0L) || !std::isfinite(squared_norm)) {
    Rcpp::stop("Feng-Zhang-Liu could not normalise a spatial direction.");
  }
  const double norm = std::sqrt(static_cast<double>(squared_norm));
  direction /= norm;
  log_radius = maximum_log + std::log(norm);
  return true;
}

bool fzl_row_difference_direction(const arma::mat& data,
                                  arma::uword first,
                                  arma::uword second,
                                  const arma::vec& log_diagonal,
                                  arma::vec& direction,
                                  double& log_radius) {
  return fzl_direction(
    data.row(first).t() - data.row(second).t(),
    log_diagonal, direction, log_radius
  );
}

bool fzl_cross_direction(const arma::mat& x, arma::uword i,
                         const arma::mat& y, arma::uword j,
                         const arma::vec& log_diagonal,
                         arma::vec& direction,
                         double& log_radius) {
  return fzl_direction(
    x.row(i).t() - y.row(j).t(), log_diagonal,
    direction, log_radius
  );
}

struct FzlScaleScore {
  arma::vec diagonal;
  double residual;
  arma::uword zero_pair_directions;
  double minimum_log_pair_radius;
};

FzlScaleScore fzl_scale_score(const arma::mat& data,
                              const std::vector<bool>& excluded,
                              const arma::vec& log_diagonal) {
  const arma::uword n = data.n_rows;
  const arma::uword p = data.n_cols;
  arma::uword m = 0;
  std::vector<arma::uword> keep;
  keep.reserve(n);
  for (arma::uword i = 0; i < n; ++i) {
    if (!excluded[i]) {
      keep.push_back(i);
      ++m;
    }
  }
  arma::mat ranks(m, p, arma::fill::zeros);
  arma::vec direction(p, arma::fill::zeros);
  arma::uword zeros = 0;
  double minimum_log_radius = std::numeric_limits<double>::infinity();
  const double m_double = static_cast<double>(m);

  for (arma::uword a = 0; a < m; ++a) {
    for (arma::uword b = a + 1; b < m; ++b) {
      double log_radius = 0.0;
      if (fzl_row_difference_direction(
            data, keep[a], keep[b], log_diagonal,
            direction, log_radius
          )) {
        ranks.row(a) += (direction / m_double).t();
        ranks.row(b) -= (direction / m_double).t();
        minimum_log_radius = std::min(minimum_log_radius, log_radius);
      } else {
        ++zeros;
      }
    }
  }

  arma::vec score = arma::mean(arma::square(ranks), 0).t();
  if (!score.is_finite() || arma::any(score <= 0.0)) {
    Rcpp::stop(
      "Feng-Zhang-Liu diagonal spatial-rank scale has non-positive "
      "coordinate rank energy; no ridge, floor, or perturbation is applied."
    );
  }
  const double score_sum = arma::sum(score);
  if (!(score_sum > 0.0) || !std::isfinite(score_sum)) {
    Rcpp::stop("Feng-Zhang-Liu spatial-rank scale score is degenerate.");
  }
  const double residual = arma::max(arma::abs(
    static_cast<double>(p) * score / score_sum - 1.0
  ));
  return FzlScaleScore{score, residual, zeros, minimum_log_radius};
}

struct FzlScaleFit {
  arma::vec log_diagonal;
  int iterations;
  double residual;
  arma::uword zero_pair_directions;
  double minimum_log_pair_radius;
};

std::string fzl_exclusion_label(const std::vector<arma::uword>& excluded) {
  if (excluded.empty()) {
    return "full sample";
  }
  std::string output = "leave-out indices ";
  for (std::size_t i = 0; i < excluded.size(); ++i) {
    if (i > 0) {
      output += ",";
    }
    output += std::to_string(excluded[i] + 1);
  }
  return output;
}

FzlScaleFit fzl_fit_scale(const arma::mat& data,
                          const std::vector<arma::uword>& exclusions,
                          double tolerance, int max_iterations,
                          bool geometric, int group) {
  const arma::uword n = data.n_rows;
  const arma::uword p = data.n_cols;
  std::vector<bool> excluded(n, false);
  for (arma::uword index : exclusions) {
    excluded[index] = true;
  }
  const arma::uword m = n - exclusions.size();
  if (m < 2) {
    Rcpp::stop(
      "Feng-Zhang-Liu scale fit for group %d (%s) needs at least two "
      "retained observations.", group,
      fzl_exclusion_label(exclusions).c_str()
    );
  }

  arma::vec mean(p, arma::fill::zeros);
  for (arma::uword i = 0; i < n; ++i) {
    if (!excluded[i]) {
      mean += data.row(i).t() / static_cast<double>(m);
    }
  }
  arma::vec variance(p, arma::fill::zeros);
  for (arma::uword i = 0; i < n; ++i) {
    if (!excluded[i]) {
      variance += arma::square(data.row(i).t() - mean);
    }
  }
  variance /= static_cast<double>(m - 1);
  if (!variance.is_finite() || arma::any(variance <= 0.0)) {
    Rcpp::stop(
      "Feng-Zhang-Liu initial marginal variance is non-positive for group "
      "%d (%s); no ridge or perturbation is applied.", group,
      fzl_exclusion_label(exclusions).c_str()
    );
  }
  arma::vec log_diagonal = arma::log(variance);
  fzl_identify(log_diagonal, geometric);

  for (int iteration = 0; iteration <= max_iterations; ++iteration) {
    const FzlScaleScore score = fzl_scale_score(
      data, excluded, log_diagonal
    );
    if (score.residual <= tolerance) {
      return FzlScaleFit{
        log_diagonal, iteration, score.residual,
        score.zero_pair_directions, score.minimum_log_pair_radius
      };
    }
    if (iteration == max_iterations) {
      Rcpp::stop(
        "Feng-Zhang-Liu diagonal spatial-rank scale failed to converge "
        "for group %d (%s) after %d updates; final residual %.17g exceeds "
        "tol %.17g.", group, fzl_exclusion_label(exclusions).c_str(),
        max_iterations, score.residual, tolerance
      );
    }
    log_diagonal += arma::log(score.diagonal);
    fzl_identify(log_diagonal, geometric);
  }
  Rcpp::stop("Unreachable Feng-Zhang-Liu scale-fit state.");
}

std::vector<std::array<arma::uword, 2>> fzl_pairs(arma::uword n) {
  std::vector<std::array<arma::uword, 2>> output;
  output.reserve(static_cast<std::size_t>(n * (n - 1) / 2));
  for (arma::uword i = 0; i < n; ++i) {
    for (arma::uword j = i + 1; j < n; ++j) {
      output.push_back({i, j});
    }
  }
  return output;
}

std::vector<std::array<arma::uword, 4>> fzl_quads(arma::uword n) {
  std::vector<std::array<arma::uword, 4>> output;
  for (arma::uword i = 0; i < n; ++i) {
    for (arma::uword j = i + 1; j < n; ++j) {
      for (arma::uword k = j + 1; k < n; ++k) {
        for (arma::uword l = k + 1; l < n; ++l) {
          output.push_back({i, j, k, l});
        }
      }
    }
  }
  return output;
}

arma::vec fzl_pooled_log_diagonal(const arma::vec& first,
                                  const arma::vec& second,
                                  double first_weight,
                                  double second_weight) {
  arma::vec output(first.n_elem, arma::fill::zeros);
  const double log_w1 = std::log(first_weight);
  const double log_w2 = std::log(second_weight);
  for (arma::uword j = 0; j < first.n_elem; ++j) {
    const double a = log_w1 + first(j);
    const double b = log_w2 + second(j);
    const double maximum = std::max(a, b);
    output(j) = maximum +
      std::log(std::exp(a - maximum) + std::exp(b - maximum));
  }
  return output;
}

arma::mat fzl_pair_indices(
    const std::vector<std::array<arma::uword, 2>>& pairs) {
  arma::mat output(pairs.size(), 2, arma::fill::zeros);
  for (std::size_t r = 0; r < pairs.size(); ++r) {
    output(r, 0) = static_cast<double>(pairs[r][0] + 1);
    output(r, 1) = static_cast<double>(pairs[r][1] + 1);
  }
  return output;
}

arma::mat fzl_quad_indices(
    const std::vector<std::array<arma::uword, 4>>& quads) {
  arma::mat output(quads.size(), 4, arma::fill::zeros);
  for (std::size_t r = 0; r < quads.size(); ++r) {
    for (arma::uword j = 0; j < 4; ++j) {
      output(r, j) = static_cast<double>(quads[r][j] + 1);
    }
  }
  return output;
}

arma::mat fzl_fit_log_matrix(const std::vector<FzlScaleFit>& fits,
                             arma::uword p) {
  arma::mat output(fits.size(), p, arma::fill::zeros);
  for (std::size_t i = 0; i < fits.size(); ++i) {
    output.row(i) = fits[i].log_diagonal.t();
  }
  return output;
}

arma::mat fzl_exp_matrix(const arma::mat& logs) {
  return arma::exp(logs);
}

arma::mat fzl_canonical_matrix(const arma::mat& logs) {
  arma::mat output(logs.n_rows, logs.n_cols, arma::fill::zeros);
  for (arma::uword i = 0; i < logs.n_rows; ++i) {
    const double maximum = logs.row(i).max();
    output.row(i) = arma::exp(logs.row(i) - maximum);
  }
  return output;
}

arma::vec fzl_fit_iterations(const std::vector<FzlScaleFit>& fits) {
  arma::vec output(fits.size(), arma::fill::zeros);
  for (std::size_t i = 0; i < fits.size(); ++i) {
    output(i) = fits[i].iterations;
  }
  return output;
}

arma::vec fzl_fit_residuals(const std::vector<FzlScaleFit>& fits) {
  arma::vec output(fits.size(), arma::fill::zeros);
  for (std::size_t i = 0; i < fits.size(); ++i) {
    output(i) = fits[i].residual;
  }
  return output;
}

arma::vec fzl_fit_zeros(const std::vector<FzlScaleFit>& fits) {
  arma::vec output(fits.size(), arma::fill::zeros);
  for (std::size_t i = 0; i < fits.size(); ++i) {
    output(i) = static_cast<double>(fits[i].zero_pair_directions);
  }
  return output;
}

arma::vec fzl_fit_min_log_radii(const std::vector<FzlScaleFit>& fits) {
  arma::vec output(fits.size(), arma::fill::zeros);
  for (std::size_t i = 0; i < fits.size(); ++i) {
    output(i) = fits[i].minimum_log_pair_radius;
  }
  return output;
}

arma::vec fzl_input_log_canonical(const arma::vec& log_diagonal,
                                  const arma::vec& log_coordinate_scale) {
  arma::vec output = log_diagonal + 2.0 * log_coordinate_scale;
  output -= output.max();
  return output;
}

arma::mat fzl_input_log_canonical_matrix(
    const arma::mat& logs, const arma::vec& log_coordinate_scale) {
  arma::mat output(logs.n_rows, logs.n_cols, arma::fill::zeros);
  for (arma::uword i = 0; i < logs.n_rows; ++i) {
    arma::vec row = logs.row(i).t() + 2.0 * log_coordinate_scale;
    row -= row.max();
    output.row(i) = row.t();
  }
  return output;
}

struct FzlWithinTrace {
  double estimate;
  double ordered_sum;
  arma::mat contributions;
  arma::uword zero_direction_evaluations;
};

FzlWithinTrace fzl_within_trace(
    const arma::mat& data,
    const std::vector<std::array<arma::uword, 4>>& quads,
    const std::vector<FzlScaleFit>& fits) {
  const arma::uword p = data.n_cols;
  arma::mat contributions(quads.size(), 24, arma::fill::zeros);
  long double total = 0.0L;
  long double correction = 0.0L;
  arma::uword zeros = 0;
  arma::vec u12(p), u34(p), u32(p), u14(p);

  for (std::size_t q = 0; q < quads.size(); ++q) {
    std::array<arma::uword, 4> order = quads[q];
    int permutation = 0;
    do {
      double log_radius = 0.0;
      if (!fzl_row_difference_direction(
            data, order[0], order[1], fits[q].log_diagonal,
            u12, log_radius)) ++zeros;
      if (!fzl_row_difference_direction(
            data, order[2], order[3], fits[q].log_diagonal,
            u34, log_radius)) ++zeros;
      if (!fzl_row_difference_direction(
            data, order[2], order[1], fits[q].log_diagonal,
            u32, log_radius)) ++zeros;
      if (!fzl_row_difference_direction(
            data, order[0], order[3], fits[q].log_diagonal,
            u14, log_radius)) ++zeros;
      const long double value = fzl_dot(u12, u34) *
        fzl_dot(u32, u14);
      contributions(q, permutation) = static_cast<double>(value);
      fzl_add(value, total, correction);
      ++permutation;
    } while (std::next_permutation(order.begin(), order.end()));
  }
  const long double ordered_sum = total + correction;
  const long double n = static_cast<long double>(data.n_rows);
  const long double falling_four = n * (n - 1.0L) *
    (n - 2.0L) * (n - 3.0L);
  const long double p_ld = static_cast<long double>(p);
  const long double estimate = 2.0L * p_ld * p_ld *
    ordered_sum / falling_four;
  return FzlWithinTrace{
    static_cast<double>(estimate), static_cast<double>(ordered_sum),
    contributions, zeros
  };
}

} // namespace


// [[Rcpp::export]]
Rcpp::List cpp_feng_zhang_liu_spatial_rank(
    const arma::mat& x, const arma::mat& y, const double tolerance,
    const int max_iterations, const std::string scale_identification) {
  if (x.n_rows < 6 || y.n_rows < 6 || x.n_cols < 1 ||
      y.n_cols != x.n_cols) {
    Rcpp::stop(
      "Feng-Zhang-Liu requires two samples with at least six rows each "
      "and the same positive number of columns."
    );
  }
  if (!x.is_finite() || !y.is_finite()) {
    Rcpp::stop("Feng-Zhang-Liu inputs must contain only finite values.");
  }
  if (!(tolerance > 0.0) || !std::isfinite(tolerance) ||
      max_iterations < 1) {
    Rcpp::stop("Invalid Feng-Zhang-Liu iteration controls.");
  }
  const bool geometric = scale_identification == "geometric";
  if (!geometric && scale_identification != "paper_trace") {
    Rcpp::stop(
      "`scale_identification` must be \"geometric\" or \"paper_trace\"."
    );
  }

  const arma::uword n1 = x.n_rows;
  const arma::uword n2 = y.n_rows;
  const arma::uword p = x.n_cols;
  const FzlPrepared prepared = fzl_prepare(x, y, geometric);
  const FzlScaleFit full1 = fzl_fit_scale(
    prepared.x, {}, tolerance, max_iterations, geometric, 1
  );
  const FzlScaleFit full2 = fzl_fit_scale(
    prepared.y, {}, tolerance, max_iterations, geometric, 2
  );

  const auto pairs1 = fzl_pairs(n1);
  const auto pairs2 = fzl_pairs(n2);
  const auto quads1 = fzl_quads(n1);
  const auto quads2 = fzl_quads(n2);
  std::vector<FzlScaleFit> pair_fits1;
  std::vector<FzlScaleFit> pair_fits2;
  std::vector<FzlScaleFit> quad_fits1;
  std::vector<FzlScaleFit> quad_fits2;
  pair_fits1.reserve(pairs1.size());
  pair_fits2.reserve(pairs2.size());
  quad_fits1.reserve(quads1.size());
  quad_fits2.reserve(quads2.size());

  for (const auto& pair : pairs1) {
    pair_fits1.push_back(fzl_fit_scale(
      prepared.x, {pair[0], pair[1]}, tolerance, max_iterations,
      geometric, 1
    ));
  }
  for (const auto& pair : pairs2) {
    pair_fits2.push_back(fzl_fit_scale(
      prepared.y, {pair[0], pair[1]}, tolerance, max_iterations,
      geometric, 2
    ));
  }
  for (const auto& quad : quads1) {
    quad_fits1.push_back(fzl_fit_scale(
      prepared.x, {quad[0], quad[1], quad[2], quad[3]},
      tolerance, max_iterations, geometric, 1
    ));
  }
  for (const auto& quad : quads2) {
    quad_fits2.push_back(fzl_fit_scale(
      prepared.y, {quad[0], quad[1], quad[2], quad[3]},
      tolerance, max_iterations, geometric, 2
    ));
  }

  const double weight1 = static_cast<double>(n1) /
    static_cast<double>(n1 + n2);
  const double weight2 = static_cast<double>(n2) /
    static_cast<double>(n1 + n2);
  arma::mat main_blocks(pairs1.size(), pairs2.size(), arma::fill::zeros);
  arma::mat trace3_blocks(pairs1.size(), pairs2.size(), arma::fill::zeros);
  long double main_total = 0.0L;
  long double main_correction = 0.0L;
  long double trace3_total = 0.0L;
  long double trace3_correction = 0.0L;
  arma::uword main_zero_evaluations = 0;
  arma::uword trace3_zero_evaluations = 0;
  arma::vec xis(p), xjt(p), xit(p), xjs(p), xpair(p), ypair(p);

  for (std::size_t a = 0; a < pairs1.size(); ++a) {
    for (std::size_t b = 0; b < pairs2.size(); ++b) {
      const arma::vec pooled = fzl_pooled_log_diagonal(
        pair_fits1[a].log_diagonal, pair_fits2[b].log_diagonal,
        weight1, weight2
      );
      const arma::uword i = pairs1[a][0];
      const arma::uword j = pairs1[a][1];
      const arma::uword s = pairs2[b][0];
      const arma::uword t = pairs2[b][1];
      double log_radius = 0.0;
      if (!fzl_cross_direction(
            prepared.x, i, prepared.y, s, pooled, xis, log_radius
          )) ++main_zero_evaluations;
      if (!fzl_cross_direction(
            prepared.x, j, prepared.y, t, pooled, xjt, log_radius
          )) ++main_zero_evaluations;
      if (!fzl_cross_direction(
            prepared.x, i, prepared.y, t, pooled, xit, log_radius
          )) ++main_zero_evaluations;
      if (!fzl_cross_direction(
            prepared.x, j, prepared.y, s, pooled, xjs, log_radius
          )) ++main_zero_evaluations;
      const long double main_value = 2.0L *
        (fzl_dot(xis, xjt) + fzl_dot(xit, xjs));
      main_blocks(a, b) = static_cast<double>(main_value);
      fzl_add(main_value, main_total, main_correction);

      if (!fzl_row_difference_direction(
            prepared.x, i, j, pooled, xpair, log_radius
          )) ++trace3_zero_evaluations;
      if (!fzl_row_difference_direction(
            prepared.y, s, t, pooled, ypair, log_radius
          )) ++trace3_zero_evaluations;
      const long double inner = fzl_dot(xpair, ypair);
      const long double trace3_value = 4.0L * inner * inner;
      trace3_blocks(a, b) = static_cast<double>(trace3_value);
      fzl_add(trace3_value, trace3_total, trace3_correction);
    }
  }

  const long double main_ordered_sum = main_total + main_correction;
  const long double main_denominator =
    static_cast<long double>(n1) * (n1 - 1.0L) *
    static_cast<long double>(n2) * (n2 - 1.0L);
  const double statistic_raw = static_cast<double>(
    main_ordered_sum / main_denominator
  );

  const FzlWithinTrace trace1 = fzl_within_trace(
    prepared.x, quads1, quad_fits1
  );
  const FzlWithinTrace trace2 = fzl_within_trace(
    prepared.y, quads2, quad_fits2
  );
  const long double trace3_ordered_sum =
    trace3_total + trace3_correction;
  const long double p_squared = static_cast<long double>(p) *
    static_cast<long double>(p);
  const long double trace3_denominator =
    static_cast<long double>(n1) * static_cast<long double>(n1) *
    static_cast<long double>(n2) * static_cast<long double>(n2);
  const double trace3 = static_cast<double>(
    p_squared * trace3_ordered_sum / trace3_denominator
  );

  const double p2 = static_cast<double>(p) * static_cast<double>(p);
  const double variance_term1 = trace1.estimate /
    (2.0 * static_cast<double>(n1) * static_cast<double>(n1 - 1) * p2);
  // The p^2 below restores the factor omitted typographically from the
  // second term of the author TeX; it is present in the symmetric oracle
  // formula and is required for consistency of the feasible estimator.
  const double variance_term2 = trace2.estimate /
    (2.0 * static_cast<double>(n2) * static_cast<double>(n2 - 1) * p2);
  const double variance_term3 = trace3 /
    (static_cast<double>(n1) * static_cast<double>(n2) * p2);
  const double sigma2 = variance_term1 + variance_term2 + variance_term3;
  if (!(sigma2 > 0.0) || !std::isfinite(sigma2)) {
    Rcpp::stop(
      "Feng-Zhang-Liu feasible variance is non-positive or non-finite "
      "(%.17g); no absolute value, ridge, floor, or perturbation is applied.",
      sigma2
    );
  }
  const double sigma = std::sqrt(sigma2);
  const double z = statistic_raw / sigma;
  if (!std::isfinite(z)) {
    Rcpp::stop("Feng-Zhang-Liu standardised statistic is non-finite.");
  }

  const arma::mat pair_logs1 = fzl_fit_log_matrix(pair_fits1, p);
  const arma::mat pair_logs2 = fzl_fit_log_matrix(pair_fits2, p);
  const arma::mat quad_logs1 = fzl_fit_log_matrix(quad_fits1, p);
  const arma::mat quad_logs2 = fzl_fit_log_matrix(quad_fits2, p);
  const arma::vec full_input_log1 = fzl_input_log_canonical(
    full1.log_diagonal, prepared.log_coordinate_scale
  );
  const arma::vec full_input_log2 = fzl_input_log_canonical(
    full2.log_diagonal, prepared.log_coordinate_scale
  );
  const arma::mat pair_input_logs1 = fzl_input_log_canonical_matrix(
    pair_logs1, prepared.log_coordinate_scale
  );
  const arma::mat pair_input_logs2 = fzl_input_log_canonical_matrix(
    pair_logs2, prepared.log_coordinate_scale
  );
  const arma::mat quad_input_logs1 = fzl_input_log_canonical_matrix(
    quad_logs1, prepared.log_coordinate_scale
  );
  const arma::mat quad_input_logs2 = fzl_input_log_canonical_matrix(
    quad_logs2, prepared.log_coordinate_scale
  );

  return Rcpp::List::create(
    Rcpp::Named("T_n") = statistic_raw,
    Rcpp::Named("z") = z,
    Rcpp::Named("sigma2_hat") = sigma2,
    Rcpp::Named("sigma_hat") = sigma,
    Rcpp::Named("trace_R1_squared_hat") = trace1.estimate,
    Rcpp::Named("trace_R2_squared_hat") = trace2.estimate,
    Rcpp::Named("trace_R1_R2_hat") = trace3,
    Rcpp::Named("variance_term1") = variance_term1,
    Rcpp::Named("variance_term2") = variance_term2,
    Rcpp::Named("variance_term3") = variance_term3,
    Rcpp::Named("main_ordered_sum") =
      static_cast<double>(main_ordered_sum),
    Rcpp::Named("main_ordered_denominator") =
      static_cast<double>(main_denominator),
    Rcpp::Named("main_pair_block_contributions") = main_blocks,
    Rcpp::Named("trace1_ordered_sum") = trace1.ordered_sum,
    Rcpp::Named("trace2_ordered_sum") = trace2.ordered_sum,
    Rcpp::Named("trace3_ordered_sum") =
      static_cast<double>(trace3_ordered_sum),
    Rcpp::Named("trace1_ordered_denominator") =
      static_cast<double>(n1) * (n1 - 1.0) * (n1 - 2.0) * (n1 - 3.0),
    Rcpp::Named("trace2_ordered_denominator") =
      static_cast<double>(n2) * (n2 - 1.0) * (n2 - 2.0) * (n2 - 3.0),
    Rcpp::Named("trace3_ordered_denominator") =
      static_cast<double>(trace3_denominator),
    Rcpp::Named("trace1_permutation_contributions") =
      trace1.contributions,
    Rcpp::Named("trace2_permutation_contributions") =
      trace2.contributions,
    Rcpp::Named("trace3_pair_block_contributions") = trace3_blocks,
    Rcpp::Named("pair_indices1") = fzl_pair_indices(pairs1),
    Rcpp::Named("pair_indices2") = fzl_pair_indices(pairs2),
    Rcpp::Named("quad_indices1") = fzl_quad_indices(quads1),
    Rcpp::Named("quad_indices2") = fzl_quad_indices(quads2),
    Rcpp::Named("full_D1_standardised") = arma::exp(full1.log_diagonal),
    Rcpp::Named("full_D2_standardised") = arma::exp(full2.log_diagonal),
    Rcpp::Named("full_log_D1_standardised") = full1.log_diagonal,
    Rcpp::Named("full_log_D2_standardised") = full2.log_diagonal,
    Rcpp::Named("full_D1_input_canonical") = arma::exp(full_input_log1),
    Rcpp::Named("full_D2_input_canonical") = arma::exp(full_input_log2),
    Rcpp::Named("full_log_D1_input_canonical") = full_input_log1,
    Rcpp::Named("full_log_D2_input_canonical") = full_input_log2,
    Rcpp::Named("full_iterations1") = full1.iterations,
    Rcpp::Named("full_iterations2") = full2.iterations,
    Rcpp::Named("full_residual1") = full1.residual,
    Rcpp::Named("full_residual2") = full2.residual,
    Rcpp::Named("full_zero_pair_directions1") =
      static_cast<double>(full1.zero_pair_directions),
    Rcpp::Named("full_zero_pair_directions2") =
      static_cast<double>(full2.zero_pair_directions),
    Rcpp::Named("full_minimum_log_pair_radius1") =
      full1.minimum_log_pair_radius,
    Rcpp::Named("full_minimum_log_pair_radius2") =
      full2.minimum_log_pair_radius,
    Rcpp::Named("leave2_D1_standardised") = fzl_exp_matrix(pair_logs1),
    Rcpp::Named("leave2_D2_standardised") = fzl_exp_matrix(pair_logs2),
    Rcpp::Named("leave2_log_D1_standardised") = pair_logs1,
    Rcpp::Named("leave2_log_D2_standardised") = pair_logs2,
    Rcpp::Named("leave2_D1_standardised_canonical") =
      fzl_canonical_matrix(pair_logs1),
    Rcpp::Named("leave2_D2_standardised_canonical") =
      fzl_canonical_matrix(pair_logs2),
    Rcpp::Named("leave2_D1_input_canonical") =
      arma::exp(pair_input_logs1),
    Rcpp::Named("leave2_D2_input_canonical") =
      arma::exp(pair_input_logs2),
    Rcpp::Named("leave2_log_D1_input_canonical") = pair_input_logs1,
    Rcpp::Named("leave2_log_D2_input_canonical") = pair_input_logs2,
    Rcpp::Named("leave2_iterations1") = fzl_fit_iterations(pair_fits1),
    Rcpp::Named("leave2_iterations2") = fzl_fit_iterations(pair_fits2),
    Rcpp::Named("leave2_residuals1") = fzl_fit_residuals(pair_fits1),
    Rcpp::Named("leave2_residuals2") = fzl_fit_residuals(pair_fits2),
    Rcpp::Named("leave2_zero_pair_directions1") = fzl_fit_zeros(pair_fits1),
    Rcpp::Named("leave2_zero_pair_directions2") = fzl_fit_zeros(pair_fits2),
    Rcpp::Named("leave2_minimum_log_pair_radii1") =
      fzl_fit_min_log_radii(pair_fits1),
    Rcpp::Named("leave2_minimum_log_pair_radii2") =
      fzl_fit_min_log_radii(pair_fits2),
    Rcpp::Named("leave4_D1_standardised") = fzl_exp_matrix(quad_logs1),
    Rcpp::Named("leave4_D2_standardised") = fzl_exp_matrix(quad_logs2),
    Rcpp::Named("leave4_log_D1_standardised") = quad_logs1,
    Rcpp::Named("leave4_log_D2_standardised") = quad_logs2,
    Rcpp::Named("leave4_D1_standardised_canonical") =
      fzl_canonical_matrix(quad_logs1),
    Rcpp::Named("leave4_D2_standardised_canonical") =
      fzl_canonical_matrix(quad_logs2),
    Rcpp::Named("leave4_D1_input_canonical") =
      arma::exp(quad_input_logs1),
    Rcpp::Named("leave4_D2_input_canonical") =
      arma::exp(quad_input_logs2),
    Rcpp::Named("leave4_log_D1_input_canonical") = quad_input_logs1,
    Rcpp::Named("leave4_log_D2_input_canonical") = quad_input_logs2,
    Rcpp::Named("leave4_iterations1") = fzl_fit_iterations(quad_fits1),
    Rcpp::Named("leave4_iterations2") = fzl_fit_iterations(quad_fits2),
    Rcpp::Named("leave4_residuals1") = fzl_fit_residuals(quad_fits1),
    Rcpp::Named("leave4_residuals2") = fzl_fit_residuals(quad_fits2),
    Rcpp::Named("leave4_zero_pair_directions1") = fzl_fit_zeros(quad_fits1),
    Rcpp::Named("leave4_zero_pair_directions2") = fzl_fit_zeros(quad_fits2),
    Rcpp::Named("leave4_minimum_log_pair_radii1") =
      fzl_fit_min_log_radii(quad_fits1),
    Rcpp::Named("leave4_minimum_log_pair_radii2") =
      fzl_fit_min_log_radii(quad_fits2),
    Rcpp::Named("main_zero_direction_evaluations") =
      static_cast<double>(main_zero_evaluations),
    Rcpp::Named("main_ordered_zero_sign_uses") =
      static_cast<double>(2 * main_zero_evaluations),
    Rcpp::Named("trace1_zero_direction_evaluations") =
      static_cast<double>(trace1.zero_direction_evaluations),
    Rcpp::Named("trace2_zero_direction_evaluations") =
      static_cast<double>(trace2.zero_direction_evaluations),
    Rcpp::Named("trace3_zero_direction_evaluations") =
      static_cast<double>(trace3_zero_evaluations),
    Rcpp::Named("preprocessed_x") = prepared.x,
    Rcpp::Named("preprocessed_y") = prepared.y,
    Rcpp::Named("preprocessing_anchor") = prepared.anchor,
    Rcpp::Named("preprocessing_log_coordinate_scale") =
      prepared.log_coordinate_scale,
    Rcpp::Named("preprocessing_scaled_coordinate_range") =
      prepared.scaled_coordinate_range,
    Rcpp::Named("preprocessing_global_operand_scale") =
      prepared.global_operand_scale,
    Rcpp::Named("sample_mean1") = prepared.mean_x,
    Rcpp::Named("sample_mean2") = prepared.mean_y,
    Rcpp::Named("scale_identification") = scale_identification,
    Rcpp::Named("tolerance") = tolerance,
    Rcpp::Named("max_iterations") = max_iterations,
    Rcpp::Named("n1") = static_cast<double>(n1),
    Rcpp::Named("n2") = static_cast<double>(n2),
    Rcpp::Named("p") = static_cast<double>(p),
    Rcpp::Named("pair_count1") = static_cast<double>(pairs1.size()),
    Rcpp::Named("pair_count2") = static_cast<double>(pairs2.size()),
    Rcpp::Named("quad_count1") = static_cast<double>(quads1.size()),
    Rcpp::Named("quad_count2") = static_cast<double>(quads2.size())
  );
}
