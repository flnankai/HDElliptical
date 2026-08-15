// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp17)]]

#include <RcppArmadillo.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <utility>
#include <vector>

namespace {

void ch4wn_require_finite_matrix(const arma::mat& x, const char* name) {
  if (!x.is_finite()) {
    Rcpp::stop("`%s` must contain only finite values.", name);
  }
}

long double ch4wn_row_dot(const arma::mat& x,
                          const arma::uword first,
                          const arma::uword second) {
  long double value = 0.0L;
  for (arma::uword j = 0; j < x.n_cols; ++j) {
    value += static_cast<long double>(x(first, j)) *
      static_cast<long double>(x(second, j));
  }
  return value;
}

long double ch4wn_lagged_column_dot(const arma::mat& x,
                                    const arma::uword lag,
                                    const arma::uword first,
                                    const arma::uword second) {
  const arma::uword effective_n = x.n_rows - lag;
  long double value = 0.0L;
  for (arma::uword t = 0; t < effective_n; ++t) {
    value += static_cast<long double>(x(t, first)) *
      static_cast<long double>(x(t + lag, second));
  }
  return value;
}

bool ch4wn_unit_direction(const arma::rowvec& value,
                          arma::rowvec& direction) {
  const double scale = arma::abs(value).max();
  if (scale == 0.0) {
    direction.zeros(value.n_elem);
    return false;
  }
  const arma::rowvec scaled = value / scale;
  const double norm = arma::norm(scaled, 2);
  if (!std::isfinite(norm) || norm <= 0.0) {
    Rcpp::stop("A spatial-sign norm is non-finite; no repair is applied.");
  }
  direction = scaled / norm;
  return true;
}

std::vector<int> ch4wn_ranks_no_ties(const arma::mat& x,
                                     const arma::uword row_start,
                                     const arma::uword length,
                                     const arma::uword column) {
  std::vector<std::pair<double, int>> ordered(length);
  for (arma::uword t = 0; t < length; ++t) {
    ordered[t] = std::make_pair(x(row_start + t, column),
                                static_cast<int>(t));
  }
  std::sort(
    ordered.begin(), ordered.end(),
    [](const std::pair<double, int>& left,
       const std::pair<double, int>& right) {
      return left.first < right.first;
    }
  );
  for (arma::uword t = 1; t < length; ++t) {
    if (ordered[t - 1].first == ordered[t].first) {
      Rcpp::stop(
        "The rank white-noise calibration requires continuous margins "
        "without ties; an exact tie was found."
      );
    }
  }
  std::vector<int> ranks(length);
  for (arma::uword t = 0; t < length; ++t) {
    ranks[ordered[t].second] = static_cast<int>(t + 1);
  }
  return ranks;
}

class ch4wn_fenwick {
 public:
  explicit ch4wn_fenwick(const int size) : tree_(size + 1, 0) {}

  void add(int index) {
    const int size = static_cast<int>(tree_.size());
    for (; index < size; index += index & -index) {
      ++tree_[index];
    }
  }

  long long prefix_sum(int index) const {
    long long result = 0;
    for (; index > 0; index -= index & -index) {
      result += tree_[index];
    }
    return result;
  }

 private:
  std::vector<int> tree_;
};

} // namespace


// [[Rcpp::export]]
Rcpp::List cpp_ch4wn_flm_core(const arma::mat& residuals,
                              const int lag,
                              const bool keep_lag) {
  ch4wn_require_finite_matrix(residuals, "residuals");
  if (residuals.n_rows < 3 || residuals.n_cols < 1) {
    Rcpp::stop("`residuals` must have at least three rows and one column.");
  }
  if (lag < 1 || static_cast<arma::uword>(lag) > residuals.n_rows - 2) {
    Rcpp::stop("`lag` must be between 1 and n - 2.");
  }

  const arma::uword n = residuals.n_rows;
  const arma::uword p = residuals.n_cols;
  const long double n_ld = static_cast<long double>(n);
  const long double ordered_denominator = n_ld *
    static_cast<long double>(n - 1);

  arma::vec marginal_variances(p, arma::fill::zeros);
  for (arma::uword j = 0; j < p; ++j) {
    long double sum = 0.0L;
    for (arma::uword t = 0; t < n; ++t) {
      const long double value = residuals(t, j);
      sum += value * value;
    }
    marginal_variances(j) = static_cast<double>(sum / n_ld);
    if (!std::isfinite(marginal_variances(j)) ||
        marginal_variances(j) <= 0.0) {
      Rcpp::stop(
        "Every coordinate must have a finite, strictly positive "
        "primary n-divisor second moment."
      );
    }
  }

  double maximum_correlation = -1.0;
  double signed_correlation_at_maximum = NA_REAL;
  int maximum_lag = 1;
  int maximum_row = 1;
  int maximum_column = 1;
  arma::vec lag_maximum_correlation(lag, arma::fill::zeros);
  arma::imat lag_maximum_indices(lag, 2, arma::fill::ones);

  for (int k = 1; k <= lag; ++k) {
    double lag_maximum = -1.0;
    for (arma::uword i = 0; i < p; ++i) {
      for (arma::uword j = 0; j < p; ++j) {
        const long double covariance = ch4wn_lagged_column_dot(
          residuals, static_cast<arma::uword>(k), i, j
        ) / n_ld;
        const long double denominator = std::sqrt(
          static_cast<long double>(marginal_variances(i)) *
          static_cast<long double>(marginal_variances(j))
        );
        const double correlation = static_cast<double>(
          covariance / denominator
        );
        if (!std::isfinite(correlation)) {
          Rcpp::stop("A lagged sample correlation is non-finite.");
        }
        const double absolute = std::abs(correlation);
        if (absolute > lag_maximum) {
          lag_maximum = absolute;
          lag_maximum_correlation(k - 1) = absolute;
          lag_maximum_indices(k - 1, 0) = static_cast<int>(i + 1);
          lag_maximum_indices(k - 1, 1) = static_cast<int>(j + 1);
        }
        if (absolute > maximum_correlation) {
          maximum_correlation = absolute;
          signed_correlation_at_maximum = correlation;
          maximum_lag = k;
          maximum_row = static_cast<int>(i + 1);
          maximum_column = static_cast<int>(j + 1);
        }
      }
    }
  }

  arma::mat gram(n, n, arma::fill::zeros);
  for (arma::uword t = 0; t < n; ++t) {
    for (arma::uword s = t; s < n; ++s) {
      const double value = static_cast<double>(
        ch4wn_row_dot(residuals, t, s)
      );
      if (!std::isfinite(value)) {
        Rcpp::stop("A row inner product is non-finite.");
      }
      gram(t, s) = value;
      gram(s, t) = value;
    }
  }

  long double trace_numerator = 0.0L;
  for (arma::uword t = 0; t < n; ++t) {
    for (arma::uword s = 0; s < n; ++s) {
      if (s == t) continue;
      const long double value = gram(t, s);
      trace_numerator += value * value;
    }
  }

  long double sum_numerator = 0.0L;
  arma::vec lag_sum_numerator(lag, arma::fill::zeros);
  arma::vec lag_ordered_pair_count(lag, arma::fill::zeros);
  for (int k = 1; k <= lag; ++k) {
    const arma::uword effective_n = n - static_cast<arma::uword>(k);
    long double lag_numerator = 0.0L;
    for (arma::uword t = 0; t < effective_n; ++t) {
      for (arma::uword s = 0; s < effective_n; ++s) {
        if (s == t) continue;
        lag_numerator += static_cast<long double>(gram(t, s)) *
          static_cast<long double>(gram(t + k, s + k));
      }
    }
    sum_numerator += lag_numerator;
    lag_sum_numerator(k - 1) = static_cast<double>(lag_numerator);
    lag_ordered_pair_count(k - 1) =
      static_cast<double>(effective_n) *
      static_cast<double>(effective_n - 1);
  }

  const double trace_estimate = static_cast<double>(
    trace_numerator / ordered_denominator
  );
  const double sum_statistic = static_cast<double>(
    sum_numerator / ordered_denominator
  );
  if (!std::isfinite(trace_estimate) || trace_estimate <= 0.0) {
    Rcpp::stop(
      "The feasible trace estimate is not finite and positive; "
      "no floor or absolute-value repair is applied."
    );
  }
  if (!std::isfinite(sum_statistic)) {
    Rcpp::stop("The FLM sum statistic is non-finite.");
  }

  const double maximum_statistic = std::sqrt(static_cast<double>(n)) *
    maximum_correlation;
  Rcpp::RObject lag_maximum_out = R_NilValue;
  Rcpp::RObject lag_indices_out = R_NilValue;
  Rcpp::RObject lag_sum_out = R_NilValue;
  if (keep_lag) {
    lag_maximum_out = Rcpp::wrap(lag_maximum_correlation);
    lag_indices_out = Rcpp::wrap(lag_maximum_indices);
    lag_sum_out = Rcpp::wrap(lag_sum_numerator);
  }

  return Rcpp::List::create(
    Rcpp::Named("maximum_statistic") = maximum_statistic,
    Rcpp::Named("maximum_absolute_correlation") = maximum_correlation,
    Rcpp::Named("signed_correlation_at_maximum") =
      signed_correlation_at_maximum,
    Rcpp::Named("maximum_lag") = maximum_lag,
    Rcpp::Named("maximum_indices") = Rcpp::IntegerVector::create(
      maximum_row, maximum_column
    ),
    Rcpp::Named("sum_statistic_scaled") = sum_statistic,
    Rcpp::Named("sum_numerator_scaled") =
      static_cast<double>(sum_numerator),
    Rcpp::Named("trace_estimate_scaled") = trace_estimate,
    Rcpp::Named("trace_numerator_scaled") =
      static_cast<double>(trace_numerator),
    Rcpp::Named("primary_ordered_denominator") =
      static_cast<double>(ordered_denominator),
    Rcpp::Named("marginal_variances_scaled") = marginal_variances,
    Rcpp::Named("lag_ordered_pair_count") = lag_ordered_pair_count,
    Rcpp::Named("lag_maximum_absolute_correlation") = lag_maximum_out,
    Rcpp::Named("lag_maximum_indices") = lag_indices_out,
    Rcpp::Named("lag_sum_numerator_scaled") = lag_sum_out
  );
}


// [[Rcpp::export]]
Rcpp::List cpp_ch4wn_spatial_sign_core(const arma::mat& residuals,
                                       const int lag,
                                       const bool keep_signs) {
  ch4wn_require_finite_matrix(residuals, "residuals");
  if (residuals.n_rows < 3 || residuals.n_cols < 1) {
    Rcpp::stop("`residuals` must have at least three rows and one column.");
  }
  if (lag < 1 || static_cast<arma::uword>(lag) > residuals.n_rows - 2) {
    Rcpp::stop("`lag` must be between 1 and n - 2.");
  }

  const arma::uword n = residuals.n_rows;
  const arma::uword p = residuals.n_cols;
  arma::mat signs(n, p, arma::fill::zeros);
  int zero_count = 0;
  for (arma::uword t = 0; t < n; ++t) {
    arma::rowvec direction(p, arma::fill::zeros);
    if (ch4wn_unit_direction(residuals.row(t), direction)) {
      signs.row(t) = direction;
    } else {
      ++zero_count;
    }
  }

  arma::mat gram(n, n, arma::fill::zeros);
  for (arma::uword t = 0; t < n; ++t) {
    for (arma::uword s = t; s < n; ++s) {
      const double value = arma::dot(signs.row(t), signs.row(s));
      gram(t, s) = value;
      gram(s, t) = value;
    }
  }

  long double trace_sum = 0.0L;
  for (arma::uword s = 0; s + 1 < n; ++s) {
    for (arma::uword t = s + 1; t < n; ++t) {
      const long double value = gram(s, t);
      trace_sum += value * value;
    }
  }
  const long double unordered_denominator =
    static_cast<long double>(n) * static_cast<long double>(n - 1);
  const double trace_estimate = static_cast<double>(
    2.0L * trace_sum / unordered_denominator
  );

  long double statistic = 0.0L;
  arma::vec lag_numerators(lag, arma::fill::zeros);
  arma::vec lag_statistics(lag, arma::fill::zeros);
  arma::vec lag_unordered_pair_count(lag, arma::fill::zeros);
  for (int h = 1; h <= lag; ++h) {
    long double numerator = 0.0L;
    const arma::uword first = static_cast<arma::uword>(h);
    for (arma::uword s = first; s + 1 < n; ++s) {
      for (arma::uword t = s + 1; t < n; ++t) {
        numerator += static_cast<long double>(gram(s - h, t - h)) *
          static_cast<long double>(gram(s, t));
      }
    }
    const long double contribution = numerator /
      static_cast<long double>(n - h);
    statistic += contribution;
    lag_numerators(h - 1) = static_cast<double>(numerator);
    lag_statistics(h - 1) = static_cast<double>(contribution);
    const double effective_n = static_cast<double>(n - h);
    lag_unordered_pair_count(h - 1) =
      effective_n * (effective_n - 1.0) / 2.0;
  }

  if (!std::isfinite(trace_estimate) || trace_estimate <= 0.0) {
    Rcpp::stop(
      "The spatial-sign trace estimate is not finite and positive; "
      "no floor or deletion is applied."
    );
  }
  const double statistic_double = static_cast<double>(statistic);
  if (!std::isfinite(statistic_double)) {
    Rcpp::stop("The spatial-sign white-noise statistic is non-finite.");
  }

  Rcpp::RObject signs_out = R_NilValue;
  if (keep_signs) signs_out = Rcpp::wrap(signs);
  return Rcpp::List::create(
    Rcpp::Named("statistic") = statistic_double,
    Rcpp::Named("trace_estimate") = trace_estimate,
    Rcpp::Named("trace_unordered_sum") = static_cast<double>(trace_sum),
    Rcpp::Named("trace_denominator") =
      static_cast<double>(unordered_denominator),
    Rcpp::Named("zero_sign_count") = zero_count,
    Rcpp::Named("lag_numerators") = lag_numerators,
    Rcpp::Named("lag_statistics") = lag_statistics,
    Rcpp::Named("lag_unordered_pair_count") = lag_unordered_pair_count,
    Rcpp::Named("signs") = signs_out
  );
}


// [[Rcpp::export]]
Rcpp::List cpp_ch4wn_rank_max_core(const arma::mat& x,
                                   const int lag,
                                   const std::string& measure,
                                   const bool keep_lag) {
  ch4wn_require_finite_matrix(x, "x");
  if (x.n_rows < 3 || x.n_cols < 1) {
    Rcpp::stop("`x` must have at least three rows and one column.");
  }
  if (lag < 1 || static_cast<arma::uword>(lag) > x.n_rows - 2) {
    Rcpp::stop("`lag` must be between 1 and n - 2.");
  }
  const bool spearman = measure == "spearman";
  const bool kendall = measure == "kendall";
  if (!spearman && !kendall) {
    Rcpp::stop("`measure` must be 'spearman' or 'kendall'.");
  }

  const arma::uword n = x.n_rows;
  const arma::uword p = x.n_cols;
  double maximum_score = -1.0;
  double measure_at_maximum = NA_REAL;
  int maximum_lag = 1;
  int maximum_row = 1;
  int maximum_column = 1;
  arma::vec lag_maximum_score(lag, arma::fill::zeros);
  arma::vec lag_measure_at_maximum(lag, arma::fill::zeros);
  arma::imat lag_maximum_indices(lag, 2, arma::fill::ones);

  for (int k = 1; k <= lag; ++k) {
    const arma::uword m = n - static_cast<arma::uword>(k);
    std::vector<std::vector<int>> first_ranks(p);
    std::vector<std::vector<int>> second_ranks(p);
    std::vector<std::vector<int>> first_order(p,
                                              std::vector<int>(m));
    for (arma::uword j = 0; j < p; ++j) {
      first_ranks[j] = ch4wn_ranks_no_ties(x, 0, m, j);
      second_ranks[j] = ch4wn_ranks_no_ties(
        x, static_cast<arma::uword>(k), m, j
      );
      for (arma::uword t = 0; t < m; ++t) {
        first_order[j][first_ranks[j][t] - 1] = static_cast<int>(t);
      }
    }

    const long double rank_mean =
      (static_cast<long double>(m) + 1.0L) / 2.0L;
    const long double spearman_denominator =
      static_cast<long double>(m) *
      (static_cast<long double>(m) * static_cast<long double>(m) - 1.0L) /
      12.0L;
    const long double kendall_pairs =
      static_cast<long double>(m) * static_cast<long double>(m - 1) / 2.0L;
    double current_lag_maximum = -1.0;

    for (arma::uword i = 0; i < p; ++i) {
      for (arma::uword j = 0; j < p; ++j) {
        double value = 0.0;
        double score = 0.0;
        if (spearman) {
          long double numerator = 0.0L;
          for (arma::uword t = 0; t < m; ++t) {
            numerator +=
              (static_cast<long double>(first_ranks[i][t]) - rank_mean) *
              (static_cast<long double>(second_ranks[j][t]) - rank_mean);
          }
          value = static_cast<double>(numerator / spearman_denominator);
          score = static_cast<double>(m) * value * value;
        } else {
          ch4wn_fenwick tree(static_cast<int>(m));
          long long inversions = 0;
          long long seen = 0;
          for (arma::uword rank = 0; rank < m; ++rank) {
            const int row = first_order[i][rank];
            const int second_rank = second_ranks[j][row];
            inversions += seen - tree.prefix_sum(second_rank);
            tree.add(second_rank);
            ++seen;
          }
          const long double concordance_difference = kendall_pairs -
            2.0L * static_cast<long double>(inversions);
          value = static_cast<double>(
            concordance_difference / kendall_pairs
          );
          const long double coefficient =
            9.0L * static_cast<long double>(m) *
            static_cast<long double>(m - 1) /
            (2.0L * (2.0L * static_cast<long double>(m) + 5.0L));
          score = static_cast<double>(coefficient) * value * value;
        }
        if (!std::isfinite(value) || !std::isfinite(score)) {
          Rcpp::stop("A rank statistic is non-finite.");
        }
        if (score > current_lag_maximum) {
          current_lag_maximum = score;
          lag_maximum_score(k - 1) = score;
          lag_measure_at_maximum(k - 1) = value;
          lag_maximum_indices(k - 1, 0) = static_cast<int>(i + 1);
          lag_maximum_indices(k - 1, 1) = static_cast<int>(j + 1);
        }
        if (score > maximum_score) {
          maximum_score = score;
          measure_at_maximum = value;
          maximum_lag = k;
          maximum_row = static_cast<int>(i + 1);
          maximum_column = static_cast<int>(j + 1);
        }
      }
    }
  }

  Rcpp::RObject lag_score_out = R_NilValue;
  Rcpp::RObject lag_measure_out = R_NilValue;
  Rcpp::RObject lag_indices_out = R_NilValue;
  if (keep_lag) {
    lag_score_out = Rcpp::wrap(lag_maximum_score);
    lag_measure_out = Rcpp::wrap(lag_measure_at_maximum);
    lag_indices_out = Rcpp::wrap(lag_maximum_indices);
  }
  return Rcpp::List::create(
    Rcpp::Named("maximum_standardized_square") = maximum_score,
    Rcpp::Named("measure_at_maximum") = measure_at_maximum,
    Rcpp::Named("maximum_lag") = maximum_lag,
    Rcpp::Named("maximum_indices") = Rcpp::IntegerVector::create(
      maximum_row, maximum_column
    ),
    Rcpp::Named("lag_maximum_standardized_square") = lag_score_out,
    Rcpp::Named("lag_measure_at_maximum") = lag_measure_out,
    Rcpp::Named("lag_maximum_indices") = lag_indices_out,
    Rcpp::Named("continuous_margin_calibration") = true,
    Rcpp::Named("ties_detected") = false
  );
}
