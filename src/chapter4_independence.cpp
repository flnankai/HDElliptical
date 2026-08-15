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

void ch4ind_require_finite(const arma::mat& x, const char* name) {
  if (!x.is_finite()) {
    Rcpp::stop("`%s` must contain only finite values.", name);
  }
}

long double ch4ind_column_dot(const arma::mat& x,
                              const arma::uword first,
                              const arma::uword second) {
  long double value = 0.0L;
  for (arma::uword row = 0; row < x.n_rows; ++row) {
    value += static_cast<long double>(x(row, first)) *
      static_cast<long double>(x(row, second));
  }
  return value;
}

std::vector<long double> ch4ind_column_norms(const arma::mat& x) {
  std::vector<long double> norms(x.n_cols, 0.0L);
  for (arma::uword column = 0; column < x.n_cols; ++column) {
    long double squared = 0.0L;
    for (arma::uword row = 0; row < x.n_rows; ++row) {
      const long double value = x(row, column);
      squared += value * value;
    }
    if (!std::isfinite(static_cast<double>(squared)) || squared <= 0.0L) {
      Rcpp::stop("Every residual column must have a finite, positive norm.");
    }
    norms[column] = std::sqrt(squared);
  }
  return norms;
}

std::vector<int> ch4ind_ranks_no_ties(const arma::mat& x,
                                      const arma::uword column,
                                      const char* name) {
  std::vector<std::pair<double, int>> ordered(x.n_rows);
  for (arma::uword row = 0; row < x.n_rows; ++row) {
    ordered[row] = std::make_pair(x(row, column), static_cast<int>(row));
  }
  std::sort(
    ordered.begin(), ordered.end(),
    [](const std::pair<double, int>& left,
       const std::pair<double, int>& right) {
      if (left.first < right.first) return true;
      if (left.first > right.first) return false;
      return left.second < right.second;
    }
  );
  for (arma::uword row = 1; row < x.n_rows; ++row) {
    if (ordered[row - 1].first == ordered[row].first) {
      Rcpp::stop(
        "The rank calibration assumes continuous margins; `%s` contains "
        "an exact tie in column %d.", name, static_cast<int>(column + 1)
      );
    }
  }
  std::vector<int> ranks(x.n_rows);
  for (arma::uword row = 0; row < x.n_rows; ++row) {
    ranks[ordered[row].second] = static_cast<int>(row + 1);
  }
  return ranks;
}

arma::imat ch4ind_rank_matrix(const arma::mat& x, const char* name) {
  arma::imat ranks(x.n_rows, x.n_cols);
  for (arma::uword column = 0; column < x.n_cols; ++column) {
    const std::vector<int> current = ch4ind_ranks_no_ties(x, column, name);
    for (arma::uword row = 0; row < x.n_rows; ++row) {
      ranks(row, column) = current[row];
    }
  }
  return ranks;
}

class ch4ind_fenwick {
 public:
  explicit ch4ind_fenwick(const int size) : tree_(size + 1, 0) {}

  void add(int index) {
    const int size = static_cast<int>(tree_.size());
    for (; index < size; index += index & -index) ++tree_[index];
  }

  long long prefix_sum(int index) const {
    long long result = 0;
    for (; index > 0; index -= index & -index) result += tree_[index];
    return result;
  }

 private:
  std::vector<int> tree_;
};

arma::mat ch4ind_rank_correlations(const arma::imat& rank_x,
                                   const arma::imat& rank_y,
                                   const int measure) {
  const arma::uword n = rank_x.n_rows;
  const arma::uword p = rank_x.n_cols;
  const arma::uword q = rank_y.n_cols;
  arma::mat answer(p, q, arma::fill::zeros);

  if (measure == 0) {
    const long double midpoint =
      (static_cast<long double>(n) + 1.0L) / 2.0L;
    const long double n_ld = static_cast<long double>(n);
    const long double coefficient = 12.0L /
      (n_ld * (n_ld * n_ld - 1.0L));
    for (arma::uword i = 0; i < p; ++i) {
      for (arma::uword j = 0; j < q; ++j) {
        long double value = 0.0L;
        for (arma::uword row = 0; row < n; ++row) {
          value +=
            (static_cast<long double>(rank_x(row, i)) - midpoint) *
            (static_cast<long double>(rank_y(row, j)) - midpoint);
        }
        answer(i, j) = static_cast<double>(coefficient * value);
      }
    }
    return answer;
  }

  if (measure != 1) Rcpp::stop("Unknown rank measure code.");
  const long long pairs = static_cast<long long>(n) *
    static_cast<long long>(n - 1) / 2;
  std::vector<int> order(n);
  for (arma::uword i = 0; i < p; ++i) {
    for (arma::uword row = 0; row < n; ++row) {
      order[rank_x(row, i) - 1] = static_cast<int>(row);
    }
    for (arma::uword j = 0; j < q; ++j) {
      ch4ind_fenwick tree(static_cast<int>(n));
      long long inversions = 0;
      for (arma::uword position = 0; position < n; ++position) {
        const int y_rank = rank_y(order[position], j);
        inversions += static_cast<long long>(position) -
          tree.prefix_sum(y_rank);
        tree.add(y_rank);
      }
      answer(i, j) = 1.0 -
        2.0 * static_cast<double>(inversions) / static_cast<double>(pairs);
    }
  }
  return answer;
}

Rcpp::List ch4ind_rank_summary(const arma::mat& correlations,
                               const double null_second_moment) {
  double maximum = -1.0;
  double signed_at_maximum = NA_REAL;
  int maximum_x = 1;
  int maximum_y = 1;
  long double sum_squares = 0.0L;
  for (arma::uword i = 0; i < correlations.n_rows; ++i) {
    for (arma::uword j = 0; j < correlations.n_cols; ++j) {
      const double current = correlations(i, j);
      const double absolute = std::abs(current);
      sum_squares += static_cast<long double>(current) * current;
      if (absolute > maximum) {
        maximum = absolute;
        signed_at_maximum = current;
        maximum_x = static_cast<int>(i + 1);
        maximum_y = static_cast<int>(j + 1);
      }
    }
  }
  const long double comparisons =
    static_cast<long double>(correlations.n_rows) * correlations.n_cols;
  return Rcpp::List::create(
    Rcpp::Named("maximum") = maximum,
    Rcpp::Named("signed_at_maximum") = signed_at_maximum,
    Rcpp::Named("maximum_index") = Rcpp::IntegerVector::create(
      maximum_x, maximum_y
    ),
    Rcpp::Named("sum_squares") = static_cast<double>(sum_squares),
    Rcpp::Named("sum_statistic") = static_cast<double>(
      sum_squares - comparisons * null_second_moment
    )
  );
}

double ch4ind_logdet_sympd(const arma::mat& x, const char* name) {
  arma::mat upper;
  if (!arma::chol(upper, x)) {
    Rcpp::stop("The centered `%s` cross-product is not positive definite.",
               name);
  }
  long double answer = 0.0L;
  for (arma::uword index = 0; index < upper.n_rows; ++index) {
    const double diagonal = upper(index, index);
    if (!std::isfinite(diagonal) || diagonal <= 0.0) {
      Rcpp::stop("The centered `%s` cross-product is numerically singular.",
                 name);
    }
    answer += 2.0L * std::log(static_cast<long double>(diagonal));
  }
  return static_cast<double>(answer);
}

arma::mat ch4ind_inverse_sqrt(const arma::mat& x, const char* name) {
  arma::vec values;
  arma::mat vectors;
  if (!arma::eig_sym(values, vectors, x)) {
    Rcpp::stop("The eigendecomposition of `%s` failed.", name);
  }
  const double largest = values.max();
  const double tolerance = 64.0 * std::numeric_limits<double>::epsilon() *
    std::max(1.0, largest);
  if (!values.is_finite() || values.min() <= tolerance) {
    Rcpp::stop("The centered `%s` cross-product is numerically singular.",
               name);
  }
  return vectors * arma::diagmat(1.0 / arma::sqrt(values)) * vectors.t();
}

}  // namespace


// [[Rcpp::export]]
Rcpp::List cpp_ch4ind_wilks_core(const arma::mat& x,
                                 const arma::mat& y) {
  ch4ind_require_finite(x, "x");
  ch4ind_require_finite(y, "y");
  if (x.n_rows != y.n_rows || x.n_rows < 3 ||
      x.n_cols < 1 || y.n_cols < 1) {
    Rcpp::stop("`x` and `y` must have the same number of rows and at least "
               "one column each.");
  }
  if (x.n_rows <= x.n_cols + y.n_cols) {
    Rcpp::stop("Classical Wilks calibration requires n > p + q.");
  }

  arma::mat centered_x = x.each_row() - arma::mean(x, 0);
  arma::mat centered_y = y.each_row() - arma::mean(y, 0);
  const double divisor = static_cast<double>(x.n_rows - 1);
  arma::mat s_xx = centered_x.t() * centered_x / divisor;
  arma::mat s_yy = centered_y.t() * centered_y / divisor;
  arma::mat s_xy = centered_x.t() * centered_y / divisor;
  arma::mat joint = arma::join_cols(
    arma::join_rows(s_xx, s_xy),
    arma::join_rows(s_xy.t(), s_yy)
  );

  const double logdet_x = ch4ind_logdet_sympd(s_xx, "x");
  const double logdet_y = ch4ind_logdet_sympd(s_yy, "y");
  const double logdet_joint = ch4ind_logdet_sympd(joint, "joint");
  const double log_lambda = logdet_joint - logdet_x - logdet_y;
  if (!std::isfinite(log_lambda) || log_lambda > 1e-8) {
    Rcpp::stop("Wilks' log Lambda is outside its numerical domain.");
  }

  const arma::mat standardized =
    ch4ind_inverse_sqrt(s_xx, "x") * s_xy *
    ch4ind_inverse_sqrt(s_yy, "y");
  arma::vec canonical;
  if (!arma::svd(canonical, standardized)) {
    Rcpp::stop("The canonical-correlation singular value decomposition failed.");
  }
  for (arma::uword index = 0; index < canonical.n_elem; ++index) {
    if (!std::isfinite(canonical[index]) || canonical[index] > 1.0 + 1e-8) {
      Rcpp::stop("A canonical correlation is outside [0, 1].");
    }
    if (canonical[index] > 1.0) canonical[index] = 1.0;
    if (canonical[index] < 0.0) canonical[index] = 0.0;
  }

  return Rcpp::List::create(
    Rcpp::Named("lambda") = std::exp(std::min(0.0, log_lambda)),
    Rcpp::Named("log_lambda") = std::min(0.0, log_lambda),
    Rcpp::Named("canonical_correlations") = canonical,
    Rcpp::Named("covariance_x") = s_xx,
    Rcpp::Named("covariance_y") = s_yy,
    Rcpp::Named("cross_covariance") = s_xy
  );
}


// [[Rcpp::export]]
Rcpp::List cpp_ch4ind_pairwise_core(const arma::mat& residuals,
                                    const bool keep_correlations) {
  ch4ind_require_finite(residuals, "residuals");
  if (residuals.n_rows < 2 || residuals.n_cols < 2) {
    Rcpp::stop("`residuals` must have at least two rows and two columns.");
  }
  const std::vector<long double> norms = ch4ind_column_norms(residuals);
  const arma::uword units = residuals.n_cols;
  arma::mat correlations;
  if (keep_correlations) {
    correlations.eye(units, units);
  }

  long double sum = 0.0L;
  long double sum_squares = 0.0L;
  double maximum = -1.0;
  double signed_at_maximum = NA_REAL;
  int maximum_first = 1;
  int maximum_second = 2;
  for (arma::uword first = 0; first + 1 < units; ++first) {
    for (arma::uword second = first + 1; second < units; ++second) {
      const long double correlation =
        ch4ind_column_dot(residuals, first, second) /
        (norms[first] * norms[second]);
      const double current = static_cast<double>(correlation);
      if (!std::isfinite(current) || std::abs(current) > 1.0 + 1e-10) {
        Rcpp::stop("A sample residual correlation is outside [-1, 1].");
      }
      const double bounded = std::max(-1.0, std::min(1.0, current));
      sum += bounded;
      sum_squares += static_cast<long double>(bounded) * bounded;
      if (keep_correlations) {
        correlations(first, second) = bounded;
        correlations(second, first) = bounded;
      }
      if (std::abs(bounded) > maximum) {
        maximum = std::abs(bounded);
        signed_at_maximum = bounded;
        maximum_first = static_cast<int>(first + 1);
        maximum_second = static_cast<int>(second + 1);
      }
    }
  }

  return Rcpp::List::create(
    Rcpp::Named("sum") = static_cast<double>(sum),
    Rcpp::Named("sum_squares") = static_cast<double>(sum_squares),
    Rcpp::Named("maximum") = maximum,
    Rcpp::Named("signed_at_maximum") = signed_at_maximum,
    Rcpp::Named("maximum_pair") = Rcpp::IntegerVector::create(
      maximum_first, maximum_second
    ),
    Rcpp::Named("correlations") = keep_correlations ?
      Rcpp::wrap(correlations) : R_NilValue
  );
}


// [[Rcpp::export]]
Rcpp::List cpp_ch4ind_serial_panel_core(
    const arma::mat& residuals,
    const double nu,
    const bool estimate_temporal,
    const arma::mat& temporal_covariance,
    const bool keep_matrices) {
  ch4ind_require_finite(residuals, "residuals");
  if (residuals.n_rows < 2 || residuals.n_cols < 3) {
    Rcpp::stop("Serial-panel calibration requires at least two time points "
               "and three panel units.");
  }
  if (!std::isfinite(nu) || nu <= std::sqrt(2.0)) {
    Rcpp::stop("`nu` must be strictly greater than sqrt(2).");
  }
  const arma::uword time = residuals.n_rows;
  const arma::uword units = residuals.n_cols;
  const std::vector<long double> norms = ch4ind_column_norms(residuals);
  arma::mat directions(time, units);
  for (arma::uword unit = 0; unit < units; ++unit) {
    for (arma::uword t = 0; t < time; ++t) {
      directions(t, unit) = static_cast<double>(
        static_cast<long double>(residuals(t, unit)) / norms[unit]
      );
    }
  }
  arma::mat correlations = directions.t() * directions;
  correlations.diag().ones();

  long double correlation_sum = 0.0L;
  double maximum_square = -1.0;
  double signed_at_maximum = NA_REAL;
  int maximum_first = 1;
  int maximum_second = 2;
  for (arma::uword first = 0; first + 1 < units; ++first) {
    for (arma::uword second = first + 1; second < units; ++second) {
      const double current = correlations(first, second);
      if (!std::isfinite(current) || std::abs(current) > 1.0 + 1e-10) {
        Rcpp::stop("A sample residual correlation is outside [-1, 1].");
      }
      correlation_sum += current;
      const double square = current * current;
      if (square > maximum_square) {
        maximum_square = square;
        signed_at_maximum = current;
        maximum_first = static_cast<int>(first + 1);
        maximum_second = static_cast<int>(second + 1);
      }
    }
  }

  const long double units_ld = static_cast<long double>(units);
  const double sum_statistic = static_cast<double>(
    std::sqrt(2.0L / (units_ld * (units_ld - 1.0L))) * correlation_sum
  );
  const arma::vec row_sums = arma::sum(correlations, 1);
  long double variance_sum = 0.0L;
  for (arma::uword first = 0; first + 1 < units; ++first) {
    for (arma::uword second = first + 1; second < units; ++second) {
      const long double current = correlations(first, second);
      const long double first_bar =
        (static_cast<long double>(row_sums[first]) - 1.0L - current) /
        static_cast<long double>(units - 2);
      const long double second_bar =
        (static_cast<long double>(row_sums[second]) - 1.0L - current) /
        static_cast<long double>(units - 2);
      variance_sum += (current - second_bar) * (current - first_bar);
    }
  }
  const double sum_variance = static_cast<double>(
    2.0L * variance_sum / (units_ld * (units_ld - 1.0L))
  );
  if (!std::isfinite(sum_variance) || sum_variance <= 0.0) {
    Rcpp::stop("The published serial-panel sum-variance estimate is not "
               "strictly positive for these residuals.");
  }

  arma::mat sigma_hat;
  arma::mat sigma_tilde;
  arma::mat u_hat;
  double p_hat = NA_REAL;
  double threshold = NA_REAL;
  double effective_dimension = NA_REAL;

  if (estimate_temporal) {
    const arma::vec time_means = arma::mean(residuals, 1);
    const arma::mat centered = residuals.each_col() - time_means;
    sigma_hat = centered * centered.t() /
      static_cast<double>(units - 1);
    const arma::vec diagonal = sigma_hat.diag();
    if (!diagonal.is_finite() || diagonal.min() <= 0.0) {
      Rcpp::stop("Every time coordinate needs positive cross-unit sample "
                 "variance for the published temporal estimator.");
    }
    const double trace_hat = arma::trace(sigma_hat);
    if (!std::isfinite(trace_hat) || trace_hat <= 0.0) {
      Rcpp::stop("The temporal sample covariance has non-positive trace.");
    }
    u_hat = residuals.t() * residuals / trace_hat;
    p_hat = (arma::accu(arma::square(u_hat)) -
      std::pow(arma::trace(u_hat), 2.0) / static_cast<double>(time)) /
      static_cast<double>(units);
    if (!std::isfinite(p_hat) || p_hat <= 0.0) {
      Rcpp::stop("The literal published threshold factor P-hat_N is not "
                 "positive; supply `temporal_covariance` rather than "
                 "replacing or projecting this estimator.");
    }
    threshold = nu * std::sqrt(
      p_hat * std::log(static_cast<double>(time)) /
      static_cast<double>(units)
    );
    sigma_tilde.zeros(time, time);
    for (arma::uword first = 0; first < time; ++first) {
      sigma_tilde(first, first) = sigma_hat(first, first);
      for (arma::uword second = first + 1; second < time; ++second) {
        const double theta = sigma_hat(first, second) /
          std::sqrt(sigma_hat(first, first) * sigma_hat(second, second));
        const double denominator = 1.0 - theta * theta;
        const double score = denominator <= 0.0 ?
          std::numeric_limits<double>::infinity() :
          std::abs(theta) / denominator;
        const double kept = score >= threshold ?
          sigma_hat(first, second) : 0.0;
        sigma_tilde(first, second) = kept;
        sigma_tilde(second, first) = kept;
      }
    }
    const double trace_tilde = arma::trace(sigma_tilde);
    const double frobenius_square = arma::accu(arma::square(sigma_tilde));
    if (!std::isfinite(trace_tilde) || !std::isfinite(frobenius_square) ||
        trace_tilde <= 0.0 || frobenius_square <= 0.0) {
      Rcpp::stop("The thresholded temporal covariance has an invalid "
                 "trace/Frobenius ratio.");
    }
    effective_dimension = trace_tilde * trace_tilde / frobenius_square;
  } else {
    ch4ind_require_finite(temporal_covariance, "temporal_covariance");
    if (temporal_covariance.n_rows != time ||
        temporal_covariance.n_cols != time ||
        !arma::approx_equal(temporal_covariance,
                            temporal_covariance.t(),
                            "absdiff", 1e-10)) {
      Rcpp::stop("`temporal_covariance` must be a symmetric T by T matrix.");
    }
    arma::mat upper;
    if (!arma::chol(upper, temporal_covariance)) {
      Rcpp::stop("`temporal_covariance` must be positive definite.");
    }
    const double trace_temporal = arma::trace(temporal_covariance);
    const double frobenius_square =
      arma::accu(arma::square(temporal_covariance));
    effective_dimension = trace_temporal * trace_temporal /
      frobenius_square;
    sigma_tilde = temporal_covariance;
  }

  if (!std::isfinite(effective_dimension) || effective_dimension <= 0.0) {
    Rcpp::stop("The temporal effective dimension is not positive and finite.");
  }

  return Rcpp::List::create(
    Rcpp::Named("sum_statistic") = sum_statistic,
    Rcpp::Named("sum_variance") = sum_variance,
    Rcpp::Named("maximum_square") = maximum_square,
    Rcpp::Named("signed_at_maximum") = signed_at_maximum,
    Rcpp::Named("maximum_pair") = Rcpp::IntegerVector::create(
      maximum_first, maximum_second
    ),
    Rcpp::Named("effective_dimension") = effective_dimension,
    Rcpp::Named("p_hat") = p_hat,
    Rcpp::Named("threshold") = threshold,
    Rcpp::Named("correlations") = keep_matrices ?
      Rcpp::wrap(correlations) : R_NilValue,
    Rcpp::Named("sigma_hat") = keep_matrices && estimate_temporal ?
      Rcpp::wrap(sigma_hat) : R_NilValue,
    Rcpp::Named("sigma_tilde") = keep_matrices ?
      Rcpp::wrap(sigma_tilde) : R_NilValue,
    Rcpp::Named("u_hat") = keep_matrices && estimate_temporal ?
      Rcpp::wrap(u_hat) : R_NilValue
  );
}


// [[Rcpp::export]]
Rcpp::List cpp_ch4ind_rank_vector_core(
    const arma::mat& x,
    const arma::mat& y,
    const int measure,
    const arma::imat& permutations,
    const bool keep_correlations,
    const bool keep_permutation) {
  ch4ind_require_finite(x, "x");
  ch4ind_require_finite(y, "y");
  if (x.n_rows != y.n_rows || x.n_rows < 3 ||
      x.n_cols < 1 || y.n_cols < 1) {
    Rcpp::stop("`x` and `y` must have the same number of rows, at least "
               "three rows, and at least one column each.");
  }
  if (measure != 0 && measure != 1) {
    Rcpp::stop("Unknown rank measure code.");
  }
  if (permutations.n_rows != x.n_rows || permutations.n_cols < 2) {
    Rcpp::stop("`permutations` must have n rows and at least two columns.");
  }

  const arma::uword n = x.n_rows;
  const arma::imat rank_x = ch4ind_rank_matrix(x, "x");
  const arma::imat rank_y = ch4ind_rank_matrix(y, "y");
  const arma::mat observed_correlations =
    ch4ind_rank_correlations(rank_x, rank_y, measure);
  const double n_double = static_cast<double>(n);
  const double null_second_moment = measure == 0 ?
    1.0 / (n_double - 1.0) :
    2.0 * (2.0 * n_double + 5.0) /
      (9.0 * n_double * (n_double - 1.0));
  const Rcpp::List observed = ch4ind_rank_summary(
    observed_correlations, null_second_moment
  );

  arma::vec permutation_statistics(permutations.n_cols);
  std::vector<int> seen(n);
  arma::imat permuted_rank_x(n, rank_x.n_cols);
  for (arma::uword draw = 0; draw < permutations.n_cols; ++draw) {
    std::fill(seen.begin(), seen.end(), 0);
    for (arma::uword row = 0; row < n; ++row) {
      const int source = permutations(row, draw);
      if (source < 1 || source > static_cast<int>(n) ||
          seen[source - 1] != 0) {
        Rcpp::stop("Every column of `permutations` must be a permutation "
                   "of 1:n.");
      }
      seen[source - 1] = 1;
      permuted_rank_x.row(row) = rank_x.row(source - 1);
    }
    const arma::mat current = ch4ind_rank_correlations(
      permuted_rank_x, rank_y, measure
    );
    const Rcpp::List summary = ch4ind_rank_summary(
      current, null_second_moment
    );
    permutation_statistics[draw] = Rcpp::as<double>(
      summary["sum_statistic"]
    );
  }
  const double permutation_mean = arma::mean(permutation_statistics);
  const double permutation_variance = arma::accu(
    arma::square(permutation_statistics - permutation_mean)
  ) / static_cast<double>(permutations.n_cols - 1);
  if (!std::isfinite(permutation_variance) ||
      permutation_variance <= 0.0) {
    Rcpp::stop("The intrinsic permutation variance estimate is not "
               "strictly positive; increase `B` or use non-degenerate data.");
  }

  return Rcpp::List::create(
    Rcpp::Named("maximum") = observed["maximum"],
    Rcpp::Named("signed_at_maximum") = observed["signed_at_maximum"],
    Rcpp::Named("maximum_index") = observed["maximum_index"],
    Rcpp::Named("sum_squares") = observed["sum_squares"],
    Rcpp::Named("sum_statistic") = observed["sum_statistic"],
    Rcpp::Named("null_second_moment") = null_second_moment,
    Rcpp::Named("permutation_variance") = permutation_variance,
    Rcpp::Named("permutation_mean") = permutation_mean,
    Rcpp::Named("correlations") = keep_correlations ?
      Rcpp::wrap(observed_correlations) : R_NilValue,
    Rcpp::Named("permutation_statistics") = keep_permutation ?
      Rcpp::wrap(permutation_statistics) : R_NilValue
  );
}
