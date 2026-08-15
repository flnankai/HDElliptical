// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp11)]]

#include <RcppArmadillo.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <vector>

namespace ch4completion {

inline double stable_norm(const arma::rowvec& x) {
  double scale = 0.0;
  long double sumsq = 1.0L;
  for (arma::uword j = 0; j < x.n_elem; ++j) {
    const double value = std::abs(x[j]);
    if (!std::isfinite(value)) {
      Rcpp::stop("A standardized residual component is not finite.");
    }
    if (value == 0.0) {
      continue;
    }
    if (scale < value) {
      const long double ratio = scale == 0.0 ? 0.0L :
        static_cast<long double>(scale) / static_cast<long double>(value);
      sumsq = 1.0L + sumsq * ratio * ratio;
      scale = value;
    } else {
      const long double ratio = static_cast<long double>(value) /
        static_cast<long double>(scale);
      sumsq += ratio * ratio;
    }
  }
  if (scale == 0.0) {
    return 0.0;
  }
  const long double answer = static_cast<long double>(scale) *
    std::sqrt(sumsq);
  if (!(answer <= static_cast<long double>(
          std::numeric_limits<double>::max()))) {
    Rcpp::stop("A standardized residual radius is not representable.");
  }
  return static_cast<double>(answer);
}

inline long double choose_long_double(const int n, const int m) {
  if (m < 0 || m > n) {
    return 0.0L;
  }
  const int k = std::min(m, n - m);
  long double answer = 1.0L;
  for (int j = 1; j <= k; ++j) {
    answer *= static_cast<long double>(n - k + j) /
      static_cast<long double>(j);
  }
  return answer;
}

inline long double factorial_long_double(const int m) {
  long double answer = 1.0L;
  for (int j = 2; j <= m; ++j) {
    answer *= static_cast<long double>(j);
  }
  return answer;
}

void combination_rec(const int n, const int m, const int start,
                     const int depth, std::vector<int>& current,
                     std::vector<std::vector<int> >& output) {
  if (depth == m) {
    output.push_back(current);
    return;
  }
  const int remaining = m - depth;
  for (int value = start; value <= n - remaining; ++value) {
    current[depth] = value;
    combination_rec(n, m, value + 1, depth + 1, current, output);
  }
}

inline std::vector<std::vector<int> > combinations(const int n,
                                                    const int m) {
  std::vector<std::vector<int> > output;
  const long double count = choose_long_double(n, m);
  if (count <= static_cast<long double>(
          std::numeric_limits<std::size_t>::max())) {
    output.reserve(static_cast<std::size_t>(count));
  }
  std::vector<int> current(m);
  combination_rec(n, m, 0, 0, current, output);
  return output;
}

inline std::vector<std::vector<int> > orderings(const int m) {
  std::vector<int> current(m);
  std::iota(current.begin(), current.end(), 0);
  std::vector<std::vector<int> > output;
  output.reserve(static_cast<std::size_t>(factorial_long_double(m)));
  do {
    output.push_back(current);
  } while (std::next_permutation(current.begin(), current.end()));
  return output;
}

inline int indicator_leq(const double left, const double right) {
  return left <= right ? 1 : 0;
}

inline int pair_contrast(const std::vector<double>& values,
                         const std::vector<int>& order,
                         const int anchor) {
  const int first = indicator_leq(values[order[0]], values[order[anchor]]) -
    indicator_leq(values[order[1]], values[order[anchor]]);
  const int second = indicator_leq(values[order[2]], values[order[anchor]]) -
    indicator_leq(values[order[3]], values[order[anchor]]);
  return first * second;
}

inline int separated_pair(const std::vector<double>& values,
                          const int a, const int b,
                          const int c, const int d) {
  return values[a] < values[c] && values[a] < values[d] &&
    values[b] < values[c] && values[b] < values[d] ? 1 : 0;
}

inline int tau_contrast(const std::vector<double>& values,
                        const std::vector<int>& order) {
  return
    separated_pair(values, order[0], order[2], order[1], order[3]) +
    separated_pair(values, order[1], order[3], order[0], order[2]) -
    separated_pair(values, order[0], order[3], order[1], order[2]) -
    separated_pair(values, order[1], order[2], order[0], order[3]);
}

inline double symmetrized_kernel(
    const std::vector<double>& x,
    const std::vector<double>& y,
    const int measure,
    const std::vector<std::vector<int> >& permutations) {
  long double total = 0.0L;
  if (measure == 0) {
    for (const std::vector<int>& order : permutations) {
      total += static_cast<long double>(pair_contrast(x, order, 4)) *
        static_cast<long double>(pair_contrast(y, order, 4));
    }
    return static_cast<double>(total / 16.0L);
  }
  if (measure == 1) {
    for (const std::vector<int>& order : permutations) {
      total += static_cast<long double>(pair_contrast(x, order, 4)) *
        static_cast<long double>(pair_contrast(y, order, 5));
    }
    return static_cast<double>(total / 32.0L);
  }
  for (const std::vector<int>& order : permutations) {
    total += static_cast<long double>(tau_contrast(x, order)) *
      static_cast<long double>(tau_contrast(y, order));
  }
  return static_cast<double>(
    total / static_cast<long double>(permutations.size())
  );
}

inline bool column_has_tie(const arma::mat& x, const arma::uword column) {
  std::vector<double> values(x.n_rows);
  for (arma::uword i = 0; i < x.n_rows; ++i) {
    values[i] = x(i, column);
  }
  std::sort(values.begin(), values.end());
  for (std::size_t i = 1; i < values.size(); ++i) {
    if (values[i] == values[i - 1]) {
      return true;
    }
  }
  return false;
}

inline arma::mat u_statistics(
    const arma::mat& x,
    const arma::mat& y,
    const int measure,
    const std::vector<std::vector<int> >& subsets,
    const std::vector<std::vector<int> >& permutations,
    const std::vector<int>* permuted_x,
    std::size_t& interrupt_counter) {
  const int order = measure == 0 ? 5 : (measure == 1 ? 6 : 4);
  arma::mat answer(x.n_cols, y.n_cols, arma::fill::zeros);
  std::vector<double> x_values(order);
  std::vector<double> y_values(order);
  for (arma::uword i = 0; i < x.n_cols; ++i) {
    for (arma::uword j = 0; j < y.n_cols; ++j) {
      long double total = 0.0L;
      for (const std::vector<int>& subset : subsets) {
        for (int k = 0; k < order; ++k) {
          const int row = subset[k];
          const int x_row = permuted_x == nullptr ? row :
            (*permuted_x)[row];
          x_values[k] = x(x_row, i);
          y_values[k] = y(row, j);
        }
        total += static_cast<long double>(symmetrized_kernel(
          x_values, y_values, measure, permutations
        ));
        ++interrupt_counter;
        if ((interrupt_counter & 4095U) == 0U) {
          Rcpp::checkUserInterrupt();
        }
      }
      answer(i, j) = static_cast<double>(
        total / static_cast<long double>(subsets.size())
      );
    }
  }
  return answer;
}

struct Summary {
  double maximum;
  arma::uword maximum_row;
  arma::uword maximum_col;
  double sum_squares;
  double sum_statistic;
};

inline Summary summarize(const arma::mat& estimates,
                         const double null_second_moment) {
  double maximum = -std::numeric_limits<double>::infinity();
  arma::uword maximum_row = 0U;
  arma::uword maximum_col = 0U;
  long double sum_squares = 0.0L;
  for (arma::uword i = 0; i < estimates.n_rows; ++i) {
    for (arma::uword j = 0; j < estimates.n_cols; ++j) {
      const double value = estimates(i, j);
      if (!std::isfinite(value)) {
        Rcpp::stop("A degenerate rank U-statistic is not finite.");
      }
      if (value > maximum) {
        maximum = value;
        maximum_row = i;
        maximum_col = j;
      }
      sum_squares += static_cast<long double>(value) * value;
    }
  }
  const long double centered = sum_squares -
    static_cast<long double>(estimates.n_elem) * null_second_moment;
  return Summary{
    maximum, maximum_row, maximum_col,
    static_cast<double>(sum_squares), static_cast<double>(centered)
  };
}

inline double null_second_moment(const int measure, const double n) {
  if (measure == 0) {
    return 2.0 * (n * n + 5.0 * n - 32.0) /
      (9.0 * n * (n - 1.0) * (n - 3.0) * (n - 4.0));
  }
  if (measure == 1) {
    return 2.0 * (n * n * n - 3.0 * n * n - 6.0 * n + 10.0) /
      (n * (n - 1.0) * (n - 2.0) * (n - 3.0) * (n - 4.0));
  }
  return 8.0 * (3.0 * n * n + 5.0 * n - 18.0) /
    (75.0 * n * (n - 1.0) * (n - 2.0) * (n - 3.0));
}

}  // namespace ch4completion


//' Standardized residual radii for the Chapter 4 completion methods
//'
//' @param residuals Observation-by-coordinate residual matrix.
//' @param diagonal Strictly positive scale diagonal.
//' @return A vector of stable Euclidean radii.
//' @keywords internal
// [[Rcpp::export]]
arma::vec cpp_ch4_completion_standardized_radii(
    const arma::mat& residuals,
    const arma::vec& diagonal) {
  if (residuals.n_cols != diagonal.n_elem || residuals.n_rows < 1U ||
      residuals.n_cols < 1U || !residuals.is_finite() ||
      !diagonal.is_finite() || arma::any(diagonal <= 0.0)) {
    Rcpp::stop("The residual matrix and positive scale diagonal are incompatible.");
  }
  arma::vec radii(residuals.n_rows, arma::fill::zeros);
  arma::rowvec standardized(residuals.n_cols);
  for (arma::uword i = 0; i < residuals.n_rows; ++i) {
    for (arma::uword j = 0; j < residuals.n_cols; ++j) {
      standardized[j] = residuals(i, j) / std::sqrt(diagonal[j]);
    }
    radii[i] = ch4completion::stable_norm(standardized);
  }
  return radii;
}


//' Weighted spatial-sign alpha quadratic form
//'
//' @param directions Observation-by-asset spatial-sign matrix.
//' @param h Residualized-intercept vector.
//' @param weights Evaluated radial weights.
//' @return The quadratic form and empirical second weight moment.
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List cpp_ch4_completion_weighted_alpha_q(
    const arma::mat& directions,
    const arma::vec& h,
    const arma::vec& weights) {
  if (directions.n_rows != h.n_elem || h.n_elem != weights.n_elem ||
      directions.n_rows < 2U || directions.n_cols < 2U ||
      !directions.is_finite() || !h.is_finite() || !weights.is_finite()) {
    Rcpp::stop("The weighted-alpha directions, h vector, and weights are incompatible.");
  }
  const double h2 = arma::dot(h, h);
  if (!std::isfinite(h2) || !(h2 > 0.0)) {
    Rcpp::stop("The residualized intercept must have positive squared norm.");
  }
  const arma::vec hw = h % weights;
  const arma::rowvec weighted_sum = hw.t() * directions;
  long double off_diagonal = static_cast<long double>(
    arma::dot(weighted_sum, weighted_sum)
  );
  long double weight_square_sum = 0.0L;
  for (arma::uword i = 0; i < directions.n_rows; ++i) {
    const long double weight_square =
      static_cast<long double>(weights[i]) * weights[i];
    weight_square_sum += weight_square;
    off_diagonal -= static_cast<long double>(h[i]) * h[i] *
      weight_square * arma::dot(directions.row(i), directions.row(i));
  }
  const long double psi2 = weight_square_sum /
    static_cast<long double>(weights.n_elem);
  const long double q = static_cast<long double>(directions.n_cols) *
    off_diagonal / static_cast<long double>(h2);
  if (!(psi2 > 0.0L) ||
      !(psi2 <= static_cast<long double>(
          std::numeric_limits<double>::max())) ||
      !std::isfinite(static_cast<double>(q))) {
    Rcpp::stop("The weighted-alpha moment or quadratic form is invalid; no repair is applied.");
  }
  return Rcpp::List::create(
    Rcpp::Named("Q") = static_cast<double>(q),
    Rcpp::Named("psi2") = static_cast<double>(psi2),
    Rcpp::Named("h2") = h2
  );
}


//' Exact degenerate rank-U core for vector independence
//'
//' @param x Observation-by-coordinate first vector block.
//' @param y Observation-by-coordinate second vector block.
//' @param measure Zero for Hoeffding D, one for BKR R, two for tau-star.
//' @param permutations Integer n-by-B matrix permuting X rows.
//' @param max_kernel_evaluations Strict upper bound on symmetrized kernel terms.
//' @param keep_estimates Whether to return the observed coordinate-pair estimates.
//' @param keep_permutation Whether to return permutation sum statistics.
//' @return Exact observed and intrinsic-permutation components.
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List cpp_ch4_completion_vector_u_core(
    const arma::mat& x,
    const arma::mat& y,
    const int measure,
    const arma::imat& permutations,
    const double max_kernel_evaluations,
    const bool keep_estimates,
    const bool keep_permutation) {
  if (x.n_rows != y.n_rows || x.n_cols < 1U || y.n_cols < 1U ||
      !x.is_finite() || !y.is_finite()) {
    Rcpp::stop("`x` and `y` must be finite matrices with matching rows and positive dimensions.");
  }
  if (measure < 0 || measure > 2) {
    Rcpp::stop("Unknown degenerate rank-U measure code.");
  }
  const int order = measure == 0 ? 5 : (measure == 1 ? 6 : 4);
  const int n = static_cast<int>(x.n_rows);
  if (n < order) {
    Rcpp::stop("The sample size is smaller than the selected U-kernel order.");
  }
  if (x.n_cols * y.n_cols < 2U) {
    Rcpp::stop("The high-dimensional max calibration requires p * q >= 2.");
  }
  if (permutations.n_rows != x.n_rows || permutations.n_cols < 2U) {
    Rcpp::stop("`permutations` must have n rows and at least two columns.");
  }
  if (!std::isfinite(max_kernel_evaluations) ||
      !(max_kernel_evaluations >= 1.0)) {
    Rcpp::stop("`max_kernel_evaluations` must be finite and at least one.");
  }
  for (arma::uword j = 0; j < x.n_cols; ++j) {
    if (ch4completion::column_has_tie(x, j)) {
      Rcpp::stop("`x` contains an exact tie; the primary continuous-margin calibration is undefined.");
    }
  }
  for (arma::uword j = 0; j < y.n_cols; ++j) {
    if (ch4completion::column_has_tie(y, j)) {
      Rcpp::stop("`y` contains an exact tie; the primary continuous-margin calibration is undefined.");
    }
  }

  const long double subset_count =
    ch4completion::choose_long_double(n, order);
  const long double ordering_count =
    ch4completion::factorial_long_double(order);
  const long double base_evaluations = subset_count * ordering_count *
    static_cast<long double>(x.n_cols) * y.n_cols;
  const long double total_evaluations = base_evaluations *
    static_cast<long double>(permutations.n_cols + 1U);
  if (total_evaluations >
      static_cast<long double>(max_kernel_evaluations)) {
    Rcpp::stop(
      "The exact high-order U-kernel workload exceeds `max_kernel_evaluations`; no approximate or silently truncated computation is used."
    );
  }

  const std::vector<std::vector<int> > subsets =
    ch4completion::combinations(n, order);
  const std::vector<std::vector<int> > kernel_orderings =
    ch4completion::orderings(order);
  const double null_moment = ch4completion::null_second_moment(
    measure, static_cast<double>(n)
  );
  if (!std::isfinite(null_moment) || !(null_moment > 0.0)) {
    Rcpp::stop("The primary finite-sample null second moment is invalid.");
  }

  std::size_t interrupt_counter = 0U;
  const arma::mat estimates = ch4completion::u_statistics(
    x, y, measure, subsets, kernel_orderings, nullptr,
    interrupt_counter
  );
  const ch4completion::Summary observed = ch4completion::summarize(
    estimates, null_moment
  );

  arma::vec permutation_statistics(permutations.n_cols);
  std::vector<int> permuted_x(n);
  std::vector<int> seen(n);
  for (arma::uword draw = 0; draw < permutations.n_cols; ++draw) {
    std::fill(seen.begin(), seen.end(), 0);
    for (int row = 0; row < n; ++row) {
      const int source = permutations(row, draw);
      if (source < 1 || source > n || seen[source - 1] != 0) {
        Rcpp::stop("Every column of `permutations` must be a permutation of 1:n.");
      }
      seen[source - 1] = 1;
      permuted_x[row] = source - 1;
    }
    const arma::mat current = ch4completion::u_statistics(
      x, y, measure, subsets, kernel_orderings, &permuted_x,
      interrupt_counter
    );
    permutation_statistics[draw] = ch4completion::summarize(
      current, null_moment
    ).sum_statistic;
  }
  const double permutation_mean = arma::mean(permutation_statistics);
  const double permutation_variance = arma::accu(arma::square(
    permutation_statistics - permutation_mean
  )) / static_cast<double>(permutations.n_cols - 1U);
  if (!std::isfinite(permutation_variance) ||
      !(permutation_variance > 0.0)) {
    Rcpp::stop(
      "The intrinsic permutation variance estimate is not strictly positive; increase `B` or use non-degenerate data."
    );
  }

  return Rcpp::List::create(
    Rcpp::Named("maximum") = observed.maximum,
    Rcpp::Named("maximum_index") = Rcpp::IntegerVector::create(
      static_cast<int>(observed.maximum_row + 1U),
      static_cast<int>(observed.maximum_col + 1U)
    ),
    Rcpp::Named("sum_squares") = observed.sum_squares,
    Rcpp::Named("sum_statistic") = observed.sum_statistic,
    Rcpp::Named("null_second_moment") = null_moment,
    Rcpp::Named("permutation_variance") = permutation_variance,
    Rcpp::Named("permutation_mean") = permutation_mean,
    Rcpp::Named("estimates") = keep_estimates ?
      Rcpp::wrap(estimates) : R_NilValue,
    Rcpp::Named("permutation_statistics") = keep_permutation ?
      Rcpp::wrap(permutation_statistics) : R_NilValue,
    Rcpp::Named("kernel_order") = order,
    Rcpp::Named("subsets") = static_cast<double>(subset_count),
    Rcpp::Named("orderings_per_subset") =
      static_cast<double>(ordering_count),
    Rcpp::Named("base_kernel_evaluations") =
      static_cast<double>(base_evaluations),
    Rcpp::Named("total_kernel_evaluations") =
      static_cast<double>(total_evaluations)
  );
}
