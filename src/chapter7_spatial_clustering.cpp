#include <RcppArmadillo.h>
#include <algorithm>
#include <array>
#include <limits>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp17)]]

namespace {

void require_matrix(const arma::mat& x, const char* name) {
  if (x.n_rows == 0 || x.n_cols == 0 || !x.is_finite()) {
    Rcpp::stop("`%s` must be a non-empty finite numeric matrix.", name);
  }
}

std::vector<arma::uword> checked_active(
    const Rcpp::IntegerVector& active, const arma::uword p) {
  if (active.size() == 0) {
    Rcpp::stop("`active` must contain at least one coordinate.");
  }
  std::vector<arma::uword> answer;
  answer.reserve(active.size());
  std::vector<bool> seen(p, false);
  for (R_xlen_t j = 0; j < active.size(); ++j) {
    if (Rcpp::IntegerVector::is_na(active[j]) ||
        active[j] < 0 || static_cast<arma::uword>(active[j]) >= p) {
      Rcpp::stop("`active` contains an invalid zero-based coordinate.");
    }
    const arma::uword index = static_cast<arma::uword>(active[j]);
    if (seen[index]) {
      Rcpp::stop("`active` must not contain duplicate coordinates.");
    }
    seen[index] = true;
    answer.push_back(index);
  }
  return answer;
}

arma::uvec checked_labels(const Rcpp::IntegerVector& labels,
                          const arma::uword n,
                          const arma::uword k) {
  if (static_cast<arma::uword>(labels.size()) != n) {
    Rcpp::stop("`labels` must contain one label per observation.");
  }
  arma::uvec answer(n);
  for (arma::uword i = 0; i < n; ++i) {
    if (Rcpp::IntegerVector::is_na(labels[i]) || labels[i] < 1 ||
        static_cast<arma::uword>(labels[i]) > k) {
      Rcpp::stop("`labels` must contain integers from 1 through K.");
    }
    answer(i) = static_cast<arma::uword>(labels[i] - 1);
  }
  return answer;
}

struct ScaledDifference {
  arma::rowvec scaled;
  double coordinate_scale;
  double difference_scale;
};

template <typename Coordinate>
ScaledDifference scaled_difference(
    const arma::rowvec& left, const arma::rowvec& right,
    const arma::uword count, const Coordinate& coordinate) {
  arma::rowvec difference(count, arma::fill::zeros);
  bool direct_finite = true;
  for (arma::uword index = 0; index < count; ++index) {
    const arma::uword j = coordinate(index);
    difference(index) = left(j) - right(j);
    if (!std::isfinite(difference(index))) {
      direct_finite = false;
    }
  }
  double coordinate_scale = 1.0;
  if (!direct_finite) {
    coordinate_scale = 0.0;
    for (arma::uword index = 0; index < count; ++index) {
      const arma::uword j = coordinate(index);
      coordinate_scale = std::max(
        coordinate_scale, std::max(std::abs(left(j)), std::abs(right(j)))
      );
    }
    if (coordinate_scale == 0.0) {
      return ScaledDifference{difference, 1.0, 0.0};
    }
    for (arma::uword index = 0; index < count; ++index) {
      const arma::uword j = coordinate(index);
      difference(index) = left(j) / coordinate_scale -
        right(j) / coordinate_scale;
    }
  }
  const double difference_scale = arma::abs(difference).max();
  if (difference_scale == 0.0) {
    difference.zeros();
    return ScaledDifference{difference, coordinate_scale, 0.0};
  }
  difference /= difference_scale;
  return ScaledDifference{difference, coordinate_scale, difference_scale};
}

double stable_product3(const double first, const double second,
                       const double third) {
  if (!std::isfinite(first) || !std::isfinite(second) ||
      !std::isfinite(third) || first < 0.0 || second < 0.0 || third < 0.0) {
    return std::numeric_limits<double>::infinity();
  }
  if (first == 0.0 || second == 0.0 || third == 0.0) {
    return 0.0;
  }
  std::array<double, 3> factors{first, second, third};
  std::sort(factors.begin(), factors.end());
  return (factors[0] * factors[2]) * factors[1];
}


double active_distance(const arma::rowvec& left,
                       const arma::rowvec& right,
                       const std::vector<arma::uword>& active) {
  const ScaledDifference residual = scaled_difference(
    left, right, static_cast<arma::uword>(active.size()),
    [&active](const arma::uword index) { return active[index]; }
  );
  if (residual.difference_scale == 0.0) {
    return 0.0;
  }
  const double answer = stable_product3(
    residual.coordinate_scale, residual.difference_scale,
    arma::norm(residual.scaled, 2)
  );
  if (!std::isfinite(answer)) {
    Rcpp::stop("A Euclidean distance overflowed the finite numeric range.");
  }
  return answer;
}

bool stable_direction(const arma::rowvec& left,
                      const arma::rowvec& right,
                      const double zero_tol,
                      arma::rowvec& direction) {
  const ScaledDifference residual = scaled_difference(
    left, right, left.n_elem,
    [](const arma::uword index) { return index; }
  );
  if (residual.difference_scale == 0.0) {
    direction.zeros(left.n_elem);
    return false;
  }
  const double scaled_norm = arma::norm(residual.scaled, 2);
  double scaled_tolerance = zero_tol / residual.coordinate_scale;
  if (scaled_tolerance != 0.0) {
    scaled_tolerance /= residual.difference_scale;
  }
  if (!(scaled_norm > scaled_tolerance)) {
    direction.zeros(left.n_elem);
    return false;
  }
  direction = residual.scaled / scaled_norm;
  return true;
}

double metric_distance(const arma::rowvec& left,
                       const arma::rowvec& right,
                       const arma::mat& inverse_cholesky) {
  const ScaledDifference residual = scaled_difference(
    left, right, left.n_elem,
    [](const arma::uword index) { return index; }
  );
  if (residual.difference_scale == 0.0) {
    return 0.0;
  }
  const arma::rowvec transformed =
    residual.scaled * inverse_cholesky.t();
  const double distance_root = stable_product3(
    residual.coordinate_scale, residual.difference_scale,
    arma::norm(transformed, 2)
  );
  if (!std::isfinite(distance_root) ||
      (distance_root != 0.0 &&
       distance_root > std::numeric_limits<double>::max() / distance_root)) {
    Rcpp::stop("A metric distance overflowed the finite numeric range.");
  }
  return distance_root * distance_root;
}

Rcpp::List assignment_result(const arma::mat& distances,
                             const bool keep_matrix) {
  const arma::uword n = distances.n_rows;
  const arma::uword k = distances.n_cols;
  Rcpp::IntegerVector labels(n);
  arma::vec assigned(n, arma::fill::zeros);
  arma::uword tie_count = 0;
  for (arma::uword i = 0; i < n; ++i) {
    arma::uword best = 0;
    double best_value = distances(i, 0);
    bool tied = false;
    for (arma::uword cluster = 1; cluster < k; ++cluster) {
      const double candidate = distances(i, cluster);
      if (candidate < best_value) {
        best = cluster;
        best_value = candidate;
        tied = false;
      } else if (candidate == best_value) {
        tied = true;
      }
    }
    labels[i] = static_cast<int>(best + 1);
    assigned(i) = best_value;
    if (tied) {
      ++tie_count;
    }
  }
  return Rcpp::List::create(
    Rcpp::Named("labels") = labels,
    Rcpp::Named("assigned_distance") = assigned,
    Rcpp::Named("distance_matrix") = keep_matrix ?
      Rcpp::wrap(distances) : R_NilValue,
    Rcpp::Named("tie_count") = static_cast<int>(tie_count),
    Rcpp::Named("tie_rule") = "first minimum cluster index"
  );
}

}  // namespace


//' Euclidean assignment kernel for Chapter 7 spatial clustering
//'
//' @param x Finite observation-by-variable matrix.
//' @param centers Finite cluster-by-variable center matrix.
//' @param active Zero-based active coordinate indices.
//' @param keep_matrix Whether to retain every observation-center distance.
//' @return Internal deterministic assignment and distance diagnostics.
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List cpp_ch7sc_assign_euclidean(
    const arma::mat& x, const arma::mat& centers,
    const Rcpp::IntegerVector& active, const bool keep_matrix) {
  require_matrix(x, "x");
  require_matrix(centers, "centers");
  if (x.n_cols != centers.n_cols) {
    Rcpp::stop("`x` and `centers` must have the same number of columns.");
  }
  const std::vector<arma::uword> coordinates = checked_active(
    active, x.n_cols
  );
  arma::mat distances(x.n_rows, centers.n_rows, arma::fill::zeros);
  for (arma::uword i = 0; i < x.n_rows; ++i) {
    for (arma::uword cluster = 0; cluster < centers.n_rows; ++cluster) {
      distances(i, cluster) = active_distance(
        x.row(i), centers.row(cluster), coordinates
      );
    }
  }
  return assignment_result(distances, keep_matrix);
}


//' SSCM and exact inverse for Chapter 7 spatial clustering
//'
//' @param x Finite observation-by-variable matrix.
//' @param centers Finite cluster-by-variable center matrix.
//' @param labels One-based cluster labels.
//' @param lambda Strictly positive ridge in the stated SSCM definition.
//' @param zero_tol Non-negative zero-residual tolerance.
//' @param keep_signs Whether to retain the residual spatial signs.
//' @return Internal SSCM, inverse, and finite-sample certificates.
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List cpp_ch7sc_sscm_metric(
    const arma::mat& x, const arma::mat& centers,
    const Rcpp::IntegerVector& labels, const double lambda,
    const double zero_tol, const bool keep_signs) {
  require_matrix(x, "x");
  require_matrix(centers, "centers");
  if (x.n_cols != centers.n_cols) {
    Rcpp::stop("`x` and `centers` must have the same number of columns.");
  }
  if (!std::isfinite(lambda) || lambda <= 0.0) {
    Rcpp::stop("`lambda` must be a finite strictly positive number.");
  }
  if (!std::isfinite(zero_tol) || zero_tol < 0.0) {
    Rcpp::stop("`zero_tol` must be finite and non-negative.");
  }
  const arma::uvec memberships = checked_labels(
    labels, x.n_rows, centers.n_rows
  );
  arma::mat signs(x.n_rows, x.n_cols, arma::fill::zeros);
  arma::uword n_zero = 0;
  for (arma::uword i = 0; i < x.n_rows; ++i) {
    arma::rowvec direction(x.n_cols, arma::fill::zeros);
    if (stable_direction(
          x.row(i), centers.row(memberships(i)), zero_tol,
      direction)) {
      signs.row(i) = direction;
    } else {
      ++n_zero;
    }
  }
  arma::mat metric = signs.t() * signs /
    static_cast<double>(x.n_rows);
  metric.diag() += lambda;
  metric = 0.5 * (metric + metric.t());
  arma::vec eigenvalues;
  if (!arma::eig_sym(eigenvalues, metric) || !eigenvalues.is_finite() ||
      !(eigenvalues.min() > 0.0)) {
    Rcpp::stop("The stated SSCM plus `lambda * I` is not numerically SPD.");
  }
  arma::mat inverse_metric;
  if (!arma::inv_sympd(inverse_metric, metric) ||
      !inverse_metric.is_finite()) {
    Rcpp::stop("Exact SPD inversion of the stated SSCM failed.");
  }
  const arma::mat identity = arma::eye(metric.n_rows, metric.n_cols);
  const double inverse_residual = arma::abs(
    metric * inverse_metric - identity
  ).max();
  const double expected_trace =
    static_cast<double>(x.n_rows - n_zero) /
      static_cast<double>(x.n_rows) +
    static_cast<double>(x.n_cols) * lambda;
  return Rcpp::List::create(
    Rcpp::Named("metric") = metric,
    Rcpp::Named("inverse") = inverse_metric,
    Rcpp::Named("signs") = keep_signs ? Rcpp::wrap(signs) : R_NilValue,
    Rcpp::Named("n_zero") = static_cast<int>(n_zero),
    Rcpp::Named("trace_expected") = expected_trace,
    Rcpp::Named("trace_error") = std::abs(arma::trace(metric) - expected_trace),
    Rcpp::Named("minimum_eigenvalue") = eigenvalues.min(),
    Rcpp::Named("inverse_residual") = inverse_residual,
    Rcpp::Named("inverse_method") = "exact SPD inverse; no pseudoinverse"
  );
}


//' SSCM metric assignment kernel for Chapter 7
//'
//' @param x Finite observation-by-variable matrix.
//' @param centers Finite cluster-by-variable center matrix.
//' @param inverse_metric Finite symmetric positive-definite inverse SSCM.
//' @param keep_matrix Whether to retain every squared metric distance.
//' @return Internal deterministic assignment and distance diagnostics.
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List cpp_ch7sc_assign_metric(
    const arma::mat& x, const arma::mat& centers,
    const arma::mat& inverse_metric, const bool keep_matrix) {
  require_matrix(x, "x");
  require_matrix(centers, "centers");
  require_matrix(inverse_metric, "inverse_metric");
  if (x.n_cols != centers.n_cols ||
      inverse_metric.n_rows != x.n_cols ||
      inverse_metric.n_cols != x.n_cols) {
    Rcpp::stop("The metric and center dimensions do not match `x`.");
  }
  const double symmetry_error = arma::abs(
    inverse_metric - inverse_metric.t()
  ).max();
  arma::vec eigenvalues;
  if (symmetry_error != 0.0 ||
      !arma::eig_sym(eigenvalues, inverse_metric) ||
      !eigenvalues.is_finite() || !(eigenvalues.min() > 0.0)) {
    Rcpp::stop("`inverse_metric` must be exactly symmetric and SPD.");
  }
  arma::mat inverse_cholesky;
  if (!arma::chol(inverse_cholesky, inverse_metric)) {
    Rcpp::stop(
      "Cholesky factorization of `inverse_metric` failed; no repair was applied."
    );
  }
  arma::mat distances(x.n_rows, centers.n_rows, arma::fill::zeros);
  for (arma::uword i = 0; i < x.n_rows; ++i) {
    for (arma::uword cluster = 0; cluster < centers.n_rows; ++cluster) {
      distances(i, cluster) = metric_distance(
        x.row(i), centers.row(cluster), inverse_cholesky
      );
    }
  }
  Rcpp::List answer = assignment_result(distances, keep_matrix);
  answer["distance_definition"] = "squared SSCM metric distance";
  return answer;
}


//' Across-center feature scores for Sparse--SM
//'
//' @param centers Finite cluster-by-variable center matrix.
//' @return Internal center average and hard-screening scores.
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List cpp_ch7sc_feature_scores(const arma::mat& centers) {
  require_matrix(centers, "centers");
  const arma::rowvec average = arma::mean(centers, 0);
  const arma::rowvec scores = arma::sum(
    arma::abs(centers.each_row() - average), 0
  );
  if (!scores.is_finite()) {
    Rcpp::stop("The across-center feature scores are non-finite.");
  }
  return Rcpp::List::create(
    Rcpp::Named("average_center") = average,
    Rcpp::Named("scores") = scores
  );
}


//' Active-subspace geometry for Sparse--SM selectors
//'
//' @param x Finite observation-by-variable matrix.
//' @param centers Finite cluster-by-variable center matrix.
//' @param labels One-based cluster labels.
//' @param overall Finite vector in the same coordinate convention as centers.
//' @param active Zero-based active coordinate indices.
//' @return Internal between- and within-spatial-median geometry.
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List cpp_ch7sc_geometry(
    const arma::mat& x, const arma::mat& centers,
    const Rcpp::IntegerVector& labels, const arma::rowvec& overall,
    const Rcpp::IntegerVector& active) {
  require_matrix(x, "x");
  require_matrix(centers, "centers");
  if (x.n_cols != centers.n_cols || overall.n_elem != x.n_cols ||
      !overall.is_finite()) {
    Rcpp::stop("The center and overall dimensions do not match `x`.");
  }
  const std::vector<arma::uword> coordinates = checked_active(
    active, x.n_cols
  );
  const arma::uvec memberships = checked_labels(
    labels, x.n_rows, centers.n_rows
  );
  arma::uvec sizes(centers.n_rows, arma::fill::zeros);
  double within_sum = 0.0;
  for (arma::uword i = 0; i < x.n_rows; ++i) {
    ++sizes(memberships(i));
    within_sum += active_distance(
      x.row(i), centers.row(memberships(i)), coordinates
    );
  }
  if (arma::any(sizes == 0)) {
    Rcpp::stop("Every cluster must be non-empty for selector geometry.");
  }
  double between_dispersion = 0.0;
  for (arma::uword cluster = 0; cluster < centers.n_rows; ++cluster) {
    const double distance = active_distance(
      centers.row(cluster), overall, coordinates
    );
    between_dispersion += static_cast<double>(sizes(cluster)) *
      distance * distance;
  }
  double pair_sum = 0.0;
  arma::uword pair_count = 0;
  for (arma::uword first = 0; first < centers.n_rows; ++first) {
    for (arma::uword second = first + 1; second < centers.n_rows; ++second) {
      pair_sum += active_distance(
        centers.row(first), centers.row(second), coordinates
      );
      ++pair_count;
    }
  }
  const double average_between = pair_count == 0 ? NA_REAL :
    pair_sum / static_cast<double>(pair_count);
  if (!std::isfinite(within_sum) || !std::isfinite(between_dispersion) ||
      (pair_count > 0 && !std::isfinite(average_between))) {
    Rcpp::stop("The selector geometry is non-finite.");
  }
  return Rcpp::List::create(
    Rcpp::Named("sizes") = sizes,
    Rcpp::Named("within_sum") = within_sum,
    Rcpp::Named("within_mean") = within_sum /
      static_cast<double>(x.n_rows),
    Rcpp::Named("between_dispersion") = between_dispersion,
    Rcpp::Named("average_between") = average_between,
    Rcpp::Named("pair_count") = static_cast<int>(pair_count),
    Rcpp::Named("coordinate_count") = static_cast<int>(coordinates.size())
  );
}
