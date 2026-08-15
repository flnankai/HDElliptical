// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp17)]]

#include <RcppArmadillo.h>

#include <algorithm>
#include <cmath>
#include <limits>

namespace {

void ch5_require_finite(const arma::mat& x, const char* name) {
  if (!x.is_finite()) {
    Rcpp::stop("`%s` must contain only finite values.", name);
  }
}


void ch5_require_finite(const arma::vec& x, const char* name) {
  if (!x.is_finite()) {
    Rcpp::stop("`%s` must contain only finite values.", name);
  }
}


void ch5_require_symmetric(const arma::mat& x, const char* name) {
  if (x.n_rows != x.n_cols) {
    Rcpp::stop("`%s` must be square.", name);
  }
  double maximum = 1.0;
  if (x.n_elem > 0) maximum = std::max(maximum, arma::abs(x).max());
  const double tolerance = 100.0 * std::numeric_limits<double>::epsilon() *
    maximum;
  if (arma::abs(x - x.t()).max() > tolerance) {
    Rcpp::stop("`%s` must be symmetric.", name);
  }
}

} // namespace


// [[Rcpp::export]]
arma::vec cpp_ch5_classical_linear_scores(const arma::mat& x,
                                           const arma::vec& coefficients,
                                           const double intercept) {
  if (x.n_cols != coefficients.n_elem) {
    Rcpp::stop("The coefficient length must equal the number of columns in `x`.");
  }
  ch5_require_finite(x, "x");
  ch5_require_finite(coefficients, "coefficients");
  if (!std::isfinite(intercept)) {
    Rcpp::stop("`intercept` must be finite.");
  }
  arma::vec score = x * coefficients + intercept;
  if (!score.is_finite()) {
    Rcpp::stop("A linear classifier score is outside the finite double range.");
  }
  return score;
}


// [[Rcpp::export]]
arma::vec cpp_ch5_classical_quadratic_scores(const arma::mat& x,
                                              const arma::mat& quadratic,
                                              const arma::vec& linear,
                                              const double intercept) {
  const arma::uword p = x.n_cols;
  if (quadratic.n_rows != p || quadratic.n_cols != p ||
      linear.n_elem != p) {
    Rcpp::stop("The quadratic and linear coefficients must match `x`.");
  }
  ch5_require_finite(x, "x");
  ch5_require_finite(quadratic, "quadratic");
  ch5_require_finite(linear, "linear");
  ch5_require_symmetric(quadratic, "quadratic");
  if (!std::isfinite(intercept)) {
    Rcpp::stop("`intercept` must be finite.");
  }

  arma::vec score(x.n_rows, arma::fill::zeros);
  for (arma::uword i = 0; i < x.n_rows; ++i) {
    const arma::vec row = x.row(i).t();
    const double value = arma::dot(row, quadratic * row) +
      arma::dot(linear, row) + intercept;
    if (!std::isfinite(value)) {
      Rcpp::stop("A quadratic classifier score is outside the finite double range.");
    }
    score(i) = value;
  }
  return score;
}


// [[Rcpp::export]]
Rcpp::List cpp_ch5_classical_mahalanobis_pairs(
    const arma::mat& x,
    const arma::vec& location1,
    const arma::mat& precision1,
    const arma::vec& location2,
    const arma::mat& precision2) {
  const arma::uword p = x.n_cols;
  if (location1.n_elem != p || location2.n_elem != p ||
      precision1.n_rows != p || precision1.n_cols != p ||
      precision2.n_rows != p || precision2.n_cols != p) {
    Rcpp::stop("Locations and precision matrices must match the columns of `x`.");
  }
  ch5_require_finite(x, "x");
  ch5_require_finite(location1, "location1");
  ch5_require_finite(location2, "location2");
  ch5_require_finite(precision1, "precision1");
  ch5_require_finite(precision2, "precision2");
  ch5_require_symmetric(precision1, "precision1");
  ch5_require_symmetric(precision2, "precision2");

  arma::mat factor1;
  arma::mat factor2;
  if (!arma::chol(factor1, precision1) || !arma::chol(factor2, precision2)) {
    Rcpp::stop("Both precision matrices must be strictly positive definite.");
  }
  arma::vec distance1(x.n_rows, arma::fill::zeros);
  arma::vec distance2(x.n_rows, arma::fill::zeros);
  for (arma::uword i = 0; i < x.n_rows; ++i) {
    const arma::vec difference1 = x.row(i).t() - location1;
    const arma::vec difference2 = x.row(i).t() - location2;
    const arma::vec whitened1 = factor1 * difference1;
    const arma::vec whitened2 = factor2 * difference2;
    distance1(i) = arma::dot(whitened1, whitened1);
    distance2(i) = arma::dot(whitened2, whitened2);
    if (!std::isfinite(distance1(i)) || !std::isfinite(distance2(i))) {
      Rcpp::stop("A squared Mahalanobis distance is outside the finite double range.");
    }
  }
  return Rcpp::List::create(
    Rcpp::Named("distance1") = distance1,
    Rcpp::Named("distance2") = distance2
  );
}


// [[Rcpp::export]]
Rcpp::List cpp_ch5_classical_two_class_moments(
    const arma::mat& x,
    const Rcpp::IntegerVector& class_id) {
  const arma::uword n = x.n_rows;
  const arma::uword p = x.n_cols;
  if (n < 2 || p < 1 || class_id.size() != static_cast<R_xlen_t>(n)) {
    Rcpp::stop("`x` and `class_id` have incompatible dimensions.");
  }
  ch5_require_finite(x, "x");
  std::vector<arma::uword> rows1;
  std::vector<arma::uword> rows2;
  rows1.reserve(n);
  rows2.reserve(n);
  for (arma::uword i = 0; i < n; ++i) {
    if (class_id[i] == NA_INTEGER) {
      Rcpp::stop("`class_id` must not contain missing values.");
    }
    if (class_id[i] == 1) rows1.push_back(i);
    else if (class_id[i] == 2) rows2.push_back(i);
    else Rcpp::stop("`class_id` entries must be one or two.");
  }
  if (rows1.empty() || rows2.empty()) {
    Rcpp::stop("Both classes must contain observations.");
  }

  arma::uvec index1(rows1.size());
  arma::uvec index2(rows2.size());
  for (arma::uword i = 0; i < index1.n_elem; ++i) index1(i) = rows1[i];
  for (arma::uword i = 0; i < index2.n_elem; ++i) index2(i) = rows2[i];
  const arma::mat x1 = x.rows(index1);
  const arma::mat x2 = x.rows(index2);
  const arma::rowvec mean1_row = arma::mean(x1, 0);
  const arma::rowvec mean2_row = arma::mean(x2, 0);
  const arma::mat centered1 = x1.each_row() - mean1_row;
  const arma::mat centered2 = x2.each_row() - mean2_row;
  const arma::mat sscp1 = centered1.t() * centered1;
  const arma::mat sscp2 = centered2.t() * centered2;
  if (!mean1_row.is_finite() || !mean2_row.is_finite() ||
      !sscp1.is_finite() || !sscp2.is_finite()) {
    Rcpp::stop("A two-class mean or within-class cross-product is outside the finite double range.");
  }
  return Rcpp::List::create(
    Rcpp::Named("mean1") = mean1_row.t(),
    Rcpp::Named("mean2") = mean2_row.t(),
    Rcpp::Named("sscp1") = sscp1,
    Rcpp::Named("sscp2") = sscp2,
    Rcpp::Named("n1") = static_cast<int>(rows1.size()),
    Rcpp::Named("n2") = static_cast<int>(rows2.size())
  );
}


// Equation (4.3) of Fan and Fan (2008), evaluated for every leading block
// after the R layer has ordered features by decreasing absolute Welch t.
// [[Rcpp::export]]
Rcpp::List cpp_ch5_classical_fair_criterion(
    const arma::vec& t_squared,
    const arma::mat& ordered_correlation,
    const int n1,
    const int n2) {
  const arma::uword p = t_squared.n_elem;
  if (p < 1 || ordered_correlation.n_rows != p ||
      ordered_correlation.n_cols != p) {
    Rcpp::stop("The FAIR statistics and correlation matrix have incompatible dimensions.");
  }
  if (n1 < 2 || n2 < 2) {
    Rcpp::stop("FAIR requires at least two observations per class.");
  }
  ch5_require_finite(t_squared, "t_squared");
  ch5_require_finite(ordered_correlation, "ordered_correlation");
  ch5_require_symmetric(ordered_correlation, "ordered_correlation");
  if (arma::any(t_squared < 0.0)) {
    Rcpp::stop("`t_squared` must be non-negative.");
  }

  arma::vec criterion(p, arma::fill::zeros);
  arma::vec lambda_max(p, arma::fill::zeros);
  long double cumulative_t_squared = 0.0L;
  const long double n = static_cast<long double>(n1) +
    static_cast<long double>(n2);
  const long double product = static_cast<long double>(n1) *
    static_cast<long double>(n2);
  for (arma::uword m = 1; m <= p; ++m) {
    cumulative_t_squared += static_cast<long double>(t_squared(m - 1));
    const arma::mat block = ordered_correlation.submat(0, 0, m - 1, m - 1);
    arma::vec eigenvalues;
    if (!arma::eig_sym(eigenvalues, block) || eigenvalues.n_elem != m) {
      Rcpp::stop("An exact FAIR truncated-correlation eigendecomposition failed.");
    }
    const double largest = eigenvalues(m - 1);
    if (!(largest > 0.0) || !std::isfinite(largest)) {
      Rcpp::stop("A FAIR truncated-correlation largest eigenvalue is not positive and finite.");
    }
    lambda_max(m - 1) = largest;
    const long double adjusted = cumulative_t_squared +
      static_cast<long double>(m) *
      (static_cast<long double>(n1) - static_cast<long double>(n2)) / n;
    const long double numerator = n * adjusted * adjusted;
    const long double denominator = static_cast<long double>(largest) *
      (static_cast<long double>(m) * product +
       product * cumulative_t_squared);
    const long double value = numerator / denominator;
    const double stored = static_cast<double>(value);
    if (!(denominator > 0.0L) || !std::isfinite(value) ||
        !std::isfinite(stored)) {
      Rcpp::stop("A primary FAIR feature-count criterion is not finite.");
    }
    criterion(m - 1) = stored;
  }
  return Rcpp::List::create(
    Rcpp::Named("criterion") = criterion,
    Rcpp::Named("lambda_max") = lambda_max
  );
}


// The Shao et al. threshold keeps diagonal entries and uses a strict
// off-diagonal comparison. No positive-definite repair occurs here.
// [[Rcpp::export]]
arma::mat cpp_ch5_classical_hard_threshold_covariance(
    const arma::mat& covariance,
    const double threshold) {
  ch5_require_finite(covariance, "covariance");
  ch5_require_symmetric(covariance, "covariance");
  if (!(threshold >= 0.0) || !std::isfinite(threshold)) {
    Rcpp::stop("`threshold` must be finite and non-negative.");
  }
  arma::mat answer = (covariance + covariance.t()) / 2.0;
  const arma::uword p = answer.n_rows;
  for (arma::uword i = 0; i < p; ++i) {
    for (arma::uword j = i + 1; j < p; ++j) {
      if (!(std::abs(answer(i, j)) > threshold)) {
        answer(i, j) = 0.0;
        answer(j, i) = 0.0;
      }
    }
  }
  return answer;
}
