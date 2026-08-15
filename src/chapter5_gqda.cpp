// Chapter 5: generalized quadratic discriminant scores.
//
// The primary GQDA rule compares the difference of two Mahalanobis
// distances with c times a signed log-determinant contrast.  This kernel
// forms both residuals in long double and uses one common row scale, so a
// large common translation or a representable common change of units does
// not introduce an avoidable subtraction or squaring overflow.

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp17)]]

#include <RcppArmadillo.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <string>
#include <vector>

namespace {

long double c5gqda_quadratic_scaled(
    const std::vector<long double>& residual,
    const arma::mat& precision) {
  const arma::uword p = precision.n_rows;
  long double answer = 0.0L;
  for (arma::uword ii = 0; ii < p; ++ii) {
    long double row = 0.0L;
    for (arma::uword jj = 0; jj < p; ++jj) {
      row += static_cast<long double>(precision(ii, jj)) * residual[jj];
    }
    answer += residual[ii] * row;
  }
  return answer;
}

double c5gqda_restore(long double scaled, long double scale) {
  if (scaled == 0.0L) return 0.0;
  const long double sign = scaled > 0.0L ? 1.0L : -1.0L;
  const long double log_absolute =
    std::log(std::fabs(scaled)) + 2.0L * std::log(scale);
  const long double log_max =
    std::log(static_cast<long double>(std::numeric_limits<double>::max()));
  const long double log_min =
    std::log(static_cast<long double>(std::numeric_limits<double>::denorm_min()));
  if (log_absolute > log_max) {
    return sign > 0.0L ? R_PosInf : R_NegInf;
  }
  if (log_absolute < log_min) {
    return sign > 0.0L ? 0.0 : -0.0;
  }
  return static_cast<double>(sign * std::exp(log_absolute));
}

}  // namespace

// [[Rcpp::export]]
Rcpp::List cpp_ch5_gqda_components(
    const arma::mat& x,
    const arma::vec& location1,
    const arma::vec& location2,
    const arma::mat& precision1,
    const arma::mat& precision2) {
  const arma::uword n = x.n_rows;
  const arma::uword p = x.n_cols;
  if (p == 0 || location1.n_elem != p || location2.n_elem != p ||
      precision1.n_rows != p || precision1.n_cols != p ||
      precision2.n_rows != p || precision2.n_cols != p) {
    Rcpp::stop("GQDA component dimensions are not conformable.");
  }
  if (!x.is_finite() || !location1.is_finite() || !location2.is_finite() ||
      !precision1.is_finite() || !precision2.is_finite()) {
    Rcpp::stop("GQDA components require finite inputs.");
  }

  arma::vec distance1(n, arma::fill::zeros);
  arma::vec distance2(n, arma::fill::zeros);
  arma::vec difference(n, arma::fill::zeros);
  arma::vec log_abs_difference(n, arma::fill::zeros);
  arma::ivec difference_sign(n, arma::fill::zeros);

  std::vector<long double> residual1(p);
  std::vector<long double> residual2(p);
  for (arma::uword ii = 0; ii < n; ++ii) {
    long double row_scale = 0.0L;
    for (arma::uword jj = 0; jj < p; ++jj) {
      const long double first = static_cast<long double>(x(ii, jj)) -
        static_cast<long double>(location1(jj));
      const long double second = static_cast<long double>(x(ii, jj)) -
        static_cast<long double>(location2(jj));
      residual1[jj] = first;
      residual2[jj] = second;
      row_scale = std::max(row_scale, std::fabs(first));
      row_scale = std::max(row_scale, std::fabs(second));
    }
    if (row_scale == 0.0L) row_scale = 1.0L;
    for (arma::uword jj = 0; jj < p; ++jj) {
      residual1[jj] /= row_scale;
      residual2[jj] /= row_scale;
    }
    const long double q1_scaled =
      c5gqda_quadratic_scaled(residual1, precision1);
    const long double q2_scaled =
      c5gqda_quadratic_scaled(residual2, precision2);
    const long double roundoff_scale = std::max(
      1.0L, std::max(std::fabs(q1_scaled), std::fabs(q2_scaled)));
    const long double tolerance =
      256.0L * static_cast<long double>(std::numeric_limits<double>::epsilon()) *
      static_cast<long double>(p) * roundoff_scale;
    if (q1_scaled < -tolerance || q2_scaled < -tolerance) {
      Rcpp::stop("A supplied GQDA precision produced a negative quadratic form.");
    }
    const long double q1 = q1_scaled < 0.0L ? 0.0L : q1_scaled;
    const long double q2 = q2_scaled < 0.0L ? 0.0L : q2_scaled;
    const long double delta = q2 - q1;
    distance1(ii) = c5gqda_restore(q1, row_scale);
    distance2(ii) = c5gqda_restore(q2, row_scale);
    difference(ii) = c5gqda_restore(delta, row_scale);
    if (delta == 0.0L) {
      difference_sign(ii) = 0;
      log_abs_difference(ii) = R_NegInf;
    } else {
      difference_sign(ii) = delta > 0.0L ? 1 : -1;
      log_abs_difference(ii) = static_cast<double>(
        std::log(std::fabs(delta)) + 2.0L * std::log(row_scale));
    }
  }

  return Rcpp::List::create(
    Rcpp::Named("distance1.squared") = distance1,
    Rcpp::Named("distance2.squared") = distance2,
    Rcpp::Named("difference") = difference,
    Rcpp::Named("difference.sign") = difference_sign,
    Rcpp::Named("difference.log.absolute") = log_abs_difference,
    Rcpp::Named("numerical.contract") = std::string(
      "long-double residual subtraction and one common row scale; no repair"));
}
