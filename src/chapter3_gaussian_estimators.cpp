// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp17)]]

#include <RcppArmadillo.h>

#include <cmath>
#include <string>

namespace {

double threshold_value(const double z, const double lambda,
                       const std::string& rule, const double shape) {
  if (!std::isfinite(z) || !std::isfinite(lambda) || lambda < 0.0) {
    Rcpp::stop("Threshold inputs must be finite and thresholds non-negative.");
  }
  const double magnitude = std::abs(z);
  if (!(magnitude > lambda)) {
    return 0.0;
  }
  if (rule == "hard") {
    return z;
  }
  const double sign = std::signbit(z) ? -1.0 : 1.0;
  if (rule == "soft") {
    return sign * (magnitude - lambda);
  }
  if (rule == "scad") {
    if (!(shape > 2.0) || !std::isfinite(shape)) {
      Rcpp::stop("The SCAD parameter must be finite and greater than two.");
    }
    if (magnitude <= 2.0 * lambda) {
      return sign * (magnitude - lambda);
    }
    if (magnitude <= shape * lambda) {
      return ((shape - 1.0) * z - sign * shape * lambda) /
        (shape - 2.0);
    }
    return z;
  }
  if (rule == "adaptive_lasso") {
    if (!(shape >= 0.0) || !std::isfinite(shape)) {
      Rcpp::stop("The adaptive-lasso exponent must be finite and non-negative.");
    }
    const double shrinkage = std::pow(lambda / magnitude, shape + 1.0);
    return z * (1.0 - shrinkage);
  }
  Rcpp::stop("Unknown generalized thresholding rule.");
  return NA_REAL;
}

} // namespace


// [[Rcpp::export]]
arma::mat cpp_ch3_gaussian_threshold_matrix(const arma::mat& covariance,
                                             const arma::mat& threshold,
                                             const std::string& rule,
                                             const double shape,
                                             const bool threshold_diagonal) {
  if (covariance.n_rows != covariance.n_cols ||
      threshold.n_rows != covariance.n_rows ||
      threshold.n_cols != covariance.n_cols) {
    Rcpp::stop("Covariance and threshold matrices must be square and conformable.");
  }
  const arma::uword p = covariance.n_rows;
  arma::mat result(p, p, arma::fill::zeros);
  for (arma::uword i = 0; i < p; ++i) {
    for (arma::uword j = i; j < p; ++j) {
      double value;
      if (i == j && !threshold_diagonal) {
        value = covariance(i, i);
      } else {
        value = threshold_value(covariance(i, j), threshold(i, j),
                                rule, shape);
      }
      result(i, j) = value;
      result(j, i) = value;
    }
  }
  return result;
}


// The Cai--Liu variability matrix.  The covariance argument is deliberately
// separate so that the R interface can expose the finite-sample consequence
// of using either divisor n (the paper) or n-1 (a common book convention).
// [[Rcpp::export]]
arma::mat cpp_ch3_gaussian_product_variability(const arma::mat& residual,
                                                const arma::mat& covariance) {
  const arma::uword n = residual.n_rows;
  const arma::uword p = residual.n_cols;
  if (covariance.n_rows != p || covariance.n_cols != p) {
    Rcpp::stop("`covariance` is not conformable with `residual`.");
  }
  arma::mat theta(p, p, arma::fill::zeros);
  for (arma::uword i = 0; i < p; ++i) {
    for (arma::uword j = i; j < p; ++j) {
      long double total = 0.0L;
      for (arma::uword k = 0; k < n; ++k) {
        const long double product =
          static_cast<long double>(residual(k, i)) * residual(k, j);
        const long double difference = product - covariance(i, j);
        total += difference * difference;
      }
      const double value = static_cast<double>(total / n);
      theta(i, j) = value;
      theta(j, i) = value;
    }
  }
  return theta;
}
