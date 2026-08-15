// Chapter 6: classical PCA/CCA and robust factor-analysis kernels.
//
// These kernels only evaluate stated cross-products and least-squares
// identities.  Matrix rank, positive definiteness, eigengaps, and all method
// contracts are certified in R.  No ridge, pseudoinverse, eigenvalue floor,
// or other numerical repair is performed here.

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp17)]]

#include <RcppArmadillo.h>

#include <cmath>
#include <limits>

namespace {

void ch6cf_check_data(const arma::mat& x, const char* name,
                      const arma::uword minimum_rows = 1) {
  if (x.n_rows < minimum_rows || x.n_cols == 0 || !x.is_finite()) {
    Rcpp::stop("`%s` must be a finite matrix with at least %d row(s) and one column.",
               name, static_cast<int>(minimum_rows));
  }
}

void ch6cf_check_center(const arma::rowvec& center, const arma::uword p,
                        const char* name) {
  if (center.n_elem != p || !center.is_finite()) {
    Rcpp::stop("`%s` must be finite and have one value per column.", name);
  }
}

void ch6cf_check_finite(const arma::mat& x, const char* description) {
  if (!x.is_finite()) {
    Rcpp::stop("%s is outside the finite double range.", description);
  }
}

}  // namespace


// [[Rcpp::export]]
Rcpp::List cpp_ch6cf_center_moments(const arma::mat& x,
                                    const arma::rowvec& center,
                                    const double divisor) {
  ch6cf_check_data(x, "x", 2);
  ch6cf_check_center(center, x.n_cols, "center");
  if (!std::isfinite(divisor) || divisor <= 0.0) {
    Rcpp::stop("`divisor` must be finite and strictly positive.");
  }

  const arma::mat centered = x.each_row() - center;
  ch6cf_check_finite(centered, "The centered data matrix");
  const arma::mat covariance = centered.t() * centered / divisor;
  ch6cf_check_finite(covariance, "The centered second-moment matrix");

  return Rcpp::List::create(
    Rcpp::Named("centered") = centered,
    Rcpp::Named("covariance") = covariance,
    Rcpp::Named("divisor") = divisor
  );
}


// [[Rcpp::export]]
Rcpp::List cpp_ch6cf_cca_moments(const arma::mat& x,
                                 const arma::mat& y,
                                 const arma::rowvec& center_x,
                                 const arma::rowvec& center_y,
                                 const double divisor) {
  ch6cf_check_data(x, "x", 2);
  ch6cf_check_data(y, "y", 2);
  if (x.n_rows != y.n_rows) {
    Rcpp::stop("`x` and `y` must have the same number of rows.");
  }
  ch6cf_check_center(center_x, x.n_cols, "center_x");
  ch6cf_check_center(center_y, y.n_cols, "center_y");
  if (!std::isfinite(divisor) || divisor <= 0.0) {
    Rcpp::stop("`divisor` must be finite and strictly positive.");
  }

  const arma::mat centered_x = x.each_row() - center_x;
  const arma::mat centered_y = y.each_row() - center_y;
  ch6cf_check_finite(centered_x, "The centered `x` matrix");
  ch6cf_check_finite(centered_y, "The centered `y` matrix");
  const arma::mat covariance_xx = centered_x.t() * centered_x / divisor;
  const arma::mat covariance_yy = centered_y.t() * centered_y / divisor;
  const arma::mat covariance_xy = centered_x.t() * centered_y / divisor;
  ch6cf_check_finite(covariance_xx, "The `x` covariance block");
  ch6cf_check_finite(covariance_yy, "The `y` covariance block");
  ch6cf_check_finite(covariance_xy, "The cross-covariance block");

  return Rcpp::List::create(
    Rcpp::Named("centered_x") = centered_x,
    Rcpp::Named("centered_y") = centered_y,
    Rcpp::Named("covariance_xx") = covariance_xx,
    Rcpp::Named("covariance_xy") = covariance_xy,
    Rcpp::Named("covariance_yy") = covariance_yy,
    Rcpp::Named("divisor") = divisor
  );
}


// [[Rcpp::export]]
Rcpp::List cpp_ch6cf_rts_components(const arma::mat& x,
                                    const arma::rowvec& center,
                                    const arma::mat& loadings) {
  ch6cf_check_data(x, "x", 2);
  ch6cf_check_center(center, x.n_cols, "center");
  if (loadings.n_rows != x.n_cols || loadings.n_cols == 0 ||
      !loadings.is_finite()) {
    Rcpp::stop("`loadings` must be a finite p by K matrix with K positive.");
  }

  const arma::mat centered = x.each_row() - center;
  ch6cf_check_finite(centered, "The centered factor data");
  const double p = static_cast<double>(x.n_cols);
  const arma::mat factor_scores = centered * loadings / p;
  const arma::mat common_component = factor_scores * loadings.t();
  const arma::mat residuals = centered - common_component;
  const arma::mat normal_equation = residuals * loadings;
  const arma::mat loading_gram = loadings.t() * loadings / p;
  ch6cf_check_finite(factor_scores, "The RTS factor-score matrix");
  ch6cf_check_finite(common_component, "The RTS common-component matrix");
  ch6cf_check_finite(residuals, "The RTS residual matrix");

  const double normal_equation_residual = normal_equation.n_elem == 0 ? 0.0 :
    arma::abs(normal_equation).max();

  return Rcpp::List::create(
    Rcpp::Named("centered") = centered,
    Rcpp::Named("factor_scores") = factor_scores,
    Rcpp::Named("common_component") = common_component,
    Rcpp::Named("residuals") = residuals,
    Rcpp::Named("loading_gram_scaled") = loading_gram,
    Rcpp::Named("normal_equation_residual") = normal_equation_residual
  );
}
