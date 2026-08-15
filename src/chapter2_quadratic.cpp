// Chapter 2 quadratic-form location tests.
//
// The kernels use the smaller of variable-space and sample-space Gram matrices
// for covariance trace functionals. They intentionally reject non-positive
// variance estimates instead of hiding degeneracy with a ridge, an absolute
// value, or a numerical floor.

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp17)]]

#include <RcppArmadillo.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include <string>

namespace {

void check_finite_matrix_quadratic(const arma::mat& x, const char* name) {
  if (!x.is_finite()) {
    Rcpp::stop("`%s` must contain only finite values.", name);
  }
}

void check_finite_vector_quadratic(const arma::vec& x, const char* name) {
  if (!x.is_finite()) {
    Rcpp::stop("`%s` must contain only finite values.", name);
  }
}

void require_positive_finite_quadratic(double value,
                                       const char* quantity,
                                       const char* test_name) {
  if (!std::isfinite(value) || value <= 0.0) {
    Rcpp::stop(
      "%s requires a finite, strictly positive %s; the input is degenerate "
      "for this calibration, and no ridge or variance repair is applied.",
      test_name,
      quantity
    );
  }
}

double squared_crossproduct_trace_quadratic(const arma::mat& x,
                                            double divisor,
                                            std::string& gram_type,
                                            double& gram_dimension) {
  if (x.n_cols <= x.n_rows) {
    const arma::mat gram = x.t() * x / divisor;
    gram_type = "primal";
    gram_dimension = static_cast<double>(x.n_cols);
    return arma::accu(arma::square(gram));
  }

  const arma::mat gram = x * x.t() / divisor;
  gram_type = "dual";
  gram_dimension = static_cast<double>(x.n_rows);
  return arma::accu(arma::square(gram));
}

arma::vec common_column_scales_quadratic(const arma::mat& x,
                                         const arma::vec& mu) {
  arma::vec scales(x.n_cols);
  for (arma::uword j = 0; j < x.n_cols; ++j) {
    double scale = std::abs(mu(j));
    for (arma::uword i = 0; i < x.n_rows; ++i) {
      scale = std::max(scale, std::abs(x(i, j)));
    }
    scales(j) = scale > 0.0 ? scale : 1.0;
  }
  return scales;
}

double common_global_scale_quadratic(const arma::mat& x,
                                     const arma::mat& y) {
  double scale = 0.0;
  for (arma::uword j = 0; j < x.n_cols; ++j) {
    for (arma::uword i = 0; i < x.n_rows; ++i) {
      scale = std::max(scale, std::abs(x(i, j)));
    }
    for (arma::uword i = 0; i < y.n_rows; ++i) {
      scale = std::max(scale, std::abs(y(i, j)));
    }
  }
  return scale > 0.0 ? scale : 1.0;
}

double restore_scaled_quantity_quadratic(double value, double multiplier) {
  if (value == 0.0) {
    return value;
  }
  if (std::isinf(multiplier)) {
    return std::copysign(std::numeric_limits<double>::infinity(), value);
  }
  return value * multiplier;
}

}  // namespace


//' Srivastava--Du one-sample statistic kernel
//'
//' @param x Numeric observation-by-variable matrix.
//' @param mu Numeric null-mean vector.
//' @return Internal list of statistic components.
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List cpp_srivastava_du_one_sample(const arma::mat& x,
                                        const arma::vec& mu) {
  check_finite_matrix_quadratic(x, "x");
  check_finite_vector_quadratic(mu, "mu");

  const arma::uword n = x.n_rows;
  const arma::uword p = x.n_cols;
  if (n < 4) {
    Rcpp::stop("Srivastava-Du requires at least four observations.");
  }
  if (p < 1) {
    Rcpp::stop("Srivastava-Du requires at least one variable.");
  }
  if (mu.n_elem != p) {
    Rcpp::stop("`mu` must have one value per column of `x`.");
  }

  const double N = static_cast<double>(n);
  const double nu = N - 1.0;
  const double p_double = static_cast<double>(p);
  const arma::vec column_scale = common_column_scales_quadratic(x, mu);
  arma::mat x_scaled = x;
  x_scaled.each_row() /= column_scale.t();
  const arma::vec mu_scaled = mu / column_scale;

  const arma::rowvec mean_scaled_row = arma::mean(x_scaled, 0);
  const arma::vec sample_mean_scaled = mean_scaled_row.t();
  const arma::vec difference_scaled = sample_mean_scaled - mu_scaled;
  arma::mat centered_scaled = x_scaled;
  centered_scaled.each_row() -= mean_scaled_row;

  const arma::vec variance_diagonal_scaled =
    arma::sum(arma::square(centered_scaled), 0).t() / nu;
  if (!variance_diagonal_scaled.is_finite() ||
      arma::any(variance_diagonal_scaled <= 0.0)) {
    Rcpp::stop(
      "Srivastava-Du requires every marginal sample variance to be finite "
      "and strictly positive; no ridge is applied."
    );
  }

  const double A = N * arma::dot(
    difference_scaled % difference_scaled,
    1.0 / variance_diagonal_scaled
  );

  arma::mat standardized = centered_scaled;
  standardized.each_row() /= arma::sqrt(variance_diagonal_scaled).t();
  std::string gram_type;
  double gram_dimension = 0.0;
  const double trace_R2 = squared_crossproduct_trace_quadratic(
    standardized,
    nu,
    gram_type,
    gram_dimension
  );
  const double correction = 1.0 + trace_R2 /
    (p_double * std::sqrt(p_double));
  const double null_center = nu * p_double / (nu - 2.0);
  const double variance = 2.0 *
    (trace_R2 - p_double * p_double / nu) * correction;

  require_positive_finite_quadratic(
    variance,
    "normalising variance",
    "Srivastava-Du"
  );
  if (!std::isfinite(A) || !std::isfinite(null_center) ||
      !std::isfinite(trace_R2) || !std::isfinite(correction)) {
    Rcpp::stop(
      "Srivastava-Du components are non-finite at the supplied numerical "
      "scale."
    );
  }

  const double z = (A - null_center) / std::sqrt(variance);
  if (!std::isfinite(z)) {
    Rcpp::stop("Srivastava-Du produced a non-finite standardised statistic.");
  }

  const arma::vec sample_mean = sample_mean_scaled % column_scale;
  const arma::vec difference = difference_scaled % column_scale;
  const arma::vec variance_diagonal = variance_diagonal_scaled %
    arma::square(column_scale);

  return Rcpp::List::create(
    Rcpp::Named("z") = z,
    Rcpp::Named("A") = A,
    Rcpp::Named("sample_mean") = sample_mean,
    Rcpp::Named("difference") = difference,
    Rcpp::Named("variance_diagonal") = variance_diagonal,
    Rcpp::Named("N") = N,
    Rcpp::Named("nu") = nu,
    Rcpp::Named("p") = p_double,
    Rcpp::Named("trace_R2") = trace_R2,
    Rcpp::Named("c") = correction,
    Rcpp::Named("null_center") = null_center,
    Rcpp::Named("variance") = variance,
    Rcpp::Named("column_scale") = column_scale,
    Rcpp::Named("gram_type") = gram_type,
    Rcpp::Named("gram_dimension") = gram_dimension
  );
}


//' Original Bai--Saranadasa two-sample statistic kernel
//'
//' @param x,y Numeric observation-by-variable matrices.
//' @return Internal list of statistic components.
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List cpp_bai_saranadasa_two_sample(const arma::mat& x,
                                         const arma::mat& y) {
  check_finite_matrix_quadratic(x, "x");
  check_finite_matrix_quadratic(y, "y");

  const arma::uword n1 = x.n_rows;
  const arma::uword n2 = y.n_rows;
  const arma::uword p = x.n_cols;
  if (n1 < 2 || n2 < 2) {
    Rcpp::stop(
      "Bai-Saranadasa requires at least two observations in each group."
    );
  }
  if (p < 1 || y.n_cols != p) {
    Rcpp::stop("`x` and `y` must have the same positive number of columns.");
  }

  const double n1_double = static_cast<double>(n1);
  const double n2_double = static_cast<double>(n2);
  const double N = n1_double + n2_double;
  const double nu = N - 2.0;
  const double p_double = static_cast<double>(p);
  const double a = N / (n1_double * n2_double);

  const double global_scale = common_global_scale_quadratic(x, y);
  const arma::mat x_scaled = x / global_scale;
  const arma::mat y_scaled = y / global_scale;

  const arma::rowvec mean_x_scaled_row = arma::mean(x_scaled, 0);
  const arma::rowvec mean_y_scaled_row = arma::mean(y_scaled, 0);
  const arma::vec mean_x_scaled = mean_x_scaled_row.t();
  const arma::vec mean_y_scaled = mean_y_scaled_row.t();
  const arma::vec difference_scaled = mean_x_scaled - mean_y_scaled;

  arma::mat centered_x_scaled = x_scaled;
  arma::mat centered_y_scaled = y_scaled;
  centered_x_scaled.each_row() -= mean_x_scaled_row;
  centered_y_scaled.each_row() -= mean_y_scaled_row;
  const arma::mat centered_scaled = arma::join_cols(
    centered_x_scaled,
    centered_y_scaled
  );

  const double trace_Sp_scaled =
    arma::accu(arma::square(centered_scaled)) / nu;
  std::string gram_type;
  double gram_dimension = 0.0;
  const double trace_Sp2_scaled = squared_crossproduct_trace_quadratic(
    centered_scaled,
    nu,
    gram_type,
    gram_dimension
  );
  const double B_scaled = trace_Sp2_scaled -
    trace_Sp_scaled * trace_Sp_scaled / nu;
  require_positive_finite_quadratic(
    B_scaled,
    "bias-correction term B",
    "Bai-Saranadasa"
  );

  const double trace_sigma2_hat_scaled =
    nu * nu / (N * (N - 3.0)) * B_scaled;
  const double variance_scaled = 2.0 * a * a * (N - 1.0) / nu *
    trace_sigma2_hat_scaled;
  const double xi_scaled = arma::dot(difference_scaled, difference_scaled) -
    a * trace_Sp_scaled;
  require_positive_finite_quadratic(
    trace_sigma2_hat_scaled,
    "estimate of tr(Sigma^2)",
    "Bai-Saranadasa"
  );
  require_positive_finite_quadratic(
    variance_scaled,
    "normalising variance",
    "Bai-Saranadasa"
  );
  if (!std::isfinite(trace_Sp_scaled) ||
      !std::isfinite(trace_Sp2_scaled) || !std::isfinite(xi_scaled)) {
    Rcpp::stop(
      "Bai-Saranadasa components are non-finite at the supplied numerical "
      "scale."
    );
  }

  const double z = xi_scaled / std::sqrt(variance_scaled);
  if (!std::isfinite(z)) {
    Rcpp::stop(
      "Bai-Saranadasa produced a non-finite standardised statistic."
    );
  }

  const double scale_squared = global_scale * global_scale;
  const double scale_fourth = scale_squared * scale_squared;
  const arma::vec mean_x = mean_x_scaled * global_scale;
  const arma::vec mean_y = mean_y_scaled * global_scale;
  const arma::vec difference = difference_scaled * global_scale;
  const double trace_Sp = restore_scaled_quantity_quadratic(
    trace_Sp_scaled,
    scale_squared
  );
  const double trace_Sp2 = restore_scaled_quantity_quadratic(
    trace_Sp2_scaled,
    scale_fourth
  );
  const double B = restore_scaled_quantity_quadratic(B_scaled, scale_fourth);
  const double trace_sigma2_hat = restore_scaled_quantity_quadratic(
    trace_sigma2_hat_scaled,
    scale_fourth
  );
  const double xi = restore_scaled_quantity_quadratic(
    xi_scaled,
    scale_squared
  );
  const double variance = restore_scaled_quantity_quadratic(
    variance_scaled,
    scale_fourth
  );

  return Rcpp::List::create(
    Rcpp::Named("z") = z,
    Rcpp::Named("xi") = xi,
    Rcpp::Named("mean_x") = mean_x,
    Rcpp::Named("mean_y") = mean_y,
    Rcpp::Named("difference") = difference,
    Rcpp::Named("n1") = n1_double,
    Rcpp::Named("n2") = n2_double,
    Rcpp::Named("N") = N,
    Rcpp::Named("nu") = nu,
    Rcpp::Named("p") = p_double,
    Rcpp::Named("a") = a,
    Rcpp::Named("trace_Sp") = trace_Sp,
    Rcpp::Named("trace_Sp2") = trace_Sp2,
    Rcpp::Named("B") = B,
    Rcpp::Named("trace_sigma2_hat") = trace_sigma2_hat,
    Rcpp::Named("variance") = variance,
    Rcpp::Named("xi_scaled") = xi_scaled,
    Rcpp::Named("variance_scaled") = variance_scaled,
    Rcpp::Named("global_scale") = global_scale,
    Rcpp::Named("gram_type") = gram_type,
    Rcpp::Named("gram_dimension") = gram_dimension
  );
}
