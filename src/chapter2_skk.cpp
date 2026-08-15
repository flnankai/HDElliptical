// Srivastava--Katayama--Kano two-sample high-dimensional mean test.
//
// The trace functionals are evaluated in variable space when p is smaller
// than n1 + n2 and in observation space otherwise.  In particular, the
// high-dimensional branch never materialises a p by p covariance matrix.

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp17)]]

#include <RcppArmadillo.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <string>

namespace {

void skk_require_finite_matrix(const arma::mat& x, const char* name) {
  if (!x.is_finite()) {
    Rcpp::stop("`%s` must contain only finite values.", name);
  }
}

void skk_require_finite_component(double value,
                                  const char* quantity) {
  if (!std::isfinite(value)) {
    Rcpp::stop(
      "Srivastava-Katayama-Kano requires a finite %s; no absolute-value "
      "or numerical-floor repair is applied.",
      quantity
    );
  }
}

void skk_require_positive_finite(double value, const char* quantity) {
  if (!std::isfinite(value) || value <= 0.0) {
    Rcpp::stop(
      "Srivastava-Katayama-Kano requires a finite, strictly positive %s; "
      "the input is degenerate for this calibration, and no ridge or "
      "variance repair is applied.",
      quantity
    );
  }
}

arma::vec skk_common_column_scales(const arma::mat& x,
                                   const arma::mat& y) {
  arma::vec scales(x.n_cols, arma::fill::ones);
  for (arma::uword j = 0; j < x.n_cols; ++j) {
    double scale = 0.0;
    for (arma::uword i = 0; i < x.n_rows; ++i) {
      scale = std::max(scale, std::abs(x(i, j)));
    }
    for (arma::uword i = 0; i < y.n_rows; ++i) {
      scale = std::max(scale, std::abs(y(i, j)));
    }
    if (scale > 0.0) {
      scales(j) = scale;
    }
  }
  return scales;
}

long double skk_sum_squares(const arma::mat& x) {
  long double total = 0.0L;
  for (arma::uword j = 0; j < x.n_cols; ++j) {
    for (arma::uword i = 0; i < x.n_rows; ++i) {
      const long double value = static_cast<long double>(x(i, j));
      total += value * value;
    }
  }
  return total;
}

long double skk_frobenius_inner(const arma::mat& x,
                                const arma::mat& y) {
  long double total = 0.0L;
  for (arma::uword j = 0; j < x.n_cols; ++j) {
    for (arma::uword i = 0; i < x.n_rows; ++i) {
      total += static_cast<long double>(x(i, j)) *
        static_cast<long double>(y(i, j));
    }
  }
  return total;
}

long double skk_trace_standardised_covariance(const arma::mat& z,
                                              long double divisor) {
  long double total = 0.0L;
  for (arma::uword j = 0; j < z.n_cols; ++j) {
    for (arma::uword i = 0; i < z.n_rows; ++i) {
      const long double value = static_cast<long double>(z(i, j));
      total += value * value;
    }
  }
  return total / divisor;
}

double skk_restore_variance(double value, double scale) {
  if (value == 0.0) {
    return 0.0;
  }
  const long double restored = static_cast<long double>(value) *
    static_cast<long double>(scale) * static_cast<long double>(scale);
  if (restored > static_cast<long double>(
        std::numeric_limits<double>::max())) {
    return std::numeric_limits<double>::infinity();
  }
  return static_cast<double>(restored);
}

arma::vec skk_restore_variance_vector(const arma::vec& value,
                                      const arma::vec& scale) {
  arma::vec restored(value.n_elem);
  for (arma::uword j = 0; j < value.n_elem; ++j) {
    restored(j) = skk_restore_variance(value(j), scale(j));
  }
  return restored;
}

}  // namespace


//' Srivastava--Katayama--Kano statistic kernel
//'
//' @param x,y Numeric observation-by-variable matrices.
//' @return Internal list of statistic components.
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List cpp_skk_two_sample(const arma::mat& x,
                              const arma::mat& y) {
  skk_require_finite_matrix(x, "x");
  skk_require_finite_matrix(y, "y");

  const arma::uword n1 = x.n_rows;
  const arma::uword n2 = y.n_rows;
  const arma::uword p = x.n_cols;
  if (n1 < 2 || n2 < 2) {
    Rcpp::stop(
      "Srivastava-Katayama-Kano requires at least two observations in "
      "each group."
    );
  }
  if (p < 1 || y.n_cols != p) {
    Rcpp::stop("`x` and `y` must have the same positive number of columns.");
  }

  const double n1_double = static_cast<double>(n1);
  const double n2_double = static_cast<double>(n2);
  const double nu1 = n1_double - 1.0;
  const double nu2 = n2_double - 1.0;
  const double p_double = static_cast<double>(p);

  // A common scale is used within each variable and across both groups.  This
  // is an algebraically neutral diagonal change of units, including when a
  // scale factor is negative in the user's data transformation.
  const arma::vec column_scale = skk_common_column_scales(x, y);
  arma::mat x_scaled = x;
  arma::mat y_scaled = y;
  x_scaled.each_row() /= column_scale.t();
  y_scaled.each_row() /= column_scale.t();

  const arma::rowvec mean_x_scaled_row = arma::mean(x_scaled, 0);
  const arma::rowvec mean_y_scaled_row = arma::mean(y_scaled, 0);
  const arma::vec mean_x_scaled = mean_x_scaled_row.t();
  const arma::vec mean_y_scaled = mean_y_scaled_row.t();
  const arma::vec difference_scaled = mean_x_scaled - mean_y_scaled;

  arma::mat centered_x = x_scaled;
  arma::mat centered_y = y_scaled;
  centered_x.each_row() -= mean_x_scaled_row;
  centered_y.each_row() -= mean_y_scaled_row;

  const arma::vec variance1_scaled =
    arma::sum(arma::square(centered_x), 0).t() / nu1;
  const arma::vec variance2_scaled =
    arma::sum(arma::square(centered_y), 0).t() / nu2;
  const arma::vec d_hat_scaled = variance1_scaled / n1_double +
    variance2_scaled / n2_double;
  if (!d_hat_scaled.is_finite() || arma::any(d_hat_scaled <= 0.0)) {
    Rcpp::stop(
      "Srivastava-Katayama-Kano requires every combined marginal variance "
      "diag(S1) / n1 + diag(S2) / n2 to be finite and strictly positive; "
      "no ridge is applied."
    );
  }

  const arma::vec inverse_sqrt_d = 1.0 / arma::sqrt(d_hat_scaled);
  arma::mat z1 = centered_x;
  arma::mat z2 = centered_y;
  z1.each_row() %= inverse_sqrt_d.t();
  z2.each_row() %= inverse_sqrt_d.t();
  if (!z1.is_finite() || !z2.is_finite()) {
    Rcpp::stop(
      "Srivastava-Katayama-Kano diagonal standardisation produced "
      "non-finite values."
    );
  }

  const long double trace_a1 = skk_trace_standardised_covariance(
    z1, static_cast<long double>(nu1)
  );
  const long double trace_a2 = skk_trace_standardised_covariance(
    z2, static_cast<long double>(nu2)
  );
  long double trace_a1_squared = 0.0L;
  long double trace_a2_squared = 0.0L;
  long double trace_a1_a2 = 0.0L;
  std::string gram_type;
  double gram_dimension = 0.0;

  const arma::uword total = n1 + n2;
  if (p <= total) {
    const arma::mat a1 = z1.t() * z1 / nu1;
    const arma::mat a2 = z2.t() * z2 / nu2;
    trace_a1_squared = skk_sum_squares(a1);
    trace_a2_squared = skk_sum_squares(a2);
    trace_a1_a2 = skk_frobenius_inner(a1, a2);
    gram_type = "primal";
    gram_dimension = p_double;
  } else {
    arma::mat standardised(total, p);
    standardised.rows(0, n1 - 1) = z1 / std::sqrt(nu1);
    standardised.rows(n1, total - 1) = z2 / std::sqrt(nu2);
    const arma::mat gram = standardised * standardised.t();
    trace_a1_squared = skk_sum_squares(
      gram.submat(0, 0, n1 - 1, n1 - 1)
    );
    trace_a2_squared = skk_sum_squares(
      gram.submat(n1, n1, total - 1, total - 1)
    );
    trace_a1_a2 = skk_sum_squares(
      gram.submat(0, n1, n1 - 1, total - 1)
    );
    gram_type = "dual";
    gram_dimension = static_cast<double>(total);
  }

  const long double p_ld = static_cast<long double>(p_double);
  const long double f1_ld = (
    trace_a1_squared - trace_a1 * trace_a1 /
      static_cast<long double>(nu1)
  ) / p_ld;
  const long double f2_ld = (
    trace_a2_squared - trace_a2 * trace_a2 /
      static_cast<long double>(nu2)
  ) / p_ld;
  const long double g_ld = trace_a1_a2 / p_ld;
  const double f1 = static_cast<double>(f1_ld);
  const double f2 = static_cast<double>(f2_ld);
  const double g = static_cast<double>(g_ld);
  // F1 and F2 are unbiased finite-sample trace corrections and can be
  // negative even when their combined variance estimate is positive. Reject
  // only non-finite components here; the actual calibration contract is the
  // strict positivity of variance_q below.
  skk_require_finite_component(f1, "F1 trace correction");
  skk_require_finite_component(f2, "F2 trace correction");
  skk_require_finite_component(g, "G cross-trace term");

  const long double variance_q_ld =
    2.0L * f1_ld /
      (static_cast<long double>(n1_double) * n1_double) +
    2.0L * f2_ld /
      (static_cast<long double>(n2_double) * n2_double) +
    4.0L * g_ld /
      (static_cast<long double>(n1_double) * n2_double);
  const long double trace_r_hat_squared_ld =
    trace_a1_squared /
      (static_cast<long double>(n1_double) * n1_double) +
    trace_a2_squared /
      (static_cast<long double>(n2_double) * n2_double) +
    2.0L * trace_a1_a2 /
      (static_cast<long double>(n1_double) * n2_double);
  const long double c_hat_ld = 1.0L + trace_r_hat_squared_ld /
    (p_ld * std::sqrt(p_ld));

  long double q_raw_ld = 0.0L;
  for (arma::uword j = 0; j < p; ++j) {
    const long double difference =
      static_cast<long double>(difference_scaled(j));
    q_raw_ld += difference * difference /
      static_cast<long double>(d_hat_scaled(j));
  }
  const long double q_scaled_ld = (q_raw_ld - p_ld) / std::sqrt(p_ld);
  const long double denominator_variance_ld = variance_q_ld * c_hat_ld;

  const double q_raw = static_cast<double>(q_raw_ld);
  const double q_scaled = static_cast<double>(q_scaled_ld);
  const double variance_q = static_cast<double>(variance_q_ld);
  const double trace_r_hat_squared =
    static_cast<double>(trace_r_hat_squared_ld);
  const double c_hat = static_cast<double>(c_hat_ld);
  const double denominator_variance =
    static_cast<double>(denominator_variance_ld);
  skk_require_positive_finite(variance_q, "variance estimate for q");
  skk_require_positive_finite(trace_r_hat_squared, "trace of R-hat squared");
  skk_require_positive_finite(c_hat, "finite-sample correction c-hat");
  skk_require_positive_finite(
    denominator_variance, "normalising variance"
  );
  if (!std::isfinite(q_raw) || !std::isfinite(q_scaled)) {
    Rcpp::stop(
      "Srivastava-Katayama-Kano produced a non-finite diagonal quadratic "
      "form."
    );
  }

  const double z = q_scaled / std::sqrt(denominator_variance);
  if (!std::isfinite(z)) {
    Rcpp::stop(
      "Srivastava-Katayama-Kano produced a non-finite standardised "
      "statistic."
    );
  }

  const arma::vec mean_x = mean_x_scaled % column_scale;
  const arma::vec mean_y = mean_y_scaled % column_scale;
  const arma::vec difference = difference_scaled % column_scale;
  const arma::vec variance1 = skk_restore_variance_vector(
    variance1_scaled, column_scale
  );
  const arma::vec variance2 = skk_restore_variance_vector(
    variance2_scaled, column_scale
  );
  const arma::vec d_hat = skk_restore_variance_vector(
    d_hat_scaled, column_scale
  );

  return Rcpp::List::create(
    Rcpp::Named("z") = z,
    Rcpp::Named("Q_raw") = q_raw,
    Rcpp::Named("q_scaled") = q_scaled,
    Rcpp::Named("F1") = f1,
    Rcpp::Named("F2") = f2,
    Rcpp::Named("G") = g,
    Rcpp::Named("variance_q") = variance_q,
    Rcpp::Named("trace_R_hat2") = trace_r_hat_squared,
    Rcpp::Named("c_hat") = c_hat,
    Rcpp::Named("denominator_variance") = denominator_variance,
    Rcpp::Named("trace_A1") = static_cast<double>(trace_a1),
    Rcpp::Named("trace_A2") = static_cast<double>(trace_a2),
    Rcpp::Named("trace_A1_squared") =
      static_cast<double>(trace_a1_squared),
    Rcpp::Named("trace_A2_squared") =
      static_cast<double>(trace_a2_squared),
    Rcpp::Named("trace_A1_A2") = static_cast<double>(trace_a1_a2),
    Rcpp::Named("mean_x") = mean_x,
    Rcpp::Named("mean_y") = mean_y,
    Rcpp::Named("difference") = difference,
    Rcpp::Named("variance1_diagonal") = variance1,
    Rcpp::Named("variance2_diagonal") = variance2,
    Rcpp::Named("D_hat_diagonal") = d_hat,
    Rcpp::Named("n1") = n1_double,
    Rcpp::Named("n2") = n2_double,
    Rcpp::Named("p") = p_double,
    Rcpp::Named("column_scale") = column_scale,
    Rcpp::Named("gram_type") = gram_type,
    Rcpp::Named("gram_dimension") = gram_dimension
  );
}
