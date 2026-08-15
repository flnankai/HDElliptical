// Zhang--Zhou--Guo normal-reference one-sample mean test.
//
// The paper's feasible statistic is the centred U-statistic
//   T = n ||xbar||^2 - tr(S)
//     = {2 / (n - 1)} sum_{i < j} x_i' x_j,
// not the motivating Gaussian oracle n ||xbar||^2.  This kernel computes T,
// the exact Gaussian/Wishart unbiased estimators of tr(Sigma^2) and
// tr(Sigma^3), and the three-cumulant shifted-scaled chi-square parameters.

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp17)]]

#include <RcppArmadillo.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <string>

namespace {

void zzg22nr_neumaier_add(long double value,
                          long double& total,
                          long double& correction) {
  const long double updated = total + value;
  if (std::abs(total) >= std::abs(value)) {
    correction += (total - updated) + value;
  } else {
    correction += (value - updated) + total;
  }
  total = updated;
}

double zzg22nr_checked_finite(long double value,
                              const std::string& label) {
  const long double maximum =
    static_cast<long double>(std::numeric_limits<double>::max());
  if (!std::isfinite(value) || std::abs(value) > maximum) {
    Rcpp::stop("Zhang-Zhou-Guo produced a non-finite %s.", label.c_str());
  }
  return static_cast<double>(value);
}

double zzg22nr_checked_positive(long double value,
                                const std::string& label) {
  const double answer = zzg22nr_checked_finite(value, label);
  if (!(answer > 0.0)) {
    Rcpp::stop(
      "Zhang-Zhou-Guo requires a strictly positive %s; no ridge, "
      "absolute-value repair, or numerical floor is applied.",
      label.c_str()
    );
  }
  return answer;
}

arma::rowvec zzg22nr_stable_column_mean(const arma::mat& x) {
  arma::rowvec answer(x.n_cols, arma::fill::zeros);
  for (arma::uword j = 0; j < x.n_cols; ++j) {
    long double total = 0.0L;
    long double correction = 0.0L;
    for (arma::uword i = 0; i < x.n_rows; ++i) {
      zzg22nr_neumaier_add(
        static_cast<long double>(x(i, j)), total, correction
      );
    }
    answer(j) = static_cast<double>(
      (total + correction) / static_cast<long double>(x.n_rows)
    );
  }
  return answer;
}

long double zzg22nr_trace(const arma::mat& matrix) {
  long double total = 0.0L;
  long double correction = 0.0L;
  for (arma::uword i = 0; i < matrix.n_rows; ++i) {
    zzg22nr_neumaier_add(
      static_cast<long double>(matrix(i, i)), total, correction
    );
  }
  return total + correction;
}

long double zzg22nr_sum_squares(const arma::mat& matrix) {
  long double total = 0.0L;
  long double correction = 0.0L;
  for (arma::uword j = 0; j < matrix.n_cols; ++j) {
    for (arma::uword i = 0; i < matrix.n_rows; ++i) {
      const long double value = static_cast<long double>(matrix(i, j));
      zzg22nr_neumaier_add(value * value, total, correction);
    }
  }
  return total + correction;
}

long double zzg22nr_trace_cube(const arma::mat& matrix) {
  const arma::mat square = matrix * matrix;
  long double total = 0.0L;
  long double correction = 0.0L;
  for (arma::uword j = 0; j < matrix.n_cols; ++j) {
    for (arma::uword i = 0; i < matrix.n_rows; ++i) {
      zzg22nr_neumaier_add(
        static_cast<long double>(square(i, j)) *
          static_cast<long double>(matrix(j, i)),
        total,
        correction
      );
    }
  }
  return total + correction;
}

long double zzg22nr_u_statistic(const arma::mat& x) {
  // Coordinatewise identity for the ordered off-diagonal sum.  Long-double
  // accumulation avoids first forming an n-by-n raw Gram matrix and is more
  // stable than subtracting two large double-precision quadratic forms.
  long double total = 0.0L;
  long double correction = 0.0L;
  for (arma::uword j = 0; j < x.n_cols; ++j) {
    long double sum = 0.0L;
    long double sum_correction = 0.0L;
    long double square_sum = 0.0L;
    long double square_correction = 0.0L;
    for (arma::uword i = 0; i < x.n_rows; ++i) {
      const long double value = static_cast<long double>(x(i, j));
      zzg22nr_neumaier_add(value, sum, sum_correction);
      zzg22nr_neumaier_add(
        value * value, square_sum, square_correction
      );
    }
    const long double column_sum = sum + sum_correction;
    const long double contribution =
      column_sum * column_sum - (square_sum + square_correction);
    zzg22nr_neumaier_add(contribution, total, correction);
  }
  return (total + correction) /
    static_cast<long double>(x.n_rows - 1U);
}

}  // namespace


//' Zhang--Zhou--Guo one-sample normal-reference kernel
//'
//' @param residual Numeric matrix of null-centred observations.  A common
//'   positive scaling of this matrix leaves the reported calibration
//'   unchanged.
//'
//' @return Internal list of the centred U-statistic, unbiased trace
//'   estimates, cumulants, and shifted-scaled chi-square parameters.
//'
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List cpp_zhang_zhou_guo_one_sample(const arma::mat& residual) {
  const arma::uword n = residual.n_rows;
  const arma::uword p = residual.n_cols;
  if (n < 4U || p < 1U) {
    Rcpp::stop(
      "Zhang-Zhou-Guo requires at least four observations and one variable."
    );
  }
  if (!residual.is_finite()) {
    Rcpp::stop(
      "Zhang-Zhou-Guo requires finite null-centred observations."
    );
  }

  double internal_scale = 0.0;
  for (arma::uword j = 0; j < p; ++j) {
    for (arma::uword i = 0; i < n; ++i) {
      internal_scale = std::max(
        internal_scale, std::abs(residual(i, j))
      );
    }
  }
  if (!(internal_scale > 0.0) || !std::isfinite(internal_scale)) {
    Rcpp::stop(
      "Zhang-Zhou-Guo requires nonzero null-centred observations and "
      "positive sample variation."
    );
  }

  const arma::mat x = residual / internal_scale;
  const arma::rowvec mean = zzg22nr_stable_column_mean(x);
  arma::mat centred = x;
  centred.each_row() -= mean;

  const long double v = static_cast<long double>(n - 1U);
  const bool primal = p <= n;
  arma::mat covariance_spectrum;
  if (primal) {
    covariance_spectrum = (centred.t() * centred) /
      static_cast<double>(v);
  } else {
    covariance_spectrum = (centred * centred.t()) /
      static_cast<double>(v);
  }

  const long double trace_s = zzg22nr_trace(covariance_spectrum);
  const long double trace_s2 = zzg22nr_sum_squares(covariance_spectrum);
  const long double trace_s3 = zzg22nr_trace_cube(covariance_spectrum);

  long double mean_square = 0.0L;
  long double mean_square_correction = 0.0L;
  for (arma::uword j = 0; j < p; ++j) {
    const long double value = static_cast<long double>(mean(j));
    zzg22nr_neumaier_add(
      value * value, mean_square, mean_square_correction
    );
  }
  mean_square += mean_square_correction;
  const long double oracle_q = static_cast<long double>(n) * mean_square;
  const long double statistic_by_subtraction = oracle_q - trace_s;
  const long double statistic = zzg22nr_u_statistic(x);

  const long double a1_hat = trace_s;
  const long double a2_bracket = trace_s2 - trace_s * trace_s / v;
  const long double a2_factor =
    v * v / ((v - 1.0L) * (v + 2.0L));
  const long double a2_hat = a2_factor * a2_bracket;

  // Appendix (A.33) of the primary paper has (v + 4), not (v + 3).
  const long double a3_bracket =
    trace_s3 - 3.0L * trace_s * trace_s2 / v +
    2.0L * trace_s * trace_s * trace_s / (v * v);
  const long double a3_factor =
    v * v * v * v /
    ((v - 1.0L) * (v + 4.0L) * (v * v - 4.0L));
  const long double a3_hat = a3_factor * a3_bracket;

  const double a2_checked = zzg22nr_checked_positive(
    a2_hat, "unbiased trace(Sigma^2) estimate"
  );
  const double a3_checked = zzg22nr_checked_positive(
    a3_hat, "unbiased trace(Sigma^3) estimate"
  );

  const long double n_ld = static_cast<long double>(n);
  const long double n_minus_two = n_ld - 2.0L;
  const long double kappa2 = 2.0L * n_ld / v * a2_hat;
  const long double kappa3 =
    8.0L * n_ld * n_minus_two / (v * v) * a3_hat;
  const long double beta0 =
    -n_ld / n_minus_two * a2_hat * a2_hat / a3_hat;
  const long double beta1 =
    n_minus_two / v * a3_hat / a2_hat;
  const long double degrees =
    n_ld * v / (n_minus_two * n_minus_two) *
    a2_hat * a2_hat * a2_hat / (a3_hat * a3_hat);

  const double beta1_checked = zzg22nr_checked_positive(
    beta1, "shifted-chi-square scale beta1"
  );
  const double degrees_checked = zzg22nr_checked_positive(
    degrees, "matched chi-square degrees of freedom"
  );
  const long double reference_argument = (statistic - beta0) / beta1;

  return Rcpp::List::create(
    Rcpp::Named("T_scaled") = zzg22nr_checked_finite(
      statistic, "centred U-statistic"
    ),
    Rcpp::Named("T_by_oracle_minus_trace_scaled") =
      zzg22nr_checked_finite(
        statistic_by_subtraction,
        "oracle-minus-trace identity statistic"
      ),
    Rcpp::Named("T_identity_error_scaled") = zzg22nr_checked_finite(
      std::abs(statistic - statistic_by_subtraction),
      "U-statistic identity error"
    ),
    Rcpp::Named("oracle_Q_scaled") = zzg22nr_checked_finite(
      oracle_q, "Gaussian oracle quadratic form"
    ),
    Rcpp::Named("mean_residual_scaled") = mean,
    Rcpp::Named("trace_S_scaled") = zzg22nr_checked_finite(
      trace_s, "trace(S)"
    ),
    Rcpp::Named("trace_S2_scaled") = zzg22nr_checked_finite(
      trace_s2, "trace(S^2)"
    ),
    Rcpp::Named("trace_S3_scaled") = zzg22nr_checked_finite(
      trace_s3, "trace(S^3)"
    ),
    Rcpp::Named("a1_hat_scaled") = zzg22nr_checked_finite(
      a1_hat, "trace(Sigma) estimate"
    ),
    Rcpp::Named("a2_bracket_scaled") = zzg22nr_checked_finite(
      a2_bracket, "trace(Sigma^2) correction bracket"
    ),
    Rcpp::Named("a2_factor") = static_cast<double>(a2_factor),
    Rcpp::Named("a2_hat_scaled") = a2_checked,
    Rcpp::Named("a3_bracket_scaled") = zzg22nr_checked_finite(
      a3_bracket, "trace(Sigma^3) correction bracket"
    ),
    Rcpp::Named("a3_factor") = static_cast<double>(a3_factor),
    Rcpp::Named("a3_hat_scaled") = a3_checked,
    Rcpp::Named("kappa1_scaled") = 0.0,
    Rcpp::Named("kappa2_scaled") = zzg22nr_checked_positive(
      kappa2, "second reference cumulant"
    ),
    Rcpp::Named("kappa3_scaled") = zzg22nr_checked_positive(
      kappa3, "third reference cumulant"
    ),
    Rcpp::Named("beta0_scaled") = zzg22nr_checked_finite(
      beta0, "shift parameter beta0"
    ),
    Rcpp::Named("beta1_scaled") = beta1_checked,
    Rcpp::Named("df") = degrees_checked,
    Rcpp::Named("reference_argument") = zzg22nr_checked_finite(
      reference_argument, "chi-square reference argument"
    ),
    Rcpp::Named("internal_scale") = internal_scale,
    Rcpp::Named("trace_computation") = primal ? "primal" : "dual",
    Rcpp::Named("constructs_p_by_p_matrix") = primal,
    Rcpp::Named("n") = static_cast<double>(n),
    Rcpp::Named("p") = static_cast<double>(p)
  );
}
