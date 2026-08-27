// Copyright (c) 2026 Long Feng
// SPDX-License-Identifier: MIT

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp17)]]

#include <RcppArmadillo.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <string>

namespace ch4afc {

inline double checked_double(const long double value,
                             const std::string& label) {
  if (!std::isfinite(value) ||
      std::abs(value) >
        static_cast<long double>(std::numeric_limits<double>::max())) {
    Rcpp::stop("The %s is not representable as a finite double.", label);
  }
  return static_cast<double>(value);
}

struct SignRow {
  arma::rowvec value;
  bool zero;
};

inline SignRow stable_sign(const arma::rowvec& x) {
  const arma::uword p = x.n_elem;
  double maximum = 0.0;
  for (arma::uword j = 0; j < p; ++j) {
    maximum = std::max(maximum, std::abs(x(j)));
  }
  SignRow answer{arma::rowvec(p, arma::fill::zeros), true};
  if (maximum == 0.0) {
    return answer;
  }
  const arma::rowvec scaled = x / maximum;
  const double norm = arma::norm(scaled, 2);
  if (!std::isfinite(norm) || !(norm > 0.0)) {
    Rcpp::stop("A spatial-sign norm is not finite and strictly positive.");
  }
  answer.value = scaled / norm;
  answer.zero = false;
  return answer;
}

inline arma::mat signs(const arma::mat& x, arma::uword& zero_rows) {
  arma::mat answer(x.n_rows, x.n_cols, arma::fill::zeros);
  zero_rows = 0U;
  for (arma::uword i = 0; i < x.n_rows; ++i) {
    const SignRow current = stable_sign(x.row(i));
    answer.row(i) = current.value;
    zero_rows += static_cast<arma::uword>(current.zero);
  }
  return answer;
}

inline long double stable_dot(const arma::rowvec& x,
                              const arma::rowvec& y) {
  long double total = 0.0L;
  long double correction = 0.0L;
  for (arma::uword j = 0; j < x.n_elem; ++j) {
    const long double term = static_cast<long double>(x(j)) *
      static_cast<long double>(y(j));
    const long double updated = total + term;
    if (std::abs(total) >= std::abs(term)) {
      correction += (total - updated) + term;
    } else {
      correction += (term - updated) + total;
    }
    total = updated;
  }
  return total + correction;
}

}  // namespace ch4afc


//' Strict least-squares projection for conditional-alpha procedures
//'
//' @param y Finite observation-by-asset matrix on a common numerical scale.
//' @param design Finite full-column-rank observation-by-regressor matrix.
//' @return Restricted residuals, residualized intercept, coefficients, and
//'   exact numerical diagnostics.
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List cpp_ch4_afc_project(const arma::mat& y,
                               const arma::mat& design) {
  if (!y.is_finite() || !design.is_finite()) {
    Rcpp::stop("`y` and `design` must contain only finite values.");
  }
  if (y.n_rows < 2U || y.n_cols < 1U || design.n_rows != y.n_rows) {
    Rcpp::stop("The projection inputs have incompatible dimensions.");
  }
  if (design.n_cols >= design.n_rows) {
    Rcpp::stop("The design must have fewer columns than observations.");
  }

  const arma::uword n = y.n_rows;
  arma::mat coefficients(design.n_cols, y.n_cols, arma::fill::zeros);
  arma::mat residuals = y;
  arma::vec h(n, arma::fill::ones);
  double reciprocal_condition = 1.0;
  double normal_equation_residual = 0.0;

  if (design.n_cols > 0U) {
    const arma::mat cross = arma::symmatu(design.t() * design);
    arma::mat root;
    if (!arma::chol(root, cross)) {
      Rcpp::stop(
        "The conditional-alpha design cross-product is not positive "
        "definite; no generalized inverse or ridge is used."
      );
    }
    reciprocal_condition = arma::rcond(cross);
    if (!std::isfinite(reciprocal_condition) ||
        !(reciprocal_condition > 0.0)) {
      Rcpp::stop("The conditional-alpha design has zero numerical reciprocal condition.");
    }
    const arma::mat rhs = design.t() * y;
    coefficients = arma::solve(
      arma::trimatu(root),
      arma::solve(arma::trimatl(root.t()), rhs)
    );
    const arma::vec one_rhs = design.t() * arma::ones<arma::vec>(n);
    const arma::vec one_coefficient = arma::solve(
      arma::trimatu(root),
      arma::solve(arma::trimatl(root.t()), one_rhs)
    );
    residuals -= design * coefficients;
    h -= design * one_coefficient;
    if (!coefficients.is_finite() || !residuals.is_finite() ||
        !h.is_finite()) {
      Rcpp::stop("The strict conditional-alpha projection produced a non-finite result.");
    }
    const double rhs_norm = std::max(1.0, arma::norm(rhs, "fro"));
    normal_equation_residual = arma::norm(
      design.t() * residuals, "fro"
    ) / rhs_norm;
  }

  const double h2 = arma::dot(h, h);
  if (!std::isfinite(h2) || !(h2 > 0.0)) {
    Rcpp::stop("The residualized intercept has no positive finite energy.");
  }
  return Rcpp::List::create(
    Rcpp::Named("residuals") = residuals,
    Rcpp::Named("h") = h,
    Rcpp::Named("h2") = h2,
    Rcpp::Named("coefficients") = coefficients,
    Rcpp::Named("reciprocal_condition") = reciprocal_condition,
    Rcpp::Named("normal_equation_residual") =
      normal_equation_residual
  );
}


//' Spatial Kendall matrix for robust latent-factor extraction
//'
//' @param x Finite observation-by-variable matrix.
//' @return Spatial Kendall matrix and the number of zero pair differences.
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List cpp_ch4_afc_spatial_kendall(const arma::mat& x) {
  if (!x.is_finite() || x.n_rows < 2U || x.n_cols < 1U) {
    Rcpp::stop("`x` must be a finite matrix with at least two rows and one column.");
  }
  arma::mat kendall(x.n_cols, x.n_cols, arma::fill::zeros);
  arma::uword zero_pairs = 0U;
  long double pair_count = 0.0L;
  for (arma::uword i = 0; i + 1U < x.n_rows; ++i) {
    for (arma::uword j = i + 1U; j < x.n_rows; ++j) {
      const ch4afc::SignRow direction =
        ch4afc::stable_sign(x.row(i) - x.row(j));
      zero_pairs += static_cast<arma::uword>(direction.zero);
      if (!direction.zero) {
        kendall += direction.value.t() * direction.value;
      }
      pair_count += 1.0L;
    }
  }
  kendall /= static_cast<double>(pair_count);
  kendall = arma::symmatu(kendall);
  return Rcpp::List::create(
    Rcpp::Named("matrix") = kendall,
    Rcpp::Named("pairs") = static_cast<double>(pair_count),
    Rcpp::Named("zero_pairs") = static_cast<double>(zero_pairs),
    Rcpp::Named("trace") = arma::trace(kendall)
  );
}


//' Light-tail conditional-alpha sum and maximum ingredients
//'
//' @param residuals Restricted residual matrix, observations by assets.
//' @param h Residualized intercept.
//' @param factor_count Number of observed factors in the primary marginal
//'   variance divisor.
//' @param design_columns Number of columns in the null sieve design.
//' @param compute_max Whether to require and compute maximum-test marginal
//'   variances.
//' @return Exact feasible sum, trace, variance, and maximum ingredients.
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List cpp_ch4_afc_light_components(
    const arma::mat& residuals,
    const arma::vec& h,
    const int factor_count,
    const int design_columns,
    const bool compute_max) {
  if (!residuals.is_finite() || !h.is_finite() ||
      residuals.n_rows < 2U || residuals.n_cols < 2U ||
      h.n_elem != residuals.n_rows ||
      (compute_max && factor_count < 0) ||
      design_columns < 0 ||
      static_cast<arma::uword>(design_columns) >= residuals.n_rows) {
    Rcpp::stop("Invalid light-tail conditional-alpha component inputs.");
  }
  const arma::uword T = residuals.n_rows;
  const arma::uword N = residuals.n_cols;
  const long double Tld = static_cast<long double>(T);
  const long double Nld = static_cast<long double>(N);
  const long double marginal_df = Tld -
    static_cast<long double>(factor_count) - 1.0L;
  if (compute_max && !(marginal_df > 0.0L)) {
    Rcpp::stop("The primary marginal variance divisor T-d-1 is not positive.");
  }

  long double sum_stat_numerator = 0.0L;
  long double mean_numerator = 0.0L;
  arma::vec t_squared(N, arma::fill::zeros);
  arma::vec marginal_variance(N, arma::fill::zeros);
  arma::rowvec column_means(N, arma::fill::zeros);
  for (arma::uword j = 0; j < N; ++j) {
    long double column_sum = 0.0L;
    long double square_sum = 0.0L;
    for (arma::uword t = 0; t < T; ++t) {
      const long double value = static_cast<long double>(residuals(t, j));
      column_sum += value;
      square_sum += value * value;
      const long double ht = static_cast<long double>(h(t));
      mean_numerator += value * value * ht * ht;
    }
    sum_stat_numerator += column_sum * column_sum;
    column_means(j) = static_cast<double>(column_sum / Tld);
    if (compute_max) {
      const long double variance = square_sum / marginal_df;
      if (!std::isfinite(variance) || !(variance > 0.0L)) {
        Rcpp::stop(
          "Every primary marginal residual variance must be strictly positive "
          "and finite; no floor is used."
        );
      }
      marginal_variance(j) = ch4afc::checked_double(
        variance, "marginal residual variance"
      );
      t_squared(j) = ch4afc::checked_double(
        column_sum * column_sum / (Tld * variance),
        "coordinatewise maximum statistic"
      );
    }
  }

  const arma::mat centered = residuals.each_row() - column_means;
  long double trace_sigma = 0.0L;
  for (arma::uword t = 0; t < T; ++t) {
    trace_sigma += ch4afc::stable_dot(centered.row(t), centered.row(t));
  }
  trace_sigma /= Tld;

  long double trace_sigma_squared = 0.0L;
  for (arma::uword t = 0; t < T; ++t) {
    for (arma::uword s = 0; s < T; ++s) {
      const long double dot = ch4afc::stable_dot(
        centered.row(t), centered.row(s)
      );
      trace_sigma_squared += dot * dot;
    }
  }
  trace_sigma_squared /= Tld * Tld;

  const long double q = static_cast<long double>(design_columns);
  const long double trace_correction_denominator =
    (Tld + q - 1.0L) * (Tld - q);
  const long double inner = trace_sigma_squared -
    trace_sigma * trace_sigma / (Tld - q);
  const long double trace_estimate = Tld * Tld * inner /
    trace_correction_denominator;

  long double h2 = 0.0L;
  long double h4 = 0.0L;
  for (arma::uword t = 0; t < T; ++t) {
    const long double square = static_cast<long double>(h(t)) *
      static_cast<long double>(h(t));
    h2 += square;
    h4 += square * square;
  }
  const long double off_diagonal_h4 = h2 * h2 - h4;
  const long double variance_estimate = 2.0L * trace_estimate *
    off_diagonal_h4 / (Nld * Nld * Tld * Tld);

  return Rcpp::List::create(
    Rcpp::Named("sum_statistic") = ch4afc::checked_double(
      sum_stat_numerator / (Nld * Tld), "conditional sum statistic"
    ),
    Rcpp::Named("mean_estimate") = ch4afc::checked_double(
      mean_numerator / (Nld * Tld), "conditional null mean estimate"
    ),
    Rcpp::Named("variance_estimate") = ch4afc::checked_double(
      variance_estimate, "conditional null variance estimate"
    ),
    Rcpp::Named("trace_sigma") = ch4afc::checked_double(
      trace_sigma, "residual covariance trace"
    ),
    Rcpp::Named("trace_sigma_squared_raw") = ch4afc::checked_double(
      trace_sigma_squared, "squared residual covariance trace"
    ),
    Rcpp::Named("trace_sigma_squared_estimate") = ch4afc::checked_double(
      trace_estimate, "bias-corrected squared covariance trace"
    ),
    Rcpp::Named("trace_inner") = ch4afc::checked_double(
      inner, "bias-corrected squared covariance trace inner term"
    ),
    Rcpp::Named("trace_correction_denominator") =
      ch4afc::checked_double(
        trace_correction_denominator, "trace correction denominator"
      ),
    Rcpp::Named("off_diagonal_h4") = ch4afc::checked_double(
      off_diagonal_h4, "off-diagonal h fourth-moment sum"
    ),
    Rcpp::Named("marginal_variance") = marginal_variance,
    Rcpp::Named("coordinate_t_squared") = t_squared,
    Rcpp::Named("maximum_statistic") = compute_max ?
      Rcpp::wrap(t_squared.max()) : R_NilValue,
    Rcpp::Named("marginal_df") = static_cast<double>(marginal_df)
  );
}


//' Robust conditional spatial-sign sum ingredients
//'
//' @param residuals Restricted residuals used in the CSS numerator.
//' @param trace_residuals Residuals from the primary uncentered-basis fit used
//'   in the off-diagonal trace estimator.
//' @param h Residualized intercept from the restricted null sieve design.
//' @return CSS numerator, trace estimator, signs, and diagnostics.
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List cpp_ch4_afc_css_components(
    const arma::mat& residuals,
    const arma::mat& trace_residuals,
    const arma::vec& h) {
  if (!residuals.is_finite() || !trace_residuals.is_finite() ||
      !h.is_finite() || residuals.n_rows < 2U ||
      residuals.n_cols < 2U ||
      trace_residuals.n_rows != residuals.n_rows ||
      trace_residuals.n_cols != residuals.n_cols ||
      h.n_elem != residuals.n_rows) {
    Rcpp::stop("Invalid conditional spatial-sign sum inputs.");
  }
  arma::uword numerator_zero_rows = 0U;
  arma::uword trace_zero_rows = 0U;
  const arma::mat numerator_signs = ch4afc::signs(
    residuals, numerator_zero_rows
  );
  const arma::mat trace_signs = ch4afc::signs(
    trace_residuals, trace_zero_rows
  );
  const long double h2 = static_cast<long double>(arma::dot(h, h));
  if (!(h2 > 1.0L) || !std::isfinite(h2)) {
    Rcpp::stop(
      "The primary CSS trace denominator requires h'h > 1."
    );
  }

  arma::rowvec weighted_sum(residuals.n_cols, arma::fill::zeros);
  for (arma::uword t = 0; t < residuals.n_rows; ++t) {
    weighted_sum += h(t) * numerator_signs.row(t);
  }
  const long double quadratic = ch4afc::stable_dot(
    weighted_sum, weighted_sum
  ) / h2;
  const long double css_numerator = quadratic - 1.0L;

  long double trace_numerator = 0.0L;
  for (arma::uword t = 0; t + 1U < residuals.n_rows; ++t) {
    const long double ht2 = static_cast<long double>(h(t)) *
      static_cast<long double>(h(t));
    for (arma::uword s = t + 1U; s < residuals.n_rows; ++s) {
      const long double hs2 = static_cast<long double>(h(s)) *
        static_cast<long double>(h(s));
      const long double dot = ch4afc::stable_dot(
        trace_signs.row(t), trace_signs.row(s)
      );
      trace_numerator += 2.0L * ht2 * hs2 * dot * dot;
    }
  }
  const long double trace_denominator = h2 * (h2 - 1.0L);
  const long double trace_estimate = trace_numerator /
    trace_denominator;

  return Rcpp::List::create(
    Rcpp::Named("quadratic") = ch4afc::checked_double(
      quadratic, "CSS sign quadratic"
    ),
    Rcpp::Named("numerator") = ch4afc::checked_double(
      css_numerator, "CSS centered numerator"
    ),
    Rcpp::Named("trace_estimate") = ch4afc::checked_double(
      trace_estimate, "CSS squared sign-scatter trace estimate"
    ),
    Rcpp::Named("trace_numerator") = ch4afc::checked_double(
      trace_numerator, "CSS trace numerator"
    ),
    Rcpp::Named("trace_denominator") = ch4afc::checked_double(
      trace_denominator, "CSS trace denominator"
    ),
    Rcpp::Named("numerator_signs") = numerator_signs,
    Rcpp::Named("trace_signs") = trace_signs,
    Rcpp::Named("numerator_zero_rows") =
      static_cast<double>(numerator_zero_rows),
    Rcpp::Named("trace_zero_rows") =
      static_cast<double>(trace_zero_rows)
  );
}
