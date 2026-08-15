// Copyright (C) 2026
// SPDX-License-Identifier: GPL-3.0-or-later

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp17)]]

#include <RcppArmadillo.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <string>
#include <vector>

namespace ch4alpha {

struct DirectionResult {
  arma::rowvec direction;
  double radius;
  bool zero;
};

inline DirectionResult direction_scaled(
    const arma::rowvec& residual,
    const arma::vec& diagonal,
    const double zero_tol) {
  const arma::uword p = residual.n_elem;
  arma::rowvec standardized(p, arma::fill::zeros);
  double maximum = 0.0;

  for (arma::uword j = 0; j < p; ++j) {
    const double root = std::sqrt(diagonal[j]);
    if (!std::isfinite(root) || !(root > 0.0)) {
      Rcpp::stop("The diagonal sign scale is not strictly positive and finite.");
    }
    const double value = residual[j] / root;
    if (!std::isfinite(value)) {
      Rcpp::stop("A standardized factor residual is not finite.");
    }
    standardized[j] = value;
    maximum = std::max(maximum, std::abs(value));
  }

  DirectionResult result{arma::rowvec(p, arma::fill::zeros), 0.0, true};
  if (maximum == 0.0) {
    return result;
  }

  const arma::rowvec scaled = standardized / maximum;
  const double scaled_norm = arma::norm(scaled, 2);
  if (!std::isfinite(scaled_norm) || !(scaled_norm > 0.0)) {
    Rcpp::stop("A standardized factor-residual norm is invalid.");
  }
  if (zero_tol > 0.0 && maximum <= zero_tol / scaled_norm) {
    return result;
  }

  const double radius = maximum * scaled_norm;
  if (!std::isfinite(radius) || !(radius > 0.0)) {
    Rcpp::stop("A standardized factor-residual radius is not representable.");
  }
  result.direction = scaled / scaled_norm;
  result.radius = radius;
  result.zero = false;
  return result;
}

inline arma::mat restricted_slopes(const arma::mat& y,
                                   const arma::mat& factors,
                                   const arma::uvec* rows,
                                   const std::string& label) {
  const arma::uword k = factors.n_cols;
  if (k == 0U) {
    return arma::mat(0U, y.n_cols);
  }

  arma::mat f_use;
  arma::mat y_use;
  if (rows == nullptr) {
    f_use = factors;
    y_use = y;
  } else {
    f_use = factors.rows(*rows);
    y_use = y.rows(*rows);
  }
  if (f_use.n_rows < k) {
    Rcpp::stop(label + " has fewer observations than factor columns.");
  }

  arma::mat answer;
  const arma::mat cross = f_use.t() * f_use;
  const arma::mat rhs = f_use.t() * y_use;
  const bool solved = arma::solve(
    answer, cross, rhs,
    arma::solve_opts::likely_sympd + arma::solve_opts::no_approx
  );
  if (!solved || !answer.is_finite()) {
    Rcpp::stop(label + " has a singular factor cross-product; no generalized inverse is used.");
  }
  return answer;
}

struct DiagonalFit {
  arma::vec diagonal;
  arma::mat directions;
  arma::vec radii;
  int iterations;
  bool converged;
  double equation_residual;
  double relative_update;
  int zero_residuals;
};

inline DiagonalFit fit_diagonal_sign_scale(
    const arma::mat& residual,
    const double tol,
    const int max_iter,
    const double zero_tol) {
  const arma::uword n = residual.n_rows;
  const arma::uword p = residual.n_cols;
  arma::vec diagonal(p, arma::fill::zeros);

  for (arma::uword j = 0; j < p; ++j) {
    long double total = 0.0L;
    for (arma::uword i = 0; i < n; ++i) {
      const long double value = static_cast<long double>(residual(i, j));
      total += value * value;
    }
    diagonal[j] = static_cast<double>(total / static_cast<long double>(n));
    if (!std::isfinite(diagonal[j]) || !(diagonal[j] > 0.0)) {
      Rcpp::stop("Every asset must have positive restricted-residual scale.");
    }
  }
  diagonal /= arma::mean(diagonal);

  arma::mat directions(n, p, arma::fill::zeros);
  arma::vec radii(n, arma::fill::zeros);
  double equation_residual = std::numeric_limits<double>::infinity();
  double relative_update = std::numeric_limits<double>::infinity();
  int zero_residuals = 0;
  bool converged = false;
  int iteration = 0;

  for (iteration = 1; iteration <= max_iter; ++iteration) {
    arma::vec moments(p, arma::fill::zeros);
    zero_residuals = 0;
    for (arma::uword i = 0; i < n; ++i) {
      const DirectionResult current = direction_scaled(
        residual.row(i), diagonal, zero_tol
      );
      directions.row(i) = current.direction;
      radii[i] = current.radius;
      zero_residuals += current.zero ? 1 : 0;
      moments += arma::square(current.direction).t();
    }
    moments *= static_cast<double>(p) / static_cast<double>(n);
    equation_residual = arma::max(arma::abs(moments - 1.0));

    arma::vec updated = diagonal % moments;
    if (!updated.is_finite() || arma::any(updated <= 0.0)) {
      Rcpp::stop(
        "The diagonal spatial-sign equation became singular; no floor or ridge is applied."
      );
    }
    const double normalizer = arma::mean(updated);
    if (!std::isfinite(normalizer) || !(normalizer > 0.0)) {
      Rcpp::stop("The diagonal spatial-sign normalization is invalid.");
    }
    updated /= normalizer;
    relative_update = arma::max(arma::abs(arma::log(updated / diagonal)));
    diagonal = updated;

    if (equation_residual <= tol && relative_update <= tol) {
      converged = true;
      break;
    }
  }

  // Diagnostics and statistics must correspond to the returned diagonal.
  arma::vec final_moments(p, arma::fill::zeros);
  zero_residuals = 0;
  for (arma::uword i = 0; i < n; ++i) {
    const DirectionResult current = direction_scaled(
      residual.row(i), diagonal, zero_tol
    );
    directions.row(i) = current.direction;
    radii[i] = current.radius;
    zero_residuals += current.zero ? 1 : 0;
    final_moments += arma::square(current.direction).t();
  }
  final_moments *= static_cast<double>(p) / static_cast<double>(n);
  equation_residual = arma::max(arma::abs(final_moments - 1.0));
  converged = converged && equation_residual <= tol;

  return DiagonalFit{
    diagonal, directions, radii, std::min(iteration, max_iter), converged,
    equation_residual, relative_update, zero_residuals
  };
}

inline arma::vec factor_intercept_residual(const arma::mat& factors) {
  const arma::uword n = factors.n_rows;
  arma::vec h(n, arma::fill::ones);
  if (factors.n_cols == 0U) {
    return h;
  }
  arma::vec coefficient;
  const bool solved = arma::solve(
    coefficient, factors.t() * factors, factors.t() * h,
    arma::solve_opts::likely_sympd + arma::solve_opts::no_approx
  );
  if (!solved || !coefficient.is_finite()) {
    Rcpp::stop("The factor cross-product is singular; no generalized inverse is used.");
  }
  h -= factors * coefficient;
  return h;
}

inline double q_statistic(const arma::mat& directions,
                          const arma::vec& h) {
  const double h2 = arma::dot(h, h);
  const arma::rowvec weighted_sum = h.t() * directions;
  long double off_diagonal = static_cast<long double>(
    arma::dot(weighted_sum, weighted_sum)
  );
  for (arma::uword i = 0; i < directions.n_rows; ++i) {
    off_diagonal -= static_cast<long double>(h[i]) *
      static_cast<long double>(h[i]) *
      static_cast<long double>(arma::dot(directions.row(i), directions.row(i)));
  }
  return static_cast<double>(
    static_cast<long double>(directions.n_cols) * off_diagonal /
      static_cast<long double>(h2)
  );
}

inline double split_trace_estimator(const arma::mat& y,
                                    const arma::mat& factors,
                                    const arma::vec& h,
                                    const arma::vec& diagonal,
                                    const double zero_tol) {
  const arma::uword n = y.n_rows;
  const arma::uword p = y.n_cols;
  const arma::uword k = factors.n_cols;
  const double h2 = arma::dot(h, h);
  if (!std::isfinite(h2) || !(h2 > 1.0)) {
    Rcpp::stop("The primary trace estimator requires h'h > 1.");
  }
  if (n < 2U * k + 2U) {
    Rcpp::stop(
      "The split leave-two-out trace estimator requires T >= 2*K + 2."
    );
  }

  long double total = 0.0L;
  std::vector<arma::uword> remaining;
  remaining.reserve(n - 2U);

  for (arma::uword first = 0; first < n; ++first) {
    for (arma::uword second = first + 1U; second < n; ++second) {
      remaining.clear();
      for (arma::uword row = 0; row < n; ++row) {
        if (row != first && row != second) {
          remaining.push_back(row);
        }
      }
      const arma::uword split = remaining.size() / 2U;
      arma::uvec rows1(split);
      arma::uvec rows2(remaining.size() - split);
      for (arma::uword index = 0; index < split; ++index) {
        rows1[index] = remaining[index];
      }
      for (arma::uword index = split; index < remaining.size(); ++index) {
        rows2[index - split] = remaining[index];
      }

      const std::string pair = "A split leave-two-out factor fit for observations " +
        std::to_string(first + 1U) + " and " + std::to_string(second + 1U);
      const arma::mat slope1 = restricted_slopes(y, factors, &rows1, pair + " (first half)");
      const arma::mat slope2 = restricted_slopes(y, factors, &rows2, pair + " (second half)");
      arma::rowvec residual_first_forward = y.row(first);
      arma::rowvec residual_second_forward = y.row(second);
      arma::rowvec residual_second_reverse = y.row(second);
      arma::rowvec residual_first_reverse = y.row(first);
      if (k > 0U) {
        residual_first_forward -= factors.row(first) * slope1;
        residual_second_forward -= factors.row(second) * slope2;
        // The primary estimator is an ordered-pair sum. Reversing the pair
        // swaps which deleted observation uses the first- and second-half
        // slope fit; it is not generally twice the forward contribution.
        residual_second_reverse -= factors.row(second) * slope1;
        residual_first_reverse -= factors.row(first) * slope2;
      }
      const DirectionResult direction_first_forward = direction_scaled(
        residual_first_forward, diagonal, zero_tol
      );
      const DirectionResult direction_second_forward = direction_scaled(
        residual_second_forward, diagonal, zero_tol
      );
      const DirectionResult direction_second_reverse = direction_scaled(
        residual_second_reverse, diagonal, zero_tol
      );
      const DirectionResult direction_first_reverse = direction_scaled(
        residual_first_reverse, diagonal, zero_tol
      );
      const long double product_forward = static_cast<long double>(
        arma::dot(
          direction_first_forward.direction,
          direction_second_forward.direction
        )
      );
      const long double product_reverse = static_cast<long double>(
        arma::dot(
          direction_second_reverse.direction,
          direction_first_reverse.direction
        )
      );
      const long double weight = static_cast<long double>(h[first]) *
        static_cast<long double>(h[first]) *
        static_cast<long double>(h[second]) *
        static_cast<long double>(h[second]);
      total += weight * (
        product_forward * product_forward +
        product_reverse * product_reverse
      );
    }
  }

  const long double multiplier = static_cast<long double>(p) *
    static_cast<long double>(p) /
    (static_cast<long double>(h2) * static_cast<long double>(h2 - 1.0));
  const double answer = static_cast<double>(multiplier * total);
  if (!std::isfinite(answer) || !(answer > 0.0)) {
    Rcpp::stop(
      "The split leave-two-out trace estimate is not strictly positive and finite; no floor is applied."
    );
  }
  return answer;
}

}  // namespace ch4alpha


//' Core OLS quantities for unconditional factor-pricing alpha tests
//'
//' @param y Observation-by-asset numeric matrix, already column scaled.
//' @param factors Observation-by-factor numeric matrix, already column scaled.
//' @return A list of OLS quantities in the scaled coordinates.
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List cpp_ch4_alpha_ols(const arma::mat& y,
                             const arma::mat& factors) {
  if (y.n_rows != factors.n_rows) {
    Rcpp::stop("`y` and `factors` must have the same number of rows.");
  }
  if (y.n_rows < 2U || y.n_cols < 1U || !y.is_finite() || !factors.is_finite()) {
    Rcpp::stop("The scaled alpha-test inputs are invalid.");
  }

  const arma::uword n = y.n_rows;
  const arma::uword assets = y.n_cols;
  const arma::uword k = factors.n_cols;
  if (n <= k + 1U) {
    Rcpp::stop("The joint intercept-factor regression has no positive residual degrees of freedom.");
  }

  arma::mat design(n, k + 1U, arma::fill::ones);
  if (k > 0U) {
    design.cols(1U, k) = factors;
  }
  arma::mat coefficients;
  const bool solved = arma::solve(
    coefficients, design.t() * design, design.t() * y,
    arma::solve_opts::likely_sympd + arma::solve_opts::no_approx
  );
  if (!solved || !coefficients.is_finite()) {
    Rcpp::stop("The joint intercept-factor design is rank deficient; no generalized inverse is used.");
  }
  const arma::mat residuals = y - design * coefficients;
  const arma::vec h = ch4alpha::factor_intercept_residual(factors);
  const double h2 = arma::dot(h, h);
  if (!std::isfinite(h2) || !(h2 > 0.0)) {
    Rcpp::stop("The intercept is in the factor span, so alpha is not identifiable.");
  }

  arma::vec sse = arma::sum(arma::square(residuals), 0).t();
  if (!sse.is_finite() || arma::any(sse <= 0.0)) {
    Rcpp::stop("Every asset must have a strictly positive OLS residual sum of squares.");
  }
  const double v = static_cast<double>(n - k - 1U);
  const arma::vec alpha = coefficients.row(0).t();
  const arma::vec t_squared = arma::square(alpha) * h2 / (sse / v);

  const arma::mat covariance_t = residuals.t() * residuals;
  arma::mat correlation(assets, assets, arma::fill::eye);
  for (arma::uword i = 0; i < assets; ++i) {
    for (arma::uword j = i + 1U; j < assets; ++j) {
      const double value = covariance_t(i, j) / std::sqrt(sse[i] * sse[j]);
      if (!std::isfinite(value)) {
        Rcpp::stop("An OLS residual correlation is not finite.");
      }
      correlation(i, j) = value;
      correlation(j, i) = value;
    }
  }

  return Rcpp::List::create(
    Rcpp::Named("alpha") = alpha,
    Rcpp::Named("coefficients") = coefficients,
    Rcpp::Named("residuals") = residuals,
    Rcpp::Named("sse") = sse,
    Rcpp::Named("t_squared") = t_squared,
    Rcpp::Named("residual_correlation") = correlation,
    Rcpp::Named("h") = h,
    Rcpp::Named("h2") = h2,
    Rcpp::Named("residual_df") = v
  );
}


//' Liu--Feng--Ma spatial-sign alpha core
//'
//' @param y Observation-by-asset numeric matrix, already column scaled.
//' @param factors Observation-by-factor numeric matrix, already column scaled.
//' @param tol Positive fixed-point tolerance.
//' @param max_iter Positive iteration limit.
//' @param zero_tol Non-negative standardized zero-radius tolerance.
//' @param compute_trace Whether to compute the primary split leave-two-out trace.
//' @return Components of the primary statistic.
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List cpp_ch4_lfm_spatial_sign_core(
    const arma::mat& y,
    const arma::mat& factors,
    const double tol,
    const int max_iter,
    const double zero_tol,
    const bool compute_trace) {
  if (y.n_rows != factors.n_rows || y.n_rows < 2U || y.n_cols < 2U ||
      !y.is_finite() || !factors.is_finite()) {
    Rcpp::stop("The scaled spatial-sign alpha inputs are invalid.");
  }
  if (!std::isfinite(tol) || !(tol > 0.0) || max_iter < 1 ||
      !std::isfinite(zero_tol) || zero_tol < 0.0) {
    Rcpp::stop("The spatial-sign iteration controls are invalid.");
  }

  const arma::mat slopes = ch4alpha::restricted_slopes(
    y, factors, nullptr, "The full-sample restricted factor fit"
  );
  arma::mat residual = y;
  if (factors.n_cols > 0U) {
    residual -= factors * slopes;
  }
  const ch4alpha::DiagonalFit fit = ch4alpha::fit_diagonal_sign_scale(
    residual, tol, max_iter, zero_tol
  );
  const arma::vec h = ch4alpha::factor_intercept_residual(factors);
  const double h2 = arma::dot(h, h);
  const double q = ch4alpha::q_statistic(fit.directions, h);
  if (!std::isfinite(q)) {
    Rcpp::stop("The spatial-sign alpha quadratic form is not finite.");
  }

  double trace = NA_REAL;
  if (compute_trace) {
    trace = ch4alpha::split_trace_estimator(
      y, factors, h, fit.diagonal, zero_tol
    );
  }

  return Rcpp::List::create(
    Rcpp::Named("Q") = q,
    Rcpp::Named("trace_R2") = trace,
    Rcpp::Named("diagonal") = fit.diagonal,
    Rcpp::Named("directions") = fit.directions,
    Rcpp::Named("radii") = fit.radii,
    Rcpp::Named("restricted_slopes") = slopes,
    Rcpp::Named("restricted_residuals") = residual,
    Rcpp::Named("h") = h,
    Rcpp::Named("h2") = h2,
    Rcpp::Named("iterations") = fit.iterations,
    Rcpp::Named("converged") = fit.converged,
    Rcpp::Named("equation_residual") = fit.equation_residual,
    Rcpp::Named("relative_update") = fit.relative_update,
    Rcpp::Named("zero_residuals") = fit.zero_residuals,
    Rcpp::Named("trace_denominator") = h2 * (h2 - 1.0),
    Rcpp::Named("trace_split") = "chronological first/second halves after deleting each pair"
  );
}
