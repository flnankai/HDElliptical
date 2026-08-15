// Wang--Peng--Li one-sample high-dimensional spatial-sign test.
//
// The feasible variance is the literal leave-two-out/cross-validation
// estimator in equation (7) of Wang, Peng, and Li (2015).  We intentionally
// do not use their faster equation (8), because that simplification uses
// ||Z_i||^2 = 1 whereas the paper defines Z_i = U(0) = 0 at a null residual.

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp17)]]

#include <RcppArmadillo.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

namespace {

void wpl_require_finite_matrix(const arma::mat& x, const char* name) {
  if (!x.is_finite()) {
    Rcpp::stop("`%s` must contain only finite values.", name);
  }
}

void wpl_neumaier_add(long double value,
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

long double wpl_compensated_dot(const arma::mat& z,
                                arma::uword i,
                                arma::uword j) {
  long double total = 0.0L;
  long double correction = 0.0L;
  for (arma::uword ell = 0; ell < z.n_cols; ++ell) {
    const long double product =
      static_cast<long double>(z(i, ell)) *
      static_cast<long double>(z(j, ell));
    wpl_neumaier_add(product, total, correction);
  }
  return total + correction;
}

struct WplDirections {
  arma::mat z;
  arma::vec squared_norm;
  arma::vec sum;
  arma::vec mean;
  arma::vec sample_mean;
  arma::uword zero_count;
  arma::uword overflow_fallback_count;
};

arma::vec wpl_stable_column_mean(const arma::mat& x) {
  arma::vec output(x.n_cols, arma::fill::zeros);
  for (arma::uword j = 0; j < x.n_cols; ++j) {
    double scale = 0.0;
    for (arma::uword i = 0; i < x.n_rows; ++i) {
      scale = std::max(scale, std::abs(x(i, j)));
    }
    if (scale == 0.0) {
      continue;
    }
    long double total = 0.0L;
    long double correction = 0.0L;
    for (arma::uword i = 0; i < x.n_rows; ++i) {
      wpl_neumaier_add(
        static_cast<long double>(x(i, j) / scale), total, correction
      );
    }
    const long double scaled_mean =
      (total + correction) / static_cast<long double>(x.n_rows);
    output(j) = static_cast<double>(scaled_mean) * scale;
  }
  return output;
}

WplDirections wpl_directions(const arma::mat& x, const arma::vec& mu) {
  const arma::uword n = x.n_rows;
  const arma::uword p = x.n_cols;
  arma::mat z(n, p, arma::fill::zeros);
  arma::vec squared_norm(n, arma::fill::zeros);
  arma::vec direction_sum(p, arma::fill::zeros);
  std::vector<long double> column_total(p, 0.0L);
  std::vector<long double> column_correction(p, 0.0L);
  std::vector<double> residual(p, 0.0);
  arma::uword zero_count = 0;
  arma::uword overflow_fallback_count = 0;

  for (arma::uword i = 0; i < n; ++i) {
    bool subtraction_overflow = false;
    double max_abs = 0.0;
    for (arma::uword j = 0; j < p; ++j) {
      residual[j] = x(i, j) - mu(j);
      if (!std::isfinite(residual[j])) {
        subtraction_overflow = true;
      } else {
        max_abs = std::max(max_abs, std::abs(residual[j]));
      }
    }

    // Use direct subtraction whenever it is representable, preserving nearby
    // large floating-point differences.  Only an actual overflow triggers
    // operand prescaling.
    if (subtraction_overflow) {
      ++overflow_fallback_count;
      double operand_scale = 0.0;
      for (arma::uword j = 0; j < p; ++j) {
        operand_scale = std::max(operand_scale, std::abs(x(i, j)));
        operand_scale = std::max(operand_scale, std::abs(mu(j)));
      }
      if (!(operand_scale > 0.0) || !std::isfinite(operand_scale)) {
        Rcpp::stop(
          "Wang-Peng-Li could not form a finite residual direction."
        );
      }
      max_abs = 0.0;
      for (arma::uword j = 0; j < p; ++j) {
        residual[j] = x(i, j) / operand_scale - mu(j) / operand_scale;
        max_abs = std::max(max_abs, std::abs(residual[j]));
      }
    }

    if (max_abs == 0.0) {
      ++zero_count;
      continue;
    }

    long double sum_squares = 0.0L;
    long double square_correction = 0.0L;
    for (arma::uword j = 0; j < p; ++j) {
      residual[j] /= max_abs;
      const long double value = static_cast<long double>(residual[j]);
      wpl_neumaier_add(
        value * value, sum_squares, square_correction
      );
    }
    const long double norm = std::sqrt(sum_squares + square_correction);
    if (!(norm > 0.0L) || !std::isfinite(norm)) {
      Rcpp::stop(
        "Wang-Peng-Li could not normalise a finite nonzero residual."
      );
    }

    long double row_norm_squared = 0.0L;
    long double row_correction = 0.0L;
    for (arma::uword j = 0; j < p; ++j) {
      const double direction = static_cast<double>(
        static_cast<long double>(residual[j]) / norm
      );
      z(i, j) = direction;
      const long double direction_ld = static_cast<long double>(direction);
      wpl_neumaier_add(
        direction_ld * direction_ld,
        row_norm_squared,
        row_correction
      );
      wpl_neumaier_add(
        direction_ld, column_total[j], column_correction[j]
      );
    }
    squared_norm(i) = static_cast<double>(
      row_norm_squared + row_correction
    );
  }

  for (arma::uword j = 0; j < p; ++j) {
    direction_sum(j) = static_cast<double>(
      column_total[j] + column_correction[j]
    );
  }

  return WplDirections{
    z,
    squared_norm,
    direction_sum,
    direction_sum / static_cast<double>(n),
    wpl_stable_column_mean(x),
    zero_count,
    overflow_fallback_count
  };
}

}  // namespace


//' Wang--Peng--Li one-sample spatial-sign statistic kernel
//'
//' @param x Numeric observation-by-variable matrix.
//' @param mu Numeric null-location vector.
//' @return Internal list of statistic and cross-validation components.
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List cpp_wang_peng_li_one_sample(const arma::mat& x,
                                       const arma::vec& mu) {
  wpl_require_finite_matrix(x, "x");
  if (!mu.is_finite()) {
    Rcpp::stop("`mu` must contain only finite values.");
  }

  const arma::uword n = x.n_rows;
  const arma::uword p = x.n_cols;
  if (n < 3) {
    Rcpp::stop(
      "Wang-Peng-Li requires at least three observations for its "
      "leave-two-out variance estimator."
    );
  }
  if (p < 1 || mu.n_elem != p) {
    Rcpp::stop("`mu` must have one entry for every column of `x`.");
  }

  const WplDirections directions = wpl_directions(x, mu);
  const arma::mat& z = directions.z;

  long double t_total = 0.0L;
  long double t_correction = 0.0L;
  for (arma::uword j = 0; j + 1 < n; ++j) {
    for (arma::uword k = j + 1; k < n; ++k) {
      wpl_neumaier_add(
        wpl_compensated_dot(z, j, k), t_total, t_correction
      );
    }
  }
  const long double t_raw_ld = t_total + t_correction;

  std::vector<long double> sign_sum(p, 0.0L);
  std::vector<long double> sign_sum_correction(p, 0.0L);
  for (arma::uword i = 0; i < n; ++i) {
    for (arma::uword ell = 0; ell < p; ++ell) {
      wpl_neumaier_add(
        static_cast<long double>(z(i, ell)),
        sign_sum[ell],
        sign_sum_correction[ell]
      );
    }
  }
  for (arma::uword ell = 0; ell < p; ++ell) {
    sign_sum[ell] += sign_sum_correction[ell];
  }

  // Equation (7) of Wang, Peng, and Li (2015).  For a rank-one product,
  // tr{(Zj-m) Zj' (Zk-m) Zk'} equals
  // [Zj'(Zk-m)] [Zk'(Zj-m)].  Terms for (j,k) and (k,j) coincide, so we
  // evaluate each unordered pair once and multiply its contribution by two.
  const long double leaveout_divisor = static_cast<long double>(n - 2);
  long double cv_total = 0.0L;
  long double cv_correction = 0.0L;
  for (arma::uword j = 0; j + 1 < n; ++j) {
    for (arma::uword k = j + 1; k < n; ++k) {
      long double first = 0.0L;
      long double first_correction = 0.0L;
      long double second = 0.0L;
      long double second_correction = 0.0L;
      for (arma::uword ell = 0; ell < p; ++ell) {
        const long double zj = static_cast<long double>(z(j, ell));
        const long double zk = static_cast<long double>(z(k, ell));
        const long double leaveout_mean =
          (sign_sum[ell] - zj - zk) / leaveout_divisor;
        wpl_neumaier_add(
          zj * (zk - leaveout_mean), first, first_correction
        );
        wpl_neumaier_add(
          zk * (zj - leaveout_mean), second, second_correction
        );
      }
      const long double contribution =
        (first + first_correction) * (second + second_correction);
      wpl_neumaier_add(
        2.0L * contribution, cv_total, cv_correction
      );
    }
  }

  const long double ordered_pair_count =
    static_cast<long double>(n) * static_cast<long double>(n - 1);
  const long double unordered_pair_count = ordered_pair_count / 2.0L;
  const long double cv_numerator_ld = cv_total + cv_correction;
  const long double trace_b2_ld = cv_numerator_ld / ordered_pair_count;
  const long double variance_ld = unordered_pair_count * trace_b2_ld;

  const double t_raw = static_cast<double>(t_raw_ld);
  const double cv_numerator = static_cast<double>(cv_numerator_ld);
  const double trace_b2_hat = static_cast<double>(trace_b2_ld);
  const double variance_hat = static_cast<double>(variance_ld);
  if (!std::isfinite(trace_b2_hat) || trace_b2_hat <= 0.0 ||
      !std::isfinite(variance_hat) || variance_hat <= 0.0) {
    Rcpp::stop(
      "Wang-Peng-Li requires its literal leave-two-out estimate of "
      "tr(B^2) to be finite and strictly positive. It was %.17g with %llu "
      "zero spatial sign(s); no absolute-value or variance-floor repair is "
      "applied.",
      trace_b2_hat,
      static_cast<unsigned long long>(directions.zero_count)
    );
  }
  if (!std::isfinite(t_raw) || !std::isfinite(cv_numerator)) {
    Rcpp::stop("Wang-Peng-Li produced a non-finite raw component.");
  }

  const double standard_error = std::sqrt(variance_hat);
  const double statistic = t_raw / standard_error;
  if (!std::isfinite(statistic)) {
    Rcpp::stop("Wang-Peng-Li produced a non-finite standardised statistic.");
  }

  long double norm_total = 0.0L;
  long double norm_correction = 0.0L;
  for (arma::uword i = 0; i < n; ++i) {
    wpl_neumaier_add(
      static_cast<long double>(directions.squared_norm(i)),
      norm_total,
      norm_correction
    );
  }

  return Rcpp::List::create(
    Rcpp::Named("z") = statistic,
    Rcpp::Named("T_raw") = t_raw,
    Rcpp::Named("trace_B2_hat") = trace_b2_hat,
    Rcpp::Named("variance_hat") = variance_hat,
    Rcpp::Named("standard_error") = standard_error,
    Rcpp::Named("cross_validation_numerator") = cv_numerator,
    Rcpp::Named("pair_count") = static_cast<double>(unordered_pair_count),
    Rcpp::Named("ordered_pair_count") =
      static_cast<double>(ordered_pair_count),
    Rcpp::Named("direction_sum") = directions.sum,
    Rcpp::Named("direction_mean") = directions.mean,
    Rcpp::Named("sum_direction_norm_squared") =
      static_cast<double>(norm_total + norm_correction),
    Rcpp::Named("sample_mean") = directions.sample_mean,
    Rcpp::Named("zero_signs") =
      static_cast<double>(directions.zero_count),
    Rcpp::Named("nonzero_signs") =
      static_cast<double>(n - directions.zero_count),
    Rcpp::Named("overflow_fallback_rows") =
      static_cast<double>(directions.overflow_fallback_count),
    Rcpp::Named("n") = static_cast<double>(n),
    Rcpp::Named("p") = static_cast<double>(p)
  );
}
