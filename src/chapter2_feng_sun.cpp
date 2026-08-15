// Feng--Sun one-sample scalar-invariant high-dimensional spatial-sign test.
//
// Each unordered observation pair gets its own diagonal HR-type location and
// scale fit from the leave-two-out sample.  The test numerator uses only the
// leave-out scale, whereas the feasible trace estimator uses both the
// leave-out location and scale, exactly as in Feng and Sun (2016).

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp17)]]

#include <RcppArmadillo.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <utility>
#include <vector>

namespace {

void fs_neumaier_add(long double value,
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

struct FsPreparedData {
  arma::mat residual;
  arma::vec anchor;
  arma::vec residual_scale_ratio;
  arma::vec log_residual_scale;
  arma::vec sample_mean;
  arma::uword subtraction_overflow_columns;
};

arma::vec fs_stable_column_mean(const arma::mat& x) {
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
      fs_neumaier_add(
        static_cast<long double>(x(i, j) / scale), total, correction
      );
    }
    const long double scaled_mean =
      (total + correction) / static_cast<long double>(x.n_rows);
    output(j) = static_cast<double>(scaled_mean) * scale;
  }
  return output;
}

FsPreparedData fs_prepare_data(const arma::mat& x, const arma::vec& mu) {
  const arma::uword n = x.n_rows;
  const arma::uword p = x.n_cols;
  arma::mat residual(n, p, arma::fill::zeros);
  arma::vec anchor(p, arma::fill::zeros);
  arma::vec ratio(p, arma::fill::zeros);
  arma::vec log_scale(p, arma::fill::zeros);
  arma::uword fallback_columns = 0;

  for (arma::uword j = 0; j < p; ++j) {
    double column_anchor = std::abs(mu(j));
    for (arma::uword i = 0; i < n; ++i) {
      column_anchor = std::max(column_anchor, std::abs(x(i, j)));
    }
    if (!(column_anchor > 0.0) || !std::isfinite(column_anchor)) {
      Rcpp::stop(
        "Feng-Sun requires every variable to have a nonzero finite scale; "
        "variable %llu is identically zero at the null location.",
        static_cast<unsigned long long>(j + 1)
      );
    }
    anchor(j) = column_anchor;

    bool subtraction_overflow = false;
    double max_abs = 0.0;
    for (arma::uword i = 0; i < n; ++i) {
      const double value = x(i, j) - mu(j);
      residual(i, j) = value;
      if (!std::isfinite(value)) {
        subtraction_overflow = true;
      } else {
        max_abs = std::max(max_abs, std::abs(value));
      }
    }

    if (subtraction_overflow) {
      ++fallback_columns;
      max_abs = 0.0;
      for (arma::uword i = 0; i < n; ++i) {
        const double value =
          x(i, j) / column_anchor - mu(j) / column_anchor;
        residual(i, j) = value;
        max_abs = std::max(max_abs, std::abs(value));
      }
      ratio(j) = max_abs;
    } else {
      ratio(j) = max_abs / column_anchor;
    }

    if (!(max_abs > 0.0) || !std::isfinite(max_abs) ||
        !(ratio(j) > 0.0) || !std::isfinite(ratio(j))) {
      Rcpp::stop(
        "Feng-Sun requires variable %llu to differ from its null location "
        "in at least one observation.",
        static_cast<unsigned long long>(j + 1)
      );
    }
    residual.col(j) /= max_abs;
    log_scale(j) = std::log(column_anchor) + std::log(ratio(j));

    double minimum = residual(0, j);
    double maximum = residual(0, j);
    for (arma::uword i = 1; i < n; ++i) {
      minimum = std::min(minimum, residual(i, j));
      maximum = std::max(maximum, residual(i, j));
    }
    if (!(maximum > minimum)) {
      Rcpp::stop(
        "Feng-Sun requires positive variation in every variable; variable "
        "%llu is constant across observations after null centering.",
        static_cast<unsigned long long>(j + 1)
      );
    }
  }

  return FsPreparedData{
    residual,
    anchor,
    ratio,
    log_scale,
    fs_stable_column_mean(x),
    fallback_columns
  };
}

bool fs_direction(const arma::mat& data,
                  arma::uword row,
                  const arma::vec& theta,
                  const arma::vec& diagonal,
                  arma::vec& direction,
                  double& radius) {
  const arma::uword p = data.n_cols;
  double max_abs = 0.0;
  for (arma::uword j = 0; j < p; ++j) {
    const double value =
      (data(row, j) - theta(j)) / std::sqrt(diagonal(j));
    if (!std::isfinite(value)) {
      Rcpp::stop(
        "Feng-Sun diagonal standardisation produced a non-finite value."
      );
    }
    direction(j) = value;
    max_abs = std::max(max_abs, std::abs(value));
  }
  if (max_abs == 0.0) {
    direction.zeros();
    radius = 0.0;
    return false;
  }

  long double sum_squares = 0.0L;
  long double correction = 0.0L;
  for (arma::uword j = 0; j < p; ++j) {
    direction(j) /= max_abs;
    const long double value = static_cast<long double>(direction(j));
    fs_neumaier_add(value * value, sum_squares, correction);
  }
  const long double norm = std::sqrt(sum_squares + correction);
  if (!(norm > 0.0L) || !std::isfinite(norm)) {
    Rcpp::stop("Feng-Sun could not normalise a nonzero spatial direction.");
  }
  radius = max_abs * static_cast<double>(norm);
  if (!(radius > 0.0) || !std::isfinite(radius)) {
    Rcpp::stop("Feng-Sun produced a non-finite standardised radius.");
  }
  direction /= static_cast<double>(norm);
  return true;
}

long double fs_dot(const arma::vec& x, const arma::vec& y) {
  long double total = 0.0L;
  long double correction = 0.0L;
  for (arma::uword j = 0; j < x.n_elem; ++j) {
    fs_neumaier_add(
      static_cast<long double>(x(j)) * static_cast<long double>(y(j)),
      total,
      correction
    );
  }
  return total + correction;
}

struct FsScore {
  arma::vec sign_sum;
  arma::vec sign_square_sum;
  double inverse_radius_sum;
  double location_residual;
  double scale_residual;
};

FsScore fs_leaveout_score(const arma::mat& data,
                          arma::uword leave_i,
                          arma::uword leave_j,
                          const arma::vec& theta,
                          const arma::vec& diagonal) {
  const arma::uword n = data.n_rows;
  const arma::uword p = data.n_cols;
  const arma::uword m = n - 2;
  arma::vec sign_sum(p, arma::fill::zeros);
  arma::vec sign_square_sum(p, arma::fill::zeros);
  arma::vec direction(p, arma::fill::zeros);
  long double inverse_radius_total = 0.0L;
  long double inverse_radius_correction = 0.0L;

  for (arma::uword row = 0; row < n; ++row) {
    if (row == leave_i || row == leave_j) {
      continue;
    }
    double radius = 0.0;
    if (!fs_direction(data, row, theta, diagonal, direction, radius)) {
      Rcpp::stop(
        "Feng-Sun diagonal HR iteration for leave-out pair (%llu, %llu) "
        "has a training observation exactly at its current location; the "
        "inverse-radius update is undefined.",
        static_cast<unsigned long long>(leave_i + 1),
        static_cast<unsigned long long>(leave_j + 1)
      );
    }
    sign_sum += direction;
    sign_square_sum += arma::square(direction);
    fs_neumaier_add(
      1.0L / static_cast<long double>(radius),
      inverse_radius_total,
      inverse_radius_correction
    );
  }

  const double inverse_radius_sum = static_cast<double>(
    inverse_radius_total + inverse_radius_correction
  );
  if (!(inverse_radius_sum > 0.0) || !std::isfinite(inverse_radius_sum)) {
    Rcpp::stop(
      "Feng-Sun diagonal HR iteration has a non-finite inverse-radius "
      "denominator for leave-out pair (%llu, %llu).",
      static_cast<unsigned long long>(leave_i + 1),
      static_cast<unsigned long long>(leave_j + 1)
    );
  }

  const double m_double = static_cast<double>(m);
  const double p_double = static_cast<double>(p);
  const arma::vec mean_sign = sign_sum / m_double;
  const arma::vec scale_equation =
    p_double * sign_square_sum / m_double;
  const double location_residual = arma::max(arma::abs(mean_sign));
  const double scale_residual =
    arma::max(arma::abs(scale_equation - 1.0));

  return FsScore{
    sign_sum,
    sign_square_sum,
    inverse_radius_sum,
    location_residual,
    scale_residual
  };
}

struct FsFit {
  arma::vec theta;
  arma::vec diagonal;
  int iterations;
  double location_residual;
  double scale_residual;
};

FsFit fs_fit_leave_two_out(const arma::mat& data,
                           arma::uword leave_i,
                           arma::uword leave_j,
                           double tolerance,
                           int max_iterations) {
  const arma::uword n = data.n_rows;
  const arma::uword p = data.n_cols;
  const arma::uword m = n - 2;
  const double m_double = static_cast<double>(m);
  const double p_double = static_cast<double>(p);

  arma::vec theta(p, arma::fill::zeros);
  for (arma::uword col = 0; col < p; ++col) {
    long double total = 0.0L;
    long double correction = 0.0L;
    for (arma::uword row = 0; row < n; ++row) {
      if (row != leave_i && row != leave_j) {
        fs_neumaier_add(
          static_cast<long double>(data(row, col)), total, correction
        );
      }
    }
    theta(col) = static_cast<double>(
      (total + correction) / static_cast<long double>(m)
    );
  }

  arma::vec diagonal(p, arma::fill::zeros);
  for (arma::uword col = 0; col < p; ++col) {
    long double total = 0.0L;
    long double correction = 0.0L;
    for (arma::uword row = 0; row < n; ++row) {
      if (row != leave_i && row != leave_j) {
        const long double value =
          static_cast<long double>(data(row, col)) -
          static_cast<long double>(theta(col));
        fs_neumaier_add(value * value, total, correction);
      }
    }
    diagonal(col) = static_cast<double>(
      (total + correction) / static_cast<long double>(m - 1)
    );
    if (!(diagonal(col) > 0.0) || !std::isfinite(diagonal(col))) {
      Rcpp::stop(
        "Feng-Sun requires every leave-two-out initial marginal variance "
        "to be finite and strictly positive; variable %llu is degenerate "
        "for pair (%llu, %llu). No ridge is applied.",
        static_cast<unsigned long long>(col + 1),
        static_cast<unsigned long long>(leave_i + 1),
        static_cast<unsigned long long>(leave_j + 1)
      );
    }
  }
  diagonal /= diagonal.max();

  for (int iteration = 0; iteration <= max_iterations; ++iteration) {
    const FsScore score = fs_leaveout_score(
      data, leave_i, leave_j, theta, diagonal
    );
    if (score.location_residual <= tolerance &&
        score.scale_residual <= tolerance) {
      return FsFit{
        theta,
        diagonal,
        iteration,
        score.location_residual,
        score.scale_residual
      };
    }
    if (iteration == max_iterations) {
      Rcpp::stop(
        "Feng-Sun diagonal HR fit did not converge for leave-out pair "
        "(%llu, %llu) in %d iterations (location residual %.6g, scale "
        "residual %.6g).",
        static_cast<unsigned long long>(leave_i + 1),
        static_cast<unsigned long long>(leave_j + 1),
        max_iterations,
        score.location_residual,
        score.scale_residual
      );
    }

    arma::vec theta_new = theta;
    for (arma::uword col = 0; col < p; ++col) {
      theta_new(col) += std::sqrt(diagonal(col)) *
        score.sign_sum(col) / score.inverse_radius_sum;
    }
    arma::vec diagonal_new = diagonal %
      (p_double * score.sign_square_sum / m_double);
    if (!theta_new.is_finite() || !diagonal_new.is_finite() ||
        arma::any(diagonal_new <= 0.0)) {
      Rcpp::stop(
        "Feng-Sun diagonal HR update became non-finite or non-positive "
        "for leave-out pair (%llu, %llu); no ridge or floor is applied.",
        static_cast<unsigned long long>(leave_i + 1),
        static_cast<unsigned long long>(leave_j + 1)
      );
    }
    diagonal_new /= diagonal_new.max();
    theta = std::move(theta_new);
    diagonal = std::move(diagonal_new);
  }

  Rcpp::stop("Internal Feng-Sun iteration error.");
}

arma::vec fs_restore_location(const arma::vec& theta_standardised,
                              const FsPreparedData& prepared,
                              const arma::vec& mu) {
  arma::vec restored(theta_standardised.n_elem, arma::fill::zeros);
  for (arma::uword j = 0; j < theta_standardised.n_elem; ++j) {
    const double inner = mu(j) / prepared.anchor(j) +
      prepared.residual_scale_ratio(j) * theta_standardised(j);
    restored(j) = prepared.anchor(j) * inner;
    if (!std::isfinite(restored(j))) {
      Rcpp::stop(
        "Feng-Sun leave-two-out location in variable %llu is outside the "
        "finite input range.",
        static_cast<unsigned long long>(j + 1)
      );
    }
  }
  return restored;
}

}  // namespace


//' Feng--Sun one-sample spatial-sign statistic kernel
//'
//' @param x Numeric observation-by-variable matrix.
//' @param mu Numeric null-location vector.
//' @param tolerance Positive equation-residual tolerance.
//' @param max_iterations Positive maximum update count.
//' @return Internal list of feasible statistic components and leave-out fits.
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List cpp_feng_sun_one_sample(const arma::mat& x,
                                   const arma::vec& mu,
                                   double tolerance,
                                   int max_iterations) {
  if (!x.is_finite()) {
    Rcpp::stop("`x` must contain only finite values.");
  }
  if (!mu.is_finite()) {
    Rcpp::stop("`mu` must contain only finite values.");
  }
  const arma::uword n = x.n_rows;
  const arma::uword p = x.n_cols;
  if (n < 4) {
    Rcpp::stop(
      "Feng-Sun requires at least four observations so every leave-two-out "
      "sample has at least two rows."
    );
  }
  if (p < 1 || mu.n_elem != p) {
    Rcpp::stop("`mu` must have one entry for every column of `x`.");
  }
  if (!std::isfinite(tolerance) || tolerance <= 0.0) {
    Rcpp::stop("`tolerance` must be finite and strictly positive.");
  }
  if (max_iterations < 1) {
    Rcpp::stop("`max_iterations` must be a positive integer.");
  }

  const FsPreparedData prepared = fs_prepare_data(x, mu);
  const arma::mat& data = prepared.residual;
  const arma::uword pair_count = n * (n - 1) / 2;
  arma::mat leaveout_theta_standardised(pair_count, p, arma::fill::zeros);
  arma::mat leaveout_theta(pair_count, p, arma::fill::zeros);
  arma::mat leaveout_diagonal_standardised(pair_count, p, arma::fill::zeros);
  arma::mat leaveout_log_diagonal_input(pair_count, p, arma::fill::zeros);
  arma::mat leaveout_diagonal_input(pair_count, p, arma::fill::zeros);
  arma::vec pair_i(pair_count, arma::fill::zeros);
  arma::vec pair_j(pair_count, arma::fill::zeros);
  arma::vec pair_test_inner(pair_count, arma::fill::zeros);
  arma::vec pair_variance_inner_squared(pair_count, arma::fill::zeros);
  arma::vec pair_iterations(pair_count, arma::fill::zeros);
  arma::vec pair_location_residual(pair_count, arma::fill::zeros);
  arma::vec pair_scale_residual(pair_count, arma::fill::zeros);
  arma::vec pair_test_zero_signs(pair_count, arma::fill::zeros);
  arma::vec pair_variance_zero_signs(pair_count, arma::fill::zeros);

  arma::vec zero_theta(p, arma::fill::zeros);
  arma::vec direction_i(p, arma::fill::zeros);
  arma::vec direction_j(p, arma::fill::zeros);
  long double test_total = 0.0L;
  long double test_correction = 0.0L;
  long double variance_square_total = 0.0L;
  long double variance_square_correction = 0.0L;
  arma::uword pair_index = 0;

  for (arma::uword i = 0; i + 1 < n; ++i) {
    for (arma::uword j = i + 1; j < n; ++j, ++pair_index) {
      const FsFit fit = fs_fit_leave_two_out(
        data, i, j, tolerance, max_iterations
      );
      pair_i(pair_index) = static_cast<double>(i + 1);
      pair_j(pair_index) = static_cast<double>(j + 1);
      pair_iterations(pair_index) = static_cast<double>(fit.iterations);
      pair_location_residual(pair_index) = fit.location_residual;
      pair_scale_residual(pair_index) = fit.scale_residual;
      leaveout_theta_standardised.row(pair_index) = fit.theta.t();
      leaveout_diagonal_standardised.row(pair_index) = fit.diagonal.t();
      leaveout_theta.row(pair_index) =
        fs_restore_location(fit.theta, prepared, mu).t();

      double maximum_log_diagonal = -std::numeric_limits<double>::infinity();
      for (arma::uword col = 0; col < p; ++col) {
        const double value = 2.0 * prepared.log_residual_scale(col) +
          std::log(fit.diagonal(col));
        leaveout_log_diagonal_input(pair_index, col) = value;
        maximum_log_diagonal = std::max(maximum_log_diagonal, value);
      }
      for (arma::uword col = 0; col < p; ++col) {
        const double canonical_log =
          leaveout_log_diagonal_input(pair_index, col) -
          maximum_log_diagonal;
        leaveout_log_diagonal_input(pair_index, col) = canonical_log;
        leaveout_diagonal_input(pair_index, col) =
          std::exp(canonical_log);
      }

      double radius_i = 0.0;
      double radius_j = 0.0;
      const bool test_nonzero_i = fs_direction(
        data, i, zero_theta, fit.diagonal, direction_i, radius_i
      );
      const bool test_nonzero_j = fs_direction(
        data, j, zero_theta, fit.diagonal, direction_j, radius_j
      );
      pair_test_zero_signs(pair_index) =
        static_cast<double>((!test_nonzero_i) + (!test_nonzero_j));
      const long double test_inner = fs_dot(direction_i, direction_j);
      pair_test_inner(pair_index) = static_cast<double>(test_inner);
      fs_neumaier_add(test_inner, test_total, test_correction);

      const bool variance_nonzero_i = fs_direction(
        data, i, fit.theta, fit.diagonal, direction_i, radius_i
      );
      const bool variance_nonzero_j = fs_direction(
        data, j, fit.theta, fit.diagonal, direction_j, radius_j
      );
      pair_variance_zero_signs(pair_index) =
        static_cast<double>((!variance_nonzero_i) + (!variance_nonzero_j));
      const long double variance_inner = fs_dot(direction_i, direction_j);
      const long double variance_inner_squared =
        variance_inner * variance_inner;
      pair_variance_inner_squared(pair_index) =
        static_cast<double>(variance_inner_squared);
      fs_neumaier_add(
        variance_inner_squared,
        variance_square_total,
        variance_square_correction
      );
    }
  }

  const long double pair_count_ld = static_cast<long double>(pair_count);
  const long double p_ld = static_cast<long double>(p);
  const long double test_sum_ld = test_total + test_correction;
  const long double variance_square_sum_ld =
    variance_square_total + variance_square_correction;
  const long double statistic_raw_ld = test_sum_ld / pair_count_ld;
  const long double trace_r2_hat_ld =
    p_ld * p_ld * variance_square_sum_ld / pair_count_ld;
  const long double sigma2_hat_ld =
    variance_square_sum_ld / (pair_count_ld * pair_count_ld);

  const double statistic_raw = static_cast<double>(statistic_raw_ld);
  const double trace_r2_hat = static_cast<double>(trace_r2_hat_ld);
  const double sigma2_hat = static_cast<double>(sigma2_hat_ld);
  if (!std::isfinite(trace_r2_hat) || trace_r2_hat <= 0.0 ||
      !std::isfinite(sigma2_hat) || sigma2_hat <= 0.0) {
    Rcpp::stop(
      "Feng-Sun requires its feasible leave-two-out variance estimate to "
      "be finite and strictly positive; no absolute-value or variance-floor "
      "repair is applied."
    );
  }
  const double sigma_hat = std::sqrt(sigma2_hat);
  const double z = statistic_raw / sigma_hat;
  if (!std::isfinite(statistic_raw) || !std::isfinite(z)) {
    Rcpp::stop("Feng-Sun produced a non-finite test statistic.");
  }

  return Rcpp::List::create(
    Rcpp::Named("z") = z,
    Rcpp::Named("T_SS") = statistic_raw,
    Rcpp::Named("test_inner_product_sum") =
      static_cast<double>(test_sum_ld),
    Rcpp::Named("variance_inner_product_squared_sum") =
      static_cast<double>(variance_square_sum_ld),
    Rcpp::Named("trace_R2_hat") = trace_r2_hat,
    Rcpp::Named("sigma2_hat") = sigma2_hat,
    Rcpp::Named("sigma_hat") = sigma_hat,
    Rcpp::Named("pair_i") = pair_i,
    Rcpp::Named("pair_j") = pair_j,
    Rcpp::Named("pair_test_inner_product") = pair_test_inner,
    Rcpp::Named("pair_variance_inner_product_squared") =
      pair_variance_inner_squared,
    Rcpp::Named("pair_iterations") = pair_iterations,
    Rcpp::Named("pair_location_residual") = pair_location_residual,
    Rcpp::Named("pair_scale_residual") = pair_scale_residual,
    Rcpp::Named("pair_test_zero_signs") = pair_test_zero_signs,
    Rcpp::Named("pair_variance_zero_signs") = pair_variance_zero_signs,
    Rcpp::Named("leaveout_theta_standardised") =
      leaveout_theta_standardised,
    Rcpp::Named("leaveout_theta") = leaveout_theta,
    Rcpp::Named("leaveout_D_standardised") =
      leaveout_diagonal_standardised,
    Rcpp::Named("leaveout_log_D_input_canonical") =
      leaveout_log_diagonal_input,
    Rcpp::Named("leaveout_D_input_canonical") =
      leaveout_diagonal_input,
    Rcpp::Named("sample_mean") = prepared.sample_mean,
    Rcpp::Named("column_log_residual_scale") =
      prepared.log_residual_scale,
    Rcpp::Named("subtraction_overflow_columns") =
      static_cast<double>(prepared.subtraction_overflow_columns),
    Rcpp::Named("pair_count") = static_cast<double>(pair_count),
    Rcpp::Named("ordered_pair_count") =
      static_cast<double>(n) * static_cast<double>(n - 1),
    Rcpp::Named("n") = static_cast<double>(n),
    Rcpp::Named("p") = static_cast<double>(p),
    Rcpp::Named("tolerance") = tolerance,
    Rcpp::Named("max_iterations") = max_iterations
  );
}
