// Feng--Zou--Wang two-sample multivariate-sign test.
//
// This is a literal implementation of the leave-one-out statistic and the
// feasible null-variance estimator in Proposition 2 of Feng, Zou and Wang
// (2016).  In particular, the numerator uses crossed leave-one-out locations,
// while the trace estimator uses own-sample leave-one-out signs, inverse
// radii, and the two full-sample diagonal scale fits.

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp17)]]

#include <RcppArmadillo.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <string>
#include <utility>

namespace {

void fzw_neumaier_add(long double value,
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

long double fzw_dot(const arma::vec& x, const arma::vec& y) {
  long double total = 0.0L;
  long double correction = 0.0L;
  for (arma::uword j = 0; j < x.n_elem; ++j) {
    fzw_neumaier_add(
      static_cast<long double>(x(j)) * static_cast<long double>(y(j)),
      total,
      correction
    );
  }
  return total + correction;
}

arma::vec fzw_stable_column_mean(const arma::mat& x) {
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
      fzw_neumaier_add(
        static_cast<long double>(x(i, j) / scale), total, correction
      );
    }
    output(j) = static_cast<double>(
      (total + correction) / static_cast<long double>(x.n_rows)
    ) * scale;
  }
  return output;
}

struct FzwPreparedData {
  arma::mat x;
  arma::mat y;
  arma::vec anchor;
  arma::vec scale_base;
  arma::vec scale_ratio;
  arma::vec log_scale;
  arma::vec mean_x;
  arma::vec mean_y;
  arma::uword subtraction_overflow_columns;
};

FzwPreparedData fzw_prepare_data(const arma::mat& x, const arma::mat& y) {
  const arma::uword n1 = x.n_rows;
  const arma::uword n2 = y.n_rows;
  const arma::uword p = x.n_cols;
  arma::mat x_out(n1, p, arma::fill::zeros);
  arma::mat y_out(n2, p, arma::fill::zeros);
  arma::vec anchor(p, arma::fill::zeros);
  arma::vec scale_base(p, arma::fill::ones);
  arma::vec scale_ratio(p, arma::fill::ones);
  arma::vec log_scale(p, arma::fill::zeros);
  arma::uword overflow_columns = 0;

  for (arma::uword j = 0; j < p; ++j) {
    anchor(j) = x(0, j);
    bool overflow = false;
    double max_difference = 0.0;
    for (arma::uword i = 0; i < n1; ++i) {
      const double value = x(i, j) - anchor(j);
      x_out(i, j) = value;
      if (!std::isfinite(value)) {
        overflow = true;
      } else {
        max_difference = std::max(max_difference, std::abs(value));
      }
    }
    for (arma::uword i = 0; i < n2; ++i) {
      const double value = y(i, j) - anchor(j);
      y_out(i, j) = value;
      if (!std::isfinite(value)) {
        overflow = true;
      } else {
        max_difference = std::max(max_difference, std::abs(value));
      }
    }

    if (overflow) {
      ++overflow_columns;
      double operand_scale = std::abs(anchor(j));
      for (arma::uword i = 0; i < n1; ++i) {
        operand_scale = std::max(operand_scale, std::abs(x(i, j)));
      }
      for (arma::uword i = 0; i < n2; ++i) {
        operand_scale = std::max(operand_scale, std::abs(y(i, j)));
      }
      if (!(operand_scale > 0.0) || !std::isfinite(operand_scale)) {
        Rcpp::stop(
          "Feng-Zou-Wang could not construct a finite common scale for "
          "variable %llu.",
          static_cast<unsigned long long>(j + 1)
        );
      }
      max_difference = 0.0;
      const double anchor_scaled = anchor(j) / operand_scale;
      for (arma::uword i = 0; i < n1; ++i) {
        const double value = x(i, j) / operand_scale - anchor_scaled;
        x_out(i, j) = value;
        max_difference = std::max(max_difference, std::abs(value));
      }
      for (arma::uword i = 0; i < n2; ++i) {
        const double value = y(i, j) / operand_scale - anchor_scaled;
        y_out(i, j) = value;
        max_difference = std::max(max_difference, std::abs(value));
      }
      scale_base(j) = operand_scale;
      scale_ratio(j) = max_difference;
    } else {
      scale_base(j) = max_difference;
      scale_ratio(j) = 1.0;
    }

    if (!(max_difference > 0.0) || !std::isfinite(max_difference) ||
        !(scale_base(j) > 0.0) || !std::isfinite(scale_base(j)) ||
        !(scale_ratio(j) > 0.0) || !std::isfinite(scale_ratio(j))) {
      Rcpp::stop(
        "Feng-Zou-Wang requires pooled positive variation in variable "
        "%llu.",
        static_cast<unsigned long long>(j + 1)
      );
    }
    x_out.col(j) /= max_difference;
    y_out.col(j) /= max_difference;
    log_scale(j) = std::log(scale_base(j)) + std::log(scale_ratio(j));
  }

  return FzwPreparedData{
    x_out,
    y_out,
    anchor,
    scale_base,
    scale_ratio,
    log_scale,
    fzw_stable_column_mean(x),
    fzw_stable_column_mean(y),
    overflow_columns
  };
}

std::string fzw_fit_label(int excluded) {
  if (excluded < 0) {
    return "full sample";
  }
  return "leave-one-out observation " + std::to_string(excluded + 1);
}

bool fzw_direction_difference(const arma::vec& value,
                              const arma::vec& theta,
                              const arma::vec& diagonal,
                              arma::vec& direction,
                              double& radius) {
  const arma::uword p = value.n_elem;
  double max_abs = 0.0;
  for (arma::uword j = 0; j < p; ++j) {
    const double standardized =
      (value(j) - theta(j)) / std::sqrt(diagonal(j));
    if (!std::isfinite(standardized)) {
      Rcpp::stop(
        "Feng-Zou-Wang diagonal standardisation produced a non-finite "
        "value; no numerical floor is applied."
      );
    }
    direction(j) = standardized;
    max_abs = std::max(max_abs, std::abs(standardized));
  }
  if (max_abs == 0.0) {
    direction.zeros();
    radius = 0.0;
    return false;
  }

  long double square_total = 0.0L;
  long double square_correction = 0.0L;
  for (arma::uword j = 0; j < p; ++j) {
    direction(j) /= max_abs;
    const long double component = static_cast<long double>(direction(j));
    fzw_neumaier_add(
      component * component, square_total, square_correction
    );
  }
  const long double norm = std::sqrt(square_total + square_correction);
  if (!(norm > 0.0L) || !std::isfinite(norm)) {
    Rcpp::stop(
      "Feng-Zou-Wang could not normalise a nonzero spatial direction."
    );
  }
  radius = max_abs * static_cast<double>(norm);
  if (!(radius > 0.0) || !std::isfinite(radius)) {
    Rcpp::stop(
      "Feng-Zou-Wang produced a non-finite standardised radius."
    );
  }
  direction /= static_cast<double>(norm);
  return true;
}

bool fzw_direction_row(const arma::mat& data,
                       arma::uword row,
                       const arma::vec& theta,
                       const arma::vec& diagonal,
                       arma::vec& direction,
                       double& radius) {
  return fzw_direction_difference(
    data.row(row).t(), theta, diagonal, direction, radius
  );
}

struct FzwScore {
  arma::vec sign_sum;
  arma::vec sign_square_sum;
  double inverse_radius_sum;
  double minimum_radius;
  double location_residual;
  double scale_residual;
};

FzwScore fzw_score(const arma::mat& data,
                   int excluded,
                   const arma::vec& theta,
                   const arma::vec& diagonal,
                   int group) {
  const arma::uword n = data.n_rows;
  const arma::uword p = data.n_cols;
  const arma::uword m = n - (excluded >= 0 ? 1 : 0);
  arma::vec sign_sum(p, arma::fill::zeros);
  arma::vec sign_square_sum(p, arma::fill::zeros);
  arma::vec direction(p, arma::fill::zeros);
  long double inverse_total = 0.0L;
  long double inverse_correction = 0.0L;
  double minimum_radius = std::numeric_limits<double>::infinity();

  for (arma::uword row = 0; row < n; ++row) {
    if (excluded >= 0 && row == static_cast<arma::uword>(excluded)) {
      continue;
    }
    double radius = 0.0;
    if (!fzw_direction_row(
          data, row, theta, diagonal, direction, radius
        )) {
      Rcpp::stop(
        "Feng-Zou-Wang diagonal HR iteration for group %d (%s) has a "
        "training observation exactly at its current location; the "
        "inverse-radius location update is undefined.",
        group,
        fzw_fit_label(excluded).c_str()
      );
    }
    sign_sum += direction;
    sign_square_sum += arma::square(direction);
    fzw_neumaier_add(
      1.0L / static_cast<long double>(radius),
      inverse_total,
      inverse_correction
    );
    minimum_radius = std::min(minimum_radius, radius);
  }

  const double inverse_radius_sum = static_cast<double>(
    inverse_total + inverse_correction
  );
  if (!(inverse_radius_sum > 0.0) ||
      !std::isfinite(inverse_radius_sum)) {
    Rcpp::stop(
      "Feng-Zou-Wang diagonal HR iteration for group %d (%s) has an "
      "invalid inverse-radius denominator.",
      group,
      fzw_fit_label(excluded).c_str()
    );
  }
  const double m_double = static_cast<double>(m);
  const double p_double = static_cast<double>(p);
  const double location_residual =
    arma::max(arma::abs(sign_sum / m_double));
  const double scale_residual = arma::max(arma::abs(
    p_double * sign_square_sum / m_double - 1.0
  ));
  return FzwScore{
    sign_sum,
    sign_square_sum,
    inverse_radius_sum,
    minimum_radius,
    location_residual,
    scale_residual
  };
}

struct FzwFit {
  arma::vec theta;
  arma::vec diagonal;
  int iterations;
  double location_residual;
  double scale_residual;
  double minimum_training_radius;
};

void fzw_normalise_diagonal(arma::vec& diagonal,
                            int group,
                            int excluded) {
  long double log_total = 0.0L;
  long double log_correction = 0.0L;
  for (arma::uword j = 0; j < diagonal.n_elem; ++j) {
    if (!(diagonal(j) > 0.0) || !std::isfinite(diagonal(j))) {
      Rcpp::stop(
        "Feng-Zou-Wang cannot normalise a non-positive or non-finite "
        "diagonal for group %d (%s).",
        group,
        fzw_fit_label(excluded).c_str()
      );
    }
    fzw_neumaier_add(
      std::log(static_cast<long double>(diagonal(j))),
      log_total,
      log_correction
    );
  }
  const long double mean_log =
    (log_total + log_correction) /
    static_cast<long double>(diagonal.n_elem);
  for (arma::uword j = 0; j < diagonal.n_elem; ++j) {
    const long double centred_log =
      std::log(static_cast<long double>(diagonal(j))) - mean_log;
    diagonal(j) = static_cast<double>(std::exp(centred_log));
    if (!(diagonal(j) > 0.0) || !std::isfinite(diagonal(j))) {
      Rcpp::stop(
        "Feng-Zou-Wang geometric-mean scale identification overflowed "
        "for group %d (%s); no floor is applied.",
        group,
        fzw_fit_label(excluded).c_str()
      );
    }
  }
}

FzwFit fzw_fit(const arma::mat& data,
               int excluded,
               double tolerance,
               int max_iterations,
               int group) {
  const arma::uword n = data.n_rows;
  const arma::uword p = data.n_cols;
  const arma::uword m = n - (excluded >= 0 ? 1 : 0);
  const double m_double = static_cast<double>(m);
  const double p_double = static_cast<double>(p);

  arma::vec theta(p, arma::fill::zeros);
  for (arma::uword col = 0; col < p; ++col) {
    long double total = 0.0L;
    long double correction = 0.0L;
    for (arma::uword row = 0; row < n; ++row) {
      if (excluded < 0 || row != static_cast<arma::uword>(excluded)) {
        fzw_neumaier_add(
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
      if (excluded < 0 || row != static_cast<arma::uword>(excluded)) {
        const long double centered =
          static_cast<long double>(data(row, col)) -
          static_cast<long double>(theta(col));
        fzw_neumaier_add(centered * centered, total, correction);
      }
    }
    diagonal(col) = static_cast<double>(
      (total + correction) / static_cast<long double>(m - 1)
    );
    if (!(diagonal(col) > 0.0) || !std::isfinite(diagonal(col))) {
      Rcpp::stop(
        "Feng-Zou-Wang requires every initial marginal variance to be "
        "finite and strictly positive; group %d variable %llu is "
        "degenerate for %s. No ridge is applied.",
        group,
        static_cast<unsigned long long>(col + 1),
        fzw_fit_label(excluded).c_str()
      );
    }
  }
  // Equation (3) identifies D only up to a common multiplier.  A unit
  // geometric mean is used for every fit.  Unlike max-diagonal
  // normalisation, this identification is equivariant with the same global
  // multiplier for all fits under a common coordinatewise change of units;
  // that coherence is essential for c-hat and the full-sample D bridges.
  fzw_normalise_diagonal(diagonal, group, excluded);

  for (int iteration = 0; iteration <= max_iterations; ++iteration) {
    const FzwScore score = fzw_score(
      data, excluded, theta, diagonal, group
    );
    if (score.location_residual <= tolerance &&
        score.scale_residual <= tolerance) {
      return FzwFit{
        theta,
        diagonal,
        iteration,
        score.location_residual,
        score.scale_residual,
        score.minimum_radius
      };
    }
    if (iteration == max_iterations) {
      Rcpp::stop(
        "Feng-Zou-Wang diagonal HR fit did not converge for group %d "
        "(%s) in %d iterations (location residual %.6g, scale residual "
        "%.6g).",
        group,
        fzw_fit_label(excluded).c_str(),
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
        "Feng-Zou-Wang diagonal HR update became non-finite or "
        "non-positive for group %d (%s); no ridge, absolute-value, or "
        "floor repair is applied.",
        group,
        fzw_fit_label(excluded).c_str()
      );
    }
    fzw_normalise_diagonal(diagonal_new, group, excluded);
    theta = std::move(theta_new);
    diagonal = std::move(diagonal_new);
  }
  Rcpp::stop("Internal Feng-Zou-Wang iteration error.");
}

arma::vec fzw_restore_location(const arma::vec& theta,
                               const FzwPreparedData& prepared) {
  arma::vec output(theta.n_elem, arma::fill::zeros);
  for (arma::uword j = 0; j < theta.n_elem; ++j) {
    const long double anchored =
      static_cast<long double>(prepared.anchor(j)) /
        static_cast<long double>(prepared.scale_base(j)) +
      static_cast<long double>(prepared.scale_ratio(j)) *
        static_cast<long double>(theta(j));
    const long double restored =
      static_cast<long double>(prepared.scale_base(j)) * anchored;
    output(j) = static_cast<double>(restored);
    if (!std::isfinite(output(j))) {
      Rcpp::stop(
        "Feng-Zou-Wang fitted location in variable %llu lies outside "
        "the finite input range.",
        static_cast<unsigned long long>(j + 1)
      );
    }
  }
  return output;
}

void fzw_store_scale(const arma::vec& diagonal,
                     const FzwPreparedData& prepared,
                     arma::rowvec& canonical,
                     arma::rowvec& canonical_log) {
  double maximum_log = -std::numeric_limits<double>::infinity();
  for (arma::uword j = 0; j < diagonal.n_elem; ++j) {
    canonical_log(j) = 2.0 * prepared.log_scale(j) +
      std::log(diagonal(j));
    maximum_log = std::max(maximum_log, canonical_log(j));
  }
  for (arma::uword j = 0; j < diagonal.n_elem; ++j) {
    canonical_log(j) -= maximum_log;
    canonical(j) = std::exp(canonical_log(j));
  }
}

struct FzwGroupFits {
  FzwFit full;
  arma::mat leave_theta;
  arma::mat leave_diagonal;
  arma::vec iterations;
  arma::vec location_residual;
  arma::vec scale_residual;
  arma::vec minimum_training_radius;
  arma::mat own_direction;
  arma::vec own_radius;
  arma::vec own_inverse_radius;
};

FzwGroupFits fzw_fit_group(const arma::mat& data,
                           double tolerance,
                           int max_iterations,
                           int group) {
  const arma::uword n = data.n_rows;
  const arma::uword p = data.n_cols;
  FzwFit full = fzw_fit(
    data, -1, tolerance, max_iterations, group
  );
  arma::mat leave_theta(n, p, arma::fill::zeros);
  arma::mat leave_diagonal(n, p, arma::fill::zeros);
  arma::vec iterations(n, arma::fill::zeros);
  arma::vec location_residual(n, arma::fill::zeros);
  arma::vec scale_residual(n, arma::fill::zeros);
  arma::vec minimum_training_radius(n, arma::fill::zeros);
  arma::mat own_direction(n, p, arma::fill::zeros);
  arma::vec own_radius(n, arma::fill::zeros);
  arma::vec own_inverse_radius(n, arma::fill::zeros);
  arma::vec direction(p, arma::fill::zeros);

  for (arma::uword i = 0; i < n; ++i) {
    const FzwFit fit = fzw_fit(
      data, static_cast<int>(i), tolerance, max_iterations, group
    );
    leave_theta.row(i) = fit.theta.t();
    leave_diagonal.row(i) = fit.diagonal.t();
    iterations(i) = static_cast<double>(fit.iterations);
    location_residual(i) = fit.location_residual;
    scale_residual(i) = fit.scale_residual;
    minimum_training_radius(i) = fit.minimum_training_radius;
    double radius = 0.0;
    if (!fzw_direction_row(
          data, i, fit.theta, fit.diagonal, direction, radius
        )) {
      Rcpp::stop(
        "Feng-Zou-Wang leave-one-out standardised radius is zero for "
        "group %d observation %llu; the feasible inverse-radius estimate "
        "of c_%d is undefined. U(0)=0 is not a repair for this divisor.",
        group,
        static_cast<unsigned long long>(i + 1),
        group
      );
    }
    own_direction.row(i) = direction.t();
    own_radius(i) = radius;
    own_inverse_radius(i) = 1.0 / radius;
    if (!std::isfinite(own_inverse_radius(i))) {
      Rcpp::stop(
        "Feng-Zou-Wang produced a non-finite inverse radius for group %d "
        "observation %llu.",
        group,
        static_cast<unsigned long long>(i + 1)
      );
    }
  }

  return FzwGroupFits{
    full,
    leave_theta,
    leave_diagonal,
    iterations,
    location_residual,
    scale_residual,
    minimum_training_radius,
    own_direction,
    own_radius,
    own_inverse_radius
  };
}

}  // namespace


//' Feng--Zou--Wang two-sample multivariate-sign kernel
//'
//' @param x First numeric observation-by-variable matrix.
//' @param y Second numeric observation-by-variable matrix.
//' @param tolerance Positive estimating-equation tolerance.
//' @param max_iterations Positive maximum update count per fit.
//' @return Internal list containing the feasible statistic and diagnostics.
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List cpp_feng_zou_wang_two_sample_sign(
    const arma::mat& x,
    const arma::mat& y,
    double tolerance,
    int max_iterations) {
  if (!x.is_finite() || !y.is_finite()) {
    Rcpp::stop("`x` and `y` must contain only finite values.");
  }
  const arma::uword n1 = x.n_rows;
  const arma::uword n2 = y.n_rows;
  const arma::uword p = x.n_cols;
  if (n1 < 3 || n2 < 3) {
    Rcpp::stop(
      "Feng-Zou-Wang requires at least three observations in each group "
      "so every leave-one-out initial variance has at least two rows."
    );
  }
  if (p < 1 || y.n_cols != p) {
    Rcpp::stop("`x` and `y` must have the same positive number of columns.");
  }
  if (!std::isfinite(tolerance) || tolerance <= 0.0) {
    Rcpp::stop("`tolerance` must be finite and strictly positive.");
  }
  if (max_iterations < 1) {
    Rcpp::stop("`max_iterations` must be a positive integer.");
  }

  const FzwPreparedData prepared = fzw_prepare_data(x, y);
  const FzwGroupFits group1 = fzw_fit_group(
    prepared.x, tolerance, max_iterations, 1
  );
  const FzwGroupFits group2 = fzw_fit_group(
    prepared.y, tolerance, max_iterations, 2
  );

  arma::mat numerator_inner(n1, n2, arma::fill::zeros);
  arma::mat numerator_contribution(n1, n2, arma::fill::zeros);
  arma::mat numerator_zero_signs(n1, n2, arma::fill::zeros);
  arma::vec direction1(p, arma::fill::zeros);
  arma::vec direction2(p, arma::fill::zeros);
  long double numerator_total = 0.0L;
  long double numerator_correction = 0.0L;
  arma::uword numerator_zero_uses = 0;

  for (arma::uword i = 0; i < n1; ++i) {
    const arma::vec d1 = group1.leave_diagonal.row(i).t();
    const arma::vec theta1 = group1.leave_theta.row(i).t();
    for (arma::uword j = 0; j < n2; ++j) {
      const arma::vec d2 = group2.leave_diagonal.row(j).t();
      const arma::vec theta2 = group2.leave_theta.row(j).t();
      double radius1 = 0.0;
      double radius2 = 0.0;
      const bool nonzero1 = fzw_direction_row(
        prepared.x, i, theta2, d1, direction1, radius1
      );
      const bool nonzero2 = fzw_direction_row(
        prepared.y, j, theta1, d2, direction2, radius2
      );
      const arma::uword zero_count =
        static_cast<arma::uword>(!nonzero1) +
        static_cast<arma::uword>(!nonzero2);
      numerator_zero_signs(i, j) = static_cast<double>(zero_count);
      numerator_zero_uses += zero_count;
      const long double inner = fzw_dot(direction1, direction2);
      const long double contribution = -inner;
      numerator_inner(i, j) = static_cast<double>(inner);
      numerator_contribution(i, j) = static_cast<double>(contribution);
      fzw_neumaier_add(
        contribution, numerator_total, numerator_correction
      );
    }
  }

  const long double cross_count =
    static_cast<long double>(n1) * static_cast<long double>(n2);
  const long double R_n_ld =
    (numerator_total + numerator_correction) / cross_count;

  const double c1_hat = arma::mean(group1.own_inverse_radius);
  const double c2_hat = arma::mean(group2.own_inverse_radius);
  if (!(c1_hat > 0.0) || !std::isfinite(c1_hat) ||
      !(c2_hat > 0.0) || !std::isfinite(c2_hat)) {
    Rcpp::stop(
      "Feng-Zou-Wang requires finite positive leave-one-out inverse-radius "
      "means; no replacement or floor is applied."
    );
  }

  arma::vec bridge1(p, arma::fill::zeros);
  arma::vec bridge2(p, arma::fill::zeros);
  for (arma::uword j = 0; j < p; ++j) {
    bridge1(j) = std::sqrt(
      group1.full.diagonal(j) / group2.full.diagonal(j)
    );
    bridge2(j) = 1.0 / bridge1(j);
  }
  if (!bridge1.is_finite() || !bridge2.is_finite()) {
    Rcpp::stop(
      "Feng-Zou-Wang full-sample diagonal bridge is non-finite; no scale "
      "repair is applied."
    );
  }

  arma::mat trace1_pair(n1, n1, arma::fill::zeros);
  arma::mat trace2_pair(n2, n2, arma::fill::zeros);
  arma::mat trace3_pair(n1, n2, arma::fill::zeros);
  long double trace1_square_total = 0.0L;
  long double trace1_square_correction = 0.0L;
  long double trace2_square_total = 0.0L;
  long double trace2_square_correction = 0.0L;
  long double trace3_square_total = 0.0L;
  long double trace3_square_correction = 0.0L;

  for (arma::uword k = 0; k < n1; ++k) {
    const arma::vec right =
      bridge1 % group1.own_direction.row(k).t();
    for (arma::uword l = 0; l < n1; ++l) {
      if (l == k) {
        continue;
      }
      const long double inner = fzw_dot(
        group1.own_direction.row(l).t(), right
      );
      const long double squared = inner * inner;
      trace1_pair(l, k) = static_cast<double>(squared);
      fzw_neumaier_add(
        squared, trace1_square_total, trace1_square_correction
      );
    }
  }
  for (arma::uword k = 0; k < n2; ++k) {
    const arma::vec right =
      bridge2 % group2.own_direction.row(k).t();
    for (arma::uword l = 0; l < n2; ++l) {
      if (l == k) {
        continue;
      }
      const long double inner = fzw_dot(
        group2.own_direction.row(l).t(), right
      );
      const long double squared = inner * inner;
      trace2_pair(l, k) = static_cast<double>(squared);
      fzw_neumaier_add(
        squared, trace2_square_total, trace2_square_correction
      );
    }
  }
  for (arma::uword l = 0; l < n1; ++l) {
    for (arma::uword k = 0; k < n2; ++k) {
      const long double inner = fzw_dot(
        group1.own_direction.row(l).t(),
        group2.own_direction.row(k).t()
      );
      const long double squared = inner * inner;
      trace3_pair(l, k) = static_cast<double>(squared);
      fzw_neumaier_add(
        squared, trace3_square_total, trace3_square_correction
      );
    }
  }

  const long double p_ld = static_cast<long double>(p);
  const long double n1_ld = static_cast<long double>(n1);
  const long double n2_ld = static_cast<long double>(n2);
  const long double ordered1 = n1_ld * static_cast<long double>(n1 - 1);
  const long double ordered2 = n2_ld * static_cast<long double>(n2 - 1);
  const long double trace1_sum =
    trace1_square_total + trace1_square_correction;
  const long double trace2_sum =
    trace2_square_total + trace2_square_correction;
  const long double trace3_sum =
    trace3_square_total + trace3_square_correction;
  const long double c1_ld = static_cast<long double>(c1_hat);
  const long double c2_ld = static_cast<long double>(c2_hat);
  const long double trace_A1_ld = p_ld * p_ld *
    (c2_ld * c2_ld) / (c1_ld * c1_ld) * trace1_sum / ordered1;
  const long double trace_A2_ld = p_ld * p_ld *
    (c1_ld * c1_ld) / (c2_ld * c2_ld) * trace2_sum / ordered2;
  const long double trace_A3_ld =
    p_ld * p_ld * trace3_sum / cross_count;
  const long double variance_term1_ld =
    2.0L * trace_A1_ld / (ordered1 * p_ld * p_ld);
  const long double variance_term2_ld =
    2.0L * trace_A2_ld / (ordered2 * p_ld * p_ld);
  const long double variance_term3_ld =
    4.0L * trace_A3_ld / (cross_count * p_ld * p_ld);
  const long double sigma2_ld = variance_term1_ld +
    variance_term2_ld + variance_term3_ld;

  const double R_n = static_cast<double>(R_n_ld);
  const double trace_A1 = static_cast<double>(trace_A1_ld);
  const double trace_A2 = static_cast<double>(trace_A2_ld);
  const double trace_A3 = static_cast<double>(trace_A3_ld);
  const double variance_term1 = static_cast<double>(variance_term1_ld);
  const double variance_term2 = static_cast<double>(variance_term2_ld);
  const double variance_term3 = static_cast<double>(variance_term3_ld);
  const double sigma2_hat = static_cast<double>(sigma2_ld);
  if (!std::isfinite(trace_A1) || trace_A1 < 0.0 ||
      !std::isfinite(trace_A2) || trace_A2 < 0.0 ||
      !std::isfinite(trace_A3) || trace_A3 < 0.0 ||
      !(sigma2_hat > 0.0) || !std::isfinite(sigma2_hat)) {
    Rcpp::stop(
      "Feng-Zou-Wang requires its feasible Proposition 2 variance "
      "estimate to be finite and strictly positive; no absolute-value, "
      "ridge, or variance-floor repair is applied."
    );
  }
  const double sigma_hat = std::sqrt(sigma2_hat);
  const double z = R_n / sigma_hat;
  if (!std::isfinite(R_n) || !std::isfinite(z)) {
    Rcpp::stop("Feng-Zou-Wang produced a non-finite test statistic.");
  }

  arma::mat leave_theta1_input(n1, p, arma::fill::zeros);
  arma::mat leave_theta2_input(n2, p, arma::fill::zeros);
  arma::mat leave_D1_input(n1, p, arma::fill::zeros);
  arma::mat leave_D2_input(n2, p, arma::fill::zeros);
  arma::mat leave_log_D1_input(n1, p, arma::fill::zeros);
  arma::mat leave_log_D2_input(n2, p, arma::fill::zeros);
  for (arma::uword i = 0; i < n1; ++i) {
    leave_theta1_input.row(i) = fzw_restore_location(
      group1.leave_theta.row(i).t(), prepared
    ).t();
    arma::rowvec canonical(p, arma::fill::zeros);
    arma::rowvec canonical_log(p, arma::fill::zeros);
    fzw_store_scale(
      group1.leave_diagonal.row(i).t(),
      prepared,
      canonical,
      canonical_log
    );
    leave_D1_input.row(i) = canonical;
    leave_log_D1_input.row(i) = canonical_log;
  }
  for (arma::uword i = 0; i < n2; ++i) {
    leave_theta2_input.row(i) = fzw_restore_location(
      group2.leave_theta.row(i).t(), prepared
    ).t();
    arma::rowvec canonical(p, arma::fill::zeros);
    arma::rowvec canonical_log(p, arma::fill::zeros);
    fzw_store_scale(
      group2.leave_diagonal.row(i).t(),
      prepared,
      canonical,
      canonical_log
    );
    leave_D2_input.row(i) = canonical;
    leave_log_D2_input.row(i) = canonical_log;
  }

  arma::rowvec full_D1_input(p, arma::fill::zeros);
  arma::rowvec full_D2_input(p, arma::fill::zeros);
  arma::rowvec full_log_D1_input(p, arma::fill::zeros);
  arma::rowvec full_log_D2_input(p, arma::fill::zeros);
  fzw_store_scale(
    group1.full.diagonal,
    prepared,
    full_D1_input,
    full_log_D1_input
  );
  fzw_store_scale(
    group2.full.diagonal,
    prepared,
    full_D2_input,
    full_log_D2_input
  );

  return Rcpp::List::create(
    Rcpp::Named("z") = z,
    Rcpp::Named("R_n") = R_n,
    Rcpp::Named("sigma2_hat") = sigma2_hat,
    Rcpp::Named("sigma_hat") = sigma_hat,
    Rcpp::Named("trace_A1_squared_hat") = trace_A1,
    Rcpp::Named("trace_A2_squared_hat") = trace_A2,
    Rcpp::Named("trace_A3tA3_hat") = trace_A3,
    Rcpp::Named("variance_term1") = variance_term1,
    Rcpp::Named("variance_term2") = variance_term2,
    Rcpp::Named("variance_term3") = variance_term3,
    Rcpp::Named("c1_hat") = c1_hat,
    Rcpp::Named("c2_hat") = c2_hat,
    Rcpp::Named("numerator_inner_product") = numerator_inner,
    Rcpp::Named("numerator_contribution") = numerator_contribution,
    Rcpp::Named("numerator_zero_signs") = numerator_zero_signs,
    Rcpp::Named("numerator_zero_sign_uses") =
      static_cast<double>(numerator_zero_uses),
    Rcpp::Named("trace1_pair_squared") = trace1_pair,
    Rcpp::Named("trace2_pair_squared") = trace2_pair,
    Rcpp::Named("trace3_pair_squared") = trace3_pair,
    Rcpp::Named("trace1_pair_squared_sum") =
      static_cast<double>(trace1_sum),
    Rcpp::Named("trace2_pair_squared_sum") =
      static_cast<double>(trace2_sum),
    Rcpp::Named("trace3_pair_squared_sum") =
      static_cast<double>(trace3_sum),
    Rcpp::Named("bridge1_diagonal") = bridge1,
    Rcpp::Named("bridge2_diagonal") = bridge2,
    Rcpp::Named("own_direction1") = group1.own_direction,
    Rcpp::Named("own_direction2") = group2.own_direction,
    Rcpp::Named("own_radius1") = group1.own_radius,
    Rcpp::Named("own_radius2") = group2.own_radius,
    Rcpp::Named("own_inverse_radius1") = group1.own_inverse_radius,
    Rcpp::Named("own_inverse_radius2") = group2.own_inverse_radius,
    Rcpp::Named("full_theta1_standardised") = group1.full.theta,
    Rcpp::Named("full_theta2_standardised") = group2.full.theta,
    Rcpp::Named("full_theta1") = fzw_restore_location(
      group1.full.theta, prepared
    ),
    Rcpp::Named("full_theta2") = fzw_restore_location(
      group2.full.theta, prepared
    ),
    Rcpp::Named("full_D1_standardised") = group1.full.diagonal,
    Rcpp::Named("full_D2_standardised") = group2.full.diagonal,
    Rcpp::Named("full_D1_input_canonical") = full_D1_input.t(),
    Rcpp::Named("full_D2_input_canonical") = full_D2_input.t(),
    Rcpp::Named("full_log_D1_input_canonical") = full_log_D1_input.t(),
    Rcpp::Named("full_log_D2_input_canonical") = full_log_D2_input.t(),
    Rcpp::Named("full_iterations1") = group1.full.iterations,
    Rcpp::Named("full_iterations2") = group2.full.iterations,
    Rcpp::Named("full_location_residual1") =
      group1.full.location_residual,
    Rcpp::Named("full_location_residual2") =
      group2.full.location_residual,
    Rcpp::Named("full_scale_residual1") = group1.full.scale_residual,
    Rcpp::Named("full_scale_residual2") = group2.full.scale_residual,
    Rcpp::Named("full_minimum_training_radius1") =
      group1.full.minimum_training_radius,
    Rcpp::Named("full_minimum_training_radius2") =
      group2.full.minimum_training_radius,
    Rcpp::Named("leaveout_theta1_standardised") = group1.leave_theta,
    Rcpp::Named("leaveout_theta2_standardised") = group2.leave_theta,
    Rcpp::Named("leaveout_theta1") = leave_theta1_input,
    Rcpp::Named("leaveout_theta2") = leave_theta2_input,
    Rcpp::Named("leaveout_D1_standardised") = group1.leave_diagonal,
    Rcpp::Named("leaveout_D2_standardised") = group2.leave_diagonal,
    Rcpp::Named("leaveout_D1_input_canonical") = leave_D1_input,
    Rcpp::Named("leaveout_D2_input_canonical") = leave_D2_input,
    Rcpp::Named("leaveout_log_D1_input_canonical") = leave_log_D1_input,
    Rcpp::Named("leaveout_log_D2_input_canonical") = leave_log_D2_input,
    Rcpp::Named("leaveout_iterations1") = group1.iterations,
    Rcpp::Named("leaveout_iterations2") = group2.iterations,
    Rcpp::Named("leaveout_location_residual1") =
      group1.location_residual,
    Rcpp::Named("leaveout_location_residual2") =
      group2.location_residual,
    Rcpp::Named("leaveout_scale_residual1") = group1.scale_residual,
    Rcpp::Named("leaveout_scale_residual2") = group2.scale_residual,
    Rcpp::Named("leaveout_minimum_training_radius1") =
      group1.minimum_training_radius,
    Rcpp::Named("leaveout_minimum_training_radius2") =
      group2.minimum_training_radius,
    Rcpp::Named("sample_mean1") = prepared.mean_x,
    Rcpp::Named("sample_mean2") = prepared.mean_y,
    Rcpp::Named("column_anchor") = prepared.anchor,
    Rcpp::Named("column_log_scale") = prepared.log_scale,
    Rcpp::Named("subtraction_overflow_columns") =
      static_cast<double>(prepared.subtraction_overflow_columns),
    Rcpp::Named("n1") = static_cast<double>(n1),
    Rcpp::Named("n2") = static_cast<double>(n2),
    Rcpp::Named("p") = static_cast<double>(p),
    Rcpp::Named("ordered_pair_count1") = static_cast<double>(ordered1),
    Rcpp::Named("ordered_pair_count2") = static_cast<double>(ordered2),
    Rcpp::Named("cross_pair_count") = static_cast<double>(cross_count),
    Rcpp::Named("tolerance") = tolerance,
    Rcpp::Named("max_iterations") = max_iterations
  );
}
