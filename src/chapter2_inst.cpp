// Feng--Liu--Ma one-sample inverse norm sign test.
//
// The primary statistic uses pair-specific leave-two-out diagonal HR fits.
// The fitted location is used only while estimating the diagonal scale: every
// endpoint and the centring vector in the feasible variance are formed from
// X - mu0.  The primary variance is the literal 2 n^{-4} ordered-pair
// estimator in Section S.3 of the paper's official supplement.

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp17)]]

#include <RcppArmadillo.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <string>
#include <vector>

namespace {

void inst_neumaier_add(const long double value,
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

double inst_checked_double(const long double value, const char* quantity) {
  const long double maximum = static_cast<long double>(
    std::numeric_limits<double>::max()
  );
  if (!std::isfinite(value) || value > maximum || value < -maximum) {
    Rcpp::stop(
      "INST produced a non-finite %s; no numerical floor, ridge, "
      "absolute-value repair, or perturbation is applied.", quantity
    );
  }
  return static_cast<double>(value);
}

struct InstPreparedData {
  arma::mat residual;
  arma::vec null_location;
  arma::vec direct_scale;
  arma::vec operand_scale;
  arma::vec null_normalized;
  arma::vec normalized_scale;
  arma::vec log_residual_scale;
  arma::uvec overflow_fallback;
  arma::vec sample_mean;
};

arma::vec inst_stable_column_mean(const arma::mat& x) {
  arma::vec answer(x.n_cols, arma::fill::zeros);
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
      inst_neumaier_add(
        static_cast<long double>(x(i, j) / scale), total, correction
      );
    }
    answer(j) = inst_checked_double(
      (total + correction) * static_cast<long double>(scale) /
        static_cast<long double>(x.n_rows),
      "sample mean"
    );
  }
  return answer;
}

InstPreparedData inst_prepare_data(const arma::mat& x,
                                   const arma::vec& mu) {
  const arma::uword n = x.n_rows;
  const arma::uword p = x.n_cols;
  const long double double_max = static_cast<long double>(
    std::numeric_limits<double>::max()
  );
  InstPreparedData prepared;
  prepared.residual.set_size(n, p);
  prepared.null_location = mu;
  prepared.direct_scale.zeros(p);
  prepared.operand_scale.ones(p);
  prepared.null_normalized.zeros(p);
  prepared.normalized_scale.ones(p);
  prepared.log_residual_scale.zeros(p);
  prepared.overflow_fallback.zeros(p);
  prepared.sample_mean = inst_stable_column_mean(x);

  for (arma::uword j = 0; j < p; ++j) {
    std::vector<long double> difference(n);
    long double maximum = 0.0L;
    bool direct_ok = true;
    for (arma::uword i = 0; i < n; ++i) {
      difference[i] = static_cast<long double>(x(i, j)) -
        static_cast<long double>(mu(j));
      direct_ok = direct_ok && std::isfinite(difference[i]);
      maximum = std::max(maximum, std::abs(difference[i]));
    }
    direct_ok = direct_ok && std::isfinite(maximum) &&
      maximum > 0.0L && maximum <= double_max;

    if (direct_ok) {
      const double scale = static_cast<double>(maximum);
      prepared.direct_scale(j) = scale;
      prepared.log_residual_scale(j) = std::log(scale);
      for (arma::uword i = 0; i < n; ++i) {
        prepared.residual(i, j) = static_cast<double>(
          difference[i] / maximum
        );
      }
    } else {
      prepared.overflow_fallback(j) = 1u;
      double operand_scale = std::abs(mu(j));
      for (arma::uword i = 0; i < n; ++i) {
        operand_scale = std::max(operand_scale, std::abs(x(i, j)));
      }
      if (!(operand_scale > 0.0) || !std::isfinite(operand_scale)) {
        Rcpp::stop(
          "INST could not construct a finite residual in variable %llu.",
          static_cast<unsigned long long>(j + 1)
        );
      }
      prepared.operand_scale(j) = operand_scale;
      prepared.null_normalized(j) = mu(j) / operand_scale;
      double normalized_maximum = 0.0;
      for (arma::uword i = 0; i < n; ++i) {
        const double value = x(i, j) / operand_scale -
          prepared.null_normalized(j);
        prepared.residual(i, j) = value;
        normalized_maximum = std::max(normalized_maximum, std::abs(value));
      }
      if (!(normalized_maximum > 0.0) ||
          !std::isfinite(normalized_maximum)) {
        Rcpp::stop(
          "INST requires variable %llu to differ from its null location "
          "in at least one observation.",
          static_cast<unsigned long long>(j + 1)
        );
      }
      prepared.normalized_scale(j) = normalized_maximum;
      prepared.log_residual_scale(j) = std::log(operand_scale) +
        std::log(normalized_maximum);
      prepared.residual.col(j) /= normalized_maximum;
    }

    double minimum = prepared.residual(0, j);
    double column_maximum = prepared.residual(0, j);
    for (arma::uword i = 1; i < n; ++i) {
      minimum = std::min(minimum, prepared.residual(i, j));
      column_maximum = std::max(column_maximum, prepared.residual(i, j));
    }
    if (!(column_maximum > minimum)) {
      Rcpp::stop(
        "INST requires positive variation in every variable; variable "
        "%llu is constant. No ridge is applied.",
        static_cast<unsigned long long>(j + 1)
      );
    }
  }
  if (!prepared.residual.is_finite()) {
    Rcpp::stop("INST internal null residual standardisation failed.");
  }
  return prepared;
}

arma::vec inst_restore_location(const arma::rowvec& theta,
                                const InstPreparedData& prepared) {
  arma::vec answer(theta.n_elem);
  for (arma::uword j = 0; j < theta.n_elem; ++j) {
    long double value;
    if (prepared.overflow_fallback(j) == 0u) {
      value = static_cast<long double>(prepared.null_location(j)) +
        static_cast<long double>(prepared.direct_scale(j)) *
          static_cast<long double>(theta(j));
    } else {
      value = static_cast<long double>(prepared.operand_scale(j)) *
        (static_cast<long double>(prepared.null_normalized(j)) +
         static_cast<long double>(prepared.normalized_scale(j)) *
           static_cast<long double>(theta(j)));
    }
    answer(j) = inst_checked_double(value, "leave-two-out location");
  }
  return answer;
}

bool inst_direction_radius(const arma::rowvec& observation,
                           const arma::rowvec& location,
                           const arma::vec& diagonal,
                           const double zero_tol,
                           arma::rowvec& direction,
                           double& radius) {
  const arma::uword p = observation.n_elem;
  direction.set_size(p);
  double maximum = 0.0;
  for (arma::uword j = 0; j < p; ++j) {
    if (!std::isfinite(diagonal(j)) || diagonal(j) <= 0.0) {
      Rcpp::stop(
        "An INST diagonal iterate is not finite and strictly positive; "
        "no ridge is applied."
      );
    }
    direction(j) = (observation(j) - location(j)) /
      std::sqrt(diagonal(j));
    if (!std::isfinite(direction(j))) {
      Rcpp::stop("INST diagonal standardisation produced a non-finite value.");
    }
    maximum = std::max(maximum, std::abs(direction(j)));
  }
  if (maximum == 0.0) {
    direction.zeros();
    radius = 0.0;
    return false;
  }

  long double squares = 0.0L;
  long double square_correction = 0.0L;
  for (arma::uword j = 0; j < p; ++j) {
    direction(j) /= maximum;
    const long double value = static_cast<long double>(direction(j));
    inst_neumaier_add(value * value, squares, square_correction);
  }
  const long double normalized_norm = std::sqrt(
    squares + square_correction
  );
  if (!(normalized_norm > 0.0L) || !std::isfinite(normalized_norm)) {
    Rcpp::stop("INST could not normalise a nonzero residual.");
  }
  radius = inst_checked_double(
    static_cast<long double>(maximum) * normalized_norm,
    "standardised radius"
  );
  if (!(radius > zero_tol)) {
    direction.zeros();
    return false;
  }
  direction /= static_cast<double>(normalized_norm);
  return true;
}

long double inst_dot(const arma::rowvec& x, const arma::rowvec& y) {
  long double total = 0.0L;
  long double correction = 0.0L;
  for (arma::uword j = 0; j < x.n_elem; ++j) {
    inst_neumaier_add(
      static_cast<long double>(x(j)) * static_cast<long double>(y(j)),
      total, correction
    );
  }
  return total + correction;
}

struct InstFit {
  arma::rowvec location;
  arma::vec diagonal;
  int iterations = 0;
  bool iteration_stable = false;
  double relative_update = R_PosInf;
  double location_update = R_PosInf;
  double diagonal_update = R_PosInf;
  double score_residual = R_PosInf;
  double location_score_residual = R_PosInf;
  double diagonal_score_residual = R_PosInf;
  double minimum_residual_distance = R_PosInf;
};

bool inst_is_omitted(const arma::uword row, const int omit_i,
                     const int omit_j) {
  return static_cast<int>(row) == omit_i ||
    static_cast<int>(row) == omit_j;
}

arma::uword inst_subset_size(const arma::mat& data, const int omit_i,
                             const int omit_j) {
  arma::uword omitted = 0;
  if (omit_i >= 0) {
    ++omitted;
  }
  if (omit_j >= 0 && omit_j != omit_i) {
    ++omitted;
  }
  return data.n_rows - omitted;
}

void inst_subset_initial(const arma::mat& data, const int omit_i,
                         const int omit_j, const std::string& label,
                         arma::rowvec& location, arma::vec& diagonal) {
  const arma::uword m = inst_subset_size(data, omit_i, omit_j);
  if (m < 2) {
    Rcpp::stop("%s requires at least two retained observations.", label.c_str());
  }
  location.zeros(data.n_cols);
  for (arma::uword j = 0; j < data.n_cols; ++j) {
    long double total = 0.0L;
    long double correction = 0.0L;
    for (arma::uword i = 0; i < data.n_rows; ++i) {
      if (!inst_is_omitted(i, omit_i, omit_j)) {
        inst_neumaier_add(
          static_cast<long double>(data(i, j)), total, correction
        );
      }
    }
    location(j) = static_cast<double>(
      (total + correction) / static_cast<long double>(m)
    );
  }

  diagonal.zeros(data.n_cols);
  for (arma::uword j = 0; j < data.n_cols; ++j) {
    long double total = 0.0L;
    long double correction = 0.0L;
    for (arma::uword i = 0; i < data.n_rows; ++i) {
      if (!inst_is_omitted(i, omit_i, omit_j)) {
        const long double residual =
          static_cast<long double>(data(i, j)) -
          static_cast<long double>(location(j));
        inst_neumaier_add(residual * residual, total, correction);
      }
    }
    diagonal(j) = static_cast<double>(
      (total + correction) / static_cast<long double>(m - 1)
    );
  }
  if (!location.is_finite() || !diagonal.is_finite() ||
      arma::any(diagonal <= 0.0)) {
    Rcpp::stop(
      "%s requires every retained marginal sample variance to be finite "
      "and strictly positive; no ridge is applied.", label.c_str()
    );
  }
}

struct InstScore {
  arma::rowvec direction_sum;
  arma::vec direction_square_sum;
  arma::vec radii;
  double minimum_radius;
  double location_residual;
  double diagonal_residual;
  double score_residual;
};

InstScore inst_subset_score(const arma::mat& data, const int omit_i,
                            const int omit_j,
                            const arma::rowvec& location,
                            const arma::vec& diagonal,
                            const double zero_tol,
                            const std::string& label) {
  const arma::uword m = inst_subset_size(data, omit_i, omit_j);
  const double p = static_cast<double>(data.n_cols);
  arma::rowvec direction_sum(data.n_cols, arma::fill::zeros);
  arma::vec direction_square_sum(data.n_cols, arma::fill::zeros);
  arma::vec radii(m);
  arma::uword position = 0;
  for (arma::uword row = 0; row < data.n_rows; ++row) {
    if (inst_is_omitted(row, omit_i, omit_j)) {
      continue;
    }
    arma::rowvec direction;
    double radius = 0.0;
    if (!inst_direction_radius(
          data.row(row), location, diagonal, zero_tol, direction, radius
        )) {
      Rcpp::stop(
        "%s is undefined because a retained observation has a "
        "standardised residual at or below `zero_tol`; no perturbation is "
        "applied.", label.c_str()
      );
    }
    direction_sum += direction;
    direction_square_sum += arma::square(direction).t();
    radii(position++) = radius;
  }
  const arma::rowvec mean_direction = direction_sum /
    static_cast<double>(m);
  const arma::vec diagonal_equation = p * direction_square_sum /
    static_cast<double>(m);
  const double location_residual = arma::norm(mean_direction, 2);
  const double diagonal_residual = arma::abs(
    diagonal_equation - 1.0
  ).max();
  return InstScore{
    direction_sum,
    direction_square_sum,
    radii,
    radii.min(),
    location_residual,
    diagonal_residual,
    std::max(location_residual, diagonal_residual)
  };
}

InstFit inst_joint_fit(const arma::mat& data, const int omit_i,
                       const int omit_j, const double tol,
                       const int max_iter, const double zero_tol,
                       const std::string& label) {
  InstFit fit;
  inst_subset_initial(
    data, omit_i, omit_j, label, fit.location, fit.diagonal
  );
  const arma::vec initial_diagonal = fit.diagonal;
  const double p = static_cast<double>(data.n_cols);
  const double m = static_cast<double>(
    inst_subset_size(data, omit_i, omit_j)
  );

  for (int update = 0; update < max_iter; ++update) {
    const InstScore score = inst_subset_score(
      data, omit_i, omit_j, fit.location, fit.diagonal,
      zero_tol, label
    );
    fit.minimum_residual_distance = std::min(
      fit.minimum_residual_distance, score.minimum_radius
    );
    const arma::vec relative_inverse_radius =
      score.minimum_radius / score.radii;
    const double denominator = arma::accu(relative_inverse_radius);
    if (!(denominator > 0.0) || !std::isfinite(denominator)) {
      Rcpp::stop(
        "%s has an invalid inverse-radius update denominator; no repair "
        "is applied.", label.c_str()
      );
    }
    const arma::rowvec step = score.direction_sum *
      (score.minimum_radius / denominator);
    const arma::rowvec next_location = fit.location +
      arma::sqrt(fit.diagonal).t() % step;
    const arma::vec next_diagonal = p * fit.diagonal %
      (score.direction_square_sum / m);
    if (!next_location.is_finite() || !next_diagonal.is_finite() ||
        arma::any(next_diagonal <= 0.0)) {
      Rcpp::stop(
        "%s produced a non-finite or non-positive iterate; no ridge or "
        "floor is applied.", label.c_str()
      );
    }
    fit.location_update = arma::abs(
      (next_location - fit.location) / arma::sqrt(initial_diagonal).t()
    ).max();
    fit.diagonal_update = arma::abs(
      arma::log(next_diagonal / fit.diagonal)
    ).max();
    fit.relative_update = std::max(
      fit.location_update, fit.diagonal_update
    );
    fit.location = next_location;
    fit.diagonal = next_diagonal;
    fit.iterations = update + 1;
    if (fit.relative_update <= tol) {
      fit.iteration_stable = true;
      break;
    }
  }

  const InstScore final_score = inst_subset_score(
    data, omit_i, omit_j, fit.location, fit.diagonal,
    zero_tol, label
  );
  fit.minimum_residual_distance = std::min(
    fit.minimum_residual_distance, final_score.minimum_radius
  );
  fit.location_score_residual = final_score.location_residual;
  fit.diagonal_score_residual = final_score.diagonal_residual;
  fit.score_residual = final_score.score_residual;
  return fit;
}

arma::rowvec inst_null_direction(const arma::mat& data,
                                 const arma::uword row,
                                 const arma::vec& diagonal,
                                 const double zero_tol,
                                 double& radius,
                                 const std::string& label) {
  arma::rowvec zero(data.n_cols, arma::fill::zeros);
  arma::rowvec direction;
  if (!inst_direction_radius(
        data.row(row), zero, diagonal, zero_tol, direction, radius
      )) {
    Rcpp::stop(
      "%s is undefined because observation %llu has a null-centred "
      "standardised radius at or below `zero_tol`; inverse-norm weighting "
      "is singular and no perturbation is applied.",
      label.c_str(), static_cast<unsigned long long>(row + 1)
    );
  }
  return direction;
}

Rcpp::List inst_pair_diagnostics(const std::vector<InstFit>& fits) {
  const R_xlen_t size = static_cast<R_xlen_t>(fits.size());
  Rcpp::IntegerVector iterations(size);
  Rcpp::LogicalVector stable(size);
  Rcpp::NumericVector relative_update(size);
  Rcpp::NumericVector location_update(size);
  Rcpp::NumericVector diagonal_update(size);
  Rcpp::NumericVector score_residual(size);
  Rcpp::NumericVector location_score(size);
  Rcpp::NumericVector diagonal_score(size);
  Rcpp::NumericVector minimum_distance(size);
  bool all_stable = true;
  int failures = 0;
  int maximum_iterations = 0;
  double worst_relative = 0.0;
  double worst_score = 0.0;
  double minimum_residual = R_PosInf;
  for (R_xlen_t i = 0; i < size; ++i) {
    iterations[i] = fits[i].iterations;
    stable[i] = fits[i].iteration_stable;
    relative_update[i] = fits[i].relative_update;
    location_update[i] = fits[i].location_update;
    diagonal_update[i] = fits[i].diagonal_update;
    score_residual[i] = fits[i].score_residual;
    location_score[i] = fits[i].location_score_residual;
    diagonal_score[i] = fits[i].diagonal_score_residual;
    minimum_distance[i] = fits[i].minimum_residual_distance;
    all_stable = all_stable && fits[i].iteration_stable;
    failures += fits[i].iteration_stable ? 0 : 1;
    maximum_iterations = std::max(maximum_iterations, fits[i].iterations);
    worst_relative = std::max(worst_relative, fits[i].relative_update);
    worst_score = std::max(worst_score, fits[i].score_residual);
    minimum_residual = std::min(
      minimum_residual, fits[i].minimum_residual_distance
    );
  }
  return Rcpp::List::create(
    Rcpp::Named("iterations") = iterations,
    Rcpp::Named("iteration.stable") = stable,
    Rcpp::Named("relative.update") = relative_update,
    Rcpp::Named("location.relative.update") = location_update,
    Rcpp::Named("log.diagonal.relative.update") = diagonal_update,
    Rcpp::Named("score.residual") = score_residual,
    Rcpp::Named("location.score.residual") = location_score,
    Rcpp::Named("diagonal.score.residual") = diagonal_score,
    Rcpp::Named("minimum.residual.distance") = minimum_distance,
    Rcpp::Named("all.iteration.stable") = all_stable,
    Rcpp::Named("stability.failures") = failures,
    Rcpp::Named("maximum.iterations") = maximum_iterations,
    Rcpp::Named("worst.relative.update") = worst_relative,
    Rcpp::Named("worst.score.residual") = worst_score,
    Rcpp::Named("smallest.residual.distance") = minimum_residual,
    Rcpp::Named("convergence.basis") =
      "relative location and log-diagonal iterate update"
  );
}

}  // namespace


//' Feng--Liu--Ma one-sample inverse norm sign statistic kernel
//'
//' @param x Numeric observation-by-variable matrix.
//' @param mu Numeric null-location vector.
//' @param tol Positive relative-update tolerance.
//' @param max_iter Positive maximum number of updates per leave-two-out fit.
//' @param zero_tol Non-negative singular-radius tolerance.
//' @return Internal list of primary statistic, direct feasible variance,
//'   cross-fitted nuisance diagnostics, and iteration diagnostics.
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List cpp_inst_one_sample(const arma::mat& x,
                               const arma::vec& mu,
                               const double tol,
                               const int max_iter,
                               const double zero_tol) {
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
      "INST requires at least four observations so every leave-two-out "
      "diagonal fit has at least two retained rows."
    );
  }
  if (p < 1 || mu.n_elem != p) {
    Rcpp::stop("`mu` must have one entry for every column of `x`.");
  }
  if (!std::isfinite(tol) || tol <= 0.0 || max_iter < 1 ||
      !std::isfinite(zero_tol) || zero_tol < 0.0) {
    Rcpp::stop("Invalid INST iteration controls.");
  }

  const InstPreparedData prepared = inst_prepare_data(x, mu);
  const arma::mat& data = prepared.residual;
  const arma::uword pair_count = n * (n - 1) / 2;
  const long double n_ld = static_cast<long double>(n);
  const long double ordered_count = n_ld * static_cast<long double>(n - 1);
  const long double n_fourth = n_ld * n_ld * n_ld * n_ld;

  arma::vec pair_i(pair_count);
  arma::vec pair_j(pair_count);
  arma::mat leaveout_location_standardized(pair_count, p);
  arma::mat leaveout_location(pair_count, p);
  arma::mat leaveout_diagonal_standardized(pair_count, p);
  arma::mat leaveout_log_diagonal_canonical(pair_count, p);
  arma::mat leaveout_diagonal_canonical(pair_count, p);
  arma::mat endpoint_radius(pair_count, 2);
  arma::mat endpoint_inverse_radius(pair_count, 2);
  arma::vec pair_test_inner(pair_count);
  arma::vec pair_statistic_kernel(pair_count);
  arma::vec pair_sign_mean_norm(pair_count);
  arma::mat pair_variance_factor(pair_count, 2);
  arma::vec pair_variance_kernel(pair_count);
  arma::vec pair_trace_inner_squared(pair_count);
  arma::vec pair_trace_zero_signs(pair_count, arma::fill::zeros);
  std::vector<InstFit> fits;
  fits.reserve(pair_count);

  long double statistic_sum = 0.0L;
  long double statistic_correction = 0.0L;
  long double variance_ordered_sum = 0.0L;
  long double variance_ordered_correction = 0.0L;
  long double trace_square_sum = 0.0L;
  long double trace_square_correction = 0.0L;
  long double inverse_radius_sum = 0.0L;
  long double inverse_radius_correction = 0.0L;
  long double inverse_radius_square_sum = 0.0L;
  long double inverse_radius_square_correction = 0.0L;

  arma::uword pair_index = 0;
  for (arma::uword i = 0; i + 1 < n; ++i) {
    for (arma::uword j = i + 1; j < n; ++j, ++pair_index) {
      const std::string label = "INST leave-two-out fit (" +
        std::to_string(i + 1) + ", " + std::to_string(j + 1) + ")";
      InstFit fit = inst_joint_fit(
        data, static_cast<int>(i), static_cast<int>(j),
        tol, max_iter, zero_tol, label
      );
      fits.push_back(fit);
      pair_i(pair_index) = static_cast<double>(i + 1);
      pair_j(pair_index) = static_cast<double>(j + 1);
      leaveout_location_standardized.row(pair_index) = fit.location;
      leaveout_location.row(pair_index) =
        inst_restore_location(fit.location, prepared).t();
      leaveout_diagonal_standardized.row(pair_index) = fit.diagonal.t();

      double maximum_log_diagonal = -std::numeric_limits<double>::infinity();
      for (arma::uword ell = 0; ell < p; ++ell) {
        const double value = 2.0 * prepared.log_residual_scale(ell) +
          std::log(fit.diagonal(ell));
        leaveout_log_diagonal_canonical(pair_index, ell) = value;
        maximum_log_diagonal = std::max(maximum_log_diagonal, value);
      }
      for (arma::uword ell = 0; ell < p; ++ell) {
        const double canonical_log =
          leaveout_log_diagonal_canonical(pair_index, ell) -
          maximum_log_diagonal;
        leaveout_log_diagonal_canonical(pair_index, ell) = canonical_log;
        leaveout_diagonal_canonical(pair_index, ell) =
          std::exp(canonical_log);
      }

      double radius_i = 0.0;
      double radius_j = 0.0;
      const arma::rowvec direction_i = inst_null_direction(
        data, i, fit.diagonal, zero_tol, radius_i, label
      );
      const arma::rowvec direction_j = inst_null_direction(
        data, j, fit.diagonal, zero_tol, radius_j, label
      );
      const double inverse_i = 1.0 / radius_i;
      const double inverse_j = 1.0 / radius_j;
      if (!std::isfinite(inverse_i) || !std::isfinite(inverse_j)) {
        Rcpp::stop(
          "%s has a non-finite inverse radius; no weight cap is applied.",
          label.c_str()
        );
      }
      endpoint_radius(pair_index, 0) = radius_i;
      endpoint_radius(pair_index, 1) = radius_j;
      endpoint_inverse_radius(pair_index, 0) = inverse_i;
      endpoint_inverse_radius(pair_index, 1) = inverse_j;
      inst_neumaier_add(
        static_cast<long double>(inverse_i), inverse_radius_sum,
        inverse_radius_correction
      );
      inst_neumaier_add(
        static_cast<long double>(inverse_j), inverse_radius_sum,
        inverse_radius_correction
      );
      inst_neumaier_add(
        static_cast<long double>(inverse_i) * inverse_i,
        inverse_radius_square_sum, inverse_radius_square_correction
      );
      inst_neumaier_add(
        static_cast<long double>(inverse_j) * inverse_j,
        inverse_radius_square_sum, inverse_radius_square_correction
      );

      const long double inner = inst_dot(direction_i, direction_j);
      const long double statistic_kernel = inner *
        static_cast<long double>(inverse_i) *
        static_cast<long double>(inverse_j);
      pair_test_inner(pair_index) = inst_checked_double(
        inner, "pair test inner product"
      );
      pair_statistic_kernel(pair_index) = inst_checked_double(
        statistic_kernel, "pair statistic kernel"
      );
      inst_neumaier_add(
        statistic_kernel, statistic_sum, statistic_correction
      );

      std::vector<long double> mean_total(p, 0.0L);
      std::vector<long double> mean_correction(p, 0.0L);
      for (arma::uword k = 0; k < n; ++k) {
        if (k == i || k == j) {
          continue;
        }
        double radius_k = 0.0;
        const arma::rowvec direction_k = inst_null_direction(
          data, k, fit.diagonal, zero_tol, radius_k, label
        );
        for (arma::uword ell = 0; ell < p; ++ell) {
          inst_neumaier_add(
            static_cast<long double>(direction_k(ell)),
            mean_total[ell], mean_correction[ell]
          );
        }
      }
      arma::rowvec leaveout_sign_mean(p);
      long double sign_mean_norm_square = 0.0L;
      long double sign_mean_norm_correction = 0.0L;
      const long double leaveout_count = static_cast<long double>(n - 2);
      for (arma::uword ell = 0; ell < p; ++ell) {
        leaveout_sign_mean(ell) = static_cast<double>(
          (mean_total[ell] + mean_correction[ell]) / leaveout_count
        );
        const long double value = leaveout_sign_mean(ell);
        inst_neumaier_add(
          value * value,
          sign_mean_norm_square, sign_mean_norm_correction
        );
      }
      pair_sign_mean_norm(pair_index) = std::sqrt(
        static_cast<double>(
          sign_mean_norm_square + sign_mean_norm_correction
        )
      );
      const long double first = inner -
        inst_dot(leaveout_sign_mean, direction_j);
      const long double second = inner -
        inst_dot(leaveout_sign_mean, direction_i);
      const long double weight_square =
        static_cast<long double>(inverse_i) * inverse_i *
        static_cast<long double>(inverse_j) * inverse_j;
      const long double variance_kernel = weight_square * first * second;
      pair_variance_factor(pair_index, 0) = inst_checked_double(
        first, "first feasible variance factor"
      );
      pair_variance_factor(pair_index, 1) = inst_checked_double(
        second, "second feasible variance factor"
      );
      pair_variance_kernel(pair_index) = inst_checked_double(
        variance_kernel, "feasible variance kernel"
      );
      inst_neumaier_add(
        2.0L * variance_kernel,
        variance_ordered_sum, variance_ordered_correction
      );

      arma::rowvec trace_direction_i;
      arma::rowvec trace_direction_j;
      double trace_radius_i = 0.0;
      double trace_radius_j = 0.0;
      const bool trace_nonzero_i = inst_direction_radius(
        data.row(i), fit.location, fit.diagonal, 0.0,
        trace_direction_i, trace_radius_i
      );
      const bool trace_nonzero_j = inst_direction_radius(
        data.row(j), fit.location, fit.diagonal, 0.0,
        trace_direction_j, trace_radius_j
      );
      if (!trace_nonzero_i) {
        trace_direction_i.zeros(p);
        pair_trace_zero_signs(pair_index) += 1.0;
      }
      if (!trace_nonzero_j) {
        trace_direction_j.zeros(p);
        pair_trace_zero_signs(pair_index) += 1.0;
      }
      const long double trace_inner = inst_dot(
        trace_direction_i, trace_direction_j
      );
      const long double trace_square = trace_inner * trace_inner;
      pair_trace_inner_squared(pair_index) = inst_checked_double(
        trace_square, "pair trace inner product squared"
      );
      inst_neumaier_add(
        trace_square, trace_square_sum, trace_square_correction
      );
    }
  }

  const long double raw_statistic_ld =
    (statistic_sum + statistic_correction) /
    static_cast<long double>(pair_count);
  const long double variance_ordered_kernel_ld =
    variance_ordered_sum + variance_ordered_correction;
  const long double variance_ld = 2.0L * variance_ordered_kernel_ld /
    n_fourth;
  const double raw_statistic = inst_checked_double(
    raw_statistic_ld, "raw statistic"
  );
  const double variance = inst_checked_double(
    variance_ld, "direct feasible variance"
  );
  if (!(variance > 0.0)) {
    Rcpp::stop(
      "INST requires the primary S.3 direct feasible variance to be "
      "strictly positive. It was %.17g; no absolute value or variance "
      "floor is applied.", variance
    );
  }
  const double standard_error = std::sqrt(variance);
  const double z = raw_statistic / standard_error;
  if (!std::isfinite(z)) {
    Rcpp::stop("INST produced a non-finite standardised statistic.");
  }

  const long double endpoint_count = ordered_count;
  const long double inverse_radius_mean_ld =
    (inverse_radius_sum + inverse_radius_correction) / endpoint_count;
  const long double nu2_ld =
    (inverse_radius_square_sum + inverse_radius_square_correction) /
    endpoint_count;
  const long double trace_r2_ld =
    static_cast<long double>(p) * static_cast<long double>(p) *
    (trace_square_sum + trace_square_correction) /
    static_cast<long double>(pair_count);
  const long double oracle_variance_ld = 2.0L * nu2_ld * nu2_ld *
    trace_r2_ld /
    (ordered_count * static_cast<long double>(p) *
     static_cast<long double>(p));
  const long double weighted_trace_ld = variance_ld * ordered_count / 2.0L;

  const Rcpp::List fit_diagnostics = inst_pair_diagnostics(fits);
  return Rcpp::List::create(
    Rcpp::Named("z") = z,
    Rcpp::Named("T") = raw_statistic,
    Rcpp::Named("variance") = variance,
    Rcpp::Named("standard_error") = standard_error,
    Rcpp::Named("ordered_variance_kernel_sum") = inst_checked_double(
      variance_ordered_kernel_ld, "ordered feasible variance kernel sum"
    ),
    Rcpp::Named("inverse_radius_mean") = inst_checked_double(
      inverse_radius_mean_ld, "cross-fitted inverse-radius mean"
    ),
    Rcpp::Named("nu2") = inst_checked_double(
      nu2_ld, "cross-fitted second inverse-radius moment"
    ),
    Rcpp::Named("c0") = inst_checked_double(
      nu2_ld, "cross-fitted inverse-norm c0"
    ),
    Rcpp::Named("trace_R2_hat") = inst_checked_double(
      trace_r2_ld, "cross-fitted trace R squared"
    ),
    Rcpp::Named("oracle_factorized_variance") = inst_checked_double(
      oracle_variance_ld, "oracle-factorized diagnostic variance"
    ),
    Rcpp::Named("weighted_score_trace_hat") = inst_checked_double(
      weighted_trace_ld, "weighted score covariance trace"
    ),
    Rcpp::Named("pair_i") = pair_i,
    Rcpp::Named("pair_j") = pair_j,
    Rcpp::Named("pair_test_inner_product") = pair_test_inner,
    Rcpp::Named("pair_statistic_kernel") = pair_statistic_kernel,
    Rcpp::Named("pair_sign_mean_norm") = pair_sign_mean_norm,
    Rcpp::Named("pair_variance_factor") = pair_variance_factor,
    Rcpp::Named("pair_variance_kernel") = pair_variance_kernel,
    Rcpp::Named("pair_trace_inner_product_squared") =
      pair_trace_inner_squared,
    Rcpp::Named("pair_trace_zero_signs") = pair_trace_zero_signs,
    Rcpp::Named("endpoint_radius") = endpoint_radius,
    Rcpp::Named("endpoint_inverse_radius") = endpoint_inverse_radius,
    Rcpp::Named("leaveout_location_standardized") =
      leaveout_location_standardized,
    Rcpp::Named("leaveout_location") = leaveout_location,
    Rcpp::Named("leaveout_diagonal_standardized") =
      leaveout_diagonal_standardized,
    Rcpp::Named("leaveout_log_diagonal_input_canonical") =
      leaveout_log_diagonal_canonical,
    Rcpp::Named("leaveout_diagonal_input_canonical") =
      leaveout_diagonal_canonical,
    Rcpp::Named("fit_diagnostics") = fit_diagnostics,
    Rcpp::Named("all_iteration_stable") =
      fit_diagnostics["all.iteration.stable"],
    Rcpp::Named("stability_failures") =
      fit_diagnostics["stability.failures"],
    Rcpp::Named("sample_mean") = prepared.sample_mean,
    Rcpp::Named("column_log_residual_scale") =
      prepared.log_residual_scale,
    Rcpp::Named("subtraction_overflow_columns") = static_cast<double>(
      arma::accu(prepared.overflow_fallback)
    ),
    Rcpp::Named("pair_count") = static_cast<double>(pair_count),
    Rcpp::Named("ordered_pair_count") = static_cast<double>(
      n * (n - 1)
    ),
    Rcpp::Named("n") = static_cast<double>(n),
    Rcpp::Named("p") = static_cast<double>(p),
    Rcpp::Named("tol") = tol,
    Rcpp::Named("max_iter") = max_iter,
    Rcpp::Named("zero_tol") = zero_tol
  );
}
