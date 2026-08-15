// Yan--Zhao--Feng inverse-norm weighted max and max-sum tests.
//
// The full-sample weighted diagonal HR recursion implements the estimating
// equations in arXiv:2501.14168 with the dimensionally coherent D^{1/2}
// premultipliers.  The sum component uses unweighted leave-two-out diagonal
// HR scales and the direct feasible variance from the official Feng--Liu--Ma
// supplement.  No ridge, radius perturbation, weight cap, absolute-value
// variance repair, or numerical floor is used.

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp17)]]

#include <RcppArmadillo.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <string>
#include <vector>

namespace {

void yzf_neumaier_add(const long double value,
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

double yzf_checked_double(const long double value, const char* quantity) {
  const long double maximum = static_cast<long double>(
    std::numeric_limits<double>::max()
  );
  if (!std::isfinite(value) || value > maximum || value < -maximum) {
    Rcpp::stop(
      "Yan--Zhao--Feng computation produced a non-finite %s; no numerical "
      "floor, ridge, absolute-value repair, weight cap, or perturbation is "
      "applied.", quantity
    );
  }
  return static_cast<double>(value);
}

long double yzf_dot(const arma::rowvec& x, const arma::rowvec& y) {
  long double total = 0.0L;
  long double correction = 0.0L;
  for (arma::uword j = 0; j < x.n_elem; ++j) {
    yzf_neumaier_add(
      static_cast<long double>(x(j)) * static_cast<long double>(y(j)),
      total, correction
    );
  }
  return total + correction;
}

arma::vec yzf_stable_column_mean(const arma::mat& x) {
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
      yzf_neumaier_add(
        static_cast<long double>(x(i, j) / scale), total, correction
      );
    }
    answer(j) = yzf_checked_double(
      (total + correction) * static_cast<long double>(scale) /
        static_cast<long double>(x.n_rows),
      "sample mean"
    );
  }
  return answer;
}

struct YzfPreparedData {
  arma::mat residual;
  arma::vec reference;
  arma::vec direct_scale;
  arma::vec operand_scale;
  arma::vec reference_normalized;
  arma::vec normalized_scale;
  arma::vec log_residual_scale;
  arma::uvec overflow_fallback;
  arma::vec sample_mean;
};

YzfPreparedData yzf_prepare_data(const arma::mat& x,
                                 const arma::vec& reference) {
  const arma::uword n = x.n_rows;
  const arma::uword p = x.n_cols;
  const long double double_max = static_cast<long double>(
    std::numeric_limits<double>::max()
  );
  YzfPreparedData prepared;
  prepared.residual.set_size(n, p);
  prepared.reference = reference;
  prepared.direct_scale.zeros(p);
  prepared.operand_scale.ones(p);
  prepared.reference_normalized.zeros(p);
  prepared.normalized_scale.ones(p);
  prepared.log_residual_scale.zeros(p);
  prepared.overflow_fallback.zeros(p);
  prepared.sample_mean = yzf_stable_column_mean(x);

  for (arma::uword j = 0; j < p; ++j) {
    std::vector<long double> difference(n);
    long double maximum = 0.0L;
    bool direct_ok = true;
    for (arma::uword i = 0; i < n; ++i) {
      difference[i] = static_cast<long double>(x(i, j)) -
        static_cast<long double>(reference(j));
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
      double operand_scale = std::abs(reference(j));
      for (arma::uword i = 0; i < n; ++i) {
        operand_scale = std::max(operand_scale, std::abs(x(i, j)));
      }
      if (!(operand_scale > 0.0) || !std::isfinite(operand_scale)) {
        Rcpp::stop(
          "Yan--Zhao--Feng computation could not construct a finite "
          "residual in variable %llu.",
          static_cast<unsigned long long>(j + 1)
        );
      }
      prepared.operand_scale(j) = operand_scale;
      prepared.reference_normalized(j) = reference(j) / operand_scale;
      double normalized_maximum = 0.0;
      for (arma::uword i = 0; i < n; ++i) {
        const double value = x(i, j) / operand_scale -
          prepared.reference_normalized(j);
        prepared.residual(i, j) = value;
        normalized_maximum = std::max(normalized_maximum, std::abs(value));
      }
      if (!(normalized_maximum > 0.0) ||
          !std::isfinite(normalized_maximum)) {
        Rcpp::stop(
          "Yan--Zhao--Feng computation requires variable %llu to differ "
          "from its centring reference in at least one observation.",
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
        "Yan--Zhao--Feng computation requires positive variation in every "
        "variable; variable %llu is constant. No ridge is applied.",
        static_cast<unsigned long long>(j + 1)
      );
    }
  }
  if (!prepared.residual.is_finite()) {
    Rcpp::stop("Yan--Zhao--Feng internal residual standardisation failed.");
  }
  return prepared;
}

arma::vec yzf_restore_location(const arma::rowvec& theta,
                               const YzfPreparedData& prepared,
                               const char* quantity) {
  arma::vec answer(theta.n_elem);
  for (arma::uword j = 0; j < theta.n_elem; ++j) {
    long double value;
    if (prepared.overflow_fallback(j) == 0u) {
      value = static_cast<long double>(prepared.reference(j)) +
        static_cast<long double>(prepared.direct_scale(j)) *
          static_cast<long double>(theta(j));
    } else {
      value = static_cast<long double>(prepared.operand_scale(j)) *
        (static_cast<long double>(prepared.reference_normalized(j)) +
         static_cast<long double>(prepared.normalized_scale(j)) *
           static_cast<long double>(theta(j)));
    }
    answer(j) = yzf_checked_double(value, quantity);
  }
  return answer;
}

bool yzf_direction_radius(const arma::rowvec& observation,
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
        "A Yan--Zhao--Feng diagonal iterate is not finite and strictly "
        "positive; no ridge is applied."
      );
    }
    direction(j) = (observation(j) - location(j)) /
      std::sqrt(diagonal(j));
    if (!std::isfinite(direction(j))) {
      Rcpp::stop(
        "Yan--Zhao--Feng diagonal standardisation produced a non-finite "
        "value."
      );
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
    yzf_neumaier_add(value * value, squares, square_correction);
  }
  const long double normalized_norm = std::sqrt(
    squares + square_correction
  );
  if (!(normalized_norm > 0.0L) || !std::isfinite(normalized_norm)) {
    Rcpp::stop("Yan--Zhao--Feng could not normalise a nonzero residual.");
  }
  radius = yzf_checked_double(
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

bool yzf_is_omitted(const arma::uword row, const int omit_i,
                    const int omit_j) {
  return static_cast<int>(row) == omit_i ||
    static_cast<int>(row) == omit_j;
}

arma::uword yzf_subset_size(const arma::mat& data, const int omit_i,
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

void yzf_subset_initial(const arma::mat& data, const int omit_i,
                        const int omit_j, const std::string& label,
                        arma::rowvec& location, arma::vec& diagonal) {
  const arma::uword retained = yzf_subset_size(data, omit_i, omit_j);
  if (retained < 2) {
    Rcpp::stop("%s requires at least two retained observations.", label.c_str());
  }
  location.zeros(data.n_cols);
  for (arma::uword j = 0; j < data.n_cols; ++j) {
    long double total = 0.0L;
    long double correction = 0.0L;
    for (arma::uword i = 0; i < data.n_rows; ++i) {
      if (!yzf_is_omitted(i, omit_i, omit_j)) {
        yzf_neumaier_add(
          static_cast<long double>(data(i, j)), total, correction
        );
      }
    }
    location(j) = static_cast<double>(
      (total + correction) / static_cast<long double>(retained)
    );
  }

  diagonal.zeros(data.n_cols);
  for (arma::uword j = 0; j < data.n_cols; ++j) {
    long double total = 0.0L;
    long double correction = 0.0L;
    for (arma::uword i = 0; i < data.n_rows; ++i) {
      if (!yzf_is_omitted(i, omit_i, omit_j)) {
        const long double residual =
          static_cast<long double>(data(i, j)) -
          static_cast<long double>(location(j));
        yzf_neumaier_add(residual * residual, total, correction);
      }
    }
    diagonal(j) = static_cast<double>(
      (total + correction) / static_cast<long double>(retained - 1)
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

struct YzfScore {
  arma::rowvec direction_sum;
  arma::vec direction_square_sum;
  arma::mat directions;
  arma::vec radii;
  double minimum_radius;
  int zero_count;
  double location_residual;
  double diagonal_residual;
  double score_residual;
};

YzfScore yzf_subset_score(const arma::mat& data, const int omit_i,
                          const int omit_j, const arma::rowvec& location,
                          const arma::vec& diagonal, const double power_m,
                          const double zero_tol,
                          const std::string& label) {
  const arma::uword retained = yzf_subset_size(data, omit_i, omit_j);
  const double p = static_cast<double>(data.n_cols);
  arma::rowvec direction_sum(data.n_cols, arma::fill::zeros);
  arma::vec direction_square_sum(data.n_cols, arma::fill::zeros);
  arma::mat directions(retained, data.n_cols, arma::fill::zeros);
  arma::vec radii(retained, arma::fill::zeros);
  arma::uword position = 0;
  int zero_count = 0;
  for (arma::uword row = 0; row < data.n_rows; ++row) {
    if (yzf_is_omitted(row, omit_i, omit_j)) {
      continue;
    }
    arma::rowvec direction;
    double radius = 0.0;
    const bool nonzero = yzf_direction_radius(
      data.row(row), location, diagonal, zero_tol, direction, radius
    );
    if (!nonzero) {
      ++zero_count;
      if (power_m < 1.0) {
        Rcpp::stop(
          "%s is undefined because a retained observation has a "
          "standardised residual at or below `zero_tol` and r^(m-1) is "
          "singular for m < 1; no perturbation is applied.", label.c_str()
        );
      }
    }
    directions.row(position) = direction;
    radii(position) = radius;
    direction_sum += direction;
    direction_square_sum += arma::square(direction).t();
    ++position;
  }

  double maximum_log_weight = -std::numeric_limits<double>::infinity();
  for (arma::uword i = 0; i < retained; ++i) {
    double log_weight;
    if (radii(i) == 0.0) {
      log_weight = power_m == 0.0 ? 0.0 :
        -std::numeric_limits<double>::infinity();
    } else {
      log_weight = power_m * std::log(radii(i));
    }
    maximum_log_weight = std::max(maximum_log_weight, log_weight);
  }
  if (!std::isfinite(maximum_log_weight)) {
    Rcpp::stop(
      "%s has no positive finite weighted-sign mass; no repair is applied.",
      label.c_str()
    );
  }
  long double weight_total = 0.0L;
  long double weight_correction = 0.0L;
  std::vector<long double> location_total(data.n_cols, 0.0L);
  std::vector<long double> location_correction(data.n_cols, 0.0L);
  for (arma::uword i = 0; i < retained; ++i) {
    long double relative_weight = 0.0L;
    if (radii(i) == 0.0) {
      relative_weight = power_m == 0.0 ?
        std::exp(-static_cast<long double>(maximum_log_weight)) : 0.0L;
    } else {
      relative_weight = std::exp(
        static_cast<long double>(power_m) * std::log(radii(i)) -
        static_cast<long double>(maximum_log_weight)
      );
    }
    yzf_neumaier_add(relative_weight, weight_total, weight_correction);
    for (arma::uword j = 0; j < data.n_cols; ++j) {
      yzf_neumaier_add(
        relative_weight * static_cast<long double>(directions(i, j)),
        location_total[j], location_correction[j]
      );
    }
  }
  const long double weight_sum = weight_total + weight_correction;
  if (!(weight_sum > 0.0L) || !std::isfinite(weight_sum)) {
    Rcpp::stop("%s has an invalid weighted-score denominator.", label.c_str());
  }
  long double location_norm_square = 0.0L;
  long double location_norm_correction = 0.0L;
  for (arma::uword j = 0; j < data.n_cols; ++j) {
    const long double component =
      (location_total[j] + location_correction[j]) / weight_sum;
    yzf_neumaier_add(
      component * component,
      location_norm_square, location_norm_correction
    );
  }
  const double location_residual = std::sqrt(static_cast<double>(
    location_norm_square + location_norm_correction
  ));
  const arma::vec diagonal_equation = p * direction_square_sum /
    static_cast<double>(retained);
  const double diagonal_residual = arma::abs(
    diagonal_equation - 1.0
  ).max();
  return YzfScore{
    direction_sum,
    direction_square_sum,
    directions,
    radii,
    radii.min(),
    zero_count,
    location_residual,
    diagonal_residual,
    std::max(location_residual, diagonal_residual)
  };
}

struct YzfFit {
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
  int zero_residual_count = 0;
};

YzfFit yzf_joint_fit(const arma::mat& data, const int omit_i,
                     const int omit_j, const double power_m,
                     const double tol, const int max_iter,
                     const double zero_tol, const std::string& label) {
  YzfFit fit;
  yzf_subset_initial(
    data, omit_i, omit_j, label, fit.location, fit.diagonal
  );
  const arma::vec initial_diagonal = fit.diagonal;
  const double p = static_cast<double>(data.n_cols);
  const double retained = static_cast<double>(
    yzf_subset_size(data, omit_i, omit_j)
  );

  for (int update = 0; update < max_iter; ++update) {
    const YzfScore score = yzf_subset_score(
      data, omit_i, omit_j, fit.location, fit.diagonal,
      power_m, zero_tol, label
    );
    fit.minimum_residual_distance = std::min(
      fit.minimum_residual_distance, score.minimum_radius
    );
    fit.zero_residual_count = std::max(
      fit.zero_residual_count, score.zero_count
    );

    arma::rowvec step(data.n_cols, arma::fill::zeros);
    if (power_m == 0.0) {
      // Preserve the Feng--Liu--Ma/INST operation order exactly for every
      // leave-two-out diagonal HR fit.
      const arma::vec relative_inverse_radius =
        score.minimum_radius / score.radii;
      const double denominator = arma::accu(relative_inverse_radius);
      if (!(denominator > 0.0) || !std::isfinite(denominator)) {
        Rcpp::stop(
          "%s has an invalid inverse-radius update denominator; no repair "
          "is applied.", label.c_str()
        );
      }
      step = score.direction_sum *
        (score.minimum_radius / denominator);
    } else if (power_m == 1.0) {
      std::vector<long double> total(data.n_cols, 0.0L);
      std::vector<long double> correction(data.n_cols, 0.0L);
      for (arma::uword i = 0; i < score.radii.n_elem; ++i) {
        for (arma::uword j = 0; j < data.n_cols; ++j) {
          yzf_neumaier_add(
            static_cast<long double>(score.radii(i)) *
              static_cast<long double>(score.directions(i, j)),
            total[j], correction[j]
          );
        }
      }
      for (arma::uword j = 0; j < data.n_cols; ++j) {
        step(j) = yzf_checked_double(
          (total[j] + correction[j]) / retained,
          "m = 1 location update"
        );
      }
    } else {
      const double minimum_radius = score.radii.min();
      if (!(minimum_radius > zero_tol)) {
        Rcpp::stop(
          "%s has a singular r^(m-1) location update; no perturbation is "
          "applied.", label.c_str()
        );
      }
      const long double anchor =
        static_cast<long double>(power_m - 1.0) *
        std::log(minimum_radius);
      long double denominator = 0.0L;
      long double denominator_correction = 0.0L;
      std::vector<long double> total(data.n_cols, 0.0L);
      std::vector<long double> correction(data.n_cols, 0.0L);
      for (arma::uword i = 0; i < score.radii.n_elem; ++i) {
        const long double relative = std::exp(
          static_cast<long double>(power_m - 1.0) *
            std::log(score.radii(i)) - anchor
        );
        yzf_neumaier_add(
          relative, denominator, denominator_correction
        );
        for (arma::uword j = 0; j < data.n_cols; ++j) {
          yzf_neumaier_add(
            relative * static_cast<long double>(score.radii(i)) *
              static_cast<long double>(score.directions(i, j)),
            total[j], correction[j]
          );
        }
      }
      const long double denominator_value =
        denominator + denominator_correction;
      if (!(denominator_value > 0.0L) ||
          !std::isfinite(denominator_value)) {
        Rcpp::stop(
          "%s has an invalid r^(m-1) update denominator; no repair is "
          "applied.", label.c_str()
        );
      }
      for (arma::uword j = 0; j < data.n_cols; ++j) {
        step(j) = yzf_checked_double(
          (total[j] + correction[j]) / denominator_value,
          "weighted location update"
        );
      }
    }

    const arma::rowvec next_location = fit.location +
      arma::sqrt(fit.diagonal).t() % step;
    const arma::vec next_diagonal = p * fit.diagonal %
      (score.direction_square_sum / retained);
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

  const YzfScore final_score = yzf_subset_score(
    data, omit_i, omit_j, fit.location, fit.diagonal,
    power_m, zero_tol, label
  );
  fit.minimum_residual_distance = std::min(
    fit.minimum_residual_distance, final_score.minimum_radius
  );
  fit.zero_residual_count = std::max(
    fit.zero_residual_count, final_score.zero_count
  );
  fit.location_score_residual = final_score.location_residual;
  fit.diagonal_score_residual = final_score.diagonal_residual;
  fit.score_residual = final_score.score_residual;
  return fit;
}

Rcpp::List yzf_fit_diagnostics(const YzfFit& fit) {
  return Rcpp::List::create(
    Rcpp::Named("iterations") = fit.iterations,
    Rcpp::Named("iteration.stable") = fit.iteration_stable,
    Rcpp::Named("relative.update") = fit.relative_update,
    Rcpp::Named("location.relative.update") = fit.location_update,
    Rcpp::Named("log.diagonal.relative.update") = fit.diagonal_update,
    Rcpp::Named("score.residual") = fit.score_residual,
    Rcpp::Named("location.score.residual") = fit.location_score_residual,
    Rcpp::Named("diagonal.score.residual") = fit.diagonal_score_residual,
    Rcpp::Named("minimum.residual.distance") =
      fit.minimum_residual_distance,
    Rcpp::Named("zero.residual.count") = fit.zero_residual_count,
    Rcpp::Named("convergence.basis") =
      "relative location and log-diagonal iterate update"
  );
}

Rcpp::List yzf_pair_diagnostics(const std::vector<YzfFit>& fits) {
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
  Rcpp::IntegerVector zero_count(size);
  bool all_stable = true;
  int failures = 0;
  int maximum_iterations = 0;
  double worst_relative = 0.0;
  double worst_score = 0.0;
  double minimum_residual = R_PosInf;
  int total_zero = 0;
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
    zero_count[i] = fits[i].zero_residual_count;
    all_stable = all_stable && fits[i].iteration_stable;
    failures += fits[i].iteration_stable ? 0 : 1;
    maximum_iterations = std::max(maximum_iterations, fits[i].iterations);
    worst_relative = std::max(worst_relative, fits[i].relative_update);
    worst_score = std::max(worst_score, fits[i].score_residual);
    minimum_residual = std::min(
      minimum_residual, fits[i].minimum_residual_distance
    );
    total_zero += fits[i].zero_residual_count;
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
    Rcpp::Named("zero.residual.count") = zero_count,
    Rcpp::Named("all.iteration.stable") = all_stable,
    Rcpp::Named("stability.failures") = failures,
    Rcpp::Named("maximum.iterations") = maximum_iterations,
    Rcpp::Named("worst.relative.update") = worst_relative,
    Rcpp::Named("worst.score.residual") = worst_score,
    Rcpp::Named("smallest.residual.distance") = minimum_residual,
    Rcpp::Named("total.zero.residual.count") = total_zero,
    Rcpp::Named("convergence.basis") =
      "relative location and log-diagonal iterate update"
  );
}

struct YzfScaleDisplay {
  arma::vec log_diagonal;
  arma::vec diagonal;
};

YzfScaleDisplay yzf_scale_display(const arma::vec& diagonal,
                                  const YzfPreparedData& prepared) {
  arma::vec log_diagonal(diagonal.n_elem);
  double maximum = -std::numeric_limits<double>::infinity();
  for (arma::uword j = 0; j < diagonal.n_elem; ++j) {
    log_diagonal(j) = 2.0 * prepared.log_residual_scale(j) +
      std::log(diagonal(j));
    if (!std::isfinite(log_diagonal(j))) {
      Rcpp::stop(
        "Yan--Zhao--Feng input-coordinate log diagonal is non-finite."
      );
    }
    maximum = std::max(maximum, log_diagonal(j));
  }
  log_diagonal -= maximum;
  return YzfScaleDisplay{log_diagonal, arma::exp(log_diagonal)};
}

double yzf_radial_weight(const double radius, const double exponent,
                         const std::string& label) {
  if (radius == 0.0) {
    if (exponent < 0.0) {
      Rcpp::stop(
        "%s is undefined because a zero radius is raised to a negative "
        "power; no perturbation is applied.", label.c_str()
      );
    }
    return exponent == 0.0 ? 1.0 : 0.0;
  }
  if (exponent == -1.0) {
    const double answer = 1.0 / radius;
    if (!std::isfinite(answer)) {
      Rcpp::stop("%s has a non-finite inverse radius.", label.c_str());
    }
    return answer;
  }
  if (exponent == 0.0) {
    return 1.0;
  }
  if (exponent == 1.0) {
    return radius;
  }
  const double answer = std::exp(exponent * std::log(radius));
  if (!std::isfinite(answer)) {
    Rcpp::stop(
      "%s produced a non-finite radial weight; no weight cap is applied.",
      label.c_str()
    );
  }
  return answer;
}

struct YzfMoment {
  double value;
  double log_value;
};

YzfMoment yzf_mean_power(const arma::vec& radii, const double exponent,
                         const std::string& label) {
  if (exponent == 0.0) {
    return YzfMoment{1.0, 0.0};
  }
  double maximum_log = -std::numeric_limits<double>::infinity();
  std::vector<double> log_terms(radii.n_elem);
  for (arma::uword i = 0; i < radii.n_elem; ++i) {
    if (radii(i) == 0.0) {
      if (exponent < 0.0) {
        Rcpp::stop(
          "%s is undefined because a zero radius enters a negative radial "
          "moment; no perturbation is applied.", label.c_str()
        );
      }
      log_terms[i] = -std::numeric_limits<double>::infinity();
    } else {
      log_terms[i] = exponent * std::log(radii(i));
      maximum_log = std::max(maximum_log, log_terms[i]);
    }
  }
  if (!std::isfinite(maximum_log)) {
    Rcpp::stop("%s is zero for every observation.", label.c_str());
  }
  long double total = 0.0L;
  long double correction = 0.0L;
  for (arma::uword i = 0; i < radii.n_elem; ++i) {
    if (std::isfinite(log_terms[i])) {
      yzf_neumaier_add(
        std::exp(static_cast<long double>(log_terms[i] - maximum_log)),
        total, correction
      );
    }
  }
  const long double log_mean =
    static_cast<long double>(maximum_log) +
    std::log((total + correction) /
             static_cast<long double>(radii.n_elem));
  if (!std::isfinite(log_mean)) {
    Rcpp::stop("%s has a non-finite log moment.", label.c_str());
  }
  const long double log_maximum = std::log(
    static_cast<long double>(std::numeric_limits<double>::max())
  );
  const double value = log_mean > log_maximum ?
    R_PosInf : static_cast<double>(std::exp(log_mean));
  return YzfMoment{value, static_cast<double>(log_mean)};
}

struct YzfMaxCore {
  double raw;
  double centred;
  double maximum_standardized_location;
  double zeta_m_minus_1;
  double zeta_2m;
  double log_zeta_m_minus_1;
  double log_zeta_2m;
  double log_moment_ratio;
  double moment_ratio;
  double finite_sample_correction;
  arma::vec radii;
  arma::vec weights;
  arma::mat directions;
};

YzfMaxCore yzf_max_core(const arma::mat& data, const YzfFit& fit,
                        const double power_m, const double zero_tol,
                        const std::string& label) {
  const arma::uword n = data.n_rows;
  const arma::uword p = data.n_cols;
  const YzfScore score = yzf_subset_score(
    data, -1, -1, fit.location, fit.diagonal,
    power_m, zero_tol, label
  );
  arma::vec weights(n);
  for (arma::uword i = 0; i < n; ++i) {
    weights(i) = yzf_radial_weight(
      score.radii(i), power_m, label
    );
  }
  const YzfMoment first = yzf_mean_power(
    score.radii, power_m - 1.0, "zeta_(m-1) estimate"
  );
  const YzfMoment second = yzf_mean_power(
    score.radii, 2.0 * power_m, "zeta_(2m) estimate"
  );
  const double log_ratio = 2.0 * first.log_value - second.log_value;
  if (!std::isfinite(log_ratio)) {
    Rcpp::stop("The weighted max radial-moment ratio is non-finite.");
  }
  const long double log_double_max = std::log(
    static_cast<long double>(std::numeric_limits<double>::max())
  );
  const double moment_ratio =
    static_cast<long double>(log_ratio) > log_double_max ?
    R_PosInf : static_cast<double>(std::exp(
      static_cast<long double>(log_ratio)
    ));
  const arma::rowvec standardized_location = fit.location /
    arma::sqrt(fit.diagonal).t();
  const double maximum = arma::abs(standardized_location).max();
  const double correction = 1.0 - 1.0 / std::sqrt(
    static_cast<double>(n)
  );
  long double raw = 0.0L;
  if (maximum > 0.0) {
    const long double log_raw =
      std::log(static_cast<long double>(n)) +
      2.0L * std::log(static_cast<long double>(maximum)) +
      static_cast<long double>(log_ratio) +
      std::log(static_cast<long double>(p)) +
      std::log(static_cast<long double>(correction));
    if (!std::isfinite(log_raw) || log_raw > log_double_max) {
      Rcpp::stop(
        "The weighted max statistic is non-finite; no rescaling cap is "
        "applied."
      );
    }
    raw = std::exp(log_raw);
  }
  const double raw_double = yzf_checked_double(raw, "weighted max statistic");
  const double centred = raw_double - 2.0 * std::log(
    static_cast<double>(p)
  ) + std::log(std::log(static_cast<double>(p)));
  if (!std::isfinite(centred)) {
    Rcpp::stop("The centred weighted max statistic is non-finite.");
  }
  return YzfMaxCore{
    raw_double,
    centred,
    maximum,
    first.value,
    second.value,
    first.log_value,
    second.log_value,
    log_ratio,
    moment_ratio,
    correction,
    score.radii,
    weights,
    score.directions
  };
}

void yzf_validate_inputs(const arma::mat& x, const double power_m,
                         const double tol, const int max_iter,
                         const double zero_tol) {
  if (!x.is_finite()) {
    Rcpp::stop("`x` must contain only finite values.");
  }
  if (!std::isfinite(power_m) || power_m > 1.0) {
    Rcpp::stop("`m` must be finite and no greater than one.");
  }
  if (!std::isfinite(tol) || tol <= 0.0 || max_iter < 1 ||
      !std::isfinite(zero_tol) || zero_tol < 0.0) {
    Rcpp::stop("Invalid Yan--Zhao--Feng iteration controls.");
  }
}

Rcpp::List yzf_full_result(const YzfPreparedData& prepared,
                           const YzfFit& fit,
                           const YzfMaxCore* maximum) {
  const YzfScaleDisplay display = yzf_scale_display(
    fit.diagonal, prepared
  );
  Rcpp::List answer = Rcpp::List::create(
    Rcpp::Named("location_standardized") = fit.location,
    Rcpp::Named("location") = yzf_restore_location(
      fit.location, prepared, "weighted location estimate"
    ),
    Rcpp::Named("diagonal_standardized") = fit.diagonal,
    Rcpp::Named("log_diagonal_input_canonical") = display.log_diagonal,
    Rcpp::Named("diagonal_input_canonical") = display.diagonal,
    Rcpp::Named("fit_diagnostics") = yzf_fit_diagnostics(fit),
    Rcpp::Named("iteration_stable") = fit.iteration_stable,
    Rcpp::Named("sample_mean") = prepared.sample_mean,
    Rcpp::Named("column_log_residual_scale") =
      prepared.log_residual_scale,
    Rcpp::Named("subtraction_overflow_columns") = static_cast<double>(
      arma::accu(prepared.overflow_fallback)
    )
  );
  if (maximum != nullptr) {
    answer["max_raw"] = maximum->raw;
    answer["max_centred"] = maximum->centred;
    answer["maximum_standardized_location"] =
      maximum->maximum_standardized_location;
    answer["zeta_m_minus_1"] = maximum->zeta_m_minus_1;
    answer["zeta_2m"] = maximum->zeta_2m;
    answer["log_zeta_m_minus_1"] = maximum->log_zeta_m_minus_1;
    answer["log_zeta_2m"] = maximum->log_zeta_2m;
    answer["log_moment_ratio"] = maximum->log_moment_ratio;
    answer["moment_ratio"] = maximum->moment_ratio;
    answer["finite_sample_correction"] =
      maximum->finite_sample_correction;
    answer["full_radius"] = maximum->radii;
    answer["full_weight"] = maximum->weights;
    answer["full_direction"] = maximum->directions;
  }
  return answer;
}

Rcpp::List yzf_sum_core(const arma::mat& data,
                        const YzfPreparedData& prepared,
                        const double power_m, const double tol,
                        const int max_iter, const double zero_tol) {
  const arma::uword n = data.n_rows;
  const arma::uword p = data.n_cols;
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
  arma::mat endpoint_weight(pair_count, 2);
  arma::vec pair_test_inner(pair_count);
  arma::vec pair_statistic_kernel(pair_count);
  arma::vec pair_sign_mean_norm(pair_count);
  arma::mat pair_variance_factor(pair_count, 2);
  arma::vec pair_variance_kernel(pair_count);
  arma::vec pair_trace_inner_squared(pair_count);
  arma::vec pair_trace_zero_signs(pair_count, arma::fill::zeros);
  arma::vec pair_null_zero_signs(pair_count, arma::fill::zeros);
  std::vector<YzfFit> fits;
  fits.reserve(pair_count);

  long double statistic_sum = 0.0L;
  long double statistic_correction = 0.0L;
  long double variance_ordered_sum = 0.0L;
  long double variance_ordered_correction = 0.0L;
  long double trace_square_sum = 0.0L;
  long double trace_square_correction = 0.0L;
  long double radial_moment_sum = 0.0L;
  long double radial_moment_correction = 0.0L;

  arma::uword pair_index = 0;
  for (arma::uword i = 0; i + 1 < n; ++i) {
    for (arma::uword j = i + 1; j < n; ++j, ++pair_index) {
      const std::string label = "weighted sum leave-two-out fit (" +
        std::to_string(i + 1) + ", " + std::to_string(j + 1) + ")";
      YzfFit fit = yzf_joint_fit(
        data, static_cast<int>(i), static_cast<int>(j), 0.0,
        tol, max_iter, zero_tol, label
      );
      fits.push_back(fit);
      pair_i(pair_index) = static_cast<double>(i + 1);
      pair_j(pair_index) = static_cast<double>(j + 1);
      leaveout_location_standardized.row(pair_index) = fit.location;
      leaveout_location.row(pair_index) = yzf_restore_location(
        fit.location, prepared, "leave-two-out location"
      ).t();
      leaveout_diagonal_standardized.row(pair_index) = fit.diagonal.t();
      const YzfScaleDisplay display = yzf_scale_display(
        fit.diagonal, prepared
      );
      leaveout_log_diagonal_canonical.row(pair_index) =
        display.log_diagonal.t();
      leaveout_diagonal_canonical.row(pair_index) = display.diagonal.t();

      arma::rowvec zero(p, arma::fill::zeros);
      arma::rowvec direction_i;
      arma::rowvec direction_j;
      double radius_i = 0.0;
      double radius_j = 0.0;
      const bool nonzero_i = yzf_direction_radius(
        data.row(i), zero, fit.diagonal, zero_tol,
        direction_i, radius_i
      );
      const bool nonzero_j = yzf_direction_radius(
        data.row(j), zero, fit.diagonal, zero_tol,
        direction_j, radius_j
      );
      if (!nonzero_i) {
        pair_null_zero_signs(pair_index) += 1.0;
      }
      if (!nonzero_j) {
        pair_null_zero_signs(pair_index) += 1.0;
      }
      if ((!nonzero_i || !nonzero_j) && power_m < 0.0) {
        Rcpp::stop(
          "%s is undefined because an endpoint null-centred radius is at "
          "or below `zero_tol` and is raised to a negative power; no "
          "perturbation is applied.", label.c_str()
        );
      }
      const double weight_i = nonzero_i ? yzf_radial_weight(
        radius_i, power_m, label + " endpoint radial weight"
      ) : (power_m == 0.0 ? 1.0 : 0.0);
      const double weight_j = nonzero_j ? yzf_radial_weight(
        radius_j, power_m, label + " endpoint radial weight"
      ) : (power_m == 0.0 ? 1.0 : 0.0);
      endpoint_radius(pair_index, 0) = radius_i;
      endpoint_radius(pair_index, 1) = radius_j;
      endpoint_weight(pair_index, 0) = weight_i;
      endpoint_weight(pair_index, 1) = weight_j;
      yzf_neumaier_add(
        static_cast<long double>(weight_i) * weight_i,
        radial_moment_sum, radial_moment_correction
      );
      yzf_neumaier_add(
        static_cast<long double>(weight_j) * weight_j,
        radial_moment_sum, radial_moment_correction
      );

      const long double inner = yzf_dot(direction_i, direction_j);
      const long double pair_weight =
        static_cast<long double>(weight_i) *
        static_cast<long double>(weight_j);
      const long double statistic_kernel = pair_weight * inner;
      pair_test_inner(pair_index) = yzf_checked_double(
        inner, "pair test inner product"
      );
      pair_statistic_kernel(pair_index) = yzf_checked_double(
        statistic_kernel, "pair statistic kernel"
      );
      yzf_neumaier_add(
        statistic_kernel, statistic_sum, statistic_correction
      );

      std::vector<long double> mean_total(p, 0.0L);
      std::vector<long double> mean_correction(p, 0.0L);
      for (arma::uword k = 0; k < n; ++k) {
        if (k == i || k == j) {
          continue;
        }
        arma::rowvec direction_k;
        double radius_k = 0.0;
        const bool nonzero_k = yzf_direction_radius(
          data.row(k), zero, fit.diagonal, zero_tol,
          direction_k, radius_k
        );
        if (!nonzero_k) {
          pair_null_zero_signs(pair_index) += 1.0;
        }
        for (arma::uword ell = 0; ell < p; ++ell) {
          yzf_neumaier_add(
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
        yzf_neumaier_add(
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
        yzf_dot(leaveout_sign_mean, direction_j);
      const long double second = inner -
        yzf_dot(leaveout_sign_mean, direction_i);
      const long double variance_kernel = pair_weight * pair_weight *
        first * second;
      pair_variance_factor(pair_index, 0) = yzf_checked_double(
        first, "first feasible variance factor"
      );
      pair_variance_factor(pair_index, 1) = yzf_checked_double(
        second, "second feasible variance factor"
      );
      pair_variance_kernel(pair_index) = yzf_checked_double(
        variance_kernel, "feasible variance kernel"
      );
      yzf_neumaier_add(
        2.0L * variance_kernel,
        variance_ordered_sum, variance_ordered_correction
      );

      arma::rowvec trace_direction_i;
      arma::rowvec trace_direction_j;
      double trace_radius_i = 0.0;
      double trace_radius_j = 0.0;
      const bool trace_nonzero_i = yzf_direction_radius(
        data.row(i), fit.location, fit.diagonal, 0.0,
        trace_direction_i, trace_radius_i
      );
      const bool trace_nonzero_j = yzf_direction_radius(
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
      const long double trace_inner = yzf_dot(
        trace_direction_i, trace_direction_j
      );
      const long double trace_square = trace_inner * trace_inner;
      pair_trace_inner_squared(pair_index) = yzf_checked_double(
        trace_square, "pair trace inner product squared"
      );
      yzf_neumaier_add(
        trace_square, trace_square_sum, trace_square_correction
      );
    }
  }

  const long double raw_statistic_ld =
    (statistic_sum + statistic_correction) /
    static_cast<long double>(pair_count);
  const long double ordered_variance_kernel_ld =
    variance_ordered_sum + variance_ordered_correction;
  const long double variance_ld = 2.0L * ordered_variance_kernel_ld /
    n_fourth;
  const double raw_statistic = yzf_checked_double(
    raw_statistic_ld, "weighted sum statistic"
  );
  const double variance = yzf_checked_double(
    variance_ld, "direct feasible variance"
  );
  if (!(variance > 0.0)) {
    Rcpp::stop(
      "The weighted sum test requires the primary direct feasible variance "
      "to be strictly positive. It was %.17g; no absolute value or variance "
      "floor is applied.", variance
    );
  }
  const double standard_error = std::sqrt(variance);
  const double z = raw_statistic / standard_error;
  if (!std::isfinite(z)) {
    Rcpp::stop("The weighted sum standardised statistic is non-finite.");
  }

  const long double zeta_2m_ld =
    (radial_moment_sum + radial_moment_correction) / ordered_count;
  const long double trace_r2_ld =
    static_cast<long double>(p) * static_cast<long double>(p) *
    (trace_square_sum + trace_square_correction) /
    static_cast<long double>(pair_count);
  const long double oracle_pair_ld = 2.0L * zeta_2m_ld * zeta_2m_ld *
    trace_r2_ld /
    (ordered_count * static_cast<long double>(p) *
     static_cast<long double>(p));
  const long double oracle_asymptotic_ld =
    2.0L * zeta_2m_ld * zeta_2m_ld * trace_r2_ld /
    (n_ld * n_ld * static_cast<long double>(p) *
     static_cast<long double>(p));
  const long double weighted_trace_ld =
    variance_ld * ordered_count / 2.0L;
  const Rcpp::List fit_diagnostics = yzf_pair_diagnostics(fits);

  return Rcpp::List::create(
    Rcpp::Named("z") = z,
    Rcpp::Named("T") = raw_statistic,
    Rcpp::Named("variance") = variance,
    Rcpp::Named("standard_error") = standard_error,
    Rcpp::Named("ordered_variance_kernel_sum") = yzf_checked_double(
      ordered_variance_kernel_ld, "ordered feasible variance kernel sum"
    ),
    Rcpp::Named("zeta_2m_crossfit") = yzf_checked_double(
      zeta_2m_ld, "cross-fitted zeta_(2m) estimate"
    ),
    Rcpp::Named("trace_R2_hat") = yzf_checked_double(
      trace_r2_ld, "cross-fitted trace R squared"
    ),
    Rcpp::Named("oracle_variance_pair_adjusted") = yzf_checked_double(
      oracle_pair_ld, "pair-adjusted oracle-factorized variance"
    ),
    Rcpp::Named("oracle_variance_asymptotic") = yzf_checked_double(
      oracle_asymptotic_ld, "asymptotic oracle-factorized variance"
    ),
    Rcpp::Named("weighted_score_trace_hat") = yzf_checked_double(
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
    Rcpp::Named("pair_null_zero_signs") = pair_null_zero_signs,
    Rcpp::Named("endpoint_radius") = endpoint_radius,
    Rcpp::Named("endpoint_weight") = endpoint_weight,
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
    Rcpp::Named("pair_count") = static_cast<double>(pair_count),
    Rcpp::Named("ordered_pair_count") = static_cast<double>(n * (n - 1))
  );
}

}  // namespace


//' Weighted scaled spatial median and diagonal scale kernel
//'
//' @param x Numeric observation-by-variable matrix.
//' @param m Finite radial power no greater than one.
//' @param tol Positive relative-update tolerance.
//' @param max_iter Positive maximum number of updates.
//' @param zero_tol Non-negative singular-radius tolerance.
//' @return Internal list containing the fitted location, diagonal scale,
//'   radial quantities, and iteration diagnostics.
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List cpp_weighted_scaled_spatial_median(
    const arma::mat& x, const double m, const double tol,
    const int max_iter, const double zero_tol) {
  yzf_validate_inputs(x, m, tol, max_iter, zero_tol);
  if (x.n_rows < 2 || x.n_cols < 1) {
    Rcpp::stop(
      "Weighted scaled spatial median fitting requires at least two rows "
      "and one variable."
    );
  }
  const arma::vec reference = yzf_stable_column_mean(x);
  const YzfPreparedData prepared = yzf_prepare_data(x, reference);
  const YzfFit fit = yzf_joint_fit(
    prepared.residual, -1, -1, m, tol, max_iter, zero_tol,
    "weighted scaled spatial median fit"
  );
  const YzfScore score = yzf_subset_score(
    prepared.residual, -1, -1, fit.location, fit.diagonal,
    m, zero_tol, "weighted scaled spatial median fit"
  );
  arma::vec weights(x.n_rows);
  for (arma::uword i = 0; i < x.n_rows; ++i) {
    weights(i) = yzf_radial_weight(
      score.radii(i), m, "weighted scaled spatial median radial weight"
    );
  }
  Rcpp::List answer = yzf_full_result(prepared, fit, nullptr);
  answer["full_radius"] = score.radii;
  answer["full_weight"] = weights;
  answer["full_direction"] = score.directions;
  answer["m"] = m;
  answer["n"] = static_cast<double>(x.n_rows);
  answer["p"] = static_cast<double>(x.n_cols);
  answer["tol"] = tol;
  answer["max_iter"] = max_iter;
  answer["zero_tol"] = zero_tol;
  return answer;
}


//' Yan--Zhao--Feng weighted max statistic kernel
//'
//' @param x Numeric observation-by-variable matrix.
//' @param mu Numeric null-location vector.
//' @param m Finite radial power no greater than one.
//' @param tol Positive relative-update tolerance.
//' @param max_iter Positive maximum number of updates.
//' @param zero_tol Non-negative singular-radius tolerance.
//' @return Internal list containing the weighted max statistic and fit.
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List cpp_yzf_weighted_max(
    const arma::mat& x, const arma::vec& mu, const double m,
    const double tol, const int max_iter, const double zero_tol) {
  yzf_validate_inputs(x, m, tol, max_iter, zero_tol);
  if (!mu.is_finite() || mu.n_elem != x.n_cols) {
    Rcpp::stop("`mu` must contain one finite value for every column of `x`.");
  }
  if (x.n_rows < 2 || x.n_cols < 2) {
    Rcpp::stop(
      "The weighted max test requires at least two observations and two "
      "variables."
    );
  }
  const YzfPreparedData prepared = yzf_prepare_data(x, mu);
  const YzfFit fit = yzf_joint_fit(
    prepared.residual, -1, -1, m, tol, max_iter, zero_tol,
    "weighted max full-sample fit"
  );
  const YzfMaxCore maximum = yzf_max_core(
    prepared.residual, fit, m, zero_tol,
    "weighted max full-sample fit"
  );
  Rcpp::List answer = yzf_full_result(prepared, fit, &maximum);
  answer["m"] = m;
  answer["n"] = static_cast<double>(x.n_rows);
  answer["p"] = static_cast<double>(x.n_cols);
  answer["tol"] = tol;
  answer["max_iter"] = max_iter;
  answer["zero_tol"] = zero_tol;
  return answer;
}


//' Yan--Zhao--Feng weighted max-sum statistic kernel
//'
//' @param x Numeric observation-by-variable matrix.
//' @param mu Numeric null-location vector.
//' @param m Finite radial power no greater than one.
//' @param tol Positive relative-update tolerance.
//' @param max_iter Positive maximum number of updates.
//' @param zero_tol Non-negative singular-radius tolerance.
//' @return Internal list containing weighted max and feasible weighted sum
//'   components and all fit diagnostics.
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List cpp_yzf_weighted_maxsum(
    const arma::mat& x, const arma::vec& mu, const double m,
    const double tol, const int max_iter, const double zero_tol) {
  yzf_validate_inputs(x, m, tol, max_iter, zero_tol);
  if (!mu.is_finite() || mu.n_elem != x.n_cols) {
    Rcpp::stop("`mu` must contain one finite value for every column of `x`.");
  }
  if (x.n_rows < 4 || x.n_cols < 2) {
    Rcpp::stop(
      "The weighted max-sum test requires at least four observations and "
      "two variables."
    );
  }
  const YzfPreparedData prepared = yzf_prepare_data(x, mu);
  const YzfFit full_fit = yzf_joint_fit(
    prepared.residual, -1, -1, m, tol, max_iter, zero_tol,
    "weighted max full-sample fit"
  );
  const YzfMaxCore maximum = yzf_max_core(
    prepared.residual, full_fit, m, zero_tol,
    "weighted max full-sample fit"
  );
  Rcpp::List answer = yzf_full_result(prepared, full_fit, &maximum);
  answer["sum"] = yzf_sum_core(
    prepared.residual, prepared, m, tol, max_iter, zero_tol
  );
  answer["m"] = m;
  answer["n"] = static_cast<double>(x.n_rows);
  answer["p"] = static_cast<double>(x.n_cols);
  answer["tol"] = tol;
  answer["max_iter"] = max_iter;
  answer["zero_tol"] = zero_tol;
  return answer;
}
