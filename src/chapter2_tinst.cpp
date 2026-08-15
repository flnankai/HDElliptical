// Huang--Liu--Zhou--Feng two-sample inverse norm sign test.
//
// This file deliberately follows the feasible statistic in Huang et al.
// (2023).  In particular, all location and diagonal estimates used by the
// statistic are observation-specific leave-one-out estimates.  The three
// trace estimators use full-sample diagonal estimates and ordered pairs.

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp17)]]

#include <RcppArmadillo.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <string>
#include <vector>

namespace {

struct TinstTransform {
  arma::mat x;
  arma::mat y;
  arma::rowvec base;
  arma::rowvec scale;
  arma::rowvec raw_scale;
  arma::rowvec midpoint_normalized;
  arma::rowvec spread_normalized;
  std::vector<bool> normalized_fallback;
};

struct TinstFit {
  arma::rowvec location;
  arma::vec diagonal;
  int iterations = 0;
  bool iteration_stable = false;
  double relative_update = R_PosInf;
  double location_change = R_PosInf;
  double diagonal_change = R_NaReal;
  double score_residual = R_PosInf;
  double location_residual = R_PosInf;
  double diagonal_residual = R_NaReal;
  double minimum_residual_distance = R_PosInf;
};

void tinst_require_finite_matrix(const arma::mat& x, const char* name) {
  if (!x.is_finite()) {
    Rcpp::stop("`%s` must contain only finite values.", name);
  }
}

void tinst_require_controls(const double tol, const int max_iter,
                            const double zero_tol) {
  if (!std::isfinite(tol) || tol <= 0.0 || max_iter < 1 ||
      !std::isfinite(zero_tol) || zero_tol < 0.0) {
    Rcpp::stop("Invalid tINST convergence controls.");
  }
}

double tinst_checked_double(const long double value, const char* quantity) {
  const long double maximum = static_cast<long double>(
    std::numeric_limits<double>::max()
  );
  if (!std::isfinite(value) || value > maximum || value < -maximum) {
    Rcpp::stop(
      "tINST produced a non-finite %s; no numerical floor, ridge, or "
      "absolute-value repair is applied.", quantity
    );
  }
  return static_cast<double>(value);
}

// A common affine transformation is algebraically neutral for tINST.  It
// keeps all subsequent residual calculations near unit scale.  Direct
// differences are preferred because they preserve small residuals under a
// large common translation.  The normalized fallback also handles an
// opposite-sign pair whose direct double-precision difference overflows.
TinstTransform tinst_common_transform(const arma::mat& x,
                                      const arma::mat& y) {
  const arma::uword p = x.n_cols;
  TinstTransform answer;
  answer.x.set_size(x.n_rows, p);
  answer.y.set_size(y.n_rows, p);
  answer.base.set_size(p);
  answer.scale.set_size(p);
  answer.raw_scale.set_size(p);
  answer.midpoint_normalized.set_size(p);
  answer.spread_normalized.set_size(p);
  answer.normalized_fallback.assign(p, false);

  const long double double_max = static_cast<long double>(
    std::numeric_limits<double>::max()
  );
  for (arma::uword j = 0; j < p; ++j) {
    const double reference = x(0, j);
    std::vector<long double> x_difference(x.n_rows);
    std::vector<long double> y_difference(y.n_rows);
    long double minimum = 0.0L;
    long double maximum = 0.0L;
    bool direct_ok = true;
    for (arma::uword i = 0; i < x.n_rows; ++i) {
      const long double difference = static_cast<long double>(x(i, j)) -
        static_cast<long double>(reference);
      x_difference[i] = difference;
      direct_ok = direct_ok && std::isfinite(difference);
      minimum = std::min(minimum, difference);
      maximum = std::max(maximum, difference);
    }
    for (arma::uword i = 0; i < y.n_rows; ++i) {
      const long double difference = static_cast<long double>(y(i, j)) -
        static_cast<long double>(reference);
      y_difference[i] = difference;
      direct_ok = direct_ok && std::isfinite(difference);
      minimum = std::min(minimum, difference);
      maximum = std::max(maximum, difference);
    }

    long double midpoint = minimum / 2.0L + maximum / 2.0L;
    long double spread = std::max(
      std::abs(minimum - midpoint), std::abs(maximum - midpoint)
    );
    direct_ok = direct_ok && std::isfinite(midpoint) &&
      std::isfinite(spread) && spread <= double_max;

    if (direct_ok) {
      if (spread == 0.0L) {
        spread = 1.0L;
      }
      for (arma::uword i = 0; i < x.n_rows; ++i) {
        answer.x(i, j) = static_cast<double>(
          (x_difference[i] - midpoint) / spread
        );
      }
      for (arma::uword i = 0; i < y.n_rows; ++i) {
        answer.y(i, j) = static_cast<double>(
          (y_difference[i] - midpoint) / spread
        );
      }
      const long double base = static_cast<long double>(reference) + midpoint;
      answer.base(j) = std::isfinite(base) && std::abs(base) <= double_max ?
        static_cast<double>(base) : reference;
      answer.scale(j) = static_cast<double>(spread);
      answer.raw_scale(j) = 1.0;
      answer.midpoint_normalized(j) = 0.0;
      answer.spread_normalized(j) = 1.0;
      continue;
    }

    answer.normalized_fallback[j] = true;
    double raw_scale = 0.0;
    for (arma::uword i = 0; i < x.n_rows; ++i) {
      raw_scale = std::max(raw_scale, std::abs(x(i, j)));
    }
    for (arma::uword i = 0; i < y.n_rows; ++i) {
      raw_scale = std::max(raw_scale, std::abs(y(i, j)));
    }
    if (!std::isfinite(raw_scale) || raw_scale <= 0.0) {
      raw_scale = 1.0;
    }
    double normalized_minimum = 1.0;
    double normalized_maximum = -1.0;
    for (arma::uword i = 0; i < x.n_rows; ++i) {
      const double value = x(i, j) / raw_scale;
      normalized_minimum = std::min(normalized_minimum, value);
      normalized_maximum = std::max(normalized_maximum, value);
    }
    for (arma::uword i = 0; i < y.n_rows; ++i) {
      const double value = y(i, j) / raw_scale;
      normalized_minimum = std::min(normalized_minimum, value);
      normalized_maximum = std::max(normalized_maximum, value);
    }
    const double normalized_midpoint = normalized_minimum / 2.0 +
      normalized_maximum / 2.0;
    double normalized_spread = std::max(
      std::abs(normalized_minimum - normalized_midpoint),
      std::abs(normalized_maximum - normalized_midpoint)
    );
    if (normalized_spread == 0.0) {
      normalized_spread = 1.0;
    }
    for (arma::uword i = 0; i < x.n_rows; ++i) {
      answer.x(i, j) = (x(i, j) / raw_scale - normalized_midpoint) /
        normalized_spread;
    }
    for (arma::uword i = 0; i < y.n_rows; ++i) {
      answer.y(i, j) = (y(i, j) / raw_scale - normalized_midpoint) /
        normalized_spread;
    }
    answer.base(j) = raw_scale * normalized_midpoint;
    const long double physical_scale =
      static_cast<long double>(raw_scale) * normalized_spread;
    answer.scale(j) = physical_scale <= double_max ?
      static_cast<double>(physical_scale) : R_PosInf;
    answer.raw_scale(j) = raw_scale;
    answer.midpoint_normalized(j) = normalized_midpoint;
    answer.spread_normalized(j) = normalized_spread;
  }

  if (!answer.x.is_finite() || !answer.y.is_finite()) {
    Rcpp::stop("tINST internal common standardisation failed.");
  }
  return answer;
}

arma::vec tinst_restore_location(const arma::rowvec& location,
                                 const TinstTransform& transform) {
  arma::vec answer(location.n_elem);
  for (arma::uword j = 0; j < location.n_elem; ++j) {
    long double restored;
    if (transform.normalized_fallback[j]) {
      restored = static_cast<long double>(transform.raw_scale(j)) *
        (static_cast<long double>(transform.midpoint_normalized(j)) +
         static_cast<long double>(transform.spread_normalized(j)) *
           static_cast<long double>(location(j)));
    } else {
      restored = static_cast<long double>(transform.base(j)) +
        static_cast<long double>(transform.scale(j)) *
          static_cast<long double>(location(j));
    }
    const long double maximum = static_cast<long double>(
      std::numeric_limits<double>::max()
    );
    answer(j) = std::isfinite(restored) && std::abs(restored) <= maximum ?
      static_cast<double>(restored) :
      (restored < 0.0L ? R_NegInf : R_PosInf);
  }
  return answer;
}

arma::vec tinst_restore_difference(const arma::rowvec& difference,
                                   const TinstTransform& transform) {
  arma::vec answer(difference.n_elem);
  for (arma::uword j = 0; j < difference.n_elem; ++j) {
    long double restored;
    if (transform.normalized_fallback[j]) {
      restored = static_cast<long double>(transform.raw_scale(j)) *
        static_cast<long double>(transform.spread_normalized(j)) *
        static_cast<long double>(difference(j));
    } else {
      restored = static_cast<long double>(transform.scale(j)) *
        static_cast<long double>(difference(j));
    }
    const long double maximum = static_cast<long double>(
      std::numeric_limits<double>::max()
    );
    answer(j) = std::isfinite(restored) && std::abs(restored) <= maximum ?
      static_cast<double>(restored) :
      (restored < 0.0L ? R_NegInf : R_PosInf);
  }
  return answer;
}

// Returns false for a residual at or below zero_tol.  Callers fitting an
// inverse-norm equation reject such a residual.  The final cross statistic,
// whose definition explicitly includes I(x != 0), instead maps it to zero.
bool tinst_direction_radius(const arma::rowvec& observation,
                            const arma::rowvec& location,
                            const arma::vec& diagonal,
                            const double zero_tol,
                            arma::rowvec& direction,
                            double& radius) {
  const arma::uword p = observation.n_elem;
  arma::rowvec epsilon(p);
  for (arma::uword j = 0; j < p; ++j) {
    if (!std::isfinite(diagonal(j)) || diagonal(j) <= 0.0) {
      Rcpp::stop(
        "A tINST diagonal iterate is not finite and strictly positive; "
        "no ridge is applied."
      );
    }
    epsilon(j) = (observation(j) - location(j)) /
      std::sqrt(diagonal(j));
  }
  if (!epsilon.is_finite()) {
    Rcpp::stop("A tINST diagonal standardisation produced non-finite values.");
  }
  const double maximum = arma::abs(epsilon).max();
  direction.zeros(p);
  if (maximum == 0.0) {
    radius = 0.0;
    return false;
  }
  const arma::rowvec normalized = epsilon / maximum;
  const double normalized_norm = std::sqrt(arma::dot(normalized, normalized));
  if (!std::isfinite(normalized_norm) || normalized_norm <= 0.0) {
    Rcpp::stop("A tINST residual has an invalid Euclidean norm.");
  }
  if (zero_tol > 0.0 && maximum <= zero_tol / normalized_norm) {
    radius = maximum * normalized_norm;
    return false;
  }
  radius = maximum * normalized_norm;
  if (!std::isfinite(radius) || radius <= 0.0) {
    Rcpp::stop("A tINST residual has a non-finite Euclidean norm.");
  }
  direction = normalized / normalized_norm;
  return true;
}

arma::uword tinst_subset_size(const arma::mat& data, const int omitted) {
  return data.n_rows - (omitted >= 0 ? 1 : 0);
}

arma::rowvec tinst_subset_mean(const arma::mat& data, const int omitted) {
  arma::rowvec total(data.n_cols, arma::fill::zeros);
  arma::uword count = 0;
  for (arma::uword i = 0; i < data.n_rows; ++i) {
    if (static_cast<int>(i) == omitted) {
      continue;
    }
    total += data.row(i);
    ++count;
  }
  return total / static_cast<double>(count);
}

arma::vec tinst_subset_variance(const arma::mat& data, const int omitted,
                                const arma::rowvec& mean,
                                const std::string& label) {
  const arma::uword count = tinst_subset_size(data, omitted);
  arma::vec variance(data.n_cols, arma::fill::zeros);
  for (arma::uword i = 0; i < data.n_rows; ++i) {
    if (static_cast<int>(i) == omitted) {
      continue;
    }
    const arma::rowvec residual = data.row(i) - mean;
    variance += arma::square(residual).t();
  }
  variance /= static_cast<double>(count - 1);
  if (!variance.is_finite() || arma::any(variance <= 0.0)) {
    Rcpp::stop(
      "%s requires every sample marginal variance to be finite and strictly "
      "positive; no ridge is applied.", label.c_str()
    );
  }
  return variance;
}

bool tinst_subset_directions(const arma::mat& data, const int omitted,
                             const arma::rowvec& location,
                             const arma::vec& diagonal,
                             const double zero_tol,
                             const std::string& label,
                             arma::mat& directions,
                             arma::vec& radii,
                             const bool fail_on_zero = true) {
  const arma::uword count = tinst_subset_size(data, omitted);
  directions.set_size(count, data.n_cols);
  radii.set_size(count);
  arma::uword position = 0;
  for (arma::uword i = 0; i < data.n_rows; ++i) {
    if (static_cast<int>(i) == omitted) {
      continue;
    }
    arma::rowvec direction;
    double radius = 0.0;
    if (!tinst_direction_radius(
          data.row(i), location, diagonal, zero_tol, direction, radius
        )) {
      if (fail_on_zero) {
        Rcpp::stop(
          "%s is undefined because an included observation has a zero "
          "standardised residual before its iteration became stable; no "
          "perturbation is applied.", label.c_str()
        );
      }
      directions.row(position).zeros();
      radii(position) = 0.0;
      return false;
    }
    directions.row(position) = direction;
    radii(position) = radius;
    ++position;
  }
  return true;
}

TinstFit tinst_joint_fit(const arma::mat& data, const int omitted,
                         const double tol, const int max_iter,
                         const double zero_tol,
                         const std::string& label) {
  const double p = static_cast<double>(data.n_cols);
  TinstFit fit;
  fit.location = tinst_subset_mean(data, omitted);
  fit.diagonal = tinst_subset_variance(
    data, omitted, fit.location, label
  );
  const arma::vec initial_diagonal = fit.diagonal;

  for (int update = 0; update < max_iter; ++update) {
    arma::mat directions;
    arma::vec radii;
    tinst_subset_directions(
      data, omitted, fit.location, fit.diagonal, zero_tol, label,
      directions, radii
    );
    const double minimum_radius = radii.min();
    fit.minimum_residual_distance = std::min(
      fit.minimum_residual_distance, minimum_radius
    );
    const arma::vec relative_inverse_radius = minimum_radius / radii;
    const double denominator = arma::accu(relative_inverse_radius);
    if (!std::isfinite(denominator) || denominator <= 0.0) {
      Rcpp::stop(
        "%s has an invalid location-update denominator; no repair is "
        "applied.", label.c_str()
      );
    }
    const arma::rowvec standardized_step = arma::sum(directions, 0) *
      (minimum_radius / denominator);
    const arma::rowvec next_location = fit.location +
      arma::sqrt(fit.diagonal).t() % standardized_step;
    const arma::vec next_diagonal = p * fit.diagonal %
      arma::mean(arma::square(directions), 0).t();
    if (!next_location.is_finite() || !next_diagonal.is_finite() ||
        arma::any(next_diagonal <= 0.0)) {
      Rcpp::stop(
        "%s produced a non-finite or non-positive iterate; no ridge is "
        "applied.", label.c_str()
      );
    }
    fit.location_change = arma::abs(
      (next_location - fit.location) / arma::sqrt(initial_diagonal).t()
    ).max();
    fit.diagonal_change = arma::abs(
      arma::log(next_diagonal / fit.diagonal)
    ).max();
    fit.relative_update = std::max(
      fit.location_change, fit.diagonal_change
    );
    fit.location = next_location;
    fit.diagonal = next_diagonal;
    fit.iterations = update + 1;
    if (fit.relative_update <= tol) {
      fit.iteration_stable = true;
      break;
    }
  }

  // Evaluate the score only as a diagnostic.  Stability is certified solely
  // by the relative iterate update, matching the paper's algorithmic wording.
  arma::mat final_directions;
  arma::vec final_radii;
  const bool final_score_defined = tinst_subset_directions(
    data, omitted, fit.location, fit.diagonal, zero_tol, label,
    final_directions, final_radii, fit.iteration_stable ? false : true
  );
  if (final_score_defined) {
    fit.minimum_residual_distance = std::min(
      fit.minimum_residual_distance, final_radii.min()
    );
    const arma::rowvec mean_direction = arma::mean(final_directions, 0);
    fit.location_residual = arma::norm(mean_direction, 2);
    const arma::vec shape_equation = p *
      arma::mean(arma::square(final_directions), 0).t();
    fit.diagonal_residual = arma::abs(shape_equation - 1.0).max();
    fit.score_residual = std::max(
      fit.location_residual, fit.diagonal_residual
    );
  } else {
    fit.minimum_residual_distance = 0.0;
    fit.location_residual = R_PosInf;
    fit.diagonal_residual = R_PosInf;
    fit.score_residual = R_PosInf;
  }
  return fit;
}

TinstFit tinst_weighted_location_fit(const arma::mat& data,
                                     const int omitted,
                                     const arma::vec& diagonal,
                                     const double tol,
                                     const int max_iter,
                                     const double zero_tol,
                                     const std::string& label) {
  TinstFit fit;
  fit.location = tinst_subset_mean(data, omitted);
  fit.diagonal = diagonal;
  const arma::vec location_scale = diagonal;

  for (int update = 0; update < max_iter; ++update) {
    arma::mat directions;
    arma::vec radii;
    tinst_subset_directions(
      data, omitted, fit.location, diagonal, zero_tol, label,
      directions, radii
    );
    const double minimum_radius = radii.min();
    fit.minimum_residual_distance = std::min(
      fit.minimum_residual_distance, minimum_radius
    );
    const arma::vec inverse_radius_weights = minimum_radius / radii;
    const double equation_denominator = arma::accu(inverse_radius_weights);
    const arma::rowvec equation_numerator =
      inverse_radius_weights.t() * directions;
    if (!std::isfinite(equation_denominator) ||
        equation_denominator <= 0.0 || !equation_numerator.is_finite()) {
      Rcpp::stop(
        "%s has an invalid inverse-norm estimating equation; no repair is "
        "applied.", label.c_str()
      );
    }
    const double update_denominator = arma::dot(
      inverse_radius_weights, inverse_radius_weights
    );
    if (!std::isfinite(update_denominator) || update_denominator <= 0.0) {
      Rcpp::stop(
        "%s has an invalid inverse-norm update denominator; no repair is "
        "applied.", label.c_str()
      );
    }
    const arma::rowvec standardized_step = minimum_radius *
      equation_numerator / update_denominator;
    const arma::rowvec next_location = fit.location +
      arma::sqrt(diagonal).t() % standardized_step;
    if (!next_location.is_finite()) {
      Rcpp::stop("%s produced a non-finite location iterate.", label.c_str());
    }
    fit.location_change = arma::abs(
      (next_location - fit.location) / arma::sqrt(location_scale).t()
    ).max();
    fit.relative_update = fit.location_change;
    fit.location = next_location;
    fit.iterations = update + 1;
    if (fit.relative_update <= tol) {
      fit.iteration_stable = true;
      break;
    }
  }


  arma::mat final_directions;
  arma::vec final_radii;
  const bool final_score_defined = tinst_subset_directions(
    data, omitted, fit.location, diagonal, zero_tol, label,
    final_directions, final_radii, fit.iteration_stable ? false : true
  );
  if (final_score_defined) {
    const double minimum_radius = final_radii.min();
    fit.minimum_residual_distance = std::min(
      fit.minimum_residual_distance, minimum_radius
    );
    const arma::vec inverse_radius_weights = minimum_radius / final_radii;
    const double denominator = arma::accu(inverse_radius_weights);
    const arma::rowvec numerator =
      inverse_radius_weights.t() * final_directions;
    fit.location_residual = arma::norm(numerator / denominator, 2);
    fit.score_residual = fit.location_residual;
  } else {
    fit.minimum_residual_distance = 0.0;
    fit.location_residual = R_PosInf;
    fit.score_residual = R_PosInf;
  }
  return fit;
}

Rcpp::List tinst_fit_diagnostic(const TinstFit& fit) {
  return Rcpp::List::create(
    Rcpp::Named("iterations") = fit.iterations,
    Rcpp::Named("iteration.stable") = fit.iteration_stable,
    Rcpp::Named("relative.update") = fit.relative_update,
    Rcpp::Named("score.residual") = fit.score_residual,
    Rcpp::Named("minimum.residual.distance") =
      fit.minimum_residual_distance,
    Rcpp::Named("location.relative.update") = fit.location_change,
    Rcpp::Named("log.diagonal.relative.update") = fit.diagonal_change,
    Rcpp::Named("location.score.residual") = fit.location_residual,
    Rcpp::Named("diagonal.score.residual") = fit.diagonal_residual,
    Rcpp::Named("convergence.basis") =
      "relative update as specified by paper"
  );
}

Rcpp::List tinst_fit_diagnostics(const std::vector<TinstFit>& fits) {
  const R_xlen_t n = static_cast<R_xlen_t>(fits.size());
  Rcpp::IntegerVector iterations(n);
  Rcpp::LogicalVector iteration_stable(n);
  Rcpp::NumericVector relative_update(n);
  Rcpp::NumericVector score_residual(n);
  Rcpp::NumericVector minimum_residual_distance(n);
  Rcpp::NumericVector location_change(n);
  Rcpp::NumericVector diagonal_change(n);
  Rcpp::NumericVector location_residual(n);
  Rcpp::NumericVector diagonal_residual(n);
  bool all_stable = true;
  double worst_relative_update = 0.0;
  double worst_score_residual = 0.0;
  double smallest_residual_distance = R_PosInf;
  int maximum_iterations = 0;
  for (R_xlen_t i = 0; i < n; ++i) {
    iterations[i] = fits[i].iterations;
    iteration_stable[i] = fits[i].iteration_stable;
    relative_update[i] = fits[i].relative_update;
    score_residual[i] = fits[i].score_residual;
    minimum_residual_distance[i] = fits[i].minimum_residual_distance;
    location_change[i] = fits[i].location_change;
    diagonal_change[i] = fits[i].diagonal_change;
    location_residual[i] = fits[i].location_residual;
    diagonal_residual[i] = fits[i].diagonal_residual;
    all_stable = all_stable && fits[i].iteration_stable;
    worst_relative_update = std::max(
      worst_relative_update, fits[i].relative_update
    );
    worst_score_residual = std::max(
      worst_score_residual, fits[i].score_residual
    );
    smallest_residual_distance = std::min(
      smallest_residual_distance, fits[i].minimum_residual_distance
    );
    maximum_iterations = std::max(maximum_iterations, fits[i].iterations);
  }
  return Rcpp::List::create(
    Rcpp::Named("iterations") = iterations,
    Rcpp::Named("iteration.stable") = iteration_stable,
    Rcpp::Named("relative.update") = relative_update,
    Rcpp::Named("score.residual") = score_residual,
    Rcpp::Named("minimum.residual.distance") = minimum_residual_distance,
    Rcpp::Named("location.relative.update") = location_change,
    Rcpp::Named("log.diagonal.relative.update") = diagonal_change,
    Rcpp::Named("location.score.residual") = location_residual,
    Rcpp::Named("diagonal.score.residual") = diagonal_residual,
    Rcpp::Named("all.iteration.stable") = all_stable,
    Rcpp::Named("worst.relative.update") = worst_relative_update,
    Rcpp::Named("worst.score.residual") = worst_score_residual,
    Rcpp::Named("smallest.residual.distance") =
      smallest_residual_distance,
    Rcpp::Named("maximum.iterations") = maximum_iterations,
    Rcpp::Named("convergence.basis") =
      "relative update as specified by paper"
  );
}

bool tinst_all_iteration_stable(const TinstFit& full1,
                                const TinstFit& full2,
                                const std::vector<TinstFit>& joint1,
                                const std::vector<TinstFit>& joint2,
                                const std::vector<TinstFit>& weighted1,
                                const std::vector<TinstFit>& weighted2,
                                int& failures) {
  failures = (!full1.iteration_stable) + (!full2.iteration_stable);
  for (const TinstFit& fit : joint1) failures += !fit.iteration_stable;
  for (const TinstFit& fit : joint2) failures += !fit.iteration_stable;
  for (const TinstFit& fit : weighted1) failures += !fit.iteration_stable;
  for (const TinstFit& fit : weighted2) failures += !fit.iteration_stable;
  return failures == 0;
}

}  // namespace


//' Huang--Liu--Zhou--Feng two-sample inverse norm sign kernel
//'
//' @param x,y Numeric observation-by-variable matrices.
//' @param tol Positive estimating-equation tolerance.
//' @param max_iter Positive maximum number of updates for every fit.
//' @param zero_tol Non-negative standardized zero-residual tolerance.
//' @return Internal list of tINST statistic and convergence components.
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List cpp_tinst_two_sample(const arma::mat& x,
                                const arma::mat& y,
                                const double tol,
                                const int max_iter,
                                const double zero_tol) {
  tinst_require_finite_matrix(x, "x");
  tinst_require_finite_matrix(y, "y");
  tinst_require_controls(tol, max_iter, zero_tol);
  const arma::uword n1 = x.n_rows;
  const arma::uword n2 = y.n_rows;
  const arma::uword p = x.n_cols;
  if (n1 < 3 || n2 < 3) {
    Rcpp::stop("tINST requires at least three observations in each group.");
  }
  if (p < 1 || y.n_cols != p) {
    Rcpp::stop("`x` and `y` must have the same positive number of columns.");
  }

  const TinstTransform transform = tinst_common_transform(x, y);
  const arma::mat& xs = transform.x;
  const arma::mat& ys = transform.y;
  const TinstFit full1 = tinst_joint_fit(
    xs, -1, tol, max_iter, zero_tol, "The group-1 full-sample tINST fit"
  );
  const TinstFit full2 = tinst_joint_fit(
    ys, -1, tol, max_iter, zero_tol, "The group-2 full-sample tINST fit"
  );

  arma::mat diagonal1(n1, p);
  arma::mat diagonal2(n2, p);
  arma::mat joint_location1(n1, p);
  arma::mat joint_location2(n2, p);
  arma::mat weighted_location1(n1, p);
  arma::mat weighted_location2(n2, p);
  arma::mat own_direction1(n1, p);
  arma::mat own_direction2(n2, p);
  arma::vec own_radius1(n1);
  arma::vec own_radius2(n2);
  std::vector<TinstFit> joint1(n1), joint2(n2), weighted1(n1), weighted2(n2);

  for (arma::uword i = 0; i < n1; ++i) {
    const std::string joint_label = "The group-1 leave-one-out tINST joint fit";
    const std::string weighted_label =
      "The group-1 leave-one-out tINST weighted-location fit";
    joint1[i] = tinst_joint_fit(
      xs, static_cast<int>(i), tol, max_iter, zero_tol, joint_label
    );
    weighted1[i] = tinst_weighted_location_fit(
      xs, static_cast<int>(i), joint1[i].diagonal,
      tol, max_iter, zero_tol, weighted_label
    );
    diagonal1.row(i) = joint1[i].diagonal.t();
    joint_location1.row(i) = joint1[i].location;
    weighted_location1.row(i) = weighted1[i].location;
    arma::rowvec direction;
    double radius = 0.0;
    if (!tinst_direction_radius(
          xs.row(i), weighted1[i].location, joint1[i].diagonal,
          zero_tol, direction, radius
        )) {
      Rcpp::stop(
        "The group-1 tINST nuisance estimates are undefined because an "
        "observation equals its leave-one-out weighted location."
      );
    }
    own_direction1.row(i) = direction;
    own_radius1(i) = radius;
  }
  for (arma::uword i = 0; i < n2; ++i) {
    const std::string joint_label = "The group-2 leave-one-out tINST joint fit";
    const std::string weighted_label =
      "The group-2 leave-one-out tINST weighted-location fit";
    joint2[i] = tinst_joint_fit(
      ys, static_cast<int>(i), tol, max_iter, zero_tol, joint_label
    );
    weighted2[i] = tinst_weighted_location_fit(
      ys, static_cast<int>(i), joint2[i].diagonal,
      tol, max_iter, zero_tol, weighted_label
    );
    diagonal2.row(i) = joint2[i].diagonal.t();
    joint_location2.row(i) = joint2[i].location;
    weighted_location2.row(i) = weighted2[i].location;
    arma::rowvec direction;
    double radius = 0.0;
    if (!tinst_direction_radius(
          ys.row(i), weighted2[i].location, joint2[i].diagonal,
          zero_tol, direction, radius
        )) {
      Rcpp::stop(
        "The group-2 tINST nuisance estimates are undefined because an "
        "observation equals its leave-one-out weighted location."
      );
    }
    own_direction2.row(i) = direction;
    own_radius2(i) = radius;
  }

  long double cross_sum = 0.0L;
  double cross_zero_vectors = 0.0;
  for (arma::uword i = 0; i < n1; ++i) {
    for (arma::uword j = 0; j < n2; ++j) {
      arma::rowvec direction1;
      arma::rowvec direction2;
      double radius1 = 0.0;
      double radius2 = 0.0;
      const bool nonzero1 = tinst_direction_radius(
        xs.row(i), weighted2[j].location, joint1[i].diagonal,
        zero_tol, direction1, radius1
      );
      const bool nonzero2 = tinst_direction_radius(
        ys.row(j), weighted1[i].location, joint2[j].diagonal,
        zero_tol, direction2, radius2
      );
      if (!nonzero1 || !nonzero2) {
        cross_zero_vectors += (!nonzero1) + (!nonzero2);
        continue;
      }
      long double directional_inner = 0.0L;
      for (arma::uword ell = 0; ell < p; ++ell) {
        directional_inner += static_cast<long double>(direction1(ell)) *
          static_cast<long double>(direction2(ell));
      }
      cross_sum += directional_inner /
        (static_cast<long double>(radius1) * radius2);
    }
  }
  const long double n1_ld = static_cast<long double>(n1);
  const long double n2_ld = static_cast<long double>(n2);
  const long double p_ld = static_cast<long double>(p);
  const double statistic = tinst_checked_double(
    -cross_sum / (n1_ld * n2_ld), "raw statistic"
  );

  long double inverse_radius_square_sum1 = 0.0L;
  long double inverse_radius_square_sum2 = 0.0L;
  for (arma::uword i = 0; i < n1; ++i) {
    const long double radius = own_radius1(i);
    inverse_radius_square_sum1 += 1.0L / (radius * radius);
  }
  for (arma::uword i = 0; i < n2; ++i) {
    const long double radius = own_radius2(i);
    inverse_radius_square_sum2 += 1.0L / (radius * radius);
  }
  const double nu12 = tinst_checked_double(
    inverse_radius_square_sum1 / n1_ld, "group-1 inverse radial moment"
  );
  const double nu22 = tinst_checked_double(
    inverse_radius_square_sum2 / n2_ld, "group-2 inverse radial moment"
  );
  // For omega(r)=1/r, nu_k2 and c_k have the same sample formula.  They are
  // kept as separately named quantities to mirror the published variance.
  const double c1 = nu12;
  const double c2 = nu22;
  if (!(c1 > 0.0) || !(c2 > 0.0)) {
    Rcpp::stop("tINST inverse radial moment estimates must be positive.");
  }

  const arma::vec ratio1 = arma::sqrt(full1.diagonal / full2.diagonal);
  const arma::vec ratio2 = arma::sqrt(full2.diagonal / full1.diagonal);
  if (!ratio1.is_finite() || !ratio2.is_finite()) {
    Rcpp::stop("tINST full-sample diagonal ratios are non-finite.");
  }
  long double a1_sum = 0.0L;
  long double a2_sum = 0.0L;
  long double a12_sum = 0.0L;
  for (arma::uword i = 0; i < n1; ++i) {
    for (arma::uword j = 0; j < n1; ++j) {
      if (i == j) continue;
      long double inner = 0.0L;
      for (arma::uword ell = 0; ell < p; ++ell) {
        inner += static_cast<long double>(own_direction1(j, ell)) *
          static_cast<long double>(ratio1(ell)) *
          static_cast<long double>(own_direction1(i, ell));
      }
      a1_sum += inner * inner;
    }
  }
  // The outer limit is n2.  The n1 limit in the submitted manuscript is a
  // typographical error and is intentionally not reproduced here.
  for (arma::uword i = 0; i < n2; ++i) {
    for (arma::uword j = 0; j < n2; ++j) {
      if (i == j) continue;
      long double inner = 0.0L;
      for (arma::uword ell = 0; ell < p; ++ell) {
        inner += static_cast<long double>(own_direction2(j, ell)) *
          static_cast<long double>(ratio2(ell)) *
          static_cast<long double>(own_direction2(i, ell));
      }
      a2_sum += inner * inner;
    }
  }
  for (arma::uword i = 0; i < n1; ++i) {
    for (arma::uword j = 0; j < n2; ++j) {
      long double inner = 0.0L;
      for (arma::uword ell = 0; ell < p; ++ell) {
        inner += static_cast<long double>(own_direction1(i, ell)) *
          static_cast<long double>(own_direction2(j, ell));
      }
      a12_sum += inner * inner;
    }
  }
  const long double p_square = p_ld * p_ld;
  const double a1 = tinst_checked_double(
    p_square * a1_sum / (n1_ld * (n1_ld - 1.0L)), "A1 trace estimate"
  );
  const double a2 = tinst_checked_double(
    p_square * a2_sum / (n2_ld * (n2_ld - 1.0L)), "A2 trace estimate"
  );
  const double a12 = tinst_checked_double(
    p_square * a12_sum / (n1_ld * n2_ld), "A12 trace estimate"
  );

  const long double nu12_ld = nu12;
  const long double nu22_ld = nu22;
  const long double c1_ld = c1;
  const long double c2_ld = c2;
  const long double variance_term1 =
    2.0L * nu12_ld * nu12_ld * c2_ld * c2_ld /
    (c1_ld * c1_ld) * static_cast<long double>(a1) /
    (n1_ld * (n1_ld - 1.0L) * p_square);
  const long double variance_term2 =
    2.0L * nu22_ld * nu22_ld * c1_ld * c1_ld /
    (c2_ld * c2_ld) * static_cast<long double>(a2) /
    (n2_ld * (n2_ld - 1.0L) * p_square);
  const long double variance_term12 =
    4.0L * nu12_ld * nu22_ld * static_cast<long double>(a12) /
    (n1_ld * n2_ld * p_square);
  const long double variance_ld = variance_term1 + variance_term2 +
    variance_term12;
  if (!std::isfinite(variance_ld) || variance_ld <= 0.0L) {
    Rcpp::stop(
      "tINST requires its published variance estimate to be finite and "
      "strictly positive; no absolute-value or floor repair is applied."
    );
  }
  const double variance = tinst_checked_double(
    variance_ld, "variance estimate"
  );
  const double z = statistic / std::sqrt(variance);
  if (!std::isfinite(z)) {
    Rcpp::stop("The standardized tINST statistic is not finite.");
  }

  int stability_failures = 0;
  const bool all_iteration_stable = tinst_all_iteration_stable(
    full1, full2, joint1, joint2, weighted1, weighted2,
    stability_failures
  );
  const arma::rowvec mean1_scaled = arma::mean(xs, 0);
  const arma::rowvec mean2_scaled = arma::mean(ys, 0);
  const arma::vec mean1 = tinst_restore_location(mean1_scaled, transform);
  const arma::vec mean2 = tinst_restore_location(mean2_scaled, transform);
  const arma::vec difference = tinst_restore_difference(
    mean1_scaled - mean2_scaled, transform
  );

  return Rcpp::List::create(
    Rcpp::Named("T") = statistic,
    Rcpp::Named("variance") = variance,
    Rcpp::Named("z") = z,
    Rcpp::Named("variance_term1") = static_cast<double>(variance_term1),
    Rcpp::Named("variance_term2") = static_cast<double>(variance_term2),
    Rcpp::Named("variance_term12") = static_cast<double>(variance_term12),
    Rcpp::Named("nu12") = nu12,
    Rcpp::Named("nu22") = nu22,
    Rcpp::Named("c1") = c1,
    Rcpp::Named("c2") = c2,
    Rcpp::Named("A1") = a1,
    Rcpp::Named("A2") = a2,
    Rcpp::Named("A12") = a12,
    Rcpp::Named("n1") = static_cast<double>(n1),
    Rcpp::Named("n2") = static_cast<double>(n2),
    Rcpp::Named("p") = static_cast<double>(p),
    Rcpp::Named("mean_x") = mean1,
    Rcpp::Named("mean_y") = mean2,
    Rcpp::Named("difference") = difference,
    Rcpp::Named("full_diagonal1_scaled") = full1.diagonal,
    Rcpp::Named("full_diagonal2_scaled") = full2.diagonal,
    Rcpp::Named("leaveout_diagonal1_scaled") = diagonal1,
    Rcpp::Named("leaveout_diagonal2_scaled") = diagonal2,
    Rcpp::Named("joint_location1_scaled") = joint_location1,
    Rcpp::Named("joint_location2_scaled") = joint_location2,
    Rcpp::Named("weighted_location1_scaled") = weighted_location1,
    Rcpp::Named("weighted_location2_scaled") = weighted_location2,
    Rcpp::Named("own_direction1") = own_direction1,
    Rcpp::Named("own_direction2") = own_direction2,
    Rcpp::Named("own_radius1") = own_radius1,
    Rcpp::Named("own_radius2") = own_radius2,
    Rcpp::Named("full_fit1") = tinst_fit_diagnostic(full1),
    Rcpp::Named("full_fit2") = tinst_fit_diagnostic(full2),
    Rcpp::Named("joint_fits1") = tinst_fit_diagnostics(joint1),
    Rcpp::Named("joint_fits2") = tinst_fit_diagnostics(joint2),
    Rcpp::Named("weighted_fits1") = tinst_fit_diagnostics(weighted1),
    Rcpp::Named("weighted_fits2") = tinst_fit_diagnostics(weighted2),
    Rcpp::Named("all_iteration_stable") = all_iteration_stable,
    Rcpp::Named("stability_failures") = stability_failures,
    Rcpp::Named("cross_zero_vectors") = cross_zero_vectors,
    Rcpp::Named("internal_anchor") = transform.base.t(),
    Rcpp::Named("internal_column_scale") = transform.scale.t(),
    Rcpp::Named("internal_overflow_fallback") =
      std::any_of(transform.normalized_fallback.begin(),
                  transform.normalized_fallback.end(),
                  [](bool value) { return value; })
  );
}
