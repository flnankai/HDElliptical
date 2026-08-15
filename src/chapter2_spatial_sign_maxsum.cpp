// Scaled spatial median and spatial-sign max statistic.
//
// This is the full-sample diagonal Hettmansperger--Randles system used by
// Liu, Feng, Zhao and Wang.  The diagonal is identified only up to a common
// positive multiplier; internally we impose unit geometric mean.  No ridge,
// floor, pseudoinverse, radius perturbation, or iterate substitution is used.

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp17)]]

#include <RcppArmadillo.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

namespace {

void ssmax_neumaier_add(const long double value,
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

double ssmax_checked_double(const long double value,
                            const char* quantity) {
  const long double maximum = static_cast<long double>(
    std::numeric_limits<double>::max()
  );
  if (!std::isfinite(value) || value > maximum || value < -maximum) {
    Rcpp::stop(
      "The scaled spatial median produced a non-finite %s; no numerical "
      "floor, ridge, perturbation, or pseudoinverse is applied.",
      quantity
    );
  }
  return static_cast<double>(value);
}

arma::vec ssmax_stable_column_mean(const arma::mat& x) {
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
      ssmax_neumaier_add(
        static_cast<long double>(x(i, j) / scale), total, correction
      );
    }
    answer(j) = ssmax_checked_double(
      (total + correction) * static_cast<long double>(scale) /
        static_cast<long double>(x.n_rows),
      "sample mean"
    );
  }
  return answer;
}

struct SsmaxPreparedData {
  arma::mat data;
  arma::vec origin;
  arma::vec sample_mean;
  arma::vec direct_scale;
  arma::vec operand_scale;
  arma::vec origin_normalized;
  arma::vec normalized_scale;
  arma::vec log_scale;
  arma::uvec overflow_fallback;
};

SsmaxPreparedData ssmax_prepare_data(const arma::mat& x,
                                     const arma::vec& target,
                                     const bool center_on_target) {
  const arma::uword n = x.n_rows;
  const arma::uword p = x.n_cols;
  const long double double_max = static_cast<long double>(
    std::numeric_limits<double>::max()
  );
  SsmaxPreparedData prepared;
  prepared.data.set_size(n, p);
  prepared.sample_mean = ssmax_stable_column_mean(x);
  prepared.origin = center_on_target ? target : prepared.sample_mean;
  prepared.direct_scale.zeros(p);
  prepared.operand_scale.ones(p);
  prepared.origin_normalized.zeros(p);
  prepared.normalized_scale.ones(p);
  prepared.log_scale.zeros(p);
  prepared.overflow_fallback.zeros(p);

  for (arma::uword j = 0; j < p; ++j) {
    std::vector<long double> difference(n);
    long double maximum = 0.0L;
    bool direct_ok = true;
    for (arma::uword i = 0; i < n; ++i) {
      difference[i] = static_cast<long double>(x(i, j)) -
        static_cast<long double>(prepared.origin(j));
      direct_ok = direct_ok && std::isfinite(difference[i]);
      maximum = std::max(maximum, std::abs(difference[i]));
    }
    direct_ok = direct_ok && std::isfinite(maximum) &&
      maximum > 0.0L && maximum <= double_max;

    if (direct_ok) {
      const double scale = static_cast<double>(maximum);
      prepared.direct_scale(j) = scale;
      prepared.log_scale(j) = std::log(scale);
      for (arma::uword i = 0; i < n; ++i) {
        prepared.data(i, j) = static_cast<double>(
          difference[i] / maximum
        );
      }
    } else {
      prepared.overflow_fallback(j) = 1u;
      double operand_scale = std::abs(prepared.origin(j));
      for (arma::uword i = 0; i < n; ++i) {
        operand_scale = std::max(operand_scale, std::abs(x(i, j)));
      }
      if (!(operand_scale > 0.0) || !std::isfinite(operand_scale)) {
        Rcpp::stop(
          "The scaled spatial median could not construct a finite common "
          "operand scale for variable %llu.",
          static_cast<unsigned long long>(j + 1)
        );
      }
      prepared.operand_scale(j) = operand_scale;
      prepared.origin_normalized(j) =
        prepared.origin(j) / operand_scale;
      double normalized_maximum = 0.0;
      for (arma::uword i = 0; i < n; ++i) {
        prepared.data(i, j) = x(i, j) / operand_scale -
          prepared.origin_normalized(j);
        normalized_maximum = std::max(
          normalized_maximum, std::abs(prepared.data(i, j))
        );
      }
      if (!(normalized_maximum > 0.0) ||
          !std::isfinite(normalized_maximum)) {
        Rcpp::stop(
          "The scaled spatial median requires positive variation in "
          "variable %llu.",
          static_cast<unsigned long long>(j + 1)
        );
      }
      prepared.normalized_scale(j) = normalized_maximum;
      prepared.data.col(j) /= normalized_maximum;
      prepared.log_scale(j) = std::log(operand_scale) +
        std::log(normalized_maximum);
    }

    double minimum = prepared.data(0, j);
    double column_maximum = prepared.data(0, j);
    for (arma::uword i = 1; i < n; ++i) {
      minimum = std::min(minimum, prepared.data(i, j));
      column_maximum = std::max(column_maximum, prepared.data(i, j));
    }
    if (!(column_maximum > minimum)) {
      Rcpp::stop(
        "The scaled spatial median requires strictly positive marginal "
        "sample variation; variable %llu is degenerate. No ridge is "
        "applied.",
        static_cast<unsigned long long>(j + 1)
      );
    }
  }
  if (!prepared.data.is_finite()) {
    Rcpp::stop("Scaled spatial median internal standardisation failed.");
  }
  return prepared;
}


arma::vec ssmax_restore_location(const arma::vec& location,
                                 const SsmaxPreparedData& prepared) {
  arma::vec answer(location.n_elem);
  for (arma::uword j = 0; j < location.n_elem; ++j) {
    long double value;
    if (prepared.direct_scale(j) > 0.0) {
      value = static_cast<long double>(prepared.origin(j)) +
        static_cast<long double>(prepared.direct_scale(j)) *
          static_cast<long double>(location(j));
    } else {
      value = static_cast<long double>(prepared.operand_scale(j)) *
        (static_cast<long double>(prepared.origin_normalized(j)) +
         static_cast<long double>(prepared.normalized_scale(j)) *
           static_cast<long double>(location(j)));
    }
    answer(j) = ssmax_checked_double(value, "restored location");
  }
  return answer;
}

arma::vec ssmax_restore_difference(const arma::vec& location,
                                   const SsmaxPreparedData& prepared) {
  arma::vec answer(location.n_elem);
  for (arma::uword j = 0; j < location.n_elem; ++j) {
    const bool direct = prepared.direct_scale(j) > 0.0;
    const double base = direct ? prepared.direct_scale(j) :
      prepared.operand_scale(j);
    const double normalized_location = direct ? location(j) :
      prepared.normalized_scale(j) * location(j);
    answer(j) = ssmax_checked_double(
      static_cast<long double>(base) *
        static_cast<long double>(normalized_location),
      "location difference from the centering origin"
    );
  }
  return answer;
}

bool ssmax_direction(const arma::rowvec& value,
                     const arma::vec& location,
                     const arma::vec& diagonal,
                     const double zero_tol,
                     arma::rowvec& direction,
                     double& radius) {
  const arma::uword p = value.n_elem;
  direction.set_size(p);
  double maximum = 0.0;
  for (arma::uword j = 0; j < p; ++j) {
    if (!(diagonal(j) > 0.0) || !std::isfinite(diagonal(j))) {
      Rcpp::stop(
        "A scaled spatial median diagonal iterate is not finite and "
        "strictly positive; no ridge or floor is applied."
      );
    }
    direction(j) = (value(j) - location(j)) /
      std::sqrt(diagonal(j));
    if (!std::isfinite(direction(j))) {
      Rcpp::stop(
        "Scaled spatial median diagonal standardisation produced a "
        "non-finite value."
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
  long double correction = 0.0L;
  for (arma::uword j = 0; j < p; ++j) {
    direction(j) /= maximum;
    const long double component = static_cast<long double>(direction(j));
    ssmax_neumaier_add(component * component, squares, correction);
  }
  const long double norm = std::sqrt(squares + correction);
  if (!(norm > 0.0L) || !std::isfinite(norm)) {
    Rcpp::stop("Scaled spatial median could not normalise a residual.");
  }
  radius = ssmax_checked_double(
    static_cast<long double>(maximum) * norm,
    "standardised residual radius"
  );
  if (!(radius > zero_tol)) {
    direction.zeros();
    return false;
  }
  direction /= static_cast<double>(norm);
  return true;
}

void ssmax_normalise_diagonal(arma::vec& diagonal) {
  long double log_total = 0.0L;
  long double correction = 0.0L;
  for (arma::uword j = 0; j < diagonal.n_elem; ++j) {
    if (!(diagonal(j) > 0.0) || !std::isfinite(diagonal(j))) {
      Rcpp::stop(
        "The scaled spatial median has a non-positive or non-finite "
        "diagonal iterate; no ridge or floor is applied."
      );
    }
    ssmax_neumaier_add(
      std::log(static_cast<long double>(diagonal(j))),
      log_total, correction
    );
  }
  const long double mean_log = (log_total + correction) /
    static_cast<long double>(diagonal.n_elem);
  for (arma::uword j = 0; j < diagonal.n_elem; ++j) {
    diagonal(j) = ssmax_checked_double(
      std::exp(
        std::log(static_cast<long double>(diagonal(j))) - mean_log
      ),
      "identified diagonal iterate"
    );
    if (!(diagonal(j) > 0.0)) {
      Rcpp::stop(
        "Scaled spatial median diagonal identification underflowed; no "
        "floor is applied."
      );
    }
  }
}

struct SsmaxScore {
  arma::rowvec direction_sum;
  arma::vec direction_square_sum;
  arma::mat directions;
  arma::vec radii;
  arma::vec inverse_radii;
  double inverse_radius_sum;
  double minimum_radius;
  double location_residual;
  double diagonal_residual;
  double score_residual;
};

SsmaxScore ssmax_score(const arma::mat& data,
                       const arma::vec& location,
                       const arma::vec& diagonal,
                       const double zero_tol) {
  const arma::uword n = data.n_rows;
  const arma::uword p = data.n_cols;
  arma::rowvec direction_sum(p, arma::fill::zeros);
  arma::vec direction_square_sum(p, arma::fill::zeros);
  arma::mat directions(n, p);
  arma::vec radii(n);
  arma::vec inverse_radii(n);
  long double inverse_total = 0.0L;
  long double inverse_correction = 0.0L;
  for (arma::uword i = 0; i < n; ++i) {
    arma::rowvec direction;
    double radius = 0.0;
    if (!ssmax_direction(
          data.row(i), location, diagonal, zero_tol,
          direction, radius
        )) {
      Rcpp::stop(
        "The scaled spatial median iteration is undefined because "
        "training observation %llu has a standardised residual at or "
        "below `zero_tol`; no perturbation is applied.",
        static_cast<unsigned long long>(i + 1)
      );
    }
    directions.row(i) = direction;
    direction_sum += direction;
    direction_square_sum += arma::square(direction).t();
    radii(i) = radius;
    const long double inverse = 1.0L /
      static_cast<long double>(radius);
    inverse_radii(i) = ssmax_checked_double(
      inverse, "inverse residual radius"
    );
    ssmax_neumaier_add(
      inverse, inverse_total, inverse_correction
    );
  }
  const double inverse_radius_sum = ssmax_checked_double(
    inverse_total + inverse_correction,
    "inverse-radius update denominator"
  );
  if (!(inverse_radius_sum > 0.0)) {
    Rcpp::stop(
      "The scaled spatial median inverse-radius denominator is invalid; "
      "no weight cap is applied."
    );
  }
  const arma::rowvec mean_direction = direction_sum /
    static_cast<double>(n);
  const arma::vec diagonal_equation =
    static_cast<double>(p) * direction_square_sum /
    static_cast<double>(n);
  const double location_residual = arma::abs(mean_direction).max();
  const double diagonal_residual = arma::abs(
    diagonal_equation - 1.0
  ).max();
  return SsmaxScore{
    direction_sum,
    direction_square_sum,
    directions,
    radii,
    inverse_radii,
    inverse_radius_sum,
    radii.min(),
    location_residual,
    diagonal_residual,
    std::max(location_residual, diagonal_residual)
  };
}

struct SsmaxFit {
  arma::vec location;
  arma::vec diagonal;
  int iterations = 0;
  bool iteration_stable = false;
  double relative_update = R_PosInf;
  double location_relative_update = R_PosInf;
  double log_diagonal_update = R_PosInf;
  double score_residual = R_PosInf;
  double location_score_residual = R_PosInf;
  double diagonal_score_residual = R_PosInf;
  double minimum_residual_distance = R_PosInf;
};

SsmaxFit ssmax_fit(const arma::mat& data,
                   const double tol,
                   const int max_iter,
                   const double zero_tol) {
  const arma::uword n = data.n_rows;
  const arma::uword p = data.n_cols;
  SsmaxFit fit;
  fit.location.zeros(p);
  for (arma::uword j = 0; j < p; ++j) {
    long double total = 0.0L;
    long double correction = 0.0L;
    for (arma::uword i = 0; i < n; ++i) {
      ssmax_neumaier_add(
        static_cast<long double>(data(i, j)), total, correction
      );
    }
    fit.location(j) = static_cast<double>(
      (total + correction) / static_cast<long double>(n)
    );
  }
  fit.diagonal.zeros(p);
  for (arma::uword j = 0; j < p; ++j) {
    long double total = 0.0L;
    long double correction = 0.0L;
    for (arma::uword i = 0; i < n; ++i) {
      const long double residual =
        static_cast<long double>(data(i, j)) -
        static_cast<long double>(fit.location(j));
      ssmax_neumaier_add(residual * residual, total, correction);
    }
    fit.diagonal(j) = ssmax_checked_double(
      (total + correction) / static_cast<long double>(n - 1),
      "initial marginal sample variance"
    );
    if (!(fit.diagonal(j) > 0.0)) {
      Rcpp::stop(
        "The scaled spatial median requires every marginal sample "
        "variance to be strictly positive; variable %llu is degenerate. "
        "No ridge is applied.",
        static_cast<unsigned long long>(j + 1)
      );
    }
  }
  ssmax_normalise_diagonal(fit.diagonal);
  const arma::vec initial_diagonal = fit.diagonal;

  for (int update = 0; update < max_iter; ++update) {
    const SsmaxScore score = ssmax_score(
      data, fit.location, fit.diagonal, zero_tol
    );
    fit.minimum_residual_distance = std::min(
      fit.minimum_residual_distance, score.minimum_radius
    );
    const arma::vec next_location = fit.location +
      arma::sqrt(fit.diagonal) %
        (score.direction_sum.t() / score.inverse_radius_sum);
    arma::vec next_diagonal = fit.diagonal %
      (static_cast<double>(p) * score.direction_square_sum /
       static_cast<double>(n));
    if (!next_location.is_finite() || !next_diagonal.is_finite() ||
        arma::any(next_diagonal <= 0.0)) {
      Rcpp::stop(
        "A scaled spatial median update became non-finite or non-positive; "
        "no ridge, floor, absolute-value repair, or pseudoinverse is "
        "applied."
      );
    }
    ssmax_normalise_diagonal(next_diagonal);
    fit.location_relative_update = arma::abs(
      (next_location - fit.location) / arma::sqrt(initial_diagonal)
    ).max();
    fit.log_diagonal_update = arma::abs(
      arma::log(next_diagonal / fit.diagonal)
    ).max();
    fit.relative_update = std::max(
      fit.location_relative_update, fit.log_diagonal_update
    );
    fit.location = next_location;
    fit.diagonal = next_diagonal;
    fit.iterations = update + 1;
    if (fit.relative_update <= tol) {
      fit.iteration_stable = true;
      break;
    }
  }

  const SsmaxScore final_score = ssmax_score(
    data, fit.location, fit.diagonal, zero_tol
  );
  fit.minimum_residual_distance = std::min(
    fit.minimum_residual_distance, final_score.minimum_radius
  );
  fit.location_score_residual = final_score.location_residual;
  fit.diagonal_score_residual = final_score.diagonal_residual;
  fit.score_residual = final_score.score_residual;
  return fit;
}

void ssmax_store_diagonal(const arma::vec& diagonal,
                          const SsmaxPreparedData& prepared,
                          arma::vec& canonical,
                          arma::vec& canonical_log) {
  canonical.set_size(diagonal.n_elem);
  canonical_log.set_size(diagonal.n_elem);
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

}  // namespace


//' Full-sample scaled spatial median kernel
//'
//' @param x Numeric observation-by-variable matrix.
//' @param target Finite centering vector. Ignored when `center_on_target` is
//'   false.
//' @param center_on_target Whether to express the fit relative to `target`.
//' @param tol Positive estimating-equation tolerance.
//' @param max_iter Positive maximum update count.
//' @param zero_tol Non-negative singular-radius tolerance.
//' @return Internal fit, max statistic ingredients, and diagnostics.
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List cpp_scaled_spatial_median(
    const arma::mat& x,
    const arma::vec& target,
    const bool center_on_target,
    const double tol,
    const int max_iter,
    const double zero_tol) {
  if (!x.is_finite() || !target.is_finite()) {
    Rcpp::stop("`x` and `target` must contain only finite values.");
  }
  if (x.n_rows < 2) {
    Rcpp::stop("The scaled spatial median requires at least two observations.");
  }
  if (x.n_cols < 1 || target.n_elem != x.n_cols) {
    Rcpp::stop("`target` must have one entry per positive matrix column.");
  }
  if (!std::isfinite(tol) || tol <= 0.0 || max_iter < 1 ||
      !std::isfinite(zero_tol) || zero_tol < 0.0) {
    Rcpp::stop("Invalid scaled spatial median iteration controls.");
  }

  const arma::uword n = x.n_rows;
  const arma::uword p = x.n_cols;
  const SsmaxPreparedData prepared = ssmax_prepare_data(
    x, target, center_on_target
  );
  const SsmaxFit fit = ssmax_fit(
    prepared.data, tol, max_iter, zero_tol
  );
  const SsmaxScore final_score = ssmax_score(
    prepared.data, fit.location, fit.diagonal, zero_tol
  );
  const double zeta1_hat = final_score.inverse_radius_sum /
    static_cast<double>(n);
  if (!(zeta1_hat > 0.0) || !std::isfinite(zeta1_hat)) {
    Rcpp::stop(
      "The inverse-radius moment is not finite and strictly positive; no "
      "weight cap is applied."
    );
  }

  arma::vec standardized_location = fit.location /
    arma::sqrt(fit.diagonal);
  arma::vec coordinate_statistic(p, arma::fill::zeros);
  long double maximum_square = 0.0L;
  for (arma::uword j = 0; j < p; ++j) {
    const long double value = static_cast<long double>(
      standardized_location(j)
    );
    maximum_square = std::max(maximum_square, value * value);
  }
  const long double n_ld = static_cast<long double>(n);
  const long double p_ld = static_cast<long double>(p);
  const long double factor_ld = 1.0L - 1.0L / std::sqrt(n_ld);
  const long double zeta_ld = static_cast<long double>(zeta1_hat);
  const long double multiplier = n_ld * p_ld * zeta_ld * zeta_ld *
    factor_ld;
  for (arma::uword j = 0; j < p; ++j) {
    const long double value = static_cast<long double>(
      standardized_location(j)
    );
    coordinate_statistic(j) = ssmax_checked_double(
      multiplier * value * value,
      "coordinatewise max statistic"
    );
  }
  const double t_max = ssmax_checked_double(
    multiplier * maximum_square,
    "spatial-sign max statistic"
  );

  arma::vec diagonal_input;
  arma::vec log_diagonal_input;
  ssmax_store_diagonal(
    fit.diagonal, prepared, diagonal_input, log_diagonal_input
  );
  const arma::vec location_input = ssmax_restore_location(
    fit.location, prepared
  );
  const arma::vec location_relative_origin = ssmax_restore_difference(
    fit.location, prepared
  );

  return Rcpp::List::create(
    Rcpp::Named("location") = location_input,
    Rcpp::Named("location_relative_origin") =
      location_relative_origin,
    Rcpp::Named("location_standardized") = fit.location,
    Rcpp::Named("standardized_location_coordinate") =
      standardized_location,
    Rcpp::Named("origin") = prepared.origin,
    Rcpp::Named("sample_mean") = prepared.sample_mean,
    Rcpp::Named("diagonal_standardized") = fit.diagonal,
    Rcpp::Named("diagonal_input_canonical") = diagonal_input,
    Rcpp::Named("log_diagonal_input_canonical") = log_diagonal_input,
    Rcpp::Named("directions") = final_score.directions,
    Rcpp::Named("radii") = final_score.radii,
    Rcpp::Named("inverse_radii") = final_score.inverse_radii,
    Rcpp::Named("zeta1_hat") = zeta1_hat,
    Rcpp::Named("T_max") = t_max,
    Rcpp::Named("coordinate_statistic") = coordinate_statistic,
    Rcpp::Named("finite_sample_factor") =
      static_cast<double>(factor_ld),
    Rcpp::Named("iterations") = fit.iterations,
    Rcpp::Named("iteration_stable") = fit.iteration_stable,
    Rcpp::Named("relative_update") = fit.relative_update,
    Rcpp::Named("location_relative_update") =
      fit.location_relative_update,
    Rcpp::Named("log_diagonal_update") = fit.log_diagonal_update,
    Rcpp::Named("score_residual") = fit.score_residual,
    Rcpp::Named("location_score_residual") =
      fit.location_score_residual,
    Rcpp::Named("diagonal_score_residual") =
      fit.diagonal_score_residual,
    Rcpp::Named("minimum_residual_distance") =
      fit.minimum_residual_distance,
    Rcpp::Named("center_on_target") = center_on_target,
    Rcpp::Named("column_log_scale") = prepared.log_scale,
    Rcpp::Named("subtraction_overflow_columns") =
      arma::accu(prepared.overflow_fallback),
    Rcpp::Named("n") = static_cast<double>(n),
    Rcpp::Named("p") = static_cast<double>(p),
    Rcpp::Named("tolerance") = tol,
    Rcpp::Named("max_iterations") = max_iter,
    Rcpp::Named("zero_tolerance") = zero_tol
  );
}
