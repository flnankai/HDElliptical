// Chapter 3: high-dimensional Hettmansperger--Randles estimation.
//
// This kernel implements Algorithm 2 of Yan, Feng, and Zhang (2025): the
// location and trace-p shape are updated jointly, and hard banding is applied
// to the standardized-sign SSCM inside every shape update.  The kernel never
// adds a ridge, floors an eigenvalue, uses a pseudoinverse, or perturbs a zero
// residual.  Every matrix square root and inverse is obtained from a checked
// symmetric eigendecomposition.

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp17)]]

#include <RcppArmadillo.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <string>
#include <vector>

namespace {

struct Ch3HdhrSpectral {
  bool valid;
  std::string message;
  arma::vec values;
  arma::mat vectors;
  double minimum;
  double maximum;
  double reciprocal_condition;
};

struct Ch3HdhrMedian {
  bool stable;
  std::string failure;
  arma::vec center;
  int iterations;
  double relative_update;
  double score_residual;
  int coincident_observations;
};

struct Ch3HdhrMap {
  bool valid;
  std::string failure_stage;
  std::string failure;
  arma::vec location_candidate;
  arma::mat shape_candidate;
  arma::mat raw_shape_candidate;
  arma::mat sscm;
  arma::mat banded_sscm;
  arma::mat signs;
  double location_relative_update;
  double shape_relative_update;
  double relative_update;
  double score_l2;
  double score_infinity;
  double minimum_radius;
  double banded_minimum_eigenvalue;
  double banded_reciprocal_condition;
  double candidate_minimum_eigenvalue;
  double candidate_reciprocal_condition;
};

long double ch3_hdhr_safe_norm(const arma::vec& x) {
  long double maximum = 0.0L;
  for (arma::uword j = 0; j < x.n_elem; ++j) {
    const long double value = std::fabs(static_cast<long double>(x(j)));
    if (value > maximum) {
      maximum = value;
    }
  }
  if (maximum == 0.0L) {
    return 0.0L;
  }
  long double sum = 0.0L;
  for (arma::uword j = 0; j < x.n_elem; ++j) {
    const long double ratio = static_cast<long double>(x(j)) / maximum;
    sum += ratio * ratio;
  }
  return maximum * std::sqrt(sum);
}

double ch3_hdhr_checked_double(const long double value,
                               const std::string& quantity) {
  const long double maximum = static_cast<long double>(
    std::numeric_limits<double>::max()
  );
  if (!std::isfinite(value) || value > maximum || value < -maximum) {
    Rcpp::stop(
      "%s is outside the finite double range; no clipping or numerical "
      "repair was applied.", quantity.c_str()
    );
  }
  const double answer = static_cast<double>(value);
  if (value != 0.0L && answer == 0.0) {
    Rcpp::stop(
      "%s underflows the positive double range; no flooring or clipping was "
      "applied.", quantity.c_str()
    );
  }
  return answer;
}

double ch3_hdhr_relative_vec(const arma::vec& candidate,
                             const arma::vec& current) {
  const double denominator = std::max(1.0, arma::norm(current, 2));
  return arma::norm(candidate - current, 2) / denominator;
}

double ch3_hdhr_relative_mat(const arma::mat& candidate,
                             const arma::mat& current) {
  const double denominator = std::max(1.0, arma::norm(current, "fro"));
  return arma::norm(candidate - current, "fro") / denominator;
}

Ch3HdhrSpectral ch3_hdhr_spectral(const arma::mat& matrix,
                                  const std::string& name) {
  Ch3HdhrSpectral result;
  result.valid = false;
  result.minimum = NA_REAL;
  result.maximum = NA_REAL;
  result.reciprocal_condition = NA_REAL;

  if (matrix.n_rows != matrix.n_cols || matrix.n_rows == 0) {
    result.message = name + " must be a non-empty square matrix.";
    return result;
  }
  if (!matrix.is_finite()) {
    result.message = name + " contains a non-finite value.";
    return result;
  }
  const double scale = arma::abs(matrix).max();
  if (!std::isfinite(scale) || scale <= 0.0) {
    result.message = name + " has no finite positive matrix scale.";
    return result;
  }
  const double asymmetry = arma::abs(matrix - matrix.t()).max();
  const long double symmetry_tolerance =
    64.0L * static_cast<long double>(
      std::numeric_limits<double>::epsilon()
    ) * static_cast<long double>(scale);
  if (!std::isfinite(asymmetry) ||
      static_cast<long double>(asymmetry) > symmetry_tolerance) {
    result.message = name + " is not symmetric at floating-point precision.";
    return result;
  }
  const arma::mat symmetric = 0.5 * (matrix + matrix.t());
  bool ok = false;
  try {
    ok = arma::eig_sym(result.values, result.vectors, symmetric);
  } catch (...) {
    ok = false;
  }
  if (!ok || result.values.n_elem == 0 || !result.values.is_finite() ||
      !result.vectors.is_finite()) {
    result.message = "The eigendecomposition of " + name + " failed.";
    return result;
  }
  result.minimum = result.values.min();
  result.maximum = result.values.max();
  if (!std::isfinite(result.minimum) || !std::isfinite(result.maximum) ||
      result.minimum <= 0.0 || result.maximum <= 0.0) {
    result.message = name + " is not positive definite.";
    return result;
  }
  result.reciprocal_condition = result.minimum / result.maximum;
  result.valid = true;
  return result;
}

arma::mat ch3_hdhr_from_eigen(const Ch3HdhrSpectral& spectral,
                              const double power,
                              const std::string& name) {
  arma::vec transformed(spectral.values.n_elem);
  for (arma::uword j = 0; j < spectral.values.n_elem; ++j) {
    const long double value = std::pow(
      static_cast<long double>(spectral.values(j)),
      static_cast<long double>(power)
    );
    transformed(j) = ch3_hdhr_checked_double(value, name);
  }
  arma::mat answer = spectral.vectors * arma::diagmat(transformed) *
    spectral.vectors.t();
  answer = 0.5 * (answer + answer.t());
  if (!answer.is_finite()) {
    Rcpp::stop("%s contains a non-finite value.", name.c_str());
  }
  return answer;
}

arma::vec ch3_hdhr_coordinate_medians(const arma::mat& x) {
  arma::vec answer(x.n_cols);
  std::vector<double> values(x.n_rows);
  for (arma::uword j = 0; j < x.n_cols; ++j) {
    for (arma::uword i = 0; i < x.n_rows; ++i) {
      values[i] = x(i, j);
    }
    std::sort(values.begin(), values.end());
    if (x.n_rows % 2 == 1) {
      answer(j) = values[x.n_rows / 2];
    } else {
      answer(j) = 0.5 * values[x.n_rows / 2 - 1] +
        0.5 * values[x.n_rows / 2];
    }
  }
  return answer;
}

Ch3HdhrMedian ch3_hdhr_spatial_median(const arma::mat& x,
                                      const double tolerance,
                                      const int maximum_iterations,
                                      const double zero_tolerance) {
  Ch3HdhrMedian result;
  result.stable = false;
  result.failure.clear();
  result.center = ch3_hdhr_coordinate_medians(x);
  result.iterations = 0;
  result.relative_update = R_PosInf;
  result.score_residual = R_PosInf;
  result.coincident_observations = 0;

  for (int iteration = 1; iteration <= maximum_iterations; ++iteration) {
    std::vector<long double> radii(x.n_rows);
    std::vector<arma::vec> differences(x.n_rows);
    long double minimum_positive = std::numeric_limits<long double>::infinity();
    arma::vec score(x.n_cols, arma::fill::zeros);
    int coincident = 0;

    for (arma::uword i = 0; i < x.n_rows; ++i) {
      differences[i] = x.row(i).t() - result.center;
      radii[i] = ch3_hdhr_safe_norm(differences[i]);
      if (radii[i] <= static_cast<long double>(zero_tolerance)) {
        ++coincident;
      } else {
        minimum_positive = std::min(minimum_positive, radii[i]);
        for (arma::uword j = 0; j < x.n_cols; ++j) {
          score(j) += static_cast<double>(
            static_cast<long double>(differences[i](j)) / radii[i]
          );
        }
      }
    }

    const long double score_norm_ld = ch3_hdhr_safe_norm(score);
    const double score_norm = ch3_hdhr_checked_double(
      score_norm_ld, "the spatial-median score norm"
    );
    result.score_residual = std::max(0.0, score_norm - coincident) /
      static_cast<double>(x.n_rows);
    result.coincident_observations = coincident;
    result.iterations = iteration;

    if (coincident == static_cast<int>(x.n_rows)) {
      result.stable = true;
      result.relative_update = 0.0;
      return result;
    }
    if (coincident > 0 && score_norm <= static_cast<double>(coincident)) {
      result.stable = true;
      result.relative_update = 0.0;
      return result;
    }
    if (!std::isfinite(minimum_positive) || minimum_positive <= 0.0L) {
      result.failure = "The spatial-median update has no positive radius.";
      return result;
    }

    arma::vec weighted_sum(x.n_cols, arma::fill::zeros);
    long double weight_sum = 0.0L;
    for (arma::uword i = 0; i < x.n_rows; ++i) {
      if (radii[i] > static_cast<long double>(zero_tolerance)) {
        const long double weight = minimum_positive / radii[i];
        weight_sum += weight;
        for (arma::uword j = 0; j < x.n_cols; ++j) {
          weighted_sum(j) += static_cast<double>(
            weight * static_cast<long double>(x(i, j))
          );
        }
      }
    }
    if (!std::isfinite(weight_sum) || weight_sum <= 0.0L) {
      result.failure = "The spatial-median weight sum is not finite positive.";
      return result;
    }
    arma::vec standard_update = weighted_sum /
      ch3_hdhr_checked_double(weight_sum, "the spatial-median weight sum");
    double retained = 0.0;
    if (coincident > 0) {
      retained = std::min(1.0, static_cast<double>(coincident) / score_norm);
    }
    arma::vec candidate = (1.0 - retained) * standard_update +
      retained * result.center;
    if (!candidate.is_finite()) {
      result.failure = "The spatial-median candidate is not finite.";
      return result;
    }
    result.relative_update = ch3_hdhr_relative_vec(candidate, result.center);
    if (!std::isfinite(result.relative_update)) {
      result.failure = "The spatial-median relative update is not finite.";
      return result;
    }
    if (result.relative_update <= tolerance) {
      result.stable = true;
      return result;
    }
    result.center = candidate;
  }

  result.failure = "The sample spatial median did not stabilize within max_iter.";
  return result;
}

arma::mat ch3_hdhr_band(const arma::mat& matrix, const int bandwidth) {
  arma::mat answer = matrix;
  for (arma::uword i = 0; i < matrix.n_rows; ++i) {
    for (arma::uword j = 0; j < matrix.n_cols; ++j) {
      const long long distance = std::llabs(
        static_cast<long long>(i) - static_cast<long long>(j)
      );
      if (distance > bandwidth) {
        answer(i, j) = 0.0;
      }
    }
  }
  return 0.5 * (answer + answer.t());
}

arma::mat ch3_hdhr_normalize_trace(const arma::mat& matrix,
                                   const arma::uword dimension,
                                   const std::string& name) {
  long double trace = 0.0L;
  for (arma::uword j = 0; j < dimension; ++j) {
    trace += static_cast<long double>(matrix(j, j));
  }
  if (!std::isfinite(trace) || trace <= 0.0L) {
    Rcpp::stop("%s has no finite positive trace.", name.c_str());
  }
  const long double factor = static_cast<long double>(dimension) / trace;
  arma::mat answer(matrix.n_rows, matrix.n_cols);
  for (arma::uword i = 0; i < matrix.n_rows; ++i) {
    for (arma::uword j = 0; j < matrix.n_cols; ++j) {
      answer(i, j) = ch3_hdhr_checked_double(
        factor * static_cast<long double>(matrix(i, j)), name
      );
    }
  }
  answer = 0.5 * (answer + answer.t());
  return answer;
}

Ch3HdhrMap ch3_hdhr_fixed_point_map(const arma::mat& x,
                                    const arma::vec& location,
                                    const arma::mat& shape,
                                    const int bandwidth,
                                    const double zero_tolerance) {
  Ch3HdhrMap result;
  result.valid = false;
  result.failure_stage.clear();
  result.failure.clear();
  result.location_relative_update = NA_REAL;
  result.shape_relative_update = NA_REAL;
  result.relative_update = NA_REAL;
  result.score_l2 = NA_REAL;
  result.score_infinity = NA_REAL;
  result.minimum_radius = NA_REAL;
  result.banded_minimum_eigenvalue = NA_REAL;
  result.banded_reciprocal_condition = NA_REAL;
  result.candidate_minimum_eigenvalue = NA_REAL;
  result.candidate_reciprocal_condition = NA_REAL;

  const Ch3HdhrSpectral shape_spectral = ch3_hdhr_spectral(
    shape, "the current HDHR shape"
  );
  if (!shape_spectral.valid) {
    result.failure_stage = "current shape";
    result.failure = shape_spectral.message;
    return result;
  }

  arma::mat square_root;
  arma::mat inverse_square_root;
  try {
    square_root = ch3_hdhr_from_eigen(
      shape_spectral, 0.5, "the current shape square root"
    );
    inverse_square_root = ch3_hdhr_from_eigen(
      shape_spectral, -0.5, "the current inverse shape square root"
    );
  } catch (std::exception& error) {
    result.failure_stage = "matrix square root";
    result.failure = error.what();
    return result;
  }

  const arma::uword n = x.n_rows;
  const arma::uword p = x.n_cols;
  result.signs.zeros(n, p);
  arma::vec score(p, arma::fill::zeros);
  std::vector<long double> radii(n);
  long double inverse_radius_sum = 0.0L;
  long double minimum_radius = std::numeric_limits<long double>::infinity();

  for (arma::uword i = 0; i < n; ++i) {
    const arma::vec difference = x.row(i).t() - location;
    const arma::vec standardized = inverse_square_root * difference;
    radii[i] = ch3_hdhr_safe_norm(standardized);
    if (!std::isfinite(radii[i]) ||
        radii[i] <= static_cast<long double>(zero_tolerance)) {
      result.failure_stage = "standardized residual";
      result.failure =
        "A standardized residual is zero or non-finite; Algorithm 2's "
        "inverse-radius update is undefined and no perturbation was applied.";
      return result;
    }
    minimum_radius = std::min(minimum_radius, radii[i]);
    inverse_radius_sum += 1.0L / radii[i];
    for (arma::uword j = 0; j < p; ++j) {
      const double sign = static_cast<double>(
        static_cast<long double>(standardized(j)) / radii[i]
      );
      result.signs(i, j) = sign;
      score(j) += sign;
    }
  }
  if (!std::isfinite(inverse_radius_sum) || inverse_radius_sum <= 0.0L) {
    result.failure_stage = "inverse-radius denominator";
    result.failure = "The inverse-radius denominator is not finite positive.";
    return result;
  }

  score /= static_cast<double>(n);
  result.score_l2 = arma::norm(score, 2);
  result.score_infinity = arma::abs(score).max();
  result.minimum_radius = ch3_hdhr_checked_double(
    minimum_radius, "the minimum standardized residual radius"
  );

  const arma::vec square_root_score = square_root * score;
  arma::vec location_increment(p);
  const long double mean_inverse_radius = inverse_radius_sum /
    static_cast<long double>(n);
  for (arma::uword j = 0; j < p; ++j) {
    location_increment(j) = ch3_hdhr_checked_double(
      static_cast<long double>(square_root_score(j)) /
        mean_inverse_radius,
      "the HDHR location increment"
    );
  }
  result.location_candidate = location + location_increment;
  if (!result.location_candidate.is_finite()) {
    result.failure_stage = "location update";
    result.failure = "The HDHR location candidate is not finite.";
    return result;
  }

  result.sscm = result.signs.t() * result.signs /
    static_cast<double>(n);
  result.sscm = 0.5 * (result.sscm + result.sscm.t());
  result.banded_sscm = ch3_hdhr_band(result.sscm, bandwidth);

  if (bandwidth == static_cast<int>(p) - 1 && n < p) {
    result.failure_stage = "banded sign SSCM";
    result.failure =
      "With bandwidth p - 1 and n < p, the unbanded sign SSCM is "
      "structurally singular; no ridge or pseudoinverse was applied.";
    return result;
  }
  const Ch3HdhrSpectral banded_spectral = ch3_hdhr_spectral(
    result.banded_sscm, "the banded standardized-sign SSCM"
  );
  result.banded_minimum_eigenvalue = banded_spectral.minimum;
  result.banded_reciprocal_condition = banded_spectral.reciprocal_condition;
  if (!banded_spectral.valid) {
    result.failure_stage = "banded sign SSCM";
    result.failure = banded_spectral.message +
      " No ridge, eigenvalue floor, or taper replacement was applied.";
    return result;
  }

  try {
    const arma::mat raw_congruence = square_root * result.sscm * square_root;
    const arma::mat structured_congruence = square_root *
      result.banded_sscm * square_root;
    result.raw_shape_candidate = ch3_hdhr_normalize_trace(
      raw_congruence, p, "the raw HDHR shape map"
    );
    result.shape_candidate = ch3_hdhr_normalize_trace(
      structured_congruence, p, "the structured HDHR shape map"
    );
  } catch (std::exception& error) {
    result.failure_stage = "shape update";
    result.failure = error.what();
    return result;
  }

  const Ch3HdhrSpectral candidate_spectral = ch3_hdhr_spectral(
    result.shape_candidate, "the structured HDHR shape candidate"
  );
  result.candidate_minimum_eigenvalue = candidate_spectral.minimum;
  result.candidate_reciprocal_condition =
    candidate_spectral.reciprocal_condition;
  if (!candidate_spectral.valid) {
    result.failure_stage = "shape candidate";
    result.failure = candidate_spectral.message +
      " No positive-definite repair was applied.";
    return result;
  }

  result.location_relative_update = ch3_hdhr_relative_vec(
    result.location_candidate, location
  );
  result.shape_relative_update = ch3_hdhr_relative_mat(
    result.shape_candidate, shape
  );
  result.relative_update = std::max(
    result.location_relative_update, result.shape_relative_update
  );
  if (!std::isfinite(result.relative_update)) {
    result.failure_stage = "relative update";
    result.failure = "The joint HDHR relative update is not finite.";
    return result;
  }
  result.valid = true;
  return result;
}

Rcpp::List ch3_hdhr_map_diagnostics(const Ch3HdhrMap& map) {
  return Rcpp::List::create(
    Rcpp::Named("location.relative.update") = map.location_relative_update,
    Rcpp::Named("shape.relative.update") = map.shape_relative_update,
    Rcpp::Named("relative.update") = map.relative_update,
    Rcpp::Named("score.residual.l2") = map.score_l2,
    Rcpp::Named("score.residual.infinity") = map.score_infinity,
    Rcpp::Named("minimum.residual.distance.scaled") = map.minimum_radius,
    Rcpp::Named("banded.sscm.minimum.eigenvalue") =
      map.banded_minimum_eigenvalue,
    Rcpp::Named("banded.sscm.rcond") =
      map.banded_reciprocal_condition,
    Rcpp::Named("shape.map.minimum.eigenvalue") =
      map.candidate_minimum_eigenvalue,
    Rcpp::Named("shape.map.rcond") =
      map.candidate_reciprocal_condition
  );
}

} // namespace


// [[Rcpp::export]]
Rcpp::List cpp_ch3_hdhr_fit(const arma::mat& x,
                            const arma::mat& pilot_precision,
                            const int bandwidth,
                            const double median_tolerance,
                            const int median_maximum_iterations,
                            const double tolerance,
                            const int maximum_iterations,
                            const double zero_tolerance) {
  if (x.n_rows < 2 || x.n_cols < 1 || !x.is_finite()) {
    Rcpp::stop("`x` must be a finite matrix with at least two rows and one column.");
  }
  if (pilot_precision.n_rows != x.n_cols ||
      pilot_precision.n_cols != x.n_cols) {
    Rcpp::stop("`pilot_precision` must be a p by p matrix matching `x`.");
  }
  if (bandwidth < 0 || bandwidth >= static_cast<int>(x.n_cols)) {
    Rcpp::stop("`bandwidth` must be an integer from zero through p - 1.");
  }

  const arma::uword n = x.n_rows;
  const arma::uword p = x.n_cols;
  arma::rowvec anchor = x.row(0);
  long double data_scale = 0.0L;
  for (arma::uword i = 0; i < n; ++i) {
    for (arma::uword j = 0; j < p; ++j) {
      const long double difference = static_cast<long double>(x(i, j)) -
        static_cast<long double>(anchor(j));
      data_scale = std::max(data_scale, std::fabs(difference));
    }
  }
  if (!std::isfinite(data_scale) || data_scale <= 0.0L) {
    return Rcpp::List::create(
      Rcpp::Named("valid") = false,
      Rcpp::Named("failure_stage") = "data scaling",
      Rcpp::Named("failure") =
        "All observations coincide, so neither shape nor inverse-radius "
        "updates are defined.",
      Rcpp::Named("iterations") = 0
    );
  }
  arma::mat scaled_x(n, p);
  for (arma::uword i = 0; i < n; ++i) {
    for (arma::uword j = 0; j < p; ++j) {
      const long double difference = static_cast<long double>(x(i, j)) -
        static_cast<long double>(anchor(j));
      const long double scaled = difference / data_scale;
      scaled_x(i, j) = static_cast<double>(scaled);
      if (difference != 0.0L && scaled_x(i, j) == 0.0) {
        return Rcpp::List::create(
          Rcpp::Named("valid") = false,
          Rcpp::Named("failure_stage") = "data scaling",
          Rcpp::Named("failure") =
            "A nonzero centered data difference underflows double precision "
            "after common scaling; no coordinate-wise rescaling was applied.",
          Rcpp::Named("iterations") = 0
        );
      }
    }
  }
  const double log_data_scale = static_cast<double>(std::log(data_scale));

  const Ch3HdhrSpectral pilot_spectral = ch3_hdhr_spectral(
    pilot_precision, "the pilot precision"
  );
  if (!pilot_spectral.valid) {
    return Rcpp::List::create(
      Rcpp::Named("valid") = false,
      Rcpp::Named("failure_stage") = "pilot precision",
      Rcpp::Named("failure") = pilot_spectral.message +
        " No ridge, eigenvalue floor, or pseudoinverse was applied.",
      Rcpp::Named("iterations") = 0,
      Rcpp::Named("pilot_minimum_eigenvalue") = pilot_spectral.minimum,
      Rcpp::Named("pilot_rcond") = pilot_spectral.reciprocal_condition,
      Rcpp::Named("log_data_scale") = log_data_scale
    );
  }

  arma::vec inverse_weights(p);
  long double inverse_weight_sum = 0.0L;
  for (arma::uword j = 0; j < p; ++j) {
    const long double weight =
      static_cast<long double>(pilot_spectral.minimum) /
      static_cast<long double>(pilot_spectral.values(j));
    inverse_weights(j) = static_cast<double>(weight);
    inverse_weight_sum += weight;
  }
  if (!std::isfinite(inverse_weight_sum) || inverse_weight_sum <= 0.0L) {
    return Rcpp::List::create(
      Rcpp::Named("valid") = false,
      Rcpp::Named("failure_stage") = "pilot inversion",
      Rcpp::Named("failure") =
        "The scale-free inverse of the pilot precision could not be formed.",
      Rcpp::Named("iterations") = 0,
      Rcpp::Named("pilot_minimum_eigenvalue") = pilot_spectral.minimum,
      Rcpp::Named("pilot_rcond") = pilot_spectral.reciprocal_condition,
      Rcpp::Named("log_data_scale") = log_data_scale
    );
  }
  arma::vec initial_shape_values = static_cast<double>(p) * inverse_weights /
    ch3_hdhr_checked_double(
      inverse_weight_sum, "the pilot inverse weight sum"
    );
  arma::mat shape = pilot_spectral.vectors *
    arma::diagmat(initial_shape_values) * pilot_spectral.vectors.t();
  shape = 0.5 * (shape + shape.t());

  Ch3HdhrMedian median = ch3_hdhr_spatial_median(
    scaled_x, median_tolerance, median_maximum_iterations, zero_tolerance
  );
  if (!median.stable) {
    return Rcpp::List::create(
      Rcpp::Named("valid") = false,
      Rcpp::Named("failure_stage") = "initial spatial median",
      Rcpp::Named("failure") = median.failure,
      Rcpp::Named("iterations") = 0,
      Rcpp::Named("location_scaled") = median.center,
      Rcpp::Named("shape_last") = shape,
      Rcpp::Named("pilot_minimum_eigenvalue") = pilot_spectral.minimum,
      Rcpp::Named("pilot_rcond") = pilot_spectral.reciprocal_condition,
      Rcpp::Named("log_data_scale") = log_data_scale,
      Rcpp::Named("median") = Rcpp::List::create(
        Rcpp::Named("stable") = false,
        Rcpp::Named("iterations") = median.iterations,
        Rcpp::Named("relative.update") = median.relative_update,
        Rcpp::Named("score.residual") = median.score_residual,
        Rcpp::Named("coincident.observations") =
          median.coincident_observations
      )
    );
  }

  arma::vec location = median.center;
  Ch3HdhrMap map;
  bool stable = false;
  int iterations = 0;
  for (int iteration = 1; iteration <= maximum_iterations; ++iteration) {
    map = ch3_hdhr_fixed_point_map(
      scaled_x, location, shape, bandwidth, zero_tolerance
    );
    iterations = iteration;
    if (!map.valid) {
      return Rcpp::List::create(
        Rcpp::Named("valid") = false,
        Rcpp::Named("failure_stage") = map.failure_stage,
        Rcpp::Named("failure") = map.failure,
        Rcpp::Named("iterations") = iterations,
        Rcpp::Named("location_scaled") = location,
        Rcpp::Named("shape_last") = shape,
        Rcpp::Named("pilot_minimum_eigenvalue") = pilot_spectral.minimum,
        Rcpp::Named("pilot_rcond") = pilot_spectral.reciprocal_condition,
        Rcpp::Named("log_data_scale") = log_data_scale,
        Rcpp::Named("median") = Rcpp::List::create(
          Rcpp::Named("stable") = true,
          Rcpp::Named("iterations") = median.iterations,
          Rcpp::Named("relative.update") = median.relative_update,
          Rcpp::Named("score.residual") = median.score_residual,
          Rcpp::Named("coincident.observations") =
            median.coincident_observations
        ),
        Rcpp::Named("map") = ch3_hdhr_map_diagnostics(map)
      );
    }
    if (map.relative_update <= tolerance) {
      stable = true;
      break;
    }
    location = map.location_candidate;
    shape = map.shape_candidate;
  }

  if (!stable) {
    return Rcpp::List::create(
      Rcpp::Named("valid") = false,
      Rcpp::Named("failure_stage") = "joint fixed-point iteration",
      Rcpp::Named("failure") =
        "The joint HDHR location/shape update did not stabilize within "
        "max_iter.",
      Rcpp::Named("iterations") = iterations,
      Rcpp::Named("location_scaled") = map.location_candidate,
      Rcpp::Named("shape_last") = map.shape_candidate,
      Rcpp::Named("raw_shape_last") = map.raw_shape_candidate,
      Rcpp::Named("sscm_last") = map.sscm,
      Rcpp::Named("banded_sscm_last") = map.banded_sscm,
      Rcpp::Named("pilot_minimum_eigenvalue") = pilot_spectral.minimum,
      Rcpp::Named("pilot_rcond") = pilot_spectral.reciprocal_condition,
      Rcpp::Named("log_data_scale") = log_data_scale,
      Rcpp::Named("median") = Rcpp::List::create(
        Rcpp::Named("stable") = true,
        Rcpp::Named("iterations") = median.iterations,
        Rcpp::Named("relative.update") = median.relative_update,
        Rcpp::Named("score.residual") = median.score_residual,
        Rcpp::Named("coincident.observations") =
          median.coincident_observations
      ),
      Rcpp::Named("map") = ch3_hdhr_map_diagnostics(map)
    );
  }

  const Ch3HdhrSpectral final_spectral = ch3_hdhr_spectral(
    shape, "the final HDHR shape"
  );
  if (!final_spectral.valid) {
    return Rcpp::List::create(
      Rcpp::Named("valid") = false,
      Rcpp::Named("failure_stage") = "final shape",
      Rcpp::Named("failure") = final_spectral.message,
      Rcpp::Named("iterations") = iterations,
      Rcpp::Named("location_scaled") = location,
      Rcpp::Named("shape_last") = shape,
      Rcpp::Named("pilot_minimum_eigenvalue") = pilot_spectral.minimum,
      Rcpp::Named("pilot_rcond") = pilot_spectral.reciprocal_condition,
      Rcpp::Named("log_data_scale") = log_data_scale
    );
  }

  arma::mat precision;
  try {
    precision = ch3_hdhr_from_eigen(
      final_spectral, -1.0, "the final HDHR shape precision"
    );
  } catch (std::exception& error) {
    return Rcpp::List::create(
      Rcpp::Named("valid") = false,
      Rcpp::Named("failure_stage") = "final shape inversion",
      Rcpp::Named("failure") = error.what(),
      Rcpp::Named("iterations") = iterations,
      Rcpp::Named("location_scaled") = location,
      Rcpp::Named("shape_last") = shape,
      Rcpp::Named("pilot_minimum_eigenvalue") = pilot_spectral.minimum,
      Rcpp::Named("pilot_rcond") = pilot_spectral.reciprocal_condition,
      Rcpp::Named("log_data_scale") = log_data_scale
    );
  }

  arma::vec location_original(p);
  for (arma::uword j = 0; j < p; ++j) {
    const long double value = static_cast<long double>(anchor(j)) +
      data_scale * static_cast<long double>(location(j));
    location_original(j) = ch3_hdhr_checked_double(
      value, "the restored HDHR location"
    );
  }

  arma::rowvec scaled_mean = arma::mean(scaled_x, 0);
  long double qda_trace_scaled = 0.0L;
  for (arma::uword i = 0; i < n; ++i) {
    for (arma::uword j = 0; j < p; ++j) {
      const long double difference =
        static_cast<long double>(scaled_x(i, j)) -
        static_cast<long double>(scaled_mean(j));
      qda_trace_scaled += difference * difference;
    }
  }
  qda_trace_scaled /= static_cast<long double>(n - 1);
  if (!std::isfinite(qda_trace_scaled) || qda_trace_scaled <= 0.0L) {
    return Rcpp::List::create(
      Rcpp::Named("valid") = false,
      Rcpp::Named("failure_stage") = "QDA scatter scale",
      Rcpp::Named("failure") =
        "The primary paper's unbiased covariance-trace add-on is not "
        "finite positive.",
      Rcpp::Named("iterations") = iterations,
      Rcpp::Named("location_scaled") = location,
      Rcpp::Named("shape_last") = shape,
      Rcpp::Named("pilot_minimum_eigenvalue") = pilot_spectral.minimum,
      Rcpp::Named("pilot_rcond") = pilot_spectral.reciprocal_condition,
      Rcpp::Named("log_data_scale") = log_data_scale
    );
  }

  return Rcpp::List::create(
    Rcpp::Named("valid") = true,
    Rcpp::Named("failure_stage") = R_NilValue,
    Rcpp::Named("failure") = R_NilValue,
    Rcpp::Named("location") = location_original,
    Rcpp::Named("location_scaled") = location,
    Rcpp::Named("shape") = shape,
    Rcpp::Named("raw_shape") = map.raw_shape_candidate,
    Rcpp::Named("precision") = precision,
    Rcpp::Named("sscm") = map.sscm,
    Rcpp::Named("banded_sscm") = map.banded_sscm,
    Rcpp::Named("signs") = map.signs,
    Rcpp::Named("iterations") = iterations,
    Rcpp::Named("iteration_stable") = true,
    Rcpp::Named("map") = ch3_hdhr_map_diagnostics(map),
    Rcpp::Named("pilot_minimum_eigenvalue") = pilot_spectral.minimum,
    Rcpp::Named("pilot_rcond") = pilot_spectral.reciprocal_condition,
    Rcpp::Named("shape_minimum_eigenvalue") = final_spectral.minimum,
    Rcpp::Named("shape_rcond") = final_spectral.reciprocal_condition,
    Rcpp::Named("log_data_scale") = log_data_scale,
    Rcpp::Named("qda_trace_scaled") =
      ch3_hdhr_checked_double(qda_trace_scaled, "the scaled QDA trace"),
    Rcpp::Named("median") = Rcpp::List::create(
      Rcpp::Named("stable") = true,
      Rcpp::Named("iterations") = median.iterations,
      Rcpp::Named("relative.update") = median.relative_update,
      Rcpp::Named("score.residual") = median.score_residual,
      Rcpp::Named("coincident.observations") =
        median.coincident_observations
    )
  );
}
