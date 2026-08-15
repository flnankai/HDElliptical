// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp17)]]

#include <RcppArmadillo.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <string>
#include <vector>

namespace ch2genericweighted {

void compensated_add(const long double value, long double& total,
                     long double& correction) {
  const long double updated = total + value;
  if (std::abs(total) >= std::abs(value)) {
    correction += (total - updated) + value;
  } else {
    correction += (value - updated) + total;
  }
  total = updated;
}

double checked_double(const long double value, const std::string& quantity) {
  const long double maximum = static_cast<long double>(
    std::numeric_limits<double>::max()
  );
  if (!std::isfinite(value) || value > maximum || value < -maximum) {
    Rcpp::stop(
      "%s is not representable as a finite double; no rescaling, cap, or "
      "floor is applied.", quantity.c_str()
    );
  }
  return static_cast<double>(value);
}

void validate_matrix(const arma::mat& x) {
  if (x.n_rows < 2U || x.n_cols < 1U || !x.is_finite()) {
    Rcpp::stop(
      "`x` must be a finite numeric matrix with at least two rows and one "
      "column."
    );
  }
}

Rcpp::NumericVector numeric_vector(const arma::vec& x) {
  Rcpp::NumericVector answer(x.n_elem);
  for (arma::uword i = 0; i < x.n_elem; ++i) {
    answer[i] = x[i];
  }
  return answer;
}

Rcpp::NumericVector numeric_vector(const arma::rowvec& x) {
  Rcpp::NumericVector answer(x.n_elem);
  for (arma::uword i = 0; i < x.n_elem; ++i) {
    answer[i] = x[i];
  }
  return answer;
}

}  // namespace ch2genericweighted


//' Strict initial values for the generic weighted HR recursion
//'
//' @param x Observation-by-coordinate finite numeric matrix.
//' @return Column means and unbiased marginal variances.
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List cpp_ch2_generic_weighted_initial(const arma::mat& x) {
  ch2genericweighted::validate_matrix(x);
  arma::vec location(x.n_cols, arma::fill::zeros);
  arma::vec diagonal(x.n_cols, arma::fill::zeros);
  for (arma::uword j = 0; j < x.n_cols; ++j) {
    long double total = 0.0L;
    long double correction = 0.0L;
    for (arma::uword i = 0; i < x.n_rows; ++i) {
      ch2genericweighted::compensated_add(
        static_cast<long double>(x(i, j)), total, correction
      );
    }
    const long double mean = (total + correction) /
      static_cast<long double>(x.n_rows);
    location[j] = ch2genericweighted::checked_double(
      mean, "A marginal initial location"
    );

    long double squares = 0.0L;
    long double square_correction = 0.0L;
    for (arma::uword i = 0; i < x.n_rows; ++i) {
      const long double residual = static_cast<long double>(x(i, j)) - mean;
      ch2genericweighted::compensated_add(
        residual * residual, squares, square_correction
      );
    }
    const long double variance = (squares + square_correction) /
      static_cast<long double>(x.n_rows - 1U);
    if (!(variance > 0.0L)) {
      Rcpp::stop(
        "Every coordinate must have strictly positive marginal sample "
        "variance; no ridge or floor is applied."
      );
    }
    diagonal[j] = ch2genericweighted::checked_double(
      variance, "A marginal initial variance"
    );
  }
  return Rcpp::List::create(
    Rcpp::Named("location") =
      ch2genericweighted::numeric_vector(location),
    Rcpp::Named("diagonal") =
      ch2genericweighted::numeric_vector(diagonal)
  );
}


//' Spatial geometry and unweighted diagonal HR update
//'
//' @param x Observation-by-coordinate finite numeric matrix.
//' @param location Current location vector.
//' @param diagonal Current strictly positive diagonal scale.
//' @param zero_tol Non-negative singular-radius threshold.
//' @return Directions, radii, and the literal unweighted HR scale update.
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List cpp_ch2_generic_weighted_geometry(
    const arma::mat& x, const arma::vec& location,
    const arma::vec& diagonal, const double zero_tol) {
  ch2genericweighted::validate_matrix(x);
  if (location.n_elem != x.n_cols || diagonal.n_elem != x.n_cols ||
      !location.is_finite() || !diagonal.is_finite() ||
      arma::any(diagonal <= 0.0)) {
    Rcpp::stop(
      "The current location and strictly positive diagonal scale are "
      "incompatible with `x`; no ridge or floor is applied."
    );
  }
  if (!std::isfinite(zero_tol) || zero_tol < 0.0) {
    Rcpp::stop("`zero_tol` must be finite and non-negative.");
  }

  arma::mat directions(x.n_rows, x.n_cols, arma::fill::zeros);
  arma::vec radii(x.n_rows, arma::fill::zeros);
  arma::vec square_sum(x.n_cols, arma::fill::zeros);
  double minimum_radius = std::numeric_limits<double>::infinity();

  for (arma::uword i = 0; i < x.n_rows; ++i) {
    arma::rowvec standardized(x.n_cols, arma::fill::zeros);
    double maximum = 0.0;
    for (arma::uword j = 0; j < x.n_cols; ++j) {
      standardized[j] = (x(i, j) - location[j]) /
        std::sqrt(diagonal[j]);
      if (!std::isfinite(standardized[j])) {
        Rcpp::stop(
          "Diagonal standardisation produced a non-finite residual in row "
          "%llu; no rescaling repair is applied.",
          static_cast<unsigned long long>(i + 1U)
        );
      }
      maximum = std::max(maximum, std::abs(standardized[j]));
    }
    if (!(maximum > 0.0)) {
      Rcpp::stop(
        "A standardized residual radius is zero in row %llu; the spatial "
        "direction and weighted update are undefined and no perturbation is "
        "applied.", static_cast<unsigned long long>(i + 1U)
      );
    }
    long double normalized_squares = 0.0L;
    long double normalized_correction = 0.0L;
    for (arma::uword j = 0; j < x.n_cols; ++j) {
      standardized[j] /= maximum;
      const long double value = static_cast<long double>(standardized[j]);
      ch2genericweighted::compensated_add(
        value * value, normalized_squares, normalized_correction
      );
    }
    const long double normalized_norm = std::sqrt(
      normalized_squares + normalized_correction
    );
    if (!(normalized_norm > 0.0L) || !std::isfinite(normalized_norm)) {
      Rcpp::stop("A nonzero standardized residual could not be normalized.");
    }
    const long double radius = static_cast<long double>(maximum) *
      normalized_norm;
    radii[i] = ch2genericweighted::checked_double(
      radius, "A standardized residual radius"
    );
    if (!(radii[i] > zero_tol)) {
      Rcpp::stop(
        "A standardized residual radius is at or below `zero_tol` in row "
        "%llu; no radius floor or perturbation is applied.",
        static_cast<unsigned long long>(i + 1U)
      );
    }
    minimum_radius = std::min(minimum_radius, radii[i]);
    for (arma::uword j = 0; j < x.n_cols; ++j) {
      directions(i, j) = standardized[j] /
        static_cast<double>(normalized_norm);
      square_sum[j] += directions(i, j) * directions(i, j);
    }
  }

  arma::vec diagonal_equation(x.n_cols, arma::fill::zeros);
  arma::vec next_diagonal(x.n_cols, arma::fill::zeros);
  const long double dimension = static_cast<long double>(x.n_cols);
  const long double observations = static_cast<long double>(x.n_rows);
  for (arma::uword j = 0; j < x.n_cols; ++j) {
    const long double equation = dimension *
      static_cast<long double>(square_sum[j]) / observations;
    diagonal_equation[j] = ch2genericweighted::checked_double(
      equation, "An unweighted HR diagonal equation component"
    );
    const long double updated = static_cast<long double>(diagonal[j]) *
      equation;
    if (!(updated > 0.0L)) {
      Rcpp::stop(
        "The unweighted HR diagonal update is non-positive in coordinate "
        "%llu; no ridge or floor is applied.",
        static_cast<unsigned long long>(j + 1U)
      );
    }
    next_diagonal[j] = ch2genericweighted::checked_double(
      updated, "An unweighted HR diagonal update"
    );
  }

  return Rcpp::List::create(
    Rcpp::Named("directions") = directions,
    Rcpp::Named("radii") = ch2genericweighted::numeric_vector(radii),
    Rcpp::Named("next_diagonal") =
      ch2genericweighted::numeric_vector(next_diagonal),
    Rcpp::Named("diagonal_equation") =
      ch2genericweighted::numeric_vector(diagonal_equation),
    Rcpp::Named("minimum_radius") = minimum_radius
  );
}


//' Literal generic weighted HR location update
//'
//' @param location Current location vector.
//' @param diagonal Current strictly positive diagonal scale.
//' @param directions Current spatial directions.
//' @param radii Current strictly positive radii.
//' @param weights Evaluated finite radial weights.
//' @return The weighted numerator, denominator, and next location.
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List cpp_ch2_generic_weighted_step(
    const arma::vec& location, const arma::vec& diagonal,
    const arma::mat& directions, const arma::vec& radii,
    const arma::vec& weights) {
  if (directions.n_rows < 2U || directions.n_cols < 1U ||
      directions.n_cols != location.n_elem ||
      diagonal.n_elem != location.n_elem ||
      radii.n_elem != directions.n_rows ||
      weights.n_elem != directions.n_rows ||
      !location.is_finite() || !diagonal.is_finite() ||
      !directions.is_finite() || !radii.is_finite() ||
      !weights.is_finite() || arma::any(diagonal <= 0.0) ||
      arma::any(radii <= 0.0)) {
    Rcpp::stop("The weighted HR update inputs are incompatible or invalid.");
  }

  long double denominator_total = 0.0L;
  long double denominator_correction = 0.0L;
  long double absolute_total = 0.0L;
  long double absolute_correction = 0.0L;
  std::vector<long double> numerator(location.n_elem, 0.0L);
  std::vector<long double> numerator_correction(location.n_elem, 0.0L);
  for (arma::uword i = 0; i < directions.n_rows; ++i) {
    const long double weight = static_cast<long double>(weights[i]);
    ch2genericweighted::compensated_add(
      weight / static_cast<long double>(radii[i]),
      denominator_total, denominator_correction
    );
    ch2genericweighted::compensated_add(
      std::abs(weight), absolute_total, absolute_correction
    );
    for (arma::uword j = 0; j < directions.n_cols; ++j) {
      ch2genericweighted::compensated_add(
        weight * static_cast<long double>(directions(i, j)),
        numerator[j], numerator_correction[j]
      );
    }
  }
  const long double denominator = denominator_total +
    denominator_correction;
  if (!(denominator > 0.0L) || !std::isfinite(denominator)) {
    Rcpp::stop(
      "The weighted location denominator sum K(r_i)/r_i is not strictly "
      "positive and finite; no sign flip, absolute value, or floor is applied."
    );
  }
  const double denominator_double = ch2genericweighted::checked_double(
    denominator, "The weighted location denominator"
  );
  const long double absolute_weight_sum = absolute_total +
    absolute_correction;
  if (!(absolute_weight_sum > 0.0L) ||
      !std::isfinite(absolute_weight_sum)) {
    Rcpp::stop("The absolute radial-weight mass is invalid.");
  }

  arma::rowvec numerator_output(location.n_elem, arma::fill::zeros);
  arma::rowvec standardized_step(location.n_elem, arma::fill::zeros);
  arma::rowvec next_location(location.n_elem, arma::fill::zeros);
  double normalized_residual = 0.0;
  for (arma::uword j = 0; j < location.n_elem; ++j) {
    const long double current_numerator = numerator[j] +
      numerator_correction[j];
    numerator_output[j] = ch2genericweighted::checked_double(
      current_numerator, "A weighted location numerator component"
    );
    const long double step = current_numerator / denominator;
    standardized_step[j] = ch2genericweighted::checked_double(
      step, "A standardized weighted location step"
    );
    const long double updated = static_cast<long double>(location[j]) +
      std::sqrt(static_cast<long double>(diagonal[j])) * step;
    next_location[j] = ch2genericweighted::checked_double(
      updated, "A weighted HR location iterate"
    );
    normalized_residual = std::max(
      normalized_residual,
      static_cast<double>(std::abs(current_numerator) /
        absolute_weight_sum)
    );
  }

  return Rcpp::List::create(
    Rcpp::Named("next_location") =
      ch2genericweighted::numeric_vector(next_location),
    Rcpp::Named("numerator") =
      ch2genericweighted::numeric_vector(numerator_output),
    Rcpp::Named("denominator") = denominator_double,
    Rcpp::Named("standardized_step") =
      ch2genericweighted::numeric_vector(standardized_step),
    Rcpp::Named("absolute_weight_sum") =
      ch2genericweighted::checked_double(
        absolute_weight_sum, "The absolute radial-weight mass"
      ),
    Rcpp::Named("normalized_location_residual") = normalized_residual
  );
}


//' Generic oracle weighted-sign quadratic U-statistic
//'
//' @param directions Observation-by-coordinate unit spatial directions.
//' @param weights Evaluated finite radial weights.
//' @return The literal pairwise score and empirical weight moment.
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List cpp_ch2_generic_weighted_quadratic(
    const arma::mat& directions, const arma::vec& weights) {
  if (directions.n_rows < 2U || directions.n_cols < 1U ||
      weights.n_elem != directions.n_rows ||
      !directions.is_finite() || !weights.is_finite()) {
    Rcpp::stop("The weighted-sign quadratic inputs are incompatible or invalid.");
  }
  long double pair_total = 0.0L;
  long double pair_correction = 0.0L;
  for (arma::uword first = 0; first + 1U < directions.n_rows; ++first) {
    for (arma::uword second = first + 1U;
         second < directions.n_rows; ++second) {
      long double inner = 0.0L;
      long double inner_correction = 0.0L;
      for (arma::uword j = 0; j < directions.n_cols; ++j) {
        ch2genericweighted::compensated_add(
          static_cast<long double>(directions(first, j)) *
            static_cast<long double>(directions(second, j)),
          inner, inner_correction
        );
      }
      const long double term = static_cast<long double>(weights[first]) *
        static_cast<long double>(weights[second]) *
        (inner + inner_correction);
      ch2genericweighted::compensated_add(
        term, pair_total, pair_correction
      );
    }
  }
  const long double n = static_cast<long double>(directions.n_rows);
  const long double score = 2.0L * (pair_total + pair_correction) /
    (n * (n - 1.0L));

  long double weight_squares = 0.0L;
  long double weight_square_correction = 0.0L;
  std::vector<long double> score_mean_total(
    directions.n_cols, 0.0L
  );
  std::vector<long double> score_mean_correction(
    directions.n_cols, 0.0L
  );
  for (arma::uword i = 0; i < directions.n_rows; ++i) {
    const long double weight = static_cast<long double>(weights[i]);
    ch2genericweighted::compensated_add(
      weight * weight, weight_squares, weight_square_correction
    );
    for (arma::uword j = 0; j < directions.n_cols; ++j) {
      ch2genericweighted::compensated_add(
        weight * static_cast<long double>(directions(i, j)),
        score_mean_total[j], score_mean_correction[j]
      );
    }
  }
  const long double empirical_nu2 = (weight_squares +
    weight_square_correction) / n;
  if (!(empirical_nu2 > 0.0L)) {
    Rcpp::stop(
      "The empirical second radial-weight moment is zero; the weighted-sign "
      "score is degenerate and no replacement weight is used."
    );
  }

  Rcpp::NumericVector score_mean(directions.n_cols);
  for (arma::uword j = 0; j < directions.n_cols; ++j) {
    score_mean[j] = ch2genericweighted::checked_double(
      (score_mean_total[j] + score_mean_correction[j]) / n,
      "A coordinate of the weighted-sign score mean"
    );
  }

  return Rcpp::List::create(
    Rcpp::Named("score") = ch2genericweighted::checked_double(
      score, "The weighted-sign quadratic score"
    ),
    Rcpp::Named("empirical_nu2") = ch2genericweighted::checked_double(
      empirical_nu2, "The empirical second radial-weight moment"
    ),
    Rcpp::Named("weighted_score_mean") = score_mean
  );
}


