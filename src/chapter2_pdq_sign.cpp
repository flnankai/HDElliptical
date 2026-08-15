#include <RcppArmadillo.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <string>
#include <vector>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp17)]]

namespace {

void pdq_require_finite(const arma::mat& x, const char* name) {
  if (!x.is_finite()) {
    Rcpp::stop("`%s` must contain only finite values.", name);
  }
}

double pdq_stable_norm(const arma::rowvec& x) {
  const double scale = arma::abs(x).max();
  if (scale == 0.0) {
    return 0.0;
  }
  const double scaled_norm = arma::norm(x / scale, 2);
  if (!std::isfinite(scaled_norm) || scaled_norm <= 0.0 ||
      scale > std::numeric_limits<double>::max() / scaled_norm) {
    return R_PosInf;
  }
  return scale * scaled_norm;
}

arma::vec pdq_row_norms(const arma::mat& x) {
  arma::vec answer(x.n_rows, arma::fill::zeros);
  for (arma::uword i = 0; i < x.n_rows; ++i) {
    answer(i) = pdq_stable_norm(x.row(i));
  }
  return answer;
}

struct PdqSigns {
  arma::mat signs;
  arma::vec radii;
  arma::vec inverse_radii;
  arma::uword zero_count;
  double minimum_positive_radius;
};

PdqSigns pdq_sign_rows(const arma::mat& x) {
  PdqSigns answer;
  answer.signs.zeros(x.n_rows, x.n_cols);
  answer.radii.zeros(x.n_rows);
  answer.inverse_radii.zeros(x.n_rows);
  answer.zero_count = 0U;
  answer.minimum_positive_radius = R_PosInf;

  for (arma::uword i = 0; i < x.n_rows; ++i) {
    const arma::rowvec row = x.row(i);
    const double scale = arma::abs(row).max();
    if (scale == 0.0) {
      ++answer.zero_count;
      continue;
    }
    const arma::rowvec scaled = row / scale;
    const double scaled_norm = arma::norm(scaled, 2);
    if (!std::isfinite(scaled_norm) || scaled_norm <= 0.0) {
      Rcpp::stop("A PDQ spatial direction could not be normalised.");
    }
    answer.signs.row(i) = scaled / scaled_norm;
    answer.inverse_radii(i) = (1.0 / scale) / scaled_norm;
    if (scale <= std::numeric_limits<double>::max() / scaled_norm) {
      answer.radii(i) = scale * scaled_norm;
      answer.minimum_positive_radius = std::min(
        answer.minimum_positive_radius, answer.radii(i)
      );
    } else {
      answer.radii(i) = R_PosInf;
    }
  }
  return answer;
}

double pdq_spatial_objective(const arma::mat& x,
                             const arma::rowvec& location) {
  const arma::vec distances = pdq_row_norms(x.each_row() - location);
  if (!distances.is_finite()) {
    return R_PosInf;
  }
  const double maximum = distances.max();
  if (maximum == 0.0) {
    return 0.0;
  }
  return maximum * arma::mean(distances / maximum);
}

double pdq_observation_subgradient(const arma::mat& x,
                                   const arma::uword candidate) {
  const PdqSigns directions = pdq_sign_rows(
    x.each_row() - x.row(candidate)
  );
  const arma::rowvec score = arma::sum(directions.signs, 0);
  const double score_norm = pdq_stable_norm(score);
  if (!std::isfinite(score_norm)) {
    return R_PosInf;
  }
  return std::max(
    0.0, score_norm - static_cast<double>(directions.zero_count)
  ) / static_cast<double>(x.n_rows);
}

struct PdqMedian {
  arma::rowvec location;
  int iterations;
  bool converged;
  double relative_change;
  double equation_residual;
  double objective;
  arma::uword zero_count;
};

PdqMedian pdq_spatial_median(const arma::mat& x, const double tol,
                             const int max_iter) {
  PdqMedian answer;
  answer.location = arma::mean(x, 0);
  answer.iterations = 0;
  answer.converged = false;
  answer.relative_change = R_PosInf;
  answer.equation_residual = R_PosInf;
  answer.objective = R_PosInf;
  answer.zero_count = 0U;

  if (!answer.location.is_finite()) {
    Rcpp::stop("The PDQ spatial-median initial value is not finite.");
  }
  const arma::vec initial_distances = pdq_row_norms(
    x.each_row() - answer.location
  );
  if (!initial_distances.is_finite()) {
    Rcpp::stop("The PDQ spatial-median initial distances are not finite.");
  }
  double convergence_scale = initial_distances.max();
  if (!std::isfinite(convergence_scale) || convergence_scale <= 0.0) {
    convergence_scale = 1.0;
  }

  for (int iteration = 1; iteration <= max_iter; ++iteration) {
    answer.iterations = iteration;
    const arma::mat residuals = x.each_row() - answer.location;
    if (!residuals.is_finite()) {
      break;
    }
    const PdqSigns directions = pdq_sign_rows(residuals);
    const arma::rowvec residual_sum = arma::sum(directions.signs, 0);
    const double residual_norm = pdq_stable_norm(residual_sum);
    const double current_residual = std::isfinite(residual_norm) ?
      std::max(
        0.0,
        residual_norm - static_cast<double>(directions.zero_count)
      ) / static_cast<double>(x.n_rows) : R_PosInf;
    if (current_residual <= tol) {
      answer.converged = true;
      answer.relative_change = 0.0;
      break;
    }

    double minimum_distance = R_PosInf;
    arma::uword nearest = 0U;
    for (arma::uword i = 0; i < x.n_rows; ++i) {
      const double distance = directions.radii(i);
      if (distance > 0.0 && distance < minimum_distance) {
        minimum_distance = distance;
        nearest = i;
      }
    }
    if (!std::isfinite(minimum_distance) || minimum_distance <= 0.0) {
      break;
    }

    const double nearest_residual = pdq_observation_subgradient(x, nearest);
    if (nearest_residual <= tol) {
      const arma::rowvec next = x.row(nearest);
      answer.relative_change = pdq_stable_norm(next - answer.location) /
        convergence_scale;
      answer.location = next;
      answer.converged = true;
      break;
    }

    double scaled_denominator = 0.0;
    for (arma::uword i = 0; i < x.n_rows; ++i) {
      if (directions.radii(i) > 0.0) {
        scaled_denominator += minimum_distance / directions.radii(i);
      }
    }
    if (!std::isfinite(scaled_denominator) || scaled_denominator <= 0.0) {
      break;
    }
    const arma::rowvec weiszfeld = answer.location +
      residual_sum * (minimum_distance / scaled_denominator);
    double move_fraction = 1.0;
    if (directions.zero_count > 0U && residual_norm > 0.0) {
      move_fraction = std::max(
        0.0,
        1.0 - static_cast<double>(directions.zero_count) / residual_norm
      );
    }
    const arma::rowvec next = move_fraction * weiszfeld +
      (1.0 - move_fraction) * answer.location;
    if (!next.is_finite()) {
      break;
    }
    answer.relative_change = pdq_stable_norm(next - answer.location) /
      convergence_scale;
    answer.location = next;
  }

  const PdqSigns final_directions = pdq_sign_rows(
    x.each_row() - answer.location
  );
  const arma::rowvec final_score = arma::sum(final_directions.signs, 0);
  const double final_norm = pdq_stable_norm(final_score);
  answer.equation_residual = std::isfinite(final_norm) ?
    std::max(
      0.0,
      final_norm - static_cast<double>(final_directions.zero_count)
    ) / static_cast<double>(x.n_rows) : R_PosInf;
  answer.converged = std::isfinite(answer.equation_residual) &&
    answer.equation_residual <= tol;
  answer.objective = pdq_spatial_objective(x, answer.location);
  answer.zero_count = final_directions.zero_count;
  return answer;
}

struct PdqPreconditioner {
  arma::mat x;
  arma::mat y;
  arma::rowvec magnitude;
  arma::rowvec midpoint_scaled;
  arma::rowvec range_scaled;
  arma::rowvec log_unit;
};

PdqPreconditioner pdq_precondition(const arma::mat& x,
                                   const arma::mat& y) {
  PdqPreconditioner answer;
  answer.x.set_size(x.n_rows, x.n_cols);
  answer.y.set_size(y.n_rows, y.n_cols);
  answer.magnitude.set_size(x.n_cols);
  answer.midpoint_scaled.set_size(x.n_cols);
  answer.range_scaled.set_size(x.n_cols);
  answer.log_unit.set_size(x.n_cols);

  for (arma::uword j = 0; j < x.n_cols; ++j) {
    const double minimum = std::min(x.col(j).min(), y.col(j).min());
    const double maximum = std::max(x.col(j).max(), y.col(j).max());
    const double magnitude = std::max(std::abs(minimum), std::abs(maximum));
    if (!std::isfinite(magnitude) || magnitude <= 0.0) {
      Rcpp::stop(
        "PDQ requires pooled variation in variable %llu; no scale floor "
        "is applied.",
        static_cast<unsigned long long>(j + 1U)
      );
    }
    const double minimum_scaled = minimum / magnitude;
    const double maximum_scaled = maximum / magnitude;
    const double range_scaled = maximum_scaled - minimum_scaled;
    if (!std::isfinite(range_scaled) || range_scaled <= 0.0) {
      Rcpp::stop(
        "PDQ requires pooled variation in variable %llu; no scale floor "
        "is applied.",
        static_cast<unsigned long long>(j + 1U)
      );
    }
    const double midpoint_scaled = minimum_scaled + 0.5 * range_scaled;
    answer.x.col(j) = (
      x.col(j) / magnitude - midpoint_scaled
    ) / range_scaled;
    answer.y.col(j) = (
      y.col(j) / magnitude - midpoint_scaled
    ) / range_scaled;
    answer.magnitude(j) = magnitude;
    answer.midpoint_scaled(j) = midpoint_scaled;
    answer.range_scaled(j) = range_scaled;
    answer.log_unit(j) = std::log(magnitude) + std::log(range_scaled);
  }
  if (!answer.x.is_finite() || !answer.y.is_finite()) {
    Rcpp::stop("The algebraically neutral PDQ preconditioning failed.");
  }
  return answer;
}

arma::uword pdq_quantile_rank(const arma::uword pair_count,
                              const double probability) {
  arma::uword rank = static_cast<arma::uword>(
    std::floor(probability * static_cast<double>(pair_count))
  );
  if (rank < 1U) {
    rank = 1U;
  }
  if (rank > pair_count) {
    rank = pair_count;
  }
  while (rank < pair_count &&
         static_cast<double>(rank) / static_cast<double>(pair_count) <
           probability) {
    ++rank;
  }
  while (rank > 1U &&
         static_cast<double>(rank - 1U) /
           static_cast<double>(pair_count) >= probability) {
    --rank;
  }
  return rank;
}

struct PdqScale {
  arma::rowvec quantile_preconditioned;
  arma::rowvec quantile_working;
  arma::rowvec log_quantile_input;
  arma::rowvec quantile_input;
  arma::rowvec diagonal_input;
  arma::rowvec diagonal_input_canonical;
  arma::rowvec diagonal_working;
  arma::uword pair_count;
  arma::uword rank;
  arma::uword nonrepresentable_quantiles;
  arma::uword nonrepresentable_diagonals;
};

double pdq_exp_display(const double logarithm, bool& representable) {
  const double log_max = std::log(std::numeric_limits<double>::max());
  const double log_min = std::log(std::numeric_limits<double>::denorm_min());
  if (logarithm > log_max) {
    representable = false;
    return R_PosInf;
  }
  if (logarithm < log_min) {
    representable = false;
    return 0.0;
  }
  const double answer = std::exp(logarithm);
  representable = std::isfinite(answer) && answer > 0.0;
  return answer;
}

PdqScale pdq_coordinate_scale(const arma::mat& sample,
                              const arma::rowvec& log_unit,
                              const double probability,
                              const int group) {
  const arma::uword n = sample.n_rows;
  const std::uint64_t count64 =
    static_cast<std::uint64_t>(n) * static_cast<std::uint64_t>(n - 1U) / 2U;
  if (count64 > static_cast<std::uint64_t>(
        std::numeric_limits<arma::uword>::max())) {
    Rcpp::stop("The PDQ pair count exceeds the addressable vector size.");
  }
  PdqScale answer;
  answer.pair_count = static_cast<arma::uword>(count64);
  answer.rank = pdq_quantile_rank(answer.pair_count, probability);
  answer.quantile_preconditioned.set_size(sample.n_cols);
  answer.nonrepresentable_quantiles = 0U;
  answer.nonrepresentable_diagonals = 0U;

  std::vector<double> differences(answer.pair_count);
  for (arma::uword j = 0; j < sample.n_cols; ++j) {
    arma::uword position = 0U;
    for (arma::uword i = 0; i + 1U < n; ++i) {
      for (arma::uword ell = i + 1U; ell < n; ++ell) {
        differences[position++] = std::abs(sample(i, j) - sample(ell, j));
      }
    }
    std::sort(differences.begin(), differences.end());
    const double quantile = differences[answer.rank - 1U];
    if (!std::isfinite(quantile) || quantile <= 0.0) {
      Rcpp::stop(
        "The PDQ U-quantile is not finite and strictly positive for group "
        "%d, variable %llu; no ridge or floor is applied.",
        group, static_cast<unsigned long long>(j + 1U)
      );
    }
    answer.quantile_preconditioned(j) = quantile;
  }

  const double maximum_quantile = answer.quantile_preconditioned.max();
  answer.quantile_working = answer.quantile_preconditioned /
    maximum_quantile;
  if (!answer.quantile_working.is_finite() ||
      arma::any(answer.quantile_working <= 0.0)) {
    Rcpp::stop(
      "PDQ within-group quantile normalisation underflowed; no scale "
      "floor is applied."
    );
  }
  answer.diagonal_working = arma::square(answer.quantile_working);
  answer.log_quantile_input = log_unit +
    arma::log(answer.quantile_preconditioned);
  answer.quantile_input.set_size(sample.n_cols);
  answer.diagonal_input.set_size(sample.n_cols);
  answer.diagonal_input_canonical.set_size(sample.n_cols);
  const double maximum_log_quantile = answer.log_quantile_input.max();

  for (arma::uword j = 0; j < sample.n_cols; ++j) {
    bool q_ok = true;
    bool d_ok = true;
    answer.quantile_input(j) = pdq_exp_display(
      answer.log_quantile_input(j), q_ok
    );
    answer.diagonal_input(j) = pdq_exp_display(
      2.0 * answer.log_quantile_input(j), d_ok
    );
    answer.diagonal_input_canonical(j) = std::exp(
      2.0 * (answer.log_quantile_input(j) - maximum_log_quantile)
    );
    if (!q_ok) {
      ++answer.nonrepresentable_quantiles;
    }
    if (!d_ok) {
      ++answer.nonrepresentable_diagonals;
    }
  }
  return answer;
}

struct PdqGroupFit {
  PdqScale scale;
  PdqMedian median;
  arma::rowvec anchor;
  arma::mat standardized_data;
  arma::rowvec location_preconditioned;
  arma::rowvec location_input;
  PdqSigns fitted;
  arma::mat omega;
  arma::mat G;
  arma::mat G_inverse;
  arma::vec G_eigenvalues;
  double G_condition;
};

PdqGroupFit pdq_fit_group(const arma::mat& sample,
                          const arma::rowvec& log_unit,
                          const arma::rowvec& magnitude,
                          const arma::rowvec& midpoint_scaled,
                          const arma::rowvec& range_scaled,
                          const double probability, const double tol,
                          const int max_iter, const int group) {
  PdqGroupFit answer;
  answer.scale = pdq_coordinate_scale(
    sample, log_unit, probability, group
  );
  answer.anchor = sample.row(0);
  answer.standardized_data = sample.each_row() - answer.anchor;
  answer.standardized_data.each_row() /= answer.scale.quantile_working;
  if (!answer.standardized_data.is_finite()) {
    Rcpp::stop(
      "PDQ standardisation is non-finite for group %d; no numerical "
      "repair is applied.", group
    );
  }
  answer.median = pdq_spatial_median(
    answer.standardized_data, tol, max_iter
  );
  answer.location_preconditioned = answer.anchor +
    answer.scale.quantile_working % answer.median.location;
  answer.location_input = magnitude % (
    midpoint_scaled + range_scaled % answer.location_preconditioned
  );
  if (!answer.location_preconditioned.is_finite() ||
      !answer.location_input.is_finite()) {
    Rcpp::stop("The PDQ spatial median is non-finite in group %d.", group);
  }

  const arma::mat fitted_residuals =
    (sample.each_row() - answer.location_preconditioned);
  arma::mat standardized_residuals = fitted_residuals;
  standardized_residuals.each_row() /= answer.scale.quantile_working;
  answer.fitted = pdq_sign_rows(standardized_residuals);
  if (answer.fitted.zero_count > 0U) {
    Rcpp::stop(
      "The PDQ fitted spatial median in group %d coincides with %llu "
      "observation(s). The inverse-radius definition of G_hat is then "
      "undefined; no deletion, ridge, or perturbation is applied.",
      group,
      static_cast<unsigned long long>(answer.fitted.zero_count)
    );
  }
  if (!answer.fitted.inverse_radii.is_finite() ||
      arma::any(answer.fitted.inverse_radii <= 0.0)) {
    Rcpp::stop(
      "The PDQ inverse fitted radii are not finite and strictly positive "
      "in group %d; no repair is applied.", group
    );
  }

  const double n = static_cast<double>(sample.n_rows);
  answer.omega = answer.fitted.signs.t() * answer.fitted.signs / n;
  answer.G.zeros(sample.n_cols, sample.n_cols);
  const arma::mat identity = arma::eye(sample.n_cols, sample.n_cols);
  for (arma::uword i = 0; i < sample.n_rows; ++i) {
    const arma::rowvec direction = answer.fitted.signs.row(i);
    answer.G += answer.fitted.inverse_radii(i) *
      (identity - direction.t() * direction);
  }
  answer.G /= n;
  answer.G = 0.5 * (answer.G + answer.G.t());
  if (!answer.G.is_finite() ||
      !arma::eig_sym(answer.G_eigenvalues, answer.G) ||
      answer.G_eigenvalues.min() <= 0.0 ||
      !arma::inv_sympd(answer.G_inverse, answer.G)) {
    Rcpp::stop(
      "The empirical PDQ G_hat matrix is not strictly positive definite "
      "in group %d; no ridge or generalized inverse is applied.", group
    );
  }
  answer.G_condition = answer.G_eigenvalues.max() /
    answer.G_eigenvalues.min();
  if (!std::isfinite(answer.G_condition)) {
    Rcpp::stop("The PDQ G_hat condition number is not finite in group %d.",
               group);
  }
  return answer;
}

arma::mat pdq_cross_signs(const arma::mat& sample,
                          const arma::rowvec& other_location,
                          const arma::rowvec& own_quantile,
                          arma::uword& zero_count) {
  arma::mat residuals = sample.each_row() - other_location;
  residuals.each_row() /= own_quantile;
  if (!residuals.is_finite()) {
    Rcpp::stop(
      "A cross-centred PDQ residual is non-finite; no clipping is applied."
    );
  }
  const PdqSigns answer = pdq_sign_rows(residuals);
  zero_count = answer.zero_count;
  return answer.signs;
}

std::uint32_t pdq_mix32(std::uint32_t value) {
  value ^= value >> 16U;
  value *= UINT32_C(0x7feb352d);
  value ^= value >> 15U;
  value *= UINT32_C(0x846ca68b);
  value ^= value >> 16U;
  return value;
}

std::uint32_t pdq_counter_word(const std::uint64_t counter,
                               const std::uint32_t seed) {
  const std::uint32_t low = static_cast<std::uint32_t>(counter);
  const std::uint32_t high = static_cast<std::uint32_t>(counter >> 32U);
  const std::uint32_t high_key = pdq_mix32(
    high + UINT32_C(0x9e3779b9)
  );
  const std::uint32_t keyed = (seed ^ high_key) +
    UINT32_C(0x9e3779b9) * (low + UINT32_C(1));
  return pdq_mix32(keyed);
}

double pdq_quadratic(const arma::rowvec& mean1,
                     const arma::rowvec& mean2,
                     const arma::mat& K1, const arma::mat& K2,
                     const arma::mat& K3) {
  return arma::as_scalar(mean1 * K1 * mean1.t()) +
    arma::as_scalar(mean2 * K2 * mean2.t()) -
    arma::as_scalar(mean1 * K3 * mean2.t());
}

}  // namespace


//' Feng--Wang PDQ spatial-sign kernel
//'
//' @param x,y Numeric matrices with observations in rows.
//' @param quantile_prob Pairwise-difference U-quantile probability.
//' @param level Test level used for the primary empirical critical value.
//' @param B Number of Rademacher draws.
//' @param seed Integer-valued counter-generator seed represented as a double.
//' @param keep_bootstrap Whether to retain bootstrap values and multipliers.
//' @param tolerance Spatial-median estimating-equation tolerance.
//' @param max_iterations Maximum spatial-median iterations.
//' @return A list of the statistic, nuisance fits, bootstrap calibration, and
//'   numerical diagnostics.
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List cpp_feng_wang_pdq_two_sample(
    const arma::mat& x, const arma::mat& y,
    const double quantile_prob, const double level, const int B,
    const double seed, const bool keep_bootstrap,
    const double tolerance, const int max_iterations) {
  pdq_require_finite(x, "x");
  pdq_require_finite(y, "y");
  if (x.n_cols != y.n_cols || x.n_cols < 2U) {
    Rcpp::stop("PDQ requires both samples to have the same p >= 2 variables.");
  }
  if (x.n_rows < 3U || y.n_rows < 3U) {
    Rcpp::stop("PDQ requires at least three observations in each group.");
  }
  if (!std::isfinite(quantile_prob) || quantile_prob <= 0.0 ||
      quantile_prob >= 1.0) {
    Rcpp::stop("`quantile_prob` must be strictly between zero and one.");
  }
  if (!std::isfinite(level) || level <= 0.0 || level >= 1.0) {
    Rcpp::stop("`level` must be strictly between zero and one.");
  }
  if (B < 1) {
    Rcpp::stop("`B` must be a positive integer.");
  }
  if (!std::isfinite(seed) || seed < 0.0 || seed > 4294967295.0 ||
      seed != std::floor(seed)) {
    Rcpp::stop("`seed` must be an integer in [0, 2^32 - 1].");
  }
  if (!std::isfinite(tolerance) || tolerance <= 0.0 ||
      max_iterations < 1) {
    Rcpp::stop("Invalid PDQ spatial-median iteration controls.");
  }

  const arma::uword n1 = x.n_rows;
  const arma::uword n2 = y.n_rows;
  const arma::uword p = x.n_cols;
  const PdqPreconditioner transformed = pdq_precondition(x, y);
  const PdqGroupFit fit1 = pdq_fit_group(
    transformed.x, transformed.log_unit, transformed.magnitude,
    transformed.midpoint_scaled, transformed.range_scaled,
    quantile_prob, tolerance, max_iterations, 1
  );
  const PdqGroupFit fit2 = pdq_fit_group(
    transformed.y, transformed.log_unit, transformed.magnitude,
    transformed.midpoint_scaled, transformed.range_scaled,
    quantile_prob, tolerance, max_iterations, 2
  );

  const arma::rowvec A12_diagonal =
    fit2.scale.quantile_working / fit1.scale.quantile_working;
  const arma::rowvec A21_diagonal =
    fit1.scale.quantile_working / fit2.scale.quantile_working;
  arma::rowvec A12_input_diagonal(p);
  arma::rowvec A21_input_diagonal(p);
  for (arma::uword j = 0; j < p; ++j) {
    bool A12_ok = true;
    bool A21_ok = true;
    A12_input_diagonal(j) = pdq_exp_display(
      fit2.scale.log_quantile_input(j) -
        fit1.scale.log_quantile_input(j),
      A12_ok
    );
    A21_input_diagonal(j) = pdq_exp_display(
      fit1.scale.log_quantile_input(j) -
        fit2.scale.log_quantile_input(j),
      A21_ok
    );
  }
  if (!A12_diagonal.is_finite() || !A21_diagonal.is_finite() ||
      arma::any(A12_diagonal <= 0.0) ||
      arma::any(A21_diagonal <= 0.0)) {
    Rcpp::stop(
      "A PDQ diagonal bridge is not finite and strictly positive; no "
      "clipping is applied."
    );
  }
  const arma::mat A12 = arma::diagmat(A12_diagonal);
  const arma::mat A21 = arma::diagmat(A21_diagonal);
  const arma::mat M1 = fit2.G * A21 * fit1.G_inverse;
  const arma::mat M2 = fit2.G_inverse * A12.t() * fit1.G;
  const arma::mat K1 = 0.5 * (M1 + M1.t());
  const arma::mat K2 = 0.5 * (M2 + M2.t());
  const arma::mat C12 = fit2.G_inverse * A12.t() * fit1.G *
    fit2.G * A21 * fit1.G_inverse;
  const arma::mat K3 = arma::eye(p, p) + C12.t();
  if (!K1.is_finite() || !K2.is_finite() || !K3.is_finite()) {
    Rcpp::stop("A PDQ K matrix is non-finite; no repair is applied.");
  }

  arma::uword cross_zero1 = 0U;
  arma::uword cross_zero2 = 0U;
  const arma::mat cross_sign1 = pdq_cross_signs(
    transformed.x, fit2.location_preconditioned,
    fit1.scale.quantile_working, cross_zero1
  );
  const arma::mat cross_sign2 = pdq_cross_signs(
    transformed.y, fit1.location_preconditioned,
    fit2.scale.quantile_working, cross_zero2
  );
  const arma::mat cross_inner = cross_sign1 * cross_sign2.t();
  const double Rhat = -arma::accu(cross_inner) /
    (static_cast<double>(n1) * static_cast<double>(n2));

  const arma::rowvec fitted_mean1 = arma::mean(fit1.fitted.signs, 0);
  const arma::rowvec fitted_mean2 = arma::mean(fit2.fitted.signs, 0);
  const double Qhat = pdq_quadratic(
    fitted_mean1, fitted_mean2, K1, K2, K3
  );
  const arma::mat within_kernel1 = fit1.fitted.signs * K1 *
    fit1.fitted.signs.t() /
    (static_cast<double>(n1) * static_cast<double>(n1));
  const arma::mat within_kernel2 = fit2.fitted.signs * K2 *
    fit2.fitted.signs.t() /
    (static_cast<double>(n2) * static_cast<double>(n2));
  const arma::mat cross_kernel = -fit1.fitted.signs * K3 *
    fit2.fitted.signs.t() /
    (2.0 * static_cast<double>(n1) * static_cast<double>(n2));
  const arma::vec bias_contributions1 = within_kernel1.diag();
  const arma::vec bias_contributions2 = within_kernel2.diag();
  const double bias1 = arma::accu(bias_contributions1);
  const double bias2 = arma::accu(bias_contributions2);
  const double bhat = bias1 + bias2;
  const double statistic = Rhat - bhat;
  const double fitted_diagonal_deleted = Qhat - bhat;
  if (!std::isfinite(Rhat) || !std::isfinite(Qhat) ||
      !std::isfinite(bhat) || !std::isfinite(statistic)) {
    Rcpp::stop("The PDQ observed statistic is non-finite; no repair is applied.");
  }

  arma::mat H1_zero = within_kernel1;
  arma::mat H2_zero = within_kernel2;
  H1_zero.diag().zeros();
  H2_zero.diag().zeros();
  const double bootstrap_variance = 2.0 * (
    arma::accu(arma::square(H1_zero)) +
    arma::accu(arma::square(H2_zero)) +
    2.0 * arma::accu(arma::square(cross_kernel))
  );
  if (!std::isfinite(bootstrap_variance) || bootstrap_variance <= 0.0) {
    Rcpp::stop(
      "The empirical diagonal-deleted PDQ bootstrap variance is not "
      "finite and strictly positive; no absolute-value or floor repair "
      "is applied."
    );
  }

  arma::vec bootstrap(B);
  Rcpp::IntegerMatrix multipliers;
  if (keep_bootstrap) {
    multipliers = Rcpp::IntegerMatrix(B, static_cast<int>(n1 + n2));
  }
  const std::uint32_t seed_u32 = static_cast<std::uint32_t>(seed);
  arma::rowvec mean_star1(p);
  arma::rowvec mean_star2(p);
  double bootstrap_sum = 0.0;
  double bootstrap_minimum = R_PosInf;
  double bootstrap_maximum = R_NegInf;
  std::uint64_t exceedances = 0U;

  for (int b = 0; b < B; ++b) {
    mean_star1.zeros();
    mean_star2.zeros();
    for (arma::uword i = 0; i < n1 + n2; ++i) {
      const std::uint64_t counter =
        static_cast<std::uint64_t>(b) *
          static_cast<std::uint64_t>(n1 + n2) +
        static_cast<std::uint64_t>(i);
      const int sign =
        (pdq_counter_word(counter, seed_u32) & UINT32_C(1)) != 0U ? 1 : -1;
      if (keep_bootstrap) {
        multipliers(b, static_cast<int>(i)) = sign;
      }
      if (i < n1) {
        mean_star1 += static_cast<double>(sign) * fit1.fitted.signs.row(i);
      } else {
        mean_star2 += static_cast<double>(sign) *
          fit2.fitted.signs.row(i - n1);
      }
    }
    mean_star1 /= static_cast<double>(n1);
    mean_star2 /= static_cast<double>(n2);
    const double value = pdq_quadratic(
      mean_star1, mean_star2, K1, K2, K3
    ) - bhat;
    if (!std::isfinite(value)) {
      Rcpp::stop("A PDQ bootstrap statistic is non-finite.");
    }
    bootstrap(b) = value;
    bootstrap_sum += value;
    bootstrap_minimum = std::min(bootstrap_minimum, value);
    bootstrap_maximum = std::max(bootstrap_maximum, value);
    if (value >= statistic) {
      ++exceedances;
    }
  }

  const double bootstrap_mean = bootstrap_sum / static_cast<double>(B);
  double bootstrap_empirical_variance = 0.0;
  for (int b = 0; b < B; ++b) {
    const double centered = bootstrap(b) - bootstrap_mean;
    bootstrap_empirical_variance += centered * centered;
  }
  bootstrap_empirical_variance /= static_cast<double>(B);
  arma::vec ordered_bootstrap = arma::sort(bootstrap);
  const arma::uword critical_rank = pdq_quantile_rank(
    static_cast<arma::uword>(B), 1.0 - level
  );
  const double critical_value = ordered_bootstrap(critical_rank - 1U);
  const double p_value = (1.0 + static_cast<double>(exceedances)) /
    (static_cast<double>(B) + 1.0);
  const bool primary_rejected = statistic > critical_value;
  const bool p_value_rejected = p_value <= level;

  Rcpp::RObject kept_bootstrap = R_NilValue;
  Rcpp::RObject kept_multipliers = R_NilValue;
  if (keep_bootstrap) {
    kept_bootstrap = Rcpp::wrap(bootstrap);
    kept_multipliers = multipliers;
  }

  return Rcpp::List::create(
    Rcpp::Named("statistic") = statistic,
    Rcpp::Named("R_hat") = Rhat,
    Rcpp::Named("b_hat") = bhat,
    Rcpp::Named("bias1") = bias1,
    Rcpp::Named("bias2") = bias2,
    Rcpp::Named("Q_hat_fitted") = Qhat,
    Rcpp::Named("fitted_diagonal_deleted") = fitted_diagonal_deleted,
    Rcpp::Named("p_value") = p_value,
    Rcpp::Named("critical_value") = critical_value,
    Rcpp::Named("critical_rank") = static_cast<double>(critical_rank),
    Rcpp::Named("primary_rejected") = primary_rejected,
    Rcpp::Named("p_value_rejected") = p_value_rejected,
    Rcpp::Named("decision_disagreement") =
      primary_rejected != p_value_rejected,
    Rcpp::Named("exceedances") = static_cast<double>(exceedances),
    Rcpp::Named("bootstrap") = kept_bootstrap,
    Rcpp::Named("multipliers") = kept_multipliers,
    Rcpp::Named("bootstrap_mean") = bootstrap_mean,
    Rcpp::Named("bootstrap_empirical_variance") =
      bootstrap_empirical_variance,
    Rcpp::Named("bootstrap_minimum") = bootstrap_minimum,
    Rcpp::Named("bootstrap_maximum") = bootstrap_maximum,
    Rcpp::Named("bootstrap_variance_formula") = bootstrap_variance,
    Rcpp::Named("seed") = seed,
    Rcpp::Named("quantile_pair_count1") =
      static_cast<double>(fit1.scale.pair_count),
    Rcpp::Named("quantile_pair_count2") =
      static_cast<double>(fit2.scale.pair_count),
    Rcpp::Named("quantile_rank1") = static_cast<double>(fit1.scale.rank),
    Rcpp::Named("quantile_rank2") = static_cast<double>(fit2.scale.rank),
    Rcpp::Named("quantile_preconditioned1") =
      fit1.scale.quantile_preconditioned.t(),
    Rcpp::Named("quantile_preconditioned2") =
      fit2.scale.quantile_preconditioned.t(),
    Rcpp::Named("quantile_working1") = fit1.scale.quantile_working.t(),
    Rcpp::Named("quantile_working2") = fit2.scale.quantile_working.t(),
    Rcpp::Named("quantile_input1") = fit1.scale.quantile_input.t(),
    Rcpp::Named("quantile_input2") = fit2.scale.quantile_input.t(),
    Rcpp::Named("log_quantile_input1") =
      fit1.scale.log_quantile_input.t(),
    Rcpp::Named("log_quantile_input2") =
      fit2.scale.log_quantile_input.t(),
    Rcpp::Named("D_input1") = fit1.scale.diagonal_input.t(),
    Rcpp::Named("D_input2") = fit2.scale.diagonal_input.t(),
    Rcpp::Named("D_input_canonical1") =
      fit1.scale.diagonal_input_canonical.t(),
    Rcpp::Named("D_input_canonical2") =
      fit2.scale.diagonal_input_canonical.t(),
    Rcpp::Named("D_working1") = fit1.scale.diagonal_working.t(),
    Rcpp::Named("D_working2") = fit2.scale.diagonal_working.t(),
    Rcpp::Named("nonrepresentable_quantiles1") =
      static_cast<double>(fit1.scale.nonrepresentable_quantiles),
    Rcpp::Named("nonrepresentable_quantiles2") =
      static_cast<double>(fit2.scale.nonrepresentable_quantiles),
    Rcpp::Named("nonrepresentable_diagonals1") =
      static_cast<double>(fit1.scale.nonrepresentable_diagonals),
    Rcpp::Named("nonrepresentable_diagonals2") =
      static_cast<double>(fit2.scale.nonrepresentable_diagonals),
    Rcpp::Named("location1") = fit1.location_input.t(),
    Rcpp::Named("location2") = fit2.location_input.t(),
    Rcpp::Named("location_preconditioned1") =
      fit1.location_preconditioned.t(),
    Rcpp::Named("location_preconditioned2") =
      fit2.location_preconditioned.t(),
    Rcpp::Named("location_standardized1") = fit1.median.location.t(),
    Rcpp::Named("location_standardized2") = fit2.median.location.t(),
    Rcpp::Named("standardized_data1") = fit1.standardized_data,
    Rcpp::Named("standardized_data2") = fit2.standardized_data,
    Rcpp::Named("fitted_signs1") = fit1.fitted.signs,
    Rcpp::Named("fitted_signs2") = fit2.fitted.signs,
    Rcpp::Named("fitted_radii1") = fit1.fitted.radii,
    Rcpp::Named("fitted_radii2") = fit2.fitted.radii,
    Rcpp::Named("fitted_inverse_radii1") = fit1.fitted.inverse_radii,
    Rcpp::Named("fitted_inverse_radii2") = fit2.fitted.inverse_radii,
    Rcpp::Named("cross_signs1") = cross_sign1,
    Rcpp::Named("cross_signs2") = cross_sign2,
    Rcpp::Named("cross_inner_products") = cross_inner,
    Rcpp::Named("cross_zero_count1") =
      static_cast<double>(cross_zero1),
    Rcpp::Named("cross_zero_count2") =
      static_cast<double>(cross_zero2),
    Rcpp::Named("Omega1") = fit1.omega,
    Rcpp::Named("Omega2") = fit2.omega,
    Rcpp::Named("G1") = fit1.G,
    Rcpp::Named("G2") = fit2.G,
    Rcpp::Named("G_inverse1") = fit1.G_inverse,
    Rcpp::Named("G_inverse2") = fit2.G_inverse,
    Rcpp::Named("G_eigenvalues1") = fit1.G_eigenvalues,
    Rcpp::Named("G_eigenvalues2") = fit2.G_eigenvalues,
    Rcpp::Named("G_condition1") = fit1.G_condition,
    Rcpp::Named("G_condition2") = fit2.G_condition,
    Rcpp::Named("A12_diagonal") = A12_diagonal.t(),
    Rcpp::Named("A21_diagonal") = A21_diagonal.t(),
    Rcpp::Named("A12_input_diagonal") = A12_input_diagonal.t(),
    Rcpp::Named("A21_input_diagonal") = A21_input_diagonal.t(),
    Rcpp::Named("K1") = K1,
    Rcpp::Named("K2") = K2,
    Rcpp::Named("K3") = K3,
    Rcpp::Named("within_kernel1") = within_kernel1,
    Rcpp::Named("within_kernel2") = within_kernel2,
    Rcpp::Named("cross_kernel") = cross_kernel,
    Rcpp::Named("bias_contributions1") = bias_contributions1,
    Rcpp::Named("bias_contributions2") = bias_contributions2,
    Rcpp::Named("median_iterations1") = fit1.median.iterations,
    Rcpp::Named("median_iterations2") = fit2.median.iterations,
    Rcpp::Named("median_converged1") = fit1.median.converged,
    Rcpp::Named("median_converged2") = fit2.median.converged,
    Rcpp::Named("median_relative_change1") =
      fit1.median.relative_change,
    Rcpp::Named("median_relative_change2") =
      fit2.median.relative_change,
    Rcpp::Named("median_equation_residual1") =
      fit1.median.equation_residual,
    Rcpp::Named("median_equation_residual2") =
      fit2.median.equation_residual,
    Rcpp::Named("median_objective1") = fit1.median.objective,
    Rcpp::Named("median_objective2") = fit2.median.objective,
    Rcpp::Named("median_zero_count1") =
      static_cast<double>(fit1.median.zero_count),
    Rcpp::Named("median_zero_count2") =
      static_cast<double>(fit2.median.zero_count),
    Rcpp::Named("preconditioned_x") = transformed.x,
    Rcpp::Named("preconditioned_y") = transformed.y,
    Rcpp::Named("precondition_magnitude") = transformed.magnitude.t(),
    Rcpp::Named("precondition_midpoint_scaled") =
      transformed.midpoint_scaled.t(),
    Rcpp::Named("precondition_range_scaled") =
      transformed.range_scaled.t(),
    Rcpp::Named("precondition_log_unit") = transformed.log_unit.t()
  );
}
