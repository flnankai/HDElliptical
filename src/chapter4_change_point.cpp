// Chapter 4: high-dimensional change-point methods.
//
// Formula kernels for Wang--Feng DMS, Liu--Feng--Peng--Wang spatial-sign
// scans, and Song--Wen--Feng ERHT.  The kernels deliberately do not repair
// zero radii, non-positive diagonal scales, or non-positive variance
// estimates.  Such cases are returned to R as explicit method failures.

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp17)]]

#include <RcppArmadillo.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <string>

namespace {

long double ch4_cp_norm_vec(const arma::vec& values) {
  long double maximum = 0.0L;
  for (arma::uword j = 0; j < values.n_elem; ++j) {
    const long double value = std::fabs(static_cast<long double>(values(j)));
    if (value > maximum) maximum = value;
  }
  if (maximum == 0.0L) return 0.0L;
  long double sum = 0.0L;
  for (arma::uword j = 0; j < values.n_elem; ++j) {
    const long double ratio = static_cast<long double>(values(j)) / maximum;
    sum += ratio * ratio;
  }
  return maximum * std::sqrt(sum);
}

double ch4_cp_double(const long double value, const std::string& name) {
  const long double largest =
    static_cast<long double>(std::numeric_limits<double>::max());
  if (!std::isfinite(value) || value > largest || value < -largest) {
    Rcpp::stop("%s is outside the finite double range.", name.c_str());
  }
  const double answer = static_cast<double>(value);
  if (value != 0.0L && answer == 0.0) {
    Rcpp::stop("%s underflows the finite double range.", name.c_str());
  }
  return answer;
}

Rcpp::List ch4_cp_failure(const std::string& failure,
                          const int iterations = 0,
                          const double update = R_NaReal) {
  return Rcpp::List::create(
    Rcpp::_["valid"] = false,
    Rcpp::_["stable"] = false,
    Rcpp::_["failure"] = failure,
    Rcpp::_["iterations"] = iterations,
    Rcpp::_["relative_update"] = update
  );
}

arma::vec ch4_cp_column_means(const arma::mat& x) {
  arma::vec result(x.n_cols, arma::fill::zeros);
  for (arma::uword j = 0; j < x.n_cols; ++j) {
    long double value = 0.0L;
    const long double anchor = static_cast<long double>(x(0, j));
    for (arma::uword i = 0; i < x.n_rows; ++i) {
      value += static_cast<long double>(x(i, j)) - anchor;
    }
    result(j) = ch4_cp_double(
      anchor + value / static_cast<long double>(x.n_rows),
      "a column mean"
    );
  }
  return result;
}

} // namespace


//' Coordinate CUSUMs and Rice difference variances
//'
//' Computes the literal gamma=0 and gamma=1/2 CUSUM arrays together with
//' `sum(diff(x)^2)/(2(n-1))`.  Columnwise anchors and long-double accumulators
//' reduce cancellation under a common translation.
//'
//' @param x Numeric observations in rows.
//' @return A list of two CUSUM matrices and the difference variance.
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List ch4_cp_cusum_cpp(const Rcpp::NumericMatrix& x) {
  const int n = x.nrow();
  const int p = x.ncol();
  if (n < 2 || p < 1) Rcpp::stop("x must have at least two rows and one column.");

  Rcpp::NumericMatrix c0(n - 1, p);
  Rcpp::NumericMatrix chalf(n - 1, p);
  Rcpp::NumericVector s2(p);

  for (int j = 0; j < p; ++j) {
    const long double anchor = static_cast<long double>(x(0, j));
    long double total = 0.0L;
    long double diff_sum = 0.0L;
    for (int i = 0; i < n; ++i) {
      total += static_cast<long double>(x(i, j)) - anchor;
      if (i > 0) {
        const long double difference =
          static_cast<long double>(x(i, j)) -
          static_cast<long double>(x(i - 1, j));
        diff_sum += difference * difference;
      }
    }
    s2[j] = ch4_cp_double(
      diff_sum / (2.0L * static_cast<long double>(n - 1)),
      "a Rice difference variance"
    );

    long double partial = 0.0L;
    for (int k = 1; k < n; ++k) {
      partial += static_cast<long double>(x(k - 1, j)) - anchor;
      const long double raw =
        (partial - static_cast<long double>(k) * total /
          static_cast<long double>(n)) /
        std::sqrt(static_cast<long double>(n));
      const long double u = static_cast<long double>(k) /
        static_cast<long double>(n);
      const long double weight = 1.0L / std::sqrt(u * (1.0L - u));
      c0(k - 1, j) = ch4_cp_double(raw, "an unweighted CUSUM");
      chalf(k - 1, j) = ch4_cp_double(raw * weight, "a weighted CUSUM");
    }
  }

  return Rcpp::List::create(
    Rcpp::_["cusum0"] = c0,
    Rcpp::_["cusum_half"] = chalf,
    Rcpp::_["difference_variance"] = s2
  );
}


//' Wang--Feng finite-difference moment estimators
//'
//' Implements the displayed leave-four and leave-three estimators in the
//' primary paper.  The exclusion set is literally
//' `{2,...,n} \\ {i1,...,im}` and no bridging difference is inserted.
//'
//' @param x Numeric observations in rows.
//' @return Trace, radial fourth-moment, and term-level diagnostics.
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List ch4_cp_dms_moments_cpp(const Rcpp::NumericMatrix& x) {
  const int n = x.nrow();
  const int p = x.ncol();
  if (n < 6) {
    return ch4_cp_failure("DMS leave-four variances require at least six observations.");
  }
  Rcpp::NumericVector trace_terms(n - 3);
  Rcpp::NumericVector fourth_terms(n - 2);

  auto leave_variance = [&](const int first_obs_zero,
                            const int last_obs_zero,
                            const int coordinate,
                            long double& answer) -> bool {
    long double sum = 0.0L;
    int count = 0;
    // d indexes the upper observation of a first difference in zero base.
    for (int d = 1; d < n; ++d) {
      if (d >= first_obs_zero && d <= last_obs_zero) continue;
      const long double difference =
        static_cast<long double>(x(d, coordinate)) -
        static_cast<long double>(x(d - 1, coordinate));
      sum += difference * difference;
      ++count;
    }
    if (count <= 0) return false;
    answer = sum / (2.0L * static_cast<long double>(count));
    return std::isfinite(answer) && answer > 0.0L;
  };

  long double trace_sum = 0.0L;
  for (int i = 0; i < n - 3; ++i) {
    long double inner = 0.0L;
    for (int j = 0; j < p; ++j) {
      long double variance = 0.0L;
      if (!leave_variance(i, i + 3, j, variance)) {
        return ch4_cp_failure(
          "A DMS leave-four diagonal variance is non-positive or non-finite."
        );
      }
      const long double left =
        static_cast<long double>(x(i, j)) -
        static_cast<long double>(x(i + 1, j));
      const long double right =
        static_cast<long double>(x(i + 2, j)) -
        static_cast<long double>(x(i + 3, j));
      inner += left * right / variance;
    }
    const long double term = inner * inner;
    trace_terms[i] = ch4_cp_double(term, "a DMS trace term");
    trace_sum += term;
  }
  const long double trace_hat =
    trace_sum / (4.0L * static_cast<long double>(n - 3));

  long double fourth_sum = 0.0L;
  for (int i = 0; i < n - 2; ++i) {
    long double inner = 0.0L;
    for (int j = 0; j < p; ++j) {
      long double variance = 0.0L;
      if (!leave_variance(i, i + 2, j, variance)) {
        return ch4_cp_failure(
          "A DMS leave-three diagonal variance is non-positive or non-finite."
        );
      }
      const long double left =
        static_cast<long double>(x(i, j)) -
        static_cast<long double>(x(i + 1, j));
      const long double right =
        static_cast<long double>(x(i + 1, j)) -
        static_cast<long double>(x(i + 2, j));
      inner += left * right / variance;
    }
    const long double term = inner * inner;
    fourth_terms[i] = ch4_cp_double(term, "a DMS fourth-moment term");
    fourth_sum += term;
  }
  const long double fourth_hat =
    fourth_sum / static_cast<long double>(n - 2) - 3.0L * trace_hat;

  return Rcpp::List::create(
    Rcpp::_["valid"] = true,
    Rcpp::_["failure"] = R_NilValue,
    Rcpp::_["trace_hat"] = ch4_cp_double(trace_hat, "the DMS trace estimate"),
    Rcpp::_["fourth_hat"] = ch4_cp_double(fourth_hat, "the DMS fourth-moment estimate"),
    Rcpp::_["trace_terms"] = trace_terms,
    Rcpp::_["fourth_terms"] = fourth_terms,
    Rcpp::_["leaveout_definition"] =
      "literal displayed A_m={2,...,n}\\{i_1,...,i_m}; no bridge"
  );
}


//' Joint scaled spatial-median and diagonal fixed point
//'
//' Implements the three updates stated by Liu, Feng, Peng, and Wang.  The
//' common diagonal scale inherited from sample-variance initialization is
//' retained; only relative diagonal scale is identified.
//'
//' @param x Numeric observations in rows.
//' @param tol Relative-update tolerance.
//' @param max_iter Maximum iterations.
//' @param zero_tol Exact/near-zero residual tolerance; zero is the paper
//'   convention and positive values are user-requested diagnostics.
//' @return Fixed point and convergence diagnostics.
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List ch4_cp_scaled_hr_cpp(const arma::mat& x,
                                const double tol = 1e-8,
                                const int max_iter = 1000,
                                const double zero_tol = 0.0) {
  const int n = static_cast<int>(x.n_rows);
  const int p = static_cast<int>(x.n_cols);
  if (n < 2 || p < 1) return ch4_cp_failure("A scaled HR fit needs at least two rows.");
  if (!(tol > 0.0) || !std::isfinite(tol) || max_iter < 1 ||
      !(zero_tol >= 0.0) || !std::isfinite(zero_tol)) {
    return ch4_cp_failure("Invalid scaled HR iteration controls.");
  }

  arma::vec theta;
  try {
    theta = ch4_cp_column_means(x);
  } catch (std::exception& error) {
    return ch4_cp_failure(error.what());
  }
  arma::vec diagonal(p, arma::fill::zeros);
  for (int j = 0; j < p; ++j) {
    long double sum = 0.0L;
    for (int i = 0; i < n; ++i) {
      const long double value =
        static_cast<long double>(x(i, j)) -
        static_cast<long double>(theta(j));
      sum += value * value;
    }
    const long double variance = sum / static_cast<long double>(n - 1);
    if (!std::isfinite(variance) || variance <= 0.0L) {
      return ch4_cp_failure(
        "A sample-variance diagonal initializer is non-positive or non-finite."
      );
    }
    try {
      diagonal(j) = ch4_cp_double(variance, "an HR diagonal initializer");
    } catch (std::exception& error) {
      return ch4_cp_failure(error.what());
    }
  }

  double relative_update = R_PosInf;
  double minimum_distance = R_PosInf;
  int iterations = 0;
  bool stable = false;
  arma::mat signs(n, p, arma::fill::zeros);

  for (int iteration = 1; iteration <= max_iter; ++iteration) {
    iterations = iteration;
    arma::vec score(p, arma::fill::zeros);
    arma::vec sign_square(p, arma::fill::zeros);
    long double inverse_distance_sum = 0.0L;
    minimum_distance = R_PosInf;

    for (int i = 0; i < n; ++i) {
      arma::vec standardized(p);
      for (int j = 0; j < p; ++j) {
        standardized(j) = (x(i, j) - theta(j)) / std::sqrt(diagonal(j));
      }
      const long double distance = ch4_cp_norm_vec(standardized);
      if (!std::isfinite(distance) || distance <= static_cast<long double>(zero_tol)) {
        return ch4_cp_failure(
          "A scaled HR residual has zero or non-finite norm; no perturbation was applied.",
          iteration,
          relative_update
        );
      }
      const double distance_d = static_cast<double>(distance);
      minimum_distance = std::min(minimum_distance, distance_d);
      inverse_distance_sum += 1.0L / distance;
      for (int j = 0; j < p; ++j) {
        const double value = standardized(j) / distance_d;
        signs(i, j) = value;
        score(j) += value;
        sign_square(j) += value * value;
      }
    }
    if (!std::isfinite(inverse_distance_sum) || inverse_distance_sum <= 0.0L) {
      return ch4_cp_failure("The scaled HR inverse-distance sum is invalid.", iteration);
    }

    arma::vec theta_new = theta;
    arma::vec diagonal_new = diagonal;
    for (int j = 0; j < p; ++j) {
      const long double increment =
        std::sqrt(static_cast<long double>(diagonal(j))) *
        static_cast<long double>(score(j)) / inverse_distance_sum;
      try {
        theta_new(j) = ch4_cp_double(
          static_cast<long double>(theta(j)) + increment,
          "a scaled HR location update"
        );
        diagonal_new(j) = ch4_cp_double(
          static_cast<long double>(p) * static_cast<long double>(diagonal(j)) *
          static_cast<long double>(sign_square(j)) /
          static_cast<long double>(n),
          "a scaled HR diagonal update"
        );
      } catch (std::exception& error) {
        return ch4_cp_failure(error.what(), iteration, relative_update);
      }
      if (!(diagonal_new(j) > 0.0) || !std::isfinite(diagonal_new(j))) {
        return ch4_cp_failure(
          "A scaled HR diagonal update is non-positive or non-finite.",
          iteration,
          relative_update
        );
      }
    }

    double theta_update = 0.0;
    for (int j = 0; j < p; ++j) {
      theta_update = std::max(
        theta_update,
        std::fabs(theta_new(j) - theta(j)) / std::sqrt(diagonal(j))
      );
    }
    double diagonal_update = 0.0;
    for (int j = 0; j < p; ++j) {
      diagonal_update = std::max(
        diagonal_update,
        std::fabs(std::log(diagonal_new(j)) - std::log(diagonal(j)))
      );
    }
    relative_update = std::max(theta_update, diagonal_update);
    theta = theta_new;
    diagonal = diagonal_new;
    if (!std::isfinite(relative_update)) {
      return ch4_cp_failure("The scaled HR relative update is non-finite.", iteration);
    }
    if (relative_update <= tol) {
      stable = true;
      break;
    }
  }

  // Diagnostics are evaluated at the returned iterate, rather than recycled
  // from the preceding map.
  arma::vec score(p, arma::fill::zeros);
  arma::vec sign_square(p, arma::fill::zeros);
  minimum_distance = R_PosInf;
  for (int i = 0; i < n; ++i) {
    arma::vec standardized(p);
    for (int j = 0; j < p; ++j) {
      standardized(j) = (x(i, j) - theta(j)) / std::sqrt(diagonal(j));
    }
    const long double distance = ch4_cp_norm_vec(standardized);
    if (!std::isfinite(distance) || distance <= static_cast<long double>(zero_tol)) {
      return ch4_cp_failure(
        "The final scaled HR iterate has a zero or non-finite residual.",
        iterations,
        relative_update
      );
    }
    const double distance_d = static_cast<double>(distance);
    minimum_distance = std::min(minimum_distance, distance_d);
    for (int j = 0; j < p; ++j) {
      const double value = standardized(j) / distance_d;
      signs(i, j) = value;
      score(j) += value;
      sign_square(j) += value * value;
    }
  }
  score /= static_cast<double>(n);
  sign_square /= static_cast<double>(n);
  const double score_l2 = arma::norm(score, 2);
  const double score_inf = arma::abs(score).max();
  const double diagonal_residual =
    arma::abs(static_cast<double>(p) * sign_square - 1.0).max();

  return Rcpp::List::create(
    Rcpp::_["valid"] = stable,
    Rcpp::_["stable"] = stable,
    Rcpp::_["failure"] = stable ? R_NilValue : Rcpp::wrap(
      "The scaled HR iteration did not stabilize within max_iter."
    ),
    Rcpp::_["location"] = theta,
    Rcpp::_["diagonal"] = diagonal,
    Rcpp::_["signs"] = signs,
    Rcpp::_["iterations"] = iterations,
    Rcpp::_["relative_update"] = relative_update,
    Rcpp::_["score_l2"] = score_l2,
    Rcpp::_["score_infinity"] = score_inf,
    Rcpp::_["diagonal_residual"] = diagonal_residual,
    Rcpp::_["minimum_distance"] = minimum_distance,
    Rcpp::_["scale_identification"] = "relative diagonal scale only"
  );
}


//' Spatial median with a coincident-point certificate
//'
//' Uses the modified Weiszfeld map.  At an iterate coinciding with observations,
//' the exact subgradient condition is checked rather than adding jitter.
//'
//' @param x Numeric observations in rows.
//' @param tol Relative-update tolerance.
//' @param max_iter Maximum iterations.
//' @param zero_tol Coincidence tolerance, with zero giving exact coincidence.
//' @return Median and convergence diagnostics.
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List ch4_cp_spatial_median_cpp(const arma::mat& x,
                                     const double tol = 1e-8,
                                     const int max_iter = 1000,
                                     const double zero_tol = 0.0) {
  const int n = static_cast<int>(x.n_rows);
  const int p = static_cast<int>(x.n_cols);
  if (n < 1 || p < 1) return ch4_cp_failure("A spatial median needs non-empty data.");
  if (!(tol > 0.0) || !std::isfinite(tol) || max_iter < 1 ||
      !(zero_tol >= 0.0) || !std::isfinite(zero_tol)) {
    return ch4_cp_failure("Invalid spatial-median controls.");
  }

  arma::vec theta;
  try {
    theta = ch4_cp_column_means(x);
  } catch (std::exception& error) {
    return ch4_cp_failure(error.what());
  }
  double relative_update = R_PosInf;
  double score_distance = R_PosInf;
  double minimum_distance = R_PosInf;
  int coincident = 0;
  int iterations = 0;
  bool stable = false;

  // A data point can itself be an exact spatial median.  Check its exact
  // subgradient condition first, so a Weiszfeld sequence approaching a
  // nonsmooth solution is returned at that observed point rather than being
  // falsely certified merely because its steps have become tiny.
  for (int candidate = 0; candidate < n && !stable; ++candidate) {
    arma::vec candidate_theta = x.row(candidate).t();
    arma::vec candidate_score(p, arma::fill::zeros);
    int candidate_coincident = 0;
    for (int i = 0; i < n; ++i) {
      arma::vec residual = x.row(i).t() - candidate_theta;
      const long double distance = ch4_cp_norm_vec(residual);
      if (!std::isfinite(distance)) {
        return ch4_cp_failure("A spatial-median residual norm is non-finite.");
      }
      if (distance <= static_cast<long double>(zero_tol)) {
        ++candidate_coincident;
      } else {
        candidate_score += residual / static_cast<double>(distance);
      }
    }
    const double candidate_score_norm = static_cast<double>(
      ch4_cp_norm_vec(candidate_score)
    );
    if (candidate_score_norm <= static_cast<double>(candidate_coincident)) {
      theta = candidate_theta;
      coincident = candidate_coincident;
      minimum_distance = 0.0;
      score_distance = 0.0;
      relative_update = 0.0;
      stable = true;
    }
  }

  for (int iteration = 1; iteration <= max_iter && !stable; ++iteration) {
    iterations = iteration;
    arma::vec weighted_sum(p, arma::fill::zeros);
    arma::vec direction_sum(p, arma::fill::zeros);
    long double weight_sum = 0.0L;
    coincident = 0;
    minimum_distance = R_PosInf;

    for (int i = 0; i < n; ++i) {
      arma::vec residual = x.row(i).t() - theta;
      const long double distance = ch4_cp_norm_vec(residual);
      if (!std::isfinite(distance)) {
        return ch4_cp_failure("A spatial-median residual norm is non-finite.", iteration);
      }
      minimum_distance = std::min(minimum_distance, static_cast<double>(distance));
      if (distance <= static_cast<long double>(zero_tol)) {
        ++coincident;
        continue;
      }
      const double distance_d = static_cast<double>(distance);
      weighted_sum += x.row(i).t() / distance_d;
      direction_sum += residual / distance_d;
      weight_sum += 1.0L / distance;
    }

    const double direction_norm = static_cast<double>(ch4_cp_norm_vec(direction_sum));
    score_distance = std::max(0.0, direction_norm - static_cast<double>(coincident)) /
      static_cast<double>(n);
    if (coincident > 0 && direction_norm <= static_cast<double>(coincident)) {
      stable = true;
      relative_update = 0.0;
      break;
    }
    if (!std::isfinite(weight_sum) || weight_sum <= 0.0L) {
      return ch4_cp_failure(
        "The spatial-median inverse-distance sum is invalid.",
        iteration,
        relative_update
      );
    }

    arma::vec target = weighted_sum / static_cast<double>(weight_sum);
    arma::vec theta_new = target;
    if (coincident > 0) {
      if (!(direction_norm > static_cast<double>(coincident))) {
        return ch4_cp_failure("The coincident-point spatial-median map is undefined.", iteration);
      }
      const double eta = static_cast<double>(coincident) / direction_norm;
      theta_new = eta * theta + (1.0 - eta) * target;
    }
    long double residual_scale_sum = 0.0L;
    for (int i = 0; i < n; ++i) {
      residual_scale_sum += ch4_cp_norm_vec(x.row(i).t() - theta);
    }
    const long double residual_scale = residual_scale_sum /
      static_cast<long double>(n);
    if (!std::isfinite(residual_scale) || residual_scale <= 0.0L) {
      return ch4_cp_failure(
        "The spatial-median residual scale is non-positive or non-finite.",
        iteration
      );
    }
    relative_update = static_cast<double>(
      ch4_cp_norm_vec(theta_new - theta) / residual_scale
    );
    theta = theta_new;
    if (!std::isfinite(relative_update) || !theta.is_finite()) {
      return ch4_cp_failure("The spatial-median update is non-finite.", iteration);
    }
    if (relative_update <= tol) {
      stable = true;
      break;
    }
  }

  // Re-evaluate the subgradient-distance certificate at the returned point.
  arma::vec direction_sum(p, arma::fill::zeros);
  coincident = 0;
  minimum_distance = R_PosInf;
  for (int i = 0; i < n; ++i) {
    arma::vec residual = x.row(i).t() - theta;
    const long double distance = ch4_cp_norm_vec(residual);
    if (!std::isfinite(distance)) {
      return ch4_cp_failure("The final spatial-median residual is non-finite.", iterations);
    }
    minimum_distance = std::min(minimum_distance, static_cast<double>(distance));
    if (distance <= static_cast<long double>(zero_tol)) {
      ++coincident;
    } else {
      direction_sum += residual / static_cast<double>(distance);
    }
  }
  score_distance = std::max(
    0.0,
    static_cast<double>(ch4_cp_norm_vec(direction_sum)) -
      static_cast<double>(coincident)
  ) / static_cast<double>(n);

  return Rcpp::List::create(
    Rcpp::_["valid"] = stable,
    Rcpp::_["stable"] = stable,
    Rcpp::_["failure"] = stable ? R_NilValue : Rcpp::wrap(
      "The spatial-median iteration did not stabilize within max_iter."
    ),
    Rcpp::_["location"] = theta,
    Rcpp::_["iterations"] = iterations,
    Rcpp::_["relative_update"] = relative_update,
    Rcpp::_["score_residual"] = score_distance,
    Rcpp::_["minimum_distance"] = minimum_distance,
    Rcpp::_["coincident_observations"] = coincident,
    Rcpp::_["convergence_basis"] =
      "exact observed-point subgradient certificate or modified Weiszfeld residual-scale-relative update"
  );
}


//' Ordered-pair squared inner-product sum
//'
//' @param signs Matrix with observations in rows.
//' @return `sum_{i != j} (u_i' u_j)^2`.
//' @keywords internal
// [[Rcpp::export]]
double ch4_cp_ordered_pair_square_sum_cpp(const arma::mat& signs) {
  const arma::uword n = signs.n_rows;
  long double answer = 0.0L;
  for (arma::uword i = 0; i < n; ++i) {
    for (arma::uword j = 0; j < n; ++j) {
      if (i == j) continue;
      long double inner = 0.0L;
      for (arma::uword k = 0; k < signs.n_cols; ++k) {
        inner += static_cast<long double>(signs(i, k)) *
          static_cast<long double>(signs(j, k));
      }
      answer += inner * inner;
    }
  }
  return ch4_cp_double(answer, "an ordered-pair sign sum");
}


//' ERHT companion center and variance
//'
//' Computes the ordered off-diagonal expression directly, avoiding
//' subtraction of two nearly equal non-negative matrix products.
//'
//' @param companion Symmetric companion matrix A.
//' @param beta Score-CUSUM weights.
//' @return Kappa, sigma squared, and the ordered off-diagonal sum.
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List ch4_cp_erht_moments_cpp(const arma::mat& companion,
                                   const arma::vec& beta) {
  const arma::uword m = companion.n_rows;
  if (companion.n_cols != m || beta.n_elem != m) {
    return ch4_cp_failure("ERHT companion dimensions do not conform.");
  }
  long double kappa = 0.0L;
  long double offdiag = 0.0L;
  for (arma::uword i = 0; i < m; ++i) {
    const long double bi2 = static_cast<long double>(beta(i)) * beta(i);
    kappa += bi2 * static_cast<long double>(companion(i, i));
    for (arma::uword j = 0; j < m; ++j) {
      if (i == j) continue;
      const long double bj2 = static_cast<long double>(beta(j)) * beta(j);
      const long double a = static_cast<long double>(companion(i, j));
      offdiag += bi2 * bj2 * a * a;
    }
  }
  const long double sigma2 = 2.0L * static_cast<long double>(m) * offdiag;
  return Rcpp::List::create(
    Rcpp::_["valid"] = std::isfinite(kappa) && std::isfinite(sigma2) &&
      sigma2 > 0.0L,
    Rcpp::_["failure"] = (std::isfinite(sigma2) && sigma2 > 0.0L) ?
      R_NilValue : Rcpp::wrap(
        "The literal ERHT ordered-pair variance is non-positive or non-finite."
      ),
    Rcpp::_["kappa"] = ch4_cp_double(kappa, "ERHT kappa"),
    Rcpp::_["sigma2"] = ch4_cp_double(sigma2, "ERHT sigma squared"),
    Rcpp::_["ordered_offdiagonal_sum"] =
      ch4_cp_double(offdiag, "ERHT ordered off-diagonal sum")
  );
}
