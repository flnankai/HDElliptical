// Chapter 3: proportionality tests and robust spatial-sign precision.
//
// The testing kernels below preserve the ordered-index normalisations in the
// primary papers.  The SCLIME kernel solves each column LP with a
// Chambolle--Pock primal-dual iteration and certifies primal feasibility,
// stationarity, dual feasibility, and the primal-dual gap.  The SGLASSO
// kernel uses a positive-definite backtracking proximal-gradient method and
// certifies the full (diagonal included) l1 subgradient KKT equation.
// No ridge, eigenvalue floor, pseudoinverse, or post-hoc feasibility repair is
// applied.

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp17)]]

#include <RcppArmadillo.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <string>
#include <vector>

namespace {

long double ch3pp_dot_ld(const arma::rowvec& x, const arma::rowvec& y) {
  long double answer = 0.0L;
  for (arma::uword j = 0; j < x.n_elem; ++j) {
    answer += static_cast<long double>(x(j)) *
      static_cast<long double>(y(j));
  }
  return answer;
}

double ch3pp_checked_double(const long double value,
                            const char* quantity) {
  const long double maximum = static_cast<long double>(
    std::numeric_limits<double>::max()
  );
  if (!std::isfinite(value) || value > maximum || value < -maximum) {
    Rcpp::stop(
      "%s is outside the finite double range; no clipping or numerical "
      "repair is applied.", quantity
    );
  }
  return static_cast<double>(value);
}

void ch3pp_validate_finite_matrix(const arma::mat& x,
                                  const char* name) {
  if (!x.is_finite()) {
    Rcpp::stop("`%s` must contain only finite values.", name);
  }
}

double ch3pp_max_abs_mat(const arma::mat& x) {
  if (x.n_elem == 0) {
    return 0.0;
  }
  return arma::abs(x).max();
}

double ch3pp_max_abs_vec(const arma::vec& x) {
  if (x.n_elem == 0) {
    return 0.0;
  }
  return arma::abs(x).max();
}

arma::vec ch3pp_soft_threshold_vec(const arma::vec& x,
                                   const double threshold) {
  arma::vec answer = arma::sign(x) %
    arma::max(arma::abs(x) - threshold, arma::zeros<arma::vec>(x.n_elem));
  return answer;
}

arma::mat ch3pp_soft_threshold_mat(const arma::mat& x,
                                   const double threshold) {
  arma::mat answer = arma::sign(x) %
    arma::max(arma::abs(x) - threshold,
              arma::zeros<arma::mat>(x.n_rows, x.n_cols));
  return answer;
}

struct Ch3ppPairSigns {
  arma::mat signs;
  std::vector<arma::uword> first;
  std::vector<arma::uword> second;
  arma::uword zero_pairs;
};

Ch3ppPairSigns ch3pp_pair_signs(const arma::mat& x,
                                const double zero_tolerance) {
  const arma::uword n = x.n_rows;
  const arma::uword p = x.n_cols;
  const arma::uword number_pairs = n * (n - 1) / 2;
  Ch3ppPairSigns result;
  result.signs.zeros(number_pairs, p);
  result.first.reserve(number_pairs);
  result.second.reserve(number_pairs);
  result.zero_pairs = 0;

  arma::uword row = 0;
  for (arma::uword i = 0; i < n; ++i) {
    for (arma::uword j = i + 1; j < n; ++j) {
      result.first.push_back(i);
      result.second.push_back(j);

      long double scale = 0.0L;
      for (arma::uword k = 0; k < p; ++k) {
        const long double difference =
          static_cast<long double>(x(i, k)) -
          static_cast<long double>(x(j, k));
        scale = std::max(scale, std::abs(difference));
      }
      if (scale == 0.0L) {
        ++result.zero_pairs;
        ++row;
        continue;
      }

      long double square_sum = 0.0L;
      for (arma::uword k = 0; k < p; ++k) {
        const long double difference =
          static_cast<long double>(x(i, k)) -
          static_cast<long double>(x(j, k));
        const long double standardized = difference / scale;
        square_sum += standardized * standardized;
      }
      const long double radius = scale * std::sqrt(square_sum);
      if (radius <= static_cast<long double>(zero_tolerance)) {
        ++result.zero_pairs;
        ++row;
        continue;
      }
      const long double normalized_radius = std::sqrt(square_sum);
      for (arma::uword k = 0; k < p; ++k) {
        const long double difference =
          static_cast<long double>(x(i, k)) -
          static_cast<long double>(x(j, k));
        result.signs(row, k) = ch3pp_checked_double(
          (difference / scale) / normalized_radius,
          "A pairwise spatial-sign coordinate"
        );
      }
      ++row;
    }
  }
  return result;
}

struct Ch3ppSclimeColumn {
  arma::vec solution;
  arma::vec dual;
  int iterations;
  bool converged;
  double relative_update;
  double primal_violation;
  double stationarity_residual;
  double dual_violation;
  double primal_objective;
  double dual_objective;
  double duality_gap;
  double relative_gap;
};

Ch3ppSclimeColumn ch3pp_sclime_column(const arma::mat& a,
                                      const arma::vec& target,
                                      const double lambda,
                                      const double tolerance,
                                      const int max_iterations,
                                      const double operator_norm) {
  const arma::uword p = a.n_cols;
  const double step = 0.99 / operator_norm;
  arma::vec primal(p, arma::fill::zeros);
  arma::vec extrapolated(p, arma::fill::zeros);
  arma::vec dual(p, arma::fill::zeros);

  Ch3ppSclimeColumn result;
  result.solution = primal;
  result.dual = dual;
  result.iterations = 0;
  result.converged = false;
  result.relative_update = std::numeric_limits<double>::infinity();
  result.primal_violation = std::numeric_limits<double>::infinity();
  result.stationarity_residual = std::numeric_limits<double>::infinity();
  result.dual_violation = std::numeric_limits<double>::infinity();
  result.primal_objective = 0.0;
  result.dual_objective = -std::numeric_limits<double>::infinity();
  result.duality_gap = std::numeric_limits<double>::infinity();
  result.relative_gap = std::numeric_limits<double>::infinity();

  for (int iteration = 1; iteration <= max_iterations; ++iteration) {
    const arma::vec dual_input = dual + step * (a * extrapolated);
    arma::vec projected = dual_input / step;
    for (arma::uword k = 0; k < p; ++k) {
      projected(k) = std::max(
        target(k) - lambda,
        std::min(target(k) + lambda, projected(k))
      );
    }
    const arma::vec dual_new = dual_input - step * projected;
    const arma::vec primal_input = primal - step * (a.t() * dual_new);
    const arma::vec primal_new = ch3pp_soft_threshold_vec(primal_input, step);
    const arma::vec extrapolated_new = 2.0 * primal_new - primal;

    const double denominator = std::max(1.0, ch3pp_max_abs_vec(primal));
    const arma::vec primal_difference = primal_new - primal;
    result.relative_update = ch3pp_max_abs_vec(primal_difference) /
      denominator;
    primal = primal_new;
    dual = dual_new;
    extrapolated = extrapolated_new;

    const arma::vec residual = a * primal - target;
    result.primal_violation = std::max(
      0.0, ch3pp_max_abs_vec(residual) - lambda
    );
    const arma::vec stationarity = a.t() * dual;
    result.dual_violation = std::max(
      0.0, ch3pp_max_abs_vec(stationarity) - 1.0
    );
    double stationarity_residual = 0.0;
    for (arma::uword k = 0; k < p; ++k) {
      double component = 0.0;
      if (std::abs(primal(k)) > tolerance) {
        component = std::abs(
          stationarity(k) + (primal(k) > 0.0 ? 1.0 : -1.0)
        );
      } else {
        component = std::max(0.0, std::abs(stationarity(k)) - 1.0);
      }
      stationarity_residual = std::max(stationarity_residual, component);
    }
    result.stationarity_residual = stationarity_residual;
    result.primal_objective = arma::accu(arma::abs(primal));
    result.dual_objective = -arma::dot(target, dual) -
      lambda * arma::accu(arma::abs(dual));
    result.duality_gap = result.primal_objective - result.dual_objective;
    const double gap_scale = std::max(
      1.0,
      std::max(std::abs(result.primal_objective),
               std::abs(result.dual_objective))
    );
    result.relative_gap = std::abs(result.duality_gap) / gap_scale;
    result.iterations = iteration;

    if (result.primal_violation <= tolerance &&
        result.dual_violation <= tolerance &&
        result.stationarity_residual <= tolerance &&
        result.relative_gap <= tolerance) {
      result.converged = true;
      break;
    }
  }

  result.solution = primal;
  result.dual = dual;
  return result;
}

double ch3pp_logdet_spd(const arma::mat& x) {
  double log_determinant = 0.0;
  double sign = 0.0;
  if (!arma::log_det(log_determinant, sign, x) || sign <= 0.0 ||
      !std::isfinite(log_determinant)) {
    return std::numeric_limits<double>::quiet_NaN();
  }
  return log_determinant;
}

double ch3pp_sglasso_objective(const arma::mat& a,
                               const arma::mat& v,
                               const double lambda) {
  const double log_determinant = ch3pp_logdet_spd(v);
  if (!std::isfinite(log_determinant)) {
    return std::numeric_limits<double>::infinity();
  }
  return arma::trace(a * v) - log_determinant +
    lambda * arma::accu(arma::abs(v));
}

double ch3pp_sglasso_kkt(const arma::mat& a,
                         const arma::mat& v,
                         const double lambda,
                         arma::mat& inverse) {
  if (!arma::inv_sympd(inverse, v)) {
    return std::numeric_limits<double>::infinity();
  }
  const arma::mat gradient = a - inverse;
  double residual = 0.0;
  for (arma::uword i = 0; i < v.n_rows; ++i) {
    for (arma::uword j = 0; j < v.n_cols; ++j) {
      double component = 0.0;
      if (v(i, j) > 0.0) {
        component = std::abs(gradient(i, j) + lambda);
      } else if (v(i, j) < 0.0) {
        component = std::abs(gradient(i, j) - lambda);
      } else {
        component = std::max(0.0, std::abs(gradient(i, j)) - lambda);
      }
      residual = std::max(residual, component);
    }
  }
  return residual;
}

} // namespace


// [[Rcpp::export]]
Rcpp::List cpp_ch3pp_sscm_components(const arma::mat& signs_x,
                                     const arma::mat& signs_y) {
  ch3pp_validate_finite_matrix(signs_x, "signs_x");
  ch3pp_validate_finite_matrix(signs_y, "signs_y");
  if (signs_x.n_cols == 0 || signs_x.n_cols != signs_y.n_cols) {
    Rcpp::stop("The two sign matrices must have the same positive width.");
  }
  if (signs_x.n_rows < 2 || signs_y.n_rows < 2) {
    Rcpp::stop("Each sign matrix must have at least two rows.");
  }

  long double sum_x = 0.0L;
  long double sum_y = 0.0L;
  long double sum_cross = 0.0L;
  for (arma::uword i = 0; i < signs_x.n_rows; ++i) {
    for (arma::uword j = 0; j < signs_x.n_rows; ++j) {
      if (i == j) continue;
      const long double dot = ch3pp_dot_ld(signs_x.row(i), signs_x.row(j));
      sum_x += dot * dot;
    }
  }
  for (arma::uword i = 0; i < signs_y.n_rows; ++i) {
    for (arma::uword j = 0; j < signs_y.n_rows; ++j) {
      if (i == j) continue;
      const long double dot = ch3pp_dot_ld(signs_y.row(i), signs_y.row(j));
      sum_y += dot * dot;
    }
  }
  for (arma::uword i = 0; i < signs_x.n_rows; ++i) {
    for (arma::uword j = 0; j < signs_y.n_rows; ++j) {
      const long double dot = ch3pp_dot_ld(signs_x.row(i), signs_y.row(j));
      sum_cross += dot * dot;
    }
  }

  const long double n1 = static_cast<long double>(signs_x.n_rows);
  const long double n2 = static_cast<long double>(signs_y.n_rows);
  return Rcpp::List::create(
    Rcpp::Named("ordered_sum_x") = ch3pp_checked_double(
      sum_x, "The first ordered-pair sum"
    ),
    Rcpp::Named("ordered_sum_y") = ch3pp_checked_double(
      sum_y, "The second ordered-pair sum"
    ),
    Rcpp::Named("cross_sum") = ch3pp_checked_double(
      sum_cross, "The cross-pair sum"
    ),
    Rcpp::Named("A") = ch3pp_checked_double(
      sum_x / (n1 * (n1 - 1.0L)), "The first SSCM component"
    ),
    Rcpp::Named("B") = ch3pp_checked_double(
      sum_y / (n2 * (n2 - 1.0L)), "The second SSCM component"
    ),
    Rcpp::Named("C") = ch3pp_checked_double(
      sum_cross / (n1 * n2), "The cross SSCM component"
    )
  );
}


// [[Rcpp::export]]
Rcpp::List cpp_ch3pp_spatial_rank_components(const arma::mat& x,
                                             const arma::mat& y,
                                             const double zero_tolerance) {
  ch3pp_validate_finite_matrix(x, "x");
  ch3pp_validate_finite_matrix(y, "y");
  if (x.n_cols == 0 || x.n_cols != y.n_cols) {
    Rcpp::stop("`x` and `y` must have the same positive number of columns.");
  }
  if (x.n_rows < 4 || y.n_rows < 4) {
    Rcpp::stop("The fourth-order spatial-rank statistic needs four rows per group.");
  }
  if (!std::isfinite(zero_tolerance) || zero_tolerance < 0.0) {
    Rcpp::stop("`zero_tolerance` must be finite and non-negative.");
  }

  const Ch3ppPairSigns sx = ch3pp_pair_signs(x, zero_tolerance);
  const Ch3ppPairSigns sy = ch3pp_pair_signs(y, zero_tolerance);
  long double within_x_unordered = 0.0L;
  long double within_y_unordered = 0.0L;
  long double cross_unordered = 0.0L;

  for (arma::uword a = 0; a < sx.signs.n_rows; ++a) {
    for (arma::uword b = a + 1; b < sx.signs.n_rows; ++b) {
      if (sx.first[a] == sx.first[b] || sx.first[a] == sx.second[b] ||
          sx.second[a] == sx.first[b] || sx.second[a] == sx.second[b]) {
        continue;
      }
      const long double dot = ch3pp_dot_ld(sx.signs.row(a), sx.signs.row(b));
      within_x_unordered += dot * dot;
    }
  }
  for (arma::uword a = 0; a < sy.signs.n_rows; ++a) {
    for (arma::uword b = a + 1; b < sy.signs.n_rows; ++b) {
      if (sy.first[a] == sy.first[b] || sy.first[a] == sy.second[b] ||
          sy.second[a] == sy.first[b] || sy.second[a] == sy.second[b]) {
        continue;
      }
      const long double dot = ch3pp_dot_ld(sy.signs.row(a), sy.signs.row(b));
      within_y_unordered += dot * dot;
    }
  }
  for (arma::uword a = 0; a < sx.signs.n_rows; ++a) {
    for (arma::uword b = 0; b < sy.signs.n_rows; ++b) {
      const long double dot = ch3pp_dot_ld(sx.signs.row(a), sy.signs.row(b));
      cross_unordered += dot * dot;
    }
  }

  // Each unordered pair of disjoint edges has 2*2*2=8 ordered
  // representations; a cross pair of unordered edges has 2*2=4.
  const long double ordered_x = 8.0L * within_x_unordered;
  const long double ordered_y = 8.0L * within_y_unordered;
  const long double ordered_cross = 4.0L * cross_unordered;
  const long double n1 = static_cast<long double>(x.n_rows);
  const long double n2 = static_cast<long double>(y.n_rows);
  const long double denom_x = n1 * (n1 - 1.0L) *
    (n1 - 2.0L) * (n1 - 3.0L);
  const long double denom_y = n2 * (n2 - 1.0L) *
    (n2 - 2.0L) * (n2 - 3.0L);
  const long double denom_cross = n1 * (n1 - 1.0L) *
    n2 * (n2 - 1.0L);

  return Rcpp::List::create(
    Rcpp::Named("A1_raw") = ch3pp_checked_double(
      ordered_x / denom_x, "The first fourth-order U-statistic"
    ),
    Rcpp::Named("A2_raw") = ch3pp_checked_double(
      ordered_y / denom_y, "The second fourth-order U-statistic"
    ),
    Rcpp::Named("C12_raw") = ch3pp_checked_double(
      ordered_cross / denom_cross, "The cross spatial-rank statistic"
    ),
    Rcpp::Named("ordered_sum_x") = ch3pp_checked_double(
      ordered_x, "The first all-distinct ordered sum"
    ),
    Rcpp::Named("ordered_sum_y") = ch3pp_checked_double(
      ordered_y, "The second all-distinct ordered sum"
    ),
    Rcpp::Named("ordered_cross_sum") = ch3pp_checked_double(
      ordered_cross, "The ordered cross sum"
    ),
    Rcpp::Named("pair_signs_x") = sx.signs,
    Rcpp::Named("pair_signs_y") = sy.signs,
    Rcpp::Named("zero_pairs_x") = static_cast<double>(sx.zero_pairs),
    Rcpp::Named("zero_pairs_y") = static_cast<double>(sy.zero_pairs)
  );
}


// [[Rcpp::export]]
Rcpp::List cpp_ch3pp_sclime(const arma::mat& a,
                            const double lambda,
                            const double tolerance,
                            const int max_iterations) {
  ch3pp_validate_finite_matrix(a, "a");
  if (a.n_rows == 0 || a.n_rows != a.n_cols) {
    Rcpp::stop("`a` must be a non-empty square matrix.");
  }
  if (!arma::approx_equal(a, a.t(), "absdiff", 1e-12)) {
    Rcpp::stop("`a` must be symmetric.");
  }
  if (!std::isfinite(lambda) || lambda < 0.0 ||
      !std::isfinite(tolerance) || tolerance <= 0.0 ||
      max_iterations < 1) {
    Rcpp::stop("Invalid SCLIME solver controls.");
  }

  arma::vec eigenvalues;
  if (!arma::eig_sym(eigenvalues, a)) {
    Rcpp::stop("The SCLIME operator norm could not be evaluated.");
  }
  const double operator_norm = arma::abs(eigenvalues).max();
  if (!std::isfinite(operator_norm) || operator_norm <= 0.0) {
    Rcpp::stop("The SCLIME coefficient matrix has zero or non-finite norm.");
  }

  const arma::uword p = a.n_cols;
  arma::mat raw_solution(p, p, arma::fill::zeros);
  arma::mat dual_solution(p, p, arma::fill::zeros);
  Rcpp::IntegerVector iterations(p);
  Rcpp::LogicalVector converged(p);
  Rcpp::NumericVector relative_update(p);
  Rcpp::NumericVector primal_violation(p);
  Rcpp::NumericVector stationarity_residual(p);
  Rcpp::NumericVector dual_violation(p);
  Rcpp::NumericVector primal_objective(p);
  Rcpp::NumericVector dual_objective(p);
  Rcpp::NumericVector duality_gap(p);
  Rcpp::NumericVector relative_gap(p);

  for (arma::uword j = 0; j < p; ++j) {
    arma::vec target(p, arma::fill::zeros);
    target(j) = 1.0;
    const Ch3ppSclimeColumn fit = ch3pp_sclime_column(
      a, target, lambda, tolerance, max_iterations, operator_norm
    );
    raw_solution.col(j) = fit.solution;
    dual_solution.col(j) = fit.dual;
    iterations[j] = fit.iterations;
    converged[j] = fit.converged;
    relative_update[j] = fit.relative_update;
    primal_violation[j] = fit.primal_violation;
    stationarity_residual[j] = fit.stationarity_residual;
    dual_violation[j] = fit.dual_violation;
    primal_objective[j] = fit.primal_objective;
    dual_objective[j] = fit.dual_objective;
    duality_gap[j] = fit.duality_gap;
    relative_gap[j] = fit.relative_gap;
  }

  arma::mat symmetrized(p, p, arma::fill::zeros);
  for (arma::uword i = 0; i < p; ++i) {
    for (arma::uword j = i; j < p; ++j) {
      const double value =
        std::abs(raw_solution(i, j)) <= std::abs(raw_solution(j, i)) ?
        raw_solution(i, j) : raw_solution(j, i);
      symmetrized(i, j) = value;
      symmetrized(j, i) = value;
    }
  }
  const double final_feasibility = std::max(
    0.0, ch3pp_max_abs_mat(
      arma::mat(a * symmetrized - arma::eye<arma::mat>(p, p))
    ) - lambda
  );
  arma::vec final_eigenvalues;
  bool eigen_ok = arma::eig_sym(final_eigenvalues, symmetrized);

  return Rcpp::List::create(
    Rcpp::Named("raw_solution") = raw_solution,
    Rcpp::Named("solution") = symmetrized,
    Rcpp::Named("dual_solution") = dual_solution,
    Rcpp::Named("iterations") = iterations,
    Rcpp::Named("converged") = converged,
    Rcpp::Named("all_converged") = Rcpp::is_true(Rcpp::all(converged)),
    Rcpp::Named("relative_update") = relative_update,
    Rcpp::Named("primal_violation") = primal_violation,
    Rcpp::Named("stationarity_residual") = stationarity_residual,
    Rcpp::Named("dual_violation") = dual_violation,
    Rcpp::Named("primal_objective") = primal_objective,
    Rcpp::Named("dual_objective") = dual_objective,
    Rcpp::Named("duality_gap") = duality_gap,
    Rcpp::Named("relative_gap") = relative_gap,
    Rcpp::Named("operator_norm") = operator_norm,
    Rcpp::Named("primal_dual_step") = 0.99 / operator_norm,
    Rcpp::Named("final_feasibility_violation") = final_feasibility,
    Rcpp::Named("final_minimum_eigenvalue") = eigen_ok ?
      final_eigenvalues.min() : NA_REAL
  );
}


// [[Rcpp::export]]
Rcpp::List cpp_ch3pp_sglasso(const arma::mat& a,
                             const double lambda,
                             const double tolerance,
                             const int max_iterations,
                             const double initial_step,
                             const int max_backtracking) {
  ch3pp_validate_finite_matrix(a, "a");
  if (a.n_rows == 0 || a.n_rows != a.n_cols) {
    Rcpp::stop("`a` must be a non-empty square matrix.");
  }
  if (!arma::approx_equal(a, a.t(), "absdiff", 1e-12)) {
    Rcpp::stop("`a` must be symmetric.");
  }
  if (!std::isfinite(lambda) || lambda <= 0.0 ||
      !std::isfinite(tolerance) || tolerance <= 0.0 ||
      max_iterations < 1 || !std::isfinite(initial_step) ||
      initial_step <= 0.0 || max_backtracking < 1) {
    Rcpp::stop("Invalid SGLASSO solver controls.");
  }

  const arma::uword p = a.n_cols;
  arma::mat current = arma::eye<arma::mat>(p, p);
  double current_objective = ch3pp_sglasso_objective(a, current, lambda);
  double step = initial_step;
  bool converged = false;
  bool backtracking_failed = false;
  int iterations = 0;
  int total_backtracking = 0;
  double relative_update = std::numeric_limits<double>::infinity();
  double kkt_residual = std::numeric_limits<double>::infinity();
  double minimum_eigenvalue = 1.0;

  for (int iteration = 1; iteration <= max_iterations; ++iteration) {
    arma::mat inverse;
    if (!arma::inv_sympd(inverse, current)) {
      Rcpp::stop(
        "The SGLASSO iterate lost positive definiteness; no ridge or "
        "eigenvalue floor is applied."
      );
    }
    const arma::mat gradient = a - inverse;
    const double smooth_current = arma::trace(a * current) -
      ch3pp_logdet_spd(current);
    bool accepted = false;
    arma::mat candidate;
    double candidate_objective = std::numeric_limits<double>::infinity();
    double candidate_step = std::min(step * 1.2, initial_step);

    for (int bt = 0; bt < max_backtracking; ++bt) {
      const arma::mat proximal_input = current - candidate_step * gradient;
      candidate = ch3pp_soft_threshold_mat(
        proximal_input, candidate_step * lambda
      );
      candidate = 0.5 * (candidate + candidate.t());
      arma::vec eigenvalues;
      const bool eigen_ok = arma::eig_sym(eigenvalues, candidate);
      if (eigen_ok && eigenvalues.min() > 0.0) {
        const arma::mat difference = candidate - current;
        const double smooth_candidate = arma::trace(a * candidate) -
          ch3pp_logdet_spd(candidate);
        const double majorizer = smooth_current +
          arma::accu(gradient % difference) +
          arma::accu(arma::square(difference)) / (2.0 * candidate_step);
        const double slack = 64.0 * std::numeric_limits<double>::epsilon() *
          std::max(1.0, std::abs(smooth_current));
        if (std::isfinite(smooth_candidate) &&
            smooth_candidate <= majorizer + slack) {
          candidate_objective = smooth_candidate +
            lambda * arma::accu(arma::abs(candidate));
          minimum_eigenvalue = eigenvalues.min();
          accepted = true;
          total_backtracking += bt;
          break;
        }
      }
      candidate_step *= 0.5;
    }

    if (!accepted) {
      backtracking_failed = true;
      iterations = iteration - 1;
      break;
    }

    const arma::mat iterate_difference = candidate - current;
    relative_update = ch3pp_max_abs_mat(iterate_difference) /
      std::max(1.0, ch3pp_max_abs_mat(current));
    current = candidate;
    current_objective = candidate_objective;
    step = candidate_step;
    arma::mat inverse_candidate;
    kkt_residual = ch3pp_sglasso_kkt(
      a, current, lambda, inverse_candidate
    );
    iterations = iteration;
    if (std::isfinite(kkt_residual) && kkt_residual <= tolerance) {
      converged = true;
      break;
    }
  }

  arma::mat inverse;
  const bool inverse_ok = arma::inv_sympd(inverse, current);
  if (inverse_ok) {
    kkt_residual = ch3pp_sglasso_kkt(a, current, lambda, inverse);
  }
  return Rcpp::List::create(
    Rcpp::Named("solution") = current,
    Rcpp::Named("inverse") = inverse_ok ? inverse :
      arma::mat(p, p, arma::fill::value(NA_REAL)),
    Rcpp::Named("objective") = current_objective,
    Rcpp::Named("iterations") = iterations,
    Rcpp::Named("converged") = converged,
    Rcpp::Named("backtracking_failed") = backtracking_failed,
    Rcpp::Named("relative_update") = relative_update,
    Rcpp::Named("kkt_residual") = kkt_residual,
    Rcpp::Named("minimum_eigenvalue") = minimum_eigenvalue,
    Rcpp::Named("accepted_step") = step,
    Rcpp::Named("total_backtracking") = total_backtracking
  );
}
