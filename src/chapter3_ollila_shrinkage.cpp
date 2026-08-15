// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp17)]]

#include <RcppArmadillo.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include <string>
#include <vector>

namespace {

struct OsGaussRule {
  std::vector<double> nodes;
  std::vector<double> weights;
};

struct OsMapValue {
  double value;
  double error;
};

void os_require_finite(const arma::mat& x, const char* name) {
  if (!x.is_finite()) {
    Rcpp::stop("`%s` must contain only finite values.", name);
  }
}

bool os_unit_row(const arma::rowvec& value, arma::rowvec& direction,
                 double& norm, double& log_norm) {
  const double scale = arma::abs(value).max();
  if (scale == 0.0) {
    direction.zeros(value.n_elem);
    norm = 0.0;
    log_norm = R_NegInf;
    return false;
  }
  const arma::rowvec scaled = value / scale;
  const double scaled_norm = arma::norm(scaled, 2);
  direction = scaled / scaled_norm;
  norm = scale * scaled_norm;
  log_norm = std::log(scale) + std::log(scaled_norm);
  return true;
}

OsGaussRule os_gauss_legendre_rule(const int order) {
  OsGaussRule rule;
  rule.nodes.resize(order);
  rule.weights.resize(order);
  const int half = (order + 1) / 2;
  const long double pi = std::acos(-1.0L);

  for (int i = 0; i < half; ++i) {
    long double z = std::cos(
      pi * (static_cast<long double>(i) + 0.75L) /
      (static_cast<long double>(order) + 0.5L)
    );
    long double previous = 0.0L;
    long double derivative = 0.0L;
    for (int iteration = 0; iteration < 100; ++iteration) {
      long double p1 = 1.0L;
      long double p2 = 0.0L;
      for (int j = 1; j <= order; ++j) {
        const long double p3 = p2;
        p2 = p1;
        p1 = (
          (2.0L * static_cast<long double>(j) - 1.0L) * z * p2 -
          (static_cast<long double>(j) - 1.0L) * p3
        ) / static_cast<long double>(j);
      }
      derivative = static_cast<long double>(order) * (z * p1 - p2) /
        (z * z - 1.0L);
      previous = z;
      z = previous - p1 / derivative;
      if (std::abs(z - previous) <= 4.0L *
          std::numeric_limits<long double>::epsilon()) {
        break;
      }
    }
    const long double weight = 2.0L /
      ((1.0L - z * z) * derivative * derivative);
    const double left = static_cast<double>((1.0L - z) / 2.0L);
    const double right = static_cast<double>((1.0L + z) / 2.0L);
    const double mapped_weight = static_cast<double>(weight / 2.0L);
    rule.nodes[i] = left;
    rule.nodes[order - 1 - i] = right;
    rule.weights[i] = mapped_weight;
    rule.weights[order - 1 - i] = mapped_weight;
  }
  return rule;
}

double os_basic_map_rule(const double lambda, const int dimension,
                         const OsGaussRule& rule) {
  if (lambda == 0.0) {
    return 0.0;
  }
  long double total = 0.0L;
  const long double lam = static_cast<long double>(lambda);
  for (std::size_t i = 0; i < rule.nodes.size(); ++i) {
    const long double v = static_cast<long double>(rule.nodes[i]);
    const long double squared = v * v;
    const long double radial_weight = std::pow(
      v, static_cast<long double>(dimension - 1)
    );
    long double integrand = 0.0L;
    if (lam >= 1.0L) {
      const long double denominator =
        (1.0L - squared) + squared / lam;
      integrand = radial_weight / denominator;
    } else {
      const long double denominator =
        squared + (1.0L - squared) * lam;
      integrand = lam * radial_weight / denominator;
    }
    total += static_cast<long double>(rule.weights[i]) * integrand;
  }
  return static_cast<double>(total);
}

OsMapValue os_basic_map(const double lambda, const int dimension,
                        const OsGaussRule& coarse,
                        const OsGaussRule& fine) {
  const double coarse_value = os_basic_map_rule(
    lambda, dimension, coarse
  );
  const double fine_value = os_basic_map_rule(lambda, dimension, fine);
  return OsMapValue{fine_value, std::abs(fine_value - coarse_value)};
}

arma::vec os_project_lower(const arma::vec& value,
                           const arma::vec& lower) {
  return arma::max(value, lower);
}

arma::vec os_project_simplex_lower(const arma::vec& value,
                                   const arma::vec& lower,
                                   const double total) {
  const double remaining = total - arma::accu(lower);
  if (remaining < 0.0) {
    Rcpp::stop("The lower bounds are infeasible for the convex constraint.");
  }
  if (remaining == 0.0) {
    return lower;
  }
  const arma::vec shifted = value - lower;
  arma::vec sorted = arma::sort(shifted, "descend");
  double cumulative = 0.0;
  double theta = 0.0;
  int rho = -1;
  for (arma::uword j = 0; j < sorted.n_elem; ++j) {
    cumulative += sorted(j);
    const double candidate = (
      cumulative - remaining
    ) / static_cast<double>(j + 1);
    if (sorted(j) - candidate > 0.0) {
      rho = static_cast<int>(j);
      theta = candidate;
    }
  }
  if (rho < 0) {
    theta = (arma::accu(sorted) - remaining) /
      static_cast<double>(sorted.n_elem);
  }
  return lower + arma::clamp(shifted - theta, 0.0, arma::datum::inf);
}

double os_qp_objective(const arma::mat& matrix, const arma::vec& linear,
                       const arma::vec& value) {
  return 0.5 * arma::as_scalar(value.t() * matrix * value) -
    arma::dot(linear, value);
}

Rcpp::List os_qp_diagnostics(const arma::mat& matrix,
                             const arma::vec& linear,
                             const arma::vec& lower,
                             const arma::vec& value,
                             const bool convex,
                             const double lipschitz,
                             const double tolerance,
                             const int iterations,
                             const bool converged,
                             const double objective_scale,
                             const double minimum_eigenvalue,
                             const int restarts) {
  const arma::vec gradient = matrix * value - linear;
  const double lower_violation = std::max(
    0.0, arma::max(lower - value)
  );
  const double equality_violation = convex ?
    std::abs(arma::accu(value) - 1.0) : 0.0;
  const double active_tolerance = 10.0 * tolerance *
    (1.0 + arma::abs(value).max());
  double kkt_residual = 0.0;
  double complementarity_gap = 0.0;
  double equality_gradient_level = NA_REAL;

  if (convex) {
    std::vector<arma::uword> free_indices;
    for (arma::uword i = 0; i < value.n_elem; ++i) {
      if (value(i) - lower(i) > active_tolerance) {
        free_indices.push_back(i);
      }
    }
    double level = gradient.min();
    if (!free_indices.empty()) {
      long double sum = 0.0L;
      for (const arma::uword index : free_indices) {
        sum += gradient(index);
      }
      level = static_cast<double>(
        sum / static_cast<long double>(free_indices.size())
      );
    }
    equality_gradient_level = level;
    for (arma::uword i = 0; i < value.n_elem; ++i) {
      const double slack = value(i) - lower(i);
      if (slack > active_tolerance) {
        kkt_residual = std::max(
          kkt_residual, std::abs(gradient(i) - level)
        );
      } else {
        kkt_residual = std::max(
          kkt_residual, std::max(0.0, level - gradient(i))
        );
      }
      complementarity_gap += std::abs(slack * (gradient(i) - level));
    }
  } else {
    for (arma::uword i = 0; i < value.n_elem; ++i) {
      const double slack = value(i) - lower(i);
      if (slack > active_tolerance) {
        kkt_residual = std::max(kkt_residual, std::abs(gradient(i)));
      } else {
        kkt_residual = std::max(
          kkt_residual, std::max(0.0, -gradient(i))
        );
      }
      complementarity_gap += std::abs(slack * gradient(i));
    }
  }

  arma::vec projected = value;
  if (lipschitz > 0.0) {
    const arma::vec trial = value - gradient / lipschitz;
    projected = convex ? os_project_simplex_lower(trial, lower, 1.0) :
      os_project_lower(trial, lower);
  }
  const double projected_gradient_residual = lipschitz > 0.0 ?
    lipschitz * arma::abs(projected - value).max() : kkt_residual;

  return Rcpp::List::create(
    Rcpp::Named("iterations") = iterations,
    Rcpp::Named("converged") = converged,
    Rcpp::Named("objective") =
      os_qp_objective(matrix, linear, value) * objective_scale,
    Rcpp::Named("objective_scaled") =
      os_qp_objective(matrix, linear, value),
    Rcpp::Named("objective_scale") = objective_scale,
    Rcpp::Named("feasibility_residual") =
      std::max(lower_violation, equality_violation),
    Rcpp::Named("lower_bound_violation") = lower_violation,
    Rcpp::Named("equality_violation") = equality_violation,
    Rcpp::Named("kkt_residual_scaled") = kkt_residual,
    Rcpp::Named("kkt_residual") = kkt_residual * objective_scale,
    Rcpp::Named("projected_gradient_residual_scaled") =
      projected_gradient_residual,
    Rcpp::Named("projected_gradient_residual") =
      projected_gradient_residual * objective_scale,
    Rcpp::Named("complementarity_gap_scaled") = complementarity_gap,
    Rcpp::Named("complementarity_gap") =
      complementarity_gap * objective_scale,
    Rcpp::Named("equality_gradient_level_scaled") =
      equality_gradient_level,
    Rcpp::Named("minimum_eigenvalue_scaled") = minimum_eigenvalue,
    Rcpp::Named("lipschitz_constant_scaled") = lipschitz,
    Rcpp::Named("restarts") = restarts,
    Rcpp::Named("tolerance") = tolerance,
    Rcpp::Named("regularization") = "none",
    Rcpp::Named("pseudoinverse") = "none",
    Rcpp::Named("ridge") = "none"
  );
}

} // namespace


// [[Rcpp::export]]
Rcpp::List cpp_ollila_sample_moments(const arma::mat& x) {
  os_require_finite(x, "x");
  if (x.n_rows < 2 || x.n_cols < 1) {
    Rcpp::stop("`x` must have at least two rows and one column.");
  }
  double data_scale = arma::abs(x).max();
  if (data_scale == 0.0) {
    data_scale = 1.0;
  }
  const arma::mat scaled = x / data_scale;
  const arma::rowvec mean_scaled = arma::mean(scaled, 0);
  const arma::mat centered = scaled.each_row() - mean_scaled;
  const arma::mat covariance_scaled = centered.t() * centered /
    static_cast<double>(x.n_rows - 1);
  arma::vec g2(x.n_cols, arma::fill::value(NA_REAL));
  arma::uvec zero_variance(x.n_cols, arma::fill::zeros);

  for (arma::uword j = 0; j < x.n_cols; ++j) {
    const arma::vec column = centered.col(j);
    const double column_scale = arma::abs(column).max();
    if (column_scale == 0.0) {
      zero_variance(j) = 1;
      continue;
    }
    const arma::vec normalized = column / column_scale;
    const double second = arma::mean(arma::square(normalized));
    const double fourth = arma::mean(arma::square(arma::square(normalized)));
    if (!(second > 0.0) || !std::isfinite(second) ||
        !std::isfinite(fourth)) {
      zero_variance(j) = 1;
      continue;
    }
    g2(j) = fourth / (second * second) - 3.0;
  }

  return Rcpp::List::create(
    Rcpp::Named("data_scale") = data_scale,
    Rcpp::Named("mean_scaled") = mean_scaled,
    Rcpp::Named("covariance_scaled") = covariance_scaled,
    Rcpp::Named("g2") = g2,
    Rcpp::Named("zero_variance") = zero_variance,
    Rcpp::Named("trace_covariance_scaled") =
      arma::trace(covariance_scaled),
    Rcpp::Named("trace_covariance_squared_scaled") =
      arma::accu(arma::square(covariance_scaled))
  );
}


// [[Rcpp::export]]
Rcpp::List cpp_ollila_sign_components(const arma::mat& residuals,
                                      const bool keep_signs) {
  os_require_finite(residuals, "residuals");
  if (residuals.n_rows < 2 || residuals.n_cols < 1) {
    Rcpp::stop("`residuals` must have at least two rows and one column.");
  }
  const arma::uword n = residuals.n_rows;
  const arma::uword p = residuals.n_cols;
  arma::mat signs(n, p, arma::fill::zeros);
  arma::vec radii(n, arma::fill::zeros);
  arma::vec log_radii(n, arma::fill::value(R_NegInf));
  arma::uword zero_count = 0;
  for (arma::uword i = 0; i < n; ++i) {
    arma::rowvec direction(p, arma::fill::zeros);
    double norm = 0.0;
    double log_norm = R_NegInf;
    if (os_unit_row(residuals.row(i), direction, norm, log_norm)) {
      signs.row(i) = direction;
      radii(i) = norm;
      log_radii(i) = log_norm;
    } else {
      ++zero_count;
    }
  }
  const arma::mat sscm = signs.t() * signs / static_cast<double>(n);
  const arma::mat shape = static_cast<double>(p) * sscm;
  const double a = arma::accu(arma::square(shape)) /
    static_cast<double>(p);

  arma::mat left;
  arma::mat eigenvectors;
  arma::vec singular_values;
  if (!arma::svd(left, singular_values, eigenvectors, signs)) {
    Rcpp::stop("The singular-value decomposition of the sign matrix failed.");
  }
  arma::vec sscm_eigenvalues(p, arma::fill::zeros);
  for (arma::uword i = 0; i < singular_values.n_elem; ++i) {
    sscm_eigenvalues(i) = singular_values(i) * singular_values(i) /
      static_cast<double>(n);
  }

  Rcpp::RObject returned_signs = R_NilValue;
  if (keep_signs) {
    returned_signs = Rcpp::wrap(signs);
  }
  return Rcpp::List::create(
    Rcpp::Named("sscm") = sscm,
    Rcpp::Named("shape") = shape,
    Rcpp::Named("a") = a,
    Rcpp::Named("radii_scaled") = radii,
    Rcpp::Named("log_radii_scaled") = log_radii,
    Rcpp::Named("n_zero") = static_cast<double>(zero_count),
    Rcpp::Named("signs") = returned_signs,
    Rcpp::Named("eigenvectors") = eigenvectors,
    Rcpp::Named("sscm_eigenvalues") = sscm_eigenvalues,
    Rcpp::Named("trace_sscm") = arma::trace(sscm)
  );
}


// [[Rcpp::export]]
Rcpp::List cpp_ollila_basic_inverse(
    const arma::vec& delta,
    const int dimension,
    const int quadrature_order,
    const double inversion_tolerance,
    const double integration_tolerance,
    const int max_iterations,
    const int max_bracket_iterations) {
  if (dimension < 2 || delta.n_elem != static_cast<arma::uword>(dimension)) {
    Rcpp::stop("`delta` must have length `dimension`, with dimension >= 2.");
  }
  if (!delta.is_finite() || arma::any(delta < 0.0) ||
      arma::any(delta > 1.0)) {
    Rcpp::stop("BASIC target eigenvalues must lie in [0, 1].");
  }
  if (quadrature_order < 8 || quadrature_order > 2048 ||
      !std::isfinite(inversion_tolerance) || inversion_tolerance <= 0.0 ||
      !std::isfinite(integration_tolerance) || integration_tolerance <= 0.0 ||
      max_iterations < 1 || max_bracket_iterations < 1) {
    Rcpp::stop("Invalid BASIC quadrature or inversion controls.");
  }
  const OsGaussRule coarse = os_gauss_legendre_rule(quadrature_order);
  const OsGaussRule fine = os_gauss_legendre_rule(2 * quadrature_order);
  arma::vec lambda(dimension, arma::fill::zeros);
  arma::vec mapped(dimension, arma::fill::zeros);
  arma::vec integration_error(dimension, arma::fill::zeros);
  arma::vec inversion_error(dimension, arma::fill::zeros);
  arma::mat brackets(dimension, 2, arma::fill::zeros);
  arma::ivec iterations(dimension, arma::fill::zeros);
  arma::uvec boundary_zero(dimension, arma::fill::zeros);

  for (int index = 0; index < dimension; ++index) {
    const double target = delta(index);
    if (target == 0.0) {
      boundary_zero(index) = 1;
      continue;
    }
    double lower = 0.0;
    double upper = std::max(1.0, static_cast<double>(dimension));
    OsMapValue upper_map = os_basic_map(
      upper, dimension, coarse, fine
    );
    int bracket_iterations = 0;
    while (upper_map.value < target &&
           bracket_iterations < max_bracket_iterations) {
      if (upper > std::numeric_limits<double>::max() / 2.0) {
        break;
      }
      upper *= 2.0;
      upper_map = os_basic_map(upper, dimension, coarse, fine);
      ++bracket_iterations;
    }
    if (!std::isfinite(upper_map.value) || upper_map.value < target) {
      Rcpp::stop(
        "BASIC inversion could not bracket eigenvalue %d (target %.17g).",
        index + 1, target
      );
    }

    bool converged = false;
    double midpoint = upper / 2.0;
    OsMapValue midpoint_map = os_basic_map(
      midpoint, dimension, coarse, fine
    );
    int iteration = 0;
    for (iteration = 1; iteration <= max_iterations; ++iteration) {
      midpoint = lower + (upper - lower) / 2.0;
      midpoint_map = os_basic_map(
        midpoint, dimension, coarse, fine
      );
      if (midpoint_map.value < target) {
        lower = midpoint;
      } else {
        upper = midpoint;
      }
      const double residual = std::abs(midpoint_map.value - target);
      const double bracket_width = upper - lower;
      if (residual <= inversion_tolerance * std::max(1.0, target) &&
          bracket_width <= inversion_tolerance * (1.0 + midpoint)) {
        converged = true;
        break;
      }
    }
    midpoint = lower + (upper - lower) / 2.0;
    midpoint_map = os_basic_map(midpoint, dimension, coarse, fine);
    const double residual = std::abs(midpoint_map.value - target);
    if (!converged || residual > inversion_tolerance * std::max(1.0, target)) {
      Rcpp::stop(
        "BASIC inversion did not converge for eigenvalue %d; residual %.6g.",
        index + 1, residual
      );
    }
    if (midpoint_map.error > integration_tolerance) {
      Rcpp::stop(
        "BASIC quadrature error %.6g exceeds tolerance for eigenvalue %d.",
        midpoint_map.error, index + 1
      );
    }
    lambda(index) = midpoint;
    mapped(index) = midpoint_map.value;
    integration_error(index) = midpoint_map.error;
    inversion_error(index) = residual;
    brackets(index, 0) = lower;
    brackets(index, 1) = upper;
    iterations(index) = iteration;
  }

  const double lambda_sum = arma::accu(lambda);
  if (!std::isfinite(lambda_sum) || lambda_sum <= 0.0) {
    Rcpp::stop("BASIC inversion produced no positive shape eigenvalue.");
  }
  const arma::vec normalized = static_cast<double>(dimension) *
    lambda / lambda_sum;
  return Rcpp::List::create(
    Rcpp::Named("lambda_raw") = lambda,
    Rcpp::Named("lambda_normalized") = normalized,
    Rcpp::Named("delta_target") = delta,
    Rcpp::Named("delta_mapped") = mapped,
    Rcpp::Named("integration_error") = integration_error,
    Rcpp::Named("inversion_error") = inversion_error,
    Rcpp::Named("brackets") = brackets,
    Rcpp::Named("iterations") = iterations,
    Rcpp::Named("boundary_zero") = boundary_zero,
    Rcpp::Named("quadrature_order_coarse") = quadrature_order,
    Rcpp::Named("quadrature_order_fine") = 2 * quadrature_order,
    Rcpp::Named("max_integration_error") = integration_error.max(),
    Rcpp::Named("max_inversion_error") = inversion_error.max(),
    Rcpp::Named("converged") = true,
    Rcpp::Named("extrapolation") = "none"
  );
}


// [[Rcpp::export]]
Rcpp::List cpp_ollila_pool_qp(
    const arma::mat& objective_matrix,
    const arma::vec& linear_term,
    const arma::vec& lower,
    const bool convex,
    const double tolerance,
    const int max_iterations) {
  os_require_finite(objective_matrix, "objective_matrix");
  if (!linear_term.is_finite() || !lower.is_finite()) {
    Rcpp::stop("QP vectors must contain only finite values.");
  }
  const arma::uword dimension = objective_matrix.n_rows;
  if (dimension < 1 || objective_matrix.n_cols != dimension ||
      linear_term.n_elem != dimension || lower.n_elem != dimension) {
    Rcpp::stop("QP dimensions are incompatible.");
  }
  if (arma::any(lower < 0.0)) {
    Rcpp::stop("QP lower bounds must be non-negative.");
  }
  if (!std::isfinite(tolerance) || tolerance <= 0.0 ||
      max_iterations < 1) {
    Rcpp::stop("Invalid QP solver controls.");
  }
  if (convex && arma::accu(lower) > 1.0) {
    Rcpp::stop("QP lower bounds are infeasible with sum(a) = 1.");
  }

  arma::mat matrix = 0.5 * (objective_matrix + objective_matrix.t());
  const double asymmetry = arma::abs(
    objective_matrix - objective_matrix.t()
  ).max();
  const double objective_scale = std::max(
    arma::abs(matrix).max(), arma::abs(linear_term).max()
  );
  if (objective_scale == 0.0) {
    arma::vec solution;
    if (convex) {
      solution = os_project_simplex_lower(lower, lower, 1.0);
    } else {
      solution = lower;
    }
    Rcpp::List diagnostics = os_qp_diagnostics(
      matrix, linear_term, lower, solution, convex, 0.0, tolerance,
      0, true, 1.0, 0.0, 0
    );
    diagnostics["objective_asymmetry"] = asymmetry;
    diagnostics["zero_objective"] = true;
    return Rcpp::List::create(
      Rcpp::Named("solution") = solution,
      Rcpp::Named("diagnostics") = diagnostics
    );
  }
  matrix /= objective_scale;
  arma::vec linear = linear_term / objective_scale;

  arma::vec eigenvalues;
  if (!arma::eig_sym(eigenvalues, matrix)) {
    Rcpp::stop("The QP objective eigendecomposition failed.");
  }
  const double minimum_eigenvalue = eigenvalues.min();
  const double lipschitz = eigenvalues.max();
  const double psd_tolerance = 100.0 * std::numeric_limits<double>::epsilon() *
    std::max(1.0, std::abs(lipschitz)) * static_cast<double>(dimension);
  if (minimum_eigenvalue < -psd_tolerance) {
    Rcpp::stop(
      "The estimated QP objective is not positive semidefinite (min eigenvalue %.6g); no ridge is applied.",
      minimum_eigenvalue * objective_scale
    );
  }

  if (!(lipschitz > 0.0)) {
    if (!convex && arma::any(linear > 0.0)) {
      Rcpp::stop("The lower-bound-only QP is unbounded for a zero Hessian.");
    }
    arma::vec solution = lower;
    if (convex) {
      const double remaining = 1.0 - arma::accu(lower);
      arma::uword best = linear.index_max();
      solution(best) += remaining;
    }
    Rcpp::List diagnostics = os_qp_diagnostics(
      matrix, linear, lower, solution, convex, 0.0, tolerance,
      0, true, objective_scale, minimum_eigenvalue, 0
    );
    diagnostics["objective_asymmetry"] = asymmetry;
    diagnostics["zero_hessian"] = true;
    return Rcpp::List::create(
      Rcpp::Named("solution") = solution,
      Rcpp::Named("diagnostics") = diagnostics
    );
  }

  auto project = [&](const arma::vec& value) {
    return convex ? os_project_simplex_lower(value, lower, 1.0) :
      os_project_lower(value, lower);
  };
  arma::vec current = project(linear / lipschitz);
  arma::vec extrapolated = current;
  double momentum = 1.0;
  bool converged = false;
  int iterations = 0;
  int restarts = 0;

  for (int iteration = 1; iteration <= max_iterations; ++iteration) {
    iterations = iteration;
    const arma::vec gradient = matrix * extrapolated - linear;
    arma::vec candidate = project(extrapolated - gradient / lipschitz);
    double candidate_objective = os_qp_objective(
      matrix, linear, candidate
    );
    const double current_objective = os_qp_objective(
      matrix, linear, current
    );
    if (candidate_objective > current_objective) {
      extrapolated = current;
      const arma::vec restart_gradient = matrix * current - linear;
      candidate = project(current - restart_gradient / lipschitz);
      candidate_objective = os_qp_objective(matrix, linear, candidate);
      momentum = 1.0;
      ++restarts;
    }
    const double next_momentum = 0.5 * (
      1.0 + std::sqrt(1.0 + 4.0 * momentum * momentum)
    );
    const arma::vec next_extrapolated = candidate +
      ((momentum - 1.0) / next_momentum) * (candidate - current);
    const double step = arma::abs(candidate - current).max();
    current = candidate;
    extrapolated = next_extrapolated;
    momentum = next_momentum;

    const arma::vec final_gradient = matrix * current - linear;
    const arma::vec projected = project(
      current - final_gradient / lipschitz
    );
    const double projected_residual = lipschitz *
      arma::abs(projected - current).max();
    if (projected_residual <= tolerance &&
        step <= tolerance * (1.0 + arma::abs(current).max())) {
      converged = true;
      break;
    }
  }

  Rcpp::List diagnostics = os_qp_diagnostics(
    matrix, linear, lower, current, convex, lipschitz, tolerance,
    iterations, converged, objective_scale, minimum_eigenvalue, restarts
  );
  diagnostics["objective_asymmetry"] = asymmetry;
  diagnostics["psd_check_tolerance_scaled"] = psd_tolerance;
  diagnostics["zero_hessian"] = false;
  return Rcpp::List::create(
    Rcpp::Named("solution") = current,
    Rcpp::Named("diagnostics") = diagnostics
  );
}
