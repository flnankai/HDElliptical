// Chapter 3: Gaussian covariance and precision benchmarks.
//
// EC2 follows the correlation-matrix formulation of Liu, Wang and Zhao:
// an off-diagonal l1 penalty, unit diagonal, and an explicit lower eigenvalue
// constraint.  Gaussian graphical lasso penalizes only off-diagonal entries.
// CLIME is solved column by column and then symmetrized by retaining the entry
// with smaller absolute value.  Every optimizer returns auditable feasibility
// and KKT or primal-dual certificates.  No ridge, eigenvalue floor other than
// EC2's method-defining tau, pseudoinverse, or post-hoc repair is used.

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp17)]]

#include <RcppArmadillo.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

namespace {

void ch3gp_validate_finite_matrix(const arma::mat& x, const char* name) {
  if (!x.is_finite()) {
    Rcpp::stop("`%s` must contain only finite values.", name);
  }
}

double ch3gp_max_abs_mat(const arma::mat& x) {
  return x.n_elem == 0 ? 0.0 : arma::abs(x).max();
}

double ch3gp_max_abs_vec(const arma::vec& x) {
  return x.n_elem == 0 ? 0.0 : arma::abs(x).max();
}

double ch3gp_symmetry_tolerance(const arma::mat& x) {
  return 1e-12 * std::max(1.0, ch3gp_max_abs_mat(x));
}

double ch3gp_soft_scalar(const double value, const double threshold) {
  if (value > threshold) return value - threshold;
  if (value < -threshold) return value + threshold;
  return 0.0;
}

arma::vec ch3gp_soft_vector(const arma::vec& x, const double threshold) {
  arma::vec answer(x.n_elem, arma::fill::zeros);
  for (arma::uword j = 0; j < x.n_elem; ++j) {
    answer(j) = ch3gp_soft_scalar(x(j), threshold);
  }
  return answer;
}

arma::mat ch3gp_soft_off_diagonal(const arma::mat& input,
                                  const double threshold,
                                  const bool unit_diagonal) {
  arma::mat answer = 0.5 * (input + input.t());
  for (arma::uword i = 0; i < answer.n_rows; ++i) {
    if (unit_diagonal) answer(i, i) = 1.0;
    for (arma::uword j = i + 1; j < answer.n_cols; ++j) {
      const double value = ch3gp_soft_scalar(answer(i, j), threshold);
      answer(i, j) = value;
      answer(j, i) = value;
    }
  }
  return answer;
}

bool ch3gp_eigenvalues(const arma::mat& x, arma::vec& values) {
  return arma::eig_sym(values, 0.5 * (x + x.t()));
}

double ch3gp_minimum_eigenvalue(const arma::mat& x) {
  arma::vec values;
  if (!ch3gp_eigenvalues(x, values) || values.n_elem == 0) {
    return NA_REAL;
  }
  return values.min();
}

arma::mat ch3gp_project_eigen_floor(const arma::mat& input,
                                    const double floor) {
  arma::vec values;
  arma::mat vectors;
  const arma::mat symmetric = 0.5 * (input + input.t());
  if (!arma::eig_sym(values, vectors, symmetric)) {
    Rcpp::stop("The EC2 spectral projection eigendecomposition failed.");
  }
  for (arma::uword j = 0; j < values.n_elem; ++j) {
    values(j) = std::max(values(j), floor);
  }
  const arma::mat projected = vectors * arma::diagmat(values) * vectors.t();
  if (!projected.is_finite()) {
    Rcpp::stop("The EC2 spectral projection produced non-finite values.");
  }
  return 0.5 * (projected + projected.t());
}

double ch3gp_off_diagonal_l1(const arma::mat& x) {
  double answer = 0.0;
  for (arma::uword i = 0; i < x.n_rows; ++i) {
    for (arma::uword j = 0; j < x.n_cols; ++j) {
      if (i != j) answer += std::abs(x(i, j));
    }
  }
  return answer;
}

double ch3gp_ec2_objective(const arma::mat& correlation,
                           const arma::mat& estimate,
                           const double lambda) {
  return 0.5 * arma::accu(arma::square(correlation - estimate)) +
    lambda * ch3gp_off_diagonal_l1(estimate);
}

struct Ch3gpEc2Certificate {
  double primal_equality;
  double sparse_map;
  double spectral_map;
  double diagonal;
  double symmetry;
  double eigenvalue_violation;
  double minimum_eigenvalue;
  double off_diagonal_stationarity;
  double spectral_dual_violation;
  double spectral_complementarity;
  double maximum;
  arma::mat spectral_multiplier;
};

Ch3gpEc2Certificate ch3gp_ec2_certificate(
    const arma::mat& correlation,
    const arma::mat& sparse,
    const arma::mat& spectral,
    const arma::mat& dual,
    const double lambda,
    const double tau,
    const double rho) {
  Ch3gpEc2Certificate certificate;
  const arma::uword p = sparse.n_rows;
  const arma::mat sparse_map = ch3gp_soft_off_diagonal(
    spectral + dual / rho, lambda / rho, true
  );
  const arma::mat spectral_map = ch3gp_project_eigen_floor(
    (correlation + rho * sparse - dual) / (1.0 + rho), tau
  );
  certificate.primal_equality = ch3gp_max_abs_mat(spectral - sparse);
  certificate.sparse_map = ch3gp_max_abs_mat(sparse_map - sparse);
  certificate.spectral_map = ch3gp_max_abs_mat(spectral_map - spectral);
  certificate.diagonal = 0.0;
  for (arma::uword i = 0; i < p; ++i) {
    certificate.diagonal = std::max(
      certificate.diagonal, std::abs(sparse(i, i) - 1.0)
    );
  }
  certificate.symmetry = ch3gp_max_abs_mat(sparse - sparse.t());
  certificate.minimum_eigenvalue = ch3gp_minimum_eigenvalue(sparse);
  certificate.eigenvalue_violation =
    std::isfinite(certificate.minimum_eigenvalue) ?
    std::max(0.0, tau - certificate.minimum_eigenvalue) :
    std::numeric_limits<double>::infinity();

  certificate.off_diagonal_stationarity = 0.0;
  for (arma::uword i = 0; i < p; ++i) {
    for (arma::uword j = i + 1; j < p; ++j) {
      double residual = 0.0;
      if (sparse(i, j) > 0.0) {
        residual = std::abs(dual(i, j) - lambda);
      } else if (sparse(i, j) < 0.0) {
        residual = std::abs(dual(i, j) + lambda);
      } else {
        residual = std::max(0.0, std::abs(dual(i, j)) - lambda);
      }
      certificate.off_diagonal_stationarity = std::max(
        certificate.off_diagonal_stationarity, residual
      );
    }
  }

  certificate.spectral_multiplier = 0.5 * (
    sparse - correlation + dual +
    (sparse - correlation + dual).t()
  );
  const double multiplier_minimum = ch3gp_minimum_eigenvalue(
    certificate.spectral_multiplier
  );
  certificate.spectral_dual_violation =
    std::isfinite(multiplier_minimum) ?
    std::max(0.0, -multiplier_minimum) :
    std::numeric_limits<double>::infinity();
  const arma::mat slack = sparse - tau * arma::eye<arma::mat>(p, p);
  certificate.spectral_complementarity = ch3gp_max_abs_mat(
    certificate.spectral_multiplier * slack
  );

  certificate.maximum = std::max({
    certificate.primal_equality,
    certificate.sparse_map,
    certificate.spectral_map,
    certificate.diagonal,
    certificate.symmetry,
    certificate.eigenvalue_violation,
    certificate.off_diagonal_stationarity,
    certificate.spectral_dual_violation,
    certificate.spectral_complementarity
  });
  return certificate;
}

Rcpp::List ch3gp_ec2_result(const arma::mat& correlation,
                            const arma::mat& sparse,
                            const arma::mat& spectral,
                            const arma::mat& dual,
                            const double lambda,
                            const double tau,
                            const double rho,
                            const int iterations,
                            const bool converged,
                            const double relative_update,
                            const double dual_residual,
                            const bool sto_feasible,
                            const std::vector<double>& history) {
  const Ch3gpEc2Certificate certificate = ch3gp_ec2_certificate(
    correlation, sparse, spectral, dual, lambda, tau, rho
  );
  return Rcpp::List::create(
    Rcpp::Named("solution") = sparse,
    Rcpp::Named("spectral_solution") = spectral,
    Rcpp::Named("dual") = dual,
    Rcpp::Named("spectral_multiplier") = certificate.spectral_multiplier,
    Rcpp::Named("objective") = ch3gp_ec2_objective(
      correlation, sparse, lambda
    ),
    Rcpp::Named("objective_history") = Rcpp::wrap(history),
    Rcpp::Named("iterations") = iterations,
    Rcpp::Named("converged") = converged,
    Rcpp::Named("relative_update") = relative_update,
    Rcpp::Named("dual_residual") = dual_residual,
    Rcpp::Named("sto_feasible") = sto_feasible,
    Rcpp::Named("primal_equality_residual") =
      certificate.primal_equality,
    Rcpp::Named("sparse_fixed_point_residual") = certificate.sparse_map,
    Rcpp::Named("spectral_fixed_point_residual") = certificate.spectral_map,
    Rcpp::Named("diagonal_residual") = certificate.diagonal,
    Rcpp::Named("symmetry_residual") = certificate.symmetry,
    Rcpp::Named("minimum_eigenvalue") = certificate.minimum_eigenvalue,
    Rcpp::Named("eigenvalue_feasibility_violation") =
      certificate.eigenvalue_violation,
    Rcpp::Named("off_diagonal_stationarity_residual") =
      certificate.off_diagonal_stationarity,
    Rcpp::Named("spectral_dual_violation") =
      certificate.spectral_dual_violation,
    Rcpp::Named("spectral_complementarity_residual") =
      certificate.spectral_complementarity,
    Rcpp::Named("kkt_maximum") = certificate.maximum
  );
}

double ch3gp_logdet_spd(const arma::mat& x) {
  double value = 0.0;
  double sign = 0.0;
  if (!arma::log_det(value, sign, x) || sign <= 0.0 ||
      !std::isfinite(value)) {
    return std::numeric_limits<double>::quiet_NaN();
  }
  return value;
}

double ch3gp_glasso_objective(const arma::mat& scatter,
                              const arma::mat& precision,
                              const double lambda) {
  const double logdet = ch3gp_logdet_spd(precision);
  if (!std::isfinite(logdet)) {
    return std::numeric_limits<double>::infinity();
  }
  return arma::trace(scatter * precision) - logdet +
    lambda * ch3gp_off_diagonal_l1(precision);
}

struct Ch3gpGlassoKkt {
  arma::mat inverse;
  double diagonal;
  double off_diagonal;
  double maximum;
};

Ch3gpGlassoKkt ch3gp_glasso_kkt(const arma::mat& scatter,
                                const arma::mat& precision,
                                const double lambda) {
  Ch3gpGlassoKkt result;
  result.inverse.set_size(precision.n_rows, precision.n_cols);
  if (!arma::inv_sympd(result.inverse, precision)) {
    result.inverse.fill(NA_REAL);
    result.diagonal = std::numeric_limits<double>::infinity();
    result.off_diagonal = std::numeric_limits<double>::infinity();
    result.maximum = std::numeric_limits<double>::infinity();
    return result;
  }
  const arma::mat gradient = scatter - result.inverse;
  result.diagonal = 0.0;
  result.off_diagonal = 0.0;
  for (arma::uword i = 0; i < precision.n_rows; ++i) {
    result.diagonal = std::max(
      result.diagonal, std::abs(gradient(i, i))
    );
    for (arma::uword j = i + 1; j < precision.n_cols; ++j) {
      double residual = 0.0;
      if (precision(i, j) > 0.0) {
        residual = std::abs(gradient(i, j) + lambda);
      } else if (precision(i, j) < 0.0) {
        residual = std::abs(gradient(i, j) - lambda);
      } else {
        residual = std::max(0.0, std::abs(gradient(i, j)) - lambda);
      }
      result.off_diagonal = std::max(result.off_diagonal, residual);
    }
  }
  result.maximum = std::max(result.diagonal, result.off_diagonal);
  return result;
}

Rcpp::List ch3gp_glasso_result(const arma::mat& scatter,
                               const arma::mat& solution,
                               const double lambda,
                               const int iterations,
                               const bool converged,
                               const bool backtracking_failed,
                               const bool objective_descent,
                               const double relative_update,
                               const double accepted_step,
                               const int total_backtracking,
                               const std::vector<double>& history) {
  const Ch3gpGlassoKkt kkt = ch3gp_glasso_kkt(
    scatter, solution, lambda
  );
  return Rcpp::List::create(
    Rcpp::Named("solution") = solution,
    Rcpp::Named("inverse") = kkt.inverse,
    Rcpp::Named("objective") = ch3gp_glasso_objective(
      scatter, solution, lambda
    ),
    Rcpp::Named("objective_history") = Rcpp::wrap(history),
    Rcpp::Named("objective_descent") = objective_descent,
    Rcpp::Named("iterations") = iterations,
    Rcpp::Named("converged") = converged,
    Rcpp::Named("backtracking_failed") = backtracking_failed,
    Rcpp::Named("relative_update") = relative_update,
    Rcpp::Named("kkt_diagonal") = kkt.diagonal,
    Rcpp::Named("kkt_off_diagonal") = kkt.off_diagonal,
    Rcpp::Named("kkt_maximum") = kkt.maximum,
    Rcpp::Named("minimum_eigenvalue") =
      ch3gp_minimum_eigenvalue(solution),
    Rcpp::Named("accepted_step") = accepted_step,
    Rcpp::Named("total_backtracking") = total_backtracking
  );
}

struct Ch3gpClimeColumn {
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

void ch3gp_clime_certificate(const arma::mat& coefficient,
                             const arma::vec& target,
                             const double lambda,
                             const arma::vec& primal,
                             const arma::vec& dual,
                             Ch3gpClimeColumn& result) {
  const arma::vec feasibility = coefficient * primal - target;
  result.primal_violation = std::max(
    0.0, ch3gp_max_abs_vec(feasibility) - lambda
  );
  const arma::vec stationarity = coefficient.t() * dual;
  result.dual_violation = std::max(
    0.0, ch3gp_max_abs_vec(stationarity) - 1.0
  );
  result.stationarity_residual = 0.0;
  for (arma::uword k = 0; k < primal.n_elem; ++k) {
    double component = 0.0;
    if (primal(k) > 0.0) {
      component = std::abs(stationarity(k) + 1.0);
    } else if (primal(k) < 0.0) {
      component = std::abs(stationarity(k) - 1.0);
    } else {
      component = std::max(0.0, std::abs(stationarity(k)) - 1.0);
    }
    result.stationarity_residual = std::max(
      result.stationarity_residual, component
    );
  }
  result.primal_objective = arma::accu(arma::abs(primal));
  result.dual_objective = -arma::dot(target, dual) -
    lambda * arma::accu(arma::abs(dual));
  result.duality_gap = result.primal_objective - result.dual_objective;
  const double scale = std::max({
    1.0, std::abs(result.primal_objective),
    std::abs(result.dual_objective)
  });
  result.relative_gap = std::abs(result.duality_gap) / scale;
}

Ch3gpClimeColumn ch3gp_clime_column(const arma::mat& coefficient,
                                    const arma::vec& target,
                                    const double lambda,
                                    const double tolerance,
                                    const int max_iterations,
                                    const double operator_norm) {
  const arma::uword p = coefficient.n_cols;
  const double step = 0.99 / operator_norm;
  arma::vec primal(p, arma::fill::zeros);
  arma::vec extrapolated(p, arma::fill::zeros);
  arma::vec dual(p, arma::fill::zeros);

  Ch3gpClimeColumn result;
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

  if (lambda >= 1.0) {
    result.relative_update = 0.0;
    ch3gp_clime_certificate(
      coefficient, target, lambda, primal, dual, result
    );
    result.converged = result.primal_violation <= tolerance &&
      result.dual_violation <= tolerance &&
      result.stationarity_residual <= tolerance &&
      result.relative_gap <= tolerance;
    return result;
  }

  if (lambda == 0.0) {
    arma::mat inverse;
    if (arma::inv_sympd(inverse, coefficient)) {
      primal = inverse * target;
      arma::vec subgradient(p, arma::fill::zeros);
      for (arma::uword k = 0; k < p; ++k) {
        if (primal(k) > 0.0) subgradient(k) = 1.0;
        if (primal(k) < 0.0) subgradient(k) = -1.0;
      }
      dual = -inverse * subgradient;
      result.solution = primal;
      result.dual = dual;
      result.relative_update = 0.0;
      ch3gp_clime_certificate(
        coefficient, target, lambda, primal, dual, result
      );
      result.converged = result.primal_violation <= tolerance &&
        result.dual_violation <= tolerance &&
        result.stationarity_residual <= tolerance &&
        result.relative_gap <= tolerance;
      return result;
    }
  }

  for (int iteration = 1; iteration <= max_iterations; ++iteration) {
    const arma::vec dual_input = dual +
      step * (coefficient * extrapolated);
    arma::vec projected = dual_input / step;
    for (arma::uword k = 0; k < p; ++k) {
      projected(k) = std::max(
        target(k) - lambda,
        std::min(target(k) + lambda, projected(k))
      );
    }
    const arma::vec dual_new = dual_input - step * projected;
    const arma::vec primal_new = ch3gp_soft_vector(
      primal - step * (coefficient.t() * dual_new), step
    );
    const arma::vec extrapolated_new = 2.0 * primal_new - primal;
    result.relative_update = ch3gp_max_abs_vec(primal_new - primal) /
      std::max(1.0, ch3gp_max_abs_vec(primal));
    primal = primal_new;
    dual = dual_new;
    extrapolated = extrapolated_new;
    result.iterations = iteration;
    ch3gp_clime_certificate(
      coefficient, target, lambda, primal, dual, result
    );
    if (result.primal_violation <= tolerance &&
        result.dual_violation <= tolerance &&
        result.stationarity_residual <= tolerance &&
        result.relative_gap <= tolerance) {
      result.converged = true;
      break;
    }
    if ((iteration & 1023) == 0) Rcpp::checkUserInterrupt();
  }
  result.solution = primal;
  result.dual = dual;
  return result;
}

}  // namespace


// [[Rcpp::export]]
Rcpp::List cpp_ch3gp_ec2_l1(const arma::mat& correlation,
                            const double lambda,
                            const double tau,
                            const double rho,
                            const double tolerance,
                            const int max_iterations) {
  ch3gp_validate_finite_matrix(correlation, "correlation");
  if (correlation.n_rows == 0 ||
      correlation.n_rows != correlation.n_cols) {
    Rcpp::stop("`correlation` must be a non-empty square matrix.");
  }
  if (ch3gp_max_abs_mat(correlation - correlation.t()) >
      ch3gp_symmetry_tolerance(correlation)) {
    Rcpp::stop("`correlation` must be symmetric.");
  }
  for (arma::uword j = 0; j < correlation.n_rows; ++j) {
    if (std::abs(correlation(j, j) - 1.0) > 1e-10) {
      Rcpp::stop("`correlation` must have unit diagonal.");
    }
  }
  if (!std::isfinite(lambda) || lambda < 0.0 ||
      !std::isfinite(tau) || tau <= 0.0 || tau > 1.0 ||
      !std::isfinite(rho) || rho <= 0.0 ||
      !std::isfinite(tolerance) || tolerance <= 0.0 ||
      max_iterations < 1) {
    Rcpp::stop("Invalid EC2 l1 solver controls.");
  }

  arma::mat sparse = ch3gp_soft_off_diagonal(
    correlation, lambda, true
  );
  const double sto_minimum = ch3gp_minimum_eigenvalue(sparse);
  if (!std::isfinite(sto_minimum)) {
    Rcpp::stop("The EC2 soft-thresholded eigenvalues could not be evaluated.");
  }
  if (sto_minimum >= tau) {
    const arma::mat dual = correlation - sparse;
    const std::vector<double> history = {
      ch3gp_ec2_objective(correlation, sparse, lambda)
    };
    return ch3gp_ec2_result(
      correlation, sparse, sparse, dual, lambda, tau, rho, 0, true,
      0.0, 0.0, true, history
    );
  }

  arma::mat spectral = ch3gp_project_eigen_floor(sparse, tau);
  arma::mat dual(sparse.n_rows, sparse.n_cols, arma::fill::zeros);
  std::vector<double> history = {
    ch3gp_ec2_objective(correlation, sparse, lambda)
  };
  bool converged = false;
  double relative_update = std::numeric_limits<double>::infinity();
  double dual_residual = std::numeric_limits<double>::infinity();
  int iterations = 0;

  for (int iteration = 1; iteration <= max_iterations; ++iteration) {
    const arma::mat sparse_previous = sparse;
    const arma::mat spectral_previous = spectral;
    sparse = ch3gp_soft_off_diagonal(
      spectral + dual / rho, lambda / rho, true
    );
    spectral = ch3gp_project_eigen_floor(
      (correlation + rho * sparse - dual) / (1.0 + rho), tau
    );
    dual += rho * (spectral - sparse);
    dual = 0.5 * (dual + dual.t());
    relative_update = ch3gp_max_abs_mat(sparse - sparse_previous) /
      std::max(1.0, ch3gp_max_abs_mat(sparse_previous));
    dual_residual = rho * ch3gp_max_abs_mat(
      spectral - spectral_previous
    );
    iterations = iteration;
    history.push_back(ch3gp_ec2_objective(
      correlation, sparse, lambda
    ));
    const Ch3gpEc2Certificate certificate = ch3gp_ec2_certificate(
      correlation, sparse, spectral, dual, lambda, tau, rho
    );
    if (certificate.maximum <= tolerance &&
        relative_update <= tolerance &&
        dual_residual <= tolerance &&
        certificate.minimum_eigenvalue > 0.0) {
      converged = true;
      break;
    }
    if ((iteration & 255) == 0) Rcpp::checkUserInterrupt();
  }

  return ch3gp_ec2_result(
    correlation, sparse, spectral, dual, lambda, tau, rho,
    iterations, converged, relative_update, dual_residual, false, history
  );
}


// [[Rcpp::export]]
Rcpp::List cpp_ch3gp_offdiag_glasso(
    const arma::mat& scatter,
    const double lambda,
    const double tolerance,
    const int max_iterations,
    const double initial_step,
    const int max_backtracking) {
  ch3gp_validate_finite_matrix(scatter, "scatter");
  if (scatter.n_rows == 0 || scatter.n_rows != scatter.n_cols) {
    Rcpp::stop("`scatter` must be a non-empty square matrix.");
  }
  if (ch3gp_max_abs_mat(scatter - scatter.t()) >
      ch3gp_symmetry_tolerance(scatter)) {
    Rcpp::stop("`scatter` must be symmetric.");
  }
  if (!std::isfinite(lambda) || lambda < 0.0 ||
      !std::isfinite(tolerance) || tolerance <= 0.0 ||
      max_iterations < 1 || !std::isfinite(initial_step) ||
      initial_step <= 0.0 || max_backtracking < 1) {
    Rcpp::stop("Invalid off-diagonal graphical-lasso solver controls.");
  }

  const arma::uword p = scatter.n_rows;
  arma::mat current(p, p, arma::fill::zeros);
  for (arma::uword i = 0; i < p; ++i) {
    if (scatter(i, i) <= 0.0 ||
        !std::isfinite(1.0 / scatter(i, i))) {
      Rcpp::stop(
        "The unpenalized diagonal objective is not coercive because a "
        "scatter diagonal is non-positive or non-invertible; no ridge "
        "was added."
      );
    }
    current(i, i) = 1.0 / scatter(i, i);
  }

  if (lambda == 0.0) {
    arma::mat exact;
    if (!arma::inv_sympd(exact, scatter)) {
      Rcpp::stop(
        "With zero penalty `scatter` must be strictly positive definite; "
        "no pseudoinverse or ridge was used."
      );
    }
    const Ch3gpGlassoKkt kkt = ch3gp_glasso_kkt(
      scatter, exact, lambda
    );
    const std::vector<double> history = {
      ch3gp_glasso_objective(scatter, exact, lambda)
    };
    return ch3gp_glasso_result(
      scatter, exact, lambda, 0,
      std::isfinite(kkt.maximum) && kkt.maximum <= tolerance,
      false, true, 0.0, NA_REAL, 0, history
    );
  }

  double current_objective = ch3gp_glasso_objective(
    scatter, current, lambda
  );
  std::vector<double> history = {current_objective};
  const Ch3gpGlassoKkt initial_kkt = ch3gp_glasso_kkt(
    scatter, current, lambda
  );
  if (std::isfinite(initial_kkt.maximum) &&
      initial_kkt.maximum <= tolerance) {
    return ch3gp_glasso_result(
      scatter, current, lambda, 0, true, false, true, 0.0,
      NA_REAL, 0, history
    );
  }

  double step = initial_step;
  double relative_update = std::numeric_limits<double>::infinity();
  bool converged = false;
  bool backtracking_failed = false;
  bool objective_descent = true;
  int iterations = 0;
  int total_backtracking = 0;

  for (int iteration = 1; iteration <= max_iterations; ++iteration) {
    arma::mat inverse;
    if (!arma::inv_sympd(inverse, current)) {
      Rcpp::stop(
        "A graphical-lasso iterate lost positive definiteness; no ridge "
        "or eigenvalue floor was applied."
      );
    }
    const arma::mat gradient = scatter - inverse;
    const double smooth_current = arma::trace(scatter * current) -
      ch3gp_logdet_spd(current);
    bool accepted = false;
    arma::mat candidate;
    double candidate_objective = std::numeric_limits<double>::infinity();
    double candidate_step = std::min(initial_step, step * 1.2);

    for (int backtracking = 0;
         backtracking < max_backtracking; ++backtracking) {
      if (!std::isfinite(candidate_step) || candidate_step <= 0.0) break;
      candidate = ch3gp_soft_off_diagonal(
        current - candidate_step * gradient,
        candidate_step * lambda, false
      );
      arma::mat upper;
      if (candidate.is_finite() && arma::chol(upper, candidate)) {
        const arma::mat difference = candidate - current;
        const double smooth_candidate = arma::trace(scatter * candidate) -
          ch3gp_logdet_spd(candidate);
        const double majorizer = smooth_current +
          arma::accu(gradient % difference) +
          arma::accu(arma::square(difference)) /
          (2.0 * candidate_step);
        candidate_objective = smooth_candidate +
          lambda * ch3gp_off_diagonal_l1(candidate);
        const double slack = 128.0 *
          std::numeric_limits<double>::epsilon() *
          std::max({1.0, std::abs(smooth_current),
                    std::abs(current_objective)});
        if (std::isfinite(smooth_candidate) &&
            std::isfinite(candidate_objective) &&
            smooth_candidate <= majorizer + slack &&
            candidate_objective <= current_objective + slack) {
          accepted = true;
          total_backtracking += backtracking;
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
    const arma::mat difference = candidate - current;
    relative_update = arma::norm(difference, "fro") /
      std::max(1.0, arma::norm(current, "fro"));
    const double monotonicity_slack = 256.0 *
      std::numeric_limits<double>::epsilon() *
      std::max(1.0, std::abs(current_objective));
    if (candidate_objective > current_objective + monotonicity_slack) {
      objective_descent = false;
    }
    current = candidate;
    current_objective = candidate_objective;
    history.push_back(current_objective);
    step = candidate_step;
    iterations = iteration;

    const Ch3gpGlassoKkt kkt = ch3gp_glasso_kkt(
      scatter, current, lambda
    );
    if (std::isfinite(kkt.maximum) && kkt.maximum <= tolerance &&
        relative_update <= tolerance) {
      converged = true;
      break;
    }
    if ((iteration & 255) == 0) Rcpp::checkUserInterrupt();
  }

  return ch3gp_glasso_result(
    scatter, current, lambda, iterations, converged,
    backtracking_failed, objective_descent, relative_update, step,
    total_backtracking, history
  );
}


// [[Rcpp::export]]
Rcpp::List cpp_ch3gp_clime(const arma::mat& scatter,
                           const double lambda,
                           const double tolerance,
                           const int max_iterations) {
  ch3gp_validate_finite_matrix(scatter, "scatter");
  if (scatter.n_rows == 0 || scatter.n_rows != scatter.n_cols) {
    Rcpp::stop("`scatter` must be a non-empty square matrix.");
  }
  if (ch3gp_max_abs_mat(scatter - scatter.t()) >
      ch3gp_symmetry_tolerance(scatter)) {
    Rcpp::stop("`scatter` must be symmetric.");
  }
  if (!std::isfinite(lambda) || lambda < 0.0 ||
      !std::isfinite(tolerance) || tolerance <= 0.0 ||
      max_iterations < 1) {
    Rcpp::stop("Invalid CLIME solver controls.");
  }
  arma::vec eigenvalues;
  if (!ch3gp_eigenvalues(scatter, eigenvalues)) {
    Rcpp::stop("The CLIME operator norm could not be evaluated.");
  }
  const double operator_norm = arma::abs(eigenvalues).max();
  if (!std::isfinite(operator_norm) || operator_norm <= 0.0) {
    Rcpp::stop(
      "The CLIME coefficient matrix has zero or non-finite norm; no "
      "constraint repair was applied."
    );
  }

  const arma::uword p = scatter.n_cols;
  arma::mat raw(p, p, arma::fill::zeros);
  arma::mat dual(p, p, arma::fill::zeros);
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
    const Ch3gpClimeColumn fit = ch3gp_clime_column(
      scatter, target, lambda, tolerance, max_iterations, operator_norm
    );
    raw.col(j) = fit.solution;
    dual.col(j) = fit.dual;
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
      // The <= tie rule is deterministic: keep raw(i,j), then mirror it.
      const double value =
        std::abs(raw(i, j)) <= std::abs(raw(j, i)) ?
        raw(i, j) : raw(j, i);
      symmetrized(i, j) = value;
      symmetrized(j, i) = value;
    }
  }
  const double final_feasibility = std::max(
    0.0, ch3gp_max_abs_mat(
      scatter * symmetrized - arma::eye<arma::mat>(p, p)
    ) - lambda
  );

  return Rcpp::List::create(
    Rcpp::Named("raw_solution") = raw,
    Rcpp::Named("solution") = symmetrized,
    Rcpp::Named("dual_solution") = dual,
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
    Rcpp::Named("symmetrized_feasibility_violation") = final_feasibility,
    Rcpp::Named("symmetrized_minimum_eigenvalue") =
      ch3gp_minimum_eigenvalue(symmetrized)
  );
}
