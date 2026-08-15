// Chapter 5: sparse quadratic discriminant analysis.
//
// The Dantzig kernels below apply their matrix operators directly.  They do
// not construct a p^2 by p^2 Kronecker matrix.  Numerical failures are
// reported to R; no eigenvalue floor, pseudoinverse, implicit ridge, or
// post-hoc feasibility projection is used.

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp17)]]

#include <RcppArmadillo.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <string>

namespace {

double c5qda_max_abs_mat(const arma::mat& x) {
  return x.n_elem == 0 ? 0.0 : arma::abs(x).max();
}

double c5qda_max_abs_vec(const arma::vec& x) {
  return x.n_elem == 0 ? 0.0 : arma::abs(x).max();
}

arma::mat c5qda_soft_mat(const arma::mat& x, const double threshold) {
  return arma::sign(x) % arma::max(
    arma::abs(x) - threshold, arma::zeros<arma::mat>(x.n_rows, x.n_cols)
  );
}

arma::vec c5qda_soft_vec(const arma::vec& x, const double threshold) {
  return arma::sign(x) % arma::max(
    arma::abs(x) - threshold, arma::zeros<arma::vec>(x.n_elem)
  );
}

void c5qda_check_square_pair(const arma::mat& s1, const arma::mat& s2) {
  if (s1.n_rows == 0 || s1.n_rows != s1.n_cols ||
      s2.n_rows != s2.n_cols || s1.n_rows != s2.n_rows) {
    Rcpp::stop("`S1` and `S2` must be non-empty square matrices of the same size.");
  }
  if (!s1.is_finite() || !s2.is_finite()) {
    Rcpp::stop("`S1` and `S2` must contain only finite values.");
  }
}

arma::mat c5qda_operator(const arma::mat& s1, const arma::mat& s2,
                         const arma::mat& d) {
  return 0.5 * (s1 * d * s2 + s2 * d * s1);
}

arma::mat c5qda_adjoint(const arma::mat& s1, const arma::mat& s2,
                        const arma::mat& y) {
  return 0.5 * (s1.t() * y * s2.t() + s2.t() * y * s1.t());
}

double c5qda_l1_stationarity(const arma::mat& primal,
                             const arma::mat& adjoint_dual,
                             const double zero_tolerance) {
  double answer = 0.0;
  for (arma::uword i = 0; i < primal.n_rows; ++i) {
    for (arma::uword j = 0; j < primal.n_cols; ++j) {
      double component;
      if (std::abs(primal(i, j)) > zero_tolerance) {
        component = std::abs(
          adjoint_dual(i, j) + (primal(i, j) > 0.0 ? 1.0 : -1.0)
        );
      } else {
        component = std::max(0.0, std::abs(adjoint_dual(i, j)) - 1.0);
      }
      answer = std::max(answer, component);
    }
  }
  return answer;
}

double c5qda_l1_stationarity(const arma::vec& primal,
                             const arma::vec& adjoint_dual,
                             const double zero_tolerance) {
  double answer = 0.0;
  for (arma::uword j = 0; j < primal.n_elem; ++j) {
    double component;
    if (std::abs(primal(j)) > zero_tolerance) {
      component = std::abs(
        adjoint_dual(j) + (primal(j) > 0.0 ? 1.0 : -1.0)
      );
    } else {
      component = std::max(0.0, std::abs(adjoint_dual(j)) - 1.0);
    }
    answer = std::max(answer, component);
  }
  return answer;
}

double c5qda_lasso_kkt(const arma::vec& beta, const arma::vec& gradient,
                       const double lambda, const double zero_tolerance) {
  double answer = 0.0;
  for (arma::uword j = 0; j < beta.n_elem; ++j) {
    double component;
    if (std::abs(beta(j)) > zero_tolerance) {
      component = std::abs(
        gradient(j) + lambda * (beta(j) > 0.0 ? 1.0 : -1.0)
      );
    } else {
      component = std::max(0.0, std::abs(gradient(j)) - lambda);
    }
    answer = std::max(answer, component);
  }
  return answer;
}

bool c5qda_logdet_positive(const arma::mat& x, double& log_determinant,
                           double& sign) {
  sign = 0.0;
  log_determinant = std::numeric_limits<double>::quiet_NaN();
  return arma::log_det(log_determinant, sign, x) && sign > 0.0 &&
    std::isfinite(log_determinant);
}

arma::mat c5qda_mle_covariance(const arma::mat& x,
                               const arma::rowvec& center) {
  arma::mat centered = x.each_row() - center;
  return centered.t() * centered / static_cast<double>(x.n_rows);
}

} // namespace


// [[Rcpp::export]]
arma::mat cpp_c5qda_operator(const arma::mat& S1, const arma::mat& S2,
                             const arma::mat& D, const bool adjoint = false) {
  c5qda_check_square_pair(S1, S2);
  if (D.n_rows != S1.n_rows || D.n_cols != S1.n_cols || !D.is_finite()) {
    Rcpp::stop("`D` must be a finite matrix conformable with `S1` and `S2`.");
  }
  return adjoint ? c5qda_adjoint(S1, S2, D) :
    c5qda_operator(S1, S2, D);
}


// [[Rcpp::export]]
Rcpp::List cpp_c5qda_li_shao(const arma::mat& x1, const arma::mat& x2,
                             const double threshold_mean,
                             const double threshold_difference,
                             const double threshold_covariance,
                             const double ridge) {
  if (x1.n_rows < 2 || x2.n_rows < 2 || x1.n_cols == 0 ||
      x1.n_cols != x2.n_cols || !x1.is_finite() || !x2.is_finite()) {
    Rcpp::stop("The two class matrices must be finite, have equal positive column counts, and at least two rows each.");
  }
  if (!std::isfinite(threshold_mean) || threshold_mean < 0.0 ||
      !std::isfinite(threshold_difference) || threshold_difference < 0.0 ||
      !std::isfinite(threshold_covariance) || threshold_covariance < 0.0 ||
      !std::isfinite(ridge) || ridge < 0.0) {
    Rcpp::stop("Thresholds and `ridge` must be finite and non-negative.");
  }

  const arma::uword p = x1.n_cols;
  const arma::rowvec mean1 = arma::mean(x1, 0);
  const arma::rowvec mean2 = arma::mean(x2, 0);
  const arma::vec difference = (mean2 - mean1).t();
  arma::vec thresholded_difference = difference;
  for (arma::uword j = 0; j < p; ++j) {
    if (!(std::abs(difference(j)) > threshold_mean)) {
      thresholded_difference(j) = 0.0;
    }
  }

  const arma::mat sample1 = c5qda_mle_covariance(x1, mean1);
  const arma::mat sample2 = c5qda_mle_covariance(x2, mean2);
  const double n1 = static_cast<double>(x1.n_rows);
  const double n2 = static_cast<double>(x2.n_rows);
  const arma::mat pooled = (n1 * sample1 + n2 * sample2) / (n1 + n2);
  arma::mat covariance1(p, p, arma::fill::zeros);
  arma::mat covariance2(p, p, arma::fill::zeros);

  for (arma::uword i = 0; i < p; ++i) {
    for (arma::uword j = 0; j < p; ++j) {
      if (std::abs(sample1(i, j) - sample2(i, j)) <=
          threshold_difference) {
        covariance1(i, j) = pooled(i, j);
        covariance2(i, j) = pooled(i, j);
      } else {
        covariance1(i, j) = sample1(i, j);
        covariance2(i, j) = sample2(i, j);
      }
      if (i != j) {
        if (!(std::abs(covariance1(i, j)) > threshold_covariance)) {
          covariance1(i, j) = 0.0;
        }
        if (!(std::abs(covariance2(i, j)) > threshold_covariance)) {
          covariance2(i, j) = 0.0;
        }
      }
    }
  }
  covariance1 = 0.5 * (covariance1 + covariance1.t());
  covariance2 = 0.5 * (covariance2 + covariance2.t());
  covariance1.diag() += ridge;
  covariance2.diag() += ridge;

  arma::mat precision1(p, p, arma::fill::value(
    std::numeric_limits<double>::quiet_NaN()
  ));
  arma::mat precision2 = precision1;
  double logdet1 = std::numeric_limits<double>::quiet_NaN();
  double logdet2 = logdet1;
  double sign1 = 0.0;
  double sign2 = 0.0;
  const bool determinant1 = c5qda_logdet_positive(covariance1, logdet1, sign1);
  const bool determinant2 = c5qda_logdet_positive(covariance2, logdet2, sign2);
  const bool inverse1 = determinant1 && arma::inv_sympd(precision1, covariance1);
  const bool inverse2 = determinant2 && arma::inv_sympd(precision2, covariance2);
  const bool valid = inverse1 && inverse2 && precision1.is_finite() &&
    precision2.is_finite();

  std::string failure_stage;
  if (!valid) {
    failure_stage = "thresholded_covariance_inversion";
  }
  return Rcpp::List::create(
    Rcpp::Named("valid") = valid,
    Rcpp::Named("failure_stage") = failure_stage,
    Rcpp::Named("mean1") = mean1,
    Rcpp::Named("mean2") = mean2,
    Rcpp::Named("difference") = difference,
    Rcpp::Named("thresholded_difference") = thresholded_difference,
    Rcpp::Named("sample_covariance1") = sample1,
    Rcpp::Named("sample_covariance2") = sample2,
    Rcpp::Named("pooled_covariance") = pooled,
    Rcpp::Named("covariance1") = covariance1,
    Rcpp::Named("covariance2") = covariance2,
    Rcpp::Named("precision1") = precision1,
    Rcpp::Named("precision2") = precision2,
    Rcpp::Named("logdet1") = logdet1,
    Rcpp::Named("logdet2") = logdet2,
    Rcpp::Named("ridge") = ridge
  );
}


// [[Rcpp::export]]
Rcpp::List cpp_c5qda_jiang_matrix(const arma::mat& S1,
                                  const arma::mat& S2,
                                  const double lambda,
                                  const double rho,
                                  const double tolerance,
                                  const int max_iterations) {
  c5qda_check_square_pair(S1, S2);
  if (!std::isfinite(lambda) || lambda < 0.0 || !std::isfinite(rho) ||
      rho <= 0.0 || !std::isfinite(tolerance) || tolerance <= 0.0 ||
      max_iterations < 1) {
    Rcpp::stop("Invalid Jiang matrix-solver control.");
  }

  arma::vec eigen1;
  arma::vec eigen2;
  arma::mat vectors1;
  arma::mat vectors2;
  const arma::mat sym1 = 0.5 * (S1 + S1.t());
  const arma::mat sym2 = 0.5 * (S2 + S2.t());
  if (!arma::eig_sym(eigen1, vectors1, sym1) ||
      !arma::eig_sym(eigen2, vectors2, sym2) ||
      eigen1.min() < -tolerance || eigen2.min() < -tolerance) {
    return Rcpp::List::create(
      Rcpp::Named("solution") = arma::mat(S1.n_rows, S1.n_cols,
                                           arma::fill::zeros),
      Rcpp::Named("converged") = false,
      Rcpp::Named("iterations") = 0,
      Rcpp::Named("failure_stage") = "sample_covariance_psd"
    );
  }

  const arma::mat target = sym1 - sym2;
  arma::mat d(S1.n_rows, S1.n_cols, arma::fill::zeros);
  arma::mat z = d;
  arma::mat scaled_dual = d;
  double primal_residual = std::numeric_limits<double>::infinity();
  double dual_residual = primal_residual;
  double relative_update = primal_residual;
  double kkt = primal_residual;
  double objective = std::numeric_limits<double>::infinity();
  int iterations = 0;
  bool converged = false;

  for (int iteration = 1; iteration <= max_iterations; ++iteration) {
    const arma::mat rhs = target + rho * (z - scaled_dual);
    arma::mat transformed = vectors1.t() * rhs * vectors2;
    for (arma::uword i = 0; i < transformed.n_rows; ++i) {
      for (arma::uword j = 0; j < transformed.n_cols; ++j) {
        const double denominator = eigen1(i) * eigen2(j) + rho;
        if (!(denominator > 0.0) || !std::isfinite(denominator)) {
          return Rcpp::List::create(
            Rcpp::Named("solution") = z,
            Rcpp::Named("converged") = false,
            Rcpp::Named("iterations") = iteration - 1,
            Rcpp::Named("failure_stage") = "jiang_matrix_linear_system"
          );
        }
        transformed(i, j) /= denominator;
      }
    }
    const arma::mat d_new = vectors1 * transformed * vectors2.t();
    const arma::mat z_old = z;
    z = c5qda_soft_mat(d_new + scaled_dual, lambda / rho);
    scaled_dual += d_new - z;

    primal_residual = c5qda_max_abs_mat(d_new - z);
    dual_residual = rho * c5qda_max_abs_mat(z - z_old);
    relative_update = c5qda_max_abs_mat(z - z_old) /
      std::max(1.0, c5qda_max_abs_mat(z_old));
    d = d_new;

    const arma::mat gradient = sym1 * z * sym2 - target;
    kkt = 0.0;
    for (arma::uword i = 0; i < z.n_rows; ++i) {
      for (arma::uword j = 0; j < z.n_cols; ++j) {
        double component;
        if (std::abs(z(i, j)) > tolerance) {
          component = std::abs(
            gradient(i, j) + lambda * (z(i, j) > 0.0 ? 1.0 : -1.0)
          );
        } else {
          component = std::max(0.0, std::abs(gradient(i, j)) - lambda);
        }
        kkt = std::max(kkt, component);
      }
    }
    objective = 0.5 * arma::accu(z % (sym1 * z * sym2)) -
      arma::accu(z % target) + lambda * arma::accu(arma::abs(z));
    iterations = iteration;
    if (primal_residual <= tolerance && dual_residual <= tolerance &&
        kkt <= tolerance) {
      converged = true;
      break;
    }
  }

  return Rcpp::List::create(
    Rcpp::Named("solution") = z,
    Rcpp::Named("converged") = converged,
    Rcpp::Named("iterations") = iterations,
    Rcpp::Named("primal_residual") = primal_residual,
    Rcpp::Named("dual_residual") = dual_residual,
    Rcpp::Named("relative_update") = relative_update,
    Rcpp::Named("kkt_residual") = kkt,
    Rcpp::Named("objective") = objective,
    Rcpp::Named("failure_stage") = converged ? "" : "jiang_matrix_solver"
  );
}


// [[Rcpp::export]]
Rcpp::List cpp_c5qda_jiang_vector(const arma::mat& H,
                                  const arma::vec& gamma,
                                  const double lambda,
                                  const double tolerance,
                                  const int max_iterations) {
  if (H.n_rows == 0 || H.n_rows != H.n_cols || gamma.n_elem != H.n_rows ||
      !H.is_finite() || !gamma.is_finite()) {
    Rcpp::stop("`H` and `gamma` must be finite and conformable.");
  }
  if (!std::isfinite(lambda) || lambda < 0.0 ||
      !std::isfinite(tolerance) || tolerance <= 0.0 || max_iterations < 1) {
    Rcpp::stop("Invalid Jiang vector-solver control.");
  }
  const arma::mat sym_h = 0.5 * (H + H.t());
  arma::vec eigenvalues;
  if (!arma::eig_sym(eigenvalues, sym_h) || eigenvalues.min() < -tolerance) {
    return Rcpp::List::create(
      Rcpp::Named("solution") = arma::vec(H.n_rows, arma::fill::zeros),
      Rcpp::Named("converged") = false,
      Rcpp::Named("iterations") = 0,
      Rcpp::Named("failure_stage") = "jiang_vector_hessian_psd"
    );
  }
  const double lipschitz = eigenvalues.max();
  if (!(lipschitz > 0.0) || !std::isfinite(lipschitz)) {
    const bool zero_optimum = c5qda_max_abs_vec(gamma) <= lambda;
    return Rcpp::List::create(
      Rcpp::Named("solution") = arma::vec(H.n_rows, arma::fill::zeros),
      Rcpp::Named("converged") = zero_optimum,
      Rcpp::Named("iterations") = 0,
      Rcpp::Named("relative_update") = 0.0,
      Rcpp::Named("kkt_residual") = std::max(0.0, c5qda_max_abs_vec(gamma) - lambda),
      Rcpp::Named("objective") = 0.0,
      Rcpp::Named("failure_stage") = zero_optimum ? "" :
        "jiang_vector_unbounded"
    );
  }

  arma::vec beta(H.n_rows, arma::fill::zeros);
  arma::vec extrapolated = beta;
  double momentum = 1.0;
  double relative_update = std::numeric_limits<double>::infinity();
  double kkt = relative_update;
  double objective = std::numeric_limits<double>::infinity();
  int iterations = 0;
  bool converged = false;
  for (int iteration = 1; iteration <= max_iterations; ++iteration) {
    const arma::vec gradient_at_extrapolated = sym_h * extrapolated - gamma;
    const arma::vec beta_new = c5qda_soft_vec(
      extrapolated - gradient_at_extrapolated / lipschitz,
      lambda / lipschitz
    );
    const double momentum_new = 0.5 *
      (1.0 + std::sqrt(1.0 + 4.0 * momentum * momentum));
    extrapolated = beta_new + ((momentum - 1.0) / momentum_new) *
      (beta_new - beta);
    relative_update = c5qda_max_abs_vec(beta_new - beta) /
      std::max(1.0, c5qda_max_abs_vec(beta));
    beta = beta_new;
    momentum = momentum_new;
    const arma::vec gradient = sym_h * beta - gamma;
    kkt = c5qda_lasso_kkt(beta, gradient, lambda, tolerance);
    objective = 0.5 * arma::dot(beta, sym_h * beta) -
      arma::dot(gamma, beta) + lambda * arma::accu(arma::abs(beta));
    iterations = iteration;
    if (relative_update <= tolerance && kkt <= tolerance) {
      converged = true;
      break;
    }
  }
  return Rcpp::List::create(
    Rcpp::Named("solution") = beta,
    Rcpp::Named("converged") = converged,
    Rcpp::Named("iterations") = iterations,
    Rcpp::Named("relative_update") = relative_update,
    Rcpp::Named("kkt_residual") = kkt,
    Rcpp::Named("objective") = objective,
    Rcpp::Named("failure_stage") = converged ? "" : "jiang_vector_solver"
  );
}


// [[Rcpp::export]]
Rcpp::List cpp_c5qda_dantzig_matrix(const arma::mat& S1,
                                    const arma::mat& S2,
                                    const double lambda,
                                    const double tolerance,
                                    const int max_iterations) {
  c5qda_check_square_pair(S1, S2);
  if (!std::isfinite(lambda) || lambda < 0.0 ||
      !std::isfinite(tolerance) || tolerance <= 0.0 || max_iterations < 1) {
    Rcpp::stop("Invalid matrix-Dantzig solver control.");
  }
  const arma::mat s1 = 0.5 * (S1 + S1.t());
  const arma::mat s2 = 0.5 * (S2 + S2.t());
  const arma::mat target = s1 - s2;
  const double operator_bound = arma::norm(s1, 2) * arma::norm(s2, 2);
  const arma::uword p = s1.n_rows;
  arma::mat primal(p, p, arma::fill::zeros);
  arma::mat extrapolated = primal;
  arma::mat dual = primal;

  if (!(operator_bound > 0.0) || !std::isfinite(operator_bound)) {
    const double feasible_maxnorm = c5qda_max_abs_mat(target);
    const bool feasible = feasible_maxnorm <= lambda + tolerance;
    return Rcpp::List::create(
      Rcpp::Named("solution") = primal,
      Rcpp::Named("dual") = dual,
      Rcpp::Named("converged") = feasible,
      Rcpp::Named("iterations") = 0,
      Rcpp::Named("operator_norm_bound") = operator_bound,
      Rcpp::Named("feasible_maxnorm") = feasible_maxnorm,
      Rcpp::Named("primal_violation") = std::max(0.0, feasible_maxnorm - lambda),
      Rcpp::Named("dual_violation") = 0.0,
      Rcpp::Named("stationarity_residual") = 0.0,
      Rcpp::Named("primal_objective") = 0.0,
      Rcpp::Named("dual_objective") = 0.0,
      Rcpp::Named("duality_gap") = 0.0,
      Rcpp::Named("relative_gap") = 0.0,
      Rcpp::Named("failure_stage") = feasible ? "" : "matrix_dantzig_infeasible"
    );
  }

  const double step = 0.99 / operator_bound;
  double primal_violation = std::numeric_limits<double>::infinity();
  double feasible_maxnorm = primal_violation;
  double dual_violation = primal_violation;
  double stationarity = primal_violation;
  double primal_objective = 0.0;
  double dual_objective = -primal_violation;
  double duality_gap = primal_violation;
  double relative_gap = primal_violation;
  double relative_update = primal_violation;
  int iterations = 0;
  bool converged = false;

  for (int iteration = 1; iteration <= max_iterations; ++iteration) {
    const arma::mat dual_input = dual + step *
      c5qda_operator(s1, s2, extrapolated);
    arma::mat projected = dual_input / step;
    for (arma::uword i = 0; i < p; ++i) {
      for (arma::uword j = 0; j < p; ++j) {
        projected(i, j) = std::max(
          target(i, j) - lambda,
          std::min(target(i, j) + lambda, projected(i, j))
        );
      }
    }
    const arma::mat dual_new = dual_input - step * projected;
    const arma::mat primal_new = c5qda_soft_mat(
      primal - step * c5qda_adjoint(s1, s2, dual_new), step
    );
    const arma::mat extrapolated_new = 2.0 * primal_new - primal;
    relative_update = c5qda_max_abs_mat(primal_new - primal) /
      std::max(1.0, c5qda_max_abs_mat(primal));
    primal = primal_new;
    dual = dual_new;
    extrapolated = extrapolated_new;

    const arma::mat residual = c5qda_operator(s1, s2, primal) - target;
    feasible_maxnorm = c5qda_max_abs_mat(residual);
    primal_violation = std::max(0.0, feasible_maxnorm - lambda);
    const arma::mat adjoint_dual = c5qda_adjoint(s1, s2, dual);
    dual_violation = std::max(0.0, c5qda_max_abs_mat(adjoint_dual) - 1.0);
    stationarity = c5qda_l1_stationarity(primal, adjoint_dual, tolerance);
    primal_objective = arma::accu(arma::abs(primal));
    dual_objective = -arma::accu(target % dual) -
      lambda * arma::accu(arma::abs(dual));
    duality_gap = primal_objective - dual_objective;
    relative_gap = std::abs(duality_gap) / std::max(
      1.0, std::max(std::abs(primal_objective), std::abs(dual_objective))
    );
    iterations = iteration;
    if (primal_violation <= tolerance && dual_violation <= tolerance &&
        stationarity <= tolerance && relative_gap <= tolerance) {
      converged = true;
      break;
    }
  }

  const arma::mat symmetrized = 0.5 * (primal + primal.t());
  const arma::mat sym_residual = c5qda_operator(s1, s2, symmetrized) - target;
  const double sym_feasible_maxnorm = c5qda_max_abs_mat(sym_residual);
  const double sym_primal_violation = std::max(
    0.0, sym_feasible_maxnorm - lambda
  );
  if (sym_primal_violation > tolerance) {
    converged = false;
  }
  return Rcpp::List::create(
    Rcpp::Named("solution") = symmetrized,
    Rcpp::Named("unsymmetrized_solution") = primal,
    Rcpp::Named("dual") = dual,
    Rcpp::Named("converged") = converged,
    Rcpp::Named("iterations") = iterations,
    Rcpp::Named("operator_norm_bound") = operator_bound,
    Rcpp::Named("relative_update") = relative_update,
    Rcpp::Named("feasible_maxnorm") = sym_feasible_maxnorm,
    Rcpp::Named("primal_violation") = sym_primal_violation,
    Rcpp::Named("unsymmetrized_feasible_maxnorm") = feasible_maxnorm,
    Rcpp::Named("dual_violation") = dual_violation,
    Rcpp::Named("stationarity_residual") = stationarity,
    Rcpp::Named("primal_objective") = arma::accu(arma::abs(symmetrized)),
    Rcpp::Named("dual_objective") = dual_objective,
    Rcpp::Named("duality_gap") = duality_gap,
    Rcpp::Named("relative_gap") = relative_gap,
    Rcpp::Named("matrix_free") = true,
    Rcpp::Named("failure_stage") = converged ? "" : "matrix_dantzig_solver"
  );
}


// [[Rcpp::export]]
Rcpp::List cpp_c5qda_dantzig_vector(const arma::mat& S,
                                    const arma::vec& target,
                                    const double lambda,
                                    const double tolerance,
                                    const int max_iterations) {
  if (S.n_rows == 0 || S.n_rows != S.n_cols || target.n_elem != S.n_rows ||
      !S.is_finite() || !target.is_finite()) {
    Rcpp::stop("`S` and `target` must be finite and conformable.");
  }
  if (!std::isfinite(lambda) || lambda < 0.0 ||
      !std::isfinite(tolerance) || tolerance <= 0.0 || max_iterations < 1) {
    Rcpp::stop("Invalid vector-Dantzig solver control.");
  }
  const arma::mat a = 0.5 * (S + S.t());
  const double operator_bound = arma::norm(a, 2);
  arma::vec primal(S.n_rows, arma::fill::zeros);
  arma::vec extrapolated = primal;
  arma::vec dual = primal;
  if (!(operator_bound > 0.0) || !std::isfinite(operator_bound)) {
    const double feasible_maxnorm = c5qda_max_abs_vec(target);
    const bool feasible = feasible_maxnorm <= lambda + tolerance;
    return Rcpp::List::create(
      Rcpp::Named("solution") = primal,
      Rcpp::Named("dual") = dual,
      Rcpp::Named("converged") = feasible,
      Rcpp::Named("iterations") = 0,
      Rcpp::Named("operator_norm_bound") = operator_bound,
      Rcpp::Named("feasible_maxnorm") = feasible_maxnorm,
      Rcpp::Named("primal_violation") = std::max(0.0, feasible_maxnorm - lambda),
      Rcpp::Named("dual_violation") = 0.0,
      Rcpp::Named("stationarity_residual") = 0.0,
      Rcpp::Named("primal_objective") = 0.0,
      Rcpp::Named("dual_objective") = 0.0,
      Rcpp::Named("duality_gap") = 0.0,
      Rcpp::Named("relative_gap") = 0.0,
      Rcpp::Named("failure_stage") = feasible ? "" : "vector_dantzig_infeasible"
    );
  }
  const double step = 0.99 / operator_bound;
  double primal_violation = std::numeric_limits<double>::infinity();
  double feasible_maxnorm = primal_violation;
  double dual_violation = primal_violation;
  double stationarity = primal_violation;
  double primal_objective = 0.0;
  double dual_objective = -primal_violation;
  double duality_gap = primal_violation;
  double relative_gap = primal_violation;
  double relative_update = primal_violation;
  int iterations = 0;
  bool converged = false;
  for (int iteration = 1; iteration <= max_iterations; ++iteration) {
    const arma::vec dual_input = dual + step * (a * extrapolated);
    arma::vec projected = dual_input / step;
    for (arma::uword j = 0; j < projected.n_elem; ++j) {
      projected(j) = std::max(
        target(j) - lambda,
        std::min(target(j) + lambda, projected(j))
      );
    }
    const arma::vec dual_new = dual_input - step * projected;
    const arma::vec primal_new = c5qda_soft_vec(
      primal - step * (a.t() * dual_new), step
    );
    const arma::vec extrapolated_new = 2.0 * primal_new - primal;
    relative_update = c5qda_max_abs_vec(primal_new - primal) /
      std::max(1.0, c5qda_max_abs_vec(primal));
    primal = primal_new;
    dual = dual_new;
    extrapolated = extrapolated_new;

    const arma::vec residual = a * primal - target;
    feasible_maxnorm = c5qda_max_abs_vec(residual);
    primal_violation = std::max(0.0, feasible_maxnorm - lambda);
    const arma::vec adjoint_dual = a.t() * dual;
    dual_violation = std::max(0.0, c5qda_max_abs_vec(adjoint_dual) - 1.0);
    stationarity = c5qda_l1_stationarity(primal, adjoint_dual, tolerance);
    primal_objective = arma::accu(arma::abs(primal));
    dual_objective = -arma::dot(target, dual) -
      lambda * arma::accu(arma::abs(dual));
    duality_gap = primal_objective - dual_objective;
    relative_gap = std::abs(duality_gap) / std::max(
      1.0, std::max(std::abs(primal_objective), std::abs(dual_objective))
    );
    iterations = iteration;
    if (primal_violation <= tolerance && dual_violation <= tolerance &&
        stationarity <= tolerance && relative_gap <= tolerance) {
      converged = true;
      break;
    }
  }
  return Rcpp::List::create(
    Rcpp::Named("solution") = primal,
    Rcpp::Named("dual") = dual,
    Rcpp::Named("converged") = converged,
    Rcpp::Named("iterations") = iterations,
    Rcpp::Named("operator_norm_bound") = operator_bound,
    Rcpp::Named("relative_update") = relative_update,
    Rcpp::Named("feasible_maxnorm") = feasible_maxnorm,
    Rcpp::Named("primal_violation") = primal_violation,
    Rcpp::Named("dual_violation") = dual_violation,
    Rcpp::Named("stationarity_residual") = stationarity,
    Rcpp::Named("primal_objective") = primal_objective,
    Rcpp::Named("dual_objective") = dual_objective,
    Rcpp::Named("duality_gap") = duality_gap,
    Rcpp::Named("relative_gap") = relative_gap,
    Rcpp::Named("failure_stage") = converged ? "" : "vector_dantzig_solver"
  );
}


// [[Rcpp::export]]
Rcpp::List cpp_c5qda_ssqda_moments(const arma::mat& x,
                                   const arma::rowvec& center,
                                   const double zero_tolerance) {
  if (x.n_rows < 3 || x.n_cols == 0 || center.n_elem != x.n_cols ||
      !x.is_finite() || !center.is_finite()) {
    Rcpp::stop("`x` must be finite with at least three rows, and `center` must be conformable.");
  }
  if (!std::isfinite(zero_tolerance) || zero_tolerance < 0.0) {
    Rcpp::stop("`zero_tolerance` must be finite and non-negative.");
  }
  const arma::uword n = x.n_rows;
  const arma::uword p = x.n_cols;
  arma::mat signs(n, p, arma::fill::zeros);
  arma::uword zero_residuals = 0;
  for (arma::uword i = 0; i < n; ++i) {
    long double squared_norm = 0.0L;
    for (arma::uword j = 0; j < p; ++j) {
      const long double value = static_cast<long double>(x(i, j)) -
        static_cast<long double>(center(j));
      squared_norm += value * value;
    }
    const long double radius = std::sqrt(squared_norm);
    if (radius <= static_cast<long double>(zero_tolerance)) {
      ++zero_residuals;
      continue;
    }
    for (arma::uword j = 0; j < p; ++j) {
      signs(i, j) = static_cast<double>(
        (static_cast<long double>(x(i, j)) -
         static_cast<long double>(center(j))) / radius
      );
    }
  }
  const arma::mat spatial_sign_covariance = signs.t() * signs /
    static_cast<double>(n);

  const arma::rowvec sample_mean = arma::mean(x, 0);
  long double centered_sum_squares = 0.0L;
  for (arma::uword i = 0; i < n; ++i) {
    for (arma::uword j = 0; j < p; ++j) {
      const long double value = static_cast<long double>(x(i, j)) -
        static_cast<long double>(sample_mean(j));
      centered_sum_squares += value * value;
    }
  }
  const long double trace_long = centered_sum_squares /
    static_cast<long double>(n - 1);
  if (!std::isfinite(trace_long) ||
      trace_long > static_cast<long double>(std::numeric_limits<double>::max())) {
    Rcpp::stop("The SSQDA trace U-statistic is outside the finite double range.");
  }
  const double trace_estimate = static_cast<double>(trace_long);
  const arma::mat covariance = trace_estimate * spatial_sign_covariance;
  return Rcpp::List::create(
    Rcpp::Named("valid") = zero_residuals == 0 &&
      spatial_sign_covariance.is_finite() && covariance.is_finite(),
    Rcpp::Named("failure_stage") = zero_residuals == 0 ? "" :
      "ssqda_zero_spatial_residual",
    Rcpp::Named("center") = center,
    Rcpp::Named("sample_mean") = sample_mean,
    Rcpp::Named("signs") = signs,
    Rcpp::Named("sscm") = spatial_sign_covariance,
    Rcpp::Named("trace_estimate") = trace_estimate,
    Rcpp::Named("covariance") = covariance,
    Rcpp::Named("zero_residuals") = static_cast<double>(zero_residuals),
    Rcpp::Named("trace_algorithm") = "O(np) triple-U identity"
  );
}


// [[Rcpp::export]]
Rcpp::List cpp_c5qda_signed_logdet(const arma::mat& x) {
  if (x.n_rows == 0 || x.n_rows != x.n_cols || !x.is_finite()) {
    Rcpp::stop("`x` must be a finite non-empty square matrix.");
  }
  double log_determinant;
  double sign;
  const bool valid = c5qda_logdet_positive(x, log_determinant, sign);
  return Rcpp::List::create(
    Rcpp::Named("valid") = valid,
    Rcpp::Named("sign") = sign,
    Rcpp::Named("log_determinant") = log_determinant,
    Rcpp::Named("failure_stage") = valid ? "" : "nonpositive_determinant"
  );
}
