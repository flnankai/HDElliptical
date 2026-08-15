// Chapter 5 direct and robust sparse linear classifiers.
//
// The Dantzig kernel solves min ||b||_1 subject to
// ||A b - d||_infinity <= lambda with a Chambolle--Pock iteration.
// The DSDA kernel minimizes n^{-1} ||y-a-Xb||_2^2 + lambda ||b||_1
// by cyclic coordinate descent. No uncertified last iterate is returned as a
// fitted method.

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp17)]]

#include <RcppArmadillo.h>

#include <algorithm>
#include <cmath>
#include <limits>

namespace {

void c5lin_finite_matrix(const arma::mat& x, const char* name) {
  if (!x.is_finite()) {
    Rcpp::stop("%s must contain only finite values.", name);
  }
}

void c5lin_finite_vector(const arma::vec& x, const char* name) {
  if (!x.is_finite()) {
    Rcpp::stop("%s must contain only finite values.", name);
  }
}

double c5lin_max_abs(const arma::vec& x) {
  return x.n_elem == 0 ? 0.0 : arma::abs(x).max();
}

arma::vec c5lin_soft(const arma::vec& x, const double threshold) {
  return arma::sign(x) %
    arma::max(arma::abs(x) - threshold,
              arma::zeros<arma::vec>(x.n_elem));
}

double c5lin_l1_stationarity(const arma::vec& coefficient,
                             const arma::vec& stationarity,
                             const double tolerance) {
  double answer = 0.0;
  for (arma::uword j = 0; j < coefficient.n_elem; ++j) {
    double component;
    if (std::abs(coefficient(j)) > tolerance) {
      component = std::abs(
        stationarity(j) + (coefficient(j) > 0.0 ? 1.0 : -1.0)
      );
    } else {
      component = std::max(0.0, std::abs(stationarity(j)) - 1.0);
    }
    answer = std::max(answer, component);
  }
  return answer;
}

struct C5linDantzig {
  arma::vec solution;
  arma::vec dual;
  int iterations;
  bool converged;
  bool operator_zero;
  bool infeasible_zero_operator;
  double relative_update;
  double primal_violation;
  double stationarity_residual;
  double dual_violation;
  double primal_objective;
  double dual_objective;
  double duality_gap;
  double relative_gap;
  double operator_norm;
  double step;
};

C5linDantzig c5lin_solve_dantzig(
    const arma::mat& a, const arma::vec& target, const double lambda,
    const double tolerance, const int maximum_iterations) {
  const arma::uword p = a.n_cols;
  C5linDantzig out;
  out.solution.zeros(p);
  out.dual.zeros(p);
  out.iterations = 0;
  out.converged = false;
  out.operator_zero = false;
  out.infeasible_zero_operator = false;
  out.relative_update = 0.0;
  out.primal_violation =
    std::max(0.0, c5lin_max_abs(target) - lambda);
  out.stationarity_residual = 0.0;
  out.dual_violation = 0.0;
  out.primal_objective = 0.0;
  out.dual_objective = 0.0;
  out.duality_gap = 0.0;
  out.relative_gap = 0.0;
  out.operator_norm = arma::norm(a, 2);
  out.step = std::numeric_limits<double>::quiet_NaN();

  if (!std::isfinite(out.operator_norm)) {
    Rcpp::stop("The Dantzig operator norm is not finite.");
  }
  if (out.operator_norm == 0.0) {
    out.operator_zero = true;
    out.infeasible_zero_operator = out.primal_violation > tolerance;
    out.converged = !out.infeasible_zero_operator;
    return out;
  }
  // Zero feasible implies zero is an exact global minimizer.
  if (c5lin_max_abs(target) <= lambda) {
    out.primal_violation = 0.0;
    out.converged = true;
    return out;
  }

  out.step = 0.99 / out.operator_norm;
  if (!std::isfinite(out.step) || out.step <= 0.0) {
    Rcpp::stop("The Dantzig primal-dual step is invalid.");
  }
  arma::vec primal(p, arma::fill::zeros);
  arma::vec extrapolated(p, arma::fill::zeros);
  arma::vec dual(p, arma::fill::zeros);

  for (int iteration = 1;
       iteration <= maximum_iterations; ++iteration) {
    const arma::vec dual_input =
      dual + out.step * (a * extrapolated);
    arma::vec projected = dual_input / out.step;
    for (arma::uword j = 0; j < p; ++j) {
      projected(j) = std::max(
        target(j) - lambda,
        std::min(target(j) + lambda, projected(j))
      );
    }
    const arma::vec dual_new =
      dual_input - out.step * projected;
    const arma::vec primal_input =
      primal - out.step * (a.t() * dual_new);
    const arma::vec primal_new =
      c5lin_soft(primal_input, out.step);
    const arma::vec extrapolated_new =
      2.0 * primal_new - primal;

    out.relative_update = c5lin_max_abs(primal_new - primal) /
      std::max(1.0, c5lin_max_abs(primal));
    primal = primal_new;
    dual = dual_new;
    extrapolated = extrapolated_new;

    const arma::vec residual = a * primal - target;
    out.primal_violation = std::max(
      0.0, c5lin_max_abs(residual) - lambda
    );
    const arma::vec stationarity = a.t() * dual;
    out.dual_violation = std::max(
      0.0, c5lin_max_abs(stationarity) - 1.0
    );
    out.stationarity_residual = c5lin_l1_stationarity(
      primal, stationarity, tolerance
    );
    out.primal_objective = arma::accu(arma::abs(primal));
    out.dual_objective = -arma::dot(target, dual) -
      lambda * arma::accu(arma::abs(dual));
    out.duality_gap =
      out.primal_objective - out.dual_objective;
    const double gap_scale = std::max(
      1.0, std::max(std::abs(out.primal_objective),
                    std::abs(out.dual_objective))
    );
    out.relative_gap = std::abs(out.duality_gap) / gap_scale;
    out.iterations = iteration;

    if (out.primal_violation <= tolerance &&
        out.dual_violation <= tolerance &&
        out.stationarity_residual <= tolerance &&
        out.relative_gap <= tolerance) {
      out.converged = true;
      break;
    }
    if ((iteration & 2047) == 0) {
      Rcpp::checkUserInterrupt();
    }
  }
  out.solution = primal;
  out.dual = dual;
  return out;
}

double c5lin_dsda_kkt(const arma::mat& centered_x,
                      const arma::vec& residual,
                      const arma::vec& coefficient,
                      const double lambda) {
  const double n = static_cast<double>(centered_x.n_rows);
  const arma::vec gradient =
    -(2.0 / n) * centered_x.t() * residual;
  double answer = std::abs(-2.0 * arma::mean(residual));
  for (arma::uword j = 0; j < coefficient.n_elem; ++j) {
    double component;
    if (coefficient(j) > 0.0) {
      component = std::abs(gradient(j) + lambda);
    } else if (coefficient(j) < 0.0) {
      component = std::abs(gradient(j) - lambda);
    } else {
      component =
        std::max(0.0, std::abs(gradient(j)) - lambda);
    }
    answer = std::max(answer, component);
  }
  return answer;
}

} // namespace


// [[Rcpp::export]]
Rcpp::List cpp_c5lin_dantzig(const arma::mat& a,
                             const arma::vec& target,
                             const double lambda,
                             const double tolerance,
                             const int maximum_iterations) {
  c5lin_finite_matrix(a, "a");
  c5lin_finite_vector(target, "target");
  if (a.n_rows == 0 || a.n_rows != a.n_cols ||
      target.n_elem != a.n_cols) {
    Rcpp::stop("a must be square and match target.");
  }
  if (!std::isfinite(lambda) || lambda <= 0.0 ||
      !std::isfinite(tolerance) || tolerance <= 0.0 ||
      maximum_iterations < 1) {
    Rcpp::stop("Invalid Dantzig solver controls.");
  }
  const C5linDantzig fit = c5lin_solve_dantzig(
    a, target, lambda, tolerance, maximum_iterations
  );
  return Rcpp::List::create(
    Rcpp::Named("solution") = fit.solution,
    Rcpp::Named("dual") = fit.dual,
    Rcpp::Named("iterations") = fit.iterations,
    Rcpp::Named("converged") = fit.converged,
    Rcpp::Named("operator_zero") = fit.operator_zero,
    Rcpp::Named("infeasible_zero_operator") =
      fit.infeasible_zero_operator,
    Rcpp::Named("relative_update") = fit.relative_update,
    Rcpp::Named("primal_violation") = fit.primal_violation,
    Rcpp::Named("stationarity_residual") =
      fit.stationarity_residual,
    Rcpp::Named("dual_violation") = fit.dual_violation,
    Rcpp::Named("primal_objective") = fit.primal_objective,
    Rcpp::Named("dual_objective") = fit.dual_objective,
    Rcpp::Named("duality_gap") = fit.duality_gap,
    Rcpp::Named("relative_gap") = fit.relative_gap,
    Rcpp::Named("operator_norm") = fit.operator_norm,
    Rcpp::Named("primal_dual_step") = fit.step
  );
}


// [[Rcpp::export]]
Rcpp::List cpp_c5lin_dsda(const arma::mat& x,
                          const arma::vec& response,
                          const double lambda,
                          const double tolerance,
                          const int maximum_iterations) {
  c5lin_finite_matrix(x, "x");
  c5lin_finite_vector(response, "response");
  if (x.n_rows < 2 || x.n_cols == 0 ||
      response.n_elem != x.n_rows) {
    Rcpp::stop("x and response have incompatible dimensions.");
  }
  if (!std::isfinite(lambda) || lambda <= 0.0 ||
      !std::isfinite(tolerance) || tolerance <= 0.0 ||
      maximum_iterations < 1) {
    Rcpp::stop("Invalid DSDA solver controls.");
  }

  const arma::uword n_rows = x.n_rows;
  const arma::uword p = x.n_cols;
  const double n = static_cast<double>(n_rows);
  const arma::rowvec x_mean = arma::mean(x, 0);
  const double response_mean = arma::mean(response);
  const arma::mat centered_x = x.each_row() - x_mean;
  const arma::vec centered_response =
    response - response_mean;
  const arma::vec coordinate_scale =
    arma::sum(arma::square(centered_x), 0).t() / n;

  arma::vec coefficient(p, arma::fill::zeros);
  arma::vec residual = centered_response;
  arma::vec dual(n_rows, arma::fill::zeros);
  bool converged = false;
  int iterations = 0;
  double relative_update = 0.0;
  double kkt_residual =
    std::numeric_limits<double>::infinity();
  double primal_objective =
    std::numeric_limits<double>::infinity();
  double dual_objective =
    -std::numeric_limits<double>::infinity();
  double duality_gap =
    std::numeric_limits<double>::infinity();
  double relative_gap =
    std::numeric_limits<double>::infinity();
  double dual_scale = 0.0;
  double dual_feasibility =
    std::numeric_limits<double>::infinity();

  for (int iteration = 1;
       iteration <= maximum_iterations; ++iteration) {
    double largest_change = 0.0;
    const double coefficient_scale =
      std::max(1.0, c5lin_max_abs(coefficient));
    for (arma::uword j = 0; j < p; ++j) {
      const double old_value = coefficient(j);
      if (coordinate_scale(j) == 0.0) {
        coefficient(j) = 0.0;
      } else {
        const double raw = arma::dot(
          centered_x.col(j),
          residual + centered_x.col(j) * old_value
        ) / n;
        const double thresholded = std::copysign(
          std::max(std::abs(raw) - lambda / 2.0, 0.0),
          raw
        );
        coefficient(j) =
          thresholded / coordinate_scale(j);
      }
      const double change = coefficient(j) - old_value;
      if (change != 0.0) {
        residual -= centered_x.col(j) * change;
      }
      largest_change =
        std::max(largest_change, std::abs(change));
    }
    residual -= arma::mean(residual);
    relative_update =
      largest_change / coefficient_scale;

    kkt_residual = c5lin_dsda_kkt(
      centered_x, residual, coefficient, lambda
    );
    primal_objective =
      arma::dot(residual, residual) / n +
      lambda * arma::accu(arma::abs(coefficient));

    const arma::vec raw_dual = -(2.0 / n) * residual;
    const double raw_constraint = c5lin_max_abs(
      centered_x.t() * raw_dual
    );
    dual_scale =
      raw_constraint > lambda ? lambda / raw_constraint : 1.0;
    dual = dual_scale * raw_dual;
    dual_feasibility = std::max(
      0.0,
      c5lin_max_abs(centered_x.t() * dual) - lambda
    );
    dual_objective =
      -arma::dot(dual, centered_response) -
      (n / 4.0) * arma::dot(dual, dual);
    duality_gap = primal_objective - dual_objective;
    const double gap_scale = std::max(
      1.0, std::max(std::abs(primal_objective),
                    std::abs(dual_objective))
    );
    relative_gap = std::abs(duality_gap) / gap_scale;
    iterations = iteration;

    if (kkt_residual <= tolerance &&
        dual_feasibility <= tolerance &&
        duality_gap >= -tolerance * gap_scale &&
        relative_gap <= tolerance) {
      converged = true;
      break;
    }
    if ((iteration & 2047) == 0) {
      Rcpp::checkUserInterrupt();
    }
  }

  const double regression_intercept =
    response_mean - arma::dot(x_mean, coefficient);
  const arma::vec fitted =
    regression_intercept + x * coefficient;
  return Rcpp::List::create(
    Rcpp::Named("coefficient") = coefficient,
    Rcpp::Named("regression_intercept") =
      regression_intercept,
    Rcpp::Named("fitted") = fitted,
    Rcpp::Named("residual") = response - fitted,
    Rcpp::Named("centered_x") = centered_x,
    Rcpp::Named("x_mean") = x_mean,
    Rcpp::Named("response_mean") = response_mean,
    Rcpp::Named("coordinate_scale") = coordinate_scale,
    Rcpp::Named("converged") = converged,
    Rcpp::Named("iterations") = iterations,
    Rcpp::Named("relative_update") = relative_update,
    Rcpp::Named("kkt_residual") = kkt_residual,
    Rcpp::Named("primal_objective") = primal_objective,
    Rcpp::Named("dual_objective") = dual_objective,
    Rcpp::Named("duality_gap") = duality_gap,
    Rcpp::Named("relative_gap") = relative_gap,
    Rcpp::Named("dual") = dual,
    Rcpp::Named("dual_scale") = dual_scale,
    Rcpp::Named("dual_feasibility") = dual_feasibility
  );
}
