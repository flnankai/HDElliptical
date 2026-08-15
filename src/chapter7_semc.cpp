// SEMC numerical kernels for HDElliptical.
//
// Algorithm provenance: Feng and Zhuang (2026), arXiv:2605.08995, and the
// authors' MIT-licensed GEMcluster implementation at commit
// 10fce04fe690fe274dd5d237cfcd3d5c6a4139f6.  This file is an independent
// rewrite of the mathematical and software contracts; it contains no code
// from the GPL-licensed huge package.
//
// Upstream reference implementation copyright (c) 2026 Dan Zhuang and
// Long Feng. HDElliptical adaptation copyright (c) 2026 HDElliptical
// contributors. SPDX-License-Identifier: MIT

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp17)]]

#include <RcppArmadillo.h>
#include <cmath>
#include <limits>

namespace {

arma::mat semc_symmetrize(const arma::mat& x) {
  return 0.5 * (x + x.t());
}

double semc_minimum_eigenvalue(const arma::mat& x) {
  arma::vec values;
  if (!arma::eig_sym(values, semc_symmetrize(x)) || values.n_elem == 0) {
    Rcpp::stop("A symmetric eigendecomposition failed.");
  }
  return values.min();
}

void semc_require_square_finite(const arma::mat& x, const char* name) {
  if (x.n_rows == 0 || x.n_rows != x.n_cols) {
    Rcpp::stop("%s must be a nonempty square matrix.", name);
  }
  if (!x.is_finite()) {
    Rcpp::stop("%s must contain only finite values.", name);
  }
}

void semc_require_spd(const arma::mat& x, double spd_tol,
                      const char* name, double* minimum = nullptr) {
  semc_require_square_finite(x, name);
  arma::mat chol_factor;
  const arma::mat sym = semc_symmetrize(x);
  const double min_eigenvalue = semc_minimum_eigenvalue(sym);
  if (!arma::chol(chol_factor, sym) || !std::isfinite(min_eigenvalue) ||
      min_eigenvalue <= spd_tol) {
    Rcpp::stop("%s is not strictly positive definite above spd_tol.", name);
  }
  if (minimum != nullptr) {
    *minimum = min_eigenvalue;
  }
}

double semc_logdet_spd(const arma::mat& x) {
  arma::mat chol_factor;
  if (!arma::chol(chol_factor, semc_symmetrize(x))) {
    return std::numeric_limits<double>::infinity();
  }
  return 2.0 * arma::sum(arma::log(chol_factor.diag()));
}

double semc_glasso_smooth_objective(const arma::mat& precision,
                                    const arma::mat& scatter) {
  const double logdet = semc_logdet_spd(precision);
  if (!std::isfinite(logdet)) {
    return std::numeric_limits<double>::infinity();
  }
  return arma::accu(scatter % precision) - logdet;
}

double semc_offdiag_l1(const arma::mat& x) {
  double value = 0.0;
  for (arma::uword j = 0; j < x.n_cols; ++j) {
    for (arma::uword i = 0; i < x.n_rows; ++i) {
      if (i != j) {
        value += std::abs(x(i, j));
      }
    }
  }
  return value;
}

arma::mat semc_offdiag_soft_threshold(const arma::mat& x, double threshold) {
  arma::mat out = x;
  for (arma::uword j = 0; j < out.n_cols; ++j) {
    for (arma::uword i = 0; i < out.n_rows; ++i) {
      if (i == j) {
        continue;
      }
      const double value = out(i, j);
      if (value > threshold) {
        out(i, j) = value - threshold;
      } else if (value < -threshold) {
        out(i, j) = value + threshold;
      } else {
        out(i, j) = 0.0;
      }
    }
  }
  return semc_symmetrize(out);
}

double semc_offdiag_glasso_kkt(const arma::mat& precision,
                               const arma::mat& scatter,
                               double lambda) {
  arma::mat inverse;
  if (!arma::inv_sympd(inverse, semc_symmetrize(precision))) {
    return std::numeric_limits<double>::infinity();
  }
  const arma::mat gradient = scatter - inverse;
  double residual = 0.0;
  for (arma::uword j = 0; j < precision.n_cols; ++j) {
    for (arma::uword i = 0; i < precision.n_rows; ++i) {
      double current = 0.0;
      if (i == j) {
        current = std::abs(gradient(i, j));
      } else if (precision(i, j) > 0.0) {
        current = std::abs(gradient(i, j) + lambda);
      } else if (precision(i, j) < 0.0) {
        current = std::abs(gradient(i, j) - lambda);
      } else {
        current = std::max(std::abs(gradient(i, j)) - lambda, 0.0);
      }
      residual = std::max(residual, current);
    }
  }
  return residual;
}

}  // namespace


// [[Rcpp::export]]
arma::mat cpp_ch7_semc_delta(const arma::mat& x,
                             const arma::mat& centers,
                             const arma::mat& precision) {
  if (x.n_rows == 0 || x.n_cols == 0 || !x.is_finite()) {
    Rcpp::stop("x must be a nonempty finite matrix.");
  }
  if (centers.n_rows == 0 || centers.n_cols != x.n_cols ||
      !centers.is_finite()) {
    Rcpp::stop("centers must be a finite K by p matrix matching x.");
  }
  if (precision.n_rows != x.n_cols || precision.n_cols != x.n_cols ||
      !precision.is_finite()) {
    Rcpp::stop("precision must be a finite p by p matrix matching x.");
  }

  arma::mat out(x.n_rows, centers.n_rows, arma::fill::zeros);
  for (arma::uword k = 0; k < centers.n_rows; ++k) {
    arma::mat residual = x.each_row() - centers.row(k);
    out.col(k) = arma::sum((residual * precision) % residual, 1);
  }
  return out;
}


// [[Rcpp::export]]
Rcpp::List cpp_ch7_semc_softmax(const arma::mat& log_scores) {
  if (log_scores.n_rows == 0 || log_scores.n_cols == 0 ||
      !log_scores.is_finite()) {
    Rcpp::stop("log_scores must be a nonempty finite matrix.");
  }
  arma::mat probabilities(log_scores.n_rows, log_scores.n_cols,
                          arma::fill::zeros);
  arma::vec log_normalizers(log_scores.n_rows, arma::fill::zeros);
  for (arma::uword i = 0; i < log_scores.n_rows; ++i) {
    const double maximum = log_scores.row(i).max();
    arma::rowvec shifted = arma::exp(log_scores.row(i) - maximum);
    const double denominator = arma::accu(shifted);
    if (!std::isfinite(denominator) || denominator <= 0.0) {
      Rcpp::stop("A posterior softmax denominator is not strictly positive.");
    }
    probabilities.row(i) = shifted / denominator;
    log_normalizers(i) = maximum + std::log(denominator);
  }
  return Rcpp::List::create(
    Rcpp::Named("probabilities") = probabilities,
    Rcpp::Named("log_normalizers") = log_normalizers
  );
}


// [[Rcpp::export]]
Rcpp::List cpp_ch7_semc_weighted_sign_scatter(
    const arma::mat& residuals,
    const arma::vec& weights,
    double radial_floor) {
  if (residuals.n_rows == 0 || residuals.n_cols == 0 ||
      !residuals.is_finite()) {
    Rcpp::stop("residuals must be a nonempty finite matrix.");
  }
  if (weights.n_elem != residuals.n_rows || !weights.is_finite() ||
      arma::any(weights < 0.0)) {
    Rcpp::stop("weights must be finite, nonnegative, and match residuals.");
  }
  if (!std::isfinite(radial_floor) || radial_floor <= 0.0) {
    Rcpp::stop("radial_floor must be finite and strictly positive.");
  }
  const double weight_sum = arma::accu(weights);
  if (!std::isfinite(weight_sum) || weight_sum <= 0.0) {
    Rcpp::stop("weights must have a strictly positive sum.");
  }

  const arma::vec radii_squared = arma::sum(arma::square(residuals), 1);
  const arma::vec denominators = arma::clamp(
    radii_squared, radial_floor, std::numeric_limits<double>::infinity()
  );
  const arma::vec coefficients = weights / denominators;
  arma::mat scatter = residuals.t() * (residuals.each_col() % coefficients);
  scatter = semc_symmetrize(scatter / weight_sum);

  return Rcpp::List::create(
    Rcpp::Named("estimate") = scatter,
    Rcpp::Named("weight_sum") = weight_sum,
    Rcpp::Named("minimum_radius_squared") = radii_squared.min(),
    Rcpp::Named("radial_floor_uses") =
      static_cast<int>(arma::accu(radii_squared < radial_floor)),
    Rcpp::Named("trace") = arma::trace(scatter)
  );
}


// [[Rcpp::export]]
Rcpp::List cpp_ch7_semc_weighted_tyler(
    const arma::mat& residuals,
    const arma::vec& weights,
    const arma::mat& initial_shape,
    double ridge,
    double radial_floor,
    double tolerance,
    int max_iterations,
    double spd_tol) {
  if (residuals.n_rows == 0 || residuals.n_cols == 0 ||
      !residuals.is_finite()) {
    Rcpp::stop("residuals must be a nonempty finite matrix.");
  }
  if (weights.n_elem != residuals.n_rows || !weights.is_finite() ||
      arma::any(weights < 0.0)) {
    Rcpp::stop("weights must be finite, nonnegative, and match residuals.");
  }
  if (initial_shape.n_rows != residuals.n_cols ||
      initial_shape.n_cols != residuals.n_cols) {
    Rcpp::stop("initial_shape must be p by p.");
  }
  if (!std::isfinite(ridge) || ridge < 0.0 || ridge >= 1.0) {
    Rcpp::stop("ridge must lie in [0, 1).");
  }
  if (!std::isfinite(radial_floor) || radial_floor <= 0.0 ||
      !std::isfinite(tolerance) || tolerance <= 0.0 ||
      max_iterations < 1 || !std::isfinite(spd_tol) || spd_tol < 0.0) {
    Rcpp::stop("Invalid Tyler iteration controls.");
  }

  const double weight_sum = arma::accu(weights);
  if (!std::isfinite(weight_sum) || weight_sum <= 0.0) {
    Rcpp::stop("weights must have a strictly positive sum.");
  }
  const arma::uword p = residuals.n_cols;
  arma::mat shape = semc_symmetrize(initial_shape);
  double minimum_eigenvalue = 0.0;
  semc_require_spd(shape, spd_tol, "initial_shape", &minimum_eigenvalue);
  const double initial_trace = arma::trace(shape);
  if (!std::isfinite(initial_trace) || initial_trace <= 0.0) {
    Rcpp::stop("initial_shape has a nonpositive trace.");
  }
  shape *= static_cast<double>(p) / initial_trace;

  bool converged = false;
  int iterations = 0;
  double relative_change = std::numeric_limits<double>::infinity();
  double minimum_quadratic = std::numeric_limits<double>::infinity();
  double maximum_quadratic = 0.0;
  int floor_uses = 0;

  for (int iteration = 1; iteration <= max_iterations; ++iteration) {
    arma::mat precision;
    if (!arma::inv_sympd(precision, shape)) {
      Rcpp::stop("The Tyler iterate could not be inverted as SPD.");
    }
    arma::vec quadratic = arma::sum((residuals * precision) % residuals, 1);
    if (!quadratic.is_finite() || arma::any(quadratic < 0.0)) {
      Rcpp::stop("The Tyler quadratic forms are not finite and nonnegative.");
    }
    minimum_quadratic = quadratic.min();
    maximum_quadratic = quadratic.max();
    floor_uses = static_cast<int>(arma::accu(quadratic < radial_floor));
    arma::vec denominators = arma::clamp(
      quadratic, radial_floor, std::numeric_limits<double>::infinity()
    );
    arma::vec coefficients =
      static_cast<double>(p) * weights / (denominators * weight_sum);
    arma::mat proposal =
      residuals.t() * (residuals.each_col() % coefficients);
    proposal = (1.0 - ridge) * semc_symmetrize(proposal) +
      ridge * arma::eye<arma::mat>(p, p);
    const double proposal_trace = arma::trace(proposal);
    if (!std::isfinite(proposal_trace) || proposal_trace <= 0.0) {
      Rcpp::stop("The Tyler proposal has a nonpositive trace.");
    }
    proposal *= static_cast<double>(p) / proposal_trace;
    semc_require_spd(proposal, spd_tol, "Tyler proposal",
                     &minimum_eigenvalue);
    relative_change = arma::norm(proposal - shape, "fro") /
      std::max(1.0, arma::norm(shape, "fro"));
    shape = proposal;
    iterations = iteration;
    if (relative_change <= tolerance) {
      converged = true;
      break;
    }
  }

  arma::mat precision;
  if (!arma::inv_sympd(precision, shape)) {
    Rcpp::stop("The final Tyler shape could not be inverted as SPD.");
  }
  semc_require_spd(shape, spd_tol, "final Tyler shape",
                   &minimum_eigenvalue);
  return Rcpp::List::create(
    Rcpp::Named("shape") = shape,
    Rcpp::Named("precision") = precision,
    Rcpp::Named("converged") = converged,
    Rcpp::Named("iterations") = iterations,
    Rcpp::Named("relative_change") = relative_change,
    Rcpp::Named("minimum_eigenvalue") = minimum_eigenvalue,
    Rcpp::Named("trace") = arma::trace(shape),
    Rcpp::Named("minimum_quadratic") = minimum_quadratic,
    Rcpp::Named("maximum_quadratic") = maximum_quadratic,
    Rcpp::Named("radial_floor_uses") = floor_uses,
    Rcpp::Named("ridge") = ridge
  );
}


// [[Rcpp::export]]
Rcpp::List cpp_ch7_semc_offdiag_glasso(
    const arma::mat& scatter,
    double lambda,
    double tolerance,
    int max_iterations,
    double initial_step,
    int max_backtracking,
    double spd_tol,
    double majorization_tolerance) {
  semc_require_square_finite(scatter, "scatter");
  if (!std::isfinite(lambda) || lambda <= 0.0 ||
      !std::isfinite(tolerance) || tolerance <= 0.0 ||
      max_iterations < 1 || !std::isfinite(initial_step) ||
      initial_step <= 0.0 || max_backtracking < 1 ||
      !std::isfinite(spd_tol) || spd_tol < 0.0 ||
      !std::isfinite(majorization_tolerance) ||
      majorization_tolerance < 0.0) {
    Rcpp::stop("Invalid graphical-lasso controls.");
  }
  const arma::mat symmetric_scatter = semc_symmetrize(scatter);
  const double mean_diagonal = arma::mean(symmetric_scatter.diag());
  if (!std::isfinite(mean_diagonal) || mean_diagonal <= 0.0) {
    Rcpp::stop("scatter must have a strictly positive mean diagonal.");
  }

  const arma::uword p = scatter.n_rows;
  arma::mat precision =
    arma::eye<arma::mat>(p, p) / mean_diagonal;
  double minimum_eigenvalue = 0.0;
  semc_require_spd(precision, spd_tol, "initial graphical-lasso precision",
                   &minimum_eigenvalue);
  bool converged = false;
  bool backtracking_failed = false;
  int iterations = 0;
  int total_backtracking = 0;
  double accepted_step = 0.0;
  double relative_update = std::numeric_limits<double>::infinity();
  double kkt_residual = semc_offdiag_glasso_kkt(
    precision, symmetric_scatter, lambda
  );

  for (int iteration = 1; iteration <= max_iterations; ++iteration) {
    arma::mat inverse;
    if (!arma::inv_sympd(inverse, precision)) {
      Rcpp::stop("The graphical-lasso iterate could not be inverted as SPD.");
    }
    const arma::mat gradient = symmetric_scatter - inverse;
    const double smooth = semc_glasso_smooth_objective(
      precision, symmetric_scatter
    );
    double step = initial_step;
    arma::mat candidate;
    bool accepted = false;
    int used_backtracking = 0;
    for (int backtrack = 0; backtrack < max_backtracking; ++backtrack) {
      candidate = semc_offdiag_soft_threshold(
        precision - step * gradient, step * lambda
      );
      arma::mat chol_factor;
      const double candidate_minimum = semc_minimum_eigenvalue(candidate);
      if (candidate_minimum > spd_tol &&
          arma::chol(chol_factor, candidate)) {
        const arma::mat difference = candidate - precision;
        const double bound = smooth + arma::accu(gradient % difference) +
          arma::accu(arma::square(difference)) / (2.0 * step) +
          majorization_tolerance;
        const double candidate_smooth = semc_glasso_smooth_objective(
          candidate, symmetric_scatter
        );
        if (std::isfinite(candidate_smooth) && candidate_smooth <= bound) {
          accepted = true;
          minimum_eigenvalue = candidate_minimum;
          used_backtracking = backtrack;
          break;
        }
      }
      step *= 0.5;
    }
    total_backtracking += used_backtracking;
    iterations = iteration;
    if (!accepted) {
      backtracking_failed = true;
      break;
    }
    relative_update = arma::norm(candidate - precision, "fro") /
      std::max(1.0, arma::norm(precision, "fro"));
    precision = candidate;
    accepted_step = step;
    kkt_residual = semc_offdiag_glasso_kkt(
      precision, symmetric_scatter, lambda
    );
    if (std::isfinite(kkt_residual) && kkt_residual <= tolerance) {
      converged = true;
      break;
    }
  }

  arma::mat inverse;
  if (!arma::inv_sympd(inverse, precision)) {
    Rcpp::stop("The final graphical-lasso precision is not invertible as SPD.");
  }
  semc_require_spd(precision, spd_tol,
                   "final graphical-lasso precision", &minimum_eigenvalue);
  const double objective = semc_glasso_smooth_objective(
    precision, symmetric_scatter
  ) + lambda * semc_offdiag_l1(precision);
  return Rcpp::List::create(
    Rcpp::Named("solution") = precision,
    Rcpp::Named("inverse") = inverse,
    Rcpp::Named("converged") = converged,
    Rcpp::Named("backtracking_failed") = backtracking_failed,
    Rcpp::Named("iterations") = iterations,
    Rcpp::Named("relative_update") = relative_update,
    Rcpp::Named("kkt_residual") = kkt_residual,
    Rcpp::Named("objective") = objective,
    Rcpp::Named("minimum_eigenvalue") = minimum_eigenvalue,
    Rcpp::Named("accepted_step") = accepted_step,
    Rcpp::Named("total_backtracking") = total_backtracking,
    Rcpp::Named("diagonal_penalty") = false
  );
}


// [[Rcpp::export]]
arma::vec cpp_ch7_semc_weighted_kde(
    const arma::vec& evaluation_grid,
    const arma::vec& observations,
    const arma::vec& weights,
    double bandwidth) {
  if (evaluation_grid.n_elem < 2 || !evaluation_grid.is_finite() ||
      observations.n_elem == 0 || !observations.is_finite() ||
      weights.n_elem != observations.n_elem || !weights.is_finite() ||
      arma::any(weights < 0.0) || !std::isfinite(bandwidth) ||
      bandwidth <= 0.0) {
    Rcpp::stop("Invalid weighted-KDE inputs.");
  }
  const double weight_sum = arma::accu(weights);
  if (!std::isfinite(weight_sum) || weight_sum <= 0.0) {
    Rcpp::stop("KDE weights must have a strictly positive sum.");
  }
  const double normalizer =
    1.0 / (std::sqrt(2.0 * arma::datum::pi) * bandwidth * weight_sum);
  arma::vec density(evaluation_grid.n_elem, arma::fill::zeros);
  for (arma::uword i = 0; i < evaluation_grid.n_elem; ++i) {
    const arma::vec standardized =
      (evaluation_grid(i) - observations) / bandwidth;
    density(i) = normalizer *
      arma::accu(weights % arma::exp(-0.5 * arma::square(standardized)));
  }
  return density;
}


