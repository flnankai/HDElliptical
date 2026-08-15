// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp17)]]

#include <RcppArmadillo.h>
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <string>
#include <vector>

namespace {

struct Ch3tegDimensions {
  std::vector<arma::uword> value;
  arma::uword product;
};

void ch3teg_validate_finite_matrix(const arma::mat& x,
                                   const std::string& name) {
  if (!x.is_finite()) {
    Rcpp::stop("`%s` must contain only finite values.", name.c_str());
  }
}

double ch3teg_max_abs(const arma::mat& x) {
  return x.is_empty() ? 0.0 : arma::abs(x).max();
}

Ch3tegDimensions ch3teg_dimensions(const Rcpp::IntegerVector& dimensions,
                                   const arma::uword expected_product) {
  if (dimensions.size() < 1) {
    Rcpp::stop("`dims` must contain at least one mode dimension.");
  }
  std::uint64_t product = 1;
  std::vector<arma::uword> value(dimensions.size());
  for (R_xlen_t k = 0; k < dimensions.size(); ++k) {
    if (dimensions[k] == NA_INTEGER || dimensions[k] < 1) {
      Rcpp::stop("Every entry of `dims` must be a positive integer.");
    }
    const std::uint64_t dimension =
      static_cast<std::uint64_t>(dimensions[k]);
    if (product > std::numeric_limits<std::uint64_t>::max() / dimension) {
      Rcpp::stop("The product of `dims` overflows the index type.");
    }
    product *= dimension;
    value[k] = static_cast<arma::uword>(dimension);
  }
  if (product != static_cast<std::uint64_t>(expected_product)) {
    Rcpp::stop("The product of `dims` must equal `ncol(signs)`.");
  }
  Ch3tegDimensions result;
  result.value = value;
  result.product = expected_product;
  return result;
}

arma::mat ch3teg_unfold(const arma::rowvec& vectorized,
                        const Ch3tegDimensions& dimensions,
                        const arma::uword mode) {
  const arma::uword rows = dimensions.value[mode];
  const arma::uword columns = dimensions.product / rows;
  arma::mat result(rows, columns, arma::fill::zeros);
  std::vector<arma::uword> coordinate(dimensions.value.size());

  for (arma::uword index = 0; index < dimensions.product; ++index) {
    arma::uword remainder = index;
    for (arma::uword k = 0; k < dimensions.value.size(); ++k) {
      coordinate[k] = remainder % dimensions.value[k];
      remainder /= dimensions.value[k];
    }
    arma::uword column = 0;
    arma::uword multiplier = 1;
    for (arma::uword k = 0; k < dimensions.value.size(); ++k) {
      if (k == mode) continue;
      column += coordinate[k] * multiplier;
      multiplier *= dimensions.value[k];
    }
    result(coordinate[mode], column) = vectorized[index];
  }
  return result;
}

arma::rowvec ch3teg_mode_product(const arma::rowvec& input,
                                 const Ch3tegDimensions& dimensions,
                                 const arma::uword mode,
                                 const arma::mat& multiplier) {
  const arma::uword mode_dimension = dimensions.value[mode];
  arma::uword stride = 1;
  for (arma::uword k = 0; k < mode; ++k) {
    stride *= dimensions.value[k];
  }
  const arma::uword block = stride * mode_dimension;
  const arma::uword outer_blocks = dimensions.product / block;
  arma::rowvec output(dimensions.product, arma::fill::zeros);

  for (arma::uword outer = 0; outer < outer_blocks; ++outer) {
    const arma::uword base = outer * block;
    for (arma::uword offset = 0; offset < stride; ++offset) {
      for (arma::uword row = 0; row < mode_dimension; ++row) {
        long double value = 0.0L;
        for (arma::uword column = 0; column < mode_dimension; ++column) {
          value += static_cast<long double>(multiplier(row, column)) *
            static_cast<long double>(input[base + column * stride + offset]);
        }
        const double converted = static_cast<double>(value);
        if (!std::isfinite(converted)) {
          Rcpp::stop("A mode product overflowed; no rescaling repair was applied.");
        }
        output[base + row * stride + offset] = converted;
      }
    }
  }
  return output;
}

double ch3teg_logdet_spd(const arma::mat& x) {
  arma::mat upper;
  if (!arma::chol(upper, x)) {
    return std::numeric_limits<double>::quiet_NaN();
  }
  const arma::vec diagonal = upper.diag();
  if (diagonal.min() <= 0.0 || !diagonal.is_finite()) {
    return std::numeric_limits<double>::quiet_NaN();
  }
  return 2.0 * arma::accu(arma::log(diagonal));
}

double ch3teg_offdiag_l1(const arma::mat& x) {
  double result = 0.0;
  for (arma::uword i = 0; i < x.n_rows; ++i) {
    for (arma::uword j = 0; j < x.n_cols; ++j) {
      if (i != j) result += std::abs(x(i, j));
    }
  }
  return result;
}

double ch3teg_objective(const arma::mat& scatter,
                        const arma::mat& precision,
                        const double penalty) {
  const double logdet = ch3teg_logdet_spd(precision);
  if (!std::isfinite(logdet)) {
    return std::numeric_limits<double>::infinity();
  }
  return arma::trace(scatter * precision) - logdet +
    penalty * ch3teg_offdiag_l1(precision);
}

arma::mat ch3teg_offdiag_prox(const arma::mat& input,
                              const double threshold) {
  arma::mat result = 0.5 * (input + input.t());
  for (arma::uword i = 0; i < result.n_rows; ++i) {
    for (arma::uword j = i + 1; j < result.n_cols; ++j) {
      const double value = result(i, j);
      double thresholded = 0.0;
      if (value > threshold) {
        thresholded = value - threshold;
      } else if (value < -threshold) {
        thresholded = value + threshold;
      }
      result(i, j) = thresholded;
      result(j, i) = thresholded;
    }
  }
  return result;
}

Rcpp::List ch3teg_kkt(const arma::mat& scatter,
                      const arma::mat& precision,
                      const double penalty) {
  arma::mat inverse;
  if (!arma::inv_sympd(inverse, precision)) {
    return Rcpp::List::create(
      Rcpp::Named("inverse") = arma::mat(
        precision.n_rows, precision.n_cols, arma::fill::value(NA_REAL)
      ),
      Rcpp::Named("diagonal") = R_PosInf,
      Rcpp::Named("off_diagonal") = R_PosInf,
      Rcpp::Named("maximum") = R_PosInf
    );
  }
  const arma::mat gradient = scatter - inverse;
  double diagonal = 0.0;
  double off_diagonal = 0.0;
  for (arma::uword i = 0; i < precision.n_rows; ++i) {
    diagonal = std::max(diagonal, std::abs(gradient(i, i)));
    for (arma::uword j = i + 1; j < precision.n_cols; ++j) {
      double residual = 0.0;
      if (precision(i, j) > 0.0) {
        residual = std::abs(gradient(i, j) + penalty);
      } else if (precision(i, j) < 0.0) {
        residual = std::abs(gradient(i, j) - penalty);
      } else {
        residual = std::max(0.0, std::abs(gradient(i, j)) - penalty);
      }
      off_diagonal = std::max(off_diagonal, residual);
    }
  }
  return Rcpp::List::create(
    Rcpp::Named("inverse") = inverse,
    Rcpp::Named("diagonal") = diagonal,
    Rcpp::Named("off_diagonal") = off_diagonal,
    Rcpp::Named("maximum") = std::max(diagonal, off_diagonal)
  );
}

Rcpp::List ch3teg_solver_result(const arma::mat& solution,
                                const arma::mat& scatter,
                                const double penalty,
                                const int iterations,
                                const bool converged,
                                const bool backtracking_failed,
                                const double relative_update,
                                const double accepted_step,
                                const int total_backtracking,
                                const std::vector<double>& history,
                                const bool objective_descent) {
  const Rcpp::List kkt = ch3teg_kkt(scatter, solution, penalty);
  arma::vec eigenvalues;
  const bool eigen_ok = arma::eig_sym(eigenvalues, solution);
  const double minimum_eigenvalue = eigen_ok ?
    eigenvalues.min() : NA_REAL;
  return Rcpp::List::create(
    Rcpp::Named("solution") = solution,
    Rcpp::Named("inverse") = kkt["inverse"],
    Rcpp::Named("objective") = ch3teg_objective(
      scatter, solution, penalty
    ),
    Rcpp::Named("objective_history") = Rcpp::wrap(history),
    Rcpp::Named("objective_descent") = objective_descent,
    Rcpp::Named("iterations") = iterations,
    Rcpp::Named("converged") = converged,
    Rcpp::Named("backtracking_failed") = backtracking_failed,
    Rcpp::Named("relative_update") = relative_update,
    Rcpp::Named("kkt_diagonal") = kkt["diagonal"],
    Rcpp::Named("kkt_off_diagonal") = kkt["off_diagonal"],
    Rcpp::Named("kkt_maximum") = kkt["maximum"],
    Rcpp::Named("minimum_eigenvalue") = minimum_eigenvalue,
    Rcpp::Named("accepted_step") = accepted_step,
    Rcpp::Named("total_backtracking") = total_backtracking
  );
}

} // namespace


// [[Rcpp::export]]
Rcpp::List cpp_ch3teg_mode_crossproducts(
    const arma::mat& signs,
    const Rcpp::IntegerVector& dims) {
  ch3teg_validate_finite_matrix(signs, "signs");
  if (signs.n_rows < 1 || signs.n_cols < 1) {
    Rcpp::stop("`signs` must be a non-empty matrix.");
  }
  const Ch3tegDimensions dimensions = ch3teg_dimensions(dims, signs.n_cols);
  Rcpp::List result(dimensions.value.size());
  const long double sample_size =
    static_cast<long double>(signs.n_rows);

  for (arma::uword mode = 0; mode < dimensions.value.size(); ++mode) {
    arma::mat crossproduct(
      dimensions.value[mode], dimensions.value[mode], arma::fill::zeros
    );
    for (arma::uword observation = 0;
         observation < signs.n_rows; ++observation) {
      const arma::mat unfolding = ch3teg_unfold(
        signs.row(observation), dimensions, mode
      );
      crossproduct += unfolding * unfolding.t();
    }
    crossproduct *= static_cast<double>(
      static_cast<long double>(dimensions.value[mode]) / sample_size
    );
    result[mode] = 0.5 * (crossproduct + crossproduct.t());
  }
  return result;
}


// [[Rcpp::export]]
Rcpp::List cpp_ch3teg_whitened_mode_scatter(
    const arma::mat& signs,
    const Rcpp::IntegerVector& dims,
    const Rcpp::List& pilot_sqrt,
    const bool keep_vectors = false) {
  ch3teg_validate_finite_matrix(signs, "signs");
  const Ch3tegDimensions dimensions = ch3teg_dimensions(dims, signs.n_cols);
  if (pilot_sqrt.size() != static_cast<R_xlen_t>(dimensions.value.size())) {
    Rcpp::stop("`pilot_sqrt` must contain one matrix per tensor mode.");
  }
  std::vector<arma::mat> square_roots(dimensions.value.size());
  for (arma::uword mode = 0; mode < dimensions.value.size(); ++mode) {
    square_roots[mode] = Rcpp::as<arma::mat>(pilot_sqrt[mode]);
    ch3teg_validate_finite_matrix(square_roots[mode], "pilot_sqrt[[k]]");
    if (square_roots[mode].n_rows != dimensions.value[mode] ||
        square_roots[mode].n_cols != dimensions.value[mode]) {
      Rcpp::stop("Every pilot square root must match its mode dimension.");
    }
    const double symmetry_error = ch3teg_max_abs(
      square_roots[mode] - square_roots[mode].t()
    );
    if (symmetry_error > 1e-12 * std::max(
          1.0, ch3teg_max_abs(square_roots[mode]))) {
      Rcpp::stop("Every pilot square root must be symmetric.");
    }
    arma::mat upper;
    if (!arma::chol(upper, square_roots[mode])) {
      Rcpp::stop(
        "Every pilot square root must be strictly positive definite; "
        "no eigenvalue repair was applied."
      );
    }
  }

  Rcpp::List scatters(dimensions.value.size());
  Rcpp::List vectors(keep_vectors ? dimensions.value.size() : 0);
  for (arma::uword target = 0; target < dimensions.value.size(); ++target) {
    arma::mat whitened = signs;
    for (arma::uword observation = 0;
         observation < signs.n_rows; ++observation) {
      arma::rowvec value = signs.row(observation);
      for (arma::uword mode = 0; mode < dimensions.value.size(); ++mode) {
        if (mode != target) {
          value = ch3teg_mode_product(
            value, dimensions, mode, square_roots[mode]
          );
        }
      }
      whitened.row(observation) = value;
    }

    arma::mat scatter(
      dimensions.value[target], dimensions.value[target], arma::fill::zeros
    );
    for (arma::uword observation = 0;
         observation < signs.n_rows; ++observation) {
      const arma::mat unfolding = ch3teg_unfold(
        whitened.row(observation), dimensions, target
      );
      scatter += unfolding * unfolding.t();
    }
    scatter *= static_cast<double>(
      static_cast<long double>(dimensions.value[target]) /
      static_cast<long double>(signs.n_rows)
    );
    scatters[target] = 0.5 * (scatter + scatter.t());
    if (keep_vectors) vectors[target] = whitened;
  }
  return Rcpp::List::create(
    Rcpp::Named("scatter") = scatters,
    Rcpp::Named("vectors") = vectors
  );
}


// [[Rcpp::export]]
Rcpp::List cpp_ch3teg_offdiag_glasso(
    const arma::mat& scatter,
    const double penalty,
    const double tolerance,
    const int max_iterations,
    const double initial_step,
    const int max_backtracking) {
  ch3teg_validate_finite_matrix(scatter, "scatter");
  if (scatter.n_rows < 1 || scatter.n_rows != scatter.n_cols) {
    Rcpp::stop("`scatter` must be a non-empty square matrix.");
  }
  const double symmetry_error = ch3teg_max_abs(scatter - scatter.t());
  if (symmetry_error > 1e-12 * std::max(1.0, ch3teg_max_abs(scatter))) {
    Rcpp::stop("`scatter` must be symmetric.");
  }
  if (!std::isfinite(penalty) || penalty < 0.0 ||
      !std::isfinite(tolerance) || tolerance <= 0.0 ||
      max_iterations < 1 || !std::isfinite(initial_step) ||
      initial_step <= 0.0 || max_backtracking < 1) {
    Rcpp::stop("Invalid off-diagonal graphical-lasso solver controls.");
  }

  const arma::uword p = scatter.n_rows;
  arma::mat current(p, p, arma::fill::zeros);
  for (arma::uword i = 0; i < p; ++i) {
    const double diagonal = scatter(i, i);
    if (!std::isfinite(diagonal) || diagonal <= 0.0 ||
        !std::isfinite(1.0 / diagonal)) {
      Rcpp::stop(
        "The unpenalized diagonal objective is not coercive because a "
        "scatter diagonal is non-positive or numerically non-invertible; "
        "no ridge was added."
      );
    }
    current(i, i) = 1.0 / diagonal;
  }

  if (penalty == 0.0) {
    arma::mat exact;
    if (!arma::inv_sympd(exact, scatter)) {
      Rcpp::stop(
        "With zero penalty the scatter matrix must be strictly positive "
        "definite; no pseudoinverse or ridge was used."
      );
    }
    const Rcpp::List kkt = ch3teg_kkt(scatter, exact, penalty);
    const double maximum = Rcpp::as<double>(kkt["maximum"]);
    const bool converged = std::isfinite(maximum) && maximum <= tolerance;
    const std::vector<double> history = {
      ch3teg_objective(scatter, exact, penalty)
    };
    return ch3teg_solver_result(
      exact, scatter, penalty, 0, converged, false, 0.0,
      NA_REAL, 0, history, true
    );
  }

  Rcpp::List initial_kkt = ch3teg_kkt(scatter, current, penalty);
  const double initial_residual = Rcpp::as<double>(initial_kkt["maximum"]);
  double current_objective = ch3teg_objective(scatter, current, penalty);
  std::vector<double> history = {current_objective};
  if (std::isfinite(initial_residual) && initial_residual <= tolerance) {
    return ch3teg_solver_result(
      current, scatter, penalty, 0, true, false, 0.0,
      NA_REAL, 0, history, true
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
      ch3teg_logdet_spd(current);
    bool accepted = false;
    arma::mat candidate;
    double candidate_objective = std::numeric_limits<double>::infinity();
    double candidate_step = std::min(initial_step, step * 1.2);

    for (int backtracking = 0;
         backtracking < max_backtracking; ++backtracking) {
      if (!std::isfinite(candidate_step) || candidate_step <= 0.0) break;
      candidate = ch3teg_offdiag_prox(
        current - candidate_step * gradient,
        candidate_step * penalty
      );
      arma::mat upper;
      if (candidate.is_finite() && arma::chol(upper, candidate)) {
        const arma::mat difference = candidate - current;
        const double smooth_candidate = arma::trace(scatter * candidate) -
          ch3teg_logdet_spd(candidate);
        const double majorizer = smooth_current +
          arma::accu(gradient % difference) +
          arma::accu(arma::square(difference)) /
            (2.0 * candidate_step);
        candidate_objective = smooth_candidate +
          penalty * ch3teg_offdiag_l1(candidate);
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
    const double denominator = std::max(1.0, arma::norm(current, "fro"));
    relative_update = arma::norm(difference, "fro") / denominator;
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

    const Rcpp::List kkt = ch3teg_kkt(scatter, current, penalty);
    const double maximum = Rcpp::as<double>(kkt["maximum"]);
    if (std::isfinite(maximum) && maximum <= tolerance &&
        relative_update <= tolerance) {
      converged = true;
      break;
    }
    if ((iteration & 255) == 0) Rcpp::checkUserInterrupt();
  }

  return ch3teg_solver_result(
    current, scatter, penalty, iterations, converged,
    backtracking_failed, relative_update, step,
    total_backtracking, history, objective_descent
  );
}
