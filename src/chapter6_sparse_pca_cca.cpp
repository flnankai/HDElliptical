// Chapter 6: sparse PCA and sparse CCA numerical kernels.
//
// These routines implement only the displayed optimization problems and
// deterministic solvers.  They never add a ridge, replace a zero norm, floor
// an eigenvalue, perturb a tie, or invent a tuning parameter.

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp17)]]

#include <RcppArmadillo.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <string>
#include <vector>

namespace {

double c6spc_max_abs(const arma::vec& x) {
  return x.n_elem == 0 ? 0.0 : arma::abs(x).max();
}


double c6spc_sign_distance(const arma::vec& x, const arma::vec& y) {
  return std::min(arma::norm(x - y, 2), arma::norm(x + y, 2));
}


void c6spc_anchor(arma::vec& x) {
  if (x.n_elem == 0) return;
  const double maximum = arma::abs(x).max();
  for (arma::uword j = 0; j < x.n_elem; ++j) {
    if (std::abs(x(j)) == maximum) {
      if (x(j) < 0.0) x *= -1.0;
      return;
    }
  }
}


void c6spc_check_matrix(const arma::mat& x, const char* name,
                        const arma::uword minimum_rows = 1) {
  if (x.n_rows < minimum_rows || x.n_cols == 0 || !x.is_finite()) {
    Rcpp::stop("`%s` must be a finite matrix with at least %d row(s) and one column.",
               name, static_cast<int>(minimum_rows));
  }
}


arma::vec c6spc_soft(const arma::vec& x, const double threshold) {
  arma::vec answer(x.n_elem, arma::fill::zeros);
  for (arma::uword j = 0; j < x.n_elem; ++j) {
    const double magnitude = std::abs(x(j));
    if (magnitude > threshold) {
      answer(j) = std::copysign(magnitude - threshold, x(j));
    }
  }
  return answer;
}


arma::mat c6spc_soft(const arma::mat& x, const double threshold) {
  arma::mat answer(x.n_rows, x.n_cols, arma::fill::zeros);
  for (arma::uword j = 0; j < x.n_cols; ++j) {
    for (arma::uword i = 0; i < x.n_rows; ++i) {
      const double magnitude = std::abs(x(i, j));
      if (magnitude > threshold) {
        answer(i, j) = std::copysign(magnitude - threshold, x(i, j));
      }
    }
  }
  return answer;
}


struct SparseUpdate {
  arma::vec vector;
  double lambda;
  double gamma;
  double l1_norm;
  double l2_norm;
  double constraint_residual;
  double kkt_residual;
  double scaled_kkt_residual;
  double ratio_residual;
  bool active;
  bool tie_branch;
  bool valid;
  std::string message;
};


double c6spc_sparse_kkt(const arma::vec& a, const arma::vec& v,
                        const double lambda, const double gamma) {
  double residual = 0.0;
  for (arma::uword j = 0; j < a.n_elem; ++j) {
    double component;
    if (v(j) != 0.0) {
      component = std::abs(
        a(j) - lambda * (v(j) > 0.0 ? 1.0 : -1.0) - gamma * v(j)
      );
    } else {
      component = std::max(0.0, std::abs(a(j)) - lambda);
    }
    residual = std::max(residual, component);
  }
  return residual;
}


SparseUpdate c6spc_sparse_update(const arma::vec& a, const double bound,
                                 const double tolerance) {
  SparseUpdate out;
  out.vector = arma::vec(a.n_elem, arma::fill::zeros);
  out.lambda = NA_REAL;
  out.gamma = NA_REAL;
  out.l1_norm = NA_REAL;
  out.l2_norm = NA_REAL;
  out.constraint_residual = std::numeric_limits<double>::infinity();
  out.kkt_residual = std::numeric_limits<double>::infinity();
  out.scaled_kkt_residual = std::numeric_limits<double>::infinity();
  out.ratio_residual = std::numeric_limits<double>::infinity();
  out.active = false;
  out.tie_branch = false;
  out.valid = false;
  out.message = "unknown sparse update failure";

  if (a.n_elem == 0 || !a.is_finite() || !std::isfinite(bound) ||
      bound < 1.0 || bound > std::sqrt(static_cast<double>(a.n_elem)) ||
      !std::isfinite(tolerance) || tolerance <= 0.0) {
    out.message = "invalid sparse-update input";
    return out;
  }
  const double max_abs = c6spc_max_abs(a);
  if (!(max_abs > 0.0) || !std::isfinite(max_abs)) {
    out.message = "the sparse-update score vector is zero";
    return out;
  }
  const double norm_a = arma::norm(a, 2);
  const double unconstrained_ratio = arma::norm(a, 1) / norm_a;

  if (unconstrained_ratio <= bound) {
    out.vector = a / norm_a;
    out.lambda = 0.0;
    out.gamma = norm_a;
    out.active = false;
  } else {
    std::vector<arma::uword> tied;
    tied.reserve(a.n_elem);
    for (arma::uword j = 0; j < a.n_elem; ++j) {
      if (std::abs(a(j)) == max_abs) tied.push_back(j);
    }
    const double tie_ratio = std::sqrt(static_cast<double>(tied.size()));
    if (bound <= tie_ratio) {
      // At a maximum-absolute-value tie, the soft-threshold ratio jumps.
      // Construct the deterministic KKT optimum on the smallest indices.
      out.tie_branch = true;
      out.active = true;
      out.lambda = max_abs;
      out.gamma = 0.0;
      arma::uword support = 1;
      while (support < tied.size() &&
             std::sqrt(static_cast<double>(support)) < bound) {
        ++support;
      }
      if (support == 1) {
        out.vector(tied[0]) = a(tied[0]) > 0.0 ? 1.0 : -1.0;
      } else {
        const double m = static_cast<double>(support - 1);
        double radicand = m * (static_cast<double>(support) - bound * bound);
        if (radicand < 0.0 &&
            std::abs(radicand) <= 32.0 * std::numeric_limits<double>::epsilon()) {
          radicand = 0.0;
        }
        if (radicand < 0.0) {
          out.message = "the deterministic tie optimum is outside double precision";
          return out;
        }
        const double large =
          (m * bound + std::sqrt(radicand)) /
          (m * static_cast<double>(support));
        const double small = bound - m * large;
        for (arma::uword k = 0; k < support - 1; ++k) {
          out.vector(tied[k]) = std::copysign(large, a(tied[k]));
        }
        out.vector(tied[support - 1]) =
          std::copysign(small, a(tied[support - 1]));
        const double constructed_norm = arma::norm(out.vector, 2);
        if (!(constructed_norm > 0.0) || !std::isfinite(constructed_norm)) {
          out.message = "the deterministic tie optimum has zero norm";
          return out;
        }
        out.vector /= constructed_norm;
      }
    } else {
      out.active = true;
      double lower = 0.0;
      double upper = std::nextafter(max_abs, 0.0);
      arma::vec candidate;
      for (int iteration = 0; iteration < 200; ++iteration) {
        const double middle = lower + (upper - lower) / 2.0;
        candidate = c6spc_soft(a, middle);
        const double norm_candidate = arma::norm(candidate, 2);
        if (!(norm_candidate > 0.0)) {
          upper = middle;
          continue;
        }
        const double ratio = arma::norm(candidate, 1) / norm_candidate;
        if (ratio > bound) {
          lower = middle;
        } else {
          upper = middle;
        }
      }
      out.lambda = lower + (upper - lower) / 2.0;
      candidate = c6spc_soft(a, out.lambda);
      out.gamma = arma::norm(candidate, 2);
      if (!(out.gamma > 0.0) || !std::isfinite(out.gamma)) {
        out.message = "soft thresholding produced a zero vector";
        return out;
      }
      out.vector = candidate / out.gamma;
    }
  }

  out.l1_norm = arma::norm(out.vector, 1);
  out.l2_norm = arma::norm(out.vector, 2);
  out.constraint_residual = std::max(
    std::abs(out.l2_norm - 1.0), std::max(0.0, out.l1_norm - bound)
  );
  out.ratio_residual = out.active ? std::abs(out.l1_norm - bound) :
    std::max(0.0, out.l1_norm - bound);
  out.kkt_residual = c6spc_sparse_kkt(
    a, out.vector, out.lambda, out.gamma
  );
  const double scale = std::max(
    1.0, std::max(max_abs, std::max(std::abs(out.lambda),
                                   std::abs(out.gamma)))
  );
  out.scaled_kkt_residual = out.kkt_residual / scale;
  out.valid = out.vector.is_finite() && std::isfinite(out.l1_norm) &&
    std::isfinite(out.l2_norm) && std::isfinite(out.kkt_residual) &&
    std::isfinite(out.scaled_kkt_residual);
  out.message = out.valid ? "ok" : "non-finite sparse-update certificate";
  return out;
}


bool c6spc_hard_step(const arma::mat& matrix, const arma::vec& vector,
                     const arma::uword sparsity, arma::vec& answer,
                     arma::uvec& support) {
  const arma::vec product = matrix * vector;
  if (!product.is_finite()) return false;
  std::vector<arma::uword> indices(product.n_elem);
  std::iota(indices.begin(), indices.end(), 0);
  std::stable_sort(indices.begin(), indices.end(),
    [&product](const arma::uword left, const arma::uword right) {
      const double a = std::abs(product(left));
      const double b = std::abs(product(right));
      if (a != b) return a > b;
      return left < right;
    }
  );
  answer.zeros(product.n_elem);
  support.set_size(sparsity);
  for (arma::uword k = 0; k < sparsity; ++k) {
    support(k) = indices[k];
    answer(indices[k]) = product(indices[k]);
  }
  const double norm = arma::norm(answer, 2);
  if (!(norm > 0.0) || !std::isfinite(norm)) return false;
  answer /= norm;
  return answer.is_finite();
}


struct FantopeProjection {
  arma::mat matrix;
  arma::vec eigenvalues;
  double trace_residual;
  bool valid;
};


FantopeProjection c6spc_fantope_projection(const arma::mat& input,
                                           const arma::uword rank) {
  FantopeProjection out;
  out.valid = false;
  out.trace_residual = std::numeric_limits<double>::infinity();
  arma::vec values;
  arma::mat vectors;
  const arma::mat symmetric = 0.5 * (input + input.t());
  if (!arma::eig_sym(values, vectors, symmetric)) return out;
  if (rank == input.n_rows) {
    out.eigenvalues = arma::vec(input.n_rows, arma::fill::ones);
    out.matrix = arma::eye<arma::mat>(input.n_rows, input.n_cols);
    out.trace_residual = 0.0;
    out.valid = true;
    return out;
  }
  double lower = values.min() - 1.0;
  double upper = values.max();
  arma::vec projected(values.n_elem);
  for (int iteration = 0; iteration < 250; ++iteration) {
    const double theta = lower + (upper - lower) / 2.0;
    for (arma::uword j = 0; j < values.n_elem; ++j) {
      projected(j) = std::min(1.0, std::max(0.0, values(j) - theta));
    }
    if (arma::accu(projected) > static_cast<double>(rank)) {
      lower = theta;
    } else {
      upper = theta;
    }
  }
  const double theta = lower + (upper - lower) / 2.0;
  for (arma::uword j = 0; j < values.n_elem; ++j) {
    projected(j) = std::min(1.0, std::max(0.0, values(j) - theta));
  }
  out.matrix = vectors * arma::diagmat(projected) * vectors.t();
  out.matrix = 0.5 * (out.matrix + out.matrix.t());
  out.eigenvalues = projected;
  out.trace_residual = std::abs(arma::accu(projected) -
                                static_cast<double>(rank));
  out.valid = out.matrix.is_finite() && projected.is_finite() &&
    std::isfinite(out.trace_residual);
  return out;
}


arma::vec c6spc_cross_apply(const arma::mat& x, const arma::mat& y,
                            const arma::mat& left, const arma::mat& right,
                            const arma::vec& values, const arma::vec& v) {
  arma::vec answer = x.t() * (y * v);
  if (values.n_elem > 0) {
    answer -= left * (values % (right.t() * v));
  }
  return answer;
}


arma::vec c6spc_cross_adjoint(const arma::mat& x, const arma::mat& y,
                              const arma::mat& left, const arma::mat& right,
                              const arma::vec& values, const arma::vec& u) {
  arma::vec answer = y.t() * (x * u);
  if (values.n_elem > 0) {
    answer -= right * (values % (left.t() * u));
  }
  return answer;
}


arma::vec c6spc_cross_start(const arma::mat& x, const arma::mat& y,
                            const arma::mat& left, const arma::mat& right,
                            const arma::vec& values) {
  const arma::uword q = y.n_cols;
  arma::uword best = 0;
  double best_norm = -1.0;
  arma::vec basis(q, arma::fill::zeros);
  for (arma::uword j = 0; j < q; ++j) {
    basis.zeros();
    basis(j) = 1.0;
    const arma::vec column = c6spc_cross_apply(
      x, y, left, right, values, basis
    );
    const double column_norm = arma::norm(column, 2);
    if (column_norm > best_norm) {
      best_norm = column_norm;
      best = j;
    }
  }
  arma::vec start(q, arma::fill::zeros);
  if (best_norm > 0.0 && std::isfinite(best_norm)) start(best) = 1.0;
  return start;
}

}  // namespace


// [[Rcpp::export]]
Rcpp::List cpp_ch6spc_tpm(const arma::mat& operator_matrix,
                          const arma::uvec& sparsity,
                          const arma::mat& initial,
                          const double tolerance,
                          const int maximum_iterations) {
  if (operator_matrix.n_rows == 0 ||
      operator_matrix.n_rows != operator_matrix.n_cols ||
      !operator_matrix.is_finite()) {
    Rcpp::stop("`operator` must be a finite non-empty square matrix.");
  }
  if (sparsity.n_elem == 0 || arma::any(sparsity < 1) ||
      arma::any(sparsity > operator_matrix.n_rows)) {
    Rcpp::stop("Every `sparsity` value must lie between one and p.");
  }
  if (!std::isfinite(tolerance) || tolerance <= 0.0 ||
      maximum_iterations < 1) {
    Rcpp::stop("Invalid truncated-power iteration controls.");
  }
  const bool supplied_initial = initial.n_elem > 0;
  if (supplied_initial &&
      (initial.n_rows != operator_matrix.n_rows ||
       initial.n_cols != sparsity.n_elem || !initial.is_finite())) {
    Rcpp::stop("`initial` must be empty or p by length(sparsity).");
  }

  const arma::uword p = operator_matrix.n_rows;
  const arma::uword components = sparsity.n_elem;
  arma::mat current = operator_matrix;
  arma::mat loadings(p, components, arma::fill::zeros);
  arma::vec values(components, arma::fill::zeros);
  Rcpp::List certificates(components);
  arma::uword completed = 0;
  int failure_component = NA_INTEGER;
  std::string failure_message = "none";

  for (arma::uword component = 0; component < components; ++component) {
    arma::vec vector;
    if (supplied_initial) {
      vector = initial.col(component);
      const double norm = arma::norm(vector, 2);
      if (!(norm > 0.0) || !std::isfinite(norm)) {
        failure_component = static_cast<int>(component + 1);
        failure_message = "the supplied initial vector has zero norm";
        break;
      }
      vector /= norm;
    } else {
      arma::vec eigenvalues;
      arma::mat eigenvectors;
      if (!arma::eig_sym(eigenvalues, eigenvectors, current) ||
          eigenvalues.n_elem == 0) {
        failure_component = static_cast<int>(component + 1);
        failure_message = "the deterministic leading-eigenvector initialization failed";
        break;
      }
      vector = eigenvectors.col(eigenvectors.n_cols - 1);
    }
    c6spc_anchor(vector);

    bool stopped = false;
    bool zero_step = false;
    int iterations = 0;
    double change = std::numeric_limits<double>::infinity();
    double previous_objective = std::numeric_limits<double>::quiet_NaN();
    double maximum_objective_decrease = 0.0;
    arma::uvec support;

    for (int iteration = 1; iteration <= maximum_iterations; ++iteration) {
      iterations = iteration;
      arma::vec candidate;
      if (!c6spc_hard_step(current, vector, sparsity(component),
                           candidate, support)) {
        zero_step = true;
        break;
      }
      change = c6spc_sign_distance(candidate, vector);
      const double objective = arma::as_scalar(candidate.t() * current * candidate);
      if (std::isfinite(previous_objective)) {
        maximum_objective_decrease = std::max(
          maximum_objective_decrease, previous_objective - objective
        );
      }
      vector = candidate;
      previous_objective = objective;
      if (change <= tolerance) {
        stopped = true;
        break;
      }
    }

    arma::vec fixed_candidate;
    arma::uvec fixed_support;
    const bool fixed_valid = !zero_step && c6spc_hard_step(
      current, vector, sparsity(component), fixed_candidate, fixed_support
    );
    const double fixed_residual = fixed_valid ?
      c6spc_sign_distance(fixed_candidate, vector) :
      std::numeric_limits<double>::infinity();
    const double unit_residual = std::abs(arma::norm(vector, 2) - 1.0);
    const arma::uword support_size = arma::accu(arma::abs(vector) > 0.0);
    const double support_violation = support_size > sparsity(component) ?
      static_cast<double>(support_size - sparsity(component)) : 0.0;
    const arma::vec product = current * vector;
    const double rayleigh = arma::dot(vector, product);
    double restricted_residual = 0.0;
    for (arma::uword j = 0; j < p; ++j) {
      if (vector(j) != 0.0) {
        restricted_residual = std::max(
          restricted_residual, std::abs(product(j) - rayleigh * vector(j))
        );
      }
    }
    const double objective_scale = std::max(1.0, std::abs(rayleigh));
    const bool converged = stopped && fixed_valid &&
      fixed_residual <= tolerance &&
      unit_residual <= 10.0 * tolerance && support_violation == 0.0 &&
      maximum_objective_decrease <= 10.0 * tolerance * objective_scale;
    const bool certified = converged && vector.is_finite() &&
      std::isfinite(rayleigh) && std::isfinite(restricted_residual);

    certificates[component] = Rcpp::List::create(
      Rcpp::Named("converged") = converged,
      Rcpp::Named("certified") = certified,
      Rcpp::Named("iterations") = iterations,
      Rcpp::Named("last_change") = change,
      Rcpp::Named("fixed_point_residual") = fixed_residual,
      Rcpp::Named("unit_norm_residual") = unit_residual,
      Rcpp::Named("support_size") = static_cast<int>(support_size),
      Rcpp::Named("support_bound") = static_cast<int>(sparsity(component)),
      Rcpp::Named("support_violation") = support_violation,
      Rcpp::Named("restricted_eigen_residual") = restricted_residual,
      Rcpp::Named("deflated_rayleigh_value") = rayleigh,
      Rcpp::Named("maximum_objective_decrease") = maximum_objective_decrease,
      Rcpp::Named("zero_matrix_vector_step") = zero_step
    );
    if (!certified) {
      failure_component = static_cast<int>(component + 1);
      failure_message = zero_step ?
        "the truncated matrix-vector product is zero" :
        "the iteration limit or fixed-point certificate was not satisfied";
      break;
    }

    c6spc_anchor(vector);
    loadings.col(component) = vector;
    values(component) = rayleigh;
    ++completed;
    const arma::mat projection = arma::eye<arma::mat>(p, p) -
      vector * vector.t();
    current = projection * current * projection;
    current = 0.5 * (current + current.t());
  }

  return Rcpp::List::create(
    Rcpp::Named("converged") = completed == components,
    Rcpp::Named("completed") = static_cast<int>(completed),
    Rcpp::Named("failure_component") = failure_component,
    Rcpp::Named("failure_message") = failure_message,
    Rcpp::Named("loadings") = loadings,
    Rcpp::Named("values") = values,
    Rcpp::Named("certificates") = certificates,
    Rcpp::Named("deflated_operator") = current
  );
}


// [[Rcpp::export]]
Rcpp::List cpp_ch6spc_fantope(const arma::mat& covariance,
                              const int rank,
                              const double tau,
                              const double rho,
                              const double tolerance,
                              const int maximum_iterations) {
  if (covariance.n_rows == 0 || covariance.n_rows != covariance.n_cols ||
      !covariance.is_finite()) {
    Rcpp::stop("`covariance` must be a finite non-empty square matrix.");
  }
  const arma::uword p = covariance.n_rows;
  if (rank < 1 || static_cast<arma::uword>(rank) > p ||
      !std::isfinite(tau) || tau < 0.0 ||
      !std::isfinite(rho) || rho <= 0.0 ||
      !std::isfinite(tolerance) || tolerance <= 0.0 ||
      maximum_iterations < 1) {
    Rcpp::stop("Invalid Fantope ADMM controls.");
  }

  arma::mat H = (static_cast<double>(rank) / static_cast<double>(p)) *
    arma::eye<arma::mat>(p, p);
  arma::mat Z = H;
  arma::mat U(p, p, arma::fill::zeros);
  int iterations = 0;
  double primal = std::numeric_limits<double>::infinity();
  double dual = std::numeric_limits<double>::infinity();
  double scaled_primal = std::numeric_limits<double>::infinity();
  double scaled_dual = std::numeric_limits<double>::infinity();
  double trace_residual = std::numeric_limits<double>::infinity();
  bool residual_converged = false;

  for (int iteration = 1; iteration <= maximum_iterations; ++iteration) {
    iterations = iteration;
    FantopeProjection projection = c6spc_fantope_projection(
      Z - U + covariance / rho, static_cast<arma::uword>(rank)
    );
    if (!projection.valid) break;
    H = projection.matrix;
    trace_residual = projection.trace_residual;
    const arma::mat previous_Z = Z;
    const arma::mat soft_argument = H + U;
    Z = c6spc_soft(soft_argument, tau / rho);
    Z = 0.5 * (Z + Z.t());
    U += H - Z;
    U = 0.5 * (U + U.t());
    primal = arma::norm(H - Z, "fro");
    dual = rho * arma::norm(Z - previous_Z, "fro");
    scaled_primal = primal / std::max(
      1.0, std::max(arma::norm(H, "fro"), arma::norm(Z, "fro"))
    );
    scaled_dual = dual / std::max(1.0, rho * arma::norm(U, "fro"));
    if (scaled_primal <= tolerance && scaled_dual <= tolerance) {
      residual_converged = true;
      break;
    }
  }

  arma::vec eigenvalues;
  arma::eig_sym(eigenvalues, H);
  const double minimum = eigenvalues.n_elem ? eigenvalues.min() : NA_REAL;
  const double maximum = eigenvalues.n_elem ? eigenvalues.max() : NA_REAL;
  const double lower_violation = std::max(0.0, -minimum);
  const double upper_violation = std::max(0.0, maximum - 1.0);
  trace_residual = std::abs(arma::trace(H) - static_cast<double>(rank));
  const double constraint_tolerance = std::max(
    10.0 * tolerance, 128.0 * std::numeric_limits<double>::epsilon() *
      static_cast<double>(p)
  );
  const bool certified = residual_converged && H.is_finite() && Z.is_finite() &&
    std::isfinite(primal) && std::isfinite(dual) &&
    trace_residual <= constraint_tolerance &&
    lower_violation <= constraint_tolerance &&
    upper_violation <= constraint_tolerance;
  const double objective = arma::accu(covariance % H) -
    tau * arma::accu(arma::abs(H));

  return Rcpp::List::create(
    Rcpp::Named("converged") = residual_converged,
    Rcpp::Named("certified") = certified,
    Rcpp::Named("iterations") = iterations,
    Rcpp::Named("H") = H,
    Rcpp::Named("Z") = Z,
    Rcpp::Named("scaled_dual_variable") = U,
    Rcpp::Named("primal_residual") = primal,
    Rcpp::Named("dual_residual") = dual,
    Rcpp::Named("scaled_primal_residual") = scaled_primal,
    Rcpp::Named("scaled_dual_residual") = scaled_dual,
    Rcpp::Named("trace_residual") = trace_residual,
    Rcpp::Named("minimum_eigenvalue") = minimum,
    Rcpp::Named("maximum_eigenvalue") = maximum,
    Rcpp::Named("lower_violation") = lower_violation,
    Rcpp::Named("upper_violation") = upper_violation,
    Rcpp::Named("objective") = objective
  );
}


// [[Rcpp::export]]
Rcpp::List cpp_ch6spc_pmd_pca(const arma::mat& x,
                              const arma::vec& l1_bounds,
                              const arma::mat& initial,
                              const double tolerance,
                              const int maximum_iterations) {
  c6spc_check_matrix(x, "x", 2);
  const arma::uword n = x.n_rows;
  const arma::uword p = x.n_cols;
  const arma::uword components = l1_bounds.n_elem;
  if (components == 0 || components > std::min(n, p) ||
      !l1_bounds.is_finite() || arma::any(l1_bounds < 1.0) ||
      arma::any(l1_bounds > std::sqrt(static_cast<double>(p))) ||
      !std::isfinite(tolerance) || tolerance <= 0.0 ||
      maximum_iterations < 1) {
    Rcpp::stop("Invalid PMD sparse-PCA inputs or controls.");
  }
  const bool supplied_initial = initial.n_elem > 0;
  if (supplied_initial &&
      (initial.n_rows != p || initial.n_cols != components ||
       !initial.is_finite())) {
    Rcpp::stop("`initial` must be empty or p by length(l1_bounds).");
  }

  arma::mat residual = x;
  arma::mat loadings(p, components, arma::fill::zeros);
  arma::mat left_vectors(n, components, arma::fill::zeros);
  arma::vec singular_values(components, arma::fill::zeros);
  Rcpp::List certificates(components);
  arma::uword completed = 0;
  int failure_component = NA_INTEGER;
  std::string failure_message = "none";

  for (arma::uword component = 0; component < components; ++component) {
    arma::vec v;
    if (supplied_initial) {
      v = initial.col(component);
      const double norm = arma::norm(v, 2);
      if (!(norm > 0.0) || !std::isfinite(norm)) {
        failure_component = static_cast<int>(component + 1);
        failure_message = "the supplied initial loading has zero norm";
        break;
      }
      v /= norm;
    } else {
      arma::mat init_u;
      arma::vec init_d;
      arma::mat init_v;
      if (!arma::svd_econ(init_u, init_d, init_v, residual) ||
          init_d.n_elem == 0 || !(init_d(0) > 0.0)) {
        failure_component = static_cast<int>(component + 1);
        failure_message = "the deterministic leading-SVD initialization failed";
        break;
      }
      v = init_v.col(0);
    }
    c6spc_anchor(v);

    int iterations = 0;
    bool stopped = false;
    bool zero_step = false;
    double last_change = std::numeric_limits<double>::infinity();
    double previous_objective = -std::numeric_limits<double>::infinity();
    double maximum_objective_decrease = 0.0;
    SparseUpdate last_update = c6spc_sparse_update(
      arma::vec(residual.n_cols, arma::fill::zeros),
      l1_bounds(component), tolerance
    );
    arma::vec u;

    for (int iteration = 1; iteration <= maximum_iterations; ++iteration) {
      iterations = iteration;
      const arma::vec xv = residual * v;
      const double xv_norm = arma::norm(xv, 2);
      if (!(xv_norm > 0.0) || !std::isfinite(xv_norm)) {
        zero_step = true;
        break;
      }
      u = xv / xv_norm;
      last_update = c6spc_sparse_update(
        residual.t() * u, l1_bounds(component), tolerance
      );
      if (!last_update.valid) {
        zero_step = true;
        break;
      }
      const arma::vec candidate = last_update.vector;
      last_change = c6spc_sign_distance(candidate, v);
      const double objective = arma::norm(residual * candidate, 2);
      if (std::isfinite(previous_objective)) {
        maximum_objective_decrease = std::max(
          maximum_objective_decrease, previous_objective - objective
        );
      }
      previous_objective = objective;
      v = candidate;
      if (last_change <= tolerance) {
        stopped = true;
        break;
      }
    }

    const arma::vec xv = residual * v;
    const double d = arma::norm(xv, 2);
    if (!(d > 0.0) || !std::isfinite(d)) zero_step = true;
    if (!zero_step) u = xv / d;
    const SparseUpdate fixed = zero_step ? last_update :
      c6spc_sparse_update(residual.t() * u, l1_bounds(component), tolerance);
    const double fixed_residual = (!zero_step && fixed.valid) ?
      c6spc_sign_distance(fixed.vector, v) :
      std::numeric_limits<double>::infinity();
    const arma::vec u_stationarity = residual * v - d * u;
    const double u_kkt = zero_step ? std::numeric_limits<double>::infinity() :
      c6spc_max_abs(u_stationarity);
    const double scaled_u_kkt = u_kkt / std::max(1.0, d);
    const double v_kkt = (!zero_step && fixed.valid) ?
      c6spc_sparse_kkt(residual.t() * u, v, fixed.lambda, fixed.gamma) :
      std::numeric_limits<double>::infinity();
    const arma::vec v_score = residual.t() * u;
    const double v_scale = (!zero_step) ?
      std::max(1.0, c6spc_max_abs(v_score)) : 1.0;
    const double scaled_v_kkt = v_kkt / v_scale;
    const double l1_norm = arma::norm(v, 1);
    const double unit_residual = std::abs(arma::norm(v, 2) - 1.0);
    const double constraint_residual = std::max(
      unit_residual, std::max(0.0, l1_norm - l1_bounds(component))
    );
    const double objective_scale = std::max(1.0, d);
    const bool converged = stopped && !zero_step && fixed.valid &&
      fixed_residual <= tolerance &&
      fixed.ratio_residual <= 10.0 * tolerance &&
      scaled_u_kkt <= 10.0 * tolerance &&
      scaled_v_kkt <= 10.0 * tolerance &&
      constraint_residual <= 10.0 * tolerance &&
      maximum_objective_decrease <= 10.0 * tolerance * objective_scale;
    const bool certified = converged && v.is_finite() && u.is_finite();

    certificates[component] = Rcpp::List::create(
      Rcpp::Named("converged") = converged,
      Rcpp::Named("certified") = certified,
      Rcpp::Named("iterations") = iterations,
      Rcpp::Named("last_change") = last_change,
      Rcpp::Named("fixed_point_residual") = fixed_residual,
      Rcpp::Named("u_kkt_residual") = u_kkt,
      Rcpp::Named("scaled_u_kkt_residual") = scaled_u_kkt,
      Rcpp::Named("v_kkt_residual") = v_kkt,
      Rcpp::Named("scaled_v_kkt_residual") = scaled_v_kkt,
      Rcpp::Named("unit_norm_residual") = unit_residual,
      Rcpp::Named("l1_norm") = l1_norm,
      Rcpp::Named("l1_bound") = l1_bounds(component),
      Rcpp::Named("constraint_residual") = constraint_residual,
      Rcpp::Named("bound_residual") = fixed.ratio_residual,
      Rcpp::Named("lambda") = fixed.lambda,
      Rcpp::Named("active_l1_constraint") = fixed.active,
      Rcpp::Named("maximum_tie_branch") = fixed.tie_branch,
      Rcpp::Named("singular_value") = d,
      Rcpp::Named("maximum_objective_decrease") = maximum_objective_decrease,
      Rcpp::Named("zero_block_step") = zero_step
    );
    if (!certified) {
      failure_component = static_cast<int>(component + 1);
      failure_message = zero_step ?
        "a PMD block matrix-vector product is zero" :
        "the iteration limit or block-KKT certificate was not satisfied";
      break;
    }

    c6spc_anchor(v);
    // Preserve d u v' when anchoring v.
    const arma::vec original_xv = residual * v;
    u = original_xv / arma::norm(original_xv, 2);
    loadings.col(component) = v;
    left_vectors.col(component) = u;
    singular_values(component) = d;
    residual -= d * u * v.t();
    ++completed;
  }

  return Rcpp::List::create(
    Rcpp::Named("converged") = completed == components,
    Rcpp::Named("completed") = static_cast<int>(completed),
    Rcpp::Named("failure_component") = failure_component,
    Rcpp::Named("failure_message") = failure_message,
    Rcpp::Named("loadings") = loadings,
    Rcpp::Named("left_vectors") = left_vectors,
    Rcpp::Named("singular_values") = singular_values,
    Rcpp::Named("certificates") = certificates,
    Rcpp::Named("residual_matrix") = residual
  );
}


// [[Rcpp::export]]
Rcpp::List cpp_ch6spc_pmd_cca(const arma::mat& x,
                              const arma::mat& y,
                              const arma::vec& l1_x,
                              const arma::vec& l1_y,
                              const arma::mat& initial_x,
                              const arma::mat& initial_y,
                              const double tolerance,
                              const int maximum_iterations) {
  c6spc_check_matrix(x, "x", 2);
  c6spc_check_matrix(y, "y", 2);
  if (x.n_rows != y.n_rows) {
    Rcpp::stop("`x` and `y` must have the same number of rows.");
  }
  const arma::uword px = x.n_cols;
  const arma::uword py = y.n_cols;
  const arma::uword components = l1_x.n_elem;
  if (components == 0 || l1_y.n_elem != components ||
      components > std::min(x.n_rows, std::min(px, py)) ||
      !l1_x.is_finite() || !l1_y.is_finite() ||
      arma::any(l1_x < 1.0) ||
      arma::any(l1_x > std::sqrt(static_cast<double>(px))) ||
      arma::any(l1_y < 1.0) ||
      arma::any(l1_y > std::sqrt(static_cast<double>(py))) ||
      !std::isfinite(tolerance) || tolerance <= 0.0 ||
      maximum_iterations < 1) {
    Rcpp::stop("Invalid PMD sparse-CCA inputs or controls.");
  }
  const bool supplied_initial = initial_x.n_elem > 0 || initial_y.n_elem > 0;
  if (supplied_initial &&
      (initial_x.n_rows != px || initial_x.n_cols != components ||
       initial_y.n_rows != py || initial_y.n_cols != components ||
       !initial_x.is_finite() || !initial_y.is_finite())) {
    Rcpp::stop("Initial CCA matrices must both be supplied with conformable dimensions.");
  }

  arma::mat left(px, 0);
  arma::mat right(py, 0);
  arma::vec deflation_values;
  arma::mat x_loadings(px, components, arma::fill::zeros);
  arma::mat y_loadings(py, components, arma::fill::zeros);
  arma::vec singular_values(components, arma::fill::zeros);
  Rcpp::List certificates(components);
  arma::uword completed = 0;
  int failure_component = NA_INTEGER;
  std::string failure_message = "none";

  for (arma::uword component = 0; component < components; ++component) {
    arma::vec v;
    arma::vec supplied_u;
    if (supplied_initial) {
      supplied_u = initial_x.col(component);
      v = initial_y.col(component);
      const double norm_u = arma::norm(supplied_u, 2);
      const double norm_v = arma::norm(v, 2);
      if (!(norm_u > 0.0) || !(norm_v > 0.0) ||
          !std::isfinite(norm_u) || !std::isfinite(norm_v)) {
        failure_component = static_cast<int>(component + 1);
        failure_message = "a supplied CCA initial vector has zero norm";
        break;
      }
      supplied_u /= norm_u;
      v /= norm_v;
      if (arma::dot(c6spc_cross_apply(x, y, left, right,
                                      deflation_values, v), supplied_u) < 0.0) {
        v *= -1.0;
      }
    } else {
      v = c6spc_cross_start(x, y, left, right, deflation_values);
      if (arma::norm(v, 2) == 0.0) {
        failure_component = static_cast<int>(component + 1);
        failure_message = "the deflated cross-covariance operator is zero";
        break;
      }
    }

    arma::vec u(px, arma::fill::zeros);
    int iterations = 0;
    bool stopped = false;
    bool zero_step = false;
    double last_change = std::numeric_limits<double>::infinity();
    double previous_objective = -std::numeric_limits<double>::infinity();
    double maximum_objective_decrease = 0.0;
    SparseUpdate update_u = c6spc_sparse_update(
      arma::vec(px, arma::fill::zeros), l1_x(component), tolerance
    );
    SparseUpdate update_v = c6spc_sparse_update(
      arma::vec(py, arma::fill::zeros), l1_y(component), tolerance
    );

    for (int iteration = 1; iteration <= maximum_iterations; ++iteration) {
      iterations = iteration;
      update_u = c6spc_sparse_update(
        c6spc_cross_apply(x, y, left, right, deflation_values, v),
        l1_x(component), tolerance
      );
      if (!update_u.valid) {
        zero_step = true;
        break;
      }
      arma::vec candidate_u = update_u.vector;
      update_v = c6spc_sparse_update(
        c6spc_cross_adjoint(x, y, left, right, deflation_values,
                            candidate_u),
        l1_y(component), tolerance
      );
      if (!update_v.valid) {
        zero_step = true;
        break;
      }
      arma::vec candidate_v = update_v.vector;
      if (u.n_elem == candidate_u.n_elem && arma::norm(u, 2) > 0.0) {
        const double common_orientation = arma::dot(candidate_u, u) +
          arma::dot(candidate_v, v);
        if (common_orientation < 0.0) {
          candidate_u *= -1.0;
          candidate_v *= -1.0;
        }
        last_change = std::max(
          c6spc_sign_distance(candidate_u, u),
          c6spc_sign_distance(candidate_v, v)
        );
      }
      u = candidate_u;
      v = candidate_v;
      const arma::vec mv = c6spc_cross_apply(
        x, y, left, right, deflation_values, v
      );
      const double objective = arma::dot(u, mv);
      if (std::isfinite(previous_objective)) {
        maximum_objective_decrease = std::max(
          maximum_objective_decrease, previous_objective - objective
        );
      }
      previous_objective = objective;
      if (iteration > 1 && last_change <= tolerance) {
        stopped = true;
        break;
      }
    }

    const arma::vec mv = c6spc_cross_apply(
      x, y, left, right, deflation_values, v
    );
    const arma::vec mtu = c6spc_cross_adjoint(
      x, y, left, right, deflation_values, u
    );
    const SparseUpdate fixed_u = zero_step ? update_u :
      c6spc_sparse_update(mv, l1_x(component), tolerance);
    const SparseUpdate fixed_v = zero_step ? update_v :
      c6spc_sparse_update(mtu, l1_y(component), tolerance);
    const double fixed_u_residual = (!zero_step && fixed_u.valid) ?
      c6spc_sign_distance(fixed_u.vector, u) :
      std::numeric_limits<double>::infinity();
    const double fixed_v_residual = (!zero_step && fixed_v.valid) ?
      c6spc_sign_distance(fixed_v.vector, v) :
      std::numeric_limits<double>::infinity();
    const double fixed_residual = std::max(
      fixed_u_residual, fixed_v_residual
    );
    const double u_kkt = (!zero_step && fixed_u.valid) ?
      c6spc_sparse_kkt(mv, u, fixed_u.lambda, fixed_u.gamma) :
      std::numeric_limits<double>::infinity();
    const double v_kkt = (!zero_step && fixed_v.valid) ?
      c6spc_sparse_kkt(mtu, v, fixed_v.lambda, fixed_v.gamma) :
      std::numeric_limits<double>::infinity();
    const double scaled_u_kkt = u_kkt / std::max(1.0, c6spc_max_abs(mv));
    const double scaled_v_kkt = v_kkt / std::max(1.0, c6spc_max_abs(mtu));
    const double u_constraint = std::max(
      std::abs(arma::norm(u, 2) - 1.0),
      std::max(0.0, arma::norm(u, 1) - l1_x(component))
    );
    const double v_constraint = std::max(
      std::abs(arma::norm(v, 2) - 1.0),
      std::max(0.0, arma::norm(v, 1) - l1_y(component))
    );
    const double d = arma::dot(u, mv);
    const double objective_scale = std::max(1.0, std::abs(d));
    const bool converged = stopped && !zero_step && fixed_u.valid &&
      fixed_v.valid && fixed_residual <= tolerance &&
      fixed_u.ratio_residual <= 10.0 * tolerance &&
      fixed_v.ratio_residual <= 10.0 * tolerance &&
      scaled_u_kkt <= 10.0 * tolerance &&
      scaled_v_kkt <= 10.0 * tolerance &&
      u_constraint <= 10.0 * tolerance &&
      v_constraint <= 10.0 * tolerance &&
      maximum_objective_decrease <= 10.0 * tolerance * objective_scale &&
      std::isfinite(d) && d > 0.0;
    const bool certified = converged && u.is_finite() && v.is_finite();

    certificates[component] = Rcpp::List::create(
      Rcpp::Named("converged") = converged,
      Rcpp::Named("certified") = certified,
      Rcpp::Named("iterations") = iterations,
      Rcpp::Named("last_change") = last_change,
      Rcpp::Named("fixed_point_residual") = fixed_residual,
      Rcpp::Named("x_kkt_residual") = u_kkt,
      Rcpp::Named("scaled_x_kkt_residual") = scaled_u_kkt,
      Rcpp::Named("y_kkt_residual") = v_kkt,
      Rcpp::Named("scaled_y_kkt_residual") = scaled_v_kkt,
      Rcpp::Named("x_constraint_residual") = u_constraint,
      Rcpp::Named("y_constraint_residual") = v_constraint,
      Rcpp::Named("x_l1_norm") = arma::norm(u, 1),
      Rcpp::Named("y_l1_norm") = arma::norm(v, 1),
      Rcpp::Named("x_l1_bound") = l1_x(component),
      Rcpp::Named("y_l1_bound") = l1_y(component),
      Rcpp::Named("x_bound_residual") = fixed_u.ratio_residual,
      Rcpp::Named("y_bound_residual") = fixed_v.ratio_residual,
      Rcpp::Named("x_lambda") = fixed_u.lambda,
      Rcpp::Named("y_lambda") = fixed_v.lambda,
      Rcpp::Named("x_maximum_tie_branch") = fixed_u.tie_branch,
      Rcpp::Named("y_maximum_tie_branch") = fixed_v.tie_branch,
      Rcpp::Named("singular_value") = d,
      Rcpp::Named("maximum_objective_decrease") = maximum_objective_decrease,
      Rcpp::Named("zero_block_step") = zero_step
    );
    if (!certified) {
      failure_component = static_cast<int>(component + 1);
      failure_message = zero_step ?
        "a PMD CCA block matrix-vector product is zero" :
        "the iteration limit or block-KKT certificate was not satisfied";
      break;
    }

    c6spc_anchor(u);
    // c6spc_anchor has fixed u's sign. If it changed the certified pair,
    // match v by checking d.
    const arma::vec anchored_mv = c6spc_cross_apply(
      x, y, left, right, deflation_values, v
    );
    if (arma::dot(u, anchored_mv) < 0.0) v *= -1.0;
    const double anchored_d = arma::dot(u, c6spc_cross_apply(
      x, y, left, right, deflation_values, v
    ));
    x_loadings.col(component) = u;
    y_loadings.col(component) = v;
    singular_values(component) = anchored_d;
    left.insert_cols(left.n_cols, u);
    right.insert_cols(right.n_cols, v);
    deflation_values.insert_rows(deflation_values.n_rows, 1);
    deflation_values(deflation_values.n_elem - 1) = anchored_d;
    ++completed;
  }

  return Rcpp::List::create(
    Rcpp::Named("converged") = completed == components,
    Rcpp::Named("completed") = static_cast<int>(completed),
    Rcpp::Named("failure_component") = failure_component,
    Rcpp::Named("failure_message") = failure_message,
    Rcpp::Named("x_loadings") = x_loadings,
    Rcpp::Named("y_loadings") = y_loadings,
    Rcpp::Named("singular_values") = singular_values,
    Rcpp::Named("certificates") = certificates
  );
}
