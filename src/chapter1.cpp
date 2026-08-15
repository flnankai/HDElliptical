#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp17)]]

namespace {

void require_finite(const arma::mat& x, const char* name) {
  if (!x.is_finite()) {
    Rcpp::stop("`%s` must contain only finite values.", name);
  }
}

void require_finite(const arma::vec& x, const char* name) {
  if (!x.is_finite()) {
    Rcpp::stop("`%s` must contain only finite values.", name);
  }
}

arma::vec row_norms(const arma::mat& x) {
  arma::vec norms(x.n_rows, arma::fill::zeros);
  for (arma::uword i = 0; i < x.n_rows; ++i) {
    const arma::rowvec row = x.row(i);
    const double row_scale = arma::abs(row).max();
    if (row_scale > 0.0) {
      norms(i) = row_scale * arma::norm(row / row_scale, 2);
    }
  }
  return norms;
}

bool difference_direction(const arma::rowvec& left,
                          const arma::rowvec& right,
                          const double zero_tol,
                          arma::rowvec& direction) {
  // Subtract first whenever that operation is finite so nearby large values
  // retain all representable low-order bits. Only opposite values whose
  // direct difference overflows use a common pre-scaling fallback.
  arma::rowvec difference = left - right;
  double coordinate_scale = 1.0;
  if (!difference.is_finite()) {
    coordinate_scale = std::max(
      arma::abs(left).max(), arma::abs(right).max()
    );
    if (coordinate_scale == 0.0) {
      direction.zeros(left.n_elem);
      return false;
    }
    difference = left / coordinate_scale - right / coordinate_scale;
  }
  const double difference_scale = arma::abs(difference).max();
  if (difference_scale == 0.0) {
    direction.zeros(left.n_elem);
    return false;
  }
  const arma::rowvec scaled = difference / difference_scale;
  const double scaled_norm = arma::norm(scaled, 2);
  double scaled_tolerance = zero_tol / coordinate_scale;
  if (scaled_tolerance != 0.0) {
    scaled_tolerance /= difference_scale;
  }
  if (!(scaled_norm > scaled_tolerance)) {
    direction.zeros(left.n_elem);
    return false;
  }
  direction = scaled / scaled_norm;
  return true;
}

arma::mat sign_rows(const arma::mat& x, const double zero_tol,
                    arma::vec* norms_out = nullptr,
                    arma::uword* zeros_out = nullptr) {
  arma::vec norms(x.n_rows, arma::fill::zeros);
  arma::mat signs(x.n_rows, x.n_cols, arma::fill::zeros);
  arma::uword zeros = 0;

  for (arma::uword i = 0; i < x.n_rows; ++i) {
    const arma::rowvec row = x.row(i);
    const double row_scale = arma::abs(row).max();
    if (row_scale == 0.0) {
      ++zeros;
      continue;
    }
    const arma::rowvec scaled = row / row_scale;
    const double scaled_norm = arma::norm(scaled, 2);
    norms(i) = row_scale * scaled_norm;
    const double scaled_tolerance = zero_tol == 0.0 ?
      0.0 : zero_tol / row_scale;
    if (scaled_norm > scaled_tolerance) {
      signs.row(i) = scaled / scaled_norm;
    } else {
      ++zeros;
    }
  }

  if (norms_out != nullptr) {
    *norms_out = norms;
  }
  if (zeros_out != nullptr) {
    *zeros_out = zeros;
  }
  return signs;
}

double observation_subgradient_residual(const arma::mat& x,
                                        const arma::vec& weights,
                                        const arma::uword candidate,
                                        const double zero_tol,
                                        const double weight_sum) {
  const arma::mat residuals = x.each_row() - x.row(candidate);
  const arma::mat directions = sign_rows(
    residuals, zero_tol, nullptr, nullptr
  );
  arma::rowvec score(x.n_cols, arma::fill::zeros);
  double coincident_weight = 0.0;
  for (arma::uword i = 0; i < x.n_rows; ++i) {
    if (arma::dot(directions.row(i), directions.row(i)) > 0.0) {
      score += weights(i) * directions.row(i);
    } else {
      coincident_weight += weights(i);
    }
  }
  const double score_norm = arma::norm(score, 2);
  return std::isfinite(score_norm) ?
    std::max(0.0, score_norm - coincident_weight) / weight_sum :
    R_PosInf;
}

double weighted_distance_objective(const arma::mat& x,
                                   const arma::rowvec& location,
                                   const arma::vec& weights) {
  const arma::mat residuals = x.each_row() - location;
  const arma::vec distances = row_norms(residuals);
  const double maximum = distances.max();
  if (maximum == 0.0) {
    return 0.0;
  }
  const arma::vec normalized_weights = weights / arma::accu(weights);
  return maximum * arma::dot(normalized_weights, distances / maximum);
}

bool spd_roots(const arma::mat& matrix, arma::mat& square_root,
               arma::mat& inverse_square_root, const double eigen_tol) {
  arma::vec eigenvalues;
  arma::mat eigenvectors;
  const arma::mat symmetric = 0.5 * (matrix + matrix.t());
  if (!arma::eig_sym(eigenvalues, eigenvectors, symmetric)) {
    return false;
  }
  const double threshold = std::max(eigen_tol, 0.0);
  if (!(eigenvalues.max() > 0.0) || eigenvalues.min() <= threshold) {
    return false;
  }
  square_root = eigenvectors * arma::diagmat(arma::sqrt(eigenvalues)) *
    eigenvectors.t();
  inverse_square_root = eigenvectors *
    arma::diagmat(1.0 / arma::sqrt(eigenvalues)) * eigenvectors.t();
  return true;
}

bool spd_inverse(const arma::mat& matrix, arma::mat& inverse,
                 const double eigen_tol) {
  arma::vec eigenvalues;
  arma::mat eigenvectors;
  const arma::mat symmetric = 0.5 * (matrix + matrix.t());
  if (!arma::eig_sym(eigenvalues, eigenvectors, symmetric)) {
    return false;
  }
  const double threshold = std::max(eigen_tol, 0.0);
  if (!(eigenvalues.max() > 0.0) || eigenvalues.min() <= threshold) {
    return false;
  }
  const arma::vec inverse_values = 1.0 / eigenvalues;
  if (!inverse_values.is_finite()) {
    return false;
  }
  inverse = eigenvectors * arma::diagmat(inverse_values) * eigenvectors.t();
  return inverse.is_finite();
}

arma::mat normalize_trace(arma::mat matrix, const double target) {
  matrix = 0.5 * (matrix + matrix.t());
  const double matrix_trace = arma::trace(matrix);
  if (!std::isfinite(matrix_trace) || matrix_trace <= 0.0) {
    Rcpp::stop("The shape update has a non-positive or non-finite trace.");
  }
  return target * matrix / matrix_trace;
}

arma::mat tyler_update(const arma::mat& centered, const arma::mat& shape,
                       const double zero_tol) {
  arma::mat inverse_shape;
  if (!spd_inverse(shape, inverse_shape, 0.0)) {
    Rcpp::stop("The current shape iterate is not positive definite.");
  }

  const arma::uword n = centered.n_rows;
  const arma::uword p = centered.n_cols;
  arma::mat update(p, p, arma::fill::zeros);
  for (arma::uword i = 0; i < n; ++i) {
    const arma::rowvec observation = centered.row(i);
    const double denominator = arma::as_scalar(
      observation * inverse_shape * observation.t());
    if (!std::isfinite(denominator) || denominator <= zero_tol * zero_tol) {
      Rcpp::stop(
        "Tyler's equation is undefined for centered observation %d; "
        "choose a different center or remove the exact center point.",
        static_cast<int>(i + 1)
      );
    }
    update += observation.t() * observation / denominator;
  }
  update *= static_cast<double>(p) / static_cast<double>(n);
  return normalize_trace(update, static_cast<double>(p));
}

}  // namespace


// [[Rcpp::export]]
Rcpp::List cpp_spatial_sign(const arma::mat& x, const double zero_tol) {
  require_finite(x, "x");
  if (!std::isfinite(zero_tol) || zero_tol < 0.0) {
    Rcpp::stop("`zero_tol` must be a finite non-negative number.");
  }
  arma::vec norms;
  arma::uword zeros = 0;
  arma::mat signs = sign_rows(x, zero_tol, &norms, &zeros);
  return Rcpp::List::create(
    Rcpp::Named("signs") = signs,
    Rcpp::Named("norms") = norms,
    Rcpp::Named("n_zero") = static_cast<double>(zeros)
  );
}


// [[Rcpp::export]]
Rcpp::List cpp_spatial_median(const arma::mat& x, const arma::vec& weights,
                              const arma::rowvec& initial,
                              const double tol, const int max_iter,
                              const double zero_tol) {
  require_finite(x, "x");
  require_finite(weights, "weights");
  if (x.n_rows == 0 || x.n_cols == 0) {
    Rcpp::stop("`x` must have at least one row and one column.");
  }
  if (weights.n_elem != x.n_rows) {
    Rcpp::stop("`weights` must have one value per row of `x`.");
  }
  if (initial.n_elem != x.n_cols || !initial.is_finite()) {
    Rcpp::stop("`initial` must be a finite vector with one value per column.");
  }
  if (arma::any(weights < 0.0) || arma::accu(weights) <= 0.0) {
    Rcpp::stop("`weights` must be non-negative and have a positive sum.");
  }
  if (!std::isfinite(tol) || tol <= 0.0 || max_iter < 1 ||
      !std::isfinite(zero_tol) || zero_tol < 0.0) {
    Rcpp::stop("Invalid convergence controls.");
  }

  const double weight_sum = arma::accu(weights);
  if (!std::isfinite(weight_sum)) {
    Rcpp::stop("The sum of `weights` must be finite.");
  }
  arma::rowvec location = initial;
  const arma::mat initial_residuals = x.each_row() - initial;
  const arma::vec initial_distances = row_norms(initial_residuals);
  const double max_initial_distance = initial_distances.max();
  double convergence_scale = 0.0;
  if (max_initial_distance > 0.0) {
    const arma::vec scaled_distances =
      initial_distances / max_initial_distance;
    convergence_scale = max_initial_distance * std::sqrt(
      arma::dot(weights, arma::square(scaled_distances)) / weight_sum
    );
  }
  if (!std::isfinite(convergence_scale) || convergence_scale <= zero_tol) {
    convergence_scale = 1.0;
  }

  bool converged = false;
  int iterations = 0;
  double change = R_PosInf;

  for (int iteration = 1; iteration <= max_iter; ++iteration) {
    iterations = iteration;
    const arma::mat residuals = x.each_row() - location;
    arma::vec distances;
    const arma::mat directions = sign_rows(
      residuals, zero_tol, &distances, nullptr
    );
    arma::rowvec residual_sum(x.n_cols, arma::fill::zeros);
    double coincident_weight = 0.0;
    double minimum_distance = R_PosInf;
    arma::uword nearest_noncoincident = 0;

    for (arma::uword i = 0; i < x.n_rows; ++i) {
      if (arma::dot(directions.row(i), directions.row(i)) > 0.0) {
        residual_sum += weights(i) * directions.row(i);
        if (distances(i) < minimum_distance) {
          minimum_distance = distances(i);
          nearest_noncoincident = i;
        }
      } else {
        coincident_weight += weights(i);
      }
    }

    const double residual_norm = arma::norm(residual_sum, 2);
    const double current_equation_residual = std::isfinite(residual_norm) ?
      std::max(0.0, residual_norm - coincident_weight) / weight_sum :
      R_PosInf;
    if (current_equation_residual <= tol) {
      converged = true;
      change = 0.0;
      break;
    }
    if (!std::isfinite(minimum_distance) || minimum_distance <= 0.0) {
      break;
    }

    // A spatial median may be an observation. Once an observation is the
    // closest noncoincident point to the iterate, test its exact subgradient
    // condition and snap to it if valid. This avoids asymptotic stalling on
    // subnormal separations without an unconditional O(n^2 p) prescan.
    const double nearest_residual = observation_subgradient_residual(
      x, weights, nearest_noncoincident, zero_tol, weight_sum
    );
    if (nearest_residual <= tol) {
      const arma::rowvec next = x.row(nearest_noncoincident);
      change = row_norms(arma::mat(next - location))(0) /
        convergence_scale;
      location = next;
      converged = true;
      break;
    }

    double scaled_denominator = 0.0;
    for (arma::uword i = 0; i < x.n_rows; ++i) {
      if (arma::dot(directions.row(i), directions.row(i)) > 0.0) {
        scaled_denominator += weights(i) *
          (minimum_distance / distances(i));
      }
    }
    if (!std::isfinite(scaled_denominator) || scaled_denominator <= 0.0) {
      break;
    }

    const arma::rowvec weiszfeld = location +
      residual_sum * (minimum_distance / scaled_denominator);
    double move_fraction = 1.0;
    if (coincident_weight > 0.0 && residual_norm > 0.0) {
      move_fraction = std::max(0.0, 1.0 - coincident_weight / residual_norm);
    }
    const arma::rowvec next =
      move_fraction * weiszfeld + (1.0 - move_fraction) * location;
    change = arma::norm(next - location, 2) / convergence_scale;
    location = next;

  }

  const arma::mat final_residuals = x.each_row() - location;
  const arma::mat final_directions = sign_rows(
    final_residuals, zero_tol, nullptr, nullptr
  );
  arma::rowvec final_sum(x.n_cols, arma::fill::zeros);
  double final_coincident_weight = 0.0;
  for (arma::uword i = 0; i < x.n_rows; ++i) {
    if (arma::dot(final_directions.row(i), final_directions.row(i)) > 0.0) {
      final_sum += weights(i) * final_directions.row(i);
    } else {
      final_coincident_weight += weights(i);
    }
  }
  const double final_norm = arma::norm(final_sum, 2);
  const double subgradient_residual = std::isfinite(final_norm) ?
    std::max(0.0, final_norm - final_coincident_weight) / weight_sum :
    R_PosInf;
  converged = std::isfinite(subgradient_residual) &&
    subgradient_residual <= tol;

  return Rcpp::List::create(
    Rcpp::Named("location") = location.t(),
    Rcpp::Named("objective") =
      weighted_distance_objective(x, location, weights),
    Rcpp::Named("iterations") = iterations,
    Rcpp::Named("converged") = converged,
    Rcpp::Named("relative_change") = change,
    Rcpp::Named("equation_residual") = subgradient_residual
  );
}


// [[Rcpp::export]]
arma::mat cpp_spatial_rank(const arma::mat& x, const arma::mat& reference,
                           const double zero_tol) {
  require_finite(x, "x");
  require_finite(reference, "reference");
  if (x.n_cols != reference.n_cols) {
    Rcpp::stop("`x` and `reference` must have the same number of columns.");
  }
  if (reference.n_rows == 0) {
    Rcpp::stop("`reference` must have at least one row.");
  }
  if (!std::isfinite(zero_tol) || zero_tol < 0.0) {
    Rcpp::stop("`zero_tol` must be a finite non-negative number.");
  }

  arma::mat ranks(x.n_rows, x.n_cols, arma::fill::zeros);
  for (arma::uword i = 0; i < x.n_rows; ++i) {
    for (arma::uword j = 0; j < reference.n_rows; ++j) {
      arma::rowvec direction(x.n_cols, arma::fill::zeros);
      if (difference_direction(
            x.row(i), reference.row(j), zero_tol, direction
          )) {
        ranks.row(i) += direction;
      }
    }
  }
  ranks /= static_cast<double>(reference.n_rows);
  return ranks;
}


// [[Rcpp::export]]
Rcpp::List cpp_spatial_kendall(const arma::mat& x, const double zero_tol) {
  require_finite(x, "x");
  if (x.n_rows < 2 || x.n_cols == 0) {
    Rcpp::stop("`x` must contain at least two observations and one variable.");
  }
  if (!std::isfinite(zero_tol) || zero_tol < 0.0) {
    Rcpp::stop("`zero_tol` must be a finite non-negative number.");
  }

  const double pair_count = static_cast<double>(x.n_rows) *
    static_cast<double>(x.n_rows - 1) / 2.0;
  double zero_pairs = 0.0;
  arma::mat estimate(x.n_cols, x.n_cols, arma::fill::zeros);
  for (arma::uword i = 0; i + 1 < x.n_rows; ++i) {
    for (arma::uword j = i + 1; j < x.n_rows; ++j) {
      arma::rowvec direction(x.n_cols, arma::fill::zeros);
      if (difference_direction(x.row(i), x.row(j), zero_tol, direction)) {
        estimate += direction.t() * direction;
      } else {
        zero_pairs += 1.0;
      }
    }
  }
  estimate /= pair_count;
  return Rcpp::List::create(
    Rcpp::Named("matrix") = estimate,
    Rcpp::Named("n_pairs") = pair_count,
    Rcpp::Named("n_zero_pairs") = zero_pairs
  );
}


// [[Rcpp::export]]
Rcpp::List cpp_tyler_shape(const arma::mat& centered,
                           const arma::mat& initial,
                           const double tol, const int max_iter,
                           const double zero_tol) {
  require_finite(centered, "centered");
  require_finite(initial, "initial");
  if (centered.n_rows == 0 || centered.n_cols == 0) {
    Rcpp::stop("`centered` must have at least one row and one column.");
  }
  if (initial.n_rows != centered.n_cols ||
      initial.n_cols != centered.n_cols) {
    Rcpp::stop("`initial` must be a square matrix matching the data dimension.");
  }
  if (!std::isfinite(tol) || tol <= 0.0 || max_iter < 1 ||
      !std::isfinite(zero_tol) || zero_tol < 0.0) {
    Rcpp::stop("Invalid convergence controls.");
  }

  const double p = static_cast<double>(centered.n_cols);
  arma::vec radii;
  arma::uword zeros = 0;
  const arma::mat directions = sign_rows(centered, zero_tol, &radii, &zeros);
  if (zeros > 0) {
    Rcpp::stop(
      "Tyler's equation is undefined because %d centered observation(s) "
      "are zero.",
      static_cast<int>(zeros)
    );
  }
  arma::mat shape = normalize_trace(initial, p);
  arma::mat inverse_shape;
  if (!spd_inverse(shape, inverse_shape, 0.0)) {
    Rcpp::stop("`initial` must be symmetric positive definite.");
  }

  bool converged = false;
  int iterations = 0;
  double relative_change = R_PosInf;
  for (int iteration = 1; iteration <= max_iter; ++iteration) {
    iterations = iteration;
    const arma::mat next = tyler_update(directions, shape, 0.0);
    relative_change = arma::norm(next - shape, "fro") /
      std::max(1.0, arma::norm(shape, "fro"));
    shape = next;
    if (relative_change <= tol) {
      converged = true;
      break;
    }
  }

  const arma::mat mapped = tyler_update(directions, shape, 0.0);
  const double equation_residual = arma::norm(mapped - shape, "fro") /
    std::max(1.0, arma::norm(shape, "fro"));
  converged = equation_residual <= tol;
  return Rcpp::List::create(
    Rcpp::Named("shape") = shape,
    Rcpp::Named("iterations") = iterations,
    Rcpp::Named("converged") = converged,
    Rcpp::Named("relative_change") = relative_change,
    Rcpp::Named("equation_residual") = equation_residual
  );
}


// [[Rcpp::export]]
double cpp_acg_loglik(const arma::mat& x, const arma::mat& shape) {
  require_finite(x, "x");
  require_finite(shape, "shape");
  if (x.n_rows == 0 || x.n_cols == 0) {
    Rcpp::stop("`x` must have at least one row and one column.");
  }
  if (shape.n_rows != x.n_cols || shape.n_cols != x.n_cols) {
    Rcpp::stop("`shape` must be a square matrix matching the data dimension.");
  }

  arma::uword zeros = 0;
  const arma::mat directions = sign_rows(x, 0.0, nullptr, &zeros);
  if (zeros > 0) {
    Rcpp::stop("Every row of `x` must be nonzero.");
  }

  arma::vec eigenvalues;
  arma::mat eigenvectors;
  const arma::mat symmetric = 0.5 * (shape + shape.t());
  if (!arma::eig_sym(eigenvalues, eigenvectors, symmetric) ||
      eigenvalues.min() <= 0.0) {
    Rcpp::stop("`shape` must be symmetric positive definite.");
  }
  const arma::mat projections = directions * eigenvectors;
  const arma::vec log_eigenvalues = arma::log(eigenvalues);
  double sum_log_quadratic = 0.0;
  for (arma::uword i = 0; i < projections.n_rows; ++i) {
    arma::vec log_terms(projections.n_cols, arma::fill::value(-arma::datum::inf));
    for (arma::uword j = 0; j < projections.n_cols; ++j) {
      const double magnitude = std::abs(projections(i, j));
      if (magnitude > 0.0) {
        log_terms(j) = 2.0 * std::log(magnitude) - log_eigenvalues(j);
      }
    }
    const double maximum = log_terms.max();
    if (!std::isfinite(maximum)) {
      Rcpp::stop("The ACG quadratic form is numerically undefined.");
    }
    const double log_quadratic = maximum +
      std::log(arma::accu(arma::exp(log_terms - maximum)));
    sum_log_quadratic += log_quadratic;
  }

  return -0.5 * static_cast<double>(x.n_rows) *
      arma::accu(log_eigenvalues) -
    0.5 * static_cast<double>(x.n_cols) * sum_log_quadratic;
}


// [[Rcpp::export]]
Rcpp::List cpp_hr_estimator(const arma::mat& x,
                            const arma::rowvec& initial_location,
                            const arma::mat& initial_shape,
                            const double tol, const int max_iter,
                            const double zero_tol) {
  require_finite(x, "x");
  require_finite(initial_shape, "initial_shape");
  if (x.n_rows == 0 || x.n_cols == 0) {
    Rcpp::stop("`x` must have at least one row and one column.");
  }
  if (initial_location.n_elem != x.n_cols ||
      !initial_location.is_finite()) {
    Rcpp::stop("`initial_location` must match the data dimension.");
  }
  if (initial_shape.n_rows != x.n_cols ||
      initial_shape.n_cols != x.n_cols) {
    Rcpp::stop("`initial_shape` must be a square matrix matching the data dimension.");
  }
  if (!std::isfinite(tol) || tol <= 0.0 || max_iter < 1 ||
      !std::isfinite(zero_tol) || zero_tol < 0.0) {
    Rcpp::stop("Invalid convergence controls.");
  }

  const double p = static_cast<double>(x.n_cols);
  arma::rowvec location = initial_location;
  arma::mat shape = normalize_trace(initial_shape, p);
  const arma::mat initial_residuals = x.each_row() - initial_location;
  const arma::vec initial_radii = row_norms(initial_residuals);
  const double maximum_initial_radius = initial_radii.max();
  double location_scale = 0.0;
  if (maximum_initial_radius > 0.0) {
    const arma::vec scaled_radii = initial_radii / maximum_initial_radius;
    location_scale = maximum_initial_radius * std::sqrt(
      arma::dot(scaled_radii, scaled_radii) /
      static_cast<double>(x.n_rows)
    );
  }
  if (!std::isfinite(location_scale) || location_scale <= zero_tol) {
    location_scale = 1.0;
  }
  bool converged = false;
  int iterations = 0;
  double location_change = R_PosInf;
  double shape_change = R_PosInf;
  arma::rowvec best_location = location;
  arma::mat best_shape = shape;
  double best_score = R_PosInf;
  double best_location_change = R_NaReal;
  double best_shape_change = R_NaReal;
  int best_iteration = 0;

  for (int iteration = 1; iteration <= max_iter; ++iteration) {
    iterations = iteration;
    arma::mat square_root;
    arma::mat inverse_square_root;
    if (!spd_roots(shape, square_root, inverse_square_root, 0.0)) {
      Rcpp::stop("The current HR shape iterate is not positive definite.");
    }

    const arma::mat residuals = x.each_row() - location;
    const arma::mat whitened = residuals * inverse_square_root;
    arma::vec norms;
    arma::uword zeros = 0;
    const arma::mat signs = sign_rows(whitened, zero_tol, &norms, &zeros);
    if (zeros > 0) {
      Rcpp::stop(
        "The HR equations are undefined because %d residual(s) are zero.",
        static_cast<int>(zeros)
      );
    }

    const arma::mat current_shape_equation = p * signs.t() * signs /
      static_cast<double>(x.n_rows);
    const double current_location_residual = arma::norm(
      arma::mean(signs, 0), 2
    );
    const double current_shape_residual = arma::norm(
      current_shape_equation -
        arma::eye<arma::mat>(x.n_cols, x.n_cols), "fro"
    );
    const double current_score = std::max(
      current_location_residual, current_shape_residual
    );
    if (current_score < best_score) {
      best_score = current_score;
      best_location = location;
      best_shape = shape;
      best_location_change = location_change;
      best_shape_change = shape_change;
      best_iteration = iteration - 1;
    }
    if (current_score <= tol) {
      converged = true;
      location_change = 0.0;
      shape_change = 0.0;
      break;
    }

    const double minimum_radius = norms.min();
    const double scaled_inverse_radius_sum = arma::accu(minimum_radius / norms);
    if (!std::isfinite(scaled_inverse_radius_sum) ||
        scaled_inverse_radius_sum <= 0.0) {
      Rcpp::stop("The HR location update has an invalid radial denominator.");
    }
    const arma::rowvec step = arma::sum(signs, 0) * square_root *
      (minimum_radius / scaled_inverse_radius_sum);
    const arma::rowvec next_location = location + step;

    arma::mat next_shape = p * square_root *
      (signs.t() * signs / static_cast<double>(x.n_rows)) * square_root;
    next_shape = normalize_trace(next_shape, p);

    location_change = arma::norm(next_location - location, 2) / location_scale;
    shape_change = arma::norm(next_shape - shape, "fro") /
      std::max(1.0, arma::norm(shape, "fro"));
    location = next_location;
    shape = next_shape;
  }

  arma::vec final_norms;
  arma::uword final_zeros = 0;
  arma::mat final_signs;
  double location_residual = R_PosInf;
  double shape_residual = R_PosInf;
  auto evaluate_final = [&]() -> bool {
    arma::mat square_root;
    arma::mat inverse_square_root;
    if (!spd_roots(shape, square_root, inverse_square_root, 0.0)) {
      return false;
    }
    const arma::mat final_whitened = (x.each_row() - location) *
      inverse_square_root;
    final_signs = sign_rows(
      final_whitened, zero_tol, &final_norms, &final_zeros
    );
    if (final_zeros > 0) {
      return false;
    }
    location_residual = arma::norm(arma::mean(final_signs, 0), 2);
    const arma::mat shape_equation = p * final_signs.t() * final_signs /
      static_cast<double>(x.n_rows);
    shape_residual = arma::norm(
      shape_equation - arma::eye<arma::mat>(x.n_cols, x.n_cols), "fro"
    );
    return std::isfinite(location_residual) && std::isfinite(shape_residual);
  };

  bool final_defined = evaluate_final();
  double final_score = final_defined ?
    std::max(location_residual, shape_residual) : R_PosInf;
  if (final_defined && final_score < best_score) {
    best_score = final_score;
    best_location = location;
    best_shape = shape;
    best_location_change = location_change;
    best_shape_change = shape_change;
    best_iteration = iterations;
  } else if (!converged && best_score < final_score) {
    location = best_location;
    shape = best_shape;
    location_change = best_location_change;
    shape_change = best_shape_change;
    final_defined = evaluate_final();
    final_score = std::max(location_residual, shape_residual);
  }
  if (!final_defined) {
    Rcpp::stop(
      "The HR equations are undefined at every retained iterate; the final "
      "candidate has %d zero residual(s).",
      static_cast<int>(final_zeros)
    );
  }
  converged = final_score <= tol;

  return Rcpp::List::create(
    Rcpp::Named("location") = location.t(),
    Rcpp::Named("shape") = shape,
    Rcpp::Named("iterations") = iterations,
    Rcpp::Named("best_iteration") = best_iteration,
    Rcpp::Named("converged") = converged,
    Rcpp::Named("location_change") = location_change,
    Rcpp::Named("shape_change") = shape_change,
    Rcpp::Named("location_equation_residual") = location_residual,
    Rcpp::Named("shape_equation_residual") = shape_residual,
    Rcpp::Named("equation_residual") = final_score,
    Rcpp::Named("n_zero_residuals") = static_cast<double>(final_zeros)
  );
}
