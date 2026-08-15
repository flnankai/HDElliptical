#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp17)]]

namespace {

bool cutoff_weight(const std::string& weight) {
  return weight == "winsor" || weight == "quadratic" ||
    weight == "ball" || weight == "shell" ||
    weight == "linear_redescending";
}

void validate_cutoffs(const arma::vec& cutoffs) {
  if (cutoffs.n_elem != 4 || !cutoffs.is_finite() ||
      arma::any(cutoffs < 0.0)) {
    Rcpp::stop(
      "`cutoffs` must contain four finite non-negative values."
    );
  }
  if (cutoffs(0) > cutoffs(1) || cutoffs(1) > cutoffs(2) ||
      cutoffs(2) > cutoffs(3)) {
    Rcpp::stop(
      "`cutoffs` must satisfy Q1 <= Q2 <= Q3 <= Q3_star."
    );
  }
}

double radial_multiplier(const std::string& weight,
                         const double radius,
                         const arma::vec& cutoffs) {
  if (weight == "identity") {
    return 1.0;
  }
  if (weight == "spatial") {
    return 1.0 / radius;
  }

  const double q1 = cutoffs(0);
  const double q2 = cutoffs(1);
  const double q3 = cutoffs(2);
  const double q3_star = cutoffs(3);
  if (weight == "winsor") {
    return radius <= q2 ? 1.0 : q2 / radius;
  }
  if (weight == "quadratic") {
    if (radius <= q2) {
      return 1.0;
    }
    const double ratio = q2 / radius;
    return ratio * ratio;
  }
  if (weight == "ball") {
    // The upper endpoint belongs to the ball.
    return radius <= q2 ? 1.0 : 0.0;
  }
  if (weight == "shell") {
    // Both endpoints belong to the shell, as in the primary definition.
    return radius >= q1 && radius <= q3 ? 1.0 : 0.0;
  }
  if (weight == "linear_redescending") {
    if (radius <= q2) {
      return 1.0;
    }
    if (radius > q3_star || q3_star == q2) {
      return 0.0;
    }
    return (q3_star - radius) / (q3_star - q2);
  }
  Rcpp::stop("Unknown radial-weight family.");
  return NA_REAL;
}

}  // namespace


// Radially transform centered observations and accumulate their outer products.
// The input is normally pre-scaled in R, so all built-in non-spatial
// multipliers lie in [0, 1]. Spatial signs are formed directly from the stable
// unit direction, avoiding multiplication by 1 / radius.
// [[Rcpp::export]]
Rcpp::List cpp_ch6rs_radial_transform(const arma::mat& centered,
                                      const std::string& weight,
                                      const arma::vec& cutoffs,
                                      const double zero_tol) {
  if (centered.n_rows == 0 || centered.n_cols == 0 ||
      !centered.is_finite()) {
    Rcpp::stop("`centered` must be a non-empty finite numeric matrix.");
  }
  if (!std::isfinite(zero_tol) || zero_tol < 0.0) {
    Rcpp::stop("`zero_tol` must be a finite non-negative number.");
  }
  if (weight != "identity" && weight != "spatial" &&
      !cutoff_weight(weight)) {
    Rcpp::stop("Unknown radial-weight family.");
  }
  if (cutoff_weight(weight)) {
    validate_cutoffs(cutoffs);
  }

  const arma::uword n = centered.n_rows;
  const arma::uword p = centered.n_cols;
  arma::mat transformed(n, p, arma::fill::zeros);
  arma::vec radii(n, arma::fill::zeros);
  arma::vec multipliers(n, arma::fill::zeros);
  arma::vec transformed_norms(n, arma::fill::zeros);
  arma::uword n_zero = 0;
  arma::uword n_nonzero_transform = 0;

  for (arma::uword i = 0; i < n; ++i) {
    const arma::rowvec row = centered.row(i);
    const double row_scale = arma::abs(row).max();
    if (row_scale == 0.0) {
      ++n_zero;
      if (weight != "spatial") {
        multipliers(i) = radial_multiplier(weight, 0.0, cutoffs);
      }
      continue;
    }
    const arma::rowvec scaled = row / row_scale;
    const double scaled_norm = arma::norm(scaled, 2);
    const double scaled_tolerance = zero_tol == 0.0 ?
      0.0 : zero_tol / row_scale;
    const double radius = row_scale * scaled_norm;
    radii(i) = radius;
    if (!(scaled_norm > scaled_tolerance)) {
      ++n_zero;
      if (weight != "spatial") {
        multipliers(i) = radial_multiplier(weight, radius, cutoffs);
      }
      continue;
    }

    if (weight == "spatial") {
      transformed.row(i) = scaled / scaled_norm;
      multipliers(i) = 1.0 / radius;
      transformed_norms(i) = 1.0;
      ++n_nonzero_transform;
      continue;
    }

    const double multiplier = radial_multiplier(weight, radius, cutoffs);
    if (!std::isfinite(multiplier) || multiplier < 0.0) {
      Rcpp::stop("The radial multiplier is non-finite or negative.");
    }
    multipliers(i) = multiplier;
    transformed.row(i) = multiplier * row;
    transformed_norms(i) = multiplier * radius;
    if (transformed_norms(i) > 0.0) {
      ++n_nonzero_transform;
    }
  }

  const arma::mat matrix_sum = transformed.t() * transformed;
  if (!matrix_sum.is_finite()) {
    Rcpp::stop("The weighted outer-product sum is non-finite.");
  }
  return Rcpp::List::create(
    Rcpp::Named("matrix_sum") = matrix_sum,
    Rcpp::Named("transformed") = transformed,
    Rcpp::Named("radii") = radii,
    Rcpp::Named("radial_multiplier") = multipliers,
    Rcpp::Named("transformed_norms") = transformed_norms,
    Rcpp::Named("n_zero") = static_cast<double>(n_zero),
    Rcpp::Named("n_nonzero_transform") =
      static_cast<double>(n_nonzero_transform)
  );
}
