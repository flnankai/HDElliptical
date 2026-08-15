// Elliptical regularized Hotelling testing of Feng, Zhou and Wang (2026).
//
// The kernel implements the fully feasible fixed-ridge statistic from the
// paper.  It works with residuals and the null displacement translated and
// scaled separately by the R wrapper.  The n by n companion representation
// avoids a persistent p by p inverse.  A row-space/null-space spectral
// decomposition evaluates the ridge quadratic without the catastrophic
// Woodbury subtraction that occurs when rho is small.

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp17)]]

#include <RcppArmadillo.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <string>
#include <vector>

namespace {

long double erht_sum(const arma::vec& x) {
  long double total = 0.0L;
  long double correction = 0.0L;
  for (arma::uword i = 0; i < x.n_elem; ++i) {
    const long double value = static_cast<long double>(x(i));
    const long double updated = total + value;
    if (std::abs(total) >= std::abs(value)) {
      correction += (total - updated) + value;
    } else {
      correction += (value - updated) + total;
    }
    total = updated;
  }
  return total + correction;
}

long double erht_dot(const arma::vec& x, const arma::vec& y) {
  long double total = 0.0L;
  long double correction = 0.0L;
  for (arma::uword i = 0; i < x.n_elem; ++i) {
    const long double value = static_cast<long double>(x(i)) *
      static_cast<long double>(y(i));
    const long double updated = total + value;
    if (std::abs(total) >= std::abs(value)) {
      correction += (total - updated) + value;
    } else {
      correction += (value - updated) + total;
    }
    total = updated;
  }
  return total + correction;
}

double erht_checked(const long double value, const char* quantity) {
  const long double maximum = static_cast<long double>(
    std::numeric_limits<double>::max()
  );
  if (!std::isfinite(value) || value > maximum || value < -maximum) {
    Rcpp::stop(
      "ERHT produced a non-finite %s on its internal scale; no ridge "
      "other than the requested rho, absolute-value repair, or numerical "
      "floor is applied.",
      quantity
    );
  }
  return static_cast<double>(value);
}

double erht_positive(const long double value, const char* quantity) {
  const double answer = erht_checked(value, quantity);
  if (!(answer > 0.0)) {
    Rcpp::stop(
      "ERHT requires a finite, strictly positive %s; no absolute-value "
      "repair or numerical floor is applied.",
      quantity
    );
  }
  return answer;
}

arma::rowvec erht_unit_row(const arma::rowvec& x,
                           double& radius) {
  const double row_scale = arma::abs(x).max();
  if (!(row_scale > 0.0) || !std::isfinite(row_scale)) {
    radius = 0.0;
    return arma::rowvec(x.n_elem, arma::fill::zeros);
  }
  const arma::rowvec scaled = x / row_scale;
  const double scaled_norm = arma::norm(scaled, 2);
  if (!(scaled_norm > 0.0) || !std::isfinite(scaled_norm)) {
    radius = 0.0;
    return arma::rowvec(x.n_elem, arma::fill::zeros);
  }
  const long double radius_ld = static_cast<long double>(row_scale) *
    static_cast<long double>(scaled_norm);
  radius = radius_ld > static_cast<long double>(
    std::numeric_limits<double>::max()
  ) ? std::numeric_limits<double>::infinity() :
    static_cast<double>(radius_ld);
  return scaled / scaled_norm;
}

arma::vec erht_power(const arma::vec& weights, const int exponent) {
  if (exponent == 0) return arma::vec(weights.n_elem, arma::fill::ones);
  if (exponent == 1) return weights;
  return arma::square(weights);
}

long double erht_weighted_off_diagonal(const arma::mat& squared_companion,
                                       const arma::vec& left,
                                       const arma::vec& right) {
  const arma::vec product = squared_companion * right;
  return erht_dot(left, product);
}

} // namespace


// [[Rcpp::export]]
Rcpp::List cpp_erht_grid(const arma::mat& residuals_scaled,
                         const arma::vec& null_difference_scaled,
                         const arma::vec& rho,
                         const bool keep_companion = false) {
  const arma::uword n = residuals_scaled.n_rows;
  const arma::uword p = residuals_scaled.n_cols;
  const arma::uword k_grid = rho.n_elem;
  if (n < 3U || p < 2U) {
    Rcpp::stop("ERHT requires at least three observations and two variables.");
  }
  if (null_difference_scaled.n_elem != p || k_grid < 1U ||
      !residuals_scaled.is_finite() ||
      !null_difference_scaled.is_finite() || !rho.is_finite()) {
    Rcpp::stop("Invalid dimensions or non-finite input supplied to ERHT.");
  }
  for (arma::uword k = 0; k < k_grid; ++k) {
    if (!(rho(k) > 0.0)) {
      Rcpp::stop("Every ERHT ridge value must be finite and positive.");
    }
  }

  arma::mat signs(n, p, arma::fill::zeros);
  arma::vec radii(n, arma::fill::zeros);
  arma::vec weights(n, arma::fill::zeros);
  const double root_p = std::sqrt(static_cast<double>(p));
  for (arma::uword i = 0; i < n; ++i) {
    double radius = 0.0;
    signs.row(i) = erht_unit_row(residuals_scaled.row(i), radius);
    if (!(radius > 0.0) || !std::isfinite(radius)) {
      Rcpp::stop(
        "ERHT inverse-distance calibration is undefined because a fitted "
        "spatial-median residual has zero or non-finite Euclidean norm."
      );
    }
    radii(i) = radius;
    weights(i) = root_p / radius;
    if (!std::isfinite(weights(i)) || !(weights(i) > 0.0)) {
      Rcpp::stop("ERHT produced a non-finite fitted inverse distance.");
    }
  }

  const double aspect = static_cast<double>(p) / static_cast<double>(n);
  arma::mat gram = aspect * (signs * signs.t());
  gram = 0.5 * (gram + gram.t());
  arma::vec gram_eigenvalues;
  if (!arma::eig_sym(gram_eigenvalues, gram)) {
    Rcpp::stop("The ERHT companion Gram eigendecomposition failed.");
  }
  const double raw_minimum_eigenvalue = gram_eigenvalues.min();
  const double spectral_scale = std::max(1.0, gram_eigenvalues.max());
  const double roundoff_tolerance =
    256.0 * std::numeric_limits<double>::epsilon() * spectral_scale;
  arma::uword clipped_eigenvalues = 0U;
  for (arma::uword i = 0; i < gram_eigenvalues.n_elem; ++i) {
    if (gram_eigenvalues(i) < -roundoff_tolerance) {
      Rcpp::stop(
        "The ERHT companion Gram matrix has a materially negative "
        "eigenvalue; the input or numerical calculation is invalid."
      );
    }
    if (gram_eigenvalues(i) < 0.0) {
      gram_eigenvalues(i) = 0.0;
      ++clipped_eigenvalues;
    }
  }

  arma::mat left_singular;
  arma::mat right_basis;
  arma::vec singular_values;
  if (!arma::svd_econ(
      left_singular, singular_values, right_basis, signs
  )) {
    Rcpp::stop("The ERHT sign-matrix economy SVD failed.");
  }
  const arma::uword spectral_terms = singular_values.n_elem;
  const arma::vec row_eigenvalues = aspect *
    arma::square(singular_values);
  const arma::vec row_coefficients = right_basis.t() *
    null_difference_scaled;
  const double singular_scale = spectral_terms > 0U ?
    singular_values.max() : 0.0;
  const double singular_tolerance =
    std::numeric_limits<double>::epsilon() *
    static_cast<double>(std::max(n, p)) * singular_scale;
  const arma::uword numerical_rank = arma::accu(
    singular_values > singular_tolerance
  );
  const arma::mat left_singular_squared = arma::square(left_singular);
  arma::vec left_null_diagonal(n, arma::fill::zeros);
  if (spectral_terms < n) {
    left_null_diagonal = 1.0 - arma::sum(left_singular_squared, 1);
    const double projector_tolerance =
      1024.0 * std::numeric_limits<double>::epsilon();
    for (arma::uword i = 0; i < n; ++i) {
      if (left_null_diagonal(i) < -projector_tolerance) {
        Rcpp::stop(
          "The ERHT left singular subspace produced a materially negative "
          "complement-projector diagonal."
        );
      }
      if (left_null_diagonal(i) < 0.0) left_null_diagonal(i) = 0.0;
    }
  }
  long double null_difference_squared = 0.0L;
  if (spectral_terms < p) {
    const arma::vec null_difference = null_difference_scaled -
      right_basis * row_coefficients;
    null_difference_squared = erht_dot(
      null_difference, null_difference
    );
  }
  const long double e_value = erht_sum(weights) /
    static_cast<long double>(n);
  const arma::vec weights_squared = arma::square(weights);
  const long double t_value = erht_sum(weights_squared) /
    static_cast<long double>(n);

  arma::vec raw_statistic(k_grid, arma::fill::zeros);
  arma::vec center(k_grid, arma::fill::zeros);
  arma::vec variance(k_grid, arma::fill::zeros);
  arma::vec z_statistic(k_grid, arma::fill::zeros);
  arma::vec kappa(k_grid, arma::fill::zeros);
  arma::vec b1(k_grid, arma::fill::zeros);
  arma::vec b2(k_grid, arma::fill::zeros);
  arma::vec denominator(k_grid, arma::fill::zeros);
  arma::vec mu_hat(k_grid, arma::fill::zeros);
  arma::vec sigma_d_squared(k_grid, arma::fill::zeros);
  arma::mat psi(k_grid, 6U, arma::fill::zeros);
  arma::mat g_values(k_grid, 3U, arma::fill::zeros);
  arma::mat companion_eigenvalues(k_grid, n, arma::fill::zeros);
  Rcpp::LogicalVector quadratic_roundoff_clipped(k_grid, false);
  Rcpp::List companion_matrices;
  if (keep_companion) companion_matrices = Rcpp::List(k_grid);

  const std::vector<std::pair<int, int>> psi_pairs = {
    {0, 0}, {0, 1}, {0, 2}, {1, 1}, {1, 2}, {2, 2}
  };

  for (arma::uword k = 0; k < k_grid; ++k) {
    const double ridge = rho(k);
    const arma::vec a_eigen = row_eigenvalues /
      (row_eigenvalues + ridge);
    // Form the complementary eigenvalues directly.  Computing 1 - a_eigen
    // loses every meaningful digit when rho is much smaller than a positive
    // row eigenvalue, precisely the regime in which the feasible centering
    // uses e-b1 and t-b2.  The B = I-A representation makes those two
    // differences direct weighted sums rather than subtractions of nearly
    // equal floating-point numbers.
    const arma::vec b_eigen = ridge /
      (row_eigenvalues + ridge);
    companion_eigenvalues.row(k).zeros();
    if (spectral_terms > 0U) {
      companion_eigenvalues.row(k).cols(0U, spectral_terms - 1U) =
        a_eigen.t();
    }

    const arma::vec b_diagonal = left_null_diagonal +
      left_singular_squared * b_eigen;
    const long double complement_mass =
      static_cast<long double>(n - spectral_terms) +
      erht_sum(b_eigen);
    const long double direct_mass = erht_sum(a_eigen);
    const bool use_complement_off_diagonal =
      complement_mass <= direct_mass;

    arma::mat off_diagonal_source;
    arma::mat companion;
    arma::vec a_diagonal;
    if (use_complement_off_diagonal) {
      arma::mat scaled_vectors = left_singular;
      scaled_vectors.each_row() %= b_eigen.t();
      arma::mat complement = scaled_vectors * left_singular.t();
      if (spectral_terms < n) {
        complement += arma::eye<arma::mat>(n, n) -
          left_singular * left_singular.t();
      }
      complement = 0.5 * (complement + complement.t());
      off_diagonal_source = complement;
      if (keep_companion) {
        companion = -complement;
        companion.diag() += 1.0;
        companion = 0.5 * (companion + companion.t());
      }
      a_diagonal = 1.0 - b_diagonal;
    } else {
      arma::mat scaled_vectors = left_singular;
      scaled_vectors.each_row() %= a_eigen.t();
      companion = scaled_vectors * left_singular.t();
      companion = 0.5 * (companion + companion.t());
      off_diagonal_source = companion;
      a_diagonal = companion.diag();
    }
    if (keep_companion) companion_matrices[k] = companion;

    long double q_form = null_difference_squared /
      static_cast<long double>(ridge);
    for (arma::uword j = 0; j < spectral_terms; ++j) {
      const long double coefficient = static_cast<long double>(
        row_coefficients(j)
      );
      q_form += coefficient * coefficient /
        static_cast<long double>(row_eigenvalues(j) + ridge);
    }
    raw_statistic(k) = erht_checked(
      static_cast<long double>(n) * q_form,
      "fixed-ridge statistic"
    );

    const long double kappa_value = direct_mass /
      static_cast<long double>(n);
    const long double e_minus_b1 = erht_dot(b_diagonal, weights) /
      static_cast<long double>(n);
    const long double t_minus_b2 = erht_dot(
      b_diagonal, weights_squared
    ) /
      static_cast<long double>(n);
    const long double b1_value = use_complement_off_diagonal ?
      e_value - e_minus_b1 :
      erht_dot(a_diagonal, weights) / static_cast<long double>(n);
    const long double b2_value = use_complement_off_diagonal ?
      t_value - t_minus_b2 :
      erht_dot(a_diagonal, weights_squared) /
        static_cast<long double>(n);
    const long double d_value = e_minus_b1 * e_minus_b1 +
      kappa_value * t_minus_b2;
    const double d_checked = erht_positive(d_value, "conditional denominator");
    const long double mu_value = kappa_value / d_value;

    // A_ij^2 = B_ij^2 off the diagonal.  Use whichever of A or B has the
    // smaller spectral mass so that forming those entries does not subtract
    // an almost-identity matrix from I.
    arma::mat squared_companion = arma::square(off_diagonal_source);
    squared_companion.diag().zeros();
    for (arma::uword q = 0; q < psi_pairs.size(); ++q) {
      const arma::vec left = erht_power(weights, psi_pairs[q].first);
      const arma::vec right = erht_power(weights, psi_pairs[q].second);
      psi(k, q) = erht_checked(
        erht_weighted_off_diagonal(squared_companion, left, right) /
          static_cast<long double>(n),
        "weighted off-diagonal companion functional"
      );
    }
    const long double psi00 = psi(k, 0);
    const long double psi01 = psi(k, 1);
    const long double psi02 = psi(k, 2);
    const long double psi11 = psi(k, 3);
    const long double psi12 = psi(k, 4);
    const long double psi22 = psi(k, 5);
    arma::mat gamma(3U, 3U, arma::fill::zeros);
    gamma(0, 0) = 2.0 * psi00;
    gamma(0, 1) = gamma(1, 0) = 2.0 * psi01;
    gamma(0, 2) = gamma(2, 0) = 2.0 * psi11;
    gamma(1, 1) = psi02 + psi11;
    gamma(1, 2) = gamma(2, 1) = 2.0 * psi12;
    gamma(2, 2) = 2.0 * psi22;

    const long double inverse_d_squared = 1.0L / (d_value * d_value);
    arma::vec g(3U);
    g(0) = erht_checked(
      e_minus_b1 * e_minus_b1 * inverse_d_squared,
      "first variance gradient component"
    );
    g(1) = erht_checked(
      2.0L * kappa_value * e_minus_b1 * inverse_d_squared,
      "second variance gradient component"
    );
    g(2) = erht_checked(
      kappa_value * kappa_value * inverse_d_squared,
      "third variance gradient component"
    );
    g_values.row(k) = g.t();
    const long double sigma_value = static_cast<long double>(
      arma::as_scalar(g.t() * gamma * g)
    );
    const double sigma_checked = erht_positive(
      sigma_value, "conditional variance functional"
    );
    const long double center_value = static_cast<long double>(n) * mu_value;
    const long double variance_value = static_cast<long double>(n) *
      sigma_value;
    const long double z_value =
      (static_cast<long double>(raw_statistic(k)) - center_value) /
      std::sqrt(variance_value);

    kappa(k) = erht_checked(kappa_value, "companion trace functional");
    b1(k) = erht_checked(b1_value, "first diagonal-weight functional");
    b2(k) = erht_checked(b2_value, "second diagonal-weight functional");
    denominator(k) = d_checked;
    mu_hat(k) = erht_checked(mu_value, "conditional center functional");
    sigma_d_squared(k) = sigma_checked;
    center(k) = erht_checked(center_value, "fixed-ridge center");
    variance(k) = erht_positive(variance_value, "fixed-ridge variance");
    z_statistic(k) = erht_checked(z_value, "standardized statistic");

    quadratic_roundoff_clipped[k] = false;
  }

  Rcpp::List answer = Rcpp::List::create(
    Rcpp::Named("raw_statistic_scaled") = raw_statistic,
    Rcpp::Named("center_scaled") = center,
    Rcpp::Named("variance_scaled") = variance,
    Rcpp::Named("z") = z_statistic,
    Rcpp::Named("rho") = rho,
    Rcpp::Named("kappa") = kappa,
    Rcpp::Named("e_scaled") = erht_checked(e_value, "mean inverse distance"),
    Rcpp::Named("t_scaled") = erht_checked(t_value, "mean squared inverse distance"),
    Rcpp::Named("b1_scaled") = b1,
    Rcpp::Named("b2_scaled") = b2,
    Rcpp::Named("D_scaled") = denominator,
    Rcpp::Named("mu_scaled") = mu_hat,
    Rcpp::Named("sigma_D_squared_scaled") = sigma_d_squared,
    Rcpp::Named("psi_scaled") = psi,
    Rcpp::Named("g_scaled") = g_values,
    Rcpp::Named("radii_scaled") = radii,
    Rcpp::Named("weights_scaled") = weights,
    Rcpp::Named("gram_eigenvalues") = gram_eigenvalues,
    Rcpp::Named("companion_eigenvalues") = companion_eigenvalues,
    Rcpp::Named("raw_minimum_gram_eigenvalue") = raw_minimum_eigenvalue,
    Rcpp::Named("roundoff_clipped_gram_eigenvalues") =
      static_cast<double>(clipped_eigenvalues),
    Rcpp::Named("quadratic_roundoff_clipped") =
      quadratic_roundoff_clipped,
    Rcpp::Named("quadratic_spectral_rank") =
      static_cast<double>(numerical_rank),
    Rcpp::Named("quadratic_spectral_terms") =
      static_cast<double>(spectral_terms),
    Rcpp::Named("quadratic_singular_tolerance") =
      singular_tolerance,
    Rcpp::Named("quadratic_route") =
      "stable economy-SVD row-space/null-space spectral decomposition",
    Rcpp::Named("aspect_ratio") = aspect,
    Rcpp::Named("n") = static_cast<int>(n),
    Rcpp::Named("p") = static_cast<int>(p),
    Rcpp::Named("keep_companion") = keep_companion
  );
  if (keep_companion) answer["companion_matrices"] = companion_matrices;
  return answer;
}
