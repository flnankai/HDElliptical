// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp17)]]

#include <RcppArmadillo.h>
#include <cmath>
#include <limits>

namespace {

bool sph_unit_direction(const arma::rowvec& value,
                        arma::rowvec& direction,
                        double& norm,
                        double& log_norm) {
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

bool sph_difference_direction(const arma::rowvec& left,
                              const arma::rowvec& right,
                              arma::rowvec& direction) {
  arma::rowvec difference = left - right;
  if (!difference.is_finite()) {
    const double scale = std::max(
      arma::abs(left).max(), arma::abs(right).max()
    );
    if (scale == 0.0) {
      direction.zeros(left.n_elem);
      return false;
    }
    difference = left / scale - right / scale;
  }
  double norm = 0.0;
  double log_norm = R_NegInf;
  return sph_unit_direction(difference, direction, norm, log_norm);
}

void sph_require_finite_matrix(const arma::mat& x, const char* name) {
  if (!x.is_finite()) {
    Rcpp::stop("`%s` must contain only finite values.", name);
  }
}

} // namespace


// [[Rcpp::export]]
Rcpp::List cpp_elliptical_sphericity_sign_core(
    const arma::mat& residuals,
    const arma::rowvec& center_scaled,
    const bool compute_second_order,
    const bool keep_sscm) {
  sph_require_finite_matrix(residuals, "residuals");
  if (residuals.n_rows < 2 || residuals.n_cols < 2) {
    Rcpp::stop("`residuals` must have at least two rows and two columns.");
  }
  if (center_scaled.n_elem != residuals.n_cols ||
      !center_scaled.is_finite()) {
    Rcpp::stop(
      "`center_scaled` must be finite and have one value per column."
    );
  }

  const arma::uword n = residuals.n_rows;
  const arma::uword p = residuals.n_cols;
  arma::mat signs(n, p, arma::fill::zeros);
  arma::vec radii(n, arma::fill::zeros);
  arma::vec log_radii(n, arma::fill::value(R_NegInf));
  arma::uword n_zero = 0;

  for (arma::uword i = 0; i < n; ++i) {
    arma::rowvec direction(p, arma::fill::zeros);
    double norm = 0.0;
    double log_norm = R_NegInf;
    if (sph_unit_direction(
          residuals.row(i), direction, norm, log_norm
        )) {
      signs.row(i) = direction;
      radii(i) = norm;
      log_radii(i) = log_norm;
    } else {
      ++n_zero;
    }
  }

  long double ordered_sign_sum = 0.0L;
  for (arma::uword i = 0; i + 1 < n; ++i) {
    for (arma::uword j = i + 1; j < n; ++j) {
      const long double inner = static_cast<long double>(
        arma::dot(signs.row(i), signs.row(j))
      );
      ordered_sign_sum += 2.0L * inner * inner;
    }
  }
  const long double denominator = static_cast<long double>(n) *
    static_cast<long double>(n - 1);
  const long double q_tilde = static_cast<long double>(p) *
    ordered_sign_sum / denominator - 1.0L;
  if (!std::isfinite(static_cast<double>(q_tilde))) {
    Rcpp::stop("The spatial-sign sum statistic is non-finite.");
  }

  arma::rowvec sscm_diagonal = arma::mean(arma::square(signs), 0);
  const double n_double = static_cast<double>(n);
  const double p_double = static_cast<double>(p);
  const double diagonal_denominator = 2.0 * (1.0 - 1.0 / p_double);
  double max_diagonal_score = -std::numeric_limits<double>::infinity();
  arma::uword max_diagonal_index = 0;
  for (arma::uword j = 0; j < p; ++j) {
    const double deviation = sscm_diagonal(j) - 1.0 / p_double;
    const double score = n_double * p_double * (p_double + 2.0) *
      deviation * deviation / diagonal_denominator;
    if (score > max_diagonal_score) {
      max_diagonal_score = score;
      max_diagonal_index = j;
    }
  }

  double max_off_diagonal_score = -std::numeric_limits<double>::infinity();
  double max_off_diagonal_value = NA_REAL;
  arma::uword max_off_diagonal_row = 0;
  arma::uword max_off_diagonal_col = 1;
  for (arma::uword j = 0; j + 1 < p; ++j) {
    for (arma::uword k = j + 1; k < p; ++k) {
      const double value = arma::dot(signs.col(j), signs.col(k)) / n_double;
      const double score = n_double * p_double * (p_double + 2.0) *
        value * value;
      if (score > max_off_diagonal_score) {
        max_off_diagonal_score = score;
        max_off_diagonal_value = value;
        max_off_diagonal_row = j;
        max_off_diagonal_col = k;
      }
    }
  }

  const bool maximum_is_diagonal =
    max_diagonal_score >= max_off_diagonal_score;
  const double max_standardized_score = maximum_is_diagonal ?
    max_diagonal_score : max_off_diagonal_score;

  Rcpp::RObject corrected_radii = R_NilValue;
  double second_order_scale = NA_REAL;
  int corrected_nonpositive = 0;
  int corrected_nonfinite = 0;
  if (compute_second_order) {
    second_order_scale = std::max(1.0, arma::abs(center_scaled).max());
    const arma::rowvec scaled_center = center_scaled / second_order_scale;
    const double center_norm_squared = arma::dot(
      scaled_center, scaled_center
    );
    arma::vec corrected(n, arma::fill::value(NA_REAL));
    for (arma::uword i = 0; i < n; ++i) {
      const double radius = radii(i) / second_order_scale;
      if (!(radius > 0.0) || !std::isfinite(radius)) {
        ++corrected_nonfinite;
        continue;
      }
      const double value = radius +
        arma::dot(scaled_center, signs.row(i)) -
        0.5 * center_norm_squared / radius;
      corrected(i) = value;
      if (!std::isfinite(value)) {
        ++corrected_nonfinite;
      } else if (!(value > 0.0)) {
        ++corrected_nonpositive;
      }
    }
    corrected_radii = Rcpp::wrap(corrected);
  }

  Rcpp::RObject sscm = R_NilValue;
  if (keep_sscm) {
    sscm = Rcpp::wrap(signs.t() * signs / n_double);
  }

  return Rcpp::List::create(
    Rcpp::Named("signs") = signs,
    Rcpp::Named("radii_scaled") = radii,
    Rcpp::Named("log_radii_scaled") = log_radii,
    Rcpp::Named("n_zero") = static_cast<double>(n_zero),
    Rcpp::Named("ordered_sign_sum") =
      static_cast<double>(ordered_sign_sum),
    Rcpp::Named("q_tilde") = static_cast<double>(q_tilde),
    Rcpp::Named("sscm_diagonal") = sscm_diagonal,
    Rcpp::Named("sscm_trace") = arma::accu(sscm_diagonal),
    Rcpp::Named("max_diagonal_score") = max_diagonal_score,
    Rcpp::Named("max_diagonal_index") =
      static_cast<double>(max_diagonal_index + 1),
    Rcpp::Named("max_off_diagonal_score") = max_off_diagonal_score,
    Rcpp::Named("max_off_diagonal_value") = max_off_diagonal_value,
    Rcpp::Named("max_off_diagonal_indices") = Rcpp::IntegerVector::create(
      static_cast<int>(max_off_diagonal_row + 1),
      static_cast<int>(max_off_diagonal_col + 1)
    ),
    Rcpp::Named("max_standardized_score") = max_standardized_score,
    Rcpp::Named("maximum_is_diagonal") = maximum_is_diagonal,
    Rcpp::Named("corrected_radii_scaled") = corrected_radii,
    Rcpp::Named("second_order_scale") = second_order_scale,
    Rcpp::Named("corrected_nonpositive") = corrected_nonpositive,
    Rcpp::Named("corrected_nonfinite") = corrected_nonfinite,
    Rcpp::Named("sscm") = sscm
  );
}


// [[Rcpp::export]]
Rcpp::List cpp_elliptical_sphericity_rank_core(
    const arma::mat& x,
    const std::string& method,
    const bool keep_pair_signs) {
  sph_require_finite_matrix(x, "x");
  if (x.n_rows < 4 || x.n_cols < 2) {
    Rcpp::stop("`x` must have at least four rows and two columns.");
  }
  if (method != "spearman" && method != "kendall") {
    Rcpp::stop("`method` must be either 'spearman' or 'kendall'.");
  }

  const arma::uword n = x.n_rows;
  const arma::uword p = x.n_cols;
  const arma::uword pair_count = n * (n - 1) / 2;
  arma::mat pair_signs(pair_count, p, arma::fill::zeros);
  arma::umat pair_index(n, n, arma::fill::zeros);
  arma::imat endpoints(pair_count, 2, arma::fill::zeros);
  arma::uword cursor = 0;
  arma::uword zero_pairs = 0;

  for (arma::uword i = 0; i + 1 < n; ++i) {
    for (arma::uword j = i + 1; j < n; ++j) {
      arma::rowvec direction(p, arma::fill::zeros);
      if (!sph_difference_direction(x.row(i), x.row(j), direction)) {
        ++zero_pairs;
      }
      pair_signs.row(cursor) = direction;
      pair_index(i, j) = cursor;
      pair_index(j, i) = cursor;
      endpoints(cursor, 0) = static_cast<int>(i + 1);
      endpoints(cursor, 1) = static_cast<int>(j + 1);
      ++cursor;
    }
  }

  long double ordered_sum = 0.0L;
  for (arma::uword a = 0; a + 3 < n; ++a) {
    for (arma::uword b = a + 1; b + 2 < n; ++b) {
      for (arma::uword c = b + 1; c + 1 < n; ++c) {
        for (arma::uword d = c + 1; d < n; ++d) {
          const long double pairing_a = static_cast<long double>(arma::dot(
            pair_signs.row(pair_index(a, b)),
            pair_signs.row(pair_index(c, d))
          ));
          const long double pairing_b = static_cast<long double>(arma::dot(
            pair_signs.row(pair_index(a, c)),
            pair_signs.row(pair_index(b, d))
          ));
          const long double pairing_c = static_cast<long double>(arma::dot(
            pair_signs.row(pair_index(a, d)),
            pair_signs.row(pair_index(b, c))
          ));
          if (method == "spearman") {
            // The 24 ordered permutations of {a,b,c,d} reduce exactly to
            // 8 * (A*B + B*C - A*C) for the three disjoint pairings.
            ordered_sum += 8.0L * (
              pairing_a * pairing_b + pairing_b * pairing_c -
              pairing_a * pairing_c
            );
          } else {
            // Each of A^2, B^2 and C^2 occurs eight times in the ordered sum.
            ordered_sum += 8.0L * (
              pairing_a * pairing_a + pairing_b * pairing_b +
              pairing_c * pairing_c
            );
          }
        }
      }
    }
  }

  const long double ordered_denominator =
    static_cast<long double>(n) * static_cast<long double>(n - 1) *
    static_cast<long double>(n - 2) * static_cast<long double>(n - 3);
  const long double trace_estimate = method == "spearman" ?
    ordered_sum / (2.0L * ordered_denominator) :
    ordered_sum / ordered_denominator;
  const long double q_tilde = method == "spearman" ?
    4.0L * static_cast<long double>(p) * trace_estimate - 1.0L :
    static_cast<long double>(p) * trace_estimate - 1.0L;
  if (!std::isfinite(static_cast<double>(ordered_sum)) ||
      !std::isfinite(static_cast<double>(trace_estimate)) ||
      !std::isfinite(static_cast<double>(q_tilde))) {
    Rcpp::stop("The rank sphericity statistic is non-finite.");
  }

  Rcpp::RObject returned_pair_signs = R_NilValue;
  Rcpp::RObject returned_endpoints = R_NilValue;
  if (keep_pair_signs) {
    returned_pair_signs = Rcpp::wrap(pair_signs);
    returned_endpoints = Rcpp::wrap(endpoints);
  }

  return Rcpp::List::create(
    Rcpp::Named("ordered_sum") = static_cast<double>(ordered_sum),
    Rcpp::Named("ordered_denominator") =
      static_cast<double>(ordered_denominator),
    Rcpp::Named("trace_estimate") = static_cast<double>(trace_estimate),
    Rcpp::Named("q_tilde") = static_cast<double>(q_tilde),
    Rcpp::Named("n_pairs") = static_cast<double>(pair_count),
    Rcpp::Named("n_zero_pairs") = static_cast<double>(zero_pairs),
    Rcpp::Named("pair_signs") = returned_pair_signs,
    Rcpp::Named("pair_endpoints") = returned_endpoints,
    Rcpp::Named("quadruple_reduction") =
      "three disjoint pairings per unordered four-set"
  );
}
