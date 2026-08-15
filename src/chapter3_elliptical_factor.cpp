// Chapter 3: elliptical factor-model matrix estimation.
//
// The POET decomposition and certified CLIME/GLASSO solvers live in the
// corresponding R wrappers and existing Chapter 3 native modules.  This file
// implements the numerically delicate one-step Tyler map used by POET-TME.
// Each residual is normalized before its quadratic form is accumulated, so
// the kernel is invariant to rowwise positive rescaling and does not overflow
// when finite coordinates have opposite signs near DBL_MAX.

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp17)]]

#include <RcppArmadillo.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

namespace {

double ch3ef_checked_double(const long double value,
                            const char* quantity) {
  const long double maximum = static_cast<long double>(
    std::numeric_limits<double>::max()
  );
  if (!std::isfinite(value) || value > maximum || value < -maximum) {
    Rcpp::stop(
      "%s is outside the finite double range; no clipping or repair is "
      "applied.", quantity
    );
  }
  return static_cast<double>(value);
}

double ch3ef_diagnostic_double(const long double value) {
  const long double maximum = static_cast<long double>(
    std::numeric_limits<double>::max()
  );
  if (value > maximum) {
    return R_PosInf;
  }
  if (value < -maximum) {
    return R_NegInf;
  }
  return static_cast<double>(value);
}

void ch3ef_validate_finite(const arma::mat& value, const char* name) {
  if (!value.is_finite()) {
    Rcpp::stop("`%s` must contain only finite values.", name);
  }
}

} // namespace


// One-step self-normalized Tyler estimator, primary-paper equation (6).
// [[Rcpp::export]]
Rcpp::List cpp_ch3ef_tyler_one_step(const arma::mat& x,
                                    const arma::vec& center,
                                    const arma::mat& precision,
                                    const double zero_tolerance) {
  ch3ef_validate_finite(x, "x");
  ch3ef_validate_finite(precision, "precision");
  if (!center.is_finite()) {
    Rcpp::stop("`center` must contain only finite values.");
  }
  const arma::uword n = x.n_rows;
  const arma::uword p = x.n_cols;
  if (n < 1 || p < 1 || center.n_elem != p ||
      precision.n_rows != p || precision.n_cols != p) {
    Rcpp::stop("Tyler one-step inputs are not conformable.");
  }
  if (!std::isfinite(zero_tolerance) || zero_tolerance < 0.0) {
    Rcpp::stop("`zero_tolerance` must be finite and non-negative.");
  }

  const double symmetry_error = arma::abs(precision - precision.t()).max();
  const double precision_scale = std::max(1.0, arma::abs(precision).max());
  if (symmetry_error > 256.0 * std::numeric_limits<double>::epsilon() *
      precision_scale) {
    Rcpp::stop("`precision` must be numerically symmetric.");
  }
  const arma::mat symmetric_precision = 0.5 * (precision + precision.t());
  arma::mat lower;
  if (!arma::chol(lower, symmetric_precision, "lower")) {
    Rcpp::stop(
      "`precision` must be positive definite; no ridge or eigenvalue floor "
      "is applied."
    );
  }

  std::vector<long double> accumulator(p * p, 0.0L);
  std::vector<long double> standardized(p, 0.0L);
  long double minimum_quadratic =
    std::numeric_limits<long double>::infinity();
  long double maximum_quadratic = 0.0L;
  long double minimum_radius =
    std::numeric_limits<long double>::infinity();
  bool subtraction_fallback = false;

  for (arma::uword i = 0; i < n; ++i) {
    long double row_scale = 0.0L;
    for (arma::uword j = 0; j < p; ++j) {
      const long double difference = static_cast<long double>(x(i, j)) -
        static_cast<long double>(center(j));
      if (!std::isfinite(difference)) {
        Rcpp::stop("A centered residual is not finite in extended precision.");
      }
      row_scale = std::max(row_scale, std::abs(difference));
      if (!std::isfinite(static_cast<double>(x(i, j) - center(j)))) {
        subtraction_fallback = true;
      }
    }
    if (!(row_scale > 0.0L)) {
      Rcpp::stop(
        "The Tyler one-step map is undefined because a fitted residual is "
        "exactly zero."
      );
    }

    long double square_sum = 0.0L;
    for (arma::uword j = 0; j < p; ++j) {
      const long double difference = static_cast<long double>(x(i, j)) -
        static_cast<long double>(center(j));
      standardized[j] = difference / row_scale;
      square_sum += standardized[j] * standardized[j];
    }
    const long double radius = row_scale * std::sqrt(square_sum);
    minimum_radius = std::min(minimum_radius, radius);
    if (radius <= static_cast<long double>(zero_tolerance)) {
      Rcpp::stop(
        "The Tyler one-step map is undefined because a fitted residual is "
        "at or below `zero_tol`."
      );
    }

    long double quadratic = 0.0L;
    for (arma::uword j = 0; j < p; ++j) {
      for (arma::uword k = 0; k < p; ++k) {
        quadratic += standardized[j] *
          static_cast<long double>(symmetric_precision(j, k)) *
          standardized[k];
      }
    }
    if (!std::isfinite(quadratic) || !(quadratic > 0.0L)) {
      Rcpp::stop(
        "A Tyler quadratic form is non-positive or non-finite; no absolute "
        "value, ridge, or floor is applied."
      );
    }
    minimum_quadratic = std::min(minimum_quadratic, quadratic);
    maximum_quadratic = std::max(maximum_quadratic, quadratic);

    for (arma::uword j = 0; j < p; ++j) {
      for (arma::uword k = 0; k < p; ++k) {
        accumulator[j + p * k] +=
          standardized[j] * standardized[k] / quadratic;
      }
    }
  }

  arma::mat estimate(p, p, arma::fill::zeros);
  const long double multiplier = static_cast<long double>(p) /
    static_cast<long double>(n);
  for (arma::uword j = 0; j < p; ++j) {
    for (arma::uword k = 0; k < p; ++k) {
      estimate(j, k) = ch3ef_checked_double(
        multiplier * accumulator[j + p * k],
        "The Tyler one-step matrix"
      );
    }
  }
  estimate = 0.5 * (estimate + estimate.t());

  return Rcpp::List::create(
    Rcpp::Named("estimate") = estimate,
    Rcpp::Named("minimum_standardized_quadratic") =
      ch3ef_diagnostic_double(minimum_quadratic),
    Rcpp::Named("maximum_standardized_quadratic") =
      ch3ef_diagnostic_double(maximum_quadratic),
    Rcpp::Named("minimum_residual_radius") =
      ch3ef_diagnostic_double(minimum_radius),
    Rcpp::Named("subtraction_overflow_fallback") = subtraction_fallback,
    Rcpp::Named("trace") = arma::trace(estimate)
  );
}
