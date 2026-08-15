// Cai--Liu--Xia precision-adjusted two-sample maximum test.
//
// This file implements both the oracle statistic and the feasible statistic
// in equations (2), (6), and (7) of Cai, Liu, and Xia (2014).  In particular,
// the feasible denominator is computed from the transformed observations; it
// is not replaced by the diagonal of an estimated precision matrix.

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp17)]]

#include <RcppArmadillo.h>
#include <algorithm>
#include <cmath>
#include <limits>

namespace {

void clx_check_finite_matrix(const arma::mat& x, const char* name) {
  if (!x.is_finite()) {
    Rcpp::stop("`%s` must contain only finite values.", name);
  }
}

void clx_check_samples(const arma::mat& x, const arma::mat& y) {
  clx_check_finite_matrix(x, "x");
  clx_check_finite_matrix(y, "y");
  if (x.n_rows < 2 || y.n_rows < 2) {
    Rcpp::stop("CLX requires at least two observations in each group.");
  }
  if (x.n_cols < 2 || y.n_cols != x.n_cols) {
    Rcpp::stop(
      "`x` and `y` must have the same number of columns, with p >= 2."
    );
  }
}

long double clx_scaled_anchor_difference(double value,
                                         double anchor,
                                         long double scale) {
  const long double raw = static_cast<long double>(value) -
    static_cast<long double>(anchor);
  if (std::isfinite(raw)) {
    return raw / scale;
  }
  return static_cast<long double>(value) / scale -
    static_cast<long double>(anchor) / scale;
}

struct ClxCenteredSample {
  arma::mat centered;
  arma::vec mean;
};

ClxCenteredSample clx_stable_center(const arma::mat& x) {
  const arma::uword n = x.n_rows;
  const arma::uword p = x.n_cols;
  const long double n_ld = static_cast<long double>(n);
  arma::mat centered(n, p);
  arma::vec mean(p);

  for (arma::uword j = 0; j < p; ++j) {
    double coordinate_scale = 0.0;
    for (arma::uword i = 0; i < n; ++i) {
      coordinate_scale = std::max(coordinate_scale, std::abs(x(i, j)));
    }
    if (coordinate_scale == 0.0) {
      coordinate_scale = 1.0;
    }
    const long double scale = static_cast<long double>(coordinate_scale);
    const double anchor = x(0, j);
    long double sum_scaled_values = 0.0L;
    long double sum_scaled_deltas = 0.0L;
    for (arma::uword i = 0; i < n; ++i) {
      sum_scaled_values += static_cast<long double>(x(i, j)) / scale;
      sum_scaled_deltas += clx_scaled_anchor_difference(
        x(i, j), anchor, scale
      );
    }
    const long double mean_scaled_delta = sum_scaled_deltas / n_ld;
    const long double mean_value = scale * sum_scaled_values / n_ld;
    const double mean_double = static_cast<double>(mean_value);
    if (!std::isfinite(mean_double)) {
      Rcpp::stop(
        "A CLX group mean is not representable as a finite double."
      );
    }
    mean(j) = mean_double;
    for (arma::uword i = 0; i < n; ++i) {
      const long double centered_value = scale * (
        clx_scaled_anchor_difference(x(i, j), anchor, scale) -
        mean_scaled_delta
      );
      const double centered_double = static_cast<double>(centered_value);
      if (!std::isfinite(centered_double)) {
        Rcpp::stop(
          "A centered CLX observation is not representable as a finite "
          "double."
        );
      }
      centered(i, j) = centered_double;
    }
  }
  return {centered, mean};
}

arma::vec clx_stable_mean_difference(const arma::mat& x,
                                      const arma::mat& y) {
  const arma::uword p = x.n_cols;
  const long double n1 = static_cast<long double>(x.n_rows);
  const long double n2 = static_cast<long double>(y.n_rows);
  arma::vec difference(p);
  for (arma::uword j = 0; j < p; ++j) {
    double coordinate_scale = 0.0;
    for (arma::uword i = 0; i < x.n_rows; ++i) {
      coordinate_scale = std::max(coordinate_scale, std::abs(x(i, j)));
    }
    for (arma::uword i = 0; i < y.n_rows; ++i) {
      coordinate_scale = std::max(coordinate_scale, std::abs(y(i, j)));
    }
    if (coordinate_scale == 0.0) {
      coordinate_scale = 1.0;
    }
    const long double scale = static_cast<long double>(coordinate_scale);
    const double anchor = x(0, j);
    long double sum_x = 0.0L;
    long double sum_y = 0.0L;
    for (arma::uword i = 0; i < x.n_rows; ++i) {
      sum_x += clx_scaled_anchor_difference(x(i, j), anchor, scale);
    }
    for (arma::uword i = 0; i < y.n_rows; ++i) {
      sum_y += clx_scaled_anchor_difference(y(i, j), anchor, scale);
    }
    const long double value = scale * (sum_x / n1 - sum_y / n2);
    const double value_double = static_cast<double>(value);
    if (!std::isfinite(value_double)) {
      Rcpp::stop(
        "A CLX mean difference is not representable as a finite double."
      );
    }
    difference(j) = value_double;
  }
  return difference;
}

arma::mat clx_checked_symmetric(const arma::mat& matrix,
                                arma::uword dimension,
                                const char* name) {
  clx_check_finite_matrix(matrix, name);
  if (matrix.n_rows != dimension || matrix.n_cols != dimension) {
    Rcpp::stop("`%s` must be a square p by p matrix.", name);
  }

  const double matrix_scale = arma::abs(matrix).max();
  const double safe_scale = std::max(
    matrix_scale, std::numeric_limits<double>::min()
  );
  const double relative_asymmetry = arma::abs(
    matrix / safe_scale - matrix.t() / safe_scale
  ).max();
  const double relative_tolerance =
    std::sqrt(std::numeric_limits<double>::epsilon());
  if (!std::isfinite(relative_asymmetry) ||
      relative_asymmetry > relative_tolerance) {
    Rcpp::stop("`%s` must be symmetric within numerical tolerance.", name);
  }
  return 0.5 * matrix + 0.5 * matrix.t();
}

void clx_require_positive_vector(const arma::vec& values,
                                 const char* quantity) {
  if (!values.is_finite() || arma::any(values <= 0.0)) {
    Rcpp::stop(
      "CLX requires every %s to be finite and strictly positive.",
      quantity
    );
  }
}

}  // namespace


//' Adaptive-threshold precision estimator for the feasible CLX test
//'
//' The pooled covariance and the variance estimates used in its entrywise
//' thresholds both have denominator N = n1 + n2.  The threshold is
//' delta * sqrt(theta_ij * log(p) / N).  A relative eigenvalue floor is
//' applied only when needed to invert the thresholded covariance matrix.
//'
//' @param x,y Numeric observation-by-variable matrices.
//' @param delta Positive adaptive-threshold multiplier.
//' @param relative_eigen_floor Non-negative relative eigenvalue floor.
//' @return Internal list containing the precision estimate and diagnostics.
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List cpp_clx_adaptive_precision(const arma::mat& x,
                                      const arma::mat& y,
                                      double delta,
                                      double relative_eigen_floor) {
  clx_check_samples(x, y);
  if (!std::isfinite(delta) || delta <= 0.0) {
    Rcpp::stop("`threshold_delta` must be a finite positive number.");
  }
  if (!std::isfinite(relative_eigen_floor) ||
      relative_eigen_floor < 0.0 || relative_eigen_floor >= 1.0) {
    Rcpp::stop(
      "`eigen_floor` must be a finite number in [0, 1)."
    );
  }

  const arma::uword n1 = x.n_rows;
  const arma::uword n2 = y.n_rows;
  const arma::uword p = x.n_cols;
  const double total = static_cast<double>(n1 + n2);
  const double log_p = std::log(static_cast<double>(p));

  const ClxCenteredSample stable_x = clx_stable_center(x);
  const ClxCenteredSample stable_y = clx_stable_center(y);
  const double maximum_x = arma::abs(stable_x.centered).max();
  const double maximum_y = arma::abs(stable_y.centered).max();
  const double centered_scale = std::max(maximum_x, maximum_y) > 0.0 ?
    std::max(maximum_x, maximum_y) : 1.0;
  const arma::mat centered_x = stable_x.centered / centered_scale;
  const arma::mat centered_y = stable_y.centered / centered_scale;

  const arma::mat pooled_scaled =
    (centered_x.t() * centered_x + centered_y.t() * centered_y) / total;
  if (!pooled_scaled.is_finite()) {
    Rcpp::stop(
      "The internally scaled CLX pooled covariance is non-finite."
    );
  }

  arma::mat theta_scaled(p, p, arma::fill::zeros);
  arma::mat thresholds_scaled(p, p, arma::fill::zeros);
  arma::mat thresholded_scaled(p, p, arma::fill::zeros);
  double minimum_threshold = std::numeric_limits<double>::infinity();
  double maximum_threshold = 0.0;
  arma::uword retained_nonzero = 0;
  arma::uword retained_off_diagonal = 0;

  for (arma::uword i = 0; i < p; ++i) {
    for (arma::uword j = i; j < p; ++j) {
      const double sigma_ij = pooled_scaled(i, j);
      long double sum_squared = 0.0L;
      for (arma::uword k = 0; k < n1; ++k) {
        const long double product =
          static_cast<long double>(centered_x(k, i)) *
          static_cast<long double>(centered_x(k, j));
        const long double residual = product -
          static_cast<long double>(sigma_ij);
        sum_squared += residual * residual;
      }
      for (arma::uword k = 0; k < n2; ++k) {
        const long double product =
          static_cast<long double>(centered_y(k, i)) *
          static_cast<long double>(centered_y(k, j));
        const long double residual = product -
          static_cast<long double>(sigma_ij);
        sum_squared += residual * residual;
      }

      const long double theta_long = sum_squared /
        static_cast<long double>(total);
      const double theta_ij = static_cast<double>(theta_long);
      if (!std::isfinite(theta_ij) || theta_ij < 0.0) {
        Rcpp::stop(
          "The CLX adaptive-threshold variance estimate is non-finite."
        );
      }
      const double lambda_ij = delta * std::sqrt(theta_ij * log_p / total);
      if (!std::isfinite(lambda_ij)) {
        Rcpp::stop("The CLX adaptive threshold is non-finite.");
      }
      const double kept = std::abs(sigma_ij) >= lambda_ij ? sigma_ij : 0.0;

      theta_scaled(i, j) = theta_scaled(j, i) = theta_ij;
      thresholds_scaled(i, j) = thresholds_scaled(j, i) = lambda_ij;
      thresholded_scaled(i, j) = thresholded_scaled(j, i) = kept;
      minimum_threshold = std::min(minimum_threshold, lambda_ij);
      maximum_threshold = std::max(maximum_threshold, lambda_ij);
      if (kept != 0.0) {
        retained_nonzero += (i == j ? 1 : 2);
        if (i != j) {
          retained_off_diagonal += 2;
        }
      }
    }
  }

  arma::vec eigenvalues;
  arma::mat eigenvectors;
  if (!arma::eig_sym(eigenvalues, eigenvectors, thresholded_scaled) ||
      !eigenvalues.is_finite()) {
    Rcpp::stop(
      "The eigendecomposition of the thresholded CLX covariance failed."
    );
  }

  const double minimum_eigenvalue_before = eigenvalues.min();
  const double maximum_eigenvalue_before = eigenvalues.max();
  const double pooled_scale = arma::abs(pooled_scaled.diag()).max();
  const double spectral_scale = std::max(
    arma::abs(eigenvalues).max(),
    pooled_scale
  );
  if (!std::isfinite(spectral_scale) || spectral_scale <= 0.0) {
    Rcpp::stop(
      "The pooled CLX covariance has no positive numerical scale; "
      "a precision matrix cannot be estimated."
    );
  }

  const double absolute_eigen_floor =
    relative_eigen_floor * spectral_scale;
  arma::vec adjusted_eigenvalues = eigenvalues;
  arma::uword adjusted_count = 0;
  for (arma::uword i = 0; i < adjusted_eigenvalues.n_elem; ++i) {
    if (adjusted_eigenvalues(i) <= 0.0 && relative_eigen_floor == 0.0) {
      Rcpp::stop(
        "The adaptive-threshold covariance is not positive definite; "
        "set `eigen_floor` to a positive value to request a diagnosed "
        "eigenvalue adjustment."
      );
    }
    if (relative_eigen_floor > 0.0 &&
        adjusted_eigenvalues(i) < absolute_eigen_floor) {
      adjusted_eigenvalues(i) = absolute_eigen_floor;
      ++adjusted_count;
    }
  }
  if (!adjusted_eigenvalues.is_finite() ||
      arma::any(adjusted_eigenvalues <= 0.0)) {
    Rcpp::stop(
      "The requested CLX eigenvalue floor did not produce a finite, "
      "positive definite covariance estimate."
    );
  }

  arma::mat adjusted_covariance_scaled = eigenvectors *
    arma::diagmat(adjusted_eigenvalues) * eigenvectors.t();
  arma::mat precision_scaled = eigenvectors *
    arma::diagmat(1.0 / adjusted_eigenvalues) * eigenvectors.t();
  adjusted_covariance_scaled = 0.5 * adjusted_covariance_scaled +
    0.5 * adjusted_covariance_scaled.t();
  precision_scaled = 0.5 * precision_scaled + 0.5 * precision_scaled.t();
  if (!adjusted_covariance_scaled.is_finite() ||
      !precision_scaled.is_finite()) {
    Rcpp::stop(
      "The adaptive-threshold CLX precision estimate is non-finite."
    );
  }

  const long double centered_scale_ld =
    static_cast<long double>(centered_scale);
  const long double scale2_ld = centered_scale_ld * centered_scale_ld;
  const double scale2 = static_cast<double>(scale2_ld);
  const double inverse_scale2 = static_cast<double>(1.0L / scale2_ld);
  if (!std::isfinite(scale2) || scale2 <= 0.0 ||
      !std::isfinite(inverse_scale2) || inverse_scale2 <= 0.0) {
    Rcpp::stop(
      "The adaptive CLX data scale is outside the supported double-precision "
      "range after squaring."
    );
  }
  const long double scale4_ld = scale2_ld * scale2_ld;
  const arma::mat pooled = pooled_scaled * scale2;
  const arma::mat thresholds = thresholds_scaled * scale2;
  const arma::mat thresholded = thresholded_scaled * scale2;
  const arma::mat adjusted_covariance = adjusted_covariance_scaled * scale2;
  const arma::mat precision = precision_scaled * inverse_scale2;
  arma::mat theta = theta_scaled;
  theta.transform([scale4_ld](double value) {
    return static_cast<double>(static_cast<long double>(value) * scale4_ld);
  });
  const arma::vec eigenvalues_before = eigenvalues * scale2;
  const arma::vec eigenvalues_after = adjusted_eigenvalues * scale2;
  const double spectral_scale_raw = spectral_scale * scale2;
  const double absolute_eigen_floor_raw = absolute_eigen_floor * scale2;
  if (!pooled.is_finite() || !thresholds.is_finite() ||
      !thresholded.is_finite() || !adjusted_covariance.is_finite() ||
      !precision.is_finite()) {
    Rcpp::stop(
      "The adaptive CLX covariance or precision is not representable at the "
      "original data scale."
    );
  }

  return Rcpp::List::create(
    Rcpp::Named("precision") = precision,
    Rcpp::Named("pooled_covariance") = pooled,
    Rcpp::Named("theta") = theta,
    Rcpp::Named("thresholds") = thresholds,
    Rcpp::Named("thresholded_covariance") = thresholded,
    Rcpp::Named("adjusted_covariance") = adjusted_covariance,
    Rcpp::Named("eigenvalues_before") = eigenvalues_before,
    Rcpp::Named("eigenvalues_after") = eigenvalues_after,
    Rcpp::Named("minimum_eigenvalue_before") =
      minimum_eigenvalue_before * scale2,
    Rcpp::Named("maximum_eigenvalue_before") =
      maximum_eigenvalue_before * scale2,
    Rcpp::Named("minimum_eigenvalue_after") = eigenvalues_after.min(),
    Rcpp::Named("maximum_eigenvalue_after") = eigenvalues_after.max(),
    Rcpp::Named("spectral_scale") = spectral_scale_raw,
    Rcpp::Named("relative_eigen_floor") = relative_eigen_floor,
    Rcpp::Named("absolute_eigen_floor") = absolute_eigen_floor_raw,
    Rcpp::Named("adjustment_applied") = adjusted_count > 0,
    Rcpp::Named("adjusted_eigenvalues") =
      static_cast<double>(adjusted_count),
    Rcpp::Named("retained_nonzero_entries") =
      static_cast<double>(retained_nonzero),
    Rcpp::Named("retained_off_diagonal_entries") =
      static_cast<double>(retained_off_diagonal),
    Rcpp::Named("minimum_threshold") = minimum_threshold * scale2,
    Rcpp::Named("maximum_threshold") = maximum_threshold * scale2,
    Rcpp::Named("threshold_delta") = delta,
    Rcpp::Named("sample_size") = total,
    Rcpp::Named("internal_data_scale") = centered_scale,
    Rcpp::Named("pooled_covariance_scaled") = pooled_scaled,
    Rcpp::Named("theta_scaled") = theta_scaled,
    Rcpp::Named("thresholds_scaled") = thresholds_scaled,
    Rcpp::Named("thresholded_covariance_scaled") = thresholded_scaled,
    Rcpp::Named("adjusted_covariance_scaled") = adjusted_covariance_scaled,
    Rcpp::Named("precision_scaled") = precision_scaled
  );
}


//' Oracle or feasible Cai--Liu--Xia statistic kernel
//'
//' @param x,y Numeric observation-by-variable matrices.
//' @param precision Symmetric precision or estimated-precision matrix.
//' @param oracle If true, use diag(precision); otherwise use equation (7).
//' @return Internal list of statistic components.
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List cpp_clx_two_sample(const arma::mat& x,
                              const arma::mat& y,
                              const arma::mat& precision,
                              bool oracle) {
  clx_check_samples(x, y);
  const arma::uword n1 = x.n_rows;
  const arma::uword n2 = y.n_rows;
  const arma::uword p = x.n_cols;
  const arma::mat omega = clx_checked_symmetric(precision, p, "precision");

  if (oracle) {
    arma::vec eigenvalues;
    if (!arma::eig_sym(eigenvalues, omega) || !eigenvalues.is_finite() ||
        arma::any(eigenvalues <= 0.0)) {
      Rcpp::stop(
        "The oracle `precision` matrix must be positive definite."
      );
    }
  }

  const double n1_double = static_cast<double>(n1);
  const double n2_double = static_cast<double>(n2);
  const double total = n1_double + n2_double;
  const double effective_n = n1_double * n2_double / total;
  const ClxCenteredSample stable_x = clx_stable_center(x);
  const ClxCenteredSample stable_y = clx_stable_center(y);
  const arma::vec mean_x = stable_x.mean;
  const arma::vec mean_y = stable_y.mean;
  const arma::vec difference = clx_stable_mean_difference(x, y);
  const arma::vec transformed_difference = omega * difference;

  arma::vec denominator;
  arma::vec variance_x;
  arma::vec variance_y;
  if (oracle) {
    denominator = omega.diag();
  } else {
    const arma::mat centered_transformed_x = stable_x.centered * omega.t();
    const arma::mat centered_transformed_y = stable_y.centered * omega.t();
    variance_x = arma::sum(arma::square(centered_transformed_x), 0).t() /
      n1_double;
    variance_y = arma::sum(arma::square(centered_transformed_y), 0).t() /
      n2_double;
    denominator = (n1_double * variance_x + n2_double * variance_y) /
      total;
  }
  clx_require_positive_vector(
    denominator,
    oracle ? "oracle diagonal precision entries" :
      "transformed empirical pooled variances"
  );

  const arma::vec scores = std::sqrt(effective_n) *
    transformed_difference / arma::sqrt(denominator);
  const arma::vec squared_scores = arma::square(scores);
  if (!scores.is_finite() || !squared_scores.is_finite()) {
    Rcpp::stop(
      "CLX produced non-finite transformed scores at the supplied "
      "numerical scale."
    );
  }
  const arma::uword argmax = squared_scores.index_max();
  const double maximum = squared_scores(argmax);
  const double centered_maximum = maximum -
    2.0 * std::log(static_cast<double>(p)) +
    std::log(std::log(static_cast<double>(p)));
  if (!std::isfinite(maximum) || !std::isfinite(centered_maximum)) {
    Rcpp::stop("CLX produced a non-finite maximum statistic.");
  }

  return Rcpp::List::create(
    Rcpp::Named("M") = maximum,
    Rcpp::Named("G") = centered_maximum,
    Rcpp::Named("argmax") = static_cast<double>(argmax + 1),
    Rcpp::Named("scores") = scores,
    Rcpp::Named("squared_scores") = squared_scores,
    Rcpp::Named("denominator") = denominator,
    Rcpp::Named("variance_x") = variance_x,
    Rcpp::Named("variance_y") = variance_y,
    Rcpp::Named("mean_x") = mean_x,
    Rcpp::Named("mean_y") = mean_y,
    Rcpp::Named("difference") = difference,
    Rcpp::Named("transformed_difference") = transformed_difference,
    Rcpp::Named("effective_n") = effective_n,
    Rcpp::Named("n1") = n1_double,
    Rcpp::Named("n2") = n2_double,
    Rcpp::Named("N") = total,
    Rcpp::Named("p") = static_cast<double>(p),
    Rcpp::Named("oracle") = oracle
  );
}
