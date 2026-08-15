// Xu--Lin--Wei--Pan (2016) analytical adaptive sum-of-powers test.
//
// The first kernel constructs inverse-variance-standardised mean contrasts
// and their estimated correlation matrix.  The second evaluates observed
// finite-power sums and their Gaussian null moments.  Bivariate Gaussian
// moments are computed from Isserlis pairings as a polynomial in rho, avoiding
// simulation and avoiding square-root formulas that are unstable near
// |rho| = 1.

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp17)]]

#include <RcppArmadillo.h>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <limits>
#include <string>
#include <vector>

namespace {

struct AspuScaledData {
  arma::mat x;
  arma::mat y;
  std::vector<long double> origin;
  std::vector<long double> scale;
};

struct AspuAccumulator {
  long double value = 0.0L;
  long double correction = 0.0L;

  void add(long double increment) {
    const long double adjusted = increment - correction;
    const long double updated = value + adjusted;
    correction = (updated - value) - adjusted;
    value = updated;
  }
};

void aspu_require_finite_matrix(const arma::mat& x, const char* name) {
  if (!x.is_finite()) {
    Rcpp::stop("`%s` must contain only finite values.", name);
  }
}

double aspu_checked_double(long double value, const char* quantity) {
  const double result = static_cast<double>(value);
  if (!std::isfinite(result)) {
    Rcpp::stop(
      "Xu-Lin-Wei-Pan aSPU produced a non-finite or nonrepresentable %s.",
      quantity
    );
  }
  return result;
}

double aspu_checked_positive_double(long double value,
                                    const char* quantity) {
  const double result = static_cast<double>(value);
  if (!std::isfinite(result) || result <= 0.0) {
    Rcpp::stop(
      "Xu-Lin-Wei-Pan aSPU requires a finite, strictly positive and "
      "double-representable %s; no variance floor is applied.",
      quantity
    );
  }
  return result;
}

double aspu_report_long_double(long double value) {
  const long double maximum = static_cast<long double>(
    std::numeric_limits<double>::max()
  );
  if (value > maximum) {
    return std::numeric_limits<double>::infinity();
  }
  if (value < -maximum) {
    return -std::numeric_limits<double>::infinity();
  }
  return static_cast<double>(value);
}

AspuScaledData aspu_scale_columns(const arma::mat& x,
                                  const arma::mat& y) {
  AspuScaledData out;
  const arma::uword n1 = x.n_rows;
  const arma::uword n2 = y.n_rows;
  const arma::uword p = x.n_cols;
  out.x.set_size(n1, p);
  out.y.set_size(n2, p);
  out.origin.resize(p);
  out.scale.resize(p);

  for (arma::uword j = 0; j < p; ++j) {
    const long double anchor = static_cast<long double>(x(0, j));
    long double maximum_difference = 0.0L;
    bool direct_finite = true;
    for (arma::uword i = 0; i < n1; ++i) {
      const long double difference =
        static_cast<long double>(x(i, j)) - anchor;
      direct_finite = direct_finite && std::isfinite(difference);
      maximum_difference = std::max(
        maximum_difference, std::abs(difference)
      );
    }
    for (arma::uword i = 0; i < n2; ++i) {
      const long double difference =
        static_cast<long double>(y(i, j)) - anchor;
      direct_finite = direct_finite && std::isfinite(difference);
      maximum_difference = std::max(
        maximum_difference, std::abs(difference)
      );
    }

    out.origin[j] = anchor;
    if (direct_finite && std::isfinite(maximum_difference)) {
      const long double scale = maximum_difference > 0.0L ?
        maximum_difference : 1.0L;
      out.scale[j] = scale;
      for (arma::uword i = 0; i < n1; ++i) {
        out.x(i, j) = static_cast<double>(
          (static_cast<long double>(x(i, j)) - anchor) / scale
        );
      }
      for (arma::uword i = 0; i < n2; ++i) {
        out.y(i, j) = static_cast<double>(
          (static_cast<long double>(y(i, j)) - anchor) / scale
        );
      }
      continue;
    }

    // Fallback for platforms where subtraction of opposite finite double
    // extremes can overflow even in long double.
    double operand_scale = std::abs(x(0, j));
    for (arma::uword i = 0; i < n1; ++i) {
      operand_scale = std::max(operand_scale, std::abs(x(i, j)));
    }
    for (arma::uword i = 0; i < n2; ++i) {
      operand_scale = std::max(operand_scale, std::abs(y(i, j)));
    }
    if (operand_scale == 0.0) {
      operand_scale = 1.0;
    }
    const long double operand_scale_ld =
      static_cast<long double>(operand_scale);
    const long double anchor_scaled = anchor / operand_scale_ld;
    maximum_difference = 0.0L;
    for (arma::uword i = 0; i < n1; ++i) {
      const long double difference =
        static_cast<long double>(x(i, j)) / operand_scale_ld -
        anchor_scaled;
      out.x(i, j) = static_cast<double>(difference);
      maximum_difference = std::max(
        maximum_difference, std::abs(difference)
      );
    }
    for (arma::uword i = 0; i < n2; ++i) {
      const long double difference =
        static_cast<long double>(y(i, j)) / operand_scale_ld -
        anchor_scaled;
      out.y(i, j) = static_cast<double>(difference);
      maximum_difference = std::max(
        maximum_difference, std::abs(difference)
      );
    }
    const long double second_scale = maximum_difference > 0.0L ?
      maximum_difference : 1.0L;
    for (arma::uword i = 0; i < n1; ++i) {
      out.x(i, j) = static_cast<double>(
        static_cast<long double>(out.x(i, j)) / second_scale
      );
    }
    for (arma::uword i = 0; i < n2; ++i) {
      out.y(i, j) = static_cast<double>(
        static_cast<long double>(out.y(i, j)) / second_scale
      );
    }
    out.scale[j] = operand_scale_ld * second_scale;
  }
  return out;
}

arma::mat aspu_sample_covariance(const arma::mat& values,
                                 const arma::rowvec& mean) {
  arma::mat centered = values;
  centered.each_row() -= mean;
  const double denominator = static_cast<double>(values.n_rows - 1);
  const arma::mat covariance = centered.t() * centered / denominator;
  if (!covariance.is_finite()) {
    Rcpp::stop("Xu-Lin-Wei-Pan aSPU produced a non-finite covariance.");
  }
  return covariance;
}

arma::mat aspu_band_covariance(const arma::mat& covariance,
                               int bandwidth) {
  if (bandwidth < 0) {
    return covariance;
  }
  arma::mat result = covariance;
  for (arma::uword j = 0; j < result.n_rows; ++j) {
    for (arma::uword k = 0; k < result.n_cols; ++k) {
      const long long distance = std::llabs(
        static_cast<long long>(j) - static_cast<long long>(k)
      );
      if (distance > static_cast<long long>(bandwidth)) {
        result(j, k) = 0.0;
      }
    }
  }
  return result;
}

arma::mat aspu_covariance_to_correlation(const arma::mat& covariance,
                                         const char* source) {
  const arma::uword p = covariance.n_rows;
  arma::mat correlation(p, p, arma::fill::eye);
  for (arma::uword j = 0; j < p; ++j) {
    if (!std::isfinite(covariance(j, j)) || covariance(j, j) <= 0.0) {
      Rcpp::stop(
        "Xu-Lin-Wei-Pan aSPU requires every %s marginal variance to be "
        "finite and strictly positive; no variance floor is applied.",
        source
      );
    }
  }
  for (arma::uword j = 0; j < p; ++j) {
    for (arma::uword k = 0; k < j; ++k) {
      const double denominator =
        std::sqrt(covariance(j, j)) * std::sqrt(covariance(k, k));
      const double value = covariance(j, k) / denominator;
      if (!std::isfinite(value)) {
        Rcpp::stop(
          "Xu-Lin-Wei-Pan aSPU produced a non-finite %s correlation.",
          source
        );
      }
      correlation(j, k) = value;
      correlation(k, j) = value;
    }
  }
  return correlation;
}

long double aspu_integer_power(long double base, int exponent) {
  long double result = 1.0L;
  long double factor = base;
  int remaining = exponent;
  while (remaining > 0) {
    if ((remaining & 1) != 0) {
      result *= factor;
    }
    remaining >>= 1;
    if (remaining > 0) {
      factor *= factor;
    }
  }
  return result;
}

long double aspu_double_factorial(int value) {
  if (value <= 0) {
    return 1.0L;
  }
  long double result = 1.0L;
  for (int current = value; current > 1; current -= 2) {
    result *= static_cast<long double>(current);
  }
  return result;
}

long double aspu_factorial(int value) {
  long double result = 1.0L;
  for (int current = 2; current <= value; ++current) {
    result *= static_cast<long double>(current);
  }
  return result;
}

long double aspu_choose(int n, int k) {
  if (k < 0 || k > n) {
    return 0.0L;
  }
  k = std::min(k, n - k);
  long double result = 1.0L;
  for (int j = 1; j <= k; ++j) {
    result *= static_cast<long double>(n - k + j) /
      static_cast<long double>(j);
  }
  return result;
}

long double aspu_univariate_normal_moment(int power,
                                          long double variance) {
  if ((power & 1) != 0) {
    return 0.0L;
  }
  return aspu_double_factorial(power - 1) *
    aspu_integer_power(variance, power / 2);
}

long double aspu_bivariate_normal_moment(int first_power,
                                         int second_power,
                                         long double first_variance,
                                         long double second_variance,
                                         long double covariance) {
  if (((first_power + second_power) & 1) != 0) {
    return 0.0L;
  }
  AspuAccumulator result;
  const int maximum_cross = std::min(first_power, second_power);
  for (int cross = 0; cross <= maximum_cross; ++cross) {
    if (((first_power - cross) & 1) != 0 ||
        ((second_power - cross) & 1) != 0) {
      continue;
    }
    const long double coefficient =
      aspu_choose(first_power, cross) *
      aspu_choose(second_power, cross) *
      aspu_factorial(cross) *
      aspu_double_factorial(first_power - cross - 1) *
      aspu_double_factorial(second_power - cross - 1);
    result.add(
      coefficient *
      aspu_integer_power(covariance, cross) *
      aspu_integer_power(
        first_variance, (first_power - cross) / 2
      ) *
      aspu_integer_power(
        second_variance, (second_power - cross) / 2
      )
    );
  }
  return result.value;
}

}  // namespace


//' Construct standardised contrasts for the analytical aSPU test
//'
//' @param x,y Numeric observation-by-variable matrices.
//' @param correlation_source One of `"common"`, `"unequal"`, or
//'   `"supplied"`.
//' @param bandwidth1,bandwidth2 Hard-band widths; `-1` means no banding.
//' @param supplied_correlation Optional supplied correlation matrix.
//' @param supplied_standard_errors Optional supplied marginal standard errors.
//' @return Internal list of contrasts and covariance diagnostics.
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List cpp_aspu_standardize(
    const arma::mat& x,
    const arma::mat& y,
    std::string correlation_source,
    int bandwidth1,
    int bandwidth2,
    const arma::mat& supplied_correlation,
    const arma::vec& supplied_standard_errors) {
  aspu_require_finite_matrix(x, "x");
  aspu_require_finite_matrix(y, "y");
  const arma::uword n1 = x.n_rows;
  const arma::uword n2 = y.n_rows;
  const arma::uword p = x.n_cols;
  if (n1 < 2 || n2 < 2) {
    Rcpp::stop(
      "Xu-Lin-Wei-Pan aSPU requires at least two observations in each "
      "group for sample variances."
    );
  }
  if (p < 1 || y.n_cols != p) {
    Rcpp::stop("`x` and `y` must have the same positive number of columns.");
  }
  if (correlation_source != "common" &&
      correlation_source != "unequal" &&
      correlation_source != "supplied") {
    Rcpp::stop(
      "`correlation_source` must be common, unequal, or supplied."
    );
  }
  if (bandwidth1 < -1 || bandwidth2 < -1 ||
      bandwidth1 >= static_cast<int>(p) ||
      bandwidth2 >= static_cast<int>(p)) {
    Rcpp::stop("Every band width must be -1 or an integer in 0, ..., p-1.");
  }

  const AspuScaledData scaled = aspu_scale_columns(x, y);
  const arma::rowvec mean1_scaled = arma::mean(scaled.x, 0);
  const arma::rowvec mean2_scaled = arma::mean(scaled.y, 0);
  const arma::rowvec difference_scaled = mean1_scaled - mean2_scaled;
  const arma::mat covariance1 = aspu_sample_covariance(
    scaled.x, mean1_scaled
  );
  const arma::mat covariance2 = aspu_sample_covariance(
    scaled.y, mean2_scaled
  );
  const arma::mat pooled_covariance = (
    static_cast<double>(n1 - 1) * covariance1 +
    static_cast<double>(n2 - 1) * covariance2
  ) / static_cast<double>(n1 + n2 - 2);

  arma::vec standard_error_scaled(p);
  arma::mat correlation;
  std::string standard_error_source;
  if (correlation_source == "common") {
    const double sample_factor = 1.0 / static_cast<double>(n1) +
      1.0 / static_cast<double>(n2);
    for (arma::uword j = 0; j < p; ++j) {
      standard_error_scaled(j) = std::sqrt(
        pooled_covariance(j, j) * sample_factor
      );
    }
    const arma::mat banded = aspu_band_covariance(
      pooled_covariance, bandwidth1
    );
    correlation = aspu_covariance_to_correlation(
      banded, "common pooled"
    );
    standard_error_source = "common pooled sample variance";
  } else if (correlation_source == "unequal") {
    const arma::mat mean_covariance =
      covariance1 / static_cast<double>(n1) +
      covariance2 / static_cast<double>(n2);
    for (arma::uword j = 0; j < p; ++j) {
      standard_error_scaled(j) = std::sqrt(mean_covariance(j, j));
    }
    const arma::mat banded_mean_covariance =
      aspu_band_covariance(covariance1, bandwidth1) /
        static_cast<double>(n1) +
      aspu_band_covariance(covariance2, bandwidth2) /
        static_cast<double>(n2);
    correlation = aspu_covariance_to_correlation(
      banded_mean_covariance, "unequal mean-difference"
    );
    standard_error_source = "unequal sample variances";
  } else {
    if (supplied_correlation.n_rows != p ||
        supplied_correlation.n_cols != p ||
        !supplied_correlation.is_finite()) {
      Rcpp::stop(
        "A finite p by p `supplied_correlation` matrix is required."
      );
    }
    correlation = supplied_correlation;
    if (supplied_standard_errors.n_elem == 0) {
      const arma::mat mean_covariance =
        covariance1 / static_cast<double>(n1) +
        covariance2 / static_cast<double>(n2);
      for (arma::uword j = 0; j < p; ++j) {
        standard_error_scaled(j) = std::sqrt(mean_covariance(j, j));
      }
      standard_error_source = "unequal sample variances";
    } else {
      if (supplied_standard_errors.n_elem != p ||
          !supplied_standard_errors.is_finite()) {
        Rcpp::stop(
          "`supplied_standard_errors` must contain p finite values."
        );
      }
      for (arma::uword j = 0; j < p; ++j) {
        const long double value = static_cast<long double>(
          supplied_standard_errors(j)
        ) / scaled.scale[j];
        standard_error_scaled(j) = aspu_checked_positive_double(
          value, "supplied standard error after internal scaling"
        );
      }
      standard_error_source = "supplied";
    }
  }

  arma::vec w(p);
  for (arma::uword j = 0; j < p; ++j) {
    const double standard_error = standard_error_scaled(j);
    if (!std::isfinite(standard_error) || standard_error <= 0.0) {
      Rcpp::stop(
        "Xu-Lin-Wei-Pan aSPU requires every marginal standard error to be "
        "finite and strictly positive; no variance floor is applied."
      );
    }
    w(j) = difference_scaled(j) / standard_error;
    if (!std::isfinite(w(j))) {
      Rcpp::stop(
        "Xu-Lin-Wei-Pan aSPU produced a non-finite standardised contrast."
      );
    }
  }

  arma::vec mean_x(p);
  arma::vec mean_y(p);
  arma::vec difference(p);
  arma::vec standard_error(p);
  arma::vec variance1_scaled = covariance1.diag();
  arma::vec variance2_scaled = covariance2.diag();
  arma::vec pooled_variance_scaled = pooled_covariance.diag();
  arma::vec column_scale(p);
  for (arma::uword j = 0; j < p; ++j) {
    const long double scale = scaled.scale[j];
    mean_x(j) = aspu_report_long_double(
      scaled.origin[j] +
      scale * static_cast<long double>(mean1_scaled(j))
    );
    mean_y(j) = aspu_report_long_double(
      scaled.origin[j] +
      scale * static_cast<long double>(mean2_scaled(j))
    );
    difference(j) = aspu_report_long_double(
      scale * static_cast<long double>(difference_scaled(j))
    );
    standard_error(j) = aspu_report_long_double(
      scale * static_cast<long double>(standard_error_scaled(j))
    );
    column_scale(j) = aspu_report_long_double(scale);
  }

  return Rcpp::List::create(
    Rcpp::Named("mean_x") = mean_x,
    Rcpp::Named("mean_y") = mean_y,
    Rcpp::Named("difference") = difference,
    Rcpp::Named("difference_scaled") = difference_scaled.t(),
    Rcpp::Named("standard_error") = standard_error,
    Rcpp::Named("standard_error_scaled") = standard_error_scaled,
    Rcpp::Named("W") = w,
    Rcpp::Named("correlation") = correlation,
    Rcpp::Named("variance1_scaled") = variance1_scaled,
    Rcpp::Named("variance2_scaled") = variance2_scaled,
    Rcpp::Named("pooled_variance_scaled") = pooled_variance_scaled,
    Rcpp::Named("column_scale") = column_scale,
    Rcpp::Named("standard_error_source") = standard_error_source,
    Rcpp::Named("correlation_source") = correlation_source,
    Rcpp::Named("bandwidth1") = bandwidth1,
    Rcpp::Named("bandwidth2") = bandwidth2,
    Rcpp::Named("n1") = static_cast<double>(n1),
    Rcpp::Named("n2") = static_cast<double>(n2),
    Rcpp::Named("p") = static_cast<double>(p)
  );
}


//' Gaussian finite-power moments for the analytical aSPU test
//'
//' @param scores Finite-power coordinate scores.
//' @param covariance Their null covariance matrix.
//' @param powers Positive finite integer powers.
//' @return Internal list of observed sums and Gaussian moment matrices.
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List cpp_aspu_power_moments(const arma::vec& scores,
                                  const arma::mat& covariance,
                                  const Rcpp::IntegerVector& powers) {
  const arma::uword p = scores.n_elem;
  const int m = powers.size();
  if (p < 1 || m < 1 || !scores.is_finite()) {
    Rcpp::stop(
      "The aSPU power kernel requires finite contrasts and at least one "
      "finite power."
    );
  }
  if (covariance.n_rows != p || covariance.n_cols != p ||
      !covariance.is_finite()) {
    Rcpp::stop("`covariance` must be a finite p by p matrix.");
  }
  const double matrix_tolerance = 1e-10;
  for (arma::uword j = 0; j < p; ++j) {
    if (covariance(j, j) <= 0.0) {
      Rcpp::stop(
        "`covariance` must have a finite, strictly positive diagonal."
      );
    }
    for (arma::uword k = 0; k < j; ++k) {
      const double scale = std::sqrt(covariance(j, j)) *
        std::sqrt(covariance(k, k));
      if (std::abs(covariance(j, k) - covariance(k, j)) >
          matrix_tolerance * scale ||
          std::abs(covariance(j, k)) >
          (1.0 + matrix_tolerance) * scale) {
        Rcpp::stop(
          "`covariance` must be symmetric and satisfy the covariance bound."
        );
      }
    }
  }

  std::vector<int> gamma(m);
  for (int r = 0; r < m; ++r) {
    gamma[r] = powers[r];
    if (gamma[r] < 1) {
      Rcpp::stop("Every finite aSPU power must be a positive integer.");
    }
  }

  arma::vec observed(m);
  arma::vec null_mean(m);
  arma::mat null_covariance(m, m, arma::fill::zeros);
  for (int r = 0; r < m; ++r) {
    AspuAccumulator observed_sum;
    for (arma::uword j = 0; j < p; ++j) {
      observed_sum.add(aspu_integer_power(
        static_cast<long double>(scores(j)), gamma[r]
      ));
    }
    observed(r) = aspu_checked_double(
      observed_sum.value, "observed finite-power statistic"
    );
    AspuAccumulator mean_sum;
    for (arma::uword j = 0; j < p; ++j) {
      mean_sum.add(aspu_univariate_normal_moment(
        gamma[r], static_cast<long double>(covariance(j, j))
      ));
    }
    null_mean(r) = aspu_checked_double(
      mean_sum.value, "finite-power null mean"
    );
  }

  for (int r = 0; r < m; ++r) {
    for (int s = r; s < m; ++s) {
      AspuAccumulator covariance_sum;
      for (arma::uword j = 0; j < p; ++j) {
        const long double variance_j = static_cast<long double>(
          covariance(j, j)
        );
        const long double mean_r_j = aspu_univariate_normal_moment(
          gamma[r], variance_j
        );
        const long double mean_s_j = aspu_univariate_normal_moment(
          gamma[s], variance_j
        );
        covariance_sum.add(
          aspu_bivariate_normal_moment(
            gamma[r], gamma[s], variance_j, variance_j, variance_j
          ) - mean_r_j * mean_s_j
        );
        for (arma::uword k = 0; k < j; ++k) {
          const long double variance_k = static_cast<long double>(
            covariance(k, k)
          );
          const long double covariance_jk = static_cast<long double>(
            covariance(j, k)
          );
          const long double mean_r_k = aspu_univariate_normal_moment(
            gamma[r], variance_k
          );
          const long double mean_s_k = aspu_univariate_normal_moment(
            gamma[s], variance_k
          );
          const long double forward = aspu_bivariate_normal_moment(
            gamma[r], gamma[s], variance_j, variance_k, covariance_jk
          ) - mean_r_j * mean_s_k;
          const long double reverse = aspu_bivariate_normal_moment(
            gamma[r], gamma[s], variance_k, variance_j, covariance_jk
          ) - mean_r_k * mean_s_j;
          covariance_sum.add(forward + reverse);
        }
      }
      const double covariance = aspu_checked_double(
        covariance_sum.value, "finite-power null covariance"
      );
      null_covariance(r, s) = covariance;
      null_covariance(s, r) = covariance;
    }
  }

  arma::vec null_variance = null_covariance.diag();
  arma::vec standardized(m);
  arma::mat null_correlation(m, m, arma::fill::eye);
  for (int r = 0; r < m; ++r) {
    if (!std::isfinite(null_variance(r)) || null_variance(r) <= 0.0) {
      Rcpp::stop(
        "Xu-Lin-Wei-Pan aSPU requires every finite-power null variance "
        "to be finite and strictly positive; no absolute-value repair or "
        "variance floor is applied."
      );
    }
    standardized(r) = (observed(r) - null_mean(r)) /
      std::sqrt(null_variance(r));
    if (!std::isfinite(standardized(r))) {
      Rcpp::stop(
        "Xu-Lin-Wei-Pan aSPU produced a non-finite standardised "
        "finite-power statistic."
      );
    }
  }
  for (int r = 0; r < m; ++r) {
    for (int s = 0; s < r; ++s) {
      const double value = null_covariance(r, s) /
        std::sqrt(null_variance(r) * null_variance(s));
      if (!std::isfinite(value)) {
        Rcpp::stop(
          "Xu-Lin-Wei-Pan aSPU produced a non-finite power correlation."
        );
      }
      null_correlation(r, s) = value;
      null_correlation(s, r) = value;
    }
  }

  return Rcpp::List::create(
    Rcpp::Named("powers") = powers,
    Rcpp::Named("observed") = observed,
    Rcpp::Named("null_mean") = null_mean,
    Rcpp::Named("null_variance") = null_variance,
    Rcpp::Named("standardized") = standardized,
    Rcpp::Named("null_covariance") = null_covariance,
    Rcpp::Named("null_correlation") = null_correlation
  );
}
