// Zhao--Feng (2026) strong-correlation one-sample spatial-sign bootstrap.
//
// The observed signs are centered at the null value.  Bootstrap signs are
// centered at the ordinary sample spatial median supplied by the R wrapper.
// The common theoretical factor sqrt(tau * choose(n, 2)) is omitted from
// both sides of the comparison, exactly as justified in the primary paper.

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp17)]]

#include <RcppArmadillo.h>

#include <algorithm>
#include <cstdint>
#include <cmath>
#include <limits>
#include <string>
#include <vector>

namespace {

void zfsc_neumaier_add(const long double value,
                       long double& total,
                       long double& correction) {
  const long double updated = total + value;
  if (std::abs(total) >= std::abs(value)) {
    correction += (total - updated) + value;
  } else {
    correction += (value - updated) + total;
  }
  total = updated;
}

double zfsc_checked_double(const long double value,
                           const char* quantity) {
  const long double maximum = static_cast<long double>(
    std::numeric_limits<double>::max()
  );
  if (!std::isfinite(value) || value > maximum || value < -maximum) {
    Rcpp::stop(
      "Zhao-Feng strong-correlation bootstrap produced a non-finite %s; "
      "no ridge, truncation, perturbation, absolute-value repair, or "
      "numerical floor is applied.",
      quantity
    );
  }
  return static_cast<double>(value);
}

struct ZfscDirections {
  arma::mat signs;
  arma::vec log_norms;
  arma::vec reported_norms;
  arma::uword zero_count;
  arma::uword overflow_fallback_count;
};

ZfscDirections zfsc_directions(const arma::mat& x,
                               const arma::vec& center) {
  const arma::uword n = x.n_rows;
  const arma::uword p = x.n_cols;
  const double maximum_log = std::log(
    std::numeric_limits<double>::max()
  );
  const double minimum_log = std::log(
    std::numeric_limits<double>::denorm_min()
  );
  ZfscDirections out;
  out.signs.zeros(n, p);
  out.log_norms.set_size(n);
  out.log_norms.fill(-std::numeric_limits<double>::infinity());
  out.reported_norms.zeros(n);
  out.zero_count = 0u;
  out.overflow_fallback_count = 0u;

  std::vector<long double> difference(p);
  std::vector<long double> scaled(p);
  for (arma::uword i = 0; i < n; ++i) {
    long double maximum = 0.0L;
    bool direct_ok = true;
    for (arma::uword j = 0; j < p; ++j) {
      difference[j] = static_cast<long double>(x(i, j)) -
        static_cast<long double>(center(j));
      direct_ok = direct_ok && std::isfinite(difference[j]);
      maximum = std::max(maximum, std::abs(difference[j]));
    }

    long double log_base_scale = 0.0L;
    if (!direct_ok || !std::isfinite(maximum)) {
      ++out.overflow_fallback_count;
      double operand_scale = 0.0;
      for (arma::uword j = 0; j < p; ++j) {
        operand_scale = std::max(operand_scale, std::abs(x(i, j)));
        operand_scale = std::max(operand_scale, std::abs(center(j)));
      }
      if (!(operand_scale > 0.0) || !std::isfinite(operand_scale)) {
        Rcpp::stop(
          "Zhao-Feng could not construct a finite residual operand scale."
        );
      }
      const long double operand_scale_ld =
        static_cast<long double>(operand_scale);
      maximum = 0.0L;
      for (arma::uword j = 0; j < p; ++j) {
        difference[j] =
          static_cast<long double>(x(i, j)) / operand_scale_ld -
          static_cast<long double>(center(j)) / operand_scale_ld;
        if (!std::isfinite(difference[j])) {
          Rcpp::stop(
            "Zhao-Feng residual subtraction remains non-finite after "
            "common-operand scaling."
          );
        }
        maximum = std::max(maximum, std::abs(difference[j]));
      }
      log_base_scale = std::log(operand_scale_ld);
    }

    if (maximum == 0.0L) {
      ++out.zero_count;
      continue;
    }
    if (!std::isfinite(maximum) || maximum < 0.0L) {
      Rcpp::stop("Zhao-Feng encountered an invalid residual scale.");
    }

    long double square_total = 0.0L;
    long double square_correction = 0.0L;
    for (arma::uword j = 0; j < p; ++j) {
      scaled[j] = difference[j] / maximum;
      zfsc_neumaier_add(
        scaled[j] * scaled[j], square_total, square_correction
      );
    }
    const long double square_norm = square_total + square_correction;
    if (!std::isfinite(square_norm) || square_norm <= 0.0L) {
      Rcpp::stop(
        "Zhao-Feng could not normalize a nonzero residual; no norm floor "
        "is applied."
      );
    }
    const long double normalized_norm = std::sqrt(square_norm);
    for (arma::uword j = 0; j < p; ++j) {
      out.signs(i, j) = zfsc_checked_double(
        scaled[j] / normalized_norm, "spatial-sign coordinate"
      );
    }
    const long double log_norm = log_base_scale + std::log(maximum) +
      std::log(normalized_norm);
    out.log_norms(i) = zfsc_checked_double(
      log_norm, "log residual norm"
    );
    if (out.log_norms(i) > maximum_log) {
      out.reported_norms(i) = R_PosInf;
    } else if (out.log_norms(i) < minimum_log) {
      out.reported_norms(i) = 0.0;
    } else {
      out.reported_norms(i) = std::exp(out.log_norms(i));
    }
  }
  return out;
}

struct ZfscPairSum {
  long double raw;
  long double diagonal;
  arma::vec sign_sum;
};

ZfscPairSum zfsc_pair_sum(const arma::mat& signs) {
  const arma::uword n = signs.n_rows;
  const arma::uword p = signs.n_cols;
  ZfscPairSum out;
  out.raw = 0.0L;
  out.diagonal = 0.0L;
  out.sign_sum.zeros(p);
  std::vector<long double> sum(p, 0.0L);
  std::vector<long double> correction(p, 0.0L);
  long double diagonal_correction = 0.0L;

  for (arma::uword i = 0; i < n; ++i) {
    long double row_square = 0.0L;
    long double row_square_correction = 0.0L;
    for (arma::uword j = 0; j < p; ++j) {
      const long double value = static_cast<long double>(signs(i, j));
      zfsc_neumaier_add(value, sum[j], correction[j]);
      zfsc_neumaier_add(
        value * value, row_square, row_square_correction
      );
    }
    zfsc_neumaier_add(
      row_square + row_square_correction,
      out.diagonal, diagonal_correction
    );
  }
  out.diagonal += diagonal_correction;

  long double sum_square = 0.0L;
  long double sum_square_correction = 0.0L;
  for (arma::uword j = 0; j < p; ++j) {
    const long double value = sum[j] + correction[j];
    out.sign_sum(j) = zfsc_checked_double(value, "summed sign coordinate");
    zfsc_neumaier_add(
      value * value, sum_square, sum_square_correction
    );
  }
  out.raw = 0.5L *
    (sum_square + sum_square_correction - out.diagonal);
  if (!std::isfinite(out.raw)) {
    Rcpp::stop("Zhao-Feng observed pair sum is non-finite.");
  }
  return out;
}

struct ZfscBootstrapSummary {
  std::uint64_t exceedances;
  long double mean;
  long double variance;
  long double minimum;
  long double maximum;
};

} // namespace


// [[Rcpp::export]]
Rcpp::List cpp_zhao_feng_strongcorr_sign_bootstrap(
    const arma::mat& x, const arma::vec& mu,
    const arma::vec& fitted_location, const std::string multiplier,
    const int B, const double alpha, const bool keep_bootstrap) {
  if (!x.is_finite() || !mu.is_finite() || !fitted_location.is_finite()) {
    Rcpp::stop("Zhao-Feng inputs and fitted location must be finite.");
  }
  if (x.n_rows < 2u || x.n_cols < 1u) {
    Rcpp::stop(
      "Zhao-Feng requires at least two observations and one variable."
    );
  }
  if (mu.n_elem != x.n_cols || fitted_location.n_elem != x.n_cols) {
    Rcpp::stop("`mu` and `fitted_location` must match the columns of `x`.");
  }
  if (multiplier != "rademacher" && multiplier != "gaussian") {
    Rcpp::stop("`multiplier` must be `\"rademacher\"` or `\"gaussian\"`.");
  }
  if (B < 1) Rcpp::stop("`B` must be a positive integer.");
  if (!std::isfinite(alpha) || alpha <= 0.0 || alpha >= 1.0) {
    Rcpp::stop("`alpha` must be strictly between zero and one.");
  }

  const arma::uword n = x.n_rows;
  const arma::uword p = x.n_cols;
  const long double pair_count =
    static_cast<long double>(n) * static_cast<long double>(n - 1u) /
    2.0L;
  const long double root_pair_count = std::sqrt(pair_count);
  const ZfscDirections observed_directions = zfsc_directions(x, mu);
  const ZfscDirections fitted_directions =
    zfsc_directions(x, fitted_location);
  const ZfscPairSum observed = zfsc_pair_sum(
    observed_directions.signs
  );
  const ZfscPairSum fitted = zfsc_pair_sum(
    fitted_directions.signs
  );

  std::vector<long double> sign_norm_square(n, 0.0L);
  for (arma::uword i = 0; i < n; ++i) {
    long double total = 0.0L;
    long double correction = 0.0L;
    for (arma::uword j = 0; j < p; ++j) {
      const long double value = static_cast<long double>(
        fitted_directions.signs(i, j)
      );
      zfsc_neumaier_add(value * value, total, correction);
    }
    sign_norm_square[i] = total + correction;
  }

  std::vector<long double> bootstrap_raw(static_cast<std::size_t>(B));
  std::vector<long double> weighted_sum(p, 0.0L);
  std::vector<long double> weighted_correction(p, 0.0L);
  std::uint64_t exceedances = 0u;
  long double bootstrap_total = 0.0L;
  long double bootstrap_total_correction = 0.0L;
  long double bootstrap_minimum =
    std::numeric_limits<long double>::infinity();
  long double bootstrap_maximum =
    -std::numeric_limits<long double>::infinity();

  for (int b = 0; b < B; ++b) {
    if ((b & 1023) == 0) Rcpp::checkUserInterrupt();
    std::fill(weighted_sum.begin(), weighted_sum.end(), 0.0L);
    std::fill(
      weighted_correction.begin(), weighted_correction.end(), 0.0L
    );
    long double diagonal = 0.0L;
    long double diagonal_correction = 0.0L;

    for (arma::uword i = 0; i < n; ++i) {
      const double draw = multiplier == "rademacher" ?
        (R::unif_rand() < 0.5 ? -1.0 : 1.0) : R::norm_rand();
      const long double draw_ld = static_cast<long double>(draw);
      for (arma::uword j = 0; j < p; ++j) {
        zfsc_neumaier_add(
          draw_ld * static_cast<long double>(
            fitted_directions.signs(i, j)
          ),
          weighted_sum[j], weighted_correction[j]
        );
      }
      zfsc_neumaier_add(
        draw_ld * draw_ld * sign_norm_square[i],
        diagonal, diagonal_correction
      );
    }

    long double weighted_norm_square = 0.0L;
    long double weighted_norm_correction = 0.0L;
    for (arma::uword j = 0; j < p; ++j) {
      const long double value =
        weighted_sum[j] + weighted_correction[j];
      zfsc_neumaier_add(
        value * value,
        weighted_norm_square, weighted_norm_correction
      );
    }
    const long double statistic = 0.5L *
      (weighted_norm_square + weighted_norm_correction -
       diagonal - diagonal_correction);
    if (!std::isfinite(statistic)) {
      Rcpp::stop(
        "Zhao-Feng bootstrap statistic is non-finite; no truncation is "
        "applied."
      );
    }
    bootstrap_raw[static_cast<std::size_t>(b)] = statistic;
    if (statistic >= observed.raw) ++exceedances;
    zfsc_neumaier_add(
      statistic, bootstrap_total, bootstrap_total_correction
    );
    bootstrap_minimum = std::min(bootstrap_minimum, statistic);
    bootstrap_maximum = std::max(bootstrap_maximum, statistic);
  }

  const long double bootstrap_mean =
    (bootstrap_total + bootstrap_total_correction) /
    static_cast<long double>(B);
  long double variance_total = 0.0L;
  long double variance_correction = 0.0L;
  for (const long double value : bootstrap_raw) {
    const long double difference = value - bootstrap_mean;
    zfsc_neumaier_add(
      difference * difference, variance_total, variance_correction
    );
  }
  const long double bootstrap_variance =
    (variance_total + variance_correction) /
    static_cast<long double>(B);

  int critical_index = static_cast<int>(std::ceil(
    (1.0 - alpha) * static_cast<double>(B)
  ));
  critical_index = std::max(1, std::min(B, critical_index));
  std::vector<long double> ordered = bootstrap_raw;
  std::nth_element(
    ordered.begin(), ordered.begin() + (critical_index - 1), ordered.end()
  );
  const long double critical =
    ordered[static_cast<std::size_t>(critical_index - 1)];
  const double p_value =
    (1.0 + static_cast<double>(exceedances)) /
    (1.0 + static_cast<double>(B));
  const double empirical_tail = static_cast<double>(exceedances) /
    static_cast<double>(B);
  const double mc_standard_error = std::sqrt(
    empirical_tail * (1.0 - empirical_tail) / static_cast<double>(B)
  );
  const bool reject_by_critical = observed.raw > critical;
  const bool reject_by_p_value = p_value <= alpha;

  Rcpp::RObject retained_raw = R_NilValue;
  Rcpp::RObject retained_normalized = R_NilValue;
  if (keep_bootstrap) {
    Rcpp::NumericVector raw(B);
    Rcpp::NumericVector normalized(B);
    for (int b = 0; b < B; ++b) {
      raw[b] = zfsc_checked_double(
        bootstrap_raw[static_cast<std::size_t>(b)],
        "retained raw bootstrap statistic"
      );
      normalized[b] = zfsc_checked_double(
        bootstrap_raw[static_cast<std::size_t>(b)] / root_pair_count,
        "retained pair-normalized bootstrap statistic"
      );
    }
    retained_raw = raw;
    retained_normalized = normalized;
  }

  return Rcpp::List::create(
    Rcpp::Named("observed_raw") = zfsc_checked_double(
      observed.raw, "observed raw pair sum"
    ),
    Rcpp::Named("observed_pair_normalized") = zfsc_checked_double(
      observed.raw / root_pair_count,
      "observed pair-normalized statistic"
    ),
    Rcpp::Named("observed_sign_diagonal") = zfsc_checked_double(
      observed.diagonal, "observed sign diagonal sum"
    ),
    Rcpp::Named("observed_sign_sum") = observed.sign_sum,
    Rcpp::Named("observed_signs") = observed_directions.signs,
    Rcpp::Named("observed_residual_norms") =
      observed_directions.reported_norms,
    Rcpp::Named("observed_log_residual_norms") =
      observed_directions.log_norms,
    Rcpp::Named("observed_zero_residuals") = static_cast<double>(
      observed_directions.zero_count
    ),
    Rcpp::Named("observed_overflow_fallback_rows") =
      static_cast<double>(
        observed_directions.overflow_fallback_count
      ),
    Rcpp::Named("fitted_raw_pair_sum") = zfsc_checked_double(
      fitted.raw, "fitted-sign raw pair sum"
    ),
    Rcpp::Named("fitted_sign_diagonal") = zfsc_checked_double(
      fitted.diagonal, "fitted sign diagonal sum"
    ),
    Rcpp::Named("fitted_sign_sum") = fitted.sign_sum,
    Rcpp::Named("fitted_signs") = fitted_directions.signs,
    Rcpp::Named("fitted_residual_norms") =
      fitted_directions.reported_norms,
    Rcpp::Named("fitted_log_residual_norms") =
      fitted_directions.log_norms,
    Rcpp::Named("fitted_zero_residuals") = static_cast<double>(
      fitted_directions.zero_count
    ),
    Rcpp::Named("fitted_overflow_fallback_rows") =
      static_cast<double>(
        fitted_directions.overflow_fallback_count
      ),
    Rcpp::Named("bootstrap_statistics_raw") = retained_raw,
    Rcpp::Named("bootstrap_statistics_pair_normalized") =
      retained_normalized,
    Rcpp::Named("bootstrap_mean_raw") = zfsc_checked_double(
      bootstrap_mean, "bootstrap mean"
    ),
    Rcpp::Named("bootstrap_variance_raw") = zfsc_checked_double(
      bootstrap_variance, "bootstrap variance"
    ),
    Rcpp::Named("bootstrap_minimum_raw") = zfsc_checked_double(
      bootstrap_minimum, "bootstrap minimum"
    ),
    Rcpp::Named("bootstrap_maximum_raw") = zfsc_checked_double(
      bootstrap_maximum, "bootstrap maximum"
    ),
    Rcpp::Named("critical_value_raw") = zfsc_checked_double(
      critical, "bootstrap critical value"
    ),
    Rcpp::Named("critical_value_pair_normalized") = zfsc_checked_double(
      critical / root_pair_count,
      "pair-normalized bootstrap critical value"
    ),
    Rcpp::Named("critical_index") = critical_index,
    Rcpp::Named("p_value_plus_one") = p_value,
    Rcpp::Named("empirical_tail") = empirical_tail,
    Rcpp::Named("exceedances") = static_cast<double>(exceedances),
    Rcpp::Named("mc_standard_error") = mc_standard_error,
    Rcpp::Named("reject_by_critical") = reject_by_critical,
    Rcpp::Named("reject_by_p_value") = reject_by_p_value,
    Rcpp::Named("decisions_differ") =
      reject_by_critical != reject_by_p_value,
    Rcpp::Named("randomization_degenerate") =
      bootstrap_minimum == bootstrap_maximum,
    Rcpp::Named("pair_count") = static_cast<double>(pair_count),
    Rcpp::Named("root_pair_count") = static_cast<double>(root_pair_count),
    Rcpp::Named("multiplier") = multiplier,
    Rcpp::Named("B") = B,
    Rcpp::Named("alpha") = alpha,
    Rcpp::Named("n") = static_cast<double>(n),
    Rcpp::Named("p") = static_cast<double>(p)
  );
}
