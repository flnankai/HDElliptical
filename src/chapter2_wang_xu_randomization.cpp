// Wang--Xu (2022) approximate randomization test for the high-dimensional
// two-sample Behrens--Fisher problem.
//
// The observed statistic is the full-sample Chen--Qin U statistic.  Its
// conditional reference distribution is built from adjacent within-group
// half-differences and Rademacher signs; it is not a permutation of pooled
// group labels.  All helpers carry a module prefix because every package C++
// translation unit is linked into one shared library.

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

struct WxPrepared {
  arma::mat x_scaled;
  arma::mat y_scaled;
  arma::mat pseudo_x_scaled;
  arma::mat pseudo_y_scaled;
  arma::vec mean_x;
  arma::vec mean_y;
  arma::vec difference;
  arma::vec anchor;
  double global_scale;
  bool subtraction_overflow_fallback;
  long double mean_difference_squared_scaled;
  long double trace_s1_scaled;
  long double trace_s2_scaled;
  long double observed_scaled;
};

struct WxRandomSummary {
  std::uint64_t count = 0U;
  std::uint64_t exceedances = 0U;
  long double mean = 0.0L;
  long double m2 = 0.0L;
  long double minimum = std::numeric_limits<long double>::infinity();
  long double maximum = -std::numeric_limits<long double>::infinity();

  void add(const long double value, const long double observed) {
    ++count;
    if (value >= observed) ++exceedances;
    const long double delta = value - mean;
    mean += delta / static_cast<long double>(count);
    m2 += delta * (value - mean);
    minimum = std::min(minimum, value);
    maximum = std::max(maximum, value);
  }

  long double variance() const {
    return count > 1U ? m2 / static_cast<long double>(count - 1U) : 0.0L;
  }
};

void wx_require_finite_matrix(const arma::mat& value, const char* name) {
  if (!value.is_finite()) {
    Rcpp::stop("`%s` must contain only finite values.", name);
  }
}

double wx_checked_double(const long double value, const char* quantity) {
  const double answer = static_cast<double>(value);
  if (!std::isfinite(answer)) {
    Rcpp::stop("%s is not representable as one finite double.", quantity);
  }
  return answer;
}

double wx_report_double(const long double value) {
  const long double largest = static_cast<long double>(
    std::numeric_limits<double>::max()
  );
  if (value > largest) return R_PosInf;
  if (value < -largest) return R_NegInf;
  return static_cast<double>(value);
}

long double wx_scaled_anchor_difference(const double value,
                                        const double anchor,
                                        const long double scale) {
  const long double raw = static_cast<long double>(value) -
    static_cast<long double>(anchor);
  if (std::isfinite(raw)) return raw / scale;
  return static_cast<long double>(value) / scale -
    static_cast<long double>(anchor) / scale;
}

double wx_choose_global_scale(const arma::mat& x, const arma::mat& y,
                              bool& fallback) {
  fallback = false;
  long double maximum_difference = 0.0L;
  bool all_direct_differences_finite = true;
  for (arma::uword j = 0; j < x.n_cols; ++j) {
    const long double anchor = static_cast<long double>(x(0, j));
    for (arma::uword i = 0; i < x.n_rows; ++i) {
      const long double difference = static_cast<long double>(x(i, j)) - anchor;
      if (!std::isfinite(difference)) {
        all_direct_differences_finite = false;
      } else {
        maximum_difference = std::max(maximum_difference, std::abs(difference));
      }
    }
    for (arma::uword i = 0; i < y.n_rows; ++i) {
      const long double difference = static_cast<long double>(y(i, j)) - anchor;
      if (!std::isfinite(difference)) {
        all_direct_differences_finite = false;
      } else {
        maximum_difference = std::max(maximum_difference, std::abs(difference));
      }
    }
  }

  if (all_direct_differences_finite &&
      maximum_difference <= static_cast<long double>(
        std::numeric_limits<double>::max()
      )) {
    const double answer = static_cast<double>(maximum_difference);
    return answer > 0.0 ? answer : 1.0;
  }

  // On ABIs where long double has double range, subtraction of opposite
  // finite extremes can overflow.  Scale operands separately in that case.
  fallback = true;
  double operand_scale = 0.0;
  for (arma::uword j = 0; j < x.n_cols; ++j) {
    for (arma::uword i = 0; i < x.n_rows; ++i) {
      operand_scale = std::max(operand_scale, std::abs(x(i, j)));
    }
    for (arma::uword i = 0; i < y.n_rows; ++i) {
      operand_scale = std::max(operand_scale, std::abs(y(i, j)));
    }
  }
  return operand_scale > 0.0 ? operand_scale : 1.0;
}

double wx_rescale_second_order(const long double value,
                               const double scale,
                               const char* quantity) {
  if (value == 0.0L) return 0.0;
  const long double scale_ld = static_cast<long double>(scale);
  const long double log_magnitude = std::log(std::abs(value)) +
    2.0L * std::log(scale_ld);
  const long double log_largest = std::log(
    static_cast<long double>(std::numeric_limits<double>::max())
  );
  if (!std::isfinite(log_magnitude) || log_magnitude > log_largest) {
    Rcpp::stop("%s is not representable as one finite double.", quantity);
  }
  const long double result = value * scale_ld * scale_ld;
  return wx_checked_double(result, quantity);
}

double wx_report_second_order(const long double value, const double scale) {
  const long double scale_ld = static_cast<long double>(scale);
  return wx_report_double(value * scale_ld * scale_ld);
}

double wx_report_fourth_order(const long double value, const double scale) {
  const long double scale_ld = static_cast<long double>(scale);
  const long double scale2 = scale_ld * scale_ld;
  return wx_report_double(value * scale2 * scale2);
}

WxPrepared wx_prepare(const arma::mat& x, const arma::mat& y) {
  WxPrepared out;
  const arma::uword n1 = x.n_rows;
  const arma::uword n2 = y.n_rows;
  const arma::uword p = x.n_cols;
  const arma::uword m1 = n1 / 2U;
  const arma::uword m2 = n2 / 2U;
  const long double n1_ld = static_cast<long double>(n1);
  const long double n2_ld = static_cast<long double>(n2);

  out.global_scale = wx_choose_global_scale(
    x, y, out.subtraction_overflow_fallback
  );
  const long double scale_ld = static_cast<long double>(out.global_scale);
  out.anchor = x.row(0).t();
  out.x_scaled.set_size(n1, p);
  out.y_scaled.set_size(n2, p);
  out.mean_x.set_size(p);
  out.mean_y.set_size(p);
  out.difference.set_size(p);

  std::vector<long double> mean_x_scaled(p, 0.0L);
  std::vector<long double> mean_y_scaled(p, 0.0L);
  for (arma::uword j = 0; j < p; ++j) {
    const double anchor = x(0, j);
    for (arma::uword i = 0; i < n1; ++i) {
      const long double value = wx_scaled_anchor_difference(
        x(i, j), anchor, scale_ld
      );
      out.x_scaled(i, j) = wx_checked_double(
        value, "a scaled first-group anchored observation"
      );
      mean_x_scaled[j] += value;
    }
    for (arma::uword i = 0; i < n2; ++i) {
      const long double value = wx_scaled_anchor_difference(
        y(i, j), anchor, scale_ld
      );
      out.y_scaled(i, j) = wx_checked_double(
        value, "a scaled second-group anchored observation"
      );
      mean_y_scaled[j] += value;
    }
    mean_x_scaled[j] /= n1_ld;
    mean_y_scaled[j] /= n2_ld;
    const long double anchor_ld = static_cast<long double>(anchor);
    out.mean_x(j) = wx_checked_double(
      anchor_ld + mean_x_scaled[j] * scale_ld,
      "a first-group sample mean"
    );
    out.mean_y(j) = wx_checked_double(
      anchor_ld + mean_y_scaled[j] * scale_ld,
      "a second-group sample mean"
    );
    out.difference(j) = wx_checked_double(
      (mean_x_scaled[j] - mean_y_scaled[j]) * scale_ld,
      "a sample-mean difference"
    );
  }

  out.mean_difference_squared_scaled = 0.0L;
  out.trace_s1_scaled = 0.0L;
  out.trace_s2_scaled = 0.0L;
  for (arma::uword j = 0; j < p; ++j) {
    const long double difference = mean_x_scaled[j] - mean_y_scaled[j];
    out.mean_difference_squared_scaled += difference * difference;
    for (arma::uword i = 0; i < n1; ++i) {
      const long double centered =
        static_cast<long double>(out.x_scaled(i, j)) - mean_x_scaled[j];
      out.trace_s1_scaled += centered * centered /
        static_cast<long double>(n1 - 1U);
    }
    for (arma::uword i = 0; i < n2; ++i) {
      const long double centered =
        static_cast<long double>(out.y_scaled(i, j)) - mean_y_scaled[j];
      out.trace_s2_scaled += centered * centered /
        static_cast<long double>(n2 - 1U);
    }
  }
  out.observed_scaled = out.mean_difference_squared_scaled -
    out.trace_s1_scaled / n1_ld - out.trace_s2_scaled / n2_ld;
  if (!std::isfinite(out.observed_scaled)) {
    Rcpp::stop("Wang-Xu produced a non-finite internally scaled CQ statistic.");
  }

  out.pseudo_x_scaled.set_size(m1, p);
  out.pseudo_y_scaled.set_size(m2, p);
  for (arma::uword i = 0; i < m1; ++i) {
    for (arma::uword j = 0; j < p; ++j) {
      out.pseudo_x_scaled(i, j) =
        (out.x_scaled(2U * i + 1U, j) - out.x_scaled(2U * i, j)) / 2.0;
    }
  }
  for (arma::uword i = 0; i < m2; ++i) {
    for (arma::uword j = 0; j < p; ++j) {
      out.pseudo_y_scaled(i, j) =
        (out.y_scaled(2U * i + 1U, j) - out.y_scaled(2U * i, j)) / 2.0;
    }
  }
  return out;
}

arma::mat wx_quadratic_kernel(const arma::mat& pseudo_x,
                              const arma::mat& pseudo_y) {
  const arma::uword m1 = pseudo_x.n_rows;
  const arma::uword m2 = pseudo_y.n_rows;
  const arma::uword p = pseudo_x.n_cols;
  const arma::uword total = m1 + m2;
  arma::mat kernel(total, total, arma::fill::zeros);
  const long double within1_denominator =
    static_cast<long double>(m1) * static_cast<long double>(m1 - 1U);
  const long double within2_denominator =
    static_cast<long double>(m2) * static_cast<long double>(m2 - 1U);
  const long double cross_denominator =
    static_cast<long double>(m1) * static_cast<long double>(m2);

  for (arma::uword i = 0; i + 1U < m1; ++i) {
    for (arma::uword j = i + 1U; j < m1; ++j) {
      long double inner = 0.0L;
      for (arma::uword q = 0; q < p; ++q) {
        inner += static_cast<long double>(pseudo_x(i, q)) *
          static_cast<long double>(pseudo_x(j, q));
      }
      const double value = wx_checked_double(
        inner / within1_denominator,
        "a first-group randomized quadratic-kernel entry"
      );
      kernel(i, j) = kernel(j, i) = value;
    }
  }
  for (arma::uword i = 0; i + 1U < m2; ++i) {
    for (arma::uword j = i + 1U; j < m2; ++j) {
      long double inner = 0.0L;
      for (arma::uword q = 0; q < p; ++q) {
        inner += static_cast<long double>(pseudo_y(i, q)) *
          static_cast<long double>(pseudo_y(j, q));
      }
      const double value = wx_checked_double(
        inner / within2_denominator,
        "a second-group randomized quadratic-kernel entry"
      );
      kernel(m1 + i, m1 + j) = kernel(m1 + j, m1 + i) = value;
    }
  }
  for (arma::uword i = 0; i < m1; ++i) {
    for (arma::uword j = 0; j < m2; ++j) {
      long double inner = 0.0L;
      for (arma::uword q = 0; q < p; ++q) {
        inner += static_cast<long double>(pseudo_x(i, q)) *
          static_cast<long double>(pseudo_y(j, q));
      }
      const double value = wx_checked_double(
        -inner / cross_denominator,
        "a cross-group randomized quadratic-kernel entry"
      );
      kernel(i, m1 + j) = kernel(m1 + j, i) = value;
    }
  }
  return kernel;
}

long double wx_evaluate_quadratic(const arma::mat& kernel,
                                  const std::vector<int>& signs) {
  long double answer = 0.0L;
  for (arma::uword i = 0; i + 1U < kernel.n_rows; ++i) {
    for (arma::uword j = i + 1U; j < kernel.n_cols; ++j) {
      answer += 2.0L * static_cast<long double>(kernel(i, j)) *
        static_cast<long double>(signs[i]) *
        static_cast<long double>(signs[j]);
    }
  }
  return answer;
}

// A counter-based 32-bit finalizer.  Each (replicate, sign coordinate) pair
// maps to its own counter, so the stream does not depend on loop chunking.
// The first release deliberately accepts only workers = 1 at the R layer;
// retaining counter semantics makes future parallelization reproducible.
std::uint32_t wx_mix32(std::uint32_t value) {
  value ^= value >> 16U;
  value *= UINT32_C(0x7feb352d);
  value ^= value >> 15U;
  value *= UINT32_C(0x846ca68b);
  value ^= value >> 16U;
  return value;
}

std::uint32_t wx_counter_word(const std::uint64_t counter,
                              const std::uint32_t seed) {
  const std::uint32_t low = static_cast<std::uint32_t>(counter);
  const std::uint32_t high = static_cast<std::uint32_t>(counter >> 32U);
  const std::uint32_t high_key = wx_mix32(
    high + UINT32_C(0x9e3779b9)
  );
  const std::uint32_t keyed = (seed ^ high_key) +
    UINT32_C(0x9e3779b9) * (low + UINT32_C(1));
  return wx_mix32(keyed);
}

Rcpp::IntegerMatrix wx_pair_indices(const arma::uword pairs) {
  Rcpp::IntegerMatrix answer(static_cast<int>(pairs), 2);
  for (arma::uword i = 0; i < pairs; ++i) {
    answer(i, 0) = static_cast<int>(2U * i + 1U);
    answer(i, 1) = static_cast<int>(2U * i + 2U);
  }
  return answer;
}

arma::mat wx_restore_pseudo_units(const arma::mat& value,
                                  const double scale) {
  arma::mat answer(value.n_rows, value.n_cols);
  const long double scale_ld = static_cast<long double>(scale);
  for (arma::uword i = 0; i < value.n_rows; ++i) {
    for (arma::uword j = 0; j < value.n_cols; ++j) {
      answer(i, j) = wx_checked_double(
        static_cast<long double>(value(i, j)) * scale_ld,
        "a Wang-Xu pseudo-observation"
      );
    }
  }
  return answer;
}

}  // namespace


//' Wang--Xu approximate-randomization kernel
//'
//' @param x,y Numeric observation-by-variable matrices.
//' @param calibration Either `"exact"` or `"monte_carlo"`; `"auto"` is
//'   resolved by the public R wrapper.
//' @param B Positive Monte Carlo draw count (ignored for exact enumeration).
//' @param seed Integer-valued counter-generator seed represented as a double.
//' @param max_exact Maximum number of global-sign-reduced exact patterns.
//' @param keep_randomized Whether to retain statistics and sign patterns.
//' @return Internal observed, pseudo-sample, reference, and diagnostic fields.
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List cpp_wang_xu_approx_randomization(
    const arma::mat& x, const arma::mat& y, const std::string calibration,
    const int B, const double seed, const double max_exact,
    const bool keep_randomized) {
  wx_require_finite_matrix(x, "x");
  wx_require_finite_matrix(y, "y");
  if (x.n_cols < 1U || y.n_cols != x.n_cols) {
    Rcpp::stop("`x` and `y` must have the same positive number of columns.");
  }
  if (x.n_rows < 4U || y.n_rows < 4U) {
    Rcpp::stop(
      "Wang-Xu calibration requires at least four observations in each "
      "group so each pseudo-sample contains at least two half-differences."
    );
  }
  if (calibration != "exact" && calibration != "monte_carlo") {
    Rcpp::stop("`calibration` must be `\"exact\"` or `\"monte_carlo\"`.");
  }
  if (B < 1) Rcpp::stop("`B` must be a positive integer.");
  if (!std::isfinite(seed) || seed < 0.0 || seed > 4294967295.0 ||
      seed != std::floor(seed)) {
    Rcpp::stop("`seed` must be an integer in [0, 2^32 - 1].");
  }
  if (!std::isfinite(max_exact) || max_exact < 1.0 ||
      max_exact != std::floor(max_exact)) {
    Rcpp::stop("`max_exact` must be a positive integer-valued number.");
  }

  const arma::uword n1 = x.n_rows;
  const arma::uword n2 = y.n_rows;
  const arma::uword p = x.n_cols;
  const arma::uword m1 = n1 / 2U;
  const arma::uword m2 = n2 / 2U;
  const arma::uword sign_dimension = m1 + m2;
  const WxPrepared prepared = wx_prepare(x, y);
  const arma::mat kernel = wx_quadratic_kernel(
    prepared.pseudo_x_scaled, prepared.pseudo_y_scaled
  );
  const double observed = wx_rescale_second_order(
    prepared.observed_scaled, prepared.global_scale,
    "the observed Wang-Xu Chen-Qin statistic"
  );
  const bool observed_underflow = prepared.observed_scaled != 0.0L &&
    observed == 0.0;

  std::uint64_t evaluations = 0U;
  if (calibration == "exact") {
    evaluations = 1U;
    const std::uint64_t maximum = static_cast<std::uint64_t>(max_exact);
    for (arma::uword j = 1U; j < sign_dimension; ++j) {
      if (evaluations > maximum / 2U) {
        Rcpp::stop(
          "Exact Wang-Xu enumeration needs more than `max_exact` "
          "global-sign-reduced patterns; use Monte Carlo calibration or "
          "increase `max_exact` explicitly."
        );
      }
      evaluations *= 2U;
    }
  } else {
    evaluations = static_cast<std::uint64_t>(B);
  }
  if (evaluations > static_cast<std::uint64_t>(
        std::numeric_limits<int>::max()
      ) && keep_randomized) {
    Rcpp::stop(
      "`keep_randomized = TRUE` cannot retain more than the R matrix row "
      "limit; reduce the reference size or do not retain all draws."
    );
  }

  Rcpp::NumericVector retained_scaled;
  Rcpp::NumericVector retained_physical;
  Rcpp::IntegerMatrix retained_signs;
  if (keep_randomized) {
    retained_scaled = Rcpp::NumericVector(static_cast<R_xlen_t>(evaluations));
    retained_physical = Rcpp::NumericVector(static_cast<R_xlen_t>(evaluations));
    retained_signs = Rcpp::IntegerMatrix(
      static_cast<int>(evaluations), static_cast<int>(sign_dimension)
    );
  }

  WxRandomSummary summary;
  std::vector<int> signs(sign_dimension, -1);
  const std::uint32_t seed_u32 = static_cast<std::uint32_t>(seed);
  for (std::uint64_t replicate = 0U; replicate < evaluations; ++replicate) {
    if ((replicate & UINT64_C(1023)) == 0U) Rcpp::checkUserInterrupt();
    if (calibration == "exact") {
      // T(E) = T(-E), so fixing the first sign to +1 enumerates one member
      // of each two-element global-sign orbit without changing its law.
      signs[0] = 1;
      for (arma::uword j = 1U; j < sign_dimension; ++j) {
        const bool bit = ((replicate >> (j - 1U)) & UINT64_C(1)) != 0U;
        signs[j] = bit ? 1 : -1;
      }
    } else {
      for (arma::uword j = 0U; j < sign_dimension; ++j) {
        const std::uint64_t counter =
          replicate * static_cast<std::uint64_t>(sign_dimension) +
          static_cast<std::uint64_t>(j);
        signs[j] = (wx_counter_word(counter, seed_u32) & UINT32_C(1)) != 0U ?
          1 : -1;
      }
    }

    const long double statistic = wx_evaluate_quadratic(kernel, signs);
    summary.add(statistic, prepared.observed_scaled);
    if (keep_randomized) {
      retained_scaled[static_cast<R_xlen_t>(replicate)] =
        wx_checked_double(statistic, "a scaled randomized statistic");
      retained_physical[static_cast<R_xlen_t>(replicate)] =
        wx_report_second_order(statistic, prepared.global_scale);
      for (arma::uword j = 0U; j < sign_dimension; ++j) {
        retained_signs(static_cast<int>(replicate), static_cast<int>(j)) =
          signs[j];
      }
    }
  }

  const bool exact = calibration == "exact";
  const double p_value = exact ?
    static_cast<double>(summary.exceedances) /
      static_cast<double>(summary.count) :
    (1.0 + static_cast<double>(summary.exceedances)) /
      (1.0 + static_cast<double>(summary.count));
  const double empirical_tail =
    static_cast<double>(summary.exceedances) /
    static_cast<double>(summary.count);
  const double mc_standard_error = exact ? NA_REAL :
    std::sqrt(empirical_tail * (1.0 - empirical_tail) /
      static_cast<double>(summary.count));
  const double total_sign_configurations = sign_dimension <= 1023U ?
    std::ldexp(1.0, static_cast<int>(sign_dimension)) : R_PosInf;
  const long double random_variance_scaled = summary.variance();

  arma::mat pseudo_x = wx_restore_pseudo_units(
    prepared.pseudo_x_scaled, prepared.global_scale
  );
  arma::mat pseudo_y = wx_restore_pseudo_units(
    prepared.pseudo_y_scaled, prepared.global_scale
  );
  Rcpp::IntegerVector discarded_x;
  Rcpp::IntegerVector discarded_y;
  if ((n1 % 2U) != 0U) {
    discarded_x = Rcpp::IntegerVector::create(static_cast<int>(n1));
  }
  if ((n2 % 2U) != 0U) {
    discarded_y = Rcpp::IntegerVector::create(static_cast<int>(n2));
  }
  Rcpp::RObject retained_physical_output = R_NilValue;
  Rcpp::RObject retained_scaled_output = R_NilValue;
  Rcpp::RObject retained_signs_output = R_NilValue;
  if (keep_randomized) {
    retained_physical_output = retained_physical;
    retained_scaled_output = retained_scaled;
    retained_signs_output = retained_signs;
  }

  return Rcpp::List::create(
    Rcpp::Named("observed") = observed,
    Rcpp::Named("observed_scaled") = wx_checked_double(
      prepared.observed_scaled, "the scaled observed statistic"
    ),
    Rcpp::Named("mean_difference_squared") = wx_report_second_order(
      prepared.mean_difference_squared_scaled, prepared.global_scale
    ),
    Rcpp::Named("mean_difference_squared_scaled") = wx_checked_double(
      prepared.mean_difference_squared_scaled,
      "the scaled squared sample-mean difference"
    ),
    Rcpp::Named("trace_S1") = wx_report_second_order(
      prepared.trace_s1_scaled, prepared.global_scale
    ),
    Rcpp::Named("trace_S2") = wx_report_second_order(
      prepared.trace_s2_scaled, prepared.global_scale
    ),
    Rcpp::Named("trace_S1_scaled") = wx_checked_double(
      prepared.trace_s1_scaled, "the scaled first covariance trace"
    ),
    Rcpp::Named("trace_S2_scaled") = wx_checked_double(
      prepared.trace_s2_scaled, "the scaled second covariance trace"
    ),
    Rcpp::Named("mean_x") = prepared.mean_x,
    Rcpp::Named("mean_y") = prepared.mean_y,
    Rcpp::Named("difference") = prepared.difference,
    Rcpp::Named("pseudo_x") = pseudo_x,
    Rcpp::Named("pseudo_y") = pseudo_y,
    Rcpp::Named("pseudo_x_scaled") = prepared.pseudo_x_scaled,
    Rcpp::Named("pseudo_y_scaled") = prepared.pseudo_y_scaled,
    Rcpp::Named("pair_indices_x") = wx_pair_indices(m1),
    Rcpp::Named("pair_indices_y") = wx_pair_indices(m2),
    Rcpp::Named("discarded_rows_x") = discarded_x,
    Rcpp::Named("discarded_rows_y") = discarded_y,
    Rcpp::Named("quadratic_kernel_scaled") = kernel,
    Rcpp::Named("p_value") = p_value,
    Rcpp::Named("empirical_tail") = empirical_tail,
    Rcpp::Named("exceedances") = static_cast<double>(summary.exceedances),
    Rcpp::Named("reference_evaluations") = static_cast<double>(summary.count),
    Rcpp::Named("total_sign_configurations") = total_sign_configurations,
    Rcpp::Named("log2_total_sign_configurations") =
      static_cast<double>(sign_dimension),
    Rcpp::Named("randomized_mean") = wx_report_second_order(
      summary.mean, prepared.global_scale
    ),
    Rcpp::Named("randomized_variance") = wx_report_fourth_order(
      random_variance_scaled, prepared.global_scale
    ),
    Rcpp::Named("randomized_minimum") = wx_report_second_order(
      summary.minimum, prepared.global_scale
    ),
    Rcpp::Named("randomized_maximum") = wx_report_second_order(
      summary.maximum, prepared.global_scale
    ),
    Rcpp::Named("randomized_mean_scaled") = wx_checked_double(
      summary.mean, "the scaled randomized mean"
    ),
    Rcpp::Named("randomized_variance_scaled") = wx_checked_double(
      random_variance_scaled, "the scaled randomized variance"
    ),
    Rcpp::Named("randomized_minimum_scaled") = wx_checked_double(
      summary.minimum, "the scaled randomized minimum"
    ),
    Rcpp::Named("randomized_maximum_scaled") = wx_checked_double(
      summary.maximum, "the scaled randomized maximum"
    ),
    Rcpp::Named("randomized_statistics") = retained_physical_output,
    Rcpp::Named("randomized_statistics_scaled") = retained_scaled_output,
    Rcpp::Named("sign_patterns") = retained_signs_output,
    Rcpp::Named("calibration") = calibration,
    Rcpp::Named("exact_reference_enumerated") = exact,
    Rcpp::Named("global_sign_symmetry_reduced") = exact,
    Rcpp::Named("plus_one_correction") = !exact,
    Rcpp::Named("mc_standard_error") = mc_standard_error,
    Rcpp::Named("minimum_attainable_p") = exact ?
      0.0 : 1.0 / (static_cast<double>(B) + 1.0),
    Rcpp::Named("seed") = exact ? NA_REAL : seed,
    Rcpp::Named("B") = exact ? NA_REAL : static_cast<double>(B),
    Rcpp::Named("n1") = static_cast<double>(n1),
    Rcpp::Named("n2") = static_cast<double>(n2),
    Rcpp::Named("m1") = static_cast<double>(m1),
    Rcpp::Named("m2") = static_cast<double>(m2),
    Rcpp::Named("p") = static_cast<double>(p),
    Rcpp::Named("sign_dimension") = static_cast<double>(sign_dimension),
    Rcpp::Named("global_scale") = prepared.global_scale,
    Rcpp::Named("subtraction_overflow_fallback") =
      prepared.subtraction_overflow_fallback,
    Rcpp::Named("observed_physical_underflow") = observed_underflow,
    Rcpp::Named("randomization_degenerate") =
      summary.minimum == summary.maximum
  );
}
