// Zhang--Zhu--Zhang normal-reference scale-invariant two-sample test.
//
// This is the feasible statistic and Welch--Satterthwaite calibration in
// Zhang, Zhu and Zhang (2023).  In particular, it is neither the raw-L2
// normal-reference statistic nor the later normal-reference F-type test.

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp17)]]

#include <RcppArmadillo.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <string>
#include <vector>

namespace {

void zzz23_neumaier_add(const long double value,
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

double zzz23_checked_double(const long double value,
                            const char* quantity) {
  const long double maximum = static_cast<long double>(
    std::numeric_limits<double>::max()
  );
  if (!std::isfinite(value) || value > maximum || value < -maximum) {
    Rcpp::stop(
      "Zhang-Zhu-Zhang produced a non-finite %s; no ridge, absolute-value "
      "repair, numerical floor, degree-of-freedom clamp, or pseudoinverse "
      "is applied.",
      quantity
    );
  }
  return static_cast<double>(value);
}

double zzz23_checked_positive_double(const long double value,
                                     const char* quantity) {
  const double answer = zzz23_checked_double(value, quantity);
  if (!(answer > 0.0)) {
    Rcpp::stop(
      "Zhang-Zhu-Zhang requires a finite, strictly positive %s; no "
      "ridge, absolute-value repair, numerical floor, degree-of-freedom "
      "clamp, or pseudoinverse is applied.",
      quantity
    );
  }
  return answer;
}

struct Zzz23PreparedData {
  arma::mat x;
  arma::mat y;
  arma::vec anchor;
  arma::vec scale_base;
  arma::vec scale_ratio;
  arma::vec log_scale;
  arma::uvec subtraction_overflow_fallback;
};

Zzz23PreparedData zzz23_prepare_data(const arma::mat& x,
                                     const arma::mat& y) {
  const arma::uword n1 = x.n_rows;
  const arma::uword n2 = y.n_rows;
  const arma::uword p = x.n_cols;
  const long double double_max = static_cast<long double>(
    std::numeric_limits<double>::max()
  );

  Zzz23PreparedData out;
  out.x.set_size(n1, p);
  out.y.set_size(n2, p);
  out.anchor.set_size(p);
  out.scale_base.ones(p);
  out.scale_ratio.ones(p);
  out.log_scale.zeros(p);
  out.subtraction_overflow_fallback.zeros(p);

  for (arma::uword j = 0; j < p; ++j) {
    const long double anchor = static_cast<long double>(x(0, j));
    out.anchor(j) = x(0, j);
    std::vector<long double> x_difference(n1);
    std::vector<long double> y_difference(n2);
    long double maximum_difference = 0.0L;
    bool direct_ok = true;

    for (arma::uword i = 0; i < n1; ++i) {
      const long double difference =
        static_cast<long double>(x(i, j)) - anchor;
      x_difference[i] = difference;
      direct_ok = direct_ok && std::isfinite(difference);
      maximum_difference = std::max(
        maximum_difference, std::abs(difference)
      );
    }
    for (arma::uword i = 0; i < n2; ++i) {
      const long double difference =
        static_cast<long double>(y(i, j)) - anchor;
      y_difference[i] = difference;
      direct_ok = direct_ok && std::isfinite(difference);
      maximum_difference = std::max(
        maximum_difference, std::abs(difference)
      );
    }
    direct_ok = direct_ok && std::isfinite(maximum_difference) &&
      maximum_difference > 0.0L && maximum_difference <= double_max;

    if (direct_ok) {
      out.scale_base(j) = static_cast<double>(maximum_difference);
      out.log_scale(j) = std::log(out.scale_base(j));
      for (arma::uword i = 0; i < n1; ++i) {
        out.x(i, j) = static_cast<double>(
          x_difference[i] / maximum_difference
        );
      }
      for (arma::uword i = 0; i < n2; ++i) {
        out.y(i, j) = static_cast<double>(
          y_difference[i] / maximum_difference
        );
      }
      continue;
    }

    // This route also covers platforms whose long double has the same
    // exponent range as double.  Dividing both operands by one common finite
    // scale before subtraction prevents overflow without changing a
    // coordinatewise scale-invariant statistic.
    out.subtraction_overflow_fallback(j) = 1u;
    double operand_scale = std::abs(x(0, j));
    for (arma::uword i = 0; i < n1; ++i) {
      operand_scale = std::max(operand_scale, std::abs(x(i, j)));
    }
    for (arma::uword i = 0; i < n2; ++i) {
      operand_scale = std::max(operand_scale, std::abs(y(i, j)));
    }
    if (!(operand_scale > 0.0) || !std::isfinite(operand_scale)) {
      Rcpp::stop(
        "Zhang-Zhu-Zhang requires pooled positive variation in variable "
        "%llu.",
        static_cast<unsigned long long>(j + 1)
      );
    }

    const long double operand_scale_ld =
      static_cast<long double>(operand_scale);
    const long double anchor_scaled = anchor / operand_scale_ld;
    maximum_difference = 0.0L;
    for (arma::uword i = 0; i < n1; ++i) {
      const long double difference =
        static_cast<long double>(x(i, j)) / operand_scale_ld -
        anchor_scaled;
      x_difference[i] = difference;
      maximum_difference = std::max(
        maximum_difference, std::abs(difference)
      );
    }
    for (arma::uword i = 0; i < n2; ++i) {
      const long double difference =
        static_cast<long double>(y(i, j)) / operand_scale_ld -
        anchor_scaled;
      y_difference[i] = difference;
      maximum_difference = std::max(
        maximum_difference, std::abs(difference)
      );
    }
    if (!(maximum_difference > 0.0L) ||
        !std::isfinite(maximum_difference) ||
        maximum_difference > double_max) {
      Rcpp::stop(
        "Zhang-Zhu-Zhang requires pooled positive variation in variable "
        "%llu.",
        static_cast<unsigned long long>(j + 1)
      );
    }

    out.scale_base(j) = operand_scale;
    out.scale_ratio(j) = static_cast<double>(maximum_difference);
    out.log_scale(j) = std::log(operand_scale) +
      std::log(out.scale_ratio(j));
    for (arma::uword i = 0; i < n1; ++i) {
      out.x(i, j) = static_cast<double>(
        x_difference[i] / maximum_difference
      );
    }
    for (arma::uword i = 0; i < n2; ++i) {
      out.y(i, j) = static_cast<double>(
        y_difference[i] / maximum_difference
      );
    }
  }

  if (!out.x.is_finite() || !out.y.is_finite()) {
    Rcpp::stop(
      "Zhang-Zhu-Zhang internal coordinate scaling produced non-finite "
      "data."
    );
  }
  return out;
}

struct Zzz23Moments {
  arma::vec mean;
  arma::vec variance;
  arma::mat centered;
};

Zzz23Moments zzz23_moments(const arma::mat& data) {
  const arma::uword n = data.n_rows;
  const arma::uword p = data.n_cols;
  const long double n_ld = static_cast<long double>(n);
  Zzz23Moments out;
  out.mean.set_size(p);
  out.variance.set_size(p);
  out.centered.set_size(n, p);

  for (arma::uword j = 0; j < p; ++j) {
    const long double anchor = static_cast<long double>(data(0, j));
    long double delta_total = 0.0L;
    long double delta_correction = 0.0L;
    for (arma::uword i = 0; i < n; ++i) {
      zzz23_neumaier_add(
        static_cast<long double>(data(i, j)) - anchor,
        delta_total, delta_correction
      );
    }
    const long double mean = anchor +
      (delta_total + delta_correction) / n_ld;
    out.mean(j) = zzz23_checked_double(mean, "scaled group mean");

    long double square_total = 0.0L;
    long double square_correction = 0.0L;
    for (arma::uword i = 0; i < n; ++i) {
      const long double residual =
        static_cast<long double>(data(i, j)) - mean;
      out.centered(i, j) = zzz23_checked_double(
        residual, "scaled centered observation"
      );
      zzz23_neumaier_add(
        residual * residual, square_total, square_correction
      );
    }
    const long double variance =
      (square_total + square_correction) /
      static_cast<long double>(n - 1u);
    if (!std::isfinite(variance) || variance < 0.0L) {
      Rcpp::stop(
        "Zhang-Zhu-Zhang encountered a non-finite or negative unbiased "
        "marginal sample variance; no numerical repair is applied."
      );
    }
    out.variance(j) = static_cast<double>(variance);
  }
  return out;
}

long double zzz23_dot_rows(const arma::mat& x,
                           const arma::uword i,
                           const arma::mat& y,
                           const arma::uword k) {
  long double total = 0.0L;
  long double correction = 0.0L;
  for (arma::uword j = 0; j < x.n_cols; ++j) {
    zzz23_neumaier_add(
      static_cast<long double>(x(i, j)) *
        static_cast<long double>(y(k, j)),
      total, correction
    );
  }
  return total + correction;
}

long double zzz23_dot_columns(const arma::mat& x,
                              const arma::uword j,
                              const arma::mat& y,
                              const arma::uword k) {
  long double total = 0.0L;
  long double correction = 0.0L;
  for (arma::uword i = 0; i < x.n_rows; ++i) {
    zzz23_neumaier_add(
      static_cast<long double>(x(i, j)) *
        static_cast<long double>(y(i, k)),
      total, correction
    );
  }
  return total + correction;
}

struct Zzz23TraceSet {
  long double trace1;
  long double trace2;
  long double square1;
  long double square2;
  long double cross;
  std::string strategy;
};

Zzz23TraceSet zzz23_trace_primal(const arma::mat& z1,
                                 const arma::mat& z2) {
  const arma::uword p = z1.n_cols;
  const long double denominator1 = static_cast<long double>(
    z1.n_rows - 1u
  );
  const long double denominator2 = static_cast<long double>(
    z2.n_rows - 1u
  );
  Zzz23TraceSet out{0.0L, 0.0L, 0.0L, 0.0L, 0.0L, "primal"};
  long double trace1_correction = 0.0L;
  long double trace2_correction = 0.0L;
  long double square1_correction = 0.0L;
  long double square2_correction = 0.0L;
  long double cross_correction = 0.0L;

  for (arma::uword j = 0; j < p; ++j) {
    for (arma::uword k = j; k < p; ++k) {
      const long double r1 =
        zzz23_dot_columns(z1, j, z1, k) / denominator1;
      const long double r2 =
        zzz23_dot_columns(z2, j, z2, k) / denominator2;
      const long double multiplicity = j == k ? 1.0L : 2.0L;
      zzz23_neumaier_add(
        multiplicity * r1 * r1,
        out.square1, square1_correction
      );
      zzz23_neumaier_add(
        multiplicity * r2 * r2,
        out.square2, square2_correction
      );
      zzz23_neumaier_add(
        multiplicity * r1 * r2,
        out.cross, cross_correction
      );
      if (j == k) {
        zzz23_neumaier_add(r1, out.trace1, trace1_correction);
        zzz23_neumaier_add(r2, out.trace2, trace2_correction);
      }
    }
  }
  out.trace1 += trace1_correction;
  out.trace2 += trace2_correction;
  out.square1 += square1_correction;
  out.square2 += square2_correction;
  out.cross += cross_correction;
  return out;
}

Zzz23TraceSet zzz23_trace_dual(const arma::mat& z1,
                               const arma::mat& z2) {
  const arma::uword n1 = z1.n_rows;
  const arma::uword n2 = z2.n_rows;
  const long double denominator1 = static_cast<long double>(n1 - 1u);
  const long double denominator2 = static_cast<long double>(n2 - 1u);
  Zzz23TraceSet out{0.0L, 0.0L, 0.0L, 0.0L, 0.0L, "dual"};
  long double trace1_correction = 0.0L;
  long double trace2_correction = 0.0L;
  long double square1_correction = 0.0L;
  long double square2_correction = 0.0L;
  long double cross_correction = 0.0L;

  for (arma::uword i = 0; i < n1; ++i) {
    for (arma::uword k = i; k < n1; ++k) {
      const long double inner = zzz23_dot_rows(z1, i, z1, k);
      const long double multiplicity = i == k ? 1.0L : 2.0L;
      zzz23_neumaier_add(
        multiplicity * inner * inner,
        out.square1, square1_correction
      );
    }
    const long double diagonal = zzz23_dot_rows(z1, i, z1, i);
    zzz23_neumaier_add(
      diagonal / denominator1,
      out.trace1, trace1_correction
    );
  }
  for (arma::uword i = 0; i < n2; ++i) {
    for (arma::uword k = i; k < n2; ++k) {
      const long double inner = zzz23_dot_rows(z2, i, z2, k);
      const long double multiplicity = i == k ? 1.0L : 2.0L;
      zzz23_neumaier_add(
        multiplicity * inner * inner,
        out.square2, square2_correction
      );
    }
    const long double diagonal = zzz23_dot_rows(z2, i, z2, i);
    zzz23_neumaier_add(
      diagonal / denominator2,
      out.trace2, trace2_correction
    );
  }
  for (arma::uword i = 0; i < n1; ++i) {
    for (arma::uword k = 0; k < n2; ++k) {
      const long double inner = zzz23_dot_rows(z1, i, z2, k);
      zzz23_neumaier_add(
        inner * inner, out.cross, cross_correction
      );
    }
  }

  out.trace1 += trace1_correction;
  out.trace2 += trace2_correction;
  out.square1 = (out.square1 + square1_correction) /
    (denominator1 * denominator1);
  out.square2 = (out.square2 + square2_correction) /
    (denominator2 * denominator2);
  out.cross = (out.cross + cross_correction) /
    (denominator1 * denominator2);
  return out;
}

arma::vec zzz23_input_mean(const arma::mat& data) {
  arma::vec answer(data.n_cols, arma::fill::zeros);
  const long double n = static_cast<long double>(data.n_rows);
  for (arma::uword j = 0; j < data.n_cols; ++j) {
    double scale = 0.0;
    for (arma::uword i = 0; i < data.n_rows; ++i) {
      scale = std::max(scale, std::abs(data(i, j)));
    }
    if (scale == 0.0) {
      answer(j) = 0.0;
      continue;
    }
    long double total = 0.0L;
    long double correction = 0.0L;
    for (arma::uword i = 0; i < data.n_rows; ++i) {
      zzz23_neumaier_add(
        static_cast<long double>(data(i, j)) /
          static_cast<long double>(scale),
        total, correction
      );
    }
    answer(j) = zzz23_checked_double(
      static_cast<long double>(scale) *
        (total + correction) / n,
      "input-scale sample mean"
    );
  }
  return answer;
}

arma::vec zzz23_report_difference(
    const arma::vec& standardized_difference,
    const Zzz23PreparedData& prepared,
    arma::uvec& representable) {
  arma::vec answer(standardized_difference.n_elem, arma::fill::zeros);
  representable.ones(standardized_difference.n_elem);
  const double maximum_log = std::log(
    std::numeric_limits<double>::max()
  );
  for (arma::uword j = 0; j < answer.n_elem; ++j) {
    const double value = standardized_difference(j);
    if (value == 0.0) {
      continue;
    }
    const double log_absolute = prepared.log_scale(j) +
      std::log(std::abs(value));
    if (log_absolute > maximum_log) {
      answer(j) = value > 0.0 ? R_PosInf : R_NegInf;
      representable(j) = 0u;
    } else {
      answer(j) = std::copysign(std::exp(log_absolute), value);
    }
  }
  return answer;
}

arma::vec zzz23_report_scale(const arma::vec& log_value) {
  arma::vec answer(log_value.n_elem);
  const double maximum_log = std::log(
    std::numeric_limits<double>::max()
  );
  const double minimum_log = std::log(
    std::numeric_limits<double>::denorm_min()
  );
  for (arma::uword j = 0; j < log_value.n_elem; ++j) {
    if (log_value(j) > maximum_log) {
      answer(j) = R_PosInf;
    } else if (log_value(j) < minimum_log) {
      answer(j) = 0.0;
    } else {
      answer(j) = std::exp(log_value(j));
    }
  }
  return answer;
}

} // namespace


// [[Rcpp::export]]
Rcpp::List cpp_zhang_zhu_zhang_two_sample(const arma::mat& x,
                                          const arma::mat& y) {
  if (x.n_rows < 3u || y.n_rows < 3u) {
    Rcpp::stop(
      "Zhang-Zhu-Zhang requires at least three observations in each "
      "sample because its trace correction contains n_i - 2."
    );
  }
  if (x.n_cols < 1u || x.n_cols != y.n_cols) {
    Rcpp::stop("`x` and `y` must have the same positive number of columns.");
  }
  if (!x.is_finite() || !y.is_finite()) {
    Rcpp::stop("`x` and `y` must contain only finite values.");
  }

  const arma::uword n1 = x.n_rows;
  const arma::uword n2 = y.n_rows;
  const arma::uword p = x.n_cols;
  const long double n1_ld = static_cast<long double>(n1);
  const long double n2_ld = static_cast<long double>(n2);
  const long double n_ld = n1_ld + n2_ld;
  const long double p_ld = static_cast<long double>(p);
  const long double weight1 = n2_ld / n_ld;
  const long double weight2 = n1_ld / n_ld;

  const Zzz23PreparedData prepared = zzz23_prepare_data(x, y);
  const Zzz23Moments moments1 = zzz23_moments(prepared.x);
  const Zzz23Moments moments2 = zzz23_moments(prepared.y);

  arma::vec pooled_diagonal(p);
  arma::vec standardized_difference(p);
  arma::vec coordinate_contribution(p);
  arma::mat z1(n1, p);
  arma::mat z2(n2, p);
  arma::uvec zero_variance1(p, arma::fill::zeros);
  arma::uvec zero_variance2(p, arma::fill::zeros);
  long double statistic_total = 0.0L;
  long double statistic_correction = 0.0L;
  long double maximum_identity_error = 0.0L;
  long double minimum_diagonal =
    std::numeric_limits<long double>::infinity();
  const long double statistic_factor = n1_ld * n2_ld /
    (n_ld * p_ld);

  for (arma::uword j = 0; j < p; ++j) {
    const long double variance1 = moments1.variance(j);
    const long double variance2 = moments2.variance(j);
    zero_variance1(j) = variance1 == 0.0L ? 1u : 0u;
    zero_variance2(j) = variance2 == 0.0L ? 1u : 0u;
    const long double diagonal =
      weight1 * variance1 + weight2 * variance2;
    if (!std::isfinite(diagonal) || diagonal <= 0.0L) {
      Rcpp::stop(
        "Zhang-Zhu-Zhang requires the crossed pooled marginal variance "
        "to be finite and strictly positive in variable %llu; no ridge "
        "or variance floor is applied.",
        static_cast<unsigned long long>(j + 1)
      );
    }
    pooled_diagonal(j) = zzz23_checked_positive_double(
      diagonal, "crossed pooled marginal variance"
    );
    minimum_diagonal = std::min(minimum_diagonal, diagonal);
    const long double identity =
      weight1 * variance1 / diagonal +
      weight2 * variance2 / diagonal;
    maximum_identity_error = std::max(
      maximum_identity_error, std::abs(identity - 1.0L)
    );

    const long double difference =
      static_cast<long double>(moments1.mean(j)) -
      static_cast<long double>(moments2.mean(j));
    standardized_difference(j) = zzz23_checked_double(
      difference, "internally scaled sample-mean difference"
    );
    const long double contribution = statistic_factor *
      difference * difference / diagonal;
    coordinate_contribution(j) = zzz23_checked_double(
      contribution, "coordinate statistic contribution"
    );
    zzz23_neumaier_add(
      contribution, statistic_total, statistic_correction
    );

    const double inverse_root = 1.0 / std::sqrt(pooled_diagonal(j));
    for (arma::uword i = 0; i < n1; ++i) {
      z1(i, j) = moments1.centered(i, j) * inverse_root;
    }
    for (arma::uword i = 0; i < n2; ++i) {
      z2(i, j) = moments2.centered(i, j) * inverse_root;
    }
  }
  if (!z1.is_finite() || !z2.is_finite()) {
    Rcpp::stop(
      "Zhang-Zhu-Zhang diagonal standardisation produced non-finite "
      "residuals; no numerical floor is applied."
    );
  }

  const long double statistic = statistic_total + statistic_correction;
  if (statistic < 0.0L) {
    Rcpp::stop(
      "Zhang-Zhu-Zhang produced a negative quadratic statistic; no "
      "absolute-value repair is applied."
    );
  }

  Zzz23TraceSet traces;
  if (p <= n1 + n2) {
    traces = zzz23_trace_primal(z1, z2);
  } else {
    traces = zzz23_trace_dual(z1, z2);
  }
  if (!std::isfinite(traces.trace1) ||
      !std::isfinite(traces.trace2) ||
      !std::isfinite(traces.square1) ||
      !std::isfinite(traces.square2) ||
      !std::isfinite(traces.cross) ||
      traces.trace1 < 0.0L || traces.trace2 < 0.0L ||
      traces.square1 < 0.0L || traces.square2 < 0.0L ||
      traces.cross < 0.0L) {
    Rcpp::stop(
      "Zhang-Zhu-Zhang produced an invalid raw covariance trace; no "
      "absolute-value repair or floor is applied."
    );
  }

  const long double factor1 =
    (n1_ld - 1.0L) * (n1_ld - 1.0L) /
    ((n1_ld - 2.0L) * (n1_ld + 1.0L));
  const long double factor2 =
    (n2_ld - 1.0L) * (n2_ld - 1.0L) /
    ((n2_ld - 2.0L) * (n2_ld + 1.0L));
  const long double bracket1 = traces.square1 -
    traces.trace1 * traces.trace1 / (n1_ld - 1.0L);
  const long double bracket2 = traces.square2 -
    traces.trace2 * traces.trace2 / (n2_ld - 1.0L);
  const long double corrected1 = factor1 * bracket1;
  const long double corrected2 = factor2 * bracket2;
  const long double trace_corrected =
    weight1 * weight1 * corrected1 +
    weight2 * weight2 * corrected2 +
    2.0L * weight1 * weight2 * traces.cross;
  if (!std::isfinite(trace_corrected) || trace_corrected <= 0.0L) {
    Rcpp::stop(
      "Zhang-Zhu-Zhang requires its bias-corrected trace estimate to be "
      "finite and strictly positive; no absolute-value repair, floor, or "
      "degree-of-freedom substitution is applied."
    );
  }

  const long double trace_raw =
    weight1 * weight1 * traces.square1 +
    weight2 * weight2 * traces.square2 +
    2.0L * weight1 * weight2 * traces.cross;
  if (!std::isfinite(trace_raw) || trace_raw <= 0.0L) {
    Rcpp::stop(
      "Zhang-Zhu-Zhang requires its raw pooled squared trace to be finite "
      "and strictly positive."
    );
  }
  const long double trace_rn =
    weight1 * traces.trace1 + weight2 * traces.trace2;
  const long double df_unadjusted =
    p_ld * p_ld / trace_corrected;
  const long double c_np = 1.0L + trace_raw /
    (p_ld * std::sqrt(p_ld));
  const bool paper_correction_applied = c_np <= 1.2L;
  const long double df_paper = paper_correction_applied ?
    df_unadjusted / c_np : df_unadjusted;

  arma::vec log_diagonal_input(p);
  for (arma::uword j = 0; j < p; ++j) {
    log_diagonal_input(j) = 2.0 * prepared.log_scale(j) +
      std::log(pooled_diagonal(j));
  }
  const double maximum_log_diagonal = log_diagonal_input.max();
  arma::vec diagonal_input_canonical =
    arma::exp(log_diagonal_input - maximum_log_diagonal);
  const arma::vec diagonal_input = zzz23_report_scale(
    log_diagonal_input
  );
  const arma::vec column_scale = zzz23_report_scale(
    prepared.log_scale
  );
  arma::uvec difference_representable;
  const arma::vec difference_input = zzz23_report_difference(
    standardized_difference, prepared, difference_representable
  );
  const arma::vec mean1_input = zzz23_input_mean(x);
  const arma::vec mean2_input = zzz23_input_mean(y);

  return Rcpp::List::create(
    Rcpp::Named("statistic") = zzz23_checked_double(
      statistic, "test statistic"
    ),
    Rcpp::Named("coordinate_contribution") = coordinate_contribution,
    Rcpp::Named("mean1") = mean1_input,
    Rcpp::Named("mean2") = mean2_input,
    Rcpp::Named("mean_difference_input") = difference_input,
    Rcpp::Named("mean_difference_input_representable") =
      difference_representable,
    Rcpp::Named("mean1_standardized") = moments1.mean,
    Rcpp::Named("mean2_standardized") = moments2.mean,
    Rcpp::Named("mean_difference_standardized") =
      standardized_difference,
    Rcpp::Named("variance1_diagonal_standardized") =
      moments1.variance,
    Rcpp::Named("variance2_diagonal_standardized") =
      moments2.variance,
    Rcpp::Named("D_hat_diagonal_standardized") = pooled_diagonal,
    Rcpp::Named("D_hat_diagonal_input") = diagonal_input,
    Rcpp::Named("D_hat_diagonal_input_canonical") =
      diagonal_input_canonical,
    Rcpp::Named("log_D_hat_diagonal_input") = log_diagonal_input,
    Rcpp::Named("trace_R1") = zzz23_checked_double(
      traces.trace1, "trace(R1 hat)"
    ),
    Rcpp::Named("trace_R2") = zzz23_checked_double(
      traces.trace2, "trace(R2 hat)"
    ),
    Rcpp::Named("trace_Rn") = zzz23_checked_double(
      trace_rn, "trace(Rn hat)"
    ),
    Rcpp::Named("trace_R1_squared_raw") = zzz23_checked_double(
      traces.square1, "raw trace(R1 hat squared)"
    ),
    Rcpp::Named("trace_R2_squared_raw") = zzz23_checked_double(
      traces.square2, "raw trace(R2 hat squared)"
    ),
    Rcpp::Named("trace_R1_R2_raw") = zzz23_checked_double(
      traces.cross, "raw trace(R1 hat R2 hat)"
    ),
    Rcpp::Named("trace_Rn_squared_raw") = zzz23_checked_positive_double(
      trace_raw, "raw trace(Rn hat squared)"
    ),
    Rcpp::Named("trace_R1_squared_bracket") = zzz23_checked_double(
      bracket1, "group-1 corrected-trace bracket"
    ),
    Rcpp::Named("trace_R2_squared_bracket") = zzz23_checked_double(
      bracket2, "group-2 corrected-trace bracket"
    ),
    Rcpp::Named("trace_R1_squared_corrected") = zzz23_checked_double(
      corrected1, "group-1 corrected squared trace"
    ),
    Rcpp::Named("trace_R2_squared_corrected") = zzz23_checked_double(
      corrected2, "group-2 corrected squared trace"
    ),
    Rcpp::Named("trace_Rn_squared_corrected") =
      zzz23_checked_positive_double(
        trace_corrected, "bias-corrected trace estimate"
      ),
    Rcpp::Named("trace_correction_factor1") = static_cast<double>(factor1),
    Rcpp::Named("trace_correction_factor2") = static_cast<double>(factor2),
    Rcpp::Named("df_unadjusted") = zzz23_checked_positive_double(
      df_unadjusted, "unadjusted reference degrees of freedom"
    ),
    Rcpp::Named("c_np") = zzz23_checked_positive_double(
      c_np, "empirical finite-sample factor"
    ),
    Rcpp::Named("df_paper") = zzz23_checked_positive_double(
      df_paper, "paper-adjusted reference degrees of freedom"
    ),
    Rcpp::Named("paper_correction_applied") =
      paper_correction_applied,
    Rcpp::Named("weight_group1_covariance") =
      static_cast<double>(weight1),
    Rcpp::Named("weight_group2_covariance") =
      static_cast<double>(weight2),
    Rcpp::Named("column_anchor") = prepared.anchor,
    Rcpp::Named("column_scale") = column_scale,
    Rcpp::Named("column_log_scale") = prepared.log_scale,
    Rcpp::Named("subtraction_overflow_fallback") =
      prepared.subtraction_overflow_fallback,
    Rcpp::Named("zero_variance_group1") = zero_variance1,
    Rcpp::Named("zero_variance_group2") = zero_variance2,
    Rcpp::Named("minimum_D_hat_standardized") =
      zzz23_checked_positive_double(
        minimum_diagonal, "minimum crossed pooled marginal variance"
      ),
    Rcpp::Named("maximum_diagonal_identity_error") =
      zzz23_checked_double(
        maximum_identity_error, "diagonal standardisation identity error"
      ),
    Rcpp::Named("trace_computation") = traces.strategy,
    Rcpp::Named("constructs_p_by_p_matrix") = false,
    Rcpp::Named("n1") = static_cast<double>(n1),
    Rcpp::Named("n2") = static_cast<double>(n2),
    Rcpp::Named("p") = static_cast<double>(p)
  );
}
