// Feng--Zou--Wang--Zhu (2015) scale-invariant Behrens--Fisher test.
//
// The variance estimator below is the estimator in the original paper: the
// two within-group traces use leave-four-out marginal variances and the cross
// trace uses leave-two-out marginal variances in both samples.  For a fixed
// unordered set of four observations, all 24 ordered kernels share the same
// diagonal denominator.  We therefore form one local 4 by 4 dual Gram matrix
// and sum the 24 kernels from it.  This is an exact symmetry reduction, not a
// replacement of the published estimator.

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp17)]]

#include <RcppArmadillo.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <string>
#include <vector>

namespace {

struct FzwzData {
  arma::uword n1;
  arma::uword n2;
  arma::uword p;
  std::vector<long double> x;
  std::vector<long double> y;
  std::vector<long double> origin;
  std::vector<long double> scale;
};

struct FzwzMoments {
  arma::uword n;
  arma::uword p;
  const std::vector<long double>* values;
  std::vector<long double> anchor;
  std::vector<long double> delta_sum;
  std::vector<long double> delta_sumsq;
  std::vector<long double> mean;
  std::vector<long double> variance;
  std::vector<long double> third;
};

struct FzwzWithinTrace {
  long double estimate;
  long double minimum_denominator;
  long double combinations;
};

struct FzwzCrossTrace {
  long double estimate;
  long double minimum_denominator;
  long double pair_combinations;
};

inline std::size_t fzwz_index(arma::uword i,
                              arma::uword j,
                              arma::uword p) {
  return static_cast<std::size_t>(i) * static_cast<std::size_t>(p) +
    static_cast<std::size_t>(j);
}

void fzwz_require_finite_matrix(const arma::mat& x, const char* name) {
  if (!x.is_finite()) {
    Rcpp::stop("`%s` must contain only finite values.", name);
  }
}

double fzwz_checked_double(long double value, const char* quantity) {
  const double result = static_cast<double>(value);
  if (!std::isfinite(result)) {
    Rcpp::stop(
      "Feng-Zou-Wang-Zhu requires a finite %s; no ridge, absolute-value "
      "repair, or numerical floor is applied.",
      quantity
    );
  }
  return result;
}

double fzwz_checked_positive_double(long double value,
                                    const char* quantity) {
  const double result = static_cast<double>(value);
  if (!std::isfinite(result) || result <= 0.0) {
    Rcpp::stop(
      "Feng-Zou-Wang-Zhu requires a finite, strictly positive and "
      "double-representable %s; no ridge, absolute-value repair, or "
      "numerical floor is applied.",
      quantity
    );
  }
  return result;
}

double fzwz_report_long_double(long double value) {
  if (value > static_cast<long double>(
        std::numeric_limits<double>::max())) {
    return std::numeric_limits<double>::infinity();
  }
  if (value < -static_cast<long double>(
        std::numeric_limits<double>::max())) {
    return -std::numeric_limits<double>::infinity();
  }
  return static_cast<double>(value);
}

long double fzwz_choose4(arma::uword n) {
  const long double n_ld = static_cast<long double>(n);
  return n_ld * (n_ld - 1.0L) * (n_ld - 2.0L) * (n_ld - 3.0L) /
    24.0L;
}

long double fzwz_choose2(arma::uword n) {
  const long double n_ld = static_cast<long double>(n);
  return n_ld * (n_ld - 1.0L) / 2.0L;
}

FzwzData fzwz_standardize_columns(const arma::mat& x,
                                  const arma::mat& y) {
  FzwzData out;
  out.n1 = x.n_rows;
  out.n2 = y.n_rows;
  out.p = x.n_cols;
  out.x.resize(static_cast<std::size_t>(out.n1) * out.p);
  out.y.resize(static_cast<std::size_t>(out.n2) * out.p);
  out.origin.resize(out.p);
  out.scale.resize(out.p);

  for (arma::uword j = 0; j < out.p; ++j) {
    const long double anchor = static_cast<long double>(x(0, j));
    long double maximum_difference = 0.0L;
    bool direct_finite = true;
    for (arma::uword i = 0; i < out.n1; ++i) {
      const long double difference =
        static_cast<long double>(x(i, j)) - anchor;
      direct_finite = direct_finite && std::isfinite(difference);
      maximum_difference = std::max(
        maximum_difference, std::abs(difference)
      );
    }
    for (arma::uword i = 0; i < out.n2; ++i) {
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
      for (arma::uword i = 0; i < out.n1; ++i) {
        out.x[fzwz_index(i, j, out.p)] =
          (static_cast<long double>(x(i, j)) - anchor) / scale;
      }
      for (arma::uword i = 0; i < out.n2; ++i) {
        out.y[fzwz_index(i, j, out.p)] =
          (static_cast<long double>(y(i, j)) - anchor) / scale;
      }
      continue;
    }

    // This fallback is relevant on platforms where long double has no wider
    // exponent range than double and subtraction of opposite finite extremes
    // overflows.  Scaling operands first retains the direction exactly enough
    // for this coordinatewise scale-invariant statistic.
    double operand_scale = std::abs(x(0, j));
    for (arma::uword i = 0; i < out.n1; ++i) {
      operand_scale = std::max(operand_scale, std::abs(x(i, j)));
    }
    for (arma::uword i = 0; i < out.n2; ++i) {
      operand_scale = std::max(operand_scale, std::abs(y(i, j)));
    }
    if (operand_scale == 0.0) {
      operand_scale = 1.0;
    }
    const long double operand_scale_ld =
      static_cast<long double>(operand_scale);
    const long double anchor_scaled = anchor / operand_scale_ld;
    maximum_difference = 0.0L;
    for (arma::uword i = 0; i < out.n1; ++i) {
      const long double difference =
        static_cast<long double>(x(i, j)) / operand_scale_ld -
        anchor_scaled;
      out.x[fzwz_index(i, j, out.p)] = difference;
      maximum_difference = std::max(
        maximum_difference, std::abs(difference)
      );
    }
    for (arma::uword i = 0; i < out.n2; ++i) {
      const long double difference =
        static_cast<long double>(y(i, j)) / operand_scale_ld -
        anchor_scaled;
      out.y[fzwz_index(i, j, out.p)] = difference;
      maximum_difference = std::max(
        maximum_difference, std::abs(difference)
      );
    }
    const long double second_scale = maximum_difference > 0.0L ?
      maximum_difference : 1.0L;
    out.scale[j] = operand_scale_ld * second_scale;
    for (arma::uword i = 0; i < out.n1; ++i) {
      out.x[fzwz_index(i, j, out.p)] /= second_scale;
    }
    for (arma::uword i = 0; i < out.n2; ++i) {
      out.y[fzwz_index(i, j, out.p)] /= second_scale;
    }
  }
  return out;
}

FzwzMoments fzwz_moments(const std::vector<long double>& values,
                         arma::uword n,
                         arma::uword p) {
  FzwzMoments out;
  out.n = n;
  out.p = p;
  out.values = &values;
  out.anchor.resize(p);
  out.delta_sum.assign(p, 0.0L);
  out.delta_sumsq.assign(p, 0.0L);
  out.mean.resize(p);
  out.variance.resize(p);
  out.third.resize(p);
  const long double n_ld = static_cast<long double>(n);

  for (arma::uword j = 0; j < p; ++j) {
    const long double anchor = values[fzwz_index(0, j, p)];
    out.anchor[j] = anchor;
    for (arma::uword i = 0; i < n; ++i) {
      const long double delta =
        values[fzwz_index(i, j, p)] - anchor;
      out.delta_sum[j] += delta;
      out.delta_sumsq[j] += delta * delta;
    }
    const long double mean = anchor + out.delta_sum[j] / n_ld;
    out.mean[j] = mean;
    long double squared_sum = 0.0L;
    long double cubed_sum = 0.0L;
    for (arma::uword i = 0; i < n; ++i) {
      const long double centered =
        values[fzwz_index(i, j, p)] - mean;
      squared_sum += centered * centered;
      cubed_sum += centered * centered * centered;
    }
    out.variance[j] = squared_sum / (n_ld - 1.0L);
    out.third[j] = cubed_sum / n_ld;
  }
  return out;
}

long double fzwz_leaveout_variance(
    const FzwzMoments& moments,
    arma::uword j,
    const arma::uword* omitted,
    arma::uword omitted_count) {
  const arma::uword remaining = moments.n - omitted_count;
  const long double remaining_ld = static_cast<long double>(remaining);
  long double sum = moments.delta_sum[j];
  long double sumsq = moments.delta_sumsq[j];
  for (arma::uword q = 0; q < omitted_count; ++q) {
    const long double delta =
      (*moments.values)[fzwz_index(omitted[q], j, moments.p)] -
      moments.anchor[j];
    sum -= delta;
    sumsq -= delta * delta;
  }
  const long double centered_sum = sumsq - sum * sum / remaining_ld;
  const long double variance = centered_sum / (remaining_ld - 1.0L);
  if (!std::isfinite(variance) || variance < 0.0L) {
    Rcpp::stop(
      "Feng-Zou-Wang-Zhu encountered a non-finite or negative "
      "leave-out marginal variance; no numerical floor is applied."
    );
  }
  return variance;
}

inline long double fzwz_gram_inner(const long double gram[4][4],
                                   int a,
                                   int b,
                                   int c,
                                   int d) {
  return gram[a][c] - gram[a][d] - gram[b][c] + gram[b][d];
}

FzwzWithinTrace fzwz_within_trace(
    const FzwzMoments& target,
    const FzwzMoments& other,
    bool target_is_group1,
    long double gamma) {
  const arma::uword n = target.n;
  const arma::uword p = target.p;
  long double kernel_sum = 0.0L;
  long double minimum_denominator =
    std::numeric_limits<long double>::infinity();
  const std::array<int, 4> initial_order = {0, 1, 2, 3};

  for (arma::uword i0 = 0; i0 + 3 < n; ++i0) {
    for (arma::uword i1 = i0 + 1; i1 + 2 < n; ++i1) {
      for (arma::uword i2 = i1 + 1; i2 + 1 < n; ++i2) {
        for (arma::uword i3 = i2 + 1; i3 < n; ++i3) {
          const arma::uword rows[4] = {i0, i1, i2, i3};
          long double gram[4][4] = {};
          for (arma::uword j = 0; j < p; ++j) {
            const long double leave_variance = fzwz_leaveout_variance(
              target, j, rows, 4
            );
            const long double denominator = target_is_group1 ?
              leave_variance + gamma * other.variance[j] :
              other.variance[j] + gamma * leave_variance;
            if (!std::isfinite(denominator) || denominator <= 0.0L) {
              Rcpp::stop(
                "Feng-Zou-Wang-Zhu requires every leave-four-out "
                "combined marginal variance to be finite and strictly "
                "positive; no ridge is applied."
              );
            }
            minimum_denominator = std::min(
              minimum_denominator, denominator
            );
            long double relative[4];
            relative[0] = 0.0L;
            const long double base =
              (*target.values)[fzwz_index(rows[0], j, p)];
            for (int a = 1; a < 4; ++a) {
              relative[a] =
                (*target.values)[fzwz_index(rows[a], j, p)] - base;
            }
            for (int a = 1; a < 4; ++a) {
              for (int b = a; b < 4; ++b) {
                gram[a][b] += relative[a] * relative[b] / denominator;
              }
            }
          }
          for (int a = 0; a < 4; ++a) {
            for (int b = a + 1; b < 4; ++b) {
              gram[b][a] = gram[a][b];
            }
          }

          std::array<int, 4> order = initial_order;
          do {
            const long double first = fzwz_gram_inner(
              gram, order[0], order[1], order[2], order[3]
            );
            const long double second = fzwz_gram_inner(
              gram, order[2], order[1], order[0], order[3]
            );
            kernel_sum += first * second;
          } while (std::next_permutation(order.begin(), order.end()));
        }
      }
    }
  }

  const long double combinations = fzwz_choose4(n);
  const long double ordered_quadruples = 24.0L * combinations;
  const long double estimate = kernel_sum /
    (2.0L * ordered_quadruples);
  if (!std::isfinite(estimate)) {
    Rcpp::stop(
      "Feng-Zou-Wang-Zhu produced a non-finite leave-four-out trace "
      "estimate."
    );
  }
  return {estimate, minimum_denominator, combinations};
}

FzwzCrossTrace fzwz_cross_trace(const FzwzMoments& group1,
                                const FzwzMoments& group2,
                                long double gamma) {
  const arma::uword p = group1.p;
  long double square_sum = 0.0L;
  long double minimum_denominator =
    std::numeric_limits<long double>::infinity();

  for (arma::uword i0 = 0; i0 + 1 < group1.n; ++i0) {
    for (arma::uword i1 = i0 + 1; i1 < group1.n; ++i1) {
      const arma::uword omitted1[2] = {i0, i1};
      for (arma::uword j0 = 0; j0 + 1 < group2.n; ++j0) {
        for (arma::uword j1 = j0 + 1; j1 < group2.n; ++j1) {
          const arma::uword omitted2[2] = {j0, j1};
          long double inner = 0.0L;
          for (arma::uword k = 0; k < p; ++k) {
            const long double variance1 = fzwz_leaveout_variance(
              group1, k, omitted1, 2
            );
            const long double variance2 = fzwz_leaveout_variance(
              group2, k, omitted2, 2
            );
            const long double denominator = variance1 + gamma * variance2;
            if (!std::isfinite(denominator) || denominator <= 0.0L) {
              Rcpp::stop(
                "Feng-Zou-Wang-Zhu requires every two-plus-two "
                "leave-out combined marginal variance to be finite and "
                "strictly positive; no ridge is applied."
              );
            }
            minimum_denominator = std::min(
              minimum_denominator, denominator
            );
            const long double difference1 =
              (*group1.values)[fzwz_index(i0, k, p)] -
              (*group1.values)[fzwz_index(i1, k, p)];
            const long double difference2 =
              (*group2.values)[fzwz_index(j0, k, p)] -
              (*group2.values)[fzwz_index(j1, k, p)];
            inner += difference1 * difference2 / denominator;
          }
          square_sum += inner * inner;
        }
      }
    }
  }

  const long double pair_combinations =
    fzwz_choose2(group1.n) * fzwz_choose2(group2.n);
  const long double estimate = square_sum /
    (4.0L * pair_combinations);
  if (!std::isfinite(estimate)) {
    Rcpp::stop(
      "Feng-Zou-Wang-Zhu produced a non-finite two-plus-two leave-out "
      "cross-trace estimate."
    );
  }
  return {estimate, minimum_denominator, pair_combinations};
}

}  // namespace


//' Feng--Zou--Wang--Zhu Behrens--Fisher statistic kernel
//'
//' @param x,y Numeric observation-by-variable matrices.
//' @return Internal list of statistic components.
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List cpp_fzwz_bf_two_sample(const arma::mat& x,
                                  const arma::mat& y) {
  fzwz_require_finite_matrix(x, "x");
  fzwz_require_finite_matrix(y, "y");
  const arma::uword n1 = x.n_rows;
  const arma::uword n2 = y.n_rows;
  const arma::uword p = x.n_cols;
  if (n1 < 6 || n2 < 6) {
    Rcpp::stop(
      "Feng-Zou-Wang-Zhu requires at least six observations in each "
      "group for its leave-four-out trace estimators."
    );
  }
  if (p < 1 || y.n_cols != p) {
    Rcpp::stop("`x` and `y` must have the same positive number of columns.");
  }

  const FzwzData data = fzwz_standardize_columns(x, y);
  const FzwzMoments moments1 = fzwz_moments(data.x, n1, p);
  const FzwzMoments moments2 = fzwz_moments(data.y, n2, p);
  const long double n1_ld = static_cast<long double>(n1);
  const long double n2_ld = static_cast<long double>(n2);
  const long double gamma = n1_ld / n2_ld;

  std::vector<long double> denominator(p);
  std::vector<long double> coordinate_a(p);
  std::vector<long double> coordinate_contribution(p);
  long double minimum_full_denominator =
    std::numeric_limits<long double>::infinity();
  long double statistic = 0.0L;
  long double bias1 = 0.0L;
  long double bias2 = 0.0L;

  for (arma::uword j = 0; j < p; ++j) {
    const long double variance1 = moments1.variance[j];
    const long double variance2 = moments2.variance[j];
    const long double d = variance1 + gamma * variance2;
    if (!std::isfinite(d) || d <= 0.0L) {
      Rcpp::stop(
        "Feng-Zou-Wang-Zhu requires every full-sample combined marginal "
        "variance sigma1.hat^2 + gamma * sigma2.hat^2 to be finite and "
        "strictly positive; no ridge is applied."
      );
    }
    denominator[j] = d;
    minimum_full_denominator = std::min(minimum_full_denominator, d);
    const long double difference = moments1.mean[j] - moments2.mean[j];
    const long double a = difference * difference -
      variance1 / n1_ld - variance2 / n2_ld;
    coordinate_a[j] = a;
    coordinate_contribution[j] = a / d;
    statistic += coordinate_contribution[j];

    const long double d2 = d * d;
    const long double d3 = d2 * d;
    bias1 += 2.0L * variance1 * variance1 /
      (n1_ld * (n1_ld - 1.0L) * d2) +
      2.0L * gamma * variance2 * variance2 /
      (n2_ld * (n2_ld - 1.0L) * d2);
    const long double skewness_combination =
      moments1.third[j] / n1_ld -
      gamma * moments2.third[j] / n2_ld;
    bias2 += 2.0L * skewness_combination * skewness_combination / d3;
  }

  const long double null_centering = bias1 + bias2;
  const FzwzWithinTrace trace1 = fzwz_within_trace(
    moments1, moments2, true, gamma
  );
  const FzwzWithinTrace trace2 = fzwz_within_trace(
    moments2, moments1, false, gamma
  );
  const FzwzCrossTrace trace12 = fzwz_cross_trace(
    moments1, moments2, gamma
  );
  const long double coefficient1 =
    2.0L / (n1_ld * (n1_ld - 1.0L));
  const long double coefficient2 =
    2.0L / (n2_ld * (n2_ld - 1.0L));
  const long double coefficient12 = 4.0L / (n1_ld * n2_ld);
  const long double variance = coefficient1 * trace1.estimate +
    coefficient2 * trace2.estimate +
    coefficient12 * trace12.estimate;
  if (!std::isfinite(variance) || variance <= 0.0L) {
    Rcpp::stop(
      "Feng-Zou-Wang-Zhu requires a finite, strictly positive published "
      "leave-out variance estimate; no absolute-value repair or variance "
      "floor is applied."
    );
  }
  const long double z = (statistic - null_centering) /
    std::sqrt(variance);

  arma::vec mean_x(p);
  arma::vec mean_y(p);
  arma::vec difference(p);
  arma::vec variance1_original(p);
  arma::vec variance2_original(p);
  arma::vec denominator_original(p);
  arma::vec coordinate_a_original(p);
  arma::vec coordinate_a_scaled(p);
  arma::vec contribution(p);
  arma::vec third1_scaled(p);
  arma::vec third2_scaled(p);
  arma::vec column_scale(p);
  for (arma::uword j = 0; j < p; ++j) {
    const long double scale = data.scale[j];
    const long double scale2 = scale * scale;
    mean_x(j) = fzwz_report_long_double(
      data.origin[j] + scale * moments1.mean[j]
    );
    mean_y(j) = fzwz_report_long_double(
      data.origin[j] + scale * moments2.mean[j]
    );
    difference(j) = fzwz_report_long_double(
      scale * (moments1.mean[j] - moments2.mean[j])
    );
    variance1_original(j) = fzwz_report_long_double(
      scale2 * moments1.variance[j]
    );
    variance2_original(j) = fzwz_report_long_double(
      scale2 * moments2.variance[j]
    );
    denominator_original(j) = fzwz_report_long_double(
      scale2 * denominator[j]
    );
    coordinate_a_original(j) = fzwz_report_long_double(
      scale2 * coordinate_a[j]
    );
    coordinate_a_scaled(j) = static_cast<double>(coordinate_a[j]);
    contribution(j) = static_cast<double>(coordinate_contribution[j]);
    third1_scaled(j) = static_cast<double>(moments1.third[j]);
    third2_scaled(j) = static_cast<double>(moments2.third[j]);
    column_scale(j) = fzwz_report_long_double(scale);
  }

  const long double minimum_leaveout_denominator = std::min(
    trace12.minimum_denominator,
    std::min(trace1.minimum_denominator, trace2.minimum_denominator)
  );
  const long double p4_n1 = 24.0L * trace1.combinations;
  const long double p4_n2 = 24.0L * trace2.combinations;
  const long double p2_n1 = n1_ld * (n1_ld - 1.0L);
  const long double p2_n2 = n2_ld * (n2_ld - 1.0L);

  return Rcpp::List::create(
    Rcpp::Named("z") = fzwz_checked_double(z, "standardised statistic"),
    Rcpp::Named("T_BF") = fzwz_checked_double(
      statistic, "initial statistic T_BF"
    ),
    Rcpp::Named("Q3") = fzwz_checked_double(
      n1_ld * statistic, "Q3 statistic"
    ),
    Rcpp::Named("null_centering") = fzwz_checked_double(
      null_centering, "estimated asymptotic null centering"
    ),
    Rcpp::Named("bias1") = fzwz_checked_double(bias1, "b1 correction"),
    Rcpp::Named("bias2") = fzwz_checked_double(bias2, "b2 correction"),
    Rcpp::Named("variance") = fzwz_checked_positive_double(
      variance, "normalising variance"
    ),
    Rcpp::Named("trace1") = fzwz_checked_double(
      trace1.estimate, "group-1 leave-four-out trace"
    ),
    Rcpp::Named("trace2") = fzwz_checked_double(
      trace2.estimate, "group-2 leave-four-out trace"
    ),
    Rcpp::Named("trace12") = fzwz_checked_double(
      trace12.estimate, "two-plus-two leave-out cross trace"
    ),
    Rcpp::Named("variance_coefficients") = Rcpp::NumericVector::create(
      Rcpp::Named("group1") = static_cast<double>(coefficient1),
      Rcpp::Named("group2") = static_cast<double>(coefficient2),
      Rcpp::Named("cross") = static_cast<double>(coefficient12)
    ),
    Rcpp::Named("mean_x") = mean_x,
    Rcpp::Named("mean_y") = mean_y,
    Rcpp::Named("difference") = difference,
    Rcpp::Named("variance1_diagonal") = variance1_original,
    Rcpp::Named("variance2_diagonal") = variance2_original,
    Rcpp::Named("D_hat_diagonal") = denominator_original,
    Rcpp::Named("A_coordinate") = coordinate_a_original,
    Rcpp::Named("A_coordinate_scaled") = coordinate_a_scaled,
    Rcpp::Named("coordinate_contribution") = contribution,
    Rcpp::Named("third1_scaled") = third1_scaled,
    Rcpp::Named("third2_scaled") = third2_scaled,
    Rcpp::Named("gamma") = static_cast<double>(gamma),
    Rcpp::Named("n1") = static_cast<double>(n1),
    Rcpp::Named("n2") = static_cast<double>(n2),
    Rcpp::Named("p") = static_cast<double>(p),
    Rcpp::Named("column_scale") = column_scale,
    Rcpp::Named("minimum_full_D_scaled") = fzwz_checked_double(
      minimum_full_denominator, "minimum full-sample combined variance"
    ),
    Rcpp::Named("minimum_leaveout_D_scaled") = fzwz_checked_double(
      minimum_leaveout_denominator, "minimum leave-out combined variance"
    ),
    Rcpp::Named("within1_combinations") = static_cast<double>(
      trace1.combinations
    ),
    Rcpp::Named("within2_combinations") = static_cast<double>(
      trace2.combinations
    ),
    Rcpp::Named("cross_pair_combinations") = static_cast<double>(
      trace12.pair_combinations
    ),
    Rcpp::Named("P4_n1") = static_cast<double>(p4_n1),
    Rcpp::Named("P4_n2") = static_cast<double>(p4_n2),
    Rcpp::Named("P2_n1") = static_cast<double>(p2_n1),
    Rcpp::Named("P2_n2") = static_cast<double>(p2_n2)
  );
}
