// Chapter 2 leave-out high-dimensional location tests.
//
// Park--Ayyala and Chen--Qin both require quadratic pair sums whose naive
// implementations repeatedly form leave-out samples. These kernels retain
// the published estimators while using precomputed first and second moments,
// symmetric pair loops, and long-double accumulation.

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp17)]]

#include <RcppArmadillo.h>
#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

namespace {

void check_finite_matrix_leaveout(const arma::mat& x, const char* name) {
  if (!x.is_finite()) {
    Rcpp::stop("`%s` must contain only finite values.", name);
  }
}

void check_finite_vector_leaveout(const arma::vec& x, const char* name) {
  if (!x.is_finite()) {
    Rcpp::stop("`%s` must contain only finite values.", name);
  }
}

double checked_double_leaveout(long double value, const char* quantity) {
  const double result = static_cast<double>(value);
  if (!std::isfinite(result)) {
    Rcpp::stop("%s is not representable as a finite double.", quantity);
  }
  return result;
}

arma::mat scaled_residuals_leaveout(const arma::mat& x,
                                    const arma::vec& mu,
                                    arma::vec& reported_scales,
                                    arma::vec& sample_mean,
                                    arma::vec& mean_difference) {
  const arma::uword n = x.n_rows;
  const arma::uword p = x.n_cols;
  const long double n_ld = static_cast<long double>(n);
  arma::mat residuals(n, p);
  reported_scales.set_size(p);
  sample_mean.set_size(p);
  mean_difference.set_size(p);

  for (arma::uword j = 0; j < p; ++j) {
    const long double mu_ld = static_cast<long double>(mu(j));
    const long double anchor = static_cast<long double>(x(0, j));
    double coordinate_scale = std::abs(mu(j));
    for (arma::uword i = 0; i < n; ++i) {
      coordinate_scale = std::max(coordinate_scale, std::abs(x(i, j)));
    }
    if (coordinate_scale == 0.0) {
      coordinate_scale = 1.0;
    }
    long double maximum_residual = 0.0L;
    long double residual_sum = 0.0L;
    long double anchored_sum = 0.0L;
    bool direct_path_is_finite = true;
    for (arma::uword i = 0; i < n; ++i) {
      const long double value = static_cast<long double>(x(i, j));
      const long double residual = value - mu_ld;
      const long double anchored = value - anchor;
      if (!std::isfinite(residual) || !std::isfinite(anchored)) {
        direct_path_is_finite = false;
      }
      maximum_residual = std::max(maximum_residual, std::abs(residual));
      residual_sum += residual;
      anchored_sum += anchored;
      if (!std::isfinite(residual_sum) || !std::isfinite(anchored_sum)) {
        direct_path_is_finite = false;
      }
    }
    direct_path_is_finite = direct_path_is_finite &&
      std::isfinite(maximum_residual) &&
      maximum_residual <=
        static_cast<long double>(std::numeric_limits<double>::max());

    long double scale = 1.0L;
    long double sample_mean_ld = 0.0L;
    long double mean_difference_ld = 0.0L;
    if (direct_path_is_finite) {
      scale = maximum_residual > 0.0L ? maximum_residual : 1.0L;
      for (arma::uword i = 0; i < n; ++i) {
        residuals(i, j) = static_cast<double>(
          (static_cast<long double>(x(i, j)) - mu_ld) / scale
        );
      }
      sample_mean_ld = anchor + anchored_sum / n_ld;
      mean_difference_ld = residual_sum / n_ld;
    } else {
      // Some ABIs implement long double with the same range as double. In
      // that case opposite finite extremes can overflow before division.
      // Scaling the operands separately retains a finite direction while
      // leaving the coordinatewise scale-invariant PA statistic unchanged.
      scale = static_cast<long double>(coordinate_scale);
      const long double mu_scaled = mu_ld / scale;
      long double scaled_x_sum = 0.0L;
      long double scaled_residual_sum = 0.0L;
      for (arma::uword i = 0; i < n; ++i) {
        const long double x_scaled =
          static_cast<long double>(x(i, j)) / scale;
        const long double residual_scaled = x_scaled - mu_scaled;
        residuals(i, j) = static_cast<double>(residual_scaled);
        scaled_x_sum += x_scaled;
        scaled_residual_sum += residual_scaled;
      }
      sample_mean_ld = scale * scaled_x_sum / n_ld;
      mean_difference_ld = scale * scaled_residual_sum / n_ld;
    }
    reported_scales(j) = static_cast<double>(scale);
    sample_mean(j) = static_cast<double>(sample_mean_ld);
    mean_difference(j) = static_cast<double>(mean_difference_ld);
  }
  return residuals;
}

double common_global_scale_leaveout(const arma::mat& x,
                                    const arma::mat& y) {
  double scale = 0.0;
  for (arma::uword j = 0; j < x.n_cols; ++j) {
    for (arma::uword i = 0; i < x.n_rows; ++i) {
      scale = std::max(scale, std::abs(x(i, j)));
    }
    for (arma::uword i = 0; i < y.n_rows; ++i) {
      scale = std::max(scale, std::abs(y(i, j)));
    }
  }
  return scale > 0.0 ? scale : 1.0;
}

std::vector<long double> column_sums_leaveout(const arma::mat& x) {
  std::vector<long double> sums(x.n_cols, 0.0L);
  for (arma::uword i = 0; i < x.n_rows; ++i) {
    for (arma::uword j = 0; j < x.n_cols; ++j) {
      sums[j] += static_cast<long double>(x(i, j));
    }
  }
  return sums;
}

std::vector<long double> column_sumsq_leaveout(const arma::mat& x) {
  std::vector<long double> sums(x.n_cols, 0.0L);
  for (arma::uword i = 0; i < x.n_rows; ++i) {
    for (arma::uword j = 0; j < x.n_cols; ++j) {
      const long double value = static_cast<long double>(x(i, j));
      sums[j] += value * value;
    }
  }
  return sums;
}

long double scaled_anchor_difference_leaveout(double value,
                                              double anchor,
                                              long double scale) {
  const long double raw_difference =
    static_cast<long double>(value) - static_cast<long double>(anchor);
  if (std::isfinite(raw_difference)) {
    return raw_difference / scale;
  }
  return static_cast<long double>(value) / scale -
    static_cast<long double>(anchor) / scale;
}

long double cq_within_trace_estimator_leaveout(const arma::mat& x,
                                               long double scale) {
  const arma::uword n = x.n_rows;
  const arma::uword p = x.n_cols;
  const long double n_ld = static_cast<long double>(n);
  const long double remaining_n = n_ld - 2.0L;
  std::vector<long double> delta_sums(p, 0.0L);
  for (arma::uword k = 0; k < p; ++k) {
    for (arma::uword i = 0; i < n; ++i) {
      delta_sums[k] += scaled_anchor_difference_leaveout(
        x(i, k), x(0, k), scale
      );
    }
  }
  long double ordered_sum = 0.0L;

  for (arma::uword i = 0; i + 1 < n; ++i) {
    for (arma::uword j = i + 1; j < n; ++j) {
      long double left = 0.0L;
      long double right = 0.0L;
      for (arma::uword k = 0; k < p; ++k) {
        const long double delta_i = scaled_anchor_difference_leaveout(
          x(i, k), x(0, k), scale
        );
        const long double delta_j = scaled_anchor_difference_leaveout(
          x(j, k), x(0, k), scale
        );
        const long double mean_delta_minus =
          (delta_sums[k] - delta_i - delta_j) / remaining_n;
        const long double xi_scaled =
          static_cast<long double>(x(i, k)) / scale;
        const long double xj_scaled =
          static_cast<long double>(x(j, k)) / scale;
        left += xi_scaled * (delta_j - mean_delta_minus);
        right += xj_scaled * (delta_i - mean_delta_minus);
      }
      ordered_sum += 2.0L * left * right;
    }
  }
  return ordered_sum / (n_ld * (n_ld - 1.0L));
}

long double cq_cross_trace_estimator_leaveout(const arma::mat& x,
                                              const arma::mat& y,
                                              long double scale) {
  const arma::uword n1 = x.n_rows;
  const arma::uword n2 = y.n_rows;
  const arma::uword p = x.n_cols;
  const long double n1_ld = static_cast<long double>(n1);
  const long double n2_ld = static_cast<long double>(n2);
  std::vector<long double> mean_x_delta(p, 0.0L);
  std::vector<long double> mean_y_delta(p, 0.0L);

  for (arma::uword k = 0; k < p; ++k) {
    for (arma::uword i = 0; i < n1; ++i) {
      mean_x_delta[k] += scaled_anchor_difference_leaveout(
        x(i, k), x(0, k), scale
      );
    }
    for (arma::uword j = 0; j < n2; ++j) {
      mean_y_delta[k] += scaled_anchor_difference_leaveout(
        y(j, k), y(0, k), scale
      );
    }
    mean_x_delta[k] /= n1_ld;
    mean_y_delta[k] /= n2_ld;
  }

  long double cross_square_sum = 0.0L;
  for (arma::uword i = 0; i < n1; ++i) {
    for (arma::uword j = 0; j < n2; ++j) {
      long double cross = 0.0L;
      for (arma::uword k = 0; k < p; ++k) {
        const long double centered_x =
          scaled_anchor_difference_leaveout(x(i, k), x(0, k), scale) -
          mean_x_delta[k];
        const long double centered_y =
          scaled_anchor_difference_leaveout(y(j, k), y(0, k), scale) -
          mean_y_delta[k];
        cross += centered_x * centered_y;
      }
      cross_square_sum += cross * cross;
    }
  }

  // The leave-one-out factors n_k/(n_k-1) and the average 1/(n1*n2)
  // simplify exactly to this denominator.
  return cross_square_sum /
    ((n1_ld - 1.0L) * (n2_ld - 1.0L));
}

}  // namespace


//' Park--Ayyala one-sample statistic kernel
//'
//' @param x Numeric observation-by-variable matrix.
//' @param mu Numeric null-mean vector.
//' @return Internal list of statistic components.
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List cpp_park_ayyala_one_sample(const arma::mat& x,
                                      const arma::vec& mu) {
  check_finite_matrix_leaveout(x, "x");
  check_finite_vector_leaveout(mu, "mu");
  const arma::uword n = x.n_rows;
  const arma::uword p = x.n_cols;
  if (n < 6) {
    Rcpp::stop("Park-Ayyala requires at least six observations.");
  }
  if (p < 1 || mu.n_elem != p) {
    Rcpp::stop("`mu` must have one value per column of `x`.");
  }

  arma::vec column_scale;
  arma::vec sample_mean;
  arma::vec difference;
  const arma::mat z = scaled_residuals_leaveout(
    x, mu, column_scale, sample_mean, difference
  );

  const std::vector<long double> sums = column_sums_leaveout(z);
  const std::vector<long double> sums_sq = column_sumsq_leaveout(z);
  const long double n_ld = static_cast<long double>(n);
  const long double remaining_n = n_ld - 2.0L;
  const long double remaining_df = n_ld - 3.0L;
  long double T1 = 0.0L;
  long double T2 = 0.0L;
  long double min_variance = std::numeric_limits<long double>::infinity();

  for (arma::uword i = 0; i + 1 < n; ++i) {
    for (arma::uword j = i + 1; j < n; ++j) {
      long double pair_product = 0.0L;
      long double left = 0.0L;
      long double right = 0.0L;
      for (arma::uword k = 0; k < p; ++k) {
        const long double zi = static_cast<long double>(z(i, k));
        const long double zj = static_cast<long double>(z(j, k));
        const long double remaining_sum = sums[k] - zi - zj;
        const long double remaining_sumsq = sums_sq[k] - zi * zi - zj * zj;
        const long double variance = (
          remaining_sumsq - remaining_sum * remaining_sum / remaining_n
        ) / remaining_df;
        if (!std::isfinite(variance) || variance <= 0.0L) {
          Rcpp::stop(
            "Park-Ayyala requires every leave-two-out marginal variance "
            "to be finite and strictly positive; no ridge is applied."
          );
        }
        min_variance = std::min(min_variance, variance);
        const long double inverse_variance = 1.0L / variance;
        const long double mean_minus = remaining_sum / remaining_n;
        pair_product += zi * zj * inverse_variance;
        left += zi * (zj - mean_minus) * inverse_variance;
        right += zj * (zi - mean_minus) * inverse_variance;
      }
      T1 += 2.0L * pair_product;
      T2 += 2.0L * left * right;
    }
  }

  const long double finite_sample_correction =
    (n_ld - 5.0L) / (n_ld - 3.0L);
  const long double bias_factor = finite_sample_correction /
    (n_ld * (n_ld - 1.0L));
  const long double U = bias_factor * T1;
  const long double ordered_pair_count = n_ld * (n_ld - 1.0L);
  // T2 / {n(n-1)} estimates tr(R^2); the U-statistic variance contributes
  // another 2 / {n(n-1)} and the same squared finite-sample correction used
  // by U. Consequently that correction cancels from the final Z statistic.
  const long double variance = finite_sample_correction *
    finite_sample_correction * 2.0L * T2 /
    (ordered_pair_count * ordered_pair_count);
  if (!std::isfinite(variance) || variance <= 0.0L) {
    Rcpp::stop(
      "Park-Ayyala requires a finite, strictly positive estimated null "
      "variance; no absolute-value repair or variance floor is applied."
    );
  }
  const long double z_statistic = U / std::sqrt(variance);

  return Rcpp::List::create(
    Rcpp::Named("z") = checked_double_leaveout(z_statistic, "Park-Ayyala Z"),
    Rcpp::Named("U") = checked_double_leaveout(U, "Park-Ayyala U"),
    Rcpp::Named("T1") = checked_double_leaveout(T1, "Park-Ayyala T1"),
    Rcpp::Named("T2") = checked_double_leaveout(T2, "Park-Ayyala T2"),
    Rcpp::Named("variance") = checked_double_leaveout(
      variance, "Park-Ayyala variance"
    ),
    Rcpp::Named("bias_factor") = static_cast<double>(bias_factor),
    Rcpp::Named("finite_sample_correction") =
      static_cast<double>(finite_sample_correction),
    Rcpp::Named("sample_mean") = sample_mean,
    Rcpp::Named("difference") = difference,
    Rcpp::Named("N") = static_cast<double>(n),
    Rcpp::Named("p") = static_cast<double>(p),
    Rcpp::Named("leaveout_pairs") = static_cast<double>(n) *
      static_cast<double>(n - 1) / 2.0,
    Rcpp::Named("min_leaveout_variance") = static_cast<double>(min_variance),
    Rcpp::Named("column_scale") = column_scale
  );
}


//' Chen--Qin two-sample statistic kernel
//'
//' @param x,y Numeric observation-by-variable matrices.
//' @return Internal list of statistic components.
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List cpp_chen_qin_two_sample(const arma::mat& x,
                                   const arma::mat& y) {
  check_finite_matrix_leaveout(x, "x");
  check_finite_matrix_leaveout(y, "y");
  const arma::uword n1 = x.n_rows;
  const arma::uword n2 = y.n_rows;
  const arma::uword p = x.n_cols;
  if (n1 < 3 || n2 < 3) {
    Rcpp::stop("Chen-Qin requires at least three observations in each group.");
  }
  if (p < 1 || y.n_cols != p) {
    Rcpp::stop("`x` and `y` must have the same positive number of columns.");
  }

  const double global_scale = common_global_scale_leaveout(x, y);
  const long double global_scale_ld =
    static_cast<long double>(global_scale);
  const long double n1_ld = static_cast<long double>(n1);
  const long double n2_ld = static_cast<long double>(n2);

  // The CQ numerator is translation invariant even though its published
  // variance estimator is not. Compute the numerator in long double after
  // subtracting one common anchor, rather than as sumsq - sum^2 / n on data
  // whose large common level can erase the within-group variation.
  long double squared_difference_scaled = 0.0L;
  long double trace_S1_scaled = 0.0L;
  long double trace_S2_scaled = 0.0L;
  arma::vec mean_x(p);
  arma::vec mean_y(p);
  arma::vec difference(p);
  for (arma::uword k = 0; k < p; ++k) {
    const double anchor = x(0, k);
    long double sum_x_delta = 0.0L;
    long double sum_y_delta = 0.0L;
    for (arma::uword i = 0; i < n1; ++i) {
      sum_x_delta += scaled_anchor_difference_leaveout(
        x(i, k), anchor, global_scale_ld
      );
    }
    for (arma::uword i = 0; i < n2; ++i) {
      sum_y_delta += scaled_anchor_difference_leaveout(
        y(i, k), anchor, global_scale_ld
      );
    }
    const long double mean_x_delta = sum_x_delta / n1_ld;
    const long double mean_y_delta = sum_y_delta / n2_ld;
    const long double mean_difference = mean_x_delta - mean_y_delta;
    const long double anchor_ld = static_cast<long double>(anchor);
    mean_x(k) = static_cast<double>(
      anchor_ld + mean_x_delta * global_scale_ld
    );
    mean_y(k) = static_cast<double>(
      anchor_ld + mean_y_delta * global_scale_ld
    );
    difference(k) = static_cast<double>(
      mean_difference * global_scale_ld
    );
    squared_difference_scaled += mean_difference * mean_difference;
    for (arma::uword i = 0; i < n1; ++i) {
      const long double centered =
        scaled_anchor_difference_leaveout(
          x(i, k), anchor, global_scale_ld
        ) - mean_x_delta;
      trace_S1_scaled += centered * centered /
        (n1_ld - 1.0L);
    }
    for (arma::uword i = 0; i < n2; ++i) {
      const long double centered =
        scaled_anchor_difference_leaveout(
          y(i, k), anchor, global_scale_ld
        ) - mean_y_delta;
      trace_S2_scaled += centered * centered /
        (n2_ld - 1.0L);
    }
  }
  const long double T_scaled = squared_difference_scaled -
    trace_S1_scaled / n1_ld - trace_S2_scaled / n2_ld;
  const long double scale2 = global_scale_ld * global_scale_ld;
  const long double scale4 = scale2 * scale2;

  const long double A1_scaled = cq_within_trace_estimator_leaveout(
    x, global_scale_ld
  );
  const long double A2_scaled = cq_within_trace_estimator_leaveout(
    y, global_scale_ld
  );
  const long double A12_scaled = cq_cross_trace_estimator_leaveout(
    x, y, global_scale_ld
  );
  const long double variance_scaled =
    2.0L * A1_scaled / (n1_ld * (n1_ld - 1.0L)) +
    2.0L * A2_scaled / (n2_ld * (n2_ld - 1.0L)) +
    4.0L * A12_scaled / (n1_ld * n2_ld);
  if (!std::isfinite(variance_scaled) || variance_scaled <= 0.0L) {
    Rcpp::stop(
      "Chen-Qin requires a finite, strictly positive original leave-out "
      "variance estimate; no absolute-value repair or variance floor is "
      "applied."
    );
  }
  const long double z_statistic = T_scaled / std::sqrt(variance_scaled);

  return Rcpp::List::create(
    Rcpp::Named("z") = checked_double_leaveout(z_statistic, "Chen-Qin Z"),
    Rcpp::Named("T") = static_cast<double>(T_scaled * scale2),
    Rcpp::Named("A1") = static_cast<double>(A1_scaled * scale4),
    Rcpp::Named("A2") = static_cast<double>(A2_scaled * scale4),
    Rcpp::Named("A12") = static_cast<double>(A12_scaled * scale4),
    Rcpp::Named("variance") = static_cast<double>(variance_scaled * scale4),
    Rcpp::Named("trace_S1") = static_cast<double>(trace_S1_scaled * scale2),
    Rcpp::Named("trace_S2") = static_cast<double>(trace_S2_scaled * scale2),
    Rcpp::Named("T_scaled") = static_cast<double>(T_scaled),
    Rcpp::Named("A1_scaled") = static_cast<double>(A1_scaled),
    Rcpp::Named("A2_scaled") = static_cast<double>(A2_scaled),
    Rcpp::Named("A12_scaled") = static_cast<double>(A12_scaled),
    Rcpp::Named("variance_scaled") = static_cast<double>(variance_scaled),
    Rcpp::Named("mean_x") = mean_x,
    Rcpp::Named("mean_y") = mean_y,
    Rcpp::Named("difference") = difference,
    Rcpp::Named("n1") = static_cast<double>(n1),
    Rcpp::Named("n2") = static_cast<double>(n2),
    Rcpp::Named("p") = static_cast<double>(p),
    Rcpp::Named("global_scale") = global_scale
  );
}
