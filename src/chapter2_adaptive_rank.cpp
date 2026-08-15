// Zhang--Feng adaptive marginal-rank tests.
//
// The score kernels below implement only the two statistics in Zhang and
// Feng (2024): the Ouyang et al. squared-rank sum component and the marginal
// rank maximum.  They deliberately do not manufacture the L_q family that is
// absent from the primary paper.  Exact ties (and one-sample zero residuals)
// are rejected because the published null moments assume continuous margins.

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp17)]]

#include <RcppArmadillo.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

namespace {

void zf_require_finite_matrix(const arma::mat& x, const char* name) {
  if (!x.is_finite()) {
    Rcpp::stop("`%s` must contain only finite values.", name);
  }
}

void zf_neumaier_add(long double value,
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

struct ZfRankEntry {
  double value;
  arma::uword row;
  bool from_x;
};

bool zf_rank_entry_less(const ZfRankEntry& lhs,
                        const ZfRankEntry& rhs) {
  if (lhs.value < rhs.value) {
    return true;
  }
  if (rhs.value < lhs.value) {
    return false;
  }
  if (lhs.from_x != rhs.from_x) {
    return lhs.from_x < rhs.from_x;
  }
  return lhs.row < rhs.row;
}

struct ZfMoments {
  double mean_u;
  double variance_u;
  double mean_squared;
  double variance_squared;
};

ZfMoments zf_one_sample_moments(arma::uword n) {
  const long double size = static_cast<long double>(n);
  const long double mean_u = size * (size + 1.0L) / 4.0L;
  const long double variance_u =
    size * (size + 1.0L) * (2.0L * size + 1.0L) / 24.0L;
  const long double variance_squared =
    (6.0L * size + 5.0L * size * size -
     30.0L * std::pow(size, 3) - 25.0L * std::pow(size, 4) +
     24.0L * std::pow(size, 5) + 20.0L * std::pow(size, 6)) /
    1440.0L;
  if (!(variance_u > 0.0L) || !(variance_squared > 0.0L) ||
      !std::isfinite(variance_u) || !std::isfinite(variance_squared)) {
    Rcpp::stop(
      "The one-sample signed-rank null moments are degenerate or non-finite."
    );
  }
  return ZfMoments{
    static_cast<double>(mean_u),
    static_cast<double>(variance_u),
    static_cast<double>(variance_u),
    static_cast<double>(variance_squared)
  };
}

ZfMoments zf_two_sample_moments(arma::uword n, arma::uword m) {
  const long double first = static_cast<long double>(n);
  const long double second = static_cast<long double>(m);
  const long double total = first + second;
  const long double product = first * second;
  const long double mean_u = product / 2.0L;
  const long double variance_u =
    product * (total + 1.0L) / 12.0L;
  const long double variance_squared =
    (product * (5.0L * total + 8.0L) -
     3.0L * total * (total + 1.0L)) *
    (total + 1.0L) * product / 360.0L;
  if (!(variance_u > 0.0L) || !(variance_squared > 0.0L) ||
      !std::isfinite(variance_u) || !std::isfinite(variance_squared)) {
    Rcpp::stop(
      "The two-sample WMW null moments are degenerate or non-finite."
    );
  }
  return ZfMoments{
    static_cast<double>(mean_u),
    static_cast<double>(variance_u),
    static_cast<double>(variance_u),
    static_cast<double>(variance_squared)
  };
}

Rcpp::List zf_standardize_scores(const arma::vec& rank_sum,
                                 const ZfMoments& moments) {
  const arma::uword p = rank_sum.n_elem;
  arma::vec centered(p, arma::fill::zeros);
  arma::vec standardized(p, arma::fill::zeros);
  arma::vec squared(p, arma::fill::zeros);
  arma::vec squared_standardized(p, arma::fill::zeros);
  const long double sd_u =
    std::sqrt(static_cast<long double>(moments.variance_u));
  const long double sd_squared =
    std::sqrt(static_cast<long double>(moments.variance_squared));
  long double squared_total = 0.0L;
  long double squared_correction = 0.0L;

  for (arma::uword j = 0; j < p; ++j) {
    const long double deviation =
      static_cast<long double>(rank_sum(j)) -
      static_cast<long double>(moments.mean_u);
    const long double square = deviation * deviation;
    centered(j) = static_cast<double>(deviation);
    standardized(j) = static_cast<double>(deviation / sd_u);
    squared(j) = static_cast<double>(square);
    squared_standardized(j) = static_cast<double>(
      (square - static_cast<long double>(moments.mean_squared)) /
      sd_squared
    );
    zf_neumaier_add(square, squared_total, squared_correction);
  }

  const long double mean_square =
    (squared_total + squared_correction) /
    static_cast<long double>(p);
  return Rcpp::List::create(
    Rcpp::Named("rank_sum") = rank_sum,
    Rcpp::Named("centered_rank_score") = centered,
    Rcpp::Named("standardized_rank_score") = standardized,
    Rcpp::Named("squared_rank_score") = squared,
    Rcpp::Named("standardized_squared_rank_score") =
      squared_standardized,
    Rcpp::Named("mean_squared_rank_score") =
      static_cast<double>(mean_square),
    Rcpp::Named("null_mean_rank_sum") = moments.mean_u,
    Rcpp::Named("null_variance_rank_sum") = moments.variance_u,
    Rcpp::Named("null_mean_squared_rank_score") =
      moments.mean_squared,
    Rcpp::Named("null_variance_squared_rank_score") =
      moments.variance_squared
  );
}

double zf_parzen_weight(double value) {
  const double magnitude = std::abs(value);
  if (magnitude < 0.5) {
    return 1.0 - 6.0 * magnitude * magnitude +
      6.0 * magnitude * magnitude * magnitude;
  }
  if (magnitude <= 1.0) {
    const double remainder = 1.0 - magnitude;
    return 2.0 * remainder * remainder * remainder;
  }
  return 0.0;
}

}  // namespace


//' One-sample Zhang--Feng marginal signed-rank scores
//'
//' @param x Numeric observation-by-variable matrix.
//' @param mu Numeric null-location vector.
//' @return Internal list of marginal scores and exact untied null moments.
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List cpp_zhang_feng_one_sample_scores(const arma::mat& x,
                                            const arma::vec& mu) {
  zf_require_finite_matrix(x, "x");
  if (!mu.is_finite()) {
    Rcpp::stop("`mu` must contain only finite values.");
  }
  const arma::uword n = x.n_rows;
  const arma::uword p = x.n_cols;
  if (n < 2 || p < 2 || mu.n_elem != p) {
    Rcpp::stop(
      "The one-sample Zhang-Feng score kernel requires n >= 2, p >= 2, "
      "and one null-location entry per variable."
    );
  }

  arma::vec rank_sum(p, arma::fill::zeros);
  arma::uword overflow_fallback_columns = 0;
  std::vector<double> residual(n, 0.0);
  std::vector<ZfRankEntry> ordering(n);

  for (arma::uword j = 0; j < p; ++j) {
    bool subtraction_overflow = false;
    for (arma::uword i = 0; i < n; ++i) {
      residual[i] = x(i, j) - mu(j);
      if (!std::isfinite(residual[i])) {
        subtraction_overflow = true;
      }
    }
    if (subtraction_overflow) {
      ++overflow_fallback_columns;
      double operand_scale = std::abs(mu(j));
      for (arma::uword i = 0; i < n; ++i) {
        operand_scale = std::max(operand_scale, std::abs(x(i, j)));
      }
      if (!(operand_scale > 0.0) || !std::isfinite(operand_scale)) {
        Rcpp::stop(
          "Could not form finite one-sample residuals in variable %llu.",
          static_cast<unsigned long long>(j + 1)
        );
      }
      for (arma::uword i = 0; i < n; ++i) {
        residual[i] = x(i, j) / operand_scale - mu(j) / operand_scale;
      }
    }

    for (arma::uword i = 0; i < n; ++i) {
      if (!std::isfinite(residual[i])) {
        Rcpp::stop(
          "Could not form a finite one-sample residual in variable %llu.",
          static_cast<unsigned long long>(j + 1)
        );
      }
      if (residual[i] == 0.0) {
        Rcpp::stop(
          "The one-sample signed-rank calibration does not allow an exact "
          "zero residual; found one in variable %llu. No zero correction "
          "is specified by Zhang and Feng (2024).",
          static_cast<unsigned long long>(j + 1)
        );
      }
      ordering[i] = ZfRankEntry{std::abs(residual[i]), i, true};
    }
    std::sort(ordering.begin(), ordering.end(), zf_rank_entry_less);
    for (arma::uword r = 1; r < n; ++r) {
      if (ordering[r - 1].value == ordering[r].value) {
        Rcpp::stop(
          "The one-sample signed-rank calibration requires untied absolute "
          "residuals; variable %llu contains an exact tie.",
          static_cast<unsigned long long>(j + 1)
        );
      }
    }
    long double positive_rank_sum = 0.0L;
    for (arma::uword r = 0; r < n; ++r) {
      if (residual[ordering[r].row] > 0.0) {
        positive_rank_sum += static_cast<long double>(r + 1);
      }
    }
    rank_sum(j) = static_cast<double>(positive_rank_sum);
  }

  Rcpp::List output = zf_standardize_scores(
    rank_sum, zf_one_sample_moments(n)
  );
  output["overflow_fallback_columns"] =
    static_cast<double>(overflow_fallback_columns);
  return output;
}


//' Two-sample Zhang--Feng marginal WMW scores
//'
//' @param x First numeric observation-by-variable matrix.
//' @param y Second numeric observation-by-variable matrix.
//' @return Internal list of marginal scores and exact untied null moments.
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List cpp_zhang_feng_two_sample_scores(const arma::mat& x,
                                            const arma::mat& y) {
  zf_require_finite_matrix(x, "x");
  zf_require_finite_matrix(y, "y");
  const arma::uword n = x.n_rows;
  const arma::uword m = y.n_rows;
  const arma::uword p = x.n_cols;
  if (n < 2 || m < 2 || p < 2 || y.n_cols != p) {
    Rcpp::stop(
      "The two-sample Zhang-Feng score kernel requires n_x >= 2, "
      "n_y >= 2, p >= 2, and equal variable counts."
    );
  }

  arma::vec rank_sum(p, arma::fill::zeros);
  std::vector<ZfRankEntry> ordering(n + m);
  for (arma::uword j = 0; j < p; ++j) {
    for (arma::uword i = 0; i < n; ++i) {
      ordering[i] = ZfRankEntry{x(i, j), i, true};
    }
    for (arma::uword i = 0; i < m; ++i) {
      ordering[n + i] = ZfRankEntry{y(i, j), i, false};
    }
    std::sort(ordering.begin(), ordering.end(), zf_rank_entry_less);
    for (arma::uword r = 1; r < n + m; ++r) {
      if (ordering[r - 1].value == ordering[r].value) {
        Rcpp::stop(
          "The two-sample WMW calibration requires untied pooled values; "
          "variable %llu contains an exact tie.",
          static_cast<unsigned long long>(j + 1)
        );
      }
    }
    long double x_rank_sum = 0.0L;
    for (arma::uword r = 0; r < n + m; ++r) {
      if (ordering[r].from_x) {
        x_rank_sum += static_cast<long double>(r + 1);
      }
    }
    const long double offset =
      static_cast<long double>(n) *
      static_cast<long double>(n + 1) / 2.0L;
    rank_sum(j) = static_cast<double>(x_rank_sum - offset);
  }

  return zf_standardize_scores(
    rank_sum, zf_two_sample_moments(n, m)
  );
}


//' Ouyang Parzen long-run variance estimate for rank-square scores
//'
//' @param score Numeric standardized squared-rank score sequence.
//' @param lag Explicit lag-window size L; lags 1 through L - 1 are used.
//' @return Internal list of sample autocovariances, weights, and tau squared.
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List cpp_zhang_feng_parzen_tau(const arma::vec& score, int lag) {
  if (!score.is_finite()) {
    Rcpp::stop("The standardized squared-rank scores must be finite.");
  }
  const arma::uword p = score.n_elem;
  if (p < 2 || lag < 1 || static_cast<arma::uword>(lag) > p) {
    Rcpp::stop("`lag` must be an integer in 1, ..., p.");
  }

  long double mean_total = 0.0L;
  long double mean_correction = 0.0L;
  for (arma::uword j = 0; j < p; ++j) {
    zf_neumaier_add(
      static_cast<long double>(score(j)), mean_total, mean_correction
    );
  }
  const long double score_mean =
    (mean_total + mean_correction) / static_cast<long double>(p);

  const int count = lag - 1;
  Rcpp::IntegerVector lags(count);
  Rcpp::NumericVector autocovariance(count);
  Rcpp::NumericVector weight(count);
  Rcpp::NumericVector weighted_autocovariance(count);
  long double weighted_total = 0.0L;
  long double weighted_correction = 0.0L;

  for (int index = 0; index < count; ++index) {
    const arma::uword k = static_cast<arma::uword>(index + 1);
    long double covariance_total = 0.0L;
    long double covariance_correction = 0.0L;
    for (arma::uword j = 0; j < p - k; ++j) {
      const long double first =
        static_cast<long double>(score(j)) - score_mean;
      const long double second =
        static_cast<long double>(score(j + k)) - score_mean;
      zf_neumaier_add(
        first * second, covariance_total, covariance_correction
      );
    }
    const long double gamma =
      (covariance_total + covariance_correction) /
      static_cast<long double>(p - k);
    const double parzen = zf_parzen_weight(
      static_cast<double>(k) / static_cast<double>(lag)
    );
    const long double contribution =
      static_cast<long double>(parzen) * gamma;
    zf_neumaier_add(
      contribution, weighted_total, weighted_correction
    );
    lags[index] = static_cast<int>(k);
    autocovariance[index] = static_cast<double>(gamma);
    weight[index] = parzen;
    weighted_autocovariance[index] = static_cast<double>(contribution);
  }

  const long double tau_squared =
    1.0L + 2.0L * (weighted_total + weighted_correction);
  return Rcpp::List::create(
    Rcpp::Named("lag_window_size") = lag,
    Rcpp::Named("lags") = lags,
    Rcpp::Named("score_mean") = static_cast<double>(score_mean),
    Rcpp::Named("autocovariance") = autocovariance,
    Rcpp::Named("weight") = weight,
    Rcpp::Named("weighted_autocovariance") = weighted_autocovariance,
    Rcpp::Named("tau_squared") = static_cast<double>(tau_squared),
    Rcpp::Named("denominator") = "p - k"
  );
}
