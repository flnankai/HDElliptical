// Li--Wang--Zou simpler spatial-sign-based two-sample test.
//
// This implements the full-sample feasible SST in Section 2.2 of Li, Wang
// and Zou (2016).  There are no leave-out fits: each group has one diagonal
// Hettmansperger--Randles fit.  The bias and variance plug-ins use the exact
// ordered-pair and cross-pair formulas displayed below Proposition 1.

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp17)]]

#include <RcppArmadillo.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <string>
#include <utility>
#include <vector>

namespace {

void lwz_neumaier_add(const long double value,
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

double lwz_checked_double(const long double value, const char* quantity) {
  const long double maximum = static_cast<long double>(
    std::numeric_limits<double>::max()
  );
  if (!std::isfinite(value) || value > maximum || value < -maximum) {
    Rcpp::stop(
      "Li-Wang-Zou SST produced a non-finite %s; no numerical floor, "
      "ridge, absolute-value repair, or pseudoinverse is applied.",
      quantity
    );
  }
  return static_cast<double>(value);
}

long double lwz_dot(const arma::rowvec& x, const arma::rowvec& y) {
  long double total = 0.0L;
  long double correction = 0.0L;
  for (arma::uword j = 0; j < x.n_elem; ++j) {
    lwz_neumaier_add(
      static_cast<long double>(x(j)) * static_cast<long double>(y(j)),
      total, correction
    );
  }
  return total + correction;
}

long double lwz_bridge_dot(const arma::rowvec& x,
                           const arma::vec& bridge,
                           const arma::rowvec& y) {
  long double total = 0.0L;
  long double correction = 0.0L;
  for (arma::uword j = 0; j < x.n_elem; ++j) {
    lwz_neumaier_add(
      static_cast<long double>(x(j)) *
        static_cast<long double>(bridge(j)) *
        static_cast<long double>(y(j)),
      total, correction
    );
  }
  return total + correction;
}

struct LwzPreparedData {
  arma::mat x;
  arma::mat y;
  arma::vec anchor;
  arma::vec scale_base;
  arma::vec scale_ratio;
  arma::vec log_scale;
  arma::uvec overflow_fallback;
};

LwzPreparedData lwz_prepare_data(const arma::mat& x,
                                 const arma::mat& y) {
  const arma::uword n1 = x.n_rows;
  const arma::uword n2 = y.n_rows;
  const arma::uword p = x.n_cols;
  const long double double_max = static_cast<long double>(
    std::numeric_limits<double>::max()
  );
  LwzPreparedData out;
  out.x.set_size(n1, p);
  out.y.set_size(n2, p);
  out.anchor.set_size(p);
  out.scale_base.ones(p);
  out.scale_ratio.ones(p);
  out.log_scale.zeros(p);
  out.overflow_fallback.zeros(p);

  for (arma::uword j = 0; j < p; ++j) {
    out.anchor(j) = x(0, j);
    std::vector<long double> x_difference(n1);
    std::vector<long double> y_difference(n2);
    long double maximum = 0.0L;
    bool direct_ok = true;
    for (arma::uword i = 0; i < n1; ++i) {
      x_difference[i] = static_cast<long double>(x(i, j)) -
        static_cast<long double>(out.anchor(j));
      direct_ok = direct_ok && std::isfinite(x_difference[i]);
      maximum = std::max(maximum, std::abs(x_difference[i]));
    }
    for (arma::uword i = 0; i < n2; ++i) {
      y_difference[i] = static_cast<long double>(y(i, j)) -
        static_cast<long double>(out.anchor(j));
      direct_ok = direct_ok && std::isfinite(y_difference[i]);
      maximum = std::max(maximum, std::abs(y_difference[i]));
    }
    direct_ok = direct_ok && std::isfinite(maximum) &&
      maximum > 0.0L && maximum <= double_max;

    if (direct_ok) {
      const double scale = static_cast<double>(maximum);
      out.scale_base(j) = scale;
      out.scale_ratio(j) = 1.0;
      out.log_scale(j) = std::log(scale);
      for (arma::uword i = 0; i < n1; ++i) {
        out.x(i, j) = static_cast<double>(x_difference[i] / maximum);
      }
      for (arma::uword i = 0; i < n2; ++i) {
        out.y(i, j) = static_cast<double>(y_difference[i] / maximum);
      }
    } else {
      out.overflow_fallback(j) = 1u;
      double operand_scale = std::abs(out.anchor(j));
      for (arma::uword i = 0; i < n1; ++i) {
        operand_scale = std::max(operand_scale, std::abs(x(i, j)));
      }
      for (arma::uword i = 0; i < n2; ++i) {
        operand_scale = std::max(operand_scale, std::abs(y(i, j)));
      }
      if (!(operand_scale > 0.0) || !std::isfinite(operand_scale)) {
        Rcpp::stop(
          "Li-Wang-Zou SST could not construct a finite common operand "
          "scale for variable %llu.",
          static_cast<unsigned long long>(j + 1)
        );
      }
      const double anchor_scaled = out.anchor(j) / operand_scale;
      double normalized_maximum = 0.0;
      for (arma::uword i = 0; i < n1; ++i) {
        out.x(i, j) = x(i, j) / operand_scale - anchor_scaled;
        normalized_maximum = std::max(
          normalized_maximum, std::abs(out.x(i, j))
        );
      }
      for (arma::uword i = 0; i < n2; ++i) {
        out.y(i, j) = y(i, j) / operand_scale - anchor_scaled;
        normalized_maximum = std::max(
          normalized_maximum, std::abs(out.y(i, j))
        );
      }
      if (!(normalized_maximum > 0.0) ||
          !std::isfinite(normalized_maximum)) {
        Rcpp::stop(
          "Li-Wang-Zou SST requires pooled positive variation in "
          "variable %llu.",
          static_cast<unsigned long long>(j + 1)
        );
      }
      out.x.col(j) /= normalized_maximum;
      out.y.col(j) /= normalized_maximum;
      out.scale_base(j) = operand_scale;
      out.scale_ratio(j) = normalized_maximum;
      out.log_scale(j) = std::log(operand_scale) +
        std::log(normalized_maximum);
    }

    double minimum = out.x(0, j);
    double column_maximum = out.x(0, j);
    for (arma::uword i = 1; i < n1; ++i) {
      minimum = std::min(minimum, out.x(i, j));
      column_maximum = std::max(column_maximum, out.x(i, j));
    }
    for (arma::uword i = 0; i < n2; ++i) {
      minimum = std::min(minimum, out.y(i, j));
      column_maximum = std::max(column_maximum, out.y(i, j));
    }
    if (!(column_maximum > minimum)) {
      Rcpp::stop(
        "Li-Wang-Zou SST requires pooled positive variation in variable "
        "%llu.", static_cast<unsigned long long>(j + 1)
      );
    }
  }
  if (!out.x.is_finite() || !out.y.is_finite()) {
    Rcpp::stop("Li-Wang-Zou SST internal data standardisation failed.");
  }
  return out;
}

arma::vec lwz_restore_location(const arma::vec& location,
                               const LwzPreparedData& prepared) {
  arma::vec answer(location.n_elem);
  for (arma::uword j = 0; j < location.n_elem; ++j) {
    const long double scale =
      static_cast<long double>(prepared.scale_base(j)) *
      static_cast<long double>(prepared.scale_ratio(j));
    const long double value =
      static_cast<long double>(prepared.anchor(j)) +
      scale * static_cast<long double>(location(j));
    answer(j) = lwz_checked_double(value, "restored location");
  }
  return answer;
}

bool lwz_direction(const arma::rowvec& value,
                   const arma::vec& location,
                   const arma::vec& diagonal,
                   const double zero_tol,
                   arma::rowvec& direction,
                   double& radius) {
  const arma::uword p = value.n_elem;
  direction.set_size(p);
  double maximum = 0.0;
  for (arma::uword j = 0; j < p; ++j) {
    if (!(diagonal(j) > 0.0) || !std::isfinite(diagonal(j))) {
      Rcpp::stop(
        "A Li-Wang-Zou SST diagonal iterate is not finite and strictly "
        "positive; no ridge or floor is applied."
      );
    }
    direction(j) = (value(j) - location(j)) /
      std::sqrt(diagonal(j));
    if (!std::isfinite(direction(j))) {
      Rcpp::stop(
        "Li-Wang-Zou SST diagonal standardisation produced a non-finite "
        "value."
      );
    }
    maximum = std::max(maximum, std::abs(direction(j)));
  }
  if (maximum == 0.0) {
    direction.zeros();
    radius = 0.0;
    return false;
  }

  long double square_total = 0.0L;
  long double square_correction = 0.0L;
  for (arma::uword j = 0; j < p; ++j) {
    direction(j) /= maximum;
    const long double component = static_cast<long double>(direction(j));
    lwz_neumaier_add(
      component * component, square_total, square_correction
    );
  }
  const long double norm = std::sqrt(square_total + square_correction);
  if (!(norm > 0.0L) || !std::isfinite(norm)) {
    Rcpp::stop("Li-Wang-Zou SST could not normalise a nonzero residual.");
  }
  radius = lwz_checked_double(
    static_cast<long double>(maximum) * norm,
    "standardised radius"
  );
  if (!(radius > zero_tol)) {
    direction.zeros();
    return false;
  }
  direction /= static_cast<double>(norm);
  return true;
}

void lwz_normalise_diagonal(arma::vec& diagonal, const int group) {
  long double log_total = 0.0L;
  long double log_correction = 0.0L;
  for (arma::uword j = 0; j < diagonal.n_elem; ++j) {
    if (!(diagonal(j) > 0.0) || !std::isfinite(diagonal(j))) {
      Rcpp::stop(
        "Li-Wang-Zou SST group %d has a non-positive or non-finite "
        "diagonal; no ridge or floor is applied.", group
      );
    }
    lwz_neumaier_add(
      std::log(static_cast<long double>(diagonal(j))),
      log_total, log_correction
    );
  }
  const long double mean_log = (log_total + log_correction) /
    static_cast<long double>(diagonal.n_elem);
  for (arma::uword j = 0; j < diagonal.n_elem; ++j) {
    const long double value = std::exp(
      std::log(static_cast<long double>(diagonal(j))) - mean_log
    );
    diagonal(j) = lwz_checked_double(value, "identified diagonal iterate");
    if (!(diagonal(j) > 0.0)) {
      Rcpp::stop(
        "Li-Wang-Zou SST group %d diagonal scale identification "
        "underflowed; no floor is applied.", group
      );
    }
  }
}

struct LwzScore {
  arma::rowvec sign_sum;
  arma::vec sign_square_sum;
  arma::vec radii;
  double inverse_radius_sum;
  double minimum_radius;
  double location_residual;
  double diagonal_residual;
  double score_residual;
};

LwzScore lwz_score(const arma::mat& data,
                   const arma::vec& location,
                   const arma::vec& diagonal,
                   const double zero_tol,
                   const int group) {
  const arma::uword n = data.n_rows;
  const arma::uword p = data.n_cols;
  arma::rowvec sign_sum(p, arma::fill::zeros);
  arma::vec sign_square_sum(p, arma::fill::zeros);
  arma::vec radii(n);
  long double inverse_total = 0.0L;
  long double inverse_correction = 0.0L;
  for (arma::uword i = 0; i < n; ++i) {
    arma::rowvec direction;
    double radius = 0.0;
    if (!lwz_direction(
          data.row(i), location, diagonal, zero_tol, direction, radius
        )) {
      Rcpp::stop(
        "Li-Wang-Zou SST group %d HR iteration is undefined because "
        "training observation %llu has a standardised residual at or "
        "below `zero_tol`; no perturbation is applied.",
        group, static_cast<unsigned long long>(i + 1)
      );
    }
    sign_sum += direction;
    sign_square_sum += arma::square(direction).t();
    radii(i) = radius;
    lwz_neumaier_add(
      1.0L / static_cast<long double>(radius),
      inverse_total, inverse_correction
    );
  }
  const double inverse_radius_sum = lwz_checked_double(
    inverse_total + inverse_correction,
    "inverse-radius update denominator"
  );
  if (!(inverse_radius_sum > 0.0)) {
    Rcpp::stop(
      "Li-Wang-Zou SST group %d has an invalid inverse-radius update "
      "denominator.", group
    );
  }
  const arma::rowvec mean_sign = sign_sum / static_cast<double>(n);
  const arma::vec diagonal_equation =
    static_cast<double>(p) * sign_square_sum / static_cast<double>(n);
  const double location_residual = arma::abs(mean_sign).max();
  const double diagonal_residual = arma::abs(
    diagonal_equation - 1.0
  ).max();
  return LwzScore{
    sign_sum,
    sign_square_sum,
    radii,
    inverse_radius_sum,
    radii.min(),
    location_residual,
    diagonal_residual,
    std::max(location_residual, diagonal_residual)
  };
}

struct LwzFit {
  arma::vec location;
  arma::vec diagonal;
  int iterations = 0;
  bool iteration_stable = false;
  double relative_update = R_PosInf;
  double location_relative_update = R_PosInf;
  double log_diagonal_update = R_PosInf;
  double score_residual = R_PosInf;
  double location_score_residual = R_PosInf;
  double diagonal_score_residual = R_PosInf;
  double minimum_residual_distance = R_PosInf;
};

LwzFit lwz_fit(const arma::mat& data,
               const double tol,
               const int max_iter,
               const double zero_tol,
               const int group) {
  const arma::uword n = data.n_rows;
  const arma::uword p = data.n_cols;
  LwzFit fit;
  fit.location.zeros(p);
  for (arma::uword j = 0; j < p; ++j) {
    long double total = 0.0L;
    long double correction = 0.0L;
    for (arma::uword i = 0; i < n; ++i) {
      lwz_neumaier_add(
        static_cast<long double>(data(i, j)), total, correction
      );
    }
    fit.location(j) = static_cast<double>(
      (total + correction) / static_cast<long double>(n)
    );
  }
  fit.diagonal.zeros(p);
  for (arma::uword j = 0; j < p; ++j) {
    long double total = 0.0L;
    long double correction = 0.0L;
    for (arma::uword i = 0; i < n; ++i) {
      const long double residual =
        static_cast<long double>(data(i, j)) -
        static_cast<long double>(fit.location(j));
      lwz_neumaier_add(residual * residual, total, correction);
    }
    fit.diagonal(j) = lwz_checked_double(
      (total + correction) / static_cast<long double>(n - 1),
      "initial marginal variance"
    );
    if (!(fit.diagonal(j) > 0.0)) {
      Rcpp::stop(
        "Li-Wang-Zou SST requires every group-specific marginal sample "
        "variance to be strictly positive; group %d variable %llu is "
        "degenerate. No ridge is applied.",
        group, static_cast<unsigned long long>(j + 1)
      );
    }
  }
  lwz_normalise_diagonal(fit.diagonal, group);
  const arma::vec initial_diagonal = fit.diagonal;

  for (int iteration = 0; iteration <= max_iter; ++iteration) {
    const LwzScore score = lwz_score(
      data, fit.location, fit.diagonal, zero_tol, group
    );
    fit.minimum_residual_distance = std::min(
      fit.minimum_residual_distance, score.minimum_radius
    );
    fit.location_score_residual = score.location_residual;
    fit.diagonal_score_residual = score.diagonal_residual;
    fit.score_residual = score.score_residual;
    fit.iterations = iteration;
    if (score.score_residual <= tol) {
      if (iteration == 0) {
        fit.relative_update = 0.0;
        fit.location_relative_update = 0.0;
        fit.log_diagonal_update = 0.0;
      }
      fit.iteration_stable = true;
      return fit;
    }
    if (iteration == max_iter) {
      return fit;
    }

    const arma::vec next_location = fit.location +
      arma::sqrt(fit.diagonal) %
        (score.sign_sum.t() / score.inverse_radius_sum);
    arma::vec next_diagonal = fit.diagonal %
      (static_cast<double>(p) * score.sign_square_sum /
       static_cast<double>(n));
    if (!next_location.is_finite() || !next_diagonal.is_finite() ||
        arma::any(next_diagonal <= 0.0)) {
      Rcpp::stop(
        "Li-Wang-Zou SST group %d HR update became non-finite or "
        "non-positive; no ridge, absolute value, floor, or pseudoinverse "
        "is applied.", group
      );
    }
    lwz_normalise_diagonal(next_diagonal, group);
    fit.location_relative_update = arma::abs(
      (next_location - fit.location) / arma::sqrt(initial_diagonal)
    ).max();
    fit.log_diagonal_update = arma::abs(
      arma::log(next_diagonal / fit.diagonal)
    ).max();
    fit.relative_update = std::max(
      fit.location_relative_update, fit.log_diagonal_update
    );
    fit.location = next_location;
    fit.diagonal = next_diagonal;
  }
  return fit;
}

struct LwzOwnSigns {
  arma::mat direction;
  arma::vec radius;
  arma::vec inverse_radius;
  double c_hat;
};

LwzOwnSigns lwz_own_signs(const arma::mat& data,
                          const LwzFit& fit,
                          const double zero_tol,
                          const int group) {
  const arma::uword n = data.n_rows;
  const arma::uword p = data.n_cols;
  arma::mat direction(n, p);
  arma::vec radius(n);
  arma::vec inverse_radius(n);
  long double inverse_total = 0.0L;
  long double inverse_correction = 0.0L;
  for (arma::uword i = 0; i < n; ++i) {
    arma::rowvec current;
    if (!lwz_direction(
          data.row(i), fit.location, fit.diagonal, zero_tol,
          current, radius(i)
        )) {
      Rcpp::stop(
        "Li-Wang-Zou SST group %d inverse-radius estimate is undefined "
        "at observation %llu; no perturbation or weight cap is applied.",
        group, static_cast<unsigned long long>(i + 1)
      );
    }
    direction.row(i) = current;
    inverse_radius(i) = 1.0 / radius(i);
    if (!std::isfinite(inverse_radius(i))) {
      Rcpp::stop(
        "Li-Wang-Zou SST group %d has a non-finite inverse radius; no "
        "weight cap is applied.", group
      );
    }
    lwz_neumaier_add(
      static_cast<long double>(inverse_radius(i)),
      inverse_total, inverse_correction
    );
  }
  return LwzOwnSigns{
    direction,
    radius,
    inverse_radius,
    lwz_checked_double(
      (inverse_total + inverse_correction) /
        static_cast<long double>(n),
      "inverse-radius moment"
    )
  };
}

void lwz_store_diagonal(const arma::vec& diagonal,
                        const LwzPreparedData& prepared,
                        arma::vec& canonical,
                        arma::vec& canonical_log) {
  canonical.set_size(diagonal.n_elem);
  canonical_log.set_size(diagonal.n_elem);
  double maximum_log = -std::numeric_limits<double>::infinity();
  for (arma::uword j = 0; j < diagonal.n_elem; ++j) {
    canonical_log(j) = 2.0 * prepared.log_scale(j) +
      std::log(diagonal(j));
    maximum_log = std::max(maximum_log, canonical_log(j));
  }
  for (arma::uword j = 0; j < diagonal.n_elem; ++j) {
    canonical_log(j) -= maximum_log;
    canonical(j) = std::exp(canonical_log(j));
  }
}

Rcpp::List lwz_fit_diagnostics(const LwzFit& fit) {
  return Rcpp::List::create(
    Rcpp::Named("iterations") = fit.iterations,
    Rcpp::Named("iteration.stable") = fit.iteration_stable,
    Rcpp::Named("relative.update") = fit.relative_update,
    Rcpp::Named("location.relative.update") =
      fit.location_relative_update,
    Rcpp::Named("log.diagonal.relative.update") =
      fit.log_diagonal_update,
    Rcpp::Named("score.residual") = fit.score_residual,
    Rcpp::Named("location.score.residual") =
      fit.location_score_residual,
    Rcpp::Named("diagonal.score.residual") =
      fit.diagonal_score_residual,
    Rcpp::Named("minimum.residual.distance") =
      fit.minimum_residual_distance,
    Rcpp::Named("convergence.basis") =
      "maximum location/diagonal estimating-equation residual"
  );
}

}  // namespace


//' Li--Wang--Zou simpler spatial-sign two-sample kernel
//'
//' @param x First numeric observation-by-variable matrix.
//' @param y Second numeric observation-by-variable matrix.
//' @param tol Positive estimating-equation tolerance.
//' @param max_iter Positive maximum update count for each full-sample fit.
//' @param zero_tol Non-negative singular-radius tolerance.
//' @return Internal list containing the feasible SST and diagnostics.
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List cpp_li_wang_zou_two_sample_sign(
    const arma::mat& x,
    const arma::mat& y,
    const double tol,
    const int max_iter,
    const double zero_tol) {
  if (!x.is_finite() || !y.is_finite()) {
    Rcpp::stop("`x` and `y` must contain only finite values.");
  }
  if (x.n_rows < 2 || y.n_rows < 2) {
    Rcpp::stop(
      "Li-Wang-Zou SST requires at least two observations in each group."
    );
  }
  if (x.n_cols < 1 || x.n_cols != y.n_cols) {
    Rcpp::stop("`x` and `y` must have the same positive number of columns.");
  }
  if (!std::isfinite(tol) || tol <= 0.0 || max_iter < 1 ||
      !std::isfinite(zero_tol) || zero_tol < 0.0) {
    Rcpp::stop("Invalid Li-Wang-Zou SST iteration controls.");
  }

  const arma::uword n1 = x.n_rows;
  const arma::uword n2 = y.n_rows;
  const arma::uword p = x.n_cols;
  const long double n1_ld = static_cast<long double>(n1);
  const long double n2_ld = static_cast<long double>(n2);
  const long double p_ld = static_cast<long double>(p);
  const LwzPreparedData prepared = lwz_prepare_data(x, y);
  const LwzFit fit1 = lwz_fit(
    prepared.x, tol, max_iter, zero_tol, 1
  );
  const LwzFit fit2 = lwz_fit(
    prepared.y, tol, max_iter, zero_tol, 2
  );
  const LwzOwnSigns own1 = lwz_own_signs(
    prepared.x, fit1, zero_tol, 1
  );
  const LwzOwnSigns own2 = lwz_own_signs(
    prepared.y, fit2, zero_tol, 2
  );

  const long double c_ratio1_ld =
    static_cast<long double>(own2.c_hat) /
    static_cast<long double>(own1.c_hat);
  const long double c_ratio2_ld = 1.0L / c_ratio1_ld;
  const double c_ratio1 = lwz_checked_double(
    c_ratio1_ld, "c2/c1 ratio"
  );
  const double c_ratio2 = lwz_checked_double(
    c_ratio2_ld, "c1/c2 ratio"
  );

  arma::vec bridge1(p);
  arma::vec bridge2(p);
  long double trace_bridge1 = 0.0L;
  long double trace_bridge1_correction = 0.0L;
  long double trace_bridge2 = 0.0L;
  long double trace_bridge2_correction = 0.0L;
  for (arma::uword j = 0; j < p; ++j) {
    const long double half_log_ratio = 0.5L * (
      std::log(static_cast<long double>(fit1.diagonal(j))) -
      std::log(static_cast<long double>(fit2.diagonal(j)))
    );
    bridge1(j) = lwz_checked_double(
      std::exp(half_log_ratio), "group-1 diagonal bridge"
    );
    bridge2(j) = lwz_checked_double(
      std::exp(-half_log_ratio), "group-2 diagonal bridge"
    );
    lwz_neumaier_add(
      static_cast<long double>(bridge1(j)),
      trace_bridge1, trace_bridge1_correction
    );
    lwz_neumaier_add(
      static_cast<long double>(bridge2(j)),
      trace_bridge2, trace_bridge2_correction
    );
  }
  const long double trace_a1_ld = c_ratio1_ld *
    (trace_bridge1 + trace_bridge1_correction);
  const long double trace_a2_ld = c_ratio2_ld *
    (trace_bridge2 + trace_bridge2_correction);

  arma::mat trace_a1_pairs(n1, n1, arma::fill::zeros);
  arma::mat trace_a2_pairs(n2, n2, arma::fill::zeros);
  arma::mat trace_a3_pairs(n1, n2);
  long double trace_a1_pair_sum = 0.0L;
  long double trace_a1_pair_correction = 0.0L;
  long double trace_a2_pair_sum = 0.0L;
  long double trace_a2_pair_correction = 0.0L;
  long double trace_a3_pair_sum = 0.0L;
  long double trace_a3_pair_correction = 0.0L;
  for (arma::uword k = 0; k < n1; ++k) {
    for (arma::uword ell = 0; ell < n1; ++ell) {
      if (ell == k) {
        continue;
      }
      const long double inner = lwz_bridge_dot(
        own1.direction.row(ell), bridge1, own1.direction.row(k)
      );
      const long double square = inner * inner;
      trace_a1_pairs(ell, k) = lwz_checked_double(
        square, "group-1 ordered trace kernel"
      );
      lwz_neumaier_add(
        square, trace_a1_pair_sum, trace_a1_pair_correction
      );
    }
  }
  for (arma::uword k = 0; k < n2; ++k) {
    for (arma::uword ell = 0; ell < n2; ++ell) {
      if (ell == k) {
        continue;
      }
      const long double inner = lwz_bridge_dot(
        own2.direction.row(ell), bridge2, own2.direction.row(k)
      );
      const long double square = inner * inner;
      trace_a2_pairs(ell, k) = lwz_checked_double(
        square, "group-2 ordered trace kernel"
      );
      lwz_neumaier_add(
        square, trace_a2_pair_sum, trace_a2_pair_correction
      );
    }
  }
  for (arma::uword ell = 0; ell < n1; ++ell) {
    for (arma::uword k = 0; k < n2; ++k) {
      const long double inner = lwz_dot(
        own1.direction.row(ell), own2.direction.row(k)
      );
      const long double square = inner * inner;
      trace_a3_pairs(ell, k) = lwz_checked_double(
        square, "cross trace kernel"
      );
      lwz_neumaier_add(
        square, trace_a3_pair_sum, trace_a3_pair_correction
      );
    }
  }

  const long double p_square = p_ld * p_ld;
  const long double trace_a1_sq_ld = p_square * c_ratio1_ld *
    c_ratio1_ld *
    (trace_a1_pair_sum + trace_a1_pair_correction) /
    (n1_ld * static_cast<long double>(n1 - 1));
  const long double trace_a2_sq_ld = p_square * c_ratio2_ld *
    c_ratio2_ld *
    (trace_a2_pair_sum + trace_a2_pair_correction) /
    (n2_ld * static_cast<long double>(n2 - 1));
  const long double trace_a3_ld = p_square *
    (trace_a3_pair_sum + trace_a3_pair_correction) /
    (n1_ld * n2_ld);

  arma::mat numerator_inner(n1, n2);
  arma::mat numerator_contribution(n1, n2);
  long double statistic_sum = 0.0L;
  long double statistic_correction = 0.0L;
  arma::uword numerator_zero_signs = 0;
  for (arma::uword i = 0; i < n1; ++i) {
    for (arma::uword j = 0; j < n2; ++j) {
      arma::rowvec direction1;
      arma::rowvec direction2;
      double radius1 = 0.0;
      double radius2 = 0.0;
      const bool nonzero1 = lwz_direction(
        prepared.x.row(i), fit2.location, fit1.diagonal,
        zero_tol, direction1, radius1
      );
      const bool nonzero2 = lwz_direction(
        prepared.y.row(j), fit1.location, fit2.diagonal,
        zero_tol, direction2, radius2
      );
      if (!nonzero1) {
        direction1.zeros(p);
        ++numerator_zero_signs;
      }
      if (!nonzero2) {
        direction2.zeros(p);
        ++numerator_zero_signs;
      }
      const long double inner = lwz_dot(direction1, direction2);
      const long double contribution = -inner;
      numerator_inner(i, j) = lwz_checked_double(
        inner, "numerator inner product"
      );
      numerator_contribution(i, j) = lwz_checked_double(
        contribution, "numerator contribution"
      );
      lwz_neumaier_add(
        contribution, statistic_sum, statistic_correction
      );
    }
  }
  const long double statistic_ld =
    (statistic_sum + statistic_correction) / (n1_ld * n2_ld);
  const long double bias_ld = trace_a1_ld / (n1_ld * p_ld) +
    trace_a2_ld / (n2_ld * p_ld);
  const long double centered_ld = statistic_ld - bias_ld;
  // The n_k^2 factors below are the literal feasible variance multipliers
  // in Li--Wang--Zou.  Only the trace estimators use n_k(n_k-1).
  const long double variance_term1_ld = 2.0L * trace_a1_sq_ld /
    (n1_ld * n1_ld * p_square);
  const long double variance_term2_ld = 2.0L * trace_a2_sq_ld /
    (n2_ld * n2_ld * p_square);
  const long double variance_term3_ld = 4.0L * trace_a3_ld /
    (n1_ld * n2_ld * p_square);
  const long double variance_ld = variance_term1_ld +
    variance_term2_ld + variance_term3_ld;
  const double statistic = lwz_checked_double(
    statistic_ld, "uncorrected statistic"
  );
  const double bias = lwz_checked_double(bias_ld, "bias estimate");
  const double centered = lwz_checked_double(
    centered_ld, "bias-corrected statistic"
  );
  const double variance = lwz_checked_double(
    variance_ld, "feasible variance"
  );
  if (!(variance > 0.0)) {
    Rcpp::stop(
      "Li-Wang-Zou SST requires its feasible variance to be strictly "
      "positive. It was %.17g; no absolute value, floor, ridge, or "
      "variance substitution is applied.", variance
    );
  }
  const double standard_error = std::sqrt(variance);
  const double z = centered / standard_error;
  if (!std::isfinite(z)) {
    Rcpp::stop("Li-Wang-Zou SST produced a non-finite Z statistic.");
  }

  arma::vec diagonal1_input;
  arma::vec log_diagonal1_input;
  arma::vec diagonal2_input;
  arma::vec log_diagonal2_input;
  lwz_store_diagonal(
    fit1.diagonal, prepared, diagonal1_input, log_diagonal1_input
  );
  lwz_store_diagonal(
    fit2.diagonal, prepared, diagonal2_input, log_diagonal2_input
  );
  const arma::vec sample_mean1 = lwz_restore_location(
    arma::mean(prepared.x, 0).t(), prepared
  );
  const arma::vec sample_mean2 = lwz_restore_location(
    arma::mean(prepared.y, 0).t(), prepared
  );
  const Rcpp::List diagnostics1 = lwz_fit_diagnostics(fit1);
  const Rcpp::List diagnostics2 = lwz_fit_diagnostics(fit2);

  return Rcpp::List::create(
    Rcpp::Named("z") = z,
    Rcpp::Named("T") = statistic,
    Rcpp::Named("bias") = bias,
    Rcpp::Named("centered") = centered,
    Rcpp::Named("variance") = variance,
    Rcpp::Named("standard_error") = standard_error,
    Rcpp::Named("c1_hat") = own1.c_hat,
    Rcpp::Named("c2_hat") = own2.c_hat,
    Rcpp::Named("c_ratio1_hat") = c_ratio1,
    Rcpp::Named("c_ratio2_hat") = c_ratio2,
    Rcpp::Named("trace_A1_hat") = lwz_checked_double(
      trace_a1_ld, "trace A1 estimate"
    ),
    Rcpp::Named("trace_A2_hat") = lwz_checked_double(
      trace_a2_ld, "trace A2 estimate"
    ),
    Rcpp::Named("trace_A1_squared_hat") = lwz_checked_double(
      trace_a1_sq_ld, "trace A1 squared estimate"
    ),
    Rcpp::Named("trace_A2_squared_hat") = lwz_checked_double(
      trace_a2_sq_ld, "trace A2 squared estimate"
    ),
    Rcpp::Named("trace_A3tA3_hat") = lwz_checked_double(
      trace_a3_ld, "trace A3 transpose A3 estimate"
    ),
    Rcpp::Named("variance_term1") = lwz_checked_double(
      variance_term1_ld, "first variance term"
    ),
    Rcpp::Named("variance_term2") = lwz_checked_double(
      variance_term2_ld, "second variance term"
    ),
    Rcpp::Named("variance_term3") = lwz_checked_double(
      variance_term3_ld, "cross variance term"
    ),
    Rcpp::Named("numerator_inner_product") = numerator_inner,
    Rcpp::Named("numerator_contribution") = numerator_contribution,
    Rcpp::Named("trace_A1_ordered_pair_squared") = trace_a1_pairs,
    Rcpp::Named("trace_A2_ordered_pair_squared") = trace_a2_pairs,
    Rcpp::Named("trace_A3_cross_pair_squared") = trace_a3_pairs,
    Rcpp::Named("bridge1_diagonal_standardized") = bridge1,
    Rcpp::Named("bridge2_diagonal_standardized") = bridge2,
    Rcpp::Named("own_direction1") = own1.direction,
    Rcpp::Named("own_direction2") = own2.direction,
    Rcpp::Named("own_radius1") = own1.radius,
    Rcpp::Named("own_radius2") = own2.radius,
    Rcpp::Named("own_inverse_radius1") = own1.inverse_radius,
    Rcpp::Named("own_inverse_radius2") = own2.inverse_radius,
    Rcpp::Named("location1_standardized") = fit1.location,
    Rcpp::Named("location2_standardized") = fit2.location,
    Rcpp::Named("location1") = lwz_restore_location(
      fit1.location, prepared
    ),
    Rcpp::Named("location2") = lwz_restore_location(
      fit2.location, prepared
    ),
    Rcpp::Named("diagonal1_standardized") = fit1.diagonal,
    Rcpp::Named("diagonal2_standardized") = fit2.diagonal,
    Rcpp::Named("diagonal1_input_canonical") = diagonal1_input,
    Rcpp::Named("diagonal2_input_canonical") = diagonal2_input,
    Rcpp::Named("log_diagonal1_input_canonical") =
      log_diagonal1_input,
    Rcpp::Named("log_diagonal2_input_canonical") =
      log_diagonal2_input,
    Rcpp::Named("fit_diagnostics1") = diagnostics1,
    Rcpp::Named("fit_diagnostics2") = diagnostics2,
    Rcpp::Named("all_iteration_stable") =
      fit1.iteration_stable && fit2.iteration_stable,
    Rcpp::Named("stability_failures") =
      static_cast<int>(!fit1.iteration_stable) +
      static_cast<int>(!fit2.iteration_stable),
    Rcpp::Named("numerator_zero_sign_uses") =
      static_cast<double>(numerator_zero_signs),
    Rcpp::Named("sample_mean1") = sample_mean1,
    Rcpp::Named("sample_mean2") = sample_mean2,
    Rcpp::Named("column_log_scale") = prepared.log_scale,
    Rcpp::Named("subtraction_overflow_columns") = static_cast<double>(
      arma::accu(prepared.overflow_fallback)
    ),
    Rcpp::Named("ordered_pair_count1") = static_cast<double>(
      n1 * (n1 - 1)
    ),
    Rcpp::Named("ordered_pair_count2") = static_cast<double>(
      n2 * (n2 - 1)
    ),
    Rcpp::Named("cross_pair_count") = static_cast<double>(n1 * n2),
    Rcpp::Named("n1") = static_cast<double>(n1),
    Rcpp::Named("n2") = static_cast<double>(n2),
    Rcpp::Named("p") = static_cast<double>(p),
    Rcpp::Named("tol") = tol,
    Rcpp::Named("max_iter") = max_iter,
    Rcpp::Named("zero_tol") = zero_tol
  );
}
