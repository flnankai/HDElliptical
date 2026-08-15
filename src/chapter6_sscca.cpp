// Chapter 6 spatial-sign sparse canonical correlation analysis.
//
// The primary SSCCA kernel alternates exact convex metric-lasso block
// problems.  Unlike the mixedCCA reference implementation, every coordinate
// update divides by the corresponding metric diagonal, so the update remains
// valid when the spatial-sign blocks are not correlation matrices.

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp17)]]

#include <RcppArmadillo.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <string>
#include <vector>

namespace {

double ch6ss_soft(const double value, const double threshold) {
  if (value > threshold) return value - threshold;
  if (value < -threshold) return value + threshold;
  return 0.0;
}

double ch6ss_max_abs(const arma::vec& value) {
  return value.n_elem ? arma::abs(value).max() : 0.0;
}

double ch6ss_lasso_kkt(const arma::mat& metric, const arma::vec& target,
                       const arma::vec& coefficient, const double lambda,
                       const double support_tolerance) {
  const arma::vec gradient = metric * coefficient - target;
  double residual = 0.0;
  for (arma::uword j = 0; j < coefficient.n_elem; ++j) {
    double component;
    if (std::abs(coefficient(j)) > support_tolerance) {
      component = std::abs(
        gradient(j) + lambda * (coefficient(j) > 0.0 ? 1.0 : -1.0)
      );
    } else {
      component = std::max(0.0, std::abs(gradient(j)) - lambda);
    }
    residual = std::max(residual, component);
  }
  return residual;
}

struct Ch6ssLassoFit {
  arma::vec coefficient;
  int iterations;
  bool converged;
  double relative_update;
  double kkt;
  double objective;
};

Ch6ssLassoFit ch6ss_metric_lasso(
    const arma::mat& metric, const arma::vec& target, arma::vec coefficient,
    const double lambda, const double tolerance, const int max_iterations) {
  arma::vec residual = target - metric * coefficient;
  double relative_update = std::numeric_limits<double>::infinity();
  double kkt = std::numeric_limits<double>::infinity();
  bool converged = false;
  int iteration = 0;
  const double kkt_scale = 1.0 + ch6ss_max_abs(target) + lambda;

  for (iteration = 1; iteration <= max_iterations; ++iteration) {
    const arma::vec old = coefficient;
    for (arma::uword j = 0; j < coefficient.n_elem; ++j) {
      const double diagonal = metric(j, j);
      const double partial = residual(j) + diagonal * coefficient(j);
      const double updated = ch6ss_soft(partial, lambda) / diagonal;
      const double change = coefficient(j) - updated;
      if (change != 0.0) {
        coefficient(j) = updated;
        residual += metric.col(j) * change;
      }
    }
    relative_update = ch6ss_max_abs(coefficient - old) /
      std::max({1.0, ch6ss_max_abs(old), ch6ss_max_abs(coefficient)});
    // Solver certification uses exact zeros; support_tol affects only BIC df.
    kkt = ch6ss_lasso_kkt(metric, target, coefficient, lambda, 0.0);
    if (relative_update <= tolerance && kkt <= tolerance * kkt_scale) {
      converged = true;
      break;
    }
  }

  const double objective = 0.5 * arma::as_scalar(
    coefficient.t() * metric * coefficient
  ) - arma::dot(target, coefficient) + lambda * arma::accu(arma::abs(coefficient));
  return {coefficient, iteration, converged, relative_update, kkt, objective};
}

struct Ch6ssBlockFit {
  arma::vec raw;
  arma::vec normalized;
  arma::vec bic_values;
  arma::vec rss_values;
  arma::vec kkt_values;
  arma::ivec iteration_values;
  arma::uvec converged_values;
  double lambda;
  double metric_scale;
  double raw_metric_squared;
  double kkt;
  double relative_update;
  int iterations;
  int selected_index;
  int degrees_freedom;
  bool converged;
  bool valid;
  std::string message;
};

Ch6ssBlockFit ch6ss_select_block(
    const arma::mat& metric, const arma::vec& target,
    const arma::vec& start, const arma::vec& lambdas,
    const int selection, const int sample_size, const double tolerance,
    const int max_iterations, const double support_tolerance) {
  const arma::uword count = lambdas.n_elem;
  arma::mat candidates(metric.n_rows, count, arma::fill::zeros);
  arma::vec bic(count, arma::fill::value(
    std::numeric_limits<double>::infinity()
  ));
  arma::vec rss(count, arma::fill::value(
    std::numeric_limits<double>::quiet_NaN()
  ));
  arma::vec kkts(count, arma::fill::value(
    std::numeric_limits<double>::infinity()
  ));
  arma::ivec iterations(count, arma::fill::zeros);
  arma::uvec converged(count, arma::fill::zeros);
  arma::vec warm = start;

  for (arma::uword index = 0; index < count; ++index) {
    Ch6ssLassoFit fit = ch6ss_metric_lasso(
      metric, target, warm, lambdas(index), tolerance, max_iterations
    );
    candidates.col(index) = fit.coefficient;
    warm = fit.coefficient;
    kkts(index) = fit.kkt;
    iterations(index) = fit.iterations;
    converged(index) = fit.converged ? 1u : 0u;
    if (!fit.converged || !fit.coefficient.is_finite()) continue;

    const double metric_squared = arma::as_scalar(
      fit.coefficient.t() * metric * fit.coefficient
    );
    if (!(metric_squared > 0.0) || !std::isfinite(metric_squared)) continue;
    const double residual_sum = metric_squared -
      2.0 * arma::dot(fit.coefficient, target) + 1.0;
    rss(index) = residual_sum;
    if (residual_sum < 0.0 || !std::isfinite(residual_sum)) continue;
    const int degrees = static_cast<int>(arma::accu(
      arma::abs(fit.coefficient) > support_tolerance
    ));

    if (selection == 0) {
      bic(index) = 0.0;
    } else if (selection == 1) {
      bic(index) = residual_sum +
        static_cast<double>(degrees) * std::log(sample_size) / sample_size;
    } else {
      if (degrees >= sample_size) continue;
      if (residual_sum == 0.0) {
        bic(index) = -std::numeric_limits<double>::infinity();
      } else {
        bic(index) = std::log(
          sample_size * residual_sum /
          static_cast<double>(sample_size - degrees)
        ) + static_cast<double>(degrees) * std::log(sample_size) / sample_size;
      }
    }
  }

  int selected = -1;
  double best = std::numeric_limits<double>::infinity();
  for (arma::uword index = 0; index < count; ++index) {
    if (std::isnan(bic(index)) || bic(index) ==
        std::numeric_limits<double>::infinity()) continue;
    if (selected < 0 || bic(index) < best) {
      selected = static_cast<int>(index);
      best = bic(index);
    }
  }

  if (selected < 0) {
    return {
      arma::vec(), arma::vec(), bic, rss, kkts, iterations, converged,
      NA_REAL, NA_REAL, NA_REAL, NA_REAL, NA_REAL, 0, -1, 0, false, false,
      "No nonzero, converged metric-lasso candidate has a valid calibration."
    };
  }

  arma::vec raw = candidates.col(selected);
  const double metric_squared = arma::as_scalar(raw.t() * metric * raw);
  if (!(metric_squared > 0.0) || !std::isfinite(metric_squared)) {
    return {
      raw, arma::vec(), bic, rss, kkts, iterations, converged,
      lambdas(selected), NA_REAL, metric_squared, kkts(selected), NA_REAL,
      iterations(selected), selected, 0, false, false,
      "The selected metric-lasso solution has zero or non-finite metric norm."
    };
  }
  const double metric_scale = std::sqrt(metric_squared);
  arma::vec normalized = raw / metric_scale;
  const int degrees = static_cast<int>(arma::accu(
    arma::abs(raw) > support_tolerance
  ));
  return {
    raw, normalized, bic, rss, kkts, iterations, converged,
    lambdas(selected), metric_scale, metric_squared, kkts(selected),
    NA_REAL, iterations(selected), selected, degrees,
    converged(selected) == 1u, true, ""
  };
}

arma::vec ch6ss_anchor_initial(const arma::mat& cross, const bool left) {
  arma::mat u;
  arma::mat v;
  arma::vec singular;
  if (!arma::svd_econ(u, singular, v, cross, "both", "std") ||
      singular.n_elem == 0 || !(singular(0) > 0.0)) {
    return arma::vec(left ? cross.n_rows : cross.n_cols, arma::fill::zeros);
  }
  return left ? u.col(0) : v.col(0);
}

double ch6ss_relative_direction(const arma::vec& current,
                                const arma::vec& previous) {
  return ch6ss_max_abs(current - previous) /
    std::max({1.0, ch6ss_max_abs(current), ch6ss_max_abs(previous)});
}

struct Ch6ssPmdUpdate {
  arma::vec direction;
  double threshold;
  double l1;
  double l2;
  bool top_tie_slack;
  bool valid;
  std::string message;
};

Ch6ssPmdUpdate ch6ss_pmd_update(const arma::vec& target,
                                const double l1_bound,
                                const double tolerance) {
  const arma::vec magnitude = arma::abs(target);
  const double maximum = magnitude.max();
  if (!(maximum > 0.0) || !std::isfinite(maximum)) {
    return {arma::vec(), NA_REAL, NA_REAL, NA_REAL, false, false,
            "A PMD block target is zero or non-finite."};
  }
  const double norm_two = arma::norm(target, 2);
  arma::vec unconstrained = target / norm_two;
  const double unconstrained_l1 = arma::accu(arma::abs(unconstrained));
  if (unconstrained_l1 <= l1_bound + tolerance) {
    return {unconstrained, 0.0, unconstrained_l1, 1.0, false, true, ""};
  }

  const double tie_tolerance = tolerance * std::max(1.0, maximum);
  const arma::uvec top = arma::find(magnitude >= maximum - tie_tolerance);
  const double root_top = std::sqrt(static_cast<double>(top.n_elem));
  if (l1_bound < root_top - tolerance) {
    arma::vec answer(target.n_elem, arma::fill::zeros);
    const double size = l1_bound / static_cast<double>(top.n_elem);
    for (arma::uword j = 0; j < top.n_elem; ++j) {
      answer(top(j)) = target(top(j)) > 0.0 ? size : -size;
    }
    return {answer, maximum, l1_bound, arma::norm(answer, 2), true, true, ""};
  }
  if (l1_bound <= 1.0 + tolerance) {
    arma::vec answer(target.n_elem, arma::fill::zeros);
    const arma::uword index = magnitude.index_max();
    answer(index) = target(index) > 0.0 ? 1.0 : -1.0;
    return {answer, maximum, 1.0, 1.0, top.n_elem > 1u, true, ""};
  }

  double lower = 0.0;
  double upper = maximum;
  arma::vec best;
  for (int iteration = 0; iteration < 240; ++iteration) {
    const double threshold = lower + 0.5 * (upper - lower);
    arma::vec soft = arma::sign(target) % arma::clamp(
      magnitude - threshold, 0.0, std::numeric_limits<double>::infinity()
    );
    const double norm = arma::norm(soft, 2);
    if (!(norm > 0.0)) {
      upper = threshold;
      continue;
    }
    soft /= norm;
    const double ratio = arma::accu(arma::abs(soft));
    if (ratio > l1_bound) {
      lower = threshold;
    } else {
      upper = threshold;
      best = soft;
    }
    if (upper - lower <= tolerance * std::max(1.0, maximum)) break;
  }
  if (best.n_elem == 0 || !best.is_finite()) {
    return {arma::vec(), NA_REAL, NA_REAL, NA_REAL, false, false,
            "The PMD soft-threshold root could not be resolved."};
  }
  return {
    best, upper, arma::accu(arma::abs(best)), arma::norm(best, 2),
    false, true, ""
  };
}

}  // namespace


// [[Rcpp::export]]
Rcpp::List cpp_ch6_sscca_primary(
    const arma::mat& metric_x, const arma::mat& metric_y,
    const arma::mat& cross, const arma::vec& lambda_x,
    const arma::vec& lambda_y, const int selection, const int sample_size,
    const double tolerance, const int max_iterations,
    const double inner_tolerance, const int inner_max_iterations,
    const double support_tolerance) {
  if (!metric_x.is_finite() || !metric_y.is_finite() || !cross.is_finite()) {
    Rcpp::stop("SSCCA matrices must be finite.");
  }
  if (metric_x.n_rows != metric_x.n_cols ||
      metric_y.n_rows != metric_y.n_cols ||
      cross.n_rows != metric_x.n_rows || cross.n_cols != metric_y.n_rows) {
    Rcpp::stop("SSCCA matrix dimensions are inconsistent.");
  }
  if (lambda_x.n_elem == 0 || lambda_y.n_elem == 0 ||
      !lambda_x.is_finite() || !lambda_y.is_finite() ||
      arma::any(lambda_x < 0.0) || arma::any(lambda_y < 0.0)) {
    Rcpp::stop("SSCCA lambda sequences must be finite and non-negative.");
  }
  if (selection < 0 || selection > 2 || sample_size < 2 ||
      !(tolerance > 0.0) || max_iterations < 1 ||
      !(inner_tolerance > 0.0) || inner_max_iterations < 1 ||
      support_tolerance < 0.0) {
    Rcpp::stop("Invalid SSCCA solver controls.");
  }
  if (arma::any(metric_x.diag() <= 0.0) ||
      arma::any(metric_y.diag() <= 0.0)) {
    Rcpp::stop("Every SSCCA marginal metric diagonal must be positive.");
  }

  arma::vec w_x = ch6ss_anchor_initial(cross, true);
  arma::vec w_y = ch6ss_anchor_initial(cross, false);
  const double mx = arma::as_scalar(w_x.t() * metric_x * w_x);
  const double my = arma::as_scalar(w_y.t() * metric_y * w_y);
  if (!(mx > 0.0) || !(my > 0.0) || !std::isfinite(mx) ||
      !std::isfinite(my)) {
    return Rcpp::List::create(
      Rcpp::Named("valid") = false,
      Rcpp::Named("converged") = false,
      Rcpp::Named("message") =
        "The deterministic leading-cross-SVD initialization has zero metric norm."
    );
  }
  w_x /= std::sqrt(mx);
  w_y /= std::sqrt(my);

  arma::vec objective_trace(max_iterations, arma::fill::zeros);
  arma::vec update_trace(max_iterations, arma::fill::zeros);
  arma::vec lambda_x_trace(max_iterations, arma::fill::zeros);
  arma::vec lambda_y_trace(max_iterations, arma::fill::zeros);
  Ch6ssBlockFit fit_x;
  Ch6ssBlockFit fit_y;
  bool valid = true;
  bool converged = false;
  std::string message;
  int completed = 0;

  for (int iteration = 0; iteration < max_iterations; ++iteration) {
    const arma::vec old_x = w_x;
    const arma::vec old_y = w_y;
    fit_x = ch6ss_select_block(
      metric_x, cross * w_y, w_x, lambda_x, selection, sample_size,
      inner_tolerance, inner_max_iterations, support_tolerance
    );
    if (!fit_x.valid) {
      valid = false;
      message = std::string("x block: ") + fit_x.message;
      break;
    }
    w_x = fit_x.normalized;
    fit_y = ch6ss_select_block(
      metric_y, cross.t() * w_x, w_y, lambda_y, selection, sample_size,
      inner_tolerance, inner_max_iterations, support_tolerance
    );
    if (!fit_y.valid) {
      valid = false;
      message = std::string("y block: ") + fit_y.message;
      break;
    }
    w_y = fit_y.normalized;

    const double association = arma::dot(w_x, cross * w_y);
    const double objective = association - fit_x.lambda *
      arma::accu(arma::abs(w_x)) - fit_y.lambda *
      arma::accu(arma::abs(w_y));
    const double update = std::max(
      ch6ss_relative_direction(w_x, old_x),
      ch6ss_relative_direction(w_y, old_y)
    );
    objective_trace(iteration) = objective;
    update_trace(iteration) = update;
    lambda_x_trace(iteration) = fit_x.lambda;
    lambda_y_trace(iteration) = fit_y.lambda;
    completed = iteration + 1;
    if (update <= tolerance) {
      converged = true;
      break;
    }
  }

  if (!valid) {
    return Rcpp::List::create(
      Rcpp::Named("valid") = false,
      Rcpp::Named("converged") = false,
      Rcpp::Named("message") = message,
      Rcpp::Named("iterations") = completed,
      Rcpp::Named("objective_trace") = objective_trace.head(completed),
      Rcpp::Named("update_trace") = update_trace.head(completed)
    );
  }

  const double association = arma::dot(w_x, cross * w_y);
  arma::uword anchor = arma::abs(w_x).index_max();
  if (w_x(anchor) < 0.0) {
    w_x = -w_x;
    w_y = -w_y;
    fit_x.raw = -fit_x.raw;
    fit_y.raw = -fit_y.raw;
  }
  const double pair_kkt_x = ch6ss_lasso_kkt(
    metric_x, cross * w_y, fit_x.raw, fit_x.lambda, 0.0
  );
  const double pair_kkt_y = ch6ss_lasso_kkt(
    metric_y, cross.t() * w_x, fit_y.raw, fit_y.lambda, 0.0
  );

  return Rcpp::List::create(
    Rcpp::Named("valid") = true,
    Rcpp::Named("converged") = converged,
    Rcpp::Named("message") = converged ? "" :
      "The alternating SSCCA iteration did not stabilize in max_iter steps.",
    Rcpp::Named("x") = w_x,
    Rcpp::Named("y") = w_y,
    Rcpp::Named("raw_x") = fit_x.raw,
    Rcpp::Named("raw_y") = fit_y.raw,
    Rcpp::Named("association") = association,
    Rcpp::Named("objective") = completed ? objective_trace(completed - 1) : NA_REAL,
    Rcpp::Named("lambda_x") = fit_x.lambda,
    Rcpp::Named("lambda_y") = fit_y.lambda,
    Rcpp::Named("selected_index_x") = fit_x.selected_index + 1,
    Rcpp::Named("selected_index_y") = fit_y.selected_index + 1,
    Rcpp::Named("degrees_freedom_x") = fit_x.degrees_freedom,
    Rcpp::Named("degrees_freedom_y") = fit_y.degrees_freedom,
    Rcpp::Named("metric_norm_x") = std::sqrt(
      arma::as_scalar(w_x.t() * metric_x * w_x)
    ),
    Rcpp::Named("metric_norm_y") = std::sqrt(
      arma::as_scalar(w_y.t() * metric_y * w_y)
    ),
    Rcpp::Named("inner_kkt_x") = fit_x.kkt,
    Rcpp::Named("inner_kkt_y") = fit_y.kkt,
    Rcpp::Named("pair_kkt_x") = pair_kkt_x,
    Rcpp::Named("pair_kkt_y") = pair_kkt_y,
    Rcpp::Named("inner_converged_x") = fit_x.converged,
    Rcpp::Named("inner_converged_y") = fit_y.converged,
    Rcpp::Named("inner_iterations_x") = fit_x.iterations,
    Rcpp::Named("inner_iterations_y") = fit_y.iterations,
    Rcpp::Named("raw_metric_scale_x") = fit_x.metric_scale,
    Rcpp::Named("raw_metric_scale_y") = fit_y.metric_scale,
    Rcpp::Named("bic_values_x") = fit_x.bic_values,
    Rcpp::Named("bic_values_y") = fit_y.bic_values,
    Rcpp::Named("rss_values_x") = fit_x.rss_values,
    Rcpp::Named("rss_values_y") = fit_y.rss_values,
    Rcpp::Named("candidate_kkt_x") = fit_x.kkt_values,
    Rcpp::Named("candidate_kkt_y") = fit_y.kkt_values,
    Rcpp::Named("candidate_iterations_x") = fit_x.iteration_values,
    Rcpp::Named("candidate_iterations_y") = fit_y.iteration_values,
    Rcpp::Named("candidate_converged_x") = fit_x.converged_values,
    Rcpp::Named("candidate_converged_y") = fit_y.converged_values,
    Rcpp::Named("iterations") = completed,
    Rcpp::Named("relative_update") = completed ? update_trace(completed - 1) : NA_REAL,
    Rcpp::Named("objective_trace") = objective_trace.head(completed),
    Rcpp::Named("update_trace") = update_trace.head(completed),
    Rcpp::Named("lambda_x_trace") = lambda_x_trace.head(completed),
    Rcpp::Named("lambda_y_trace") = lambda_y_trace.head(completed),
    Rcpp::Named("sign_anchor_x") = static_cast<int>(anchor + 1)
  );
}


// [[Rcpp::export]]
Rcpp::List cpp_ch6_sign_whitened_pmd(
    const arma::mat& cross, const double l1_x, const double l1_y,
    const double tolerance, const int max_iterations) {
  if (!cross.is_finite() || cross.n_rows == 0 || cross.n_cols == 0 ||
      !(l1_x >= 1.0) || !(l1_y >= 1.0) || !(tolerance > 0.0) ||
      max_iterations < 1) {
    Rcpp::stop("Invalid sign-whitened PMD inputs.");
  }
  arma::vec x = ch6ss_anchor_initial(cross, true);
  arma::vec y = ch6ss_anchor_initial(cross, false);
  if (arma::norm(x, 2) == 0.0 || arma::norm(y, 2) == 0.0) {
    return Rcpp::List::create(
      Rcpp::Named("valid") = false,
      Rcpp::Named("converged") = false,
      Rcpp::Named("message") = "The whitened cross block is zero."
    );
  }
  arma::vec objective_trace(max_iterations, arma::fill::zeros);
  arma::vec update_trace(max_iterations, arma::fill::zeros);
  Ch6ssPmdUpdate update_x;
  Ch6ssPmdUpdate update_y;
  bool converged = false;
  int completed = 0;

  for (int iteration = 0; iteration < max_iterations; ++iteration) {
    const arma::vec old_x = x;
    const arma::vec old_y = y;
    update_x = ch6ss_pmd_update(cross * y, l1_x, tolerance);
    if (!update_x.valid) {
      return Rcpp::List::create(
        Rcpp::Named("valid") = false,
        Rcpp::Named("converged") = false,
        Rcpp::Named("message") = std::string("x block: ") + update_x.message
      );
    }
    x = update_x.direction;
    update_y = ch6ss_pmd_update(cross.t() * x, l1_y, tolerance);
    if (!update_y.valid) {
      return Rcpp::List::create(
        Rcpp::Named("valid") = false,
        Rcpp::Named("converged") = false,
        Rcpp::Named("message") = std::string("y block: ") + update_y.message
      );
    }
    y = update_y.direction;
    const double objective = arma::dot(x, cross * y);
    const double update = std::max(
      ch6ss_relative_direction(x, old_x),
      ch6ss_relative_direction(y, old_y)
    );
    objective_trace(iteration) = objective;
    update_trace(iteration) = update;
    completed = iteration + 1;
    if (update <= tolerance) {
      converged = true;
      break;
    }
  }

  arma::uword anchor = arma::abs(x).index_max();
  if (x(anchor) < 0.0) {
    x = -x;
    y = -y;
  }
  return Rcpp::List::create(
    Rcpp::Named("valid") = true,
    Rcpp::Named("converged") = converged,
    Rcpp::Named("message") = converged ? "" :
      "The alternating sign-whitened PMD iteration did not stabilize in max_iter steps.",
    Rcpp::Named("x") = x,
    Rcpp::Named("y") = y,
    Rcpp::Named("association") = arma::dot(x, cross * y),
    Rcpp::Named("l1_x") = arma::accu(arma::abs(x)),
    Rcpp::Named("l1_y") = arma::accu(arma::abs(y)),
    Rcpp::Named("l2_x") = arma::norm(x, 2),
    Rcpp::Named("l2_y") = arma::norm(y, 2),
    Rcpp::Named("threshold_x") = update_x.threshold,
    Rcpp::Named("threshold_y") = update_y.threshold,
    Rcpp::Named("top_tie_slack_x") = update_x.top_tie_slack,
    Rcpp::Named("top_tie_slack_y") = update_y.top_tie_slack,
    Rcpp::Named("constraint_violation_x") = std::max({
      0.0, arma::norm(x, 2) - 1.0, arma::accu(arma::abs(x)) - l1_x
    }),
    Rcpp::Named("constraint_violation_y") = std::max({
      0.0, arma::norm(y, 2) - 1.0, arma::accu(arma::abs(y)) - l1_y
    }),
    Rcpp::Named("iterations") = completed,
    Rcpp::Named("relative_update") = completed ? update_trace(completed - 1) : NA_REAL,
    Rcpp::Named("objective_trace") = objective_trace.head(completed),
    Rcpp::Named("update_trace") = update_trace.head(completed),
    Rcpp::Named("sign_anchor_x") = static_cast<int>(anchor + 1)
  );
}
