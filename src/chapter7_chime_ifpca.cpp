// Chapter 7 CHIME and influential-features PCA kernels.
//
// The CHIME implementation below solves the printed convex quadratic-lasso
// update directly and certifies its KKT equations.  It intentionally contains
// no covariance ridge, inverse, condition-number repair, or label-based tuning.

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp17)]]

#include <RcppArmadillo.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <string>
#include <vector>

namespace {

double c7_soft(const double value, const double penalty) {
  if (value > penalty) return value - penalty;
  if (value < -penalty) return value + penalty;
  return 0.0;
}

double c7_max_abs(const arma::vec& value) {
  return value.n_elem ? arma::abs(value).max() : 0.0;
}

double c7_max_abs_matrix(const arma::mat& value) {
  return value.n_elem ? arma::abs(value).max() : 0.0;
}

double c7_chime_kkt(const arma::mat& covariance, const arma::vec& difference,
                    const arma::vec& coefficient, const double lambda) {
  const arma::vec gradient = covariance * coefficient - difference;
  double residual = 0.0;
  for (arma::uword j = 0; j < coefficient.n_elem; ++j) {
    double component;
    if (coefficient(j) != 0.0) {
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

struct C7ChimeLasso {
  arma::vec coefficient;
  int iterations;
  bool converged;
  double relative_update;
  double kkt;
  double kkt_scale;
  double objective;
};

C7ChimeLasso c7_chime_lasso(
    const arma::mat& covariance, const arma::vec& difference,
    arma::vec coefficient, const double lambda, const double tolerance,
    const int max_iterations) {
  for (arma::uword j = 0; j < covariance.n_rows; ++j) {
    if (!(covariance(j, j) > 0.0) || !std::isfinite(covariance(j, j))) {
      Rcpp::stop(
        "A CHIME quadratic-lasso metric has a non-positive diagonal; "
        "no ridge or diagonal floor was applied."
      );
    }
  }

  arma::vec residual = difference - covariance * coefficient;
  double relative_update = std::numeric_limits<double>::infinity();
  double kkt = std::numeric_limits<double>::infinity();
  const double kkt_scale = 1.0 + c7_max_abs(difference) + lambda;
  bool converged = false;
  int iterations = 0;

  for (int sweep = 1; sweep <= max_iterations; ++sweep) {
    const arma::vec previous = coefficient;
    for (arma::uword j = 0; j < coefficient.n_elem; ++j) {
      const double diagonal = covariance(j, j);
      const double partial = residual(j) + diagonal * coefficient(j);
      const double updated = c7_soft(partial, lambda) / diagonal;
      const double change = updated - coefficient(j);
      if (change != 0.0) {
        coefficient(j) = updated;
        residual -= covariance.col(j) * change;
      }
    }
    iterations = sweep;
    relative_update = c7_max_abs(coefficient - previous) /
      std::max({1.0, c7_max_abs(previous), c7_max_abs(coefficient)});
    // Exact zeros define the lasso subgradient.  Reporting tolerances in R do
    // not alter the numerical KKT certificate.
    kkt = c7_chime_kkt(covariance, difference, coefficient, lambda);
    if (relative_update <= tolerance && kkt <= tolerance * kkt_scale) {
      converged = true;
      break;
    }
  }

  const double objective = 0.5 * arma::dot(
    coefficient, covariance * coefficient
  ) - arma::dot(coefficient, difference) +
    lambda * arma::accu(arma::abs(coefficient));
  return {
    coefficient, iterations, converged, relative_update, kkt, kkt_scale,
    objective
  };
}

arma::vec c7_chime_responsibility(
    const arma::mat& x, const double omega, const arma::vec& mu1,
    const arma::vec& mu2, const arma::vec& beta) {
  if (!(omega > 0.0 && omega < 1.0) || !std::isfinite(omega)) {
    Rcpp::stop(
      "A CHIME mixing weight reached the boundary; no probability clipping "
      "was applied."
    );
  }
  const arma::rowvec midpoint = 0.5 * (mu1 + mu2).t();
  arma::vec answer(x.n_rows, arma::fill::zeros);
  const double log_odds = std::log(omega) - std::log1p(-omega);
  for (arma::uword i = 0; i < x.n_rows; ++i) {
    const double eta = arma::dot(x.row(i) - midpoint, beta.t());
    if (!std::isfinite(eta)) {
      Rcpp::stop(
        "A CHIME discriminant score is non-finite; no clipping or rescaling "
        "was applied."
      );
    }
    const double value = log_odds - eta;
    if (value >= 0.0) {
      answer(i) = 1.0 / (1.0 + std::exp(-value));
    } else {
      const double exponential = std::exp(value);
      answer(i) = exponential / (1.0 + exponential);
    }
  }
  return answer;
}

arma::mat c7_chime_covariance(
    const arma::mat& x, const arma::vec& responsibility,
    const arma::vec& mu1, const arma::vec& mu2) {
  arma::mat answer(x.n_cols, x.n_cols, arma::fill::zeros);
  for (arma::uword i = 0; i < x.n_rows; ++i) {
    const arma::vec first = x.row(i).t() - mu1;
    const arma::vec second = x.row(i).t() - mu2;
    answer += (1.0 - responsibility(i)) * (first * first.t()) +
      responsibility(i) * (second * second.t());
  }
  answer /= static_cast<double>(x.n_rows);
  if (!answer.is_finite()) {
    Rcpp::stop(
      "The CHIME weighted covariance is non-finite; no clipping or scaling "
      "fallback was applied."
    );
  }
  return 0.5 * (answer + answer.t());
}

double c7_normal_cdf(const double value) {
  return R::pnorm(value, 0.0, 1.0, true, false);
}

struct C7KsScore {
  double score;
  double d_plus;
  double d_minus;
};

C7KsScore c7_ks_score(arma::vec value) {
  value = arma::sort(value);
  const double n = static_cast<double>(value.n_elem);
  double d_plus = 0.0;
  double d_minus = 0.0;
  for (arma::uword i = 0; i < value.n_elem; ++i) {
    const double distribution = c7_normal_cdf(value(i));
    d_plus = std::max(
      d_plus, (static_cast<double>(i) + 1.0) / n - distribution
    );
    d_minus = std::max(
      d_minus, distribution - static_cast<double>(i) / n
    );
  }
  return {std::sqrt(n) * std::max(d_plus, d_minus), d_plus, d_minus};
}

arma::vec c7_standardize_vector(const arma::vec& value) {
  const double n = static_cast<double>(value.n_elem);
  const double center = arma::mean(value);
  const arma::vec centered = value - center;
  const double scale = std::sqrt(arma::dot(centered, centered) / (n - 1.0));
  if (!(scale > 0.0) || !std::isfinite(scale)) {
    Rcpp::stop(
      "An IF-PCA feature has zero or non-finite n-1 sample standard "
      "deviation; no feature was silently dropped."
    );
  }
  return centered / scale;
}

double c7_squared_distance(const arma::rowvec& first,
                           const arma::rowvec& second) {
  const arma::rowvec difference = first - second;
  return arma::dot(difference, difference);
}

}  // namespace


// [[Rcpp::export]]
Rcpp::List cpp_ch7_chime_fit(
    const arma::mat& x, const double initial_omega,
    const arma::vec& initial_mu1, const arma::vec& initial_mu2,
    const arma::mat& initial_covariance, const arma::vec& lambda,
    const double tolerance, const int max_sweeps) {
  if (x.n_rows < 2 || x.n_cols < 1 || lambda.n_elem < 1) {
    Rcpp::stop("Invalid CHIME core dimensions.");
  }
  if (initial_mu1.n_elem != x.n_cols || initial_mu2.n_elem != x.n_cols ||
      initial_covariance.n_rows != x.n_cols ||
      initial_covariance.n_cols != x.n_cols) {
    Rcpp::stop("Initial CHIME dimensions do not match the data.");
  }
  if (!(initial_omega > 0.0 && initial_omega < 1.0) ||
      !std::isfinite(initial_omega)) {
    Rcpp::stop("The initial CHIME mixing weight must lie strictly in (0, 1).");
  }

  const arma::uword stages = lambda.n_elem;
  const arma::uword updates = stages - 1;
  arma::vec omega(stages, arma::fill::zeros);
  arma::mat mu1(stages, x.n_cols, arma::fill::zeros);
  arma::mat mu2(stages, x.n_cols, arma::fill::zeros);
  arma::mat beta(stages, x.n_cols, arma::fill::zeros);
  arma::mat gamma_history(x.n_rows, updates, arma::fill::zeros);
  arma::vec parameter_change(updates, arma::fill::zeros);
  arma::ivec lasso_iterations(stages, arma::fill::zeros);
  arma::uvec lasso_converged(stages, arma::fill::zeros);
  arma::vec lasso_update(stages, arma::fill::zeros);
  arma::vec lasso_kkt(stages, arma::fill::zeros);
  arma::vec lasso_kkt_scale(stages, arma::fill::zeros);
  arma::vec lasso_objective(stages, arma::fill::zeros);
  arma::uvec label_swapped(stages, arma::fill::zeros);
  Rcpp::List covariance_history(stages);

  double canonical_omega = initial_omega;
  arma::vec canonical_mu1 = initial_mu1;
  arma::vec canonical_mu2 = initial_mu2;
  if (canonical_omega > 0.5) {
    canonical_omega = 1.0 - canonical_omega;
    std::swap(canonical_mu1, canonical_mu2);
    label_swapped(0) = 1u;
  }
  omega(0) = canonical_omega;
  mu1.row(0) = canonical_mu1.t();
  mu2.row(0) = canonical_mu2.t();
  covariance_history[0] = initial_covariance;
  C7ChimeLasso initial_fit = c7_chime_lasso(
    initial_covariance, canonical_mu1 - canonical_mu2,
    arma::vec(x.n_cols, arma::fill::zeros), lambda(0), tolerance, max_sweeps
  );
  if (!initial_fit.converged || !initial_fit.coefficient.is_finite()) {
    Rcpp::stop(
      "The initial CHIME quadratic-lasso update did not meet its KKT and "
      "relative-update tolerances."
    );
  }
  beta.row(0) = initial_fit.coefficient.t();
  lasso_iterations(0) = initial_fit.iterations;
  lasso_converged(0) = initial_fit.converged ? 1u : 0u;
  lasso_update(0) = initial_fit.relative_update;
  lasso_kkt(0) = initial_fit.kkt;
  lasso_kkt_scale(0) = initial_fit.kkt_scale;
  lasso_objective(0) = initial_fit.objective;

  arma::mat current_covariance = initial_covariance;
  for (arma::uword stage = 1; stage < stages; ++stage) {
    const arma::vec old_mu1 = mu1.row(stage - 1).t();
    const arma::vec old_mu2 = mu2.row(stage - 1).t();
    const arma::vec old_beta = beta.row(stage - 1).t();
    const double old_omega = omega(stage - 1);
    const arma::vec responsibility = c7_chime_responsibility(
      x, old_omega, old_mu1, old_mu2, old_beta
    );

    const double second_weight = arma::accu(responsibility);
    const double first_weight = static_cast<double>(x.n_rows) - second_weight;
    if (!(first_weight > 0.0 && second_weight > 0.0) ||
        !std::isfinite(first_weight) || !std::isfinite(second_weight)) {
      Rcpp::stop(
        "A CHIME responsibility update emptied a component; no probability "
        "floor or component reset was applied."
      );
    }

    double new_omega = second_weight / static_cast<double>(x.n_rows);
    arma::vec new_mu1 = x.t() * (1.0 - responsibility) / first_weight;
    arma::vec new_mu2 = x.t() * responsibility / second_weight;
    arma::vec canonical_responsibility = responsibility;
    arma::vec lasso_start = old_beta;
    if (new_omega > 0.5) {
      new_omega = 1.0 - new_omega;
      std::swap(new_mu1, new_mu2);
      canonical_responsibility = 1.0 - responsibility;
      lasso_start = -old_beta;
      label_swapped(stage) = 1u;
    }
    gamma_history.col(stage - 1) = canonical_responsibility;
    current_covariance = c7_chime_covariance(
      x, canonical_responsibility, new_mu1, new_mu2
    );
    C7ChimeLasso fit = c7_chime_lasso(
      current_covariance, new_mu1 - new_mu2, lasso_start, lambda(stage),
      tolerance, max_sweeps
    );
    if (!fit.converged || !fit.coefficient.is_finite()) {
      Rcpp::stop(
        "A CHIME quadratic-lasso update did not meet its KKT and "
        "relative-update tolerances."
      );
    }

    omega(stage) = new_omega;
    mu1.row(stage) = new_mu1.t();
    mu2.row(stage) = new_mu2.t();
    beta.row(stage) = fit.coefficient.t();
    covariance_history[stage] = current_covariance;
    lasso_iterations(stage) = fit.iterations;
    lasso_converged(stage) = fit.converged ? 1u : 0u;
    lasso_update(stage) = fit.relative_update;
    lasso_kkt(stage) = fit.kkt;
    lasso_kkt_scale(stage) = fit.kkt_scale;
    lasso_objective(stage) = fit.objective;
    parameter_change(stage - 1) = std::abs(new_omega - old_omega) +
      c7_max_abs(new_mu1 - old_mu1) + c7_max_abs(new_mu2 - old_mu2) +
      c7_max_abs(fit.coefficient - old_beta);
  }

  const arma::vec final_mu1 = mu1.row(stages - 1).t();
  const arma::vec final_mu2 = mu2.row(stages - 1).t();
  const arma::vec final_beta = beta.row(stages - 1).t();
  const double final_omega = omega(stages - 1);
  const arma::vec final_responsibility = c7_chime_responsibility(
    x, final_omega, final_mu1, final_mu2, final_beta
  );
  const arma::rowvec midpoint = 0.5 * (final_mu1 + final_mu2).t();
  arma::vec score(x.n_rows, arma::fill::zeros);
  arma::ivec cluster(x.n_rows, arma::fill::zeros);
  const double threshold = std::log(final_omega) - std::log1p(-final_omega);
  for (arma::uword i = 0; i < x.n_rows; ++i) {
    score(i) = arma::dot(x.row(i) - midpoint, final_beta.t());
    if (!std::isfinite(score(i))) {
      Rcpp::stop("The final CHIME score is non-finite.");
    }
    cluster(i) = score(i) >= threshold ? 1 : 2;
  }

  return Rcpp::List::create(
    Rcpp::Named("omega") = omega,
    Rcpp::Named("mu1") = mu1,
    Rcpp::Named("mu2") = mu2,
    Rcpp::Named("beta") = beta,
    Rcpp::Named("covariance") = current_covariance,
    Rcpp::Named("covariance_history") = covariance_history,
    Rcpp::Named("responsibility_history") = gamma_history,
    Rcpp::Named("responsibility") = final_responsibility,
    Rcpp::Named("score") = score,
    Rcpp::Named("threshold") = threshold,
    Rcpp::Named("cluster") = cluster,
    Rcpp::Named("parameter_change") = parameter_change,
    Rcpp::Named("lasso_iterations") = lasso_iterations,
    Rcpp::Named("lasso_converged") = lasso_converged,
    Rcpp::Named("lasso_relative_update") = lasso_update,
    Rcpp::Named("lasso_kkt") = lasso_kkt,
    Rcpp::Named("lasso_kkt_scale") = lasso_kkt_scale,
    Rcpp::Named("label_swapped") = label_swapped,
    Rcpp::Named("lasso_objective") = lasso_objective
  );
}


// [[Rcpp::export]]
Rcpp::List cpp_ch7_if_scores(const arma::mat& x,
                             const bool software_convention) {
  if (x.n_rows < 3 || x.n_cols < 1) {
    Rcpp::stop("IF-PCA requires at least three rows and one feature.");
  }
  arma::rowvec centers = arma::mean(x, 0);
  arma::rowvec scales(x.n_cols, arma::fill::zeros);
  arma::mat standardized(x.n_rows, x.n_cols, arma::fill::zeros);
  arma::vec score(x.n_cols, arma::fill::zeros);
  arma::vec d_plus(x.n_cols, arma::fill::zeros);
  arma::vec d_minus(x.n_cols, arma::fill::zeros);
  const double software_scale = std::sqrt(
    1.0 - 1.0 / static_cast<double>(x.n_rows)
  );

  for (arma::uword j = 0; j < x.n_cols; ++j) {
    const arma::vec centered = x.col(j) - centers(j);
    const double scale = std::sqrt(
      arma::dot(centered, centered) /
      (static_cast<double>(x.n_rows) - 1.0)
    );
    if (!(scale > 0.0) || !std::isfinite(scale)) {
      Rcpp::stop(
        "An IF-PCA feature has zero or non-finite n-1 sample standard "
        "deviation; no feature was silently dropped."
      );
    }
    scales(j) = scale;
    standardized.col(j) = centered / scale;
    arma::vec ks_input = standardized.col(j);
    if (software_convention) ks_input /= software_scale;
    const C7KsScore fit = c7_ks_score(ks_input);
    score(j) = fit.score;
    d_plus(j) = fit.d_plus;
    d_minus(j) = fit.d_minus;
  }
  return Rcpp::List::create(
    Rcpp::Named("center") = centers,
    Rcpp::Named("scale") = scales,
    Rcpp::Named("standardized") = standardized,
    Rcpp::Named("score") = score,
    Rcpp::Named("d_plus") = d_plus,
    Rcpp::Named("d_minus") = d_minus,
    Rcpp::Named("software_ks_scale") = software_scale
  );
}


// [[Rcpp::export]]
arma::vec cpp_ch7_if_null_scores(const int sample_size,
                                 const int repetitions,
                                 const bool software_convention) {
  if (sample_size < 3 || repetitions < 1) {
    Rcpp::stop("Invalid IF-PCA null calibration dimensions.");
  }
  arma::vec answer(repetitions, arma::fill::zeros);
  const double software_scale = std::sqrt(
    1.0 - 1.0 / static_cast<double>(sample_size)
  );
  for (int repetition = 0; repetition < repetitions; ++repetition) {
    arma::vec value(sample_size, arma::fill::zeros);
    for (int i = 0; i < sample_size; ++i) value(i) = R::rnorm(0.0, 1.0);
    value = c7_standardize_vector(value);
    if (software_convention) value /= software_scale;
    answer(repetition) = c7_ks_score(value).score;
  }
  return answer;
}


// [[Rcpp::export]]
Rcpp::List cpp_ch7_if_deterministic_kmeans(
    const arma::mat& embedding, const int clusters,
    const double tolerance, const int max_iterations) {
  if (embedding.n_rows < static_cast<arma::uword>(clusters) ||
      embedding.n_cols < 1 || clusters < 2) {
    return Rcpp::List::create(
      Rcpp::Named("valid") = false,
      Rcpp::Named("message") = "The IF-PCA embedding cannot support K centers."
    );
  }

  arma::mat centers(clusters, embedding.n_cols, arma::fill::zeros);
  std::vector<bool> chosen(embedding.n_rows, false);
  const arma::rowvec global_center = arma::mean(embedding, 0);
  arma::uword first = 0;
  double farthest = -1.0;
  for (arma::uword i = 0; i < embedding.n_rows; ++i) {
    const double distance = c7_squared_distance(embedding.row(i), global_center);
    if (distance > farthest) {
      farthest = distance;
      first = i;
    }
  }
  centers.row(0) = embedding.row(first);
  chosen[first] = true;

  for (int center = 1; center < clusters; ++center) {
    arma::uword selected = 0;
    double best = -1.0;
    bool found = false;
    for (arma::uword i = 0; i < embedding.n_rows; ++i) {
      if (chosen[i]) continue;
      double nearest = std::numeric_limits<double>::infinity();
      for (int previous = 0; previous < center; ++previous) {
        nearest = std::min(
          nearest,
          c7_squared_distance(embedding.row(i), centers.row(previous))
        );
      }
      if (!found || nearest > best) {
        best = nearest;
        selected = i;
        found = true;
      }
    }
    if (!found || !(best > 0.0)) {
      return Rcpp::List::create(
        Rcpp::Named("valid") = false,
        Rcpp::Named("message") =
          "The IF-PCA embedding has fewer than K distinct rows."
      );
    }
    centers.row(center) = embedding.row(selected);
    chosen[selected] = true;
  }

  arma::ivec assignment(embedding.n_rows, arma::fill::value(-1));
  arma::ivec previous(embedding.n_rows, arma::fill::value(-2));
  bool converged = false;
  int iterations = 0;
  double relative_update = std::numeric_limits<double>::infinity();

  for (int iteration = 1; iteration <= max_iterations; ++iteration) {
    previous = assignment;
    for (arma::uword i = 0; i < embedding.n_rows; ++i) {
      int selected = 0;
      double best = c7_squared_distance(embedding.row(i), centers.row(0));
      for (int center = 1; center < clusters; ++center) {
        const double candidate = c7_squared_distance(
          embedding.row(i), centers.row(center)
        );
        // Strict inequality makes exact distance ties select the smaller label.
        if (candidate < best) {
          best = candidate;
          selected = center;
        }
      }
      assignment(i) = selected;
    }

    arma::mat updated(clusters, embedding.n_cols, arma::fill::zeros);
    arma::ivec counts(clusters, arma::fill::zeros);
    for (arma::uword i = 0; i < embedding.n_rows; ++i) {
      updated.row(assignment(i)) += embedding.row(i);
      counts(assignment(i)) += 1;
    }
    for (int center = 0; center < clusters; ++center) {
      if (counts(center) == 0) {
        return Rcpp::List::create(
          Rcpp::Named("valid") = false,
          Rcpp::Named("message") =
            "A deterministic IF-PCA k-means update produced an empty cluster; no center reset was applied."
        );
      }
      updated.row(center) /= static_cast<double>(counts(center));
    }
    relative_update = c7_max_abs_matrix(updated - centers) /
      std::max({1.0, c7_max_abs_matrix(centers), c7_max_abs_matrix(updated)});
    centers = updated;
    iterations = iteration;
    if (arma::all(assignment == previous) || relative_update <= tolerance) {
      converged = true;
      break;
    }
  }

  // Assign once more to the returned centers, retaining the same tie rule.
  double objective = 0.0;
  for (arma::uword i = 0; i < embedding.n_rows; ++i) {
    int selected = 0;
    double best = c7_squared_distance(embedding.row(i), centers.row(0));
    for (int center = 1; center < clusters; ++center) {
      const double candidate = c7_squared_distance(
        embedding.row(i), centers.row(center)
      );
      if (candidate < best) {
        best = candidate;
        selected = center;
      }
    }
    assignment(i) = selected;
    objective += best;
  }

  return Rcpp::List::create(
    Rcpp::Named("valid") = true,
    Rcpp::Named("cluster") = assignment + 1,
    Rcpp::Named("centers") = centers,
    Rcpp::Named("objective") = objective,
    Rcpp::Named("iterations") = iterations,
    Rcpp::Named("converged") = converged,
    Rcpp::Named("relative_update") = relative_update,
    Rcpp::Named("initialization") =
      "farthest from the global center, then farthest from prior centers",
    Rcpp::Named("tie_rule") = "smallest cluster index"
  );
}
