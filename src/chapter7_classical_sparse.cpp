// Chapter 7: classical and sparse prototype-clustering kernels.
//
// The routines below implement the displayed objectives directly.  They do
// not invent a ridge, pseudo-inverse, eigenvalue floor, random start, or empty
// cluster repair.  Every such choice is made by the R caller and recorded.

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp17)]]

#include <RcppArmadillo.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <sstream>
#include <string>
#include <unordered_set>
#include <vector>

namespace {

struct Assignment {
  arma::uvec labels;
  arma::mat distances;
  int ties;
};

struct RepairResult {
  arma::uvec labels;
  arma::mat centers;
  int repairs;
  int ties;
};

struct WeightResult {
  arma::vec weights;
  double threshold;
  double gamma;
  double l1_norm;
  double l2_norm;
  double constraint_residual;
  double stationarity_residual;
  int iterations;
  bool active;
  bool tie_branch;
  bool converged;
};

void ch7cs_check_matrix(const arma::mat& x, const char* name,
                        const arma::uword minimum_rows = 1) {
  if (x.n_rows < minimum_rows || x.n_cols == 0 || !x.is_finite()) {
    Rcpp::stop("`%s` must be a finite matrix with at least %d row(s) and one column.",
               name, static_cast<int>(minimum_rows));
  }
}

arma::uvec ch7cs_sizes(const arma::uvec& labels, const arma::uword clusters) {
  arma::uvec sizes(clusters, arma::fill::zeros);
  for (arma::uword i = 0; i < labels.n_elem; ++i) {
    if (labels(i) >= clusters) Rcpp::stop("Internal cluster label is invalid.");
    ++sizes(labels(i));
  }
  return sizes;
}

arma::mat ch7cs_mean_centers(const arma::mat& x, const arma::uvec& labels,
                             const arma::uword clusters) {
  arma::mat centers(clusters, x.n_cols, arma::fill::zeros);
  arma::uvec sizes = ch7cs_sizes(labels, clusters);
  if (arma::any(sizes == 0)) Rcpp::stop("A cluster is empty.");
  for (arma::uword i = 0; i < x.n_rows; ++i) centers.row(labels(i)) += x.row(i);
  for (arma::uword k = 0; k < clusters; ++k) centers.row(k) /= sizes(k);
  return centers;
}

double ch7cs_median(std::vector<double> values) {
  if (values.empty()) Rcpp::stop("A cluster is empty.");
  std::sort(values.begin(), values.end());
  const std::size_t n = values.size();
  if (n % 2 == 1) return values[n / 2];
  // The midpoint convention is part of the public sparse-K-median contract.
  return values[n / 2 - 1] / 2.0 + values[n / 2] / 2.0;
}

arma::mat ch7cs_median_centers(const arma::mat& x, const arma::uvec& labels,
                               const arma::uword clusters) {
  arma::uvec sizes = ch7cs_sizes(labels, clusters);
  if (arma::any(sizes == 0)) Rcpp::stop("A cluster is empty.");
  arma::mat centers(clusters, x.n_cols, arma::fill::zeros);
  for (arma::uword k = 0; k < clusters; ++k) {
    for (arma::uword j = 0; j < x.n_cols; ++j) {
      std::vector<double> values;
      values.reserve(sizes(k));
      for (arma::uword i = 0; i < x.n_rows; ++i) {
        if (labels(i) == k) values.push_back(x(i, j));
      }
      centers(k, j) = ch7cs_median(values);
    }
  }
  return centers;
}

arma::mat ch7cs_partial_centers(const arma::mat& x,
                                const arma::uvec& labels,
                                arma::mat centers,
                                const bool medians) {
  const arma::uvec sizes = ch7cs_sizes(labels, centers.n_rows);
  for (arma::uword k = 0; k < centers.n_rows; ++k) {
    if (sizes(k) == 0) continue;
    if (!medians) {
      centers.row(k).zeros();
      for (arma::uword i = 0; i < x.n_rows; ++i) {
        if (labels(i) == k) centers.row(k) += x.row(i);
      }
      centers.row(k) /= sizes(k);
    } else {
      for (arma::uword j = 0; j < x.n_cols; ++j) {
        std::vector<double> values;
        values.reserve(sizes(k));
        for (arma::uword i = 0; i < x.n_rows; ++i) {
          if (labels(i) == k) values.push_back(x(i, j));
        }
        centers(k, j) = ch7cs_median(values);
      }
    }
  }
  return centers;
}

Assignment ch7cs_assign(const arma::mat& x, const arma::mat& centers,
                        const arma::vec& weights, const bool squared) {
  const bool weighted = weights.n_elem != 0;
  if (weighted && weights.n_elem != x.n_cols) {
    Rcpp::stop("Internal feature-weight vector has the wrong length.");
  }
  Assignment answer;
  answer.labels = arma::uvec(x.n_rows, arma::fill::zeros);
  answer.distances = arma::mat(x.n_rows, centers.n_rows, arma::fill::zeros);
  answer.ties = 0;
  for (arma::uword i = 0; i < x.n_rows; ++i) {
    double best = std::numeric_limits<double>::infinity();
    arma::uword best_k = 0;
    int tied = 0;
    for (arma::uword k = 0; k < centers.n_rows; ++k) {
      double distance = 0.0;
      for (arma::uword j = 0; j < x.n_cols; ++j) {
        const double difference = x(i, j) - centers(k, j);
        const double contribution = squared ? difference * difference :
          std::abs(difference);
        distance += weighted ? weights(j) * contribution : contribution;
      }
      answer.distances(i, k) = distance;
      if (distance < best) {
        best = distance;
        best_k = k;
        tied = 0;
      } else if (distance == best) {
        // Iteration is in increasing k, so exact ties retain the first label.
        ++tied;
      }
    }
    answer.labels(i) = best_k;
    answer.ties += tied;
  }
  return answer;
}

RepairResult ch7cs_repair(const arma::mat& x, arma::uvec labels,
                          arma::mat centers, const arma::vec& weights,
                          const bool squared, const bool medians,
                          const int empty_action) {
  RepairResult answer;
  answer.repairs = 0;
  answer.ties = 0;
  const arma::uword clusters = centers.n_rows;
  while (true) {
    arma::uvec sizes = ch7cs_sizes(labels, clusters);
    arma::uword empty = clusters;
    for (arma::uword k = 0; k < clusters; ++k) {
      if (sizes(k) == 0) {
        empty = k;
        break;
      }
    }
    if (empty == clusters) break;
    if (empty_action == 0) {
      Rcpp::stop("An assignment produced an empty cluster; `empty_action = \"error\"` performs no repair.");
    }

    double farthest = -1.0;
    arma::uword chosen = x.n_rows;
    int tied = 0;
    for (arma::uword i = 0; i < x.n_rows; ++i) {
      const arma::uword donor = labels(i);
      if (sizes(donor) < 2) continue;
      double distance = 0.0;
      for (arma::uword j = 0; j < x.n_cols; ++j) {
        const double difference = x(i, j) - centers(donor, j);
        const double contribution = squared ? difference * difference :
          std::abs(difference);
        distance += weights.n_elem == 0 ? contribution :
          weights(j) * contribution;
      }
      if (distance > farthest) {
        farthest = distance;
        chosen = i;
        tied = 0;
      } else if (distance == farthest) {
        // Increasing row order makes the first row the deterministic winner.
        ++tied;
      }
    }
    if (chosen == x.n_rows) {
      Rcpp::stop("No cluster of size at least two can donate an observation to repair the empty cluster.");
    }
    labels(chosen) = empty;
    ++answer.repairs;
    answer.ties += tied;
    // Multiple empties are filled one at a time; non-empty centers are updated.
    centers = ch7cs_partial_centers(x, labels, centers, medians);
  }
  answer.labels = labels;
  answer.centers = centers;
  return answer;
}

std::string ch7cs_state_key(const arma::uvec& labels) {
  std::ostringstream stream;
  for (arma::uword i = 0; i < labels.n_elem; ++i) {
    if (i != 0) stream << ',';
    stream << labels(i);
  }
  return stream.str();
}

double ch7cs_wcss(const arma::mat& x, const arma::uvec& labels,
                  const arma::mat& centers, const arma::vec& weights,
                  const bool squared) {
  double answer = 0.0;
  for (arma::uword i = 0; i < x.n_rows; ++i) {
    for (arma::uword j = 0; j < x.n_cols; ++j) {
      const double difference = x(i, j) - centers(labels(i), j);
      const double contribution = squared ? difference * difference :
        std::abs(difference);
      answer += weights.n_elem == 0 ? contribution : weights(j) * contribution;
    }
  }
  return answer;
}

arma::vec ch7cs_sparse_kmeans_scores(const arma::mat& x,
                                     const arma::uvec& labels,
                                     const arma::mat& centers) {
  const arma::uvec sizes = ch7cs_sizes(labels, centers.n_rows);
  const arma::rowvec overall = arma::mean(x, 0);
  arma::vec scores(x.n_cols, arma::fill::zeros);
  for (arma::uword j = 0; j < x.n_cols; ++j) {
    long double between = 0.0L;
    for (arma::uword k = 0; k < centers.n_rows; ++k) {
      const long double difference = static_cast<long double>(centers(k, j)) -
        static_cast<long double>(overall(j));
      between += static_cast<long double>(sizes(k)) * difference * difference;
    }
    // The displayed ordered-pair BCSS is exactly twice TSS minus WCSS.
    scores(j) = static_cast<double>(2.0L * between);
  }
  return scores;
}

arma::rowvec ch7cs_global_median(const arma::mat& x) {
  arma::rowvec answer(x.n_cols, arma::fill::zeros);
  for (arma::uword j = 0; j < x.n_cols; ++j) {
    std::vector<double> values(x.n_rows);
    for (arma::uword i = 0; i < x.n_rows; ++i) values[i] = x(i, j);
    answer(j) = ch7cs_median(values);
  }
  return answer;
}

arma::vec ch7cs_sparse_kmedian_scores(const arma::mat& x,
                                      const arma::uvec& labels,
                                      const arma::mat& centers,
                                      const arma::rowvec& global_median,
                                      const double score_tolerance,
                                      int& roundoff_adjustments) {
  arma::vec scores(x.n_cols, arma::fill::zeros);
  roundoff_adjustments = 0;
  for (arma::uword j = 0; j < x.n_cols; ++j) {
    long double total = 0.0L;
    long double within = 0.0L;
    for (arma::uword i = 0; i < x.n_rows; ++i) {
      total += std::abs(static_cast<long double>(x(i, j)) -
                        static_cast<long double>(global_median(j)));
      within += std::abs(static_cast<long double>(x(i, j)) -
                         static_cast<long double>(centers(labels(i), j)));
    }
    double score = static_cast<double>(total - within);
    const double scale = std::max(1.0, static_cast<double>(total));
    if (score < -score_tolerance * scale) {
      Rcpp::stop("A sparse-K-median feature improvement is negative beyond `score_tol`.");
    }
    if (score < 0.0) {
      score = 0.0;
      ++roundoff_adjustments;
    }
    scores(j) = score;
  }
  return scores;
}

WeightResult ch7cs_weight_update(const arma::vec& scores, const double bound,
                                 const double tolerance,
                                 const int maximum_iterations) {
  WeightResult answer;
  answer.weights = arma::vec(scores.n_elem, arma::fill::zeros);
  answer.threshold = NA_REAL;
  answer.gamma = NA_REAL;
  answer.l1_norm = NA_REAL;
  answer.l2_norm = NA_REAL;
  answer.constraint_residual = std::numeric_limits<double>::infinity();
  answer.stationarity_residual = std::numeric_limits<double>::infinity();
  answer.iterations = 0;
  answer.active = false;
  answer.tie_branch = false;
  answer.converged = false;
  if (scores.n_elem == 0 || !scores.is_finite() || arma::any(scores < 0.0)) {
    Rcpp::stop("Sparse feature scores must be finite and non-negative.");
  }
  const double maximum = scores.max();
  if (!(maximum > 0.0)) {
    Rcpp::stop("All sparse-clustering feature scores are zero; no feature weight is defined.");
  }
  if (bound == 1.0) {
    arma::uword chosen = 0;
    for (arma::uword j = 1; j < scores.n_elem; ++j) {
      if (scores(j) > scores(chosen)) chosen = j;
    }
    answer.weights(chosen) = 1.0;
    double threshold = 0.0;
    for (arma::uword j = 0; j < scores.n_elem; ++j) {
      if (j != chosen) threshold = std::max(threshold, scores(j));
    }
    answer.threshold = threshold;
    answer.gamma = maximum - threshold;
    answer.l1_norm = 1.0;
    answer.l2_norm = 1.0;
    answer.constraint_residual = 0.0;
    answer.stationarity_residual = 0.0;
    answer.active = true;
    answer.tie_branch = arma::accu(scores == maximum) > 1;
    answer.converged = true;
    return answer;
  }

  const double norm0 = arma::norm(scores, 2);
  const double ratio0 = arma::accu(scores) / norm0;
  double threshold = 0.0;
  arma::vec soft = scores;
  if (ratio0 > bound + tolerance) {
    answer.active = true;
    double lower = 0.0;
    double upper = maximum;
    for (int iteration = 1; iteration <= maximum_iterations; ++iteration) {
      threshold = lower / 2.0 + upper / 2.0;
      soft = arma::clamp(scores - threshold, 0.0,
                         std::numeric_limits<double>::infinity());
      const double norm = arma::norm(soft, 2);
      const double ratio = norm > 0.0 ? arma::accu(soft) / norm : 1.0;
      answer.iterations = iteration;
      if (std::abs(ratio - bound) <= tolerance) break;
      if (ratio > bound) lower = threshold;
      else upper = threshold;
    }
    soft = arma::clamp(scores - threshold, 0.0,
                       std::numeric_limits<double>::infinity());
  }
  const double gamma = arma::norm(soft, 2);
  if (!(gamma > 0.0) || !std::isfinite(gamma)) {
    Rcpp::stop("The sparse weight normalization is zero or non-finite.");
  }
  answer.weights = soft / gamma;
  answer.threshold = threshold;
  answer.gamma = gamma;
  answer.l1_norm = arma::accu(answer.weights);
  answer.l2_norm = arma::norm(answer.weights, 2);
  answer.constraint_residual = answer.active ?
    std::abs(answer.l1_norm - bound) :
    std::max(0.0, answer.l1_norm - bound);
  double stationarity = 0.0;
  for (arma::uword j = 0; j < scores.n_elem; ++j) {
    const double residual = answer.weights(j) > 0.0 ?
      std::abs(scores(j) - threshold - gamma * answer.weights(j)) :
      std::max(0.0, scores(j) - threshold);
    stationarity = std::max(stationarity, residual);
  }
  answer.stationarity_residual = stationarity;
  answer.converged = answer.constraint_residual <= 10.0 * tolerance &&
    stationarity <= 10.0 * tolerance * std::max(1.0, maximum) &&
    std::abs(answer.l2_norm - 1.0) <= 10.0 * tolerance;
  return answer;
}

Rcpp::List ch7cs_weight_list(const WeightResult& weight) {
  return Rcpp::List::create(
    Rcpp::Named("threshold") = weight.threshold,
    Rcpp::Named("gamma") = weight.gamma,
    Rcpp::Named("l1_norm") = weight.l1_norm,
    Rcpp::Named("l2_norm") = weight.l2_norm,
    Rcpp::Named("constraint_residual") = weight.constraint_residual,
    Rcpp::Named("stationarity_residual") = weight.stationarity_residual,
    Rcpp::Named("iterations") = weight.iterations,
    Rcpp::Named("constraint_active") = weight.active,
    Rcpp::Named("s_one_tie_branch") = weight.tie_branch,
    Rcpp::Named("certified") = weight.converged
  );
}

struct GmmEstep {
  arma::mat responsibilities;
  double log_likelihood;
  double row_sum_residual;
};

GmmEstep ch7cs_gmm_estep(const arma::mat& x, const arma::mat& means,
                         const arma::cube& covariances,
                         const arma::vec& proportions) {
  const double log_two_pi = std::log(2.0 * arma::datum::pi);
  arma::mat log_joint(x.n_rows, means.n_rows, arma::fill::zeros);
  for (arma::uword k = 0; k < means.n_rows; ++k) {
    arma::mat lower;
    if (!arma::chol(lower, covariances.slice(k), "lower")) {
      Rcpp::stop("A Gaussian component covariance is not positive definite; no ridge or eigenvalue floor was inserted.");
    }
    const double log_determinant = 2.0 * arma::accu(arma::log(lower.diag()));
    for (arma::uword i = 0; i < x.n_rows; ++i) {
      const arma::vec residual = (x.row(i) - means.row(k)).t();
      arma::vec standardized;
      if (!arma::solve(standardized, arma::trimatl(lower), residual,
                       arma::solve_opts::fast)) {
        Rcpp::stop("Triangular solve failed in the Gaussian E-step.");
      }
      log_joint(i, k) = std::log(proportions(k)) -
        0.5 * (static_cast<double>(x.n_cols) * log_two_pi +
               log_determinant + arma::dot(standardized, standardized));
    }
  }
  GmmEstep answer;
  answer.responsibilities = arma::mat(x.n_rows, means.n_rows,
                                      arma::fill::zeros);
  answer.log_likelihood = 0.0;
  answer.row_sum_residual = 0.0;
  for (arma::uword i = 0; i < x.n_rows; ++i) {
    const double maximum = log_joint.row(i).max();
    const double sum_exp = arma::accu(arma::exp(log_joint.row(i) - maximum));
    if (!(sum_exp > 0.0) || !std::isfinite(sum_exp)) {
      Rcpp::stop("Log-sum-exp normalization failed in the Gaussian E-step.");
    }
    const double normalizer = maximum + std::log(sum_exp);
    answer.log_likelihood += normalizer;
    answer.responsibilities.row(i) = arma::exp(log_joint.row(i) - normalizer);
    answer.row_sum_residual = std::max(
      answer.row_sum_residual,
      std::abs(arma::accu(answer.responsibilities.row(i)) - 1.0)
    );
  }
  return answer;
}

struct GmmMstep {
  arma::mat means;
  arma::cube covariances;
  arma::vec proportions;
  arma::vec masses;
  arma::ivec raw_ranks;
  arma::vec minimum_eigenvalues;
};

GmmMstep ch7cs_gmm_mstep(const arma::mat& x,
                         const arma::mat& responsibilities,
                         const double ridge, const double rank_tolerance) {
  const arma::uword clusters = responsibilities.n_cols;
  GmmMstep answer;
  answer.means = arma::mat(clusters, x.n_cols, arma::fill::zeros);
  answer.covariances = arma::cube(x.n_cols, x.n_cols, clusters,
                                  arma::fill::zeros);
  answer.proportions = arma::vec(clusters, arma::fill::zeros);
  answer.masses = arma::sum(responsibilities, 0).t();
  answer.raw_ranks = arma::ivec(clusters, arma::fill::zeros);
  answer.minimum_eigenvalues = arma::vec(clusters, arma::fill::zeros);
  for (arma::uword k = 0; k < clusters; ++k) {
    const double mass = answer.masses(k);
    if (!(mass > 0.0) || !std::isfinite(mass)) {
      Rcpp::stop("A Gaussian component has zero or non-finite effective mass; no component repair was applied.");
    }
    answer.proportions(k) = mass / static_cast<double>(x.n_rows);
    answer.means.row(k) = (responsibilities.col(k).t() * x) / mass;
    arma::mat raw(x.n_cols, x.n_cols, arma::fill::zeros);
    for (arma::uword i = 0; i < x.n_rows; ++i) {
      const arma::vec residual = (x.row(i) - answer.means.row(k)).t();
      raw += responsibilities(i, k) * (residual * residual.t());
    }
    raw /= mass;
    raw = arma::symmatu(raw);
    arma::vec eigenvalues;
    if (!arma::eig_sym(eigenvalues, raw)) {
      Rcpp::stop("Eigenvalue diagnostics failed for a Gaussian covariance update.");
    }
    const double largest = eigenvalues.max();
    const double cutoff = largest > 0.0 ? rank_tolerance * largest : 0.0;
    answer.raw_ranks(k) = largest > 0.0 ?
      static_cast<int>(arma::accu(eigenvalues > cutoff)) : 0;
    arma::mat updated = raw;
    if (ridge > 0.0) updated.diag() += ridge;
    arma::mat lower;
    if (!arma::chol(lower, updated, "lower")) {
      Rcpp::stop("A Gaussian covariance M-step is singular; set an explicit positive `ridge` only if regularized EM is intended.");
    }
    if (!arma::eig_sym(eigenvalues, updated)) {
      Rcpp::stop("Eigenvalue diagnostics failed for an updated Gaussian covariance.");
    }
    answer.minimum_eigenvalues(k) = eigenvalues.min();
    answer.covariances.slice(k) = updated;
  }
  return answer;
}

double ch7cs_gmm_parameter_residual(const arma::mat& means,
                                    const arma::cube& covariances,
                                    const arma::vec& proportions,
                                    const GmmMstep& candidate) {
  double residual = arma::abs(means - candidate.means).max();
  residual = std::max(residual,
                      arma::abs(proportions - candidate.proportions).max());
  for (arma::uword k = 0; k < covariances.n_slices; ++k) {
    residual = std::max(
      residual,
      arma::abs(covariances.slice(k) - candidate.covariances.slice(k)).max()
    );
  }
  double scale = std::max(1.0, arma::abs(means).max());
  scale = std::max(scale, arma::abs(proportions).max());
  for (arma::uword k = 0; k < covariances.n_slices; ++k) {
    scale = std::max(scale, arma::abs(covariances.slice(k)).max());
  }
  return residual / scale;
}

} // namespace


// [[Rcpp::export]]
Rcpp::List cpp_ch7cs_lloyd_core(const arma::mat& x,
                                const arma::mat& initial_centers,
                                const int empty_action,
                                const double solver_tolerance,
                                const int maximum_iterations,
                                const bool keep_distances) {
  ch7cs_check_matrix(x, "x", 2);
  ch7cs_check_matrix(initial_centers, "initial_centers");
  if (initial_centers.n_cols != x.n_cols || initial_centers.n_rows > x.n_rows) {
    Rcpp::stop("`initial_centers` has incompatible dimensions.");
  }
  arma::vec no_weights;
  arma::mat centers = initial_centers;
  Assignment assignment = ch7cs_assign(x, centers, no_weights, true);
  int assignment_ties = assignment.ties;
  RepairResult repaired = ch7cs_repair(
    x, assignment.labels, centers, no_weights, true, false, empty_action
  );
  int repairs = repaired.repairs;
  int repair_ties = repaired.ties;
  arma::uvec labels = repaired.labels;
  centers = ch7cs_mean_centers(x, labels, centers.n_rows);
  std::vector<double> objective_trace;
  objective_trace.push_back(ch7cs_wcss(x, labels, centers, no_weights, true));
  std::unordered_set<std::string> states;
  states.insert(ch7cs_state_key(labels));
  bool converged = false;
  bool cycle = false;
  double maximum_increase = 0.0;
  int iterations = 0;
  for (int iteration = 1; iteration <= maximum_iterations; ++iteration) {
    assignment = ch7cs_assign(x, centers, no_weights, true);
    assignment_ties += assignment.ties;
    repaired = ch7cs_repair(
      x, assignment.labels, centers, no_weights, true, false, empty_action
    );
    repairs += repaired.repairs;
    repair_ties += repaired.ties;
    arma::uvec next_labels = repaired.labels;
    arma::mat next_centers = ch7cs_mean_centers(x, next_labels, centers.n_rows);
    const double objective = ch7cs_wcss(
      x, next_labels, next_centers, no_weights, true
    );
    const double previous = objective_trace.back();
    maximum_increase = std::max(maximum_increase, objective - previous);
    objective_trace.push_back(objective);
    iterations = iteration;
    if (arma::all(next_labels == labels)) {
      labels = next_labels;
      centers = next_centers;
      converged = true;
      break;
    }
    const std::string key = ch7cs_state_key(next_labels);
    if (states.find(key) != states.end()) {
      labels = next_labels;
      centers = next_centers;
      cycle = true;
      break;
    }
    states.insert(key);
    labels = next_labels;
    centers = next_centers;
  }
  assignment = ch7cs_assign(x, centers, no_weights, true);
  const double assignment_residual = arma::accu(assignment.labels != labels);
  const bool monotone = maximum_increase <= solver_tolerance *
    std::max(1.0, std::abs(objective_trace.front()));
  const bool certified = converged && !cycle && monotone &&
    assignment_residual == 0.0;
  return Rcpp::List::create(
    Rcpp::Named("cluster") = labels + 1,
    Rcpp::Named("centers") = centers,
    Rcpp::Named("size") = ch7cs_sizes(labels, centers.n_rows),
    Rcpp::Named("objective") = objective_trace.back(),
    Rcpp::Named("objective_trace") = objective_trace,
    Rcpp::Named("iterations") = iterations,
    Rcpp::Named("converged") = converged,
    Rcpp::Named("cycle_detected") = cycle,
    Rcpp::Named("assignment_ties") = assignment_ties,
    Rcpp::Named("empty_repairs") = repairs,
    Rcpp::Named("repair_ties") = repair_ties,
    Rcpp::Named("maximum_objective_increase") = maximum_increase,
    Rcpp::Named("assignment_residual") = assignment_residual,
    Rcpp::Named("certified") = certified,
    Rcpp::Named("distances") = keep_distances ?
      Rcpp::wrap(assignment.distances) : R_NilValue
  );
}


// [[Rcpp::export]]
Rcpp::List cpp_ch7cs_gmm_core(const arma::mat& x, arma::mat means,
                              arma::cube covariances,
                              arma::vec proportions, const double ridge,
                              const double solver_tolerance,
                              const double monotone_tolerance,
                              const double rank_tolerance,
                              const int maximum_iterations,
                              const bool keep_responsibilities) {
  ch7cs_check_matrix(x, "x", 2);
  if (means.n_rows == 0 || means.n_cols != x.n_cols || !means.is_finite() ||
      covariances.n_rows != x.n_cols || covariances.n_cols != x.n_cols ||
      covariances.n_slices != means.n_rows || !covariances.is_finite() ||
      proportions.n_elem != means.n_rows || !proportions.is_finite() ||
      arma::any(proportions <= 0.0)) {
    Rcpp::stop("Gaussian mixture inputs have incompatible dimensions or non-finite values.");
  }
  GmmEstep estep = ch7cs_gmm_estep(x, means, covariances, proportions);
  std::vector<double> log_likelihood_trace;
  log_likelihood_trace.push_back(estep.log_likelihood);
  bool converged = false;
  bool monotone = true;
  double maximum_decrease = 0.0;
  double update_residual = std::numeric_limits<double>::infinity();
  int iterations = 0;
  GmmMstep mstep;
  for (int iteration = 1; iteration <= maximum_iterations; ++iteration) {
    mstep = ch7cs_gmm_mstep(x, estep.responsibilities, ridge, rank_tolerance);
    update_residual = ch7cs_gmm_parameter_residual(
      means, covariances, proportions, mstep
    );
    means = mstep.means;
    covariances = mstep.covariances;
    proportions = mstep.proportions;
    GmmEstep next = ch7cs_gmm_estep(x, means, covariances, proportions);
    const double gain = next.log_likelihood - estep.log_likelihood;
    maximum_decrease = std::max(maximum_decrease, -gain);
    if (ridge == 0.0 &&
        gain < -monotone_tolerance *
          std::max(1.0, std::abs(estep.log_likelihood))) {
      monotone = false;
    }
    log_likelihood_trace.push_back(next.log_likelihood);
    estep = next;
    iterations = iteration;
    const double likelihood_residual = std::abs(gain) /
      std::max(1.0, std::abs(estep.log_likelihood));
    if (update_residual <= solver_tolerance &&
        likelihood_residual <= solver_tolerance) {
      converged = true;
      break;
    }
  }
  GmmMstep fixed = ch7cs_gmm_mstep(
    x, estep.responsibilities, ridge, rank_tolerance
  );
  const double fixed_point_residual = ch7cs_gmm_parameter_residual(
    means, covariances, proportions, fixed
  );
  arma::uvec labels(x.n_rows, arma::fill::zeros);
  int label_ties = 0;
  for (arma::uword i = 0; i < x.n_rows; ++i) {
    arma::uword best = 0;
    for (arma::uword k = 1; k < means.n_rows; ++k) {
      if (estep.responsibilities(i, k) > estep.responsibilities(i, best)) {
        best = k;
      } else if (estep.responsibilities(i, k) ==
                 estep.responsibilities(i, best)) {
        ++label_ties;
      }
    }
    labels(i) = best;
  }
  const bool certified = converged && (ridge > 0.0 || monotone) &&
    fixed_point_residual <= 10.0 * solver_tolerance &&
    estep.row_sum_residual <= 100.0 * std::numeric_limits<double>::epsilon();
  return Rcpp::List::create(
    Rcpp::Named("means") = means,
    Rcpp::Named("covariances") = covariances,
    Rcpp::Named("proportions") = proportions,
    Rcpp::Named("responsibilities") = keep_responsibilities ?
      Rcpp::wrap(estep.responsibilities) : R_NilValue,
    Rcpp::Named("cluster") = labels + 1,
    Rcpp::Named("log_likelihood") = estep.log_likelihood,
    Rcpp::Named("log_likelihood_trace") = log_likelihood_trace,
    Rcpp::Named("iterations") = iterations,
    Rcpp::Named("converged") = converged,
    Rcpp::Named("effective_masses") = fixed.masses,
    Rcpp::Named("raw_covariance_ranks") = fixed.raw_ranks,
    Rcpp::Named("minimum_eigenvalues") = fixed.minimum_eigenvalues,
    Rcpp::Named("responsibility_row_sum_residual") =
      estep.row_sum_residual,
    Rcpp::Named("fixed_point_residual") = fixed_point_residual,
    Rcpp::Named("maximum_log_likelihood_decrease") = maximum_decrease,
    Rcpp::Named("monotone") = ridge == 0.0 ? monotone : false,
    Rcpp::Named("label_ties") = label_ties,
    Rcpp::Named("certified") = certified
  );
}


// [[Rcpp::export]]
Rcpp::List cpp_ch7cs_sparse_kmeans_core(
    const arma::mat& x, const arma::mat& initial_centers,
    const double l1_bound, const int empty_action,
    const double weight_tolerance, const int weight_maximum_iterations,
    const double solver_tolerance, const int maximum_iterations,
    const bool keep_distances) {
  ch7cs_check_matrix(x, "x", 2);
  ch7cs_check_matrix(initial_centers, "initial_centers");
  if (initial_centers.n_cols != x.n_cols || initial_centers.n_rows > x.n_rows) {
    Rcpp::stop("`initial_centers` has incompatible dimensions.");
  }
  arma::vec no_weights;
  arma::mat centers = initial_centers;
  Assignment assignment = ch7cs_assign(x, centers, no_weights, true);
  int assignment_ties = assignment.ties;
  RepairResult repaired = ch7cs_repair(
    x, assignment.labels, centers, no_weights, true, false, empty_action
  );
  int repairs = repaired.repairs;
  int repair_ties = repaired.ties;
  arma::uvec labels = repaired.labels;
  centers = ch7cs_mean_centers(x, labels, centers.n_rows);
  std::unordered_set<std::string> states;
  states.insert(ch7cs_state_key(labels));
  std::vector<double> objective_trace;
  bool converged = false;
  bool cycle = false;
  double maximum_decrease = 0.0;
  int iterations = 0;
  WeightResult weight;
  arma::vec scores;
  for (int iteration = 1; iteration <= maximum_iterations; ++iteration) {
    scores = ch7cs_sparse_kmeans_scores(x, labels, centers);
    weight = ch7cs_weight_update(
      scores, l1_bound, weight_tolerance, weight_maximum_iterations
    );
    if (!weight.converged) {
      Rcpp::stop("The sparse-K-means feature-weight update is uncertified.");
    }
    const double objective = arma::dot(weight.weights, scores);
    if (!objective_trace.empty()) {
      maximum_decrease = std::max(
        maximum_decrease, objective_trace.back() - objective
      );
    }
    objective_trace.push_back(objective);
    assignment = ch7cs_assign(x, centers, weight.weights, true);
    assignment_ties += assignment.ties;
    repaired = ch7cs_repair(
      x, assignment.labels, centers, weight.weights, true, false, empty_action
    );
    repairs += repaired.repairs;
    repair_ties += repaired.ties;
    arma::uvec next_labels = repaired.labels;
    arma::mat next_centers = ch7cs_mean_centers(x, next_labels, centers.n_rows);
    iterations = iteration;
    if (arma::all(next_labels == labels)) {
      labels = next_labels;
      centers = next_centers;
      converged = true;
      break;
    }
    const std::string key = ch7cs_state_key(next_labels);
    if (states.find(key) != states.end()) {
      labels = next_labels;
      centers = next_centers;
      cycle = true;
      break;
    }
    states.insert(key);
    labels = next_labels;
    centers = next_centers;
  }
  // Recompute the weight block for the returned partition.
  scores = ch7cs_sparse_kmeans_scores(x, labels, centers);
  weight = ch7cs_weight_update(
    scores, l1_bound, weight_tolerance, weight_maximum_iterations
  );
  const double objective = arma::dot(weight.weights, scores);
  if (objective_trace.empty() || objective != objective_trace.back()) {
    if (!objective_trace.empty()) {
      maximum_decrease = std::max(
        maximum_decrease, objective_trace.back() - objective
      );
    }
    objective_trace.push_back(objective);
  }
  assignment = ch7cs_assign(x, centers, weight.weights, true);
  const double assignment_residual = arma::accu(assignment.labels != labels);
  const bool monotone = maximum_decrease <= solver_tolerance *
    std::max(1.0, std::abs(objective_trace.front()));
  const bool certified = converged && !cycle && monotone && weight.converged &&
    assignment_residual == 0.0;
  return Rcpp::List::create(
    Rcpp::Named("cluster") = labels + 1,
    Rcpp::Named("centers") = centers,
    Rcpp::Named("size") = ch7cs_sizes(labels, centers.n_rows),
    Rcpp::Named("weights") = weight.weights,
    Rcpp::Named("feature_scores") = scores,
    Rcpp::Named("half_scaled_feature_scores") = scores / 2.0,
    Rcpp::Named("objective") = objective,
    Rcpp::Named("objective_trace") = objective_trace,
    Rcpp::Named("iterations") = iterations,
    Rcpp::Named("converged") = converged,
    Rcpp::Named("cycle_detected") = cycle,
    Rcpp::Named("assignment_ties") = assignment_ties,
    Rcpp::Named("empty_repairs") = repairs,
    Rcpp::Named("repair_ties") = repair_ties,
    Rcpp::Named("maximum_objective_decrease") = maximum_decrease,
    Rcpp::Named("assignment_residual") = assignment_residual,
    Rcpp::Named("weight_certificate") = ch7cs_weight_list(weight),
    Rcpp::Named("weighted_within_loss") =
      ch7cs_wcss(x, labels, centers, weight.weights, true),
    Rcpp::Named("certified") = certified,
    Rcpp::Named("distances") = keep_distances ?
      Rcpp::wrap(assignment.distances) : R_NilValue
  );
}


// [[Rcpp::export]]
Rcpp::List cpp_ch7cs_sparse_kmedian_core(
    const arma::mat& x, const arma::mat& initial_centers,
    const double l1_bound, const int empty_action,
    const double score_tolerance, const double weight_tolerance,
    const int weight_maximum_iterations, const double solver_tolerance,
    const int maximum_iterations, const bool keep_distances) {
  ch7cs_check_matrix(x, "x", 2);
  ch7cs_check_matrix(initial_centers, "initial_centers");
  if (initial_centers.n_cols != x.n_cols || initial_centers.n_rows > x.n_rows) {
    Rcpp::stop("`initial_centers` has incompatible dimensions.");
  }
  arma::vec no_weights;
  arma::mat centers = initial_centers;
  Assignment assignment = ch7cs_assign(x, centers, no_weights, false);
  int assignment_ties = assignment.ties;
  RepairResult repaired = ch7cs_repair(
    x, assignment.labels, centers, no_weights, false, true, empty_action
  );
  int repairs = repaired.repairs;
  int repair_ties = repaired.ties;
  arma::uvec labels = repaired.labels;
  centers = ch7cs_median_centers(x, labels, centers.n_rows);
  const arma::rowvec global_median = ch7cs_global_median(x);
  std::unordered_set<std::string> states;
  states.insert(ch7cs_state_key(labels));
  std::vector<double> objective_trace;
  bool converged = false;
  bool cycle = false;
  double maximum_decrease = 0.0;
  int iterations = 0;
  int total_roundoff_adjustments = 0;
  int adjustments = 0;
  WeightResult weight;
  arma::vec scores;
  for (int iteration = 1; iteration <= maximum_iterations; ++iteration) {
    scores = ch7cs_sparse_kmedian_scores(
      x, labels, centers, global_median, score_tolerance, adjustments
    );
    total_roundoff_adjustments += adjustments;
    weight = ch7cs_weight_update(
      scores, l1_bound, weight_tolerance, weight_maximum_iterations
    );
    if (!weight.converged) {
      Rcpp::stop("The sparse-K-median feature-weight update is uncertified.");
    }
    const double objective = arma::dot(weight.weights, scores);
    if (!objective_trace.empty()) {
      maximum_decrease = std::max(
        maximum_decrease, objective_trace.back() - objective
      );
    }
    objective_trace.push_back(objective);
    assignment = ch7cs_assign(x, centers, weight.weights, false);
    assignment_ties += assignment.ties;
    repaired = ch7cs_repair(
      x, assignment.labels, centers, weight.weights, false, true, empty_action
    );
    repairs += repaired.repairs;
    repair_ties += repaired.ties;
    arma::uvec next_labels = repaired.labels;
    arma::mat next_centers = ch7cs_median_centers(x, next_labels,
                                                  centers.n_rows);
    iterations = iteration;
    if (arma::all(next_labels == labels)) {
      labels = next_labels;
      centers = next_centers;
      converged = true;
      break;
    }
    const std::string key = ch7cs_state_key(next_labels);
    if (states.find(key) != states.end()) {
      labels = next_labels;
      centers = next_centers;
      cycle = true;
      break;
    }
    states.insert(key);
    labels = next_labels;
    centers = next_centers;
  }
  scores = ch7cs_sparse_kmedian_scores(
    x, labels, centers, global_median, score_tolerance, adjustments
  );
  total_roundoff_adjustments += adjustments;
  weight = ch7cs_weight_update(
    scores, l1_bound, weight_tolerance, weight_maximum_iterations
  );
  const double objective = arma::dot(weight.weights, scores);
  if (objective_trace.empty() || objective != objective_trace.back()) {
    if (!objective_trace.empty()) {
      maximum_decrease = std::max(
        maximum_decrease, objective_trace.back() - objective
      );
    }
    objective_trace.push_back(objective);
  }
  assignment = ch7cs_assign(x, centers, weight.weights, false);
  const double assignment_residual = arma::accu(assignment.labels != labels);
  const bool monotone = maximum_decrease <= solver_tolerance *
    std::max(1.0, std::abs(objective_trace.front()));
  const bool certified = converged && !cycle && monotone && weight.converged &&
    assignment_residual == 0.0;
  return Rcpp::List::create(
    Rcpp::Named("cluster") = labels + 1,
    Rcpp::Named("centers") = centers,
    Rcpp::Named("global_median") = global_median,
    Rcpp::Named("size") = ch7cs_sizes(labels, centers.n_rows),
    Rcpp::Named("weights") = weight.weights,
    Rcpp::Named("feature_scores") = scores,
    Rcpp::Named("objective") = objective,
    Rcpp::Named("objective_trace") = objective_trace,
    Rcpp::Named("iterations") = iterations,
    Rcpp::Named("converged") = converged,
    Rcpp::Named("cycle_detected") = cycle,
    Rcpp::Named("assignment_ties") = assignment_ties,
    Rcpp::Named("empty_repairs") = repairs,
    Rcpp::Named("repair_ties") = repair_ties,
    Rcpp::Named("score_roundoff_adjustments") =
      total_roundoff_adjustments,
    Rcpp::Named("maximum_objective_decrease") = maximum_decrease,
    Rcpp::Named("assignment_residual") = assignment_residual,
    Rcpp::Named("weight_certificate") = ch7cs_weight_list(weight),
    Rcpp::Named("weighted_within_loss") =
      ch7cs_wcss(x, labels, centers, weight.weights, false),
    Rcpp::Named("certified") = certified,
    Rcpp::Named("distances") = keep_distances ?
      Rcpp::wrap(assignment.distances) : R_NilValue
  );
}
