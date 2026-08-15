// Zhang--Feng (2026) radial--directional tests for elliptical models.
//
// The kernel treats the supplied location and positive-definite shape as the
// complete standardisation contract.  It computes log radii and directions
// without materialising an inverse square root with extreme eigenvalue
// powers, and optionally performs the paper's intrinsic radial--directional
// bootstrap.  It does not estimate, ridge, floor, project, or otherwise repair
// the supplied shape.

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp17)]]

#include <RcppArmadillo.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <string>

namespace {

struct ZfrdStandardised {
  arma::mat directions;
  arma::vec log_radii;
  arma::vec eigenvalues_scaled;
  double shape_scale;
  double data_scale;
  double minimum_log_radius;
  double maximum_log_radius;
  arma::uword subtraction_fallbacks;
};

struct ZfrdCorrelations {
  arma::vec correlations;
  double log_variance_sum;
  arma::vec direction_variance_sums;
  arma::uword roundoff_clips;
};

long double zfrd_safe_norm(const arma::vec& x) {
  long double maximum = 0.0L;
  for (arma::uword j = 0; j < x.n_elem; ++j) {
    maximum = std::max(
      maximum, std::fabs(static_cast<long double>(x(j)))
    );
  }
  if (maximum == 0.0L) return 0.0L;
  long double sum = 0.0L;
  for (arma::uword j = 0; j < x.n_elem; ++j) {
    const long double ratio = static_cast<long double>(x(j)) / maximum;
    sum += ratio * ratio;
  }
  return maximum * std::sqrt(sum);
}

double zfrd_checked_double(const long double value,
                           const char* quantity) {
  const long double limit = static_cast<long double>(
    std::numeric_limits<double>::max()
  );
  if (!std::isfinite(value) || value > limit || value < -limit) {
    Rcpp::stop(
      "The radial--directional computation produced a non-finite %s; no "
      "ridge, floor, projection, clipping, or pseudoinverse is applied.",
      quantity
    );
  }
  return static_cast<double>(value);
}

ZfrdStandardised zfrd_standardise(const arma::mat& x,
                                  const arma::vec& location,
                                  const arma::mat& shape,
                                  const double relative_zero_tol) {
  const arma::uword n = x.n_rows;
  const arma::uword p = x.n_cols;

  double shape_scale = 0.0;
  for (arma::uword j = 0; j < p; ++j) {
    for (arma::uword k = 0; k < p; ++k) {
      shape_scale = std::max(shape_scale, std::fabs(shape(j, k)));
    }
  }
  if (!(shape_scale > 0.0) || !std::isfinite(shape_scale)) {
    Rcpp::stop("`shape` must have a finite positive numerical scale.");
  }
  const arma::mat shape_scaled = shape / shape_scale;
  arma::vec eigenvalues;
  arma::mat eigenvectors;
  if (!arma::eig_sym(eigenvalues, eigenvectors, shape_scaled)) {
    Rcpp::stop("The supplied shape eigendecomposition failed.");
  }
  if (eigenvalues.n_elem != p || !eigenvalues.is_finite() ||
      !eigenvectors.is_finite() || arma::any(eigenvalues <= 0.0)) {
    Rcpp::stop(
      "`shape` must be numerically symmetric positive definite; no "
      "eigenvalue floor or positive-definite projection is applied."
    );
  }

  double data_scale = 0.0;
  for (arma::uword i = 0; i < n; ++i) {
    for (arma::uword j = 0; j < p; ++j) {
      data_scale = std::max(data_scale, std::fabs(x(i, j)));
    }
  }
  for (arma::uword j = 0; j < p; ++j) {
    data_scale = std::max(data_scale, std::fabs(location(j)));
  }
  if (!(data_scale > 0.0)) data_scale = 1.0;

  ZfrdStandardised out;
  out.directions.zeros(n, p);
  out.log_radii.set_size(n);
  out.eigenvalues_scaled = eigenvalues;
  out.shape_scale = shape_scale;
  out.data_scale = data_scale;
  out.minimum_log_radius = std::numeric_limits<double>::infinity();
  out.maximum_log_radius = -std::numeric_limits<double>::infinity();
  out.subtraction_fallbacks = 0;
  const long double log_constant =
    std::log(static_cast<long double>(data_scale)) -
    0.5L * std::log(static_cast<long double>(shape_scale));
  for (arma::uword i = 0; i < n; ++i) {
    arma::vec residual(p);
    for (arma::uword j = 0; j < p; ++j) {
      const long double difference =
        static_cast<long double>(x(i, j)) -
        static_cast<long double>(location(j));
      double scaled = static_cast<double>(
        difference / static_cast<long double>(data_scale)
      );
      if (!std::isfinite(scaled)) {
        scaled = x(i, j) / data_scale - location(j) / data_scale;
        ++out.subtraction_fallbacks;
      }
      if (!std::isfinite(scaled)) {
        Rcpp::stop(
          "A finite observation-minus-location residual cannot be "
          "represented after safe scaling."
        );
      }
      residual(j) = scaled;
    }
    if (zfrd_safe_norm(residual) == 0.0L) {
      Rcpp::stop(
        "A fitted residual is exactly zero, so its log radius and direction "
        "are undefined; no perturbation is applied."
      );
    }

    const arma::vec projected = eigenvectors.t() * residual;
    long double maximum_log_component =
      -std::numeric_limits<long double>::infinity();
    arma::vec log_absolute(p);
    for (arma::uword j = 0; j < p; ++j) {
      if (projected(j) == 0.0) {
        log_absolute(j) = -std::numeric_limits<double>::infinity();
      } else {
        const long double value =
          std::log(std::fabs(static_cast<long double>(projected(j)))) -
          0.5L * std::log(static_cast<long double>(eigenvalues(j)));
        log_absolute(j) = static_cast<double>(value);
        maximum_log_component = std::max(maximum_log_component, value);
      }
    }
    if (!std::isfinite(maximum_log_component)) {
      Rcpp::stop("A nonzero residual vanished in the shape eigenbasis.");
    }

    arma::vec eigen_direction(p, arma::fill::zeros);
    long double norm_squared = 0.0L;
    for (arma::uword j = 0; j < p; ++j) {
      if (projected(j) == 0.0) continue;
      const long double magnitude = std::exp(
        static_cast<long double>(log_absolute(j)) - maximum_log_component
      );
      const long double signed_value = projected(j) < 0.0 ?
        -magnitude : magnitude;
      eigen_direction(j) = static_cast<double>(signed_value);
      norm_squared += signed_value * signed_value;
    }
    if (!(norm_squared > 0.0L) || !std::isfinite(norm_squared)) {
      Rcpp::stop("A standardized residual norm is not finite and positive.");
    }
    const long double scaled_norm = std::sqrt(norm_squared);
    eigen_direction /= static_cast<double>(scaled_norm);
    const arma::vec direction = eigenvectors * eigen_direction;
    if (!direction.is_finite()) {
      Rcpp::stop("A standardized direction is non-finite.");
    }
    out.directions.row(i) = direction.t();

    const long double log_radius = log_constant + maximum_log_component +
      std::log(scaled_norm);
    if (!std::isfinite(log_radius)) {
      Rcpp::stop("A standardized log radius is non-finite.");
    }
    out.log_radii(i) = static_cast<double>(log_radius);
    out.minimum_log_radius = std::min(
      out.minimum_log_radius, out.log_radii(i)
    );
    out.maximum_log_radius = std::max(
      out.maximum_log_radius, out.log_radii(i)
    );
  }
  if (relative_zero_tol > 0.0) {
    const long double log_ratio =
      static_cast<long double>(out.minimum_log_radius) -
      static_cast<long double>(out.maximum_log_radius);
    const long double log_tolerance = std::log(
      static_cast<long double>(relative_zero_tol)
    );
    if (log_ratio <= log_tolerance) {
      Rcpp::stop(
        "A standardized radius is too small relative to the largest fitted "
        "radius under `radius_zero_tol`; its log-radius correlation is "
        "treated as undefined."
      );
    }
  }
  return out;
}

ZfrdCorrelations zfrd_correlations(const arma::vec& log_radii,
                                   const arma::mat& directions) {
  const arma::uword n = directions.n_rows;
  const arma::uword p = directions.n_cols;
  long double mean_log = 0.0L;
  for (arma::uword i = 0; i < n; ++i) {
    mean_log += static_cast<long double>(log_radii(i));
  }
  mean_log /= static_cast<long double>(n);
  long double log_ss = 0.0L;
  long double maximum_centered_log = 0.0L;
  for (arma::uword i = 0; i < n; ++i) {
    const long double centered =
      static_cast<long double>(log_radii(i)) - mean_log;
    log_ss += centered * centered;
    maximum_centered_log = std::max(
      maximum_centered_log, std::fabs(centered)
    );
  }
  const long double epsilon = static_cast<long double>(
    std::numeric_limits<double>::epsilon()
  );
  const long double log_variance_tolerance = 256.0L * epsilon * epsilon *
    std::max(
      1.0L,
      static_cast<long double>(n) * maximum_centered_log *
        maximum_centered_log
    );
  if (!(log_ss > log_variance_tolerance) || !std::isfinite(log_ss)) {
    Rcpp::stop(
      "The fitted log radii have zero, numerically unresolved, or "
      "non-finite empirical variance."
    );
  }

  ZfrdCorrelations out;
  out.correlations.set_size(p);
  out.direction_variance_sums.set_size(p);
  out.log_variance_sum = zfrd_checked_double(log_ss, "log-radius variance");
  out.roundoff_clips = 0;
  for (arma::uword j = 0; j < p; ++j) {
    long double mean_u = 0.0L;
    for (arma::uword i = 0; i < n; ++i) {
      mean_u += static_cast<long double>(directions(i, j));
    }
    mean_u /= static_cast<long double>(n);
    long double u_ss = 0.0L;
    long double cross = 0.0L;
    long double maximum_centered_u = 0.0L;
    for (arma::uword i = 0; i < n; ++i) {
      const long double l_centered =
        static_cast<long double>(log_radii(i)) - mean_log;
      const long double u_centered =
        static_cast<long double>(directions(i, j)) - mean_u;
      u_ss += u_centered * u_centered;
      cross += l_centered * u_centered;
      maximum_centered_u = std::max(
        maximum_centered_u, std::fabs(u_centered)
      );
    }
    const long double u_variance_tolerance = 256.0L * epsilon * epsilon *
      std::max(
        1.0L,
        static_cast<long double>(n) * maximum_centered_u *
          maximum_centered_u
      );
    if (!(u_ss > u_variance_tolerance) || !std::isfinite(u_ss)) {
      Rcpp::stop(
        "At least one fitted direction coordinate has zero or non-finite "
        "empirical variance."
      );
    }
    long double correlation = cross / std::sqrt(log_ss * u_ss);
    if (!std::isfinite(correlation)) {
      Rcpp::stop("A radial--directional correlation is non-finite.");
    }
    if (correlation > 1.0L || correlation < -1.0L) {
      if (std::fabs(std::fabs(correlation) - 1.0L) > 1e-12L) {
        Rcpp::stop("A radial--directional correlation materially exceeds one.");
      }
      correlation = correlation > 0.0L ? 1.0L : -1.0L;
      ++out.roundoff_clips;
    }
    out.correlations(j) = static_cast<double>(correlation);
    out.direction_variance_sums(j) = zfrd_checked_double(
      u_ss, "direction-coordinate variance"
    );
  }
  return out;
}

Rcpp::List zfrd_statistics(const arma::vec& correlations,
                           const arma::uword n,
                           const arma::uword p) {
  long double sum_squared = 0.0L;
  long double maximum_squared = 0.0L;
  for (arma::uword j = 0; j < p; ++j) {
    const long double value = static_cast<long double>(correlations(j));
    const long double squared = value * value;
    sum_squared += squared;
    maximum_squared = std::max(maximum_squared, squared);
  }
  const long double t_sum = static_cast<long double>(n) * sum_squared;
  const long double t_max_raw =
    static_cast<long double>(n) * maximum_squared;
  const long double t_max = t_max_raw -
    2.0L * std::log(static_cast<long double>(p)) +
    std::log(std::log(static_cast<long double>(p)));
  return Rcpp::List::create(
    Rcpp::Named("T_sum") = zfrd_checked_double(t_sum, "sum statistic"),
    Rcpp::Named("T_max_raw") = zfrd_checked_double(
      t_max_raw, "raw max statistic"
    ),
    Rcpp::Named("T_max") = zfrd_checked_double(
      t_max, "centered max statistic"
    )
  );
}

} // anonymous namespace

// [[Rcpp::export]]
Rcpp::List cpp_zhang_feng_radial_directional(
    const arma::mat& x,
    const arma::vec& location,
    const arma::mat& shape,
    const double relative_zero_tol,
    const int bootstrap_replicates) {
  const arma::uword n = x.n_rows;
  const arma::uword p = x.n_cols;
  if (n < 3 || p < 2 || location.n_elem != p ||
      shape.n_rows != p || shape.n_cols != p) {
    Rcpp::stop("Invalid dimensions in the radial--directional native kernel.");
  }
  if (!x.is_finite() || !location.is_finite() || !shape.is_finite() ||
      !std::isfinite(relative_zero_tol) || relative_zero_tol < 0.0 ||
      relative_zero_tol >= 1.0 ||
      bootstrap_replicates < 0) {
    Rcpp::stop("Invalid finite-input contract in the native kernel.");
  }

  const ZfrdStandardised standardised = zfrd_standardise(
    x, location, shape, relative_zero_tol
  );
  const ZfrdCorrelations observed = zfrd_correlations(
    standardised.log_radii, standardised.directions
  );
  const Rcpp::List observed_statistics = zfrd_statistics(
    observed.correlations, n, p
  );

  Rcpp::List bootstrap = R_NilValue;
  if (bootstrap_replicates > 0) {
    if (bootstrap_replicates < 2) {
      Rcpp::stop("At least two bootstrap replicates are required.");
    }
    Rcpp::RNGScope scope;
    arma::vec t_sum(bootstrap_replicates);
    arma::vec t_max(bootstrap_replicates);
    arma::vec sampled_log_radii(n);
    arma::mat random_directions(n, p);
    for (int b = 0; b < bootstrap_replicates; ++b) {
      for (arma::uword i = 0; i < n; ++i) {
        arma::uword index = static_cast<arma::uword>(
          std::floor(R::runif(0.0, 1.0) * static_cast<double>(n))
        );
        if (index >= n) index = n - 1;
        sampled_log_radii(i) = standardised.log_radii(index);

        arma::vec row(p);
        for (arma::uword j = 0; j < p; ++j) {
          row(j) = R::rnorm(0.0, 1.0);
        }
        const long double row_norm = zfrd_safe_norm(row);
        if (!(row_norm > 0.0L) || !std::isfinite(row_norm)) {
          Rcpp::stop("A bootstrap Gaussian direction has invalid norm.");
        }
        row /= static_cast<double>(row_norm);
        random_directions.row(i) = row.t();
      }
      const ZfrdCorrelations replicate = zfrd_correlations(
        sampled_log_radii, random_directions
      );
      const Rcpp::List statistics = zfrd_statistics(
        replicate.correlations, n, p
      );
      t_sum(b) = Rcpp::as<double>(statistics["T_sum"]);
      t_max(b) = Rcpp::as<double>(statistics["T_max"]);
    }
    bootstrap = Rcpp::List::create(
      Rcpp::Named("T_sum") = t_sum,
      Rcpp::Named("T_max") = t_max,
      Rcpp::Named("mean_sum") = arma::mean(t_sum),
      Rcpp::Named("sd_sum") = arma::stddev(t_sum, 0),
      Rcpp::Named("mean_max") = arma::mean(t_max),
      Rcpp::Named("sd_max") = arma::stddev(t_max, 0)
    );
  }

  return Rcpp::List::create(
    Rcpp::Named("correlations") = observed.correlations,
    Rcpp::Named("directions") = standardised.directions,
    Rcpp::Named("log_radii") = standardised.log_radii,
    Rcpp::Named("T_sum") = observed_statistics["T_sum"],
    Rcpp::Named("T_max_raw") = observed_statistics["T_max_raw"],
    Rcpp::Named("T_max") = observed_statistics["T_max"],
    Rcpp::Named("log_radius_variance_sum") =
      observed.log_variance_sum,
    Rcpp::Named("direction_variance_sums") =
      observed.direction_variance_sums,
    Rcpp::Named("shape_eigenvalues_scaled") =
      standardised.eigenvalues_scaled,
    Rcpp::Named("shape_scale") = standardised.shape_scale,
    Rcpp::Named("data_scale") = standardised.data_scale,
    Rcpp::Named("minimum_log_radius") =
      standardised.minimum_log_radius,
    Rcpp::Named("maximum_log_radius") =
      standardised.maximum_log_radius,
    Rcpp::Named("subtraction_fallbacks") =
      static_cast<double>(standardised.subtraction_fallbacks),
    Rcpp::Named("correlation_roundoff_clips") =
      static_cast<double>(observed.roundoff_clips),
    Rcpp::Named("bootstrap") = bootstrap
  );
}
