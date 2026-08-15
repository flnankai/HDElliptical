#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp17)]]

namespace {

void require_finite(const arma::mat& x, const char* name) {
  if (!x.is_finite()) {
    Rcpp::stop("`%s` must contain only finite values.", name);
  }
}

void require_finite(const arma::vec& x, const char* name) {
  if (!x.is_finite()) {
    Rcpp::stop("`%s` must contain only finite values.", name);
  }
}

void require_condition_tolerance(const double tol) {
  if (!std::isfinite(tol) || tol <= 0.0 || tol >= 1.0) {
    Rcpp::stop("`tol` must be a finite number strictly between zero and one.");
  }
}

Rcpp::List chol_quadratic(const arma::mat& covariance,
                          const arma::vec& difference,
                          const double multiplier,
                          const double tol,
                          const char* covariance_name) {
  arma::mat lower;
  if (!arma::chol(lower, covariance, "lower")) {
    Rcpp::stop(
      "The %s covariance matrix is not positive definite; "
      "Hotelling's test does not use a generalized inverse.",
      covariance_name
    );
  }

  const double reciprocal_condition = arma::rcond(covariance);
  if (!std::isfinite(reciprocal_condition) || reciprocal_condition <= tol) {
    Rcpp::stop(
      "The %s covariance matrix is singular or numerically ill-conditioned "
      "at `tol = %.3g` (reciprocal condition number %.3g).",
      covariance_name, tol, reciprocal_condition
    );
  }

  arma::vec whitened;
  const bool solved = arma::solve(
    whitened, arma::trimatl(lower), difference, arma::solve_opts::fast
  );
  if (!solved || !whitened.is_finite()) {
    Rcpp::stop("The Cholesky triangular solve failed for the %s covariance matrix.",
               covariance_name);
  }

  const double quadratic = multiplier * arma::dot(whitened, whitened);
  if (!std::isfinite(quadratic)) {
    Rcpp::stop("Hotelling's quadratic form is not finite.");
  }

  return Rcpp::List::create(
    Rcpp::Named("t_squared") = quadratic,
    Rcpp::Named("rcond") = reciprocal_condition,
    Rcpp::Named("rank") = static_cast<int>(covariance.n_rows),
    Rcpp::Named("solver") = "Cholesky"
  );
}

void require_zero_tolerance(const double zero_tol) {
  if (!std::isfinite(zero_tol) || zero_tol < 0.0) {
    Rcpp::stop("`zero_tol` must be a finite non-negative number.");
  }
}

bool unit_direction(const arma::rowvec& value, const double zero_tol,
                    arma::rowvec& direction) {
  const double scale = arma::abs(value).max();
  direction.zeros(value.n_elem);
  if (!std::isfinite(scale)) {
    Rcpp::stop("A direction vector has a non-finite norm.");
  }
  if (scale == 0.0) {
    return false;
  }

  const arma::rowvec scaled = value / scale;
  const double scaled_norm = std::sqrt(arma::dot(scaled, scaled));
  if (!std::isfinite(scaled_norm)) {
    Rcpp::stop("A direction vector has a non-finite norm.");
  }
  if (scaled_norm <= zero_tol / scale) {
    return false;
  }
  direction = scaled / scaled_norm;
  return true;
}

bool difference_direction(const arma::rowvec& first,
                          const arma::rowvec& second,
                          const double zero_tol,
                          arma::rowvec& direction) {
  const double scale = std::max(
    arma::abs(first).max(), arma::abs(second).max()
  );
  direction.zeros(first.n_elem);
  if (scale == 0.0) {
    return false;
  }
  const arma::rowvec normalized = first / scale - second / scale;
  const double normalized_tolerance = zero_tol == 0.0 ?
    0.0 : zero_tol / scale;
  if (!std::isfinite(normalized_tolerance)) {
    return false;
  }
  return unit_direction(normalized, normalized_tolerance, direction);
}

bool signed_rank_direction(const arma::rowvec& first,
                           const arma::rowvec& second,
                           const arma::rowvec& location,
                           const double zero_tol,
                           arma::rowvec& direction) {
  const double scale = std::max(
    std::max(arma::abs(first).max(), arma::abs(second).max()),
    arma::abs(location).max()
  );
  direction.zeros(first.n_elem);
  if (scale == 0.0) {
    return false;
  }
  const arma::rowvec normalized = first / scale + second / scale -
    2.0 * (location / scale);
  const double normalized_tolerance = zero_tol == 0.0 ?
    0.0 : zero_tol / scale;
  if (!std::isfinite(normalized_tolerance)) {
    return false;
  }
  return unit_direction(normalized, normalized_tolerance, direction);
}

Rcpp::List chi_square_quadratic(const arma::mat& moment,
                                const arma::vec& difference,
                                const double multiplier,
                                const double tol,
                                const char* moment_name) {
  arma::mat lower;
  if (!arma::chol(lower, moment, "lower")) {
    Rcpp::stop(
      "The %s matrix is not positive definite; the classical chi-square "
      "calibration is unavailable.",
      moment_name
    );
  }

  const double reciprocal_condition = arma::rcond(moment);
  if (!std::isfinite(reciprocal_condition) || reciprocal_condition <= tol) {
    Rcpp::stop(
      "The %s matrix is singular or numerically ill-conditioned at "
      "`tol = %.3g` (reciprocal condition number %.3g).",
      moment_name, tol, reciprocal_condition
    );
  }

  arma::vec whitened;
  const bool solved = arma::solve(
    whitened, arma::trimatl(lower), difference, arma::solve_opts::fast
  );
  if (!solved || !whitened.is_finite()) {
    Rcpp::stop("The Cholesky triangular solve failed for the %s matrix.",
               moment_name);
  }

  const double statistic = multiplier * arma::dot(whitened, whitened);
  if (!std::isfinite(statistic)) {
    Rcpp::stop("The classical location-test quadratic form is not finite.");
  }

  return Rcpp::List::create(
    Rcpp::Named("statistic") = statistic,
    Rcpp::Named("rcond") = reciprocal_condition,
    Rcpp::Named("rank") = static_cast<int>(moment.n_rows),
    Rcpp::Named("solver") = "Cholesky"
  );
}

}  // namespace


// [[Rcpp::export]]
Rcpp::List cpp_hotelling_one_sample(const arma::mat& x,
                                    const arma::vec& mu,
                                    const double tol) {
  require_finite(x, "x");
  require_finite(mu, "mu");
  require_condition_tolerance(tol);

  const arma::uword n = x.n_rows;
  const arma::uword p = x.n_cols;
  if (p == 0 || n < 2) {
    Rcpp::stop("`x` must have at least two rows and one column.");
  }
  if (mu.n_elem != p) {
    Rcpp::stop("`mu` must have one value per column of `x`.");
  }
  if (n <= p) {
    Rcpp::stop("One-sample Hotelling's test requires `n > p`.");
  }

  const arma::vec sample_mean = arma::trans(arma::mean(x, 0));
  const arma::mat centered = x.each_row() - sample_mean.t();
  const arma::mat covariance = centered.t() * centered /
    static_cast<double>(n - 1);
  const arma::vec difference = sample_mean - mu;

  Rcpp::List result = chol_quadratic(
    covariance, difference, static_cast<double>(n), tol, "sample"
  );
  result["mean"] = sample_mean;
  result["difference"] = difference;
  result["covariance"] = covariance;
  return result;
}


// [[Rcpp::export]]
Rcpp::List cpp_hotelling_two_sample(const arma::mat& x,
                                    const arma::mat& y,
                                    const double tol) {
  require_finite(x, "x");
  require_finite(y, "y");
  require_condition_tolerance(tol);

  const arma::uword n1 = x.n_rows;
  const arma::uword n2 = y.n_rows;
  const arma::uword p = x.n_cols;
  if (p == 0 || y.n_cols != p) {
    Rcpp::stop("`x` and `y` must have the same positive number of columns.");
  }
  if (n1 < 2 || n2 < 2) {
    Rcpp::stop("Two-sample Hotelling's test requires at least two rows per sample.");
  }
  const arma::uword total = n1 + n2;
  if (total <= p + 1) {
    Rcpp::stop("Two-sample Hotelling's test requires `n1 + n2 > p + 1`.");
  }

  const arma::vec mean_x = arma::trans(arma::mean(x, 0));
  const arma::vec mean_y = arma::trans(arma::mean(y, 0));
  const arma::mat centered_x = x.each_row() - mean_x.t();
  const arma::mat centered_y = y.each_row() - mean_y.t();
  const arma::mat pooled_covariance =
    (centered_x.t() * centered_x + centered_y.t() * centered_y) /
    static_cast<double>(total - 2);
  const arma::vec difference = mean_x - mean_y;
  const double multiplier = static_cast<double>(n1) *
    static_cast<double>(n2) / static_cast<double>(total);

  Rcpp::List result = chol_quadratic(
    pooled_covariance, difference, multiplier, tol, "pooled"
  );
  result["mean_x"] = mean_x;
  result["mean_y"] = mean_y;
  result["difference"] = difference;
  result["covariance"] = pooled_covariance;
  return result;
}


// [[Rcpp::export]]
Rcpp::List cpp_spatial_sign_test(const arma::mat& x,
                                 const arma::vec& mu,
                                 const double tol,
                                 const double zero_tol) {
  require_finite(x, "x");
  require_finite(mu, "mu");
  require_condition_tolerance(tol);
  require_zero_tolerance(zero_tol);

  const arma::uword n = x.n_rows;
  const arma::uword p = x.n_cols;
  if (p == 0 || n < 2) {
    Rcpp::stop("`x` must have at least two rows and one column.");
  }
  if (mu.n_elem != p) {
    Rcpp::stop("`mu` must have one value per column of `x`.");
  }
  if (n <= p) {
    Rcpp::stop("The classical spatial-sign test requires `n > p`.");
  }

  arma::mat signs(n, p, arma::fill::zeros);
  arma::uword n_zero = 0;
  const arma::rowvec location = mu.t();
  for (arma::uword i = 0; i < n; ++i) {
    arma::rowvec direction(p);
    if (difference_direction(x.row(i), location, zero_tol, direction)) {
      signs.row(i) = direction;
    } else {
      ++n_zero;
    }
  }

  const arma::vec mean_sign = arma::trans(arma::mean(signs, 0));
  const arma::mat second_moment = signs.t() * signs /
    static_cast<double>(n);
  Rcpp::List result = chi_square_quadratic(
    second_moment, mean_sign, static_cast<double>(n), tol,
    "spatial-sign second-moment"
  );
  result["mean_sign"] = mean_sign;
  result["second_moment"] = second_moment;
  result["n_zero"] = static_cast<double>(n_zero);
  return result;
}


// [[Rcpp::export]]
Rcpp::List cpp_spatial_signed_rank_test(const arma::mat& x,
                                        const arma::vec& mu,
                                        const double tol,
                                        const double zero_tol) {
  require_finite(x, "x");
  require_finite(mu, "mu");
  require_condition_tolerance(tol);
  require_zero_tolerance(zero_tol);

  const arma::uword n = x.n_rows;
  const arma::uword p = x.n_cols;
  if (p == 0 || n < 2) {
    Rcpp::stop("`x` must have at least two rows and one column.");
  }
  if (mu.n_elem != p) {
    Rcpp::stop("`mu` must have one value per column of `x`.");
  }
  if (n <= p) {
    Rcpp::stop("The classical spatial signed-rank test requires `n > p`.");
  }

  arma::mat ranks(n, p, arma::fill::zeros);
  double n_zero_ordered_pairs = 0.0;
  const arma::rowvec location = mu.t();
  for (arma::uword i = 0; i < n; ++i) {
    for (arma::uword j = i; j < n; ++j) {
      arma::rowvec direction(p);
      const bool nonzero = signed_rank_direction(
        x.row(i), x.row(j), location, zero_tol, direction
      );
      if (nonzero) {
        ranks.row(i) += direction;
        if (j != i) {
          ranks.row(j) += direction;
        }
      } else {
        n_zero_ordered_pairs += (j == i) ? 1.0 : 2.0;
      }
    }
  }
  ranks /= static_cast<double>(n);

  const arma::vec mean_rank = arma::trans(arma::mean(ranks, 0));
  const arma::mat rank_second_moment = ranks.t() * ranks /
    static_cast<double>(n);
  // sqrt(n) * mean_rank has covariance 4 B_R, hence the required 1/4.
  Rcpp::List result = chi_square_quadratic(
    rank_second_moment, mean_rank, static_cast<double>(n) / 4.0, tol,
    "signed-rank second-moment"
  );
  result["mean_rank"] = mean_rank;
  result["rank_second_moment"] = rank_second_moment;
  result["n_zero_ordered_pairs"] = n_zero_ordered_pairs;
  return result;
}


// [[Rcpp::export]]
Rcpp::List cpp_spatial_rank_test(const arma::mat& x,
                                 const arma::mat& y,
                                 const double tol,
                                 const double zero_tol) {
  require_finite(x, "x");
  require_finite(y, "y");
  require_condition_tolerance(tol);
  require_zero_tolerance(zero_tol);

  const arma::uword n1 = x.n_rows;
  const arma::uword n2 = y.n_rows;
  const arma::uword p = x.n_cols;
  const arma::uword total = n1 + n2;
  if (p == 0 || y.n_cols != p) {
    Rcpp::stop("`x` and `y` must have the same positive number of columns.");
  }
  if (n1 < 2 || n2 < 2) {
    Rcpp::stop("The pooled spatial-rank test requires at least two rows per sample.");
  }
  if (total <= p) {
    Rcpp::stop("The pooled spatial-rank test requires `n1 + n2 > p`.");
  }

  const arma::mat pooled = arma::join_cols(x, y);
  arma::mat ranks(total, p, arma::fill::zeros);
  double n_zero_ordered_pairs = static_cast<double>(total);
  for (arma::uword i = 0; i < total; ++i) {
    for (arma::uword j = i + 1; j < total; ++j) {
      arma::rowvec direction(p);
      if (difference_direction(
            pooled.row(i), pooled.row(j), zero_tol, direction
          )) {
        ranks.row(i) += direction;
        ranks.row(j) -= direction;
      } else {
        n_zero_ordered_pairs += 2.0;
      }
    }
  }
  ranks /= static_cast<double>(total);

  const arma::vec mean_rank_x = arma::trans(
    arma::mean(ranks.rows(0, n1 - 1), 0)
  );
  const arma::vec mean_rank_y = arma::trans(
    arma::mean(ranks.rows(n1, total - 1), 0)
  );
  const arma::vec difference = mean_rank_x - mean_rank_y;
  const arma::mat rank_covariance = ranks.t() * ranks /
    static_cast<double>(total - 1);
  const double multiplier = static_cast<double>(n1) *
    static_cast<double>(n2) / static_cast<double>(total);

  Rcpp::List result = chi_square_quadratic(
    rank_covariance, difference, multiplier, tol, "pooled-rank covariance"
  );
  result["mean_rank_x"] = mean_rank_x;
  result["mean_rank_y"] = mean_rank_y;
  result["difference"] = difference;
  result["rank_covariance"] = rank_covariance;
  result["n_zero_ordered_pairs"] = n_zero_ordered_pairs;
  return result;
}
