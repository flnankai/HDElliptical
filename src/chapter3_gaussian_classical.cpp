// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp17)]]

#include <RcppArmadillo.h>

#include <cmath>
#include <limits>

// Stable row directions and log-radii.  Computing a direction after first
// dividing by the largest absolute coordinate avoids overflow in x'x and
// retains subnormal coordinates whenever double precision permits it.
// [[Rcpp::export]]
Rcpp::List cpp_ch3_gaussian_sign_geometry(const arma::mat& residual) {
  const arma::uword n = residual.n_rows;
  const arma::uword p = residual.n_cols;
  arma::mat directions(n, p, arma::fill::zeros);
  arma::vec log_radii(n, arma::fill::zeros);

  for (arma::uword i = 0; i < n; ++i) {
    double row_max = 0.0;
    for (arma::uword j = 0; j < p; ++j) {
      const double value = residual(i, j);
      if (!std::isfinite(value)) {
        Rcpp::stop("The sphericized residual matrix contains a non-finite value.");
      }
      row_max = std::max(row_max, std::abs(value));
    }
    if (!(row_max > 0.0)) {
      Rcpp::stop("A residual is exactly zero; its multivariate sign is undefined.");
    }

    long double sumsq = 0.0L;
    for (arma::uword j = 0; j < p; ++j) {
      const long double scaled =
        static_cast<long double>(residual(i, j)) /
        static_cast<long double>(row_max);
      sumsq += scaled * scaled;
    }
    if (!(sumsq > 0.0L) || !std::isfinite(sumsq)) {
      Rcpp::stop("A residual norm is numerically undefined.");
    }
    const long double scaled_norm = std::sqrt(sumsq);
    log_radii(i) = std::log(row_max) +
      static_cast<double>(std::log(scaled_norm));
    for (arma::uword j = 0; j < p; ++j) {
      directions(i, j) = static_cast<double>(
        (static_cast<long double>(residual(i, j)) /
         static_cast<long double>(row_max)) / scaled_norm
      );
    }
  }

  return Rcpp::List::create(
    Rcpp::Named("directions") = directions,
    Rcpp::Named("log_radii") = log_radii
  );
}


// [[Rcpp::export]]
Rcpp::List cpp_ch3_gaussian_rank_shape(const arma::mat& directions,
                                        const arma::vec& scores,
                                        const double score_second_moment) {
  const arma::uword n = directions.n_rows;
  const arma::uword p = directions.n_cols;
  if (scores.n_elem != n) {
    Rcpp::stop("`scores` must have one entry per observation.");
  }
  if (!(score_second_moment > 0.0) ||
      !std::isfinite(score_second_moment)) {
    Rcpp::stop("`score_second_moment` must be finite and positive.");
  }

  arma::mat scatter(p, p, arma::fill::zeros);
  for (arma::uword i = 0; i < n; ++i) {
    const double score = scores(i);
    if (!std::isfinite(score)) {
      Rcpp::stop("Every score must be finite.");
    }
    scatter += score * directions.row(i).t() * directions.row(i);
  }
  scatter /= static_cast<double>(n);

  const double trace_scatter = arma::trace(scatter);
  const double trace_square = arma::accu(scatter % scatter);
  const double contrast = trace_square -
    trace_scatter * trace_scatter / static_cast<double>(p);
  // A tiny negative value can arise only from the final subtraction.  It is
  // returned unchanged: the R layer decides whether the statistic is usable.
  const double statistic =
    static_cast<double>(n) * static_cast<double>(p) *
    (static_cast<double>(p) + 2.0) * contrast /
    (2.0 * score_second_moment);

  return Rcpp::List::create(
    Rcpp::Named("scatter") = scatter,
    Rcpp::Named("trace") = trace_scatter,
    Rcpp::Named("trace_square") = trace_square,
    Rcpp::Named("contrast") = contrast,
    Rcpp::Named("statistic") = statistic
  );
}
