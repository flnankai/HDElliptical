// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp17)]]

#include <RcppArmadillo.h>

#include <cmath>

namespace {

struct TraceUStatistics {
  long double y1;
  long double y2;
  long double y3;
  long double y4;
  long double y5;
  long double t1;
  long double t2;
};

TraceUStatistics trace_u_statistics(const arma::mat& x) {
  const arma::uword n = x.n_rows;
  if (n < 4) {
    Rcpp::stop("At least four observations are required for fourth-order trace U-statistics.");
  }
  arma::mat gram = x * x.t();
  arma::vec off_sum(n, arma::fill::zeros);
  long double diagonal_sum = 0.0L;
  long double off_sum_total = 0.0L;
  long double squared_off_sum = 0.0L;

  for (arma::uword i = 0; i < n; ++i) {
    diagonal_sum += gram(i, i);
    long double row = 0.0L;
    for (arma::uword j = 0; j < n; ++j) {
      if (i == j) continue;
      const long double value = gram(i, j);
      row += value;
      squared_off_sum += value * value;
    }
    off_sum(i) = static_cast<double>(row);
    off_sum_total += row;
  }

  long double row_square_sum = 0.0L;
  for (arma::uword i = 0; i < n; ++i) {
    const long double value = off_sum(i);
    row_square_sum += value * value;
  }
  const long double triple_sum = row_square_sum - squared_off_sum;
  const long double quadruple_sum =
    off_sum_total * off_sum_total - 4.0L * row_square_sum +
    2.0L * squared_off_sum;

  const long double nd = static_cast<long double>(n);
  const long double p2 = nd * (nd - 1.0L);
  const long double p3 = p2 * (nd - 2.0L);
  const long double p4 = p3 * (nd - 3.0L);
  const long double y1 = diagonal_sum / nd;
  const long double y2 = squared_off_sum / p2;
  const long double y3 = off_sum_total / p2;
  const long double y4 = triple_sum / p3;
  const long double y5 = quadruple_sum / p4;

  TraceUStatistics result;
  result.y1 = y1;
  result.y2 = y2;
  result.y3 = y3;
  result.y4 = y4;
  result.y5 = y5;
  result.t1 = y1 - y3;
  result.t2 = y2 - 2.0L * y4 + y5;
  return result;
}

long double cross_trace_u_statistic(const arma::mat& x,
                                    const arma::mat& y) {
  const arma::uword n1 = x.n_rows;
  const arma::uword n2 = y.n_rows;
  arma::mat cross = x * y.t();
  long double square_sum = 0.0L;
  long double total = 0.0L;
  long double row_square_sum = 0.0L;
  long double column_square_sum = 0.0L;

  for (arma::uword i = 0; i < n1; ++i) {
    long double row = 0.0L;
    for (arma::uword j = 0; j < n2; ++j) {
      const long double value = cross(i, j);
      square_sum += value * value;
      row += value;
      total += value;
    }
    row_square_sum += row * row;
  }
  for (arma::uword j = 0; j < n2; ++j) {
    long double column = 0.0L;
    for (arma::uword i = 0; i < n1; ++i) {
      column += cross(i, j);
    }
    column_square_sum += column * column;
  }

  const long double a = square_sum;
  const long double b1 = column_square_sum - a;
  const long double b2 = row_square_sum - a;
  const long double d = total * total - row_square_sum -
    column_square_sum + a;
  const long double n1d = static_cast<long double>(n1);
  const long double n2d = static_cast<long double>(n2);
  return a / (n1d * n2d) -
    b1 / (n1d * n2d * (n1d - 1.0L)) -
    b2 / (n1d * n2d * (n2d - 1.0L)) +
    d / (n1d * n2d * (n1d - 1.0L) * (n2d - 1.0L));
}

} // namespace


// [[Rcpp::export]]
Rcpp::List cpp_ch3_gaussian_trace_u_statistics(const arma::mat& x) {
  const TraceUStatistics value = trace_u_statistics(x);
  return Rcpp::List::create(
    Rcpp::Named("Y1") = static_cast<double>(value.y1),
    Rcpp::Named("Y2") = static_cast<double>(value.y2),
    Rcpp::Named("Y3") = static_cast<double>(value.y3),
    Rcpp::Named("Y4") = static_cast<double>(value.y4),
    Rcpp::Named("Y5") = static_cast<double>(value.y5),
    Rcpp::Named("T1") = static_cast<double>(value.t1),
    Rcpp::Named("T2") = static_cast<double>(value.t2)
  );
}


// [[Rcpp::export]]
Rcpp::List cpp_ch3_gaussian_li_chen(const arma::mat& x,
                                    const arma::mat& y) {
  if (x.n_rows < 4 || y.n_rows < 4) {
    Rcpp::stop("Each sample requires at least four observations.");
  }
  if (x.n_cols != y.n_cols) {
    Rcpp::stop("The samples must have the same number of variables.");
  }
  const TraceUStatistics tx = trace_u_statistics(x);
  const TraceUStatistics ty = trace_u_statistics(y);
  const long double cross = cross_trace_u_statistic(x, y);
  const long double statistic = tx.t2 + ty.t2 - 2.0L * cross;
  return Rcpp::List::create(
    Rcpp::Named("A1") = static_cast<double>(tx.t2),
    Rcpp::Named("A2") = static_cast<double>(ty.t2),
    Rcpp::Named("C") = static_cast<double>(cross),
    Rcpp::Named("T") = static_cast<double>(statistic)
  );
}
