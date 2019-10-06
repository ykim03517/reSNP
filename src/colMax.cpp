#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector colMax(NumericMatrix m) {
  int ncol = m.ncol();
  NumericVector max(ncol);
  for (int i = 0; i < ncol; i++)
    // Get col i with m(_, i).
    max[i] = Rcpp::max( m(_, i) );
  return max;
}


