#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
NumericVector rowMax(NumericMatrix m) {
  int nrow = m.nrow();
  NumericVector max(nrow);
  for (int i = 0; i < nrow; i++)
    // Get row i with m(i, _).
    max[i] = Rcpp::max( m(i, _) );
  return max;
}
