#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
NumericVector trueG(NumericVector g, double c) {
  int n = g.length();
  NumericVector out(n);
  for (int i = 0; i < n; i++) {
    if ( g(i) == 0 ) {
      out[i] = g(i);
    } else if ( g(i) == 1 ) {
      out[i] = c;
    } else  {
      out[i] = g(i) - 1;
    }
  }
  return out;
}
