#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]

                     NumericVector arrayC(NumericVector input, IntegerVector dim) {
                  input.attr("dim") = dim;
                  return input;
                     }

