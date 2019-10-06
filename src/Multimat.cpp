#include <RcppEigen.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix Multmat(Rcpp::NumericMatrix tmm,  Rcpp::NumericMatrix tm22) {
const Eigen::Map<Eigen::MatrixXd> ttm(Rcpp::as<Eigen::Map<Eigen::MatrixXd> >(tmm));
const Eigen::Map<Eigen::MatrixXd> ttm2(Rcpp::as<Eigen::Map<Eigen::MatrixXd> >(tm22));
                       
Eigen::MatrixXd prod = ttm*ttm2;
return(wrap(prod));
}

                       
