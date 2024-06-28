#include <RcppArmadillo.h>
#include <fstream>
using namespace Rcpp;
using namespace std;
// [[Rcpp::plugins("cpp11")]]
// [[Rcpp::depends(RcppArmadillo)]]

double rnb(double r, double p){
  double y = R::rnbinom(r, p);
  return y;
}

// [[Rcpp::export]]
NumericVector apply_rnbinom(NumericVector counts, NumericVector proportions, double size) {

  return mapply(counts, proportions, rnb);
}


