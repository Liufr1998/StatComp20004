#include <Rcpp.h>
#include <cmath>
using namespace Rcpp;

//' @title A random walk Metropolis sampler for generating the standard Laplace distribution using Rcpp.
//' @description A random walk Metropolis sampler for generating the standard Laplace distribution using Rcpp.
//' @param sigma the standard deviation of the increment
//' @param x0 the initial value
//' @param N sample size
//' @return a random sample of size \code{N} for the standard Laplace distribution
//' @examples
//' \dontrun{
//' CrwM(sigma= 2, x0 = 25, N = 2000)
//' }
//' @export
// [[Rcpp::export]]
NumericVector CrwM(double sigma, double x0, int N){
  NumericVector x(N);
  NumericVector u = runif(N);
  x[0] = x0;
  
  
  for (int i = 1; i < N; ++i) {
    NumericVector z = rnorm(1, x[i-1], sigma);
    if (u[i] <= (exp(-fabs(z[0])) / exp(-fabs(x[(i-1)])))) {
      x[i] = z[0];
    } else {
      x[i] = x[(i-1)];
    }
  }
  
  return x;
}
