# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#' @title A random walk Metropolis sampler for generating the standard Laplace distribution using Rcpp.
#' @description A random walk Metropolis sampler for generating the standard Laplace distribution using Rcpp.
#' @param sigma the standard deviation of the increment
#' @param x0 the initial value
#' @param N sample size
#' @return a random sample of size \code{N} for the standard Laplace distribution
#' @examples
#' \dontrun{
#' CrwM(sigma= 2, x0 = 25, N = 2000)
#' }
#' @export
CrwM <- function(sigma, x0, N) {
    .Call('_StatComp20004_CrwM', PACKAGE = 'StatComp20004', sigma, x0, N)
}

