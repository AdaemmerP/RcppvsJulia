#include <Rcpp.h>
#include <RcppParallel.h>
// [[Rcpp::depends(RcppParallel)]]
using namespace Rcpp;
using namespace RcppParallel;

struct CIworker : public Worker {
  /* input to read from */
  const RVector<double> X;
  const double a1, a2, z;
  const int n;
  /* output matrix to write to */
  RMatrix<double> mat;
  /* initialize from Rcpp input variables and output matrix */
  CIworker(const NumericVector X, double a1, double a2, double z, int n,
           NumericMatrix mat) : X(X), a1(a1), a2(a2), z(z), n(n), mat(mat) {}
  /* function call operator that works for the specified range (begin/end) */
  void operator()(std::size_t begin, std::size_t end) {
    for (std::size_t i = begin; i < end; i++) {
      double c0 = X[i]/n;
      double c1 = z*sqrt(c0*(1-c0)/n);
      double b0 = (1/(1 + z*z/n)) * (c0 + z*z/(2*n));
      double b1 = z/(1 + z*z/n) * sqrt( (c0*(1-c0)/n) + z*z/(4*n*n) );
      /* calculate intervals and write to output matrix */
      mat(i, 0) = c0-c1; mat(i, 1) = c0+c1;
      mat(i, 2) = b0-b1; mat(i, 3) = b0+b1;
      mat(i, 4) = R::qbeta(a1, X[i], n-X[i]+1, true, false);
      mat(i, 5) = R::qbeta(a2, X[i]+1, n-X[i], true, false);
    }
  }
};

// [[Rcpp::export]]
NumericVector CIrcpp2(int n, double p, int N, double alpha) {

  const double& a1 = alpha/2, a2 = 1-alpha/2;
  const double&  z = R::qnorm(a2, 0, 1, true, false);
  const NumericVector& X(Rcpp::rbinom(N, n, p));

  /* allocate the output matrix */
  NumericMatrix CIV(N, 6);
  /* create a worker */
  CIworker f1(X, a1, a2, z, n, CIV);
  /* call it with parallelFor */
  parallelFor(0, N, f1);

  // return vector with coverage rates
  return NumericVector::create(
    mean( (CIV(_, 0) <= p) & (p <= CIV(_, 1)) ),
    mean( (CIV(_, 2) <= p) & (p <= CIV(_, 3)) ),
    mean( (CIV(_, 4) <= p) & (p <= CIV(_, 5)) ));
}
