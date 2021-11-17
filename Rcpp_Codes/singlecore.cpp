#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector CIrcpp1(const int& n, const double& p, const int& N, const double& alpha) {

  const double& a1 = alpha/2, a2 = 1-alpha/2,
  z = R::qnorm(a2, 0, 1, true, false), zz = z+z;

  double x, c0, c1, b0, b1;
  NumericMatrix CIV(N, 6);

  for (int i = 0; i < N; i++) {
    x = R::rbinom(n, p);

    /* Wald CI */
    c0 = x/n;
    c1 = z*sqrt(c0*(1-c0)/n);
    CIV(i, 0) = c0-c1;
    CIV(i, 1) = c0+c1;

    /* Wilson CI */
    b0 = (1/(1 + zz/n)) * (c0 + zz/(2*n));
    b1 = z/(1 + zz/n) * sqrt( (c0*(1-c0)/n) + zz/(4*n*n) );
    CIV(i, 2) = b0 - b1;
    CIV(i, 3) = b0 + b1;

    /* Clopper-Pearson CI */
    CIV(i, 4) = R::qbeta(a1, x, n-x+1, true, false);
    CIV(i, 5) = R::qbeta(a2, x+1, n-x, true, false);
  }

  /* compute coverage rates */
  return NumericVector::create(
    mean( (CIV(_, 0) <= p) & (p <= CIV(_, 1)) ),
    mean( (CIV(_, 2) <= p) & (p <= CIV(_, 3)) ),
    mean( (CIV(_, 4) <= p) & (p <= CIV(_, 5)) ));
}
