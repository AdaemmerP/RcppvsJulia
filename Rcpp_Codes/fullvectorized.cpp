#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::plugins(cpp11)]]

NumericVector f2(double j, int n, double zz, double z, double a1, double a2) {
    double c0, c1, b0, b1;
    c0 = j/n;
    c1 = z*sqrt(c0*(1-c0)/n);
    b0 = (1/(1 + zz/n)) * (c0 + zz/(2*n));
    b1 = z/(1 + zz/n) * sqrt( (c0*(1-c0)/n) + zz/(4*n*n) );
    return NumericVector::create(
      c0-c1, c0+c1, b0-b1, b0+b1,
      R::qbeta(a1, j, n-j+1, true, false),
      R::qbeta(a2, j+1, n-j, true, false));
}

// [[Rcpp::export]]
NumericVector CIrcpp4(const int& n, const double& p, const int& N, const double& alpha) {
  double L1, L2, L3;
  const double& a1 = alpha/2, a2 = 1-alpha/2;
  const double&  z = R::qnorm(a2, 0, 1, true, false), zz = z*z;
  NumericVector::iterator i;
  NumericVector X(Rcpp::rbinom(N, n, p));
  NumericMatrix WWC(n+1, 6);

  for (int j = 0; j <= n; j++) WWC.row(j) = f2(j, n, zz, z, a1, a2);

  for (i = X.begin(); i != X.end(); ++i) {
    if ((WWC(*i, 0) <= p) && (p <= WWC(*i, 1))) L1++; /* Wald CI */
    if ((WWC(*i, 2) <= p) && (p <= WWC(*i, 3))) L2++; /* Wilson CI */
    if ((WWC(*i, 4) <= p) && (p <= WWC(*i, 5))) L3++; /* Clopper-Pearson CI */
  }
  return NumericVector::create(L1/N, L2/N, L3/N);
}
