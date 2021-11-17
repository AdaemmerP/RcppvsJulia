#include <Rcpp.h>
using namespace Rcpp;
// [[Rcpp::plugins(cpp11)]]

// [[Rcpp::export]]
NumericVector CIrcpp3(const int& n, const double& p, const int& N, const double& alpha) {

  double c0, c1, b0, b1;
  const double& a1 = alpha/2, a2 = 1-alpha/2;
  const double&  z = R::qnorm(a2, 0, 1, true, false), zz = z*z;
  double WaL, WaU, WiL, WiU, ClL, ClU, L1, L2, L3;
  std::unordered_map<int, double> CPL, CPU;

  NumericVector X(Rcpp::rbinom(N, n, p));
  NumericVector::iterator i;

  /* unordered map for Clopper-Pearson CI values */
  for (int j = 0; j <= n; j++) {
    CPL[j] = R::qbeta(a1, j, n-j+1, true, false);
    CPU[j] = R::qbeta(a2, j+1, n-j, true, false);
  }

  /* get intervals */
  for (i = X.begin(); i != X.end(); ++i) {
    
    /* Wald CI */
    c0  = *i/n;
    c1  = z*sqrt(c0*(1-c0)/n);
    WaL = c0-c1;
    WaU = c0+c1;

    /* Wilson CI */
    b0  = (1/(1 + zz/n)) * (c0 + zz/(2*n));
    b1  = z/(1 + zz/n) * sqrt( (c0*(1-c0)/n) + zz/(4*n*n) );
    WiL = b0 - b1;
    WiU = b0 + b1;

    if (    (WaL <= p) && (p <= WaU))     L1++;
    if (    (WiL <= p) && (p <= WiU))     L2++;
    if ((CPL[*i] <= p) && (p <= CPU[*i])) L3++;
  }

  /* compute coverage rates */
  return NumericVector::create(L1/N, L2/N, L3/N);
}