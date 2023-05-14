#include <Rcpp.h>
using namespace Rcpp;


//
//[[Rcpp::export]]
std::vector<double> weird1compute(NumericVector m, NumericVector v, NumericVector h) {
  int n = m.size();
  std::vector<double> x(n); //KLD increase for each point
  //double sum = 0;
  for(int i = 0; i < n; ++i) {
    x[i] = (pow(3*pow(h[i],2) - 6*h[i]*pow(m[i],2) + pow(m[i],4),2) -
      4*(9*pow(h[i],3) - 63*pow(h[i],2)*pow(m[i],2) + 45*h[i]*pow(m[i],4) - 7*pow(m[i],6))*v[i] +
      6*(21*pow(h[i],2) - 90*h[i]*pow(m[i],2) + 35*pow(m[i],4))*pow(v[i],2) + 60*(-3*h[i] + 7*pow(m[i],2))*pow(v[i],3) +
      105*pow(v[i],4))/(64.*pow(h[i],6)) ;
  }
  return x;
}

//[[Rcpp::export]]
double weird2compute(NumericVector m, NumericMatrix v, NumericVector h) {
  int n = m.size();
  double sum = 0;

  double sumij;
  double hi;
  double hj;
  double m1;
  double m2;
  double v1;
  double v2;
  double v12;

  for(int i = 0; i < (n-1); ++i) {
    for(int j = (i+1); j < n; ++j){
      hi = h[i];
      hj = h[j];
      m1 = m[i];
      m2 = m[j];
      v1 = v(i,i);
      v2 = v(j,j);
      v12 = v(i,j);
      if(v(i,j)!=0){
        sumij = (pow(m1,4)*pow(m2,4) + 6*pow(m1,2)*pow(m2,4)*v1 + 3*pow(m2,4)*pow(v1,2) +
          3*pow(hj,2)*(pow(m1,4) + 6*pow(m1,2)*v1 + 3*pow(v1,2)) + 16*pow(m1,3)*pow(m2,3)*v12 +
          48*m1*pow(m2,3)*v1*v12 + 72*pow(m1,2)*pow(m2,2)*pow(v12,2) + 72*pow(m2,2)*v1*pow(v12,2) +
          96*m1*m2*pow(v12,3) + 24*pow(v12,4) + 6*
          (pow(m2,2)*(pow(m1,4) + 6*pow(m1,2)*v1 + 3*pow(v1,2)) + 8*m1*m2*(pow(m1,2) + 3*v1)*v12 +
          12*(pow(m1,2) + v1)*pow(v12,2))*v2 + 3*(pow(m1,4) + 6*pow(m1,2)*v1 + 3*pow(v1,2))*pow(v2,2) +
          3*pow(hi,2)*(3*pow(hj,2) + pow(m2,4) + 6*pow(m2,2)*v2 + 3*pow(v2,2) - 6*hj*(pow(m2,2) + v2)) -
          6*hj*(8*pow(m1,3)*m2*v12 + 24*m1*m2*v1*v12 + pow(m1,4)*(pow(m2,2) + v2) +
          6*pow(m1,2)*(2*pow(v12,2) + v1*(pow(m2,2) + v2)) + 3*v1*(4*pow(v12,2) + v1*(pow(m2,2) + v2))) +
          6*hi*(-3*pow(hj,2)*(pow(m1,2) + v1) - pow(m2,2)*
          (pow(m2,2)*(pow(m1,2) + v1) + 8*m1*m2*v12 + 12*pow(v12,2)) -
          6*(pow(m2,2)*(pow(m1,2) + v1) + 4*m1*m2*v12 + 2*pow(v12,2))*v2 - 3*(pow(m1,2) + v1)*pow(v2,2) +
          6*hj*(4*m1*m2*v12 + 2*pow(v12,2) + pow(m1,2)*(pow(m2,2) + v2) + v1*(pow(m2,2) + v2))))/
            (64.*pow(hi,3)*pow(hj,3));
      }
      else{
        sumij = ((3*pow(hi,2) + pow(m1,4) + 6*pow(m1,2)*v1 + 3*pow(v1,2) - 6*hi*(pow(m1,2) + v1))*
          (3*pow(hj,2) + pow(m2,4) + 6*pow(m2,2)*v2 + 3*pow(v2,2) - 6*hj*(pow(m2,2) + v2)))/
            (64.*pow(hi,3)*pow(hj,3));
      }
      sum = sum + sumij;
    }
  }
  return sum;
}
