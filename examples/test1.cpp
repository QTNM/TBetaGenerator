//
// test aspects of TBeta
//

// std libs
#include <iostream>
#include <complex>

// us
#include "TBetaGenerator.hpp"

double testfunc1(double x, double par1, double par2) {
  // par1,2 not used but required for interface to Simpson
  // in library; in interval (0,1) should be answer=0.25.
  return x*x*x;
}

double testfunc2(double x, double par1, double par2) {
  // par1,2 not used but required for interface to Simpson
  // in library; in interval (1,100) should be answer=ln(100).
  if (x<=0.0) return 0.0;
  return 1.0/x;
}

int main() {
  std::cout.precision(10);
  std::cout << "tbeta pi = " << TBeta::pi << std::endl;

  // test maths functions
  std::cout << "TBeta heavyside(1): "  << TBeta::heavyside(1.0)  << std::endl;
  std::cout << "TBeta heavyside(0): "  << TBeta::heavyside(0.0)  << std::endl;
  std::cout << "TBeta heavyside(-2): " << TBeta::heavyside(-2.0) << std::endl;

  std::complex<double> z(0.5, 0); // real only; answer = sqrt(pi)
  std::cout << "TBeta Gamma(0.5): " << TBeta::Gamma(z) << std::endl;
  z = std::complex<double>(2.5, 1); // complex; answer see scipy.special.gamma manual.
  std::cout << "TBeta Gamma(2.5+i): " << TBeta::Gamma(z) << std::endl;

  // check integration with testfunc from above
  std::function<double(double,double,double)> fn  = testfunc1;
  TBeta::SimpsonIntegral sint(fn, 0.0, 0.0); // set dummy par1,2 to 0
  std::cout << "check integral (0.25 analytically): " << sint.integrate(0.0, 1.0, 100) << std::endl;

  std::function<double(double,double,double)> fn2 = testfunc2;
  TBeta::SimpsonIntegral sint2(fn2, 0.0, 0.0); // set dummy par1,2 to 0
  std::cout << "check integral (ln(100) analytically): " << sint2.integrate(1.0, 100.0, 1000) << std::endl;

  return 0;
}
