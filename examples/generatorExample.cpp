//
// test building electron distribution
//

// std libs
#include <iostream>
#include <random>

// us
#include "TBetaGenerator.hpp"

// build generation functor
// operator calls eactly one 
// distribution function.
// Make different functor for 
// alternative distribution.
class betaGenerator {
  // data parameters
private:
  bool order_;
  double munu_;
  double mNu_;
  double mixing_;

public:
  // fix parameter at construction
  betaGenerator(bool o, double ml, double mn, double eta) :
    order_(o),
    munu_(ml),
    mNu_(mn),
    mixing_(eta) {} // with parameter constructor

  double operator() (double x) {
    return TBeta::dGammadE(order_, munu_, mNu_, mixing_, x);
  }
};

int main() {
  typedef std::piecewise_linear_distribution<double> pld_type;
  std::cout.precision(10);
  
  // beta decay parameter from G4 macro, for instance
  bool order = true;    // normal order
  double mnu = 1.0e-4;  // 0.1 eV neutrino mass [keV]
  double mN  = 0.0;     // no sterile neutrino
  double eta = 0.0;     // no mixing

  // distribution parameter
  int nw = 100; // nw - number of bins
  double lbound = 0.2; // lower energy bound [keV]
  double ubound = TBeta::endAt(mnu, 1); // max energy
  
  pld_type ed(nw, lbound, ubound, betaGenerator(order, mnu, mN, eta));

  std::cout<< "energy intervals:" << std::endl;
  for (double val : ed.intervals()) std::cout << val << " ";
  std::cout << std::endl;

  std::cout<< "energy densities:" << std::endl;
  for (double val : ed.densities()) std::cout << val << " ";
  std::cout << std::endl;
  return 0;
}
