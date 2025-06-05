//
// test building electron distribution
//

// std libs
#include <iostream>
#include <random>
#include <chrono>

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

  std::cout << "endpoint energy mnu=0: " << TBeta::endpoint(0.0) << std::endl;
  std::cout << "endAt energy mnu=0: " << TBeta::endAt(0.0,1) << std::endl;
  std::cout << "energy intervals:" << std::endl;
  for (double val : ed.intervals()) std::cout << val << " ";
  std::cout << std::endl;

  std::cout<< "energy densities:" << std::endl;
  for (double val : ed.densities()) std::cout << val << " ";
  std::cout << std::endl;
  std::cout << std::endl;

  std::cout<< "draw 2 random numbers from the distribution:" << std::endl;
  std::cout<< "These are electron energies in [keV]." << std::endl;
  unsigned int seed = std::chrono::system_clock::now().time_since_epoch().count();
  std::default_random_engine generator(seed);
  std::cout << "First rnd=" << ed(generator) 
	    << "; 2nd=" << ed(generator) << std::endl;

  return 0;
}
