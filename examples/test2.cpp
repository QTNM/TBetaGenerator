//
// test aspects of TBeta
//

// std libs
#include <iostream>

// us
#include "TBetaGenerator.hpp"

int main() {
  double mnu = 1.e-4; // 0.1 eV
  double en  = 18.5; // keV
  double endp = TBeta::endAt(mnu, 1);

  std::cout.precision(10);
  // endpoints
  std::cout << "tbeta end at mnu=0.1, n=1: " << endp << std::endl;
  std::cout << "tbeta end at mnu=0.1, n=2: " << TBeta::endAt(mnu, 2) << std::endl;
  std::cout << "tbeta end at mnu=0.1, n=3: " << TBeta::endAt(mnu, 3) << std::endl;
  std::cout << std::endl;

  // corrections
  std::cout << "S-1 ppm at e=0.04; mnu=0.1: " << (TBeta::S(2, 0.04)-1.0)*1e6 << std::endl;
  std::cout << "S-1 ppm at e=0.1; mnu=0.1: " << (TBeta::S(2, 0.1)-1.0)*1e6 << std::endl;
  std::cout << "S-1 ppm at e=0.2; mnu=0.1: " << (TBeta::S(2, 0.2)-1.0)*1e6 << std::endl;
  std::cout << "S-1 ppm at e=0.3; mnu=0.1: " << (TBeta::S(2, 0.3)-1.0)*1e6 << std::endl;
  std::cout << "S-1 ppm at e=1; mnu=0.1: " << (TBeta::S(2, 1.0)-1.0)*1e6 << std::endl;
  std::cout << "S-1 ppm at e=5; mnu=0.1: " << (TBeta::S(2, 5.0)-1.0)*1e6 << std::endl;
  std::cout << "S-1 ppm at e=endp; mnu=0.1: " << (TBeta::S(2, endp)-1.0)*1e6 << std::endl;
  std::cout << std::endl;

  bool order = true;
  double mN  = 0.0;
  double eta = 0.0;

  std::cout << "tbeta dGdE at 0.2keV: " << TBeta::dGammadE(order, mnu, mN, eta, 0.2) << std::endl;
  std::cout << "tbeta dGdE at 1keV: " << TBeta::dGammadE(order, mnu, mN, eta, 1.0) << std::endl;
  std::cout << "tbeta dGdE at 5keV: " << TBeta::dGammadE(order, mnu, mN, eta, 5.0) << std::endl;
  std::cout << "tbeta dGdE at 18.5keV: " << TBeta::dGammadE(order, mnu, mN, eta, en) << std::endl;
  std::cout << "tbeta dGdE at endpoint: " << TBeta::dGammadE(order, mnu, mN, eta, endp) << std::endl;
  std::cout << "tbeta dGdE 1 meV beyond endpoint: " << TBeta::dGammadE(order, mnu, mN, eta, endp+1.e-6) << std::endl;
  std::cout << std::endl;

  std::cout << "tbeta full dGdE at 0.2keV: " << TBeta::dGammadEFull(order, mnu, mN, eta, 0.2) << std::endl;
  std::cout << "tbeta full dGdE at 1keV: " << TBeta::dGammadEFull(order, mnu, mN, eta, 1.0) << std::endl;
  std::cout << "tbeta full dGdE at 5keV: " << TBeta::dGammadEFull(order, mnu, mN, eta, 5.0) << std::endl;
  std::cout << "tbeta full dGdE at 18.5keV: " << TBeta::dGammadEFull(order, mnu, mN, eta, en) << std::endl;
  std::cout << "tbeta full dGdE at endpoint: " << TBeta::dGammadEFull(order, mnu, mN, eta, endp) << std::endl;
  std::cout << "tbeta full dGdE 1meV beyond endpoint: " << TBeta::dGammadEFull(order, mnu, mN, eta, endp+1.e-6) << std::endl;

  return 0;
}
