# Tritium Beta Decay Event Generator

Authors: J. Canning, F. Deppisch and W. Pei (UCL, 2022) in Mathematica for the QTNM project.

Translation to C++ for Geant4 applications: Y. Ramachers (University of Warwick, 2022)

## Notes 

This is a header-only library of simple functions. Examples testing the functions
can be found in the examples folder.

## Build instruction

Requirement: CMake 3.1x, C++11 compiler

At Warwick, SCRTP, use cvmfs as the easiest environment setup (with bash):

source /cvmfs/sft.cern.ch/lcg/views/LCG_99/x86_64-centos7-gcc10-opt/setup.sh

which sets up Geant4 10.7 and GCC10 on a CentOS7 background. Just create a 'build' 
directory, then 

cd build; cmake ..; make

and run in the build directory, for instance as 

./test1.exe
