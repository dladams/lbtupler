// reducedPDG.cxx

#include "reducedPDG.h"
#include <cmath>

using std::abs;

int reducedPDG(int pdg) {
  int apdg = abs(pdg);
  if ( pdg ==    11 ) return  1;
  if ( pdg ==   -11 ) return -1;
  if ( pdg ==    13 ) return  2;
  if ( pdg ==   -13 ) return -2;
  if ( pdg ==  2212 ) return  3;
  if ( pdg == -2212 ) return -3;
  if ( pdg == 211 ) return 4;
  if ( pdg == -211 ) return -4;
  if ( pdg == 22 ) return 6;
  if ( pdg ==  2112 ) return 7;
  if ( apdg == 12 ) return 8;
  if ( apdg == 14 ) return 8;
  if ( apdg == 16 ) return 8;
  if ( pdg ==  111 ) return 9;
  if ( pdg > 1000000 ) return 11;
  return -12;
}
