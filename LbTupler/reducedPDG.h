// reducedPDG.h

#ifndef reducedPDG_H
#define reducedPDG_H

// David Adams
// Reduced PDG codes.
//
//   1 - electron
//   2 - muon
//   3 - proton
//   4 - charged pion
//   6 - gamma
//   7 - neutron
//   8 - neutrino
//   9 - pi0
//  11 - PDG > 1.e9 (nucleus)
//  12 - Everything else
//
// Returns opposite sign for charged antiparticle.

int reducedPDG(int pdg);

#endif
