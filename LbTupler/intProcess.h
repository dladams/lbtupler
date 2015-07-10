// intProcess.h

#ifndef intProcess_H
#define intProcess_H

// David Adams
// March 2015
//
// Converts a NUTOOLS McParticle process to an integer.
//
//  0 - primary
//  1 - Decay
//  3 - muMinusCaptureAtRest
//  4 - muPlusCaptureAtRest
//  5 - nCapture
// 10 - CoulombScat
// 11 - ProtonInelastic
// 12 - NeutronInelastic
// 13 - hadElastic
// 14 - PositronNuclear
// -1 - Anything else

#include <string>

int intProcess(std::string sproc);

#endif
