// intProcess.cxx

#include "intProcess.h"

//**********************************************************************

int intProcess(std::string sproc) {
  if ( sproc == "primary"              ) return  0;
  if ( sproc == "Decay"                ) return  1;
  if ( sproc == "muMinusCaptureAtRest" ) return  3;
  if ( sproc == "muPlusCaptureAtRest"  ) return  4;
  if ( sproc == "nCapture"             ) return  5;
  if ( sproc == "CoulombScat"          ) return 10;
  if ( sproc == "ProtonInelastic"      ) return 11;
  if ( sproc == "NeutronInelastic"     ) return 12;
  if ( sproc == "hadElastic"           ) return 13;
  if ( sproc == "PositronNuclear"      ) return 14;
  return -1;
}

//**********************************************************************
