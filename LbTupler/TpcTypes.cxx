// TpcTypes.cxx

#include "TpcTypes.h"
#include <limits>

//**********************************************************************

tpc::Channel tpc::badChannel() {
  return std::numeric_limits<Channel>::max();
}

//**********************************************************************

tpc::Tick tpc::badTick() {
  return std::numeric_limits<Tick>::max();
}

//**********************************************************************

tpc::Index tpc::badIndex() {
  return std::numeric_limits<Index>::max();
}

//**********************************************************************

tpc::Index tpc::badIndex2() {
  return std::numeric_limits<Index>::max() - 1;
}

//**********************************************************************

