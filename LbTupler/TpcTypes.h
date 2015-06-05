// TpcTypes.h

#ifndef TpcTypes_H
#define TpcTypes_H

#include "Range.h"

namespace tpc {

typedef unsigned int Channel;
typedef unsigned int Index;
typedef int Tick;

typedef Range<Tick> TickRange;

// Value used for an invalid or undefined channel or tick.
Channel badChannel();
Tick badTick();
Index badIndex();

}

#endif
