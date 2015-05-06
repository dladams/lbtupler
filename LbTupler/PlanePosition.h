// PlanePosition.h

#ifndef PlanePosition_H
#define PlanePosition_H

// David Adams
// May 2015
//
// Struct that holds a position on a LAr TPC readout plane in LAr detector coordinates.

#include <vector>

class PlanePosition {

public:

  typedef unsigned int Index;

public:

  // Data.
  Index plane;          // plane # in the TPC
  Index rop;            // Global readout plane index.
  Index channel;        // Global channel.
  Index ropchannel;     // Channel in the readout plane.
  double tick;          // TDC tick float (t + x/v_drift)/t_bin
  int itick;            // TDC tick.
  bool valid;           // True if this is a valid plane position

  // Default ctor.
  PlanePosition() : plane(0), rop(0), channel(0), ropchannel(0), tick(0), valid(false) { }

};

// Vector of plane positions, e.g. positions for all planes in a TPC.
typedef std::vector<PlanePosition> PlanePositionVector;

#endif
