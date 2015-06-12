// TpcSegment.cxx

#include "TpcSegment.h"
#include <cmath>

//**********************************************************************

TpcSegment::TpcSegment()
: tpc(-1),
  enter(-1), exit(-1),
  x1(-1.e20), y1(-1.e20), z1(-1.e20), e1(1.e20),
  x2( 1.e20), y2( 1.e20), z2( 1.e20), e2(-1e20) { }

//**********************************************************************

TpcSegment::
TpcSegment(int atpc, float ax1, float ay1, float az1, float ae1, int aenter)
: tpc(atpc),
  enter(aenter), exit(0),
  x1(ax1), y1(ay1), z1(az1), e1(ae1),
  x2(ax1), y2(ay1), z2(az1), e2(ae1) { }

//**********************************************************************

void TpcSegment::addPoint(float ax, float ay, float az, float ae) {
  x2 = ax;
  y2 = ay;
  z2 = az;
  e2 = ae;
  double dx = x2 - x1;
  double dy = y2 - y1;
  double dz = z2 - z1;
  float dlen = sqrt(dx*dx + dy*dy + dz*dz);
  length += dlen;
}

//**********************************************************************
