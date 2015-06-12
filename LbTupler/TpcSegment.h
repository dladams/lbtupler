// TpcSegment.h

#ifndef TpcSegment_H
#define TpcSegment_H

// David Adams
// June 2015
//
// Class to describe the segment of an MC particle track
// passing through a TPC.

#include <memory>
#include <vector>

class TpcSegment {

public:  // methods

  // Default ctor.
  TpcSegment();

  // Ctor for the first point.
  TpcSegment(int atpc, float ax1, float ay1, float az1, float ae1, int aenter);

  // Add a point to the segment. This becomes the last.
  void addPoint(float ax, float ay, float az, float ae);

public:  // data

  int tpc;     // TPC index
  int enter;   // 0/1 if the particle doesn't/does enter this TPC.
  int exit;    // 0/1 if the particle doesn't/does exit this TPC.
  float x1;    // First point
  float y1;
  float z1;
  float e1;
  float x2;    // Last point
  float y2;
  float z2;
  float e2;
  float length;

};

typedef std::shared_ptr<TpcSegment> TpcSegmentPtr;
typedef std::vector<TpcSegmentPtr> TpcSegmentVector;

#endif
