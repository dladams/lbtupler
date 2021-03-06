// TpcSignalMatcher

#ifndef TpcSignalMatcher_H
#define TpcSignalMatcher_H

#include <vector>
#include <iosfwd>
#include "TpcSignalMap.h"

// David Adams
// May 2015
//
// Class to match entries in two vectors of TpcSignalMap objects.

class TpcSignalMatcher {

public:

  typedef unsigned int Index;
  typedef std::vector<Index> IndexVector;
  typedef std::vector<float> FloatVector;
  typedef TpcSignalMap T1;
  typedef TpcSignalMap T2;
  typedef std::shared_ptr<T1> P1;
  typedef std::vector<P1> C1;
  typedef std::shared_ptr<T2> P2;
  typedef std::vector<P2> C2;
  typedef double (*Distance)(const T1&, const T2&);

public:

  // Ctor.
  //  c1 - reference collection
  //  c2 - match collection
  //  matchByRop - If true, only objects with the same ROP are matched. If so,
  //               all objects must have an assigned ROP.
  TpcSignalMatcher(const C1& c1, const C2& c2, bool matchByRop, int dbg =0);

  // Vector of reference objects.
  const C1& referenceVector() const { return m_c1; }

  // Vector of matched objects.
  const C2& matchVector() const { return m_c2; }

  // Distance metric.
  Distance distance() const;

  // Maximum allowed distance.
  double maxDistance() const;

  // # of entries (ref-match pairs)
  unsigned int size() const;

  // Rec index for each ref index.
  int matchIndex(Index iref) const;
  float matchDistance(Index iref) const;

  // Return a string describing the match.
  //   opt = 0: Index pairs
  //         1: simple graphical
  std::string show(int opt) const;

  // Print. Same options as show.
  std::ostream& print(std::ostream& out, int iopt) const;

private:

  const C1& m_c1;
  const C2& m_c2;
  IndexVector m_matchIndex;
  FloatVector m_matchDistance;

};

#endif
