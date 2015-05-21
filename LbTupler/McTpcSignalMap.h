// McTpcSignalMap.h

#ifndef McTpcSignalMap_H
#define McTpcSignalMap_H

// David Adams
// May 2015
//
// Class to hold a TPC signal map associated with an MC particle.

#include "TpcSignalMap.h"
namespace simb {
class MCParticle;
}
namespace sim {
class Cluster;
}

class McTpcSignalMap : public TpcSignalMap {

public:

  typedef TpcSignalMap::Index   Index;
  typedef TpcSignalMap::Channel Channel;
  typedef TpcSignalMap::Tick    Tick;
  typedef TpcSignalMap::Signal  Signal;

public:

  // Default ctor.
  McTpcSignalMap();

  // Ctor for a MC track.
  McTpcSignalMap(const simb::MCParticle& par, const GeoHelper* pgh =nullptr);

  // Add contributions from a SimChannel.
  int addSimChannel(const sim::SimChannel& sch);

  // Getters.
  unsigned int trackID() const;
  int pdg() const;
  int rpdg() const;

  // Output stream.
  //   out - stream to insert output
  //   detail - 0 = single line
  //            1 = line for each channel
  //            2 = line for each hit
  //            3 = line for each tick
  std::ostream& print(std::ostream& out, int detail =1, std::string prefix ="") const;

private:

  int m_trackID;
  int m_pdg;
  int m_rpdg;

};

std::ostream& operator<<(const McTpcSignalMap& rhs, std::ostream& out);

#endif
