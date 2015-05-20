// MCTrackPerf.h

#ifndef MCTrackPerf_H
#define MCTrackPerf_H

// David Adams
// May 2015
//
// Class to hold the signals and hits associated with a MC particle or track.

#include "TpcSignalMap.h"
namespace simb {
class MCParticle;
}
namespace sim {
class Cluster;
}

class MCTrackPerf {

public:

  typedef TpcSignalMap::Index   Index;
  typedef TpcSignalMap::Channel Channel;
  typedef TpcSignalMap::Tick    Tick;
  typedef TpcSignalMap::Signal  Signal;

public:

  // Default ctor.
  MCTrackPerf();

  // Ctor for a MC track.
  MCTrackPerf(const simb::MCParticle& par, const GeoHelper* pgh =nullptr);

  // Add an energy deposit in a bin.
  int addSignal(Channel chan, Tick tick, Signal signal);

  // Add contributions from a SimChannel.
  int addSimChannel(const sim::SimChannel& sch);

  // Build hits from the tick energy deposits in each channel.
  // A hit is a contiguous set of ticks.
  int buildHits();

  // Getters.
  unsigned int trackID() const;
  int pdg() const;
  int rpdg() const;
  const TpcSignalMap& hits() const;

  // Output stream.
  //   out - stream to insert output
  //   detail - 0 = single line
  //            1 = line for each channel
  //            2 = line for each hit
  //            3 = line for each tick
  std::ostream& print(std::ostream& out, int detail =1, std::string prefix ="") const;

  // Fill a Channel vs. tick histogram.
  int fillRopChannelTickHist(TH2* ph, Index irop) const;

  // Return the distance to a cluster.
  //double clusterDistance(const recob::Cluster& clu) const;

private:

  int m_trackID;
  int m_pdg;
  int m_rpdg;
  TpcSignalMap m_hits;

};

std::ostream& operator<<(const MCTrackPerf& rhs, std::ostream& out);

#endif
