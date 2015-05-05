// MCTrackPerf.h

// Class to evaluate the performance of cluster finding. There are three steps:
// 1. Declare the MC particle(s).
// 2. Loop over SimHits to find the charge in each channel-tick bin.
// 3. Find the closest cluster for each MC track.

#include "ChannelHits.h"
namespace simb {
class MCParticle;
}
namespace sim {
class Cluster;
}

class MCTrackPerf {

public:

  typedef ChannelHits::Channel Channel;
  typedef ChannelHits::Tick    Tick;
  typedef ChannelHits::Signal  Signal;

public:

  // Default ctor.
  MCTrackPerf();

  // Ctor for a MC track.
  MCTrackPerf(const simb::MCParticle& par);

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
  const ChannelHits& hits() const;

  // Output stream.
  //   out - stream to insert output
  //   detail - 0 = single line
  //            1 = line for each channel
  //            2 = line for each hit
  //            3 = line for each tick
  std::ostream& print(std::ostream& out, int detail =1, std::string prefix ="") const;

  // Fill a Channel vs. tick histogram.
  //int fillChannelTickHist() const;

  // Return the distance to a cluster.
  //double clusterDistance(const recob::Cluster& clu) const;

private:

  int m_trackID;
  int m_pdg;
  int m_rpdg;
  ChannelHits m_hits;

};

std::ostream& operator<<(const MCTrackPerf& rhs, std::ostream& out);
