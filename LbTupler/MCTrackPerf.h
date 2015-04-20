// MCTrackPerf.h

// Class to evaluate the performance of cluster finding. There are three steps:
// 1. Declare the MC particle(s).
// 2. Loop over SimHits to find the charge in each channel-tick bin.
// 3. Find the closest cluster for each MC track.

#include <string>
#include <vector>
#include <map>
#include <iostream>
namespace simb {
class MCParticle;
}
namespace sim {
class SimChannel;
}
namespace recob {
class Cluster;
}

class MCTrackPerf {

public:

  struct Hit {
    unsigned int tick1;
    unsigned int tick2;
    double charge;
    Hit();
    Hit(unsigned int atick1, unsigned int atick2, double acharge);
  };

public:

  // Value used for an invalid or undefined channel or tick.
  static unsigned int badIndex();

public:

  typedef std::map<unsigned int, double> TickMap;
  typedef std::map<unsigned int, TickMap> TickChannelMap;
  typedef std::vector<Hit> HitVector;
  typedef std::map<unsigned int, HitVector> HitChannelMap;

  // Default ctor.
  MCTrackPerf();

  // Ctor for a MC track.
  MCTrackPerf(const simb::MCParticle& par);

  // Add an energy deposit in a bin.
  int add(unsigned int chan, unsigned int tick, double signal);

  // Add contributions from a SimChannel.
  int addSimChannel(const sim::SimChannel& sch);

  // Build hits from the tick energy deposits in each channel.
  // A hit is a contiguous set of ticks.
  int buildHits();

  // Getters.
  unsigned int trackID() const;
  int pdg() const;
  int rpdg() const;
  const TickChannelMap& tickChargeMap() const;
  const HitChannelMap& hitChargeMap() const;

  // Range of data.
  unsigned int channelMin() const;
  unsigned int channelMax() const;
  unsigned int tickMin() const;
  unsigned int tickMax() const;

  // The size is the number of included channel-tick bins.
  unsigned int size() const;

  // Output stream.
  //   out - stream to insert output
  //   detail - 0 = no printing
  //            1 = single line
  //            2 = line for each channel
  std::ostream& print(std::ostream& out, int detail =1, std::string prefix ="") const;

  // Fill a Channel vs. tick histogram.
  //int fillChannelTickHist() const;

  // Return the distance to a cluster.
  //double clusterDistance(const recob::Cluster& clu) const;

private:

  int m_trackID;
  int m_pdg;
  int m_rpdg;
  TickChannelMap m_tickchg;  // m_tickchg[chan][tick] is the charge for (chan, tick)
  HitChannelMap m_hitchg;    // m_hitchg[chan][tick] is the charge for (chan, hit)
  unsigned int m_tickMin;
  unsigned int m_tickMax;

};

std::ostream& operator<<(const MCTrackPerf& rhs, std::ostream& out);
