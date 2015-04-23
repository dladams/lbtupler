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

  typedef unsigned int Channel;
  typedef int Tick;
  typedef double Signal;

public:

  struct Hit {
    Tick tick1;
    Tick tick2;
    Signal signal;
    Hit();
    Hit(unsigned int atick1, unsigned int atick2, double asignal);
  };

public:

  // Value used for an invalid or undefined channel or tick.
  static Channel badChannel();
  static Tick badTick();

public:

  typedef std::map<Tick, double> TickMap;
  typedef std::map<Channel, TickMap> TickChannelMap;
  typedef std::vector<Hit> HitVector;
  typedef std::map<Channel, HitVector> HitChannelMap;

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
  const TickChannelMap& tickSignalMap() const;
  const HitChannelMap& hitSignalMap() const;

  // Range of data.
  Channel channelMin() const;
  Channel channelMax() const;
  Tick tickMin() const;
  Tick tickMax() const;

  // The size is the number of included channel-tick bins.
  unsigned int size() const;

  // The number of included channels.
  unsigned int channelCount() const;

  // Return the total signal in the channel.
  Signal tickSignal() const;
  Signal hitSignal() const;

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
  TickChannelMap m_ticksig;  // m_ticksig[chan][tick] is the signal for (chan, tick)
  HitChannelMap m_hitsig;    // m_hitsig[chan][hit] is the hit for (chan, hit number)
  Tick m_tickMin;
  Tick m_tickMax;

};

std::ostream& operator<<(const MCTrackPerf& rhs, std::ostream& out);
