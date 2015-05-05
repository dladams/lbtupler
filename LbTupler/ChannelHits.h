// ChannelHits.h

// Map of hits associated with channels. A hit is defined as a contiguous
// set of ticks.

#include <string>
#include <vector>
#include <map>
#include <iostream>

namespace sim {
class SimChannel;
}

class ChannelHits {

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

  typedef std::map<Tick, double> SignalTickMap;
  typedef std::map<Channel, SignalTickMap> TickChannelMap;
  typedef std::vector<Hit> HitVector;
  typedef std::map<Channel, HitVector> HitChannelMap;

  // Default ctor.
  ChannelHits();

  // Add an energy deposit in a bin.
  int addSignal(Channel chan, Tick tick, Signal signal);

  // Add contributions from a SimChannel for track tid.
  int addSimChannel(const sim::SimChannel& sch, unsigned int tid);

  // Build hits from the tick energy deposits in each channel.
  // A hit is a contiguous set of ticks.
  int buildHits();

  // Getters.
  const TickChannelMap& tickSignalMap() const;
  const HitChannelMap& hitSignalMap() const;

  // Range of data.
  Channel channelMin() const;
  Channel channelMax() const;
  Tick tickMin() const;
  Tick tickMax() const;

  // The number of included channels.
  unsigned int channelCount() const;

  // The number of included ticks.
  unsigned int size() const;
  unsigned int tickCount() const;

  // The number of hits.
  unsigned int hitCount() const;

  // Return the total signal summed over ticks and channels.
  Signal tickSignal() const;

  // Return the total signal summed over hits and channels.
  Signal hitSignal() const;

  // Output stream.
  //   out - stream to insert output
  //   detail - 0 = single line
  //            1 = line for each channel
  //            2 = line for each hit
  //            3 = line for each tick
  std::ostream& print(std::ostream& out, int detail =1,
                      std::string hdrprefix ="", std::string prefix ="") const;

  // Fill a Channel vs. tick histogram.
  //int fillChannelTickHist() const;

  // Return the distance to a cluster.
  //double clusterDistance(const recob::Cluster& clu) const;

private:

  TickChannelMap m_ticksig;  // m_ticksig[chan][tick] is the signal for (chan, tick)
  HitChannelMap m_hitsig;    // m_hitsig[chan][hit] is the hit for (chan, hit number)
  Tick m_tickMin;
  Tick m_tickMax;

};

std::ostream& operator<<(const ChannelHits& rhs, std::ostream& out);
