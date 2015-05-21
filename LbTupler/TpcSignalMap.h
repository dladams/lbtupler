// TpcSignalMap.h

// Map of hits associated with channels. A hit is defined as a contiguous
// set of ticks.

#include <string>
#include <vector>
#include <map>
#include <iostream>
#include "art/Framework/Core/FindManyP.h"
#include "Range.h"

namespace sim {
class SimChannel;
}
namespace recob {
class Hit;
}
class TH2;
class GeoHelper;

class TpcSignalMap {

public:

  typedef unsigned int Channel;
  typedef int Tick;
  typedef Range<Tick> TickRange;
  typedef unsigned int Index;
  typedef double Signal;
  typedef std::vector<art::Ptr<recob::Hit>> AssociatedHits;

public:

  struct Hit {
    TickRange ticks;
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
  typedef std::vector<Index> IndexVector;

  // Default ctor.
  TpcSignalMap();

  // Ctor from a geometry helper.
  // If this is used, then each channel (and thus each signal and hit)
  // is assigned to a ROP (APA readout plane).
  TpcSignalMap(const GeoHelper* pgh);

  // Add a signal (energy deposit or ADC count) in a bin.
  int addSignal(Channel chan, Tick tick, Signal signal);

  // Add contributions from a SimChannel for track tid.
  int addSimChannel(const sim::SimChannel& sch, unsigned int tid);

  // Add a recob::Hit and its signals.
  int addHit(const recob::Hit& hit, int verbose =0);

  // Add the hits and signals associated with a recob::Cluster.
  int addCluster(const AssociatedHits& hits, int verbose =0);

  // Build hits from the tick signals in each channel.
  // A hit is a contiguous set of ticks.
  int buildHits();

  // Getters.
  const GeoHelper* geometryHelper() const { return m_pgh; }
  const TickChannelMap& tickSignalMap() const;
  const HitChannelMap& hitSignalMap() const;

  // Range of data.
  Channel channelMin() const;
  Channel channelMax() const;
  Tick tickMin() const;
  Tick tickMax() const;
  TickRange tickRange() const;

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

  // Read the number of filled channel-tick bins for a ROP.
  Index ropNbin(Index irop) const;

  // Output stream.
  //   out - stream to insert output
  //   detail -  0 = single line
  //             1 = line for each channel
  //             2 = line for each hit
  //             3 = line for each tick
  //            10+X = also print ROP counts
  std::ostream& print(std::ostream& out, int detail =1,
                      std::string hdrprefix ="", std::string prefix ="") const;

  // Fill a Channel vs. tick histogram.
  int fillChannelTickHist(TH2* ph) const;

  // Fill a Channel vs. tick histogram for a given ROP.
  int fillRopChannelTickHist(TH2* ph, Index irop) const;

  // Return the distance to a cluster.
  //double clusterDistance(const recob::Cluster& clu) const;

private:

  const GeoHelper* m_pgh;    // Geometry helper maps channels to ROPs.
  TickChannelMap m_ticksig;  // m_ticksig[chan][tick] is the signal for (chan, tick)
  HitChannelMap m_hitsig;    // m_hitsig[chan][hit] is the hit for (chan, hit number)
  TickRange m_tickRange;
  IndexVector m_ropnbin;

};

std::ostream& operator<<(const TpcSignalMap& rhs, std::ostream& out);