// TpcSignalMap.h

#ifndef TpcSignalMap_H
#define TpcSignalMap_H

// David Adams
// June 2015
//
// Objects of this class hold a map of float signals indexed by TPC channel
// and TDC tick. They also hold hits indexed by channel.

#include <string>
#include <vector>
#include <map>
#include <iostream>
#include <memory>
#include "art/Framework/Core/FindManyP.h"
#include "TpcTypes.h"
#include "TpcSegment.h"
#include "GeoHelper.h"

namespace simb {
class MCParticle;
}
namespace sim {
class SimChannel;
}
namespace recob {
class Hit;
}
class TH2;

class TpcSignalMap {

public:

  typedef std::string Name;
  typedef tpc::Index Index;
  typedef tpc::Channel Channel;
  typedef tpc::Tick Tick;
  typedef tpc::TickRange TickRange;
  typedef double Signal;
  typedef std::vector<art::Ptr<recob::Hit>> AssociatedHits;
  typedef std::shared_ptr<TpcSignalMap> TpcSignalMapPtr;
  typedef std::vector<TpcSignalMapPtr> TpcSignalMapVector;

public:

  struct McInfo {
    const simb::MCParticle* ppar;
    int trackID;
    int pdg;
    int rpdg;
    McInfo(const simb::MCParticle& par);
  };

  struct Hit {
    TickRange ticks;
    Tick tick2;
    Signal signal;
    Hit();
    Hit(unsigned int atick1, unsigned int atick2, double asignal);
  };

public:

  typedef std::shared_ptr<McInfo> McInfoPtr;
  typedef std::map<Tick, double> SignalTickMap;
  typedef std::map<Channel, SignalTickMap> TickChannelMap;
  typedef std::map<Index, TickChannelMap> TpcTickChannelMap;
  typedef std::vector<Hit> HitVector;
  typedef std::map<Channel, HitVector> HitChannelMap;
  typedef std::vector<Index> IndexVector;

  // Default ctor.
  TpcSignalMap();

  // Ctor from name and geometry helper.
  // If this is used, then each channel (and thus each signal and hit)
  // is assigned to a ROP (APA readout plane).
  TpcSignalMap(Name name, const GeoHelper* pgh);

  // Ctor adding an MC particle.
  TpcSignalMap(std::string name, const simb::MCParticle& par, const GeoHelper* pgh =nullptr);

  // Copy keeping only the signals for a given range of channels.
  // Dtor.
  virtual ~TpcSignalMap();

  // Add a signal (energy deposit or ADC count) in a bin.
  // If itpc is valid, then the singal is assigned to that tpc.
  int addSignal(Channel chan, Tick tick, Signal signal, Index itpc =GeoHelper::badIndex());

  // Add contributions from a SimChannel for track tid.
  int addSimChannel(const sim::SimChannel& sch, unsigned int tid, bool usetpc);

  // Add contributions from a SimChannel taking track ID from the MC info.
  // If usetpc is true, then the contributions are assigned to a TPC.
  int addSimChannel(const sim::SimChannel& sch, bool usetpc);

  // Add a recob::Hit and its signals.
  int addHit(const recob::Hit& hit, int verbose =0);

  // Add the hits and signals associated with a recob::Cluster.
  int addCluster(const AssociatedHits& hits, int verbose =0);

  // Attach a segment to this object.
  int addSegment(TpcSegmentPtr pseg);

  // Add segments.
  int addSegments(const TpcSegmentVector& segs);

  // Copy the segments from another signal map;
  int copySegments(const TpcSignalMap& tsm);

  // Copy the segments from signal maps with the same MC track as this.
  int copySegments(const TpcSignalMapVector& tsms);

  // Build hits from the tick signals in each channel.
  // A hit is a contiguous set of ticks.
  int buildHits();

  // Getters.
  std::string name() const { return m_name; }
  const GeoHelper* geometryHelper() const { return m_pgh; }
  IndexVector tpcs() const;
  const TickChannelMap& tickSignalMap(Index itpc) const;
  const HitChannelMap& hitSignalMap() const;

  // Return the TPC indices that are common whith ain input list.
  IndexVector sharedTpcs(const IndexVector& intpcs) const;

  // Range of data.
  Channel channelMin() const;
  Channel channelMax() const;
  Tick tickMin() const;
  Tick tickMax() const;
  TickRange tickRange() const;

  // The number of included channels, summing over all TPCs.
  unsigned int channelCount() const;

  // The number of included ticks.
  unsigned int tickCount() const;

  // The number of channel-tick bins.
  unsigned int binCount() const;
  unsigned int size() const;

  // The number of ROPs with signals.
  unsigned int ropCount() const;

  // The number of hits.
  unsigned int hitCount() const;

  // Return the total signal summed over ticks and channels.
  Signal tickSignal() const;

  // Return the total signal summed over hits and channels.
  Signal hitSignal() const;

  // Read the number of filled channel-tick bins for a ROP.
  Index ropNbin(Index irop) const;

  // Return if there is MC info associated with this object.
  bool haveMcinfo() const;

  // Return the MC info associated with this object. Null for none.
  const McInfo* mcinfo() const;

  // Set the ROP for this set of signals.
  void setRop(Index rop);

  // Return if this set of signals is assigned to a ROP.
  bool haveRop() const;

  // Return the ROP to which this set of signals is assigned.
  Index rop() const;

  // Return the segments.
  const TpcSegmentVector& segments() const;

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

  // Split this object along ROP boundaries, i.e. one object for each ROP
  // with signals. The new object are appended to tsms.
  // If splitByTpc is true, then a separate map is created for each TPC
  // in each ROP.
  int splitByRop(TpcSignalMapVector& tsms, bool splitByTpc =false) const;

private:

  std::string m_name;             // Name.
  const GeoHelper* m_pgh;         // Geometry helper maps channels to ROPs.
  TpcTickChannelMap m_tpcticksig; // m_ticksig[itpc][chan][tick] is the signal for (itpc, chan, tick)
  HitChannelMap m_hitsig;         // m_hitsig[chan][hit] is the hit for (chan, hit number)
  TickRange m_tickRange;          // Range of ticks covered by this signal map.
  IndexVector m_ropnbin;          // Number of filled channel-tick bins for each ROP
  McInfoPtr m_pmci;               // Managing pointer to MC info.
  Index m_rop;                    // ROP for this object (badIndex() if not defined)
  TpcSegmentVector m_segments;    // Attached segments.

};

std::ostream& operator<<(const TpcSignalMap& rhs, std::ostream& out);

#endif
