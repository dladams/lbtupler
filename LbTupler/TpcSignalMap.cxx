// TpcSignalMap.cxx

#include "TpcSignalMap.h"
#include <iomanip>
#include <sstream>
#include "TH2.h"
#include "SimulationBase/MCParticle.h"
#include "Simulation/SimChannel.h"
#include "RecoBase/Hit.h"
#include "Geometry/Geometry.h"
#include "reducedPDG.h"

using std::string;
using std::cout;
using std::endl;
using std::setw;
using std::ostream;
using std::ostringstream;
using simb::MCParticle;
using sim::SimChannel;
using tpc::badIndex;
using tpc::badChannel;
using tpc::badTick;

typedef TpcSignalMap::Tick             Tick;
typedef TpcSignalMap::TickRange        TickRange;
typedef TpcSignalMap::Channel          Channel;
typedef TpcSignalMap::Signal           Signal;
typedef TpcSignalMap::Index            Index;
typedef TpcSignalMap::IndexVector      IndexVector;
typedef TpcSignalMap::IndexPairVector  IndexPairVector;

//**********************************************************************
// Local definitions.
//**********************************************************************

namespace {
int dbg() { return 0; }
}

//**********************************************************************
// Sub class methods.
//**********************************************************************

TpcSignalMap::McInfo::McInfo(const simb::MCParticle& par)
: ppar(&par),
  trackID(par.TrackId()),
  pdg(par.PdgCode()),
  rpdg(reducedPDG(pdg)) { }

//**********************************************************************

TpcSignalMap::Hit::Hit()
: ticks(0, -1) { }

//**********************************************************************

TpcSignalMap::Hit::Hit(unsigned int atick1, unsigned int atick2, double asignal)
: ticks(atick1, atick2), signal(asignal) { }

//**********************************************************************
// Main class.
//**********************************************************************

TpcSignalMap::TpcSignalMap()
: m_pgh(nullptr),
  m_usetpc(false),
  m_tickRange(std::numeric_limits<Tick>::max(), std::numeric_limits<Tick>::min()),
  m_rop(badIndex()) { }

//**********************************************************************

TpcSignalMap::TpcSignalMap(string name, const GeoHelper* pgh, bool ausetpc)
: m_name(name),
  m_pgh(pgh),
  m_usetpc(ausetpc),
  m_tickRange(std::numeric_limits<Tick>::max(), std::numeric_limits<Tick>::min()),
  m_rop(badIndex()) {
  if ( m_pgh != nullptr ) {
    m_ropnbin.resize(m_pgh->nrop(), 0);
  }
}

//**********************************************************************

TpcSignalMap::TpcSignalMap(string name, const simb::MCParticle& par, const GeoHelper* pgh, bool ausetpc)
: m_name(name),
  m_pgh(pgh),
  m_usetpc(ausetpc),
  m_tickRange(std::numeric_limits<Tick>::max(), std::numeric_limits<Tick>::min()),
  m_pmci(new McInfo(par)),
  m_rop(badIndex()) {
  if ( m_pgh != nullptr ) {
    m_ropnbin.resize(m_pgh->nrop(), 0);
  }
}

//**********************************************************************

TpcSignalMap::~TpcSignalMap() { }

//**********************************************************************

int TpcSignalMap::addSignal(Channel chan, Tick tick, Signal signal, Index itpcin) {
  const string myname = "TpcSignalMap::add: ";
  if ( dbg() ) std::cout << myname
                         << "Channel " << chan << ", tick " << tick
                         << " has signal " << signal << endl;
  Index itpc = itpcin;
  if ( usetpc() && itpc == badIndex() ) {
    cout << myname << "ERROR: Ignoring request to add TPC signal with invalid TPC." << endl;
    return 1;
  }
  // Ignore the TPC if object is not in the usetpc state.
  if ( ! usetpc() ) itpc = badIndex();
  TickChannelMap& ticksigmap = m_tpcticksig[itpc];
  if ( ticksigmap.find(chan) == ticksigmap.end() ) ticksigmap[chan] = SignalTickMap();
  SignalTickMap& ticksig = ticksigmap[chan];
  if ( ticksig.find(tick) == ticksig.end() ) {
    ticksig[tick] = signal;
    if ( m_pgh != nullptr ) {
      Index irop = m_pgh->channelRop(chan);
      if ( irop < m_ropnbin.size() ) {
        ++m_ropnbin[irop];
      } else {
        cout << myname << "Channel " << chan << " does not have a ROP." << endl;
      }
    }
  } else {
    ticksig[tick] += signal;
  }
  Tick& tick1 = m_tickRange.first();
  Tick& tick2 = m_tickRange.last();
  if ( tick < tick1 ) tick1 = tick;
  if ( tick > tick2 ) tick2 = tick;
  return 0;
}

//**********************************************************************

int TpcSignalMap::addSimChannel(const SimChannel& simchan, unsigned int tid) {
  const string myname = "TpcSignalMap::addSimChannel: ";
  if ( usetpc() && geometryHelper() == nullptr ) {
    cout << myname << "ERROR: Geometry is required to find the TPC." << endl;
    return 1;
  }
  Channel chan = simchan.Channel();
  auto const& tickides = simchan.TDCIDEMap();
  if ( dbg() ) std::cout << myname << "Track " << tid << ", channel " << chan << endl;
  for ( auto const& tickide : tickides ) {
    Tick tick = tickide.first;
    // Protect against negative ticks.
    if ( tick > 63000 ) continue;
    auto& ides = tickide.second;
    if ( dbg() ) std::cout << myname << "  Adding Tick=" << tick << " with IDE count " << ides.size() << endl;
    for ( auto& ide : ides ) {
      if ( dbg() > 1 ) std::cout << myname << "    Adding TrackID=" << ide.trackID << endl;
      if ( abs(ide.trackID) == tid ) {
        //Signal sig = ide.numElectrons;
        Signal sig = ide.energy;
        if ( dbg() > 1 ) std::cout << myname << "    Signal=" << sig << endl;
        Index itpc = badIndex();
        if ( usetpc() ) {
          double pos[3];
          pos[0] = ide.x;
          pos[1] = ide.y;
          pos[2] = ide.z;
          geo::TPCID tpcid = geometryHelper()->geometry()->FindTPCAtPosition(pos);
          if ( ! tpcid.isValid ) {
            std::cout << myname << "WARNING: IDE is not inside a TPC!" << endl;
            continue;
          }
          if ( tpcid.Cryostat != 0 ) {
            std::cout << myname << "WARNING: IDE is not in cryostat 0!" << endl;
            continue;
          }
          itpc = tpcid.TPC;
        }
        addSignal(chan, tick, sig, itpc);
      } else {
        if ( dbg() > 1 ) std::cout << myname << "  Skipping track " << ide.trackID << endl;
      }
    }  // End loop over IDE's for this tick
  }  // End loop over ticks for this sim channel
  return 0;
}

//**********************************************************************

int TpcSignalMap::addSimChannel(const sim::SimChannel& sch) {
  if ( mcinfo() == nullptr ) return -1;
  return addSimChannel(sch, mcinfo()->trackID);
}

//**********************************************************************

int TpcSignalMap::addHit(const recob::Hit& rhit, int dbg) {
  const string myname = "TpcSignalMap::addHit: ";
  Index chan = rhit.Channel();
  Tick tick1 = rhit.StartTick();
  Tick tick2 = rhit.EndTick();
  Signal sig = rhit.Integral();
  double tick0 = rhit.PeakTime();
  double sigma = rhit.RMS();
  if ( dbg > 0 ) cout << myname << chan << ": " << tick1 << "-" << tick2
                      << " peak=" << tick0
                      << " +/- " << sigma
                      << " signal=" << sig
                      << endl;
  Index itpc = badIndex();
  if ( usetpc() ) {
    cout << myname << "ERROR: Unable to extract TPC from hit." << endl;
    abort();
  }
  HitChannelMap& hitsigmap = m_tpchitsig[itpc];
  HitVector& hits = hitsigmap[chan];
  hits.emplace(hits.end(), tick1, tick2, sig);
  // Find the signal fraction in each tick bin and add it to the signal map.
  double sigsum = 0.0;
  Tick ntick = tick2 - tick1 + 1;
  if ( ntick == 1 ) {
    // All signal in one bin.
    sigsum += sig;
    addSignal(chan, tick1, sig);
  } else if ( ntick == 2 ) {
    // Linear division between bins.
    double sigbin = (tick0 - tick1 -0.5)*sig;
    sigsum += sigbin;
    addSignal(chan, tick1, sigbin);
    sigbin = sig - sigbin;
    sigsum += sigbin;
    addSignal(chan, tick2, sigbin);
  } else {
    // Gaussian distribution of signal.
    if ( sigma <= 0.0 ) {
      cout << myname << "Input hit has invalid sigma: " << sigma << endl;
      return 1;
    }
    double fac = 1.0/(sqrt(2.0)*sigma);
    for ( Tick tick=tick1; tick<=tick2; ++tick ) {
      double x1 = fac*fabs(tick-tick0);
      double erf1 = TMath::Erf(x1);
      double x2 = fac*fabs(tick+1.0-tick0);
      double erf2 = TMath::Erf(x2);
      bool inpeakbin = tick0 > tick && tick0 < tick+1.0;
      double erfbin = inpeakbin ? erf1 + erf2 : fabs(erf2 - erf1);
      double sigbin = 0.5*erfbin*sig;
      sigsum += sigbin;
      if ( dbg > 2 ) cout << myname << "  " << tick << ": " << sigbin << endl;
      addSignal(chan, tick, sigbin);
    }
  }
  if ( dbg > 1 ) cout << myname << "  (Bin sum)/integral = " << sigsum << "/" << sig << " = " << sigsum/sig << endl;
  return 0;
}

//**********************************************************************

int TpcSignalMap::addCluster(const AssociatedHits& hits, int dbg) {
  const string myname = "TpcSignalMap::addCluster: ";
  for ( const auto& phit : hits ) {
    if ( phit.isNull() ) {
      cout << myname << "WARNING: Hit is missing." << endl;
      continue;
    }
    const recob::Hit& hit = *phit;
    addHit(hit, dbg);
  }
  return 0;
}

//**********************************************************************

int TpcSignalMap::addSegment(TpcSegmentPtr pseg) {
  const string myname = "TpcSignalMap::addSegment: ";
  if ( dbg() ) cout << myname << "Adding a segment." << endl;
  m_segments.push_back(pseg);
  return 0;
}

//**********************************************************************

int TpcSignalMap::addSegments(const TpcSegmentVector& segs) {
  const string myname = "TpcSignalMap::addSegments: ";
  m_segments.insert(m_segments.end(), segs.begin(), segs.end());
  return 0;
}

//**********************************************************************

int TpcSignalMap::copySegments(const TpcSignalMap& tsm) {
  const string myname = "TpcSignalMap::copySegments: ";
  return addSegments(tsm.segments());
}

//**********************************************************************

int TpcSignalMap::copySegments(const TpcSignalMapVector& tsms) {
  const string myname = "TpcSignalMap::copySegments: ";
  if ( ! haveMcinfo() ) return 1;
  int itrk = mcinfo()->trackID;
  for ( const TpcSignalMapPtr& ptsm : tsms ) {
    if ( ptsm->haveMcinfo() ) {
      if ( ptsm->mcinfo()->trackID == itrk ) {
        copySegments(*ptsm);
      }
    }
  }
  return 0;
}

//**********************************************************************

int TpcSignalMap::buildHits() {
  const string myname = "TpcSignalMap::buildHits: ";
  // Loop over TPCs.
  for ( auto itpcticksig : m_tpcticksig ) {
    Index itpc = itpcticksig.first;
    const TickChannelMap& ticksig = itpcticksig.second;
    // Loop over channels.
    for ( const auto ent : ticksig ) {
      Channel chan = ent.first;
      //const auto& ticks = ent.second;
      const SignalTickMap& ticks = ent.second;
      Signal hitsig = 0.0;
      Tick tick1 = badTick();
      Tick tick2 = badTick();
      HitChannelMap& hitsigmap = m_tpchitsig[itpc];
      HitVector& hits = hitsigmap[chan];
      if ( hits.size() ) {
        cout << myname << "ERROR: Channel hits already defined. Size = " << hits.size() << endl;
        return 1;
      }
      for ( auto ticksig : ticks ) {
        Tick tick = ticksig.first;
        Signal sig = ticksig.second;
        if ( tick1 != badTick() && tick != tick2+1 ) {
          hits.push_back(Hit(tick1, tick2, hitsig));
          tick1 = badTick();
          tick2 = badTick();
          hitsig = 0.0;
        }
        if ( tick1 == badTick() ) {
          tick1 = tick;
          tick2 = tick1;
        } else {
          tick2 = tick;
        }
        hitsig += sig;
      }
      hits.push_back(Hit(tick1, tick2, hitsig));
    }  // End loop over channels.
  }  // End loop over TPCs.
  return 0;
}

//**********************************************************************

IndexVector TpcSignalMap::tpcs() const {
  IndexVector out;
  for ( auto tpcticksig : m_tpcticksig ) {
    Index itpc = tpcticksig.first;
    out.push_back(itpc);
  }
  return out;
}

//**********************************************************************

const TpcSignalMap::TickChannelMap&
TpcSignalMap::tickSignalMap(Index itpc) const {
  static const TickChannelMap empty;
  auto its = m_tpcticksig.find(itpc);
  if ( its == m_tpcticksig.end() ) return empty;
  return its->second;
}

//**********************************************************************

const TpcSignalMap::HitChannelMap& TpcSignalMap::hitSignalMap(Index itpc) const {
  TpcHitChannelMap::const_iterator ithcm = m_tpchitsig.find(itpc);
  if ( ithcm == m_tpchitsig.end() ) {
    static HitChannelMap empty;
    return empty;
  }
  return ithcm->second;
}

//**********************************************************************

Index TpcSignalMap::ropNbin(Index irop) const {
  if ( irop >= m_ropnbin.size() ) return 0;
  return m_ropnbin[irop];
}

//**********************************************************************

bool TpcSignalMap::haveMcinfo() const {
  return m_pmci.get() != nullptr;
}

//**********************************************************************

const TpcSignalMap::McInfo* TpcSignalMap::mcinfo() const {
  return m_pmci.get();
}

//**********************************************************************

void TpcSignalMap::setRop(Index rop) {
  m_rop = rop;
}

//**********************************************************************

bool TpcSignalMap::haveRop() const {
  return m_rop != badIndex();
}

//**********************************************************************

Index TpcSignalMap::rop() const {
  return m_rop;
}

//**********************************************************************

int TpcSignalMap::check() const {
  bool havebadtpc_sig = m_tpcticksig.find(badIndex()) != m_tpcticksig.end();
  bool havebadtpc_hit = m_tpchitsig.find(badIndex()) != m_tpchitsig.end();
  unsigned int ntpc_sig = m_tpcticksig.size();
  if ( havebadtpc_sig ) --ntpc_sig;
  unsigned int ntpc_hit = m_tpchitsig.size();
  if ( havebadtpc_hit ) --ntpc_hit;
  if ( usetpc() ) {
    if ( havebadtpc_sig ) return 1;
    if ( havebadtpc_hit ) return 2;
  } else {
    if ( ntpc_sig ) return 3;
    if ( ntpc_hit ) return 4;
  }
  return 0;
}

//**********************************************************************

IndexVector TpcSignalMap::sharedTpcs(const IndexVector& intpcs) const {
  IndexVector out;
  IndexVector mytpcs = tpcs();
  for ( Index itpc : intpcs ) {
    if ( find(mytpcs.begin(), mytpcs.end(), itpc) != mytpcs.end() ) {
      out.push_back(itpc);
    }
  }
  return out;
}

//**********************************************************************

IndexPairVector TpcSignalMap::
sharedTpcPairs(const IndexVector& intpcs, bool same) const {
  IndexPairVector out;
  IndexVector mytpcs = tpcs();
  for ( Index itpc : mytpcs ) {
    for ( Index jtpc : intpcs ) {
      bool match = (itpc == jtpc) ||
                   ( !same && (itpc == badIndex() || jtpc == badIndex()) );
      if ( match ) out.push_back(IndexPair(itpc,jtpc));
    }
  }
  return out;
}

//**********************************************************************

Channel TpcSignalMap::channelMin() const {
  Channel chmin = badChannel();
  for ( auto tpcticksig : m_tpcticksig ) {
    TickChannelMap& ticksig = tpcticksig.second;
    if ( ticksig.size() != 0 ) {
      Channel newchmin = ticksig.begin()->first;
      if ( chmin == badChannel() || newchmin < chmin ) chmin = newchmin;
    }
  }
  return chmin;
}

//**********************************************************************

Channel TpcSignalMap::channelMax() const {
  Channel chmax = badChannel();
  for ( auto tpcticksig : m_tpcticksig ) {
    TickChannelMap& ticksig = tpcticksig.second;
    if ( ticksig.size() != 0 ) {
      Channel newchmax = ticksig.rbegin()->first;
      if ( chmax == badChannel() || newchmax > chmax ) chmax = newchmax;
    }
  }
  return chmax;
}

//**********************************************************************

Tick TpcSignalMap::tickMin() const {
  return m_tickRange.first();
}

//**********************************************************************

Tick TpcSignalMap::tickMax() const {
  return m_tickRange.last();
}

//**********************************************************************

TickRange TpcSignalMap::tickRange() const {
  return m_tickRange;
}

//**********************************************************************

unsigned int TpcSignalMap::channelCount() const {
  unsigned int count = 0;
  for ( auto tpcticksig : m_tpcticksig ) {
    TickChannelMap& ticksig = tpcticksig.second;
    count += ticksig.size();
  }
  return count;
}

//**********************************************************************

unsigned int TpcSignalMap::tickCount() const {
  return tickRange().size();
}

//**********************************************************************

unsigned int TpcSignalMap::binCount() const {
  unsigned int nbin = 0;
  for ( auto tpcticksig : m_tpcticksig ) {
    TickChannelMap& ticksig = tpcticksig.second;
    for ( const auto& echan : ticksig ) {
      nbin += echan.second.size();
    }
  }
  return nbin;
}

//**********************************************************************

unsigned int TpcSignalMap::size() const {
  return binCount();
}

//**********************************************************************

unsigned int TpcSignalMap::ropCount() const {
  const GeoHelper& geohelp = *geometryHelper();
  unsigned int nrop = 0;
  for ( Index irop=0; irop<geohelp.nrop(); ++irop ) {
    if ( ropNbin(irop) ) ++nrop;
  }
  return nrop;
}

//**********************************************************************

unsigned int TpcSignalMap::hitCount() const {
  unsigned int count = 0;
  for ( TpcHitChannelMap::value_type tpchitvecs : m_tpchitsig ) {
    for ( HitChannelMap::value_type chahitvec : tpchitvecs.second ) {
      const HitVector& hitvec = chahitvec.second;
      count += hitvec.size();
    }
  }
  return count;
}

//**********************************************************************

Signal TpcSignalMap::tickSignal() const {
  Signal sigtot = 0.0;
  for ( TpcTickChannelMap::value_type tpctickmap : m_tpcticksig ) {
    for ( TickChannelMap::value_type chticksigmap : tpctickmap.second ) {
      for ( SignalTickMap::value_type ticksig : chticksigmap.second ) {
        sigtot += ticksig.second;
      }
    }
  }
  return sigtot;
}

//**********************************************************************

Signal TpcSignalMap::hitSignal() const {
  Signal sig = 0.0;
  for ( TpcHitChannelMap::value_type tpchitvecs : m_tpchitsig ) {
    for ( HitChannelMap::value_type chahitvec : tpchitvecs.second ) {
      for ( const Hit& hit : chahitvec.second ) {
        sig += hit.signal;
      }
    }
  }
  return sig;
}

//**********************************************************************

int TpcSignalMap::fillChannelTickHist(TH2* ph) const {
  for ( auto tpcticksig : m_tpcticksig ) {
    TickChannelMap& ticksigmap = tpcticksig.second;
    for ( const auto& chanticksigs : ticksigmap ) {
      unsigned int chan = chanticksigs.first;
      for ( const auto& ticksig : chanticksigs.second ) {
        unsigned int tick = ticksig.first;
        double sig = ticksig.second;
        ph->Fill(tick, chan, sig);
      }
    }
  }
  return 0;
}

//**********************************************************************

int TpcSignalMap::fillRopChannelTickHist(TH2* ph, Index irop) const {
  const int dbg = 0;
  const string myname = "TpcSignalMap::fillRopChannelTickHist: ";
  if ( m_pgh == nullptr ) return -1;
  if ( dbg ) cout << myname << "Filling histogram " << ph->GetName() << endl;
  for ( auto tpcticksig : m_tpcticksig ) {
    TickChannelMap& ticksigmap = tpcticksig.second;
    for ( const auto& chanticksigs : ticksigmap ) {
      unsigned int chan = chanticksigs.first;
      if ( m_pgh->channelRop(chan)  == irop ) {
        Index ropchan = chan - m_pgh->ropFirstChannel(irop);
        for ( const auto& ticksig : chanticksigs.second ) {
          unsigned int tick = ticksig.first;
          double sig = ticksig.second;
          if ( dbg ) cout << myname << "ROPchan, tick, sig = " << ropchan
                          << ", " << tick << ", " << sig << endl;
          ph->Fill(tick, ropchan, sig);
        }
      }
    }
  }
  return 0;
}

//**********************************************************************

ostream& TpcSignalMap::print(ostream& out, int fulldetail, string hdrprefix, string prefix) const {
  const string myname = "TpcSignalMap::print: ";
  if ( int cstat = check() ) {
    out << myname << "ERROR: Object " << name() << " is invalid. Error " << cstat << "." << endl;
    return out;
  }
  int detail2 = fulldetail/10;
  int detail1 = fulldetail - detail2*10;
  // Line for each channel displaying the time range(s) and signal
  if ( detail1 == 0 ) {
    out << hdrprefix;
    if ( haveMcinfo() ) out << "MC track " << setw(3) << mcinfo()->trackID << " ";
    out << "TpcSignalMap " << setw(14) << name() << " ";
    if ( rop() == badIndex() ) {
      unsigned int nrop = ropCount();
      out << " uses " << setw(2) << nrop << " ROPs and";
    } else {
      out << " with ROP " << setw(2) << rop();
    }
    IndexVector mytpcs = tpcs();
    if ( usetpc() ) {
      unsigned int ntpc = mytpcs.size();
      out << setw(3) << ntpc;
      out << " TPC";
      if ( ntpc != 1 ) out << "s";
      else out << " ";
    } else {
      out << " no TPCs";
    }
    out << " and has " << setw(4) << channelCount() << " channels with "
        << setw(6) << hitCount() << " hits and " << setw(5) << tickCount() << " ticks"
        << " and " << setw(2) << segments().size() << " segments." << endl;
  } else if ( detail1 == 1 ) {
    out << hdrprefix << "TpcSignalMap map has " << channelCount() << " channels:" << endl;
    Tick badtick = badTick();
    for ( auto tpcticksig : m_tpcticksig ) {
      Index itpc = tpcticksig.first;
      TickChannelMap& ticksigmap = tpcticksig.second;
      if ( itpc == badIndex() ) {
        out << prefix << "No TPC." << endl;
      } else {
        out << prefix << "TPC " << itpc << endl;
      }
      for ( const auto& echan : ticksigmap ) {
        Channel chan = echan.first;
        out << prefix << "  " << setw(5) << chan << ": ";
        Signal channelSignal = 0;
        Tick tick1 = badtick;
        Tick tick2 = badtick;
        for ( const auto& etick : echan.second ) {
          Tick tick = etick.first;
          Signal signal = etick.second;
          channelSignal += signal;
          if ( tick1 == badtick ) {
            // Start the first range.
            tick1 = tick;
            tick2 = tick1;
          } else if ( tick == tick2+1 ) {
            // Extend the current range.
            tick2 = tick;
          } else {
            // Write the current range and start a new one.
            out << tick1;
            if ( tick2 > tick1 ) out << ":" << tick2;
            out << ", ";
            tick1 = tick;
            tick2 = tick1;
          }
        }  // end loop over ticks
        if ( tick1 != badtick ) {
          out << tick1;
          if ( tick2 > tick1 ) out << ":" << tick2;
          out << ", ADC = " << channelSignal << " MeV" << endl;
        } else {
          out << "No ticks filled." << endl;
        }
      }  // End loop over channels
    }  // End loop over TPCs
  // Line for each hit displaying the time range(s) and signal
  } else if ( detail1 == 2 ) {
    out << prefix << "TpcSignalMap map has " << hitCount() << " hits:" << endl;
    if ( m_tpchitsig.size() == 0 ) {
      out << prefix << "  Hit signal map is empty." << endl;
    }
    // Loop over channels.
    for ( TpcHitChannelMap::value_type tpchits : m_tpchitsig ) {
      for ( const auto& echan : tpchits.second ) {
        Channel chan = echan.first;
        // Loop over hits in the channel.
        for ( const auto& ehit : echan.second ) {
          out << prefix << "  " << setw(5) << chan << ": " << ehit.ticks.first();
          if ( ehit.ticks.last() > ehit.ticks.first() ) out << ":" << ehit.ticks.last();
          out << ", ADC = " << ehit.signal << " MeV" << endl;
        }  // End loop over hits.
      }  // End loop over channels.
    }  // End loop over TPCs.
  // Line for each channel displaying the time range(s) and signal
  } else if ( detail1 == 3 ) {
    out << prefix << "TpcSignalMap map has " << tickCount() << " ticks:" << endl;
    for ( auto tpcticksig : m_tpcticksig ) {
      Index itpc = tpcticksig.first;
      TickChannelMap& ticksigmap = tpcticksig.second;
      if ( itpc == badIndex() ) {
        out << prefix << "No TPC." << endl;
      } else {
        out << prefix << "TPC " << itpc << endl;
      }
      for ( const auto& echan : ticksigmap ) {
        Channel chan = echan.first;
        for ( const auto& etick : echan.second ) {
          Tick tick = etick.first;
          Signal signal = etick.second;
          out << prefix << "  " << setw(5) << chan << ": " << setw(6) << tick
              << " ADC = " << signal << " MeV" << endl;
        }
      }  // End loop over channels
    }  // End loop over TPCs
  } else {
    out << prefix << "  Invalid print option: " << fulldetail << endl;
    return out;
  }
  // Print ROP counts.
  if ( detail2 == 1 ) {
    if ( m_ropnbin.size() ) {
      out << prefix << "  ROP      Name    Nbin" << endl;
      for ( Index irop=0; irop<m_ropnbin.size(); ++irop ) {
        out << prefix << setw(5) << irop << setw(10) << m_pgh->ropName(irop) << setw(8) << m_ropnbin[irop] << endl;
      }
    }
  }
  return out;
}

//**********************************************************************

int TpcSignalMap::splitByRop(TpcSignalMapVector& tsms, bool splitByTpc) const {
  const string myname = "TpcSignalMap::splitByRop: ";
  if ( geometryHelper() == nullptr ) {
    cout << myname << "ERROR: Geometry helper not found. " << endl;
    return 1;
  }
  if ( haveRop() ) {
    cout << myname << "ERROR: ROP is already assigned." << endl;
    return 2;
  }
  if ( int cstat = check() ) {
    cout << myname << "ERROR: Object " << name() << " is in invalid state " << cstat << "." << endl;
    return 3;
  }
  const GeoHelper& geohelp = *geometryHelper();
  // Loop over ROPs in the detector.
  for ( Index irop=0; irop<geohelp.nrop(); ++irop ) {
    IndexVector roptpcs = geohelp.ropTpcs(irop);
    // Find the TPCs to cover.
    IndexVector tpcs;
    if ( splitByTpc ) {
      tpcs = sharedTpcs(roptpcs);
    } else {
      tpcs.push_back(GeoHelper::badIndex());
    }
    for ( Index itpc : tpcs ) {
      // Suffix to make names unique for ROPs with mutiple TPCs.
      string namesuf;
      if ( splitByTpc && itpc!=badIndex() && roptpcs.size()>1 ) {
        for ( unsigned int ktpc=0; ktpc<roptpcs.size(); ++ktpc ) {
          if ( roptpcs[ktpc] == itpc ) {
            ostringstream ssnamesuf;
            ssnamesuf << ktpc+1;
            //ssnamesuf << "tpc" << itpc;
            namesuf = ssnamesuf.str();
          }
        }
      }
      // Find the signals for the ROP and TPC.
      auto itts = m_tpcticksig.find(itpc);
      if ( itts == m_tpcticksig.end() ) continue;
      const TickChannelMap& ticksig = itts->second;
      auto ithv = m_tpchitsig.find(itpc);
      if ( ithv == m_tpchitsig.end() ) continue;
      const HitChannelMap& hitmap = ithv->second;
      Index ch1 = geohelp.ropFirstChannel(irop);
      Index ch2 = ch1 + geohelp.ropNChannel(irop);
      if ( dbg() ) cout << myname << "Creating map for ROP " << irop << endl;
      TpcSignalMapPtr psm(new TpcSignalMap("tmp", m_pgh, true));
      for ( Channel ich=ch1; ich<ch2; ++ich ) {
        TickChannelMap::const_iterator its = ticksig.find(ich);
        HitChannelMap::const_iterator ihs = hitmap.find(ich);
        if ( its != ticksig.end() ) {
          for ( auto ticksig : its->second ) {
            psm->addSignal(ich, ticksig.first, ticksig.second, itpc);
          }
        }
        if ( ihs != hitmap.end() ) psm->m_tpchitsig[itpc][ich] = ihs->second;
      }
      if ( dbg() ) cout << myname << "  ROP " << irop << " has "
                        << psm->tickCount() << " ticks and "
                        << psm->hitCount() << " hits." << endl;
      if ( psm->tickCount() || psm->m_tpchitsig.size() ) {
        if ( dbg() ) cout << myname << "  Creating new map." << endl;
        psm->m_name = m_name + geohelp.ropName(irop) + namesuf;
        psm->m_pmci = m_pmci;
        psm->m_rop = irop;
        // Find segments that use this ROP.
        for ( TpcSegmentPtr pseg : m_segments ) {
          for ( Index itpc : geohelp.ropTpcs(irop) ) {
            if ( int(itpc) == pseg->tpc ) {
              psm->m_segments.push_back(pseg);
              break;
            }
          }
        }
        if ( dbg() ) cout << myname << "Kept " << psm->segments().size()
                          << " of " << m_segments.size() << endl;
        tsms.push_back(psm);
      }
    }  // End loop over TPCs
  }  // End loop over ROPs
  return 0;
}

//**********************************************************************

const TpcSegmentVector& TpcSignalMap::segments() const {
  return m_segments;
}

//**********************************************************************

ostream& operator<<(const TpcSignalMap& rhs, ostream& out) {
  return rhs.print(out);
}

//**********************************************************************
