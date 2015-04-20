// MCTrackPerf.cxx

#include "MCTrackPerf.h"
#include <iomanip>
#include "SimulationBase/MCParticle.h"
#include "Simulation/SimChannel.h"
#include "reducedPDG.h"

using std::string;
using std::cout;
using std::endl;
using std::setw;
using std::ostream;
using simb::MCParticle;
using sim::SimChannel;

//**********************************************************************
// Local definitions.
//**********************************************************************

namespace {
int dbg() { return 1; }
}

//**********************************************************************
// Sub class.
//**********************************************************************

MCTrackPerf::Hit::Hit()
: tick1(badIndex()), tick2(badIndex()), charge(0.0) { }

//**********************************************************************

MCTrackPerf::Hit::Hit(unsigned int atick1, unsigned int atick2, double acharge)
: tick1(atick1), tick2(atick2), charge(acharge) { }

//**********************************************************************

unsigned int MCTrackPerf::badIndex() {
  return -1;
}

//**********************************************************************

MCTrackPerf::MCTrackPerf()
: m_trackID(-1),
  m_pdg(0),
  m_rpdg(0),
  m_tickMin(badIndex()),
  m_tickMax(badIndex()) { }

//**********************************************************************

MCTrackPerf::MCTrackPerf(const MCParticle& par)
: m_trackID(par.TrackId()),
  m_pdg(par.PdgCode()),
  m_rpdg(reducedPDG(m_pdg)) { }

//**********************************************************************

int MCTrackPerf::add(unsigned int chan, unsigned int tick, double signal) {
  const string myname = "MCTrackPerf::add: ";
  if ( dbg() ) std::cout << myname << "MC Track " << m_trackID
                         << ": channel " << chan << ", tick " << tick
                         << " has signal " << signal << endl;
  if ( m_tickchg.find(chan) == m_tickchg.end() ) m_tickchg[chan] = TickMap();
  TickMap& tickchg = m_tickchg[chan];
  if ( tickchg.find(tick) == tickchg.end() ) {
    tickchg[tick] = signal;
  } else {
    tickchg[tick] += signal;
  }
  return 0;
}

//**********************************************************************

int MCTrackPerf::addSimChannel(const SimChannel& simchan) {
  const string myname = "MCTrackPerf::addSimChannel: ";
  unsigned int chan = simchan.Channel();
  bool havechan = m_tickchg.find(chan) != m_tickchg.end();
  TickMap newtickchg;
  TickMap& tickchg = havechan ? m_tickchg[chan] : newtickchg;
  auto const& tickides = simchan.TDCIDEMap();
  if ( dbg() ) std::cout << myname << "Track " << m_trackID << ", channel " << chan << endl;
  for ( auto const& tickide : tickides ) {
    unsigned int tick = tickide.first;
    auto& ides = tickide.second;
    for ( auto& ide : ides ) {
      if ( abs(ide.trackID) == m_trackID ) {
        double chg = ide.numElectrons;
        if ( dbg() ) std::cout << myname << "  Adding Tick=" << tick << ", Chg=" << chg << endl;
        if ( m_tickMin == badIndex() || tick < m_tickMin ) m_tickMin = tick;
        if ( m_tickMax == badIndex() || tick > m_tickMax ) m_tickMax = tick;
        if ( tickchg.find(tick) == tickchg.end() ) {
          tickchg[tick] = chg;
        } else {
          tickchg[tick] += chg;
        }
      } else {
        if ( dbg() > 1 ) std::cout << myname << "  Skipping track " << ide.trackID << endl;
      }
    }  // End loop over IDE's for this tick
  }  // End loop over ticks for this sim channel
  if ( !havechan && tickchg.size() ) {
    m_tickchg[chan] = tickchg;
  }
  return 0;
}

//**********************************************************************

int MCTrackPerf::buildHits() {
  const string myname = "MCTrackPerf::buildHits: ";
  // Loop over channels.
  for ( const auto ent : m_tickchg ) {
    unsigned int chan = ent.first;
    const auto& ticks = ent.second;
    double hitchg = 0.0;
    unsigned int tick1 = badIndex();
    unsigned int tick2 = badIndex();
    HitVector& hits = m_hitchg[chan];
    if ( hits.size() ) {
      cout << myname << "ERROR: Channel hits already defined." << endl;
      return 1;
    }
    for ( auto tickchg : ticks ) {
      unsigned int tick = tickchg.first;
      double chg = tickchg.first;
      if ( tick1 != badIndex() && tick != tick2+1 ) {
        hits.push_back(Hit(tick1, tick2, chg));
        tick1 = badIndex();
        tick2 = badIndex();
        hitchg = 0.0;
      }
      if ( tick1 == badIndex() ) {
        tick1 = tick;
        tick2 = tick1;
      } else {
        tick2 = tick;
      }
      hitchg += chg;
    }
    hits.push_back(Hit(tick1, tick2, hitchg));
  }  // End loop over channels.
  return 0;
}

//**********************************************************************

unsigned int MCTrackPerf::trackID() const {
  return m_trackID;
}

//**********************************************************************

int MCTrackPerf::pdg() const {
  return m_pdg;
}

//**********************************************************************

int MCTrackPerf::rpdg() const {
  return m_rpdg;
}

//**********************************************************************

const MCTrackPerf::TickChannelMap& MCTrackPerf::tickChargeMap() const {
  return m_tickchg;
}

//**********************************************************************

const MCTrackPerf::HitChannelMap& MCTrackPerf::hitChargeMap() const {
  return m_hitchg;
}

//**********************************************************************

unsigned int MCTrackPerf::channelMin() const {
  if ( m_tickchg.size() == 0 ) return badIndex();
  return m_tickchg.begin()->first;
}

//**********************************************************************

unsigned int MCTrackPerf::channelMax() const {
  if ( m_tickchg.size() == 0 ) return badIndex();
  return m_tickchg.rbegin()->first;
}

//**********************************************************************

unsigned int MCTrackPerf::tickMin() const {
  return m_tickMin;
}

//**********************************************************************

unsigned int MCTrackPerf::tickMax() const {
  return m_tickMax;
}

//**********************************************************************

unsigned int MCTrackPerf::size() const {
  unsigned int nbin = 0;
  for ( const auto& echan : tickChargeMap() ) {
    nbin += echan.second.size();
  }
  return nbin;
}

//**********************************************************************

ostream& MCTrackPerf::print(ostream& out, int detail, string prefix) const {
  // Header line describing the MC particle.
  if ( detail > 0 ) {
    out << prefix << "MCParticle " << trackID()
        << ", PDG=" << pdg() << ", RPDG=" << rpdg()
        << ", #bin=" << size()
        << endl;
  }
  // Line for each channel displaying the time range(s) and charge
  if ( detail == 1 ) {
    const unsigned int badtick = badIndex();
    for ( const auto& echan : tickChargeMap() ) {
      unsigned int chan = echan.first;
      out << prefix << "  " << setw(5) << chan << ": ";
      double channelCharge = 0;
      unsigned int tick1 = badtick;
      unsigned int tick2 = badtick;
      for ( const auto& etick : echan.second ) {
        unsigned int tick = etick.first;
        double charge = etick.second;
        channelCharge += charge;
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
          if ( tick2 > tick1 ) out << "-" << tick2;
          out << ", ";
          tick1 = tick;
          tick2 = tick1;
        }
      }  // end loop over ticks
      if ( tick1 != badtick ) {
        out << tick1;
        if ( tick2 > tick1 ) out << "-" << tick2;
        out << ", ADC = " << channelCharge << " MeV" << endl;
      } else {
        out << "No ticks filled." << endl;
      }
    }  // End loop over channels
  // Line for each hit displaying the time range(s) and charge
  } else if ( detail == 2 ) {
    if ( hitChargeMap().size() == 0 ) {
      out << prefix << "  Hit charge map is empty." << endl;
    }
    // Loop over channels.
    for ( const auto& echan : hitChargeMap() ) {
      unsigned int chan = echan.first;
      // Loop over hits in the channel.
      for ( const auto& ehit : echan.second ) {
        out << prefix << "  " << setw(5) << chan << ": " << ehit.tick1;
        if ( ehit.tick2 > ehit.tick1 ) out << "-" << ehit.tick2;
        out << ", ADC = " << ehit.charge << " MeV" << endl;
      }  // End loop over hits.
    }  // End loop over channels.
  } else {
    out << prefix << "  Invalid print option: " << detail << endl;
  }
  return out;
}

//**********************************************************************

ostream& operator<<(const MCTrackPerf& rhs, ostream& out) {
  return rhs.print(out);
}

//**********************************************************************
