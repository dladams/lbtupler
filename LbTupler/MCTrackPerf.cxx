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

typedef MCTrackPerf::Tick    Tick;
typedef MCTrackPerf::Channel Channel;
typedef MCTrackPerf::Signal  Signal;

//**********************************************************************
// Local definitions.
//**********************************************************************

namespace {
int dbg() { return 0; }
}

//**********************************************************************
// Sub class.
//**********************************************************************

MCTrackPerf::Hit::Hit()
: tick1(badTick()), tick2(badTick()), signal(0.0) { }

//**********************************************************************

MCTrackPerf::Hit::Hit(unsigned int atick1, unsigned int atick2, double asignal)
: tick1(atick1), tick2(atick2), signal(asignal) { }

//**********************************************************************

Channel MCTrackPerf::badChannel() {
  return std::numeric_limits<Channel>::max();
}

//**********************************************************************

Tick MCTrackPerf::badTick() {
  return std::numeric_limits<Tick>::min();
}

//**********************************************************************

MCTrackPerf::MCTrackPerf()
: m_trackID(-1),
  m_pdg(0),
  m_rpdg(0),
  m_tickMin(badTick()),
  m_tickMax(badTick()) { }

//**********************************************************************

MCTrackPerf::MCTrackPerf(const MCParticle& par)
: m_trackID(par.TrackId()),
  m_pdg(par.PdgCode()),
  m_rpdg(reducedPDG(m_pdg)) { }

//**********************************************************************

int MCTrackPerf::addSignal(Channel chan, Tick tick, Signal signal) {
  const string myname = "MCTrackPerf::add: ";
  if ( dbg() ) std::cout << myname << "MC Track " << m_trackID
                         << ": channel " << chan << ", tick " << tick
                         << " has signal " << signal << endl;
  if ( m_ticksig.find(chan) == m_ticksig.end() ) m_ticksig[chan] = TickMap();
  TickMap& ticksig = m_ticksig[chan];
  if ( ticksig.find(tick) == ticksig.end() ) {
    ticksig[tick] = signal;
  } else {
    ticksig[tick] += signal;
  }
  return 0;
}

//**********************************************************************

int MCTrackPerf::addSimChannel(const SimChannel& simchan) {
  const string myname = "MCTrackPerf::addSimChannel: ";
  Channel chan = simchan.Channel();
  bool havechan = m_ticksig.find(chan) != m_ticksig.end();
  TickMap newticksig;
  TickMap& ticksig = havechan ? m_ticksig[chan] : newticksig;
  auto const& tickides = simchan.TDCIDEMap();
  if ( dbg() ) std::cout << myname << "Track " << m_trackID << ", channel " << chan << endl;
  for ( auto const& tickide : tickides ) {
    Tick tick = tickide.first;
    auto& ides = tickide.second;
    for ( auto& ide : ides ) {
      if ( abs(ide.trackID) == m_trackID ) {
        //Signal sig = ide.numElectrons;
        Signal sig = ide.energy;
        if ( dbg() ) std::cout << myname << "  Adding Tick=" << tick << ", Sig=" << sig << endl;
        if ( m_tickMin == badTick() || tick < m_tickMin ) m_tickMin = tick;
        if ( m_tickMax == badTick() || tick > m_tickMax ) m_tickMax = tick;
        if ( ticksig.find(tick) == ticksig.end() ) {
          ticksig[tick] = sig;
        } else {
          ticksig[tick] += sig;
        }
      } else {
        if ( dbg() > 1 ) std::cout << myname << "  Skipping track " << ide.trackID << endl;
      }
    }  // End loop over IDE's for this tick
  }  // End loop over ticks for this sim channel
  if ( !havechan && ticksig.size() ) {
    m_ticksig[chan] = ticksig;
  }
  return 0;
}

//**********************************************************************

int MCTrackPerf::buildHits() {
  const string myname = "MCTrackPerf::buildHits: ";
  // Loop over channels.
  for ( const auto ent : m_ticksig ) {
    Channel chan = ent.first;
    const auto& ticks = ent.second;
    Signal hitsig = 0.0;
    Tick tick1 = badTick();
    Tick tick2 = badTick();
    HitVector& hits = m_hitsig[chan];
    if ( hits.size() ) {
      cout << myname << "ERROR: Channel hits already defined." << endl;
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

const MCTrackPerf::TickChannelMap& MCTrackPerf::tickSignalMap() const {
  return m_ticksig;
}

//**********************************************************************

const MCTrackPerf::HitChannelMap& MCTrackPerf::hitSignalMap() const {
  return m_hitsig;
}

//**********************************************************************

Channel MCTrackPerf::channelMin() const {
  if ( m_ticksig.size() == 0 ) return badChannel();
  return m_ticksig.begin()->first;
}

//**********************************************************************

Channel MCTrackPerf::channelMax() const {
  if ( m_ticksig.size() == 0 ) return badChannel();
  return m_ticksig.rbegin()->first;
}

//**********************************************************************

Tick MCTrackPerf::tickMin() const {
  return m_tickMin;
}

//**********************************************************************

Tick MCTrackPerf::tickMax() const {
  return m_tickMax;
}

//**********************************************************************

unsigned int MCTrackPerf::size() const {
  unsigned int nbin = 0;
  for ( const auto& echan : tickSignalMap() ) {
    nbin += echan.second.size();
  }
  return nbin;
}

//**********************************************************************

unsigned int MCTrackPerf::channelCount() const {
  return tickSignalMap().size();
}

//**********************************************************************

Signal MCTrackPerf::tickSignal() const {
  Signal sig = 0.0;
  for ( const auto& echan : tickSignalMap() ) {
    for ( const auto& etick : echan.second ) {
      sig += etick.second;
    }
  }
  return sig;
}

//**********************************************************************

Signal MCTrackPerf::hitSignal() const {
  Signal sig = 0.0;
  for ( const auto& echan : hitSignalMap() ) {
    for ( const auto& ehit : echan.second ) {
      sig += ehit.signal;
    }
  }
  return sig;
}

//**********************************************************************

ostream& MCTrackPerf::print(ostream& out, int detail, string prefix) const {
  // Header line describing the MC particle.
  if ( detail > 0 ) {
    out << prefix << "MCParticle " << trackID()
        << ", PDG=" << pdg() << ", RPDG=" << rpdg()
        << ", #chan=" << channelCount()
        << ", #bin=" << size()
        << endl;
  }
  // Line for each channel displaying the time range(s) and signal
  if ( detail == 1 ) {
    Tick badtick = badTick();
    for ( const auto& echan : tickSignalMap() ) {
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
  // Line for each hit displaying the time range(s) and signal
  } else if ( detail == 2 ) {
    if ( hitSignalMap().size() == 0 ) {
      out << prefix << "  Hit signal map is empty." << endl;
    }
    // Loop over channels.
    for ( const auto& echan : hitSignalMap() ) {
      Channel chan = echan.first;
      // Loop over hits in the channel.
      for ( const auto& ehit : echan.second ) {
        out << prefix << "  " << setw(5) << chan << ": " << ehit.tick1;
        if ( ehit.tick2 > ehit.tick1 ) out << ":" << ehit.tick2;
        out << ", ADC = " << ehit.signal << " MeV" << endl;
      }  // End loop over hits.
    }  // End loop over channels.
  // Line for each channel displaying the time range(s) and signal
  } else if ( detail == 3 ) {
    for ( const auto& echan : tickSignalMap() ) {
      Channel chan = echan.first;
      for ( const auto& etick : echan.second ) {
        Tick tick = etick.first;
        Signal signal = etick.second;
        out << prefix << "  " << setw(5) << chan << ": " << setw(6) << tick
            << " ADC = " << signal << " MeV" << endl;
      }
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
