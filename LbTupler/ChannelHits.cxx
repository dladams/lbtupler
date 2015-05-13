// ChannelHits.cxx

#include "ChannelHits.h"
#include <iomanip>
#include "TH2.h"
#include "SimulationBase/MCParticle.h"
#include "Simulation/SimChannel.h"
#include "RecoBase/Hit.h"
#include "GeoHelper.h"

using std::string;
using std::cout;
using std::endl;
using std::setw;
using std::ostream;
using simb::MCParticle;
using sim::SimChannel;

typedef ChannelHits::Tick    Tick;
typedef ChannelHits::Channel Channel;
typedef ChannelHits::Signal  Signal;
typedef ChannelHits::Index   Index;

//**********************************************************************
// Local definitions.
//**********************************************************************

namespace {
int dbg() { return 0; }
}

//**********************************************************************
// Sub class.
//**********************************************************************

ChannelHits::Hit::Hit()
: tick1(badTick()), tick2(badTick()), signal(0.0) { }

//**********************************************************************

ChannelHits::Hit::Hit(unsigned int atick1, unsigned int atick2, double asignal)
: tick1(atick1), tick2(atick2), signal(asignal) { }

//**********************************************************************
// Main class.
//**********************************************************************

Channel ChannelHits::badChannel() {
  return std::numeric_limits<Channel>::max();
}

//**********************************************************************

Tick ChannelHits::badTick() {
  return std::numeric_limits<Tick>::min();
}

//**********************************************************************

ChannelHits::ChannelHits()
: m_pgh(nullptr),
  m_tickMin(badTick()),
  m_tickMax(badTick()) { }

//**********************************************************************

ChannelHits::ChannelHits(const GeoHelper* pgh)
: m_pgh(pgh),
  m_tickMin(badTick()),
  m_tickMax(badTick()) {
  if ( m_pgh != nullptr ) {
    m_ropnbin.resize(m_pgh->nrop(), 0);
  }
}

//**********************************************************************

int ChannelHits::addSignal(Channel chan, Tick tick, Signal signal) {
  const string myname = "ChannelHits::add: ";
  if ( dbg() ) std::cout << myname
                         << "Channel " << chan << ", tick " << tick
                         << " has signal " << signal << endl;
  if ( m_ticksig.find(chan) == m_ticksig.end() ) m_ticksig[chan] = SignalTickMap();
  SignalTickMap& ticksig = m_ticksig[chan];
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
  return 0;
}

//**********************************************************************

int ChannelHits::addHit(const recob::Hit& rhit, int dbg) {
  const string myname = "ChannelHits::addHit: ";
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
  HitVector& hits = m_hitsig[chan];
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

int ChannelHits::addSimChannel(const SimChannel& simchan, unsigned int tid) {
  const string myname = "ChannelHits::addSimChannel: ";
  Channel chan = simchan.Channel();
  auto const& tickides = simchan.TDCIDEMap();
  if ( dbg() ) std::cout << myname << "Track " << tid << ", channel " << chan << endl;
  for ( auto const& tickide : tickides ) {
    Tick tick = tickide.first;
    auto& ides = tickide.second;
    for ( auto& ide : ides ) {
      if ( abs(ide.trackID) == tid ) {
        //Signal sig = ide.numElectrons;
        Signal sig = ide.energy;
        if ( dbg() ) std::cout << myname << "  Adding Tick=" << tick << ", Sig=" << sig << endl;
        if ( m_tickMin == badTick() || tick < m_tickMin ) m_tickMin = tick;
        if ( m_tickMax == badTick() || tick > m_tickMax ) m_tickMax = tick;
        addSignal(chan, tick, sig);
      } else {
        if ( dbg() > 1 ) std::cout << myname << "  Skipping track " << ide.trackID << endl;
      }
    }  // End loop over IDE's for this tick
  }  // End loop over ticks for this sim channel
  return 0;
}

//**********************************************************************

int ChannelHits::buildHits() {
  const string myname = "ChannelHits::buildHits: ";
  // Loop over channels.
  for ( const auto ent : m_ticksig ) {
    Channel chan = ent.first;
    //const auto& ticks = ent.second;
    const SignalTickMap& ticks = ent.second;
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

const ChannelHits::TickChannelMap& ChannelHits::tickSignalMap() const {
  return m_ticksig;
}

//**********************************************************************

const ChannelHits::HitChannelMap& ChannelHits::hitSignalMap() const {
  return m_hitsig;
}

//**********************************************************************

Index ChannelHits::ropNbin(Index irop) const {
  if ( irop >=  m_ropnbin.size() ) return 0;
  return m_ropnbin[irop];
}

//**********************************************************************

Channel ChannelHits::channelMin() const {
  if ( m_ticksig.size() == 0 ) return badChannel();
  return m_ticksig.begin()->first;
}

//**********************************************************************

Channel ChannelHits::channelMax() const {
  if ( m_ticksig.size() == 0 ) return badChannel();
  return m_ticksig.rbegin()->first;
}

//**********************************************************************

Tick ChannelHits::tickMin() const {
  return m_tickMin;
}

//**********************************************************************

Tick ChannelHits::tickMax() const {
  return m_tickMax;
}

//**********************************************************************

unsigned int ChannelHits::size() const {
  return tickCount();
}

//**********************************************************************

unsigned int ChannelHits::channelCount() const {
  return tickSignalMap().size();
}

//**********************************************************************

unsigned int ChannelHits::tickCount() const {
  unsigned int nbin = 0;
  for ( const auto& echan : tickSignalMap() ) {
    nbin += echan.second.size();
  }
  return nbin;
}

//**********************************************************************

unsigned int ChannelHits::hitCount() const {
  unsigned int nbin = 0;
  for ( const auto& echan : hitSignalMap() ) {
    nbin += echan.second.size();
  }
  return nbin;
}

//**********************************************************************

Signal ChannelHits::tickSignal() const {
  Signal sig = 0.0;
  for ( const auto& echan : tickSignalMap() ) {
    for ( const auto& etick : echan.second ) {
      sig += etick.second;
    }
  }
  return sig;
}

//**********************************************************************

Signal ChannelHits::hitSignal() const {
  Signal sig = 0.0;
  for ( const auto& echan : hitSignalMap() ) {
    for ( const auto& ehit : echan.second ) {
      sig += ehit.signal;
    }
  }
  return sig;
}

//**********************************************************************

int ChannelHits::fillChannelTickHist(TH2* ph) const {
  for ( const auto& chanticksigs : tickSignalMap() ) {
    unsigned int chan = chanticksigs.first;
    for ( const auto& ticksig : chanticksigs.second ) {
      unsigned int tick = ticksig.first;
      double sig = ticksig.second;
      ph->Fill(tick, chan, sig);
    }
  }
  return 0;
}

//**********************************************************************

int ChannelHits::fillRopChannelTickHist(TH2* ph, Index irop) const {
  const int dbg = 1;
  const string myname = "ChannelHits::fillRopChannelTickHist: ";
  if ( m_pgh == nullptr ) return -1;
  if ( dbg ) cout << myname << "Filling histogram " << ph->GetName() << endl;
  for ( const auto& chanticksigs : tickSignalMap() ) {
    unsigned int chan = chanticksigs.first;
    if ( m_pgh->channelRop(chan)  == irop ) {
      Index ropchan = chan - m_pgh->ropFirstChannel(irop);
      for ( const auto& ticksig : chanticksigs.second ) {
        unsigned int tick = ticksig.first;
        double sig = ticksig.second;
        cout << myname << "ROPchan, tick, sig = " << ropchan << ", " << tick << ", " << sig << endl;
        ph->Fill(tick, ropchan, sig);
      }
    }
  }
  return 0;
}

//**********************************************************************

ostream& ChannelHits::print(ostream& out, int fulldetail, string hdrprefix, string prefix) const {
  int detail2 = fulldetail/10;
  int detail1 = fulldetail - detail2*10;
  cout << "DETAIL: " << fulldetail << ", " << detail1 << ", " << detail2 << endl;
  // Line for each channel displaying the time range(s) and signal
  if ( detail1 == 0 ) {
    out << hdrprefix << "ChannelHits map has " << channelCount() << " channels with "
        << hitCount() << " hits and " << tickCount() << " ticks." << endl;
  } else if ( detail1 == 1 ) {
    out << hdrprefix << "ChannelHits map has " << channelCount() << " channels:" << endl;
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
  } else if ( detail1 == 2 ) {
    out << prefix << "ChannelHits map has " << hitCount() << " channels:" << endl;
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
  } else if ( detail1 == 3 ) {
    out << prefix << "ChannelHits map has " << tickCount() << " ticks:" << endl;
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

ostream& operator<<(const ChannelHits& rhs, ostream& out) {
  return rhs.print(out);
}

//**********************************************************************
