// MCTrackPerf.cxx

#include "MCTrackPerf.h"
#include <iomanip>
#include <sstream>
#include "SimulationBase/MCParticle.h"
#include "Simulation/SimChannel.h"
#include "reducedPDG.h"

using std::string;
using std::cout;
using std::endl;
using std::setw;
using std::ostream;
using std::ostringstream;
using simb::MCParticle;
using sim::SimChannel;

typedef MCTrackPerf::Tick    Tick;
typedef MCTrackPerf::Channel Channel;
typedef MCTrackPerf::Signal  Signal;
typedef MCTrackPerf::Index   Index;

//**********************************************************************
// Local definitions.
//**********************************************************************

namespace {
}

//**********************************************************************

MCTrackPerf::MCTrackPerf()
: m_trackID(-1),
  m_pdg(0),
  m_rpdg(0) { }

//**********************************************************************

MCTrackPerf::MCTrackPerf(const MCParticle& par, const GeoHelper* pgh)
: m_trackID(par.TrackId()),
  m_pdg(par.PdgCode()),
  m_rpdg(reducedPDG(m_pdg)),
  m_hits(pgh) { }

//**********************************************************************

int MCTrackPerf::addSignal(Channel chan, Tick tick, Signal signal) {
  return m_hits.addSignal(chan, tick, signal);
}

//**********************************************************************

int MCTrackPerf::addSimChannel(const SimChannel& simchan) {
  return m_hits.addSimChannel(simchan, trackID());
}

//**********************************************************************

int MCTrackPerf::buildHits() {
  return m_hits.buildHits();
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

const TpcSignalMap& MCTrackPerf::hits() const {
  return m_hits;
}

//**********************************************************************

int MCTrackPerf::fillRopChannelTickHist(TH2* ph, Index irop) const {
  return m_hits.fillRopChannelTickHist(ph, irop);
}

//**********************************************************************

ostream& MCTrackPerf::print(ostream& out, int detail, string prefix) const {
  ostringstream sout;
  sout << prefix << "MCParticle " << setw(3) << trackID()
        << ", PDG=" << setw(6) << pdg() << ", RPDG=" << setw(4) << rpdg() << ", ";
  return m_hits.print(out, detail, sout.str(), prefix);
}

//**********************************************************************

ostream& operator<<(const MCTrackPerf& rhs, ostream& out) {
  return rhs.print(out);
}

//**********************************************************************
