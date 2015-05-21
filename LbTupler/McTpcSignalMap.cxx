// McTpcSignalMap.cxx

#include "McTpcSignalMap.h"
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

typedef McTpcSignalMap::Tick    Tick;
typedef McTpcSignalMap::Channel Channel;
typedef McTpcSignalMap::Signal  Signal;
typedef McTpcSignalMap::Index   Index;

//**********************************************************************
// Local definitions.
//**********************************************************************

namespace {
}

//**********************************************************************

McTpcSignalMap::McTpcSignalMap()
: m_trackID(-1),
  m_pdg(0),
  m_rpdg(0) { }

//**********************************************************************

McTpcSignalMap::McTpcSignalMap(const MCParticle& par, const GeoHelper* pgh)
: TpcSignalMap(pgh),
  m_trackID(par.TrackId()),
  m_pdg(par.PdgCode()),
  m_rpdg(reducedPDG(m_pdg)) { }

//**********************************************************************

int McTpcSignalMap::addSimChannel(const sim::SimChannel& sch) {
  return TpcSignalMap::addSimChannel(sch, trackID());
}

//**********************************************************************

unsigned int McTpcSignalMap::trackID() const {
  return m_trackID;
}

//**********************************************************************

int McTpcSignalMap::pdg() const {
  return m_pdg;
}

//**********************************************************************

int McTpcSignalMap::rpdg() const {
  return m_rpdg;
}

//**********************************************************************

ostream& McTpcSignalMap::print(ostream& out, int detail, string prefix) const {
  ostringstream sout;
  sout << prefix << "MCParticle " << setw(3) << trackID()
        << ", PDG=" << setw(6) << pdg() << ", RPDG=" << setw(4) << rpdg() << ", ";
  return TpcSignalMap::print(out, detail, sout.str(), prefix);
}

//**********************************************************************

ostream& operator<<(const McTpcSignalMap& rhs, ostream& out) {
  return rhs.print(out);
}

//**********************************************************************
