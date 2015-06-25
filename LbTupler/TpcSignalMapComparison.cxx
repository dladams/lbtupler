// TpcSignalMapComparison.cxx

#include "TpcSignalMapComparison.h"
#include <iomanip>
#include "Geometry/Geometry.h"
#include "TpcSignalMap.h"
#include "GeoHelper.h"

using std::string;

typedef TpcSignalMap::Channel Channel;
typedef TpcSignalMap::Index Index;
typedef TpcSignalMap::IndexPair IndexPair;

#include <iostream>
using std::cout;
using std::endl;
using std::setw;
using tpc::Tick;
using tpc::badTick;
typedef TpcSignalMap::TickChannelMap TickChannelMap;
typedef TpcSignalMap::SignalTickMap SignalTickMap;

namespace {
  int dbg() { return 0; }
}

//**********************************************************************

TpcSignalMapComparison::
TpcSignalMapComparison(const TpcSignalMap& tsm1, const TpcSignalMap& tsm2)
: m_tsm1(tsm1), m_tsm2(tsm2),
  m_rop(tpc::badIndex()),
  m_chbegin(0), m_chend(tpc::badChannel()),
  m_nchanref(tsm1.channelCount()), m_nchanmat(tsm2.channelCount()) {
  const string myname = "TpcSignalMapComparison::ctor: ";
  if ( tsm1.channelCount() == 0 ) return;
  if ( tsm2.channelCount() == 0 ) return;
  // Find the geometry.
  const GeoHelper* pgeohelp = tsm1.geometryHelper();
  if ( pgeohelp == nullptr ) pgeohelp = tsm2.geometryHelper();
  if ( pgeohelp == nullptr ) {
    cout << myname << "ERROR: Geometry not found." << endl;
    return;
  }
  // Find the ROP.
  Index irop1 = tsm1.rop();
  Index irop2 = tsm2.rop();
  Index badrop = tpc::badIndex();
  if ( irop1 != badrop ) {
    m_rop = irop1;
    if ( irop2 != badrop &&  irop2 != irop1 ) {
      cout << myname << "ERROR: Inconsistent ROPs: " << irop1 << ", " << irop2 << endl;
      return;
    }
  } else {
    if ( irop2 != badrop ) m_rop = irop2;
  }
  // Find the allowed channel range.
  m_chbegin = 0;
  m_chend = pgeohelp->geometry()->Nchannels();
  if ( m_rop != badrop ) {
    m_chbegin = pgeohelp->ropFirstChannel(m_rop);
    m_chend = m_chbegin + pgeohelp->ropNChannel(m_rop);
    m_nchanref = 0;
    m_nchanmat = 0;
    for ( Index itpc : m_tsm1.tpcs() ) {
      for ( const auto& chanticksigs1 : m_tsm1.tickSignalMap(itpc) ) {
        Channel ch1 = chanticksigs1.first;
        if ( ch1 < m_chbegin ) continue;
        if ( ch1 >= m_chend ) continue;
        ++m_nchanref;
      }
    }
    for ( Index itpc : m_tsm2.tpcs() ) {
      for ( const auto& chanticksigs2 : m_tsm2.tickSignalMap(itpc) ) {
        Channel ch2 = chanticksigs2.first;
        if ( ch2 < m_chbegin ) continue;
        if ( ch2 >= m_chend ) continue;
        ++m_nchanmat;
      }
    }
  }
}

//**********************************************************************

double TpcSignalMapComparison::channelFraction() const {
  const string myname = "TpcSignalMapComparison::channelFraction: ";
  double chanfrac = 0.0;
  // Find the fraction of overlapping channels (in the allowed range).
  unsigned int nchanNum = 0;
  unsigned int nchanDen = 0;
  if ( dbg() ) cout << myname << "Checking " << referenceChannelCount() << " and "
                    << matchChannelCount() << " channels"
                    << " over range [" << m_chbegin << "," << m_chend << "):" << endl;
  if ( dbg() > 1 ) {
    cout << myname << "Reference TPCs:";
    for ( Index itpc : m_tsm1.tpcs() ) cout << " " << itpc;
    cout << endl;
    cout << myname << "    Match TPCs:";
    for ( Index itpc : m_tsm2.tpcs() ) cout << " " << itpc;
    cout << endl;
  }
  for ( IndexPair ijtpc : m_tsm1.sharedTpcPairs(m_tsm2.tpcs()) ) {
    Index itpc1 = ijtpc.first;
    Index itpc2 = ijtpc.second;
    if ( dbg() > 1 ) cout << myname << "  TPCs (" << itpc1 << ", " << itpc2 << ")" << endl;;
    for ( const auto& ent1 : m_tsm1.tickSignalMap(itpc1) ) {
      Channel ch1 = ent1.first;
      if ( dbg() > 2 ) cout << myname << "    Channel " << ch1 << endl;;
      if ( ch1 < m_chbegin ) continue;
      if ( ch1 >= m_chend ) continue;
      const TickChannelMap& tsmap2 = m_tsm2.tickSignalMap(itpc2);
      bool match = tsmap2.find(ch1) != tsmap2.end();
      if ( dbg() ) cout << (match ? "+" : ".");
      if ( match ) ++nchanNum;
      ++nchanDen;
    }
  }
  if ( nchanDen > 0 ) chanfrac = double(nchanNum)/double(nchanDen);
  if ( dbg() ) cout << myname << "--chanfrac=" << chanfrac << endl;
  return chanfrac;
}

//**********************************************************************

double TpcSignalMapComparison::binFraction() const {
  const string myname = "TpcSignalMapComparison::binFraction: ";
  double binfrac = 0.0;
  if ( referenceChannelCount() == 0 ) {
    if ( dbg() ) cout << myname << "No reference channels." << endl;
    return binfrac;
  }
  if ( matchChannelCount() == 0 ) {
    if ( dbg() ) cout << myname << "No match channels." << endl;
    return binfrac;
  }
  // Find the fraction of overlapping bins (in the allowed channel range).
  unsigned int nbinNum = 0;
  unsigned int nbinDen = 0;
  if ( dbg() ) cout << myname << "Checking " << referenceChannelCount() << " and "
                    << matchChannelCount() << " channels for ROP " << rop() << endl;
  // Loop over channels in the reference.
  for ( IndexPair ijtpc : m_tsm1.sharedTpcPairs(m_tsm2.tpcs()) ) {
    Index itpc1 = ijtpc.first;
    Index itpc2 = ijtpc.second;
    for ( const auto& chanticksigs1 : m_tsm1.tickSignalMap(itpc1) ) {
      Channel ch1 = chanticksigs1.first;
      const SignalTickMap& ticksigs1 = chanticksigs1.second;
      if ( ch1 < m_chbegin ) continue;
      if ( ch1 >= m_chend ) continue;
      if ( dbg() ) cout << myname << setw(6) << ch1;
      const auto& ichanticksigs2 = m_tsm2.tickSignalMap(itpc2).find(ch1);
      const Tick badtick = badTick();
      Tick ticklast = badtick;
      // Find the matching channel in the match.
      const SignalTickMap* pticksigs2 = nullptr;
      if ( ichanticksigs2 != m_tsm2.tickSignalMap(itpc2).end() ) {
        pticksigs2 = &ichanticksigs2->second;
      }
      // Loop over ticks in the reference channel.
      for ( auto ticksig1 : ticksigs1 ) {
        Tick tick = ticksig1.first;
        bool match = (pticksigs2 != nullptr) && (pticksigs2->find(tick) != pticksigs2->end());
        if ( match ) ++nbinNum;
        ++nbinDen;
        if ( dbg() ) {
          if ( ticklast == badtick || tick != ticklast+1 ) cout << " ";
          cout << (match ? "+" : ".");
        }
        ticklast = tick;
      }
      if ( dbg() ) cout << endl;
    }  // End loop over channels.
  }  // End loop over TPCs.
  if ( nbinDen > 0 ) binfrac = double(nbinNum)/double(nbinDen);
  if ( dbg() ) cout << myname << " binfrac=" << binfrac << endl;
  return binfrac;
}

//**********************************************************************
