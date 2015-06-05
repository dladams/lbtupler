// TpcSignalMatcher.cxx

#include "TpcSignalMatcher.h"
#include <sstream>
#include <ostream>
#include <iomanip>
#include "TpcSignalMapComparison.h"

using std::string;
using std::ostringstream;
using std::setw;
using std::cout;
using std::ostream;
using std::endl;
using std::vector;
using tpc::badIndex;

typedef TpcSignalMatcher::Index Index;
typedef TpcSignalMatcher::Distance Distance;
typedef TpcSignalMatcher::T1 T1;
typedef TpcSignalMatcher::T2 T2;

//**********************************************************************

namespace {

double channelFractionDistance(const T1& o1, const T2& o2) {
  TpcSignalMapComparison com(o1, o2);
  return 1.0 - com.channelFraction();
}

double binFractionDistance(const T1& o1, const T2& o2) {
  TpcSignalMapComparison com(o1, o2);
  return 1.0 - com.binFraction();
}

}

//**********************************************************************

TpcSignalMatcher::TpcSignalMatcher(const C1& c1, const C2& c2, bool ropMatch, int dbg)
: m_c1(c1), m_c2(c2) {
  const string myname = "TpcSignalMatcher::ctor: ";
  if ( dbg ) {
    cout << myname << "C1 size: " <<  m_c1.size() << endl;
    cout << myname << "C2 size: " <<  m_c2.size() << endl;
  }
  unsigned int i1 = 0;
  // Create vector that holds the indices of c2 indexed by ROP.
  // If ropMatch is not true, all are recorded with ROP=0.
  vector<IndexVector> indicesByRop(1);
  if ( ropMatch ) {
    Index maxrop = 0;
    for ( const auto& p1 : m_c1 ) {
      Index irop = p1->rop();
      if ( dbg > 1 ) cout << myname << "  irop1: " << irop << endl;
      if ( irop == badIndex() ) {
        cout << "WARNING: Found reference " << p1->name() << " without ROP." << endl;
        return;
      }
      if ( irop > maxrop ) maxrop = irop;
    }
    for ( const auto& p2 : m_c2 ) {
      Index irop = p2->rop();
      if ( dbg > 1 ) cout << myname << "  irop2: " << irop << endl;
      if ( irop == badIndex() ) {
        cout << "WARNING: Found match " << p2->name() << " without ROP." << endl;
        return;
      }
      if ( irop > maxrop ) maxrop = irop;
    }
    if ( dbg > 1 ) cout << myname << "Setting index vector size to " << maxrop << endl;
    indicesByRop.resize(maxrop+1);
    if ( dbg ) {
      cout << "Index counts before filling:" << endl;
      for ( Index irop=0; irop<indicesByRop.size(); ++irop ) {
        cout << setw(4) << irop << ": " << indicesByRop[irop].size() << endl;
      }
    }
    for ( Index i2=0; i2<m_c2.size(); ++i2 ) {
      P2 p2 = m_c2[i2];
      Index irop = p2->rop();
      indicesByRop[irop].push_back(i2);
    }
    if ( dbg ) {
      cout << "Index counts:" << endl;
      for ( Index irop=0; irop<indicesByRop.size(); ++irop ) {
        cout << setw(4) << irop << ": " << indicesByRop[irop].size() << endl;
      }
    }
  } else {
    for ( Index i2=0; i2<m_c2.size(); ++i2 ) {
      indicesByRop[0].push_back(i2);
    }
  }
  for ( const auto& p1 : m_c1 ) {
    if ( dbg > 1 ) cout << myname << "  Reference candidate " << i1 << endl;
    double dmin = maxDistance();
    Index imin = m_c2.size();
    Index irop = ropMatch ? p1->rop() : 0;
    if ( dbg > 1 ) cout << myname << "  Match vector size is " << indicesByRop[irop].size() << endl;
    for ( Index i2 : indicesByRop[irop] ) {
      P2 p2 = m_c2[i2];
      if ( dbg > 1 ) cout << myname << "    Matching candidate " << i2 << endl;
      double dis = distance()(*p1, *p2);
      if ( dbg > 1 ) cout << myname << "    Matching distance = " << dis << endl;
      if ( dis < dmin ) {
        imin = i2;
        dmin = dis;
      }
    }
    m_matchIndex.push_back(imin);
    m_matchDistance.push_back(dmin);
    if ( dbg ) cout << myname << "  Match to " << imin << " with distance " << dmin << endl;
    ++i1;
  }
}

//**********************************************************************

Distance TpcSignalMatcher::distance() const {
  if ( 1 ) return &binFractionDistance;
  else return &channelFractionDistance;
}

//**********************************************************************

double TpcSignalMatcher::maxDistance() const {
  return 1.0;
}

//**********************************************************************

int TpcSignalMatcher::matchIndex(Index iref) const {
  if ( iref >= m_matchIndex.size() ) return -1;
  Index i2 = m_matchIndex[iref];
  if ( i2 >= m_c2.size() ) return -1;
  return i2;
}

//**********************************************************************

float TpcSignalMatcher::matchDistance(Index iref) const {
  if ( iref >= m_matchIndex.size() ) return -1.0;
  Index i2 = m_matchIndex[iref];
  if ( i2 >= m_c2.size() ) return -1.0;
  return m_matchDistance[iref];
}

//**********************************************************************

string TpcSignalMatcher::show(int opt) const {
  ostringstream ssout;
  int wnam = 12;
  int wnch = 5;
  int wdis = 9;
  if ( opt == 0 ) {
    ssout << setw(wnam) << "Reference"
          << " " << setw(wnch) << "Nchan"
          << " " << setw(wnam) << "Matched"
          << " " << setw(wdis) << "Distance";
    for ( Index i1=0; i1<m_c1.size(); ++i1 ) {
      ssout << "\n";
      int i2 = matchIndex(i1);
      ssout << setw(wnam) << m_c1[i1]->name()
            << " " << setw(wnch) << m_c1[i1]->tickSignalMap().size();
      if ( i2 >= 0 ) {
        ssout << " " << setw(wnam) << m_c2[i2]->name()
              << " " << setw(wdis) << matchDistance(i1);
      }
    }
  } else {
    ssout << "Invalid option to TpcSignalMatcher::show";
  }
  return ssout.str();
}

//**********************************************************************

ostream& TpcSignalMatcher::print(ostream& out, int opt) const {
  out << show(opt) << endl;
  return out;
}

//**********************************************************************
