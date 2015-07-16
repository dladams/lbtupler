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
using tpc::badIndex2;

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
  // Do the matching.
  for ( const auto& p1 : m_c1 ) {
    if ( dbg > 1 ) cout << myname << "  Reference candidate " << i1 << endl;
    double dmin = maxDistance();
    Status stat = UNMATCHED;
    Index imin = m_c2.size();
    Index irop = ropMatch ? p1->rop() : 0;
    if ( dbg > 1 ) cout << myname << "  Match vector size for ROP " << irop << " is "
                        << indicesByRop[irop].size() << endl;
    for ( Index i2 : indicesByRop[irop] ) {
      P2 p2 = m_c2[i2];
      double dis = distance()(*p1, *p2);
      if ( dbg > 1 ) cout << myname << "    Matching candidate " << i2 << " has distance " << dis << endl;
      if ( dis < dmin ) {
        stat = MATCHED;
        imin = i2;
        dmin = dis;
      }
    }
    m_matchStatus.push_back(stat);
    m_matchIndex.push_back(imin);
    m_matchDistance.push_back(dmin);
    if ( dbg ) {
      if ( stat == MATCHED ) cout << myname << "  Match to " << imin << " with distance " << dmin << endl;
      else cout << myname << "  No match found" << endl;
    }
    ++i1;
  }
  // Remove duplicate matches.
  // For now remove the match with fewer bins.
  // Loop over all matches.
  for ( Index i1=0; i1<m_c1.size(); ++i1 ) {
    Status stati = matchStatus(i1);
    if ( stati != MATCHED ) continue;
    unsigned int i2 = matchIndex(i1);
    if ( i2 != badIndex() ) {
      // Loop over later matches.
      for ( Index j1=i1+1; j1<m_c1.size(); ++j1 ) {
        Status statj = matchStatus(i1);
        if ( statj != MATCHED ) continue;
        unsigned int j2 = matchIndex(j1);
        if ( j2 == i2 ) {
          unsigned int ni2 = c2.at(i2)->size();
          unsigned int nj2 = c2.at(j2)->size();
          bool dropj = nj2 < ni2;
          if ( dropj ) {
            m_matchStatus[j1] = DUPLICATE;
          } else {
            m_matchStatus[i1] = DUPLICATE;
            break;
          }
        }
      }
    }
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

unsigned int TpcSignalMatcher::size() const {
  return  m_matchIndex.size();
}

//**********************************************************************

TpcSignalMatcher::Status TpcSignalMatcher::matchStatus(Index iref) const {
  if ( iref >= m_matchStatus.size() ) return UNDEFINEDSTATUS;
  return m_matchStatus[iref];
}

//**********************************************************************

Index TpcSignalMatcher::matchIndex(Index iref) const {
  if ( iref >= m_matchIndex.size() ) return -1;
  return m_matchIndex[iref];
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
  int widx = 4;
  int wnam = 14;
  int wnch = 5;
  int wnti = 5;
  int wnbi = 7;
  int wdis = 9;
  if ( opt == 0 ) {
    ssout << setw(widx) << "idx"
          << " " << setw(wnam) << "Reference"
          << " " << setw(wnch) << "Nchan"
          << " " << setw(wnch) << "Ntick"
          << " " << setw(wnbi) << "Nbin"
          << " " << setw(wnam) << "Matched"
          << " " << setw(wdis) << "Distance";
    for ( Index i1=0; i1<m_c1.size(); ++i1 ) {
      Status stat = matchStatus(i1);
      ssout << "\n";
      ssout << setw(widx) << i1;
      ssout << " " << setw(wnam) << m_c1[i1]->name()
            << " " << setw(wnch) << m_c1[i1]->channelCount()
            << " " << setw(wnti) << m_c1[i1]->tickCount()
            << " " << setw(wnbi) << m_c1[i1]->binCount();
      int i2 = matchIndex(i1);
      if ( stat==MATCHED || stat==DUPLICATE ) {
        ssout << " " << setw(wnam) << m_c2[i2]->name()
              << " " << setw(wdis) << matchDistance(i1);
      }
      if ( stat == DUPLICATE ) ssout << " (duplicate)";
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
