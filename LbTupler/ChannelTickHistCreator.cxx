// ChannelTickHistCreator.cxx

#include "ChannelTickHistCreator.h"
#include <iostream>
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "TH2.h"

using std::string;
using std::cout;
using std::endl;
using art::TFileDirectory;

//**********************************************************************

ChannelTickHistCreator::
ChannelTickHistCreator(TFileDirectory& tfs, string sevt, int tick1, int tick2, 
                       string zlab, double zmin, double zmax, int ncontour)
: m_tfs(tfs),
  m_sevt(sevt),
  m_tickRange(tick1, tick2),
  m_zlab(zlab),
  m_zmin(zmin), m_zmax(zmax),
  m_ncontour(ncontour) { }

//**********************************************************************

TH2* ChannelTickHistCreator::
create(string slab, unsigned int chan1, unsigned int chan2, string stitle,
       string sevtNameSuffix, string sevtTitleSuffix, TickRange atickRange) const {
  const string myname = "ChannelTickHistCreator::create: ";
  const int dbg = 0;    // 0 for normal running
  TH2* ph = nullptr;
  if ( chan2 <= chan1 ) return nullptr;
  int nchan = chan2 - chan1;
  int mtick1 = m_tickRange.first();
  int mtick2 = m_tickRange.last();
  int tick1 = mtick1;
  int tick2 = mtick2;
  int atick1 = atickRange.first();
  int atick2 = atickRange.last();
  // If caller has supplied a tick range, use it to find the x-axis range.
  if ( atick2 >= atick1 ) {
    tick1 = atick1 - 5;
    tick2 = atick2 + 6;
    int fac = 50;
    int dtick1 = tick1/fac;
    int dtick2 = tick2/fac;
    if ( tick1 < 0 ) dtick1 -= 1;
    if ( tick2 > 0 ) dtick2 += 1;
    tick1 = fac*dtick1;
    tick2 = fac*dtick2;
    if ( tick1 < mtick1 ) tick1 = mtick1;
    if ( tick2 > mtick2 ) tick2 = mtick2;
    if ( dbg ) cout << myname << "X-axis: " << atick1 << "-" << atick2 << " ==> "
                    << tick1 << "-" << tick2 << endl;
    // Check if ther is no overlap between object and requested ranges.
    if ( tick2 <= tick1 ) return nullptr;
  }
  int ntick = tick2 - tick1;
  string hname = "h" + m_sevt + sevtNameSuffix + "_" + slab;
  string title = stitle + " event " + m_sevt;
  if ( sevtTitleSuffix.size() ) title += " " + sevtTitleSuffix;
  title += ";TDC tick;Channel;" + m_zlab;
  if ( dbg > 0 ) cout << myname << "Creating hit histo " << hname << " with " << ntick
                      << " TDC bins and " << nchan << " channel bins" << endl;
  if ( dbg > 1 ) cout << myname << "  Using TFileService." << endl;
  ph = m_tfs.make<TH2F>(hname.c_str(), title.c_str(),
                          ntick, tick1, tick2, nchan, chan1, chan2);
  if ( ph == nullptr ) {
    cout << myname << "Unable to create histogram " << hname << endl;
  } else {
    if ( dbg > 1 ) cout << myname << "  Created histogram " << ph->GetName() << endl;
    ph->GetZaxis()->SetRangeUser(m_zmin, m_zmax);
    ph->SetContour(m_ncontour);
    ph->SetStats(0);
  }
  return ph;
}

//**********************************************************************
