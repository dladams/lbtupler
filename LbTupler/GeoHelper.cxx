// GeoHelper.cxx

#include "GeoHelper.h"
#include "Geometry/Geometry.h"
#include "Utilities/DetectorProperties.h"
#include "Geometry/TPCGeo.h"
#include "Geometry/PlaneGeo.h"
#include <iostream>
#include <sstream>
#include <iomanip>
#include <set>

using std::ostream;
using std::cout;
using std::endl;
using std::ostringstream;
using std::setw;
using std::vector;
using std::set;
using std::string;
using geo::TPCGeo;
using geo::PlaneID;
using geo::PlaneGeo;
using geo::kU;
using geo::kV;
using geo::kZ;

typedef GeoHelper::Status Status;
typedef GeoHelper::Index Index;

//**********************************************************************

Index GeoHelper::badIndex() {
  return std::numeric_limits<Index>::max();
}

//**********************************************************************

GeoHelper::GeoHelper(const geo::Geometry* pgeo, const util::DetectorProperties* pdetp, Status dbg)
: m_pgeo(pgeo), m_pdetp(pdetp), m_dbg(dbg), m_ntpc(0), m_ntpp(0), m_napa(0), m_nrop(0) {
  // Find the total number of TPC.
  for ( Index icry=0; icry<ncryostat(); ++icry ) m_ntpc += m_pgeo->NTPC(icry);
  // Fill the TPC info w/o ROP or APA.
  for ( Index icry=0; icry<ncryostat(); ++icry ) {
    for ( Index itpc=0; itpc<m_pgeo->NTPC(icry); ++itpc ) {
      m_ntpp += m_pgeo->Nplanes(itpc, icry);
      m_tpccry.push_back(icry);
      ostringstream ssname;
      ssname << "TPC";
      if ( m_ntpc > 9 && itpc < 10 ) ssname << "0";
      ssname << itpc;
      m_tpcname.push_back(ssname.str());
    }
  }
  fillStandardApaMapping();
}

//**********************************************************************

unsigned int GeoHelper::ncryostat() const {
  return m_pgeo->Ncryostats();
}

//**********************************************************************

int GeoHelper::tpcCorners(unsigned int icry, unsigned int itpc, double* pos1, double* pos2) const {
  const TPCGeo& tpcgeo = m_pgeo->TPC(itpc, icry);
  double origin[3] = {0.0, 0.0, 0.0};
  double tpcpos[3] = {0.0, 0.0, 0.0};
  tpcgeo.LocalToWorld(origin, tpcpos);
  double haflen[3] = {tpcgeo.ActiveHalfWidth(), tpcgeo.ActiveHalfHeight(), 0.5*tpcgeo.ActiveLength()};
  for ( unsigned int ixyz=0; ixyz<3; ++ixyz ) {
    pos1[ixyz] = tpcpos[ixyz] - haflen[ixyz];
    pos2[ixyz] = tpcpos[ixyz] + haflen[ixyz];
  }
  return 0;
}

//**********************************************************************

Index GeoHelper::rop(geo::PlaneID pid) const {
  const string myname = "GeoHelper::Rop: ";
  auto itrplane = m_tpprop.find(pid);
  if ( itrplane == m_tpprop.end() ) {
    cout << myname << "WARNING: Unknown plane." << endl;
    return badIndex();
  }
  return itrplane->second;
}

//**********************************************************************

ostream& GeoHelper::print(ostream& out, int iopt, std::string prefix) const {
  if ( m_pgeo == nullptr ) {
    out << prefix << "Geometry is not defined." << endl;
    return out;
  }
  unsigned int wlab = 30;
  unsigned int wcry =  4;
  unsigned int wtpc =  4;
  unsigned int wdim = 12;
  unsigned int ncry = ncryostat();
  out << prefix << setw(wlab) << "Name: " << m_pgeo->DetectorName() << endl;
  out << prefix << setw(wlab) << "Cryostat count: " << ncry << endl;
  for ( unsigned int icry=0; icry<ncry; ++icry ) {
    ostringstream sslab;
    sslab << "Cryostat " << icry << " TPC count: ";
    out << prefix << setw(wlab) << sslab.str() << m_pgeo->NTPC(icry) << endl;
  }
  out << prefix << setw(wlab) <<"TPC channel count: " << m_pgeo->Nchannels() << endl;
  out << prefix << setw(wlab) <<"Optical channel count: " << m_pgeo->NOpDet() << endl;
  out << prefix << setw(wlab) <<"Scintillator channel count: " << m_pgeo->NAuxDets() << endl;
  out << prefix << "TPC sizes [cm]:" << endl;
  out << prefix << setw(wcry) << "Cry" << setw(wtpc) << "TPC"
      << setw(wdim) << "Width" << setw(wdim) << "Height"<< setw(wdim) << "Length"
      << setw(wdim) << "   xc" << "    yc"<< setw(wdim) << "    zc"
      << endl;
  // Check plane counts.
  for ( unsigned int icry=0; icry<ncry; ++icry ) {
    for ( unsigned int itpc=0; itpc<m_pgeo->NTPC(icry); ++itpc ) {
      const TPCGeo& tpcgeo = m_pgeo->TPC(itpc, icry);
      unsigned int nplanes = tpcgeo.Nplanes();
      if ( nplanes != 3 ) {
        out << prefix << "WARNING: Cryostat " << icry << " TPC " << itpc
            << " has plane count " << nplanes << " where 3 is expected." << endl;
      }
    }  // end loop over TPCs
  }  // end loop over cryostats
  // Display size and position of each TPC.
  for ( unsigned int icry=0; icry<ncry; ++icry ) {
    for ( unsigned int itpc=0; itpc<m_pgeo->NTPC(icry); ++itpc ) {
      const TPCGeo& tpcgeo = m_pgeo->TPC(itpc, icry);
      double origin[3] = {0.0, 0.0, 0.0};
      double tpcpos[3] = {0.0, 0.0, 0.0};
      tpcgeo.LocalToWorld(origin, tpcpos);
      out << prefix << setw(wcry) << icry << setw(wtpc) << itpc
          << setw(wdim) << 2.0*tpcgeo.ActiveHalfWidth()
          << setw(wdim) << 2.0*tpcgeo.ActiveHalfHeight()
          << setw(wdim) <<     tpcgeo.ActiveLength()
          << setw(wdim) << tpcpos[0]
          << setw(wdim) << tpcpos[1]
          << setw(wdim) << tpcpos[2]
          << endl;
    }  // end loop over TPCs
  }  // end loop over cryostats
  out << prefix << "TPC corners [cm]:" << endl;
  out << prefix << setw(wcry) << "Cry" << setw(wtpc) << "TPC"
      << setw(wdim) << "x1" << setw(wdim) << "y1"<< setw(wdim) << "z1"
      << setw(wdim) << "x2" << setw(wdim) << "y2"<< setw(wdim) << "z2"
      << endl;
  // Display corners of each TPC.
  for ( unsigned int icry=0; icry<ncry; ++icry ) {
    for ( unsigned int itpc=0; itpc<m_pgeo->NTPC(icry); ++itpc ) {
      const TPCGeo& tpcgeo = m_pgeo->TPC(itpc, icry);
      unsigned int nplanes = tpcgeo.Nplanes();
      if ( nplanes != 3 ) {
        out << prefix << "WARNING: Cryostat " << icry << " TPC " << itpc
            << " has plane count " << nplanes << " where 3 is expected." << endl;
        continue;
      }
      double pos1[3];
      double pos2[3];
      tpcCorners(icry, itpc, pos1, pos2);
      out << prefix << setw(wcry) << icry << setw(wtpc) << itpc
          << setw(wdim) << pos1[0]
          << setw(wdim) << pos1[1]
          << setw(wdim) << pos1[2]
          << setw(wdim) << pos2[0]
          << setw(wdim) << pos2[1]
          << setw(wdim) << pos2[2]
          << endl;
    }  // end loop over TPCs
  }  // end loop over cryostats
  // Display ROPs.
  cout << prefix << "Detector ROP count: " << nrop() << endl;
  cout << prefix << setw(18) << "Name"
       << setw(5) << "View" << setw(8) << "First ch" << setw(6) << "Nchan" << endl;
  for ( Index irop=0; irop<nrop(); ++irop ) {
    cout << prefix << setw(4) << irop << ". " << setw(12) << ropName(irop)
         << setw(5) << ropView(irop) << setw(8) << ropFirstChannel(irop)
         << setw(6) << ropNChannel(irop) << endl;
  }
/*
      const double* pos0 = tpcgeo.PlaneLocation(0);
      const double* pos1 = tpcgeo.PlaneLocation(1);
      const double* pos2 = tpcgeo.PlaneLocation(2);
      if ( (pos1[1] != pos0[1]) || (pos2[1] != pos0[1]) || (pos1[2] != pos0[2]) ||(pos2[2] != pos0[2]) ) {
        out << prefix << "WARNING: Cryostat " << icry << " TPC " << itpc
            << " 2nd or 3rd plane has unexpected position:" << endl;
        out << prefix << setw(wdim) << pos0[0] << setw(wdim) << pos0[1] << setw(wdim) << pos0[2] << endl;
        out << prefix << setw(wdim) << pos1[0] << setw(wdim) << pos1[1] << setw(wdim) << pos1[2] << endl;
        out << prefix << setw(wdim) << pos2[0] << setw(wdim) << pos2[1] << setw(wdim) << pos2[2] << endl;
        continue;
      }
*/
  //out << prefix << setw(wlab) <<"TPC plane count: " << m_pgeo->Nplanes(itpc, icry) << endl;
  //out << prefix << setw(wlab) <<"TPC wire count: " << m_pgeo->Nwires(ipla, itpc, icry) << endl;
  //out << prefix << setw(wlab) <<"TPC view count: " << m_pgeo->Nviews() << endl;
  return out;
}
  
//**********************************************************************

Status GeoHelper::fillStandardApaMapping() {
  const string myname = "GeoHelper::fillStandardApaMapping: ";
  // Fill the TPC arrays.
  Index itpc = 0;  // Global index for the current TPC
  Index itpp = 0;  // Global index for the current TPC plane
  vector<int> firsttpcplanewire;
  firsttpcplanewire.push_back(0);
  Index nzplane = 0;
  // Loop over cryostats.
  for ( Index icry=0; icry<ncryostat(); ++icry ) {
    if ( m_dbg > 1 ) cout << myname << "Begin loop over cryostat " << icry << endl;
    // Loop over TPCs in the current cryostat.
    for ( Index icrytpc=0; icrytpc<m_pgeo->NTPC(icry); ++icrytpc ) {
      if ( m_dbg > 1 ) cout << myname << "Begin loop over TPC " << icrytpc << endl;
      const TPCGeo& tpcgeo = m_pgeo->TPC(icrytpc, icry);
      unsigned int iapa = itpc/2;   // Stanard mapping: one APA for each adjacent pair of TPCs
      int nplane = m_pgeo->Nplanes(itpc, icry);
      int itdcrop = 0;   // # readouts for this TDC
      m_tpcapa.push_back(iapa);
      // Loop over planes in the TPC.
      for ( int ipla=0; ipla<nplane; ++ipla ) {
        if ( m_dbg > 1 ) cout << myname << "Begin loop over plane " << ipla << endl;
        const PlaneGeo& plageo = tpcgeo.Plane(ipla);
        int nwire = m_pgeo->Nwires(ipla, icrytpc, icry);
        int firstwire = firsttpcplanewire[itpp];
        int lastwire = firstwire + nwire - 1;
        if ( itpc<ntpc()-1 || ipla<nplane-1 ) {
          firsttpcplanewire.push_back(lastwire + 1);
        }
        // Loop over wires and find the channels.
        set<int> chans;
        for ( int iwir=0; iwir<nwire; ++iwir ) {
          int icha = m_pgeo->PlaneWireToChannel(ipla, iwir, itpc, icry);
          chans.insert(icha);
        }
        int nchan = chans.size();
        int firstchan = -1;
        int lastchan = -1;
        if ( nchan ) {
          firstchan = *chans.cbegin();
          lastchan = *chans.crbegin();
        }
        if ( lastchan <= firstchan ) {
          cout << myname << "ERROR: Invalid channel range." << endl;
          abort();
        }
        if ( m_dbg > 0 ) {
          ostringstream ssplane;
          cout << myname << "Plane" << ipla <<  " has " << nwire
               << " wires: [" << setw(4) << firstwire << "," << setw(4) << lastwire << "]"
               << " and " << nchan << "/" << lastchan - firstchan + 1 << " channels: ["
               << setw(4) << firstchan << "," << setw(4) << lastchan << "]" << endl;
        }
        // Check if the range for the current readout plane (ROP) is already covered.
        Index irop;
        if ( m_dbg > 1 ) cout << myname << "Checking if ROP is covered." << endl;
        for ( irop=0; irop<m_nrop; ++irop ) {
          int ropfirst = m_ropfirstchan[irop];
          int roplast = ropfirst + m_ropnchan[irop] - 1;
          if ( m_dbg > 2 ) cout << myname << "  Checking range: " << ropfirst << "..." << roplast << endl;
          if ( firstchan > roplast ) continue;               // New range is after the ROP
          if ( lastchan < ropfirst ) continue;               // New range is before the ROP
          if ( firstchan==ropfirst && lastchan==roplast ) break;  // Exact match
          // We could extend the range here but current geometry does not require this.
          cout << myname << myname << "ERROR: Invalid channel range overlap." << endl;
          abort();
        }
        // If this is a new ROP, then record its channel range and TPC, and assign it to
        // an APA. For now, the latter assumes APA ordering follows TPC and is
        if ( irop == m_nrop ) {  // We have a new readout plane.
          if ( m_dbg > 1 ) cout << myname << "Adding new ROP: " << irop << endl;
          if ( m_dbg > 2 ) cout << myname << "Extending APA range: " << m_napa << "..." << iapa << endl;
          for ( unsigned int japa=m_napa; japa<=iapa; ++japa ) {
            nzplane = 0;               // There are not yet any z-planes in this APA
            m_apanrop.push_back(0);
            ++m_napa;
          }
          if ( m_dbg > 2 ) cout << myname << "Recording APA ROP " << iapa << "/" << m_apanrop.size() << endl;
          ++m_apanrop[iapa];
          if ( m_dbg > 2 ) cout << myname << "Recording ROP channel info." << endl;
          m_ropfirstchan.push_back(firstchan);
          m_ropnchan.push_back(lastchan - firstchan + 1 );
          if ( m_dbg > 2 ) cout << myname << "Recording ROP TPC." << endl;
          m_roptpc.push_back(itpc);
          m_ropapa.push_back(iapa);
          ostringstream ssrop;
          if ( m_dbg > 2 ) cout << myname << "Finding view." << endl;
          geo::View_t view = plageo.View();
          string sview = "?";
          string svcount = "";
          if ( view == kU ) sview = "u";
          if ( view == kV ) sview = "v";
          if ( view == kZ ) {
            sview = "z";
            ++nzplane;
            ostringstream ssvcount;
            ssvcount << nzplane;
            svcount = ssvcount.str();
          }
          m_ropview.push_back(view);
          ssrop << "apa" << iapa << sview;
          if ( view == kZ ) ssrop << nzplane;
          if ( m_dbg > 2 ) cout << myname << "Recording view " << ssrop.str() << endl;
          m_ropname.push_back(ssrop.str());
          ++m_nrop;
          ++itdcrop;
        } else {
          if ( m_dbg > 2 ) cout << myname << "No action for existing ROP." << endl;
        }
        // Add this TPC plane to the TPC-plane-to-ROP map.
        if ( m_dbg > 1 ) cout << myname << "Adding plane to ROP map." << endl;
        PlaneID pid(icry, icrytpc, ipla);
        if ( m_tpprop.find(pid) != m_tpprop.end() ) {
          cout << myname << "ERROR: Duplicate TPC plane." << endl;
          abort();
        }
        m_tpprop[pid] = irop;
        ++itpp;
        if ( m_dbg > 1 ) cout << myname << "End loop over planes." << endl;
      }  // End loop over planes in the TPC
      ++itpc;
      if ( m_dbg > 1 ) cout << myname << "End loop over TPCs." << endl;
    }  // End loop over TPCs in the cryostat
    if ( m_dbg > 1 ) cout << myname << "End loop over cryostats." << endl;
  }  // End loop over cryostats

  return 0;
}

//**********************************************************************
