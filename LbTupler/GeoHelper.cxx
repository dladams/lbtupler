// GeoHelper.cxx

#include "GeoHelper.h"
#include "Geometry/Geometry.h"
#include "Utilities/DetectorProperties.h"
#include "Utilities/LArProperties.h"
#include "Geometry/TPCGeo.h"
#include "Geometry/PlaneGeo.h"
#include <iostream>
#include <sstream>
#include <iomanip>
#include <set>

#undef LARSOFT0401

using std::ostream;
using std::cout;
using std::endl;
using std::ostringstream;
using std::setw;
using std::vector;
using std::set;
using std::string;
using util::LArProperties;
using util::DetectorProperties;
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

GeoHelper::GeoHelper(const geo::Geometry* pgeo, Status dbg)
: m_pgeo(pgeo), m_dbg(dbg), m_ntpc(0), m_ntpp(0), m_napa(0), m_nrop(0) {
  // Find the total number of TPC.
  for ( Index icry=0; icry<ncryostat(); ++icry ) m_ntpc += m_pgeo->NTPC(icry);
  // Fill the TPC info w/o ROP or APA.
  unsigned int iglotpc = 0;
  for ( Index icry=0; icry<ncryostat(); ++icry ) {
    for ( Index icrytpc=0; icrytpc<m_pgeo->NTPC(icry); ++icrytpc ) {
      m_ntpp += m_pgeo->Nplanes(icrytpc, icry);
      m_tpccry.push_back(icry);
      ostringstream ssname;
      ssname << "TPC";
      if ( m_ntpc > 99 && iglotpc < 100 ) ssname << "0";
      if ( m_ntpc > 9 && iglotpc < 10 ) ssname << "0";
      ssname << iglotpc;
      m_tpcname.push_back(ssname.str());
      ++iglotpc;
    }
  }
  fillStandardApaMapping();
}

//**********************************************************************

DetectorProperties& GeoHelper::detectorProperties() const {
  static DetectorProperties& ref = *art::ServiceHandle<DetectorProperties>();
  return ref;
}

//**********************************************************************

LArProperties& GeoHelper::larProperties() const {
  static LArProperties& ref = *art::ServiceHandle<LArProperties>();
  return ref;
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

PlanePositionVector GeoHelper::planePositions(const double postim[]) const {
  const string myname = "LbTupler::planePositions: ";
  PlanePositionVector pps;
  if ( geometry() == nullptr ) {
    cout << myname << "ERROR: Geometry is null." << endl;
    return pps;
  }
  double samplingrate = detectorProperties().SamplingRate();
  if ( samplingrate <= 0 ) {
    cout << myname << "ERROR: Invalid sampling rate." << endl;
    return pps;
  }
  geo::TPCID tpcid = geometry()->FindTPCAtPosition(postim);
  if ( ! tpcid.isValid ) return pps;
  unsigned int nplane = geometry()->Nplanes(tpcid.TPC, tpcid.Cryostat);
  for ( unsigned int ipla=0; ipla<nplane; ++ipla ) {
    PlaneID tpp(tpcid, ipla);
    unsigned int irop = rop(tpp);
    double tick = detectorProperties().ConvertXToTicks(postim[0], ipla, tpcid.TPC, tpcid.Cryostat);
    double tickoff = postim[3]/samplingrate;
    tick += tickoff;
    if ( m_dbg > 2 ) cout << myname << "Tick offset: " << tickoff << " = " << postim[3] << "/" << samplingrate << endl;
    unsigned int ichan = geometry()->NearestChannel(postim, ipla, tpcid.TPC, tpcid.Cryostat);
    unsigned int iropchan = ichan - ropFirstChannel(irop);
    int itick = int(tick);
    if ( tick < 0.0 ) itick += -1;
    PlanePosition pp;
    pp.plane = ipla;
    pp.rop = irop;
    pp.channel = ichan;
    pp.tick = tick;
    pp.itick = itick;
    pp.ropchannel = iropchan;
    pp.valid = true;
    pps.push_back(pp);
  }
  return pps;
}

//**********************************************************************

ostream& GeoHelper::print(ostream& out, int iopt, std::string prefix) const {
  if ( m_pgeo == nullptr ) {
    out << prefix << "Geometry is not defined." << endl;
    return out;
  }
  cout << prefix << "             Total # TPCs: " << ntpc() << endl;
  cout << prefix << "     Total # TPC channels: " << geometry()->Nchannels() << endl;
#ifdef LARSOFT0401
  cout << prefix << "Total # optical detectors: " << geometry()->NOpDet() << endl;
#else
  cout << prefix << "Total # optical detectors: " << geometry()->NOpDets() << endl;
#endif
  cout << prefix << " Total # optical channels: " << geometry()->NOpChannels() << endl;
  cout << prefix << endl;
  cout << prefix << "There are " << ntpc() << " TPCs:" << endl;
  cout << prefix << "      name       APA" << endl;
  for ( unsigned int itpc=0; itpc<ntpc(); ++itpc ) {
    cout << prefix << setw(10) << tpcName(itpc) << setw(10) << tpcApa(itpc) << endl;
  }
  cout << prefix << endl;
  cout << prefix << "There are " << nrop() << " ROPs (readout planes):" << endl;
  cout << prefix << "      name  1st chan     #chan  orient" << endl;
  for ( unsigned int irop=0; irop<nrop(); ++irop ) {
    cout << prefix << setw(10) << ropName(irop) << setw(10) << ropFirstChannel(irop)
         << setw(10) << ropNChannel(irop) << setw(8) << ropView(irop) << endl;
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
#ifdef LARSOFT0401
  out << prefix << setw(wlab) <<"Optical channel count: " << m_pgeo->NOpDet() << endl;
#else
  out << prefix << setw(wlab) <<"Optical channel count: " << m_pgeo->NOpDets() << endl;
#endif
  out << prefix << setw(wlab) <<"Scintillator channel count: " << m_pgeo->NAuxDets() << endl;
  out << prefix << "TPC sizes [cm]:" << endl;
  out << prefix << setw(wcry) << "Cry" << setw(wtpc) << "TPC"
      << setw(wdim) << "Width" << setw(wdim) << "Height"<< setw(wdim) << "Length"
      << setw(wdim) << "   xc" << setw(wdim) << "    yc"<< setw(wdim) << "    zc"
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
  cout << prefix << setw(4) << "ROP" << setw(12) << "Name"
       << setw(5) << "View" << setw(8) << "1st ch" << setw(6) << "Nchan" << setw(6) << "TPCs" << endl;
  for ( Index irop=0; irop<nrop(); ++irop ) {
    ostringstream ssrops;
    for ( Index itpc : ropTpcs(irop) ) {
      if ( ssrops.str().size() ) ssrops << ",";
      ssrops << itpc;
    }
    cout << prefix << setw(4) << irop << setw(12) << ropName(irop)
         << setw(5) << ropView(irop) << setw(8) << ropFirstChannel(irop)
         << setw(6) << ropNChannel(irop) << setw(6) << ssrops.str();
    cout << endl;
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
  // Show LAr properties.
  double driftspeed = larProperties().DriftVelocity();
  double samplingrate = detectorProperties().SamplingRate();
  double samplingdriftspeed = driftspeed*samplingrate/1000.0;
  if ( m_dbg > 0 ) cout << prefix << "LAr drift speed = " << driftspeed << " cm/us" << endl;
  if ( m_dbg > 0 ) cout << prefix << "Sampling period = " << samplingrate << " ns" << endl;
  if ( m_dbg > 0 ) cout << prefix << "Sampling speed = " << samplingdriftspeed << " cm/tick" << endl;

  return out;
}
  
//**********************************************************************

// This geometry dump code is from Michelle Stancari.
ostream& GeoHelper::dump(ostream& out) const {
  const string myname = "GeoHelper::dump: ";
  double xyz[3]; 
  double abc[3];
  int chan;
  int cryo = geometry()->Ncryostats();
  for (int c=0; c<cryo; ++c){
    int tpc = geometry()->NTPC(c);
    for (int t=0; t<tpc; ++t){
      int Nplanes = geometry()->Nplanes(t,c);
      for (int p=0;p<Nplanes;++p) {
        int Nwires = geometry()->Nwires(p,t,c);
        out << "FLAG " << endl;
        for (int w=0; w<Nwires; ++w){
          geometry()->WireEndPoints(c,t,p,w,xyz,abc);
          chan = geometry()->PlaneWireToChannel(p,w,t,c);
          out << myname << setw(4) << chan << " " << c << " " << t << " " << p << setw(4) << w << " "
                    << xyz[0] << " " << xyz[1] << " " << xyz[2] <<  " "
                    << abc[0] << " " << abc[1] << " " << abc[2] << endl;
        }
      }
    }
  }
  return out;
}

//**********************************************************************

Index GeoHelper::channelRop(Index chan) const {
  const string myname = "GeoHelper::channelRop: ";
  Index irop = nrop();
  for ( irop=0; irop<nrop(); ++irop ) {
    unsigned int lastchan = ropFirstChannel(irop) + ropNChannel(irop) - 1;
    if ( chan <= lastchan ) return irop;
  }
  return irop;
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
    if ( m_dbg > 1 ) cout << myname << "Begin cryostat " << icry << endl;
    // Loop over TPCs in the current cryostat.
    for ( Index icrytpc=0; icrytpc<m_pgeo->NTPC(icry); ++icrytpc ) {
      if ( m_dbg > 1 ) cout << myname << "Begin cryostat " << icry << ",  TPC " << icrytpc << endl;
      const TPCGeo& tpcgeo = m_pgeo->TPC(icrytpc, icry);
      unsigned int iapa = itpc/2;   // Standard mapping: one APA for each adjacent pair of TPCs
      int nplane = m_pgeo->Nplanes(icrytpc, icry);
      if ( nplane == 0 ) {
        cout << myname << "WARNING: Cryostat " << icry << ", TPC " << icrytpc << " has no planes." << endl;
      }
      int itdcrop = 0;   // # readouts for this TDC
      m_tpcapa.push_back(iapa);
      // Loop over planes in the TPC.
      for ( int ipla=0; ipla<nplane; ++ipla ) {
        if ( m_dbg > 1 ) cout << myname << "Begin TPC plane " << ipla << endl;
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
          int icha = m_pgeo->PlaneWireToChannel(ipla, iwir, icrytpc, icry);
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
          if ( m_dbg > 0 ) {
            cout << myname << "Plane" << ipla <<  " has " << nwire
                 << " wires: [" << setw(4) << firstwire << "," << setw(4) << lastwire << "]"
                 << " and " << nchan << "/" << lastchan - firstchan + 1 << " channels: ["
                 << setw(4) << firstchan << "," << setw(4) << lastchan << "]" << endl;
          }
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
          IndexVector tpcs;
          tpcs.push_back(itpc);
          m_roptpc.push_back(tpcs);
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
        // ROP already exists--find the ROP for this plane.
        } else {
          Index irop = channelRop(firstchan);
          Index irop2 = channelRop(lastchan);
          if ( irop2 != irop ) {
            cout << myname << "ERROR: TPC plane covers multiple ROPs." << endl;
          }
          IndexVector& tpcs = m_roptpc[irop];
          if ( std::find(tpcs.begin(), tpcs.end(), itpc) == tpcs.end() ) {
            if ( m_dbg > 2 ) cout << myname << "Add TPC " << itpc << " to ROP " << irop << endl;
            tpcs.push_back(itpc);
          }
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
      }  // End loop over planes in the TPC
      ++itpc;
      if ( m_dbg > 1 ) cout << myname << "End TPC " << itpc << " (cryostat " << icry
                            << ", TPC " << icrytpc << ")" << endl;
    }  // End loop over TPCs in the cryostat
    if ( m_dbg > 1 ) cout << myname << "End cryostat " << icry << endl;
  }  // End loop over cryostats
  if ( m_dbg > 1 ) cout << myname << "End loop over cryostats." << endl;

  return 0;
}

//**********************************************************************
