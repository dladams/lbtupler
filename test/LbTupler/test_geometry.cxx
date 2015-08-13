// test_geometry.cxx

// David Adams
// August 2015
//
// Test script for DUNE geometries.

#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/make_ParameterSet.h"
#include "Geometry/GeometryCore.h"

#include <string>
#include <iostream>
#include <sstream>
#include <vector>
#include <cassert>

using std::string;
using std::cout;
using std::endl;
using std::ostringstream;
using std::vector;
using geo::TPCID;
using geo::TPCGeo;

string spos(double pos[3]) {
  ostringstream ssout;
  ssout << "(" << pos[0] << ", " << pos[1] << ", " << pos[2] << ")";
  return ssout.str();
}

bool inrange(double x, double x1, double x2) {
  if ( x2 > x1 ) {
    if ( x < x1 ) return false;
    if ( x > x2 ) return false;
  } else {
    if ( x < x2 ) return false;
    if ( x > x1 ) return false;
  }
  return true;
}

int main(int narg, char** argv) {
  //string geopath = "services.Geometry";
  const string myname = "test_geometry: ";
  string gname = "lbne10kt";
  if ( narg > 1 ) {
    gname = argv[1];
  }
  vector<string> gnames = {"lbne10kt", "lbne35t"};
  // (xref,yref,zref) is a point inside a TPC
  vector<double> xref = {-500.0, 0.0};
  vector<double> yref = {   0.0, 0.0};
  vector<double> zref = {   0.0, 0.0};
  vector<int> ncryExp = {2, 1};
  vector<int> ntpcExp = {120, 6};
  cout << "Checking geometry name: " << gname << endl;
  unsigned int ndet = gnames.size();
  unsigned int idet = 0;
  if ( gname == "-h" || gname == "help" ) {
    cout << "Usage: " << argv[0] << " [GEONAME]" << endl;
    return 0;
  }
  for ( ; idet<ndet; ++idet ) {
    if ( gnames[idet] == gname ) break;
  }
  assert( idet != ndet );
  string spar = R"(Name: ")" + gname + R"("
DisableWiresInG4: true
SurfaceY: 0
)";
  cout << myname << "Config: \n" << spar << endl;
  fhicl::ParameterSet parset;
  fhicl::make_ParameterSet(spar, parset);
  geo::GeometryCore* pgeo = new geo::GeometryCore(parset);
  // Load the geometry.
  string gdmlfile = gname + ".gdml";
  string rootfile = gdmlfile;
  string fullgdmlfile = gdmlfile;
  string fullrootfile;
  cet::search_path sp("FW_SEARCH_PATH");
  if ( ! sp.find_file(rootfile, fullrootfile) ) {
    cout << myname << "Unable to find the root geometry file: " << rootfile << endl;
    return 2;
  }
  cout << myname << "GDML file: " << fullgdmlfile << endl;
  cout << myname << "ROOT file: " << fullrootfile << endl;
  pgeo->LoadGeometryFile(fullgdmlfile, fullrootfile);
  // Test the geometry.
  cout << myname << "Retrieving basic geometry info" << endl;
  int ncry = pgeo->Ncryostats();
  int ntpc = pgeo->NTPC(0);
  string detname = pgeo->DetectorName();
  double ysurf = pgeo->SurfaceY();
  double totmass = pgeo->TotalMass();
  double pos[3] = {xref[idet], yref[idet], zref[idet]};
  TPCID tid = pgeo->FindTPCAtPosition(pos);
  assert(tid.isValid);
  TPCGeo tgeo = pgeo->TPC(tid);
  cout << myname << "Displaying basic geometry info" << endl;
  cout << myname << "   Geometry name: " << detname << endl;
  cout << myname << "     # cryostats: " << ncry << endl;
  cout << myname << "    # TPC for C0: " << ntpc << endl;
  cout << myname << "       Surface Y: " << ysurf << " cm" << endl;
  cout << myname << "   Detector mass: " << totmass << " kg" << endl;
  cout << myname << "   Reference pos: " << "(" << pos[0] << ", " << pos[1]
                                         << ", " << pos[2] << ")" << endl;
  cout << myname << "   Reference cry: " << tid.Cryostat << endl;
  cout << myname << "   Reference TPC: " << tid.TPC << endl;
  cout << myname << "Checking basic geometry info" << endl;
  assert(detname==gname);
  assert(ncry==ncryExp[idet]);
  assert(ntpc==ntpcExp[idet]);
  cout << myname << "Ended check of basic geometry" << endl;
  // Check TPC positions.
  cout << myname << "Checking TPC position" << endl;
  //double pos[3] = {-303.606, 587.567, 10.9497};
  double locpos2[3] = {tgeo.HalfWidth(), tgeo.HalfHeight(), 0.5*tgeo.Length()};
  double locpos1[3] = {-locpos2[0], -locpos2[1], -locpos2[2]};
  double pos1[3], pos2[3];
  double locpos[3];
  tgeo.LocalToWorld(locpos1, pos1);
  tgeo.LocalToWorld(locpos2, pos2);
  tgeo.WorldToLocal(pos, locpos);
  cout << myname << "World coordinates:" << endl;
  cout << myname << "         point: " << spos(pos) << endl;
  cout << myname << "  TPC corner 1: " << spos(pos1) << endl;
  cout << myname << "  TPC corner 2: " << spos(pos2) << endl;
  cout << myname << "TPC coordinates:" << endl;
  cout << myname << "         point: " << spos(locpos) << endl;
  cout << myname << "  TPC corner 1: " << spos(locpos1) << endl;
  cout << myname << "  TPC corner 2: " << spos(locpos2) << endl;
  assert( inrange(locpos[0], locpos1[0], locpos2[0]) );
  assert( inrange(locpos[1], locpos1[1], locpos2[1]) );
  assert( inrange(locpos[2], locpos1[2], locpos2[2]) );
  cout << myname << "Done." << endl;
  return 0;
}
