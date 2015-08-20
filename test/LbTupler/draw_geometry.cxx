// draw_geometry.cxx

// David Adams
// August 2015
//
// Draw DUNE geometries.

#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/make_ParameterSet.h"
#include "Geometry/GeometryCore.h"
#include "TROOT.h"
#include "TApplication.h"
#include "TRint.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TGeoManager.h"
#include "TEveManager.h"
#include "TEveGeoNode.h"
#include "LbTupler/LArGeoManager.h"

#include <string>
#include <iostream>
#include <sstream>
#include <vector>
#include <cassert>

using std::string;
using std::cout;
using std::endl;
using std::cin;
using std::ostringstream;
using std::vector;
using geo::TPCID;
using geo::TPCGeo;

class TpcBuilder {
public:
  double width;
  double height;
  double length;
  int ncopy = 1;
  TpcBuilder(double awidth, double aheight, double alength)
  : width(awidth), height(aheight), length(alength) { }
  void addTpc(TGeoVolume* ptop, int itpc, double xc, double yc, double zc) const {
    static TGeoMaterial* pmatLAr = new TGeoMaterial("LAr", 38.95,18,1.4);
    static TGeoMedium* pmed = new TGeoMedium("LAr" , 3, pmatLAr);
    ostringstream ssname;
    ssname << "TPC" << itpc;
    TGeoVolume* pvol = gGeoManager->MakeBox(ssname.str().c_str(), pmed, 0.5*width, 0.5*height, 0.5*length);
    TGeoTranslation* ptrans = new TGeoTranslation(xc, yc, zc);
    pvol->SetTransparency(75);
    pvol->SetLineColor(33);
    ptop->AddNode(pvol, ncopy, ptrans);
  }
};

int drawgeom(string gname) {
  const string myname = "drawgeom: ";
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
    return 1;
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
  cout << myname << "Displaying basic geometry info" << endl;
  cout << myname << "   Geometry name: " << detname << endl;
  cout << myname << "     # cryostats: " << ncry << endl;
  cout << myname << "    # TPC for C0: " << ntpc << endl;
  cout << myname << "       Surface Y: " << ysurf << " cm" << endl;
  return 0;
}

int main(int narg, char** argv) {
  string gname = "dune10kt_v1";
  if ( narg > 1 ) gname = argv[1];
  TRint* ptapp = new TRint("myapp", &narg, argv, nullptr, -1);
  drawgeom(gname);
  LArGeoManager* plgm = new LArGeoManager;
  plgm->draw();
  ptapp->Run(true);
  cout << "Terminating Root windowing." << endl;
  ptapp->Terminate(0);
  cout << "Exiting." << endl;
  return 0;
}
