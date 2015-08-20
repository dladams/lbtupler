// det35t.cxx

// David Adams
// August, 2015
//
// Displays 35-ton detector using hardwired geometry.

#include <iostream>
#include <sstream>
#include <string>
#include <memory>
#include "TGeoManager.h"
#include "TEveManager.h"
#include "TEveGeoNode.h"
#ifndef __CINT__
#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/make_ParameterSet.h"
#include "Geometry/GeometryCore.h"
#include "Geometry/ChannelMapAlg.h"
#include "lbne/Geometry/ChannelMap35OptAlg.h"
#endif
#include "LbTupler/GeoHelper.h"

using std::cout;
using std::endl;
using std::string;
using std::ostringstream;
using std::shared_ptr;
using geo::TPCID;
using geo::TPCGeo;

class LArDetector {
public:
  TGeoManager* pgeo = nullptr;
  LArDetector() {
    pgeo = gGeoManager;
  }
  LArDetector(string name, string title) {
    pgeo = new TGeoManager(name.c_str(), title.c_str());
    // Create materials.
    TGeoMaterial* pmatVac = new TGeoMaterial("Vacuum", 0,0,0);
    TGeoMaterial* pmatAl = new TGeoMaterial("Al", 26.98,13,2.7);
    TGeoMaterial* pmatLAr = new TGeoMaterial("LAr", 38.95,18,1.4);
    // Create media.
    TGeoMedium* pvac = new TGeoMedium("Vacuum",1, pmatVac);
    TGeoMedium* pAl = new TGeoMedium("Aluminium",2, pmatAl);
    TGeoMedium* pLAr = new TGeoMedium("LAr" , 3, pmatLAr);
    string vname = "TOP_" + name;
    pgeo->SetTopVolume(pgeo->MakeBox(vname.c_str(), pvac, 1000., 1000., 1000.));
    pgeo->GetTopVolume()->SetTransparency(75);
  }
  void addTpc(int itpc, double dx, double dy, double dz, double xc, double yc, double zc) {
    ostringstream ssname;
    ssname << "TPC" << itpc;
    TGeoMedium* pmed = pgeo->GetMedium("LAr");
    TGeoVolume* pvol = pgeo->MakeBox(ssname.str().c_str(), pmed, 0.5*dx, 0.5*dy, 0.5*dz);
    TGeoTranslation* ptrans = new TGeoTranslation(xc, yc, zc);
    pvol->SetTransparency(75);
    pvol->SetLineColor(33);
    pgeo->GetTopVolume()->AddNode(pvol, 0, ptrans);
  }
  void addTpc(const TGeoVolume* apvol, double xc, double yc, double zc) {
    string oldname = apvol->GetName();
    string newname = oldname + "Copy";
    TGeoVolume* pvol = dynamic_cast<TGeoVolume*>(apvol->Clone(newname.c_str()));
    TGeoTranslation* ptrans = new TGeoTranslation(xc, yc, zc);
    //pvol->SetTransparency(75);
    //pvol->SetLineColor(33);
    pgeo->GetTopVolume()->AddNode(pvol, 0, ptrans);
  }
  void draw() const {
    // Close and display the detector.
    pgeo->CloseGeometry();
    delete gEve;
    gEve = nullptr;
    TEveManager::Create();
    TGeoNode* ptopnode = pgeo->GetTopNode();
    TEveGeoTopNode* pevetopnode = new TEveGeoTopNode(pgeo, pgeo->GetTopNode());
    gEve->AddGlobalElement(pevetopnode);
    gEve->Redraw3D(kTRUE);
  }
};

int detlar(string gname ="lbne35t4apa_v5", bool useDefault =false) {

  const string myname = "detlar: ";

  // Configure the geometry.
  string spar = "Name: \"" + gname + "\"\nDisableWiresInG4: true\nSurfaceY: 0";
  cout << myname << "Config: \n" << spar << endl;
  fhicl::ParameterSet parset;
  fhicl::make_ParameterSet(spar, parset);
  geo::GeometryCore* pdetgeo = new geo::GeometryCore(parset);

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
  pdetgeo->LoadGeometryFile(fullgdmlfile, fullrootfile);

  // Add the geometry channel map.
  bool useChannels = false;
  if ( useChannels ) {
    cout << myname << "Adding channel map." << endl;
    spar = "SortingParameters: {}  # to use default";
    fhicl::make_ParameterSet(spar, parset);
    shared_ptr<geo::ChannelMapAlg> pChannelMap(new geo::ChannelMap35OptAlg(parset));
    pdetgeo->ApplyChannelMap(pChannelMap);
  }

  // Create geometry helper.
  cout << myname << "Creating geometry helper." << endl;
  GeoHelper gh(pdetgeo, useChannels);
  gh.print();

  if ( useDefault ) {
    LArDetector defdet;
    defdet.draw();
    return 0;
  }

  // Make the top detector volume.
  LArDetector det("mydet", "mydet");

  // Loop over TPCs
  unsigned int ncry = pdetgeo->Ncryostats();
  cout << myname << "   Geometry name: " << pdetgeo->DetectorName() << endl;
  cout << myname << "     # cryostats: " << ncry << endl;
  unsigned int itpc = 0;
  for ( unsigned int icry=0; icry<ncry; ++icry ) {
    unsigned int ntpc = pdetgeo->NTPC(icry);
    cout << myname << " Cryostat " << icry << " # TPC: " << ntpc << endl;
    for ( unsigned int icrytpc=0; icrytpc<ntpc; ++icrytpc ) {
      //cout << "  icry=" << icry << ", icrytpc=" << icrytpc << endl;
      TPCID tpcid(icry, icrytpc);
      const TPCGeo* ptpcgeo = pdetgeo->TPCPtr(tpcid);
      const TGeoVolume* pvol = ptpcgeo->TotalVolume();
      double pos[3];
      gh.tpcCenter(icry, icrytpc, pos);
      det.addTpc(itpc, gh.width(icry, icrytpc), gh.height(icry, icrytpc),  gh.length(icry, icrytpc),
                 pos[0], pos[1], pos[2]);
      ++itpc;
    }
  }

  det.draw();

  return 0;
}
