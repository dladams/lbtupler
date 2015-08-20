// LArGeoManager.cxx

#include "LArGeoManager.h"

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

//**********************************************************************

LArGeoManager::LArGeoManager() {
  m_pgeo = gGeoManager;
}

//**********************************************************************

LArGeoManager::LArGeoManager(const GeoHelper& gh, bool addTpcs) {
  m_pgh = &gh;
  if ( m_pgh->geometry() == nullptr ) return;
  string name = "largeo_" + m_pgh->geometry()->DetectorName();
  string title = "LAr geometry " + m_pgh->geometry()->DetectorName();
  m_pgeo = new TGeoManager(name.c_str(), title.c_str());
  // Create materials.
  TGeoMaterial* pmatVac = new TGeoMaterial("Vacuum", 0,0,0);
  TGeoMaterial* pmatAl = new TGeoMaterial("Al", 26.98,13,2.7);
  TGeoMaterial* pmatLAr = new TGeoMaterial("LAr", 38.95,18,1.4);
  // Create media.
  TGeoMedium* pvac = new TGeoMedium("Vacuum",1, pmatVac);
  new TGeoMedium("Aluminium",2, pmatAl);
  new TGeoMedium("LAr" , 3, pmatLAr);
  string vname = "TOP_" + name;
  m_pgeo->SetTopVolume(m_pgeo->MakeBox(vname.c_str(), pvac, 1000., 1000., 1000.));
  m_pgeo->GetTopVolume()->SetTransparency(75);
  if ( addTpcs ) {
    unsigned int ncry = m_pgh->geometry()->Ncryostats();
    for ( unsigned int icry=0; icry<ncry; ++icry ) {
      unsigned int ntpc = m_pgh->geometry()->NTPC(icry);
      for ( unsigned int icrytpc=0; icrytpc<ntpc; ++icrytpc ) {
        addTpc(icry, icrytpc);
      }
    }
  }
}

//**********************************************************************

void LArGeoManager::
addTpc(int itpc, double dx, double dy, double dz, double xc, double yc, double zc) const {
  if ( m_pgh == nullptr ) return;
  if ( m_pgeo == nullptr ) return;
  ostringstream ssname;
  ssname << "TPC" << itpc;
  TGeoMedium* pmed = m_pgeo->GetMedium("LAr");
  TGeoVolume* pvol = m_pgeo->MakeBox(ssname.str().c_str(), pmed, 0.5*dx, 0.5*dy, 0.5*dz);
  TGeoTranslation* ptrans = new TGeoTranslation(xc, yc, zc);
  pvol->SetTransparency(75);
  pvol->SetLineColor(33);
  m_pgeo->GetTopVolume()->AddNode(pvol, 0, ptrans);
}

//**********************************************************************

void LArGeoManager::addTpc(int icrytpc, int icry) const {
  if ( m_pgh == nullptr ) return;
  if ( m_pgh->geometry() == nullptr ) return;
  const GeoHelper& gh = *m_pgh;
  int itpc = icrytpc;
  int ncry = m_pgh->geometry()->Ncryostats();
  for ( int icry=0; icry<ncry; ++icry ) {
    itpc += m_pgh->geometry()->NTPC(icry);
  }
  double dx = gh.width(icry, icrytpc);
  double dy = gh.height(icry, icrytpc);
  double dz = gh.length(icry, icrytpc);
  double pos[3];
  gh.tpcCenter(icry, icrytpc, pos);
  addTpc(itpc, dx, dy, dz, pos[0], pos[1], pos[2]);
}

//**********************************************************************

void LArGeoManager::draw() const {
  // Close and display the detector.
  m_pgeo->CloseGeometry();
  delete gEve;
  gEve = nullptr;
  TEveManager::Create();
  TEveGeoTopNode* pevetopnode = new TEveGeoTopNode(m_pgeo, m_pgeo->GetTopNode());
  gEve->AddGlobalElement(pevetopnode);
  gEve->Redraw3D(kTRUE);
}

//**********************************************************************
