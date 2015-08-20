// LArGeoManager.h

// David Adams
// August, 2015
//
// Class to display LArSoft detectors.

class TGeoManager;

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

class LArGeoManager {

public:

  // Default ctor uses the current TGeoManager.
  // Presumably that created when the LAr geometry is loaded.
  LArGeoManager();

  // Ctor from GeoHelper.
  // If addTpcs is true, thena volume is created for each TPC.
  LArGeoManager(const GeoHelper& gh, bool addTpcs =true);

  // Add a TPC volume with the specified index, size and center position.
  void addTpc(int itpc, double dx, double dy, double dz, double xc, double yc, double zc) const;

  // Add a volume for the specified TPC.
  void addTpc(int itpc, int icry =0) const;

  // Display the detector.
  void draw() const;

private:

  TGeoManager* m_pgeo = nullptr;
  const GeoHelper* m_pgh = nullptr;
  char m_tranparency = 75;
  int m_color = 33;

};
