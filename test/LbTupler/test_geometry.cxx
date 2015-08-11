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
#include <cassert>

using std::string;
using std::cout;
using std::endl;

int main() {
  //string geopath = "services.Geometry";
  const string myname = "test_geometry: ";
  string spar;
  string sep = "\n";
  spar += "{";
  spar += sep;
  spar += "  DisableWiresInG4: true";
  spar += sep;
  spar += "  GDML: \"lbne10kt.gdml\"";
  spar += sep;
  spar += "  Name: \"lbne10kt\"";
  spar += sep;
  spar += "  ROOT: \"lbne10kt.gdml\"";
  spar += sep;
  spar += "  SurfaceY: 0";
  spar += sep;
  spar += "  service_type: \"Geometry\"";
  spar += sep;
  spar += "}";
  cout << myname << "Config: \n" << spar << endl;
  fhicl::ParameterSet parset;
  fhicl::make_ParameterSet(spar, parset);
  geo::GeometryCore* pgeo = new geo::GeometryCore(parset);
  cout << myname << "Geometry name: " << pgeo->DetectorName() << endl;
  assert(false);
  cout << myname << "Ending test" << endl;
  return 0;
}
