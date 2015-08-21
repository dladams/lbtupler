// draw_geometry.cxx

// David Adams
// August 2015
//
// Draw DUNE geometries.

#include "TApplication.h"
#include "TRint.h"
#include "LbTupler/LArGeoManager.h"

#include <string>
#include <iostream>

using std::string;
using std::cout;
using std::endl;

int main(int narg, char** argv) {
  string gname = "dune10kt_v1";
  bool useDefault = false;
  if ( narg > 1 ) gname = argv[1];
  if ( narg > 2 ) useDefault = true;
  //TRint* ptapp = new TRint("myapp", &narg, argv, nullptr, -1);
  TApplication* ptapp = new TApplication("myapp", &narg, argv, nullptr, -1);
  GeoHelper gh(gname);
  if ( useDefault ) {
    LArGeoManager lgm;
    lgm.draw();
  } else {
    LArGeoManager lgm(gh, 1);
    lgm.draw();
  }
  cout << "To quit, use Browser->\"Quit Root\" quot or Command: .q" << endl;
  ptapp->Run(true);
  ptapp->Terminate(0);
  cout << endl;   // Root 5.34 leaves an unflushed line.
  return 0;
}
