// det35t.cxx

// David Adams
// August, 2015
//
// Displays 35-ton detector using hardwired geometry.

#include <iostream>
#include <string>
#include "LbTupler/GeoHelper.h"
#include "LbTupler/LArGeoManager.h"

using std::cout;
using std::endl;
using std::string;

int detlar(string gname ="lbne35t4apa_v5", bool useDefault =false, int print =false) {

  const string myname = "detlar: ";

  // Create geometry helper.
  if ( print ) cout << myname << "Creating geometry helper." << endl;
  GeoHelper gh(gname);
  if ( print ) gh.print();

  if ( useDefault ) {
    LArGeoManager lgm;
    lgm.draw();
  } else {
    LArGeoManager lgm(gh, 1);
    lgm.draw();
  }

  return 0;
}
