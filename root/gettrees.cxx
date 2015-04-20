// gettrees.cxx

#include "gettrees.h"
#include <iostream>
#include <iomanip>
#include "TDirectory.h"
#include "TTree.h"

using std::string;
using std::cout;
using std::endl;
using std::hex;
using std::dec;

TTree* simtree(string tname) {
  const string myname = "simtree: ";
  static TTree* ptree = nullptr;
  if ( tname.size() ) {
    if ( tname == "null" ) ptree = nullptr;
    else {
      TObject* pobj = nullptr;
      gDirectory->GetObject(tname.c_str(), pobj);
      if ( pobj == nullptr ) {
        cout << myname << "Object " << tname << " not found in " << gDirectory->GetName() << endl;
        ptree = nullptr;
      } else {
        ptree = dynamic_cast<TTree*>(pobj);
        if ( pobj == nullptr ) {
          cout << myname << "Object " << tname << " is not a TTree." << endl;
        }
      }
    }
    cout << myname << "Tree set to " << hex << long(ptree) << dec << endl;
  }
  return ptree;
}

