// TpcSignalMatchTree.cxx

#include "TpcSignalMatchTree.h"

#include <iostream>

// Art includes
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"

// ROOT includes.
#include "TTree.h"

// Local includes.
#include "TpcSignalMatcher.h"
#include "TpcTypes.h"

using std::cout;
using std::endl;
using std::string;
using tpc::Index;
using tpc::badIndex;

//************************************************************************

TpcSignalMatchTree::TpcSignalMatchTree(string tname)
: m_tname(tname) {

  const string myname = "TpcSignalMatchTree: ";
  const int dbg = 0;

  // Access ART's TFileService, which will handle creating and writing
  // histograms and n-tuples for us. 
  art::ServiceHandle<art::TFileService> tfs;

  // Define the tree.
  m_ptree = tfs->make<TTree>(tname.c_str(), m_tname.c_str());

  // Define the branches (columns) of the tree.
  m_ptree->Branch("event",       &fevent,          "event/I");
  m_ptree->Branch("run",         &fRun,            "run/I");
  m_ptree->Branch("subrun",      &fSubRun,         "subrun/I");
  m_ptree->Branch("ref",         &fref,            "ref/I");       // Ref index
  m_ptree->Branch("match",       &fmatch,          "match/I");     // Match index for
  m_ptree->Branch("distance",    &fdistance,       "distance/F");  // distance for match
  m_ptree->Branch("rop",         &frop,            "rop/I");       // read out plane index
  m_ptree->Branch("rnbin",       &frnbin,          "rnbin/I");     // # bins in ref
  m_ptree->Branch("mnbin",       &fmnbin,          "mnbin/I");     // # bins in ref
  m_ptree->Branch("rnseg",       &frnseg,          "rnseg/I");     // # bins in ref

  if ( dbg > 0 ) {
    cout << myname << "Initialization complete." << endl;
    cout << myname << "   Tree: " << m_ptree->GetName() << endl;
  }
}
 
//************************************************************************

// Destructor
TpcSignalMatchTree::~TpcSignalMatchTree() { }
   
//************************************************************************

int TpcSignalMatchTree::fill(const art::Event& evt, const TpcSignalMatcher& match) {
  const string myname = "TpcSignalMatchTreevent: ";

  fevent  = evt.id().event(); 
  fRun    = evt.run();
  fSubRun = evt.subRun();
  cout << myname << "            Match size: " << match.size() << endl;
  cout << myname << " Reference vector size: " << match.referenceVector().size() << endl;
  cout << myname << "     Match vector size: " << match.matchVector().size() << endl;
  for ( unsigned int ient =0; ient<match.size(); ++ ient ) {
    fref = ient;
    fmatch = match.matchIndex(ient);
    fdistance = match.matchDistance(ient);
    Index imat = match.matchIndex(ient);
    if ( imat == badIndex() ) {
      frop = badIndex();
      fmnbin = 0;
      frnbin = 0;
      frnseg = 0;
    } else {
      const TpcSignalMap& rtsm = *match.referenceVector().at(ient);
      const TpcSignalMap& mtsm = *match.matchVector().at(imat);
      frop = rtsm.rop();
      fmnbin = mtsm.binCount();
      frnbin = rtsm.binCount();
      frnseg = rtsm.segments().size();
    }
    m_ptree->Fill();
  }


  return 0;
}

//**********************************************************************
