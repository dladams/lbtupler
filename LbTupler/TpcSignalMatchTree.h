// TpcSignalMatchTree.h

// David Adams
// June 2015
//
// Defines and fills a tree describng a TpcSignalMatch.

#ifndef TpcSignalMatchTree_Module
#define TpcSignalMatchTree_Module

#include <string>

namespace art {
class Event;
}
class TTree;
class TpcSignalMatcher;

class TpcSignalMatchTree {

public:

  // Ctor.
  //   tname - Name for the Root tree.
  TpcSignalMatchTree(std::string tname);

  // Dtor.
  ~TpcSignalMatchTree();

  // Fill the tree with a match.
  int fill(const art::Event& evt, const TpcSignalMatcher& match);

private:

  // Control parameters.
  std::string m_tname;           // Tree name.

  // The tree.
  TTree* m_ptree;

  static const int maxent = 200;

  // The variables that will go into the n-tuple.
  int fevent;
  int fRun;
  int fSubRun;
  int fref;           // Reference index (e.g. MC index)
  int fmatch;         // Match index (e.g. cluster number)
  float fdistance;    // Math-reference distance
  int frop;           // ROP
  int frnbin;         // # channel-tick bins in the reference.
  int fmnbin;         // # channel-tick bins in the match.
  int frnseg;         // # segments attached to the reference.

}; // class TpcSignalMatchTree

#endif
