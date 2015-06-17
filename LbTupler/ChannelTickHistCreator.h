// ChannelTickHistCreator.h

#ifndef ChannelTickHistCreator_H
#define ChannelTickHistCreator_h

// David Adams
// May 2015
//
// Class to create channel vs. tick histograms.

#include <string>
#include "Range.h"

class TH2;
namespace art {
class TFileDirectory;
}

class ChannelTickHistCreator {

public:  // typedefs

  typedef int Tick;
  typedef Range<Tick> TickRange;

public:  // methods

  // Ctor.
  //   ptfs = pointer to art TFile service. Used if non-null.
  //   sevt = Event ID string
  //   tick1, tick2 = Range for the x-axis
  //   zlab = z-axis lable, e.g. "ADC counts" or "Energy [MeV]"
  //   zmin, zmax = z-axis range
  //   ncontour = # contours for contour plotting
  ChannelTickHistCreator(art::TFileDirectory& tfs, std::string sevt, int tick1, int tick2,
                         std::string zlab, double zmin, double zmax, int ncontour);

  // Create a histogram for a given channel range.
  //   slab = Unique label, e.g. "hitapa2v" or "rawall"
  //   chan1, chan2 = y-axis range
  //   stitle Unique prefix for histogram title, e.g. "Hits for APA plane 2u"
  //   sevtNameSuffix = field appended to event ID in histogram name, e.g. "trk123"
  //   sevtTitleSuffix = field appended to event ID in histogram title, e.g. "track 123"
  //   tickRange - used to devine range of x-axis iff range has at least one tick
  TH2* create(std::string slab, unsigned int chan1, unsigned int chan2, std::string stitle,
              std::string sevtNameSuffix ="", std::string sevtTitleSuffix ="",
              TickRange tickRange =TickRange(0,-1));

private:  // data

  art::TFileDirectory& m_tfs;
  std::string m_sevt;
  TickRange m_tickRange;
  std::string m_zlab;
  double m_zmin;
  double m_zmax;
  int m_ncontour;

};

#endif
