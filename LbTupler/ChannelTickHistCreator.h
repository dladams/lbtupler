// ChannelTickHistCreator.h

#ifndef ChannelTickHistCreator_H
#define ChannelTickHistCreator_h

// David Adams
// May 2015
//
// Class to create channel vs. tick histograms.

#include <string>

class TH2;
namespace art {
class TFileService;
}

class ChannelTickHistCreator {

public:  // methods

  // Ctor.
  //   ptfs = pointer to art TFile service. Used if non-null.
  //   sevt = Event ID string
  //   tick1, tick2 = Range for the x-axis
  //   zlab = z-axis lable, e.g. "ADC counts" or "Energy [MeV]"
  //   zmin, zmax = z-axis range
  //   ncontour = # contours for contour plotting
  ChannelTickHistCreator(art::TFileService* ptfs, std::string sevt, int tick1, int tick2,
                         std::string zlab, double zmin, double zmax, int ncontour);

  // Create a histogram for a given channel range.
  //   slab = Unique label, e.g. "hitapa2v" or "rawall"
  //   chan1, chan2 = y-axis range
  //   stitle Unique prefix for histogram title, e.g. "Hits for APA plane 2u"
  //   sevtNameSuffix = field appended to event ID in histogram name, e.g. "trk123"
  //   sevtTitleSuffix = field appended to event ID in histogram title, e.g. "track 123"
  TH2* create(std::string slab, unsigned int chan1, unsigned int chan2, std::string stitle,
              std::string sevtNameSuffix ="", std::string sevtTitleSuffix ="");

private:  // data

  art::TFileService* m_ptfs;
  std::string m_sevt;
  int m_tick1;
  int m_tick2;
  std::string m_zlab;
  double m_zmin;
  double m_zmax;
  int m_ncontour;

};

#endif
