// SimChannelTupler_module.cc

#ifndef SimChannelTupler_Module
#define SimChannelTupler_Module

// David Adams
// May 2015
//
// Defines a Root tree that describes SimChannel data.

#include <vector>

class GeoHelper;
namespace art {
class Event;
class TFileService;
}
namespace sim {
class SimChannel;
}
class TTree;

class SimChannelTupler {

public:

  typedef std::vector<sim::SimChannel> SimChannelVector;
 
public:
 
  // Standard constructor and destructor for an ART module.
  SimChannelTupler(GeoHelper& geohelp, art::TFileService& tfs, unsigned int nscMax);

  // Add the sim channels for an event.
  void fill(const art::Event& evt, const SimChannelVector& scs); 

private:

  // The n-tuples we'll create.
  TTree* m_ptree;

  // The variables that will go into the n-tuple.

  // Tree data.
  int m_event;
  int m_run;
  int m_subrun;
  unsigned int m_nchan;                 // # channels with SimChannel info
  std::vector<unsigned int> m_chan;     // Readout channel number
  std::vector<float> m_energy;          // Energy summed over all TDCs
  std::vector<float> m_charge;          // Charge summed over all TDCs
  std::vector<float> m_meantick;        // Energy-weighted mean of the ticks for the channel.
  std::vector<int> m_peaktick;          // Tick sample with largest energy for the channel.
  std::vector<float> m_peakenergy;      // Energy of TDC sample with largest energy
  std::vector<float> m_peakcharge;      // Charge of TDC sample with largest energy
  std::vector<unsigned int> m_peaknide; // # energy deposits contributing to the TDC.

  // Geometry service.
  GeoHelper& m_geohelp;

  // Maximum size for tree arrays.
  unsigned int m_nscMax;

};

#endif
