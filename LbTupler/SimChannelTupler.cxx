// SimChannelTupler.cxx

// David Adams
// May 2015

#include "SimChannelTupler.h"
#include <iostream>
#include <iomanip>
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "Simulation/SimChannel.h"
#include "TTree.h"
#include "GeoHelper.h"

using std::string;
using std::cout;
using std::endl;
using std::setw;

//************************************************************************

SimChannelTupler::
SimChannelTupler(GeoHelper& geohelp, art::TFileService& tfs, unsigned int nscMax)
: m_geohelp(geohelp), m_nscMax(nscMax) {

  const string myname = "SimChannelTupler::ctor: ";

  m_chan.reserve(nscMax);
  m_energy.reserve(nscMax);
  m_charge.reserve(nscMax);
  m_peaktick.reserve(nscMax);
  m_meantick.reserve(nscMax);
  m_peakenergy.reserve(nscMax);
  m_peakcharge.reserve(nscMax);
  m_peaknide.reserve(nscMax);

  // Sim channel tree.
  m_ptree = tfs.make<TTree>("LbTuplerSimChannel", "LbTuplerSimChannel");
  m_ptree->Branch("event",     &m_event,             "event/I");
  m_ptree->Branch("subrun",    &m_subrun,            "subrun/I");
  m_ptree->Branch("run",       &m_run,               "run/I");
  m_ptree->Branch("nchan",     &m_nchan,             "nchan/i");
  m_ptree->Branch("chan",       m_chan.data(),       "chan[nchan]/i");        // Channel
  m_ptree->Branch("energy",     m_energy.data(),     "energy[nchan]/F");      // Total energy in channel
  m_ptree->Branch("charge",     m_charge.data(),     "charge[nchan]/F");      // Total charge in channel
  m_ptree->Branch("meantick",   m_meantick.data(),   "meantick[nchan]/F");    // Mean tick
  m_ptree->Branch("peaktick",   m_peaktick.data(),   "peaktick[nchan]/I");    // Tick with the peak charge
  m_ptree->Branch("peakenergy", m_peakenergy.data(), "peakenergy[nchan]/F");  // Energy in peak tick
  m_ptree->Branch("peakcharge", m_peakcharge.data(), "peakcharge[nchan]/F");  // Chargy in peak tick
  m_ptree->Branch("peaknide",   m_peaknide.data(),   "peaknide[nchan]/i");    // # IDE in peak tick
}
 
//************************************************************************

void SimChannelTupler::fill(const art::Event& evt, const SimChannelVector& scs) {

  const string myname = "SimChannelTupler::fill: ";
  int dbg = 0;

  // Record event info.
  m_event  = evt.id().event();
  m_run    = evt.run();
  m_subrun = evt.subRun();

  // Loop over the SimChannel objects in the event.
  if ( dbg > 0 ) cout << myname << "Looping over sim channels (size = " << scs.size() << ")"
                         << " to fill SimChannel tree." << endl;
  m_nchan = 0;
  for ( auto& val : m_chan ) val = 0;
  for ( auto& val : m_charge ) val = 0.0;
  for ( auto& val : m_energy ) val = 0.0;
  for ( auto const& simchan : scs ) {
    auto ichan = simchan.Channel();
    unsigned int irop = m_geohelp.channelRop(ichan);
    if ( irop == m_geohelp.nrop() ) {
      cout << myname << "ERROR: SimChannel channel " << ichan << " is not in a readout plane." << endl;
      abort();
    }
    int iropchan = ichan - int(m_geohelp.ropFirstChannel(irop));
    if ( iropchan < 0 || iropchan >= int(m_geohelp.ropNChannel(irop)) ) {
      cout << myname << "ERROR: ROP channel " << iropchan << " is out of range [0, "
           << m_geohelp.ropNChannel(irop) << ")" << endl;
      abort();
    }
    if ( dbg > 2 ) cout << myname << "SimChannel " << setw(4) << m_nchan
                        << ": ROP " << setw(2) << irop
                        << ", Global/ROP channel " << setw(4) << ichan
                        << "/" << setw(3) << iropchan;
    if ( m_nchan < m_nscMax  ) {
      m_chan[m_nchan] = ichan;
      auto const& idemap = simchan.TDCIDEMap();
      double energy = 0.0;  // Total energy in the channel
      double charge = 0.0;  // Total charge in the channel
      double energytick = 0.0;  // Total energy*tick in the channel
      // Loop over TDC samples.
      short tickPeak = 0;              // Tick with the largest energy for this channel
      double tickEnergyMax = 0.0;      // Energy for that tick
      double tickChargeMax = 0.0;      // Charge for that tick
      unsigned int tickNideMax = 0.0;  // # IDE for that tick
      for ( auto const& ideent : idemap ) {
        short tick = ideent.first;
        auto idevec = ideent.second;
        double tickEnergy = 0.0;
        double tickCharge = 0.0;
        unsigned int tickNide = idevec.size();
        // Sum over the tracks contributing to this channel sample.
        for ( auto& ide : idevec ) {
          tickEnergy += ide.energy;
          tickCharge += ide.numElectrons;
        }
        // Update the peak channel data.
        if ( tickPeak == 0 || tickEnergy > tickEnergyMax ) {
          tickPeak = tick;
          tickEnergyMax = tickEnergy;
          tickChargeMax = tickCharge;
          tickNideMax = tickNide;
        }
        energy += tickEnergy;
        charge += tickCharge;
        energytick += tick*tickEnergy;
      }
      if ( dbg > 2 ) {
        double ren = 0.01*int(100.0*energy+0.49999);
        double rchg = int(charge+0.499999);
        cout << "E=" << ren << " MeV; Q=" << rchg << endl;
      }
      m_energy[m_nchan] = energy;
      m_charge[m_nchan] = charge;
      m_peaktick[m_nchan] = tickPeak;
      m_meantick[m_nchan] = energytick/energy;
      m_peakenergy[m_nchan] = tickEnergyMax;
      m_peakcharge[m_nchan] = tickChargeMax;
      m_peaknide[m_nchan] = tickNideMax;
      ++m_nchan;
    } else {
      cout << myname << "WARNING: Skipping sim channel " << ichan << endl;
    }
  } // end loop over sim channels in the event. 

  // Fill SimChannel tree.
  m_ptree->Fill();

}
