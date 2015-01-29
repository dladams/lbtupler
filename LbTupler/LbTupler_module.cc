// LbTupler_module.cc
// A basic "skeleton" to read in art::Event records from a file,
// access their information, and do something with them. 

// See
// <https://cdcvs.fnal.gov/redmine/projects/larsoftsvn/wiki/Using_the_Framework>
// for a description of the ART classes used here.

// Almost everything you see in the code below may have to be changed
// by you to suit your task. The example task is to make histograms
// and n-tuples related to dE/dx of particle tracks in the detector.

// As you try to understand why things are done a certain way in this
// example ("What's all this stuff about 'auto const&'?"), it will help
// to read ADDITIONAL_NOTES.txt in the same directory as this file.

#ifndef LbTupler_Module
#define LbTupler_Module

#include <iostream>
#include <iomanip>

// LArSoft includes
#include "Simulation/SimChannel.h"
#include "Simulation/LArG4Parameters.h"
#include "RecoBase/Hit.h"
#include "RecoBase/Cluster.h"
#include "Geometry/Geometry.h"
#include "SimulationBase/MCParticle.h"
#include "SimulationBase/MCTruth.h"
#include "SimpleTypesAndConstants/geo_types.h"

// Framework includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/FindManyP.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/ParameterSet.h"
#include "cetlib/exception.h"

// ROOT includes. Note: To look up the properties of the ROOT classes,
// use the ROOT web site; e.g.,
// <http://root.cern.ch/root/html532/ClassIndex.html>
#include "TH1.h"
#include "TH2.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include "TVector3.h"

// C++ Includes
#include <map>
#include <vector>
#include <algorithm>
#include <iostream>
#include <string>
#include <cmath>

using std::cout;
using std::endl;
using std::setw;

namespace LbTupler {

  //-----------------------------------------------------------------------
  //-----------------------------------------------------------------------
  // class definition

  class LbTupler : public art::EDAnalyzer 
  {
  public:
 
    // Standard constructor and destructor for an ART module.
    explicit LbTupler(fhicl::ParameterSet const& pset);
    virtual ~LbTupler();

    // This method is called once, at the start of the job. In this
    // example, it will define the histograms and n-tuples we'll write.
    void beginJob();

    // This method is called once, at the start of each run. It's a
    // good place to read databases or files that may have
    // run-dependent information.
    void beginRun(const art::Run& run);

    // This method reads in any parameters from the .fcl files. This
    // method is called 'reconfigure' because it might be called in the
    // middle of a job; e.g., if the user changes parameter values in an
    // interactive event display.
    void reconfigure(fhicl::ParameterSet const& pset);

    // The analysis routine, called once per event. 
    void analyze (const art::Event& evt); 

  private:

    // The stuff below is the part you'll most likely have to change to
    // go from this custom example to your own task.

    // The parameters we'll read from the .fcl file.
    int fdbg;                             // Debug level. Larger for more log noise.
    std::string fSimulationProducerLabel; // The name of the producer that tracked simulated particles through the detector
    std::string fHitProducerLabel;        // The name of the producer that created hits
    std::string fClusterProducerLabel;    // The name of the producer that created clusters
    std::string fRawDigitLabel;           // The name of the producer that created the raw digits.
    int fSelectedPDG;                     // PDG code of particle we'll focus on
    double fBinSize;                      // For dE/dx work: the value of dx. 

    // Pointers to the histograms we'll create. 
    TH1D* fPDGCodeHist;
    TH1D* fMomentumHist;
    TH1D* fTrackLengthHist;

    // The n-tuples we'll create.
    TTree* fSimulationNtuple;
    TTree* fReconstructionNtuple;
    TTree* fSimChannelNtuple;

    // The variables that will go into the n-tuple.
    int fEvent;
    int fRun;
    int fSubRun;
    int fPDG;
    int ftrackid;
    // Arrays for 4-vectors: (x,y,z,t) and (Px,Py,Pz,E).
    // Note: old-style C++ arrays are considered obsolete. However,
    // to create simple n-tuples, we still need to use them. 
    float fStartXYZT[4];
    float fEndXYZT[4];
    float fStartPE[4];
    float fEndPE[4];
    static const unsigned int maxpt = 20000;
    unsigned int fnpt;
    float fptx[maxpt];
    float fpty[maxpt];
    float fptz[maxpt];
    float fptt[maxpt];
    float fpte[maxpt];
    // Number of dE/dx bins in a given track. 
    int fndedxBins;
    // The vector that will be used to accumulate dE/dx values.
    std::vector<double> fdedxBins;

    // Other variables that will be shared between different methods.
    double                            fElectronsToGeV; // conversion factor

    // The maximum size of fdedxBins; in other words, it's the
    // capacity of the vector. If we ever perform an operation that
    // causes it to exceed this limit, it would mean fdedxBins would
    // shift in memory, the n-tuple branch would point to the wrong
    // location in memory, and everything would go wonky.
    unsigned int fMaxCapacity;

    // Sim channel info.
    unsigned int fscCapacity;
    unsigned int fscCount;
    std::vector<unsigned int> fscChannel;  // Readout channel number
    std::vector<float> fscEnergy;          // Energy summed over all TDCs
    std::vector<float> fscCharge;          // Charge summed over all TDCs
    std::vector<int> fscTdcPeak;           // TDC sample with largest energy for the channel.
    std::vector<float> fscTdcMean;         // Energy-weighted mean of the TDC samples for the channel.
    std::vector<float> fscTdcEnergy;       // Energy of TDC sample with largest energy
    std::vector<float> fscTdcCharge;       // Charge of TDC sample with largest energy
    std::vector<unsigned int> fscTdcNide;  // # energy deposits contributing to the TDC.

    // Geometry service.
    art::ServiceHandle<geo::Geometry> fGeometry;       // pointer to Geometry service

  }; // class LbTupler


  //-----------------------------------------------------------------------
  //-----------------------------------------------------------------------
  // class implementation

  //-----------------------------------------------------------------------
  // Constructor
  LbTupler::LbTupler(fhicl::ParameterSet const& parameterSet)
    : EDAnalyzer(parameterSet), fdbg(0) {
    // Read in the parameters from the .fcl file.
    this->reconfigure(parameterSet);
  }

  //-----------------------------------------------------------------------
  // Destructor
  LbTupler::~LbTupler() 
  {}
   
  //-----------------------------------------------------------------------
  void LbTupler::beginJob() {
    const std::string myname = "LbTupler::beginJob: ";

    // Access ART's TFileService, which will handle creating and writing
    // histograms and n-tuples for us. 
    art::ServiceHandle<art::TFileService> tfs;

    // Fetch geometry.
    int ncryo = fGeometry->Ncryostats();
    int ntpc = fGeometry->NTPC();
    int ntpcchan = fGeometry->Nchannels();
    std::vector<int> ntpcplane;
    std::vector<int> ntpcplanewire;
    std::vector<int> firsttpcplanewire;
    firsttpcplanewire.push_back(0);
    int icry = 0;
    int itpcplane = 0;
    // Loop over TPCs.
    for ( int itpc =0; itpc<ntpc; ++itpc ) {
      int nplane = fGeometry->Nplanes(itpc, icry);
      ntpcplane.push_back(nplane);
      // Loop over planes in the TPC.
      for ( int ipla=0; ipla<nplane; ++ipla ) {
        int nwire = fGeometry->Nwires(ipla, itpc, icry);
        ntpcplanewire.push_back(nwire);
        int firstwire = firsttpcplanewire[itpcplane];
        int lastwire = firstwire + nwire - 1;
        if ( itpc<ntpc-1 || ipla<nplane-1 ) {
          firsttpcplanewire.push_back(lastwire + 1);
        }
        // Loop over wires and find the channels.
        std::set<int> chans;
        for ( int iwir=firstwire; iwir<=lastwire; ++iwir ) {
          int icha = fGeometry->PlaneWireToChannel(ipla, iwir, itpc, icry);
          chans.insert(icha);
        }
        int nchan = chans.size();
        int firstchan = -1;
        int lastchan = -1;
        if ( nchan ) {
          firstchan = *chans.cbegin();
          lastchan = *chans.crbegin();
        }
        if ( fdbg > 0 ) {
          cout << myname << "TPC " << itpc << ", plane " << ipla <<  " has " << nwire
               << " wires: [" << setw(4) << firstwire << "," << setw(4) << lastwire << "]"
               << " and " << nchan << "/" << lastchan - firstchan + 1 << " channels: ["
               << setw(4) << firstchan << "," << setw(4) << lastchan << "]" << endl;
        }
        ++itpcplane;
      }
    }
    if ( fdbg > 0 ) {
      cout << myname << "        Total # cryostats: " << ncryo << endl;
      cout << myname << "             Total # TPCs: " << ntpc << endl;
      cout << myname << "     Total # TPC channels: " << ntpcchan << endl;
      cout << myname << "Total # optical detectors: " << fGeometry->NOpDet(icry) << endl;
      cout << myname << " Total # optical channels: " << fGeometry->NOpChannels() << endl;
    }
    double detLength = fGeometry->DetLength(); 
    double detWidth  = fGeometry->DetHalfWidth()  * 2.;
    double detHeight = fGeometry->DetHalfHeight() * 2.;
    if ( fdbg > 0 ) {
      cout << myname << "Detector length: " << detLength << endl;
      cout << myname << "Detector width:  " << detWidth  << endl;
      cout << myname << "Detector height: " << detHeight << endl;
    }
    double detSize = std::sqrt( detLength*detLength + detWidth*detWidth + detHeight*detHeight );

    // Set up the vector that will be used to accumulate the dE/dx
    // values into bins. We want the capacity of this vector to be big
    // enough to hold the bins of the longest possible track.  The
    // following code will work if the detector has just one TPC and
    // all tracks are contained within it; if there's a possibility
    // that a single track could go through or outside multiple TPCs
    // then you'll have to think some more.
    // The reason we're doing this is that it's important that the
    // memory location of fdEdxBins.data() not change throughout the
    // execution of the program; see the TTree::Branch() calls below.
    fMaxCapacity = 1;
    if ( fBinSize > 0.0 ) fMaxCapacity = (unsigned int)(detSize/fBinSize ) + 1;
    fdedxBins.reserve(fMaxCapacity);
    fscChannel.reserve(fscCapacity);
    fscEnergy.reserve(fscCapacity);
    fscCharge.reserve(fscCapacity);
    fscTdcPeak.reserve(fscCapacity);
    fscTdcMean.reserve(fscCapacity);
    fscTdcEnergy.reserve(fscCapacity);
    fscTdcCharge.reserve(fscCapacity);
    fscTdcNide.reserve(fscCapacity);

    // Define the histograms. Putting semi-colons around the title
    // causes it to be displayed as the x-axis label if the histogram
    // is drawn.
    fPDGCodeHist     = tfs->make<TH1D>("pdgcodes",";PDG Code;",                  5000, -2500, 2500);
    fMomentumHist    = tfs->make<TH1D>("mom",     ";particle Momentum (GeV);",    100, 0.,    10.);
    fTrackLengthHist = tfs->make<TH1D>("length",  ";particle track length (cm);", 200, 0, detLength);

    // Define our n-tuples, which are limited forms of ROOT
    // TTrees. Start with the TTree itself.
    fSimulationNtuple     = tfs->make<TTree>("LbTuplerSimulation",    "LbTuplerSimulation");
    fReconstructionNtuple = tfs->make<TTree>("LbTuplerReconstruction","LbTuplerReconstruction");
    fSimChannelNtuple = tfs->make<TTree>("LbTuplerSimChannel","LbTuplerSimChannel");

    // Define the branches (columns) of our simulation n-tuple. When
    // we write a variable, we give the address of the variable to
    // TTree::Branch.
    fSimulationNtuple->Branch("event",       &fEvent,          "event/I");
    fSimulationNtuple->Branch("subrun",      &fSubRun,         "subrun/I");
    fSimulationNtuple->Branch("run",         &fRun,            "run/I");
    fSimulationNtuple->Branch("trackid",     &ftrackid,        "trackid/I");
    fSimulationNtuple->Branch("pdg",         &fPDG,            "pdg/I");
    // When we write arrays, we give the address of the array to
    // TTree::Branch; in C++ this is simply the array name.
    fSimulationNtuple->Branch("StartXYZT",   fStartXYZT,       "StartXYZT[4]/F");
    fSimulationNtuple->Branch("EndXYZT",     fEndXYZT,         "EndXYZT[4]/F");
    fSimulationNtuple->Branch("StartPE",     fStartPE,         "StartPE[4]/F");
    fSimulationNtuple->Branch("EndPE",       fEndPE,           "EndPE[4]/F");
    // Trajectory points.
    fSimulationNtuple->Branch("npt",       &fnpt,          "npt/i");
    fSimulationNtuple->Branch("ptx",       fptx,           "ptx[npt]/F");
    fSimulationNtuple->Branch("pty",       fpty,           "pty[npt]/F");
    fSimulationNtuple->Branch("ptz",       fptz,           "ptz[npt]/F");
    fSimulationNtuple->Branch("ptt",       fptt,           "ptt[npt]/F");
    fSimulationNtuple->Branch("pte",       fpte,           "pte[npt]/F");
    // For a variable-length array: include the number of bins.
    fSimulationNtuple->Branch("ndedx",       &fndedxBins,      "ndedx/I");
    // We're using a memory trick here: the data() method returns the
    // address of the array inside the vector. Note after we call this
    // method, the address of fdedxBins must not change.
    fSimulationNtuple->Branch("dedx",        fdedxBins.data(), "dedx[ndedx]/D");

    // A similar definition for the reconstruction n-tuple. Note that we
    // use some of the same variables in both n-tuples.
    fReconstructionNtuple->Branch("Event",   &fEvent,          "Event/I");
    fReconstructionNtuple->Branch("SubRun",  &fSubRun,         "SubRun/I");
    fReconstructionNtuple->Branch("Run",     &fRun,            "Run/I");
    fReconstructionNtuple->Branch("TrackID", &ftrackid,        "TrackID/I");
    fReconstructionNtuple->Branch("PDG",     &fPDG,            "PDG/I");
    fReconstructionNtuple->Branch("ndedx",   &fndedxBins,      "ndedx/I");
    fReconstructionNtuple->Branch("dedx",    fdedxBins.data(), "dedx[ndedx]/D");

    // Sim channel tree.
    fSimChannelNtuple->Branch("event",   &fEvent,          "event/I");
    fSimChannelNtuple->Branch("subrun",  &fSubRun,         "subrun/I");
    fSimChannelNtuple->Branch("run",     &fRun,            "run/I");
    fSimChannelNtuple->Branch("nchan", &fscCount, "nchan/i");
    fSimChannelNtuple->Branch("chan", fscChannel.data(), "channnel[nchan]/i");
    fSimChannelNtuple->Branch("energy", fscEnergy.data(), "energy[nchan]/F");
    fSimChannelNtuple->Branch("charge", fscCharge.data(), "charge[nchan]/F");
    fSimChannelNtuple->Branch("tdcpeak", fscTdcPeak.data(), "tdcpeak[nchan]/I");
    fSimChannelNtuple->Branch("tdcmean", fscTdcMean.data(), "tdcmean[nchan]/F");
    fSimChannelNtuple->Branch("tdcenergy", fscTdcEnergy.data(), "tdcenergy[nchan]/F");
    fSimChannelNtuple->Branch("tdccharge", fscTdcCharge.data(), "tdccharge[nchan]/F");
    fSimChannelNtuple->Branch("tdcnide",   fscTdcNide.data(),   "tdcnide[nchan]/i");
  }
   
  //-----------------------------------------------------------------------
  void LbTupler::beginRun(const art::Run& /*run*/)
  {
    // How to convert from number of electrons to GeV.  The ultimate
    // source of this conversion factor is
    // ${LARSIM_DIR}/include/SimpleTypesAndConstants/PhysicalConstants.h.
    fElectronsToGeV = 1./4.2e7;
  }

  //-----------------------------------------------------------------------
  void LbTupler::reconfigure(fhicl::ParameterSet const& p)
  {
    // Read parameters from the .fcl file. The names in the arguments
    // to p.get<TYPE> must match names in the .fcl file.
    fdbg                     = p.get<int>        ("DebugLevel");
    fSimulationProducerLabel = p.get<std::string>("SimulationLabel");
    fRawDigitLabel           = p.get<std::string>("RawDigitLabel");
    fHitProducerLabel        = p.get<std::string>("HitLabel");
    fClusterProducerLabel    = p.get<std::string>("ClusterLabel");
    fSelectedPDG             = p.get<int        >("PDGcode");
    fBinSize                 = p.get<double     >("BinSize");
    fscCapacity              = p.get<double     >("SimChannelSize");
    return;
  }

  //-----------------------------------------------------------------------
  void LbTupler::analyze(const art::Event& event) 
  {
    int dbg = fdbg;
    const std::string myname = "LbTupler::analyze: ";
    // Start by fetching some basic event information for our n-tuple.
    fEvent  = event.id().event(); 
    fRun    = event.run();
    fSubRun = event.subRun();
    if ( dbg > 0 ) cout << myname << "Processing run " << fRun << "-" << fSubRun
         << ", event " << fEvent << endl;

    // This is the standard method of reading multiple objects
    // associated with the same event; see
    // <https://cdcvs.fnal.gov/redmine/projects/larsoftsvn/wiki/Using_the_Framework>
    // for more information. Define a "handle" to point to a vector of
    // the objects.
    art::Handle< std::vector<simb::MCParticle> > particleHandle;
    // Then tell the event to fill the vector with all the objects of
    // that type produced by a particular producer.
    event.getByLabel(fSimulationProducerLabel, particleHandle);
    if ( dbg > 1 ) cout << myname << "Particle count: " << particleHandle->size() << endl;

    // Get all the simulated channels for the event. These channels
    // include the energy deposited for each track.
    art::Handle< std::vector<sim::SimChannel> > simChannelHandle;
    event.getByLabel(fSimulationProducerLabel, simChannelHandle);
    if ( dbg > 1 ) cout << myname << "Sim channel count: " << simChannelHandle->size() << endl;
    if ( simChannelHandle->size() > fscCapacity ) {
      cout << myname << "WARNING: Sim channel count exceeds TTree capacity." << endl;
    }

    // Get the raw digits for the event.
    art::Handle< std::vector<raw::RawDigit> > rawDigitHandle;
    event.getByLabel(fRawDigitLabel, rawDigitHandle);
    if ( dbg > 1 ) cout << myname << "Raw digit count: " << rawDigitHandle->size() << endl;

    // Loop over particles.
    for ( auto const& particle : (*particleHandle) ) {
      // For the methods you can call to get particle information,
      // see ${NUTOOLS_DIR}/include/SimulationBase/MCParticle.h.
      ftrackid = particle.TrackId();

      // Add the address of the MCParticle to the map, with the track ID as the key.
      //dla particleMap[ftrackid] = &particle;

      // Histogram the PDG code of every particle in the event.
      fPDG = particle.PdgCode();
      fPDGCodeHist->Fill( fPDG );


      // For this example, we want to fill the n-tuples and histograms
      // only with information from the primary particles in the
      // event, whose PDG codes match a value supplied in the .fcl file.
      if ( true || (particle.Process() == "primary"  &&  fPDG == fSelectedPDG) ) {
        // A particle has a trajectory, consisting of a set of
        // 4-positions and 4-mommenta.
        size_t numberTrajectoryPoints = particle.NumberTrajectoryPoints();

        // For trajectories, as for vectors and arrays, the
        // first point is #0, not #1.
        int last = numberTrajectoryPoints - 1;
        const TLorentzVector& positionStart = particle.Position(0);
        const TLorentzVector& positionEnd   = particle.Position(last);
        const TLorentzVector& momentumStart = particle.Momentum(0);
        const TLorentzVector& momentumEnd   = particle.Momentum(last);

        // Make a histogram of the starting momentum.
        fMomentumHist->Fill( momentumStart.P() );

        // Fill arrays with the 4-values. (Don't be fooled by
        // the name of the method; it just puts the numbers from
        // the 4-vector into the array.)
        positionStart.GetXYZT( fStartXYZT );
        positionEnd.GetXYZT( fEndXYZT );
        momentumStart.GetXYZT( fStartPE );
        momentumEnd.GetXYZT( fEndPE );

        // Fill trajectory.
        fnpt = 0.0;
        for ( unsigned int ipt=0; ipt<numberTrajectoryPoints; ++ipt ) {
          if ( ipt >= maxpt ) {
            cout << myname << "Found more than " << maxpt << " trajectory points."
                 << " The remainder will be skipped." << endl;
            break;
          }
          const auto& pos = particle.Position(ipt);
          const auto& mom = particle.Momentum(ipt);
          fptx[ipt] = pos.X();
          fpty[ipt] = pos.Y();
          fptz[ipt] = pos.Z();
          fptt[ipt] = pos.T();
          fpte[ipt] = pos.E();
          ++fnpt;
        }

        // Use a polar-coordinate view of the 4-vectors to
        // get the track length.
        double trackLength = ( positionEnd - positionStart ).Rho();
        if ( dbg > 2 ) cout << myname << "Track process: " << particle.Process() << endl;
        if ( dbg > 2 ) cout << myname << " Track length: " << trackLength << " mm" << endl;

        // Make a histogram of the track length.
        fTrackLengthHist->Fill( trackLength );

        // Determine the number of dE/dx bins for the n-tuple.
        fndedxBins = 1;
        if ( fBinSize > 0.0 ) fndedxBins = int( trackLength / fBinSize ) + 1;
        // Initialize the vector of dE/dx bins to zero.
        for ( auto& val : fdedxBins ) val = 0.0;

        // To look at the energy deposited by this particle's track,
        // we loop over the SimChannel objects in the event.
        for ( auto const& channel : (*simChannelHandle) ) {
          // Get the numeric ID associated with this channel.
          // (The channel number is a 32-bit unsigned int, which
          // normally requires a special data type. Let's use
          // "auto" so we don't have to remember that.)
          //dla auto channelNumber = channel.Channel();

          // A little care: There is more than one plane that
          // accumulates charge in the TPC. We only want to
          // include the energy from the collection plane.
          // (geo::kCollection is defined in
          // ${LARSIM_DIR}/include/SimpleTypesAndConstants/geo_types.h)
          if ( true ) {
          //dla if ( fGeometry->SignalType( channelNumber ) == geo::kCollection ) {
            // Each channel has a map inside it that connects
            // a time slice to energy deposits in the
            // detector.  We'll use "auto", but it's worth
            // noting that the full type of this map is
            // std::map<unsigned short, std::vector<sim::IDE>>;
            auto const& timeSlices = channel.TDCIDEMap();

            // For every time slice in this channel:
            for ( auto const& timeSlice : timeSlices ) {
              // Each entry in a map is a pair<first,second>.
              // For the timeSlices map, the 'first' is a time
              // slice number, which we don't care about in this
              // example. The 'second' is a vector of IDE
              // objects.
              auto const& energyDeposits = timeSlice.second;
              // Loop over the energy deposits. The type of
              // 'energyDeposit' will be sim::IDE, which is
              // defined in ${LARSIM_DIR}/include/Simulation/SimChannel.h.
              for ( auto const& energyDeposit : energyDeposits ) {
                // Check if the track that deposited the
                // energy matches the track of the particle.
                if ( energyDeposit.trackID == ftrackid ) {
                  // Get the (x,y,z) of the energy deposit.
                  TVector3 location( energyDeposit.x,
                                     energyDeposit.y,
                                     energyDeposit.z );

                  // Distance from the start of the track.
                  double distance = ( location - positionStart.Vect() ).Mag();

                  // Into which bin of the dE/dx array do we add the energy?
                  unsigned int bin = 0;
                  if ( fBinSize > 0.0 ) bin = (unsigned int)( distance / fBinSize );

                  // Is the dE/dx array big enough to include this bin?
                  // And do we have the capacity for it?
                  if ( fdedxBins.size() < bin+1 ) {
                    if ( bin+1 > fMaxCapacity ) {
                      throw cet::exception("LbTupler")
                        << " Exceeded capacity of dedx vector; "
                        << " you need to revise maximum track length calculation\n";
                    }
                    //  Increase the array size, padding it with zeros. 
                    fdedxBins.resize( bin+1 , 0 );
                  }

                  // Add the energy deposit to that bin. (If you look at the
                  // definition of sim::IDE, you'll see that there's another
                  // way to get the energy. Are the two methods equivalent?
                  // Compare the results and see!)
                  double nel = energyDeposit.numElectrons;
                  double de = nel * fElectronsToGeV;
                  fdedxBins[bin] += de;
                  if ( dbg > 4 ) cout << myname << "  Nel, dE: " << nel << ", " << de << endl;
                }
              } // For each energy deposit
            } // For each time slice
          } // If channel is in a collection plane
        } // For each SimChannel

        // At this point we've filled in the values of all the
        // variables and arrays we want to write to the
        // n-tuple. The following command actually writes those
        // values.
        fSimulationNtuple->Fill();

      } // if selected
    } // end loop over all particles in the event. 

    // Loop over the SimChannel objects in the event.
    fscCount = 0;
    for ( auto& val : fscChannel ) val = 0;
    for ( auto& val : fscCharge ) val = 0.0;
    for ( auto& val : fscEnergy ) val = 0.0;
    for ( auto const& chan : (*simChannelHandle) ) {
      auto ichan = chan.Channel();
      if ( dbg > 3 ) cout << myname << " " << fscCount << ": Sim channel: " << ichan << endl;
      if ( fscCount < fscCapacity ) {
        fscChannel[fscCount] = ichan;
        auto const& idemap = chan.TDCIDEMap();
        double energy = 0.0;  // Total energy in the channel
        double charge = 0.0;  // Total charge in the channel
        double energytdc = 0.0;  // Total energy*tdc in the channel
        // Loop over TDC samples.
        short tdcPeak = 0;     // TDC sample with the largest energy for this channel
        double tdcEnergyMax = 0.0;  // Energy for that TDC
        double tdcChargeMax = 0.0;  // Charge for that TDC
        unsigned int tdcNideMax = 0.0;  // Charge for that TDC
        for ( auto const& ideent : idemap ) {
          short tdc = ideent.first;
          auto idevec = ideent.second;
          double tdcEnergy = 0.0;
          double tdcCharge = 0.0;
          unsigned int tdcNide = idevec.size();
          // Sum over the tracks contributing to this channel sample.
          for ( auto& ide : idevec ) {
            tdcEnergy += ide.energy;
            tdcCharge += ide.numElectrons;
          }
          if ( tdcPeak == 0 || tdcEnergy > tdcEnergyMax ) {
            tdcPeak = tdc;
            tdcEnergyMax = tdcEnergy;
            tdcChargeMax = tdcCharge;
            tdcNideMax = tdcNide;
          }
          energy += tdcEnergy;
          charge += tdcCharge;
          energytdc += tdcEnergy*tdc;
        }
        if ( dbg > 3 ) cout << myname << "  Deposited energy, charge: "
          << energy << ", " << charge << endl;
        fscEnergy[fscCount] = energy;
        fscCharge[fscCount] = charge;
        fscTdcPeak[fscCount] = tdcPeak;
        fscTdcMean[fscCount] = energytdc/energy;
        fscTdcEnergy[fscCount] = tdcEnergyMax;
        fscTdcCharge[fscCount] = tdcChargeMax;
        fscTdcNide[fscCount] = tdcNideMax;
        ++fscCount;
      } else {
        if ( dbg > 3 ) cout << myname << "Sim channel: " << ichan << endl;
      }
    } // end loop over sim channels in the event. 
    fSimChannelNtuple->Fill();

    return;
  }

  // This macro has to be defined for this module to be invoked from a
  // .fcl file; see LbTupler.fcl for more information.
  DEFINE_ART_MODULE(LbTupler)

} // namespace LbTupler

#endif // LbTupler_Module
