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
#include <sstream>

// LArSoft includes
#include "Simulation/SimChannel.h"
#include "Simulation/LArG4Parameters.h"
#include "RecoBase/Hit.h"
#include "RecoBase/Cluster.h"
#include "Geometry/Geometry.h"
//#include "HitFinderLBNE/APAGeometryAlg.h"
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
using std::ostringstream;
using std::setw;
using std::string;
using std::vector;
using std::set;
using std::sqrt;

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
    string fSimulationProducerLabel; // The name of the producer that tracked simulated particles through the detector
    string fHitProducerLabel;        // The name of the producer that created hits
    string fClusterProducerLabel;    // The name of the producer that created clusters
    string fRawDigitLabel;           // The name of the producer that created the raw digits.
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
    int fevent;
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
    vector<double> fdedxBins;

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
    unsigned int ftdcTickMin;         // First TDC bin to draw.
    unsigned int ftdcTickMax;         // Last+1 TDC bin to draw.
    vector<unsigned int> fscChannel;  // Readout channel number
    vector<float> fscEnergy;          // Energy summed over all TDCs
    vector<float> fscCharge;          // Charge summed over all TDCs
    vector<int> fscTdcPeak;           // TDC sample with largest energy for the channel.
    vector<float> fscTdcMean;         // Energy-weighted mean of the TDC samples for the channel.
    vector<float> fscTdcEnergy;       // Energy of TDC sample with largest energy
    vector<float> fscTdcCharge;       // Charge of TDC sample with largest energy
    vector<unsigned int> fscTdcNide;  // # energy deposits contributing to the TDC.

    // Geometry service.
    art::ServiceHandle<geo::Geometry> fGeometry;       // pointer to Geometry service

    // Geometry information.
    unsigned int fncryo;         // # cryostats
    unsigned int fntpc;          // # TPC
    unsigned int fntpp;          // Total # TPC planes
    vector<int> fntpcplane;      // # planes in each TPC
    vector<int> fntpcplanewire;  // # wires in each plane
    vector<string> ftpcname;     // Names for the TPCs.
    vector<string> fplanename;   // Names for the TPC planes.
    unsigned int fnrop;          // Total # readout planes (ROPs)
    vector<int> fropfirstchan;   // first channel for each readout plane
    vector<int> fropnchan;       // # channels for each readout plane
    vector<int> froptpc;         // First TPC for each readout plane
    vector<string> fropname;     // Names for the TPC readout planes, e.g. rop1x1.
    unsigned int fnapa;          // Total # APA
    vector<int> fapanrop;        // # ROP for each APA.


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
    const string myname = "LbTupler::beginJob: ";

    // Access ART's TFileService, which will handle creating and writing
    // histograms and n-tuples for us. 
    art::ServiceHandle<art::TFileService> tfs;

    // Fetch geometry.
    fncryo = fGeometry->Ncryostats();
    fntpc = fGeometry->NTPC();
    int ntpcchan = fGeometry->Nchannels();
    vector<int> firsttpcplanewire;
    firsttpcplanewire.push_back(0);
    int icry = 0;
    fntpp = 0;
    // Loop over TPCs.
    for ( unsigned int itpc =0; itpc<fntpc; ++itpc ) {
      int nplane = fGeometry->Nplanes(itpc, icry);
      ostringstream sstpc;
      sstpc << "TPC" << itpc;
      string tpcname = sstpc.str();
      fntpcplane.push_back(nplane);
      ftpcname.push_back(tpcname);
      int itdcrop = 0;   // # readouts for this TDC
      // Loop over planes in the TPC.
      for ( int ipla=0; ipla<nplane; ++ipla ) {
        int nwire = fGeometry->Nwires(ipla, itpc, icry);
        fntpcplanewire.push_back(nwire);
        ostringstream ssplane;
        ssplane << tpcname << "Plane" << ipla;
        string planename = ssplane.str();
        fplanename.push_back(planename);
        int firstwire = firsttpcplanewire[fntpp];
        int lastwire = firstwire + nwire - 1;
        if ( itpc<fntpc-1 || ipla<nplane-1 ) {
          firsttpcplanewire.push_back(lastwire + 1);
        }
        // Loop over wires and find the channels.
        set<int> chans;
        //for ( int iwir=firstwire; iwir<=lastwire; ++iwir ) {
        for ( int iwir=0; iwir<nwire; ++iwir ) {
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
        if ( lastchan <= firstchan ) {
          cout << myname << "ERROR: Invalid channel range." << endl;
          abort();
        }
        if ( fdbg > 0 ) {
          if ( fntpp == 0 ) cout << myname << "TPC wires:" << endl;
          cout << myname << planename <<  " has " << nwire
               << " wires: [" << setw(4) << firstwire << "," << setw(4) << lastwire << "]"
               << " and " << nchan << "/" << lastchan - firstchan + 1 << " channels: ["
               << setw(4) << firstchan << "," << setw(4) << lastchan << "]" << endl;
        }
        // Check if the range for the current readout plane (ROP) is already covered.
        unsigned int irop = fnrop;
        for ( irop=0; irop<fnrop; ++irop ) {
          int ropfirst = fropfirstchan[irop];
          int roplast = ropfirst + fropnchan[irop] - 1;
          if ( firstchan > roplast ) continue;               // New range is after the ROP
          if ( lastchan < ropfirst ) continue;               // New range is before the ROP
          if ( firstchan==ropfirst && lastchan==roplast ) break;  // Exact match
          // We could extend the range here but current geometry does not require this.
          cout << myname << "ERROR: Invalid channel range overlap." << endl;
          abort();
        }
        // If this is a new ROP, then record its channel range and TPC, and assign it to
        // an APA. For now, the latter assumes APA ordering follows TPC and is
        if ( irop == fnrop ) {  // We have a new readout plane.
          unsigned int iapa = itpc/2;
          for ( unsigned int japa=fnapa; japa<=iapa; ++japa ) {
            fapanrop.push_back(0);
            ++fnapa;
          }
          ++fapanrop[iapa];
          fropfirstchan.push_back(firstchan);
          fropnchan.push_back(lastchan - firstchan + 1 );
          froptpc.push_back(itpc);
          ostringstream ssrop;
          const vector<string> pnames = {"u", "v", "x1", "x2", "a", "b", "c", "d", "e"};
          ssrop << "apa" << iapa << pnames[fapanrop[iapa]-1];
          fropname.push_back(ssrop.str());
          ++fnrop;
          ++itdcrop;
        }
        ++fntpp;
      }
    }
    if ( fdbg > 0 ) {
      cout << myname << "        Total # cryostats: " << fncryo << endl;
      cout << myname << "             Total # TPCs: " << fntpc << endl;
      cout << myname << "     Total # TPC channels: " << ntpcchan << endl;
      cout << myname << "Total # optical detectors: " << fGeometry->NOpDet(icry) << endl;
      cout << myname << " Total # optical channels: " << fGeometry->NOpChannels() << endl;
      cout << myname << "There are " << fnrop << " ROPs (readout planes):" << endl;
      cout << myname << "      name     first chan     #chan" << endl;
      for ( unsigned int irop=0; irop<fnrop; ++irop ) {
        cout << myname << setw(10) << fropname[irop] << setw(10) << fropfirstchan[irop]
             << setw(15) << fropnchan[irop] << endl;
      }
    }
    double detLength = fGeometry->DetLength(); 
    double detWidth  = fGeometry->DetHalfWidth()  * 2.;
    double detHeight = fGeometry->DetHalfHeight() * 2.;
    if ( fdbg > 0 ) {
      cout << myname << "Detector length: " << detLength << " cm" << endl;
      cout << myname << "Detector width:  " << detWidth  << " cm" << endl;
      cout << myname << "Detector height: " << detHeight << " cm" << endl;
    }
    double detSize = sqrt( detLength*detLength + detWidth*detWidth + detHeight*detHeight );

    // Geometry dump from Michelle.
    if ( fdbg > 3 ) {
      double xyz[3];                                                                                                                                   
      double abc[3];                                                                 
      int chan;
      int cryo = fGeometry->Ncryostats();                                                             
      for (int c=0; c<cryo;++c){                   
        int tpc =fGeometry->NTPC(c);                                               
        for (int t=0; t<tpc; ++t){                                            
          int Nplanes=fGeometry->Nplanes(t,c);                                      
          for (int p=0;p<Nplanes;++p) {                                        
            int Nwires = fGeometry->Nwires(p,t,c);                                  
            cout << "FLAG " << endl;
            for (int w=0;w<Nwires;++w){
              fGeometry->WireEndPoints(c,t,p,w,xyz,abc);
              chan=fGeometry->PlaneWireToChannel(p,w,t,c);  
              cout << "FLAG " << setw(4) << chan << " " << c << " " << t << " " << p << setw(4) << w << " "
                        << xyz[0] << " " << xyz[1] << " " << xyz[2] <<  " "
                        << abc[0] << " " << abc[1] << " " << abc[2] << endl;                                         
            }
          }
        }
      }
    }

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
    fSimChannelNtuple     = tfs->make<TTree>("LbTuplerSimChannel","LbTuplerSimChannel");

    // Define the branches (columns) of our simulation n-tuple. When
    // we write a variable, we give the address of the variable to
    // TTree::Branch.
    fSimulationNtuple->Branch("event",       &fevent,          "event/I");
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
    fReconstructionNtuple->Branch("Event",   &fevent,          "Event/I");
    fReconstructionNtuple->Branch("SubRun",  &fSubRun,         "SubRun/I");
    fReconstructionNtuple->Branch("Run",     &fRun,            "Run/I");
    fReconstructionNtuple->Branch("TrackID", &ftrackid,        "TrackID/I");
    fReconstructionNtuple->Branch("PDG",     &fPDG,            "PDG/I");
    fReconstructionNtuple->Branch("ndedx",   &fndedxBins,      "ndedx/I");
    fReconstructionNtuple->Branch("dedx",    fdedxBins.data(), "dedx[ndedx]/D");

    // Sim channel tree.
    fSimChannelNtuple->Branch("event",   &fevent,          "event/I");
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
    fSimulationProducerLabel = p.get<string>("SimulationLabel");
    fRawDigitLabel           = p.get<string>("RawDigitLabel");
    fHitProducerLabel        = p.get<string>("HitLabel");
    fClusterProducerLabel    = p.get<string>("ClusterLabel");
    fSelectedPDG             = p.get<int>("PDGcode");
    fBinSize                 = p.get<double>("BinSize");
    fscCapacity              = p.get<double>("SimChannelSize");
    ftdcTickMin              = p.get<int>("TdcTickMin");
    ftdcTickMax              = p.get<int>("TdcTickMax");
    return;
  }

  //-----------------------------------------------------------------------
  void LbTupler::analyze(const art::Event& event) 
  {
    // Access ART's TFileService, which will handle creating and writing
    // histograms and n-tuples for us. 
    art::ServiceHandle<art::TFileService> tfs;

    int dbg = fdbg;
    const string myname = "LbTupler::analyze: ";
    // Start by fetching some basic event information for our n-tuple.
    fevent  = event.id().event(); 
    fRun    = event.run();
    fSubRun = event.subRun();
    if ( dbg > 0 ) cout << myname << "Processing run " << fRun << "-" << fSubRun
         << ", event " << fevent << endl;

    // This is the standard method of reading multiple objects
    // associated with the same event; see
    // <https://cdcvs.fnal.gov/redmine/projects/larsoftsvn/wiki/Using_the_Framework>
    // for more information. Define a "handle" to point to a vector of
    // the objects.
    art::Handle< vector<simb::MCParticle> > particleHandle;
    // Then tell the event to fill the vector with all the objects of
    // that type produced by a particular producer.
    event.getByLabel(fSimulationProducerLabel, particleHandle);
    if ( dbg > 1 ) cout << myname << "Particle count: " << particleHandle->size() << endl;

    // Get all the simulated channels for the event. These channels
    // include the energy deposited for each track.
    art::Handle< vector<sim::SimChannel> > simChannelHandle;
    event.getByLabel(fSimulationProducerLabel, simChannelHandle);
    if ( dbg > 1 ) cout << myname << "Sim channel count: " << simChannelHandle->size() << endl;
    if ( simChannelHandle->size() > fscCapacity ) {
      cout << myname << "WARNING: Sim channel count exceeds TTree capacity." << endl;
    }

    // Get the raw digits for the event.
    art::Handle< vector<raw::RawDigit> > rawDigitHandle;
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
          // DLA: I include all planes.
          if ( true ) {
          //dla if ( fGeometry->SignalType( channelNumber ) == geo::kCollection ) {
            // Each channel has a map inside it that connects
            // a time slice to energy deposits in the
            // detector.  We'll use "auto", but it's worth
            // noting that the full type of this map is
            // map<unsigned short, vector<sim::IDE>>;
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

    // Create sim plane histograms.
    vector<TH2*> sphists;
    ostringstream ssevt;
    ssevt << fevent;
    string sevt = ssevt.str();
    while ( sevt.size() < 4 ) sevt = "0" + sevt;
    for ( unsigned int irop=0; irop<fnrop; ++irop ) {
      int nchan = fropnchan[irop];
      int ntick = ftdcTickMax-ftdcTickMin;
      string hname = "h" + sevt + "sim" + fropname[irop];
      string title = fropname[irop] + ";TDC tick;Channel;Energy [MeV]";
      if ( dbg > 1 ) cout << myname << "Creating sim histo " << hname << " with " << ntick
                          << " TDC bins " << " and " << nchan << " channel bins" << endl;
      TH2* ph = tfs->make<TH2F>(hname.c_str(), title.c_str(),
                                ntick, ftdcTickMin, ftdcTickMax,
                                nchan, 0, nchan);
      sphists.push_back(ph);
    }

    // Loop over the SimChannel objects in the event.
    fscCount = 0;
    for ( auto& val : fscChannel ) val = 0;
    for ( auto& val : fscCharge ) val = 0.0;
    for ( auto& val : fscEnergy ) val = 0.0;
    for ( auto const& chan : (*simChannelHandle) ) {
      auto ichan = chan.Channel();
      // Find the channel number in the ROP.
      unsigned int irop = fnrop;
      for ( irop=0; irop<fnrop; ++irop ) {
        unsigned int lastchan = fropfirstchan[irop] + fropnchan[irop] - 1;
        if ( ichan <= lastchan ) break;
      }
      if ( irop == fnrop ) {
        cout << myname << "ERROR: SimChannel channel " << ichan << " is not in a readout plane." << endl;
        abort();
      }
      int iropchan = ichan - fropfirstchan[irop];
      if ( iropchan < 0 || iropchan >= fropnchan[irop] ) {
        cout << myname << "ERROR: ROP channel " << iropchan << " is out of range [0, "
             << fropnchan[irop] << ")" << endl;
        abort();
      }
      // Fetch the histogram for the current ROP.
      TH2* psphist = sphists[irop];
      if ( dbg > 3 ) cout << myname << " " << fscCount << ": Sim channel=" << ichan
                          << ", ROP=" << irop << ", ROP channel=" << iropchan << endl;
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
          // Update the peak channel data.
          if ( tdcPeak == 0 || tdcEnergy > tdcEnergyMax ) {
            tdcPeak = tdc;
            tdcEnergyMax = tdcEnergy;
            tdcChargeMax = tdcCharge;
            tdcNideMax = tdcNide;
          }
          // Increment the sums for this channel.
          energy += tdcEnergy;
          charge += tdcCharge;
          energytdc += tdcEnergy*tdc;
          // Fill the histogram.
          psphist->Fill(tdc, iropchan, energy);
          if ( dbg > 3 ) cout << myname << "Sim/ROP channel " << ichan << "/" << iropchan
                              << ", tdc " << tdc
                              << " has " << int(100.0*tdcEnergy+0.4999)/100.0 << " MeV" << endl;
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
        cout << myname << "WARNING: Skipping sim channel " << ichan << endl;
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
