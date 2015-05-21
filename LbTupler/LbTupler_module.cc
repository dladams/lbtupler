// LbTupler_module.cc

// David Adams
// February 2015

#ifndef LbTupler_Module
#define LbTupler_Module

#include <iostream>
#include <iomanip>
#include <sstream>

// LArSoft includes
#include "SimpleTypesAndConstants/geo_types.h" // geo::View_t, geo::SignalType, geo::WireID
#include "Geometry/Geometry.h"
#include "Simulation/SimChannel.h"
#include "Simulation/LArG4Parameters.h"
#include "RecoBase/Hit.h"
#include "RecoBase/Wire.h"
#include "RecoBase/Cluster.h"
#include "Utilities/DetectorProperties.h"
#include "Utilities/LArProperties.h"
//#include "HitFinderLBNE/APAGeometryAlg.h"
#include "SimulationBase/MCParticle.h"
#include "SimulationBase/MCTruth.h"
#include "SimpleTypesAndConstants/geo_types.h"
#include "RawData/raw.h"
#include "RawData/RawDigit.h"

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
#include <vector>
#include <map>
#include <algorithm>
#include <iostream>
#include <string>
#include <cmath>

// Local includes.
#include "reducedPDG.h"
#include "intProcess.h"
#include "McTpcSignalMap.h"
#include "PlanePosition.h"
#include "GeoHelper.h"
#include "ChannelTickHistCreator.h"
#include "MCTrajectoryFollower.h"
#include "SimChannelTupler.h"

using std::cout;
using std::endl;
using std::ostringstream;
using std::setw;
using std::string;
using std::vector;
using std::map;
using std::set;
using std::sqrt;
using geo::WireID;
using geo::PlaneID;
using geo::View_t;
using geo::kU;
using geo::kV;
using geo::kZ;

typedef vector<TpcSignalMap> TpcSignalMapVector;

namespace LbTupler {

//**********************************************************************
// Primary class.
//**********************************************************************

class LbTupler : public art::EDAnalyzer {
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

  // Return the ADC-to-energy conversion factor for a channel.
  double adc2de(unsigned int ichan) const;

  // Display a line summarizing a 2D histogram.
  void summarize2dHist(TH2* ph, string prefix, unsigned int wnam, unsigned int wbin, unsigned int went);

private:

  // The stuff below is the part you'll most likely have to change to
  // go from this custom example to your own task.

  // The parameters we'll read from the .fcl file.
  int fdbg;                        // Debug level. Larger for more log noise.
  bool fDoTruth;                   // Read truth container.
  bool fDoMCParticles;             // Create MC particle tree and histograms.
  bool fDoMcTpcSignalMap;             // Evalueate MC track performance.
  bool fDoSimChannels;             // Create the SimChannel performance objects and histograms.
  bool fDoSimChannelTree;          // Create the SimChannel tree.
  bool fDoRaw;                     // Create the RawDigits tree and histograms.
  bool fDoWires;                   // Create the Wire histograms.
  bool fDoHits;                    // Create the Hit histograms.
  bool fDoClusters;                // Create the Cluster tree.
  string fTruthProducerLabel;      // The name of the producer that tracked simulated particles through the detector
  string fParticleProducerLabel;   // The name of the producer that tracked simulated particles through the detector
  string fSimulationProducerLabel; // The name of the producer that tracked simulated particles through the detector
  string fHitProducerLabel;        // The name of the producer that created hits
  string fWireProducerLabel;       // The name of the producer that created wires
  string fClusterProducerLabel;    // The name of the producer that created clusters
  string fRawDigitProducerLabel;   // The name of the producer that created the raw digits.
  int fSelectedPDG;                     // PDG code of particle we'll focus on
  double fBinSize;                      // For dE/dx work: the value of dx. 

  // Pointers to the histograms we'll create. 
  TH1D* fpdgCodeHist;
  TH1D* fMomentumHist;
  TH1D* fTrackLengthHist;

  // The MCParticle trajectory manager.
  MCTrajectoryFollower* m_pmctraj;

  // The n-tuples we'll create.
  SimChannelTupler* m_sctupler;
  TTree* fReconstructionNtuple;

  // The variables that will go into the n-tuple.
  int fevent;
  int fRun;
  int fSubRun;
  int fpdg;                  // PDG ID
  int ftrackid;              // Track ID

  // Other variables that will be shared between different methods.
  double                            fElectronsToGeV; // conversion factor
  double fmcpdsmax;  // Maximum step size for filling the MC particle trajectory hists
  double fadcmevu;   // MeV to ADC conversion factor for U-planes.
  double fadcmevv;   // MeV to ADC conversion factor for V-planes.
  double fadcmevz;   // MeV to ADC conversion factor for X-planes.
  bool fhistusede;   // If true, raw and wire spectra are converted to MeV.
  double fdemaxmcp;  // Max energy deposit for histogram ranges for McParticle
  double fdemax;     // Max energy deposit for histogram ranges

  // Maximum # of sim channels (for tree).
  unsigned int fscCapacity;

  // Tick range for histograms.
  unsigned int ftdcTickMin;         // First TDC bin to draw.
  unsigned int ftdcTickMax;         // Last+1 TDC bin to draw.

  // Geometry service.
  GeoHelper* fgeohelp;
  art::ServiceHandle<geo::Geometry> fGeometry;       // pointer to Geometry service

  // LArProperties service (for drift speed)
  art::ServiceHandle<util::DetectorProperties> fdetprop;

}; // class LbTupler


//************************************************************************
// class implementation
//************************************************************************

//************************************************************************

// Constructor
LbTupler::LbTupler(fhicl::ParameterSet const& parameterSet)
: EDAnalyzer(parameterSet), fdbg(0),
  m_pmctraj(nullptr), m_sctupler(nullptr), fgeohelp(nullptr) {
  // Read in the parameters from the .fcl file.
  this->reconfigure(parameterSet);
}

//************************************************************************

// Destructor
LbTupler::~LbTupler() { }
   
//************************************************************************

void LbTupler::beginJob() {
  const string myname = "LbTupler::beginJob: ";

  if ( fgeohelp == nullptr ) {
    cout << myname << "ERROR: Geometry helper is absent." << endl;
    return;
  }
  const GeoHelper& geohelp = *fgeohelp;

  // Access ART's TFileService, which will handle creating and writing
  // histograms and n-tuples for us. 
  art::ServiceHandle<art::TFileService> tfs;

  // Fetch geometry.
  unsigned int icry = 0;
  if ( fdbg > 0 ) {
    cout << myname << "             Total # TPCs: " << geohelp.ntpc() << endl;
    cout << myname << "     Total # TPC channels: " << fGeometry->Nchannels() << endl;
    cout << myname << "Total # optical detectors: " << fGeometry->NOpDet(icry) << endl;
    cout << myname << " Total # optical channels: " << fGeometry->NOpChannels() << endl;
    cout << myname << endl;
    cout << myname << "There are " << geohelp.ntpc() << " TPCs:" << endl;
    cout << myname << "      name       APA" << endl;
    for ( unsigned int itpc=0; itpc<geohelp.ntpc(); ++itpc ) {
      cout << myname << setw(10) << geohelp.tpcName(itpc) << setw(10) << geohelp.tpcApa(itpc) << endl;
    }
    cout << myname << endl;
    cout << myname << "There are " << geohelp.nrop() << " ROPs (readout planes):" << endl;
    cout << myname << "      name  1st chan     #chan  orient" << endl;
    for ( unsigned int irop=0; irop<geohelp.nrop(); ++irop ) {
      cout << myname << setw(10) << geohelp.ropName(irop) << setw(10) << geohelp.ropFirstChannel(irop)
           << setw(10) << geohelp.ropNChannel(irop) << setw(8) << geohelp.ropView(irop) << endl;
    }
  }
  double detLength = fGeometry->DetLength(); 
  double detWidth  = fGeometry->DetHalfWidth()  * 2.;
  double detHeight = fGeometry->DetHalfHeight() * 2.;
  double detVolume = sqrt( detLength*detLength + detWidth*detWidth + detHeight*detHeight );
  if ( fdbg > 0 ) {
    cout << myname << "Detector length: " << detLength << " cm" << endl;
    cout << myname << "Detector width:  " << detWidth  << " cm" << endl;
    cout << myname << "Detector height: " << detHeight << " cm" << endl;
    cout << myname << "Detector volume: " << detVolume << " cm^3" << endl;
  }

  // Geometry dump from Michelle.
  if ( fdbg > 4 ) {
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

  // Define the histograms. Putting semi-colons around the title
  // causes it to be displayed as the x-axis label if the histogram
  // is drawn.
  fpdgCodeHist     = tfs->make<TH1D>("pdgcodes",";PDG Code;",                  5000, -2500, 2500);
  fMomentumHist    = tfs->make<TH1D>("mom",     ";particle Momentum (GeV);",    100, 0.,    10.);
  fTrackLengthHist = tfs->make<TH1D>("length",  ";particle track length (cm);", 200, 0, detLength);

  // Define our n-tuples, which are limited forms of ROOT
  // TTrees. Start with the TTree itself.

  // Create the MCParticle trajectory manager. It builds the simulation tree and fills
  // the signal map for each selected MC particle.
  m_pmctraj = new MCTrajectoryFollower(fmcpdsmax, true, fgeohelp, 0);

  // Reconstruction tree.
  fReconstructionNtuple = tfs->make<TTree>("LbTuplerReconstruction","LbTuplerReconstruction");
  fReconstructionNtuple->Branch("Event",   &fevent,          "Event/I");
  fReconstructionNtuple->Branch("SubRun",  &fSubRun,         "SubRun/I");
  fReconstructionNtuple->Branch("Run",     &fRun,            "Run/I");
  fReconstructionNtuple->Branch("TrackID", &ftrackid,        "TrackID/I");
  fReconstructionNtuple->Branch("PDG",     &fpdg,            "PDG/I");

  // Sim channel tree.
  if ( fDoSimChannels && fDoSimChannelTree ) {
    m_sctupler = new SimChannelTupler(*fgeohelp, *tfs, fscCapacity);
  }
}
 
//************************************************************************

void LbTupler::beginRun(const art::Run& /*run*/) {
  // How to convert from number of electrons to GeV.  The ultimate
  // source of this conversion factor is
  // ${LARSIM_DIR}/include/SimpleTypesAndConstants/PhysicalConstants.h.
  fElectronsToGeV = 1./4.2e7;
}

//************************************************************************

void LbTupler::reconfigure(fhicl::ParameterSet const& p) {
  const string myname = "LbTupler::reconfigure: ";
  // Read parameters from the .fcl file. The names in the arguments
  // to p.get<TYPE> must match names in the .fcl file.
  fdbg                     = p.get<int>        ("DebugLevel");
  fDoTruth                 = p.get<bool>("DoTruth");
  fDoMCParticles           = p.get<bool>("DoMCParticles");
  fDoMcTpcSignalMap           = p.get<bool>("DoMcTpcSignalMap");
  fDoSimChannels           = p.get<bool>("DoSimChannels");
  fDoSimChannelTree        = p.get<bool>("DoSimChannelTree");
  fDoRaw                   = p.get<bool>("DoRaw");
  fDoWires                 = p.get<bool>("DoWires");
  fDoHits                  = p.get<bool>("DoHits");
  fDoClusters              = p.get<bool>("DoClusters");
  fTruthProducerLabel      = p.get<string>("TruthLabel");
  fParticleProducerLabel   = p.get<string>("ParticleLabel");
  fSimulationProducerLabel = p.get<string>("SimulationLabel");
  fRawDigitProducerLabel   = p.get<string>("RawDigitLabel");
  fHitProducerLabel        = p.get<string>("HitLabel");
  fWireProducerLabel       = p.get<string>("WireLabel");
  fClusterProducerLabel    = p.get<string>("ClusterLabel");
  fSelectedPDG             = p.get<int>("PDGcode");
  fBinSize                 = p.get<double>("BinSize");
  fscCapacity              = p.get<double>("SimChannelSize");
  ftdcTickMin              = p.get<int>("TdcTickMin");
  ftdcTickMax              = p.get<int>("TdcTickMax");
  fmcpdsmax                = p.get<double>("McParticleDsMax");
  fadcmevu                 = p.get<double>("AdcToMeVConversionU");
  fadcmevv                 = p.get<double>("AdcToMeVConversionV");
  fadcmevz                 = p.get<double>("AdcToMeVConversionZ");
  fdemaxmcp                = p.get<double>("HistDEMaxMcParticle");
  fdemax                   = p.get<double>("HistDEMax");
  fhistusede               = p.get<bool>("HistUseDE");
  string sep = ": ";
  int wlab = 20;
  if ( fdbg > 0 ) {
    string prefix = myname + "  ";
    cout << myname << setw(wlab) << "Module properties:" << endl;
    cout << prefix << setw(wlab) << "DebugLevel" << sep << fdbg << endl;
    cout << prefix << setw(wlab) << "DoMCParticles" << sep << fDoMCParticles << endl;
    cout << prefix << setw(wlab) << "DoMcTpcSignalMap" << sep << fDoMcTpcSignalMap << endl;
    cout << prefix << setw(wlab) << "DoSimChannels" << sep << fDoSimChannels << endl;
    cout << prefix << setw(wlab) << "DoSimChannelTree" << sep << fDoSimChannelTree << endl;
    cout << prefix << setw(wlab) << "DoRaw" << sep << fDoRaw << endl;
    cout << prefix << setw(wlab) << "DoWires" << sep << fDoWires << endl;
    cout << prefix << setw(wlab) << "DoHits" << sep << fDoHits << endl;
    cout << prefix << setw(wlab) << "DoClusters" << sep << fDoClusters << endl;
    cout << prefix << setw(wlab) << "TruthLabel" << sep << fTruthProducerLabel << endl;
    cout << prefix << setw(wlab) << "ParticleLabel" << sep << fParticleProducerLabel << endl;
    cout << prefix << setw(wlab) << "SimulationLabel" << sep << fSimulationProducerLabel << endl;
    cout << prefix << setw(wlab) << "RawDigitLabel" << sep << fRawDigitProducerLabel << endl;
    cout << prefix << setw(wlab) << "HitLabel" << sep << fHitProducerLabel << endl;
    cout << prefix << setw(wlab) << "WireLabel" << sep << fWireProducerLabel << endl;
    cout << prefix << setw(wlab) << "ClusterLabel" << sep << fClusterProducerLabel << endl;
    cout << prefix << setw(wlab) << "PDGcode" << sep << fSelectedPDG << endl;
    cout << prefix << setw(wlab) << "BinSize" << sep << fBinSize << endl;
    cout << prefix << setw(wlab) << "SimChannelSize" << sep << fscCapacity << endl;
    cout << prefix << setw(wlab) << "TdcTickMin" << sep << ftdcTickMin << endl;
    cout << prefix << setw(wlab) << "TdcTickMax" << sep << ftdcTickMax << endl;
    cout << prefix << setw(wlab) << "McParticleDsMax" << sep << fmcpdsmax << endl;
    cout << prefix << setw(wlab) << "AdcToMeVConversionU" << sep << fadcmevu << endl;
    cout << prefix << setw(wlab) << "AdcToMeVConversionV" << sep << fadcmevv << endl;
    cout << prefix << setw(wlab) << "AdcToMeVConversionZ" << sep << fadcmevz << endl;
    cout << prefix << setw(wlab) << "HistDEMaxMcParticle" << sep << fdemaxmcp << endl;
    cout << prefix << setw(wlab) << "HistDEMax" << sep << fdemax << endl;
    cout << prefix << setw(wlab) << "HistUseDE" << sep << fhistusede << endl;
  }

  cout << myname << endl;
  cout << myname << "Summary from geometry helper:" << endl;
  fgeohelp = new GeoHelper(&*fGeometry, 1);
  fgeohelp->print(cout, 0, myname);
  return;
}

//************************************************************************

void LbTupler::analyze(const art::Event& event) {
  const string myname = "LbTupler::analyze: ";

  // Access ART's TFileService, which will handle creating and writing
  // histograms and trees.
  art::ServiceHandle<art::TFileService> tfs;
  art::TFileService* ptfs = &*tfs;

  // Start by fetching some basic event information for trees and histogram labels.
  fevent  = event.id().event(); 
  fRun    = event.run();
  fSubRun = event.subRun();
  if ( fdbg > 0 ) cout << myname << "Processing run " << fRun << "-" << fSubRun
                       << ", event " << fevent << endl;

  // Create string representation of the event number.
  ostringstream ssevt;
  ssevt << fevent;
  string sevt = ssevt.str();
  string sevtf = sevt;
  while ( sevtf.size() < 4 ) sevtf = "0" + sevtf;

  // Channel-tick histogram for the simulation data products.
  ChannelTickHistCreator hcreateSim(ptfs, sevt, ftdcTickMin, ftdcTickMax, "Energy [MeV]", 0, 1.0, 20);
  ChannelTickHistCreator hcreateSimPeak(ptfs, sevt, ftdcTickMin, ftdcTickMax, "Energy [MeV]", 0, 5.0, 20);

  // Channel-tick histogram creators for the reconstructed data products.
  string ztitle = "ADC counts";
  double zmax = 150;
  int ncontour = 30;
  if ( fhistusede ) {
    ztitle = "Energy [MeV]";
    zmax = fdemax;
    ncontour = 40;
  }
  ChannelTickHistCreator hcreateReco(ptfs, sevt, ftdcTickMin, ftdcTickMax, ztitle, 0, zmax, ncontour);
  ChannelTickHistCreator hcreateRecoNeg(ptfs, sevt, ftdcTickMin, ftdcTickMax, ztitle, -zmax, zmax, ncontour);
  ChannelTickHistCreator hcreateRecoPeak(ptfs, sevt, ftdcTickMin, ftdcTickMax, ztitle, 0, 5*zmax, ncontour);

  // Formatting.
  int wnam = 12 + sevtf.size();                  // Base width for a histogram name.

  // Check gemetry helper.
  if ( fgeohelp == nullptr ) {
    cout << myname << "ERROR: Geometry helper is absent." << endl;
    return;
  }
  const GeoHelper& geohelp = *fgeohelp;

  //************************************************************************
  // MC Truth
  //************************************************************************

  if ( fDoTruth ) {
    // Get the MC Truth for the event.
    // See $NUTOOLS_DIR/include/SimulationBase/MCTruth.h
    art::Handle< vector<simb::MCTruth> > truthHandle;
    event.getByLabel(fTruthProducerLabel, truthHandle);
    if ( fdbg > 1 ) cout << myname << "Truth count: " << truthHandle->size() << endl;
  }

  //************************************************************************
  // MC particles
  //************************************************************************

  // Get all the MC particles for the event.
  art::Handle< vector<simb::MCParticle> > particleHandle;
  if ( fDoMCParticles || fDoMcTpcSignalMap ) {
    event.getByLabel(fParticleProducerLabel, particleHandle);
    if ( fdbg > 1 ) cout << myname << "MCParticle count: " << particleHandle->size() << endl;
  }

  // Initialize the trajectory follower for this event.
  int rstat = m_pmctraj->beginEvent(event, *particleHandle);
  if ( rstat ) {
    cout << myname << "ERROR: Trajectory begin event returned " << rstat << endl;
    return;
  }

  // Create vector of selected MC particles for analysis.
  vector<McTpcSignalMap> selectedMcTpcSignalMapsMC;    // Filled with MCParticle hits.
  vector<McTpcSignalMap> selectedMcTpcSignalMapsSC;    // Filled with SimChannel.
  if ( fDoMcTpcSignalMap ) {
    if ( fdbg > 1 ) cout << myname << "Selecting MC particles." << endl;
    for ( auto const& particle : (*particleHandle) ) {
      int trackid = particle.TrackId();
      int rpdg = reducedPDG(particle.PdgCode());
      int proc = intProcess(particle.Process());
      // Select particles.
      // 21apr2015: Keep also gammas
      if ( proc == 0 && rpdg < 7 ) {
        if ( fdbg > 1 ) cout << myname << "  Selecting";
        // Add selected track to the MCParticle performance list.
        selectedMcTpcSignalMapsMC.emplace(selectedMcTpcSignalMapsMC.end(), *&particle, fgeohelp);
        // Fill the MC performance information including signal map.
        McTpcSignalMap& mctp = selectedMcTpcSignalMapsMC.back();
        m_pmctraj->addMCParticle(particle, &mctp);
        mctp.buildHits();
        // Add selected track to the SimHists performance list.
        selectedMcTpcSignalMapsSC.emplace(selectedMcTpcSignalMapsSC.end(), *&particle, fgeohelp);
      } else {
        if ( fdbg > 1 ) cout << myname << "  Rejecting";
      }
      if ( fdbg > 1 ) cout << " MC particle " << setw(6) << trackid 
                          << " with RPDG=" << setw(2) << rpdg
                          << " and PROC=" << setw(2) << proc
                          << " and PDG=" << setw(10) << particle.PdgCode()
                          << endl;
    }
    // Display the selected track performance objects.
    int flag = 0;
    if ( fdbg > 2 ) flag = 11;
    if ( fdbg > 0 ) cout << myname << "Summary of selected MC track performance objects (size = "
                         << selectedMcTpcSignalMapsMC.size() << "):" << endl;
    for ( auto& mctrackperf : selectedMcTpcSignalMapsMC ) {
      mctrackperf.print(cout, flag, myname + "  ");
    }  // End loop over selected MC tracks
  }  // end DoMcTpcSignalMap

  // MCParticles histograms and tree.
  if ( fDoMCParticles ) {

    // Create the event channel-tick bin histograms for all selected MC particles.
    if ( fdbg > 1 ) cout << myname << "Create and fill MCParticle histograms for all selected tracks." << endl;
    for ( unsigned int irop=0; irop<geohelp.nrop(); ++irop ) {
      TH2* ph = hcreateSimPeak.create("mcp" + geohelp.ropName(irop), 0, geohelp.ropNChannel(irop),
                                      "MC particle signals for " + geohelp.ropName(irop));
      for ( auto& mctrackperf : selectedMcTpcSignalMapsMC ) {
        mctrackperf.fillRopChannelTickHist(ph, irop);
      }
      if ( fdbg > 1 ) summarize2dHist(ph, myname, wnam, 4, 4);
    }

    // Create the channel-tick bin histograms for each selected MC particle.
    if ( fdbg > 1 ) cout << myname << "Create and fill MCParticle histograms for each selected track." << endl;
    for ( auto& mctrackperf : selectedMcTpcSignalMapsMC ) {
      ostringstream ssmcp;
      ssmcp << mctrackperf.trackID();
      string smcp = ssmcp.str();
      for ( unsigned int irop=0; irop<geohelp.nrop(); ++irop ) {
        if ( mctrackperf.ropNbin(irop) == 0 ) continue;
        TH2* ph = hcreateSimPeak.create("mcp" + geohelp.ropName(irop), 0, geohelp.ropNChannel(irop),
                                        "MC particle signals for " + geohelp.ropName(irop),
                                        "mct" + smcp, "particle " + smcp, mctrackperf.tickRange());
        mctrackperf.fillRopChannelTickHist(ph, irop);
        if ( fdbg > 1 ) summarize2dHist(ph, myname, wnam+8, 4, 4);
      }  // End loop over ROPs
    }  // End loop over selected MC particles

  }  // end Do MCParticle

  //************************************************************************
  // Sim channels.
  //************************************************************************

  if ( fDoSimChannels ) {

    // Get all the simulated channels for the event. These channels
    // include the energy deposited for each track.
    art::Handle<vector<sim::SimChannel>> simChannelHandle;
    if ( fDoMCParticles || fDoSimChannels || fDoMcTpcSignalMap ) {
      event.getByLabel(fSimulationProducerLabel, simChannelHandle);
      if ( fdbg > 1 ) cout << myname << "Sim channel count: " << simChannelHandle->size() << endl;
    }

    // Check array sizes.
    if ( simChannelHandle->size() > fscCapacity ) {
      cout << myname << "WARNING: Sim channel count exceeds TTree capacity." << endl;
    }

    // Add sim channel info and hits to the sim channel performance objects.
    if ( fDoMcTpcSignalMap ) {
      if ( fdbg > 0 ) cout << myname << "Adding sim channels and hits to SimChannel performance objects (size = "
                          << simChannelHandle->size() << ")" << endl;
      for ( auto& mctp : selectedMcTpcSignalMapsSC ) {
        // Add the sim channel info to the selected tracks.
        for ( auto const& simchan : (*simChannelHandle) ) {
          mctp.addSimChannel(*&simchan);
        }  // End loop over sim channels in the event. 
        mctp.buildHits();
      }  // End loop over selected SimChannel MC tracks
      // Display the selected track performance objects.
      int flag = 0;
      if ( fdbg > 2 ) flag = 11;
      if ( fdbg > 0 ) {
        cout << myname << "Summary of selected-track SimChannel performance objects (size = "
             << selectedMcTpcSignalMapsSC.size() << ")" << endl;
        for ( auto& mctp : selectedMcTpcSignalMapsSC ) {
          mctp.print(cout, flag, myname + "  ");
        }  // End loop over selected MC tracks
      }
    }  // end DoMcTpcSignalMap

    // Create the Sim channel histograms: one for each plane with all selected particles
    // and one for each plane and selected particle.
    cout << myname << "Summary of SimChannel histograms for each selected particle:" << endl;
    vector<TH2*> sphists;
    for ( unsigned int irop=0; irop<geohelp.nrop(); ++irop ) {
      TH2* ph = hcreateSim.create("sim" + geohelp.ropName(irop), 0, geohelp.ropNChannel(irop),
                                  "Sim channels for " + geohelp.ropName(irop));
      sphists.push_back(ph);
      for ( const auto& mctp : selectedMcTpcSignalMapsSC ) {
        if ( mctp.ropNbin(irop) == 0 ) continue;
        mctp.fillRopChannelTickHist(ph, irop);
        unsigned int itrk = mctp.trackID();
        ostringstream sstrk;
        sstrk << itrk;
        string strk = sstrk.str();
        TH2* pht = hcreateSim.create("sim" + geohelp.ropName(irop), 0, geohelp.ropNChannel(irop),
                                     "Sim channels for " + geohelp.ropName(irop),
                                     "mct"+strk, "MC particle " + strk, mctp.tickRange());
        mctp.fillRopChannelTickHist(pht, irop);
        summarize2dHist(pht, myname, wnam+8, 4, 4);
      }
    }

    // Display the contents of each SimChannel histogram.
    if ( fdbg > 1 ) {
      cout << myname << "Summary of SimChannel histograms:" << endl;
      for ( unsigned int irop=0; irop<geohelp.nrop(); ++irop ) {
        summarize2dHist(sphists[irop], myname, wnam, 4, 4);
      }
    }

    // Fill tree.
    if ( fDoSimChannelTree ) {
      m_sctupler->fill(event, *simChannelHandle);
      if ( fdbg > 1 ) cout << myname << "Filling SimChannel tree." << endl;
    }

  }  // end DoSimChannel

  // Display the performance objects.
  if ( fDoMcTpcSignalMap && fdbg > 2 ) {
    for ( unsigned int imcs=0; imcs<selectedMcTpcSignalMapsSC.size(); ++imcs ) {
      const auto mctsc = selectedMcTpcSignalMapsSC.at(imcs);
      const auto mctmc = selectedMcTpcSignalMapsMC.at(imcs);
      cout << myname << "Dumping McTpcSignalMap sim channels for event " << fevent
           << " track " << mctsc.trackID()
           << "\n" << myname << "  Total tick/hit signal: "
           << mctsc.tickSignal() << "/" << mctsc.hitSignal() << " MeV:" << endl;
      mctsc.print(cout,  0, myname + "  ");
       mctsc.print(cout, 12, myname + "  ");
      cout << myname << "Dumping McTpcSignalMap MC hits for event " << fevent
           << " track " << mctmc.trackID()
          << "\n" << myname << "  Total tick/hit signal: "
           << mctmc.tickSignal() << "/" << mctmc.hitSignal() << " MeV:" << endl;
      mctmc.print(cout,  0, myname + "  ");
      mctmc.print(cout, 11, myname + "  ");
    }  // End loop over selected MC tracks
  }  // end DoMcTpcSignalMap

  //************************************************************************
  // Raw digits.
  //************************************************************************

  if ( fDoRaw ) {
    // Get the raw digits for the event.
    art::Handle< vector<raw::RawDigit> > rawDigitHandle;
    event.getByLabel(fRawDigitProducerLabel, rawDigitHandle);
    if ( fdbg > 1 ) cout << myname << "Raw channel count: " << rawDigitHandle->size() << endl;

    // Create the Raw digit histograms.
    vector<TH2*> rawhists;
    for ( unsigned int irop=0; irop<geohelp.nrop(); ++irop ) {
      TH2* ph = hcreateRecoNeg.create("raw" + geohelp.ropName(irop), 0, geohelp.ropNChannel(irop),
                                   "Raw signals for " + geohelp.ropName(irop));
      rawhists.push_back(ph);
    }

    for ( auto const& digit : (*rawDigitHandle) ) {
      int ichan = digit.Channel();
      unsigned int irop = geohelp.channelRop(ichan);
      TH2* ph = rawhists[irop];
      unsigned int iropchan = ichan - geohelp.ropFirstChannel(irop);
      int nadc = digit.NADC();
      vector<short> adcs;
      raw::Uncompress(digit.ADCs(), adcs, digit.Compression());
      unsigned int nzero = 0;
      for ( auto adc : adcs ) if ( adc == 0.0 ) ++nzero;
      if ( fdbg > 3 ) cout << myname << "Digit channel " << ichan
                          << " (ROP-chan = " << irop << "-" << iropchan
                          << ") has " << nadc << " ADCs and "
                          << digit.Samples() << " samples. Uncompressed size is " << adcs.size()
                          << ". Number filled is " << adcs.size()-nzero << endl;
      for ( unsigned int tick=0; tick<adcs.size(); ++tick ) {
        double wt = adcs[tick];
        if ( wt == 0 ) continue;
        if ( fhistusede ) wt *= adc2de(ichan);
        ph->Fill(tick, iropchan, wt);
      }
    }

    // Display the contents of each raw data histogram.
    if ( fdbg > 1 ) {
      cout << myname << "Summary of raw data histograms:" << endl;
      for ( unsigned int irop=0; irop<geohelp.nrop(); ++irop ) {
        summarize2dHist(rawhists[irop], myname, wnam, 4, 7);
      }
    }

  }  // end DoRawDigit

  //************************************************************************
  // Deconvoluted signals (aka wires).
  //************************************************************************

  if ( fDoWires ) {
    // See $LARDATA_DIR/include/RecoBase/Wire.h
    art::Handle< vector<recob::Wire> > wiresHandle;
    event.getByLabel(fWireProducerLabel, wiresHandle);
    if ( fdbg > 1 ) cout << myname << "Deconvoluted channel count: " << wiresHandle->size() << endl;

    // Create the wire histograms.
    vector<TH2*> dcohists;
    for ( unsigned int irop=0; irop<geohelp.nrop(); ++irop ) {
      TH2* ph = hcreateReco.create("dco" + geohelp.ropName(irop), 0, geohelp.ropNChannel(irop),
                                   "Deconvoluted signals for " + geohelp.ropName(irop));
      dcohists.push_back(ph);
    }

    for ( auto const& wire : (*wiresHandle) ) {
      int ichan = wire.Channel();
      unsigned int irop = geohelp.channelRop(ichan);
      unsigned int iropchan = ichan - geohelp.ropFirstChannel(irop);
      auto sigs = wire.Signal();
      const auto& roisigs = wire.SignalROI();
      TH2* ph = dcohists[irop];
      if ( fdbg > 3 ) cout << myname << "Deconvoluted channel " << ichan
                          << " (ROP-chan = " << irop << "-" << iropchan << ")"
                          << " with view " << wire.View()
                          << " has " << sigs.size() << " signals"
                          << " and " << roisigs.size() << " ROIs"
                          << "." << endl;
      for ( unsigned int tick=0; tick<sigs.size(); ++tick ) {
        double wt = roisigs[tick];
        if ( wt == 0 ) continue;
        if ( fhistusede ) wt *= adc2de(ichan);
        ph->Fill(tick, iropchan, wt);
      }
    }

    // Display the contents of each deconvoluted signal histogram.
    if ( fdbg > 1 ) {
      cout << myname << "Summary of deconvoluted data histograms:" << endl;
      for ( unsigned int irop=0; irop<geohelp.nrop(); ++irop ) {
        summarize2dHist(dcohists[irop], myname, wnam, 4, 7);
      }
    }

  }  // end DoWires

  //************************************************************************
  // Hits.
  //************************************************************************

  if ( fDoHits ) {
    // Get the hits for the event.
    // See $LARDATA_DIR/include/RecoBase/Hit.h
    art::Handle< vector<recob::Hit> > hitsHandle;
    event.getByLabel(fHitProducerLabel, hitsHandle);
    if ( fdbg > 1 ) cout << myname << "Hit count: " << hitsHandle->size() << endl;

    // Create the hit histograms.
    vector<TH2*> hithists;
    for ( unsigned int irop=0; irop<geohelp.nrop(); ++irop ) {
      TH2* ph = hcreateRecoPeak.create("hip" + geohelp.ropName(irop), 0, geohelp.ropNChannel(irop),
                                       "Hit peaks for " + geohelp.ropName(irop));
      hithists.push_back(ph);
    }

    for ( auto const& hit : (*hitsHandle) ) {
      int ichan = hit.Channel();
      unsigned int irop = geohelp.channelRop(ichan);
      unsigned int iropchan = ichan - geohelp.ropFirstChannel(irop);
      TH2* ph = hithists[irop];
      if ( fdbg > 3 ) cout << myname << "Hit channel " << ichan
                          << " (ROP-chan = " << irop << "-" << iropchan << ")"
                          << " with view " << hit.View()
                          << " has charge " << hit.SummedADC() 
                          << " at TDC " << hit.PeakTime()
                          << "." << endl;
      double wt = hit.SummedADC();
      if ( wt == 0 ) continue;
      if ( fhistusede ) wt *= adc2de(ichan);
      if ( fdbg > 3 ) cout << myname << "    Hit histo " << ph->GetName() << " time/channel/wt = "
                           << hit.PeakTime() << "/" << iropchan << "/" << wt << endl;
      ph->Fill(hit.PeakTime(), iropchan, wt);
    }

    if ( fdbg > 1 ) cout << myname << "Summary of hit peak histograms:" << endl;
    for ( unsigned int irop=0; irop<geohelp.nrop(); ++irop ) {
      TH2* ph = hithists[irop];
      if ( fdbg > 1 ) summarize2dHist(ph, myname, wnam, 4, 7);
    }

    // Put all hits into a signal map.
    TpcSignalMap hitsSignalMap(fgeohelp);
    for ( auto const& hit : (*hitsHandle) ) {
      hitsSignalMap.addHit(hit, 0);
    }

    // Display the hits performance
    int flag = 0;
    if ( fdbg > 2 ) flag = 1;
    if ( fdbg > 0 ) {
      cout << myname << "Summary of hit performance" << endl;
      hitsSignalMap.print(cout, flag, myname + "  ");
    }

    // Create the new hit histograms.
    TH2* phallhits = hcreateReco.create("hsgall", 0, fGeometry->Nchannels(), "Hits for ");
    hitsSignalMap.fillChannelTickHist(phallhits);
    if ( fdbg > 1 ) cout << myname << "Summary of hit histograms:" << endl;
    for ( unsigned int irop=0; irop<geohelp.nrop(); ++irop ) {
      TH2* ph = hcreateReco.create("hit" + geohelp.ropName(irop), 0, geohelp.ropNChannel(irop),
                                   "Hits for " + geohelp.ropName(irop));
      hitsSignalMap.fillRopChannelTickHist(ph,irop);
      if ( fdbg > 1 ) summarize2dHist(ph, myname, wnam, 4, 7);
    }
    if ( fdbg > 1 ) summarize2dHist(phallhits, myname, wnam, 4, 7);

  }  // end DoHits

  //************************************************************************
  // Clusters.
  //************************************************************************

  if ( fDoClusters ) {

    // Get the clusters for the event.
    // See $LARDATA_DIR/include/RecoBase/Cluster.h
    art::Handle<vector<recob::Cluster>> clustersHandle;
    event.getByLabel(fClusterProducerLabel, clustersHandle);
    std::vector<art::Ptr<recob::Cluster>> clusters;
    art::fill_ptr_vector(clusters, clustersHandle);
    if ( fdbg > 1 ) cout << myname << "Cluster count: " << clusters.size() << endl;
    // Get the cluster hit associations for the event.
    art::FindManyP<recob::Hit> clusterHits(clustersHandle, event, fClusterProducerLabel);
    if ( ! clusterHits.isValid() ) {
      cout << myname << "ERROR: Cluster hit association not found." << endl;
      abort();
    }

    // Create performance object for all clusters and for each individual cluster.
    TpcSignalMap allClusterSignalMap(&geohelp);
    TpcSignalMapVector clusterSignalMap(clusters.size(), &geohelp);
    // Loop over clusters.
    cout << myname << "Looping over clusters (size = " << clusters.size() << ")" << endl;
    for ( unsigned int iclu=0; iclu<clusters.size(); ++iclu ) {
      //art::Ptr<recob::Cluster> pclu = clusters[iclu];
      std::vector<art::Ptr<recob::Hit>> hits = clusterHits.at(iclu);
      allClusterSignalMap.addCluster(hits);
      clusterSignalMap[iclu].addCluster(hits);
    }

    // Display the hits performance
    int flag = fdbg > 2 ? 1 : 0;
    if ( fdbg > 0 ) {
      cout << myname << "Summary of hit performance" << endl;
      for ( unsigned int iclu=0; iclu<clusters.size(); ++iclu ) {
        ostringstream ssclu;
        ssclu << setw(6) << iclu << ":";
        clusterSignalMap[iclu].print(cout, flag, myname + ssclu.str());
      }
      allClusterSignalMap.print(cout, flag, myname + "   all:");
    }

    // Create the channel-tick histograms for all clusters.
    if ( fdbg > 1 ) cout << myname << "Summary of cluster hit histograms:" << endl;
    for ( unsigned int irop=0; irop<geohelp.nrop(); ++irop ) {
      TH2* ph = hcreateReco.create("clu" + geohelp.ropName(irop), 0, geohelp.ropNChannel(irop),
                                   "Cluster hits for " + geohelp.ropName(irop));
      allClusterSignalMap.fillRopChannelTickHist(ph,irop);
      if ( fdbg > 1 ) summarize2dHist(ph, myname, wnam, 4, 7);
    }

    // Create the channel-tick histograms for each cluster.
    if ( fdbg > 1 ) cout << myname << "Summary of per-cluster hit histograms:" << endl;
    for ( unsigned int iclu=0; iclu<clusters.size(); ++iclu ) {
      ostringstream ssclu;
      ssclu << iclu;
      string sclu = ssclu.str();
      TpcSignalMap& ch = clusterSignalMap[iclu];
      for ( unsigned int irop=0; irop<geohelp.nrop(); ++irop ) {
        if ( ch.ropNbin(irop) == 0 ) {
          if ( fdbg > 2 ) cout << myname << "  Skipping " << irop << endl;
          continue;
        }
        TH2* ph = hcreateReco.create("clu" + geohelp.ropName(irop), 0, geohelp.ropNChannel(irop),
                                     "Cluster hits for " + geohelp.ropName(irop),
                                     "clu" + sclu, "cluster " + sclu, ch.tickRange());
        ch.fillRopChannelTickHist(ph,irop);
        if ( fdbg > 1 ) summarize2dHist(ph, myname, wnam+10, 4, 7);
      }
    }

  }  // end DoClusters

  //************************************************************************
  // Done.
  //************************************************************************

  return;
}

//************************************************************************
// Return the ADC-to-energy conversion factor for a channel.
//************************************************************************

double LbTupler::adc2de(unsigned int ichan) const {
  string myname = "LbTupler::adc2de: ";
  const GeoHelper& geohelp = *fgeohelp;
  double cfac = 1.0;
  unsigned int irop = geohelp.channelRop(ichan);
  View_t view = geohelp.ropView(irop);
  if      ( view == kU ) cfac = fadcmevu;
  else if ( view == kV ) cfac = fadcmevv;
  else if ( view == kZ ) cfac = fadcmevz;
  else {
    cout << myname << "ERROR: plane does not have specified orientation: " << irop << endl;
    abort();
  }
  return cfac;
}

//************************************************************************

void LbTupler::
summarize2dHist(TH2* ph, string prefix,
                unsigned int wnam, unsigned int wbin, unsigned int went) {
  cout << prefix << "  " << setw(wnam) << std::left << ph->GetName()
       << std::right << " bins=" << setw(wbin) << ph->GetNbinsY() << "x" << ph->GetNbinsX()
       << ", entries=" << setw(went) << ph->GetEntries() << endl;
}

//************************************************************************

  // This macro has to be defined for this module to be invoked from a
  // .fcl file; see LbTupler.fcl for more information.
  DEFINE_ART_MODULE(LbTupler)

//************************************************************************

} // namespace LbTupler

#endif // LbTupler_Module
