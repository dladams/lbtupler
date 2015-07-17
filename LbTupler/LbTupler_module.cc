// LbTupler_module.cc

// David Adams
// February 2015

#ifndef LbTupler_Module
#define LbTupler_Module

#include <iomanip>
#include <sstream>
#include <utility>

// LArSoft includes
#include "SimpleTypesAndConstants/geo_types.h" // geo::View_t, geo::SignalType, geo::WireID
#include "Geometry/Geometry.h"
#include "Simulation/SimChannel.h"
#include "Simulation/LArG4Parameters.h"
#include "RecoBase/Hit.h"
#include "RecoBase/Wire.h"
#include "RecoBase/Cluster.h"
#include "RecoBase/Track.h"
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
#include <memory>

// Local includes.
#include "reducedPDG.h"
#include "intProcess.h"
#include "PlanePosition.h"
#include "GeoHelper.h"
#include "ChannelTickHistCreator.h"
#include "MCTrajectoryFollower.h"
#include "SimChannelTupler.h"
#include "TpcSignalMatcher.h"
#include "TpcSignalMatchTree.h"

using std::cout;
using std::endl;
using std::ostringstream;
using std::setw;
using std::string;
using std::vector;
using std::map;
using std::set;
using std::sqrt;
using std::shared_ptr;
using std::pair;
using geo::WireID;
using geo::PlaneID;
using geo::View_t;
using geo::kU;
using geo::kV;
using geo::kZ;
using simb::MCParticle;
using recob::Cluster;
using recob::Track;
using tpc::IndexVector;

typedef shared_ptr<TpcSignalMap> TpcSignalMapPtr;
typedef vector<TpcSignalMapPtr> TpcSignalMapVector;
typedef shared_ptr<TpcSignalMapVector> TpcSignalMapVectorPtr;
typedef TpcSignalMap::Index Index;
typedef map<int, const MCParticle*> ParticleMap;
typedef map<Index, IndexVector> IndexVectorMap;
typedef pair<TpcSignalMapPtr, TpcSignalMapVectorPtr> ClusterResult;
typedef shared_ptr<TpcSignalMatchTree> TpcSignalMatchTreePtr;
typedef art::Ptr<recob::Cluster> ClusterPtr;
typedef vector<ClusterPtr> ClusterPtrVector;

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
  void summarize2dHist(TH2* ph, string prefix, unsigned int wnam,
                       unsigned int wbin, unsigned int went) const;

  // Process a cluster container. A signal map is created for all clusters and for each cluster.
  //   event - event ID
  //   conname - name of the cluster container (producer)
  //   label - label for naming maps and histograms, e.g. "clu"
  //   phcreate - if non-null histograms are created with this
  //   wnam - width for displaying histogram names
  ClusterResult processClusters(const art::Event& evt, string conname, string label,
                                const ChannelTickHistCreator* phcreate, unsigned int wnam) const;

  // Process a track contaienr.
  void processTracks(const art::Event& evt, string conname, string label) const;

private:

  // The stuff below is the part you'll most likely have to change to
  // go from this custom example to your own task.

  // The parameters we'll read from the fcl file.
  int fdbg;                            // Debug level. Larger for more log noise.
  bool fDoTruth;                       // Read truth container.
  bool fDoMcParticleTree;              // Create the McParticle tree.
  bool fDoSimChannelTree;              // Create the SimChannel tree.
  bool fDoMcParticleSignalHists;       // Create signal histograms for McParticles
  bool fDoMcDescendantSignalHists;     // Create signal histograms for McParticle descendants
  bool fDoSimChannelSignalHists;       // Create signal histograms for SimChannels
  bool fDoRawSignalHists;              // Create signal histograms for the RawDigits
  bool fDoDeconvolutedSignalHists;     // Create signal histograms for the Wires (deconvoluted signals)
  bool fDoHitSignalHists;              // Create signal histograms for the the Hits.
  bool fDoClusterSignalHists;          // Create signal histograms for the Clusters.
  bool fDoRefClusterSignalHists;       // Create signal histograms for the reference Clusters.
  bool fDoMcParticleClusterMatching;   // Match clusters to McParticle signals.
  bool fDoMcDescendantClusterMatching; // Match clusters to McParticle descendant signals.
  bool fDoSimChannelClusterMatching;   // Match clusters to SimChannel signals.
  bool fDoRefClusterClusterMatching;   // Match clusters to reference clusters.
  string fTruthProducerLabel;          // The name of the producer that tracked simulated particles through the detector
  string fParticleProducerLabel;       // The name of the producer that tracked simulated particles through the detector
  string fSimulationProducerLabel;     // The name of the producer that tracked simulated particles through the detector
  string fHitProducerLabel;            // The name of the producer that created hits
  string fWireProducerLabel;           // The name of the producer that created wires
  string fClusterProducerLabel;        // The name of the producer that created clusters
  string fRefClusterProducerLabel;     // The name of the producer that created reference clusters
  string fTrackProducerLabel;          // The name of the producer that created tracks
  string fRawDigitProducerLabel;       // The name of the producer that created the raw digits.
  bool fUseGammaNotPi0;                // Flag to select MCParticle gamma from pi0 instead of pi0
  double fBinSize;                     // For dE/dx work: the value of dx. 

  // Derived control parameters.
  bool fDoMcParticles;             // Read MC particles.
  bool fDoSimChannels;             // Read SimChannels
  bool fDoRaw;                     // Read raw data
  bool fDoWires;                   // Read wire (deconvoluted) data
  bool fDoHits;                    // Read hit data
  bool fDoClusters;                // Read cluster data
  bool fDoTracks;                  // Read track data
  bool fDoMcParticleSelection;     // Select MC particles.
  bool fDoMcParticleSignalMaps;    // Fill MC particle signal maps.
  bool fDoMcDescendantSignalMaps;  // Fill MC particle descendant signal maps.
  bool fDoSimChannelSignalMaps;    // Fill SimChannel signal maps.
  bool fDoClusterSignalMaps;       // Fill Cluster signal maps.
  bool fDoRefClusterSignalMaps;    // Fill Cluster signal maps.

  // Pointers to the histograms we'll create. 
  TH1D* fpdgCodeHist;
  TH1D* fMomentumHist;

  // The MCParticle trajectory managers.
  MCTrajectoryFollower* m_pmctrajmc;
  MCTrajectoryFollower* m_pmctrajmd;

  // The match trees.
  TpcSignalMatchTreePtr m_ptsmtSimChannelCluster;
  TpcSignalMatchTreePtr m_ptsmtRefClusterCluster;

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
  m_pmctrajmc(nullptr), m_pmctrajmd(nullptr),
  m_sctupler(nullptr),
  fgeohelp(nullptr) {
  // Read in the parameters from the .fcl file.
  this->reconfigure(parameterSet);
}

//************************************************************************

// Destructor
LbTupler::~LbTupler() { }
   
//************************************************************************

void LbTupler::beginJob() {
  const string myname = "LbTupler::beginJob: ";
  if ( fdbg > 0 ) cout << myname << "Begin begin job." << endl;

  if ( fgeohelp == nullptr ) {
    cout << myname << "ERROR: Geometry helper is absent." << endl;
    return;
  }
  const GeoHelper& geohelp = *fgeohelp;

  // Access ART's TFileService, which will handle creating and writing
  // histograms and n-tuples for us. 
  art::ServiceHandle<art::TFileService> tfs;

  // Define the histograms. Putting semi-colons around the title
  // causes it to be displayed as the x-axis label if the histogram
  // is drawn.
  fpdgCodeHist     = tfs->make<TH1D>("pdgcodes",";PDG Code;",                  5000, -2500, 2500);
  fMomentumHist    = tfs->make<TH1D>("mom",     ";particle Momentum (GeV);",    100, 0.,    10.);

  // Define our n-tuples, which are limited forms of ROOT
  // TTrees. Start with the TTree itself.

  // Create the MCParticle trajectory manager. It builds the simulation tree and fills
  // the signal map for each selected MC particle.
  unsigned int minNptdet = 1;
  string mcptreename;
  if ( fDoMcParticleTree ) mcptreename = "McParticleTree";
  m_pmctrajmc = new MCTrajectoryFollower(fmcpdsmax, mcptreename, fgeohelp, minNptdet, 0);
  m_pmctrajmd = new MCTrajectoryFollower(fmcpdsmax, "", fgeohelp, 0, 0);

  // Match trees.
  m_ptsmtSimChannelCluster.reset(new TpcSignalMatchTree("SimChannelClusterMatch"));
  m_ptsmtRefClusterCluster.reset(new TpcSignalMatchTree("RefClusterClusterMatch"));

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

  if ( fdbg > 0 ) cout << myname << "End begin job." << endl;
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
  fdbg                           = p.get<int>("DebugLevel");
  fDoTruth                       = p.get<bool>("DoTruth");
  fDoMcParticleTree              = p.get<bool>("DoMcParticleTree");
  fDoSimChannelTree              = p.get<bool>("DoSimChannelTree");
  fDoMcParticleSignalHists       = p.get<bool>("DoMcParticleSignalHists");
  fDoMcDescendantSignalHists     = p.get<bool>("DoMcDescendantSignalHists");
  fDoSimChannelSignalHists       = p.get<bool>("DoSimChannelSignalHists");
  fDoRawSignalHists              = p.get<bool>("DoRawSignalHists");
  fDoDeconvolutedSignalHists     = p.get<bool>("DoDeconvolutedSignalHists");
  fDoHitSignalHists              = p.get<bool>("DoHitSignalHists");
  fDoClusterSignalHists          = p.get<bool>("DoClusterSignalHists");
  fDoRefClusterSignalHists       = p.get<bool>("DoRefClusterSignalHists");
  fDoTracks                      = p.get<bool>("DoTracks");
  fDoMcParticleClusterMatching   = p.get<bool>("DoMcParticleClusterMatching");
  fDoMcDescendantClusterMatching = p.get<bool>("DoMcDescendantClusterMatching");
  fDoSimChannelClusterMatching   = p.get<bool>("DoSimChannelClusterMatching");
  fDoRefClusterClusterMatching   = p.get<bool>("DoSimChannelClusterMatching");
  fTruthProducerLabel            = p.get<string>("TruthLabel");
  fParticleProducerLabel         = p.get<string>("ParticleLabel");
  fSimulationProducerLabel       = p.get<string>("SimulationLabel");
  fRawDigitProducerLabel         = p.get<string>("RawDigitLabel");
  fHitProducerLabel              = p.get<string>("HitLabel");
  fWireProducerLabel             = p.get<string>("WireLabel");
  fClusterProducerLabel          = p.get<string>("ClusterLabel");
  fRefClusterProducerLabel       = p.get<string>("RefClusterLabel");
  fTrackProducerLabel            = p.get<string>("TrackLabel");
  fUseGammaNotPi0                = p.get<bool>("UseGammaNotPi0");
  fBinSize                       = p.get<double>("BinSize");
  fscCapacity                    = p.get<double>("SimChannelSize");
  ftdcTickMin                    = p.get<int>("TdcTickMin");
  ftdcTickMax                    = p.get<int>("TdcTickMax");
  fmcpdsmax                      = p.get<double>("McParticleDsMax");
  fadcmevu                       = p.get<double>("AdcToMeVConversionU");
  fadcmevv                       = p.get<double>("AdcToMeVConversionV");
  fadcmevz                       = p.get<double>("AdcToMeVConversionZ");
  fdemaxmcp                      = p.get<double>("HistDEMaxMcParticle");
  fdemax                         = p.get<double>("HistDEMax");
  fhistusede                     = p.get<bool>("HistUseDE");

  // Derived control flags.
  fDoMcParticleSignalMaps   = fDoMcParticleSignalHists   || fDoMcParticleClusterMatching;
  fDoMcDescendantSignalMaps = fDoMcDescendantSignalHists || fDoMcDescendantClusterMatching;
  fDoSimChannelSignalMaps   = fDoSimChannelSignalHists   || fDoSimChannelClusterMatching;
  fDoClusterSignalMaps = fDoClusterSignalHists ||
                         fDoMcParticleClusterMatching || fDoMcDescendantClusterMatching ||
                         fDoSimChannelClusterMatching || fDoRefClusterClusterMatching;
  fDoRefClusterSignalMaps = fDoRefClusterSignalHists || fDoRefClusterClusterMatching;
  fDoMcParticleSelection = fDoMcParticleSignalMaps || fDoMcDescendantSignalMaps || fDoSimChannelSignalMaps ||
                           fDoMcParticleTree;
  fDoMcParticles = fDoMcParticleSelection;
  fDoSimChannels = fDoSimChannelSignalMaps || fDoSimChannelTree;
  fDoRaw = fDoRawSignalHists;
  fDoWires = fDoDeconvolutedSignalHists;
  fDoHits = fDoHitSignalHists;
  fDoClusters = fDoClusterSignalMaps;

  // Display properties.
  string sep = ": ";
  int wlab = 30;
  if ( fdbg > 0 ) {
    string prefix = myname + "  ";
    cout << endl;
    cout << myname << "Module properties:" << endl;
    cout << prefix << setw(wlab) << "DebugLevel" << sep << fdbg << endl;
    cout << prefix << setw(wlab) << "DoTruth" << sep << fDoTruth << endl;
    cout << prefix << setw(wlab) << "DoMcParticleTree" << sep << fDoMcParticleTree << endl;
    cout << prefix << setw(wlab) << "DoSimChannelTree" << sep << fDoSimChannelTree << endl;
    cout << prefix << setw(wlab) << "DoMcParticleSignalHists" << sep << fDoMcParticleSignalHists << endl;
    cout << prefix << setw(wlab) << "DoMcDescendantSignalHists" << sep << fDoMcDescendantSignalHists << endl;
    cout << prefix << setw(wlab) << "DoSimChannelSignalHists" << sep << fDoSimChannelSignalHists << endl;
    cout << prefix << setw(wlab) << "DoDeconvolutedSignalHists" << sep << fDoDeconvolutedSignalHists << endl;
    cout << prefix << setw(wlab) << "DoRawSignalHists" << sep << fDoRawSignalHists << endl;
    cout << prefix << setw(wlab) << "DoHitSignalHists" << sep << fDoHitSignalHists << endl;
    cout << prefix << setw(wlab) << "DoClusterSignalHists" << sep << fDoClusterSignalHists << endl;
    cout << prefix << setw(wlab) << "DoRefClusterSignalHists" << sep << fDoRefClusterSignalHists << endl;
    cout << prefix << setw(wlab) << "DoMcParticleClusterMatching" << sep << fDoMcParticleClusterMatching << endl;
    cout << prefix << setw(wlab) << "DoMcDescendantClusterMatching" << sep << fDoMcDescendantClusterMatching << endl;
    cout << prefix << setw(wlab) << "DoSimChannelClusterMatching" << sep << fDoSimChannelClusterMatching << endl;
    cout << prefix << setw(wlab) << "DoRefClusterClusterMatching" << sep << fDoRefClusterClusterMatching << endl;
    cout << prefix << setw(wlab) << "TruthLabel" << sep << fTruthProducerLabel << endl;
    cout << prefix << setw(wlab) << "ParticleLabel" << sep << fParticleProducerLabel << endl;
    cout << prefix << setw(wlab) << "SimulationLabel" << sep << fSimulationProducerLabel << endl;
    cout << prefix << setw(wlab) << "RawDigitLabel" << sep << fRawDigitProducerLabel << endl;
    cout << prefix << setw(wlab) << "HitLabel" << sep << fHitProducerLabel << endl;
    cout << prefix << setw(wlab) << "WireLabel" << sep << fWireProducerLabel << endl;
    cout << prefix << setw(wlab) << "ClusterLabel" << sep << fClusterProducerLabel << endl;
    cout << prefix << setw(wlab) << "RefClusterLabel" << sep << fRefClusterProducerLabel << endl;
    cout << prefix << setw(wlab) << "UseGammaNotPi0" << sep << fUseGammaNotPi0 << endl;
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

  if ( fdbg > 1 ) {
    string prefix = myname + "  ";
    cout << endl;
    cout << myname << "Derived properties:" << endl;
    cout << prefix << setw(wlab) << "DoMcParticles" << sep << fDoMcParticles << endl;
    cout << prefix << setw(wlab) << "DoSimChannels" << sep << fDoSimChannels << endl;
    cout << prefix << setw(wlab) << "DoRaw" << sep << fDoRaw << endl;
    cout << prefix << setw(wlab) << "DoWires" << sep << fDoWires << endl;
    cout << prefix << setw(wlab) << "DoHits" << sep << fDoHits << endl;
    cout << prefix << setw(wlab) << "DoClusters" << sep << fDoClusters << endl;
    cout << prefix << setw(wlab) << "DoMcParticleSelection" << sep << fDoMcParticleSelection << endl;
    cout << prefix << setw(wlab) << "DoMcParticleSignalMaps" << sep << fDoMcParticleSignalMaps << endl;
    cout << prefix << setw(wlab) << "DoMcDescendantSignalMaps" << sep << fDoMcDescendantSignalMaps << endl;
    cout << prefix << setw(wlab) << "DoSimChannelSignalMaps" << sep << fDoSimChannelSignalMaps << endl;
    cout << prefix << setw(wlab) << "DoClusterSignalMaps" << sep << fDoClusterSignalMaps << endl;
    cout << prefix << setw(wlab) << "DoRefClusterSignalMaps" << sep << fDoRefClusterSignalMaps << endl;
  }
  cout << myname << endl;
  fgeohelp = new GeoHelper(&*fGeometry, 0);
  cout << myname << "Summary from geometry helper:" << endl;
  fgeohelp->print(cout, 0, myname);
  // Geometry dump from Michelle.
  if ( fdbg > 4 ) fgeohelp->dump(cout);
  return;
}

//************************************************************************

void LbTupler::analyze(const art::Event& event) {
  const string myname = "LbTupler::analyze: ";

  // Access ART's TFileService, which will handle creating and writing
  // histograms and trees.
  art::ServiceHandle<art::TFileService> tfsHandle;
  art::TFileService* ptfs = &*tfsHandle;

  // Start by fetching some basic event information for trees and histogram labels.
  fevent  = event.id().event(); 
  fRun    = event.run();
  fSubRun = event.subRun();
  if ( fdbg > 0 ) cout << myname << "Processing run " << fRun << "-" << fSubRun
                       << ", event " << fevent << endl;

  // Create string representations of the event number.
  ostringstream ssevt;
  ssevt << fevent;
  string sevt = ssevt.str();
  string sevtf = sevt;
  while ( sevtf.size() < 4 ) sevtf = "0" + sevtf;

  // Create directory for event-level histograms.
  art::TFileDirectory htfs = ptfs->mkdir("event" + sevt);

  // Channel-tick histogram creators for the simulation data products.
  ChannelTickHistCreator hcreateSim(htfs, sevt, ftdcTickMin, ftdcTickMax, "Energy [MeV]", 0, 1.0, 20);
  ChannelTickHistCreator hcreateSimPeak(htfs, sevt, ftdcTickMin, ftdcTickMax, "Energy [MeV]", 0, 5.0, 20);

  // Channel-tick histogram creators for the reconstructed data products.
  string ztitle = "ADC counts";
  double zmax = 200;
  int ncontour = 20;
  if ( fhistusede ) {
    ztitle = "Energy [MeV]";
    zmax = fdemax;
    ncontour = 40;
  }
  ChannelTickHistCreator hcreateReco(htfs, sevt, ftdcTickMin, ftdcTickMax, ztitle, 0, zmax, ncontour);
  ChannelTickHistCreator hcreateRecoNeg(htfs, sevt, ftdcTickMin, ftdcTickMax, ztitle, -zmax, zmax, ncontour);
  ChannelTickHistCreator hcreateRecoPeak(htfs, sevt, ftdcTickMin, ftdcTickMax, ztitle, 0, 5*zmax, ncontour);

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

  // Vectors of signal maps.
  TpcSignalMapVector selectedMcTpcSignalMapsMC;    // Filled with MCParticle hits.
  TpcSignalMapVector selectedMcTpcSignalMapsMD;    // Filled with MCParticle and descendant hits.
  TpcSignalMapVector selectedMcTpcSignalMapsSC;    // Filled with SimChannel.
  TpcSignalMapVector selectedMcTpcSignalMapsMCbyROP;
  TpcSignalMapVector selectedMcTpcSignalMapsMDbyROP;

  // Particle descendant vector indexed by trackid.
  IndexVectorMap descendants;

  // Get all the MC particles for the event.
  art::Handle< vector<MCParticle> > particleHandle;
  if ( fDoMcParticles ) {
    event.getByLabel(fParticleProducerLabel, particleHandle);
    if ( fdbg > 1 ) cout << myname << "MCParticle count: " << particleHandle->size() << endl;

    // Initialize the trajectory follower for this event.
    if ( int rstatmc = m_pmctrajmc->beginEvent(event, *particleHandle) ) {
      cout << myname << "ERROR: Trajectory begin event returned " << rstatmc << endl;
      return;
    }
    if ( int rstatmd = m_pmctrajmd->beginEvent(event, *particleHandle) ) {
      cout << myname << "ERROR: Trajectory begin event returned " << rstatmd << endl;
      return;
    }

    // Create vector of selected MC particles for analysis.
    // Note that descendants typically adds nothing. The descendants are not saved for showers.
    // Use the SimChannel map to see showers.
    ParticleMap pars;   // Particle map indexed by trackid so we can find ancestors.
    if ( fDoMcParticleSelection ) {
      if ( fdbg > 1 ) cout << myname << "Selecting MC particles." << endl;
      for ( auto const& particle : (*particleHandle) ) {
        int trackid = particle.TrackId();
        pars[trackid] = &particle;
        int rpdg = reducedPDG(particle.PdgCode());
        int proc = intProcess(particle.Process());
        // Select particles.
        // 21apr2015: Keep also gammas
        // 08jul2015: Keep also gamma from initial state pi0
        // 089ul2015: Keep all with PROC=0 except pi0
        bool pi0gamma = false;  // True iff this is a pi0 from an initial-state pi0
        if ( rpdg == 6 && proc == 1 ) {
          int itrkmom = particle.Mother();
          ParticleMap::const_iterator iparmom = pars.find(itrkmom);
          if ( iparmom != pars.end() ) {
            const MCParticle& parmom = *iparmom->second;
            int procmom = intProcess(parmom.Process());
            int pdgmom = parmom.PdgCode();
            pi0gamma = pdgmom == 111 && procmom == 0;
            if ( fdbg > 2 ) {
              cout << myname << "Checked pi0 gamma: parent " << itrkmom << " has proc=" << proc
                   << " and pdg=" << pdgmom << " ==> select=" << pi0gamma << endl;
            }
          } else {
            cout << myname << "WARNING: Parent missing in particle map." << endl;
          }
        }
        // Select initial-state particles.
        bool select = proc == 0;
        // Except, use decay gammas for pi0.
        if ( fUseGammaNotPi0 ) {
          if ( rpdg == 9 ) select = false;
          else select |= pi0gamma;
        }
        if ( select ) {
          // Create an McParticle map to decide if this particle is selected.
          ostringstream ssnam;
          ssnam << "mcp";
          if ( trackid < 100 ) ssnam << 0;
          if ( trackid < 10 ) ssnam << 0;
          ssnam << trackid;
          bool usetpc = true;
          TpcSignalMapPtr pmctpmc(new TpcSignalMap(ssnam.str(), *&particle, fgeohelp, usetpc));
          int keepstat = m_pmctrajmc->addMCParticle(particle, pmctpmc.get(), false);
          // Keep tracks inside detector.
          if ( keepstat == 0 ) {
            string snam = ssnam.str();
            if ( fDoMcParticleSignalMaps ) {
              pmctpmc->buildHits();
              selectedMcTpcSignalMapsMC.push_back(pmctpmc);
            }
            TpcSignalMapPtr pmctpmd;
            if ( fDoMcDescendantSignalMaps ) {
              snam[2] = 'd';  // Use "mcd" instead of "mcp" for map with descendants
              pmctpmd.reset(new TpcSignalMap(snam, *&particle, fgeohelp, usetpc));
              m_pmctrajmd->addMCParticle(particle, pmctpmd.get(), true, &descendants[trackid]);
              pmctpmd->buildHits();
              selectedMcTpcSignalMapsMD.push_back(pmctpmd);
            }
            if ( fDoSimChannelSignalMaps ) {
              snam[2] = 's';  // Use "mcs" instead of "mcp" for SimHits.
              TpcSignalMapPtr pmctpsc(new TpcSignalMap(snam, *&particle, fgeohelp, usetpc));
              selectedMcTpcSignalMapsSC.push_back(pmctpsc);
              if ( fDoMcDescendantSignalMaps ) pmctpsc->copySegments(*pmctpmd);
            }
            if ( fdbg > 1 ) cout << myname << "  Selected";
          } else {
            if ( fdbg > 1 ) cout << myname << "  Dropped ";
          }
        } else {
          if ( fdbg > 1 ) cout << myname << "  Rejected";
        }
        if ( fdbg > 1 ) {
          cout << " MC particle " << setw(6) << trackid  << " with"
               << " PDG=" << setw(10) << particle.PdgCode()
               << " RPDG=" << setw(2) << rpdg
               << " Status=" << setw(2) << particle.StatusCode()
               << " PROC=" << setw(2) << proc
               << " parent=" << setw(3) << particle.Mother()
               << " nchild=" << setw(3) << particle.NumberDaughters()
               << " at (" << setw(3) << int(particle.Vx())
               <<    ", " << setw(3) << int(particle.Vy())
               <<    ", " << setw(3) << int(particle.Vz()) << ") cm";
          double efac = 1000.0;
          if ( abs(rpdg)==3 || rpdg==7 || rpdg==11 ) {
            cout << " Ekin=" << setw(5) << int(efac*(particle.E() - particle.Mass()));
          } else {
            cout << " Etot=" << setw(5) << int(efac*particle.E());
          }
          cout << " MeV";
          cout << endl;
        }
      }
      // Display the MC particle signal maps.
      int flag = 0;
      if ( fdbg > 2 ) flag = 11;
      if ( fDoMcParticleSignalMaps ) {
        if ( fdbg > 0 ) cout << myname << "Summary of selected MC particle signal maps (size = "
                             << selectedMcTpcSignalMapsMC.size() << "):" << endl;
        for ( auto pmctp : selectedMcTpcSignalMapsMC ) {
          pmctp->print(cout, flag, myname + "  ");
        }  // End loop over selected MC tracks
      }
      if ( fDoMcDescendantSignalMaps ) {
        if ( fdbg > 0 ) cout << myname << "Summary of selected MC particle plus descendant signal maps (size = "
                             << selectedMcTpcSignalMapsMD.size() << "):" << endl;
        for ( auto pmctp : selectedMcTpcSignalMapsMD ) {
          pmctp->print(cout, flag, myname + "  ");
        }  // End loop over selected MC tracks
      }
    }  // end DoMcParticleSelection
  
    // Split the MC signal maps by ROP.
    if ( fDoMcParticleSignalMaps ) {
      for ( const TpcSignalMapPtr& ptsm : selectedMcTpcSignalMapsMC ) {
        ptsm->splitByRop(selectedMcTpcSignalMapsMCbyROP, true);
      }
      if ( fdbg > 0 ) {
        cout << myname << "Summary of selected MC particle by ROP maps (size = "
             << selectedMcTpcSignalMapsMCbyROP.size() << "):" << endl;
        for ( auto pmctp : selectedMcTpcSignalMapsMCbyROP ) pmctp->print(cout, 0, myname + "  ");
      }
    }  // end DoMcParticleSignalMaps

    // Split the MD signal maps by ROP.
    if ( fDoMcDescendantSignalMaps ) {
      for ( const TpcSignalMapPtr& ptsm : selectedMcTpcSignalMapsMD ) {
        ptsm->splitByRop(selectedMcTpcSignalMapsMDbyROP, true);
      }
      if ( fdbg > 0 ) {
        cout << myname << "Summary of selected MD particle by ROP maps (size = "
             << selectedMcTpcSignalMapsMDbyROP.size() << "):" << endl;
        for ( auto pmctp : selectedMcTpcSignalMapsMDbyROP ) pmctp->print(cout, 0, myname + "  ");
      }
    }  // end DoMcDescendantSignalMaps
  
    // Create the event channel-tick bin histograms for all selected MC particles.
    if ( fDoMcParticleSignalHists ) {
      if ( fdbg > 1 ) cout << myname << "Create and fill MC particle histograms for all selected tracks." << endl;
      for ( unsigned int irop=0; irop<geohelp.nrop(); ++irop ) {
        TH2* ph = hcreateSim.create("mcp" + geohelp.ropName(irop), 0, geohelp.ropNChannel(irop),
                                        "MC particle signals for " + geohelp.ropName(irop));
        for ( auto pmctp : selectedMcTpcSignalMapsMC ) pmctp->fillRopChannelTickHist(ph, irop);
        if ( fdbg > 1 ) summarize2dHist(ph, myname, wnam, 4, 4);
      }
    }

    // Create the event channel-tick bin histograms for all selected MC particles with descendants.
    if ( fDoMcDescendantSignalHists ) {
      if ( fdbg > 1 ) cout << myname << "Create and fill MC descendant histograms for all selected tracks." << endl;
      for ( unsigned int irop=0; irop<geohelp.nrop(); ++irop ) {
        TH2* ph = hcreateSim.create("mcd" + geohelp.ropName(irop), 0, geohelp.ropNChannel(irop),
                                          "MC par+desc signals for " + geohelp.ropName(irop));
        for ( auto pmctp : selectedMcTpcSignalMapsMD ) pmctp->fillRopChannelTickHist(ph, irop);
        if ( fdbg > 1 ) summarize2dHist(ph, myname, wnam, 4, 4);
      }
    }

    // Create the channel-tick bin histograms for each selected MC particle by ROP.
    if ( fDoMcParticleSignalHists ) {
      if ( fdbg > 1 ) cout << myname << "Create and fill MC particle histograms for each selected track." << endl;
      for ( auto pmctp : selectedMcTpcSignalMapsMCbyROP ) {
        ostringstream ssmcp;
        ssmcp << pmctp->mcinfo()->trackID;
        string smcp = ssmcp.str();
        Index irop = pmctp->rop();
        TH2* ph = hcreateSim.create(pmctp->name(), 0, geohelp.ropNChannel(irop),
                                        "MC particle signals for " + geohelp.ropName(irop),
                                        "", "particle " + smcp, pmctp->tickRange());
        if ( ph != nullptr ) {
          pmctp->fillRopChannelTickHist(ph, irop);
          if ( fdbg > 1 ) summarize2dHist(ph, myname, wnam+8, 4, 4);
        }
      }
    }

    if ( fDoMcDescendantSignalHists ) {
      if ( fdbg > 1 ) cout << myname << "Create and fill MC descendant histograms for each selected track." << endl;
      for ( auto pmctp : selectedMcTpcSignalMapsMDbyROP ) {
        ostringstream ssmcd;
        ssmcd << pmctp->mcinfo()->trackID;
        string smcd = ssmcd.str();
          Index irop = pmctp->rop();
          TH2* ph = hcreateSim.create(pmctp->name(), 0, geohelp.ropNChannel(irop),
                                    "MD particle signals for " + geohelp.ropName(irop),
                                    "", "particle " + smcd, pmctp->tickRange());
        if ( ph != nullptr ) {
          pmctp->fillRopChannelTickHist(ph, irop);
          if ( fdbg > 1 ) summarize2dHist(ph, myname, wnam+8, 4, 4);
        }
      }  // End loop over selected MC particles by ROP
    }

  }  // end DoMcParticles

  //************************************************************************
  // Sim channels.
  //************************************************************************

  TpcSignalMapVector selectedMcTpcSignalMapsSCbyROP;
  if ( fDoSimChannels ) {

    // Get all the simulated channels for the event. These channels
    // include the energy deposited for each track.
    art::Handle<vector<sim::SimChannel>> simChannelHandle;
    event.getByLabel(fSimulationProducerLabel, simChannelHandle);
    if ( fdbg > 1 ) cout << myname << "Sim channel count: " << simChannelHandle->size() << endl;

    // Check array sizes.
    if ( simChannelHandle->size() > fscCapacity ) {
      cout << myname << "WARNING: Sim channel count exceeds TTree capacity." << endl;
    }

    // Signal map for all sim hits and same split by ROP. 
    bool usetpc = true;
    TpcSignalMapPtr ptpsim(new TpcSignalMap("allsim", fgeohelp, usetpc));
    TpcSignalMapVector tpsimByRop;
   
    // Add sim channel info and hits to the sim channel signal maps.
    // We should but don't include contributions from descendants, e.g. mu->e.
    if ( fDoSimChannelSignalMaps ) {
      if ( fdbg > 0 ) cout << myname << "Adding sim channels and hits to SimChannel signal maps (size = "
                          << simChannelHandle->size() << ")" << endl;

      // Add SimChannels to the complete signal map.
      for ( auto const& simchan : (*simChannelHandle) ) {
        ptpsim->addSimChannel(*&simchan, -1);
      }
      ptpsim->buildHits();

      // Split the complete SimChannel signal map by ROP.
      ptpsim->splitByRop(tpsimByRop, true);

      // Add SimChannels to selected-track SimChannel signal maps.
      for ( auto pmctp : selectedMcTpcSignalMapsSC ) {
        // Add the sim channel info to the selected tracks.
        for ( auto const& simchan : (*simChannelHandle) ) {
          Index tid = pmctp->mcinfo()->trackID;
          pmctp->addSimChannel(*&simchan, descendants[tid]);
        }  // End loop over sim channels in the event. 
        pmctp->buildHits();
      }  // End loop over selected SimChannel MC tracks

      // Split the selected-track SimChannel signal maps by ROP.
      for ( const TpcSignalMapPtr& ptsm : selectedMcTpcSignalMapsSC ) {
        ptsm->splitByRop(selectedMcTpcSignalMapsSCbyROP, true);
      }  

      // Display the selected-track SimChannel signal maps.
      int flag = 0;
      if ( fdbg > 2 ) flag = 11;
      if ( fdbg > 0 ) {
        cout << myname << "Summary of complete SimChannel signal map:" << endl;
        ptpsim->print(cout, flag, myname + "  ");
        cout << myname << "Summary of complete SimChannel signal maps (size = "
             << tpsimByRop.size() << ") after splitting" << endl;
        for ( auto ptsm : tpsimByRop ) {
          ptsm->print(cout, flag, myname + "  ");
        }  // End loop over signal maps
        cout << myname << "Summary of selected-track SimChannel signal maps (size = "
             << selectedMcTpcSignalMapsSC.size() << ")" << endl;
        for ( auto ptsm : selectedMcTpcSignalMapsSC ) {
          ptsm->print(cout, flag, myname + "  ");
        }  // End loop over signal maps
        cout << myname << "Summary of selected-track SimChannel signal maps (size = "
             << selectedMcTpcSignalMapsSCbyROP.size() << ") after splitting" << endl;
        for ( auto ptsm : selectedMcTpcSignalMapsSCbyROP ) {
          ptsm->print(cout, flag, myname + "  ");
        }  // End loop over signal maps
      }  // End dbg

    }  // end DoMcTpcSignalMaps

    // Create the Sim channel histograms:
    // one for each plane with all sim hits
    // one for each plane with all selected particles.
    // and one for each plane and selected particle.
    if ( fDoSimChannelSignalHists ) {
  
      // All sim hits.
      vector<TH2*> spahists;
      for ( unsigned int irop=0; irop<geohelp.nrop(); ++irop ) {
        TH2* ph = hcreateSim.create("sim" + geohelp.ropName(irop), 0, geohelp.ropNChannel(irop),
                                    "Sim channels for " + geohelp.ropName(irop));
        spahists.push_back(ph);
        for ( const auto pmctp : tpsimByRop ) {
          if ( pmctp->ropNbin(irop) == 0 ) continue;
          pmctp->fillRopChannelTickHist(ph, irop);
        }
      }
      // Display the contents of each SimChannel histogram.
      if ( fdbg > 1 ) {
        cout << myname << "Summary of complete SimChannel histograms:" << endl;
        for ( unsigned int irop=0; irop<geohelp.nrop(); ++irop ) {
          summarize2dHist(spahists[irop], myname, wnam, 4, 4);
        }
      }
  
      // All selected particles.
      vector<TH2*> sphists;
      for ( unsigned int irop=0; irop<geohelp.nrop(); ++irop ) {
        TH2* ph = hcreateSim.create("ssi" + geohelp.ropName(irop), 0, geohelp.ropNChannel(irop),
                                    "Sim channels for " + geohelp.ropName(irop));
        sphists.push_back(ph);
        for ( const auto pmctp : selectedMcTpcSignalMapsSC ) {
          if ( pmctp->ropNbin(irop) == 0 ) continue;
          pmctp->fillRopChannelTickHist(ph, irop);
        }
      }
      // Display the contents of each SimChannel histogram.
      if ( fdbg > 1 ) {
        cout << myname << "Summary of SimChannel histograms for all selected particles:" << endl;
        for ( unsigned int irop=0; irop<geohelp.nrop(); ++irop ) {
          summarize2dHist(sphists[irop], myname, wnam, 4, 4);
        }
      }
  
      // Create the Sim channel histograms: one for each plane and selected particle.
      cout << myname << "Summary of SimChannel histograms for each selected particle:" << endl;
      for ( const auto pmctp : selectedMcTpcSignalMapsSCbyROP ) {
        unsigned int itrk = pmctp->mcinfo()->trackID;
        ostringstream sstrk;
        sstrk << itrk;
        string strk = sstrk.str();
        Index irop = pmctp->rop();
        TH2* pht = hcreateSim.create(pmctp->name(), 0, geohelp.ropNChannel(irop),
                                     "Sim channels for " + geohelp.ropName(irop),
                                     "", "MC particle " + strk, pmctp->tickRange());
        if ( pht != nullptr ) {
          pmctp->fillRopChannelTickHist(pht, irop);
          summarize2dHist(pht, myname, wnam+8, 4, 4);
        }
      }
  
    }  // end DoSimChannelSignalHists

    // Fill tree.
    if ( fDoSimChannelTree ) {
      if ( fdbg > 1 ) cout << myname << "Filling SimChannel tree." << endl;
      m_sctupler->fill(event, *simChannelHandle);
    }

  }  // end DoSimChannel

  // Display the signal maps.
  if ( fDoSimChannelSignalMaps && fdbg > 2 ) {
    for ( unsigned int imcs=0; imcs<selectedMcTpcSignalMapsSC.size(); ++imcs ) {
      const auto mctsc = *selectedMcTpcSignalMapsSC.at(imcs);
      const auto mctmc = *selectedMcTpcSignalMapsMC.at(imcs);
      cout << myname << "Dumping McTpcSignalMap sim channels for event " << fevent
           << " track " << mctsc.mcinfo()->trackID
           << "\n" << myname << "  Total tick/hit signal: "
           << mctsc.tickSignal() << "/" << mctsc.hitSignal() << " MeV:" << endl;
      mctsc.print(cout,  0, myname + "  ");
       mctsc.print(cout, 12, myname + "  ");
      cout << myname << "Dumping McTpcSignalMap MC hits for event " << fevent
           << " track " << mctmc.mcinfo()->trackID
           << "\n" << myname << "  Total tick/hit signal: "
           << mctmc.tickSignal() << "/" << mctmc.hitSignal() << " MeV:" << endl;
      mctmc.print(cout,  0, myname + "  ");
      mctmc.print(cout, 11, myname + "  ");
    }  // End loop over selected MC tracks
  }  // end DoMcTpcSignalMaps

  //************************************************************************
  // Raw digits.
  //************************************************************************

  if ( fDoRaw ) {
    // Get the raw digits for the event.
    art::Handle< vector<raw::RawDigit> > rawDigitHandle;
    event.getByLabel(fRawDigitProducerLabel, rawDigitHandle);
    if ( fdbg > 1 ) cout << myname << "Raw channel count: " << rawDigitHandle->size() << endl;

    // Create the Raw digit histograms.
    if ( fDoRawSignalHists ) {
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
    }  // end DoRawSignalHists

  }  // end DoRawDigit

  //************************************************************************
  // Deconvoluted signals (aka wires).
  //************************************************************************

  if ( fDoWires ) {
    // See $LARDATA_DIR/include/RecoBase/Wire.h
    art::Handle< vector<recob::Wire> > wiresHandle;
    event.getByLabel(fWireProducerLabel, wiresHandle);
    if ( fdbg > 1 ) cout << myname << "Deconvoluted channel count: " << wiresHandle->size() << endl;

    // Create the deconvoluted signal histograms.
    if ( fDoDeconvolutedSignalHists ) {
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
    }  // end DoDeconvolutedSignalHists

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
    if ( fDoHitSignalHists ) {
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
      TpcSignalMap hitsSignalMap("allhits", fgeohelp);
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

    }  // end DoHitSignalHists

  }  // end DoHits

  //************************************************************************
  // Clusters.
  //************************************************************************

  // Reference clusters.
  TpcSignalMapVectorPtr prefClusterSignalMaps;
  if ( fDoRefClusterSignalMaps ) {
    const ChannelTickHistCreator* phcreate = fDoClusterSignalHists ? &hcreateReco : nullptr;
    TpcSignalMapVectorPtr pclusterSignalMaps;
    ClusterResult clures = processClusters(event, fRefClusterProducerLabel, "rcl", phcreate, wnam);
    prefClusterSignalMaps = clures.second;
  }

  // Clusters to evaluate.
  if ( fDoClusters ) {

    TpcSignalMapVectorPtr pclusterSignalMaps;
    const ChannelTickHistCreator* phcreate = fDoClusterSignalHists ? &hcreateReco : nullptr;
    if ( fDoClusterSignalMaps ) {
      ClusterResult clures = processClusters(event, fClusterProducerLabel, "clu", phcreate, wnam);
      pclusterSignalMaps = clures.second;
    }  // end DoClusterSignalMaps

    // Evaluate the MC cluster-finding performance.
    if ( fDoMcParticleClusterMatching ) {
      if ( fdbg > 1 ) cout << myname << "Matching MC and clusters." << endl;
      TpcSignalMatcher clumatchmc(selectedMcTpcSignalMapsMCbyROP, *pclusterSignalMaps, true, 0);
      clumatchmc.print(cout, 0);
    }

    // Evaluate the MD cluster-finding performance.
    if ( fDoMcDescendantClusterMatching ) {
      if ( fdbg > 1 ) cout << myname << "Matching MD and clusters." << endl;
      TpcSignalMatcher clumatchmd(selectedMcTpcSignalMapsMDbyROP, *pclusterSignalMaps, true, 0);
      clumatchmd.print(cout, 0);
    }

    // Evaluate the SC cluster-finding performance.
    if ( fDoSimChannelClusterMatching ) {
      if ( fdbg > 1 ) cout << myname << "Matching SC and clusters." << endl;
      TpcSignalMatcher match(selectedMcTpcSignalMapsSCbyROP, *pclusterSignalMaps, true, 0);
      match.print(cout, 0);
      if ( fdbg > 1 ) cout << myname << "Filling match tree." << endl;
      if ( m_ptsmtSimChannelCluster ) m_ptsmtSimChannelCluster->fill(event, match);
    }

    // Evaluate the reference cluster cluster-finding performance.
    if ( fDoRefClusterClusterMatching ) {
      if ( fdbg > 1 ) cout << myname << "Matching reference clusters and clusters." << endl;
      TpcSignalMatcher match(*prefClusterSignalMaps, *pclusterSignalMaps, true, 0);
      match.print(cout, 0);
      if ( fdbg > 1 ) cout << myname << "Filling match tree." << endl;
      if ( m_ptsmtRefClusterCluster ) m_ptsmtRefClusterCluster->fill(event, match);
    }

  }  // end DoClusters

  //************************************************************************
  // Tracks.
  //************************************************************************

  if ( fDoTracks ) {
    processTracks(event, fTrackProducerLabel, "trk");
  }

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
                unsigned int wnam, unsigned int wbin, unsigned int went) const {
  cout << prefix << "  " << setw(wnam) << std::left << ph->GetName()
       << std::right << " bins=" << setw(wbin) << ph->GetNbinsY() << "x" << ph->GetNbinsX()
       << ", entries=" << setw(went) << ph->GetEntries() << endl;
}

//************************************************************************

ClusterResult LbTupler::
processClusters(const art::Event& event, string conname, string label,
                const ChannelTickHistCreator* phcreateReco, unsigned int wnam) const {
  const string myname = "LbTupler::processClusters: ";
  ClusterResult out;
  GeoHelper& geohelp = *fgeohelp;

  // Get the clusters for the event.
  // See $LARDATA_DIR/include/RecoBase/Cluster.h
  art::Handle<vector<recob::Cluster>> clustersHandle;
  event.getByLabel(conname, clustersHandle);
  std::vector<art::Ptr<recob::Cluster>> clusters;
  art::fill_ptr_vector(clusters, clustersHandle);
  if ( fdbg > 1 ) cout << myname << "Cluster count: " << clusters.size() << endl;
  // Get the cluster hit associations for the event.
  art::FindManyP<recob::Hit> clusterHits(clustersHandle, event, conname);
  if ( ! clusterHits.isValid() ) {
    cout << myname << "ERROR: Cluster hit association not found." << endl;
    abort();
  }

  // Create signal maps for all clusters and for each individual cluster.
  out.first.reset(new  TpcSignalMap("all" + label, &geohelp));
  TpcSignalMap& allClusterSignalMap = *out.first;
  out.second.reset(new TpcSignalMapVector);
  TpcSignalMapVector& clusterSignalMap = *out.second;
  for ( unsigned int iclu=0; iclu<clusters.size(); ++iclu ) {
    ostringstream ssnam;
    ssnam << label;
    if ( iclu < 100 ) ssnam << 0;
    if ( iclu < 10 ) ssnam << 0;
    ssnam << iclu;
    clusterSignalMap.push_back(TpcSignalMapPtr(new TpcSignalMap(ssnam.str(), &geohelp)));
  }
  // Loop over clusters.
  cout << myname << "Looping over clusters (size = " << clusters.size() << ")" << endl;
  for ( unsigned int iclu=0; iclu<clusters.size(); ++iclu ) {
    art::Ptr<recob::Cluster> pclu = clusters[iclu];
    Index irop = geohelp.rop(pclu->Plane());
    std::vector<art::Ptr<recob::Hit>> hits = clusterHits.at(iclu);
    allClusterSignalMap.addCluster(hits);
    clusterSignalMap[iclu]->addCluster(hits);
    clusterSignalMap[iclu]->setRop(irop);
  }

  // Display the cluster signal map.
  int flag = fdbg > 2 ? 1 : 0;
  if ( fdbg > 0 ) {
    cout << myname << "Summary of cluster signal maps" << endl;
    allClusterSignalMap.print(cout, flag, myname + "   all: ");
    for ( unsigned int iclu=0; iclu<clusters.size(); ++iclu ) {
      ostringstream ssclu;
      ssclu << setw(6) << iclu << ": ";
      clusterSignalMap[iclu]->print(cout, flag, myname + ssclu.str());
    }
  }

  if ( phcreateReco != nullptr ) {
    const ChannelTickHistCreator& hcreateReco = *phcreateReco;
    // Create the channel-tick histograms for all clusters.
    if ( fdbg > 1 ) cout << myname << "Summary of cluster hit histograms:" << endl;
    for ( unsigned int irop=0; irop<geohelp.nrop(); ++irop ) {
      TH2* ph = hcreateReco.create(label + geohelp.ropName(irop), 0, geohelp.ropNChannel(irop),
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
      TpcSignalMapPtr pch = clusterSignalMap[iclu];
      for ( unsigned int irop=0; irop<geohelp.nrop(); ++irop ) {
        if ( pch->ropNbin(irop) == 0 ) {
          if ( fdbg > 2 ) cout << myname << "  Skipping " << irop << endl;
          continue;
        }
        TH2* ph = hcreateReco.create(pch->name(), 0, geohelp.ropNChannel(irop),
                                     "Cluster hits for " + geohelp.ropName(irop),
                                     //"", "cluster " + sclu, pch->tickRange());
                                     "", pch->name(), pch->tickRange());
        if ( ph != nullptr ) {
          pch->fillRopChannelTickHist(ph,irop);
          if ( fdbg > 1 ) summarize2dHist(ph, myname, wnam+10, 4, 7);
        }
      }
    }
  }
  return out;
}

//************************************************************************

// For now, this method simply dumps the track-cluster associations to the log.
// Why do I see so many repeated clusters associated with a track?

void LbTupler::
processTracks(const art::Event& event, string conname, string label) const {
  const string myname = "LbTupler::processTracks: ";
  GeoHelper& geohelp = *fgeohelp;

  // Get the tracks for the event.
  // See $LARDATA_DIR/include/RecoBase/Track.h
  art::Handle<vector<recob::Track>> tracksHandle;
  event.getByLabel(conname, tracksHandle);
  std::vector<art::Ptr<recob::Track>> tracks;
  art::fill_ptr_vector(tracks, tracksHandle);
  if ( fdbg > 1 ) cout << myname << "Track count: " << tracks.size() << endl;
  // Get the cluster hit associations for the event.
  art::FindManyP<recob::Cluster> trackClusters(tracksHandle, event, conname);
  if ( ! trackClusters.isValid() ) {
    cout << myname << "ERROR: Track-cluster association not found." << endl;
    abort();
  }
  if ( trackClusters.size() != tracks.size() ) {
    cout << myname << "ERROR: Track-cluster association size differs from tracks: "
         << trackClusters.size() << " != " << tracks.size() << endl;
    abort();
  }
  if ( fdbg > 1 ) cout << myname << "Track-cluster association count: "
                       << trackClusters.size() << endl;
  for ( unsigned int itrk=0; itrk<tracks.size(); ++itrk ) {
    //art::Ptr<recob::Track> ptrk = tracks[itrk];
    //Index irop = geohelp.rop(pclu->Plane());
    ClusterPtrVector clusters = trackClusters.at(itrk);
    cout << myname << "Track " << itrk << " cluster multiplicity: " << clusters.size() << endl;
    for ( ClusterPtr pclu : clusters ) {
      cout << myname << *pclu << endl;
    }
  }
}

//************************************************************************

  // This macro has to be defined for this module to be invoked from a
  // .fcl file; see LbTupler.fcl for more information.
  DEFINE_ART_MODULE(LbTupler)

//************************************************************************

} // namespace LbTupler

#endif // LbTupler_Module
