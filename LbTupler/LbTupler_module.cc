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

#include "reducedPDG.h"
#include "intProcess.h"
#include "MCTrackPerf.h"
#include "PlanePosition.h"

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

#include "GeoHelper.h"

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

private:

  // The stuff below is the part you'll most likely have to change to
  // go from this custom example to your own task.

  // The parameters we'll read from the .fcl file.
  int fdbg;                        // Debug level. Larger for more log noise.
  bool fDoTruth;                   // Read truth container.
  bool fDoMCParticles;             // Create MC particle tree and histograms.
  bool fDoMCTrackPerf;             // Evalueate MC track performance.
  bool fDoSimChannels;             // Create the SimChannel tree and histograms.
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

  // The n-tuples we'll create.
  TTree* fSimulationNtuple;
  TTree* fSimChannelNtuple;
  TTree* fReconstructionNtuple;

  // The variables that will go into the n-tuple.
  int fevent;
  int fRun;
  int fSubRun;
  int fpdg;                  // PDG ID
  int frpdg;                 // reduced PDG ID
  int fproc;                 // Process
  int ftrackid;              // Track ID
  int fparent;               // Parent track ID
  unsigned int fnchild;      // # children
  unsigned int fndetchild;   // # children in detector
  int fndetin;               // # TPC entries from non-TPC
  int fndetout;              // # TPC exits to non-TPC
  int fntpcin;               // # TPC entries from non-TPC or other TPC
  int fntpcout;              // # TPC exits to not-TPC or other TPC
  int fncryin;               // # cryo entries from non-cryo or other cryo
  int fncryout;              // # cryo exits to non-cryo or other cryo
  static const unsigned int fmaxchild = 20;
  int fchild[fmaxchild];
  int fdetchild[fmaxchild];

  // Arrays for 4-vectors: (x,y,z,t) and (Px,Py,Pz,E).
  // Note: old-style C++ arrays are considered obsolete. However,
  // to create simple n-tuples, we still need to use them. 
  float fStartXYZT[4];
  float fEndXYZT[4];
  float fStartPE[4];
  float fEndPE[4];
  static const unsigned int maxpt = 20000;
  unsigned int fnpt;     // # points in trajectory
  unsigned int fnptdet;  // # trajectory points in any TPC
  unsigned int fnptcry;  // # trajectory points in any cryostat
  unsigned int fntpc;    // # TPC
  unsigned int fnapa;    // Total # APA
  unsigned int fnrop;    // Total # readout planes (ROPs)
  vector<int> fnpttpc;   // # trajectory point in each TPC
  vector<int> fnptapa;   // # trajectory point in each APA
  vector<int> fnptrop;   // # trajectory point in each ROP
  float fptx[maxpt];     // point x
  float fpty[maxpt];     // point y
  float fptz[maxpt];     // point z
  float fptt[maxpt];     // point t
  float fpte[maxpt];     // point E
  int fpttpc[maxpt];     // point TPC
  int fptapa[maxpt];     // point APA
  int fptuchan[maxpt];   // point channel in U-plane
  int fptvchan[maxpt];   // point channel in V-plane
  int fptzchan[maxpt];   // point channel in Z-plane
  float fptutick[maxpt]; // point tick in U-plane
  float fptvtick[maxpt]; // point tick in V-plane
  float fptztick[maxpt]; // point tick in Z-plane
  // Length of track in detector.
  float fdetlen;        // Length of track in detector.
  float fdettickmin;    // Smallest TDC tick in detector.
  float fdettickmax;    // Largest TDC tick in detector.
  float fdetx1;         // Detector entry x 
  float fdety1;         // Detector entry y 
  float fdetz1;         // Detector entry z 
  float fdetx2;         // Detector exit x 
  float fdety2;         // Detector exit y 
  float fdetz2;         // Detector exit z 
  // dE/dx bins on given track. 
  int fndedxBins;
  vector<double> fdedxBins;

  // Other variables that will be shared between different methods.
  double                            fElectronsToGeV; // conversion factor
  double fmcpdsmax;  // Maximum step size for filling the MC particle trajectory hists
  double fadcmevu;   // MeV to ADC conversion factor for U-planes.
  double fadcmevv;   // MeV to ADC conversion factor for V-planes.
  double fadcmevz;   // MeV to ADC conversion factor for X-planes.
  bool fhistusede;   // If true, raw and wire spectra are converted to MeV.
  double fdemaxmcp;  // Max energy deposit for histogram ranges for McParticle
  double fdemax;     // Max energy deposit for histogram ranges

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
  GeoHelper* fgeohelp;
  art::ServiceHandle<geo::Geometry> fGeometry;       // pointer to Geometry service

/*
  // Geometry information.
  vector<int> fropfirstchan;    // first channel for each readout plane
  vector<int> fropnchan;        // # channels for each readout plane
  vector<int> froptpc;          // First TPC for each readout plane
  vector<int> fropapa;          // the APA for each readout plane
  vector<string> fropname;      // Names for the TPC readout planes, e.g. rop1x1.
  vector<string> froporient;    // Wire orientation for each ROP plane: u, v or z.
  vector<int> fapanrop;         // # ROP for each APA.
  map<PlaneID, unsigned int> ftpcplanerop; // # ROP for each TPC plane.
*/

  // LArProperties service (for drift speed)
  art::ServiceHandle<util::DetectorProperties> fdetprop;
  double fdriftspeed;
  double fsamplingrate;
  double fsamplingdriftspeed;
  double fsamplingoffsetx;
  double fsamplingoffsetu;
  double fsamplingoffsetv;

}; // class LbTupler


//************************************************************************
// class implementation
//************************************************************************

//************************************************************************

// Constructor
LbTupler::LbTupler(fhicl::ParameterSet const& parameterSet)
: EDAnalyzer(parameterSet), fdbg(0), fgeohelp(0) {
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
  fntpc = geohelp.ntpc();
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
  if ( fdbg > 0 ) {
    cout << myname << "Detector length: " << detLength << " cm" << endl;
    cout << myname << "Detector width:  " << detWidth  << " cm" << endl;
    cout << myname << "Detector height: " << detHeight << " cm" << endl;
  }
  double detSize = sqrt( detLength*detLength + detWidth*detWidth + detHeight*detHeight );

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

  // Fetch LAr properties.
  art::ServiceHandle<util::LArProperties> larprop;
  fdriftspeed = larprop->DriftVelocity();
  fsamplingrate = fdetprop->SamplingRate();
  fsamplingdriftspeed = fdriftspeed*fsamplingrate/1000.0;
  if ( fdbg > 0 ) cout << myname << "LAr drift speed = " << fdriftspeed << " cm/us" << endl;
  if ( fdbg > 0 ) cout << myname << "Sampling period = " << fsamplingrate << " ns" << endl;
  if ( fdbg > 0 ) cout << myname << "Sampling speed = " << fsamplingdriftspeed << " cm/tick" << endl;

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
  fpdgCodeHist     = tfs->make<TH1D>("pdgcodes",";PDG Code;",                  5000, -2500, 2500);
  fMomentumHist    = tfs->make<TH1D>("mom",     ";particle Momentum (GeV);",    100, 0.,    10.);
  fTrackLengthHist = tfs->make<TH1D>("length",  ";particle track length (cm);", 200, 0, detLength);

  // Define our n-tuples, which are limited forms of ROOT
  // TTrees. Start with the TTree itself.
  fSimulationNtuple     = tfs->make<TTree>("LbTuplerSimulation",    "LbTuplerSimulation");
  fReconstructionNtuple = tfs->make<TTree>("LbTuplerReconstruction","LbTuplerReconstruction");
  fSimChannelNtuple     = tfs->make<TTree>("LbTuplerSimChannel","LbTuplerSimChannel");

  // Set array sizes.
  fnpttpc.resize(geohelp.ntpc());
  fnptapa.resize(geohelp.napa());
  fnptrop.resize(geohelp.nrop());

  // Define the branches (columns) of our simulation n-tuple. When
  // we write a variable, we give the address of the variable to
  // TTree::Branch.
  fSimulationNtuple->Branch("event",       &fevent,          "event/I");
  fSimulationNtuple->Branch("run",         &fRun,            "run/I");
  fSimulationNtuple->Branch("subrun",      &fSubRun,         "subrun/I");
  fSimulationNtuple->Branch("pdg",         &fpdg,            "pdg/I");              // PDG ID
  fSimulationNtuple->Branch("rpdg",        &frpdg,           "rpdg/I");             // reduced PDG ID
  fSimulationNtuple->Branch("proc",        &fproc,           "proc/I");             // reduced PDG ID
  fSimulationNtuple->Branch("trackid",     &ftrackid,        "trackid/I");          // Track ID
  fSimulationNtuple->Branch("parent",      &fparent,         "parent/I");           // Parent 
  fSimulationNtuple->Branch("nchild",      &fnchild,         "nchild/I");           // # children
  fSimulationNtuple->Branch("child",       fchild,           "child[nchild]/I");    // children
  fSimulationNtuple->Branch("ndetchild",   &fndetchild,      "ndetchild/I");        // # children in det
  fSimulationNtuple->Branch("detchild",    fdetchild,        "detchild[ndetchild]/I"); // children in det
  fSimulationNtuple->Branch("ndetin",      &fndetin,         "ndetin/I");           // # detector entries
  fSimulationNtuple->Branch("ndetout",     &fndetout,        "ndetout/I");          // # detector exits
  fSimulationNtuple->Branch("ntpcin",      &fntpcin,         "ntpcin/I");           // # detector entries
  fSimulationNtuple->Branch("ntpcout",     &fntpcout,        "ntpcout/I");          // # detector exits
  fSimulationNtuple->Branch("ncryin",      &fncryin,         "ncryin/I");           // # detector entries
  fSimulationNtuple->Branch("ncryout",     &fncryout,        "ncryout/I");          // # detector exits
  // When we write arrays, we give the address of the array to
  // TTree::Branch; in C++ this is simply the array name.
  fSimulationNtuple->Branch("StartXYZT",   fStartXYZT,       "StartXYZT[4]/F");     // Starting point
  fSimulationNtuple->Branch("EndXYZT",     fEndXYZT,         "EndXYZT[4]/F");       // Ending point
  fSimulationNtuple->Branch("StartPE",     fStartPE,         "StartPE[4]/F");       // Starting momentum
  fSimulationNtuple->Branch("EndPE",       fEndPE,           "EndPE[4]/F");         // Ending momentum
  // Trajectory points.
  fSimulationNtuple->Branch("npt",          &fnpt,          "npt/i");             // # points
  fSimulationNtuple->Branch("nptdet",       &fnptdet,       "nptdet/i");          // # points in detector
  fSimulationNtuple->Branch("nptcry",       &fnptcry,       "nptcry/i");          // # points in detector
  fSimulationNtuple->Branch("ntpc",         &fntpc,         "ntpc/i");            // # TPC
  fSimulationNtuple->Branch("npttpc",       fnpttpc.data(), "npttpc[ntpc]/i");    // # points in each TPC
  fSimulationNtuple->Branch("napa",         &fnapa,         "napa/i");            // # APA
  fSimulationNtuple->Branch("nptapa",       fnptapa.data(), "nptapa[napa]/i");    // # points in each APA
  fSimulationNtuple->Branch("nrop",         &fnrop,         "nrop/i");            // # ROP
  fSimulationNtuple->Branch("nptrop",       fnptrop.data(), "nptrop[nrop]/i");    // # points in each ROP
  fSimulationNtuple->Branch("ptx",          fptx,           "ptx[npt]/F");
  fSimulationNtuple->Branch("pty",          fpty,           "pty[npt]/F");
  fSimulationNtuple->Branch("ptz",          fptz,           "ptz[npt]/F");
  fSimulationNtuple->Branch("ptt",          fptt,           "ptt[npt]/F");
  fSimulationNtuple->Branch("pte",          fpte,           "pte[npt]/F");
  fSimulationNtuple->Branch("pttpc",        fpttpc,         "pttpc[npt]/I");
  fSimulationNtuple->Branch("ptapa",        fptapa,         "ptapa[npt]/I");
  fSimulationNtuple->Branch("ptuchan",      fptuchan,       "ptuchan[npt]/I");
  fSimulationNtuple->Branch("ptvchan",      fptvchan,       "ptvchan[npt]/I");
  fSimulationNtuple->Branch("ptzchan",      fptzchan,       "ptzchan[npt]/I");
  fSimulationNtuple->Branch("ptutick",      fptutick,       "ptutick[npt]/F");
  fSimulationNtuple->Branch("ptvtick",      fptvtick,       "ptvtick[npt]/F");
  fSimulationNtuple->Branch("ptztick",      fptztick,       "ptztick[npt]/F");
  // length of track in detector
  fSimulationNtuple->Branch("detlen",       &fdetlen,       "detlen/F");
  fSimulationNtuple->Branch("dettickmin",   &fdettickmin,   "dettickmin/F");
  fSimulationNtuple->Branch("dettickmax",   &fdettickmax,   "dettickmax/F");
  fSimulationNtuple->Branch("detx1",        &fdetx1,        "detx1/F");
  fSimulationNtuple->Branch("dety1",        &fdety1,        "dety1/F");
  fSimulationNtuple->Branch("detz1",        &fdetz1,        "detz1/F");
  fSimulationNtuple->Branch("detx2",        &fdetx2,        "detx2/F");
  fSimulationNtuple->Branch("dety2",        &fdety2,        "dety2/F");
  fSimulationNtuple->Branch("detz2",        &fdetz2,        "detz2/F");
  // dE/dx bins
  fSimulationNtuple->Branch("ndedx",        &fndedxBins,      "ndedx/I");
  fSimulationNtuple->Branch("dedx",         fdedxBins.data(), "dedx[ndedx]/D");

  // A similar definition for the reconstruction n-tuple. Note that we
  // use some of the same variables in both n-tuples.
  fReconstructionNtuple->Branch("Event",   &fevent,          "Event/I");
  fReconstructionNtuple->Branch("SubRun",  &fSubRun,         "SubRun/I");
  fReconstructionNtuple->Branch("Run",     &fRun,            "Run/I");
  fReconstructionNtuple->Branch("TrackID", &ftrackid,        "TrackID/I");
  fReconstructionNtuple->Branch("PDG",     &fpdg,            "PDG/I");
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
  fDoMCTrackPerf           = p.get<bool>("DoMCTrackPerf");
  fDoSimChannels           = p.get<bool>("DoSimChannels");
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
    cout << prefix << setw(wlab) << "DoMCTrackPerf" << sep << fDoMCTrackPerf << endl;
    cout << prefix << setw(wlab) << "DoSimChannels" << sep << fDoSimChannels << endl;
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
  fgeohelp = new GeoHelper(&*fGeometry, &*fdetprop, 1);
  fgeohelp->print(cout, 0, myname);
  fnrop = fgeohelp->nrop();
  return;
}

//************************************************************************

void LbTupler::analyze(const art::Event& event) {
  const string myname = "LbTupler::analyze: ";

  // Access ART's TFileService, which will handle creating and writing
  // histograms and n-tuples for us. 
  art::ServiceHandle<art::TFileService> tfs;

  // Check gemetry helper.
  if ( fgeohelp == nullptr ) {
    cout << myname << "ERROR: Geometry helper is absent." << endl;
    return;
  }
  const GeoHelper& geohelp = *fgeohelp;

  // Start by fetching some basic event information for our n-tuple.
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
  if ( fDoMCParticles || fDoMCTrackPerf ) {
    event.getByLabel(fParticleProducerLabel, particleHandle);
    if ( fdbg > 1 ) cout << myname << "MCParticle count: " << particleHandle->size() << endl;
  }

  // Get all the simulated channels for the event. These channels
  // include the energy deposited for each track.
  art::Handle<vector<sim::SimChannel>> simChannelHandle;
  if ( fDoMCParticles || fDoSimChannels || fDoMCTrackPerf ) {
    event.getByLabel(fSimulationProducerLabel, simChannelHandle);
    if ( fdbg > 1 ) cout << myname << "Sim channel count: " << simChannelHandle->size() << endl;
  }

  // Create vector of selected MC particles for analysis.
  vector<MCTrackPerf> selectedMCTrackPerfsMC;    // Filled with MCParticle hits.
  vector<MCTrackPerf> selectedMCTrackPerfsSC;    // Filled with SimChannel.
  if ( fDoMCTrackPerf ) {
    for ( auto const& particle : (*particleHandle) ) {
      int trackid = particle.TrackId();
      int rpdg = reducedPDG(particle.PdgCode());
      int proc = intProcess(particle.Process());
      // Select particles.
      // 21apr2015: Keep also gammas
      if ( proc == 0 && rpdg < 7 ) {
        if ( fdbg > 1 ) cout << myname << "Selecting";
        selectedMCTrackPerfsMC.emplace(selectedMCTrackPerfsMC.end(), *&particle, fgeohelp);
        selectedMCTrackPerfsSC.emplace(selectedMCTrackPerfsSC.end(), *&particle, fgeohelp);
      } else {
        if ( fdbg > 1 ) cout << myname << "Rejecting";
      }
      if ( fdbg > 1 ) cout << " MC particle " << setw(6) << trackid 
                          << " with RPDG=" << setw(2) << rpdg
                          << " and PROC=" << setw(2) << proc
                          << " and PDG=" << setw(10) << particle.PdgCode()
                          << endl;
    }
  }  // end DoMCTrackPerf

  // MCParticles histograms and tree.
  if ( fDoMCParticles ) {

    // Create the event channel-tick bin histograms for all selected MC particles.
    int ntick = ftdcTickMax - ftdcTickMin;
    string hnprefix = "h" + sevtf + "_mcp";
    string htprefix = "Event " + sevt + " MC particles for ";
    string htsuffix = ";TDC tick;Channel;Energy [MeV]";
    vector<TH2*> mcphists;
    for ( unsigned int irop=0; irop<geohelp.nrop(); ++irop ) {
      int nchan = geohelp.ropNChannel(irop);
      string hname = hnprefix + geohelp.ropName(irop);
      string title = htprefix + geohelp.ropName(irop) + htsuffix;
      if ( fdbg > 1 ) cout << myname << "Creating sim histo " << hname << " with " << ntick
                          << " TDC bins " << " and " << nchan << " channel bins" << endl;
      TH2* ph = tfs->make<TH2F>(hname.c_str(), title.c_str(),
                                ntick, ftdcTickMin, ftdcTickMax,
                                nchan, 0, nchan);
      ph->GetZaxis()->SetRangeUser(0.0, fdemaxmcp);
      ph->SetContour(20);
      ph->SetStats(0);
      mcphists.push_back(ph);
    }

    // Loop over particles and fetch the # points inside the detector for each.
    // Might later want to add the # descendants with points inside the detector.
    map<unsigned int, unsigned int> ndetptmap;
    for ( auto const& particle : (*particleHandle) ) {
      unsigned int npt = 0;
      unsigned int tid = particle.TrackId();
      size_t numberTrajectoryPoints = particle.NumberTrajectoryPoints();
      for ( unsigned int ipt=0; ipt<numberTrajectoryPoints; ++ipt ) {
        const auto& pos = particle.Position(ipt);
        double xyzt[4] = {pos.X(), pos.Y(), pos.Z(), pos.T()};
        geo::TPCID tpcid = fGeometry->FindTPCAtPosition(xyzt);
        if ( tpcid.isValid ) ++npt;
      }
      ndetptmap[tid] = npt;
    }

    // Loop over particles.
    // See ${NUTOOLS_DIR}/include/SimulationBase/MCParticle.h.
    for ( auto const& particle : (*particleHandle) ) {
      ftrackid = particle.TrackId();
      fpdg = particle.PdgCode();
      frpdg = reducedPDG(fpdg);
      fparent = particle.Mother();
      fproc = intProcess(particle.Process());
      if ( fproc < 0 ) {
        cout << myname << "WARNING: Unknown process: " << particle.Process() << endl;
      }
      fnchild = particle.NumberDaughters();
      if ( fnchild > fmaxchild ) {
        cout << myname << "WARNING: Too many child particles: " << fnchild << endl;
        fnchild = fmaxchild;
      }
      fndetchild = 0;
      fndetin = 0;
      fndetout = 0;
      fntpcin = 0;
      fntpcout = 0;
      fncryin = 0;
      fncryout = 0;
      for ( unsigned int ichi=0; ichi<fnchild; ++ichi ) {
        unsigned int tid = particle.Daughter(ichi);
        fchild[ichi] = tid;
        if ( ndetptmap[tid] ) fdetchild[fndetchild++] = tid;
      }
      if ( fdbg > 2 ) {
        cout << myname << "ID=" << ftrackid << ", PDG=" << fpdg
             << ", RPDG=" << frpdg
             << ", status=" << particle.StatusCode()
             << ", Process: " << particle.Process() << "(" << fproc << ")"
             << ", Parent: " << fparent
             << endl;
      }

      // Find matching selected track.
      MCTrackPerf* pmctp = nullptr;
      for ( auto& mctp : selectedMCTrackPerfsMC ) {
        if ( int(mctp.trackID()) == ftrackid ) {
          pmctp = &mctp;
          if ( fdbg > 1 ) cout << myname << "  This is a selected track." << endl;
          break;
        }
      }

      // Fill histograms.
      fpdgCodeHist->Fill( fpdg );

      // For this example, we want to fill the n-tuples and histograms
      // only with information from the primary particles in the
      // event, whose PDG codes match a value supplied in the .fcl file.
      if ( true || ((particle.Process() == "primary"  &&  fpdg == fSelectedPDG)) ) {
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
        fnpt = 0;
        fnptdet = 0;
        fnptcry = 0;
        for ( auto& count : fnpttpc ) count = 0;
        for ( auto& count : fnptapa ) count = 0;
        for ( auto& count : fnptrop ) count = 0;
      
        double x0 = 0.0;
        double y0 = 0.0;
        double z0 = 0.0;
        double t0 = 0.0;
        double e0 = 0.0;
        bool indet0 = false;
        const unsigned int notcry = UINT_MAX;
        unsigned int icry0 = notcry;
        geo::TPCID tpcid0;
        if ( tpcid0.isValid ) {
          cout << myname << "ERROR: Initial TPCID is valid." << endl;
          abort();
        }
        fdetlen = 0.0;
        fdettickmin =  1000000.0;
        fdettickmax = -1000000.0;
        fdetx1 = 1.e6;
        fdety1 = 1.e6;
        fdetz1 = 1.e6;
        fdetx2 = 1.e6;
        fdety2 = 1.e6;
        fdetz2 = 1.e6;
        for ( unsigned int ipt=0; ipt<numberTrajectoryPoints; ++ipt ) {
          if ( ipt >= maxpt ) {
            cout << myname << "Found more than " << maxpt << " trajectory points."
                 << " The remainder will be skipped." << endl;
            break;
          }
          bool lastpoint = (ipt+1 == numberTrajectoryPoints);
          const auto& pos = particle.Position(ipt);
          const auto& mom = particle.Momentum(ipt);
          fptx[ipt] = pos.X();
          fpty[ipt] = pos.Y();
          fptz[ipt] = pos.Z();
          fptt[ipt] = pos.T();
          fpte[ipt] = mom.E();
          double xyzt[4] = {pos.X(), pos.Y(), pos.Z(), pos.T()};
          geo::TPCID tpcid = fGeometry->FindTPCAtPosition(xyzt);
          unsigned int icry = fGeometry->FindCryostatAtPosition(xyzt);
          if ( icry != notcry ) ++fnptcry;
          fptuchan[ipt] = -1;
          fptvchan[ipt] = -1;
          fptzchan[ipt] = -1;
          fptutick[ipt] = -1.0;
          fptvtick[ipt] = -1.0;
          fptztick[ipt] = -1.0;
          bool indet = false;
          ++fnpt;
          double x = pos.X();
          double y = pos.Y();
          double z = pos.Z();
          double t = pos.T();
          double e = particle.E(ipt);
          double detot = 1000.0*(e0 - e);
          if ( ! tpcid.isValid ) {
            fpttpc[ipt] = -1;
            fptapa[ipt] = -1;
          } else if ( tpcid.Cryostat != 0 ) {
            fpttpc[ipt] = -2;
            fptapa[ipt] = -2;
          } else {
            indet = true;
            unsigned int itpc = tpcid.TPC;
            unsigned int iapa = geohelp.tpcApa(itpc);
            fpttpc[ipt] = itpc;
            fptapa[ipt] = iapa;
            if ( fnptdet == 0 ) {
              fdetx1 = x;
              fdety1 = y;
              fdetz1 = z;
            }
            fdetx2 = x;
            fdety2 = y;
            fdetz2 = z;
            // Find the (channel, tick) for each plane in this TPC.
            PlanePositionVector pps = geohelp.planePositions(xyzt);
            if ( pps.size() ) {
              ++fnptdet;
              ++fnpttpc[itpc];
              ++fnptapa[iapa];
            }
            for ( const auto& pp : pps ) {
              unsigned int irop = pp.rop;
              View_t orient = geohelp.ropView(irop);
              ++fnptrop[irop];
              if ( orient == kU ) {
                fptuchan[ipt] = pp.ropchannel;
                fptutick[ipt] = pp.tick;
              } else if ( orient == kV ) {
                fptvchan[ipt] = pp.ropchannel;
                fptvtick[ipt] = pp.tick;
              } else if ( orient == kZ ) {
                fptzchan[ipt] = pp.ropchannel;
                fptztick[ipt] = pp.tick;
              }
            }
          }
          if ( fdbg > 3 ) cout << myname << "  MC Particle " << ftrackid
                               << " point " << fnpt << ": xyzt=("
                               << x << ", " << y << ", " << z << ", " << t << ")"
                               << ", E=" << 1000*e << " MeV"
                               << ", TPC=" << tpcid.TPC
                               << endl;
          if ( ipt ) {
            if ( !indet0 && indet ) ++fndetin;
            if ( indet0 && !indet ) ++fndetout;
            if ( tpcid.isValid && (!tpcid0.isValid or tpcid.TPC != tpcid0.TPC ) ) {
              ++fntpcin;
              if ( fdbg > 3 ) {
                cout << myname << "    Entering TPC " << tpcid.TPC;
                if ( tpcid0.isValid ) cout << ", exiting TPC " << tpcid0.TPC;
                cout << endl;
              }
            }
            if ( tpcid0.isValid && (!tpcid.isValid or tpcid.TPC != tpcid0.TPC ) ) {
              ++fntpcout;
              if ( fdbg > 3 ) {
                cout << myname << "    Exiting TPC " << tpcid0.TPC;
                if ( tpcid.isValid ) cout << ", entering TPC " << tpcid.TPC;
                cout << endl;
              }
            }
            if ( icry!=notcry && icry!=icry0 ) {
              ++fncryin;
              if ( fdbg > 3 ) {
                cout << myname << "    Entering cryostat " << icry;
                if ( icry0 != notcry ) cout << ", exiting cryostat " << icry0;
                cout << endl;
              }
            }
            if ( icry0!=notcry && icry!=icry0 ) {
              ++fncryout;
              if ( fdbg > 3 ) {
                cout << myname << "    Exiting cryostat " << icry0;
                if ( icry != notcry ) cout << ", entering cryostat " << icry;
                cout << endl;
              }
            }
            if ( lastpoint && fdbg > 3 ) {
              cout << myname << "    Last trajectory point; # children is " << fnchild << endl; 
            }
          }
          // Fill time-channel histogram with dE for each pair of adjacent tracjectory points.
          // Increase the number of points to ensure granularity less than fmcpdsmax;
          if ( ipt != 0 && indet0 && indet) {
            double dx = x - x0;
            double dy = y - y0;
            double dz = z - z0;
            double dt = t - t0;
            double ds = sqrt(dx*dx + dy*dy + dz*dz);
            unsigned int nstep = ds/fmcpdsmax + 1;
            double invstep = 1.0/nstep;
            double destep = invstep*detot;
            if ( fdbg > 3 ) cout << myname << "    # steps: " << nstep << endl;
            for ( unsigned int istp=0; istp<nstep; ++istp ) {
              double xa = x0 + (istp+0.5)*invstep*dx;
              double ya = y0 + (istp+0.5)*invstep*dy;
              double za = z0 + (istp+0.5)*invstep*dz;
              double ta = t0 + (istp+0.5)*invstep*dt;
              double postim[4] = {xa, ya, za, ta};
              PlanePositionVector pps = geohelp.planePositions(postim);
              for ( const auto& pp : pps ) {
                if ( ! pp.valid ) cout << myname << "    Invalid plane position!" << endl;
                TH2* ph = mcphists[pp.rop];
                if ( fdbg > 3 ) {
                  cout << myname << "    Filling " << ph->GetName()
                       << ": tick=" << pp.tick
                       << ", chan=" << pp.ropchannel
                       << ", DE=" << destep << " MeV" << endl;
                }
                ph->Fill(pp.tick, pp.ropchannel, destep);
                if ( pp.tick < fdettickmin ) fdettickmin = pp.tick;
                if ( pp.tick > fdettickmax ) fdettickmax = pp.tick;
                if ( pmctp != nullptr ) pmctp->addSignal(pp.channel, pp.tick, destep);
              }
            }
            // If this and the last point are in the detector, increment the detector path length.
            if ( indet0 && indet ) fdetlen += ds;
          }
          x0 = x,
          y0 = y;
          z0 = z;
          t0 = t;
          e0 = e;
          indet0 = indet;
          tpcid0 = tpcid;
          icry0 = icry;
        }
        if ( fnptdet != ndetptmap[ftrackid] ) {
          cout << myname << "WARNING: Inconsistent # detector points: "
               << fnptdet << " != " << ndetptmap[ftrackid] << endl;
        }

        // Use a polar-coordinate view of the 4-vectors to
        // get the track length.
        double trackLength = ( positionEnd - positionStart ).Rho();
        if ( fdbg > 2 ) cout << myname << "  Track length: " << trackLength << " mm" << endl;

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
                  if ( fdbg > 4 ) cout << myname << "  Nel, dE: " << nel << ", " << de << endl;
                }
              } // For each energy deposit
            } // For each time slice
          } // If channel is in a collection plane
        } // For each SimChannel

        // At this point we've filled in the values of all the
        // variables and arrays we want to write to the
        // n-tuple. The following command actually writes those
        // values.
      } // if selected
      fSimulationNtuple->Fill();
    } // end loop over all particles in the event. 
    if ( fdbg > 2 ) {
      cout << myname << "Tree " << fSimulationNtuple->GetName()
                     << " entry count  is " << fSimulationNtuple->GetEntries() << endl;
      for ( auto ph : mcphists ) {
        cout << myname << "Histogram " << setw(20) << ph->GetName()
             << " entry count is " << ph->GetEntries() << endl;
      }
    }

    // Create the channel-tick bin histograms for each selected MC particle.
    vector<vector<TH2*>> mcphists_mcp;
    for ( auto& mctrackperf : selectedMCTrackPerfsMC ) {
      vector<TH2*>& hists = *mcphists_mcp.emplace(mcphists_mcp.end(), geohelp.nrop(), nullptr);
      ostringstream ssmcp;
      ssmcp << mctrackperf.trackID();
      string smcp = ssmcp.str();
      hnprefix = "h" + sevtf + "mcp" + smcp + "_mcp";
      htprefix = "Event " + sevt + " MC particle " + smcp + " for ";
      for ( unsigned int irop=0; irop<geohelp.nrop(); ++irop ) {
        if ( mctrackperf.hits().ropNbin(irop) == 0 ) continue;
        int nchan = geohelp.ropNChannel(irop);
        string hname = hnprefix + geohelp.ropName(irop);
        string title = htprefix + geohelp.ropName(irop) + htsuffix;
        TH2* ph = tfs->make<TH2F>(hname.c_str(), title.c_str(),
                                  ntick, ftdcTickMin, ftdcTickMax,
                                  nchan, 0, nchan);
        ph->GetZaxis()->SetRangeUser(0.0, fdemaxmcp);
        ph->SetContour(20);
        ph->SetStats(0);
        mctrackperf.fillRopChannelTickHist(ph, irop);
        hists[irop] = ph;
      }  // End loop over ROPs
    }  // End loop over selected MC particles

  }  // end Do MCParticle

  if ( fDoMCTrackPerf ) {
    int flag = 0;
    if ( fdbg > 0 ) flag = 11;
    if ( fdbg > 0 ) cout << myname << "Selected MC track list with MCParticle fill (size = "
                         << selectedMCTrackPerfsMC.size() << ")" << endl;
    for ( auto& mctrackperf : selectedMCTrackPerfsMC ) {
      mctrackperf.print(cout, flag, myname + "  ");
    }  // End loop over selected MC tracks
  }  // end DoMCTrackPerf

  //************************************************************************
  // Sim channels.
  //************************************************************************

  if ( fDoSimChannels ) {

    // Check array sizes.
    if ( simChannelHandle->size() > fscCapacity ) {
      cout << myname << "WARNING: Sim channel count exceeds TTree capacity." << endl;
    }

    // Create the Sim channel histograms.
    vector<TH2*> sphists;
    for ( unsigned int irop=0; irop<geohelp.nrop(); ++irop ) {
      int nchan = geohelp.ropNChannel(irop);
      int ntick = ftdcTickMax-ftdcTickMin;
      string hname = "h" + sevtf + "sim" + geohelp.ropName(irop);
      string title = "Sim channels for TPC plane " + geohelp.ropName(irop) + " event " + sevt + ";TDC tick;Channel;Energy [MeV]";
      if ( fdbg > 1 ) cout << myname << "Creating sim histo " << hname << " with " << ntick
                          << " TDC bins " << " and " << nchan << " channel bins" << endl;
      TH2* ph = tfs->make<TH2F>(hname.c_str(), title.c_str(),
                                ntick, ftdcTickMin, ftdcTickMax,
                                nchan, 0, nchan);
      ph->GetZaxis()->SetRangeUser(0.0, fdemax);
      ph->SetContour(20);
      ph->SetStats(0);
      sphists.push_back(ph);
    }

    // Loop over the SimChannel objects in the event.
    fscCount = 0;
    for ( auto& val : fscChannel ) val = 0;
    for ( auto& val : fscCharge ) val = 0.0;
    for ( auto& val : fscEnergy ) val = 0.0;
    for ( auto const& simchan : (*simChannelHandle) ) {
      auto ichan = simchan.Channel();
      unsigned int irop = geohelp.channelRop(ichan);
      if ( irop == geohelp.nrop() ) {
        cout << myname << "ERROR: SimChannel channel " << ichan << " is not in a readout plane." << endl;
        abort();
      }
      int iropchan = ichan - int(geohelp.ropFirstChannel(irop));
      if ( iropchan < 0 || iropchan >= int(geohelp.ropNChannel(irop)) ) {
        cout << myname << "ERROR: ROP channel " << iropchan << " is out of range [0, "
             << geohelp.ropNChannel(irop) << ")" << endl;
        abort();
      }
      // Fetch the histogram for the current ROP.
      TH2* psphist = sphists[irop];
      if ( fdbg > 3 ) cout << myname << "SimChannel " << setw(4) << fscCount
                          << ": ROP " << setw(2) << irop
                          << ", Global/ROP channel " << setw(4) << ichan
                          << "/" << setw(3) << iropchan;
      if ( fscCount < fscCapacity ) {
        fscChannel[fscCount] = ichan;
        auto const& idemap = simchan.TDCIDEMap();
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
          psphist->Fill(tdc, iropchan, tdcEnergy);
          if ( fdbg > 4 ) cout << "\n" << myname << "    " << psphist->GetName()
                              << ": Global/ROP channel " << setw(4) << ichan
                              << "/" << setw(3) << iropchan
                              << ", tdc " << tdc
                              << " has " << int(100.0*tdcEnergy+0.4999)/100.0 << " MeV";
        }
        if ( fdbg > 3 ) {
          if ( fdbg > 4 ) cout << "\n" << myname << "    Total ";
          else cout << ": ";
          double ren = 0.01*int(100.0*energy+0.49999);
          double rchg = int(charge+0.499999);
          cout << "E=" << ren << " MeV; Q=" << rchg << endl;
        }
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

    // Fill SimChannel tree.
    fSimChannelNtuple->Fill();


  }  // end DoSimChannel

  // Add sim channel info to the sim channel performance objects.
  if ( fDoMCTrackPerf ) {
    if ( fdbg > 0 ) cout << myname << "Selected MC track list with SimChannel fill (size = "
                        << selectedMCTrackPerfsSC.size() << ")" << endl;
    for ( auto& mctrackperf : selectedMCTrackPerfsSC ) {
      // Add the sim channel info to the selected tracks.
      for ( auto const& simchan : (*simChannelHandle) ) {
        mctrackperf.addSimChannel(*&simchan);
      }  // End loop over sim channels in the event. 
    }  // End loop over selected SimChannel MC tracks
  }  // end DoMCTrackPerf

  // Add hits to the performance objects.
  if ( fDoMCTrackPerf ) {
    for ( auto& mctrackperf : selectedMCTrackPerfsSC ) mctrackperf.buildHits();
    for ( auto& mctrackperf : selectedMCTrackPerfsMC ) mctrackperf.buildHits();
  }

  // Display the MC tracks.
  if ( fDoMCTrackPerf && fdbg > 3 ) {
    for ( unsigned int imcs=0; imcs<selectedMCTrackPerfsSC.size(); ++imcs ) {
      const auto mctsc = selectedMCTrackPerfsSC.at(imcs);
      const auto mctmc = selectedMCTrackPerfsMC.at(imcs);
      if ( fdbg > 3 ) {
        cout << myname << "Dumping MCTrackPerf sim channels for event " << fevent
             << " track " << mctsc.trackID()
             << "\n" << myname << "  Total tick/hit signal: "
             << mctsc.hits().tickSignal() << "/" << mctsc.hits().hitSignal() << " MeV:" << endl;
        mctsc.print(cout,  0, myname + "  ");
        mctsc.print(cout, 12, myname + "  ");
      }
      if ( fdbg > 3 ) {
        cout << myname << "Dumping MCTrackPerf MC hits for event " << fevent
             << " track " << mctmc.trackID()
             << "\n" << myname << "  Total tick/hit signal: "
             << mctmc.hits().tickSignal() << "/" << mctmc.hits().hitSignal() << " MeV:" << endl;
        mctmc.print(cout, 11, myname + "  ");
      }
    }  // End loop over selected MC tracks
  }  // end DoMCTrackPerf

  //************************************************************************
  // Raw digits.
  //************************************************************************

  // Assign histogram titles.
  string ztitle = "ADC counts";
  double zmax = 150;
  int ncontour = 30;
  if ( fhistusede ) {
    ztitle = "Energy [MeV]";
    zmax = fdemax;
    ncontour = 40;
  }

  if ( fDoRaw ) {
    // Get the raw digits for the event.
    art::Handle< vector<raw::RawDigit> > rawDigitHandle;
    event.getByLabel(fRawDigitProducerLabel, rawDigitHandle);
    if ( fdbg > 1 ) cout << myname << "Raw digit count: " << rawDigitHandle->size() << endl;

    // Create the Raw digit histograms.
    vector<TH2*> rawhists;
    for ( unsigned int irop=0; irop<geohelp.nrop(); ++irop ) {
      int nchan = geohelp.ropNChannel(irop);
      int ntick = ftdcTickMax-ftdcTickMin;
      string hname = "h" + sevtf + "raw" + geohelp.ropName(irop);
      string title = "Raw digits for TPC plane " + geohelp.ropName(irop) + " event " + sevt
                     + ";TDC tick;Channel;" + ztitle;
      if ( fdbg > 1 ) cout << myname << "Creating sim histo " << hname << " with " << ntick
                          << " TDC bins " << " and " << nchan << " channel bins" << endl;
      TH2* ph = tfs->make<TH2F>(hname.c_str(), title.c_str(),
                                ntick, ftdcTickMin, ftdcTickMax,
                                nchan, 0, nchan);
      ph->GetZaxis()->SetRangeUser(-zmax, zmax);
      ph->SetContour(ncontour);
      ph->SetStats(0);
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
  }  // end DoRawDigit

  //************************************************************************
  // Wires.
  //************************************************************************

  if ( fDoWires ) {
    // See $LARDATA_DIR/include/RecoBase/Wire.h
    art::Handle< vector<recob::Wire> > wiresHandle;
    event.getByLabel(fWireProducerLabel, wiresHandle);
    if ( fdbg > 1 ) cout << myname << "Wire count: " << wiresHandle->size() << endl;

    // Create the wire histograms.
    vector<TH2*> wirehists;
    for ( unsigned int irop=0; irop<geohelp.nrop(); ++irop ) {
      int nchan = geohelp.ropNChannel(irop);
      int ntick = ftdcTickMax-ftdcTickMin;
      string hname = "h" + sevtf + "wir" + geohelp.ropName(irop);
      string title = "Wire signals for TPC plane " + geohelp.ropName(irop) + " event " + sevt
                     + ";TDC tick;Channel;" + ztitle;
      if ( fdbg > 1 ) cout << myname << "Creating wire histo " << hname << " with " << ntick
                          << " TDC bins " << " and " << nchan << " channel bins" << endl;
      TH2* ph = tfs->make<TH2F>(hname.c_str(), title.c_str(),
                                ntick, ftdcTickMin, ftdcTickMax,
                                nchan, 0, nchan);
      ph->GetZaxis()->SetRangeUser(-zmax, zmax);
      ph->SetContour(ncontour);
      ph->SetStats(0);
      wirehists.push_back(ph);
    }

    for ( auto const& wire : (*wiresHandle) ) {
      int ichan = wire.Channel();
      unsigned int irop = geohelp.channelRop(ichan);
      unsigned int iropchan = ichan - geohelp.ropFirstChannel(irop);
      auto sigs = wire.Signal();
      const auto& roisigs = wire.SignalROI();
      TH2* ph = wirehists[irop];
      if ( fdbg > 3 ) cout << myname << "Wire channel " << ichan
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
      int nchan = geohelp.ropNChannel(irop);
      int ntick = ftdcTickMax-ftdcTickMin;
      string hname = "h" + sevtf + "hit" + geohelp.ropName(irop);
      string title = "Hits for TPC plane " + geohelp.ropName(irop) + " event " + sevt
                     + ";TDC tick;Channel;" + ztitle;
      if ( fdbg > 1 ) cout << myname << "Creating hit histo " << hname << " with " << ntick
                          << " TDC bins " << " and " << nchan << " channel bins" << endl;
      TH2* ph = tfs->make<TH2F>(hname.c_str(), title.c_str(),
                                ntick, ftdcTickMin, ftdcTickMax,
                                nchan, 0, nchan);
      ph->GetZaxis()->SetRangeUser(-zmax, zmax);
      ph->SetContour(ncontour);
      ph->SetStats(0);
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
      ph->Fill(hit.PeakTime(), hit.Channel(), wt);
    }
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
    // Loop over clusters.
    int wclus = 3;
    int wtpc = 1;
    int wapa = 1;
    int wchan = 4;
    int whit = 3;
    int wtick = 5;
    int wadc = 4;
    int wsadc = 6;
    int wwid = 18;
    for ( unsigned int iclu=0; iclu<clusters.size(); ++iclu ) {
      art::Ptr<recob::Cluster> pclu = clusters[iclu];
      std::vector<art::Ptr<recob::Hit>> hits = clusterHits.at(iclu);
      if ( fdbg > 2 ) {
        int cltpc = pclu->Plane().TPC;
        int clapa = geohelp.tpcApa(cltpc);
        float cladc = pclu->SummedADC();
        float clchg = pclu->Integral();
        cout << myname << "  Cluster " << setw(wclus) << iclu
             << " TPC: " << setw(wtpc) << cltpc
             << " APA: " << setw(wapa) << clapa
             << " hit count: " << setw(whit) << hits.size()
             << ", ADC=" << setw(wsadc) << int(cladc)
             << ", Q=" << clchg
             << endl;
        }
      // Loop over hits on the cluster.
      for ( auto& phit : hits ) {
        if ( phit.isNull() ) {
          cout << myname << "ERROR: Cluster hit is mising." << endl;
          abort();
        }
        WireID wid = phit->WireID();
        ostringstream sswid;
        sswid << wid;
        unsigned int chan = phit->Channel();
        int lid = phit->LocalIndex();
        unsigned int tick1 = phit->StartTick();
        unsigned int tick2 = phit->EndTick();
        if ( fdbg > 3 ) cout << myname << "    "
                            << setw(wchan) << chan << "-" << lid
                            << " " << setw(wwid) << sswid.str() << " I:" << lid
                            << "  (" << setw(wtick) << tick1 << ", " << setw(wtick) << tick2 << ")"
                            << ", ADC=" << setw(wadc) << int(phit->SummedADC())
                            << endl;
      }  // end loop over cluster hits
      ++iclu;
    }  // end loop over event clusters
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

  // This macro has to be defined for this module to be invoked from a
  // .fcl file; see LbTupler.fcl for more information.
  DEFINE_ART_MODULE(LbTupler)

//************************************************************************

} // namespace LbTupler

#endif // LbTupler_Module
