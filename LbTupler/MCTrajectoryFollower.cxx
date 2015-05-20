// MCTrajectoryFollower_module.cxx

#include "MCTrajectoryFollower.h"

#include <iostream>
#include <iomanip>
#include <sstream>

// Art includes
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"

// LArSoft includes
#include "SimulationBase/MCParticle.h"
#include "Geometry/Geometry.h"

// ROOT includes.
#include "TTree.h"

// Local includes.
#include "intProcess.h"
#include "reducedPDG.h"
#include "GeoHelper.h"
#include "MCTrackPerf.h"

using std::cout;
using std::endl;
using std::string;
using std::ostringstream;
using simb::MCParticle;
using geo::View_t;
using geo::kU;
using geo::kV;
using geo::kZ;

//************************************************************************

// Constructor
MCTrajectoryFollower::MCTrajectoryFollower(double dsmax, bool filltree, const GeoHelper* pgeohelp, int dbg)
: m_dbg(dbg), m_filltree(filltree), m_dsmax(dsmax), m_geohelp(pgeohelp)  {

  const string myname = "MCTrajectory:ctor: ";

  if ( m_geohelp == nullptr ) {
    cout << myname << "ERROR: Geometry helper is absent." << endl;
    return;
  }
  const GeoHelper& geohelp = *m_geohelp;

  // Access ART's TFileService, which will handle creating and writing
  // histograms and n-tuples for us. 
  art::ServiceHandle<art::TFileService> tfs;

  // Define the tree.
  m_ptree = tfs->make<TTree>("LbTuplerSimulation",    "LbTuplerSimulation");

  // Set array sizes.
  fnpttpc.resize(geohelp.ntpc());
  fnptapa.resize(geohelp.napa());
  fnptrop.resize(geohelp.nrop());

  // Define the branches (columns) of the tree.
  if ( m_filltree ) {
    m_ptree->Branch("event",       &fevent,          "event/I");
    m_ptree->Branch("run",         &fRun,            "run/I");
    m_ptree->Branch("subrun",      &fSubRun,         "subrun/I");
    m_ptree->Branch("pdg",         &fpdg,            "pdg/I");              // PDG ID
    m_ptree->Branch("rpdg",        &frpdg,           "rpdg/I");             // reduced PDG ID
    m_ptree->Branch("proc",        &fproc,           "proc/I");             // reduced PDG ID
    m_ptree->Branch("trackid",     &ftrackid,        "trackid/I");          // Track ID
    m_ptree->Branch("parent",      &fparent,         "parent/I");           // Parent 
    m_ptree->Branch("nchild",      &fnchild,         "nchild/I");           // # children
    m_ptree->Branch("child",       fchild,           "child[nchild]/I");    // children
    m_ptree->Branch("ndetchild",   &fndetchild,      "ndetchild/I");        // # children in det
    m_ptree->Branch("detchild",    fdetchild,        "detchild[ndetchild]/I"); // children in det
    m_ptree->Branch("ndetin",      &fndetin,         "ndetin/I");           // # detector entries
    m_ptree->Branch("ndetout",     &fndetout,        "ndetout/I");          // # detector exits
    m_ptree->Branch("ntpcin",      &fntpcin,         "ntpcin/I");           // # detector entries
    m_ptree->Branch("ntpcout",     &fntpcout,        "ntpcout/I");          // # detector exits
    m_ptree->Branch("ncryin",      &fncryin,         "ncryin/I");           // # detector entries
    m_ptree->Branch("ncryout",     &fncryout,        "ncryout/I");          // # detector exits
    m_ptree->Branch("StartXYZT",   fStartXYZT,       "StartXYZT[4]/F");     // Starting point
    m_ptree->Branch("EndXYZT",     fEndXYZT,         "EndXYZT[4]/F");       // Ending point
    m_ptree->Branch("StartPE",     fStartPE,         "StartPE[4]/F");       // Starting momentum
    m_ptree->Branch("EndPE",       fEndPE,           "EndPE[4]/F");         // Ending momentum
    // Trajectory points.
    m_ptree->Branch("npt",          &fnpt,          "npt/i");             // # points
    m_ptree->Branch("nptdet",       &fnptdet,       "nptdet/i");          // # points in detector
    m_ptree->Branch("nptcry",       &fnptcry,       "nptcry/i");          // # points in detector
    m_ptree->Branch("ntpc",         &fntpc,         "ntpc/i");            // # TPC
    m_ptree->Branch("npttpc",       fnpttpc.data(), "npttpc[ntpc]/i");    // # points in each TPC
    m_ptree->Branch("napa",         &fnapa,         "napa/i");            // # APA
    m_ptree->Branch("nptapa",       fnptapa.data(), "nptapa[napa]/i");    // # points in each APA
    m_ptree->Branch("nrop",         &fnrop,         "nrop/i");            // # ROP
    m_ptree->Branch("nptrop",       fnptrop.data(), "nptrop[nrop]/i");    // # points in each ROP
    m_ptree->Branch("ptx",          fptx,           "ptx[npt]/F");
    m_ptree->Branch("pty",          fpty,           "pty[npt]/F");
    m_ptree->Branch("ptz",          fptz,           "ptz[npt]/F");
    m_ptree->Branch("ptt",          fptt,           "ptt[npt]/F");
    m_ptree->Branch("pte",          fpte,           "pte[npt]/F");
    m_ptree->Branch("pttpc",        fpttpc,         "pttpc[npt]/I");
    m_ptree->Branch("ptapa",        fptapa,         "ptapa[npt]/I");
    m_ptree->Branch("ptuchan",      fptuchan,       "ptuchan[npt]/I");
    m_ptree->Branch("ptvchan",      fptvchan,       "ptvchan[npt]/I");
    m_ptree->Branch("ptzchan",      fptzchan,       "ptzchan[npt]/I");
    m_ptree->Branch("ptutick",      fptutick,       "ptutick[npt]/F");
    m_ptree->Branch("ptvtick",      fptvtick,       "ptvtick[npt]/F");
    m_ptree->Branch("ptztick",      fptztick,       "ptztick[npt]/F");
    // length of track in detector
    m_ptree->Branch("detlen",       &fdetlen,       "detlen/F");
    m_ptree->Branch("dettickmin",   &fdettickmin,   "dettickmin/F");
    m_ptree->Branch("dettickmax",   &fdettickmax,   "dettickmax/F");
    m_ptree->Branch("detx1",        &fdetx1,        "detx1/F");
    m_ptree->Branch("dety1",        &fdety1,        "dety1/F");
    m_ptree->Branch("detz1",        &fdetz1,        "detz1/F");
    m_ptree->Branch("detx2",        &fdetx2,        "detx2/F");
    m_ptree->Branch("dety2",        &fdety2,        "dety2/F");
    m_ptree->Branch("detz2",        &fdetz2,        "detz2/F");
  }
}
 
//************************************************************************

// Destructor
MCTrajectoryFollower::~MCTrajectoryFollower() { }
   
//************************************************************************

int MCTrajectoryFollower::beginEvent(const art::Event& event, const MCParticleVector& pars) {
  const string myname = "MCTrajectory::beginEvent: ";

  // Check gemetry helper.
  if ( m_geohelp == nullptr ) {
    cout << myname << "ERROR: Geometry helper is absent." << endl;
    return 1;
  }

  // Start by fetching some event information.
  fevent  = event.id().event(); 
  fRun    = event.run();
  fSubRun = event.subRun();
  if ( m_dbg > 0 ) cout << myname << "Processing run " << fRun << "-" << fSubRun
       << ", event " << fevent << endl;

  // Create string representation of the event number.
  ostringstream ssevt;
  ssevt << fevent;
  string sevt = ssevt.str();
  string sevtf = sevt;
  while ( sevtf.size() < 4 ) sevtf = "0" + sevtf;

  m_ndetptmap.clear();
  // Loop over particles and fetch the # points inside the detector for each.
  // Might later want to add the # descendants with points inside the detector.
  for ( auto const& par : pars ) {
    unsigned int npt = 0;
    unsigned int tid = par.TrackId();
    size_t numberTrajectoryPoints = par.NumberTrajectoryPoints();
      for ( unsigned int ipt=0; ipt<numberTrajectoryPoints; ++ipt ) {
      const auto& pos = par.Position(ipt);
      double xyzt[4] = {pos.X(), pos.Y(), pos.Z(), pos.T()};
      geo::TPCID tpcid = m_geohelp->geometry()->FindTPCAtPosition(xyzt);
      if ( tpcid.isValid ) ++npt;
    }
    m_ndetptmap[tid] = npt;
  }

  return 0;
}

//************************************************************************

int MCTrajectoryFollower::endEvent() {
  const string myname = "MCTrajectory::endEvent: ";
  if ( m_filltree ) {
    m_ptree->Fill();
  }
  return 0;
}

//************************************************************************

int MCTrajectoryFollower::addMCParticle(const MCParticle& particle, MCTrackPerf* pmctp) {
  const string myname = "MCTrajectory::addMCParticle: ";
  if ( m_geohelp == nullptr ) {
    cout << myname << "ERROR: Geometry helper is absent." << endl;
    return 1;
  }
  const GeoHelper& geohelp = *m_geohelp;
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
    if ( m_ndetptmap[tid] ) fdetchild[fndetchild++] = tid;
  }
  if ( m_dbg > 2 ) {
    cout << myname << "ID=" << ftrackid << ", PDG=" << fpdg
         << ", RPDG=" << frpdg
         << ", status=" << particle.StatusCode()
         << ", Process: " << particle.Process() << "(" << fproc << ")"
         << ", Parent: " << fparent
         << endl;
  }

  size_t numberTrajectoryPoints = particle.NumberTrajectoryPoints();
  int last = numberTrajectoryPoints - 1;
  const TLorentzVector& positionStart = particle.Position(0);
  const TLorentzVector& positionEnd   = particle.Position(last);
  const TLorentzVector& momentumStart = particle.Momentum(0);
  const TLorentzVector& momentumEnd   = particle.Momentum(last);

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
    geo::TPCID tpcid = geohelp.geometry()->FindTPCAtPosition(xyzt);
    unsigned int icry = geohelp.geometry()->FindCryostatAtPosition(xyzt);
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
    if ( m_dbg > 3 ) cout << myname << "  MC Particle " << ftrackid
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
        if ( m_dbg > 3 ) {
          cout << myname << "    Entering TPC " << tpcid.TPC;
          if ( tpcid0.isValid ) cout << ", exiting TPC " << tpcid0.TPC;
          cout << endl;
        }
      }
      if ( tpcid0.isValid && (!tpcid.isValid or tpcid.TPC != tpcid0.TPC ) ) {
        ++fntpcout;
        if ( m_dbg > 3 ) {
          cout << myname << "    Exiting TPC " << tpcid0.TPC;
          if ( tpcid.isValid ) cout << ", entering TPC " << tpcid.TPC;
          cout << endl;
        }
      }
      if ( icry!=notcry && icry!=icry0 ) {
        ++fncryin;
        if ( m_dbg > 3 ) {
          cout << myname << "    Entering cryostat " << icry;
          if ( icry0 != notcry ) cout << ", exiting cryostat " << icry0;
          cout << endl;
        }
      }
      if ( icry0!=notcry && icry!=icry0 ) {
        ++fncryout;
        if ( m_dbg > 3 ) {
          cout << myname << "    Exiting cryostat " << icry0;
          if ( icry != notcry ) cout << ", entering cryostat " << icry;
          cout << endl;
        }
      }
      if ( lastpoint && m_dbg > 3 ) {
        cout << myname << "    Last trajectory point; # children is " << fnchild << endl; 
      }
    }
    // Fill MCTrackPerf with dE for each pair of adjacent tracjectory points.
    // Increase the number of points to ensure granularity less than fmcpdsmax;
    if ( ipt != 0 && indet0 && indet) {
      double dx = x - x0;
      double dy = y - y0;
      double dz = z - z0;
      double dt = t - t0;
      double ds = sqrt(dx*dx + dy*dy + dz*dz);
      unsigned int nstep = ds/m_dsmax + 1;
      double invstep = 1.0/nstep;
      double destep = invstep*detot;
      if ( m_dbg > 3 ) cout << myname << "    # steps: " << nstep << endl;
      for ( unsigned int istp=0; istp<nstep; ++istp ) {
        double xa = x0 + (istp+0.5)*invstep*dx;
        double ya = y0 + (istp+0.5)*invstep*dy;
        double za = z0 + (istp+0.5)*invstep*dz;
        double ta = t0 + (istp+0.5)*invstep*dt;
        double postim[4] = {xa, ya, za, ta};
        PlanePositionVector pps = geohelp.planePositions(postim);
        for ( const auto& pp : pps ) {
          if ( ! pp.valid ) cout << myname << "    Invalid plane position!" << endl;
          if ( m_dbg > 3 ) {
            cout << myname << "    Filling: "
                 << "tick=" << pp.tick
                 << ", chan=" << pp.ropchannel
                 << ", DE=" << destep << " MeV" << endl;
          }
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
  }  // End loop over trajectory points.

  return 0;

}

//**********************************************************************
