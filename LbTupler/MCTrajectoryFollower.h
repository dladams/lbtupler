// MCTrajectoryFollower.h

// David Adams
// May 2015
//
// Uses MCParticles to build a Root tree and to add hits to MCTrackPerf objects.
// Interpolates between the points on the MCParticle trajectory.

#ifndef MCTracjectoryFollower_Module
#define MCTracjectoryFollower_Module

#include <vector>
#include <map>

namespace art {
class Event;
}
namespace simb {
class MCParticle;
}
class TTree;
class GeoHelper;
class MCTrackPerf;

class MCTrajectoryFollower {

public:

  typedef std::vector<simb::MCParticle> MCParticleVector; 

  // Ctor.
  //   filltree - If true the Root tree is created
  //   geohelp - GeoHelper used to access geometry information
  MCTrajectoryFollower(double dsmax, bool filltree =true, const GeoHelper* geohelp =nullptr, int dbg =0);

  // Dtor.
  ~MCTrajectoryFollower();

  // This method is called at the start of each event.
  int beginEvent(const art::Event& event, const MCParticleVector& pars =MCParticleVector());

  // This method is called at the end of each event.
  int endEvent();

  // Call this to add an MCParticle to the current event.
  //   par - The input MCParticle.
  //   pmctp - If non-null, the MCparticle is used to fill this performance object.
  int addMCParticle(const simb::MCParticle& par, MCTrackPerf* pmctp =nullptr);

private:

  // Control parameters.
  int m_dbg;                        // Debug level. Larger for more log noise.
  bool m_filltree;                  // Create and fill the MCParticle tree.
  double m_dsmax;                   // Maximum step size for interpolation of trajectory

  // The tree.
  TTree* m_ptree;

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
  std::vector<int> fnpttpc;   // # trajectory point in each TPC
  std::vector<int> fnptapa;   // # trajectory point in each APA
  std::vector<int> fnptrop;   // # trajectory point in each ROP
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

  // Geometry.
  const GeoHelper* m_geohelp;

  // # detector poings for each MC particle.
  std::map<unsigned int, unsigned int> m_ndetptmap;

}; // class MCTrajectory

#endif
