// MCTrajectoryFollower.h

// David Adams
// May 2015
//
// Uses MCParticles to build a Root tree and to add hits to TpcSignalMap objects.
// Interpolates between the points on the MCParticle trajectory.

#ifndef MCTracjectoryFollower_Module
#define MCTracjectoryFollower_Module

#include <vector>
#include <map>
#include "TpcTypes.h"

namespace art {
class Event;
}
namespace simb {
class MCParticle;
}
class TTree;
class GeoHelper;
class TpcSignalMap;

class MCTrajectoryFollower {

public:

  typedef std::vector<simb::MCParticle> MCParticleVector; 

  // Ctor.
  //   dsmax  - Maximum step size used inf following trajectory.
  //            Interpolation is used when the steps are larger.
  //   tname - Name for the Root tree. If blank, no tree is filled.
  //   geohelp - GeoHelper used to access geometry information
  MCTrajectoryFollower(double dsmax, std::string tname, const GeoHelper* geohelp =nullptr,
                       unsigned int minNptdet =0, int dbg =0);

  // Dtor.
  ~MCTrajectoryFollower();

  // This method is called at the start of each event.
  int beginEvent(const art::Event& event, const MCParticleVector& pars =MCParticleVector());

  // Call this to add an MCParticle to the current event.
  //   par - The input MCParticle.
  //   pmctp - If non-null, the MCparticle is used to fill this performance object.
  //   useDescendants - If true, descendants are also added to the signal map.
  //   ptids - if non-null, the ID for this track and descendants (if used) are added to
  //           this vector
  // Returns 0 if particle is accepted, >0 if rejected, <0 for error.
  int addMCParticle(const simb::MCParticle& par, TpcSignalMap* pmtsm =nullptr,
                    bool useDescendants =false, tpc::IndexVector* ptids =nullptr);

private:

  // Control parameters.
  int m_dbg;                     // Debug level. Larger for more log noise.
  std::string m_tname;           // Tree name.
  bool m_filltree;               // Create and fill the MCParticle tree.
  double m_dsmax;                // Maximum step size for interpolation of trajectory
  unsigned int m_minNptdet;      // Min # points in detector to keep the particle
  int m_generation;              // Depth in the descendant tree.

  // The tree.
  TTree* m_ptree;

  // The variables that will go into the n-tuple.
  int fevent;
  int fRun;
  int fSubRun;
  int fpdg;                  // PDG ID
  int frpdg;                 // reduced PDG ID
  int fproc;                 // Process
  int fitrk;                 // Track index (0, 1, ...)
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
  static const unsigned int fmaxchild = 50;
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

  // # detector points for each MC particle.
  std::map<unsigned int, unsigned int> m_ndetptmap;

  const art::Event* m_pevt;
  const MCParticleVector* m_ppars;

}; // class MCTrajectory

#endif
