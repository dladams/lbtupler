// GeoHelper.h
//
// David Adams
// April 2015
//
// Helper for accessing LArSoft geometry information.

#include <string>
#include <iosfwd>
#include <vector>
#include <map>
// We need the folowing because simple enums cannot be forward declared.
#include "SimpleTypesAndConstants/geo_types.h"

namespace geo {
class Geometry;
}
namespace util {
class DetectorProperties;
}
class GeoHelper {

public:
  typedef int Status;
  typedef unsigned int Index;
  typedef std::vector<Index> IndexVector;
  typedef std::vector<IndexVector> IndexVectorVector;
  typedef std::string Name;
  typedef std::vector<Name> NameVector;
  typedef std::vector<geo::View_t> ViewVector;
  typedef std::map<geo::PlaneID, Index> PlaneIDIndexMap;

public:

  // Value for an invalid index.
  static Index badIndex();

public:

  // Ctor from LArSoft geometry service and detector properties.
  GeoHelper(const geo::Geometry* pgeo, const util::DetectorProperties* pdetp =nullptr, Status dbg =0);

  // Return the corners of the active TPC volume.
  Status tpcCorners(unsigned int icry, unsigned int itpc, double* pos1, double* pos2) const;

  // Total number of cryostats.
  Index ncryostat() const;

  // Total number of APAs.
  Index napa() const { return m_napa; }

  // Total number of TPCs.
  Index ntpc() const { return m_ntpc; }

  // Total number of TPC planes.
  Index ntpp() const { return m_ntpp; }

  // Total number of readout planes.
  Index nrop() const { return m_nrop; }

  // TPC index (0, ..., ntpc()-1)
  Index tpcIndex(unsigned int icry, unsigned int itpcplane);

  // TPC properties.
  Index tpcCryostat(Index itpc) const { return m_tpccry.at(itpc); }
  Index tpcApa(Index itpc) const { return m_tpcapa.at(itpc); }
  Name tpcName(Index itpc) const { return m_tpcname.at(itpc); }

  // ROP properties.
  Name ropName(Index irop) const { return m_ropname.at(irop); }
  Index ropFirstChannel(Index irop) const { return m_ropfirstchan.at(irop); }
  Index ropNChannel(Index irop) const { return m_ropnchan.at(irop); }
  geo::View_t ropView(Index irop) const { return m_ropview.at(irop); }

  // Return the ROP for a TPC plane.
  Index rop(geo::PlaneID pid) const;

  // Return a global plane index for a plane in a TPC
  // Print detector description.
  std::ostream& print(std::ostream& out, int iopt =0, std::string prefix ="") const;

private:

  // Fill info pertaining to APAs and ROPs assuming standard TPC-APA mapping.
  Status fillStandardApaMapping();

private:

  const geo::Geometry* m_pgeo;
  const util::DetectorProperties* m_pdetp;
  Status m_dbg;
  Index m_ntpc;                 // Total # TPCs in the detector
  Index m_ntpp;                 // Total # TPC planes in the detector
  Index m_napa;                 // Total # APAs in the detector
  Index m_nrop;                 // Total # readout planes in the detector
  IndexVector m_tpccry;         // cryostat for each TPC
  IndexVector m_tpcapa;         // APA for each TPC
  NameVector m_tpcname;         // Name for each TPC
  IndexVector m_apanrop;        // # ROP for each APA
  IndexVector m_ropfirstchan;   // First channel for each ROP;
  IndexVector m_ropnchan;       // # channels for each ROP;
  IndexVector m_roptpc;         // Global TPC index for each ROP.
  IndexVector m_ropapa;         // Global APA index for each ROP.
  ViewVector m_ropview;         // View (kU, kV, kZ) for each ROP
  NameVector m_ropname;         // Name for each ROP
  PlaneIDIndexMap m_tpprop;     // Global ROP index for each TPC plane.
};
