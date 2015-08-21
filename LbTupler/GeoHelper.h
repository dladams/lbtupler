// GeoHelper.h
//
// David Adams
// April 2015
//
// Helper for accessing LArSoft geometry information.

#ifndef GeoHelper_H
#define GeoHelper_H

#include <string>
#include <iostream>
#include <vector>
#include <map>
// We need the folowing because simple enums cannot be forward declared.
#include "SimpleTypesAndConstants/geo_types.h"
#include "PlanePosition.h"

namespace geo {
class GeometryCore;
}
namespace util {
class LArProperties;
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
  // Note: DetectorProperties has non-const methods.
  GeoHelper(const geo::GeometryCore* pgeo, bool useChannels, Status dbg =0);

  // Ctor from geometry name, e.g. "dune10kt_v1".
  // This is for use outside the art framework.
  // if useChannels is true, then the apprpriate channel map is loaded.
  GeoHelper(std::string gname, bool useChannels =false, Status dbg =0);

  // Return the geometry.
  const geo::GeometryCore* geometry() const { return m_pgeo; }

  // Do we have channel mapping.
  bool haveChannelMap() const { return m_haveChannelMap; }

  // Return the LAr properties.
  util::LArProperties* larProperties() const;

  // Return the detector properties.
  util::DetectorProperties* detectorProperties() const;

  // Return dimensions for a TPC.
  // useActive = true for the "active" dimensions.
  double  width(unsigned int icry, unsigned int itpc, bool useActive =false) const;
  double height(unsigned int icry, unsigned int itpc, bool useActive =false) const;
  double length(unsigned int icry, unsigned int itpc, bool useActive =false) const;

  // Return the center of a TPC volume.
  Status tpcCenter(unsigned int icry, unsigned int itpc, double* pos) const;

  // Return the corners of a TPC volume.
  Status tpcCorners(unsigned int icry, unsigned int itpc, double* pos1, double* pos2) const;

  // Total number of cryostats.
  Index ncryostat() const;

  // Total number of APAs.
  Index napa() const { return m_napa; }

  // Total number of TPCs.
  Index ntpc();
  Index ntpc() const;

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
  const IndexVector& ropTpcs(Index irop) const { return m_roptpc.at(irop); }

  // Return the ROP for a TPC plane.
  Index rop(geo::PlaneID pid) const;

  // Print detector description.
  std::ostream& print(std::ostream& out =std::cout, int iopt =0, std::string prefix ="") const;

  // Geometry dump from Michelle.
  std::ostream& dump(std::ostream& out =std::cout) const;

  // Return the Rop for a channel.
  // Returns nrop() if the channel is not valid.
  Index channelRop(Index chan) const;

  // Return all TPC plane postitions a spacetime point.
  // xyzt = {x, y, z, t} [cm,ns]
  PlanePositionVector planePositions(const double xyzt[]) const;

private:

  // Fill info pertaining to APAs and ROPs assuming standard TPC-APA mapping.
  Status fillStandardApaMapping();

private:

  const geo::GeometryCore* m_pgeo;
  bool m_haveChannelMap;        // Does the geometry have a channel map.
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
  IndexVectorVector m_roptpc;   // Global TPC indices for each ROP.
  IndexVector m_ropapa;         // Global APA index for each ROP.
  ViewVector m_ropview;         // View (kU, kV, kZ) for each ROP
  NameVector m_ropname;         // Name for each ROP
  PlaneIDIndexMap m_tpprop;     // Global ROP index for each TPC plane.
};

#endif
