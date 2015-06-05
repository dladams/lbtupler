// TpcSignalMapComparison.h

#ifndef TpcSignalMapComparison_H
#define TpcSignalMapComparison_H

// David Adams
// May 2015
//
// Compares two TpcSignalMap objects.

#include "TpcTypes.h"

class TpcSignalMap;

class TpcSignalMapComparison {

public:

  typedef tpc::Index Index;
  typedef tpc::Channel Channel;

  // Ctor.
  //  tsm1 - signal map for the reference
  //  tsm2 - signal map for the match
  TpcSignalMapComparison(const TpcSignalMap& tsm1, const TpcSignalMap& tsm2);

  // Return reference and match objects.
  const TpcSignalMap& reference() const { return m_tsm1; }
  const TpcSignalMap& match() const { return m_tsm2; }

  // Return the ROP found for this pair of objects.
  // Returns tpc::badIndex() if none.
  Index rop() const { return m_rop; }

  // Return the range of channels for this object's ROP.
  Channel beginChannel() const { return m_chbegin; }
  Channel endChannel() const { return m_chend; }

  // Return the number of channels in the ROP range.
  Index referenceChannelCount() const { return m_nchanref; }
  Index matchChannelCount() const { return m_nchanmat; }

  // Metrics.
  // channelFraction = fraction of reference channels included in the match
  double channelFraction() const;
  // binFraction = fraction of reference channel-tick bins included in match
  double binFraction() const;

private:

  const TpcSignalMap& m_tsm1;
  const TpcSignalMap& m_tsm2;
  Index m_rop;
  Channel m_chbegin;
  Channel m_chend;
  Index m_nchanref;
  Index m_nchanmat;

};

#endif
