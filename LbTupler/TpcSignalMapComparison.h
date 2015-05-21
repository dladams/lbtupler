// TpcSignalMapComparison.h

#ifndef TpcSignalMapComparator_H
#define TpcSignalMapComparator_H

class TpcSignalMap;

// David Adams
// May 2015
//
// Compares two TpcSignalMap objects.

class TpcSignalMapComparison {

public:

  TpcSignalMapComparison(const TpcSignalMap& tsm1, const TpcSignalMap& tsm2);

  // Metrics.
  double channelFraction() const { return m_chanfrac; }

private:

  const TpcSignalMap& m_tsm1;
  const TpcSignalMap& m_tsm2;
  double m_chanfrac;

};

#endif
