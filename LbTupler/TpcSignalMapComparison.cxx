// TpcSignalMapComparison.cxx

#include "TpcSignalMapComparison.h"
#include "TpcSignalMap.h"

typedef TpcSignalMap::Channel Channel;

//**********************************************************************

TpcSignalMapComparison::
TpcSignalMapComparison(const TpcSignalMap& tsm1, const TpcSignalMap& tsm2)
: m_tsm1(tsm1), m_tsm2(tsm2), m_chanfrac(-1.0) {
  Channel ch2min = tsm2.channelMin();
  Channel ch2max = tsm2.channelMax();
  unsigned int nchanNum = 0;
  unsigned int nchanDen = 0;
  for ( Channel ch1=tsm1.channelMin(); ch1<=tsm1.channelMax(); ++ch1 ) {
    if ( ch1 >= ch2min && ch1 <= ch2max ) ++nchanNum;
    ++nchanDen;
  }
  if ( nchanDen > 0 ) m_chanfrac = double(nchanNum)/double(nchanDen);
}

//**********************************************************************
