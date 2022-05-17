#pragma once

#include <IFace\PointOfInterest.h>

class CPyPointOfInterest
{
public:
   void Init(IBroker* pBroker);
   void Reset();

   std::vector<pgsPointOfInterest> GetPointsOfInterest(const CSegmentKey& segmentKey) const;
   std::vector<pgsPointOfInterest> GetSpanPointsOfInterest(const CSpanKey& spanKey) const;

private:
   CComPtr<IPointOfInterest> m_pPoi;
}; 
