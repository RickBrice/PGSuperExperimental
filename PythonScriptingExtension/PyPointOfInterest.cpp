#include "stdafx.h"
#include "PyPointOfInterest.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

void CPyPointOfInterest::Init(IBroker* pBroker)
{
   pBroker->GetInterface(IID_IPointOfInterest, (IUnknown**)&m_pPoi);
}

void CPyPointOfInterest::Reset()
{
   m_pPoi.Release();
}

std::vector<pgsPointOfInterest> CPyPointOfInterest::GetPointsOfInterest(const CSegmentKey& segmentKey) const
{
   PoiList poiList;
   m_pPoi->GetPointsOfInterest(segmentKey,&poiList);
   std::vector<pgsPointOfInterest> vPoi;
   MakePoiVector(poiList, &vPoi);
   return vPoi;
}

std::vector<pgsPointOfInterest> CPyPointOfInterest::GetSpanPointsOfInterest(const CSpanKey& spanKey) const
{
   PoiList poiList;
   m_pPoi->GetPointsOfInterest(spanKey,&poiList);
   std::vector<pgsPointOfInterest> vPoi;
   MakePoiVector(poiList, &vPoi);
   return vPoi;
}
