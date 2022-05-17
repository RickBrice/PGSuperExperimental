#pragma once

#include <IFace\MomentCapacity.h>

class CPyMomentCapacity
{
public:
   void Init(IBroker* pBroker);
   void Reset();

   const MINMOMENTCAPDETAILS& GetMinMomentCapacityDetails(IntervalIndexType intervalIdx, const pgsPointOfInterest& poi, bool bPositiveMoment);

private:
   CComPtr<IMomentCapacity> m_pMomentCapacity;
}; 
