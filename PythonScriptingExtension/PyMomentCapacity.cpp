#include "stdafx.h"
#include "PyMomentCapacity.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

void CPyMomentCapacity::Init(IBroker* pBroker)
{
   pBroker->GetInterface(IID_IMomentCapacity, (IUnknown**)&m_pMomentCapacity);
}

void CPyMomentCapacity::Reset()
{
   m_pMomentCapacity.Release();
}

const MINMOMENTCAPDETAILS& CPyMomentCapacity::GetMinMomentCapacityDetails(IntervalIndexType intervalIdx, const pgsPointOfInterest& poi, bool bPositiveMoment)
{
   return *(m_pMomentCapacity->GetMinMomentCapacityDetails(intervalIdx, poi, bPositiveMoment));
}
