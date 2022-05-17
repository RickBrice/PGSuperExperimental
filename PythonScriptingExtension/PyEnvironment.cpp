#include "stdafx.h"
#include "PyEnvironment.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

void CPyEnvironment::Init(IBroker* pBroker)
{
   pBroker->GetInterface(IID_IEnvironment, (IUnknown**)&m_pEnvironment);
}

void CPyEnvironment::Reset()
{
   m_pEnvironment.Release();
}

enumExposureCondition CPyEnvironment::GetExposureCondition() const
{
   return m_pEnvironment->GetExposureCondition();
}

void CPyEnvironment::SetExposureCondition(enumExposureCondition newVal)
{
   m_pEnvironment->SetExposureCondition(newVal);
}

Float64 CPyEnvironment::GetRelHumidity() const
{
   return m_pEnvironment->GetRelHumidity();
}

void CPyEnvironment::SetRelHumidity(Float64 newVal)
{
   m_pEnvironment->SetRelHumidity(newVal);
}
