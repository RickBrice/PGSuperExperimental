#include "stdafx.h"
#include "PyLimits.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

void CPyLimits::Init(IBroker* pBroker)
{
   pBroker->GetInterface(IID_ILimits, (IUnknown**)&m_pLimits);
}

void CPyLimits::Reset()
{
   m_pLimits.Release();
}

Float64 CPyLimits::GetMaxSlabFc(pgsTypes::ConcreteType concType) const
{
   return m_pLimits->GetMaxSlabFc(concType);
}

Float64 CPyLimits::GetMaxSegmentFci(pgsTypes::ConcreteType concType) const
{
   return m_pLimits->GetMaxSegmentFci(concType);
}

Float64 CPyLimits::GetMaxSegmentFc(pgsTypes::ConcreteType concType) const
{
   return m_pLimits->GetMaxSegmentFc(concType);
}

Float64 CPyLimits::GetMaxClosureFci(pgsTypes::ConcreteType concType) const
{
   return m_pLimits->GetMaxClosureFci(concType);
}

Float64 CPyLimits::GetMaxClosureFc(pgsTypes::ConcreteType concType) const
{
   return m_pLimits->GetMaxClosureFc(concType);
}

Float64 CPyLimits::GetMaxConcreteUnitWeight(pgsTypes::ConcreteType concType) const
{
   return m_pLimits->GetMaxConcreteUnitWeight(concType);
}

Float64 CPyLimits::GetMaxConcreteAggSize(pgsTypes::ConcreteType concType) const
{
   return m_pLimits->GetMaxConcreteAggSize(concType);
}
