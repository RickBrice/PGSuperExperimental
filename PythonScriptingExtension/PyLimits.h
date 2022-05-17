#pragma once

#include <IFace\Project.h>

class CPyLimits
{
public:
   void Init(IBroker* pBroker);
   void Reset();

   Float64 GetMaxSlabFc(pgsTypes::ConcreteType concType) const;
   Float64 GetMaxSegmentFci(pgsTypes::ConcreteType concType) const;
   Float64 GetMaxSegmentFc(pgsTypes::ConcreteType concType) const;
   Float64 GetMaxClosureFci(pgsTypes::ConcreteType concType) const;
   Float64 GetMaxClosureFc(pgsTypes::ConcreteType concType) const;
   Float64 GetMaxConcreteUnitWeight(pgsTypes::ConcreteType concType) const;
   Float64 GetMaxConcreteAggSize(pgsTypes::ConcreteType concType) const;

private:
   CComPtr<ILimits> m_pLimits;
}; 
