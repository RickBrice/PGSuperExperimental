#pragma once

#include <IFace\Project.h>

class CPyEnvironment
{
public:
   void Init(IBroker* pBroker);
   void Reset();

   enumExposureCondition GetExposureCondition() const;
   void SetExposureCondition(enumExposureCondition newVal);
   Float64 GetRelHumidity() const;
   void SetRelHumidity(Float64 newVal);

private:
   CComPtr<IEnvironment> m_pEnvironment;
}; 
