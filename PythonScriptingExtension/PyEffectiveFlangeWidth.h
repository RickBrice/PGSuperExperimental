#pragma once

#include <IFace\Project.h>

class CPyEffectiveFlangeWidth
{
public:
   void Init(IBroker* pBroker);
   void Reset();

   bool GetIgnoreEffectiveFlangeWidthLimits() const;
   void SetIgnoreEffectiveFlangeWidthLimits(bool bIgnore);

private:
   CComPtr<IEffectiveFlangeWidth> m_pEffectiveFlangeWidth;
}; 
