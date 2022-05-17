#pragma once

#include <IFace\Project.h>

class CPyLoadModifiers
{
public:
   void Init(IBroker* pBroker);
   void Reset();

   void SetDuctilityModifier(ILoadModifiers::Level level, Float64 value);
   Float64 GetDuctilityFactor() const;
   ILoadModifiers::Level GetDuctilityLevel() const;

   void SetImportanceModifier(ILoadModifiers::Level level, Float64 value);
   Float64 GetImportanceFactor() const;
   ILoadModifiers::Level GetImportanceLevel() const;

   void SetRedundancyModifier(ILoadModifiers::Level level, Float64 value);
   Float64 GetRedundancyFactor() const;
   ILoadModifiers::Level GetRedundancyLevel() const;

private:
   CComPtr<ILoadModifiers> m_pLoadModifiers;
}; 
