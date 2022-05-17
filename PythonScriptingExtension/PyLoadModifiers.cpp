#include "stdafx.h"
#include "PyLoadModifiers.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

void CPyLoadModifiers::Init(IBroker* pBroker)
{
   pBroker->GetInterface(IID_ILoadModifiers, (IUnknown**)&m_pLoadModifiers);
}

void CPyLoadModifiers::Reset()
{
   m_pLoadModifiers.Release();
}

void CPyLoadModifiers::SetDuctilityModifier(ILoadModifiers::Level level, Float64 value)
{
   m_pLoadModifiers->SetDuctilityFactor(level, value);
}

ILoadModifiers::Level CPyLoadModifiers::GetDuctilityLevel() const
{
   return m_pLoadModifiers->GetDuctilityLevel();
}

Float64 CPyLoadModifiers::GetDuctilityFactor() const
{
   return m_pLoadModifiers->GetDuctilityFactor();
}

void CPyLoadModifiers::SetImportanceModifier(ILoadModifiers::Level level, Float64 value)
{
   m_pLoadModifiers->SetImportanceFactor(level, value);
}

ILoadModifiers::Level CPyLoadModifiers::GetImportanceLevel() const
{
   return m_pLoadModifiers->GetImportanceLevel();
}

Float64 CPyLoadModifiers::GetImportanceFactor() const
{
   return m_pLoadModifiers->GetImportanceFactor();
}

void CPyLoadModifiers::SetRedundancyModifier(ILoadModifiers::Level level, Float64 value)
{
   m_pLoadModifiers->SetRedundancyFactor(level, value);
}

ILoadModifiers::Level CPyLoadModifiers::GetRedundancyLevel() const
{
   return m_pLoadModifiers->GetRedundancyLevel();
}

Float64 CPyLoadModifiers::GetRedundancyFactor() const
{
   return m_pLoadModifiers->GetRedundancyFactor();
}
