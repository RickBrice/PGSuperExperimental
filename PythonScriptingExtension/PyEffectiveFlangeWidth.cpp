#include "stdafx.h"
#include "PyEffectiveFlangeWidth.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

void CPyEffectiveFlangeWidth::Init(IBroker* pBroker)
{
   pBroker->GetInterface(IID_IEffectiveFlangeWidth, (IUnknown**)&m_pEffectiveFlangeWidth);
}

void CPyEffectiveFlangeWidth::Reset()
{
   m_pEffectiveFlangeWidth.Release();
}

bool CPyEffectiveFlangeWidth::GetIgnoreEffectiveFlangeWidthLimits() const
{
   return m_pEffectiveFlangeWidth->IgnoreEffectiveFlangeWidthLimits();
}

void CPyEffectiveFlangeWidth::SetIgnoreEffectiveFlangeWidthLimits(bool bIgnore)
{
   m_pEffectiveFlangeWidth->IgnoreEffectiveFlangeWidthLimits(bIgnore);
}
