#include "stdafx.h"
#include "PyGirderModelView.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

void CPyGirderModelView::Init(IGirderModelViewController* pViewController)
{
   m_pMyViewController = pViewController;
   CPyViewControllerBase::Init(pViewController);
}

bool CPyGirderModelView::GetSyncWithBridgeModelView() const
{
   return m_pMyViewController->SyncWithBridgeModelView();
}

void CPyGirderModelView::SetSyncWithBridgeModelView(bool bSync)
{
   m_pMyViewController->SyncWithBridgeModelView();
}

const CGirderKey& CPyGirderModelView::GetGirder() const
{
   return m_pMyViewController->GetGirder();
}

void CPyGirderModelView::SelectGirder(const CGirderKey& girderKey)
{
   m_pMyViewController->SelectGirder(girderKey);
}

EventIndexType CPyGirderModelView::GetEvent() const
{
   return m_pMyViewController->GetEvent();
}

bool CPyGirderModelView::SetEvent(EventIndexType eventIdx)
{
   return m_pMyViewController->SetEvent(eventIdx);
}

Float64 CPyGirderModelView::GetCutLocation() const
{
   return m_pMyViewController->GetCutLocation();
}

void CPyGirderModelView::CutAt(Float64 Xg)
{
   m_pMyViewController->CutAt(Xg);
}

void CPyGirderModelView::CutAtNext()
{
   m_pMyViewController->CutAtNext();
}

void CPyGirderModelView::CutAtPrev()
{
   m_pMyViewController->CutAtPrev();
}

Float64 CPyGirderModelView::GetMinCutLocation() const
{
   Float64 min, max;
   m_pMyViewController->GetCutRange(&min, &max);
   return min;
}

Float64 CPyGirderModelView::GetMaxCutLocation() const
{
   Float64 min, max;
   m_pMyViewController->GetCutRange(&min, &max);
   return max;
}

void CPyGirderModelView::SetShowStrands(bool bShow)
{
   m_pMyViewController->ShowStrands(bShow);
}

bool CPyGirderModelView::GetShowStrands() const
{
   return m_pMyViewController->ShowStrands();
}

void CPyGirderModelView::SetShowStrandCG(bool bShow)
{
   m_pMyViewController->ShowStrandCG(bShow);
}

bool CPyGirderModelView::GetShowStrandCG() const
{
   return m_pMyViewController->ShowStrandCG();
}

void CPyGirderModelView::SetShowCG(bool bShow)
{
   m_pMyViewController->ShowCG(bShow);
}

bool CPyGirderModelView::GetShowCG() const
{
   return m_pMyViewController->ShowCG();
}

void CPyGirderModelView::SetShowSectionProperties(bool bShow)
{
   m_pMyViewController->ShowSectionProperties(bShow);
}

bool CPyGirderModelView::GetShowSectionProperties() const
{
   return m_pMyViewController->ShowSectionProperties();
}

void CPyGirderModelView::SetShowDimensions(bool bShow)
{
   m_pMyViewController->ShowDimensions(bShow);
}

bool CPyGirderModelView::GetShowDimensions() const
{
   return m_pMyViewController->ShowDimensions();
}

void CPyGirderModelView::SetShowLongitudinalReinforcement(bool bShow)
{
   m_pMyViewController->ShowLongitudinalReinforcement(bShow);
}

bool CPyGirderModelView::GetShowLongitudinalReinforcement() const
{
   return m_pMyViewController->ShowLongitudinalReinforcement();
}

void CPyGirderModelView::SetShowTransverseReinforcement(bool bShow)
{
   m_pMyViewController->ShowTransverseReinforcement(bShow);
}

bool CPyGirderModelView::GetShowTransverseReinforcement() const
{
   return m_pMyViewController->ShowTransverseReinforcement();
}

void CPyGirderModelView::SetShowLoads(bool bShow)
{
   m_pMyViewController->ShowLoads(bShow);
}

bool CPyGirderModelView::GetShowLoads() const
{
   return m_pMyViewController->ShowLoads();
}

void CPyGirderModelView::SetSchematic(bool bSchematic)
{
   m_pMyViewController->Schematic(bSchematic);
}

bool CPyGirderModelView::GetSchematic() const
{
   return m_pMyViewController->Schematic();
}
