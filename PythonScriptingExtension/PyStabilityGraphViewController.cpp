#include "stdafx.h"
#include "PyStabilityGraphViewController.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

void CPyStabilityGraphViewController::Init(IStabilityGraphViewController* pViewController)
{
   m_pMyViewController = pViewController;
   CPyViewControllerBase::Init(pViewController);
}

const CSegmentKey& CPyStabilityGraphViewController::GetSegment() const
{
   return m_pMyViewController->GetSegment();
}

void CPyStabilityGraphViewController::SelectSegment(const CSegmentKey& segmentKey)
{
   m_pMyViewController->SelectSegment(segmentKey);
}

void CPyStabilityGraphViewController::SetViewMode(IStabilityGraphViewController::ViewMode mode)
{
   m_pMyViewController->SetViewMode(mode);
}

IStabilityGraphViewController::ViewMode CPyStabilityGraphViewController::GetViewMode() const
{
   return m_pMyViewController->GetViewMode();
}

void CPyStabilityGraphViewController::SetShowGrid(bool bShow)
{
   m_pMyViewController->ShowGrid(bShow);
}

bool CPyStabilityGraphViewController::GetShowGrid() const
{
   return m_pMyViewController->ShowGrid();
}

CPyStabilityGraphViewController* CPyStabilityGraphViewController::GetGraphController(CPyViewControllerBase* pController)
{
   auto pMyController = static_cast<CPyStabilityGraphViewController*>(pController);
   return pMyController;
}
