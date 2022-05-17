#include "stdafx.h"
#include "PyBridgeModelView.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

void CPyBridgeModelView::Init(IBridgeModelViewController* pViewController,ISelection* pSelection)
{
   m_pViewController = pViewController;
   m_pSelection = pSelection;
   CPyViewControllerBase::Init(pViewController);
}

void CPyBridgeModelView::GetGroupRange(GroupIndexType* pStartGroupIdx, GroupIndexType* pEndGroupIdx) const
{
   m_pViewController->GetGroupRange(pStartGroupIdx, pEndGroupIdx);
}

void CPyBridgeModelView::SetGroupRange(GroupIndexType startGroupIdx, GroupIndexType endGroupIdx)
{
   m_pViewController->SetGroupRange(startGroupIdx, endGroupIdx);
}

Float64 CPyBridgeModelView::GetCutStation() const
{
   return m_pViewController->GetCutStation();
}

void CPyBridgeModelView::SetCutStation(Float64 station)
{
   m_pViewController->SetCutStation(station);
}

void CPyBridgeModelView::SetViewMode(IBridgeModelViewController::ViewMode mode)
{
   m_pViewController->SetViewMode(mode);
}

IBridgeModelViewController::ViewMode CPyBridgeModelView::GetViewMode() const
{
   return m_pViewController->GetViewMode();
}

void CPyBridgeModelView::SetNorthUp(bool bNorthUp)
{
   m_pViewController->NorthUp(bNorthUp);
}

bool CPyBridgeModelView::GetNorthUp() const
{
   return m_pViewController->NorthUp();
}

void CPyBridgeModelView::SetShowLabels(bool bShowLabels)
{
   m_pViewController->ShowLabels(bShowLabels);
}

bool CPyBridgeModelView::GetShowLabels() const
{
   return m_pViewController->ShowLabels();
}

void CPyBridgeModelView::SetShowDimensions(bool bShowDimensions)
{
   m_pViewController->ShowDimensions(bShowDimensions);
}

bool CPyBridgeModelView::GetShowDimensions() const
{
   return m_pViewController->ShowDimensions();
}

void CPyBridgeModelView::SetShowBridge(bool bShowBridge)
{
   m_pViewController->ShowBridge(bShowBridge);
}

bool CPyBridgeModelView::GetShowBridge() const
{
   return m_pViewController->ShowBridge();
}

void CPyBridgeModelView::SetSchematic(bool bSchematic)
{
   m_pViewController->Schematic(bSchematic);
}

bool CPyBridgeModelView::GetSchematic() const
{
   return m_pViewController->Schematic();
}

void CPyBridgeModelView::ClearSelection()
{
   m_pSelection->ClearSelection();
}

PierIndexType CPyBridgeModelView::GetSelectedPier()
{
   return m_pSelection->GetSelectedPier();
}

SpanIndexType CPyBridgeModelView::GetSelectedSpan()
{
   return m_pSelection->GetSelectedSpan();
}

CGirderKey CPyBridgeModelView::GetSelectedGirder()
{
   return m_pSelection->GetSelectedGirder();
}

CSegmentKey CPyBridgeModelView::GetSelectedSegment()
{
   return m_pSelection->GetSelectedSegment();
}

bool CPyBridgeModelView::IsDeckSelected()
{
   return m_pSelection->IsDeckSelected();
}

bool CPyBridgeModelView::IsAlignmentSelected()
{
   return m_pSelection->IsAlignmentSelected();
}

bool CPyBridgeModelView::IsRailingSystemSelected(pgsTypes::TrafficBarrierOrientation orientation)
{
   return m_pSelection->IsRailingSystemSelected(orientation);
}

void CPyBridgeModelView::SelectPier(PierIndexType pierIdx)
{
   m_pSelection->SelectPier(pierIdx);
}

void CPyBridgeModelView::SelectSpan(SpanIndexType spanIdx)
{
   m_pSelection->SelectSpan(spanIdx);
}

void CPyBridgeModelView::SelectGirder(const CGirderKey& girderKey)
{
   m_pSelection->SelectGirder(girderKey);
}

void CPyBridgeModelView::SelectSegment(const CSegmentKey& segmentKey)
{
   m_pSelection->SelectSegment(segmentKey);
}

void CPyBridgeModelView::SelectDeck()
{
   m_pSelection->SelectDeck();
}

void CPyBridgeModelView::SelectAlignment()
{
   m_pSelection->SelectAlignment();
}

void CPyBridgeModelView::SelectRailingSystem(pgsTypes::TrafficBarrierOrientation orientation)
{
   m_pSelection->SelectRailingSystem(orientation);
}
