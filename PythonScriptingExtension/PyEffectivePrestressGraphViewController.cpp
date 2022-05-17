#include "stdafx.h"
#include "PyEffectivePrestressGraphViewController.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

void CPyEffectivePrestressGraphViewController::Init(IEffectivePrestressGraphViewController* pViewController)
{
   m_pMyViewController = pViewController;
   CPyViewControllerBase::Init(pViewController);
}

IntervalIndexType CPyEffectivePrestressGraphViewController::GetMinimumInterval() const
{
   IntervalIndexType minIntervalIdx, maxIntervalIdx;
   m_pMyViewController->GetIntervalRange(&minIntervalIdx, &maxIntervalIdx);
   return minIntervalIdx;
}

IntervalIndexType CPyEffectivePrestressGraphViewController::GetMaximumInterval() const
{
   IntervalIndexType minIntervalIdx, maxIntervalIdx;
   m_pMyViewController->GetIntervalRange(&minIntervalIdx, &maxIntervalIdx);
   return maxIntervalIdx;
}

void CPyEffectivePrestressGraphViewController::SelectInterval(IntervalIndexType intervalIdx)
{
   m_pMyViewController->SelectInterval(intervalIdx);
}

void CPyEffectivePrestressGraphViewController::SelectIntervals(const boost::python::list& listIntervals)
{
   boost::python::stl_input_iterator<IntervalIndexType> begin(listIntervals), end;
   std::vector<IntervalIndexType> vIntervals(begin, end);
   m_pMyViewController->SelectIntervals(vIntervals);
}

boost::python::list CPyEffectivePrestressGraphViewController::GetSelectedIntervals() const
{
   std::vector<IntervalIndexType> vIntervals = m_pMyViewController->GetSelectedIntervals();
   boost::python::list list;
   std::for_each(vIntervals.cbegin(), vIntervals.cend(), [&list](const auto& intervalIdx) {list.append(intervalIdx);});
   return list;
}

const CGirderKey& CPyEffectivePrestressGraphViewController::GetGirder() const
{
   return m_pMyViewController->GetGirder();
}

void CPyEffectivePrestressGraphViewController::SelectGirder(const CGirderKey& girderKey)
{
   m_pMyViewController->SelectGirder(girderKey);
}

void CPyEffectivePrestressGraphViewController::SetViewMode(IEffectivePrestressGraphViewController::ViewMode mode)
{
   m_pMyViewController->SetViewMode(mode);
}

IEffectivePrestressGraphViewController::ViewMode CPyEffectivePrestressGraphViewController::GetViewMode() const
{
   return m_pMyViewController->GetViewMode();
}

void CPyEffectivePrestressGraphViewController::SetStrandType(IEffectivePrestressGraphViewController::StrandType strandType)
{
   m_pMyViewController->SetStrandType(strandType);
}

IEffectivePrestressGraphViewController::StrandType CPyEffectivePrestressGraphViewController::GetStrandType() const
{
   return m_pMyViewController->GetStrandType();
}

void CPyEffectivePrestressGraphViewController::SetDuct(IEffectivePrestressGraphViewController::DuctType ductType, DuctIndexType ductIdx)
{
   m_pMyViewController->SetDuct(ductType,ductIdx);
}

DuctIndexType CPyEffectivePrestressGraphViewController::GetDuct() const
{
   return m_pMyViewController->GetDuct();
}

IEffectivePrestressGraphViewController::DuctType CPyEffectivePrestressGraphViewController::GetDuctType() const
{
   return m_pMyViewController->GetDuctType();
}

void CPyEffectivePrestressGraphViewController::SetShowGrid(bool bShow)
{
   m_pMyViewController->ShowGrid(bShow);
}

bool CPyEffectivePrestressGraphViewController::GetShowGrid() const
{
   return m_pMyViewController->ShowGrid();
}

void CPyEffectivePrestressGraphViewController::SetShowGirder(bool bShow)
{
   m_pMyViewController->ShowGirder(bShow);
}

bool CPyEffectivePrestressGraphViewController::GetShowGirder() const
{
   return m_pMyViewController->ShowGirder();
}

CPyEffectivePrestressGraphViewController* CPyEffectivePrestressGraphViewController::GetGraphController(CPyViewControllerBase* pController)
{
   auto pMyController = static_cast<CPyEffectivePrestressGraphViewController*>(pController);
   return pMyController;
}
