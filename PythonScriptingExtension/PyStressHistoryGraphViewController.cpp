#include "stdafx.h"
#include "PyStressHistoryGraphViewController.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

void CPyStressHistoryGraphViewController::Init(IStressHistoryGraphViewController* pViewController)
{
   m_pMyViewController = pViewController;
   CPyViewControllerBase::Init(pViewController);
}

void CPyStressHistoryGraphViewController::SelectLocation(const pgsPointOfInterest& poi)
{
   m_pMyViewController->SelectLocation(poi);
}

const pgsPointOfInterest& CPyStressHistoryGraphViewController::GetLocation() const
{
   return m_pMyViewController->GetLocation();
}

void CPyStressHistoryGraphViewController::SetXAxisType(IStressHistoryGraphViewController::XAxisType type)
{
   m_pMyViewController->SetXAxisType(type);
}

IStressHistoryGraphViewController::XAxisType CPyStressHistoryGraphViewController::GetXAxisType() const
{
   return m_pMyViewController->GetXAxisType();
}

void CPyStressHistoryGraphViewController::SelectStressLocation(pgsTypes::StressLocation stressLocation, bool bEnable)
{
   m_pMyViewController->Stresses(stressLocation, bEnable);
}

bool CPyStressHistoryGraphViewController::IsStressLocationSelected(pgsTypes::StressLocation stressLocation) const
{
   return m_pMyViewController->Stresses(stressLocation);
}

void CPyStressHistoryGraphViewController::SetShowGrid(bool bShow)
{
   m_pMyViewController->ShowGrid(bShow);
}

bool CPyStressHistoryGraphViewController::GetShowGrid() const
{
   return m_pMyViewController->ShowGrid();
}

CPyStressHistoryGraphViewController* CPyStressHistoryGraphViewController::GetGraphController(CPyViewControllerBase* pController)
{
   auto pMyController = static_cast<CPyStressHistoryGraphViewController*>(pController);
   return pMyController;
}
