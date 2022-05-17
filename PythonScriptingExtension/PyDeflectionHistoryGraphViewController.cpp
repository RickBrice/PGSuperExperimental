#include "stdafx.h"
#include "PyDeflectionHistoryGraphViewController.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

void CPyDeflectionHistoryGraphViewController::Init(IDeflectionHistoryGraphViewController* pViewController)
{
   m_pMyViewController = pViewController;
   CPyViewControllerBase::Init(pViewController);
}

void CPyDeflectionHistoryGraphViewController::SelectLocation(const pgsPointOfInterest& poi)
{
   m_pMyViewController->SelectLocation(poi);
}

const pgsPointOfInterest& CPyDeflectionHistoryGraphViewController::GetLocation() const
{
   return m_pMyViewController->GetLocation();
}

void CPyDeflectionHistoryGraphViewController::SetXAxisType(IDeflectionHistoryGraphViewController::XAxisType type)
{
   m_pMyViewController->SetXAxisType(type);
}

IDeflectionHistoryGraphViewController::XAxisType CPyDeflectionHistoryGraphViewController::GetXAxisType() const
{
   return m_pMyViewController->GetXAxisType();
}

void CPyDeflectionHistoryGraphViewController::SetIncludeElevationAdjustment(bool bAdjust)
{
   m_pMyViewController->IncludeElevationAdjustment(bAdjust);
}

bool CPyDeflectionHistoryGraphViewController::GetIncludeElevationAdjustment() const
{
   return m_pMyViewController->IncludeElevationAdjustment();
}

void CPyDeflectionHistoryGraphViewController::SetIncludeUnrecoverableDefl(bool bInclude)
{
   m_pMyViewController->IncludeUnrecoverableDefl(bInclude);
}

bool CPyDeflectionHistoryGraphViewController::GetIncludeUnrecoverableDefl() const
{
   return m_pMyViewController->IncludeUnrecoverableDefl();
}

void CPyDeflectionHistoryGraphViewController::SetShowGrid(bool bShow)
{
   m_pMyViewController->ShowGrid(bShow);
}

bool CPyDeflectionHistoryGraphViewController::GetShowGrid() const
{
   return m_pMyViewController->ShowGrid();
}

CPyDeflectionHistoryGraphViewController* CPyDeflectionHistoryGraphViewController::GetGraphController(CPyViewControllerBase* pController)
{
   auto pMyController = static_cast<CPyDeflectionHistoryGraphViewController*>(pController);
   return pMyController;
}
