#pragma once

#include "PyViewControllerBase.h"
#include <EAF\EAFViewController.h>
#include <Graphing\StressHistoryGraphViewController.h>

class CPyStressHistoryGraphViewController : public CPyViewControllerBase
{
public:
   void Init(IStressHistoryGraphViewController* pViewController);

   void SelectLocation(const pgsPointOfInterest& poi);
   const pgsPointOfInterest& GetLocation() const;

   void SetXAxisType(IStressHistoryGraphViewController::XAxisType type);
   IStressHistoryGraphViewController::XAxisType GetXAxisType() const;

   void SelectStressLocation(pgsTypes::StressLocation stressLocation, bool bEnable);
   bool IsStressLocationSelected(pgsTypes::StressLocation stressLocation) const;

   void SetShowGrid(bool bShow);
   bool GetShowGrid() const;

   static CPyStressHistoryGraphViewController* GetGraphController(CPyViewControllerBase* pController);

private:
   CComPtr<IStressHistoryGraphViewController> m_pMyViewController;
};

