#pragma once

#include "PyViewControllerBase.h"
#include <EAF\EAFViewController.h>
#include <Graphs\DeflectionHistoryGraphViewController.h>

class CPyDeflectionHistoryGraphViewController : public CPyViewControllerBase
{
public:
   void Init(IDeflectionHistoryGraphViewController* pViewController);

   void SelectLocation(const pgsPointOfInterest& poi);
   const pgsPointOfInterest& GetLocation() const;

   void SetXAxisType(IDeflectionHistoryGraphViewController::XAxisType type);
   IDeflectionHistoryGraphViewController::XAxisType GetXAxisType() const;

   void SetIncludeElevationAdjustment(bool bAdjust);
   bool GetIncludeElevationAdjustment() const;

   void SetIncludeUnrecoverableDefl(bool bInclude);
   bool GetIncludeUnrecoverableDefl() const;

   void SetShowGrid(bool bShow);
   bool GetShowGrid() const;

   static CPyDeflectionHistoryGraphViewController* GetGraphController(CPyViewControllerBase* pController);

private:
   CComPtr<IDeflectionHistoryGraphViewController> m_pMyViewController;
};

