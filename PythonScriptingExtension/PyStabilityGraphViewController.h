#pragma once

#include "PyViewControllerBase.h"
#include <EAF\EAFViewController.h>
#include <Graphing\StabilityGraphViewController.h>

class CPyStabilityGraphViewController : public CPyViewControllerBase
{
public:
   void Init(IStabilityGraphViewController* pViewController);

   const CSegmentKey& GetSegment() const;
   void SelectSegment(const CSegmentKey& segmentKey);

   void SetViewMode(IStabilityGraphViewController::ViewMode mode);
   IStabilityGraphViewController::ViewMode GetViewMode() const;

   void SetShowGrid(bool bShow);
   bool GetShowGrid() const;

   static CPyStabilityGraphViewController* GetGraphController(CPyViewControllerBase* pController);

private:
   CComPtr<IStabilityGraphViewController> m_pMyViewController;
};

