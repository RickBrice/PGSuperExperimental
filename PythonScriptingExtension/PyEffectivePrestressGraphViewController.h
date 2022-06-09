#pragma once

#include "PyViewControllerBase.h"
#include <EAF\EAFViewController.h>
#include <Graphing\EffectivePrestressGraphViewController.h>

#include <boost\python\list.hpp>

class CPyEffectivePrestressGraphViewController : public CPyViewControllerBase
{
public:
   void Init(IEffectivePrestressGraphViewController* pViewController);

   IntervalIndexType GetMinimumInterval() const;
   IntervalIndexType GetMaximumInterval() const;

   void SelectInterval(IntervalIndexType intervalIdx);

   void SelectIntervals(const boost::python::list& listIntervals);
   boost::python::list GetSelectedIntervals() const;

   const CGirderKey& GetGirder() const;
   void SelectGirder(const CGirderKey& girderKey);

   void SetViewMode(IEffectivePrestressGraphViewController::ViewMode mode);
   IEffectivePrestressGraphViewController::ViewMode GetViewMode() const;

   void SetStrandType(IEffectivePrestressGraphViewController::StrandType strandType);
   IEffectivePrestressGraphViewController::StrandType GetStrandType() const;

   // set the duct index to -1 to show effective prestress of pretensioned strands
   void SetDuct(IEffectivePrestressGraphViewController::DuctType ductType, DuctIndexType ductIdx);
   DuctIndexType GetDuct() const;
   IEffectivePrestressGraphViewController::DuctType GetDuctType() const;

   void SetShowGrid(bool bShow);
   bool GetShowGrid() const;

   void SetShowGirder(bool bShow);
   bool GetShowGirder() const;

   static CPyEffectivePrestressGraphViewController* GetGraphController(CPyViewControllerBase* pController);

private:
   CComPtr<IEffectivePrestressGraphViewController> m_pMyViewController;
};

