#pragma once

#include "PyViewControllerBase.h"
#include <EAF\EAFViewController.h>
#include <Graphs\GirderPropertiesGraphViewController.h>

class CPyGirderPropertiesGraphViewController : public CPyViewControllerBase
{
public:
   void Init(IGirderPropertiesGraphViewController* pViewController);

   void SetPropertyType(CGirderPropertiesGraphBuilder::PropertyType propertyType);
   CGirderPropertiesGraphBuilder::PropertyType GetPropertyType() const;
   bool IsInvariantProperty(CGirderPropertiesGraphBuilder::PropertyType propertyType) const;

   void SetSectionPropertyType(pgsTypes::SectionPropertyType type);
   pgsTypes::SectionPropertyType GetSectionPropertyType() const;

   void SelectGirder(const CGirderKey& girderKey);
   const CGirderKey& GetGirder() const;

   void SetInterval(IntervalIndexType intervalIdx);
   IntervalIndexType GetInterval() const;
   IntervalIndexType GetFirstInterval() const;
   IntervalIndexType GetLastInterval() const;

   void SetShowGrid(bool bShow);
   bool GetShowGrid() const;

   void SetShowGirder(bool bShow);
   bool GetShowGirder() const;

   static CPyGirderPropertiesGraphViewController* GetGraphController(CPyViewControllerBase* pController);

private:
   CComPtr<IGirderPropertiesGraphViewController> m_pMyViewController;
};

