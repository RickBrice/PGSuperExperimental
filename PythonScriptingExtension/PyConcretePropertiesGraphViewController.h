#pragma once

#include "PyViewControllerBase.h"
#include <EAF\EAFViewController.h>
#include <Graphing\ConcretePropertiesGraphViewController.h>

class CPyConcretePropertiesGraphViewController : public CPyViewControllerBase
{
public:
   void Init(IConcretePropertiesGraphViewController* pViewController);

   void SetGraphType(IConcretePropertiesGraphViewController::GraphType type);
   IConcretePropertiesGraphViewController::GraphType GetGraphType() const;
   void SetElementType(IConcretePropertiesGraphViewController::ElementType type);
   IConcretePropertiesGraphViewController::ElementType GetElementType() const;
   void SetXAxisType(IConcretePropertiesGraphViewController::XAxisType type);
   IConcretePropertiesGraphViewController::XAxisType GetXAxisType() const;

   void SetSegment(const CSegmentKey& segmentKey);
   const CSegmentKey& GetSegment() const;
   void SetClosureJoint(const CClosureKey& closureKey);
   const CClosureKey& GetClosureJoint() const;

   void SetShowGrid(bool bShow);
   bool GetShowGrid() const;

   static CPyConcretePropertiesGraphViewController* GetGraphController(CPyViewControllerBase* pController);

private:
   CComPtr<IConcretePropertiesGraphViewController> m_pMyViewController;
};


