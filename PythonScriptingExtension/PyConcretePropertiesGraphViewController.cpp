#include "stdafx.h"
#include "PyConcretePropertiesGraphViewController.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

void CPyConcretePropertiesGraphViewController::Init(IConcretePropertiesGraphViewController* pViewController)
{
   m_pMyViewController = pViewController;
   CPyViewControllerBase::Init(pViewController);
}

void CPyConcretePropertiesGraphViewController::SetGraphType(IConcretePropertiesGraphViewController::GraphType type)
{
   m_pMyViewController->SetGraphType(type);
}

IConcretePropertiesGraphViewController::GraphType CPyConcretePropertiesGraphViewController::GetGraphType() const
{
   return m_pMyViewController->GetGraphType();
}

void CPyConcretePropertiesGraphViewController::SetElementType(IConcretePropertiesGraphViewController::ElementType type)
{
   m_pMyViewController->SetElementType(type);
}

IConcretePropertiesGraphViewController::ElementType CPyConcretePropertiesGraphViewController::GetElementType() const
{
   return m_pMyViewController->GetElementType();
}

void CPyConcretePropertiesGraphViewController::SetXAxisType(IConcretePropertiesGraphViewController::XAxisType type)
{
   m_pMyViewController->SetXAxisType(type);
}

IConcretePropertiesGraphViewController::XAxisType CPyConcretePropertiesGraphViewController::GetXAxisType() const
{
   return m_pMyViewController->GetXAxisType();
}

void CPyConcretePropertiesGraphViewController::SetSegment(const CSegmentKey& segmentKey)
{
   m_pMyViewController->SetSegment(segmentKey);
}

const CSegmentKey& CPyConcretePropertiesGraphViewController::GetSegment() const
{
   return m_pMyViewController->GetSegment();
}

void CPyConcretePropertiesGraphViewController::SetClosureJoint(const CClosureKey& closureKey)
{
   m_pMyViewController->SetClosureJoint(closureKey);
}

const CClosureKey& CPyConcretePropertiesGraphViewController::GetClosureJoint() const
{
   return m_pMyViewController->GetClosureJoint();
}

void CPyConcretePropertiesGraphViewController::SetShowGrid(bool bShow)
{
   m_pMyViewController->ShowGrid(bShow);
}

bool CPyConcretePropertiesGraphViewController::GetShowGrid() const
{
   return m_pMyViewController->ShowGrid();
}

CPyConcretePropertiesGraphViewController* CPyConcretePropertiesGraphViewController::GetGraphController(CPyViewControllerBase* pController)
{
   auto pMyController = static_cast<CPyConcretePropertiesGraphViewController*>(pController);
   return pMyController;
}
