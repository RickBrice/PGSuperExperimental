#include "stdafx.h"
#include "PyGirderPropertiesGraphViewController.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

void CPyGirderPropertiesGraphViewController::Init(IGirderPropertiesGraphViewController* pViewController)
{
   m_pMyViewController = pViewController;
   CPyViewControllerBase::Init(pViewController);
}

void CPyGirderPropertiesGraphViewController::SetPropertyType(CGirderPropertiesGraphBuilder::PropertyType propertyType)
{
   if (!m_pMyViewController->SetPropertyType(propertyType))
   {
      PyErr_SetString(PyExc_TypeError, "Invalid girder property type");
   }
}

CGirderPropertiesGraphBuilder::PropertyType CPyGirderPropertiesGraphViewController::GetPropertyType() const
{
   return m_pMyViewController->GetPropertyType();
}

bool CPyGirderPropertiesGraphViewController::IsInvariantProperty(CGirderPropertiesGraphBuilder::PropertyType propertyType) const
{
   return m_pMyViewController->IsInvariantProperty(propertyType);
}

void CPyGirderPropertiesGraphViewController::SetSectionPropertyType(pgsTypes::SectionPropertyType type)
{
   if (!m_pMyViewController->SetSectionPropertyType(type))
   {
      PyErr_SetString(PyExc_TypeError, "Invalid section property type");
   }
}

pgsTypes::SectionPropertyType CPyGirderPropertiesGraphViewController::GetSectionPropertyType() const
{
   return m_pMyViewController->GetSectionPropertyType();
}

const CGirderKey& CPyGirderPropertiesGraphViewController::GetGirder() const
{
   return m_pMyViewController->GetGirder();
}

void CPyGirderPropertiesGraphViewController::SelectGirder(const CGirderKey& girderKey)
{
   m_pMyViewController->SelectGirder(girderKey);
}

void CPyGirderPropertiesGraphViewController::SetInterval(IntervalIndexType intervalIdx)
{
   m_pMyViewController->SetInterval(intervalIdx);
}

IntervalIndexType CPyGirderPropertiesGraphViewController::GetInterval() const
{
   return m_pMyViewController->GetInterval();
}

IntervalIndexType CPyGirderPropertiesGraphViewController::GetFirstInterval() const
{
   return m_pMyViewController->GetFirstInterval();
}

IntervalIndexType CPyGirderPropertiesGraphViewController::GetLastInterval() const
{
   return m_pMyViewController->GetLastInterval();
}

void CPyGirderPropertiesGraphViewController::SetShowGrid(bool bShow)
{
   m_pMyViewController->ShowGrid(bShow);
}

bool CPyGirderPropertiesGraphViewController::GetShowGrid() const
{
   return m_pMyViewController->ShowGrid();
}

void CPyGirderPropertiesGraphViewController::SetShowGirder(bool bShow)
{
   m_pMyViewController->ShowGirder(bShow);
}

bool CPyGirderPropertiesGraphViewController::GetShowGirder() const
{
   return m_pMyViewController->ShowGirder();
}

CPyGirderPropertiesGraphViewController* CPyGirderPropertiesGraphViewController::GetGraphController(CPyViewControllerBase* pController)
{
   auto pMyController = static_cast<CPyGirderPropertiesGraphViewController*>(pController);
   return pMyController;
}
