#include "stdafx.h"
#include "PyLoadsView.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

void CPyLoadsView::Init(ILoadsViewController* pViewController)
{
   m_pViewController = pViewController;
   CPyViewControllerBase::Init(pViewController);
}

void CPyLoadsView::SortBy(ILoadsViewController::Field field, ILoadsViewController::Direction direction)
{
   m_pViewController->SortBy(field, direction);
}

IndexType CPyLoadsView::GetLoadCount() const
{
   return m_pViewController->GetLoadCount();
}

std::_tstring CPyLoadsView::GetFieldValue(IndexType idx, ILoadsViewController::Field field) const
{
   return m_pViewController->GetFieldValue(idx, field);
}

void CPyLoadsView::DeleteLoad(IndexType idx)
{
   m_pViewController->DeleteLoad(idx);
}
