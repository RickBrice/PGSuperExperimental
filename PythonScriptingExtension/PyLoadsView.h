#pragma once

#include <LoadsViewController.h>
#include "PyViewControllerBase.h"

class CPyLoadsView : public CPyViewControllerBase
{
public:
   void Init(ILoadsViewController* pViewController);

   void SortBy(ILoadsViewController::Field field, ILoadsViewController::Direction direction);
   IndexType GetLoadCount() const;
   std::_tstring GetFieldValue(IndexType idx, ILoadsViewController::Field field) const;
   void DeleteLoad(IndexType idx);

private:
   CComPtr<ILoadsViewController> m_pViewController;
};

