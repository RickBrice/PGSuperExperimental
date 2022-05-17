#pragma once

#include <IReportManager.h>

class CPyReportManager
{
public:
   void Init(IBroker* pBroker);
   void Reset();

   IndexType GetReportBuilderCount();
   boost::python::list GetReportNames();

private:
   CComPtr<IReportManager> m_pReportManager;
};