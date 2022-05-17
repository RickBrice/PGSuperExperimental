#include "stdafx.h"
#include "PyReportManager.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

void CPyReportManager::Init(IBroker* pBroker)
{
   pBroker->GetInterface(IID_IReportManager, (IUnknown**)&m_pReportManager);
}

void CPyReportManager::Reset()
{
   m_pReportManager.Release();
}

IndexType CPyReportManager::GetReportBuilderCount()
{
   return m_pReportManager->GetReportBuilderCount();
}

boost::python::list CPyReportManager::GetReportNames()
{
   USES_CONVERSION;
   boost::python::list list;
   std::vector<std::_tstring> vNames = m_pReportManager->GetReportNames();
   for (const auto& name : vNames)
   {
#if defined _UNICODE
      boost::python::handle<> handle(PyUnicode_FromWideChar(name.c_str(),-1));
#else
      boost::python::handle<> handle(PyUnicode_FromString(T2A(name.c_str())));
#endif
      list.append(boost::python::object(handle));
   }

   return list;
}