#include "stdafx.h"
#include "PyGraphManager.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

void CPyGraphManager::Init(IBroker* pBroker)
{
   pBroker->GetInterface(IID_IGraphManager, (IUnknown**)&m_pGraphManager);
}

void CPyGraphManager::Reset()
{
   m_pGraphManager.Release();
}

IndexType CPyGraphManager::GetGraphBuilderCount()
{
   return m_pGraphManager->GetGraphBuilderCount();
}

boost::python::list CPyGraphManager::GetGraphNames()
{
   USES_CONVERSION;
   boost::python::list list;
   std::vector<std::_tstring> vNames = m_pGraphManager->GetGraphNames();
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
