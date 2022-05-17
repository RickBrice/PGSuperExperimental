#pragma once

#include <IGraphManager.h>

class CPyGraphManager
{
public:
   void Init(IBroker* pBroker);
   void Reset();

   IndexType GetGraphBuilderCount();
   boost::python::list GetGraphNames();

private:
   CComPtr<IGraphManager> m_pGraphManager;
};