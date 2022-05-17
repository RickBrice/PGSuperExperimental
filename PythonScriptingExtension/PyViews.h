#pragma once

#include <IFace\Views.h>

#include "PyBridgeModelView.h"
#include "PyGirderModelView.h"
#include "PyLoadsView.h"

class CPyViews
{
public:
   void Init(IBroker* pBroker);
   void Reset();

   CPyBridgeModelView CreateBridgeModelView();
   CPyGirderModelView CreateGirderView();
   CPyLoadsView CreateLoadsView();
   void CreateLibraryEditorView();
   void CreateGraphView(IndexType idx);
   CPyViewControllerBase* CreateGraphViewByName(LPCSTR strName);
   void CreateReport(IndexType idx,bool bPrompt);

private:
   IBroker* m_pBroker; // weak reference
   CComPtr<IViews> m_pViews;

   CPyViewControllerBase* CreateGraphContollerWrapper(IEAFViewController* pController) const;
};
