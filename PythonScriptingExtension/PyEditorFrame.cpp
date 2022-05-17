#include "stdafx.h"
#include "Resource.h"
#include "PyEditorFrame.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

IMPLEMENT_DYNCREATE(CPyEditorFrame, CMDIChildWnd)

CPyEditorFrame::CPyEditorFrame()
{
}

CPyEditorFrame::~CPyEditorFrame()
{
}

BEGIN_MESSAGE_MAP(CPyEditorFrame, CMDIChildWnd)
END_MESSAGE_MAP()

BOOL CPyEditorFrame::LoadFrame(UINT nIDResource, DWORD dwDefaultStyle,
   CWnd* pParentWnd, CCreateContext* pContext)
{
   if (!__super::LoadFrame(nIDResource, dwDefaultStyle, pParentWnd, pContext))
   {
      return FALSE;
   }

   AFX_MANAGE_STATE(AfxGetStaticModuleState());
   CString strFullString;
   if (strFullString.LoadString(nIDResource))
   {
      AfxExtractSubString(m_strTitle, strFullString, 0);    // first sub-string
      SetWindowText(m_strTitle);
   }


   return TRUE;
}
