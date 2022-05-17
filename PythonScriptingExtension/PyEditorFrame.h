#pragma once

class CPyEditorFrame : public CMDIChildWnd
{
public:
   DECLARE_DYNCREATE(CPyEditorFrame)

   virtual BOOL LoadFrame(UINT nIDResource, DWORD dwDefaultStyle,
      CWnd* pParentWnd, CCreateContext* pContext = NULL);

protected:
   CPyEditorFrame(void);
   ~CPyEditorFrame(void);

   DECLARE_MESSAGE_MAP()
};
