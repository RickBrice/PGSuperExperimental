#pragma once


class CPyEditorDocTemplate : public CMultiDocTemplate
{
public:
   CPyEditorDocTemplate(UINT nIDResource, CRuntimeClass* pDocClass, CRuntimeClass* pFrameClass, CRuntimeClass* pViewClass);

   virtual CFrameWnd* CreateNewFrame(CDocument* pDoc, CFrameWnd* pOther);
   virtual CDocument* OpenDocumentFile(LPCTSTR lpszPathName, BOOL bMakeVisible = TRUE);
};
