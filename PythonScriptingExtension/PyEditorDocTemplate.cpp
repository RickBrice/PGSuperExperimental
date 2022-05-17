
#include "stdafx.h"
#include "PyEditorDocTemplate.h"
#include <EAF\EAFUtilities.h>

#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

//////////////////////////////// Macros / Defines /////////////////////////////

#ifdef _DEBUG
#define new DEBUG_NEW
#endif #ifdef _DEBUG


//////////////////////////////// Implementation ///////////////////////////////

CPyEditorDocTemplate::CPyEditorDocTemplate(UINT nIDResource, CRuntimeClass* pDocClass, CRuntimeClass* pFrameClass, CRuntimeClass* pViewClass) :
   CMultiDocTemplate(nIDResource, pDocClass, pFrameClass, pViewClass)
{
}

CFrameWnd* CPyEditorDocTemplate::CreateNewFrame(CDocument* pDoc, CFrameWnd* pOther)
{
   //Validate our parameters
   if (pDoc != NULL)
      ASSERT_VALID(pDoc);
   // create a frame wired to the specified document

   ASSERT(m_nIDResource != 0); // must have a resource ID to load from
   CCreateContext context;
   context.m_pCurrentFrame = pOther;
   context.m_pCurrentDoc = pDoc;
   context.m_pNewViewClass = m_pViewClass;
   context.m_pNewDocTemplate = this;

   if (m_pFrameClass == NULL)
   {
      ASSERT(FALSE);
      return NULL;
   }
   CFrameWnd* pFrame = static_cast<CFrameWnd*>(m_pFrameClass->CreateObject());
   if (pFrame == NULL)
   {
      return NULL;
   }
   ASSERT_KINDOF(CFrameWnd, pFrame);


   //create new from resource
   if (!pFrame->LoadFrame(m_nIDResource,
      WS_OVERLAPPEDWINDOW | FWS_ADDTOTITLE,  //default frame styles
      EAFGetMainFrame(), &context))                //This is the critical change compared to CMultiDocTemplate::CreateNewFrame where we pass "m_pFrameWnd" instead of NULL
   {
      //frame will be deleted in PostNcDestroy cleanup
      return NULL;
   }

   //it worked !
   return pFrame;
}

CDocument* CPyEditorDocTemplate::OpenDocumentFile(LPCTSTR lpszPathName, BOOL bMakeVisible)
{
   return __super::OpenDocumentFile(lpszPathName, FALSE/*dont add to MRU*/, bMakeVisible);
}