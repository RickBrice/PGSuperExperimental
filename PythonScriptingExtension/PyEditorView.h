#pragma once

#include "ScintillaDocView.h"

class CPyEditorView : public CScintillaView
{
public:
   CPyEditorView();
   ~CPyEditorView();
   afx_msg int OnCreate(LPCREATESTRUCT lpCreateStruct);

   DECLARE_MESSAGE_MAP()

   DECLARE_DYNCREATE(CPyEditorView)

protected:
   void SetAStyle(CScintillaCtrl& rCtrl, int style, COLORREF fore, COLORREF back = RGB(0xFF, 0xFF, 0xFF), int size = -1, const char* face = nullptr);
};
