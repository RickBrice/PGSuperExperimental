#include "stdafx.h"
#include "PyEditorView.h"
#include "SciLexer.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

//Taken from "lib\keyword.py" in the Python distribution
const TCHAR pyKeyWords[] = _T("False None True and as assert break class continue def del elif else except finally for from global if import in is lambda nonlocal not or pass raise return try while with yield");

#define SC_PYTHON_MARK_ENABLED_BREAKPOINT    0 
#define SC_PYTHON_MARK_DISABLED_BREAKPOINT   1
#define SC_PYTHON_MARK_INSTRUCTION_POINTER   2


IMPLEMENT_DYNCREATE(CPyEditorView,CScintillaView)

BEGIN_MESSAGE_MAP(CPyEditorView,CScintillaView)
   ON_WM_CREATE()
END_MESSAGE_MAP()

CPyEditorView::CPyEditorView()
{

}

CPyEditorView::~CPyEditorView()
{

}

int CPyEditorView::OnCreate(LPCREATESTRUCT lpCreateStruct)
{
   if (CScintillaView::OnCreate(lpCreateStruct) == -1)
      return -1;

   CScintillaCtrl& rCtrl = GetCtrl();

   //Setup the Lexer
   rCtrl.SetLexer(SCLEX_PYTHON);
   rCtrl.SetKeyWords(0, pyKeyWords);
   SetAStyle(rCtrl, STYLE_DEFAULT, RGB(0, 0, 0), RGB(0xff, 0xff, 0xff));
   rCtrl.StyleClearAll();
   SetAStyle(rCtrl, SCE_P_DEFAULT, RGB(0, 0, 0));
   SetAStyle(rCtrl, SCE_P_COMMENTLINE, RGB(0, 0x80, 0));
   SetAStyle(rCtrl, SCE_P_COMMENTBLOCK, RGB(0, 0x80, 0));
   SetAStyle(rCtrl, SCE_P_NUMBER, RGB(0, 0x80, 0x80));
   SetAStyle(rCtrl, SCE_P_WORD, RGB(0, 0, 0x80));
   rCtrl.StyleSetBold(SCE_P_WORD, 1);
   SetAStyle(rCtrl, SCE_P_WORD2, RGB(0, 0, 0x80));
   rCtrl.StyleSetBold(SCE_P_WORD2, 1);
   SetAStyle(rCtrl, SCE_P_STRING, RGB(0x80, 0, 0x80));
   SetAStyle(rCtrl, SCE_P_CHARACTER, RGB(0x80, 0, 0x80));
   SetAStyle(rCtrl, SCE_P_TRIPLE, RGB(0x80, 0, 0x80));
   SetAStyle(rCtrl, SCE_P_TRIPLEDOUBLE, RGB(0x80, 0, 0x80));
   SetAStyle(rCtrl, SCE_P_IDENTIFIER, RGB(0, 0, 0));
   SetAStyle(rCtrl, SCE_P_CLASSNAME, RGB(0, 0, 0));
   SetAStyle(rCtrl, SCE_P_DEFNAME, RGB(0, 0, 0));
   SetAStyle(rCtrl, SCE_P_OPERATOR, RGB(0x80, 0, 0));

   //Setup our markers
   rCtrl.SetMarginWidthN(1, 16);
   rCtrl.SetMarginTypeN(1, SC_MARGIN_SYMBOL);
   rCtrl.MarkerDefine(SC_PYTHON_MARK_ENABLED_BREAKPOINT, SC_MARK_CIRCLE);
   rCtrl.MarkerSetFore(SC_PYTHON_MARK_ENABLED_BREAKPOINT, RGB(255, 0, 0));
   rCtrl.MarkerSetBack(SC_PYTHON_MARK_ENABLED_BREAKPOINT, RGB(255, 0, 0));
   rCtrl.MarkerDefine(SC_PYTHON_MARK_DISABLED_BREAKPOINT, SC_MARK_CIRCLE);
   rCtrl.MarkerSetFore(SC_PYTHON_MARK_DISABLED_BREAKPOINT, RGB(255, 0, 0));
   rCtrl.MarkerDefine(SC_PYTHON_MARK_INSTRUCTION_POINTER, SC_MARK_SHORTARROW);
   rCtrl.MarkerSetFore(SC_PYTHON_MARK_INSTRUCTION_POINTER, RGB(0, 0, 0));
   rCtrl.MarkerSetBack(SC_PYTHON_MARK_INSTRUCTION_POINTER, RGB(255, 255, 0));

   //Setup the notification for modified text
   rCtrl.SetModEventMask(SC_MOD_INSERTTEXT | SC_MOD_DELETETEXT);

   return 0;
}


void CPyEditorView::SetAStyle(CScintillaCtrl& rCtrl, int style, COLORREF fore, COLORREF back, int size, const char* face)
{
   rCtrl.StyleSetFore(style, fore);
   rCtrl.StyleSetBack(style, back);
   if (size >= 1)
      rCtrl.StyleSetSize(style, size);
   if (face != nullptr)
      rCtrl.StyleSetFont(style, face);
}
