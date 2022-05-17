// StreamEdit.cpp : implementation file
//

#include "stdafx.h"
#include "StreamEdit.h"


// CStreamEdit

IMPLEMENT_DYNAMIC(CStreamEdit, CEdit)

CStreamEdit::CStreamEdit()
{

}

CStreamEdit::~CStreamEdit()
{
}


BEGIN_MESSAGE_MAP(CStreamEdit, CEdit)
END_MESSAGE_MAP()



// CStreamEdit message handlers


std::streamsize CStreamEdit::xsputn(const char *_Ptr, std::streamsize _Count)
{
   CString strText;
   GetWindowText(strText);
   CString newText(_Ptr);
   newText.Replace(_T("\n"), _T("\r\n"));
   strText += newText;
   SetWindowText(strText);
   UpdateWindow();
   return _Count;
}