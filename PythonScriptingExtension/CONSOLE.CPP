// Console.cpp : implementation file
//

#include "stdafx.h"
#include "resource.h"
#include "Console.h"

// CConsole

IMPLEMENT_DYNAMIC(CConsole, CDialogBar)

CConsole::CConsole() :
   CDialogBar(),
   m_strCommandPrompt(_T(">>> ")),
   m_strContinuePrompt(_T("... "))
{
   m_strCurrentPrompt = m_strCommandPrompt;
}

CConsole::~CConsole()
{
}

void CConsole::DoDataExchange(CDataExchange* pDX)
{
   CDialogBar::DoDataExchange(pDX);
   DDX_Control(pDX, IDC_COMMAND, m_cbCommand);
   DDX_Control(pDX, IDC_OUTPUT, m_wndOutput);
}

BEGIN_MESSAGE_MAP(CConsole, CDialogBar)
   ON_COMMAND(IDC_RUN, &CConsole::OnRun)
   ON_UPDATE_COMMAND_UI(IDC_RUN, &CConsole::OnUpdateRun)
   ON_MESSAGE(WM_INITDIALOG,OnInitDialog)
END_MESSAGE_MAP()



// CConsole message handlers
LRESULT CConsole::OnInitDialog(WPARAM wParam, LPARAM lParam)
{
   HandleInitDialog(wParam, lParam);
   UpdateData(FALSE);
   return FALSE;
}


void CConsole::OnUpdateRun(CCmdUI *pCmdUI)
{
   pCmdUI->Enable(TRUE);
}

void CConsole::OnRun()
{
   CString strCommand;
   m_cbCommand.GetWindowText(strCommand);

   HistoryAdd(strCommand);
   WriteCommand(m_strCurrentPrompt);
   WriteCommand(strCommand);
   WriteCommand(CString("\r\n"));

   //auto* old = std::cout.rdbuf();
   //std::cout.rdbuf(&m_wndOutput);

   USES_CONVERSION;
   int result = PyRun_SimpleString(OLE2A(strCommand));
   //m_pPython->RunCommand(strCommand);

   //std::cout.rdbuf(old);

   m_cbCommand.SetEditSel(0, -1);
   m_cbCommand.Clear();
}

void CConsole::HistoryAdd(const CString& strCommand)
{
   auto i = m_cbCommand.FindStringExact(-1, strCommand);
   if (i != CB_ERR)
   {
      m_cbCommand.DeleteString(i);
   }
   i = m_cbCommand.InsertString(0, strCommand);
   m_cbCommand.SetCurSel(i);
}

void CConsole::WriteCommand(const CString& strCommand)
{
   CString str;
   m_wndOutput.GetWindowText(str);
   str += strCommand;
   m_wndOutput.SetWindowText(str);
   m_wndOutput.UpdateWindow();
}
