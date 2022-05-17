// MacroDlg.cpp : implementation file
//

#include "stdafx.h"
#include "MacroDlg.h"
#include "afxdialogex.h"
#include "resource.h"

#include <EAF\EAFUtilities.h>
#include <EAF\EAFApp.h>
#include "Scripting.h"

#include "PGSuperScriptingAgentImp.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif


// CMacroDlg dialog

IMPLEMENT_DYNAMIC(CMacroDlg, CDialog)

CMacroDlg::CMacroDlg(CWnd* pParent /*=NULL*/)
	: CDialog(IDD_MACROS, pParent)
{
   m_pScriptEngine = nullptr;
}

CMacroDlg::~CMacroDlg()
{
}

void CMacroDlg::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);
}

void CMacroDlg::SetScriptEngine(CPyScriptEngine* pScriptEngine)
{
   m_pScriptEngine = pScriptEngine;
}

BEGIN_MESSAGE_MAP(CMacroDlg, CDialog)
   ON_BN_CLICKED(IDC_RUN, &CMacroDlg::OnRun)
   ON_BN_CLICKED(IDC_DEBUG, &CMacroDlg::OnDebug)
   ON_BN_CLICKED(IDC_EDIT, &CMacroDlg::OnEdit)
   ON_BN_CLICKED(IDC_CREATE, &CMacroDlg::OnCreate)
   ON_BN_CLICKED(IDC_DELETE, &CMacroDlg::OnDelete)
   ON_CBN_DBLCLK(IDC_MACROS, &CMacroDlg::OnCbnDblclkMacros)
   ON_WM_CLOSE()
   ON_WM_SHOWWINDOW()
END_MESSAGE_MAP()


// CMacroDlg message handlers


BOOL CMacroDlg::OnInitDialog()
{
   CDialog::OnInitDialog();

   // TODO:  Add extra initialization here
   CComboBox* pCB = (CComboBox*)GetDlgItem(IDC_MACROS);

   CComPtr<IBroker> pBroker;
   EAFGetBroker(&pBroker);
   GET_IFACE2(pBroker, IScripting, pScripting);
   std::vector<CString> vScripts = pScripting->GetScriptNames();
   for (const auto& name : vScripts)
   {
      pCB->AddString(name);
   }

   return TRUE;  // return TRUE unless you set the focus to a control
                 // EXCEPTION: OCX Property Pages should return FALSE
}

void CMacroDlg::OnRun()
{
   // TODO: Add your control notification handler code here
   CComboBox* pCB = (CComboBox*)GetDlgItem(IDC_MACROS);
   CString strScript;
   pCB->GetLBText(pCB->GetCurSel(), strScript);

   CComPtr<IBroker> pBroker;
   EAFGetBroker(&pBroker);
   GET_IFACE2(pBroker, IScripting, pScripting);
   strScript = pScripting->GetScriptFile(strScript);

   m_pScriptEngine->Run(strScript);
}


void CMacroDlg::OnDebug()
{
   CComboBox* pCB = (CComboBox*)GetDlgItem(IDC_MACROS);
   int curSel = pCB->GetCurSel();

   CString strScript;
   pCB->GetLBText(curSel, strScript);

   CComPtr<IBroker> pBroker;
   EAFGetBroker(&pBroker);
   GET_IFACE2(pBroker, IScripting, pScripting);
   CString strScriptFile = pScripting->GetScriptFile(strScript);
   m_pScriptEngine->Debug(strScriptFile);
}


void CMacroDlg::OnEdit()
{
   CComboBox* pCB = (CComboBox*)GetDlgItem(IDC_MACROS);
   int curSel = pCB->GetCurSel();

   CString strScript;
   pCB->GetLBText(curSel, strScript);

   CComPtr<IBroker> pBroker;
   EAFGetBroker(&pBroker);
   GET_IFACE2(pBroker, IScripting, pScripting);
   CString strScriptFile = pScripting->GetScriptFile(strScript);
   m_pScriptEngine->Edit(strScriptFile);
}


void CMacroDlg::OnCreate()
{
   // TODO: Add your control notification handler code here
   CComboBox* pCB = (CComboBox*)GetDlgItem(IDC_MACROS);
   CString strScript;
   pCB->GetWindowText(strScript);

   CComPtr<IBroker> pBroker;
   EAFGetBroker(&pBroker);
   GET_IFACE2(pBroker, IScripting, pScripting);
   CString strScriptFile = pScripting->GetScriptFile(strScript);

   CFileFind finder;
   BOOL bFound = finder.FindFile(strScriptFile);
   if (bFound)
   {
      CString strError;
      strError.Format(_T("Could not create script file %s. The file already exists."), strScriptFile);
      AfxMessageBox(strError);
      return;
   }

   CFile f;
   CFileException e;
   if (!f.Open(strScriptFile, CFile::modeCreate | CFile::modeWrite, &e))
   {
      CString strError;
      TCHAR szError[255];
      e.GetErrorMessage(szError, 255);
      strError.Format(_T("Could not create script file %s (%s)"), strScriptFile, szError);
      AfxMessageBox(strError);
   }
   else
   {
      CString strBoilerPlate(_T("#Script\r\nimport BridgeLink\r\npgsuper = BridgeLink.GetPGSuper()\r\n#your script goes here"));
      f.Write(strBoilerPlate,strBoilerPlate.GetLength());
      f.Close();
      pCB->AddString(strScript);
      m_pScriptEngine->Edit(strScriptFile);
   }
}


void CMacroDlg::OnDelete()
{
   CComboBox* pCB = (CComboBox*)GetDlgItem(IDC_MACROS);
   int curSel = pCB->GetCurSel();

   CString strScript;
   pCB->GetLBText(curSel, strScript);

   CComPtr<IBroker> pBroker;
   EAFGetBroker(&pBroker);
   GET_IFACE2(pBroker, IScripting, pScripting);
   CString strScriptFile = pScripting->GetScriptFile(strScript);

   CString strMsg;
   strMsg.Format(_T("Are you sure you want to delete %s.\nThis action cannot be undone."), strScript);
   if (AfxMessageBox(strMsg,MB_YESNOCANCEL) == IDYES)
   {
      try
      {
         CFile::Remove(strScriptFile);
         pCB->DeleteString(curSel);
      }
      catch (CFileException* pEx)
      {
         TCHAR szError[255];
         pEx->GetErrorMessage(szError, 255);

         CString strMsg;
         strMsg.Format(_T("Cannot delete %s (%s)"), strScriptFile, szError);
         AfxMessageBox(strMsg);
         pEx->Delete();
      };
   }

}

void CMacroDlg::OnCbnDblclkMacros()
{
   OnRun();
}

void CMacroDlg::OnCancel()
{
   SaveWindowState();
   __super::OnCancel();
}

void CMacroDlg::OnClose()
{
   SaveWindowState();
   CDialog::OnClose();
}

void CMacroDlg::SaveWindowState()
{
   CEAFApp* pApp = EAFGetApp();

   // Save the layout of the application window
   WINDOWPLACEMENT wp;
   wp.length = sizeof wp;
   if (GetWindowPlacement(&wp))
   {
      pApp->WriteWindowPlacement(_T("Window Positions"), _T("Macros"), &wp);
   }
}

void CMacroDlg::OnShowWindow(BOOL bShow, UINT nStatus)
{
   CDialog::OnShowWindow(bShow, nStatus);

   // Restore the layout of the application window
   if (bShow)
   {
      WINDOWPLACEMENT wp;
      CEAFApp* pApp = EAFGetApp();
      if (pApp->ReadWindowPlacement(_T("Window Positions"), _T("Macros"), &wp))
      {
         SetWindowPos(NULL, wp.rcNormalPosition.left, wp.rcNormalPosition.top, wp.rcNormalPosition.right - wp.rcNormalPosition.left, wp.rcNormalPosition.bottom - wp.rcNormalPosition.top, 0);
      }
   }
}
