#pragma once

#include "StreamEdit.h"

// CConsole

class CConsole : public CDialogBar
{
	DECLARE_DYNAMIC(CConsole)

public:
	CConsole();
	virtual ~CConsole();

protected:
	DECLARE_MESSAGE_MAP()

   CComboBox m_cbCommand;
   CStreamEdit m_wndOutput;


   CString m_strCurrentPrompt;
   CString m_strCommandPrompt;
   CString m_strContinuePrompt;

   void HistoryAdd(const CString& strCommand);

public:
   afx_msg void OnRun();
   afx_msg void OnUpdateRun(CCmdUI *pCmdUI);
   afx_msg LRESULT OnInitDialog(WPARAM, LPARAM);
   virtual void DoDataExchange(CDataExchange* pDX);

   void WriteCommand(const CString& strCommand);
};


