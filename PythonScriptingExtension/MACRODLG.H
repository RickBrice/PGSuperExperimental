#pragma once

class CPyScriptEngine;

// CMacroDlg dialog

class CMacroDlg : public CDialog
{
	DECLARE_DYNAMIC(CMacroDlg)

public:
	CMacroDlg(CWnd* pParent = NULL);   // standard constructor
	virtual ~CMacroDlg();

   void SetScriptEngine(CPyScriptEngine* pScriptEngine);

// Dialog Data
#ifdef AFX_DESIGN_TIME
	enum { IDD = IDD_MACROS };
#endif

protected:
   CPyScriptEngine* m_pScriptEngine;

	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV support

   virtual void OnCancel();

   void SaveWindowState();

	DECLARE_MESSAGE_MAP()
public:
   afx_msg void OnRun();
   afx_msg void OnDebug();
   afx_msg void OnEdit();
   afx_msg void OnCreate();
   afx_msg void OnDelete();
   virtual BOOL OnInitDialog();
   afx_msg void OnCbnDblclkMacros();
   afx_msg void OnClose();
   afx_msg void OnShowWindow(BOOL bShow, UINT nStatus);
};
