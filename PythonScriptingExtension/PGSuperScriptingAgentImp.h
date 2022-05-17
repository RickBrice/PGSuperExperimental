///////////////////////////////////////////////////////////////////////
// ExtensionAgentExample - Extension Agent Example Project for PGSuper
// Copyright © 1999-2017  Washington State Department of Transportation
//                        Bridge and Structures Office
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the Alternate Route Open Source License as 
// published by the Washington State Department of Transportation, 
// Bridge and Structures Office.
//
// This program is distributed in the hope that it will be useful, but 
// distribution is AS IS, WITHOUT ANY WARRANTY; without even the implied 
// warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See 
// the Alternate Route Open Source License for more details.
//
// You should have received a copy of the Alternate Route Open Source 
// License along with this program; if not, write to the Washington 
// State Department of Transportation, Bridge and Structures Office, 
// P.O. Box  47340, Olympia, WA 98503, USA or e-mail 
// Bridge_Support@wsdot.wa.gov
///////////////////////////////////////////////////////////////////////

// ExampleExtensionAgent.h : Declaration of the CExampleExtensionAgent

#pragma once
#include "resource.h"       // main symbols
#include <EAF\EAFUIIntegration.h>
#include <EAF\EAFInterfaceCache.h>

#include "Scripting.h"
#include "MacroDlg.h"

#include "Console.h"

#include "PyEditorDoc.h"
#include "PyEditorFrame.h"
#include "PyEditorView.h"

#include "PyScriptEngine.h"



// {9AD09682-03D8-42CF-B17A-6E679ABFD880}
DEFINE_GUID(CLSID_PGSuperScriptingAgent,
   0x9AD09682, 0x03D8, 0x42CF, 0xB1, 0x7a, 0x6e, 0x67, 0x9a, 0xbf, 0xd8, 0x80);

#if defined(_WIN32_WCE) && !defined(_CE_DCOM) && !defined(_CE_ALLOW_SINGLE_THREADED_OBJECTS_IN_MTA)
#error "Single-threaded COM objects are not properly supported on Windows CE platform, such as the Windows Mobile platforms that do not include full DCOM support. Define _CE_ALLOW_SINGLE_THREADED_OBJECTS_IN_MTA to force ATL to support creating single-thread COM object's and allow use of it's single-threaded COM object implementations. The threading model in your rgs file was set to 'Free' as that is the only threading model supported in non DCOM Windows CE platforms."
#endif

class CPGSuperScriptingAgent;
class CScintillaCtrl;

class CMyCmdTarget : public CCmdTarget
{
public:
   CMyCmdTarget() { m_pMyAgent = nullptr; };

   void Init(CPGSuperScriptingAgent* pAgent);
   void Reset();

   afx_msg void OnMacros();
   afx_msg void OnInteractive();
   afx_msg void OnShowConsole();
   afx_msg void OnUpdateShowConsole(CCmdUI *pCmdUI);

   CPGSuperScriptingAgent* m_pMyAgent;

   CMacroDlg m_MacroDlg;

   DECLARE_MESSAGE_MAP()
};


// CExtensionAgent

class ATL_NO_VTABLE CPGSuperScriptingAgent :
	public CComObjectRootEx<CComSingleThreadModel>,
	public CComCoClass<CPGSuperScriptingAgent, &CLSID_PGSuperScriptingAgent>,
   public IEAFCommandCallback,
	public IAgentEx,
   //public IAgentPersist,
   public IAgentUIIntegration,
   public IEAFProcessCommandLine, // we add some command line options for launching scripts
   public IScripting
{
public:
   CPGSuperScriptingAgent()
	{
	}

DECLARE_REGISTRY_RESOURCEID(IDR_PGSUPERSCRIPTINGAGENT)

DECLARE_NOT_AGGREGATABLE(CPGSuperScriptingAgent)

BEGIN_COM_MAP(CPGSuperScriptingAgent)
	COM_INTERFACE_ENTRY(IAgent)
	COM_INTERFACE_ENTRY(IAgentEx)
	//COM_INTERFACE_ENTRY(IAgentPersist)
   COM_INTERFACE_ENTRY(IAgentUIIntegration)
   COM_INTERFACE_ENTRY(IEAFProcessCommandLine)
   COM_INTERFACE_ENTRY(IScripting)
END_COM_MAP()



	DECLARE_PROTECT_FINAL_CONSTRUCT()
   DECLARE_EAF_AGENT_DATA;

	HRESULT FinalConstruct();
	void FinalRelease()
	{
	}

   CMyCmdTarget m_MyCommandTarget;

   DWORD m_dwExtendUICookie;

   CBitmap m_bmpMenu;

   CEAFMenu* m_pMyMenu;
   void CreateMenus();
   void RemoveMenus();

   UINT m_ToolBarID;
   void CreateToolBar();
   void RemoveToolBar();

   void RegisterViews();
   void UnregisterViews();

   CConsole m_Console;
   CPyScriptEngine ScriptEngine;

   void EditScript(LPCTSTR lpszScriptFile);

// IAgentEx
public:
   STDMETHOD(SetBroker)(IBroker* pBroker) override;
   STDMETHOD(RegInterfaces)() override;
   STDMETHOD(Init)() override;
   STDMETHOD(Init2)() override;
   STDMETHOD(Reset)() override;
   STDMETHOD(ShutDown)() override;
   STDMETHOD(GetClassID)(CLSID* pCLSID) override;

//// IAgentPersist
//public:
//   STDMETHOD(Load)(/*[in]*/ IStructuredLoad* pStrLoad) override;
//   STDMETHOD(Save)(/*[in]*/ IStructuredSave* pStrSave) override;

// IAgentUIIntegration
public:
   STDMETHOD(IntegrateWithUI)(BOOL bIntegrate) override;

// IEAFCommandCallback
public:
   virtual BOOL OnCommandMessage(UINT nID,int nCode,void* pExtra,AFX_CMDHANDLERINFO* pHandlerInfo) override;
   virtual BOOL GetStatusBarMessageString(UINT nID, CString& rMessage) const override;
   virtual BOOL GetToolTipMessageString(UINT nID, CString& rMessage) const override;

// IEAFProcessCommandLine
public:
   virtual BOOL ProcessCommandLineOptions(CEAFCommandLineInfo& cmdInfo) override;

// IScripting
public:
   virtual LPCTSTR GetScriptFolder() const override;
   virtual std::vector<CString> GetScriptNames() const override;
   virtual std::vector<CString> GetScriptFiles() const override;
   virtual CString GetScriptFile(LPCTSTR lpszScriptName) const override;

private:
   void RegisterUIExtensions();
   void UnregisterUIExtensions();

   // the script editor view is not a view of the PGSuper document
   //long m_SourceEditorViewKey;
   //void CreateScriptEditorView();

   CPyEditorDoc m_ScriptDocument;
   CPyEditorFrame* m_pFrame;
   CPyEditorView* m_pView;
};

OBJECT_ENTRY_AUTO(CLSID_PGSuperScriptingAgent, CPGSuperScriptingAgent)
