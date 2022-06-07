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

// dllmain.h : Declaration of module class.

class CPGSuperScriptingExtensionModule : public CAtlDllModuleT< CPGSuperScriptingExtensionModule >
{
public :
//	DECLARE_LIBID(LIBID_PGSuperScriptingExtensionLib)
	DECLARE_REGISTRY_APPID_RESOURCEID(IDR_PGSUPERSCRIPTINGEXTENSION, "{4E7C50C2-59D4-4E54-892C-FB0428A54100}")
};

extern class CPGSuperScriptingExtensionModule _AtlModule;

#include "Bindings.h"
class CExtensionApp : public CWinApp
{
public:
   CPyPGSuper m_pyPGSuper;

   // Overrides
   virtual BOOL InitInstance() override;
   virtual int ExitInstance() override;

   CString GetScriptPath() const;

   DECLARE_MESSAGE_MAP()

private:
   HMODULE LoadLibraryFromApplicationDirectory(LPCTSTR lpFileName);
   void ReadSettings();
   void WriteSettings();

   HINSTANCE m_hSciDLL;
   CString m_ScriptPath;
};
