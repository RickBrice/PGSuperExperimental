#pragma once

#include <EAF\EAFCommandLineInfo.h>

class CScriptingCommandLineInfo : public CEAFCommandLineInfo
{
public:
   CScriptingCommandLineInfo();

   UINT m_nParams;
   BOOL m_bRunScript;
   CString m_StartupScript;

   virtual void ParseParam(LPCTSTR lpszParam, BOOL bFlag, BOOL bLast) override;
   virtual CString GetUsageMessage() override;
};
