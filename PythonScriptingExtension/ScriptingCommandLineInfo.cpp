#include "stdafx.h"
#include "ScriptingCommandLineInfo.h"


#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

CScriptingCommandLineInfo::CScriptingCommandLineInfo() :
   CEAFCommandLineInfo(),
   m_nParams(0),
   m_bRunScript(FALSE)
{
}

void CScriptingCommandLineInfo::ParseParam(LPCTSTR lpszParam, BOOL bFlag, BOOL bLast)
{
   CString strParam(lpszParam);

   CCommandLineInfo::ParseParam(lpszParam, bFlag, bLast);

   m_nParams++;

   if (strParam.Left(6).CompareNoCase(_T("script")) == 0)
   {
      m_bRunScript = TRUE;
      m_StartupScript = strParam.Right(strParam.GetLength() - 7);
   }
}

CString CScriptingCommandLineInfo::GetUsageMessage()
{
   return CString(_T("/script=<script file>"));
}
