#pragma once
#include "PyIORedirector.h"

class CPyScriptEngine
{
public:
   CPyScriptEngine();
   ~CPyScriptEngine();

   void Interactive();
   void Run(LPCTSTR lpszScriptFile);
   void Edit(LPCTSTR lpszScriptFile);
   void Debug(LPCTSTR lpszScriptFile);

private:
   std::unique_ptr<CPyIORedirector> ioRedirect;

   void ShowConsole(BOOL bShow);
   HWND m_hwndConsole;
};
