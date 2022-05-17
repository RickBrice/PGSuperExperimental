#include "stdafx.h"
#include "PyScriptEngine.h"

#include <Python.h>

#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif



CPyScriptEngine::CPyScriptEngine()
{
   AllocConsole();

   SetConsoleTitle(_T("Python"));

   m_hwndConsole = GetConsoleWindow();
   ShowConsole(FALSE);

   ioRedirect = std::make_unique<CPyIORedirector>();
}

CPyScriptEngine::~CPyScriptEngine()
{
   ioRedirect.reset();
   FreeConsole();
}

void CPyScriptEngine::ShowConsole(BOOL bShow)
{
   ::ShowWindow(m_hwndConsole, bShow ? SW_SHOW : SW_HIDE);
}

void CPyScriptEngine::Interactive()
{
   // still haven't figured out how to do an interactive console yet
   ShowConsole(TRUE);
   int result = PyRun_InteractiveLoop(stdin, "stdin");
   //m_PyIoRedirector.ShowConsole(FALSE);
   /*
   AllocConsole();

   //boost::python::object stdOut;
   //boost::python::object stdIn;
   //boost::python::object stdErr;
   //SetupRedirectIO(stdOut, stdIn, stdErr);

   HANDLE h = GetStdHandle(STD_INPUT_HANDLE);

   int result = PyRun_InteractiveLoopFlags((FILE*)h, "<stdin>", nullptr);
   //int result = Py_Main(__argc, __targv);

   //CloseIOFiles(stdErr, stdOut, stdIn);
   FreeConsole();

   //PyThreadState* pThreadState = Py_NewInterpreter();
   */
}

void CPyScriptEngine::Run(LPCTSTR lpszScriptFile)
{
   USES_CONVERSION;

   bool bKeepConsoleOpenOnError = true;

   FILE* pFile;
   if (_tfopen_s(&pFile, lpszScriptFile, _T("r")) == 0)
   {
      ShowConsole(TRUE);

      int result = PyRun_SimpleFile(pFile, T2A(lpszScriptFile));
      fclose(pFile);

      if (result == 0 || bKeepConsoleOpenOnError == false)
      {
         ShowConsole(FALSE);
      }
   }
   else
   {
#pragma Reminder("WORKING HERE - need error handling if script file cannot be opened")
   }
}

void CPyScriptEngine::Edit(LPCTSTR lpszScriptFile)
{
   AfxMessageBox(_T("Edit.... not implemented"));
}

void CPyScriptEngine::Debug(LPCTSTR lpszScriptFile)
{
   AfxMessageBox(_T("Debug.... not implemented"));
}
