#include "stdafx.h"
#include "PyIORedirector.h"

//#include <fcntl.h> // for _O_TEXT

#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

static const WORD MAX_CONSOLE_LINES = 500;
CPyIORedirector::CPyIORedirector()
{
   SetupRedirectIO(m_pyStdOut, m_pyStdIn, m_pyStdErr);

   //// Get STDOUT handle
   //HANDLE ConsoleOutput = GetStdHandle(STD_OUTPUT_HANDLE);
   //int SystemOutput = _open_osfhandle(intptr_t(ConsoleOutput), _O_TEXT);
   //FILE *COutputHandle = _fdopen(SystemOutput, "w");

   //// Get STDERR handle
   //HANDLE ConsoleError = GetStdHandle(STD_ERROR_HANDLE);
   //int SystemError = _open_osfhandle(intptr_t(ConsoleError), _O_TEXT);
   //FILE *CErrorHandle = _fdopen(SystemError, "w");

   //// Get STDIN handle
   //HANDLE ConsoleInput = GetStdHandle(STD_INPUT_HANDLE);
   //int SystemInput = _open_osfhandle(intptr_t(ConsoleInput), _O_TEXT);
   //FILE *CInputHandle = _fdopen(SystemInput, "r");

   ////make cout, wcout, cin, wcin, wcerr, cerr, wclog and clog point to console as well
   //std::ios::sync_with_stdio(true);

   //// Redirect the CRT standard input, output, and error handles to the console
   //freopen_s(&CInputHandle, "CONIN$", "r", stdin);
   //freopen_s(&COutputHandle, "CONOUT$", "w", stdout);
   //freopen_s(&CErrorHandle, "CONOUT$", "w", stderr);

   ////Clear the error state for each of the C++ standard stream objects. We need to do this, as
   ////attempts to access the standard streams before they refer to a valid target will cause the
   ////iostream objects to enter an error state. In versions of Visual Studio after 2005, this seems
   ////to always occur during startup regardless of whether anything has been read from or written to
   ////the console or not.
   //std::wcout.clear();
   //std::cout.clear();
   //std::wcerr.clear();
   //std::cerr.clear();
   //std::wcin.clear();
   //std::cin.clear();

}

CPyIORedirector::~CPyIORedirector()
{
   CloseIOFiles(m_pyStdErr, m_pyStdOut, m_pyStdIn);
}

void CPyIORedirector::SetupRedirectIO(boost::python::object& stdOut, boost::python::object& stdIn, boost::python::object &stdErr)
{
   boost::python::object io = boost::python::import(boost::python::str("io"));
   boost::python::object sys = boost::python::import(boost::python::str("sys"));

   RedirectIO(io, sys, "stdout", "CONOUT$", "wt", stdOut);
   RedirectIO(io, sys, "stdin", "CONIN$", "rt", stdIn);
   //RedirectIO(io, sys, "stderr", "CONERR$", "wt", stdErr); // this just goes to a file named CONERR$
   RedirectIO(io, sys, "stderr", "CONOUT$", "wt", stdErr);
}

bool CPyIORedirector::RedirectIO(boost::python::object& io, boost::python::object& sys, LPCSTR pszPythonName, LPCSTR pszName, LPCSTR pszMode, boost::python::object& fileObject)
{
   //Validate our parameters
   ASSERT(io.ptr() != nullptr);
   ASSERT(sys.ptr() != nullptr);

   fileObject = io.attr("open")(pszName, pszMode);
   if (fileObject.ptr() == nullptr)
   {
      //   AddOutputText(IDS_PYDEBUG_FAIL_OPEN_FILE_HANDLE, GetPythonErrorString());
      return false;
   }

   boost::python::setattr(sys, boost::python::str(pszPythonName), fileObject);
   return true;
}


void CPyIORedirector::ClosePythonFile(boost::python::object& file)
{
   //Validate our parameters
   ASSERT(file.ptr() != nullptr);

   boost::python::handle<> result = boost::python::handle<>(PyObject_CallMethod(file.ptr(), "close", ""));
}

void CPyIORedirector::CloseIOFiles(boost::python::object& stdErr, boost::python::object& stdOut, boost::python::object& stdIn)
{
   if (stdErr.ptr() != nullptr)
      ClosePythonFile(stdErr);

   if (stdIn.ptr() != nullptr)
      ClosePythonFile(stdIn);

   if (stdOut.ptr() != nullptr)
      ClosePythonFile(stdOut);
}

