#pragma once

class CPyIORedirector
{
public:
   CPyIORedirector();
   ~CPyIORedirector();

private:
   //FILE pstdout, pstdin, pstderr;
   //intptr_t stdInHandle, stdOutHandle, stdErrHandle;
   //int stdInConHandle, stdOutConHandle, stdErrConHandle;

   boost::python::object m_pyStdOut;
   boost::python::object m_pyStdIn;
   boost::python::object m_pyStdErr;

   bool RedirectIO(boost::python::object& ioModule, boost::python::object& sysModule, LPCSTR pszPythonName, LPCSTR pszName, LPCSTR pszMode, boost::python::object& fileObject);
   void SetupRedirectIO(boost::python::object& stdOut, boost::python::object& stdIn, boost::python::object &stdErr);
   void ClosePythonFile(boost::python::object& file);
   void CloseIOFiles(boost::python::object& stdErr, boost::python::object& stdOut, boost::python::object& stdIn);
};

