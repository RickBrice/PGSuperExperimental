#include "stdafx.h"
#include "ExceptionHandlers.h"

void handle_unwind_exception(CXUnwind* pException)
{
   USES_CONVERSION;
   std::_tstring strMsg;
   pException->GetErrorMessage(&strMsg);
   PyErr_SetString(PyExc_RuntimeError, T2A(strMsg.c_str()));
   //PyErr_Print();
   //throw pException;
}
