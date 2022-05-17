#pragma once

#include <EAF\EAFViewController.h>

class CPyViewControllerBase
{
public:
   void Init(IEAFViewController* pViewController);

   void Close();
   bool IsOpen() const;

   void Minimize();
   void Maximize();
   void Restore();

private:
   CComPtr<IEAFViewController> m_pViewController;
};

