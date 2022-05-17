#include "stdafx.h"
#include "PyViewControllerBase.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

void CPyViewControllerBase::Init(IEAFViewController* pViewController)
{
   m_pViewController = pViewController;
}

void CPyViewControllerBase::Close()
{
   m_pViewController->Close();
}

bool CPyViewControllerBase::IsOpen() const
{
   return m_pViewController->IsOpen();
}

void CPyViewControllerBase::Minimize()
{
   m_pViewController->Minimize();
}

void CPyViewControllerBase::Maximize()
{
   m_pViewController->Maximize();
}

void CPyViewControllerBase::Restore()
{
   m_pViewController->Restore();
}

