#include "stdafx.h"
#include "PyProjectProperties.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

void CPyProjectProperties::Init(IBroker* pBroker)
{
   pBroker->GetInterface(IID_IProjectProperties, (IUnknown**)&m_pProjectProperties);
}

void CPyProjectProperties::Reset()
{
   m_pProjectProperties.Release();
}

std::string CPyProjectProperties::GetBridgeName() const
{
   USES_CONVERSION;
   return T2A(m_pProjectProperties->GetBridgeName());
}

void CPyProjectProperties::SetBridgeName(LPCSTR strName)
{
   USES_CONVERSION;
   m_pProjectProperties->SetBridgeName(A2T(strName));
}

std::string CPyProjectProperties::GetBridgeID() const
{
   USES_CONVERSION;
   return T2A(m_pProjectProperties->GetBridgeID());
}

void CPyProjectProperties::SetBridgeID(LPCSTR strID)
{
   USES_CONVERSION;
   m_pProjectProperties->SetBridgeID(A2T(strID));
}

std::string CPyProjectProperties::GetJobNumber() const
{
   USES_CONVERSION;
   return T2A(m_pProjectProperties->GetJobNumber());
}

void CPyProjectProperties::SetJobNumber(LPCSTR strJobNumber)
{
   USES_CONVERSION;
   m_pProjectProperties->SetJobNumber(A2T(strJobNumber));
}

std::string CPyProjectProperties::GetEngineer() const
{
   USES_CONVERSION;
   return T2A(m_pProjectProperties->GetEngineer());
}

void CPyProjectProperties::SetEngineer(LPCSTR strEngineer)
{
   USES_CONVERSION;
   m_pProjectProperties->SetEngineer(A2T(strEngineer));
}

std::string CPyProjectProperties::GetCompany() const
{
   USES_CONVERSION;
   return T2A(m_pProjectProperties->GetCompany());
}

void CPyProjectProperties::SetCompany(LPCSTR strCompany)
{
   USES_CONVERSION;
   m_pProjectProperties->SetCompany(A2T(strCompany));
}

std::string CPyProjectProperties::GetComments() const
{
   USES_CONVERSION;
   return T2A(m_pProjectProperties->GetComments());
}

void CPyProjectProperties::SetComments(LPCSTR strComments)
{
   USES_CONVERSION;
   m_pProjectProperties->SetComments(A2T(strComments));
}
