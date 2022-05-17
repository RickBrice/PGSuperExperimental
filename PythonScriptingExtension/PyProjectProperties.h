#pragma once

#include <IFace\Project.h>

class CPyProjectProperties
{
public:
   void Init(IBroker* pBroker);
   void Reset();

   std::string GetBridgeName() const;
   void SetBridgeName(LPCSTR strName);

   std::string GetBridgeID() const;
   void SetBridgeID(LPCSTR strID);

   std::string GetJobNumber() const;
   void SetJobNumber(LPCSTR strJobNumber);

   std::string GetEngineer() const;
   void SetEngineer(LPCSTR strEngineer);

   std::string GetCompany() const;
   void SetCompany(LPCSTR strCompany);

   std::string GetComments() const;
   void SetComments(LPCSTR strComments);

private:
   CComPtr<IProjectProperties> m_pProjectProperties;
}; 
