#include "stdafx.h"
#include "PyRoadway.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

void CPyRoadway::Init(IBroker* pBroker)
{
   pBroker->GetInterface(IID_IRoadwayData, (IUnknown**)&m_pAlignment);
}

void CPyRoadway::Reset()
{
   m_pAlignment.Release();
}

void CPyRoadway::SetAlignmentData(const AlignmentData2& data)
{
   m_pAlignment->SetAlignmentData2(data);
}

AlignmentData2 CPyRoadway::GetAlignmentData() const
{
   return m_pAlignment->GetAlignmentData2();
}

void CPyRoadway::SetProfileData(const ProfileData2& data)
{
   m_pAlignment->SetProfileData2(data);
}

ProfileData2 CPyRoadway::GetProfileData() const
{
   return m_pAlignment->GetProfileData2();
}

void CPyRoadway::SetRoadwaySectionData(const RoadwaySectionData& data)
{
   m_pAlignment->SetRoadwaySectionData(data);
}

RoadwaySectionData CPyRoadway::GetRoadwaySectionData() const
{
   return m_pAlignment->GetRoadwaySectionData();
}

