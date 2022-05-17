#pragma once

#include <IFace\Project.h>

class CPyRoadway
{
public:
   void Init(IBroker* pBroker);
   void Reset();

   void SetAlignmentData(const AlignmentData2& data);
   AlignmentData2 GetAlignmentData() const;

   void SetProfileData(const ProfileData2& data);
   ProfileData2 GetProfileData() const;

   void SetRoadwaySectionData(const RoadwaySectionData& data);
   RoadwaySectionData GetRoadwaySectionData() const;

private:
   CComPtr<IRoadwayData> m_pAlignment;
}; 
