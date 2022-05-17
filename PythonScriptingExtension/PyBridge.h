#pragma once

#include <IFace\Bridge.h>

class CPyBridge
{
public:
   void Init(IBroker* pBroker);
   void Reset();

   bool HasAsymmetricGirders() const;
   bool HasTiltedGirders() const;
   Float64 GetAlignmentOffset() const;
   Float64 GetLength() const;
   Float64 GetNominalSpanLength(SpanIndexType spanIdx) const;
   SpanIndexType GetSpanCount() const;
   PierIndexType GetPierCount() const;
   SupportIndexType GetTemporarySupportCount() const;
   GroupIndexType GetGirderGroupCount() const;
   GirderIndexType GetGirderlineCount() const;
   GirderIndexType GetGirderCount(GroupIndexType grpIdx) const;
   GirderIndexType GetGirderCountBySpan(SpanIndexType spanIdx) const;
   PierIndexType GetGirderGroupStartPier(GroupIndexType grpIdx) const;
   PierIndexType GetGirderGroupEndPier(GroupIndexType grpIdx) const;
   SpanIndexType GetGirderGroupStartSpan(GroupIndexType grpIdx) const;
   SpanIndexType GetGirderGroupEndSpan(GroupIndexType grpIdx) const;
   GroupIndexType GetGirderGroupIndex(SpanIndexType spanIdx) const;

   Float64 GetSpanLength(const CSpanKey& spanKey) const;
   Float64 GetFullSpanLength(const CSpanKey& spanKey) const;
   Float64 GetGirderlineLength(GirderIndexType gdrLineIdx) const;
   Float64 GetGirderLayoutLength(const CGirderKey& girderKey) const;
   Float64 GetGirderSpanLength(const CGirderKey& girderKey) const;
   Float64 GetGirderLength(const CGirderKey& girderKey) const;
   Float64 GetCantileverLength(const CSpanKey& spanKey, pgsTypes::MemberEndType endType) const;

   SegmentIndexType GetSegmentCount(const CGirderKey& girderKey) const;
   Float64 GetSegmentLength(const CSegmentKey& segmentKey) const;
   Float64 GetSegmentSpanLength(const CSegmentKey& segmentKey) const;
   Float64 GetSegmentLayoutLength(const CSegmentKey& segmentKey) const;
   Float64 GetSegmentFramingLength(const CSegmentKey& segmentKey) const;
   Float64 GetSegmentPlanLength(const CSegmentKey& segmentKey) const;
   Float64 GetSegmentStartEndDistance(const CSegmentKey& segmentKey) const;
   Float64 GetSegmentEndEndDistance(const CSegmentKey& segmentKey) const;
   Float64 GetSegmentStartBearingOffset(const CSegmentKey& segmentKey) const;
   Float64 GetSegmentEndBearingOffset(const CSegmentKey& segmentKey) const;
   Float64 GetSegmentStartSupportWidth(const CSegmentKey& segmentKey) const;
   Float64 GetSegmentEndSupportWidth(const CSegmentKey& segmentKey) const;
   bool ModelCantilever(const CSegmentKey& segmentKey,pgsTypes::MemberEndType endType) const;
   Float64 GetSegmentSlope(const CSegmentKey& segmentKey) const;
   Float64 GetSlabOffset(const CSegmentKey& segmentKey, pgsTypes::MemberEndType end) const;
   Float64 GetElevationAdjustment(IntervalIndexType intervalIdx, const pgsPointOfInterest& poi) const;
   Float64 GetRotationAdjustment(IntervalIndexType intervalIdx, const pgsPointOfInterest& poi) const;
   Float64 GetCLPierToCLBearingDistance(const CSegmentKey& segmentKey, pgsTypes::MemberEndType endType, pgsTypes::MeasurementType measure) const;
   Float64 GetCLPierToSegmentEndDistance(const CSegmentKey& segmentKey, pgsTypes::MemberEndType endType, pgsTypes::MeasurementType measure) const;
   //void GetSegmentBearing(const CSegmentKey& segmentKey, IDirection** ppBearing) const;
   //void GetSegmentNormal(const CSegmentKey& segmentKey, IDirection** ppNormal) const;
   //void GetSegmentAngle(const CSegmentKey& segmentKey, pgsTypes::MemberEndType endType, IAngle** ppAngle) const;
   //void GetSegmentSkewAngle(const CSegmentKey& segmentKey, pgsTypes::MemberEndType endType, IAngle** ppAngle) const;
   //CSegmentKey GetSegmentAtPier(PierIndexType pierIdx, const CGirderKey& girderKey) const;
   //void GetSpansForSegment(const CSegmentKey& segmentKey, SpanIndexType* pStartSpanIdx, SpanIndexType* pEndSpanIdx) const;
   //Float64 GetSegmentOffset(const CSegmentKey& segmentKey, Float64 station) const;
   //void GetPoint(const CSegmentKey& segmentKey, Float64 Xpoi, pgsTypes::PlanCoordinateType pcType, IPoint2d** ppPoint) const;
   //void GetPoint(const pgsPointOfInterest& poi, pgsTypes::PlanCoordinateType pcType, IPoint2d** ppPoint) const;
   //void GetStationAndOffset(const CSegmentKey& segmentKey, Float64 Xpoi, Float64* pStation, Float64* pOffset) const;
   //void GetStationAndOffset(const pgsPointOfInterest& poi, Float64* pStation, Float64* pOffset) const;
   //bool GetSegmentPierIntersection(const CSegmentKey& segmentKey, PierIndexType pierIdx, pgsTypes::PlanCoordinateType pcType, IPoint2d** ppPoint) const;
   //bool GetSegmentTempSupportIntersection(const CSegmentKey& segmentKey, SupportIndexType tsIdx, pgsTypes::PlanCoordinateType pcType, IPoint2d** ppPoint) const;
   //bool IsInteriorGirder(const CGirderKey& girderKey) const;
   //bool IsExteriorGirder(const CGirderKey& girderKey) const;
   //bool IsLeftExteriorGirder(const CGirderKey& girderKey) const;
   //bool IsRightExteriorGirder(const CGirderKey& girderKey) const;
   //bool IsObtuseCorner(const CSpanKey& spanKey, pgsTypes::MemberEndType endType) const;
   //bool AreGirderTopFlangesRoughened(const CSegmentKey& segmentKey) const;
   //bool GetSpan(Float64 station, SpanIndexType* pSpanIdx) const;
   //void GetDistanceBetweenGirders(const pgsPointOfInterest& poi, Float64 *pLeft, Float64* pRight) const;
   //void GetBottomFlangeClearance(const pgsPointOfInterest& poi, Float64* pLeft, Float64* pRight) const;
   //std::vector<SpaceBetweenGirder> GetGirderSpacing(Float64 station) const;
   //std::vector<Float64> GetGirderSpacing(SpanIndexType spanIdx, Float64 Xspan) const;
   //std::vector<Float64> GetGirderSpacing(PierIndexType pierIdx, pgsTypes::PierFaceType pierFace, pgsTypes::MeasurementLocation measureLocation, pgsTypes::MeasurementType measureType) const;
   //Float64 GetGirderOffset(GirderIndexType gdrIdx, PierIndexType pierIdx, pgsTypes::PierFaceType pierFace, pgsTypes::OffsetMeasurementType offsetMeasureDatum) const;
   //void GetSpacingAlongGirder(const CGirderKey& girderKey, Float64 Xg, Float64* leftSpacing, Float64* rightSpacing) const;
   //void GetSpacingAlongGirder(const pgsPointOfInterest& poi, Float64* leftSpacing, Float64* rightSpacing) const;
   //GDRCONFIG GetSegmentConfiguration(const CSegmentKey& segmentKey) const;


private:
   CComPtr<IBridge> m_pBridge;
}; 
