#include "stdafx.h"
#include "PyBridgeTypes.h"

#include "PyBridge.h"

using namespace boost::python;

void export_bridge_types()
{
   class_<CPyBridge>("Bridge")
      .def("HasAsymmetricGirders",&CPyBridge::HasAsymmetricGirders)
      .def("HasTiltedGirders",&CPyBridge::HasTiltedGirders)
      .def("GetLength",&CPyBridge::GetLength)
      .def("GetNominalSpanLength",&CPyBridge::GetNominalSpanLength)
      .def("GetAlignmentOffset", &CPyBridge::GetAlignmentOffset)
      .def("GetSpanCount",&CPyBridge::GetSpanCount)
      .def("GetPierCount",&CPyBridge::GetPierCount)
      .def("GetTemporarySupportCount",&CPyBridge::GetTemporarySupportCount)
      .def("GetGirderGroupCount",&CPyBridge::GetGirderGroupCount)
      .def("GetGirderCount", &CPyBridge::GetGirderCount)
      .def("GetGirderlineCount",&CPyBridge::GetGirderlineCount)
      .def("GetGirderCountBySpan",&CPyBridge::GetGirderCountBySpan)
      .def("GetGirderGroupStartPier",&CPyBridge::GetGirderGroupStartPier)
      .def("GetGirderGroupEndPier",&CPyBridge::GetGirderGroupEndPier)
      .def("GetGirderGroupStartSpan",&CPyBridge::GetGirderGroupStartSpan)
      .def("GetGirderGroupEndSpan",&CPyBridge::GetGirderGroupEndSpan)
      .def("GetGirderGroupIndex",&CPyBridge::GetGirderGroupIndex)
      .def("GetSpanLength",&CPyBridge::GetSpanLength)
      .def("GetFullSpanLength",&CPyBridge::GetFullSpanLength)
      .def("GetGirderlineLength",&CPyBridge::GetGirderlineLength)
      .def("GetGirderLayoutLength",&CPyBridge::GetGirderLayoutLength)
      .def("GetGirderSpanLength",&CPyBridge::GetGirderSpanLength)
      .def("GetGirderLength",&CPyBridge::GetGirderLength)
      .def("GetCantileverLength",&CPyBridge::GetCantileverLength)

      .def("GetSegmentCount", &CPyBridge::GetSegmentCount)
      .def("GetSegmentLength",&CPyBridge::GetSegmentLength)
      .def("GetSegmentSpanLength",&CPyBridge::GetSegmentSpanLength)
      .def("GetSegmentLayoutLength", &CPyBridge::GetSegmentLayoutLength)
      .def("GetSegmentFramingLength",&CPyBridge::GetSegmentFramingLength)
      .def("GetSegmentPlanLength",&CPyBridge::GetSegmentPlanLength)
      .def("GetSegmentStartEndDistance",&CPyBridge::GetSegmentStartEndDistance)
      .def("GetSegmentEndEndDistance",&CPyBridge::GetSegmentEndEndDistance)
      .def("GetSegmentStartBearingOffset",&CPyBridge::GetSegmentStartBearingOffset)
      .def("GetSegmentEndBearingOffset",&CPyBridge::GetSegmentEndBearingOffset)
      .def("GetSegmentStartSupportWidth",&CPyBridge::GetSegmentStartSupportWidth)
      .def("GetSegmentEndSupportWidth",&CPyBridge::GetSegmentEndSupportWidth)
      .def("ModelCantilever",&CPyBridge::ModelCantilever)
      .def("GetSegmentSlope",&CPyBridge::GetSegmentSlope)
 
      .def("GetElevationAdjustment",&CPyBridge::GetElevationAdjustment)
      .def("GetRotationAdjustment",&CPyBridge::GetRotationAdjustment)
      .def("GetCLPierToCLBearingDistance",&CPyBridge::GetCLPierToCLBearingDistance)
      .def("GetCLPierToSegmentEndDistance",&CPyBridge::GetCLPierToSegmentEndDistance)
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
   ;
}