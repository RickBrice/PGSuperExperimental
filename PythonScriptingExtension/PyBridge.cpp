#include "stdafx.h"
#include "PyBridge.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

void CPyBridge::Init(IBroker* pBroker)
{
   pBroker->GetInterface(IID_IBridge, (IUnknown**)&m_pBridge);
}

void CPyBridge::Reset()
{
   m_pBridge.Release();
}

bool CPyBridge::HasAsymmetricGirders() const
{
   return m_pBridge->HasAsymmetricGirders();
}

bool CPyBridge::HasTiltedGirders() const
{
   return m_pBridge->HasTiltedGirders();
}

Float64 CPyBridge::GetAlignmentOffset() const
{
   return m_pBridge->GetAlignmentOffset();
}

Float64 CPyBridge::GetLength() const
{
   return m_pBridge->GetLength();
}

Float64 CPyBridge::GetNominalSpanLength(SpanIndexType spanIdx) const
{
   return m_pBridge->GetSpanLength(spanIdx);
}

SpanIndexType CPyBridge::GetSpanCount() const
{
   return m_pBridge->GetSpanCount();
}

PierIndexType CPyBridge::GetPierCount() const
{
   return m_pBridge->GetPierCount();
}

SupportIndexType CPyBridge::GetTemporarySupportCount() const
{
   return m_pBridge->GetTemporarySupportCount();
}

GroupIndexType CPyBridge::GetGirderGroupCount() const
{
   return m_pBridge->GetGirderGroupCount();
}

GirderIndexType CPyBridge::GetGirderlineCount() const
{
   return m_pBridge->GetGirderlineCount();
}

GirderIndexType CPyBridge::GetGirderCount(GroupIndexType grpIdx) const
{
   return m_pBridge->GetGirderCount(grpIdx);
}

GirderIndexType CPyBridge::GetGirderCountBySpan(SpanIndexType spanIdx) const
{
   return m_pBridge->GetGirderCountBySpan(spanIdx);
}

PierIndexType CPyBridge::GetGirderGroupStartPier(GroupIndexType grpIdx) const
{
   return m_pBridge->GetGirderGroupStartPier(grpIdx);
}

PierIndexType CPyBridge::GetGirderGroupEndPier(GroupIndexType grpIdx) const
{
   return m_pBridge->GetGirderGroupEndPier(grpIdx);
}

SpanIndexType CPyBridge::GetGirderGroupStartSpan(GroupIndexType grpIdx) const
{
   return m_pBridge->GetGirderGroupStartSpan(grpIdx);
}

SpanIndexType CPyBridge::GetGirderGroupEndSpan(GroupIndexType grpIdx) const
{
   return m_pBridge->GetGirderGroupEndSpan(grpIdx);
}

GroupIndexType CPyBridge::GetGirderGroupIndex(SpanIndexType spanIdx) const
{
   return m_pBridge->GetGirderGroupIndex(spanIdx);
}

Float64 CPyBridge::GetSpanLength(const CSpanKey& spanKey) const
{
   return m_pBridge->GetSpanLength(spanKey);
}

Float64 CPyBridge::GetFullSpanLength(const CSpanKey& spanKey) const
{
   return m_pBridge->GetFullSpanLength(spanKey);
}

Float64 CPyBridge::GetGirderlineLength(GirderIndexType gdrLineIdx) const
{
   return m_pBridge->GetGirderlineLength(gdrLineIdx);
}

Float64 CPyBridge::GetGirderLayoutLength(const CGirderKey& girderKey) const
{
   return m_pBridge->GetGirderLayoutLength(girderKey);
}

Float64 CPyBridge::GetGirderSpanLength(const CGirderKey& girderKey) const
{
   return m_pBridge->GetGirderSpanLength(girderKey);
}

Float64 CPyBridge::GetGirderLength(const CGirderKey& girderKey) const
{
   return m_pBridge->GetGirderLength(girderKey);
}

Float64 CPyBridge::GetCantileverLength(const CSpanKey& spanKey, pgsTypes::MemberEndType endType) const
{
   return m_pBridge->GetCantileverLength(spanKey, endType);
}

SegmentIndexType CPyBridge::GetSegmentCount(const CGirderKey& girderKey) const
{
	return m_pBridge->GetSegmentCount(girderKey);
}

Float64 CPyBridge::GetSegmentLength(const CSegmentKey& segmentKey) const
{
   return m_pBridge->GetSegmentLength(segmentKey);
}

Float64 CPyBridge::GetSegmentSpanLength(const CSegmentKey& segmentKey) const
{
   return m_pBridge->GetSegmentSpanLength(segmentKey);
}

Float64 CPyBridge::GetSegmentLayoutLength(const CSegmentKey& segmentKey) const
{
   return m_pBridge->GetSegmentLayoutLength(segmentKey);
}

Float64 CPyBridge::GetSegmentFramingLength(const CSegmentKey& segmentKey) const
{
   return m_pBridge->GetSegmentFramingLength(segmentKey);
}

Float64 CPyBridge::GetSegmentPlanLength(const CSegmentKey& segmentKey) const
{
   return m_pBridge->GetSegmentPlanLength(segmentKey);
}

Float64 CPyBridge::GetSegmentStartEndDistance(const CSegmentKey& segmentKey) const
{
   return m_pBridge->GetSegmentStartEndDistance(segmentKey);
}

Float64 CPyBridge::GetSegmentEndEndDistance(const CSegmentKey& segmentKey) const
{
   return m_pBridge->GetSegmentEndEndDistance(segmentKey);
}

Float64 CPyBridge::GetSegmentStartBearingOffset(const CSegmentKey& segmentKey) const
{
   return m_pBridge->GetSegmentStartBearingOffset(segmentKey);
}

Float64 CPyBridge::GetSegmentEndBearingOffset(const CSegmentKey& segmentKey) const
{
   return m_pBridge->GetSegmentEndBearingOffset(segmentKey);
}

Float64 CPyBridge::GetSegmentStartSupportWidth(const CSegmentKey& segmentKey) const
{
   return m_pBridge->GetSegmentStartSupportWidth(segmentKey);
}

Float64 CPyBridge::GetSegmentEndSupportWidth(const CSegmentKey& segmentKey) const
{
   return m_pBridge->GetSegmentEndSupportWidth(segmentKey);
}

bool CPyBridge::ModelCantilever(const CSegmentKey& segmentKey, pgsTypes::MemberEndType endType) const
{
   bool bLeft, bRight;
   m_pBridge->ModelCantilevers(segmentKey, &bLeft, &bRight);
   return (endType == pgsTypes::metStart ? bLeft : bRight);
}

Float64 CPyBridge::GetSegmentSlope(const CSegmentKey& segmentKey) const
{
   return m_pBridge->GetSegmentSlope(segmentKey);
}

Float64 CPyBridge::GetSlabOffset(const CSegmentKey& segmentKey, pgsTypes::MemberEndType end) const
{
   return m_pBridge->GetSlabOffset(segmentKey,end);
}

Float64 CPyBridge::GetElevationAdjustment(IntervalIndexType intervalIdx, const pgsPointOfInterest& poi) const
{
   return m_pBridge->GetElevationAdjustment(intervalIdx, poi);
}

Float64 CPyBridge::GetRotationAdjustment(IntervalIndexType intervalIdx, const pgsPointOfInterest& poi) const
{
   return m_pBridge->GetRotationAdjustment(intervalIdx, poi);
}

Float64 CPyBridge::GetCLPierToCLBearingDistance(const CSegmentKey& segmentKey, pgsTypes::MemberEndType endType, pgsTypes::MeasurementType measure) const
{
   return m_pBridge->GetCLPierToCLBearingDistance(segmentKey, endType, measure);
}

Float64 CPyBridge::GetCLPierToSegmentEndDistance(const CSegmentKey& segmentKey, pgsTypes::MemberEndType endType, pgsTypes::MeasurementType measure) const
{
   return m_pBridge->GetCLPierToSegmentEndDistance(segmentKey, endType, measure);
}

