#include "stdafx.h"

#include <IFace\Project.h>
#include <IFace\VersionInfo.h>
#include <IFace\Alignment.h>


#include "IfcAlignmentConverter.h"
#include "IfcAlignmentConverterException.h"


#if defined _DEBUG
#pragma comment(lib,"F:/IfcOpenShell/_build-vs2015-x64/Debug/IfcParse.lib")
#else
#pragma comment(lib,"F:/IfcOpenShell/_build-vs2015-x64/Release/IfcParse.lib")
#endif

#define ALIGNMENT_NAME _T("Alignment1")

// Constants for tracking the state of converting the profile data
#define PROFILE_NOT_STARTED -1
#define PROFILE_ESTABLISHED  1

// Use this throw macro when the data conversion cannot continue
// The catcher, or other, is responsible for deleting it
#define IFC_THROW(_s_) throw new CIfcAlignmentConverterException(_s_);

Float64 CIfcAlignmentConverter::m_Precision = 0.001;

HRESULT SameLocation(IPoint2d* pnt1, IPoint2d* pnt2,Float64 tolerance)
{
   Float64 x1, y1, x2, y2;
   pnt1->get_X(&x1);
   pnt1->get_Y(&y1);
   pnt2->get_X(&x2);
   pnt2->get_Y(&y2);

   return IsEqual(x1, x2, tolerance) && IsEqual(y1, y2, tolerance) ? S_OK : S_FALSE;
}

template <typename Schema>
typename Schema::IfcAlignment2DHorizontalSegment* CreateLineSegment(typename Schema::IfcCartesianPoint* pStartPoint, Float64 angle, Float64 length)
{
   return new Schema::IfcAlignment2DHorizontalSegment(boost::none, boost::none, boost::none, new Schema::IfcLineSegment2D(pStartPoint, angle, length));
}

template <typename Schema>
typename Schema::IfcAlignmentHorizontalSegment* CreateLineSegment(typename Schema::IfcCartesianPoint* pStartPoint, Float64 angle, Float64 length)
{
   return new Schema::IfcAlignmentHorizontalSegment(boost::none, boost::none, pStartPoint, angle, 0.0, 0.0, length, boost::none, Schema::IfcAlignmentHorizontalSegmentTypeEnum::IfcAlignmentHorizontalSegmentType_LINE);
}

template <typename Schema>
typename Schema::IfcAlignment2DHorizontalSegment* CreateTransitionCurve(typename Schema::IfcCartesianPoint* pStartPoint, Float64 startDirection, Float64 Ls, Float64 R1, Float64 R2, bool bIsCCW1, bool bIsCCW2)
{
   return new Schema::IfcAlignment2DHorizontalSegment(boost::none, boost::none, boost::none, new Schema::IfcTransitionCurveSegment2D(pStartPoint, startDirection, Ls, IsZero(R1) ? boost::none : boost::optional<Float64>(R1), IsZero(R2) ? boost::none : boost::optional<Float64>(R2), bIsCCW1, bIsCCW2, Schema::IfcTransitionCurveType::IfcTransitionCurveType_CLOTHOIDCURVE));
}

template <typename Schema>
typename Schema::IfcAlignmentHorizontalSegment* CreateTransitionCurve(typename Schema::IfcCartesianPoint* pStartPoint, Float64 startDirection, Float64 Ls, Float64 R1, Float64 R2, bool bIsCCW1, bool bIsCCW2)
{
   return new Schema::IfcAlignmentHorizontalSegment(boost::none, boost::none, pStartPoint, startDirection, (bIsCCW1 ? 1.0 : -1.0)*R1, (bIsCCW2 ? 1.0 : -1.0)*R2, Ls, boost::none, Schema::IfcAlignmentHorizontalSegmentTypeEnum::IfcAlignmentHorizontalSegmentType_CLOTHOID);
}

template <typename Schema>
typename Schema::IfcAlignment2DHorizontalSegment* CreateCircularCurve(typename Schema::IfcCartesianPoint* pStartPoint, Float64 startDirection, Float64 Lc, Float64 R, bool bIsCCW)
{
   return new Schema::IfcAlignment2DHorizontalSegment(boost::none, boost::none, boost::none, new Schema::IfcCircularArcSegment2D(pStartPoint, startDirection, Lc, R, bIsCCW));
}

template <typename Schema>
typename Schema::IfcAlignmentHorizontalSegment* CreateCircularCurve(typename Schema::IfcCartesianPoint* pStartPoint, Float64 startDirection, Float64 Lc, Float64 R, bool bIsCCW)
{
   return new Schema::IfcAlignmentHorizontalSegment(boost::none, boost::none, pStartPoint, startDirection, (bIsCCW ? 1.0 : -1.0)*R, (bIsCCW ? 1.0 : -1.0)*R, Lc, boost::none, Schema::IfcAlignmentHorizontalSegmentTypeEnum::IfcAlignmentHorizontalSegmentType_CIRCULARARC);
}

template <typename Schema>
typename Schema::IfcAlignment2DHorizontal* CreateAlignment(std::string guid, typename Schema::IfcOwnerHistory* ownerHistory, typename Schema::IfcLocalPlacement* placement, typename Schema::IfcProductRepresentation* representation, Float64 startStation, boost::shared_ptr<IfcTemplatedEntityList<typename Schema::IfcAlignment2DHorizontalSegment>>& alignment_segments)
{
   return new Schema::IfcAlignment2DHorizontal(startStation, alignment_segments);
}

template <typename Schema>
typename Schema::IfcAlignmentHorizontal* CreateAlignment(std::string  guid, typename Schema::IfcOwnerHistory* ownerHistory, typename Schema::IfcLocalPlacement* placement, typename Schema::IfcProductRepresentation* representation, Float64 startStation, boost::shared_ptr<IfcTemplatedEntityList<typename Schema::IfcAlignmentHorizontalSegment>>& alignment_segments)
{
   return new Schema::IfcAlignmentHorizontal(guid, ownerHistory, boost::none, boost::none, boost::none, placement, representation, startStation, alignment_segments);
}

template <typename Schema>
typename boost::shared_ptr<IfcTemplatedEntityList<typename Schema::IfcAlignment2DHorizontalSegment>> CreateHorizontalSegments()
{
   return boost::shared_ptr<IfcTemplatedEntityList<typename Schema::IfcAlignment2DHorizontalSegment>>(new IfcTemplatedEntityList<typename Schema::IfcAlignment2DHorizontalSegment>());
}

template <typename Schema>
typename boost::shared_ptr<IfcTemplatedEntityList<typename Schema::IfcAlignmentHorizontalSegment>> CreateHorizontalSegments()
{
   return boost::shared_ptr<IfcTemplatedEntityList<typename Schema::IfcAlignmentHorizontalSegment>>(new IfcTemplatedEntityList<typename Schema::IfcAlignmentHorizontalSegment>());
}

template <typename Schema,typename AlignmentType>
typename AlignmentType* BuildAlignment(IBroker* pBroker, IfcHierarchyHelper<Schema>& file)
{
   auto alignment_segments = CreateHorizontalSegments<Schema>();

   GET_IFACE2(pBroker, IGeometry, pGeometry);
   GET_IFACE2(pBroker, IRoadway, pAlignment);
   Float64 startStation, startElevation, startGrade;
   CComPtr<IPoint2d> startPoint;
   pAlignment->GetStartPoint(2, &startStation, &startElevation, &startGrade, &startPoint);

   // create the start point
   std::vector<Float64> vCoordinates;
   vCoordinates.resize(2);
   startPoint->Location(&vCoordinates[0], &vCoordinates[1]);
   auto ifc_start_point = new Schema::IfcCartesianPoint(vCoordinates);

   // loop over all the horizontal curves
   CComPtr<IPoint2d> prevPoint = startPoint;
   auto ifc_prev_point = ifc_start_point;
   IndexType nHCurves = pAlignment->GetCurveCount();
   for (IndexType i = 0; i < nHCurves; i++)
   {
      CComPtr<IHorzCurve> curve;
      pAlignment->GetCurve(i, &curve);

      CComPtr<IPoint2d> pntTS;
      pAlignment->GetCurvePoint(i, cptTS, pgsTypes::pcGlobal, &pntTS);

      // create a line segment from end of previous alignment segment to the start of this curve
      if (SameLocation(pntTS,prevPoint, CIfcAlignmentConverter::GetPrecision()) == S_FALSE)
      {
         Float64 dist;
         CComPtr<IDirection> direction;
         pGeometry->Inverse(prevPoint, pntTS, &dist, &direction);
         Float64 angle;
         direction->get_Value(&angle);

         auto ifc_line_segment = CreateLineSegment<Schema>(ifc_prev_point, angle, dist);
         alignment_segments->push(ifc_line_segment);
      }

      // create this horizontal curve

      std::array<Float64, 2> Lspiral;
      curve->get_SpiralLength(spEntry, &Lspiral[spEntry]);
      curve->get_SpiralLength(spExit, &Lspiral[spExit]);

      Float64 Lc;
      curve->get_CurveLength(&Lc);

      Float64 R;
      curve->get_Radius(&R);

      ATLASSERT(0 < R); // need to deal with zero radius curves, which gives us an angle point in the alignment

      CurveDirectionType curve_direction;
      curve->get_Direction(&curve_direction);
      bool bIsCCW = (curve_direction == cdLeft) ? true : false;

      // curve starts at the Tangent to Spiral point (TS)
      pntTS->Location(&vCoordinates[0], &vCoordinates[1]);
      auto ifc_start = new Schema::IfcCartesianPoint(vCoordinates);

      if (0.0 < Lspiral[spEntry])
      {
         // there is an entry spiral
         CComPtr<IDirection> bkTangentBrg;
         curve->get_BkTangentBrg(&bkTangentBrg);
         Float64 bk_tangent_direction;
         bkTangentBrg->get_Value(&bk_tangent_direction);

         auto ifc_entry_spiral = CreateTransitionCurve<Schema>(ifc_start, bk_tangent_direction, Lspiral[spEntry], 0, R, bIsCCW, bIsCCW);
         alignment_segments->push(ifc_entry_spiral);

         // spiral ends at the Spiral to Curve point (CS)
         CComPtr<IPoint2d> pntSC;
         pAlignment->GetCurvePoint(i, cptSC, pgsTypes::pcGlobal, &pntSC);
         pntSC->Location(&vCoordinates[0], &vCoordinates[1]);
         ifc_start = new Schema::IfcCartesianPoint(vCoordinates);
      }

      // build the horizontal curve
      CComPtr<IDirection> bkTangentBrgCurve;
      curve->get_CurveBkTangentBrg(&bkTangentBrgCurve);
      Float64 bk_tangent_direction_curve;
      bkTangentBrgCurve->get_Value(&bk_tangent_direction_curve);
      auto ifc_hcurve = CreateCircularCurve<Schema>(ifc_start, bk_tangent_direction_curve, Lc, R, bIsCCW);
      alignment_segments->push(ifc_hcurve);

      if (0.0 < Lspiral[spExit])
      {
         // there is an exit spiral

         CComPtr<IDirection> fwdTangentBrgCurve;
         curve->get_CurveFwdTangentBrg(&fwdTangentBrgCurve); // forward tangent of curve is start tangent to exit spiral
         Float64 fwd_tangent_direction_curve;
         fwdTangentBrgCurve->get_Value(&fwd_tangent_direction_curve);

         // spiral starts at the Curve to Spiral point (CS)
         CComPtr<IPoint2d> pntCS;
         pAlignment->GetCurvePoint(i, cptCS, pgsTypes::pcGlobal, &pntCS);
         pntCS->Location(&vCoordinates[0], &vCoordinates[1]);
         ifc_start = new Schema::IfcCartesianPoint(vCoordinates);

         auto ifc_exit_spiral = CreateTransitionCurve<Schema>(ifc_start, fwd_tangent_direction_curve, Lspiral[spExit], R, 0, bIsCCW, bIsCCW);
         alignment_segments->push(ifc_exit_spiral);
      }

      // end of this curve (Spiral to Tangent, ST) becomes previous point for next alignment segment
      prevPoint.Release();
      pAlignment->GetCurvePoint(i, cptST, pgsTypes::pcGlobal, &prevPoint);
      prevPoint->Location(&vCoordinates[0], &vCoordinates[1]);
      ifc_prev_point = new Schema::IfcCartesianPoint(vCoordinates);
   }

   // build a linear segment from end of previous alignment segment to the end of the alignment
   Float64 endStation, endElevation, endGrade;
   CComPtr<IPoint2d> endPoint;
   pAlignment->GetEndPoint(2, &endStation, &endElevation, &endGrade, &endPoint);

   if (SameLocation(prevPoint,endPoint,CIfcAlignmentConverter::GetPrecision()) == S_FALSE)
   {
      // end the alignment with a line segment
      Float64 dist;
      CComPtr<IDirection> direction;
      pGeometry->Inverse(prevPoint, endPoint, &dist, &direction);
      Float64 angle;
      direction->get_Value(&angle);

      auto ifc_line_segment = CreateLineSegment<Schema>(ifc_start_point, angle, dist);
      alignment_segments->push(ifc_line_segment);
   }

   // create a horizontal alignment from all the alignment segments
   auto horizontal_alignment = CreateAlignment<Schema>(IfcParse::IfcGlobalId(), file.getSingle<Schema::IfcOwnerHistory>(), file.getSingle<Schema::IfcLocalPlacement>(), file.getSingle<Schema::IfcProductRepresentation>(), startStation, alignment_segments);

   return horizontal_alignment;
}

template <typename Schema>
typename Schema::IfcAlignment2DHorizontal* BuildAlignment(IBroker* pBroker, IfcHierarchyHelper<Schema>& file)
{
   return BuildAlignment<Schema, Schema::IfcAlignment2DHorizontal>(pBroker, file);
}

template <typename Schema>
typename Schema::IfcAlignmentHorizontal* BuildAlignment(IBroker* pBroker, IfcHierarchyHelper<Schema>& file)
{
   return BuildAlignment<Schema, Schema::IfcAlignmentHorizontal>(pBroker, file);
}

template <typename Schema>
void CreateLinearSegment(Float64 startDistAlong, Float64 L, Float64 startHeight, Float64 startGradient, typename Schema::IfcAlignment2DVerSegLine** ppSegment)
{
   *ppSegment = new Schema::IfcAlignment2DVerSegLine(boost::none, boost::none, boost::none, startDistAlong, L, startHeight, startGradient);
}

template <typename Schema>
void CreateLinearSegment(Float64 startDistAlong, Float64 L, Float64 startHeight, Float64 startGradient, typename Schema::IfcAlignmentVerticalSegment** ppSegment)
{
   *ppSegment = new Schema::IfcAlignmentVerticalSegment(boost::none, boost::none, startDistAlong, L, startHeight, startGradient, startGradient, boost::none, Schema::IfcAlignmentVerticalSegmentTypeEnum::IfcAlignmentVerticalSegmentType_CONSTANTGRADIENT);
}

template <typename Schema>
void CreateParabolicSegment(Float64 startDistAlong, Float64 L, Float64 startHeight, Float64 startGradient, Float64 endGradient, Float64 k, typename Schema::IfcAlignment2DVerSegParabolicArc** ppSegment)
{
   *ppSegment = new Schema::IfcAlignment2DVerSegParabolicArc(boost::none, boost::none, boost::none, startDistAlong, L, startHeight, startGradient, fabs(1 / k), ::BinarySign(k) < 0 ? true : false);
}

template <typename Schema>
void CreateParabolicSegment(Float64 startDistAlong, Float64 L, Float64 startHeight, Float64 startGradient, Float64 endGradient, Float64 k, typename Schema::IfcAlignmentVerticalSegment** ppSegment)
{
   *ppSegment = new Schema::IfcAlignmentVerticalSegment(boost::none, boost::none, startDistAlong, L, startHeight, startGradient, endGradient, k, Schema::IfcAlignmentVerticalSegmentTypeEnum::IfcAlignmentVerticalSegmentType_PARABOLICARC);
}

template <typename Schema>
typename Schema::IfcAlignment2DVertical* CreateProfile(std::string guid, typename Schema::IfcOwnerHistory* ownerHistory, typename Schema::IfcLocalPlacement* placement, typename Schema::IfcProductRepresentation* representation, boost::shared_ptr<IfcTemplatedEntityList<typename Schema::IfcAlignment2DVerticalSegment>>& profile_segments)
{
   return new Schema::IfcAlignment2DVertical(profile_segments);
}


template <typename Schema>
typename Schema::IfcAlignmentVertical* CreateProfile(std::string guid, typename Schema::IfcOwnerHistory* ownerHistory, typename Schema::IfcLocalPlacement* placement, typename Schema::IfcProductRepresentation* representation, boost::shared_ptr<IfcTemplatedEntityList<typename Schema::IfcAlignmentVerticalSegment>>& profile_segments)
{
   return new Schema::IfcAlignmentVertical(guid, ownerHistory, boost::none, boost::none, boost::none, placement, representation, profile_segments);
}

template <typename Schema>
typename boost::shared_ptr<IfcTemplatedEntityList<typename Schema::IfcAlignment2DVerticalSegment>> CreateProfileSegments()
{
   return boost::shared_ptr<IfcTemplatedEntityList<typename Schema::IfcAlignment2DVerticalSegment>>(new IfcTemplatedEntityList<typename Schema::IfcAlignment2DVerticalSegment>());
}

template <typename Schema>
typename boost::shared_ptr<IfcTemplatedEntityList<typename Schema::IfcAlignmentVerticalSegment>> CreateProfileSegments()
{
   return boost::shared_ptr<IfcTemplatedEntityList<typename Schema::IfcAlignmentVerticalSegment>>(new IfcTemplatedEntityList<typename Schema::IfcAlignmentVerticalSegment>());
}

template <typename Schema,typename ProfileType,typename LinearSegment,typename ParabolicSegment>
typename ProfileType* BuildProfile(IBroker* pBroker, IfcHierarchyHelper<Schema>& file)
{
   auto profile_segments = CreateProfileSegments<Schema>();

   // Profile is defined by profile segments located at "distance from start" of the alignment and "length".
   // We can't use stations to define the profile.
   // Distance from start is taken to be Station - Start Station

   GET_IFACE2(pBroker, IRoadway, pAlignment);

   Float64 startStation, startElevation, startGrade;
   CComPtr<IPoint2d> startPoint;
   pAlignment->GetStartPoint(2, &startStation, &startElevation, &startGrade, &startPoint);

   Float64 prev_end_dist_along = 0; // startStation; // this is distance along alignment, not station
   Float64 prev_end_gradiant = startGrade;
   Float64 prev_end_height = startElevation;

   IndexType nVCurves = pAlignment->GetVertCurveCount();
   for (IndexType i = 0; i < nVCurves; i++)
   {
      CComPtr<IVertCurve> curve;
      pAlignment->GetVertCurve(i, &curve);

      CComPtr<IProfilePoint> startPoint;
      curve->get_BVC(&startPoint);
      CComPtr<IStation> station;
      startPoint->get_Station(&station);
      Float64 start_height;
      startPoint->get_Elevation(&start_height);
      ZoneIndexType zoneIdx;
      Float64 start_dist_along;
      station->GetStation(&zoneIdx, &start_dist_along);
      start_dist_along -= startStation;
#pragma Reminder("WORKING HERE - How to deal with station equations?") // see IfcReferent

      if (!IsEqual(prev_end_dist_along, start_dist_along))
      {
         // create a linear segment between the last profile element and this curve
         Float64 length = start_dist_along - prev_end_dist_along;
         LinearSegment* linear_segment;
         CreateLinearSegment<Schema>(prev_end_dist_along, length, prev_end_height, prev_end_gradiant, &linear_segment);
         profile_segments->push(linear_segment);
      }

      Float64 l1, l2;
      curve->get_L1(&l1);
      curve->get_L2(&l2);
      if (!IsEqual(l1, l2) && !IsZero(l2))
      {
         // compound curve
         CComPtr<IProfilePoint> pviPoint;
         curve->get_PVI(&pviPoint);
         CComPtr<IStation> pviStation;
         pviPoint->get_Station(&pviStation);
         Float64 pviElevation;
         curve->Elevation(CComVariant(pviStation), &pviElevation);
         Float64 pviGrade;
         curve->Grade(CComVariant(pviStation), &pviGrade);

         Float64 start_gradient, end_gradient;
         curve->get_EntryGrade(&start_gradient);
         curve->get_ExitGrade(&end_gradient);

         Float64 k1, k2;
         curve->get_K1(&k1);
         curve->get_K2(&k2);

         ParabolicSegment* parabolic_segment1;
         CreateParabolicSegment<Schema>(start_dist_along, l1, start_height, start_gradient, pviGrade, k1, &parabolic_segment1);
         profile_segments->push(parabolic_segment1);

         ParabolicSegment* parabolic_segment2;
         CreateParabolicSegment<Schema>(start_dist_along + l1, l2, pviElevation, pviGrade, end_gradient, k2, &parabolic_segment2);
         profile_segments->push(parabolic_segment2);
      }
      else
      {
         Float64 horizontal_length;
         curve->get_Length(&horizontal_length);
         Float64 start_gradient, end_gradient;
         curve->get_EntryGrade(&start_gradient);
         curve->get_ExitGrade(&end_gradient);

         if (IsEqual(start_gradient, end_gradient))
         {
            // this is just a straight line
            LinearSegment* linear_segment;
            CreateLinearSegment<Schema>(prev_end_dist_along, l1, prev_end_height, start_gradient, &linear_segment);
            profile_segments->push(linear_segment);
         }
         else
         {
            Float64 k;
            curve->get_K1(&k);
            ParabolicSegment* parabolic_segment;
            CreateParabolicSegment<Schema>(start_dist_along, horizontal_length, start_height, start_gradient, end_gradient, k, &parabolic_segment);
            profile_segments->push(parabolic_segment);
         }
      }

      // setup parameters for next loop
      CComPtr<IProfilePoint> evc;
      curve->get_EVC(&evc);
      CComPtr<IStation> evcStation;
      evc->get_Station(&evcStation);
      evcStation->GetStation(&zoneIdx, &prev_end_dist_along);
      prev_end_dist_along -= startStation;
      evc->get_Elevation(&prev_end_height);
      curve->get_ExitGrade(&prev_end_gradiant);
#pragma Reminder("WORKING HERE - How to deal with station equations?") // see IfcReferent
   }

   Float64 endStation, endElevation, endGrade;
   CComPtr<IPoint2d> endPoint;
   pAlignment->GetEndPoint(2, &endStation, &endElevation, &endGrade, &endPoint);
   if (!IsEqual(prev_end_dist_along, endStation))
   {
      // create a linear segment between the last profile element and the end of the alignment
      ATLASSERT(IsEqual(prev_end_gradiant, endGrade));
      Float64 length = endStation - startStation - prev_end_dist_along;
      LinearSegment* linear_segment;
      CreateLinearSegment<Schema>(prev_end_dist_along, length, prev_end_height, prev_end_gradiant, &linear_segment);
      profile_segments->push(linear_segment);

      // check elevation
      ATLASSERT(IsEqual(endElevation, prev_end_height + length*prev_end_gradiant));
   }

   auto vertical_profile = CreateProfile<Schema>(IfcParse::IfcGlobalId(), file.getSingle<Schema::IfcOwnerHistory>(), file.getSingle<Schema::IfcLocalPlacement>(), file.getSingle<Schema::IfcProductRepresentation>(), profile_segments);

   return vertical_profile;
}

template <typename Schema>
typename Schema::IfcAlignment2DVertical* BuildProfile(IBroker* pBroker, IfcHierarchyHelper<Schema>& file)
{
   return BuildProfile<Schema, Schema::IfcAlignment2DVertical, Schema::IfcAlignment2DVerSegLine,Schema::IfcAlignment2DVerSegParabolicArc>(pBroker, file);
}

template <typename Schema>
typename Schema::IfcAlignmentVertical* BuildProfile(IBroker* pBroker, IfcHierarchyHelper<Schema>& file)
{
   return BuildProfile<Schema, Schema::IfcAlignmentVertical, Schema::IfcAlignmentVerticalSegment,Schema::IfcAlignmentVerticalSegment>(pBroker, file);
}

CIfcAlignmentConverter::CIfcAlignmentConverter(void)
{
   m_pLengthUnit = nullptr;

   m_bAlignmentStarted = false;
   m_ProfileState = PROFILE_NOT_STARTED;

   m_CogoEngine.CoCreateInstance(CLSID_CogoEngine);
   m_GeomUtil.CoCreateInstance(CLSID_GeomUtil);

   m_LastAlignmentType = Unknown;
}

CIfcAlignmentConverter::~CIfcAlignmentConverter(void)
{
}

template <typename Schema>
void CIfcAlignmentConverter::InitUnits(IfcParse::IfcFile& file)
{
   auto geometric_representation_contexts = file.instances_by_type<Schema::IfcGeometricRepresentationContext>();
   ATLASSERT(geometric_representation_contexts->size() == 1);// From my understanding there should only be the one
   auto geometric_representation_context = *(geometric_representation_contexts->begin());
   if (geometric_representation_context->hasPrecision())
   {
      m_Precision = geometric_representation_context->Precision();
   }


#pragma Reminder("WORKING HERE - UNITS - THERE ARE MANY CASES THIS DOESN'T DEAL WITH")
   auto unit_assignment_instances = file.instances_by_type<Schema::IfcUnitAssignment>();
   ATLASSERT(unit_assignment_instances->size() == 1);
   auto unit_assignment = *(unit_assignment_instances->begin());
   auto units = unit_assignment->Units();
   for (auto unit : *units)
   {
      auto derived_unit = unit->as<Schema::IfcDerivedUnit>();
      auto monitary_unit = unit->as<Schema::IfcMonetaryUnit>();
      auto si_unit = unit->as<Schema::IfcSIUnit>();
      auto conversion_based_unit = unit->as<Schema::IfcConversionBasedUnit>();
      auto conversion_based_unit_with_offset = unit->as<Schema::IfcConversionBasedUnitWithOffset>();

      if (si_unit && si_unit->Name() == Schema::IfcSIUnitName::IfcSIUnitName_METRE)
      {
         if (si_unit->hasPrefix())
         {
            switch (si_unit->Prefix())
            {
            case Schema::IfcSIPrefix::IfcSIPrefix_KILO:
               m_pLengthUnit = &unitMeasure::Kilometer;
               break;

            case Schema::IfcSIPrefix::IfcSIPrefix_CENTI:
               m_pLengthUnit = &unitMeasure::Centimeter;
               break;

            case Schema::IfcSIPrefix::IfcSIPrefix_MILLI:
               m_pLengthUnit = &unitMeasure::Millimeter;
               break;

            default:
               ATLASSERT(false); // unit prefix isn't supported
            }
         }
         else
         {
            m_pLengthUnit = &unitMeasure::Meter;
         }
         break;
      }

      if (conversion_based_unit && conversion_based_unit->UnitType() == Schema::IfcUnitEnum::IfcUnit_LENGTHUNIT)
      {
         if (conversion_based_unit->UnitType() == Schema::IfcUnitEnum::IfcUnit_LENGTHUNIT)
         {
            auto measure_with_unit = conversion_based_unit->ConversionFactor();
            auto unit_component = measure_with_unit->UnitComponent()->as<Schema::IfcSIUnit>();
            ATLASSERT(unit_component->Name() == Schema::IfcSIUnitName::IfcSIUnitName_METRE);
            ATLASSERT(unit_component->hasPrefix() == false); // not dealing with conversion factors to anything but meter

            Float64 conversion_factor;
            try
            {
               auto value_component = measure_with_unit->ValueComponent()->as<Schema::IfcLengthMeasure>();
               ATLASSERT(value_component); // not dealing with anything but simple conversion factors
               conversion_factor = 1 / (*value_component);
            }
            catch (IfcParse::IfcInvalidTokenException& e)
            {
               // Was expecting something like 
               // #15 = IFCMEASUREWITHUNIT(IFCLENGTHMEASURE(3.28083333333333), #16);
               // where the expected token is IFCLENGTHMEASURE, but instead found something like
               // #15=IFCMEASUREWITHUNIT(3.28083333333333,#16);
               // we'll just get the value and keep going
               TRACE(e.what());
               Argument* pArgument = measure_with_unit->get("ValueComponent");
               ATLASSERT(pArgument->type() == IfcUtil::Argument_DOUBLE);
               double value = double(*pArgument);
               conversion_factor = 1 / value;
            }

            if (IsEqual(conversion_factor, unitMeasure::Feet.GetConvFactor()))
            {
               m_pLengthUnit = &unitMeasure::Feet;
            }
            else if (IsEqual(conversion_factor, unitMeasure::USSurveyFoot.GetConvFactor()))
            {
               m_pLengthUnit = &unitMeasure::USSurveyFoot;
            }
            else if (IsEqual(conversion_factor, unitMeasure::Inch.GetConvFactor()))
            {
               m_pLengthUnit = &unitMeasure::Inch;
            }
            else if (IsEqual(conversion_factor, unitMeasure::Mile.GetConvFactor()))
            {
               m_pLengthUnit = &unitMeasure::Mile;
            }
            else if (IsEqual(conversion_factor, unitMeasure::Yard.GetConvFactor()))
            {
               m_pLengthUnit = &unitMeasure::Yard;
            }
            else if (IsEqual(conversion_factor, unitMeasure::USSurveyYard.GetConvFactor()))
            {
               m_pLengthUnit = &unitMeasure::USSurveyYard;
            }
            else
            {
               ATLASSERT(false); // we don't have a unit of measure for this
            }
         }
         break;
      }
   }
}

void CIfcAlignmentConverter::ConvertToIfc(IBroker* pBroker, const CString& strFilePath, CIfcAlignmentConverter::SchemaType schemaType)
{
   switch (schemaType)
   {
   case Schema_4x1: ConvertToIfc<Ifc4x1>(pBroker, strFilePath); break;
   case Schema_4x2: ConvertToIfc<Ifc4x2>(pBroker, strFilePath); break;
   case Schema_4x3_rc1: ConvertToIfc<Ifc4x3_rc1>(pBroker, strFilePath); break;
   case Schema_4x3_rc2: ConvertToIfc<Ifc4x3_rc2>(pBroker, strFilePath); break;
   default:
      ATLASSERT(false); // is there a new schema type
   }
}


template <typename Schema>
void CIfcAlignmentConverter::ConvertToIfc(IBroker* pBroker, const CString& strFilePath)
{
   USES_CONVERSION;

   IfcHierarchyHelper<Schema> file;
   file.header().file_name().name(T2A(strFilePath));

   //auto project = file.addProject(); // Don't like the default units in IfcOpenShell so we have do build our own
   /////////////////////////// The following is copied from addProject and tweeked
   IfcEntityList::ptr units(new IfcEntityList);
   Schema::IfcDimensionalExponents* dimexp = new Schema::IfcDimensionalExponents(0, 0, 0, 0, 0, 0, 0);
   Schema::IfcSIUnit* unit1 = new Schema::IfcSIUnit(Schema::IfcUnitEnum::IfcUnit_LENGTHUNIT, boost::none, Schema::IfcSIUnitName::IfcSIUnitName_METRE);
   Schema::IfcSIUnit* unit2 = new Schema::IfcSIUnit(Schema::IfcUnitEnum::IfcUnit_PLANEANGLEUNIT, boost::none, Schema::IfcSIUnitName::IfcSIUnitName_RADIAN);

   units->push(unit1);
   units->push(unit2);

   Schema::IfcUnitAssignment* unit_assignment = new Schema::IfcUnitAssignment(units);

   Schema::IfcRepresentationContext::list::ptr rep_contexts(new Schema::IfcRepresentationContext::list);
   Schema::IfcProject* project = new Schema::IfcProject(IfcParse::IfcGlobalId(), nullptr, boost::none, boost::none, boost::none, boost::none, boost::none, rep_contexts, unit_assignment);

   file.addEntity(dimexp);
   file.addEntity(unit1);
   file.addEntity(unit2);
   file.addEntity(unit_assignment);
   file.addEntity(project);
   ///////////////////////////////////////// end of copy from addProject

   auto site = file.addSite(project);

   GET_IFACE2(pBroker, IProjectProperties, pProjectProperties);
   project->setName(T2A(pProjectProperties->GetBridgeName()));
   site->setName(T2A(pProjectProperties->GetBridgeName()));

   auto owner_history = file.getSingle<Schema::IfcOwnerHistory>();
   GET_IFACE2(pBroker, IVersionInfo, pVersionInfo);
   owner_history->OwningApplication()->setApplicationFullName("BridgeLink::PGSuper");
   owner_history->OwningApplication()->setApplicationIdentifier("PGSuper");
   owner_history->OwningApplication()->setVersion(T2A(pVersionInfo->GetVersion(true)));
   owner_history->OwningApplication()->ApplicationDeveloper()->setIdentification("Washington State Department of Transportation, Bridge and Structures Office");
   owner_history->OwningApplication()->ApplicationDeveloper()->setName("Richard Brice, PE");

   auto local_placement = file.addLocalPlacement();

   auto geometric_representation_context = new Schema::IfcGeometricRepresentationContext(boost::none, boost::none, 3, boost::none, local_placement, nullptr);
   file.addEntity(geometric_representation_context); // ADD THE CONTEXT TO THE FILE!!!!
   auto contexts = project->RepresentationContexts(); // get the old context
   contexts->push(geometric_representation_context); // add the new context
   project->setRepresentationContexts(contexts); // set the context back into the project

   auto horizontal_alignment = BuildAlignment<Schema>(pBroker, file);
   auto vertical_profile = BuildProfile<Schema>(pBroker, file);

   //////////////////////////////////////////////////////////////////////////////////

   // create an alignment curve in 3D using the horizontal alignment and vertical profile (nullptr for profile since it's optional)
   // we could also use an IfcOffsetCurveByDistances or an IfcPolyline to define the reference system curve
#pragma Reminder("WORKING HERE - IfcAlignmentCurve is marked as depreciated in 4x3rc2 so what is the alternative?")
   auto axis = new Schema::IfcAlignmentCurve(horizontal_alignment, vertical_profile, boost::none);

   // create an alignment... this is the linear positioning element. it could be an alignment curve (as used above) or IfcOffsetCurveByDistances or IfcPolyline
   auto alignment = new Schema::IfcAlignment(IfcParse::IfcGlobalId(), owner_history, std::string("Alignment"), boost::none, boost::none, local_placement, file.getSingle<Schema::IfcProductRepresentation>(), axis, boost::none);
   file.addEntity(alignment);

   boost::shared_ptr<IfcTemplatedEntityList<Schema::IfcProduct>> related_elements(new IfcTemplatedEntityList<Schema::IfcProduct>());
   related_elements->push(alignment);

   auto spatial_structure = new Schema::IfcRelContainedInSpatialStructure(IfcParse::IfcGlobalId(), owner_history, boost::none, boost::none, related_elements, site);
   file.addEntity(spatial_structure);

   site->ContainsElements()->push(spatial_structure);

   std::ofstream ofs(T2A(strFilePath));
   ofs << file;
}

//LX::CrossSects* CIfcAlignmentConverter::CreateCrossSections(IBroker* pBroker, LX::IFactory* pFactory)
//{
//   // This method creates
//   // Alignment + CrossSects + CrossSect(station) + DesignCrossSectSurface + CrossSectPnt
//   //                                                                      + CrossSectPnt
//   //                                                                      + CrossSectPnt
//   //                                                                      + CrossSectPnt
//   //
//   // For each cross section defined in the PGSuper alignment
//
//   // NOTE: The converter cannot yet read the cross section from a LandXML file. So that
//   // there is consistent functionality between the importer and exporter, this data
//   // will not be exported in this version
//
//   ATLASSERT(false); // SHOULD NEVER GET HERE (see note above)
//
//
//   USES_CONVERSION;
//
//
//   LX::CrossSects* pCrossSections = pFactory->createCrossSects();
//   pCrossSections->setDesc(_T("Cross sections exported from PGSuper"));
//   /*
//
//   GET_IFACE2(pBroker,IRoadway,pPGSuperAlignment);
//   GET_IFACE2(pBroker,IRoadwayData,pRoadwayData);
//   GET_IFACE2(pBroker,IEAFDisplayUnits,pDisplayUnits);
//   GET_IFACE2(pBroker,IBridge,pBridge);
//
//   // get offset to left and right side of the bridge
//   Float64 left_offset  = 0;
//   Float64 right_offset = 0;
//   PierIndexType nPiers = pBridge->GetPierCount();
//   for (PierIndexType pierIdx = 0; pierIdx < nPiers; pierIdx++ )
//   {
//   Float64 left_edge_offset  = pBridge->GetLeftSlabEdgeOffset(pierIdx);
//   Float64 right_edge_offset = pBridge->GetRightSlabEdgeOffset(pierIdx);
//
//   left_offset  = Min(left_offset, left_edge_offset);
//   right_offset = Max(right_offset,right_edge_offset);
//   }
//
//   // make cross section a little wider than the bridge
//   left_offset  *= 1.05;
//   right_offset *= 1.05;
//
//   const RoadwaySectionData& sectionData = pRoadwayData->GetRoadwaySectionData();
//   auto iter = std::cbegin(sectionData.Superelevations);
//   auto end = std::cend(sectionData.Superelevations);
//   for ( ; iter != end; iter++)
//   {
//   const auto& crownData = *iter;
//
//   LX::CrossSect* pCrossSection = pFactory->createCrossSect();
//   Float64 station = ::ConvertFromSysUnits(crownData.Station,unitMeasure::Feet);
//   pCrossSection->setSta(station);
//
//   CString strName;
//   strName.Format(_T("%s"),FormatStation(pDisplayUnits->GetStationFormat(),crownData.Station));
//   pCrossSection->setDesc( strName.GetBuffer() );
//
//
//   LX::DesignCrossSectSurf* pDesignCrossSectionSurface = pFactory->createDesignCrossSectSurf();
//   pDesignCrossSectionSurface->setName( LX::String(_T("Finished roadway surface")) );
//   pDesignCrossSectionSurface->setSide(LX::EnumSideofRoadType::Values::k_both);
//
//   LX::CrossSectPnt* pLeftPnt   = pFactory->createCrossSectPnt();
//   LX::CrossSectPnt* pCenterPnt = pFactory->createCrossSectPnt();
//   LX::CrossSectPnt* pRightPnt  = pFactory->createCrossSectPnt();
//
//   Float64 offset, elevation;
//
//   // left edge of deck
//   offset = ::ConvertFromSysUnits(left_offset,unitMeasure::Feet);
//   elevation = pPGSuperAlignment->GetElevation(crownData.Station, left_offset );
//   elevation = ::ConvertFromSysUnits(elevation,unitMeasure::Feet);
//   pLeftPnt->setDataFormat(LX::EnumDataFormatType::Values::k_Offset_Elevation);
//   pLeftPnt->addItem(offset);
//   pLeftPnt->addItem(elevation);
//   pLeftPnt->setAlignRef(ALIGNMENT_NAME);
//   pDesignCrossSectionSurface->CrossSectPnt().addItem(pLeftPnt);
//
//   if ( crownData.CrownPointOffset < 0 )
//   {
//   // profile grade point is left of crown point
//   LX::CrossSectPnt* pPGLPnt = pFactory->createCrossSectPnt();
//   offset = 0.0;
//   elevation = pPGSuperAlignment->GetElevation(crownData.Station,offset);
//   elevation = ::ConvertFromSysUnits(elevation,unitMeasure::Feet);
//   pPGLPnt->setDataFormat(LX::EnumDataFormatType::Values::k_Offset_Elevation);
//   pPGLPnt->addItem(offset);
//   pPGLPnt->addItem(elevation);
//   pPGLPnt->setAlignRef(ALIGNMENT_NAME);
//   pDesignCrossSectionSurface->CrossSectPnt().addItem(pPGLPnt);
//   }
//
//   // crown point
//   offset = ::ConvertFromSysUnits(crownData.CrownPointOffset,unitMeasure::Feet);
//   elevation = pPGSuperAlignment->GetElevation(crownData.Station, crownData.CrownPointOffset );
//   elevation = ::ConvertFromSysUnits(elevation,unitMeasure::Feet);
//   pCenterPnt->setDataFormat(LX::EnumDataFormatType::Values::k_Offset_Elevation);
//   pCenterPnt->addItem(offset);
//   pCenterPnt->addItem(elevation);
//   pCenterPnt->setAlignRef(ALIGNMENT_NAME);
//   pDesignCrossSectionSurface->CrossSectPnt().addItem(pCenterPnt);
//
//   if ( 0 < crownData.CrownPointOffset )
//   {
//   // profile grade point is right of crown point
//   LX::CrossSectPnt* pPGLPnt = pFactory->createCrossSectPnt();
//   offset = 0.0;
//   elevation = pPGSuperAlignment->GetElevation(crownData.Station,offset);
//   elevation = ::ConvertFromSysUnits(elevation,unitMeasure::Feet);
//   pPGLPnt->setDataFormat(LX::EnumDataFormatType::Values::k_Offset_Elevation);
//   pPGLPnt->addItem(offset);
//   pPGLPnt->addItem(elevation);
//   pPGLPnt->setAlignRef(ALIGNMENT_NAME);
//   pDesignCrossSectionSurface->CrossSectPnt().addItem(pPGLPnt);
//   }
//
//   // right edge of deck
//   offset = ::ConvertFromSysUnits(right_offset,unitMeasure::Feet);
//   elevation = pPGSuperAlignment->GetElevation(crownData.Station, right_offset );
//   elevation = ::ConvertFromSysUnits(elevation,unitMeasure::Feet);
//   pRightPnt->setDataFormat(LX::EnumDataFormatType::Values::k_Offset_Elevation);
//   pRightPnt->addItem(offset);
//   pRightPnt->addItem(elevation);
//   pRightPnt->setAlignRef(ALIGNMENT_NAME);
//   pDesignCrossSectionSurface->CrossSectPnt().addItem(pRightPnt);
//
//   pCrossSection->DesignCrossSectSurf().addItem(pDesignCrossSectionSurface);
//   pCrossSections->CrossSect().addItem(pCrossSection);
//   }
//   */
//   return pCrossSections;
//}
//
//LX::Roadway* CIfcAlignmentConverter::CreateRoadway(IBroker* pBroker, LX::IFactory* pFactory)
//{
//   LX::Roadway* pRoadway = pFactory->createRoadway();
//   pRoadway->setName(_T("Roadway"));
//   LX::StringCollection* pAlignRefs = pFactory->createStringCollection();
//   pAlignRefs->addItem(ALIGNMENT_NAME);
//   pRoadway->setAlignmentRefs(pAlignRefs);
//
//   GET_IFACE2(pBroker, IBridge, pBridge);
//   PierIndexType nPiers = pBridge->GetPierCount();
//   Float64 startStation = pBridge->GetPierStation(0);
//   Float64 endStation = pBridge->GetPierStation(nPiers - 1);
//
//   startStation = ::ConvertFromSysUnits(startStation, unitMeasure::Feet);
//   endStation = ::ConvertFromSysUnits(endStation, unitMeasure::Feet);
//
//   LX::BridgeElement* pBridgeElement = pFactory->createBridgeElement();
//   pBridgeElement->setStaStart(startStation);
//   pBridgeElement->setStaEnd(endStation);
//   pRoadway->BridgeElement().addItem(pBridgeElement);
//
//   return pRoadway;
//}
//

HRESULT CIfcAlignmentConverter::ConvertToPGSuper(IBroker* pBroker, CString& strFilePath)
{
   USES_CONVERSION;
   IfcParse::IfcFile file(T2A(strFilePath.GetBuffer()));

   if (!file.good())
   {
      AfxMessageBox(_T("Unable to parse .ifc file"));
      return S_OK;
   }

   AlignmentData2 alignment_data;
   ProfileData2 profile_data;
   RoadwaySectionData section_data;

   m_Notes.clear();

   auto strSchemaName = file.schema()->name();
   bool bResult = false;
   if (strSchemaName == std::string("IFC4X1"))
   {
      bResult = ConvertToPGSuper<Ifc4x1>(file, &alignment_data, &profile_data, &section_data);
   }
   else if (strSchemaName == std::string("IFC4X2"))
   {
      bResult = ConvertToPGSuper<Ifc4x2>(file, &alignment_data, &profile_data, &section_data);
   }
   else if (strSchemaName == std::string("IFC4X3_RC1"))
   {
      bResult = ConvertToPGSuper<Ifc4x3_rc1>(file, &alignment_data, &profile_data, &section_data);
   }
   //else if (strSchemaName == std::string("IFC4X3_RC2"))
   //{
#pragma Reminder("WORKING HERE - Still can't read 4x3 rc2 files")
      // ConvertToPGSuper has pre-rc2 data types and constructs - need to generalize
      //bResult = ConvertToPGSuper<Ifc4x3_rc2>(file, &alignment_data, &profile_data, &section_data);
   //}
   else
   {
      ATLASSERT(false); // is there a new schema?
   }

   auto notes = GetNotes();
   std::_tstring strNotes;
   for (auto note : notes)
   {
      strNotes += note + _T("\n");
   }
   AfxMessageBox(strNotes.c_str(), MB_OK);

   if (bResult)
   {
      GET_IFACE2(pBroker, IEvents, pEvents);
      pEvents->HoldEvents();

      GET_IFACE2(pBroker, IRoadwayData, pRoadwayData);
      pRoadwayData->SetAlignmentData2(alignment_data);
      pRoadwayData->SetProfileData2(profile_data);
      pRoadwayData->SetRoadwaySectionData(section_data);

      pEvents->FirePendingEvents();
   }

   return S_OK;
}

template <typename Schema>
bool CIfcAlignmentConverter::ConvertToPGSuper(IfcParse::IfcFile& file, AlignmentData2* pAlignmentData, ProfileData2* pProfileData, RoadwaySectionData* pRoadwaySectionData)
{
   InitUnits<Schema>(file);

   auto alignment = GetAlignment<Schema>(file);
   if (alignment == nullptr)
      return false;

   LoadAlignment<Schema>(alignment);
   *pAlignmentData = m_AlignmentData;

   LoadProfile<Schema>(alignment);
   *pProfileData = m_ProfileData;

#pragma Reminder("WORKING HERE - Roadway Section Data")
   // this is dummy data
   RoadwaySectionTemplate roadway_template;
   roadway_template.LeftSlope = -0.02;
   roadway_template.RightSlope = -0.02;
   roadway_template.Station = 0;
   m_RoadwaySectionData.NumberOfSegmentsPerSection = 2;
   m_RoadwaySectionData.ControllingRidgePointIdx = 1;
   m_RoadwaySectionData.RoadwaySectionTemplates.push_back(roadway_template);
   *pRoadwaySectionData = m_RoadwaySectionData;

   return true;
}

std::vector<std::_tstring> CIfcAlignmentConverter::GetNotes()
{
   return m_Notes;
}

template <typename Schema>
typename Schema::IfcAlignment* CIfcAlignmentConverter::GetAlignment(IfcParse::IfcFile& file)
{
   USES_CONVERSION;

   Schema::IfcAlignment::list::ptr alignments = file.instances_by_type<Schema::IfcAlignment>();
   std::vector<Schema::IfcAlignment*> valid_alignments;

   for (auto alignment : *alignments)
   {
      if (IsValidAlignment<Schema>(alignment))
         valid_alignments.push_back(alignment);
   }

   if (valid_alignments.size() == 0)
   {
      AfxMessageBox(_T("File does not contain alignments that are compatible with this software."), MB_OK);
   }
   else
   {
      std::ostringstream os;
      for (auto alignment : valid_alignments)
      {
         auto strLabel = (alignment->hasName() ? alignment->Name() : alignment->hasDescription() ? alignment->Description() : "Unnamed");
         os << strLabel << std::endl;
      }
      int result = AfxChoose(_T("Select Alignment"), _T("Select alignment to import"), A2T(os.str().c_str()), 0, TRUE);
      if (result < 0)
         return nullptr; // dialog was canceled
      else
         return valid_alignments[result];
   }


   return nullptr;
}

template <typename Schema>
void CIfcAlignmentConverter::LoadAlignment(typename Schema::IfcAlignment* pAlignment)
{
   m_bAlignmentStarted = false; // the alignment datablock has not yet been started

   // initalize the alignment data
   m_AlignmentData.Direction = 0.00;
   m_AlignmentData.xRefPoint = 0.00;
   m_AlignmentData.yRefPoint = 0.00;
   m_AlignmentData.HorzCurves.clear();

   auto axis = pAlignment->Axis();
   auto curve = axis->as<Schema::IfcAlignmentCurve>();
   auto horizontal = curve->Horizontal();

   Float64 current_station; // station at the start of the current element
   if(horizontal->hasStartDistAlong())
   {
      // as I understand IFC 8.7.3.1, StartDistAlong is the value of the distance along at the start of the alignment... that seems like a starting station
      current_station = ::ConvertToSysUnits(horizontal->StartDistAlong(), *m_pLengthUnit);
   }
   else
   {
      current_station = 0.0;
   }

   // alignment is made up of Line, Spiral, and/or Curve elements
   auto segments = horizontal->Segments();
   auto begin = segments->begin();
   auto iter = begin;
   auto end = segments->end();
   for (; iter != end; iter++)
   {
      auto segment(*iter);
      bool bIsThereANextSegment = ((iter+1) != end);

      auto linear = segment->CurveGeometry()->as<Schema::IfcLineSegment2D>();
      auto transition = segment->CurveGeometry()->as<Schema::IfcTransitionCurveSegment2D>();
      auto curve = segment->CurveGeometry()->as<Schema::IfcCircularArcSegment2D>();

      Schema::IfcTransitionCurveSegment2D* entrySpiral = nullptr;
      Schema::IfcTransitionCurveSegment2D* exitSpiral = nullptr;

      Float64 end_station = current_station;

      if (linear)
      {
         end_station = OnLine<Schema,Schema::IfcLineSegment2D>(current_station, linear);
      }
      else if (transition)
      {
         // PGSuper can only handle
         // Spiral-Curve
         // Spiral-Curve-Spiral
         // Curve-Spiral
         //
         // Curve-Spiral-Curve, where Spiral is a transition spiral with the start and end radius equal
         // to the curve radii... PGSuper cannot do this case


         entrySpiral = transition;
         if (bIsThereANextSegment)
         {
            // if there is a next segment, check to see if it is a curve
            iter++; // advance to next segment
            curve = (*iter)->CurveGeometry()->as<Schema::IfcCircularArcSegment2D>();
            if (curve)
            {
               // it's a curve... is there an element that follows the curve?
               bIsThereANextSegment = ((iter+1) != end);
               if (bIsThereANextSegment)
               {
                  // if there is a next segment, see if it is a spiral
                  exitSpiral = (*(iter+1))->CurveGeometry()->as<Schema::IfcTransitionCurveSegment2D>();

                  // if not a spiral, pExitSpiral will be nullptr
                  // this is OK, it just means we have a Spiral-Curve situation
                  if (exitSpiral)
                  {
                     // it is an exit spiral, so advance the iterator
                     iter++;
                     bIsThereANextSegment = ((iter+1) != end);
                  }
               }

               end_station = OnCurve<Schema, Schema::IfcTransitionCurveSegment2D, Schema::IfcCircularArcSegment2D>(current_station, entrySpiral, curve, exitSpiral);
            }
            else
            {
               IFC_THROW(_T("A curve must follow a spiral")); // because PGSuper can't handle it otherwise
               ATLASSERT(false); // a curve must follow a spiral
            }
         }
         else
         {
            m_Notes.push_back(std::_tstring(_T("Element ignored: The last element in the alignment cannot be a Spiral."))); // PGSuper can't model a lone spiral
         }
      }
      else if (curve)
      {
         // looking for Curve-Spiral case
         if (bIsThereANextSegment)
         {
            // check to see if the next element is a spiral
            exitSpiral = (*(iter + 1))->CurveGeometry()->as<Schema::IfcTransitionCurveSegment2D>();
            
            // if not a spiral, pExitSpiral will be nullptr
            // this is OK, it just means we have a Spiral-Curve situation
            if (exitSpiral)
            {
               iter++;
               bIsThereANextSegment = ((iter+1) != end);
            }

            // Check if the next object is a Curve and if
            // the exit spiral and the curve are touching
            if (bIsThereANextSegment)
            {
               auto next_curve = (*(iter))->CurveGeometry()->as<Schema::IfcCircularArcSegment2D>();
               if (next_curve && exitSpiral)
               {
                  CComPtr<IPoint2d> pntSpiralStart, pntSpiralPI, pntSpiralEnd;
                  GetSpiralPoints<Schema>(exitSpiral, &pntSpiralStart, &pntSpiralPI, &pntSpiralEnd);

                  CComPtr<IPoint2d> pntNextCurveStart, pntNextCurvePI, pntNextCurveEnd, pntNextCurveCenter;
                  GetCurvePoints<Schema,Schema::IfcCircularArcSegment2D>(next_curve, &pntNextCurveStart, &pntNextCurvePI, &pntNextCurveEnd, &pntNextCurveCenter);

                  if (SameLocation(pntSpiralEnd,pntNextCurveStart,m_Precision) == S_FALSE)
                  {
                     IFC_THROW(_T("PGSuper cannot model a transition spiral between two circular curves."));
                  }
               }
            }
         }

         ATLASSERT(entrySpiral == nullptr); // this must be the case
         end_station = OnCurve<Schema, Schema::IfcTransitionCurveSegment2D, Schema::IfcCircularArcSegment2D>(current_station, entrySpiral, curve, exitSpiral);
      }
      else
      {
         ATLASSERT(false);
      }
      current_station = end_station;
   }
}

template <typename Schema>
void CIfcAlignmentConverter::LoadProfile(typename Schema::IfcAlignment* pAlignment)
{
   m_ProfileState = PROFILE_NOT_STARTED;
   m_ProfileData.VertCurves.clear();

   auto axis = pAlignment->Axis();
   auto curve = axis->as<Schema::IfcAlignmentCurve>();
   auto vertical = curve->hasVertical() ? curve->Vertical() : nullptr;

   auto segments = vertical ? vertical->Segments() : nullptr;

  if (vertical == nullptr || segments->size() == 0)
   {
      // the profile geometry list is empty so assume a flat grade
      m_Notes.push_back(std::_tstring(_T("A profile was not found or the profile does not contain segments. Assuming a default profile.")));
      m_ProfileData.Station = 0.0;
      m_ProfileData.Elevation = 0.0;
      m_ProfileData.Grade = 0.0;
      return;
   }

   for (auto segment : *segments)
   {
      auto linear_segment = segment->as<Schema::IfcAlignment2DVerSegLine>();
      auto parabolic_arc = segment->as<Schema::IfcAlignment2DVerSegParabolicArc>();
      auto circular_arc = segment->as<Schema::IfcAlignment2DVerSegCircularArc>();

      if (linear_segment)
      {
         LinearSegment<Schema, Schema::IfcAlignment2DVerSegLine>(linear_segment);
      }
      else if (parabolic_arc)
      {
         ParabolicSegment<Schema, Schema::IfcAlignment2DVerSegParabolicArc>(parabolic_arc);
      }
      else if (circular_arc)
      {
#pragma Reminder("WORKING HERE - Need to deal with circular arcs") // treat it as a parabola for now
         parabolic_arc = (Schema::IfcAlignment2DVerSegParabolicArc*)circular_arc;
         ParabolicSegment<Schema, Schema::IfcAlignment2DVerSegParabolicArc>(parabolic_arc);
         // TODO: provide a better exception
         //IFC_THROW(_T("Circular curve element was found in the profile definition. PGSuper does not support circular vertical curves."));
      }
      else
      {
         ATLASSERT(false); // is there a new type ???
                           // TODO: provide a better exception
         IFC_THROW(_T("An unknown profile element was encountered"));
      }
   }

   if (m_ProfileState != PROFILE_ESTABLISHED)
   {
      // we are out of elements and the profile definition is not finished

      // THIS IS A GOOD PLACE TO GENERATE INFORMATION MESSAGES ABOUT ANY ASSUMPTIONS
      // THIS IMPORTER HAD TO MAKE
      if (m_ProfileData.VertCurves.size() == 0)
      {
         // A second element was not provided so the main grade
         // could not be established... use the default value of 0
         m_Notes.push_back(std::_tstring(_T("More elements are needed to determine the starting grade of the profile. Assuming a grade of 0.0%")));
      }
      else
      {
         // The exit profile of the last vertical curve could not be established...
         // We could either pop the last vertical curve out of the list or use the
         // default exit grade of 0.
         //
         // Use the default exit grade of 0.
         m_Notes.push_back(std::_tstring(_T("More elements are needed to determine the exit grade of the last vertical curve. Assuming a grade of 0.0%")));
      }
   }
}

//void CIfcAlignmentConverter::LoadCrossSections(LX::CrossSects* pCrossSects, LX::String& strSurfaceName)
//{
//   // NOTE: Not sure how to read superelevation and cross section data from LandXML files
//   // This version does not support importing data for PGSuper cross sections
//   //
//   // This is prototype code and doesn't work.
//
//   ATLASSERT(false); // SHOULD NEVER GET HERE (see note above)
//                     /*
//                     m_RoadwaySectionData.Superelevations.clear();
//
//                     std::vector<DesignCrossSectData> vSectionData; // collect the section data so that it can be massaged into PGSuper format
//
//                     LX::CrossSectCollectionIterator* pCrossSectCollectionIter = pCrossSects->CrossSect().iterator();
//                     while ( !pCrossSectCollectionIter->atEnd() )
//                     {
//                     LX::CrossSect* pCrossSect = pCrossSectCollectionIter->current();
//
//                     LX::DesignCrossSectSurfCollectionIterator* pDesignCrossSectSurfIter = pCrossSect->DesignCrossSectSurf().iterator();
//                     while ( !pDesignCrossSectSurfIter->atEnd() )
//                     {
//                     LX::DesignCrossSectSurf* pDesignCrossSectSurf = pDesignCrossSectSurfIter->current();
//                     if ( pDesignCrossSectSurf->hasValue_Name() && pDesignCrossSectSurf->getName() == strSurfaceName )
//                     {
//                     // this is our surface.. get the cross section points and turn them into PGSuper input
//                     //pDesignCrossSectSurf->CrossSectPnt(); // these are the points we want
//                     // need to find the crown point and compute slope for left and right side and crown point offset
//                     // need to deal with left, right, or both sides
//
//                     DesignCrossSectData sectionData;
//                     sectionData.Station = pCrossSect->getSta();
//                     GetSlopes(pDesignCrossSectSurf,&sectionData);
//                     vSectionData.push_back(sectionData);
//                     }
//                     pDesignCrossSectSurfIter->next();
//                     }
//                     pCrossSectCollectionIter->next();
//                     }
//
//                     // put the sections in order
//                     std::sort(vSectionData.begin(),vSectionData.end(),CIfcAlignmentConverter::Compare);
//                     MergeSections(vSectionData);
//
//                     std::vector<DesignCrossSectData>::iterator iter;
//                     for ( iter = vSectionData.begin(); iter != vSectionData.end(); iter++ )
//                     {
//                     DesignCrossSectData& sd = *iter;
//                     ATLASSERT(sd.Side == LX::EnumSideofRoadType::Values::k_both);
//                     CrownData2 cd;
//                     cd.Station = sd.Station;
//                     cd.Left = sd.LeftSlope;
//                     cd.Right = sd.RightSlope;
//                     cd.CrownPointOffset = sd.CrownPointOffset;
//
//                     m_RoadwaySectionData.Superelevations.push_back(cd);
//                     }
//                     */
//}
//

template <typename Schema,typename Segment>
Float64 CIfcAlignmentConverter::OnLine(Float64 startStation, typename Segment* pLine)
{
   Float64 sx, sy;
   GetPoint<Schema>(pLine->StartPoint(), &sx, &sy);

   Float64 length = pLine->SegmentLength();
   Float64 startDirection = pLine->StartDirection();
   return OnLine(sx, sy, startStation, startDirection, length);
}

Float64 CIfcAlignmentConverter::OnLine(Float64 sx,Float64 sy,Float64 startStation,Float64 startDirection, Float64 length)
{
   Float64 end_station = startStation + length;

   if (!m_bAlignmentStarted)
   {
      // the bridge starts somewhere in this line segment
      m_AlignmentData.RefStation = startStation;
      m_AlignmentData.xRefPoint = sx;
      m_AlignmentData.yRefPoint = sy;
      m_AlignmentData.Direction = startDirection;

      if (IsZero(m_AlignmentData.Direction))
         m_AlignmentData.Direction = 0.0;
      else if (m_AlignmentData.Direction < 0)
         m_AlignmentData.Direction += TWO_PI;

      m_bAlignmentStarted = true;
   }
   else
   {
      if (m_LastAlignmentType != Curve)
      {
         // add an angle point
         HorzCurveData hcData;
         hcData.PIStation = startStation;
         hcData.Radius = 0;
         hcData.EntrySpiral = 0;
         hcData.ExitSpiral = 0;
         hcData.bFwdTangent = true;
         hcData.FwdTangent = startDirection;

         m_AlignmentData.HorzCurves.push_back(hcData);
      }
   }

   m_LastAlignmentType = Line;

   return end_station;
}

template <typename Schema,typename SpiralType, typename CurveType>
Float64 CIfcAlignmentConverter::OnCurve(Float64 startStation, typename SpiralType* pEntrySpiral, typename CurveType* pCurve, typename SpiralType* pExitSpiral)
{
   ATLASSERT(pCurve != nullptr);

   Float64 radius = pCurve->Radius();

   // Get all the construction points
   CComPtr<IPoint2d> pntEntryStart, pntEntryPI, pntEntryEnd;
   CComPtr<IPoint2d> pntExitStart, pntExitPI, pntExitEnd;
   CComPtr<IPoint2d> pntCurveStart, pntCurvePI, pntCurveEnd, pntCurveCenter;

   Float64 entry_spiral_length = 0;
   Float64 exit_spiral_length = 0;
   if (pEntrySpiral)
   {
      GetSpiralPoints<Schema>(pEntrySpiral, &pntEntryStart, &pntEntryPI, &pntEntryEnd);
      entry_spiral_length = ::ConvertToSysUnits(pEntrySpiral->SegmentLength(), *m_pLengthUnit);

      if (pntEntryStart == nullptr || pntEntryPI == nullptr || pntEntryEnd == nullptr)
      {
         m_Notes.push_back(std::_tstring(_T("Entry spiral ignored.")));
         pntEntryStart.Release();
         pntEntryPI.Release();
         pntEntryEnd.Release();
         pEntrySpiral = nullptr;
         entry_spiral_length = 0;
      }
      else
      {
         if(pEntrySpiral->hasStartRadius() && !IsZero(pEntrySpiral->StartRadius()))
         {
            m_Notes.push_back(std::_tstring(_T("Start radius of entry spiral taken to be infinite")));
         }
         CheckSpiralType<Schema,SpiralType>(pEntrySpiral);
      }
   }

   GetCurvePoints<Schema,CurveType>(pCurve, &pntCurveStart, &pntCurvePI, &pntCurveEnd, &pntCurveCenter);

   if (pntCurveStart == nullptr || pntCurvePI == nullptr || pntCurveEnd == nullptr || pntCurveCenter == nullptr)
   {
      m_Notes.push_back(std::_tstring(_T("Zero radius curve could not be constructed.")));
      return startStation;
   }

   if (pExitSpiral)
   {
      GetSpiralPoints<Schema>(pExitSpiral, &pntExitStart, &pntExitPI, &pntExitEnd);
      exit_spiral_length = ::ConvertToSysUnits(pExitSpiral->SegmentLength(), *m_pLengthUnit);

      if (pntExitStart == nullptr || pntExitPI == nullptr || pntExitEnd == nullptr)
      {
         m_Notes.push_back(std::_tstring(_T("Exit spiral ignored.")));
         pntExitStart.Release();
         pntExitPI.Release();
         pntExitEnd.Release();
         pExitSpiral = nullptr;
         exit_spiral_length = 0;
      }
      else
      {
         if (pExitSpiral->hasEndRadius() && !IsZero(pExitSpiral->EndRadius()))
         {
            m_Notes.push_back(std::_tstring(_T("End radius of exit spiral taken to be infinite")));
         }
         CheckSpiralType<Schema,SpiralType>(pExitSpiral);
      }
   }

   Float64 sx, sy;
   pntCurveStart->get_X(&sx);
   pntCurveStart->get_Y(&sy);

#if defined _DEBUG
   // check radius based on the circular curve parameters
   Float64 ex, ey;
   Float64 cx, cy;
   Float64 px, py;

   pntCurvePI->get_X(&px);
   pntCurvePI->get_Y(&py);

   pntCurveEnd->get_X(&ex);
   pntCurveEnd->get_Y(&ey);

   pntCurveCenter->get_X(&cx);
   pntCurveCenter->get_Y(&cy);

   Float64 dx = sx - cx;
   Float64 dy = sy - cy;
   ATLASSERT(IsEqual(radius, sqrt(dx*dx + dy*dy)));
#endif // _DEBUG

   // Determine the control points
   CComPtr<IPoint2d> pntStart, pntPI, pntEnd;
   if (pEntrySpiral && !pExitSpiral)
   {
      // Spiral-Curve
      pntStart = pntEntryStart;

      // PI is at the intersection of the forward and back tangents
      CComPtr<IIntersect2> intersect;
      m_CogoEngine->get_Intersect(&intersect);
      intersect->LinesByPoints(pntEntryStart, pntEntryPI, 0.0, pntCurvePI, pntCurveEnd, 0.0, &pntPI);

      pntEnd = pntCurveEnd;

      if (!IsEqual(::ConvertToSysUnits(pEntrySpiral->EndRadius(), *m_pLengthUnit), radius))
      {
         m_Notes.push_back(std::_tstring(_T("End radius of the entry sprial does not match the radius of the circular curve. The entry spiral end radius will be ignored.")));
      }

      if (SameLocation(pntEntryEnd,pntCurveStart,m_Precision) == S_FALSE)
      {
         m_Notes.push_back(std::_tstring(_T("The end of the entry spiral does not coincide with the start of the circular curve. The end of the entry spiral has been adjusted.")));
      }
   }
   else if (!pEntrySpiral && pExitSpiral)
   {
      // Curve-Spiral
      pntStart = pntCurveStart;

      // PI is at the intersection of the forward and back tangents
      CComPtr<IIntersect2> intersect;
      m_CogoEngine->get_Intersect(&intersect);
      intersect->LinesByPoints(pntStart, pntCurvePI, 0.0, pntExitPI, pntExitEnd, 0.0, &pntPI);

      pntEnd = pntExitEnd;

      if (!IsEqual(::ConvertToSysUnits(pExitSpiral->StartRadius(), *m_pLengthUnit), radius))
      {
         m_Notes.push_back(std::_tstring(_T("Start radius of the exit sprial does not match the radius of the circular curve. The exit spiral start radius will be ignored.")));
      }

      if (SameLocation(pntCurveEnd,pntExitStart,m_Precision) == S_FALSE)
      {
         m_Notes.push_back(std::_tstring(_T("The start of the exit spiral does not coincide with the end of the circular curve. The exit spiral has been adjusted.")));
      }
   }
   else if (pEntrySpiral && pExitSpiral)
   {
      // Spiral-Curve-Spiral
      pntStart = pntEntryStart;

      CComPtr<IIntersect2> intersect;
      m_CogoEngine->get_Intersect(&intersect);
      intersect->LinesByPoints(pntEntryStart, pntEntryPI, 0.0, pntExitPI, pntExitEnd, 0.0, &pntPI);

      ATLASSERT(pntPI->SameLocation(pntCurvePI));

      pntEnd = pntExitEnd;


      if (!IsEqual(::ConvertToSysUnits(pEntrySpiral->EndRadius(), *m_pLengthUnit), radius))
      {
         m_Notes.push_back(std::_tstring(_T("End radius of the entry sprial does not match the radius of the circular curve. The entry spiral end radius will be ignored.")));
      }

      if (!IsEqual(::ConvertToSysUnits(pExitSpiral->StartRadius(), *m_pLengthUnit), radius))
      {
         m_Notes.push_back(std::_tstring(_T("Start radius of the exit sprial does not match the radius of the circular curve. The exit spiral start radius will be ignored.")));
      }

      if (SameLocation(pntEntryEnd,pntCurveStart,m_Precision) == S_FALSE)
      {
         m_Notes.push_back(std::_tstring(_T("The end of the entry spiral does not coincide with the start of the circular curve. The entry spiral has been adjusted.")));
      }

      if (SameLocation(pntCurveEnd,pntExitStart,m_Precision) == S_FALSE)
      {
         m_Notes.push_back(std::_tstring(_T("The start of the exit spiral does not coincide with the end of the circular curve. The exit spiral has been adjusted.")));
      }
   }
   else
   {
      // Curve
      pntStart = pntCurveStart;
      pntPI = pntCurvePI;
      pntEnd = pntCurveEnd;
   }

   // create a horizontal curve object so that we can get some information from it
   CComPtr<IHorzCurve> hc;
   hc.CoCreateInstance(CLSID_HorzCurve);
   hc->putref_PBT(pntStart);
   hc->putref_PI(pntPI);
   hc->putref_PFT(pntEnd);
   hc->put_Radius(radius);
   hc->put_SpiralLength(spEntry, entry_spiral_length);
   hc->put_SpiralLength(spExit, exit_spiral_length);

   Float64 length;
   hc->get_TotalLength(&length);

   CComPtr<IAngle> objAngle;
   hc->get_CircularCurveAngle(&objAngle);

   Float64 end_station = startStation + length;

   if (!m_bAlignmentStarted)
   {
      // this is the first element... start the alignment data
      m_AlignmentData.RefStation = startStation;
      m_AlignmentData.xRefPoint = sx;
      m_AlignmentData.yRefPoint = sy;

      CComPtr<IDirection> dirBkTangent;
      hc->get_BkTangentBrg(&dirBkTangent);
      dirBkTangent->get_Value(&m_AlignmentData.Direction);

      m_bAlignmentStarted = true;
   }

   HorzCurveData hcData;
   hcData.EntrySpiral = entry_spiral_length;
   hcData.ExitSpiral = exit_spiral_length;
   hcData.Radius = radius;

   Float64 tangent;
   hc->get_BkTangentLength(&tangent);
   hcData.PIStation = startStation + tangent;

   CComPtr<IAngle> curve_angle;
   hc->get_CurveAngle(&curve_angle);

   Float64 delta;
   curve_angle->get_Value(&delta);

   CurveDirectionType dir;
   hc->get_Direction(&dir);

   if (dir == cdRight)
      delta *= -1;

   hcData.bFwdTangent = false;
   hcData.FwdTangent = delta;

   m_AlignmentData.HorzCurves.push_back(hcData);

   m_LastAlignmentType = Curve;

   return end_station;
}

template <typename Schema, typename CurveType>
void CIfcAlignmentConverter::GetCurvePoints(typename CurveType* pCurve, IPoint2d** ppStart, IPoint2d** ppPI, IPoint2d** ppEnd, IPoint2d** ppCenter)
{
   auto pStart = pCurve->StartPoint();
   auto bkTangentBrg = pCurve->StartDirection();
   auto L = pCurve->SegmentLength();
   auto R = pCurve->Radius();

   Float64 delta = L / R;
   Float64 T = R*tan(delta/2);

   Float64 sx, sy;
   GetPoint<Schema>(pStart, &sx, &sy);
   CComPtr<IPoint2d> pntStart;
   pntStart.CoCreateInstance(CLSID_Point2d);
   pntStart->Move(sx, sy);
   pntStart.CopyTo(ppStart);

   CComPtr<ILocate2> locate;
   m_CogoEngine->get_Locate(&locate);

   locate->ByDistDir(*ppStart, T, CComVariant(bkTangentBrg), 0.0, ppPI);

   Float64 fwdTangentBrg = bkTangentBrg + (pCurve->IsCCW() ? 1 : -1)*delta;

   locate->ByDistDir(*ppPI, T, CComVariant(fwdTangentBrg), 0.0, ppEnd);

   locate->ByDistDir(*ppStart, R, CComVariant(bkTangentBrg + (pCurve->IsCCW() ? 1 : -1)*PI_OVER_2), 0.0, ppCenter);
}

Float64 SpiralX(Float64 ls, Float64 angle)
{
   return ls*(1 - pow(angle, 2) / 10 + pow(angle, 4) / 216 - pow(angle, 6) / 9360);
}

Float64 SpiralY(Float64 ls, Float64 angle)
{
   return ls*(angle / 3 - pow(angle, 3) / 42 + pow(angle, 5) / 1320 - pow(angle, 7) / 75600);
}

template <typename Schema,typename SpiralType>
void CIfcAlignmentConverter::GetSpiralPoints(typename SpiralType* pSpiral, IPoint2d** ppStart, IPoint2d** ppPI, IPoint2d** ppEnd)
{
   auto pStart = pSpiral->StartPoint();
   auto bkTangentBrg = pSpiral->StartDirection();
   auto L = pSpiral->SegmentLength();
   auto R = (pSpiral->hasStartRadius() ? pSpiral->StartRadius() : pSpiral->EndRadius());
   bool bIsCCW = (pSpiral->hasStartRadius() ? pSpiral->IsStartRadiusCCW() : pSpiral->IsEndRadiusCCW());

   Float64 sx, sy;
   GetPoint<Schema>(pStart, &sx, &sy);
   CComPtr<IPoint2d> pntStart;
   pntStart.CoCreateInstance(CLSID_Point2d);
   pntStart->Move(sx, sy);
   pntStart.CopyTo(ppStart);

   Float64 DE = L / (2 * R); // deflection angle
   Float64 X = SpiralX(L, DE);
   Float64 Y = SpiralY(L, DE);
   Float64 v = Y / sin(DE); // short tangent
   Float64 u = X - Y / tan(DE); // long tangent

   CComPtr<ILocate2> locate;
   m_CogoEngine->get_Locate(&locate);

   if (pSpiral->hasStartRadius())
   {
      locate->ByDistDir(*ppStart, v, CComVariant(bkTangentBrg), 0.0, ppPI);
      locate->ByDistDir(*ppPI, u, CComVariant(bkTangentBrg + (bIsCCW ? 1 : -1)*DE), 0.0, ppEnd);
   }
   else
   {
      locate->ByDistDir(*ppStart, u, CComVariant(bkTangentBrg), 0.0, ppPI);
      locate->ByDistDir(*ppPI, v, CComVariant(bkTangentBrg + (bIsCCW ? 1 : -1)*DE), 0.0, ppEnd);
   }
}

template <typename Schema,typename LineSegmentType>
void CIfcAlignmentConverter::LinearSegment(typename LineSegmentType* pLinearSegment)
{
   Float64 length = ::ConvertToSysUnits(pLinearSegment->HorizontalLength(),*m_pLengthUnit);
   Float64 start_gradient = pLinearSegment->StartGradient();
   Float64 start_dist = ::ConvertToSysUnits(pLinearSegment->StartDistAlong(),*m_pLengthUnit);

   Float64 start_height = ::ConvertToSysUnits(pLinearSegment->StartHeight(),*m_pLengthUnit);

   if (m_ProfileState == PROFILE_NOT_STARTED)
   {
      m_ProfileData.Station = start_dist;
      m_ProfileData.Elevation = start_height;
      m_ProfileData.Grade = start_gradient;
      m_ProfileState = PROFILE_ESTABLISHED;
   }
   else
   {
      // PGSuper models linear segments as zero length vertical curves
      VertCurveData vcData;
      vcData.PVIStation = start_dist + length/2;
      vcData.L1 = length;
      vcData.L2 = 0;
      vcData.ExitGrade = start_gradient;

      m_ProfileData.VertCurves.push_back(vcData);

      m_ProfileState = PROFILE_ESTABLISHED;
   }
}

template <typename Schema, typename ParabolicSegmentType>
void CIfcAlignmentConverter::ParabolicSegment(typename ParabolicSegmentType* pParaCurve)
{
   // finish any open profile element
   Float64 start_gradient = pParaCurve->StartGradient();
   Float64 start_dist = ::ConvertToSysUnits(pParaCurve->StartDistAlong(), *m_pLengthUnit);
   Float64 start_height = ::ConvertToSysUnits(pParaCurve->StartHeight(), *m_pLengthUnit);
   Float64 length = ::ConvertToSysUnits(pParaCurve->HorizontalLength(), *m_pLengthUnit);
   Float64 R = ::ConvertToSysUnits(pParaCurve->ParabolaConstant(), *m_pLengthUnit);
   if (pParaCurve->IsConvex())
      R *= -1;

   Float64 exit_gradient = length/R + start_gradient;

   if (m_ProfileState == PROFILE_NOT_STARTED)
   {
      m_ProfileData.Station = start_dist;
      m_ProfileData.Elevation = start_height;
      m_ProfileData.Grade = start_gradient;
      m_ProfileState = PROFILE_ESTABLISHED;
   }

   // add this vertical curve
   VertCurveData vcData;
   vcData.PVIStation = start_dist + length/2;
   vcData.L1 = length;
   vcData.L2 = 0;
   vcData.ExitGrade = exit_gradient;

   m_ProfileData.VertCurves.push_back(vcData);
}

template <typename Schema, typename SpiralType>
void CIfcAlignmentConverter::CheckSpiralType(typename SpiralType* pSpiral)
{
   switch (pSpiral->TransitionCurveType())
   {
   case Schema::IfcTransitionCurveType::IfcTransitionCurveType_BIQUADRATICPARABOLA:
   case Schema::IfcTransitionCurveType::IfcTransitionCurveType_BLOSSCURVE:
   case Schema::IfcTransitionCurveType::IfcTransitionCurveType_COSINECURVE:
   case Schema::IfcTransitionCurveType::IfcTransitionCurveType_CUBICPARABOLA:
   case Schema::IfcTransitionCurveType::IfcTransitionCurveType_SINECURVE:
      m_Notes.push_back(std::_tstring(_T("Spiral type not supported. Assuming clothoid")));
      break;

   case Schema::IfcTransitionCurveType::IfcTransitionCurveType_CLOTHOIDCURVE:
      // this is ok... it is what we were expecting
      break;

   default:
      ATLASSERT(false); // is there a new spiral type???
      m_Notes.push_back(std::_tstring(_T("Spiral type not defined. Assuming clothoid")));
      break;
   }
}
//
//
//
//void CIfcAlignmentConverter::GetSlopes(LX::DesignCrossSectSurf* pDesignSurf, DesignCrossSectData* pSectionData)
//{
//   pSectionData->Side = pDesignSurf->hasValue_Side() ? pDesignSurf->getSide() : LX::EnumSideofRoadType::Values::k_both;
//
//   std::map<Float64, Float64> offset_elevation; // sorted list of offset/elevation values
//                                                // what about closed shapes??? There could be multiple elevations at the same offset.
//                                                // how do you know you've got the roadway surface??
//
//                                                // NOTE: This loop doesn't work yet... but it is where I'm working...
//                                                // need to convert this data to single left and right slopes
//   LX::CrossSectPntCollectionIterator* pIter = pDesignSurf->CrossSectPnt().iterator();
//   while (!pIter->atEnd())
//   {
//      LX::CrossSectPnt* pPnt = pIter->current();
//
//      if (!pPnt->hasValue_DataFormat() || pPnt->getDataFormat() == LX::EnumDataFormatType::Values::k_Offset_Elevation)
//      {
//         Float64 offset = pPnt->at(0); // offset from PGL  < 0 = left of PGL
//         Float64 elevation = pPnt->at(1);
//         offset_elevation.insert(std::make_pair(offset, elevation));
//      }
//      else
//      {
//         Float64 slope = pPnt->at(0); // ft/ft? or %???
//         Float64 distance = pPnt->at(1); // cummulate from PGL?
//      }
//
//      pIter->next();
//   }
//
//   // assuming everything is offset/elevation...
//   // find the zero point, this is the PGL
//
//   // lower_bound returns the iterator that is at the position where the key value is equal to or greater than 0.0
//   std::map<Float64, Float64>::iterator found = offset_elevation.lower_bound(0.0);
//   if (found == offset_elevation.end())
//   {
//      ATLASSERT(false); // there isn't a zero point... or anything greater than zero... now what?
//   }
//
//   std::map<Float64, Float64>::iterator left_iter = found;
//   std::map<Float64, Float64>::iterator right_iter = found;
//
//   left_iter--;
//   right_iter++;
//
//   std::pair<Float64, Float64> PGL = *found;
//   std::pair<Float64, Float64> left = *left_iter;
//   std::pair<Float64, Float64> right = *right_iter;
//
//   Float64 left_slope = (left.second - PGL.second) / (PGL.first - left.first);
//   Float64 right_slope = (right.second - PGL.second) / (right.first - PGL.first);
//
//   pSectionData->CrownPointOffset = PGL.first;
//   pSectionData->LeftSlope = left_slope;
//   pSectionData->RightSlope = right_slope;
//}
//
//bool CIfcAlignmentConverter::Compare(const DesignCrossSectData& a, const DesignCrossSectData& b)
//{
//   if (a.Station < b.Station)
//      return true;
//
//   if (b.Station < a.Station)
//      return false;
//
//   return (a.Side < b.Side); // if stations are same, sorts right, left, then both
//}
//
//void CIfcAlignmentConverter::MergeSections(std::vector<DesignCrossSectData>& vSectionData)
//{
//   // there could be several records for a single station... one for left and one for right...
//   // merge these into a single record... when this method is done, all sections in the vector
//   // have Side = both
//}
