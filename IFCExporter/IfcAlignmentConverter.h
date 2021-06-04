#pragma once

#include "F:/IfcOpenShell/_installed-vs2015-x64/include/ifcparse/IfcHierarchyHelper.h"

#include "F:/IfcOpenShell/_installed-vs2015-x64/include/ifcparse/Ifc4x1.h"
#include "F:/IfcOpenShell/_installed-vs2015-x64/include/ifcparse/Ifc4x2.h"
#include "F:/IfcOpenShell/_installed-vs2015-x64/include/ifcparse/Ifc4x3_rc1.h"
#include "F:/IfcOpenShell/_installed-vs2015-x64/include/ifcparse/Ifc4x3_rc2.h"

#include <IFace\Project.h>
#include <MfcTools\Prompts.h>

///////////////////////////////////////////////////////////////////////////
// CIfcAlignmentConverter
//
// Converts data between IFC and PGSuper data structures
class CIfcAlignmentConverter
{
public:
   CIfcAlignmentConverter(void);
   ~CIfcAlignmentConverter(void);

   enum SchemaType
   {
      Schema_4x1,
      Schema_4x2,
      Schema_4x3_rc1,
      Schema_4x3_rc2
   };

   // Converts PGSuper data to Ifc data
   void ConvertToIfc(IBroker* pBroker, const CString& strFilePath,SchemaType schemaType);

   // Converts Ifc data to PGSuper data
   HRESULT ConvertToPGSuper(IBroker* pBroker, CString& strFilePath);

   // Returns a list of notes that were generated during the LandXML to PGSuper conversion process
   std::vector<std::_tstring> GetNotes();

   static Float64 GetPrecision() { return m_Precision; }

private:
   static Float64 m_Precision;
   const unitLength* m_pLengthUnit;
   std::vector<std::_tstring> m_Notes;
   CComPtr<ICogoEngine> m_CogoEngine;
   CComPtr<IGeomUtil2d> m_GeomUtil;

   AlignmentData2 m_AlignmentData;
   ProfileData2   m_ProfileData;
   RoadwaySectionData m_RoadwaySectionData;

   bool m_bAlignmentStarted;
   int m_ProfileState; // -1 = not yet started, 0 = started, but grade not determined, 1 = first point established

   //LX::CrossSects* CreateCrossSections(IBroker* pBroker, LX::IFactory* pFactory);
   //LX::Roadway*    CreateRoadway(IBroker* pBroker, LX::IFactory* pFactory);


   template <typename Schema>
   bool ConvertToPGSuper(IfcParse::IfcFile& file, AlignmentData2* pAlignmentData, ProfileData2* pProfileData, RoadwaySectionData* pRoadwaySectionData);

   template <typename Schema>
   typename Schema::IfcAlignment* GetAlignment(IfcParse::IfcFile& file);

   template <typename Schema>
   void InitUnits(IfcParse::IfcFile& file);

   template <typename Schema>
   void ConvertToIfc(IBroker* pBroker, const CString& strFilePath);

   template <typename Schema>
   void LoadAlignment(typename Schema::IfcAlignment* pAlignment);

   template <typename Schema>
   void LoadProfile(typename Schema::IfcAlignment* pAlignment);
   //void LoadCrossSections(LX::CrossSects* pCrossSects, LX::String& strSurfaceName);

   template <typename Schema>
   void GetPoint(typename Schema::IfcCartesianPoint* pPoint, Float64* pX, Float64* pY)
   {
      auto coordinates = pPoint->Coordinates();
      ATLASSERT(2 <= coordinates.size());
      *pX = ::ConvertToSysUnits(coordinates[0],*m_pLengthUnit);
      *pY = ::ConvertToSysUnits(coordinates[1],*m_pLengthUnit);
   }

   enum LastAlignmentType { Unknown, Line, Curve } m_LastAlignmentType;

   // adds a line to the alignment. returns the station at the end of the line
   template <typename Schema, typename Segment> Float64 OnLine(Float64 startStation, typename Segment* pLine);
   Float64 OnLine(Float64 sx, Float64 sy, Float64 startStation, Float64 startDirection, Float64 length);

   // adds a curve to the alignment. returns the station at the end of the curve
   template <typename Schema,typename Spiral,typename Curve>
   Float64 OnCurve(Float64 startStation, typename Spiral* pEntrySpiral, typename Curve* pCurve, typename Spiral* pExitSpiral);

   // adds linear segment to the profile
   template <typename Schema, typename LineSegmentType>
   void LinearSegment(typename LineSegmentType* pLinearSegment);

   // adds a parabolic curve to the proflie
   template <typename Schema, typename ParabolicSegmentType>
   void ParabolicSegment(typename ParabolicSegmentType* pParaCurve);

   template <typename Schema,typename CurveType>
   void GetCurvePoints(typename CurveType* pCurve, IPoint2d** ppStart, IPoint2d** ppPI, IPoint2d** ppEnd, IPoint2d** ppCenter);

   template <typename Schema,typename SpiralType>
   void GetSpiralPoints(typename SpiralType* pSpiral, IPoint2d** ppStart, IPoint2d** ppPI, IPoint2d** ppEnd);

   template <typename Schema, typename SpiralType>
   void CheckSpiralType(typename SpiralType* pSpiral);

   //// NOTE: Not converting LandXML to cross slope data. This struct and the functions that follow are
   //// trial implementations, but they don't work.
   //struct DesignCrossSectData
   //{
   //   Float64 Station;
   //   LX::EnumSideofRoadType::Values Side;

   //   Float64 CrownPointOffset;
   //   Float64 LeftSlope;
   //   Float64 RightSlope;
   //};
   //void GetSlopes(LX::DesignCrossSectSurf* pDesignSurf, DesignCrossSectData* pSectionData);
   //static bool Compare(const DesignCrossSectData& a, const DesignCrossSectData& b);
   //void MergeSections(std::vector<DesignCrossSectData>& vSectionData);


   template <typename Schema>
   bool IsValidAlignment(typename Schema::IfcAlignment* pAlignment)
   {
      auto axis = pAlignment->Axis()->as<Schema::IfcAlignmentCurve>();
      auto horizontal = axis->Horizontal();
      auto segments = horizontal->Segments();
      auto nSegments = segments->size();
      if (nSegments == 0)
      {
         // alignment doesn't have any segments
         return false;
      }
      else if (nSegments == 1)
      {
         // our model doesn't support isolated transition segments
         // transition curves must be adjacent to circular curves
         auto segment = (*segments->begin());
         auto transition = segment->CurveGeometry()->as<Schema::IfcTransitionCurveSegment2D>();
         return (transition == nullptr ? true : false);
      }
      else if (nSegments == 2)
      {
         // our model doesn't support isolated transition segments
         // transition curves must be adjacent to circular curves
         // can't have two transitions adjacent to each other either
         auto segment1 = (*segments->begin());
         auto transition1 = segment1->CurveGeometry()->as<Schema::IfcTransitionCurveSegment2D>();

         auto segment2 = (*segments->begin() + 1);
         auto transition2 = segment2->CurveGeometry()->as<Schema::IfcTransitionCurveSegment2D>();
         return (transition1 == nullptr || transition2 == nullptr ? true : false);
      }
      else
      {
         // walk the segments - if we run into a transition curve we have to check the following
         // * transition is adjacent to a circular curve
         // * common radius with circular curve and transition curve are equal
         // * radius of transition curve away from circular curve is infinite
         auto begin = segments->begin();
         auto iter = begin;
         auto end = segments->end();
         for (; iter != end; iter++)
         {
            auto segment(*iter);
            auto transition = segment->CurveGeometry()->as<Schema::IfcTransitionCurveSegment2D>();
            if (transition)
            {
               Float64 start_radius = (transition->hasStartRadius() ? transition->StartRadius() : 0.0);
               Float64 end_radius = (transition->hasEndRadius() ? transition->EndRadius() : 0.0);
               if (!IsZero(start_radius) && !IsZero(end_radius))
                  return false; // one radius must be zero

               if ((iter == begin && !IsZero(start_radius)) || (iter == end - 1 && !IsZero(end_radius)))
                  return false; // if starting with a transition curve, start radius must be zero or if ending with a transition curve, end radius must be zero

               if (iter != begin && !IsZero(start_radius))
               {
                  // transition starts with a radius so a circular curve must preceed this transition curve
                  auto arc = (*(iter - 1))->CurveGeometry()->as<Schema::IfcCircularArcSegment2D>();
                  if (!arc) return false; // previous is not an arc

                  if (!IsEqual(arc->Radius(), start_radius)) return false; // common radii must be equal

                  if (arc->IsCCW() != transition->IsStartRadiusCCW()) return false; // curves must be same direction
               }

               if (iter != end - 1 && !IsZero(end_radius))
               {
                  // transition ends with a radius so a circular curve most come after this transition curve
                  auto arc = (*(iter + 1))->CurveGeometry()->as<Schema::IfcCircularArcSegment2D>();
                  if (!arc) return false; // previous is not an arc

                  if (!IsEqual(arc->Radius(), end_radius)) return false; // common radii must be equal

                  if (arc->IsCCW() != transition->IsEndRadiusCCW()) return false; // curves must be same direction
               }
            }
         }
      }


      return true;
   }
};

