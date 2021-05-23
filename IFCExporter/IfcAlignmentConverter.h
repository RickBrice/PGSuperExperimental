#pragma once

#include "F:/IfcOpenShell/_installed-vs2015-x64/include/ifcparse/IfcHierarchyHelper.h"

#if defined USE_IFC4X1
#define IfcSchema Ifc4x1
#include "F:/IfcOpenShell/_installed-vs2015-x64/include/ifcparse/Ifc4x1.h"
#elif defined USE_IFC4X3_RC2
#define IfcSchema Ifc4x3_rc2
#include "F:/IfcOpenShell/_installed-vs2015-x64/include/ifcparse/Ifc4x3_rc2.h"
#endif

#include <IFace\Project.h>

///////////////////////////////////////////////////////////////////////////
// CIfcAlignmentConverter
//
// Converts data between IFC and PGSuper data structures
class CIfcAlignmentConverter
{
public:
   CIfcAlignmentConverter(void);
   ~CIfcAlignmentConverter(void);

   // Converts PGSuper data to Ifc
   void ConvertToIfc(IBroker* pBroker, const CString& strFilePath);

   // Converts Ifc data to PGSuper
   //void ConvertToPGSuper(LX::LandXML* pLandXml, LX::Alignment* pAlignment, LX::ProfAlign* pProfAlign, LX::CrossSects* pCrossSects, LX::String& strSurfaceName, AlignmentData2* pAlignmentData, ProfileData2* pProfileData, RoadwaySectionData* pRoadwaySectionData);
   void ConvertToPGSuper(Ifc4x1::IfcAlignment* pAlignment, AlignmentData2* pAlignmentData, ProfileData2* pProfileData, RoadwaySectionData* pRoadwaySectionData);

   void SetLengthUnit(const unitLength* pUnitLength) { m_pLengthUnit = pUnitLength; }

   Ifc4x1::IfcAlignment* GetAlignment(IfcParse::IfcFile& file);

   // Returns a list of notes that were generated during the LandXML to PGSuper conversion process
   std::vector<std::_tstring> GetNotes();

private:
   const unitLength* m_pLengthUnit;
   std::vector<std::_tstring> m_Notes;
   CComPtr<ICogoEngine> m_CogoEngine;
   CComPtr<IGeomUtil2d> m_GeomUtil;

   AlignmentData2 m_AlignmentData;
   ProfileData2   m_ProfileData;
   RoadwaySectionData m_RoadwaySectionData;

   bool m_bAlignmentStarted;
   int m_ProfileState; // -1 = not yet started, 0 = started, but grade not determined, 1 = first point established

   //LX::Alignment*  CreateAlignment(IBroker* pBroker, LX::IFactory* pFactory);
   //LX::Profile*    CreateProfile(IBroker* pBroker, LX::IFactory* pFactory);
   //LX::CrossSects* CreateCrossSections(IBroker* pBroker, LX::IFactory* pFactory);
   //LX::Roadway*    CreateRoadway(IBroker* pBroker, LX::IFactory* pFactory);

   void SetupUnits(Ifc4x1::IfcAlignment* pAlignment);
   void LoadAlignment(Ifc4x1::IfcAlignment* pAlignment);
   void LoadProfile(Ifc4x1::IfcAlignment* pAlignment);
   //void LoadCrossSections(LX::CrossSects* pCrossSects, LX::String& strSurfaceName);

   void GetPoint(Ifc4x1::IfcCartesianPoint* pPoint, Float64* pX, Float64* pY);

   enum LastAlignmentType { Unknown, Line, Curve } m_LastAlignmentType;

   //// adds a line to the alignment. returns the station at the end of the line
   Float64 OnLine(Float64 startStation, Ifc4x1::IfcLineSegment2D* pLine);

   //// adds a curve to the alignment. returns the station at the end of the curve
   Float64 OnCurve(Float64 startStation, Ifc4x1::IfcTransitionCurveSegment2D* pEntrySpiral, Ifc4x1::IfcCircularArcSegment2D* pCurve, Ifc4x1::IfcTransitionCurveSegment2D* pExitSpiral);

   //// adds a PVI to the proflie. returns the elevation at the PVI
   //std::pair<Float64, Float64> OnPVI(std::pair<Float64, Float64> prevPI, const LX::PVI* pPVI);
   void LinearSegment(Ifc4x1::IfcAlignment2DVerSegLine* pLinearSegment);

   //// adds a parabolic curve to the proflie. returns the elevation at the PVI
   //std::pair<Float64, Float64> OnParaCurve(std::pair<Float64, Float64> prevPI, const LX::ParaCurve* pParaCurve);
   void ParabolicSegment(Ifc4x1::IfcAlignment2DVerSegParabolicArc* pParaCurve);

   //// adds a unsymmetrical parabolic curve to the proflie. returns the elevation at the PVI
   //std::pair<Float64, Float64> OnUnsymParaCurve(std::pair<Float64, Float64> prevPI, const LX::UnsymParaCurve* pUnsymParaCurve);

   void GetCurvePoints(Ifc4x1::IfcCircularArcSegment2D* pCurve, IPoint2d** ppStart, IPoint2d** ppPI, IPoint2d** ppEnd, IPoint2d** ppCenter);
   void GetSpiralPoints(Ifc4x1::IfcTransitionCurveSegment2D* pSpiral, IPoint2d** ppStart, IPoint2d** ppPI, IPoint2d** ppEnd);

   void CheckSpiralType(Ifc4x1::IfcTransitionCurveSegment2D* pSpiral);

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

   bool IsValidAlignment(Ifc4x1::IfcAlignment* pAlignment);
};

