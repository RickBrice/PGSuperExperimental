#include "stdafx.h"

#include <MfcTools\Prompts.h>

#include <IFace\Project.h>
#include <IFace\VersionInfo.h>
#include <IFace\Alignment.h>

// This is defined as a compiler setting
//#define USE_IFC4X1
//#define USE_IFC4X3_RC2

#include "IfcAlignmentConverter.h"
#include "IfcAlignmentConverterException.h"


#if defined _DEBUG
#pragma comment(lib,"F:/IfcOpenShell/_build-vs2015-x64/Debug/IfcParse.lib")
#else
#pragma comment(lib,"F:/IfcOpenShell/_build-vs2015-x64/Release/IfcParse.lib")
#endif


#if defined USE_IFC4X1
IfcSchema::IfcAlignment2DHorizontal* BuildAlignment(IBroker* pBroker, IfcHierarchyHelper<IfcSchema>& file);
IfcSchema::IfcAlignment2DVertical* BuildProfile(IBroker* pBroker, IfcHierarchyHelper<IfcSchema>& file);
#elif defined USE_IFC4X3_RC2
IfcSchema::IfcAlignmentHorizontal* BuildAlignment(IBroker* pBroker, IfcHierarchyHelper<IfcSchema>& file);
IfcSchema::IfcAlignmentVertical* BuildProfile(IBroker* pBroker, IfcHierarchyHelper<IfcSchema>& file);
#endif


#define ALIGNMENT_NAME _T("Alignment1")

// Constants for tracking the state of converting the profile data
#define PROFILE_NOT_STARTED -1
#define PROFILE_ESTABLISHED  1

HRESULT SameLocation(IPoint2d* pnt1, IPoint2d* pnt2)
{
   Float64 x1, y1, x2, y2;
   pnt1->get_X(&x1);
   pnt1->get_Y(&y1);
   pnt2->get_X(&x2);
   pnt2->get_Y(&y2);

   return IsEqual(x1, x2) && IsEqual(y1, y2) ? S_OK : S_FALSE;
}

// Use this throw macro when the data conversion cannot continue
// The catcher, or other, is responsible for deleting it
#define IFC_THROW(_s_) throw new CIfcAlignmentConverterException(_s_);

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

void CIfcAlignmentConverter::ConvertToIfc(IBroker* pBroker, const CString& strFilePath)
{
   USES_CONVERSION;

   IfcHierarchyHelper<IfcSchema> file;
//   file.header().file_name().name(T2A(dlg.GetFileName().GetBuffer()));
   file.header().file_name().name(T2A(strFilePath));

   //auto project = file.addProject(); // Don't like the default units in IfcOpenShell so we have do build our own
   /////////////////////////// The following is copied from addProject and tweeked
   IfcEntityList::ptr units(new IfcEntityList);
   IfcSchema::IfcDimensionalExponents* dimexp = new IfcSchema::IfcDimensionalExponents(0, 0, 0, 0, 0, 0, 0);
   IfcSchema::IfcSIUnit* unit1 = new IfcSchema::IfcSIUnit(IfcSchema::IfcUnitEnum::IfcUnit_LENGTHUNIT, boost::none, IfcSchema::IfcSIUnitName::IfcSIUnitName_METRE);
   IfcSchema::IfcSIUnit* unit2 = new IfcSchema::IfcSIUnit(IfcSchema::IfcUnitEnum::IfcUnit_PLANEANGLEUNIT, boost::none, IfcSchema::IfcSIUnitName::IfcSIUnitName_RADIAN);

   units->push(unit1);
   units->push(unit2);

   IfcSchema::IfcUnitAssignment* unit_assignment = new IfcSchema::IfcUnitAssignment(units);

   IfcSchema::IfcRepresentationContext::list::ptr rep_contexts(new IfcSchema::IfcRepresentationContext::list);
   IfcSchema::IfcProject* project = new IfcSchema::IfcProject(IfcParse::IfcGlobalId(), nullptr, boost::none, boost::none, boost::none, boost::none, boost::none, rep_contexts, unit_assignment);

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

   auto owner_history = file.getSingle<IfcSchema::IfcOwnerHistory>();
   GET_IFACE2(pBroker, IVersionInfo, pVersionInfo);
   owner_history->OwningApplication()->setApplicationFullName("BridgeLink::PGSuper");
   owner_history->OwningApplication()->setApplicationIdentifier("PGSuper");
   owner_history->OwningApplication()->setVersion(T2A(pVersionInfo->GetVersion(true)));
   owner_history->OwningApplication()->ApplicationDeveloper()->setIdentification("Washington State Department of Transportation, Bridge and Structures Office");
   owner_history->OwningApplication()->ApplicationDeveloper()->setName("Richard Brice, PE");

   auto local_placement = file.addLocalPlacement();

   auto geometric_representation_context = new IfcSchema::IfcGeometricRepresentationContext(boost::none, boost::none, 3, boost::none, local_placement, nullptr);
   file.addEntity(geometric_representation_context); // ADD THE CONTEXT TO THE FILE!!!!
   auto contexts = project->RepresentationContexts(); // get the old context
   contexts->push(geometric_representation_context); // add the new context
   project->setRepresentationContexts(contexts); // set the context back into the project

   auto horizontal_alignment = BuildAlignment(pBroker, file);
   auto vertical_profile = BuildProfile(pBroker, file);

   //////////////////////////////////////////////////////////////////////////////////

   // create an alignment curve in 3D using the horizontal alignment and vertical profile (nullptr for profile since it's optional)
   // we could also use an IfcOffsetCurveByDistances or an IfcPolyline to define the reference system curve
#if defined USE_IFC4X3_RC2
#pragma Reminder("WORKING HERE - IfcAlignmentCurve is marked as depreciated so what is the alternative?")
#endif
   auto axis = new IfcSchema::IfcAlignmentCurve(horizontal_alignment, vertical_profile, boost::none);

   // create an alignment... this is the linear positioning element. it could be an alignment curve (as used above) or IfcOffsetCurveByDistances or IfcPolyline
   auto alignment = new IfcSchema::IfcAlignment(IfcParse::IfcGlobalId(), owner_history, std::string("Alignment"), boost::none,
      boost::none /*object type*/, local_placement /*object placement*/, file.getSingle<IfcSchema::IfcProductRepresentation>()/*representation*/, axis, boost::none /*IfcSchema::IfcAlignmentTypeEnum::IfcAlignmentType_USERDEFINED*/);
   file.addEntity(alignment);

   boost::shared_ptr<IfcTemplatedEntityList<IfcSchema::IfcProduct>> related_elements(new IfcTemplatedEntityList<IfcSchema::IfcProduct>());
   related_elements->push(alignment);

   auto spatial_structure = new IfcSchema::IfcRelContainedInSpatialStructure(IfcParse::IfcGlobalId(), owner_history, boost::none, boost::none, related_elements, site);
   file.addEntity(spatial_structure);

   site->ContainsElements()->push(spatial_structure);

   std::ofstream ofs(T2A(strFilePath));
   ofs << file;
}


//void CIfcAlignmentConverter::ConvertToLandXML(IBroker* pBroker, LX::Document* pDocument)
//{
//   LX::IFactory* pFactory = pDocument->getFactory();
//
//   // Create LandXML node and set its attributes.
//   LX::LandXML* pLandXml = pFactory->createLandXML();
//   pLandXml->setLanguage(_T("English"));
//   pLandXml->setVersion(_T("1.2")); // LandXML version
//
//   CTime now = CTime::GetCurrentTime();
//   CString strTime = now.Format(_T("%H:%M:%S"));
//   CString strDate = now.Format(_T("%Y-%m-%d"));
//
//   pLandXml->setTime(strTime.GetBuffer());
//   pLandXml->setDate(strDate.GetBuffer());
//
//
//   LX::Application* pApplication = pFactory->createApplication();
//   pApplication->setName(_T("PGSuper Extension Pack"));
//   pApplication->setManufacturer(_T("BridgeSight Software"));
//   pApplication->setManufacturerURL(_T("http://www.bridgesight.com"));
//
//
//   GET_IFACE2(pBroker, IProjectProperties, pProjectProperties);
//   LX::Author* pAuthor = pFactory->createAuthor();
//   pAuthor->setCompany(pProjectProperties->GetCompany());
//   pAuthor->setCreatedBy(pProjectProperties->GetEngineer());
//   pApplication->Author().addItem(pAuthor);
//
//   pLandXml->setApplication(pApplication);
//
//   // Create Project tag.
//   LX::Project* pProject = pFactory->createProject();
//   LX::String strName(pProjectProperties->GetBridgeName());
//   pProject->setName(strName);
//   pLandXml->setProject(pProject);
//
//   // Create Units tag
//   LX::Units* pUnits = pFactory->createUnits();
//
//   // Use Imperial Units
//   LX::Imperial* pImperial = pFactory->createImperial();
//
//   pImperial->setAngularUnit(LX::EnumAngularType::Values::k_decimal_degrees);
//   pImperial->setAreaUnit(LX::EnumImpArea::Values::k_squareFoot);
//   pImperial->setDiameterUnit(LX::EnumImpDiameter::Values::k_foot);
//   pImperial->setDirectionUnit(LX::EnumAngularType::Values::k_decimal_degrees);
//   pImperial->setHeightUnit(LX::EnumImpHeight::Values::k_foot);
//   pImperial->setLinearUnit(LX::EnumImpLinear::Values::k_foot);
//   pImperial->setVolumeUnit(LX::EnumImpVolume::Values::k_cubicYard);
//   pImperial->setTemperatureUnit(LX::EnumImpTemperature::k_fahrenheit);
//   pImperial->setPressureUnit(LX::EnumImpPressure::k_inHG);
//
//   // Set the Imperial units in the Units object.
//   pUnits->setSelectedUnits(pImperial);
//   pLandXml->setUnits(pUnits);
//
//   // Create Alignments
//   LX::Alignments* pAlignments = pFactory->createAlignments();
//   LX::Alignment*  pAlignment = CreateAlignment(pBroker, pFactory);
//   LX::Profile*    pProfile = CreateProfile(pBroker, pFactory);
//   //LX::CrossSects* pCrossSects = CreateCrossSections(pBroker,pFactory);
//   pAlignment->Profile().addItem(pProfile);
//   //pAlignment->CrossSects().addItem(pCrossSects);
//   pAlignments->Alignment().addItem(pAlignment);
//   pLandXml->Alignments().addItem(pAlignments);
//
//   LX::Roadways* pRoadways = pFactory->createRoadways();
//   LX::Roadway*  pRoadway = CreateRoadway(pBroker, pFactory);
//   pRoadways->Roadway().addItem(pRoadway);
//   pLandXml->Roadways().addItem(pRoadways);
//
//   // Set the root node of the document.
//   pDocument->setRootObject(pLandXml);
//}
//
//LX::Alignment* CIfcAlignmentConverter::CreateAlignment(IBroker* pBroker, LX::IFactory* pFactory)
//{
//   LX::Alignment*  pAlignment = pFactory->createAlignment();
//   pAlignment->setName(ALIGNMENT_NAME);
//   pAlignment->setDesc(_T("Alignment exported from PGSuper"));
//
//   // create geometry of the alignment
//   LX::CoordGeom* pCoordGeom = pFactory->createCoordGeom();
//   pCoordGeom->setName(pAlignment->getName());
//   pAlignment->setCoordGeom(pCoordGeom); // attach them together
//
//
//   GET_IFACE2(pBroker, IBridge, pBridge);
//   PierIndexType nPiers = pBridge->GetPierCount();
//   SpanIndexType nSpans = pBridge->GetSpanCount();
//
//   Float64 first_pier_station = pBridge->GetPierStation(0);
//   Float64 last_pier_station = pBridge->GetPierStation(nPiers - 1);
//
//   Float64 first_span_length = pBridge->GetPierStation(1) - pBridge->GetPierStation(0);
//   Float64 last_span_length = pBridge->GetPierStation(nPiers - 1) - pBridge->GetPierStation(nPiers - 2);
//
//
//   // figure out parameters for the alignment
//   GET_IFACE2(pBroker, IRoadway, pPGSuperAlignment);
//   GET_IFACE2(pBroker, IRoadwayData, pRoadwayData);
//
//   const AlignmentData2& alignmentData = pRoadwayData->GetAlignmentData2();
//
//   // get beginning of alignment station and coordinate
//   Float64 TS_station = first_pier_station - first_span_length / 2;
//   Float64 ST_station = last_pier_station + last_span_length / 2;
//
//   if (alignmentData.HorzCurves.size() != 0)
//   {
//      // alignment has curves
//      HorzCurveData hcData = alignmentData.HorzCurves.front();
//      if (!IsZero(hcData.Radius))
//      {
//         // alignment starts and ends with a real horizontal curve... get the TS and ST stations
//
//         CComPtr<IHorzCurve> first_curve;
//         pPGSuperAlignment->GetCurve(0, &first_curve);
//         Float64 tangent;
//         first_curve->get_Tangent(&tangent);
//         TS_station = hcData.PIStation - tangent;
//      }
//
//      hcData = alignmentData.HorzCurves.back();
//      if (!IsZero(hcData.Radius))
//      {
//         CComPtr<IHorzCurve> last_curve;
//         pPGSuperAlignment->GetCurve(alignmentData.HorzCurves.size() - 1, &last_curve);
//         Float64 tangent;
//         last_curve->get_Tangent(&tangent);
//         Float64 length;
//         last_curve->get_TotalLength(&length);
//         ST_station = hcData.PIStation - tangent + length;
//      }
//   }
//
//   Float64 start_station = Min(alignmentData.RefStation, first_pier_station - first_span_length / 2, TS_station);
//   Float64 start_x, start_y;
//
//   CComPtr<IPoint2d> pntStart;
//   pPGSuperAlignment->GetPoint(start_station, 0.0, nullptr, pgsTypes::pcGlobal, &pntStart);
//   pntStart->get_X(&start_x);
//   pntStart->get_Y(&start_y);
//
//   Float64 end_station = Max(last_pier_station + last_span_length / 2, ST_station);
//   Float64 end_x, end_y;
//   CComPtr<IPoint2d> pntEnd;
//   pPGSuperAlignment->GetPoint(end_station, 0.0, nullptr, pgsTypes::pcGlobal, &pntEnd);
//   pntEnd->get_X(&end_x);
//   pntEnd->get_Y(&end_y);
//
//   start_x = ::ConvertFromSysUnits(start_x, unitMeasure::Feet);
//   start_y = ::ConvertFromSysUnits(start_y, unitMeasure::Feet);
//   end_x = ::ConvertFromSysUnits(end_x, unitMeasure::Feet);
//   end_y = ::ConvertFromSysUnits(end_y, unitMeasure::Feet);
//
//   // lay out alignment.... start point to end point if no curves
//   // else start point, line to first ST, CURVE, line from TS to ST, Curve, etc... finish with line from ST to end
//   // Check for zero length lines and skip them
//   pAlignment->setStaStart(::ConvertFromSysUnits(start_station, unitMeasure::Feet));
//   pAlignment->setLength(::ConvertFromSysUnits(end_station - start_station, unitMeasure::Feet));
//   if (alignmentData.HorzCurves.size() == 0)
//   {
//      // Use a straight line to represent the alignment
//
//      // Create line.
//      LX::Line* pLine = pFactory->createLine();
//      LX::Start* pStart = nullptr;
//      LX::End*   pEnd = nullptr;
//
//      pStart = pFactory->createStart();
//      pStart->addItem(start_y); // northing
//      pStart->addItem(start_x); // easting
//
//      pEnd = pFactory->createEnd();
//      pEnd->addItem(end_y); // northing
//      pEnd->addItem(end_x); // easting
//
//      pLine->setStart(pStart);
//      pLine->setEnd(pEnd);
//
//      pCoordGeom->GeomList().addItem(pLine);
//   }
//   else
//   {
//      // use a series of lines and curves to represent the alignment
//      Float64 prev_x = start_x;
//      Float64 prev_y = start_y;
//      auto begin = std::cbegin(alignmentData.HorzCurves);
//      auto iter = begin;
//      auto end = std::cend(alignmentData.HorzCurves);
//      for (; iter != end; iter++)
//      {
//         const auto& hcData = *iter;
//         if (IsZero(hcData.Radius))
//         {
//            CComPtr<IPoint2d> pntPI;
//            pPGSuperAlignment->GetPoint(hcData.PIStation, 0.00, nullptr, pgsTypes::pcGlobal, &pntPI);
//
//            LX::Start* pStart = nullptr;
//            LX::POI*   pPI = nullptr;
//            LX::End*   pEnd = nullptr;
//            if (iter == begin)
//            {
//               // need a point on back tangent
//               ATLASSERT(start_station < hcData.PIStation);
//               CComPtr<IPoint2d> pntPBT;
//               pPGSuperAlignment->GetPoint(start_station, 0.00, nullptr, pgsTypes::pcGlobal, &pntPBT);
//               Float64 x, y;
//               pntPBT->get_X(&x);
//               pntPBT->get_Y(&y);
//               x = ::ConvertFromSysUnits(x, unitMeasure::Feet);
//               y = ::ConvertFromSysUnits(y, unitMeasure::Feet);
//               prev_x = x;
//               prev_y = y;
//               pStart = pFactory->createStart();
//               pStart->addItem(y);
//               pStart->addItem(x);
//            }
//
//            // PI point
//            pPI = pFactory->createPOI();
//            Float64 x, y;
//            pntPI->get_X(&x);
//            pntPI->get_Y(&y);
//            x = ::ConvertFromSysUnits(x, unitMeasure::Feet);
//            y = ::ConvertFromSysUnits(y, unitMeasure::Feet);
//            prev_x = x;
//            prev_y = y;
//            pPI->addItem(y); // Northing
//            pPI->addItem(x); // Easting
//
//            if (iter == alignmentData.HorzCurves.end() - 1)
//            {
//               // need a point on forward tangent
//               ATLASSERT(hcData.PIStation < end_station);
//               CComPtr<IPoint2d> pntPFT;
//               pPGSuperAlignment->GetPoint(end_station, 0.00, nullptr, pgsTypes::pcGlobal, &pntPFT);
//               Float64 x, y;
//               pntPFT->get_X(&x);
//               pntPFT->get_Y(&y);
//               x = ::ConvertFromSysUnits(x, unitMeasure::Feet);
//               y = ::ConvertFromSysUnits(y, unitMeasure::Feet);
//               prev_x = x;
//               prev_y = y;
//               pEnd = pFactory->createEnd();
//               pEnd->addItem(y);
//               pEnd->addItem(x);
//            }
//
//            if (pStart == nullptr && pEnd == nullptr)
//            {
//               pCoordGeom->GeomList().addItem(pPI);
//            }
//
//            if (pStart)
//            {
//               LX::Line* pLine = pFactory->createLine();
//               pLine->setStart(pStart);
//
//               LX::End* pPIEnd = pFactory->createEnd();
//               pPIEnd->addItem(pPI->at(0));
//               pPIEnd->addItem(pPI->at(1));
//               pLine->setEnd(pPIEnd);
//
//               pCoordGeom->GeomList().addItem(pLine);
//            }
//
//            if (pEnd)
//            {
//               LX::Line* pLine = pFactory->createLine();
//
//               LX::Start* pPIStart = pFactory->createStart();
//               pPIStart->addItem(pPI->at(0));
//               pPIStart->addItem(pPI->at(1));
//               pLine->setStart(pPIStart);
//               pLine->setEnd(pEnd);
//
//               pCoordGeom->GeomList().addItem(pLine);
//            }
//         }
//         else
//         {
//            // horizontal curve
//            CComPtr<IHorzCurve> hCurve;
//            CollectionIndexType curveIdx = CollectionIndexType(iter - alignmentData.HorzCurves.begin());
//            pPGSuperAlignment->GetCurve(curveIdx, &hCurve);
//            CComPtr<IPoint2d> pntTS, pntSC, pntPI, pntCS, pntST, pntCC;
//            pPGSuperAlignment->GetCurvePoint(curveIdx, cptTS, pgsTypes::pcGlobal, &pntTS);
//            pPGSuperAlignment->GetCurvePoint(curveIdx, cptSC, pgsTypes::pcGlobal, &pntSC);
//            pPGSuperAlignment->GetCurvePoint(curveIdx, cptPI, pgsTypes::pcGlobal, &pntPI);
//            pPGSuperAlignment->GetCurvePoint(curveIdx, cptCS, pgsTypes::pcGlobal, &pntCS);
//            pPGSuperAlignment->GetCurvePoint(curveIdx, cptST, pgsTypes::pcGlobal, &pntST);
//            pPGSuperAlignment->GetCurvePoint(curveIdx, cptCC, pgsTypes::pcGlobal, &pntCC);
//
//            Float64 radius;
//            hCurve->get_Radius(&radius);
//
//            Float64 Lin, Lout;
//            hCurve->get_SpiralLength(spEntry, &Lin);
//            hCurve->get_SpiralLength(spExit, &Lout);
//
//            CurveDirectionType dirType;
//            hCurve->get_Direction(&dirType);
//
//            Float64 x, y;
//            pntTS->get_X(&x);
//            pntTS->get_Y(&y);
//            x = ::ConvertFromSysUnits(x, unitMeasure::Feet);
//            y = ::ConvertFromSysUnits(y, unitMeasure::Feet);
//
//            // line before curve
//            Float64 dx = x - prev_x;
//            Float64 dy = y - prev_y;
//            Float64 length = sqrt(dx*dx + dy*dy);
//            if (!IsZero(length))
//            {
//               LX::Line* pLine = pFactory->createLine();
//               LX::Start* pStart = nullptr;
//               LX::End*   pEnd = nullptr;
//
//               pStart = pFactory->createStart();
//               pStart->addItem(prev_y); // northing
//               pStart->addItem(prev_x); // easting
//
//               pEnd = pFactory->createEnd();
//               pEnd->addItem(y); // northing
//               pEnd->addItem(x); // easting
//
//               pLine->setStart(pStart);
//               pLine->setEnd(pEnd);
//
//               pCoordGeom->GeomList().addItem(pLine);
//            }
//
//            // Entry spiral
//            if (!IsZero(Lin))
//            {
//               LX::Spiral* pSpiral = pFactory->createSpiral();
//               pSpiral->setSpiType(LX::EnumSpiralType::k_clothoid); // only kind of spiral PGSuper supports
//               pSpiral->setRot(dirType == cdLeft ? LX::EnumClockwise::k_ccw : LX::EnumClockwise::k_cw);
//
//               pSpiral->setLength(::ConvertFromSysUnits(Lin, unitMeasure::Feet));
//               pSpiral->setRadiusStart(std::numeric_limits<Float64>::infinity());
//               pSpiral->setRadiusEnd(::ConvertFromSysUnits(radius, unitMeasure::Feet));
//
//               CComPtr<IPoint2d> pntSPI;
//               hCurve->get_SPI(spEntry, &pntSPI);
//
//               pntTS->get_X(&x);
//               pntTS->get_Y(&y);
//               LX::Start* pStart = pFactory->createStart();
//               x = ::ConvertFromSysUnits(x, unitMeasure::Feet);
//               y = ::ConvertFromSysUnits(y, unitMeasure::Feet);
//               pStart->addItem(y); // northing
//               pStart->addItem(x); // easting
//               pSpiral->setStart(pStart);
//
//               pntSPI->get_X(&x);
//               pntSPI->get_Y(&y);
//               LX::POI* pPI = pFactory->createPOI();
//               x = ::ConvertFromSysUnits(x, unitMeasure::Feet);
//               y = ::ConvertFromSysUnits(y, unitMeasure::Feet);
//               pPI->addItem(y); // northing
//               pPI->addItem(x); // easting
//               pSpiral->setPI(pPI);
//
//               pntSC->get_X(&x);
//               pntSC->get_Y(&y);
//               LX::End* pEnd = pFactory->createEnd();
//               x = ::ConvertFromSysUnits(x, unitMeasure::Feet);
//               y = ::ConvertFromSysUnits(y, unitMeasure::Feet);
//               pEnd->addItem(y); // northing
//               pEnd->addItem(x); // easting
//               pSpiral->setEnd(pEnd);
//
//               pCoordGeom->GeomList().addItem(pSpiral);
//            }
//
//            // Horizontal curve
//            LX::Curve* pCurve = pFactory->createCurve();
//            pCurve->setRot(dirType == cdLeft ? LX::EnumClockwise::k_ccw : LX::EnumClockwise::k_cw);
//
//            LX::Start* pStart = pFactory->createStart();
//            pntSC->get_X(&x);
//            pntSC->get_Y(&y);
//            x = ::ConvertFromSysUnits(x, unitMeasure::Feet);
//            y = ::ConvertFromSysUnits(y, unitMeasure::Feet);
//            pStart->addItem(y);
//            pStart->addItem(x);
//            pCurve->setStart(pStart);
//
//            LX::POI* pPI = pFactory->createPOI();
//            pntPI->get_X(&x);
//            pntPI->get_Y(&y);
//            x = ::ConvertFromSysUnits(x, unitMeasure::Feet);
//            y = ::ConvertFromSysUnits(y, unitMeasure::Feet);
//            pPI->addItem(y);
//            pPI->addItem(x);
//            pCurve->setPI(pPI);
//
//            LX::Center* pCenter = pFactory->createCenter();
//            pntCC->get_X(&x);
//            pntCC->get_Y(&y);
//            x = ::ConvertFromSysUnits(x, unitMeasure::Feet);
//            y = ::ConvertFromSysUnits(y, unitMeasure::Feet);
//            pCenter->addItem(y);
//            pCenter->addItem(x);
//            pCurve->setCenter(pCenter);
//
//            LX::End* pEnd = pFactory->createEnd();
//            pntCS->get_X(&x);
//            pntCS->get_Y(&y);
//            x = ::ConvertFromSysUnits(x, unitMeasure::Feet);
//            y = ::ConvertFromSysUnits(y, unitMeasure::Feet);
//            pEnd->addItem(y);
//            pEnd->addItem(x);
//            pCurve->setEnd(pEnd);
//            prev_x = x;
//            prev_y = y;
//
//            pCoordGeom->GeomList().addItem(pCurve);
//
//            // Exit spiral
//            if (!IsZero(Lout))
//            {
//               LX::Spiral* pSpiral = pFactory->createSpiral();
//               pSpiral->setSpiType(LX::EnumSpiralType::k_clothoid); // only kind of spiral PGSuper supports
//               pSpiral->setRot(dirType == cdLeft ? LX::EnumClockwise::k_ccw : LX::EnumClockwise::k_cw);
//
//               pSpiral->setLength(::ConvertFromSysUnits(Lout, unitMeasure::Feet));
//               pSpiral->setRadiusStart(::ConvertFromSysUnits(radius, unitMeasure::Feet));
//               pSpiral->setRadiusEnd(std::numeric_limits<Float64>::infinity());
//
//               CComPtr<IPoint2d> pntSPI;
//               hCurve->get_SPI(spExit, &pntSPI);
//
//               pntCS->get_X(&x);
//               pntCS->get_Y(&y);
//               LX::Start* pStart = pFactory->createStart();
//               x = ::ConvertFromSysUnits(x, unitMeasure::Feet);
//               y = ::ConvertFromSysUnits(y, unitMeasure::Feet);
//               pStart->addItem(y); // northing
//               pStart->addItem(x); // easting
//               pSpiral->setStart(pStart);
//
//               pntSPI->get_X(&x);
//               pntSPI->get_Y(&y);
//               LX::POI* pPI = pFactory->createPOI();
//               x = ::ConvertFromSysUnits(x, unitMeasure::Feet);
//               y = ::ConvertFromSysUnits(y, unitMeasure::Feet);
//               pPI->addItem(y); // northing
//               pPI->addItem(x); // easting
//               pSpiral->setPI(pPI);
//
//               pntST->get_X(&x);
//               pntST->get_Y(&y);
//               LX::End* pEnd = pFactory->createEnd();
//               x = ::ConvertFromSysUnits(x, unitMeasure::Feet);
//               y = ::ConvertFromSysUnits(y, unitMeasure::Feet);
//               pEnd->addItem(y); // northing
//               pEnd->addItem(x); // easting
//               pSpiral->setEnd(pEnd);
//
//               pCoordGeom->GeomList().addItem(pSpiral);
//
//               prev_x = x;
//               prev_y = y;
//            }
//         }
//      }
//
//      // last line
//      Float64 dx = end_x - prev_x;
//      Float64 dy = end_y - prev_y;
//      Float64 length = sqrt(dx*dx + dy*dy);
//      if (!IsZero(length))
//      {
//         LX::Line* pLine = pFactory->createLine();
//         LX::Start* pStart = nullptr;
//         LX::End*   pEnd = nullptr;
//
//         pStart = pFactory->createStart();
//         pStart->addItem(prev_y); // northing
//         pStart->addItem(prev_x); // easting
//
//         pEnd = pFactory->createEnd();
//         pEnd->addItem(end_y); // northing
//         pEnd->addItem(end_x); // easting
//
//         pLine->setStart(pStart);
//         pLine->setEnd(pEnd);
//
//         pCoordGeom->GeomList().addItem(pLine);
//      }
//   }
//
//   return pAlignment;
//}
//
//LX::Profile* CIfcAlignmentConverter::CreateProfile(IBroker* pBroker, LX::IFactory* pFactory)
//{
//   LX::Profile* pProfile = pFactory->createProfile();
//   pProfile->setDesc(_T("Profile exported from PGSuper"));
//
//   // ProfAlign
//   LX::ProfAlign* pProfAlign = pFactory->createProfAlign();
//   pProfAlign->setName(_T("Profile Grade Line"));
//
//   GET_IFACE2(pBroker, IRoadwayData, pRoadwayData);
//
//   const ProfileData2& profileData = pRoadwayData->GetProfileData2();
//   if (profileData.VertCurves.size() == 0)
//   {
//      LX::PVI* pBeginVI = pFactory->createPVI();
//      Float64 station = profileData.Station;
//      Float64 elevation = profileData.Elevation;
//      station = ::ConvertFromSysUnits(station, unitMeasure::Feet);
//      elevation = ::ConvertFromSysUnits(elevation, unitMeasure::Feet);
//
//      pBeginVI->addItem(station);
//      pBeginVI->addItem(elevation);
//
//      LX::PVI* pEndVI = pFactory->createPVI();
//      station += 100;
//      elevation += 100 * profileData.Grade;
//
//      pEndVI->addItem(station);
//      pEndVI->addItem(elevation);
//
//      pProfAlign->VertGeomList().addItem(pBeginVI);
//      pProfAlign->VertGeomList().addItem(pEndVI);
//   }
//   else
//   {
//      GET_IFACE2(pBroker, IRoadway, pPGSuperAlignment);
//
//      auto begin = std::cbegin(profileData.VertCurves);
//      auto iter = begin;
//      auto end = std::cend(profileData.VertCurves);
//      for (; iter != end; iter++)
//      {
//         const auto& vcData = *iter;
//         if (IsZero(vcData.L1))
//         {
//            // this is a PVI
//            if (iter == begin)
//            {
//               LX::PVI* pPVI = pFactory->createPVI();
//               Float64 station = (IsEqual(profileData.Station, vcData.PVIStation) ? Min(0.0, profileData.Station) : profileData.Station);
//               Float64 elevation = pPGSuperAlignment->GetElevation(station, 0.0);
//               station = ::ConvertFromSysUnits(station, unitMeasure::Feet);
//               elevation = ::ConvertFromSysUnits(elevation, unitMeasure::Feet);
//               pPVI->addItem(station);
//               pPVI->addItem(elevation);
//               pProfAlign->VertGeomList().addItem(pPVI);
//            }
//
//            LX::PVI* pPVI = pFactory->createPVI();
//            Float64 elevation = pPGSuperAlignment->GetElevation(vcData.PVIStation, 0.0);
//            elevation = ::ConvertFromSysUnits(elevation, unitMeasure::Feet);
//            pPVI->addItem(::ConvertFromSysUnits(vcData.PVIStation, unitMeasure::Feet));
//            pPVI->addItem(elevation);
//            pProfAlign->VertGeomList().addItem(pPVI);
//
//            if (iter == end - 1)
//            {
//               LX::PVI* pPVI = pFactory->createPVI();
//               Float64 elevation = pPGSuperAlignment->GetElevation(vcData.PVIStation + 100, 0.0);
//               elevation = ::ConvertFromSysUnits(elevation, unitMeasure::Feet);
//               pPVI->addItem(::ConvertFromSysUnits(vcData.PVIStation + 100.0, unitMeasure::Feet));
//               pPVI->addItem(elevation);
//               pProfAlign->VertGeomList().addItem(pPVI);
//            }
//         }
//         else
//         {
//            // this is a vertical curve
//
//            CComPtr<IVertCurve> vertCurve;
//            CollectionIndexType idx = CollectionIndexType(iter - begin);
//            pPGSuperAlignment->GetVertCurve(idx, &vertCurve);
//
//            // first VI point is the begin of the first vert curve
//            if (iter == begin)
//            {
//               CComPtr<IProfilePoint> pntBVC;
//               vertCurve->get_BVC(&pntBVC);
//               LX::PVI* pBeginVI = pFactory->createPVI();
//               CComPtr<IStation> objStation;
//               Float64 station, elevation;
//               pntBVC->get_Station(&objStation);
//               objStation->get_Value(&station);
//               pntBVC->get_Elevation(&elevation);
//               station = ::ConvertFromSysUnits(station, unitMeasure::Feet);
//               elevation = ::ConvertFromSysUnits(elevation, unitMeasure::Feet);
//
//               pBeginVI->addItem(station);
//               pBeginVI->addItem(elevation);
//               pProfAlign->VertGeomList().addItem(pBeginVI);
//            }
//
//            // do all the vert curves
//            CComPtr<IProfilePoint> pntPVI;
//            vertCurve->get_PVI(&pntPVI);
//            CComPtr<IStation> objStation;
//            Float64 station, elevation;
//            pntPVI->get_Station(&objStation);
//            objStation->get_Value(&station);
//            pntPVI->get_Elevation(&elevation);
//            station = ::ConvertFromSysUnits(station, unitMeasure::Feet);
//            elevation = ::ConvertFromSysUnits(elevation, unitMeasure::Feet);
//
//            Float64 l1, l2;
//            vertCurve->get_L1(&l1);
//            vertCurve->get_L2(&l2);
//            l1 = ::ConvertFromSysUnits(l1, unitMeasure::Feet);
//            l2 = ::ConvertFromSysUnits(l2, unitMeasure::Feet);
//
//            if (IsEqual(l1, l2))
//            {
//               // symetrical VC
//               LX::ParaCurve* pParaCurve = pFactory->createParaCurve();
//               pParaCurve->addItem(station);
//               pParaCurve->addItem(elevation);
//               pParaCurve->setLength(l1 + l2);
//               pProfAlign->VertGeomList().addItem(pParaCurve);
//            }
//            else
//            {
//               // unsymetrical VC
//               LX::UnsymParaCurve* pParaCurve = pFactory->createUnsymParaCurve();
//               pParaCurve->addItem(station);
//               pParaCurve->addItem(elevation);
//               pParaCurve->setLengthIn(l1);
//               pParaCurve->setLengthOut(l2);
//               pProfAlign->VertGeomList().addItem(pParaCurve);
//            }
//
//            // end with the EVC of the last vert curve
//            if (iter == end - 1)
//            {
//               CComPtr<IProfilePoint> pntEVC;
//               vertCurve->get_EVC(&pntEVC);
//               LX::PVI* pEndVI = pFactory->createPVI();
//               CComPtr<IStation> objStation;
//               Float64 station, elevation;
//               pntEVC->get_Station(&objStation);
//               objStation->get_Value(&station);
//               pntEVC->get_Elevation(&elevation);
//               station = ::ConvertFromSysUnits(station, unitMeasure::Feet);
//               elevation = ::ConvertFromSysUnits(elevation, unitMeasure::Feet);
//
//               pEndVI->addItem(station);
//               pEndVI->addItem(elevation);
//               pProfAlign->VertGeomList().addItem(pEndVI);
//            }
//         }
//      }
//   }
//
//   pProfile->ProfAlign().addItem(pProfAlign);
//   return pProfile;
//}
//
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

void CIfcAlignmentConverter::ConvertToPGSuper(Ifc4x1::IfcAlignment* pAlignment, AlignmentData2* pAlignmentData, ProfileData2* pProfileData, RoadwaySectionData* pRoadwaySectionData)
{
   SetupUnits(pAlignment);
   LoadAlignment(pAlignment);
   *pAlignmentData = m_AlignmentData;

   LoadProfile(pAlignment);
   *pProfileData = m_ProfileData;

#pragma Reminder("WORKING HERE - Roadway Section Data")
}

std::vector<std::_tstring> CIfcAlignmentConverter::GetNotes()
{
   return m_Notes;
}

void CIfcAlignmentConverter::SetupUnits(Ifc4x1::IfcAlignment* pAlignment)
{
#pragma Reminder("WORKING HERE - Units")
   // need to work back from pAlignment to the project to get the unit of measure in the file
   // for now, we will assume meter because that is what PGSuper writes and it is what is used in the IFC example file
   m_pLengthUnit = &unitMeasure::Meter;

//   LX::Units* pUnits = pLandXml->getUnits();
//   LX::Object* pSelectedUnits = pUnits->getSelectedUnits();
//   LX::Metric* pMetric = dynamic_cast<LX::Metric*>(pSelectedUnits);
//   LX::Imperial* pImperial = dynamic_cast<LX::Imperial*>(pSelectedUnits);
//   if (pMetric)
//   {
//      switch (pMetric->getLinearUnit())
//      {
//      case LX::EnumMetLinear::k_centimeter:
//         m_pLengthUnit = &unitMeasure::Centimeter;
//         break;
//
//      case LX::EnumMetLinear::k_kilometer:
//         m_pLengthUnit = &unitMeasure::Kilometer;
//         break;
//
//      case LX::EnumMetLinear::k_meter:
//         m_pLengthUnit = &unitMeasure::Meter;
//         break;
//
//      case LX::EnumMetLinear::k_millimeter:
//         m_pLengthUnit = &unitMeasure::Millimeter;
//         break;
//
//      default:
//         XML_THROW(_T("An unknown unit of measure of encountered"));
//      }
//   }
//   else if (pImperial)
//   {
//      switch (pImperial->getLinearUnit())
//      {
//      case LX::EnumImpLinear::k_foot:
//         m_pLengthUnit = &unitMeasure::Feet;
//         break;
//
//      case LX::EnumImpLinear::k_USSurveyFoot:
//         m_pLengthUnit = &unitMeasure::USSurveyFoot;
//         break;
//
//      case LX::EnumImpLinear::k_inch:
//         m_pLengthUnit = &unitMeasure::Inch;
//         break;
//
//      case LX::EnumImpLinear::k_mile:
//         m_pLengthUnit = &unitMeasure::Mile;
//         break;
//
//      default:
//         XML_THROW(_T("An unknown unit of measure of encountered"));
//      }
//   }
//   else
//   {
//      XML_THROW(_T("An unknown unit type of encountered"));
//   }
}

void CIfcAlignmentConverter::LoadAlignment(Ifc4x1::IfcAlignment* pAlignment)
{
   m_bAlignmentStarted = false; // the alignment datablock has not yet been started

   // initalize the alignment data
   m_AlignmentData.Direction = 0.00;
   m_AlignmentData.xRefPoint = 0.00;
   m_AlignmentData.yRefPoint = 0.00;
   m_AlignmentData.HorzCurves.clear();

   auto axis = pAlignment->Axis();
   auto curve = axis->as<Ifc4x1::IfcAlignmentCurve>();
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

      auto linear = segment->CurveGeometry()->as<Ifc4x1::IfcLineSegment2D>();
      auto transition = segment->CurveGeometry()->as<Ifc4x1::IfcTransitionCurveSegment2D>();
      auto curve = segment->CurveGeometry()->as<Ifc4x1::IfcCircularArcSegment2D>();

      Ifc4x1::IfcTransitionCurveSegment2D* entrySpiral = nullptr;
      Ifc4x1::IfcTransitionCurveSegment2D* exitSpiral = nullptr;

      Float64 end_station = current_station;

      if (linear)
      {
         end_station = OnLine(current_station, linear);
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
            curve = (*iter)->CurveGeometry()->as<Ifc4x1::IfcCircularArcSegment2D>();
            if (curve)
            {
               // it's a curve... is there an element that follows the curve?
               bIsThereANextSegment = ((iter+1) != end);
               if (bIsThereANextSegment)
               {
                  // if there is a next segment, see if it is a spiral
                  exitSpiral = (*(iter+1))->CurveGeometry()->as<Ifc4x1::IfcTransitionCurveSegment2D>();

                  // if not a spiral, pExitSpiral will be nullptr
                  // this is OK, it just means we have a Spiral-Curve situation
                  if (exitSpiral)
                  {
                     // it is an exit spiral, so advance the iterator
                     iter++;
                     bIsThereANextSegment = ((iter+1) != end);
                  }
               }

               end_station = OnCurve(current_station, entrySpiral, curve, exitSpiral);
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
            exitSpiral = (*(iter + 1))->CurveGeometry()->as<Ifc4x1::IfcTransitionCurveSegment2D>();
            
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
               auto next_curve = (*(iter))->CurveGeometry()->as<Ifc4x1::IfcCircularArcSegment2D>();
               if (next_curve && exitSpiral)
               {
                  CComPtr<IPoint2d> pntSpiralStart, pntSpiralPI, pntSpiralEnd;
                  GetSpiralPoints(exitSpiral, &pntSpiralStart, &pntSpiralPI, &pntSpiralEnd);

                  CComPtr<IPoint2d> pntNextCurveStart, pntNextCurvePI, pntNextCurveEnd, pntNextCurveCenter;
                  GetCurvePoints(next_curve, &pntNextCurveStart, &pntNextCurvePI, &pntNextCurveEnd, &pntNextCurveCenter);

                  if (SameLocation(pntSpiralEnd, pntNextCurveStart) == S_FALSE)
                  {
                     IFC_THROW(_T("PGSuper cannot model a transition spiral between two circular curves."));
                  }
               }
            }
         }

         ATLASSERT(entrySpiral == nullptr); // this must be the case
         end_station = OnCurve(current_station, entrySpiral, curve, exitSpiral);
      }
      else
      {
         ATLASSERT(false);
      }
      current_station = end_station;
   }
}

void CIfcAlignmentConverter::LoadProfile(Ifc4x1::IfcAlignment* pAlignment)
{
   m_ProfileState = PROFILE_NOT_STARTED;
   m_ProfileData.VertCurves.clear();

   auto axis = pAlignment->Axis();
   auto curve = axis->as<Ifc4x1::IfcAlignmentCurve>();
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
      auto linear_segment = segment->as<Ifc4x1::IfcAlignment2DVerSegLine>();
      auto parabolic_arc = segment->as<Ifc4x1::IfcAlignment2DVerSegParabolicArc>();
      auto circular_arc = segment->as<Ifc4x1::IfcAlignment2DVerSegCircularArc>();

      if (linear_segment)
      {
         LinearSegment(linear_segment);
      }
      else if (parabolic_arc)
      {
         ParabolicSegment(parabolic_arc);
      }
      else if (circular_arc)
      {
#pragma Reminder("WORKING HERE - Need to deal with circular arcs") // treat it as a parabola for now
         parabolic_arc = (Ifc4x1::IfcAlignment2DVerSegParabolicArc*)circular_arc;
         ParabolicSegment(parabolic_arc);
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
void CIfcAlignmentConverter::GetPoint(Ifc4x1::IfcCartesianPoint* pPoint, Float64* pX, Float64* pY)
{
   auto coordinates = pPoint->Coordinates();
   ATLASSERT(2 <= coordinates.size());
   *pX = coordinates[0];
   *pY = coordinates[1];
}

Float64 CIfcAlignmentConverter::OnLine(Float64 startStation, Ifc4x1::IfcLineSegment2D* pLine)
{
   Float64 sx, sy;
   GetPoint(pLine->StartPoint(), &sx, &sy);

   Float64 end_station = startStation + pLine->SegmentLength();

   if (!m_bAlignmentStarted)
   {
      // the bridge starts somewhere in this line segment
      m_AlignmentData.RefStation = startStation;
      m_AlignmentData.xRefPoint = sx;
      m_AlignmentData.yRefPoint = sy;
      m_AlignmentData.Direction = pLine->StartDirection();

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
         hcData.FwdTangent = pLine->StartDirection();

         m_AlignmentData.HorzCurves.push_back(hcData);
      }
   }

   m_LastAlignmentType = Line;

   return end_station;
}

Float64 CIfcAlignmentConverter::OnCurve(Float64 startStation, Ifc4x1::IfcTransitionCurveSegment2D* pEntrySpiral, Ifc4x1::IfcCircularArcSegment2D* pCurve, Ifc4x1::IfcTransitionCurveSegment2D* pExitSpiral)
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
      GetSpiralPoints(pEntrySpiral, &pntEntryStart, &pntEntryPI, &pntEntryEnd);
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
         CheckSpiralType(pEntrySpiral);
      }
   }

   GetCurvePoints(pCurve, &pntCurveStart, &pntCurvePI, &pntCurveEnd, &pntCurveCenter);

   if (pntCurveStart == nullptr || pntCurvePI == nullptr || pntCurveEnd == nullptr || pntCurveCenter == nullptr)
   {
      m_Notes.push_back(std::_tstring(_T("Zero radius curve could not be constructed.")));
      return startStation;
   }

   if (pExitSpiral)
   {
      GetSpiralPoints(pExitSpiral, &pntExitStart, &pntExitPI, &pntExitEnd);
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
         CheckSpiralType(pExitSpiral);
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

      if (SameLocation(pntEntryEnd, pntCurveStart) == S_FALSE)
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

      if (SameLocation(pntCurveEnd, pntExitStart) == S_FALSE)
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

      if (SameLocation(pntEntryEnd, pntCurveStart) == S_FALSE)
      {
         m_Notes.push_back(std::_tstring(_T("The end of the entry spiral does not coincide with the start of the circular curve. The entry spiral has been adjusted.")));
      }

      if (SameLocation(pntCurveEnd, pntExitStart) == S_FALSE)
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

void CIfcAlignmentConverter::GetCurvePoints(Ifc4x1::IfcCircularArcSegment2D* pCurve, IPoint2d** ppStart, IPoint2d** ppPI, IPoint2d** ppEnd, IPoint2d** ppCenter)
{
   auto pStart = pCurve->StartPoint();
   auto bkTangentBrg = pCurve->StartDirection();
   auto L = pCurve->SegmentLength();
   auto R = pCurve->Radius();

   Float64 delta = L / R;
   Float64 T = R*tan(delta/2);

   Float64 sx, sy;
   GetPoint(pStart, &sx, &sy);
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

void CIfcAlignmentConverter::GetSpiralPoints(Ifc4x1::IfcTransitionCurveSegment2D* pSpiral, IPoint2d** ppStart, IPoint2d** ppPI, IPoint2d** ppEnd)
{
   auto pStart = pSpiral->StartPoint();
   auto bkTangentBrg = pSpiral->StartDirection();
   auto L = pSpiral->SegmentLength();
   auto R = (pSpiral->hasStartRadius() ? pSpiral->StartRadius() : pSpiral->EndRadius());
   bool bIsCCW = (pSpiral->hasStartRadius() ? pSpiral->IsStartRadiusCCW() : pSpiral->IsEndRadiusCCW());

   Float64 sx, sy;
   GetPoint(pStart, &sx, &sy);
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

void CIfcAlignmentConverter::LinearSegment(Ifc4x1::IfcAlignment2DVerSegLine* pLinearSegment)
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

void CIfcAlignmentConverter::ParabolicSegment(Ifc4x1::IfcAlignment2DVerSegParabolicArc* pParaCurve)
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

void CIfcAlignmentConverter::CheckSpiralType(Ifc4x1::IfcTransitionCurveSegment2D* pSpiral)
{
   switch (pSpiral->TransitionCurveType())
   {
   case Ifc4x1::IfcTransitionCurveType::IfcTransitionCurveType_BIQUADRATICPARABOLA:
   case Ifc4x1::IfcTransitionCurveType::IfcTransitionCurveType_BLOSSCURVE:
   case Ifc4x1::IfcTransitionCurveType::IfcTransitionCurveType_COSINECURVE:
   case Ifc4x1::IfcTransitionCurveType::IfcTransitionCurveType_CUBICPARABOLA:
   case Ifc4x1::IfcTransitionCurveType::IfcTransitionCurveType_SINECURVE:
      m_Notes.push_back(std::_tstring(_T("Spiral type not supported. Assuming clothoid")));
      break;

   case Ifc4x1::IfcTransitionCurveType::IfcTransitionCurveType_CLOTHOIDCURVE:
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



#if defined USE_IFC4X1
IfcSchema::IfcAlignment2DHorizontal* BuildAlignment(IBroker* pBroker, IfcHierarchyHelper<IfcSchema>& file)
#elif defined USE_IFC4X3_RC2
IfcSchema::IfcAlignmentHorizontal* BuildAlignment(IBroker* pBroker, IfcHierarchyHelper<IfcSchema>& file)
#endif
{
#if defined USE_IFC4X1
   boost::shared_ptr<IfcTemplatedEntityList<IfcSchema::IfcAlignment2DHorizontalSegment>> alignment_segments(new IfcTemplatedEntityList<IfcSchema::IfcAlignment2DHorizontalSegment>());
#elif defined USE_IFC4X3_RC2
   boost::shared_ptr<IfcTemplatedEntityList<IfcSchema::IfcAlignmentHorizontalSegment>> alignment_segments(new IfcTemplatedEntityList<IfcSchema::IfcAlignmentHorizontalSegment>());
#endif

   GET_IFACE2(pBroker, IGeometry, pGeometry);
   GET_IFACE2(pBroker, IRoadway, pAlignment);
   Float64 startStation, startElevation, startGrade;
   CComPtr<IPoint2d> startPoint;
   pAlignment->GetStartPoint(2, &startStation, &startElevation, &startGrade, &startPoint);

   // create the start point
   std::vector<Float64> vCoordinates;
   vCoordinates.resize(2);
   startPoint->Location(&vCoordinates[0], &vCoordinates[1]);
   auto ifc_start_point = new IfcSchema::IfcCartesianPoint(vCoordinates);

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
      if (pntTS->SameLocation(prevPoint) == S_FALSE)
      {
         Float64 dist;
         CComPtr<IDirection> direction;
         pGeometry->Inverse(prevPoint, pntTS, &dist, &direction);
         Float64 angle;
         direction->get_Value(&angle);

#if defined USE_IFC4X1
         auto ifc_line_segment = new IfcSchema::IfcAlignment2DHorizontalSegment(boost::none, boost::none, boost::none, new IfcSchema::IfcLineSegment2D(ifc_prev_point, angle, dist));
#elif defined USE_IFC4X3_RC2
         auto ifc_line_segment = new IfcSchema::IfcAlignmentHorizontalSegment(boost::none, boost::none, ifc_prev_point, angle, 0.0, 0.0, dist, boost::none, IfcSchema::IfcAlignmentHorizontalSegmentTypeEnum::IfcAlignmentHorizontalSegmentType_LINE);
#endif
         alignment_segments->push(ifc_line_segment);
      }

      // create this horizontal curve
      CComPtr<IDirection> bkTangentBrg;
      curve->get_BkTangentBrg(&bkTangentBrg);
      Float64 bk_tangent_direction;
      bkTangentBrg->get_Value(&bk_tangent_direction);

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
      auto ifc_start = new IfcSchema::IfcCartesianPoint(vCoordinates);

      if (0.0 < Lspiral[spEntry])
      {
         // there is an entry spiral
#if defined USE_IFC4X1
         auto ifc_entry_spiral = new IfcSchema::IfcAlignment2DHorizontalSegment(boost::none, std::string("TS"), std::string("SC"), new IfcSchema::IfcTransitionCurveSegment2D(ifc_start, bk_tangent_direction, Lspiral[spEntry], boost::none, R, bIsCCW, bIsCCW, IfcSchema::IfcTransitionCurveType::IfcTransitionCurveType_CLOTHOIDCURVE));
#elif defined USE_IFC4X3_RC2
         auto ifc_entry_spiral = new IfcSchema::IfcAlignmentHorizontalSegment(std::string("TS"), std::string("SC"), ifc_start, bk_tangent_direction, 0.0, (bIsCCW ? 1.0 : -1.0)*R, Lspiral[spEntry], boost::none, IfcSchema::IfcAlignmentHorizontalSegmentTypeEnum::IfcAlignmentHorizontalSegmentType_CLOTHOID);
#endif
         alignment_segments->push(ifc_entry_spiral);

         // spiral ends at the Spiral to Curve point (CS)
         CComPtr<IPoint2d> pntSC;
         pAlignment->GetCurvePoint(i, cptSC, pgsTypes::pcGlobal, &pntSC);
         pntSC->Location(&vCoordinates[0], &vCoordinates[1]);
         ifc_start = new IfcSchema::IfcCartesianPoint(vCoordinates);
      }

      // build the horizontal curve
#if defined USE_IFC4X1
      auto ifc_hcurve = new IfcSchema::IfcAlignment2DHorizontalSegment(boost::none, std::string("SC"), std::string("CS"), new IfcSchema::IfcCircularArcSegment2D(ifc_start, bk_tangent_direction, Lc, R, bIsCCW));
#elif defined USE_IFC4X3_RC2
      auto ifc_hcurve = new IfcSchema::IfcAlignmentHorizontalSegment(std::string("TS"), std::string("SC"), ifc_start, bk_tangent_direction, (bIsCCW ? 1.0 : -1.0)*R, (bIsCCW ? 1.0 : -1.0)*R, Lc, boost::none, IfcSchema::IfcAlignmentHorizontalSegmentTypeEnum::IfcAlignmentHorizontalSegmentType_CIRCULARARC);
#endif
      alignment_segments->push(ifc_hcurve);

      if (0.0 < Lspiral[spExit])
      {
         // there is an exit spiral

         // spiral starts at the Curve to Spiral point (CS)
         CComPtr<IPoint2d> pntCS;
         pAlignment->GetCurvePoint(i, cptCS, pgsTypes::pcGlobal, &pntCS);
         pntCS->Location(&vCoordinates[0], &vCoordinates[1]);
         ifc_start = new IfcSchema::IfcCartesianPoint(vCoordinates);

#if defined USE_IFC4X1
         auto ifc_exit_spiral = new IfcSchema::IfcAlignment2DHorizontalSegment(boost::none, std::string("CS"), std::string("ST"), new IfcSchema::IfcTransitionCurveSegment2D(ifc_start, bk_tangent_direction, Lspiral[spExit], R, boost::none, bIsCCW, bIsCCW, IfcSchema::IfcTransitionCurveType::IfcTransitionCurveType_CLOTHOIDCURVE));
#elif defined USE_IFC4X3_RC2
         auto ifc_exit_spiral = new IfcSchema::IfcAlignmentHorizontalSegment(std::string("TS"), std::string("SC"), ifc_start, bk_tangent_direction, (bIsCCW ? 1.0 : -1.0)*R, 0.0, Lspiral[spExit], boost::none, IfcSchema::IfcAlignmentHorizontalSegmentTypeEnum::IfcAlignmentHorizontalSegmentType_CLOTHOID);
#endif
         alignment_segments->push(ifc_exit_spiral);
      }

      // end of this curve (Spiral to Tangent, ST) becomes previous point for next alignment segment
      prevPoint.Release();
      pAlignment->GetCurvePoint(i, cptST, pgsTypes::pcGlobal, &prevPoint);
      prevPoint->Location(&vCoordinates[0], &vCoordinates[1]);
      ifc_prev_point = new IfcSchema::IfcCartesianPoint(vCoordinates);
   }

   // build a linear segment from end of previous alignment segment to the end of the alignment
   Float64 endStation, endElevation, endGrade;
   CComPtr<IPoint2d> endPoint;
   pAlignment->GetEndPoint(2, &endStation, &endElevation, &endGrade, &endPoint);

   if (prevPoint->SameLocation(endPoint) == S_FALSE)
   {
      // end the alignment with a line segment
      Float64 dist;
      CComPtr<IDirection> direction;
      pGeometry->Inverse(prevPoint, endPoint, &dist, &direction);
      Float64 angle;
      direction->get_Value(&angle);

#if defined USE_IFC4X1
      auto ifc_line_segment = new IfcSchema::IfcAlignment2DHorizontalSegment(boost::none, boost::none, std::string("End"), new IfcSchema::IfcLineSegment2D(ifc_prev_point, angle, dist));
#elif defined USE_IFC4X3_RC2
      auto ifc_line_segment = new IfcSchema::IfcAlignmentHorizontalSegment(boost::none, std::string("End"), ifc_start_point, angle, 0.0, 0.0, dist, boost::none, IfcSchema::IfcAlignmentHorizontalSegmentTypeEnum::IfcAlignmentHorizontalSegmentType_LINE);
#endif
      alignment_segments->push(ifc_line_segment);
   }

   // create a horizontal alignment from all the alignment segments
#if defined USE_IFC4X1
   auto horizontal_alignment = new IfcSchema::IfcAlignment2DHorizontal(startStation, alignment_segments);
#elif defined USE_IFC4X3_RC2
   auto horizontal_alignment = new IfcSchema::IfcAlignmentHorizontal(IfcParse::IfcGlobalId(), file.getSingle<IfcSchema::IfcOwnerHistory>(), boost::none, boost::none, boost::none, file.getSingle<IfcSchema::IfcLocalPlacement>(), file.getSingle<IfcSchema::IfcProductRepresentation>(), startStation, alignment_segments);
#endif

   return horizontal_alignment;
}

#if defined USE_IFC4X1
IfcSchema::IfcAlignment2DVertical* BuildProfile(IBroker* pBroker, IfcHierarchyHelper<IfcSchema>& file)
#elif defined USE_IFC4X3_RC2
IfcSchema::IfcAlignmentVertical* BuildProfile(IBroker* pBroker, IfcHierarchyHelper<IfcSchema>& file)
#endif
{
#if defined USE_IFC4X1
   boost::shared_ptr<IfcTemplatedEntityList<IfcSchema::IfcAlignment2DVerticalSegment>> profile_segments(new IfcTemplatedEntityList<IfcSchema::IfcAlignment2DVerticalSegment>());
#elif defined USE_IFC4X3_RC2
   boost::shared_ptr<IfcTemplatedEntityList<IfcSchema::IfcAlignmentVerticalSegment>> profile_segments(new IfcTemplatedEntityList<IfcSchema::IfcAlignmentVerticalSegment>());
#endif

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
#if defined USE_IFC4X1
         profile_segments->push(new IfcSchema::IfcAlignment2DVerSegLine(boost::none, std::string("Start"), boost::none, prev_end_dist_along, length, prev_end_height, prev_end_gradiant));
#elif defined USE_IFC4X3_RC2
         profile_segments->push(new IfcSchema::IfcAlignmentVerticalSegment(std::string("Start"), boost::none, prev_end_dist_along, length, prev_end_height, prev_end_gradiant, prev_end_gradiant, 0, IfcSchema::IfcAlignmentVerticalSegmentTypeEnum::IfcAlignmentVerticalSegmentType_CONSTANTGRADIENT));
#endif
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

         Float64 start_gradient;
         curve->get_EntryGrade(&start_gradient);

         Float64 k1, k2;
         curve->get_K1(&k1);
         curve->get_K2(&k2);

#if defined USE_IFC4X1
         profile_segments->push(new IfcSchema::IfcAlignment2DVerSegParabolicArc(boost::none, std::string("BVC"), boost::none, start_dist_along, l1, start_height, start_gradient, fabs(1 / k1), ::BinarySign(k1) < 0 ? true : false));
         profile_segments->push(new IfcSchema::IfcAlignment2DVerSegParabolicArc(boost::none, boost::none, std::string("EVC"), start_dist_along + l1, l2, pviElevation, pviGrade, fabs(1 / k2), ::BinarySign(k2) < 0 ? true : false));
#elif defined USE_IFC4X3_RC2
         profile_segments->push(new IfcSchema::IfcAlignmentVerticalSegment(boost::none, boost::none, start_dist_along, l1, start_height, start_gradient, pviGrade, k1, IfcSchema::IfcAlignmentVerticalSegmentTypeEnum::IfcAlignmentVerticalSegmentType_PARABOLICARC));
         profile_segments->push(new IfcSchema::IfcAlignmentVerticalSegment(boost::none, boost::none, start_dist_along + l1, l2, pviElevation, pviGrade, end_gradient, k2, IfcSchema::IfcAlignmentVerticalSegmentTypeEnum::IfcAlignmentVerticalSegmentType_PARABOLICARC));
#endif
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
#if defined USE_IFC4X1
            profile_segments->push(new IfcSchema::IfcAlignment2DVerSegLine(boost::none, boost::none, boost::none, prev_end_dist_along, l1, prev_end_height, start_gradient));
#elif defined USE_IFC4X3_RC2
            profile_segments->push(new IfcSchema::IfcAlignmentVerticalSegment(std::string("Start"), boost::none, prev_end_dist_along, l1, prev_end_height, prev_end_gradiant, prev_end_gradiant, 0, IfcSchema::IfcAlignmentVerticalSegmentTypeEnum::IfcAlignmentVerticalSegmentType_CONSTANTGRADIENT));
#endif
         }
         else
         {
            Float64 k;
            curve->get_K1(&k);
#if defined USE_IFC4X1
            profile_segments->push(new IfcSchema::IfcAlignment2DVerSegParabolicArc(boost::none, std::string("BVC"), std::string("EVC"), start_dist_along, horizontal_length, start_height, start_gradient, fabs(1 / k), ::BinarySign(k) < 0 ? true : false));
#elif defined USE_IFC4X3_RC2
            profile_segments->push(new IfcSchema::IfcAlignmentVerticalSegment(std::string("BVC"), std::string("EVC"), start_dist_along, horizontal_length, start_height, start_gradient, end_gradient, radius_of_curvature, IfcSchema::IfcAlignmentVerticalSegmentTypeEnum::IfcAlignmentVerticalSegmentType_PARABOLICARC));
#endif
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
#if defined USE_IFC4X1
      profile_segments->push(new IfcSchema::IfcAlignment2DVerSegLine(boost::none, boost::none, std::string("End"), prev_end_dist_along, length, prev_end_height, prev_end_gradiant));
#elif defined USE_IFC4X3_RC2
      profile_segments->push(new IfcSchema::IfcAlignmentVerticalSegment(boost::none, std::string("End"), prev_end_dist_along, length, prev_end_height, prev_end_gradiant, prev_end_gradiant, 0, IfcSchema::IfcAlignmentVerticalSegmentTypeEnum::IfcAlignmentVerticalSegmentType_CONSTANTGRADIENT));
#endif

      // check elevation
      ATLASSERT(IsEqual(endElevation, prev_end_height + length*prev_end_gradiant));
   }

#if defined USE_IFC4X1
   auto vertical_profile = new IfcSchema::IfcAlignment2DVertical(profile_segments);
#elif defined USE_IFC4X3_RC2
   auto vertical_profile = new IfcSchema::IfcAlignmentVertical(IfcParse::IfcGlobalId(), file.getSingle<IfcSchema::IfcOwnerHistory>(), boost::none, boost::none, boost::none, file.getSingle<IfcSchema::IfcLocalPlacement>(), file.getSingle<IfcSchema::IfcProductRepresentation>(), profile_segments);
#endif

   return vertical_profile;
}

bool CIfcAlignmentConverter::IsValidAlignment(Ifc4x1::IfcAlignment* pAlignment)
{
   auto axis = pAlignment->Axis()->as<Ifc4x1::IfcAlignmentCurve>();
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
      auto transition = segment->CurveGeometry()->as<Ifc4x1::IfcTransitionCurveSegment2D>();
      return (transition == nullptr ? true : false);
   }
   else if (nSegments == 2)
   {
      // our model doesn't support isolated transition segments
      // transition curves must be adjacent to circular curves
      // can't have two transitions adjacent to each other either
      auto segment1 = (*segments->begin());
      auto transition1 = segment1->CurveGeometry()->as<Ifc4x1::IfcTransitionCurveSegment2D>();

      auto segment2 = (*segments->begin() + 1);
      auto transition2 = segment2->CurveGeometry()->as<Ifc4x1::IfcTransitionCurveSegment2D>();
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
         auto transition = segment->CurveGeometry()->as<Ifc4x1::IfcTransitionCurveSegment2D>();
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
               auto arc = (*(iter - 1))->CurveGeometry()->as<Ifc4x1::IfcCircularArcSegment2D>();
               if (!arc) return false; // previous is not an arc

               if (!IsEqual(arc->Radius(), start_radius)) return false; // common radii must be equal

               if (arc->IsCCW() != transition->IsStartRadiusCCW()) return false; // curves must be same direction
            }

            if (iter != end - 1 && !IsZero(end_radius))
            {
               // transition ends with a radius so a circular curve most come after this transition curve
               auto arc = (*(iter + 1))->CurveGeometry()->as<Ifc4x1::IfcCircularArcSegment2D>();
               if (!arc) return false; // previous is not an arc

               if (!IsEqual(arc->Radius(), end_radius)) return false; // common radii must be equal

               if (arc->IsCCW() != transition->IsEndRadiusCCW()) return false; // curves must be same direction
            }
         }
      }
   }


   return true;
}

Ifc4x1::IfcAlignment* CIfcAlignmentConverter::GetAlignment(IfcParse::IfcFile& file)
{
   USES_CONVERSION;



   //std::ostringstream os;
   //os << "Schema: " << file.schema()->name() << std::endl;
   //AfxMessageBox(A2T(os.str().c_str()));

   Ifc4x1::IfcAlignment::list::ptr alignments = file.instances_by_type<Ifc4x1::IfcAlignment>();
   std::vector<Ifc4x1::IfcAlignment*> valid_alignments;

   for (auto alignment : *alignments)
   {
      if (IsValidAlignment(alignment))
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
