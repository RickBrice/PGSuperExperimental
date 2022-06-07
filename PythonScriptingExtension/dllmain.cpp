#include "stdafx.h"

// Experimental code for creating Python bindings to PGSuper
// Currently using Boost.Python, however it is starting to become out of date with Python
// Another option is pybind (https://github.com/pybind, https://pybind11.readthedocs.io/en/latest/)
#pragma warning (disable : 4267)

#include "resource.h"

#include "PyEditorDocTemplate.h"
#include "PyEditorDoc.h"
#include "PyEditorFrame.h"
#include "PyEditorView.h"

#include <initguid.h>
#include "dllmain.h"
#include <WBFLCore_i.c>
#include "PGSuperScriptingAgentImp.h"
#include "PGSuperCatCom.h"
#include "PGSpliceCatCom.h"
#include <EAF\EAFUIIntegration.h>
#include <EAF\EAFStatusCenter.h>
#include <IGraphManager.h>
#include <IFace\Selection.h>
#include <IFace\AnalysisResults.h>
#include "Scripting.h"

#include "Bindings.h"

#include "ExceptionHandlers.h"

#include <boost\python\suite\indexing\vector_indexing_suite.hpp>

// include files for the "export" functions
#include "PyPGSuperTypes.h"
#include "PyViewTypes.h"
#include "PyGraphTypes.h"
#include "PyBridgeTypes.h"
#include "PyPointOfInterestTypes.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

CPGSuperScriptingExtensionModule _AtlModule;

CPyPGSuper& GetApp() 
{ 
   AFX_MANAGE_STATE(AfxGetStaticModuleState());
   CExtensionApp* pApp = (CExtensionApp*)AfxGetApp();
   return pApp->m_pyPGSuper;  
}

inline GroupIndexType AllGroups() { return ALL_GROUPS; }
inline GirderIndexType AllGirders() { return ALL_GIRDERS; }
inline SegmentIndexType AllSegments() { return ALL_SEGMENTS; }
inline SpanIndexType AllSpans() { return ALL_SPANS; }
inline PierIndexType AllPiers() { return ALL_PIERS; }
inline StageIndexType AllStages() { return ALL_STAGES; }
inline DuctIndexType AllDucts() { return ALL_DUCTS; }

BOOST_PYTHON_MODULE(BridgeLink)
{
   using namespace boost::python;

   register_exception_translator<CXUnwind*>(handle_unwind_exception);

   def("GetPGSuper", GetApp, return_value_policy<reference_existing_object>());

   export_pgsuper_types();
   export_view_types();
   export_graph_types();
   export_bridge_types();
   export_point_of_interest_types();

   // need a better way of dealing with C++ #define constants
   def("AllGroups", AllGroups);
   def("AllGirders", AllGirders);
   def("AllSegments", AllSegments);
   def("AllSpans", AllSpans);
   def("AllPiers", AllPiers);
   def("AllStages", AllStages);
   def("AllDucts", AllDucts);

   class_<MINMOMENTCAPDETAILS>("MinimumMomentCapacityDetails")
      .def_readonly("Mcr", &MINMOMENTCAPDETAILS::Mcr)
      .def_readonly("Mu", &MINMOMENTCAPDETAILS::Mu)
      ;

   class_<CGirderKey>("GirderKey")
      .def(init<GroupIndexType,GirderIndexType>())
      .def(init<const CGirderKey&>())
      .def_readwrite("groupIndex", &CGirderKey::groupIndex)
      .def_readwrite("girderIndex", &CGirderKey::girderIndex)
      .def("IsEqual",&CGirderKey::IsEqual)
      .def(self == self) // __eq__
      .def(self != self) // __ne__
      .def(self < self) // __lt__
      ;

   class_<CSegmentKey, bases<CGirderKey>>("SegmentKey")
      .def(init<GroupIndexType, GirderIndexType, SegmentIndexType>())
      .def(init<const CSegmentKey&>())
      .def_readwrite("segmentIndex",&CSegmentKey::segmentIndex)
      .def("IsEqual", &CSegmentKey::IsEqual)
      .def(self == self) // __eq__
      .def(self != self) // __ne__
      .def(self < self) // __lt__
      ;

   class_<CSpanKey>("SpanKey")
      .def(init<SpanIndexType, GirderIndexType>())
      .def(init<const CSpanKey&>())
      .def_readwrite("spanIndex", &CSpanKey::spanIndex)
      .def_readwrite("girderIndex", &CSpanKey::girderIndex)
      .def("IsEqual", &CSpanKey::IsEqual)
      .def(self == self) // __eq__
      .def(self != self) // __ne__
      .def(self < self) // __lt__
      ;

#pragma Reminder("WORKING HERE - need to flesh out pgsPointOfInterest")
   class_<pgsPointOfInterest>("PointOfInterest")
      .def(init<const pgsPointOfInterest&>())
      .def(init<const CSegmentKey&, Float64, PoiAttributeType>())
      .def("GetSegmentKey",&pgsPointOfInterest::GetSegmentKey, return_internal_reference<>())
      .def("GetDistFromStart",&pgsPointOfInterest::GetDistFromStart)
      ;
   
   class_<std::vector<pgsPointOfInterest>>("PointsOfInterest")
      .def(vector_indexing_suite<std::vector<pgsPointOfInterest>>())
      ;


   class_<CPyProjectProperties>("ProjectProperties")
      .add_property("BridgeName", &CPyProjectProperties::GetBridgeName, &CPyProjectProperties::SetBridgeName)
      .add_property("BridgeID", &CPyProjectProperties::GetBridgeID, &CPyProjectProperties::SetBridgeID)
      .add_property("JobNumber", &CPyProjectProperties::GetJobNumber, &CPyProjectProperties::SetJobNumber)
      .add_property("Engineer", &CPyProjectProperties::GetEngineer, &CPyProjectProperties::SetEngineer)
      .add_property("Company", &CPyProjectProperties::GetCompany, &CPyProjectProperties::SetCompany)
      .add_property("Comments", &CPyProjectProperties::GetComments, &CPyProjectProperties::SetComments)
      ;


   enum_<enumExposureCondition>("ExposureCondition")
      .value("Normal", expNormal)
      .value("Severe", expSevere)
      ;

   class_<CPyEnvironment>("Environment")
      .add_property("ExposureCondition", &CPyEnvironment::GetExposureCondition, &CPyEnvironment::SetExposureCondition)
      .add_property("RelHumidity", &CPyEnvironment::GetRelHumidity, &CPyEnvironment::SetRelHumidity)
      ;

   class_<CPyGraphManager>("Graphs")
      .def("Count", &CPyGraphManager::GetGraphBuilderCount)
      .def("Names", &CPyGraphManager::GetGraphNames)
      ;

   class_<CPyReportManager>("Reports")
      .def("Count", &CPyReportManager::GetReportBuilderCount)
      .def("Names", &CPyReportManager::GetReportNames)
      ;

   class_<CompoundCurveData>("HorzCurve")
      .def_readwrite("PIStation", &CompoundCurveData::PIStation)
      .def_readwrite("Radius", &CompoundCurveData::Radius)
      .def_readwrite("FwdTangent", &CompoundCurveData::FwdTangent)
      .def_readwrite("IsFwdTangent", &CompoundCurveData::bFwdTangent)
      .def_readwrite("EntrySpiral", &CompoundCurveData::EntrySpiral)
      .def_readwrite("ExitSpiral", &CompoundCurveData::ExitSpiral)
      ;

   class_<std::vector<CompoundCurveData>>("HorzCurves")
      .def(vector_indexing_suite<std::vector<CompoundCurveData>>())
      ;

   class_<AlignmentData2>("Alignment")
      .def_readwrite("Direction", &AlignmentData2::Direction)
      .def_readwrite("RefStation", &AlignmentData2::RefStation)
      .def_readwrite("RefPointX", &AlignmentData2::xRefPoint)
      .def_readwrite("RefPointY", &AlignmentData2::yRefPoint)
      .def_readwrite("HorzCurves", &AlignmentData2::CompoundCurves)
      ;

   class_<VertCurveData>("VertCurve")
      .def_readwrite("PVIStation", &VertCurveData::PVIStation)
      .def_readwrite("ExitGrade", &VertCurveData::ExitGrade)
      .def_readwrite("L1", &VertCurveData::L1)
      .def_readwrite("L2", &VertCurveData::L2)
      ;

   class_<std::vector<VertCurveData>>("VertCurves")
      .def(vector_indexing_suite<std::vector<VertCurveData>>())
      ;

   class_<ProfileData2>("Profile")
      .def_readwrite("Station", &ProfileData2::Station)
      .def_readwrite("Elevation", &ProfileData2::Elevation)
      .def_readwrite("Grade", &ProfileData2::Grade)
      .def_readwrite("VertCurves", &ProfileData2::VertCurves)
      ;
   
   class_<std::vector<RoadwaySectionTemplate>>("Templates")
      .def(vector_indexing_suite<std::vector<RoadwaySectionTemplate>>())
      ;
   
   class_<std::vector<RoadwaySegmentData>>("TemplateSegments")
      .def(vector_indexing_suite<std::vector<RoadwaySegmentData>>())
      ;

   class_<RoadwaySegmentData>("TemplateSegment")
      .def_readwrite("Length", &RoadwaySegmentData::Length)
      .def_readwrite("Slope", &RoadwaySegmentData::Slope)
      ;

   class_<RoadwaySectionData>("RoadwaySections")
      .def_readwrite("SegmentsPerSection", &RoadwaySectionData::NumberOfSegmentsPerSection)
      .def_readwrite("AlignmentPointIndex", &RoadwaySectionData::AlignmentPointIdx)
      .def_readwrite("ProfileGradePointIndex",&RoadwaySectionData::ProfileGradePointIdx)
      .def_readwrite("Templates", &RoadwaySectionData::RoadwaySectionTemplates)
      ;

   class_<RoadwaySectionTemplate>("SectionTemplate")
      .def_readwrite("Station", &RoadwaySectionTemplate::Station)
      .def_readwrite("LeftSlope", &RoadwaySectionTemplate::LeftSlope)
      .def_readwrite("RightSlope", &RoadwaySectionTemplate::RightSlope)
      ;

   class_<CPyRoadway>("Roadway")
      .add_property("Alignment", &CPyRoadway::GetAlignmentData, &CPyRoadway::SetAlignmentData)
      .add_property("Profile", &CPyRoadway::GetProfileData, &CPyRoadway::SetProfileData)
      .add_property("RoadwaySections", &CPyRoadway::GetRoadwaySectionData, &CPyRoadway::SetRoadwaySectionData)
      ;

   enum_<ILoadModifiers::Level>("LoadModifierLevel")
      .value("Low", ILoadModifiers::Low)
      .value("Normal", ILoadModifiers::Normal)
      .value("High", ILoadModifiers::High)
      ;

   class_<CPyLoadModifiers>("LoadModifiers")
      .def("SetDuctilityModifier", &CPyLoadModifiers::SetDuctilityModifier)
      .def("GetDuctilityLevel", &CPyLoadModifiers::GetDuctilityLevel)
      .def("GetDuctilityFactor", &CPyLoadModifiers::GetDuctilityFactor)
      .def("SetImportanceModifier", &CPyLoadModifiers::SetImportanceModifier)
      .def("GetImportanceLevel", &CPyLoadModifiers::GetImportanceLevel)
      .def("GetImportanceFactor", &CPyLoadModifiers::GetImportanceFactor)
      .def("SetRedundancyModifier", &CPyLoadModifiers::SetRedundancyModifier)
      .def("GetRedundancyLevel", &CPyLoadModifiers::GetRedundancyLevel)
      .def("GetRedundancyFactor", &CPyLoadModifiers::GetRedundancyFactor)
      ;

   class_<CPyLimits>("Limits")
      .def("GetMaxSlabFc", &CPyLimits::GetMaxSlabFc)
      .def("GetMaxSegmentFci", &CPyLimits::GetMaxSegmentFci)
      .def("GetMaxSegmentFc", &CPyLimits::GetMaxSegmentFc)
      .def("GetMaxClosureFci", &CPyLimits::GetMaxClosureFci)
      .def("GetMaxClosureFc", &CPyLimits::GetMaxClosureFc)
      .def("GetMaxConcreteUnitWeight", &CPyLimits::GetMaxConcreteUnitWeight)
      .def("GetMaxConcreteAggSize", &CPyLimits::GetMaxConcreteAggSize)
      ;

   class_<CPyEffectiveFlangeWidth>("EffectiveFlangeWidth")
      .def("GetIgnoreEffectiveFlangeWidthLimits", &CPyEffectiveFlangeWidth::GetIgnoreEffectiveFlangeWidthLimits)
      .def("SetIgnoreEffectiveFlangeWidthLimits", &CPyEffectiveFlangeWidth::SetIgnoreEffectiveFlangeWidthLimits)
      ;

   class_<CPyLossParameters>("LossParameters")
      .def("GetLossMethod", &CPyLossParameters::GetLossMethod)
      .def("GetTimeDependentModel", &CPyLossParameters::GetTimeDependentModel)
      .def("SetIgnoreCreepEffects", &CPyLossParameters::SetIgnoreCreepEffects)
      .def("GetIgnoreCreepEffects", &CPyLossParameters::GetIgnoreCreepEffects)
      .def("SetIgnoreShrinkageEffects", &CPyLossParameters::SetIgnoreShrinkageEffects)
      .def("GetIgnoreShrinkageEffects", &CPyLossParameters::GetIgnoreShrinkageEffects)
      .def("SetIgnoreRelaxationEffects", &CPyLossParameters::SetIgnoreRelaxationEffects)
      .def("GetIgnoreRelaxationEffects", &CPyLossParameters::GetIgnoreRelaxationEffects)
      .def("SetIgnoreTimeDependentEffects", &CPyLossParameters::SetIgnoreTimeDependentEffects)
      .def("SetTendonPostTensionParameters", &CPyLossParameters::SetTendonPostTensionParameters)
      .def("GetTendonPostTensionParameters", &CPyLossParameters::GetTendonPostTensionParameters)
      .def("SetTemporaryStrandPostTensionParameters", &CPyLossParameters::SetTemporaryStrandPostTensionParameters)
      .def("GetTemporaryStrandPostTensionParameters", &CPyLossParameters::GetTemporaryStrandPostTensionParameters)
      .def("SetUseGeneralLumpSumLosses", &CPyLossParameters::SetUseGeneralLumpSumLosses)
      .def("GetUseGeneralLumpSumLosses", &CPyLossParameters::GetUseGeneralLumpSumLosses)
      .def("GetBeforeXferLosses", &CPyLossParameters::GetBeforeXferLosses)
      .def("SetBeforeXferLosses", &CPyLossParameters::SetBeforeXferLosses)
      .def("GetAfterXferLosses", &CPyLossParameters::GetAfterXferLosses)
      .def("SetAfterXferLosses", &CPyLossParameters::SetAfterXferLosses)
      .def("GetLiftingLosses", &CPyLossParameters::GetLiftingLosses)
      .def("SetLiftingLosses", &CPyLossParameters::SetLiftingLosses)
      .def("GetShippingLosses", &CPyLossParameters::GetShippingLosses)
      .def("SetShippingLosses", &CPyLossParameters::SetShippingLosses)
      .def("GetBeforeTempStrandRemovalLosses", &CPyLossParameters::GetBeforeTempStrandRemovalLosses)
      .def("SetBeforeTempStrandRemovalLosses", &CPyLossParameters::SetBeforeTempStrandRemovalLosses)
      .def("GetAfterTempStrandRemovalLosses", &CPyLossParameters::GetAfterTempStrandRemovalLosses)
      .def("SetAfterTempStrandRemovalLosses", &CPyLossParameters::SetAfterTempStrandRemovalLosses)
      .def("GetAfterDeckPlacementLosses", &CPyLossParameters::GetAfterDeckPlacementLosses)
      .def("SetAfterDeckPlacementLosses", &CPyLossParameters::SetAfterDeckPlacementLosses)
      .def("GetAfterSIDLLosses", &CPyLossParameters::GetAfterSIDLLosses)
      .def("SetAfterSIDLLosses", &CPyLossParameters::SetAfterSIDLLosses)
      .def("GetFinalLosses", &CPyLossParameters::GetFinalLosses)
      .def("SetFinalLosses", &CPyLossParameters::SetFinalLosses)
   ;

   class_<CPyPGSuper>("PGSuper")
      .add_property("VersionString", &CPyPGSuper::GetVersionString)
      .add_property("VersionStringWithBuildNumber", &CPyPGSuper::GetVersionStringWithBuildNumber)
      .add_property("Version", &CPyPGSuper::GetVersion)
      .add_property("VersionWithBuildNumber", &CPyPGSuper::GetVersionWithBuildNumber)

      .add_property("Properties", &CPyPGSuper::Properties)
      .add_property("Environment", &CPyPGSuper::Environment)
      .add_property("Roadway", &CPyPGSuper::Roadway)
      .add_property("Views",&CPyPGSuper::Views)
      .add_property("Graphs", &CPyPGSuper::Graphs)
      .add_property("Reports", &CPyPGSuper::Reports)
      .add_property("LoadModifiers",&CPyPGSuper::LoadModifiers)
      .add_property("Limits",&CPyPGSuper::Limits)
      .add_property("EffectiveFlangeWidth",&CPyPGSuper::EffectiveFlangeWidth)
      .add_property("LossParameters",&CPyPGSuper::LossParameters)

      .add_property("Bridge",&CPyPGSuper::Bridge)
      .add_property("PointOfInterest",&CPyPGSuper::PointOfInterest)

      .add_property("MomentCapacity",&CPyPGSuper::MomentCapacity)
      ;
}

BEGIN_MESSAGE_MAP(CExtensionApp, CWinApp)
END_MESSAGE_MAP()

CExtensionApp theApp;

BOOL CExtensionApp::InitInstance()
{
   SetRegistryKey(_T("Washington State Department of Transportation"));
   ReadSettings();

   if (m_pDocManager == nullptr)
   {
      m_pDocManager = new CDocManager();
   }
   CPyEditorDocTemplate* pDocTemplate = new CPyEditorDocTemplate(IDR_PYEDITOR, RUNTIME_CLASS(CPyEditorDoc), RUNTIME_CLASS(CPyEditorFrame), RUNTIME_CLASS(CPyEditorView));
   m_pDocManager->AddDocTemplate(pDocTemplate);


   //Load the scintilla dll
   //m_hSciDLL = LoadLibraryFromApplicationDirectory(_T("SciLexer.dll"));
   //if (m_hSciDLL == nullptr)
   //{
   //   //AfxMessageBox(_T("Scintilla DLL is not installed, Please download the SciTE editor and copy the SciLexer.dll into this application's directory"));
   //   //return FALSE;
   //}

   //HINSTANCE hDLL = LoadLibraryFromApplicationDirectory(_T("PyDebugPP.dll"));

   PyImport_AppendInittab("BridgeLink", &PyInit_BridgeLink);
   Py_Initialize();
   //PyEval_InitThreads(); // initialize threading and aquire the GIL

   ATLASSERT(Py_IsInitialized());

   return CWinApp::InitInstance();
}

int CExtensionApp::ExitInstance()
{
   WriteSettings();

   Py_FinalizeEx();

   if (m_pDocManager)
   {
      delete m_pDocManager;
      m_pDocManager = nullptr;
   }
   
   return CWinApp::ExitInstance();
}

HMODULE CExtensionApp::LoadLibraryFromApplicationDirectory(LPCTSTR lpFileName)
{
   //Get the Application diretory
   TCHAR szFullPath[_MAX_PATH];
   szFullPath[0] = _T('\0');
   if (GetModuleFileName(nullptr, szFullPath, _countof(szFullPath)) == 0)
      return nullptr;

   //Form the new path
   TCHAR szDrive[_MAX_DRIVE];
   szDrive[0] = _T('\0');
   TCHAR szDir[_MAX_DIR];
   szDir[0] = _T('\0');
   _tsplitpath_s(szFullPath, szDrive, sizeof(szDrive) / sizeof(TCHAR), szDir, sizeof(szDir) / sizeof(TCHAR), nullptr, 0, nullptr, 0);
   TCHAR szFname[_MAX_FNAME];
   szFname[0] = _T('\0');
   TCHAR szExt[_MAX_EXT];
   szExt[0] = _T('\0');
   _tsplitpath_s(lpFileName, nullptr, 0, nullptr, 0, szFname, sizeof(szFname) / sizeof(TCHAR), szExt, sizeof(szExt) / sizeof(TCHAR));
   _tmakepath_s(szFullPath, sizeof(szFullPath) / sizeof(TCHAR), szDrive, szDir, szFname, szExt);

   //Delegate to LoadLibrary    
   return LoadLibrary(szFullPath);
}

CString CExtensionApp::GetScriptPath() const
{
   return m_ScriptPath;
}

void CExtensionApp::ReadSettings()
{
   AFX_MANAGE_STATE(AfxGetStaticModuleState());

   LPWSTR path;
   HRESULT hr = ::SHGetKnownFolderPath(FOLDERID_RoamingAppData, KF_FLAG_DEFAULT, NULL, &path);
   ATLASSERT(SUCCEEDED(hr));

   CString strDefaultPath;
   strDefaultPath.Format(_T("%s\\PGSuper\\Scripts"), path);
   m_ScriptPath = GetProfileString(_T("Options"), _T("ScriptPath"), strDefaultPath);
}

void CExtensionApp::WriteSettings()
{
   AFX_MANAGE_STATE(AfxGetStaticModuleState());
   WriteProfileString(_T("Options"), _T("ScriptPath"), m_ScriptPath);
}