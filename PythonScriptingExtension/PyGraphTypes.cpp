#include "stdafx.h"

#pragma warning (disable : 4267)

#include "PyGraphTypes.h"

#include <Graphs\StabilityGraphViewController.h>
#include <Graphs\EffectivePrestressGraphViewController.h>
#include <Graphs\ConcretePropertiesGraphViewController.h>
#include <Graphs\DeflectionHistoryGraphViewController.h>
#include <Graphs\StressHistoryGraphViewController.h>
#include <Graphs\AnalysisResultsGraphViewController.h>

#include "PyStabilityGraphViewController.h"
#include "PyEffectivePrestressGraphViewController.h"
#include "PyConcretePropertiesGraphViewController.h"
#include "PyDeflectionHistoryGraphViewController.h"
#include "PyStressHistoryGraphViewController.h"
#include "PyAnalysisResultsGraphViewController.h"
#include "PyGirderPropertiesGraphViewController.h"

#include <boost\python\suite\indexing\vector_indexing_suite.hpp>


using namespace boost::python;

void export_graph_types()
{
   /////////////////////////////////////////////////////
   // concrete properties graph
   /////////////////////////////////////////////////////
   enum_<IConcretePropertiesGraphViewController::GraphType>("ConcreteProperty")
      .value("Fc", IConcretePropertiesGraphViewController::Fc)
      .value("Ec", IConcretePropertiesGraphViewController::Ec)
      .value("Shrinkage", IConcretePropertiesGraphViewController::Shrinkage)
      .value("Creep", IConcretePropertiesGraphViewController::Creep)
      ;

   enum_<IConcretePropertiesGraphViewController::ElementType>("ConcretePropertiesElement")
      .value("Segment", IConcretePropertiesGraphViewController::Segment)
      //.value("Closure", IConcretePropertiesGraphViewController::Closure) // PGSplice
      .value("Deck", IConcretePropertiesGraphViewController::Deck)
      ;

#pragma Reminder("WORKING HERE - need to consolodate axis type so we don't have a different enum for each graph type")
   enum_<IConcretePropertiesGraphViewController::XAxisType>("ConcretePropertiesXAxis")
      .value("TimeLinear", IConcretePropertiesGraphViewController::TimeLinear)
      .value("TimeLog", IConcretePropertiesGraphViewController::TimeLog)
      .value("AgeLinear", IConcretePropertiesGraphViewController::AgeLinear)
      .value("AgeLog", IConcretePropertiesGraphViewController::AgeLog)
      .value("Interval", IConcretePropertiesGraphViewController::Interval)
      ;

   class_<CPyConcretePropertiesGraphViewController, bases<CPyViewControllerBase>>("ConcretePropertiesGraph")
      .def("GetGraphController", &CPyConcretePropertiesGraphViewController::GetGraphController, return_value_policy<reference_existing_object>())
      .staticmethod("GetGraphController")
      .add_property("Property", &CPyConcretePropertiesGraphViewController::GetGraphType, &CPyConcretePropertiesGraphViewController::SetGraphType)
      .add_property("Element", &CPyConcretePropertiesGraphViewController::GetElementType, &CPyConcretePropertiesGraphViewController::SetElementType)
      .add_property("XAxis", &CPyConcretePropertiesGraphViewController::GetXAxisType, &CPyConcretePropertiesGraphViewController::SetXAxisType)
      .def("GetSegment", &CPyConcretePropertiesGraphViewController::GetSegment, return_internal_reference<>())
      .def("SelectSegment", &CPyConcretePropertiesGraphViewController::SetSegment)
      .def("GetClosureJoint", &CPyConcretePropertiesGraphViewController::GetClosureJoint, return_internal_reference<>())
      .def("SelectClosureJoint", &CPyConcretePropertiesGraphViewController::SetClosureJoint)
      .add_property("ShowGrid", &CPyConcretePropertiesGraphViewController::GetShowGrid, &CPyConcretePropertiesGraphViewController::SetShowGrid)
      ;

   /////////////////////////////////////////////////////
   // effective prestress graph
   /////////////////////////////////////////////////////
   enum_<IEffectivePrestressGraphViewController::ViewMode>("EffectivePrestressGraphMode")
      .value("Stress", IEffectivePrestressGraphViewController::Stress)
      .value("Force", IEffectivePrestressGraphViewController::Force)
      ;

   enum_<IEffectivePrestressGraphViewController::StrandType>("EffectivePrestressStrandType")
      .value("Permanent", IEffectivePrestressGraphViewController::Permanent)
      .value("Temporary", IEffectivePrestressGraphViewController::Temporary)
      ;

   class_<CPyEffectivePrestressGraphViewController, bases<CPyViewControllerBase>>("EffectivePrestressGraph")
      .def("GetGraphController", &CPyEffectivePrestressGraphViewController::GetGraphController, return_value_policy<reference_existing_object>())
      .staticmethod("GetGraphController")
      .def("GetGirder", &CPyEffectivePrestressGraphViewController::GetGirder, return_internal_reference<>())
      .def("SelectGirder", &CPyEffectivePrestressGraphViewController::SelectGirder)
      .def("GetMinimumInterval", &CPyEffectivePrestressGraphViewController::GetMinimumInterval)
      .def("GetMaximumInterval", &CPyEffectivePrestressGraphViewController::GetMaximumInterval)
      .def("SelectInterval", &CPyEffectivePrestressGraphViewController::SelectInterval)
      .def("SelectIntervals", &CPyEffectivePrestressGraphViewController::SelectIntervals)
      .def("GetSelectedIntervals", &CPyEffectivePrestressGraphViewController::GetSelectedIntervals)
      .add_property("Mode", &CPyEffectivePrestressGraphViewController::GetViewMode, &CPyEffectivePrestressGraphViewController::SetViewMode)
      .add_property("StrandType", &CPyEffectivePrestressGraphViewController::GetStrandType, &CPyEffectivePrestressGraphViewController::SetStrandType)
      .add_property("Duct", &CPyEffectivePrestressGraphViewController::GetDuct, &CPyEffectivePrestressGraphViewController::SetDuct)
      .add_property("ShowGrid", &CPyEffectivePrestressGraphViewController::GetShowGrid, &CPyEffectivePrestressGraphViewController::SetShowGrid)
      .add_property("ShowGirder", &CPyEffectivePrestressGraphViewController::GetShowGirder, &CPyEffectivePrestressGraphViewController::SetShowGirder)
      ;

   /////////////////////////////////////////////////////
   // stability graph
   /////////////////////////////////////////////////////

   enum_<IStabilityGraphViewController::ViewMode>("StabilityGraphMode")
      .value("Lifting", IStabilityGraphViewController::Lifting)
      .value("Hauling", IStabilityGraphViewController::Hauling)
      ;

   class_<CPyStabilityGraphViewController, bases<CPyViewControllerBase>>("StabilityGraph")
      .def("GetGraphController", &CPyStabilityGraphViewController::GetGraphController, return_value_policy<reference_existing_object>())
      .staticmethod("GetGraphController")
      .def("GetSegment", &CPyStabilityGraphViewController::GetSegment, return_internal_reference<>())
      .def("SelectSegment", &CPyStabilityGraphViewController::SelectSegment)
      .add_property("Mode", &CPyStabilityGraphViewController::GetViewMode, &CPyStabilityGraphViewController::SetViewMode)
      .add_property("ShowGrid", &CPyStabilityGraphViewController::GetShowGrid, &CPyStabilityGraphViewController::SetShowGrid)
      ;

   /////////////////////////////////////////////////////
   // deflection history graph
   /////////////////////////////////////////////////////
#pragma Reminder("WORKING HERE - need to consolodate axis type so we don't have a different enum for each graph type")
   enum_<IDeflectionHistoryGraphViewController::XAxisType>("DeflectionHistoryXAxis")
      .value("TimeLinear", IDeflectionHistoryGraphViewController::TimeLinear)
      .value("TimeLog", IDeflectionHistoryGraphViewController::TimeLog)
      .value("Interval", IDeflectionHistoryGraphViewController::Interval)
      ;

   class_<CPyDeflectionHistoryGraphViewController, bases<CPyViewControllerBase>>("DeflectionHistoryGraph")
      .def("GetGraphController", &CPyDeflectionHistoryGraphViewController::GetGraphController, return_value_policy<reference_existing_object>())
      .staticmethod("GetGraphController")
      .def("GetLocation",&CPyDeflectionHistoryGraphViewController::GetLocation, return_internal_reference<>())
      .def("SelectLocation", &CPyDeflectionHistoryGraphViewController::SelectLocation)
      .add_property("XAxis", &CPyDeflectionHistoryGraphViewController::GetXAxisType, &CPyDeflectionHistoryGraphViewController::SetXAxisType)
      .add_property("IncludeElevationAdjustment", &CPyDeflectionHistoryGraphViewController::GetIncludeElevationAdjustment, &CPyDeflectionHistoryGraphViewController::SetIncludeElevationAdjustment)
      .add_property("IncludeUnrecoverableDefl", &CPyDeflectionHistoryGraphViewController::GetIncludeUnrecoverableDefl, &CPyDeflectionHistoryGraphViewController::SetIncludeUnrecoverableDefl)
      .add_property("ShowGrid", &CPyDeflectionHistoryGraphViewController::GetShowGrid, &CPyDeflectionHistoryGraphViewController::SetShowGrid)
      ;

   /////////////////////////////////////////////////////
   // stress history graph
   /////////////////////////////////////////////////////

#pragma Reminder("WORKING HERE - need to consolodate axis type so we don't have a different enum for each graph type")
   enum_<IStressHistoryGraphViewController::XAxisType>("StressHistoryXAxis")
      .value("TimeLinear", IStressHistoryGraphViewController::TimeLinear)
      .value("TimeLog", IStressHistoryGraphViewController::TimeLog)
      .value("Interval", IStressHistoryGraphViewController::Interval)
      ;

   class_<CPyStressHistoryGraphViewController, bases<CPyViewControllerBase>>("StressHistoryGraph")
      .def("GetGraphController", &CPyStressHistoryGraphViewController::GetGraphController, return_value_policy<reference_existing_object>())
      .staticmethod("GetGraphController")
      .def("GetLocation",&CPyStressHistoryGraphViewController::GetLocation, return_internal_reference<>())
      .def("SelectLocation", &CPyStressHistoryGraphViewController::SelectLocation)
      .add_property("XAxis", &CPyStressHistoryGraphViewController::GetXAxisType, &CPyStressHistoryGraphViewController::SetXAxisType)
      .def("SelectStressLocation", &CPyStressHistoryGraphViewController::SelectStressLocation)
      .def("IsStressLocationSelected", &CPyStressHistoryGraphViewController::IsStressLocationSelected)
      .add_property("ShowGrid", &CPyStressHistoryGraphViewController::GetShowGrid, &CPyStressHistoryGraphViewController::SetShowGrid)
      ;

   /////////////////////////////////////////////////////
   // analysis results graph
   /////////////////////////////////////////////////////

   enum_<IAnalysisResultsGraphViewController::GraphMode>("AnalysisResultsGraphMode")
      .value("Interval", IAnalysisResultsGraphViewController::Interval)
      .value("Loading", IAnalysisResultsGraphViewController::Loading)
      ;

   enum_<ActionType>("ActionType")
      .value("Axial", actionAxial)
      .value("Shear", actionShear)
      .value("Moment", actionMoment)
      .value("Reaction", actionReaction)
      .value("Deflection", actionDeflection)
      .value("XDeflection", actionXDeflection)
      .value("Rotation", actionRotation)
      .value("Stress", actionStress)
      ;

   class_<std::vector<ActionType>>("ActionTypes")
      .def(vector_indexing_suite<std::vector<ActionType>>())
      ;

   enum_<ResultsType>("AnalysisResultsType")
      .value("Cummulative", rtCumulative)
      .value("Incremental", rtIncremental)
      ;


   class_<CPyAnalysisResultsGraphViewController, bases<CPyViewControllerBase>>("AnalysisResultsGraph")
      .def("GetGraphController", &CPyAnalysisResultsGraphViewController::GetGraphController, return_value_policy<reference_existing_object>())
      .staticmethod("GetGraphController")
      .add_property("GraphMode", &CPyAnalysisResultsGraphViewController::GetGraphMode, &CPyAnalysisResultsGraphViewController::SetGraphMode)
      .def("GetGirder", &CPyAnalysisResultsGraphViewController::GetGirder, return_internal_reference<>())
      .def("SelectGirder", &CPyAnalysisResultsGraphViewController::SelectGirder)
      .add_property("ResultsType", &CPyAnalysisResultsGraphViewController::GetResultsType, &CPyAnalysisResultsGraphViewController::SetResultsType)
      .def("GetActionTypes", &CPyAnalysisResultsGraphViewController::GetActionTypes)
      .def("GetActionName", &CPyAnalysisResultsGraphViewController::GetActionName)
      .add_property("ActionType", &CPyAnalysisResultsGraphViewController::GetActionType, &CPyAnalysisResultsGraphViewController::SetActionType)
      .add_property("AnalysisType", &CPyAnalysisResultsGraphViewController::GetAnalysisType, &CPyAnalysisResultsGraphViewController::SetAnalysisType)
      .def("GetGraphTypeCount",&CPyAnalysisResultsGraphViewController::GetGraphTypeCount)
      .def("GetGraphType",&CPyAnalysisResultsGraphViewController::GetGraphType)
      .def("SelectGraphTypeByIndex",&CPyAnalysisResultsGraphViewController::SelectGraphTypeByIndex)
      .def("SelectGraphTypeByName",&CPyAnalysisResultsGraphViewController::SelectGraphTypeByName)
      .def("GetGraphCount",&CPyAnalysisResultsGraphViewController::GetGraphCount)
      .def("GetSelectedGraphCount",&CPyAnalysisResultsGraphViewController::GetSelectedGraphCount)
      .def("GetSelectedGraphs",&CPyAnalysisResultsGraphViewController::GetSelectedGraphs)
      .def("GetGraphName",&CPyAnalysisResultsGraphViewController::GetGraphName)
      .def("SelectGraphByIndex", &CPyAnalysisResultsGraphViewController::SelectGraphByIndex)
      .def("SelectGraphByName", &CPyAnalysisResultsGraphViewController::SelectGraphByName)
      .def("SelectGraphsByIndex", &CPyAnalysisResultsGraphViewController::SelectGraphsByIndex)
      .def("SelectGraphsByName", &CPyAnalysisResultsGraphViewController::SelectGraphsByName)
      .def("SelectStressLocation", &CPyAnalysisResultsGraphViewController::SelectStressLocation)
      .def("IsStressLocationSelected", &CPyAnalysisResultsGraphViewController::IsStressLocationSelected)
      .add_property("ShowGrid", &CPyAnalysisResultsGraphViewController::GetShowGrid, &CPyAnalysisResultsGraphViewController::SetShowGrid)
      .add_property("ShowGirder", &CPyAnalysisResultsGraphViewController::GetShowGirder, &CPyAnalysisResultsGraphViewController::SetShowGirder)
      ;


   /////////////////////////////////////////////////////
   // girder properties graph
   /////////////////////////////////////////////////////
   enum_<CGirderPropertiesGraphBuilder::PropertyType>("GirderPropertyType")
      .value("Height", CGirderPropertiesGraphBuilder::Height)
      .value("Area", CGirderPropertiesGraphBuilder::Area)
      .value("MomentOfInertia", CGirderPropertiesGraphBuilder::MomentOfInertia)
      .value("Centroid", CGirderPropertiesGraphBuilder::Centroid)
      .value("SectionModulus", CGirderPropertiesGraphBuilder::SectionModulus)
      .value("KernPoint", CGirderPropertiesGraphBuilder::KernPoint)
      .value("StrandEccentricity", CGirderPropertiesGraphBuilder::StrandEccentricity)
      .value("TendonEccentricity", CGirderPropertiesGraphBuilder::TendonEccentricity) // for PGSplice
      .value("TendonProfile", CGirderPropertiesGraphBuilder::TendonProfile) // for PGSplice
      .value("EffectiveFlangeWidth", CGirderPropertiesGraphBuilder::EffectiveFlangeWidth)
      .value("Fc", CGirderPropertiesGraphBuilder::Fc)
      .value("Ec", CGirderPropertiesGraphBuilder::Ec)
      ;


   class_<CPyGirderPropertiesGraphViewController, bases<CPyViewControllerBase>>("GirderPropertiesGraph")
      .def("GetGraphController", &CPyGirderPropertiesGraphViewController::GetGraphController, return_value_policy<reference_existing_object>())
      .staticmethod("GetGraphController")
      .add_property("GirderProperty",&CPyGirderPropertiesGraphViewController::GetPropertyType,&CPyGirderPropertiesGraphViewController::SetPropertyType)
      .def("IsInvariantProperty", &CPyGirderPropertiesGraphViewController::IsInvariantProperty)
      .add_property("SectionProperty",&CPyGirderPropertiesGraphViewController::GetSectionPropertyType,&CPyGirderPropertiesGraphViewController::SetSectionPropertyType)
      .def("GetGirder", &CPyGirderPropertiesGraphViewController::GetGirder, return_internal_reference<>())
      .def("SelectGirder", &CPyGirderPropertiesGraphViewController::SelectGirder)
      .add_property("Interval", &CPyGirderPropertiesGraphViewController::GetInterval,&CPyGirderPropertiesGraphViewController::SetInterval)
      .def("GetFirstInterval", &CPyGirderPropertiesGraphViewController::GetFirstInterval)
      .def("GetLastInterval", &CPyGirderPropertiesGraphViewController::GetLastInterval)
      .add_property("ShowGrid", &CPyGirderPropertiesGraphViewController::GetShowGrid, &CPyGirderPropertiesGraphViewController::SetShowGrid)
      .add_property("ShowGirder", &CPyGirderPropertiesGraphViewController::GetShowGirder, &CPyGirderPropertiesGraphViewController::SetShowGirder)
      ;
}