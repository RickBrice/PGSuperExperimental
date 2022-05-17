#include "stdafx.h"
#include "PyViewTypes.h"

#include "Bindings.h"

using namespace boost::python;

void export_view_types()
{
   // view creation methods
   class_<CPyViews>("Views", "Manipulate views")
      .def("CreateBridgeModelView", &CPyViews::CreateBridgeModelView)
      .def("CreateGirderView", &CPyViews::CreateGirderView)
      .def("CreateLoadsView", &CPyViews::CreateLoadsView)
      .def("CreateLibraryEditorView", &CPyViews::CreateLibraryEditorView)
      .def("CreateGraphView", &CPyViews::CreateGraphView)
      .def("CreateGraphViewByName", &CPyViews::CreateGraphViewByName, return_value_policy<manage_new_object>())
      .def("CreateReport", &CPyViews::CreateReport)
      ;

   // base class controller for all views
   class_<CPyViewControllerBase>("ViewController")
      .def("Close", &CPyViewControllerBase::Close)
      .def("IsOpen", &CPyViewControllerBase::IsOpen)
      .def("Minimize", &CPyViewControllerBase::Minimize)
      .def("Maximize", &CPyViewControllerBase::Maximize)
      .def("Restore", &CPyViewControllerBase::Restore)
      ;

   // bridge model view
   enum_<IBridgeModelViewController::ViewMode>("BridgeModelViewMode")
      .value("Bridge", IBridgeModelViewController::Bridge)
      .value("Alignment", IBridgeModelViewController::Alignment)
      ;

   class_<CPyBridgeModelView, bases<CPyViewControllerBase>>("BridgeModelView")
      .def("SetGroupRange", &CPyBridgeModelView::SetGroupRange)
      .def("GetGroupRange", &CPyBridgeModelView::GetGroupRange)
      .add_property("CutStation", &CPyBridgeModelView::GetCutStation, &CPyBridgeModelView::SetCutStation)
      .add_property("ViewMode", &CPyBridgeModelView::GetViewMode, &CPyBridgeModelView::SetViewMode)
      .add_property("NorthUp", &CPyBridgeModelView::GetNorthUp, &CPyBridgeModelView::SetNorthUp)
      .add_property("ShowLabels", &CPyBridgeModelView::GetShowLabels, &CPyBridgeModelView::SetShowLabels)
      .add_property("ShowDimensions", &CPyBridgeModelView::GetShowDimensions, &CPyBridgeModelView::SetShowDimensions)
      .add_property("ShowBridge", &CPyBridgeModelView::GetShowBridge, &CPyBridgeModelView::SetShowBridge)
      .add_property("Schematic", &CPyBridgeModelView::GetSchematic, &CPyBridgeModelView::SetSchematic)

      .def("ClearSelection", &CPyBridgeModelView::ClearSelection)

      .def("GetSelectPier", &CPyBridgeModelView::GetSelectedPier)
      .def("SelectPier", &CPyBridgeModelView::SelectPier)

      .def("GetSelectSpan", &CPyBridgeModelView::GetSelectedSpan)
      .def("SelectSpan", &CPyBridgeModelView::SelectSpan)

      .def("GetSelectGirder", &CPyBridgeModelView::GetSelectedGirder)
      .def("SelectGirder", &CPyBridgeModelView::SelectGirder)

      .def("GetSelectSegment", &CPyBridgeModelView::GetSelectedSegment)
      .def("SelectSegment", &CPyBridgeModelView::SelectSegment)

      .def("IsDeckSelected", &CPyBridgeModelView::IsDeckSelected)
      .def("SelectDeck", &CPyBridgeModelView::SelectDeck)

      .def("IsAlignmentSelected", &CPyBridgeModelView::IsAlignmentSelected)
      .def("SelectAlignment", &CPyBridgeModelView::SelectAlignment)

      .def("IsRailingSystemSelected", &CPyBridgeModelView::IsRailingSystemSelected)
      .def("SelectRailingSystem", &CPyBridgeModelView::SelectRailingSystem)
      ;

   // girder model view
   class_<CPyGirderModelView, bases<CPyViewControllerBase>>("GirderModelView")
      .add_property("SyncWithBridgeModelView", &CPyGirderModelView::GetSyncWithBridgeModelView, &CPyGirderModelView::SetSyncWithBridgeModelView)
      .def("GetGirder", &CPyGirderModelView::GetGirder, return_internal_reference<>())
      .def("SelectGirder", &CPyGirderModelView::SelectGirder)
      .add_property("Event", &CPyGirderModelView::GetEvent, &CPyGirderModelView::SetEvent)
      .def("GetCutLocation", &CPyGirderModelView::GetCutLocation)
      .def("CutAt", &CPyGirderModelView::CutAt)
      .def("CutAtNext", &CPyGirderModelView::CutAtNext)
      .def("CutAtPrev", &CPyGirderModelView::CutAtPrev)
      .def("GetMinCutLocation", &CPyGirderModelView::GetMinCutLocation)
      .def("GetMaxCutLocation", &CPyGirderModelView::GetMaxCutLocation)
      .add_property("ShowStrands", &CPyGirderModelView::GetShowStrands, &CPyGirderModelView::SetShowStrands)
      .add_property("ShowStrandCG", &CPyGirderModelView::GetShowStrandCG, &CPyGirderModelView::SetShowStrandCG)
      .add_property("ShowCG", &CPyGirderModelView::GetShowCG, &CPyGirderModelView::SetShowCG)
      .add_property("ShowSectionProperties", &CPyGirderModelView::GetShowSectionProperties, &CPyGirderModelView::SetShowSectionProperties)
      .add_property("ShowDimensions", &CPyGirderModelView::GetShowDimensions, &CPyGirderModelView::SetShowDimensions)
      .add_property("ShowLongitindalReinforcement", &CPyGirderModelView::GetShowLongitudinalReinforcement, &CPyGirderModelView::SetShowLongitudinalReinforcement)
      .add_property("ShowTransverseReinforcement", &CPyGirderModelView::GetShowTransverseReinforcement, &CPyGirderModelView::SetShowTransverseReinforcement)
      .add_property("ShowLoads", &CPyGirderModelView::GetShowLoads, &CPyGirderModelView::SetShowLoads)
      .add_property("Schematic", &CPyGirderModelView::GetSchematic, &CPyGirderModelView::SetSchematic)
      ;

   // loads view
   enum_<ILoadsViewController::Field>("LoadsViewField")
      .value("Type", ILoadsViewController::Type)
      .value("Event", ILoadsViewController::Event)
      .value("LoadCase", ILoadsViewController::LoadCase)
      .value("Location", ILoadsViewController::Location)
      .value("Magnitude", ILoadsViewController::Magnitude)
      .value("Description", ILoadsViewController::Description)
      ;

   enum_<ILoadsViewController::Direction>("LoadsViewSortDirection")
      .value("Ascending", ILoadsViewController::Ascending)
      .value("Descending", ILoadsViewController::Descending)
      ;

   class_<CPyLoadsView, bases<CPyViewControllerBase>>("LoadsView")
      .def("SortBy", &CPyLoadsView::SortBy)
      .add_property("LoadCount", &CPyLoadsView::GetLoadCount)
      .add_property("FieldValue", &CPyLoadsView::GetFieldValue)
      .def("DeleteLoad", &CPyLoadsView::DeleteLoad)
      ;

}
