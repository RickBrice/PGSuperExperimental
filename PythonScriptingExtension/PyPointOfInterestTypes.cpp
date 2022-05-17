#include "stdafx.h"
#include "PyBridgeTypes.h"

#include "PyPointOfInterest.h"

using namespace boost::python;

void export_point_of_interest_types()
{
   class_<CPyPointOfInterest>("PointOfInterest")
      .def("GetPointsOfInterest", &CPyPointOfInterest::GetPointsOfInterest)
      .def("GetSpanPointsOfInterest", &CPyPointOfInterest::GetSpanPointsOfInterest)
      ;
}