///////////////////////////////////////////////////////////////////////
// IEPluginExample
// Copyright © 1999-2020  Washington State Department of Transportation
//                        Bridge and Structures Office
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the Alternate Route Open Source License as 
// published by the Washington State Department of Transportation, 
// Bridge and Structures Office.
//
// This program is distributed in the hope that it will be useful, but 
// distribution is AS IS, WITHOUT ANY WARRANTY; without even the implied 
// warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See 
// the Alternate Route Open Source License for more details.
//
// You should have received a copy of the Alternate Route Open Source 
// License along with this program; if not, write to the Washington 
// State Department of Transportation, Bridge and Structures Office, 
// P.O. Box  47340, Olympia, WA 98503, USA or e-mail 
// Bridge_Support@wsdot.wa.gov
///////////////////////////////////////////////////////////////////////

// PGSuperDataImporter.cpp : Implementation of CPGSuperDataImporter
#include "stdafx.h"
#include "IFCExporter.h"
#include "PGSuperDataImporter.h"
#include <EAF\EAFAutoProgress.h>
#include <IFace\Project.h>

/////////////////////////////////////////////////////////////////////////////
// CPGSuperDataImporter
STDMETHODIMP CPGSuperDataImporter::Init(UINT nCmdID)
{
   return S_OK;
}

STDMETHODIMP CPGSuperDataImporter::GetMenuText(BSTR*  bstrText) const
{
   *bstrText = CComBSTR("Alignment from IFC File");
   return S_OK;
}

STDMETHODIMP CPGSuperDataImporter::GetBitmapHandle(HBITMAP* phBmp) const
{
   *phBmp = nullptr;
   return S_OK;
}

STDMETHODIMP CPGSuperDataImporter::GetCommandHintText(BSTR*  bstrText) const
{
   *bstrText = CComBSTR("Status line hint text\nTool tip text");
   return S_OK;   
}

STDMETHODIMP CPGSuperDataImporter::Import(IBroker* pBroker)
{
   AFX_MANAGE_STATE(AfxGetStaticModuleState());
   USES_CONVERSION;

   CFileDialog dlg(TRUE, _T("ifc"),NULL,OFN_FILEMUSTEXIST | OFN_PATHMUSTEXIST,_T("IFC Files (*.ifc)|*.ifc||"));
   if (dlg.DoModal() == IDOK)
   {
      CString fileName = dlg.GetPathName();

      IfcParse::IfcFile file(T2A(fileName.GetBuffer()));

#pragma Reminder("WORKING HERE - UNITS - THERE ARE MANY CASES THIS DOESN'T DEAL WITH")
      auto unit_assignment_instances = file.instances_by_type<Ifc4x1::IfcUnitAssignment>();
      ATLASSERT(unit_assignment_instances->size() == 1);
      auto unit_assignment = *(unit_assignment_instances->begin());
      auto units = unit_assignment->Units();
      for (auto unit : *units)
      {
         auto derived_unit = unit->as<Ifc4x1::IfcDerivedUnit>();
         auto monitary_unit = unit->as<Ifc4x1::IfcMonetaryUnit>();
         auto si_unit = unit->as<Ifc4x1::IfcSIUnit>();
         auto conversion_based_unit = unit->as<Ifc4x1::IfcConversionBasedUnit>();
         auto conversion_based_unit_with_offset = unit->as<Ifc4x1::IfcConversionBasedUnitWithOffset>();

         if (si_unit && si_unit->Name() == Ifc4x1::IfcSIUnitName::IfcSIUnitName_METRE)
         {
            if (si_unit->hasPrefix())
            {
               switch (si_unit->Prefix())
               {
               case Ifc4x1::IfcSIPrefix::IfcSIPrefix_KILO:
                  m_IfcConverter.SetLengthUnit(&unitMeasure::Kilometer);
                  break;

               case Ifc4x1::IfcSIPrefix::IfcSIPrefix_CENTI:
                  m_IfcConverter.SetLengthUnit(&unitMeasure::Centimeter);
                  break;

               case Ifc4x1::IfcSIPrefix::IfcSIPrefix_MILLI:
                  m_IfcConverter.SetLengthUnit(&unitMeasure::Millimeter);
                  break;

               default:
                  ATLASSERT(false); // unit prefix isn't supported
               }
            }
            else
            {
               m_IfcConverter.SetLengthUnit(&unitMeasure::Meter);
            }
            break;
         }

         if (conversion_based_unit && conversion_based_unit->UnitType() == Ifc4x1::IfcUnitEnum::IfcUnit_LENGTHUNIT)
         {
            if (conversion_based_unit->UnitType() == Ifc4x1::IfcUnitEnum::IfcUnit_LENGTHUNIT)
            {
               auto measure_with_unit = conversion_based_unit->ConversionFactor();
               auto unit_component = measure_with_unit->UnitComponent()->as<Ifc4x1::IfcSIUnit>();
               ATLASSERT(unit_component->Name() == Ifc4x1::IfcSIUnitName::IfcSIUnitName_METRE);
               ATLASSERT(unit_component->hasPrefix() == false); // not dealing with conversion factors to anything but meter
               auto value_component = measure_with_unit->ValueComponent()->as<Ifc4x1::IfcReal>();
               ATLASSERT(value_component); // not dealing with anything but simple conversion factors
               if (IsEqual((Float64)(*value_component), unitMeasure::Feet.GetConvFactor()))
               {
                  m_IfcConverter.SetLengthUnit(&unitMeasure::Feet);
               }
               else if (IsEqual((Float64)(*value_component), unitMeasure::Inch.GetConvFactor()))
               {
                  m_IfcConverter.SetLengthUnit(&unitMeasure::Inch);
               }
               else if (IsEqual((Float64)(*value_component), unitMeasure::Mile.GetConvFactor()))
               {
                  m_IfcConverter.SetLengthUnit(&unitMeasure::Mile);
               }

            }
            break;
         }
      }

      if (!file.good())
      {
         AfxMessageBox(_T("Unable to parse .ifc file"));
         return S_OK;
      }

      auto alignment = m_IfcConverter.GetAlignment(file);
      if (alignment)
      {
         AlignmentData2 alignment_data;
         ProfileData2 profile_data;
         RoadwaySectionData section_data;
         m_IfcConverter.ConvertToPGSuper(alignment, &alignment_data,&profile_data,&section_data);

         GET_IFACE2(pBroker, IEvents, pEvents);
         pEvents->HoldEvents();

         GET_IFACE2(pBroker, IRoadwayData, pRoadwayData);
         pRoadwayData->SetAlignmentData2(alignment_data);
         pRoadwayData->SetProfileData2(profile_data);
         pRoadwayData->SetRoadwaySectionData(section_data);

         pEvents->FirePendingEvents();
      }
   }

   return S_OK;
}

