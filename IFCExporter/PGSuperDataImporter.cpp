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
         m_IfcConverter.ConvertToPGSuper(alignment, &alignment_data,&profile_data);

         GET_IFACE2(pBroker, IEvents, pEvents);
         pEvents->HoldEvents();

         GET_IFACE2(pBroker, IRoadwayData, pRoadwayData);
         pRoadwayData->SetAlignmentData2(alignment_data);
         pRoadwayData->SetProfileData2(profile_data);

         pEvents->FirePendingEvents();
      }
   }

   return S_OK;
}

