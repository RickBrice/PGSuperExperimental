#pragma once

#include "PyViewControllerBase.h"
#include <EAF\EAFViewController.h>
#include <Graphs\AnalysisResultsGraphViewController.h>

class CPyAnalysisResultsGraphViewController : public CPyViewControllerBase
{
public:
   void Init(IAnalysisResultsGraphViewController* pViewController);

   void SetGraphMode(IAnalysisResultsGraphViewController::GraphMode mode);
   IAnalysisResultsGraphViewController::GraphMode GetGraphMode() const;

   const CGirderKey& GetGirder() const;
   void SelectGirder(const CGirderKey& girderKey);

   std::vector<ActionType> GetActionTypes() const;
   boost::python::str GetActionName(ActionType action) const;

   void SetResultsType(ResultsType resultsType);
   ResultsType GetResultsType() const;

   void SetActionType(ActionType actionType);
   ActionType GetActionType() const;

   void SetAnalysisType(pgsTypes::AnalysisType analysisType);
   pgsTypes::AnalysisType GetAnalysisType() const;

   IndexType GetGraphTypeCount() const;
   boost::python::str GetGraphType(IndexType idx) const;
   void SelectGraphTypeByIndex(IndexType idx);
   void SelectGraphTypeByName(LPCSTR lpszType);

   IndexType GetGraphCount() const;
   IndexType GetSelectedGraphCount() const;
   boost::python::list GetSelectedGraphs() const;
   boost::python::str GetGraphName(IndexType graphIdx) const;
   void SelectGraphByIndex(IndexType graphIdx);
   void SelectGraphByName(LPCSTR lpszGraphName);
   void SelectGraphsByIndex(const boost::python::list& lGraphs);
   void SelectGraphsByName(const boost::python::list& lGraphs);

   void SetIncludeElevationAdjustment(bool bInclude);
   bool GetIncludeElevationAdjustment() const;

   void SetIncludeUnrecoverableDefl(bool bInclude);
   bool GetIncludeUnrecoverableDefl() const;


   void SetShowGrid(bool bShow);
   bool GetShowGrid() const;

   void SetShowGirder(bool bShow);
   bool GetShowGirder() const;

   void SelectStressLocation(pgsTypes::StressLocation location, bool bSelect);
   bool IsStressLocationSelected(pgsTypes::StressLocation location) const;

   static CPyAnalysisResultsGraphViewController* GetGraphController(CPyViewControllerBase* pController);

private:
   CComPtr<IAnalysisResultsGraphViewController> m_pMyViewController;
};

