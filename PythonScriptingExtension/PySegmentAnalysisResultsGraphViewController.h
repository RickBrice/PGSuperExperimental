#pragma once

#include "PyViewControllerBase.h"
#include <EAF\EAFViewController.h>
#include <Graphs\SegmentAnalysisResultsGraphViewController.h>

class CPySegmentAnalysisResultsGraphViewController : public CPyViewControllerBase
{
public:
   void Init(ISegmentAnalysisResultsGraphViewController* pViewController);

   void SetGraphMode(ISegmentAnalysisResultsGraphViewController::GraphMode mode);
   ISegmentAnalysisResultsGraphViewController::GraphMode GetGraphMode() const;

   CSegmentKey GetSegment() const;
   void SelectSegment(const CSegmentKey& segment);

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

   void SetShowGrid(bool bShow);
   bool GetShowGrid() const;

   void SetShowGirder(bool bShow);
   bool GetShowGirder() const;

   void SelectStressLocation(pgsTypes::StressLocation location, bool bSelect);
   bool IsStressLocationSelected(pgsTypes::StressLocation location) const;

   static CPySegmentAnalysisResultsGraphViewController* GetGraphController(CPyViewControllerBase* pController);

private:
   CComPtr<ISegmentAnalysisResultsGraphViewController> m_pMyViewController;
};

