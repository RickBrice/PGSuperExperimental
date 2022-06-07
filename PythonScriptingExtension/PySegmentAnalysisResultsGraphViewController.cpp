#include "stdafx.h"
#include "PySegmentAnalysisResultsGraphViewController.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

void CPySegmentAnalysisResultsGraphViewController::Init(ISegmentAnalysisResultsGraphViewController* pViewController)
{
   m_pMyViewController = pViewController;
   CPyViewControllerBase::Init(pViewController);
}

void CPySegmentAnalysisResultsGraphViewController::SetGraphMode(ISegmentAnalysisResultsGraphViewController::GraphMode mode)
{
   m_pMyViewController->SetGraphMode(mode);
}

ISegmentAnalysisResultsGraphViewController::GraphMode CPySegmentAnalysisResultsGraphViewController::GetGraphMode() const
{
   return m_pMyViewController->GetGraphMode();
}

CSegmentKey CPySegmentAnalysisResultsGraphViewController::GetSegment() const
{
   return m_pMyViewController->GetSegment();
}

void CPySegmentAnalysisResultsGraphViewController::SelectSegment(const CSegmentKey& segmentKey)
{
   m_pMyViewController->SelectSegment(segmentKey);
}

void CPySegmentAnalysisResultsGraphViewController::SetResultsType(ResultsType resultsType)
{
   m_pMyViewController->SetResultsType(resultsType);
}

ResultsType CPySegmentAnalysisResultsGraphViewController::GetResultsType() const
{
   return m_pMyViewController->GetResultsType();
}

std::vector<ActionType> CPySegmentAnalysisResultsGraphViewController::GetActionTypes() const
{
   return m_pMyViewController->GetActionTypes();
}

boost::python::str CPySegmentAnalysisResultsGraphViewController::GetActionName(ActionType action) const
{
   USES_CONVERSION;
   const char* pstrName = T2A(m_pMyViewController->GetActionName(action));
   boost::python::str strName(pstrName);
   return strName;
}

void CPySegmentAnalysisResultsGraphViewController::SetActionType(ActionType actionType)
{
   m_pMyViewController->SetActionType(actionType);
}

ActionType CPySegmentAnalysisResultsGraphViewController::GetActionType() const
{
   return m_pMyViewController->GetActionType();
}

void CPySegmentAnalysisResultsGraphViewController::SetAnalysisType(pgsTypes::AnalysisType analysisType)
{
   m_pMyViewController->SetAnalysisType(analysisType);
}

pgsTypes::AnalysisType CPySegmentAnalysisResultsGraphViewController::GetAnalysisType() const
{
   return m_pMyViewController->GetAnalysisType();
}

void CPySegmentAnalysisResultsGraphViewController::SelectStressLocation(pgsTypes::StressLocation location, bool bSelect)
{
   m_pMyViewController->SelectStressLocation(location, bSelect);
}

bool CPySegmentAnalysisResultsGraphViewController::IsStressLocationSelected(pgsTypes::StressLocation location) const
{
   return m_pMyViewController->IsStressLocationSelected(location);
}

IndexType CPySegmentAnalysisResultsGraphViewController::GetGraphTypeCount() const
{
   return m_pMyViewController->GetGraphTypeCount();
}

boost::python::str CPySegmentAnalysisResultsGraphViewController::GetGraphType(IndexType idx) const
{
   USES_CONVERSION;
   const char* pstrName = T2A(m_pMyViewController->GetGraphType(idx));
   boost::python::str strName(pstrName);
   return strName;
}

void CPySegmentAnalysisResultsGraphViewController::SelectGraphTypeByIndex(IndexType idx)
{
   m_pMyViewController->SelectGraphType(idx);
}

void CPySegmentAnalysisResultsGraphViewController::SelectGraphTypeByName(LPCSTR lpszType)
{
   USES_CONVERSION;
   m_pMyViewController->SelectGraphType(A2T(lpszType));
}

IndexType CPySegmentAnalysisResultsGraphViewController::GetGraphCount() const
{
   return m_pMyViewController->GetGraphCount();
}

IndexType CPySegmentAnalysisResultsGraphViewController::GetSelectedGraphCount() const
{
   return m_pMyViewController->GetSelectedGraphCount();
}

boost::python::list CPySegmentAnalysisResultsGraphViewController::GetSelectedGraphs() const
{
   boost::python::list list;
   std::vector<IndexType> vSelectedGraphs = m_pMyViewController->GetSelectedGraphs();
   std::for_each(vSelectedGraphs.cbegin(), vSelectedGraphs.cend(), [&list](const auto& index) {list.append(index);});
   return list;
}

boost::python::str CPySegmentAnalysisResultsGraphViewController::GetGraphName(IndexType graphIdx) const
{
   USES_CONVERSION;
   const char* pstrName = T2A(m_pMyViewController->GetGraphName(graphIdx));
   boost::python::str strGraphName(pstrName);
   return strGraphName;
}

void CPySegmentAnalysisResultsGraphViewController::SelectGraphByIndex(IndexType graphIdx)
{
   m_pMyViewController->SelectGraph(graphIdx);
}

void CPySegmentAnalysisResultsGraphViewController::SelectGraphByName(LPCSTR lpszGraphName)
{
   USES_CONVERSION;
   m_pMyViewController->SelectGraph(A2T(lpszGraphName));
}

void CPySegmentAnalysisResultsGraphViewController::SelectGraphsByIndex(const boost::python::list& lGraphs)
{
   boost::python::stl_input_iterator<IndexType> begin(lGraphs), end;
   std::vector<IndexType> vGraphs(begin, end);
   m_pMyViewController->SelectGraphs(vGraphs);
}

void CPySegmentAnalysisResultsGraphViewController::SelectGraphsByName(const boost::python::list& lGraphs)
{
   boost::python::stl_input_iterator<const char*> begin(lGraphs), end;
   std::vector<CString> vGraphs(begin, end);
   m_pMyViewController->SelectGraphs(vGraphs);
}

void CPySegmentAnalysisResultsGraphViewController::SetShowGrid(bool bShow)
{
   m_pMyViewController->ShowGrid(bShow);
}

bool CPySegmentAnalysisResultsGraphViewController::GetShowGrid() const
{
   return m_pMyViewController->ShowGrid();
}

void CPySegmentAnalysisResultsGraphViewController::SetShowGirder(bool bShow)
{
   m_pMyViewController->ShowGirder(bShow);
}

bool CPySegmentAnalysisResultsGraphViewController::GetShowGirder() const
{
   return m_pMyViewController->ShowGirder();
}

CPySegmentAnalysisResultsGraphViewController* CPySegmentAnalysisResultsGraphViewController::GetGraphController(CPyViewControllerBase* pController)
{
   auto pMyController = static_cast<CPySegmentAnalysisResultsGraphViewController*>(pController);
   return pMyController;
}
