#include "stdafx.h"
#include "PyAnalysisResultsGraphViewController.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

void CPyAnalysisResultsGraphViewController::Init(IAnalysisResultsGraphViewController* pViewController)
{
   m_pMyViewController = pViewController;
   CPyViewControllerBase::Init(pViewController);
}

void CPyAnalysisResultsGraphViewController::SetGraphMode(IAnalysisResultsGraphViewController::GraphMode mode)
{
   m_pMyViewController->SetGraphMode(mode);
}

IAnalysisResultsGraphViewController::GraphMode CPyAnalysisResultsGraphViewController::GetGraphMode() const
{
   return m_pMyViewController->GetGraphMode();
}

const CGirderKey& CPyAnalysisResultsGraphViewController::GetGirder() const
{
   return m_pMyViewController->GetGirder();
}

void CPyAnalysisResultsGraphViewController::SelectGirder(const CGirderKey& girderKey)
{
   m_pMyViewController->SelectGirder(girderKey);
}

void CPyAnalysisResultsGraphViewController::SetResultsType(ResultsType resultsType)
{
   m_pMyViewController->SetResultsType(resultsType);
}

ResultsType CPyAnalysisResultsGraphViewController::GetResultsType() const
{
   return m_pMyViewController->GetResultsType();
}

std::vector<ActionType> CPyAnalysisResultsGraphViewController::GetActionTypes() const
{
   return m_pMyViewController->GetActionTypes();
}

boost::python::str CPyAnalysisResultsGraphViewController::GetActionName(ActionType action) const
{
   USES_CONVERSION;
   const char* pstrName = T2A(m_pMyViewController->GetActionName(action));
   boost::python::str strName(pstrName);
   return strName;
}

void CPyAnalysisResultsGraphViewController::SetActionType(ActionType actionType)
{
   m_pMyViewController->SetActionType(actionType);
}

ActionType CPyAnalysisResultsGraphViewController::GetActionType() const
{
   return m_pMyViewController->GetActionType();
}

void CPyAnalysisResultsGraphViewController::SetAnalysisType(pgsTypes::AnalysisType analysisType)
{
   m_pMyViewController->SetAnalysisType(analysisType);
}

pgsTypes::AnalysisType CPyAnalysisResultsGraphViewController::GetAnalysisType() const
{
   return m_pMyViewController->GetAnalysisType();
}

void CPyAnalysisResultsGraphViewController::SelectStressLocation(pgsTypes::StressLocation location, bool bSelect)
{
   m_pMyViewController->SelectStressLocation(location, bSelect);
}

bool CPyAnalysisResultsGraphViewController::IsStressLocationSelected(pgsTypes::StressLocation location) const
{
   return m_pMyViewController->IsStressLocationSelected(location);
}

IndexType CPyAnalysisResultsGraphViewController::GetGraphTypeCount() const
{
   return m_pMyViewController->GetGraphTypeCount();
}

boost::python::str CPyAnalysisResultsGraphViewController::GetGraphType(IndexType idx) const
{
   USES_CONVERSION;
   const char* pstrName = T2A(m_pMyViewController->GetGraphType(idx));
   boost::python::str strName(pstrName);
   return strName;
}

void CPyAnalysisResultsGraphViewController::SelectGraphTypeByIndex(IndexType idx)
{
   m_pMyViewController->SelectGraphType(idx);
}

void CPyAnalysisResultsGraphViewController::SelectGraphTypeByName(LPCSTR lpszType)
{
   USES_CONVERSION;
   m_pMyViewController->SelectGraphType(A2T(lpszType));
}

IndexType CPyAnalysisResultsGraphViewController::GetGraphCount() const
{
   return m_pMyViewController->GetGraphCount();
}

IndexType CPyAnalysisResultsGraphViewController::GetSelectedGraphCount() const
{
   return m_pMyViewController->GetSelectedGraphCount();
}

boost::python::list CPyAnalysisResultsGraphViewController::GetSelectedGraphs() const
{
   boost::python::list list;
   std::vector<IndexType> vSelectedGraphs = m_pMyViewController->GetSelectedGraphs();
   std::for_each(vSelectedGraphs.cbegin(), vSelectedGraphs.cend(), [&list](const auto& index) {list.append(index);});
   return list;
}

boost::python::str CPyAnalysisResultsGraphViewController::GetGraphName(IndexType graphIdx) const
{
   USES_CONVERSION;
   const char* pstrName = T2A(m_pMyViewController->GetGraphName(graphIdx));
   boost::python::str strGraphName(pstrName);
   return strGraphName;
}

void CPyAnalysisResultsGraphViewController::SelectGraphByIndex(IndexType graphIdx)
{
   m_pMyViewController->SelectGraph(graphIdx);
}

void CPyAnalysisResultsGraphViewController::SelectGraphByName(LPCSTR lpszGraphName)
{
   USES_CONVERSION;
   m_pMyViewController->SelectGraph(A2T(lpszGraphName));
}

void CPyAnalysisResultsGraphViewController::SelectGraphsByIndex(const boost::python::list& lGraphs)
{
   boost::python::stl_input_iterator<IndexType> begin(lGraphs), end;
   std::vector<IndexType> vGraphs(begin, end);
   m_pMyViewController->SelectGraphs(vGraphs);
}

void CPyAnalysisResultsGraphViewController::SelectGraphsByName(const boost::python::list& lGraphs)
{
   boost::python::stl_input_iterator<const char*> begin(lGraphs), end;
   std::vector<CString> vGraphs(begin, end);
   m_pMyViewController->SelectGraphs(vGraphs);
}

void CPyAnalysisResultsGraphViewController::SetIncludeElevationAdjustment(bool bInclude)
{
   m_pMyViewController->IncludeElevationAdjustment(bInclude);
}

bool CPyAnalysisResultsGraphViewController::GetIncludeElevationAdjustment() const
{
   return m_pMyViewController->IncludeElevationAdjustment();
}

void CPyAnalysisResultsGraphViewController::SetIncludeUnrecoverableDefl(bool bInclude)
{
   m_pMyViewController->IncludeUnrecoverableDefl(bInclude);
}

bool CPyAnalysisResultsGraphViewController::GetIncludeUnrecoverableDefl() const
{
   return m_pMyViewController->IncludeUnrecoverableDefl();
}

void CPyAnalysisResultsGraphViewController::SetShowGrid(bool bShow)
{
   m_pMyViewController->ShowGrid(bShow);
}

bool CPyAnalysisResultsGraphViewController::GetShowGrid() const
{
   return m_pMyViewController->ShowGrid();
}

void CPyAnalysisResultsGraphViewController::SetShowGirder(bool bShow)
{
   m_pMyViewController->ShowGirder(bShow);
}

bool CPyAnalysisResultsGraphViewController::GetShowGirder() const
{
   return m_pMyViewController->ShowGirder();
}

CPyAnalysisResultsGraphViewController* CPyAnalysisResultsGraphViewController::GetGraphController(CPyViewControllerBase* pController)
{
   auto pMyController = static_cast<CPyAnalysisResultsGraphViewController*>(pController);
   return pMyController;
}
