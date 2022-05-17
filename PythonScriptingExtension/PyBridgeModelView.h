#pragma once

#include <BridgeModelViewController.h>
#include <IFace\Selection.h>
#include "PyViewControllerBase.h"

class CPyBridgeModelView : public CPyViewControllerBase
{
public:
   void Init(IBridgeModelViewController* pViewController,ISelection* pSelection);

   void GetGroupRange(GroupIndexType* pStartGroupIdx, GroupIndexType* pEndGroupIdx) const;
   void SetGroupRange(GroupIndexType startGroupIdx, GroupIndexType endGroupIdx);
   Float64 GetCutStation() const;
   void SetCutStation(Float64 station);
   void SetViewMode(IBridgeModelViewController::ViewMode mode);
   IBridgeModelViewController::ViewMode GetViewMode() const;

   void SetNorthUp(bool bNorthUp);
   bool GetNorthUp() const;
   void SetShowLabels(bool bShowLabels);
   bool GetShowLabels() const;
   void SetShowDimensions(bool bShowDimensions);
   bool GetShowDimensions() const;
   void SetShowBridge(bool bShowBridge);
   bool GetShowBridge() const;
   void SetSchematic(bool bSchematic);
   bool GetSchematic() const;

   void ClearSelection();

   PierIndexType GetSelectedPier();
   SpanIndexType GetSelectedSpan();
   CGirderKey GetSelectedGirder();
   CSegmentKey GetSelectedSegment();
   CClosureKey GetSelectedClosureJoint();
   SupportIDType GetSelectedTemporarySupport();
   bool IsDeckSelected();
   bool IsAlignmentSelected();
   bool IsRailingSystemSelected(pgsTypes::TrafficBarrierOrientation orientation);

   void SelectPier(PierIndexType pierIdx);
   void SelectSpan(SpanIndexType spanIdx);
   void SelectGirder(const CGirderKey& girderKey);
   void SelectSegment(const CSegmentKey& segmentKey);
   void SelectClosureJoint(const CClosureKey& closureKey);
   void SelectTemporarySupport(SupportIDType tsID);
   void SelectDeck();
   void SelectAlignment();
   void SelectRailingSystem(pgsTypes::TrafficBarrierOrientation orientation);

private:
   CComPtr<IBridgeModelViewController> m_pViewController;
   CComPtr<ISelection> m_pSelection;
};

