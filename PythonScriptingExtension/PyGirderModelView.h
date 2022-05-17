#pragma once

#include <GirderModelViewController.h>
#include "PyViewControllerBase.h"

class CPyGirderModelView : public CPyViewControllerBase
{
public:
   void Init(IGirderModelViewController* pViewController);

   bool GetSyncWithBridgeModelView() const;
   void SetSyncWithBridgeModelView(bool bSync);
   const CGirderKey& GetGirder() const;
   void SelectGirder(const CGirderKey& girderKey);
   EventIndexType GetEvent() const;
   bool SetEvent(EventIndexType eventIdx);
   Float64 GetCutLocation() const;
   void CutAt(Float64 Xg);
   void CutAtNext();
   void CutAtPrev();
   Float64 GetMinCutLocation() const;
   Float64 GetMaxCutLocation() const;
   void SetShowStrands(bool bShow);
   bool GetShowStrands() const;
   void SetShowStrandCG(bool bShow);
   bool GetShowStrandCG() const;
   void SetShowCG(bool bShow);
   bool GetShowCG() const;
   void SetShowSectionProperties(bool bShow);
   bool GetShowSectionProperties() const;
   void SetShowDimensions(bool bShow);
   bool GetShowDimensions() const;
   void SetShowLongitudinalReinforcement(bool bShow);
   bool GetShowLongitudinalReinforcement() const;
   void SetShowTransverseReinforcement(bool bShow);
   bool GetShowTransverseReinforcement() const;
   void SetShowLoads(bool bShow);
   bool GetShowLoads() const;
   void SetSchematic(bool bSchematic);
   bool GetSchematic() const;

private:
   CComPtr<IGirderModelViewController> m_pMyViewController;
};

