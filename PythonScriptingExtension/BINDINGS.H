#pragma once

#include "PyProjectProperties.h"
#include "PyEnvironment.h"
#include "PyRoadway.h"
#include "PyViews.h"
#include "PyGraphManager.h"
#include "PyReportManager.h"
#include "PyLoadModifiers.h"
#include "PyBridge.h"
#include "PyPointOfInterest.h"
#include "PyLimits.h"
#include "PyEffectiveFlangeWidth.h"
#include "PyLossParameters.h"
#include "PyMomentCapacity.h"

#include <IFace\VersionInfo.h>


class CPyPGSuper
{
public:
   void Init(IBroker* pBroker);
   void Reset();

   std::string GetVersionString() const;
   std::string GetVersionStringWithBuildNumber() const;

   std::string GetVersion() const;
   std::string GetVersionWithBuildNumber() const;

   CPyProjectProperties Properties();
   CPyEnvironment Environment();
   CPyRoadway Roadway();
   CPyViews Views();
   CPyGraphManager Graphs();
   CPyReportManager Reports();
   CPyLoadModifiers LoadModifiers();
   CPyBridge Bridge();
   CPyPointOfInterest PointOfInterest();
   CPyLimits Limits();
   CPyEffectiveFlangeWidth EffectiveFlangeWidth();
   CPyLossParameters LossParameters();
   CPyMomentCapacity MomentCapacity();

private:
   CComPtr<IVersionInfo> m_pVersionInfo;

   CPyProjectProperties m_Properties;
   CPyEnvironment m_Environment;
   CPyRoadway m_Roadway;
   CPyViews m_Views;
   CPyGraphManager m_GraphManager;
   CPyReportManager m_ReportManager;
   CPyLoadModifiers m_LoadModifiers;
   CPyBridge m_Bridge;
   CPyPointOfInterest m_PointOfInterest;
   CPyLimits m_Limits;
   CPyEffectiveFlangeWidth m_EffectiveFlangeWidth;
   CPyLossParameters m_LossParameters;
   CPyMomentCapacity m_MomentCapacity;
};
