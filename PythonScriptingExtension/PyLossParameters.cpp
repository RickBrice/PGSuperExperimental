#include "stdafx.h"
#include "PyLossParameters.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif

void CPyLossParameters::Init(IBroker* pBroker)
{
   pBroker->GetInterface(IID_ILossParameters, (IUnknown**)&m_pLossParameters);
}

void CPyLossParameters::Reset()
{
   m_pLossParameters.Release();
}

pgsTypes::LossMethod CPyLossParameters::GetLossMethod() const
{
   return m_pLossParameters->GetLossMethod();
}

pgsTypes::TimeDependentModel CPyLossParameters::GetTimeDependentModel() const
{
   return m_pLossParameters->GetTimeDependentModel();
}

void CPyLossParameters::SetIgnoreCreepEffects(bool bIgnore)
{
   m_pLossParameters->IgnoreCreepEffects(bIgnore);
}

bool CPyLossParameters::GetIgnoreCreepEffects() const
{
   return m_pLossParameters->IgnoreCreepEffects();
}

void CPyLossParameters::SetIgnoreShrinkageEffects(bool bIgnore)
{
   m_pLossParameters->IgnoreShrinkageEffects(bIgnore);
}

bool CPyLossParameters::GetIgnoreShrinkageEffects() const
{
   return m_pLossParameters->IgnoreShrinkageEffects();
}

void CPyLossParameters::SetIgnoreRelaxationEffects(bool bIgnore)
{
   m_pLossParameters->IgnoreRelaxationEffects(bIgnore);
}

bool CPyLossParameters::GetIgnoreRelaxationEffects() const
{
   return m_pLossParameters->IgnoreRelaxationEffects();
}

void CPyLossParameters::SetIgnoreTimeDependentEffects(bool bIgnoreCreep, bool bIgnoreShrinkage, bool bIgnoreRelaxation)
{
   m_pLossParameters->IgnoreTimeDependentEffects(bIgnoreCreep, bIgnoreShrinkage, bIgnoreRelaxation);
}

void CPyLossParameters::SetTendonPostTensionParameters(Float64 Dset, Float64 wobble, Float64 friction)
{
   m_pLossParameters->SetTendonPostTensionParameters(Dset, wobble, friction);
}

void CPyLossParameters::GetTendonPostTensionParameters(Float64* Dset, Float64* wobble, Float64* friction) const
{
   m_pLossParameters->GetTendonPostTensionParameters(Dset, wobble, friction);
}

void CPyLossParameters::SetTemporaryStrandPostTensionParameters(Float64 Dset, Float64 wobble, Float64 friction)
{
   m_pLossParameters->SetTemporaryStrandPostTensionParameters(Dset, wobble, friction);
}

void CPyLossParameters::GetTemporaryStrandPostTensionParameters(Float64* Dset, Float64* wobble, Float64* friction) const
{
   m_pLossParameters->GetTemporaryStrandPostTensionParameters(Dset, wobble, friction);
}

void CPyLossParameters::SetUseGeneralLumpSumLosses(bool bLumpSum)
{
   m_pLossParameters->UseGeneralLumpSumLosses(bLumpSum);
}

bool CPyLossParameters::GetUseGeneralLumpSumLosses() const
{
   return m_pLossParameters->UseGeneralLumpSumLosses();
}

Float64 CPyLossParameters::GetBeforeXferLosses() const
{
   return m_pLossParameters->GetBeforeXferLosses();
}

void CPyLossParameters::SetBeforeXferLosses(Float64 loss)
{
   m_pLossParameters->SetBeforeXferLosses(loss);
}

Float64 CPyLossParameters::GetAfterXferLosses() const
{
   return m_pLossParameters->GetAfterXferLosses();
}

void CPyLossParameters::SetAfterXferLosses(Float64 loss)
{
   m_pLossParameters->SetAfterXferLosses(loss);
}

Float64 CPyLossParameters::GetLiftingLosses() const
{
   return m_pLossParameters->GetLiftingLosses();
}

void CPyLossParameters::SetLiftingLosses(Float64 loss)
{
   m_pLossParameters->SetLiftingLosses(loss);
}

Float64 CPyLossParameters::GetShippingLosses() const
{
   return m_pLossParameters->GetShippingLosses();
}

void CPyLossParameters::SetShippingLosses(Float64 loss)
{
   m_pLossParameters->SetShippingLosses(loss);
}

Float64 CPyLossParameters::GetBeforeTempStrandRemovalLosses() const
{
   return m_pLossParameters->GetBeforeTempStrandRemovalLosses();
}

void CPyLossParameters::SetBeforeTempStrandRemovalLosses(Float64 loss)
{
   m_pLossParameters->SetBeforeTempStrandRemovalLosses(loss);
}

Float64 CPyLossParameters::GetAfterTempStrandRemovalLosses() const
{
   return m_pLossParameters->GetAfterTempStrandRemovalLosses();
}

void CPyLossParameters::SetAfterTempStrandRemovalLosses(Float64 loss)
{
   m_pLossParameters->SetAfterTempStrandRemovalLosses(loss);
}

Float64 CPyLossParameters::GetAfterDeckPlacementLosses() const
{
   return m_pLossParameters->GetAfterDeckPlacementLosses();
}

void CPyLossParameters::SetAfterDeckPlacementLosses(Float64 loss)
{
   m_pLossParameters->SetAfterDeckPlacementLosses(loss);
}

Float64 CPyLossParameters::GetAfterSIDLLosses() const
{
   return m_pLossParameters->GetAfterSIDLLosses();
}

void CPyLossParameters::SetAfterSIDLLosses(Float64 loss)
{
   m_pLossParameters->SetAfterSIDLLosses(loss);
}

Float64 CPyLossParameters::GetFinalLosses() const
{
   return m_pLossParameters->GetFinalLosses();
}

void CPyLossParameters::SetFinalLosses(Float64 loss)
{
   m_pLossParameters->SetFinalLosses(loss);
}
