#pragma once

#include <IFace\Project.h>

class CPyLossParameters
{
public:
   void Init(IBroker* pBroker);
   void Reset();

   pgsTypes::LossMethod GetLossMethod() const;
   pgsTypes::TimeDependentModel GetTimeDependentModel() const;
   void SetIgnoreCreepEffects(bool bIgnore);
   bool GetIgnoreCreepEffects() const;
   void SetIgnoreShrinkageEffects(bool bIgnore);
   bool GetIgnoreShrinkageEffects() const;
   void SetIgnoreRelaxationEffects(bool bIgnore);
   bool GetIgnoreRelaxationEffects() const;
   void SetIgnoreTimeDependentEffects(bool bIgnoreCreep, bool bIgnoreShrinkage, bool bIgnoreRelaxation);
   void SetTendonPostTensionParameters(Float64 Dset, Float64 wobble, Float64 friction);
   void GetTendonPostTensionParameters(Float64* Dset, Float64* wobble, Float64* friction) const;
   void SetTemporaryStrandPostTensionParameters(Float64 Dset, Float64 wobble, Float64 friction);
   void GetTemporaryStrandPostTensionParameters(Float64* Dset, Float64* wobble, Float64* friction) const;
   void SetUseGeneralLumpSumLosses(bool bLumpSum);
   bool GetUseGeneralLumpSumLosses() const;
   Float64 GetBeforeXferLosses() const;
   void SetBeforeXferLosses(Float64 loss);
   Float64 GetAfterXferLosses() const;
   void SetAfterXferLosses(Float64 loss);
   Float64 GetLiftingLosses() const;
   void SetLiftingLosses(Float64 loss);
   Float64 GetShippingLosses() const;
   void SetShippingLosses(Float64 loss);
   Float64 GetBeforeTempStrandRemovalLosses() const;
   void SetBeforeTempStrandRemovalLosses(Float64 loss);
   Float64 GetAfterTempStrandRemovalLosses() const;
   void SetAfterTempStrandRemovalLosses(Float64 loss);
   Float64 GetAfterDeckPlacementLosses() const;
   void SetAfterDeckPlacementLosses(Float64 loss);
   Float64 GetAfterSIDLLosses() const;
   void SetAfterSIDLLosses(Float64 loss);
   Float64 GetFinalLosses() const;
   void SetFinalLosses(Float64 loss);

private:
   CComPtr<ILossParameters> m_pLossParameters;
}; 
