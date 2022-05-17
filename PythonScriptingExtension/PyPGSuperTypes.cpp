#include "stdafx.h"
#include "PyPGSuperTypes.h"

using namespace boost::python;

void export_pgsuper_types()
{
   enum_<pgsTypes::AnalysisType>("AnalysisType")
      .value("Simple", pgsTypes::Simple)
      .value("Continuous", pgsTypes::Continuous)
      .value("Envelope", pgsTypes::Envelope)
      ;

   enum_<pgsTypes::StressLocation>("StressLocation")
      .value("TopGirder", pgsTypes::TopGirder)
      .value("BottomGirder", pgsTypes::BottomGirder)
      .value("TopDeck", pgsTypes::TopDeck)
      .value("BottomDeck", pgsTypes::BottomDeck)
      ;


   enum_<pgsTypes::TrafficBarrierOrientation>("TrafficBarrierOrientation")
      .value("Left", pgsTypes::tboLeft)
      .value("Right", pgsTypes::tboRight)
      ;

   enum_<pgsTypes::SectionPropertyType>("SectionPropertyType")
      .value("GrossNoncomposite", pgsTypes::sptGrossNoncomposite)
      .value("Gross", pgsTypes::sptGross)
      .value("TransformedNoncomposite", pgsTypes::sptTransformedNoncomposite)
      .value("Transformed", pgsTypes::sptTransformed)
      .value("NetGirder", pgsTypes::sptNetGirder)
      .value("NetDeck", pgsTypes::sptNetDeck)
      ;

   enum_<pgsTypes::MemberEndType>("MemberEndType")
      .value("Start", pgsTypes::metStart)
      .value("End", pgsTypes::metEnd)
      ;

   enum_<pgsTypes::ConcreteType>("ConcreteType")
      .value("Normal", pgsTypes::Normal)
      .value("AllLightweight", pgsTypes::AllLightweight)
      .value("SandLightweight", pgsTypes::SandLightweight)
      ;

   enum_<pgsTypes::LossMethod>("LossMethod")
      .value("AASHTO_REFINED", pgsTypes::AASHTO_REFINED)
      .value("AASHTO_LUMPSUM", pgsTypes::AASHTO_LUMPSUM)
      .value("GENERAL_LUMPSUM",pgsTypes::GENERAL_LUMPSUM)
      .value("WSDOT_LUMPSUM",pgsTypes::WSDOT_LUMPSUM)
      .value("AASHTO_LUMPSUM_2005",pgsTypes::AASHTO_LUMPSUM_2005)
      .value("AASHTO_REFINED_2005",pgsTypes::AASHTO_REFINED_2005)
      .value("WSDOT_LUMPSUM_2005",pgsTypes::WSDOT_LUMPSUM_2005)
      .value("WSDOT_REFINED_2005",pgsTypes::WSDOT_REFINED_2005)
      .value("WSDOT_REFINED",pgsTypes::WSDOT_REFINED)
      .value("TXDOT_REFINED_2004",pgsTypes::TXDOT_REFINED_2004)
      .value("TXDOT_REFINED_2013",pgsTypes::TXDOT_REFINED_2013)
      .value("TIME_STEP",pgsTypes::TIME_STEP)
      ;

   enum_<pgsTypes::TimeDependentModel>("TimeDependentModel")
      .value("AASHTO", pgsTypes::tdmAASHTO)
      .value("ACI209", pgsTypes::tdmACI209)
      .value("CEBFIP", pgsTypes::tdmCEBFIP)
      ;
}