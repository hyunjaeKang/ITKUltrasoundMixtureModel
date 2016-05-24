/*=========================================================================
 *
 *  Copyright Insight Software Consortium
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *         http://www.apache.org/licenses/LICENSE-2.0.txt
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *
 *=========================================================================*/
#ifndef __itkGeneralizeGammaExpectationMaximizationMixtureModelEstimator_hxx
#define __itkGeneralizeGammaExpectationMaximizationMixtureModelEstimator_hxx

//#include "itkGeneralizeGammaExpectationMaximizationMixtureModelEstimator.h"
#include "itkGeneralizedGammaExpectationMaximizationMixtureModelEstimator.h"


namespace itk
{
namespace Statistics
{


template< class TSample >
GeneralizeGammaExpectationMaximizationMixtureModelEstimator< TSample >
::GeneralizeGammaExpectationMaximizationMixtureModelEstimator()
{}

template< class TSample >
const typename GeneralizeGammaExpectationMaximizationMixtureModelEstimator< TSample >::ComponentVectorType &
GeneralizeGammaExpectationMaximizationMixtureModelEstimator< TSample >
::GetComponents() const
{
  return Superclass::m_ComponentVector;
}





template< class TSample >
const typename
  GeneralizeGammaExpectationMaximizationMixtureModelEstimator< TSample >::MembershipFunctionVectorObjectType *
GeneralizeGammaExpectationMaximizationMixtureModelEstimator< TSample >
::GetOutput() const
{

  size_t                   numberOfComponents = Superclass::m_ComponentVector.size();
  MembershipFunctionVectorType & membershipFunctionsVector =
    Superclass::m_MembershipFunctionsObject->Get();

  MeasurementVectorSizeType measurementVectorSize =
    Superclass::m_Sample->GetMeasurementVectorSize();


  if(measurementVectorSize != 1)
    itkExceptionMacro( << "Measurement Vector dimension must be 1" );


  ParameterVectorRealType l;
  NumericTraits<ParameterVectorRealType>::SetLength(l, measurementVectorSize);


  ParameterVectorRealType a;
  NumericTraits<ParameterVectorRealType>::SetLength(a, measurementVectorSize);

  ParameterVectorRealType c;
  NumericTraits<ParameterVectorRealType>::SetLength(c, measurementVectorSize);

  typename Superclass::ComponentType::ParametersType parameters;

  for ( size_t i = 0; i < numberOfComponents; ++i )
    {

    parameters = Superclass::m_ComponentVector[i]->GetFullParameters();
    typename GeneralizeGammaMembershipFunctionType::Pointer membershipFunction =
    GeneralizeGammaMembershipFunctionType::New();
    membershipFunction->SetMeasurementVectorSize(measurementVectorSize);
    unsigned int parameterIndex = 0;
    for ( unsigned int j = 0; j < measurementVectorSize; j++ )
      {
      l[j] = parameters[parameterIndex];
      ++parameterIndex;
      }

    for ( unsigned int j = 0; j < measurementVectorSize; j++ )
      {
      a[j] = parameters[parameterIndex];
      ++parameterIndex;
      }

    for ( unsigned int j = 0; j < measurementVectorSize; j++ )
      {
      c[j] = parameters[parameterIndex];
      ++parameterIndex;
      }

    membershipFunction->SetL(l);
    membershipFunction->SetA(a);
    membershipFunction->SetC(c);
    membershipFunctionsVector.push_back( membershipFunction.GetPointer() );
    }

  return static_cast< const MembershipFunctionVectorObjectType * >( Superclass::m_MembershipFunctionsObject );



}

} // end of namespace Statistics
} // end of namespace itk

#endif
