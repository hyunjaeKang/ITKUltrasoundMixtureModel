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
#ifndef __itkExpectationMaximizationMixtureModelEstimator_hxx
#define __itkExpectationMaximizationMixtureModelEstimator_hxx

#include "itkExpectationMaximizationUSMixtureModelEstimator.h"
#include "itkNumericTraits.h"

namespace itk
{
namespace Statistics
{
template< class TSample >
ExpectationMaximizationUSMixtureModelEstimator< TSample >
::ExpectationMaximizationUSMixtureModelEstimator()
{
  m_TerminationCode = NOT_CONVERGED;

  m_MembershipFunctionsObject            = MembershipFunctionVectorObjectType::New();
  m_MembershipFunctionsWeightArrayObject =
    MembershipFunctionsWeightsArrayObjectType::New();
  m_Sample = 0;
  m_MaxIteration = 100;
}

template< class TSample >
void
ExpectationMaximizationUSMixtureModelEstimator< TSample >
::PrintSelf(std::ostream & os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);
  os << indent << "Maximum Iteration: "
     << this->GetMaximumIteration() << std::endl;
  os << indent << "Sample: "
     << this->GetSample() << std::endl;
  os << indent << "Number Of Components: "
     << this->GetNumberOfComponents() << std::endl;
  for ( unsigned int i = 0; i < this->GetNumberOfComponents(); i++ )
    {
    os << indent << "Component Membership Function[" << i << "]: "
       << this->GetComponentMembershipFunction(i) << std::endl;
    }
  os << indent << "Termination Code: "
     << this->GetTerminationCode() << std::endl;
  os << indent << "Initial Proportions: "
     << this->GetInitialProportions() << std::endl;
  os << indent << "Proportions: "
     << this->GetProportions() << std::endl;
  os << indent << "Calculated Expectation: " << this->CalculateExpectation() << std::endl;
}

template< class TSample >
void
ExpectationMaximizationUSMixtureModelEstimator< TSample >
::SetMaximumIteration(int numberOfIterations)
{
  m_MaxIteration = numberOfIterations;
}

template< class TSample >
int
ExpectationMaximizationUSMixtureModelEstimator< TSample >
::GetMaximumIteration() const
{
  return m_MaxIteration;
}

template< class TSample >
void
ExpectationMaximizationUSMixtureModelEstimator< TSample >
::SetInitialProportions(ProportionVectorType & proportions)
{
  m_InitialProportions = proportions;
}

template< class TSample >
const typename ExpectationMaximizationUSMixtureModelEstimator< TSample >::ProportionVectorType &
ExpectationMaximizationUSMixtureModelEstimator< TSample >
::GetInitialProportions() const
{
  return m_InitialProportions;
}

template< class TSample >
const typename ExpectationMaximizationUSMixtureModelEstimator< TSample >::ProportionVectorType &
ExpectationMaximizationUSMixtureModelEstimator< TSample >
::GetProportions() const
{
  return m_Proportions;
}

template< class TSample >
void
ExpectationMaximizationUSMixtureModelEstimator< TSample >
::SetSample(const TSample *sample)
{
  m_Sample = sample;
}

template< class TSample >
const TSample *
ExpectationMaximizationUSMixtureModelEstimator< TSample >
::GetSample() const
{
  return m_Sample;
}

template< class TSample >
int
ExpectationMaximizationUSMixtureModelEstimator< TSample >
::AddComponent(ComponentType *component)
{
  m_ComponentVector.push_back(component);
  return static_cast< int >( m_ComponentVector.size() );
}


template< class TSample >
const typename ExpectationMaximizationUSMixtureModelEstimator< TSample >::ComponentVectorType &
ExpectationMaximizationUSMixtureModelEstimator< TSample >
::GetComponents() const
{
return m_ComponentVector;
}

template< class TSample >
typename ExpectationMaximizationUSMixtureModelEstimator< TSample >::ComponentType *
ExpectationMaximizationUSMixtureModelEstimator< TSample >
::GetComponent(unsigned int i) const
{
  if (i < m_ComponentVector.size())
    return m_ComponentVector[i];
  else
    return NULL;
}


template< class TSample >
unsigned int
ExpectationMaximizationUSMixtureModelEstimator< TSample >
::GetNumberOfComponents() const
{
  return m_ComponentVector.size();
}

template< class TSample >
typename ExpectationMaximizationUSMixtureModelEstimator< TSample >::TERMINATION_CODE
ExpectationMaximizationUSMixtureModelEstimator< TSample >
::GetTerminationCode() const
{
  return m_TerminationCode;
}

template< class TSample >
typename ExpectationMaximizationUSMixtureModelEstimator< TSample >::ComponentMembershipFunctionType *
ExpectationMaximizationUSMixtureModelEstimator< TSample >
::GetComponentMembershipFunction(int componentIndex) const
{
  return ( m_ComponentVector[componentIndex] )->GetMembershipFunction();
}

template< class TSample >
bool
ExpectationMaximizationUSMixtureModelEstimator< TSample >
::CalculateDensities()
{
  bool componentModified = false;

  for ( size_t i = 0; i < m_ComponentVector.size(); i++ )
    {
    if ( ( m_ComponentVector[i] )->AreParametersModified() )
      {
      componentModified = true;
      break;
      }
    }

  if ( !componentModified )
    {
    return false;
    }

  double                temp;
  size_t                numberOfComponents = m_ComponentVector.size();
  std::vector< double > tempWeights(numberOfComponents, 0. );

  typename TSample::ConstIterator iter = m_Sample->Begin();
  typename TSample::ConstIterator last = m_Sample->End();

  size_t componentIndex;

  typedef typename TSample::AbsoluteFrequencyType FrequencyType;
  FrequencyType frequency;
  FrequencyType zeroFrequency = NumericTraits< FrequencyType >::Zero;
  typename TSample::MeasurementVectorType mvector;
  double density;
  double densitySum;
  double minDouble = NumericTraits<double>::epsilon();

  SizeValueType measurementVectorIndex = 0;

  while ( iter != last )
    {
    mvector = iter.GetMeasurementVector();
    frequency = iter.GetFrequency();
    densitySum = 0.0;
    if ( frequency > zeroFrequency )
      {
      for ( componentIndex = 0; componentIndex < numberOfComponents;
            ++componentIndex )
        {
        double t_prop = m_Proportions[componentIndex];
        double t_value = m_ComponentVector[componentIndex]->Evaluate(mvector);
        density = t_prop * t_value;
        tempWeights[componentIndex] = density;
        densitySum += density;
        }

      for ( componentIndex = 0; componentIndex < numberOfComponents;
            ++componentIndex )
        {
        temp = tempWeights[componentIndex];

        // just to make sure temp does not blow up!
        if ( densitySum > NumericTraits<double>::epsilon() )
          {
          temp /= densitySum;
          }
        m_ComponentVector[componentIndex]->SetWeight(measurementVectorIndex,
                                                     temp);
        }
      }
    else
      {
      for ( componentIndex = 0; componentIndex < numberOfComponents;
            ++componentIndex )
        {
        m_ComponentVector[componentIndex]->SetWeight(measurementVectorIndex,
                                                     minDouble);
        }
      }

    ++iter;
    ++measurementVectorIndex;
    }

  return true;
}

template< class TSample >
double
ExpectationMaximizationUSMixtureModelEstimator< TSample >
::CalculateExpectation() const
{
  double sum = 0.0;

  if ( m_Sample )
    {
    unsigned int  measurementVectorIndex;
    SizeValueType size = m_Sample->Size();
    double        logProportion;
    double        temp;
    for ( size_t componentIndex = 0;
          componentIndex < m_ComponentVector.size();
          ++componentIndex )
      {
      temp = m_Proportions[componentIndex];

      // if temp is below the smallest positive double number
      // the log may blow up
      if( temp > NumericTraits<double>::epsilon() )
        {
        logProportion = vcl_log( temp );
        }
      else
        {
        logProportion = NumericTraits< double >::NonpositiveMin();
        }
      for ( measurementVectorIndex = 0; measurementVectorIndex < size;
            measurementVectorIndex++ )
        {
        temp = m_ComponentVector[componentIndex]->
               GetWeight(measurementVectorIndex);
        if( temp > NumericTraits<double>::epsilon() )
          {
          sum += temp * ( logProportion + vcl_log( temp ) );
          }
        else
          {
          // let's throw an exception
          itkExceptionMacro( << "temp is null" );
          }
        }
      }
    }
  return sum;
}

template< class TSample >
bool
ExpectationMaximizationUSMixtureModelEstimator< TSample >
::UpdateComponentParameters()
{
  bool           updated = false;
  ComponentType *component;

  for ( size_t componentIndex = 0; componentIndex < m_ComponentVector.size();
        ++componentIndex )
    {
    component = m_ComponentVector[componentIndex];
    component->Update();
      // Update all the componets
//    if ( component->AreParametersModified() )
//      {
//      return true;
//      }
    }

  return updated;
}

template< class TSample >
bool
ExpectationMaximizationUSMixtureModelEstimator< TSample >
::UpdateProportions()
{
  size_t numberOfComponents = m_ComponentVector.size();
  size_t sampleSize = m_Sample->Size();
  double totalFrequency = static_cast< double >( m_Sample->GetTotalFrequency() );
  size_t   i, j;
  double tempSum;
  bool   updated = false;

  for ( i = 0; i < numberOfComponents; ++i )
    {
    tempSum = 0.;

    if( totalFrequency > NumericTraits<double>::epsilon() )
      {
      for ( j = 0; j < sampleSize; ++j )
        {
        tempSum += ( m_ComponentVector[i]->GetWeight(j)
                     * m_Sample->GetFrequency(j) );
        }

      tempSum /= totalFrequency;
      }

    if ( tempSum != m_Proportions[i] )
      {
      m_Proportions[i] = tempSum;
      updated = true;
      }
    }

  return updated;
}

template< class TSample >
void
ExpectationMaximizationUSMixtureModelEstimator< TSample >
::GenerateData()
{
  m_Proportions = m_InitialProportions;

  int iteration = 0;
  m_CurrentIteration = 0;
  while ( iteration < m_MaxIteration )
    {
    m_CurrentIteration = iteration;

    this->InvokeEvent(IterationEvent());


    if ( this->CalculateDensities() )
      {
      this->UpdateComponentParameters();
      this->UpdateProportions();
      }
    else
      {
      m_TerminationCode = CONVERGED;
      break;
      }
    ++iteration;
    }
}



template< class TSample >
const typename ExpectationMaximizationUSMixtureModelEstimator< TSample
                                                             >::MembershipFunctionsWeightsArrayObjectType *
ExpectationMaximizationUSMixtureModelEstimator< TSample >
::GetMembershipFunctionsWeightsArray() const
{
  size_t           numberOfComponents = m_ComponentVector.size();
  ProportionVectorType & membershipFunctionsWeightVector =
    m_MembershipFunctionsWeightArrayObject->Get();

  membershipFunctionsWeightVector.SetSize(numberOfComponents);
  for ( size_t i = 0; i < numberOfComponents; ++i )
    {
    membershipFunctionsWeightVector[i] = m_Proportions[i];
    }

  return static_cast< const MembershipFunctionsWeightsArrayObjectType * >( m_MembershipFunctionsWeightArrayObject );
}

template< class TSample >
void
ExpectationMaximizationUSMixtureModelEstimator< TSample >
::Update()
{
  this->GenerateData();
}



template< class TSample >
const typename
ExpectationMaximizationUSMixtureModelEstimator< TSample >::MembershipFunctionVectorObjectType *
ExpectationMaximizationUSMixtureModelEstimator< TSample >
::GetOutput() const
{

  size_t numberOfComponents = m_ComponentVector.size();

  MembershipFunctionVectorType & membershipFunctionsVector = m_MembershipFunctionsObject->Get();


  for ( size_t i = 0; i < numberOfComponents; ++i )
    {
    LightObject::Pointer loPtr = m_ComponentVector[i]->GetMembershipFunction()->Clone();

    MembershipFunctionType * membershipFunction =
    dynamic_cast<MembershipFunctionType *>(loPtr.GetPointer());

    membershipFunctionsVector.push_back( membershipFunction);
    }

return static_cast< const MembershipFunctionVectorObjectType * >( m_MembershipFunctionsObject );
}

} // end of namespace Statistics
} // end of namespace itk

#endif