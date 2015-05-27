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
#ifndef __itkRayleighMixtureModelComponent_hxx
#define __itkRayleighMixtureModelComponent_hxx

#include <iostream>


namespace itk
{
namespace Statistics
{
template< class TSample >
RayleighMixtureModelComponent< TSample >
::RayleighMixtureModelComponent()
{
  m_RayleighMembershipFunction = NativeMembershipFunctionType::New();
  this->SetMembershipFunction( (MembershipFunctionType *)
                               m_RayleighMembershipFunction.GetPointer() );
  m_Sigma.Fill(0.0);
  m_SigmaEstimate.Fill(0.0);
}

template< class TSample >
void
RayleighMixtureModelComponent< TSample >
::PrintSelf(std::ostream & os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);
  os << indent << "Sigma: " << m_Sigma << std::endl;
  os << indent << "Sigma Estimated: " << m_SigmaEstimate << std::endl;
  os << indent << "RayleighMembershipFunction: " << m_RayleighMembershipFunction << std::endl;
}

template< class TSample >
void
RayleighMixtureModelComponent< TSample >
::SetSample(const TSample *sample)
{
  Superclass::SetSample(sample);

  const MeasurementVectorSizeType measurementVectorLength =
  sample->GetMeasurementVectorSize();

  if(measurementVectorLength != 1)
    itkExceptionMacro(<< "Measurement Vector dimension must be 1");

  m_RayleighMembershipFunction->SetMeasurementVectorSize(measurementVectorLength);

  NumericTraits<MeasurementVectorType>::SetLength(m_Sigma, measurementVectorLength);
  m_Sigma.Fill(NumericTraits< double >::One);

  NumericTraits<MeasurementVectorType>::SetLength(m_SigmaEstimate, measurementVectorLength);
  m_SigmaEstimate.Fill(NumericTraits< double >::One);

  typename NativeMembershipFunctionType::SigmaVectorType sigma;
  NumericTraits<typename NativeMembershipFunctionType::SigmaVectorType>::SetLength(sigma,
                                                                                  measurementVectorLength);


  for ( unsigned int i = 0; i < measurementVectorLength; ++i )
    {
    sigma[i] = m_Sigma[i];
    }
  

  m_RayleighMembershipFunction->SetSigma(sigma);
}

template< class TSample >
void
RayleighMixtureModelComponent< TSample >
::SetParameters(const ParametersType & parameters)
{
  Superclass::SetParameters(parameters);

  unsigned int paramIndex = 0;
  bool changed = false;

  unsigned int i;

  MeasurementVectorSizeType measurementVectorSize =
  this->GetSample()->GetMeasurementVectorSize();


  // Sigma parameter
  for ( i = 0; i < measurementVectorSize; i++ )
    {
    if ( m_Sigma[i] != parameters[paramIndex] )
      {
      m_Sigma[i] = parameters[paramIndex];

      m_SigmaEstimate[i] = parameters[paramIndex];
      
      changed = true;
      }

    ++paramIndex;
    }

  typename NativeMembershipFunctionType::SigmaVectorType sigma;
  NumericTraits<typename NativeMembershipFunctionType::SigmaVectorType>::SetLength(sigma,
                                                                                  measurementVectorSize);
  for ( i = 0; i < measurementVectorSize; ++i )
    {
    sigma[i] = m_Sigma[i];
    }

  m_RayleighMembershipFunction->SetSigma(sigma);

  this->AreParametersModified(changed);
}

template< class TSample >
double
RayleighMixtureModelComponent< TSample >
::CalculateParametersChange()
{
  unsigned int i;


  double                    temp;
  double                    changes = 0.0;
  MeasurementVectorSizeType measurementVectorSize =
  this->GetSample()->GetMeasurementVectorSize();

  for ( i = 0; i < measurementVectorSize; i++ )
    {
    temp = m_Sigma[i] - m_SigmaEstimate[i];
    changes += temp * temp;
    }


  changes = vcl_sqrt(changes);
  return changes;

}

template< class TSample >
void
RayleighMixtureModelComponent< TSample >
::GenerateData()
{

  MeasurementVectorSizeType measurementVectorSize =
  this->GetSample()->GetMeasurementVectorSize();

  this->AreParametersModified(false);


  typename TSample::MeasurementVectorType measurements;


  this->CalculateSigma();

  MeasurementVectorSizeType   i;
  double         changes;
  bool           changed = false;
  ParametersType parameters = this->GetFullParameters(); // parameters = [sigma]
  MeasurementVectorSizeType            paramIndex  = 0;


  for ( i = 0; i < measurementVectorSize; i++ )
    {
    changes = vnl_math_abs( m_Sigma[i] - m_SigmaEstimate[i] );

    if ( changes > this->GetMinimalParametersChange() )
      {
      changed = true;
      break;
      }
    }

  if ( changed )
    {
    for ( i = 0; i < measurementVectorSize; i++ )
      {
      m_Sigma[i] = m_SigmaEstimate[i];
      }

    
    for ( paramIndex = 0; paramIndex < measurementVectorSize; paramIndex++ )
      {
      parameters[paramIndex] = m_SigmaEstimate[paramIndex];
      }
    this->AreParametersModified(true);
    }
  else
    {
    paramIndex = measurementVectorSize;
    }

  
  //THIS IS NEEDED TO update m_Sigma. SHOULD BE REMOVED

  // Alpha update
  paramIndex = 0;
  for ( i = 0; i < measurementVectorSize; i++ )
    {
    m_Sigma[i] = parameters[paramIndex];
    ++paramIndex;
    }

  typename NativeMembershipFunctionType::SigmaVectorType sigma;
  NumericTraits<typename NativeMembershipFunctionType::SigmaVectorType>::SetLength(sigma,
                                                                                  measurementVectorSize);

  for ( i = 0; i < measurementVectorSize; ++i )
    {
    sigma[i] = m_Sigma[i];
    }
  m_RayleighMembershipFunction->SetSigma(sigma);

  
  
  Superclass::SetParameters(parameters);
}

template< class TSample >
void
RayleighMixtureModelComponent< TSample >
::CalculateSigma()
{
  const WeightArrayType & weights = this->GetWeights();
  unsigned long int idW = 0;


  typename TSample::ConstIterator iter = this->GetSample()->Begin();
  typename TSample::ConstIterator end =  this->GetSample()->End();
  typename TSample::MeasurementVectorType measurements;


  MeasurementVectorSizeType measurementVectorSize =
  this->GetSample()->GetMeasurementVectorSize();




//  ke(n) = sum(w_t(n,:))/(sum(w_t(:)));
//  fe(n) = sum(w_t(n,:).*y.^2/2)/(sum(w_t(n,:))+eps);
//
//  if isnan(fe(n)), keyboard, end
//
//    p(n,:) = (y./(fe(n)+eps)).*exp(-y.^2/(2*(fe(n)+eps))) + eps;
//  K = K + p(n,:).*ke(n);


  typename NativeMembershipFunctionType::SigmaVectorType cte1;
  NumericTraits<typename NativeMembershipFunctionType::SigmaVectorType>::SetLength(cte1,
                                                                                   measurementVectorSize);
  cte1.Fill(NumericTraits< double >::Zero);


  double tmpWeights = 0;

  while ( iter != end )
    {
    measurements = iter.GetMeasurementVector();
    double w = weights[idW];
    tmpWeights += w;
    for ( unsigned int dim = 0; dim < measurementVectorSize; ++dim )
      {
      cte1[dim] += vcl_pow(measurements[dim],2) * w;
      }

    
    ++iter;
    ++idW;
    }

  tmpWeights *= 2;

  tmpWeights += vnl_math::eps;


  for ( unsigned int dim = 0; dim < measurementVectorSize; ++dim )
    {
    double tmp = cte1[dim] / tmpWeights;
    if (tmp != tmp)
      {
      std::cout << "RayleighMixtureModelComponent:: Nan is reached" << std::endl;
      itkExceptionMacro(<<"Nan is reached");
      }
    m_SigmaEstimate[dim] = vcl_sqrt(tmp);
    }
}



template< class TSample >
double
RayleighMixtureModelComponent< TSample >
::SortValue() const
{
    // The sort order is specified by the mean value E{} = sigma * sqrt(pi/2)
  double mean = m_Sigma[0] * vcl_sqrt(vnl_math::pi/2);

  return mean;
}

} // end of namespace Statistics
} // end of namespace itk

#endif
