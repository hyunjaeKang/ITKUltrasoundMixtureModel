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
#ifndef __itkRayleighMembershipFunction_hxx
#define __itkRayleighMembershipFunction_hxx

#include "itkRayleighMembershipFunction.h"

namespace itk
{
namespace Statistics
{
template< class TMeasurementVector >
RayleighMembershipFunction< TMeasurementVector >
::RayleighMembershipFunction()
{
  NumericTraits<SigmaVectorType>::SetLength(m_Sigma, this->GetMeasurementVectorSize());
  m_Sigma.Fill( NumericTraits< double >::One );
}

template< class TMeasurementVector >
void
RayleighMembershipFunction< TMeasurementVector >
::PrintSelf(std::ostream & os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);

  os << indent << "Sigma: " << m_Sigma << std::endl;
}

template< class TMeasurementVector >
void
RayleighMembershipFunction< TMeasurementVector >
::SetSigma(const SigmaVectorType & sigma)
{
  if ( this->GetMeasurementVectorSize() )
    {
    MeasurementVectorTraits::Assert(sigma,
                                    this->GetMeasurementVectorSize(),
                                    "RayleighMembershipFunction::SetSigma(): Size of mean vector specified does not match the size of a measurement vector.");
    }
  else
    {
    // not already set, cache the size
    this->SetMeasurementVectorSize( sigma.Size() );
    }


  if(sigma.Size() != 1)
    itkExceptionMacro( << "Measurement Vector dimension must be 1" );

  if(sigma[0] <= 0)
    itkExceptionMacro( << "Sigma must be > 0" );

  if ( m_Sigma != sigma )
    {
    m_Sigma = sigma;
    this->Modified();
    }
}



template< class TMeasurementVector >
void
RayleighMembershipFunction< TMeasurementVector >
::SetParameters(const ParameterVectorRealType & param)
{
this->SetSigma(param[0]);
}


template< class TMeasurementVector >
inline double
RayleighMembershipFunction< TMeasurementVector >
::Evaluate(const MeasurementVectorType & measurement) const
{
  if(this->GetMeasurementVectorSize() != 1)
    itkExceptionMacro( << "Measurement Vector dimension must be 1" );

  double xx = vcl_pow(measurement[0],2);
  double sigma = vcl_pow(m_Sigma[0],2) + vnl_math::eps;
  double result = measurement[0]/sigma * vcl_exp(-xx/(2*sigma)) + vnl_math::eps;

  return result;
}


template< class TVector >
typename LightObject::Pointer
RayleighMembershipFunction< TVector >
::InternalClone() const
{
  LightObject::Pointer loPtr = Superclass::InternalClone();
  typename Self::Pointer membershipFunction =
    dynamic_cast<Self *>(loPtr.GetPointer());
  if(membershipFunction.IsNull())
    {
    itkExceptionMacro(<< "downcast to type "
                      << this->GetNameOfClass()
                      << " failed.");
    }

  membershipFunction->SetMeasurementVectorSize( this->GetMeasurementVectorSize() );
  membershipFunction->SetSigma(this->GetSigma());


  return loPtr;
}
} // end namespace Statistics
} // end of namespace itk

#endif
