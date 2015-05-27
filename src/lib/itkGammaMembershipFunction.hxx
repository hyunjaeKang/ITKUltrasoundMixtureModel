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
#ifndef __itkGammaMembershipFunction_hxx
#define __itkGammaMembershipFunction_hxx

#include "itkGammaMembershipFunction.h"

#include "vnl/vnl_gamma.h"


namespace itk
{
namespace Statistics
{
template< class TMeasurementVector >
GammaMembershipFunction< TMeasurementVector >
::GammaMembershipFunction()
{
  NumericTraits<AlphaVectorType>::SetLength(m_Alpha,
                                            this->GetMeasurementVectorSize());
  m_Alpha.Fill( NumericTraits< double >::One );

  NumericTraits<BetaVectorType>::SetLength(m_Beta,
                                           this->GetMeasurementVectorSize());
  m_Beta.Fill( NumericTraits< double >::One );

}

template< class TMeasurementVector >
void
GammaMembershipFunction< TMeasurementVector >
::PrintSelf(std::ostream & os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);

  os << indent << "Alpha: " << m_Alpha << std::endl;
  os << indent << "Beta: " << m_Beta << std::endl;
}

template< class TMeasurementVector >
void
GammaMembershipFunction< TMeasurementVector >
::SetAlpha(const AlphaVectorType & alpha)
{
  if ( this->GetMeasurementVectorSize() )
    {
    MeasurementVectorTraits::Assert(alpha,
                                    this->GetMeasurementVectorSize(),
                                    "GammaMembershipFunction::SetMean(): Size of mean vector specified does not match the size of a measurement vector.");
    }
  else
    {
    // not already set, cache the size
    this->SetMeasurementVectorSize( alpha.Size() );
    }


  if(alpha.Size() != 1)
    itkExceptionMacro( << "Measurement Vector dimension must be 1" );

  if(alpha[0] <= 0)
    itkExceptionMacro( << "Alpha must be > 0" );

  if ( m_Alpha != alpha )
    {
    m_Alpha = alpha;
    this->Modified();
    }
}



template< class TMeasurementVector >
void
GammaMembershipFunction< TMeasurementVector >
::SetBeta(const BetaVectorType & beta)
{
  if ( this->GetMeasurementVectorSize() )
    {
    MeasurementVectorTraits::Assert(beta,
                                    this->GetMeasurementVectorSize(),
                                    "GammaMembershipFunction::SetMean(): Size of mean vector specified does not match the size of a measurement vector.");
    }
  else
    {
    // not already set, cache the size
    this->SetMeasurementVectorSize( beta.Size() );
    }


  if(beta.Size() != 1)
    itkExceptionMacro( << "Measurement Vector dimension must be 1" );

  if(beta[0] <= 0)
    itkExceptionMacro( << "Beta must be > 0" );

  if ( m_Alpha != beta )
    {
    m_Beta = beta;
    this->Modified();
    }
}


template< class TMeasurementVector >
void
GammaMembershipFunction< TMeasurementVector >
::SetParameters(const ParameterVectorRealType & param)
{
this->SetAlpha(param[0]);
this->SetBeta(param[1]);

}


template< class TMeasurementVector >
inline double
GammaMembershipFunction< TMeasurementVector >
::Evaluate(const MeasurementVectorType & measurement) const
{
  if(this->GetMeasurementVectorSize() != 1)
    itkExceptionMacro( << "Measurement Vector dimension must be 1" );



  double result = 1/(vcl_pow(m_Beta[0],m_Alpha[0]) * vnl_gamma(m_Alpha[0])) *
    vcl_pow(measurement[0],m_Alpha[0]-1) * vcl_exp(-measurement[0]/m_Beta[0]);

  return result;
}

template< class TMeasurementVector >
inline double
GammaMembershipFunction< TMeasurementVector >
::DerivativeEvaluation(const MeasurementVectorType & measurement) const
{
  if(this->GetMeasurementVectorSize() != 1)
    itkExceptionMacro( << "Measurement Vector dimension must be 1" );

  double result = (m_Alpha[0] - 1 -  measurement[0]/m_Beta[0]) *
  vcl_pow(measurement[0],m_Alpha[0]-2) * vcl_exp(-measurement[0]/m_Beta[0]);

  return result/(vcl_pow(m_Beta[0],m_Alpha[0]) * vnl_gamma(m_Alpha[0]));

}

template< class TVector >
typename LightObject::Pointer
GammaMembershipFunction< TVector >
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

  membershipFunction->SetAlpha(this->GetAlpha());
  membershipFunction->SetBeta(this->GetBeta());


  return loPtr;
}
} // end namespace Statistics
} // end of namespace itk

#endif
