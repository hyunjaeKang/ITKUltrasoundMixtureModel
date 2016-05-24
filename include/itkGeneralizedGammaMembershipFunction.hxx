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
#ifndef __itkGeneralizedGammaMembershipFunction_hxx
#define __itkGeneralizedGammaMembershipFunction_hxx

#include "itkGeneralizedGammaMembershipFunction.h"

#include "vnl/vnl_gamma.h"


namespace itk
{
namespace Statistics
{
template< class TMeasurementVector >
GeneralizedGammaMembershipFunction< TMeasurementVector >
::GeneralizedGammaMembershipFunction()
{
  NumericTraits<ParameterVectorRealType>::SetLength(m_L, this->GetMeasurementVectorSize());
  m_L.Fill( NumericTraits< double >::One );

  NumericTraits<ParameterVectorRealType>::SetLength(m_A, this->GetMeasurementVectorSize());
  m_A.Fill( NumericTraits< double >::One );

  NumericTraits<ParameterVectorRealType>::SetLength(m_C, this->GetMeasurementVectorSize());
  m_C.Fill( NumericTraits< double >::One );

}

template< class TMeasurementVector >
void
GeneralizedGammaMembershipFunction< TMeasurementVector >
::PrintSelf(std::ostream & os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);

  os << indent << "L: " << m_L << std::endl;
  os << indent << "A: " << m_A << std::endl;
  os << indent << "C: " << m_C << std::endl;
}

template< class TMeasurementVector >
void
GeneralizedGammaMembershipFunction< TMeasurementVector >
::SetA(const ParameterVectorRealType & a)
{
  if ( this->GetMeasurementVectorSize() )
    {
    MeasurementVectorTraits::Assert(a, this->GetMeasurementVectorSize(),
                                    "GeneralizedGammaMembershipFunction::SetA(): Size of mean vector specified does not match the size of a measurement vector.");
    }
  else
    {
    // not already set, cache the size
    this->SetMeasurementVectorSize( a.Size() );
    }


  if(a.Size() != 1)
    itkExceptionMacro( << "Measurement Vector dimension must be 1" );

  if(a[0] <= 0)
    itkExceptionMacro( << "A must be > 0" );

  if ( m_A != a )
    {
    m_A = a;
    this->Modified();
    }
}

template< class TMeasurementVector >
void
GeneralizedGammaMembershipFunction< TMeasurementVector >
::SetL(const ParameterVectorRealType & l)
{
if ( this->GetMeasurementVectorSize() )
  {
  MeasurementVectorTraits::Assert(l, this->GetMeasurementVectorSize(),
                                  "GeneralizedGammaMembershipFunction::SetL(): Size of mean vector specified does not match the size of a measurement vector.");
  }
else
  {
    // not already set, cache the size
  this->SetMeasurementVectorSize( l.Size() );
  }


if(l.Size() != 1)
  itkExceptionMacro( << "Measurement Vector dimension must be 1" );

if(l[0] <= 0)
  itkExceptionMacro( << "L must be > 0" );

if ( m_L != l )
  {
  m_L = l;
  this->Modified();
  }
}


template< class TMeasurementVector >
void
GeneralizedGammaMembershipFunction< TMeasurementVector >
::SetC(const ParameterVectorRealType & c)
{
if ( this->GetMeasurementVectorSize() )
  {
  MeasurementVectorTraits::Assert(c, this->GetMeasurementVectorSize(),
                                  "GeneralizedGammaMembershipFunction::SetC(): Size of mean vector specified does not match the size of a measurement vector.");
  }
else
  {
    // not already set, cache the size
  this->SetMeasurementVectorSize( c.Size() );
  }


if(c.Size() != 1)
  itkExceptionMacro( << "Measurement Vector dimension must be 1" );

if(c[0] <= 0)
  itkExceptionMacro( << "C must be > 0" );

if ( m_C != c )
  {
  m_C = c;
  this->Modified();
  }
}


template< class TMeasurementVector >
void
GeneralizedGammaMembershipFunction< TMeasurementVector >
::SetParameters(const ParameterVectorRealType & param)
{
  this->SetL(param[0]);
  this->SetA(param[1]);
  this->SetC(param[2]);
}



template< class TMeasurementVector >
inline double
GeneralizedGammaMembershipFunction< TMeasurementVector >
::Evaluate(const MeasurementVectorType & measurement) const
{
  if(this->GetMeasurementVectorSize() != 1)
    itkExceptionMacro( << "Measurement Vector dimension must be 1" );


//  f(x|l,a,c) = abs(c)* a^(l*c)* x^(l*c-1) * exp(-(a*x)^c)/Gamma(l)

  double result = vcl_abs(m_C[0]) * vcl_pow(m_A[0], m_L[0]*m_C[0])
  * vcl_pow(measurement[0], m_L[0]*m_C[0]-1)
  * vcl_exp(-vcl_pow(measurement[0]*m_A[0], m_C[0]));

  result /= vnl_gamma(m_L[0]);

  return result;
}

template< class TMeasurementVector >
inline double
GeneralizedGammaMembershipFunction< TMeasurementVector >
::DerivativeEvaluation(const MeasurementVectorType & measurement) const
{
  if(this->GetMeasurementVectorSize() != 1)
    itkExceptionMacro( << "Measurement Vector dimension must be 1" );


  itkExceptionMacro( << "DerivativeEvaluateion() is not implemented !" );

  return 0;
}

template< class TVector >
typename LightObject::Pointer
GeneralizedGammaMembershipFunction< TVector >
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

  membershipFunction->SetL(this->GetL());
  membershipFunction->SetA(this->GetA());
  membershipFunction->SetC(this->GetC());


  return loPtr;
}
} // end namespace Statistics
} // end of namespace itk

#endif
