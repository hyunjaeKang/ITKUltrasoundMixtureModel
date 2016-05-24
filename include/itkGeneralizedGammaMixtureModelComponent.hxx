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
#ifndef __itkGeneralizedGammaMixtureModelComponent_hxx
#define __itkGeneralizedGammaMixtureModelComponent_hxx

#include <iostream>

#include <vnl/algo/vnl_amoeba.h>

#include "itkGeneralizedGammaMixtureModelComponent.h"
#include "itkGeneralizedGammaMLEFunction.h"




namespace itk
{
namespace Statistics
{
template< class TSample >
GeneralizedGammaMixtureModelComponent< TSample >
::GeneralizedGammaMixtureModelComponent()
{
  m_GeneralizedGammaMembershipFunction = NativeMembershipFunctionType::New();
  this->SetMembershipFunction( (MembershipFunctionType *)
                               m_GeneralizedGammaMembershipFunction.GetPointer() );


  m_L.Fill(0.0);
  m_A.Fill(0.0);
  m_C.Fill(0.0);

  m_LEstimate.Fill(0.0);
  m_AEstimate.Fill(0.0);
  m_CEstimate.Fill(0.0);

}

template< class TSample >
void
GeneralizedGammaMixtureModelComponent< TSample >
::PrintSelf(std::ostream & os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);

  os << indent << "L: " << m_L << std::endl;
  os << indent << "A: " << m_A << std::endl;
  os << indent << "C: " << m_C << std::endl;
  os << indent << "L Estimated: " << m_LEstimate << std::endl;
  os << indent << "A Estimated" << m_AEstimate << std::endl;
  os << indent << "C Estimated" << m_CEstimate << std::endl;
  os << indent << "GeneralizedGammaMembershipFunction: "
     << m_GeneralizedGammaMembershipFunction << std::endl;
}

template< class TSample >
void
GeneralizedGammaMixtureModelComponent< TSample >
::SetSample(const TSample *sample)
{
  Superclass::SetSample(sample);

  const MeasurementVectorSizeType measurementVectorLength =
  sample->GetMeasurementVectorSize();

  if(measurementVectorLength != 1)
    itkExceptionMacro(<< "Measurement Vector dimension must be 1");


  m_GeneralizedGammaMembershipFunction->SetMeasurementVectorSize(measurementVectorLength);

  NumericTraits<MeasurementVectorType>::SetLength(m_L, measurementVectorLength);
  m_L.Fill(NumericTraits< double >::One);

  NumericTraits<MeasurementVectorType>::SetLength(m_A, measurementVectorLength);
  m_A.Fill(NumericTraits< double >::One);

  NumericTraits<MeasurementVectorType>::SetLength(m_C, measurementVectorLength);
  m_C.Fill(NumericTraits< double >::One);

  NumericTraits<MeasurementVectorType>::SetLength(m_LEstimate, measurementVectorLength);
  m_LEstimate.Fill(NumericTraits< double >::One);

  NumericTraits<MeasurementVectorType>::SetLength(m_AEstimate, measurementVectorLength);
  m_AEstimate.Fill(NumericTraits< double >::One);

  NumericTraits<MeasurementVectorType>::SetLength(m_CEstimate, measurementVectorLength);
  m_CEstimate.Fill(NumericTraits< double >::One);


  ParameterVectorRealType l;
  ParameterVectorRealType a;
  ParameterVectorRealType c;

  NumericTraits<ParameterVectorRealType>::SetLength(l, measurementVectorLength);

  NumericTraits<ParameterVectorRealType>::SetLength(a, measurementVectorLength);

  NumericTraits<ParameterVectorRealType>::SetLength(c, measurementVectorLength);

  for ( unsigned int i = 0; i < measurementVectorLength; ++i )
    {
    l[i] = m_L[i];
    a[i] = m_A[i];
    c[i] = m_C[i];
    }


  m_GeneralizedGammaMembershipFunction->SetL(l);
  m_GeneralizedGammaMembershipFunction->SetA(a);
  m_GeneralizedGammaMembershipFunction->SetC(c);

}

template< class TSample >
void
GeneralizedGammaMixtureModelComponent< TSample >
::SetParameters(const ParametersType & parameters)
{
  Superclass::SetParameters(parameters);

  unsigned int paramIndex = 0;
  bool changed = false;

  unsigned int i;

  MeasurementVectorSizeType measurementVectorSize =
  this->GetSample()->GetMeasurementVectorSize();


  // L parameter
  for ( i = 0; i < measurementVectorSize; i++ )
    {
    if ( m_L[i] != parameters[paramIndex] )
      {
      m_L[i] = parameters[paramIndex];

      m_LEstimate[i] = parameters[paramIndex];

      changed = true;
      }

    ++paramIndex;
    }

  ParameterVectorRealType l;
  NumericTraits<ParameterVectorRealType>::SetLength(l, measurementVectorSize);
  for ( i = 0; i < measurementVectorSize; ++i )
    {
    l[i] = m_L[i];
    }

  m_GeneralizedGammaMembershipFunction->SetL(l);


  // A parameter
  for ( i = 0; i < measurementVectorSize; i++ )
    {
    if ( m_A[i] != parameters[paramIndex] )
      {
      m_A[i] = parameters[paramIndex];

      m_AEstimate[i] = parameters[paramIndex];

      changed = true;
      }

    ++paramIndex;
    }

  ParameterVectorRealType a;
  NumericTraits<ParameterVectorRealType>::SetLength(a, measurementVectorSize);
  for ( i = 0; i < measurementVectorSize; ++i )
    {
    a[i] = m_A[i];
    }

  m_GeneralizedGammaMembershipFunction->SetA(a);


  // C parameter
  for ( i = 0; i < measurementVectorSize; i++ )
    {
    if ( m_C[i] != parameters[paramIndex] )
      {
      m_C[i] = parameters[paramIndex];

      m_CEstimate[i] = parameters[paramIndex];

      changed = true;
      }

    ++paramIndex;
    }

  ParameterVectorRealType c;
  NumericTraits<ParameterVectorRealType>::SetLength(c, measurementVectorSize);
  for ( i = 0; i < measurementVectorSize; ++i )
    {
    c[i] = m_C[i];
    }

  m_GeneralizedGammaMembershipFunction->SetC(c);


  this->AreParametersModified(changed);
}

template< class TSample >
double
GeneralizedGammaMixtureModelComponent< TSample >
::CalculateParametersChange()
{
  unsigned int i;


  double                    temp;
  double                    changes = 0.0;
  MeasurementVectorSizeType measurementVectorSize =
  this->GetSample()->GetMeasurementVectorSize();

  for ( i = 0; i < measurementVectorSize; i++ )
    {
    temp = m_L[i] - m_LEstimate[i];
    changes += temp * temp;

    temp = m_A[i] - m_AEstimate[i];
    changes += temp * temp;

    temp = m_C[i] - m_CEstimate[i];
    changes += temp * temp;
    }


  changes = vcl_sqrt(changes);
  return changes;
}

template< class TSample >
void
GeneralizedGammaMixtureModelComponent< TSample >
::GenerateData()
{

  MeasurementVectorSizeType measurementVectorSize =
  this->GetSample()->GetMeasurementVectorSize();

  this->AreParametersModified(false);


  typename TSample::MeasurementVectorType measurements;


  // Parameter estimation (L, A and C)
  this->CalculateParameters();

  MeasurementVectorSizeType   i;
  double         changes;
  bool           changed = false;
  ParametersType parameters = this->GetFullParameters(); // parameters=[l,a,c]
  MeasurementVectorSizeType            paramIndex  = 0;



  for ( i = 0; i < measurementVectorSize; i++ )
    {
    changes = vnl_math_abs( m_L[i] - m_LEstimate[i] );

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
      m_L[i] = m_LEstimate[i];
      }


    for ( paramIndex = 0; paramIndex < measurementVectorSize; paramIndex++ )
      {
      parameters[paramIndex] = m_LEstimate[paramIndex];
      }
    this->AreParametersModified(true);
    }
  else
    {
    paramIndex = measurementVectorSize;
    }

  for ( i = 0; i < measurementVectorSize; i++ )
    {
    changes = vnl_math_abs( m_A[i] - m_AEstimate[i] );

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
      m_A[i] = m_AEstimate[i];
      }


    for ( i = 0; i < measurementVectorSize; i++ )
      {
      parameters[paramIndex] = m_AEstimate[i];
      ++paramIndex;
      }
    this->AreParametersModified(true);
    }
  else
    {
    paramIndex += measurementVectorSize;
    }



  for ( i = 0; i < measurementVectorSize; i++ )
    {
    changes = vnl_math_abs( m_C[i] - m_CEstimate[i] );

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
      m_C[i] = m_CEstimate[i];
      }


    for ( i = 0; i < measurementVectorSize; i++ )
      {
      parameters[paramIndex] = m_CEstimate[i];
      ++paramIndex;
      }
    this->AreParametersModified(true);
    }



  //THIS IS NEEDED TO update m_L, m_A and m_C. SHOULD BE REMOVED


  paramIndex = 0;
  for ( i = 0; i < measurementVectorSize; i++ )
    {
    m_L[i] = parameters[paramIndex];
    ++paramIndex;
    }

  ParameterVectorRealType l;
  NumericTraits<ParameterVectorRealType>::SetLength(l, measurementVectorSize);

  for ( i = 0; i < measurementVectorSize; ++i )
    {
    l[i] = m_L[i];
    }
  m_GeneralizedGammaMembershipFunction->SetL(l);





  for ( i = 0; i < measurementVectorSize; i++ )
    {
    m_A[i] = parameters[paramIndex];
    ++paramIndex;
    }

  ParameterVectorRealType a;
  NumericTraits<ParameterVectorRealType>::SetLength(a, measurementVectorSize);

  for ( i = 0; i < measurementVectorSize; ++i )
    {
    a[i] = m_A[i];
    }
  m_GeneralizedGammaMembershipFunction->SetA(a);


  for ( i = 0; i < measurementVectorSize; i++ )
    {
    m_C[i] = parameters[paramIndex];
    ++paramIndex;
    }

  ParameterVectorRealType c;
  NumericTraits<ParameterVectorRealType>::SetLength(c, measurementVectorSize);

  for ( i = 0; i < measurementVectorSize; ++i )
    {
    c[i] = m_C[i];
    }
  m_GeneralizedGammaMembershipFunction->SetC(c);




  Superclass::SetParameters(parameters);
}

template< class TSample >
void
GeneralizedGammaMixtureModelComponent< TSample >
::CalculateParameters()
{
  const WeightArrayType  *weights = &(this->GetWeights());
  const SampleType *sample = this->GetSample();

  MeasurementVectorSizeType measurementVectorSize =
  sample->GetMeasurementVectorSize();

  typedef GeneralizedGammaMLEFunction<SampleType, WeightArrayType, ParameterVectorRealType> GeneralizedGammaMLEFunctionType;


  for ( unsigned int dim = 0; dim < measurementVectorSize; ++dim )
    {

    // C optimization
    GeneralizedGammaMLEFunctionType func(sample, m_L, m_A, weights, dim);


      // Creamos el optimizador
    vnl_amoeba optim(func);
      // Set upt initial guest
    vnl_vector<double> c(1);
    c[0]= m_C[dim];

    optim.set_max_iterations(30);
    optim.set_x_tolerance(1e-3);
    optim.set_f_tolerance(1e-3);


//    optim.verbose = 1;


    optim.minimize(c);





      // Once the c value was estimated, we estimate the l and a parameters

    std::vector<double> mlparameters = func.Evaluate(c[0]); // parameters [l, a, c]
    m_LEstimate[dim] = mlparameters[1];
    m_AEstimate[dim] = mlparameters[2];
    m_CEstimate[dim] = c[0];
    }

}



template< class TSample >
double
GeneralizedGammaMixtureModelComponent< TSample >
::SortValue() const
{
    // The sort order is specified by the mean value:
//  E{}= 1/a * special.gamma(l + 1/c) / special.gamma(l)
  double mean = 1./m_A[0] * vnl_gamma(m_L[0] + 1./m_C[0]) / vnl_gamma(m_L[0]);
  return mean;
}

} // end of namespace Statistics
} // end of namespace itk

#endif
