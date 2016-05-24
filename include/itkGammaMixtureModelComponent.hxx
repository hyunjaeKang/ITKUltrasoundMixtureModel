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
#ifndef __itkGammaMixtureModelComponent_hxx
#define __itkGammaMixtureModelComponent_hxx

#include "itkGammaMixtureModelComponent.h"

#include <iostream>
#include <vnl/algo/vnl_brent_minimizer.h>

#include <boost/math/tools/roots.hpp>



namespace itk
{
namespace Statistics
{
template< class TSample >
GammaMixtureModelComponent< TSample >
::GammaMixtureModelComponent()
{
  m_GammaMembershipFunction = NativeMembershipFunctionType::New();
  this->SetMembershipFunction( (MembershipFunctionType *)
                               m_GammaMembershipFunction.GetPointer() );
  m_Alpha.Fill(0.0);
  m_Beta.Fill(0.0);

  m_AlphaEstimate.Fill(0.0);
  m_BetaEstimate.Fill(0.0);
}

template< class TSample >
void
GammaMixtureModelComponent< TSample >
::PrintSelf(std::ostream & os, Indent indent) const
{
  Superclass::PrintSelf(os, indent);
  os << indent << "Alpha: " << m_Alpha << std::endl;
  os << indent << "Beta: " << m_Beta << std::endl;
  os << indent << "Alpha Estimated: " << m_AlphaEstimate << std::endl;
  os << indent << "Beta: Estimated" << m_BetaEstimate << std::endl;
  os << indent << "GammaMembershipFunction: " << m_GammaMembershipFunction << std::endl;
}

template< class TSample >
void
GammaMixtureModelComponent< TSample >
::SetSample(const TSample *sample)
{
  Superclass::SetSample(sample);

  const MeasurementVectorSizeType measurementVectorLength =
  sample->GetMeasurementVectorSize();

  if(measurementVectorLength != 1)
    itkExceptionMacro(<< "Measurement Vector dimension must be 1");


  m_GammaMembershipFunction->SetMeasurementVectorSize(measurementVectorLength);




  NumericTraits<MeasurementVectorType>::SetLength(m_Alpha, measurementVectorLength);
  m_Alpha.Fill(NumericTraits< double >::One);

  NumericTraits<MeasurementVectorType>::SetLength(m_Beta, measurementVectorLength);
  m_Beta.Fill(NumericTraits< double >::One);


  NumericTraits<MeasurementVectorType>::SetLength(m_AlphaEstimate, measurementVectorLength);
  m_AlphaEstimate.Fill(NumericTraits< double >::One);

  NumericTraits<MeasurementVectorType>::SetLength(m_BetaEstimate, measurementVectorLength);
  m_BetaEstimate.Fill(NumericTraits< double >::One);



  typename NativeMembershipFunctionType::AlphaVectorType alpha;
  typename NativeMembershipFunctionType::BetaVectorType beta;

  NumericTraits<typename NativeMembershipFunctionType::AlphaVectorType>::SetLength(alpha,
                                                                                  measurementVectorLength);

  NumericTraits<typename NativeMembershipFunctionType::BetaVectorType>::SetLength(beta,
                                                                                  measurementVectorLength);

  for ( unsigned int i = 0; i < measurementVectorLength; ++i )
    {
    alpha[i] = m_Alpha[i];
    beta[i] = m_Beta[i];
    }


  m_GammaMembershipFunction->SetAlpha(alpha);
  m_GammaMembershipFunction->SetBeta(beta);
}

template< class TSample >
void
GammaMixtureModelComponent< TSample >
::SetParameters(const ParametersType & parameters)
{
  Superclass::SetParameters(parameters);

  unsigned int paramIndex = 0;
  bool changed = false;

  unsigned int i;

  MeasurementVectorSizeType measurementVectorSize =
  this->GetSample()->GetMeasurementVectorSize();


  // Alpha parameter
  for ( i = 0; i < measurementVectorSize; i++ )
    {
    if ( m_Alpha[i] != parameters[paramIndex] )
      {
      m_Alpha[i] = parameters[paramIndex];

      m_AlphaEstimate[i] = parameters[paramIndex];

      changed = true;
      }

    ++paramIndex;
    }

  typename NativeMembershipFunctionType::AlphaVectorType alpha;
  NumericTraits<typename NativeMembershipFunctionType::AlphaVectorType>::SetLength(alpha,
                                                                                  measurementVectorSize);
  for ( i = 0; i < measurementVectorSize; ++i )
    {
    alpha[i] = m_Alpha[i];
    }

  m_GammaMembershipFunction->SetAlpha(alpha);



  // Beta parameter

  for ( i = 0; i < measurementVectorSize; i++ )
    {
    if ( m_Beta[i] != parameters[paramIndex] )
      {
      m_Beta[i] = parameters[paramIndex];

      m_BetaEstimate[i] = parameters[paramIndex];

      changed = true;
      }

    ++paramIndex;
    }


  typename NativeMembershipFunctionType::BetaVectorType beta;
  NumericTraits<typename NativeMembershipFunctionType::BetaVectorType>::SetLength(beta,
                                                                                  measurementVectorSize);

  for ( i = 0; i < measurementVectorSize; ++i )
    {
    beta[i] = m_Beta[i];
    }

  m_GammaMembershipFunction->SetBeta(beta);

  this->AreParametersModified(changed);
}

template< class TSample >
double
GammaMixtureModelComponent< TSample >
::CalculateParametersChange()
{
  unsigned int i;


  double                    temp;
  double                    changes = 0.0;
  MeasurementVectorSizeType measurementVectorSize =
  this->GetSample()->GetMeasurementVectorSize();

  for ( i = 0; i < measurementVectorSize; i++ )
    {
    temp = m_Alpha[i] - m_AlphaEstimate[i];
    changes += temp * temp;

    temp = m_Beta[i] - m_BetaEstimate[i];
    changes += temp * temp;
    }


  changes = vcl_sqrt(changes);
  return changes;

}

template< class TSample >
void
GammaMixtureModelComponent< TSample >
::GenerateData()
{

  MeasurementVectorSizeType measurementVectorSize =
  this->GetSample()->GetMeasurementVectorSize();

  this->AreParametersModified(false);


  typename TSample::MeasurementVectorType measurements;


  // Estimate Alpha and Beta
  this->CalculateAlphaBeta();

  MeasurementVectorSizeType   i;
  double         changes;
  bool           changed = false;
  ParametersType parameters = this->GetFullParameters(); // parameters = [alpha,beta]
  MeasurementVectorSizeType            paramIndex  = 0;


  for ( i = 0; i < measurementVectorSize; i++ )
    {
    changes = vnl_math_abs( m_Alpha[i] - m_AlphaEstimate[i] );

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
      m_Alpha[i] = m_AlphaEstimate[i];
      }


    for ( paramIndex = 0; paramIndex < measurementVectorSize; paramIndex++ )
      {
      parameters[paramIndex] = m_AlphaEstimate[paramIndex];
      }
    this->AreParametersModified(true);
    }
  else
    {
    paramIndex = measurementVectorSize;
    }


  for ( i = 0; i < measurementVectorSize; i++ )
    {
    changes = vnl_math_abs( m_Beta[i] - m_BetaEstimate[i] );

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
      m_Beta[i] = m_BetaEstimate[i];
      }


    for ( i = 0; i < measurementVectorSize; i++ )
      {
      parameters[paramIndex] = m_BetaEstimate[i];
      ++paramIndex;
      }
    this->AreParametersModified(true);
    }





  //THIS IS NEEDED TO update m_alpha and m_Beta.SHOULD BE REMOVED

  // Update Alpha
  paramIndex = 0;
  for ( i = 0; i < measurementVectorSize; i++ )
    {
    m_Alpha[i] = parameters[paramIndex];
    ++paramIndex;
    }

  typename NativeMembershipFunctionType::AlphaVectorType alpha;
  NumericTraits<typename NativeMembershipFunctionType::AlphaVectorType>::SetLength(alpha,
                                                                                  measurementVectorSize);

  for ( i = 0; i < measurementVectorSize; ++i )
    {
    alpha[i] = m_Alpha[i];
    }
  m_GammaMembershipFunction->SetAlpha(alpha);




  // Update Beta
  for ( i = 0; i < measurementVectorSize; i++ )
    {
    m_Beta[i] = parameters[paramIndex];
    ++paramIndex;
    }

  typename NativeMembershipFunctionType::BetaVectorType beta;
  NumericTraits<typename NativeMembershipFunctionType::BetaVectorType>::SetLength(beta,
                                                                                  measurementVectorSize);

  for ( i = 0; i < measurementVectorSize; ++i )
    {
    beta[i] = m_Beta[i];
    }
  m_GammaMembershipFunction->SetBeta(beta);

  Superclass::SetParameters(parameters);
}

template< class TSample >
void
GammaMixtureModelComponent< TSample >
::CalculateAlphaBeta()
{
  const WeightArrayType & weights = this->GetWeights();
  unsigned long int idW = 0;


  typename TSample::ConstIterator iter = this->GetSample()->Begin();
  typename TSample::ConstIterator end =  this->GetSample()->End();
  typename TSample::MeasurementVectorType measurements;


  MeasurementVectorSizeType measurementVectorSize =
  this->GetSample()->GetMeasurementVectorSize();



  typename NativeMembershipFunctionType::AlphaVectorType cte1;
  NumericTraits<typename NativeMembershipFunctionType::AlphaVectorType>::SetLength(cte1,
                                                                                   measurementVectorSize);
  cte1.Fill(NumericTraits< double >::Zero);

  typename NativeMembershipFunctionType::AlphaVectorType cte2;
  NumericTraits<typename NativeMembershipFunctionType::AlphaVectorType>::SetLength(cte2,
                                                                                   measurementVectorSize);
  cte2.Fill(NumericTraits< double >::Zero);



  double tmpWeights = 0;




  while ( iter != end )
    {
    measurements = iter.GetMeasurementVector();

    double w = weights[idW];
    tmpWeights += w;

    for ( unsigned int dim = 0; dim < measurementVectorSize; ++dim )
      {
      cte1[dim] += measurements[dim] * w;
      cte2[dim] += vcl_log(measurements[dim] + vnl_math::eps) * w;

      if (cte1[dim] != cte1[dim] || cte2[dim] != cte2[dim])
        {
        std::cout << "GammaMixtureModelComponent:: Nan is reached" << std::endl;
        itkExceptionMacro(<<"Nan is reached");
        }
      }


    ++iter;
    ++idW;
    }


  tmpWeights += vnl_math::eps;

  GammaMLEFunction func;


  for ( unsigned int dim = 0; dim < measurementVectorSize; ++dim )
    {
    // Compute Alpha
    double cte = vcl_log(cte1[dim]/tmpWeights + vnl_math::eps) - cte2[dim]/tmpWeights;

    func.SetConstant(cte);

    double upper;
    double lower;


    if(func(m_Alpha[dim]) > 0)
      {
      upper = m_Alpha[dim];
      lower = 0.5 * upper;

      while(func(lower) >0)
        {
        upper = lower;
        lower = 0.5 * upper;
        if(lower < NumericTraits< double >::min())
          itkExceptionMacro(<<"Unable to reach a maximum likelihood solution");
        }
      }
    else
      {
      lower = m_Alpha[dim];
      upper = 2 * lower;
      while(func(upper) <0)
        {
        lower = upper;
        upper = 2 * lower;
        if(upper > NumericTraits< double >::max())
          itkExceptionMacro(<<"Unable to reach a maximum likelihood solution");
        }
      }



    // Using Boost
    typedef std::pair<double, double> Result;
    boost::uintmax_t max_iter=500;
    boost::math::tools::eps_tolerance<double> tol(30);

    Result r = boost::math::tools::toms748_solve(func, lower, upper, tol, max_iter);

    m_AlphaEstimate[dim] = (r.first + r.second) /2;


      // There are some minimization issues with VNL
//    vnl_brent_minimizer optim(func);
//    optim.set_max_function_evals(500);
//    optim.set_x_tolerance(1e-8);
//    optim.set_f_tolerance(1e-8);
//    m_AlphaEstimate[dim] = optim.minimize_given_bounds(lower, (lower+upper)/2, upper);


    // Compute Beta
    m_BetaEstimate[dim] = cte1[dim]/(tmpWeights * m_AlphaEstimate[dim] + vnl_math::eps);

    }
}



template< class TSample >
double
GammaMixtureModelComponent< TSample >
::SortValue() const
{
    // The sort order is specified by the mean value (E{} = alpha * beta)
  double mean = m_Alpha[0] * m_Beta[0];

  return mean;
}

} // end of namespace Statistics
} // end of namespace itk

#endif
