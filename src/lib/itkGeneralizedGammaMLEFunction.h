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
#ifndef __itkGeneralizedGammaMLEFunction_h
#define __itkGeneralizedGammaMLEFunction_h

#include <vnl/vnl_cost_function.h>
#include <vnl/algo/vnl_brent_minimizer.h>
#include <boost/math/tools/roots.hpp>



#include "itkGammaMLEFunction.h"


#include <exception>
class MyException: public std::exception
{
public:
  MyException(const char *text){m_Text = text;}
  
  virtual const char* what() const throw()
  {
  return m_Text;
  }

private:
  const char *m_Text;
};



namespace itk
{
namespace Statistics
{
/** \class GeneralizedGammaMLEFunction
 * \brief GeneralizedGammaMLEFunction corresponds to the Gamma cost function 
 * required for solving the expectation maximization strategy
 *
 * GeneralizedGammaMLEFunction is a subclass of vnl_cost_function
 * that models the Eq. 18 and 19.
 *
 * This class makes use of the following internal parametrization:
 * f(x|a,v,p) = abs(p) x^(p*v-1) * exp(-(x/a)^p)/a^(p*v)*Gamma(v);
 *
 * and returns the parameters according to the following parametrization:
 * f(x|l,a,c) =  abs(c) * a^(l*c) * x^(l*c-1) *exp(-(a*x)^c)/Gamma(l);
 *
 * where:
 * [a,v,p] vs [a,l,c]
 * a = 1/a
 * v = l
 * p = c
 *
 * \ingroup ITKStatistics
 */

template< class TSample, class WeightArrayType, class ParameterVectorRealType>
  class GeneralizedGammaMLEFunction: public vnl_cost_function
{
public:
  typedef typename TSample::MeasurementVectorType     MeasurementVectorType;
  typedef typename TSample::MeasurementVectorSizeType MeasurementVectorSizeType;

  

  
  GeneralizedGammaMLEFunction() : m_Sample(NULL), m_L(1), m_A(1), m_Dim(0), m_Weights(NULL){}
  GeneralizedGammaMLEFunction(const TSample *sample, ParameterVectorRealType l,
                             ParameterVectorRealType a, const WeightArrayType *weights,
                             unsigned long int dim)
  {
  m_Sample = sample;
  m_L = l;
  m_A = a;
  m_Weights = weights;
  m_Dim = dim;
  }

    // Method to use Boost
  double operator()(const double x)
  {
  vnl_vector<double> val(1);
  val[0] = x;
  return f(val);
  }

  double f(const vnl_vector<double>& c)
  {
  std::vector<double> param = Evaluate(c[0]); // param=[value, l, a]
  return param[0];
  }


  

  std::vector<double> Evaluate(double c)
  {
  std::vector<double> param(3);
  param[0] = 0; param[1] = 0; param[2] = 0;

  

  MeasurementVectorSizeType measurementVectorSize =
  m_Sample->GetMeasurementVectorSize();

  MeasurementVectorType mv;
  

  /*
   * Internal parametrizatino
   * f(x|a,v,p) = abs(p) x^(p*v-1) * exp(-(x/a)^p)/a^(p*v)*Gamma(v);
   * Donde:
   * a = 1/a
   * v = l
   * p = c
   */

// Unused
//  double a = 1./m_A[m_Dim];
    
    
  double v = m_L[m_Dim];
  double p = c;


  // Transform the GG to a Gamma
  double alpha = v;



  typename TSample::Pointer sample = TSample::New();
  sample->SetMeasurementVectorSize( measurementVectorSize );

  typename TSample::ConstIterator iter = m_Sample->Begin();
  typename TSample::ConstIterator end =  m_Sample->End();

  
  while(iter != end)
    {
    mv = iter.GetMeasurementVector();
    mv[m_Dim] = vcl_pow(mv[m_Dim], p);

    sample->PushBack( mv );
    ++iter;
    }

  

  iter = sample->Begin();
  end =  sample->End();

  double cte1 = 0;
  double cte2 = 0;
  double alphaEstimate = 0;
  double betaEstimate = 0;


  double tmpWeights = 0;
  unsigned long int idW = 0;
  
  while ( iter != end )
    {
    mv = iter.GetMeasurementVector();

    double w = (*m_Weights)[idW];
    tmpWeights += w;

    cte1 += mv[m_Dim] * w;
    cte2 += vcl_log(mv[m_Dim] + vnl_math::eps) * w;

    ++iter;
    ++idW;
    }



  tmpWeights += vnl_math::eps;

  GammaMLEFunction func;


  // Compute Alpha
  double cte = vcl_log(cte1/tmpWeights + vnl_math::eps) - cte2/tmpWeights;

  func.SetConstant(cte);

  double upper;
  double lower;


  if(func(alpha) > 0)
    {
    upper = alpha;
    lower = 0.5 * upper;

    while(func(lower) >0)
      {
      upper = lower;
      lower = 0.5 * upper;
      if(lower < NumericTraits< double >::min())
        {
        MyException e("Unable to reach a maximum likelihood solution");
        std::cerr << e.what() << std::endl;
        throw e;
        }
      }
    }
  else
    {
    lower = alpha;
    upper = 2 * lower;
    while(func(upper) <0)
      {
      lower = upper;
      upper = 2 * lower;
      if(upper > NumericTraits< double >::max())
        {
        MyException e("Unable to reach a maximum likelihood solution");
        std::cerr << e.what() << std::endl;
        throw e;
        }
      
      }
    }


  // Using Boost
  typedef std::pair<double, double> Result;
  boost::uintmax_t max_iter=500;
  boost::math::tools::eps_tolerance<double> tol(30);

  Result r = boost::math::tools::toms748_solve(func, lower, upper, tol, max_iter);

  alphaEstimate = (r.first + r.second) /2;
  

//  // Using VNL
//  vnl_brent_minimizer optim(func);
//  optim.set_max_function_evals(100);
//  optim.set_x_tolerance(1e-5);
//  optim.set_f_tolerance(1e-5);
//  optim.set_verbose(true);
//
//  alphaEstimate = optim.minimize_given_bounds(lower, (lower+upper)/2, upper);

//  std::cout << "Return error: " << optim.get_failure_code() << std::endl;

  

  // Compute Beta
  betaEstimate = cte1/(tmpWeights * alphaEstimate + vnl_math::eps);


  // Compute the Generalized Gamma paramters [p, a , v]
  double v_final = alphaEstimate;
  double a_final = vcl_pow(betaEstimate, 1./p) + vnl_math::eps;


  double value = 0;
  
  iter = m_Sample->Begin();
  end =  m_Sample->End();
  idW = 0;

  while ( iter != end )
    {
    mv = iter.GetMeasurementVector();

    double w = (*m_Weights)[idW];

    value += w * (1. / p + p * v_final * vcl_log(mv[m_Dim] / a_final) -
                     vcl_pow(mv[m_Dim] / a_final, p) * vcl_log(mv[m_Dim] / a_final));

    ++iter;
    ++idW;
    }

  value = vcl_abs(value);
  

  // Change from internal to output parametrization (l, a)
  param[0] = value;
  param[1] = v_final;
  param[2] = 1./a_final;

  return param;
  }


private:
  const TSample           *m_Sample;
  const WeightArrayType   *m_Weights;
  ParameterVectorRealType  m_L;
  ParameterVectorRealType  m_A;
  unsigned long int        m_Dim;

}; // end class GeneralizedGammaMLEFunction

} // end of namespace Statistics
} // end of namespace itk
#endif