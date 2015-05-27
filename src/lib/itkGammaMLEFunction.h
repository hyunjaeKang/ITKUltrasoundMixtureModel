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
#ifndef __itkGammaMLEFunction_h
#define __itkGammaMLEFunction_h

#include <vnl/vnl_cost_function.h>
#include <vnl/vnl_gamma.h>

namespace itk
{
namespace Statistics
{


/** \class GammaMLEFunction
 * \brief GammaMLEFunction corresponds to the Gamma cost function required for
 * solving the expectation maximization strategy
 *
 * GammaMembershipFunction is a subclass of vnl_cost_function
 * that models the Eq. 14
 *
 * \ingroup ITKStatistics
 */
  class GammaMLEFunction: public vnl_cost_function
{
public:
  GammaMLEFunction() : m_Constant(0){}
  GammaMLEFunction(double cte) {m_Constant = cte;}

    // Method to use Boost
  double operator()(const double x)
  {
  vnl_vector<double> val(1);
  val[0] = x;
  return f(val);
  }

  double f(const vnl_vector<double>& x)
  {
    // log-verosimilitud
  double result = m_Constant - vcl_log(x[0]+vnl_math::eps) + vnl_digamma(x[0]+vnl_math::eps);
  return result;
  }


  void SetConstant(double k){m_Constant = k;}
  double GetConstant(void) {return m_Constant;}


private:
  double m_Constant;

}; // end class GammaMLEFunction

} // end of namespace Statistics
} // end of namespace itk
#endif
