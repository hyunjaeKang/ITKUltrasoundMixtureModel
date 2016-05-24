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
#ifndef __itkGeneralizedGammaInitHelper_h
#define __itkGeneralizedGammaInitHelper_h


#include <vnl/vnl_cost_function.h>

#include "itkMixtureModelInitHelper.h"
#include "itkGammaMLEFunction.h"


namespace itk
{
namespace Statistics
{


/** \class GeneralizedGammaInitHelper
 * \brief GeneralizedGammaInitHelper implements the initialization strategy 
 * defined in eq. 2 and 17.
 *
 * \ingroup ITKStatistics
 */
  
template< class TSample >
class ITK_EXPORT GeneralizedGammaInitHelper:
  public MixtureModelInitHelper<TSample>
{
public:
  /** Standard class typedef */
  typedef GeneralizedGammaInitHelper                    Self;
  typedef MixtureModelInitHelper<TSample>              Superclass;
  typedef SmartPointer< Self >                         Pointer;
  typedef SmartPointer< const Self >                   ConstPointer;

  /** Standard macros */
  itkTypeMacro(GeneralizedGammaInitHelper, MixtureModelInitHelper);
  itkNewMacro(Self);


  /** TSample template argument related typedefs */
  typedef TSample                                 SampleType;
  typedef typename TSample::MeasurementType       MeasurementType;
  typedef typename TSample::MeasurementVectorType MeasurementVectorType;


  /** Type of the array of the proportion values */
  typedef typename Superclass::ProportionVectorType ProportionVectorType;

  /** Type of the array of the proportion values */
  typedef typename Superclass::ParameterType ParameterType;
  typedef typename Superclass::ParameterVectorType ParameterVectorType;


protected:
  GeneralizedGammaInitHelper();
  virtual ~GeneralizedGammaInitHelper() {}

  /** Starts the estimation process */
  void GenerateData();

  void PrintSelf(std::ostream & os, Indent indent) const;


  static std::vector<double> GamFit(SampleType * sample);


  static typename SampleType::Pointer LogSample(SampleType * sample);

  static typename SampleType::Pointer PowSample(SampleType * sample, double p);


  /*
   * Generalize Gamma MLE function to minimize using vnl_amoeba 
   */
  class GeneralizedGammaMLEInitFunction: public vnl_cost_function
  {
  public:
  GeneralizedGammaMLEInitFunction() : m_Sample(NULL){}
  GeneralizedGammaMLEInitFunction(SampleType *sample) {m_Sample = sample;}

  double f(const vnl_vector<double>& x)
    {

    double c = x[0];
    
    // Transformamos la GG en una Gamma  sample^c
    typename SampleType::Pointer powsample = PowSample(m_Sample, c);

      // Estimamos los parametros de la gamma usando MLE
    std::vector<double> param = GamFit(powsample);

      // Minimizamos la log-verosimilitud
    double l = param[0];
    double a = vcl_pow(param[1], 1./c);
    
    typename SampleType::Pointer pdf = GeneralGamma(m_Sample, l, a, c);


    int measurementVectorSize = pdf->GetMeasurementVectorSize();
    MeasurementVectorType c_out;
    c_out.Fill(0);
    

    typename SampleType::ConstIterator it = pdf->Begin();

    while ( it != pdf->End() )
      {
      MeasurementVectorType mv = it.GetMeasurementVector();

      for ( unsigned int j = 0 ; j < measurementVectorSize; j++ )
        {
        double logmv = -vcl_log(mv[j]);
        c_out[j] += logmv;
        }
      ++it;
      }

    return c_out[0];
    }


    //TODO: GeneralGamma deberia usar la GeneralizedGammaMembershipFunction
  typename SampleType::Pointer GeneralGamma(SampleType *x, double l,
                                            double a, double c)
  {
  typename SampleType::Pointer pdf = SampleType::New();


  int measurementVectorSize =x->GetMeasurementVectorSize();
  pdf->SetMeasurementVectorSize( measurementVectorSize );
  typename SampleType::ConstIterator it = x->Begin();

  double cte = vcl_abs(c) * vcl_pow(a, l*c) / vnl_gamma(l);


  while ( it != x->End() )
    {
    MeasurementVectorType mv = it.GetMeasurementVector();
    MeasurementVectorType val;


    for ( unsigned int j = 0 ; j < measurementVectorSize; j++ )
      {
      val[j] = cte * vcl_pow(mv[j], l*c-1) * vcl_exp(-vcl_pow(a*mv[j], c));
      }


    pdf->PushBack(val);
    ++it;
    }
  
  return pdf;
  }


  private:
  SampleType            *m_Sample;
  
  }; // end class ggammaMLEFunction


};  // end of class
} // end of namespace Statistics
} // end of namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkGeneralizedGammaInitHelper.hxx"
#endif

#endif
