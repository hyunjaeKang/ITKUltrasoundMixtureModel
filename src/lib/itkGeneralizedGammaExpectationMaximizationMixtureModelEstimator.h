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
#ifndef __itkGeneralizeGammaExpectationMaximizationMixtureModelEstimator_h
#define __itkGeneralizeGammaExpectationMaximizationMixtureModelEstimator_h


#include "itkExpectationMaximizationMixtureModelEstimatorBase.h"

#include "itkGeneralizeGammaMembershipFunction.h"


namespace itk
{
namespace Statistics
{

  /*
   * Se utiliza la parametrizacion interna
   * f(x|l,a,c) = abs(c)* a^(l*c)* x^(l*c-1) * exp(-(a*x)^c)/Gamma(l)
   */

template< class TSample >
class ITK_EXPORT GeneralizeGammaExpectationMaximizationMixtureModelEstimator:
  public ExpectationMaximizationMixtureModelEstimatorBase< TSample >
{
public:
  /** Standard class typedef */
  typedef GeneralizeGammaExpectationMaximizationMixtureModelEstimator Self;
  typedef ExpectationMaximizationMixtureModelEstimatorBase< TSample > Superclass;
  typedef SmartPointer< Self >                         Pointer;
  typedef SmartPointer< const Self >                   ConstPointer;

  /** Standard macros */
  itkTypeMacro(GeneralizeGammaExpectationMaximizationMixtureModelEstimator,
               ExpectationMaximizationMixtureModelEstimator);
  itkNewMacro(Self);

  /** TSample template argument related typedefs */
  typedef TSample                                 SampleType;
  typedef typename TSample::MeasurementType       MeasurementType;
  typedef typename TSample::MeasurementVectorType MeasurementVectorType;

  /** Typedef requried to generate dataobject decorated output that can
    * be plugged into SampleClassifierFilter */
  typedef GeneralizeGammaMembershipFunction< MeasurementVectorType >
  GeneralizeGammaMembershipFunctionType;

  typedef typename GeneralizeGammaMembershipFunctionType::Pointer
  GeneralizeGammaMembershipFunctionPointer;


  typedef typename GeneralizeGammaMembershipFunctionType::ParameterVectorRealType
  ParameterVectorRealType;


  typedef typename Superclass::MembershipFunctionVectorObjectType
  MembershipFunctionVectorObjectType;

  typedef typename Superclass::ComponentVectorType ComponentVectorType;


  typedef typename Superclass::MembershipFunctionVectorType
  MembershipFunctionVectorType;

  typedef typename SampleType::MeasurementVectorSizeType MeasurementVectorSizeType;



  const  ComponentVectorType & GetComponents() const;


  /** Output Membership function vector containing the membership functions with
    * the final optimized parameters */
  const MembershipFunctionVectorObjectType * GetOutput() const;

protected:
  GeneralizeGammaExpectationMaximizationMixtureModelEstimator();
  virtual ~GeneralizeGammaExpectationMaximizationMixtureModelEstimator() {}

};  // end of class
} // end of namespace Statistics
} // end of namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkGeneralizeGammaExpectationMaximizationMixtureModelEstimator.hxx"
#endif

#endif
