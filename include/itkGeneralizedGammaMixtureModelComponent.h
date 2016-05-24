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
#ifndef __itkGeneralizedGammaMixtureModelComponent_h
#define __itkGeneralizedGammaMixtureModelComponent_h

#include "itkUSMixtureModelComponentBase.h"
#include "itkGeneralizedGammaMembershipFunction.h"



namespace itk
{
namespace Statistics
{
/** \class GeneralizedGammaMixtureModelComponent
 * \brief is a component (derived from MixtureModelComponentBase) for
 * GeneralizedGamma class. This class is used in
 * ExpectationMaximizationMixtureModelEstimator.
 *
 * On every iteration of EM estimation, this class's GenerateData
 * method is called to compute the new distribution parameters.
 *
 *
 * This class makes use of the following parametrization:
 * f(x|l,a,c) =  abs(c) * a^(l*c) * x^(l*c-1) *exp(-(a*x)^c)/Gamma(l);
 */

template< class TSample >
class ITK_EXPORT GeneralizedGammaMixtureModelComponent:
  public USMixtureModelComponentBase< TSample >
{
public:
  /**Standard class typedefs. */
  typedef GeneralizedGammaMixtureModelComponent        Self;
  typedef USMixtureModelComponentBase< TSample > Superclass;
  typedef SmartPointer< Self >                 Pointer;
  typedef SmartPointer< const Self >           ConstPointer;

  /**Standard Macros */
  itkTypeMacro(GeneralizedGammaMixtureModelComponent, USMixtureModelComponentBase);
  itkNewMacro(Self);

  /** Typedefs from the superclass */
  typedef TSample                                 SampleType;
  typedef typename Superclass::MeasurementVectorType     MeasurementVectorType;
  typedef typename Superclass::MeasurementVectorSizeType MeasurementVectorSizeType;
  typedef typename Superclass::MembershipFunctionType    MembershipFunctionType;
  typedef typename Superclass::WeightArrayType           WeightArrayType;
  typedef typename Superclass::ParametersType            ParametersType;


  /** Type of the membership function. GeneralizedGamma density function */
  typedef GeneralizedGammaMembershipFunction< MeasurementVectorType >
  NativeMembershipFunctionType;

  typedef typename NativeMembershipFunctionType::MeasurementVectorRealType ParameterVectorRealType;



  /** Sets the input sample */
  void SetSample(const TSample *sample);

  /** Sets the component's distribution parameters. */
  void SetParameters(const ParametersType & parameters);

  double SortValue() const;


protected:
  GeneralizedGammaMixtureModelComponent();
  virtual ~GeneralizedGammaMixtureModelComponent() {}
  void PrintSelf(std::ostream & os, Indent indent) const;

  /** Returns the sum of squared changes in parameters between
   * iterations */
  double CalculateParametersChange();

  /** Computes the new distribution parameters */
  void GenerateData();


  /** Estimate the new Alpha and Beta **/
  void CalculateParameters();




private:

  typename NativeMembershipFunctionType::Pointer m_GeneralizedGammaMembershipFunction;

  ParameterVectorRealType m_L;
  ParameterVectorRealType m_A;
  ParameterVectorRealType m_C;

  ParameterVectorRealType m_LEstimate;
  ParameterVectorRealType m_AEstimate;
  ParameterVectorRealType m_CEstimate;


};  // end of class



} // end of namespace Statistics
} // end of namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkGeneralizedGammaMixtureModelComponent.hxx"
#endif

#endif
