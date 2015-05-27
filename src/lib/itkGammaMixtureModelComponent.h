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
#ifndef __itkGammaMixtureModelComponent_h
#define __itkGammaMixtureModelComponent_h

#include "itkMixtureModelComponentBase.h"
#include "itkGammaMembershipFunction.h"

#include "itkGammaMLEFunction.h"



namespace itk
{
namespace Statistics
{
/** \class GammaMixtureModelComponent
 * \brief is a component (derived from MixtureModelComponentBase) for
 * Gamma class. This class is used in
 * ExpectationMaximizationMixtureModelEstimator.
 *
 * On every iteration of EM estimation, this class's GenerateData
 * method is called to compute the new distribution parameters.
 *
 *
 * \sa MixtureModelComponentBase, ExpectationMaximizationMixtureModelEstimator
 * \ingroup ITKStatistics
 */

template< class TSample >
class ITK_EXPORT GammaMixtureModelComponent:
  public MixtureModelComponentBase< TSample >
{
public:
  /**Standard class typedefs. */
  typedef GammaMixtureModelComponent        Self;
  typedef MixtureModelComponentBase< TSample > Superclass;
  typedef SmartPointer< Self >                 Pointer;
  typedef SmartPointer< const Self >           ConstPointer;

  /**Standard Macros */
  itkTypeMacro(GammaMixtureModelComponent, MixtureModelComponentBase);
  itkNewMacro(Self);

  /** Typedefs from the superclass */
  typedef typename Superclass::MeasurementVectorType     MeasurementVectorType;
  typedef typename Superclass::MeasurementVectorSizeType MeasurementVectorSizeType;
  typedef typename Superclass::MembershipFunctionType    MembershipFunctionType;
  typedef typename Superclass::WeightArrayType           WeightArrayType;
  typedef typename Superclass::ParametersType            ParametersType;

  /** Type of the membership function. Gamma density function */
  typedef GammaMembershipFunction< MeasurementVectorType >
  NativeMembershipFunctionType;


  /** Sets the input sample */
  void SetSample(const TSample *sample);

  /** Sets the component's distribution parameters. */
  void SetParameters(const ParametersType & parameters);


  double SortValue() const;

protected:
  GammaMixtureModelComponent();
  virtual ~GammaMixtureModelComponent() {}
  void PrintSelf(std::ostream & os, Indent indent) const;

  /** Returns the sum of squared changes in parameters between
   * iterations */
  double CalculateParametersChange();

  /** Computes the new distribution parameters */
  void GenerateData();


  /** Estimate the new Alpha and Beta **/
  void CalculateAlphaBeta();


  


private:
  typename NativeMembershipFunctionType::Pointer m_GammaMembershipFunction;

  typename NativeMembershipFunctionType::AlphaVectorType m_Alpha;
  typename NativeMembershipFunctionType::BetaVectorType m_Beta;

  typename NativeMembershipFunctionType::AlphaVectorType m_AlphaEstimate;
  typename NativeMembershipFunctionType::BetaVectorType m_BetaEstimate;

};  // end of class


  
} // end of namespace Statistics
} // end of namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkGammaMixtureModelComponent.hxx"
#endif

#endif
