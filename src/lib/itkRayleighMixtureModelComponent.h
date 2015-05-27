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
#ifndef __itkRayleighMixtureModelComponent_h
#define __itkRayleighMixtureModelComponent_h

#include "itkMixtureModelComponentBase.h"
#include "itkRayleighMembershipFunction.h"


namespace itk
{
namespace Statistics
{
/** \class RayleighMixtureModelComponent
 * \brief is a component (derived from MixtureModelComponentBase) for
 * Rayleigh class. This class is used in
 * ExpectationMaximizationMixtureModelEstimator.
 *
 * On every iteration of EM estimation, this class's GenerateData
 * method is called to compute the new distribution parameters.
 *
 * \sa MixtureModelComponentBase, ExpectationMaximizationMixtureModelEstimator
 * \ingroup ITKStatistics
 */

template< class TSample >
class ITK_EXPORT RayleighMixtureModelComponent:
  public MixtureModelComponentBase< TSample >
{
public:
  /**Standard class typedefs. */
  typedef RayleighMixtureModelComponent        Self;
  typedef MixtureModelComponentBase< TSample > Superclass;
  typedef SmartPointer< Self >                 Pointer;
  typedef SmartPointer< const Self >           ConstPointer;

  /**Standard Macros */
  itkTypeMacro(RayleighMixtureModelComponent, MixtureModelComponentBase);
  itkNewMacro(Self);

  /** Typedefs from the superclass */
  typedef typename Superclass::MeasurementVectorType     MeasurementVectorType;
  typedef typename Superclass::MeasurementVectorSizeType MeasurementVectorSizeType;
  typedef typename Superclass::MembershipFunctionType    MembershipFunctionType;
  typedef typename Superclass::WeightArrayType           WeightArrayType;
  typedef typename Superclass::ParametersType            ParametersType;

  /** Type of the membership function. Rayleigh density function */
  typedef RayleighMembershipFunction< MeasurementVectorType >
  NativeMembershipFunctionType;


  /** Sets the input sample */
  void SetSample(const TSample *sample);

  /** Sets the component's distribution parameters. */
  void SetParameters(const ParametersType & parameters);


  double SortValue() const;

protected:
  RayleighMixtureModelComponent();
  virtual ~RayleighMixtureModelComponent() {}
  void PrintSelf(std::ostream & os, Indent indent) const;

  /** Returns the sum of squared changes in parameters between
   * iterations */
  double CalculateParametersChange();

  /** Computes the new distribution parameters */
  void GenerateData();


  /** Estimate the new Sigma **/
  void CalculateSigma();



private:
  typename NativeMembershipFunctionType::Pointer m_RayleighMembershipFunction;
  typename NativeMembershipFunctionType::SigmaVectorType m_Sigma;
  typename NativeMembershipFunctionType::SigmaVectorType m_SigmaEstimate;


};  // end of class


  
} // end of namespace Statistics
} // end of namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkRayleighMixtureModelComponent.hxx"
#endif

#endif
