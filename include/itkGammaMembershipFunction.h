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
#ifndef __itkGammaMembershipFunction_h
#define __itkGammaMembershipFunction_h

#include "itkMixtureModelMembershipFunctionBase.h"

namespace itk
{
namespace Statistics
{
/** \class GammaMembershipFunction
 * \brief GammaMembershipFunction models class membership through a
 * Gamma function.
 *
 * GammaMembershipFunction is a subclass of MembershipFunctionBase
 * that models class membership (or likelihood) using a 
 * Gamma function. The parameters of the Gamma are established using 
 * the method SetParameters() or the methods SetAlpha() and SetBeta().
 *
 * \ingroup ITKStatistics
 */

template< class TMeasurementVector >
class ITK_EXPORT GammaMembershipFunction:
  public MixtureModelMembershipFunctionBase< TMeasurementVector >
{
public:
  /** Standard class typedefs */
  typedef GammaMembershipFunction                                  Self;
  typedef MixtureModelMembershipFunctionBase< TMeasurementVector > Superclass;
  typedef SmartPointer< Self >                                     Pointer;
  typedef SmartPointer< const Self >                               ConstPointer;

  /** Standard macros */
  itkTypeMacro(GammaMembershipFunction, MixtureModelMembershipFunctionBase);
  itkNewMacro(Self);

  /** SmartPointer class for superclass */
  typedef typename Superclass::Pointer MembershipFunctionPointer;

  /** Typedef alias for the measurement vectors */
  typedef TMeasurementVector MeasurementVectorType;

  /** Length of each measurement vector */
  typedef typename Superclass::MeasurementVectorSizeType MeasurementVectorSizeType;

  /** Type of the mean vector. RealType on a vector-type is the same
   * vector-type but with a real element type.  */
  typedef typename itk::NumericTraits< MeasurementVectorType >::RealType MeasurementVectorRealType;
  typedef MeasurementVectorRealType AlphaVectorType;
  typedef MeasurementVectorRealType BetaVectorType;


  typedef typename Superclass::ParameterVectorRealType ParameterVectorRealType;


  void SetParameters(const ParameterVectorRealType & param);

  /** Set the mean of the Gamma distribution. Mean is a vector type
   * similar to the measurement type but with a real element type. */
  void SetAlpha(const AlphaVectorType & alpha);

  /** Get the mean of the Gamma distribution. Mean is a vector type
   * similar to the measurement type but with a real element type. */
  itkGetConstReferenceMacro(Alpha, AlphaVectorType);

  /** Set the mean of the Gamma distribution. Mean is a vector type
   * similar to the measurement type but with a real element type. */
  void SetBeta(const BetaVectorType & beta);

  /** Get the mean of the Gamma distribution. Mean is a vector type
   * similar to the measurement type but with a real element type. */
  itkGetConstReferenceMacro(Beta, BetaVectorType);

  
  /** Evaluate the probability density of a measurement vector. */
  double Evaluate(const MeasurementVectorType & measurement) const;

  /** Method to clone a membership function, i.e. create a new instance of
   * the same type of membership function and configure its ivars to
   * match. */
  virtual typename LightObject::Pointer InternalClone() const;


  /** Evaluamos la derivada de la funcion gamma
   */
  double DerivativeEvaluation(const MeasurementVectorType & measurement) const;

//  static const unsigned int nParametersByClass = 2;// alpha, beta

protected:
  GammaMembershipFunction(void);
  virtual ~GammaMembershipFunction(void) {}
  void PrintSelf(std::ostream & os, Indent indent) const;

private:
  GammaMembershipFunction(const Self &);   //purposely not implemented
  void operator=(const Self &); //purposely not implemented

  AlphaVectorType m_Alpha;
  BetaVectorType  m_Beta;
  
};
} // end of namespace Statistics
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkGammaMembershipFunction.hxx"
#endif

#endif
