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
#ifndef __itkRayleighMembershipFunction_h
#define __itkRayleighMembershipFunction_h

#include "itkMixtureModelMembershipFunctionBase.h"

namespace itk
{
namespace Statistics
{
/** \class RayleighMembershipFunction
 * \brief RayleighMembershipFunction models class membership through a
 * Rayleigh function.
 *
 * RayleighMembershipFunction is a subclass of MembershipFunctionBase
 * that models class membership (or likelihood) using a
 * Rayleigh function. The parameters of the Rayleigh are established using
 * the method SetParameters() or  SetSigma().
 * sharply near the mean.
 *
 * \ingroup ITKStatistics
 */

template< class TMeasurementVector >
class ITK_EXPORT RayleighMembershipFunction:
  public MixtureModelMembershipFunctionBase< TMeasurementVector >
{
public:
  /** Standard class typedefs */
  typedef RayleighMembershipFunction                                  Self;
  typedef MixtureModelMembershipFunctionBase< TMeasurementVector > Superclass;
  typedef SmartPointer< Self >                                     Pointer;
  typedef SmartPointer< const Self >                               ConstPointer;

  /** Standard macros */
  itkTypeMacro(RayleighMembershipFunction, MixtureModelMembershipFunctionBase);
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
  typedef MeasurementVectorRealType SigmaVectorType;



  typedef typename Superclass::ParameterVectorRealType ParameterVectorRealType;


  void SetParameters(const ParameterVectorRealType & param);

  /** Set the mean of the Rayleigh distribution. Mean is a vector type
   * similar to the measurement type but with a real element type. */
  void SetSigma(const SigmaVectorType & alpha);

  /** Get the mean of the Rayleigh distribution. Mean is a vector type
   * similar to the measurement type but with a real element type. */
  itkGetConstReferenceMacro(Sigma, SigmaVectorType);

  
  /** Evaluate the probability density of a measurement vector. */
  double Evaluate(const MeasurementVectorType & measurement) const;

  /** Method to clone a membership function, i.e. create a new instance of
   * the same type of membership function and configure its ivars to
   * match. */
  virtual typename LightObject::Pointer InternalClone() const;


protected:
  RayleighMembershipFunction(void);
  virtual ~RayleighMembershipFunction(void) {}
  void PrintSelf(std::ostream & os, Indent indent) const;

private:
  RayleighMembershipFunction(const Self &);   //purposely not implemented
  void operator=(const Self &); //purposely not implemented

  SigmaVectorType m_Sigma;

  
};
} // end of namespace Statistics
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkRayleighMembershipFunction.hxx"
#endif

#endif
