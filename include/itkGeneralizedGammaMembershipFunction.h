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
#ifndef __itkGeneralizedGammaMembershipFunction_h
#define __itkGeneralizedGammaMembershipFunction_h


#include "itkMixtureModelMembershipFunctionBase.h"

namespace itk
{
namespace Statistics
{
/** \class GeneralizedGammaMembershipFunction
 * \brief GeneralizedGammaMembershipFunction models class membership through a
 * GeneralizedGamma function.
 *
 * GeneralizedGammaMembershipFunction is a subclass of MembershipFunctionBase
 * that models class membership (or likelihood) using a 
 * Gamma function. The parameters of the Gamma are established using
 * the method SetParameters() or the methods SetL(), SetA() and SetC().
 *
 *
 * The parametrization used is defined as follows:
 * f(x|l,a,c) = abs(c)* a^(l*c)* x^(l*c-1) * exp(-(a*x)^c)/Gamma(l)
 *
 * \ingroup ITKStatistics
 */

template< class TMeasurementVector >
class ITK_EXPORT GeneralizedGammaMembershipFunction:
  public MixtureModelMembershipFunctionBase< TMeasurementVector >
{
public:
  /** Standard class typedefs */
  typedef GeneralizedGammaMembershipFunction                        Self;
  typedef MixtureModelMembershipFunctionBase< TMeasurementVector > Superclass;
  typedef SmartPointer< Self >                                     Pointer;
  typedef SmartPointer< const Self >                               ConstPointer;

  /** Standard macros */
  itkTypeMacro(GeneralizedGammaMembershipFunction, MixtureModelMembershipFunctionBase);
  itkNewMacro(Self);

  /** SmartPointer class for superclass */
  typedef typename Superclass::Pointer MembershipFunctionPointer;

  /** Typedef alias for the measurement vectors */
  typedef TMeasurementVector MeasurementVectorType;

  /** Length of each measurement vector */
  typedef typename Superclass::MeasurementVectorSizeType MeasurementVectorSizeType;

  /** Type of the mean vector. RealType on a vector-type is the same
   * vector-type but with a real element type.  */

  typedef typename Superclass::ParameterVectorRealType ParameterVectorRealType;



  void SetParameters(const ParameterVectorRealType & param);


  /** Set the mean of the GeneralizedGamma distribution. Mean is a vector type
   * similar to the measurement type but with a real element type. */
  void SetL(const ParameterVectorRealType & l);

  /** Get the mean of the GeneralizedGamma distribution. Mean is a vector type
   * similar to the measurement type but with a real element type. */
  itkGetConstReferenceMacro(L, ParameterVectorRealType);

  /** Set the mean of the GeneralizedGamma distribution. Mean is a vector type
   * similar to the measurement type but with a real element type. */
  void SetA(const ParameterVectorRealType & a);

  

  /** Get the mean of the GeneralizedGamma distribution. Mean is a vector type
   * similar to the measurement type but with a real element type. */
  itkGetConstReferenceMacro(A, ParameterVectorRealType);



  void SetC(const ParameterVectorRealType & c);
  /** Get the mean of the GeneralizedGamma distribution. Mean is a vector type
   * similar to the measurement type but with a real element type. */
  itkGetConstReferenceMacro(C, ParameterVectorRealType);

  
  /** Evaluate the probability density of a measurement vector. */
  double Evaluate(const MeasurementVectorType & measurement) const;

  /** Method to clone a membership function, i.e. create a new instance of
   * the same type of membership function and configure its ivars to
   * match. */
  virtual typename LightObject::Pointer InternalClone() const;


  /** Evaluamos la derivada de la funcion GeneralizedGamma
   */
  double DerivativeEvaluation(const MeasurementVectorType & measurement) const;

//  static const unsigned int nParameters = 3;// l, a, c
//  unsigned int GetNumberOfParameters() const {return nParameters;};

protected:
  GeneralizedGammaMembershipFunction(void);
  virtual ~GeneralizedGammaMembershipFunction(void) {}
  void PrintSelf(std::ostream & os, Indent indent) const;

private:
  GeneralizedGammaMembershipFunction(const Self &);   //purposely not implemented
  void operator=(const Self &); //purposely not implemented

  ParameterVectorRealType m_L;
  ParameterVectorRealType m_A;
  ParameterVectorRealType m_C;
  
};
} // end of namespace Statistics
} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkGeneralizedGammaMembershipFunction.hxx"
#endif

#endif
