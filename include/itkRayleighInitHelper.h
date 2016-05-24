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
#ifndef __itkRayleighInitHelper_h
#define __itkRayleighInitHelper_h

#include "itkMixtureModelInitHelper.h"

namespace itk
{
namespace Statistics
{
/** \class RayleighInitHelper
 * \brief RayleighInitHelper implements the initialization strategy defined in
 * eq. 2 and 10.
 *
 * \ingroup ITKStatistics
 */
  
template< class TSample >
class ITK_EXPORT RayleighInitHelper:
  public MixtureModelInitHelper<TSample>
{
public:
  /** Standard class typedef */
  typedef RayleighInitHelper                              Self;
  typedef MixtureModelInitHelper<TSample>              Superclass;
  typedef SmartPointer< Self >                         Pointer;
  typedef SmartPointer< const Self >                   ConstPointer;

  /** Standard macros */
  itkTypeMacro(RayleighInitHelper, MixtureModelInitHelper);
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
  RayleighInitHelper();
  virtual ~RayleighInitHelper() {}
  void PrintSelf(std::ostream & os, Indent indent) const;


  /** Starts the estimation process */
  void GenerateData();



};  // end of class
} // end of namespace Statistics
} // end of namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkRayleighInitHelper.hxx"
#endif

#endif
