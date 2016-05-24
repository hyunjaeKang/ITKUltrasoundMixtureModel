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
#ifndef __itkMixtureModelInitHelper_hxx
#define __itkMixtureModelInitHelper_hxx


namespace itk
{
namespace Statistics
{
template< class TSample >
MixtureModelInitHelper< TSample >::MixtureModelInitHelper()
  {
  m_Sample = NULL;
  m_NumberOfParameters = 2;
  m_NumberOfClasses = 2;
  m_Verbosity= false;
  m_Iterations = 30;
  }

template< class TSample >
const typename MixtureModelInitHelper< TSample >::ProportionVectorType &
MixtureModelInitHelper< TSample >::GetProportions() const
  {
  return m_Proportions;
  }


template< class TSample >
const typename MixtureModelInitHelper< TSample >::ParameterVectorType &
MixtureModelInitHelper< TSample >::GetParameters() const
  {
  return m_Parameters;
  }



template< class TSample >
void MixtureModelInitHelper< TSample >::Update()
  {
  this->GenerateData();
  }


template< class TSample >
void
MixtureModelInitHelper< TSample >
::PrintSelf(std::ostream & os, Indent indent) const
  {
  Superclass::PrintSelf(os, indent);

  os << indent << "Maximum Iteration: "
  << this->GetMaximumIteration() << std::endl;
  os << indent << "Sample: "
  << this->GetSample() << std::endl;
  os << indent << "Number Of Classes: "
  << this->GetNumberOfClasses() << std::endl;
  os << indent << "Number Of Parameters: "
  << this->GetNumberOfParameters() << std::endl;
  os << indent << "Verbosity: ";
  if(this->GetVerbosity())
    os << "On";
  else
    os << "Off";

  os << std::endl;

  }

} // end of namespace Statistics
} // end of namespace itk

#endif