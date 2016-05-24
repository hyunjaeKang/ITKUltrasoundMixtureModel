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
#ifndef __itkLikelihoodGeneralizeGammaMapImageFilter_hxx
#define __itkLikelihoodGeneralizeGammaMapImageFilter_hxx

#include "itkLikelihoodMapImageFilter.h"

namespace itk
{

template< class TInputImage >
void
LikelihoodMapImageFilter< TInputImage >
::GenerateData()
{
itkDebugMacro(<< "LikelihoodGeneralizeGammaMapImageFilter generating data ");

if(m_Proportions.Size() == 0 ||
   m_MembershipFunctionsObject->Get().size() == 0)
  itkExceptionMacro(<<"Parameters must be set before call Update");

typename UnaryFunctorType::Pointer unaryFilter = UnaryFunctorType::New();


unaryFilter->GetFunctor().SetMembershipFunctionsObject(m_MembershipFunctionsObject);
  unaryFilter->GetFunctor().SetProportions(m_Proportions);
unaryFilter->SetInput(this->GetInput());
unaryFilter->GraftOutput(this->GetOutput());
unaryFilter->Update();
this->GraftOutput(unaryFilter->GetOutput());
}

template< class TInputImage >
void
LikelihoodMapImageFilter< TInputImage >
::PrintSelf(std::ostream & os, Indent indent) const
{
Superclass::PrintSelf(os, indent);


for ( unsigned int i = 0 ; i < m_Proportions.Size() ; i++ )
  {
  os << indent << "Cluster[" << i << "]" << std::endl;
  os << indent << "    Parameters:" << std::endl;
  os << indent << "         " << m_MembershipFunctionsObject->Get()[i]<< std::endl;
  os << indent << "    Proportion: ";
  os << indent << "         " << m_Proportions[i] << std::endl;
  }
}

} // end namespace itk

#endif
