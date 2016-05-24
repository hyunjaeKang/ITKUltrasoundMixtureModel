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
#ifndef __itkTissueCharacterizationFilter_hxx
#define __itkTissueCharacterizationFilter_hxx

#include "itkTissueCharacterizationFilter.h"

#include <itkCastImageFilter.h>
#include "itkVector.h"

#include "itkImageRegionConstIterator.h"
#include "itkImageRegionIterator.h"
#include "itkCastImageFilter.h"



#include "itkLikelihoodMapImageFilter.h"


#include <algorithm>    // std::sort

namespace itk
{

  template< class TInputImage>
  TissueCharacterizationFilter< TInputImage >::TissueCharacterizationFilter():
  m_MinValue(1), m_BloodCharacterisation(false), m_Verbosity(false),
  m_Sampling(5), m_Iterations(400), m_MinimalParametersChange(1.0e-06)
  {
  m_Estimator = EstimatorType::New();
  m_InitHelper = NULL;
  m_MembershipFunctionsObject = MembershipFunctionVectorObjectType::New();
  }


  template< class TInputImage>
  int
  TissueCharacterizationFilter< TInputImage >::AddComponent(ComponentType *component)
  {
  return m_Estimator->AddComponent(component);
  }

  template< class TInputImage>
  void
  TissueCharacterizationFilter< TInputImage >::GenerateData()
  {
  itkDebugMacro(<< "GeneralizedGammaTissueCharacterization_hxx generating data ");


  typename TInputImage::ConstPointer inputImage = this->GetInput();

  double nClasses = m_Estimator->GetNumberOfComponents();


  ParameterVectorType finalParameter( nClasses );

  ProportionVectorType finalProportion( nClasses );


  finalProportion.Fill(- itk::NumericTraits< double >::One);


  /*
   * 1) Create the sample
   */

  typename SampleType::Pointer sample = SampleType::New();
  sample->SetMeasurementVectorSize( Self::MeasurementVectorSize ); // length of measurement vectors
                                                             // in the sample. 1-dimension.


  typedef itk::ImageRegionConstIterator< TInputImage > ConstIteratorType;
  ConstIteratorType it( inputImage, inputImage->GetRequestedRegion() );


  it.GoToBegin();

  while ( !it.IsAtEnd())
    {
    MeasurementVectorType mv;
    mv[0] = it.Get();

    if(mv[0]>=m_MinValue)
      {
      sample->PushBack( mv );
      }

    for (unsigned int i=0; !it.IsAtEnd() && i<m_Sampling;i++,++it);

    }


  /*
   * 2) Estimate the initial parameters
   */


  if (m_InitHelper.IsNotNull())
    {
    m_InitHelper->SetSample(sample);
    m_InitHelper->SetNumberOfClasses(nClasses);
    m_InitHelper->SetMaximumIteration(m_Iterations);
    if (m_Verbosity)
      m_InitHelper->SetVerbosityOn();

    m_InitHelper->Update();

    m_InitProportions = m_InitHelper->GetProportions();
    m_InitParameters = m_InitHelper->GetParameters();
    }



  double numberOfParameters = m_InitParameters[0].GetSize();

  if(m_Verbosity)
    {
    std::cout<<std::endl
    <<"-----\tInitial parameters \t-----"<<std::endl;

    for(unsigned int i=0; i<nClasses;i++)
      {
      ParametersType param = m_InitParameters[i];
      std::cout<<"Cluster [ "<<i<<"]"<<std::endl;
      std::cout<<"\tParameters: "<<std::endl;

      std::cout<<"\t\t[";
      for (unsigned int j=0; j< numberOfParameters -1; j++)
        std::cout<<param[j]<<",";
      std::cout<<param[numberOfParameters -1]<<"]"<<std::endl;

      std::cout<<"\tProportion: "<<m_InitProportions[i]<<std::endl;
      }
    }



  /*
   * 3) Expectation-Maximization algorithm
   */
  for ( unsigned int i = 0 ; i < nClasses ; i++ )
    {
    (m_Estimator->GetComponent(i))->SetSample( sample );
    (m_Estimator->GetComponent(i))->SetParameters( m_InitParameters[i] );
    (m_Estimator->GetComponent(i))->SetMinimalParametersChange( m_MinimalParametersChange );
    }



  m_Estimator->SetSample( sample );
  m_Estimator->SetMaximumIteration( m_Iterations );
  m_Estimator->SetInitialProportions( m_InitProportions );


  typedef CommandIterationUpdate<SampleType> TObserver;

  typename TObserver::Pointer observer = TObserver::New();
  m_Estimator->AddObserver( itk::IterationEvent(), observer );

  if(m_Verbosity)
    std::cout<<std::endl
    <<"EM-Algorithm"<<std::endl;


  m_Estimator->Update();

  /*
   * 4) Identified the tissue
   */

  const ComponentVectorType components = m_Estimator->GetComponents();

  std::vector< std::pair< double, uint > > idx(nClasses);

  for(unsigned int i=0; i<nClasses;i++)
    {
    double sort_val  = components[i]->SortValue();
    idx[i].first = sort_val;
    idx[i].second = i;

    }



  if (m_BloodCharacterisation)
    std::sort(idx.begin(),idx.end(),TissueCharacterizationFilter< TInputImage >::compare_descending);
  else
    std::sort(idx.begin(),idx.end(),TissueCharacterizationFilter< TInputImage >::compare_ascending);




  MembershipFunctionVectorType & membershipFunctionsVector = m_MembershipFunctionsObject->Get();


  for(unsigned int i=0; i<nClasses;i++)
    {

    ParametersType param = (components[idx[i].second])->GetFullParameters();
    finalParameter[i] = param;
    finalProportion[i] = m_Estimator->GetProportions()[idx[i].second];

    LightObject::Pointer loPtr = (components[idx[i].second])->GetMembershipFunction()->Clone();

    MembershipFunctionType * membershipFunction =
    dynamic_cast<MembershipFunctionType *>(loPtr.GetPointer());

    membershipFunctionsVector.push_back( membershipFunction);
    }



  if(m_Verbosity)
    {
    std::cout<<"-----\tFinal parameters\t-----"<<std::endl;

    for ( unsigned int i = 0 ; i < nClasses ; i++ )
      {
      std::cout << "Cluster[" << i << "]" << std::endl;
      std::cout << "    Parameters:" << std::endl;
      std::cout << "         " << finalParameter[i]<< std::endl;
      std::cout << "    Proportion: ";
      std::cout << "         " << finalProportion[i] << std::endl;
      }


      // Termination code


    std::cout<<std::endl
    <<"Termination code: ";
    switch (m_Estimator->GetTerminationCode())
      {
        case EstimatorType::CONVERGED:
        std::cout<<"Converged"<<std::endl;
        break;
        case EstimatorType::NOT_CONVERGED:
        std::cout<<"Max Iteration achieve (Not Converged)"<<std::endl;
        break;
      }


    std::cout << "Building the likelihood map" << std::endl;
    }

  m_Parameters = finalParameter;
  m_Proportions = finalProportion;



  /*
   * 5) Create the likelihood map for the tissue
   */

  typedef itk::LikelihoodMapImageFilter<TInputImage> LikelihoodMapType;
  typename LikelihoodMapType::Pointer lmapFilter = LikelihoodMapType::New();

  lmapFilter->SetInput(inputImage);
  lmapFilter->SetProportions(m_Proportions);
  lmapFilter->SetMembershipFunctionsObject(this->GetOutputMembershipFunctions());
  lmapFilter->Update();





  typedef itk::CastImageFilter< typename LikelihoodMapType::OutputImageType,
  OutputImageType > TOutputCastFilter;

  typename TOutputCastFilter::Pointer caster = TOutputCastFilter::New();
  caster->SetInput(lmapFilter->GetOutput());
  caster->Update();


  this->GraftOutput(caster->GetOutput());


  }



template< class TInputImage>
const typename TissueCharacterizationFilter< TInputImage >::MembershipFunctionVectorObjectType *
TissueCharacterizationFilter< TInputImage >::GetOutputMembershipFunctions() const
{
  return static_cast< const MembershipFunctionVectorObjectType * >( m_MembershipFunctionsObject );
}


} // end namespace itk

#endif
