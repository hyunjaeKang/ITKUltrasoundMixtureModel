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
#ifndef __itkGaussianInitHelper_hxx
#define __itkGaussianInitHelper_hxx



#include "itkKdTree.h"
#include "itkWeightedCentroidKdTreeGenerator.h"
#include "itkKdTreeBasedKmeansEstimator.h"
#include "itkMinimumDecisionRule.h"
#include "itkDistanceToCentroidMembershipFunction.h"
#include "itkSampleClassifierFilter.h"

#include "itkCovarianceSampleFilter.h"


namespace itk
{
namespace Statistics
{
template< class TSample >
GaussianInitHelper< TSample >
::GaussianInitHelper()
{
  this->SetNumberOfParameters(2); // mean(1) + covariance(1) 1D
  this->SetNumberOfClasses(2);
}



template< class TSample >
void
GaussianInitHelper< TSample >
::GenerateData()
{

  const unsigned int measurementVectorSize = 1; // Vector 1-dimensional


  if(this->GetVerbosity())
    std::cout<<"Gaussian Initialization: K-means"<<std::endl;

  /*
   * k-means clustering initialization
   */

  typedef itk::Statistics::WeightedCentroidKdTreeGenerator< SampleType >
  TreeGeneratorType;
  typename TreeGeneratorType::Pointer treeGenerator = TreeGeneratorType::New();

  treeGenerator->SetSample( this->GetSample() );
  treeGenerator->SetBucketSize( 16 );
  treeGenerator->Update();

  typedef typename TreeGeneratorType::KdTreeType TreeType;
  typedef itk::Statistics::KdTreeBasedKmeansEstimator<TreeType> KmeanEstimatorType;
  typename KmeanEstimatorType::Pointer kmeanEstimator = KmeanEstimatorType::New();

  typename KmeanEstimatorType::ParametersType initialMeans(this->GetNumberOfClasses());
  initialMeans.Fill(itk::NumericTraits< double >::Zero);


  kmeanEstimator->SetParameters( initialMeans );
  kmeanEstimator->SetKdTree( treeGenerator->GetOutput() );
  kmeanEstimator->SetMaximumIteration( this->GetMaximumIteration() );
  kmeanEstimator->SetCentroidPositionChangesThreshold(0.0);
  kmeanEstimator->StartOptimization();

  typename KmeanEstimatorType::ParametersType estimatedMeans = kmeanEstimator->GetParameters();


  typedef itk::Statistics::DistanceToCentroidMembershipFunction< MeasurementVectorType >
  MembershipFunctionType;
  typedef typename MembershipFunctionType::Pointer                      MembershipFunctionPointer;

  typedef itk::Statistics::MinimumDecisionRule DecisionRuleType;

  DecisionRuleType::Pointer decisionRule = DecisionRuleType::New();

  typedef itk::Statistics::SampleClassifierFilter< SampleType > ClassifierType;
  typename ClassifierType::Pointer classifier = ClassifierType::New();

  classifier->SetDecisionRule(decisionRule);
  classifier->SetInput( this->GetSample() );
  classifier->SetNumberOfClasses( this->GetNumberOfClasses() );

  typedef typename ClassifierType::ClassLabelVectorObjectType               ClassLabelVectorObjectType;
  typedef typename ClassifierType::ClassLabelVectorType                     ClassLabelVectorType;
  typedef typename ClassifierType::MembershipFunctionVectorObjectType       MembershipFunctionVectorObjectType;
  typedef typename ClassifierType::MembershipFunctionVectorType             MembershipFunctionVectorType;

  typename ClassLabelVectorObjectType::Pointer  classLabelsObject = ClassLabelVectorObjectType::New();
  classifier->SetClassLabels( classLabelsObject );



  ClassLabelVectorType &  classLabelsVector = classLabelsObject->Get();


  /*
   * Label the cluster's samples
   */
  for (unsigned int i=0; i<this->GetNumberOfClasses();i++)
    classLabelsVector.push_back( i+1 );



  typename MembershipFunctionVectorObjectType::Pointer membershipFunctionsObject =
  MembershipFunctionVectorObjectType::New();
  classifier->SetMembershipFunctions( membershipFunctionsObject );

  MembershipFunctionVectorType &  membershipFunctionsVector = membershipFunctionsObject->Get();

  typename MembershipFunctionType::CentroidType origin( measurementVectorSize );
  int index = 0;
  for ( unsigned int i = 0 ; i < this->GetNumberOfClasses() ; i++ )
    {
    MembershipFunctionPointer membershipFunction = MembershipFunctionType::New();
    for ( unsigned int j = 0 ; j < measurementVectorSize; j++ )
      {
      origin[j] = estimatedMeans[index++];
      }
    membershipFunction->SetCentroid( origin );
    membershipFunctionsVector.push_back( membershipFunction.GetPointer() );
    }

  classifier->Update();




  /*
   *  Mean and variance are computed from each cluster
   */
  std::vector< typename SampleType::Pointer > samplesCluster;


  for(unsigned int i=0; i<this->GetNumberOfClasses();i++)
    {
    samplesCluster.push_back( SampleType::New() );
    (samplesCluster[i])->SetMeasurementVectorSize( 1 );
    }



  const typename ClassifierType::MembershipSampleType* membershipSample = classifier->GetOutput();
  typename ClassifierType::MembershipSampleType::ConstIterator it3 = membershipSample->Begin();

  while ( it3 != membershipSample->End() )
    {
    MeasurementVectorType mv = it3.GetMeasurementVector();
    unsigned int idClass = it3.GetClassLabel();

    (samplesCluster[idClass-1])->PushBack(mv);
    ++it3;
    }



  ProportionVectorType proportions(this->GetNumberOfClasses());

  ParameterVectorType parameters( this->GetNumberOfClasses() );

  typedef itk::Statistics::CovarianceSampleFilter<SampleType> CovarianceFilterType;

  double nSample = this->GetSample()->Size();


  for(unsigned int i=0; i<this->GetNumberOfClasses();i++)
    {
    ParameterType params( this->GetNumberOfParameters() );
    
    typename CovarianceFilterType::Pointer cfilter = CovarianceFilterType::New();
    cfilter->SetInput(samplesCluster[i]);
    cfilter->Update();

    proportions[i] = (samplesCluster[i])->Size() / nSample;

    MeasurementVectorType mean = cfilter->GetMean();
    typename CovarianceFilterType::MatrixType cov = cfilter->GetCovarianceMatrix();


    unsigned int paramIndex = 0;
    unsigned int ii, jj;

    for ( ii = 0; ii < measurementVectorSize; ++ii )
      {
      params[paramIndex] = mean[ii];
      ++paramIndex;
      }


    for ( ii = 0; ii < cov.Rows(); ii++ )
      {
      for ( jj = 0; jj < cov.Cols(); jj++ )
        {
        params[paramIndex] = cov.GetVnlMatrix().get(ii, jj);
        ++paramIndex;
        }
      }

    
    parameters[i]=params;

    }

  this->SetProportions(proportions);
  this->SetParameters(parameters);
  
}


template< class TSample >
void
GaussianInitHelper< TSample >
::PrintSelf(std::ostream & os, Indent indent) const
{
Superclass::PrintSelf(os, indent);
  
  os << indent << "Gamma Init Helper " << std::endl;

}

} // end of namespace Statistics
} // end of namespace itk

#endif