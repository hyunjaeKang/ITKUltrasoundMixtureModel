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
#ifndef __itkGeneralizedGammaInitHelper_hxx
#define __itkGeneralizedGammaInitHelper_hxx


#include "itkKdTree.h"
#include "itkWeightedCentroidKdTreeGenerator.h"
#include "itkKdTreeBasedKmeansEstimator.h"
#include "itkMinimumDecisionRule.h"
#include "itkDistanceToCentroidMembershipFunction.h"
#include "itkSampleClassifierFilter.h"

#include "itkCovarianceSampleFilter.h"

#include <vnl/algo/vnl_amoeba.h>
#include <vnl/algo/vnl_brent_minimizer.h>
#include <boost/math/tools/roots.hpp>



namespace itk
{
namespace Statistics
{

template< class TSample >
GeneralizedGammaInitHelper< TSample >::GeneralizedGammaInitHelper()
{
  this->SetNumberOfParameters(3);
  this->SetNumberOfClasses(2);
}

template< class TSample >
void
GeneralizedGammaInitHelper< TSample >
::GenerateData()
{

  const unsigned int measurementVectorSize = 1; // Vector 1-dimensional

  if(this->GetVerbosity())
    std::cout<<"Generalized Gamma Initialization: K-means"<<std::endl;

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
  initialMeans.Fill(itk::NumericTraits< double >::One);


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

    // creamos el vector de muestras por clase
  for(unsigned int i=0; i<this->GetNumberOfClasses();i++)
    {
    samplesCluster.push_back( SampleType::New() );
    (samplesCluster[i])->SetMeasurementVectorSize( 1 );
    }



    // segun la clasificacion rellenamos el vector de muestras
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

  double nSample = this->GetSample()->Size();



  for(unsigned int i=0; i<this->GetNumberOfClasses();i++)
    {
    ParameterType params( this->GetNumberOfParameters() );
      //Pesos
    proportions[i] = (samplesCluster[i])->Size() / nSample;



    // p parameter is obtained by minimizing the ggammaMLE function
    GeneralizedGammaMLEInitFunction func(samplesCluster[i]);

    vnl_amoeba optim(func);

    // Set upt initial guest
    vnl_vector<double> x(1);
    x[0]= 1;

    optim.set_max_iterations(30);
    optim.set_x_tolerance(1e-8);
    optim.set_f_tolerance(1e-8);
    optim.minimize(x);

    double c = x[0];


      // Transforma GG to Gamma  (x^c)
    typename SampleType::Pointer powsample = PowSample(samplesCluster[i], c);

      // Gamma parameters estimation
    std::vector<double> param = GamFit(powsample);

      // Transform the Gamma parameters to Generalized Gamma
    double a = vcl_pow(param[1], -1/c);
    double l = param[0];

      // params = [l, a, c]
    params[0] = l;
    params[1] = a;
    params[2] = c;
    
    parameters[i]=params;
    }


  this->SetProportions(proportions);
  this->SetParameters(parameters);
}


template< class TSample >
typename TSample::Pointer
GeneralizedGammaInitHelper< TSample >
::PowSample(TSample *sample, double p) 
{

typename SampleType::Pointer powsample = SampleType::New();


int measurementVectorSize =sample->GetMeasurementVectorSize();
powsample->SetMeasurementVectorSize( measurementVectorSize );
typename SampleType::ConstIterator it = sample->Begin();

while ( it != sample->End() )
  {
  MeasurementVectorType mv = it.GetMeasurementVector();
  MeasurementVectorType mvlog;

  for ( unsigned int j = 0 ; j < measurementVectorSize; j++ )
    {
    mvlog[j] = vcl_pow(mv[j], p);
    }


  powsample->PushBack(mvlog);
  ++it;
  }

return powsample;
}
  

template< class TSample >
typename TSample::Pointer
GeneralizedGammaInitHelper< TSample >
::LogSample(TSample *sample) 
{
typename SampleType::Pointer logsample = SampleType::New();


int measurementVectorSize =sample->GetMeasurementVectorSize();
logsample->SetMeasurementVectorSize( measurementVectorSize );
typename SampleType::ConstIterator it = sample->Begin();

while ( it != sample->End() )
  {
  MeasurementVectorType mv = it.GetMeasurementVector();
  MeasurementVectorType mvlog;

  for ( unsigned int j = 0 ; j < measurementVectorSize; j++ )
    {
    mvlog[j] = vcl_log(mv[j]);
    }


  logsample->PushBack(mvlog);
  ++it;
  }

return logsample;
}


  
template< class TSample >
std::vector<double>
GeneralizedGammaInitHelper< TSample >
::GamFit(TSample *sample)
{
  typedef typename itk::Statistics::CovarianceSampleFilter<SampleType> CovarianceFilterType;

  std::vector<double> param(2);


  typename CovarianceFilterType::Pointer cfilter = CovarianceFilterType::New();
  cfilter->SetInput(sample);
  cfilter->Update();
  MeasurementVectorType mean = cfilter->GetMean();


  typename SampleType::Pointer logsample = LogSample(sample);
  int measurementVectorSize =sample->GetMeasurementVectorSize();


  typename CovarianceFilterType::Pointer cLogfilter = CovarianceFilterType::New();
  cLogfilter->SetInput(logsample);
  cLogfilter->Update();
  MeasurementVectorType meanLog = cLogfilter->GetMean();

  MeasurementVectorType cte;

  for ( unsigned int j = 0 ; j < measurementVectorSize; j++ )
    {
    cte[j]= vcl_log(mean[j]) - meanLog[j];
    }
  


  GammaMLEFunction func(cte[0]);

  double aest = (3-cte[0] + vcl_sqrt(vcl_pow(cte[0]-3, 2) + 24*cte[0])) / (12*cte[0]);
  double xa = aest*(1-0.4);
  double xb = aest*(1+0.4);


  // Boost
  typedef std::pair<double, double> Result;
  boost::uintmax_t max_iter=500;
  boost::math::tools::eps_tolerance<double> tol(30);

  Result r = boost::math::tools::toms748_solve(func, xa, xb, tol, max_iter);

  param[0] = (r.first + r.second) /2;



//  Using VNL
//  vnl_brent_minimizer optim(func);
//  optim.set_max_function_evals(500);
//  optim.set_x_tolerance(1e-8);
//  optim.set_f_tolerance(1e-8);
//  param[0] = optim.minimize_given_bounds(xa, (xa+xb)/2, xb);

  

  MeasurementVectorType theta = mean / param[0];
  param[1] = theta[0];

    // param = [shape, scale]

  return param;
}


template< class TSample >
void
GeneralizedGammaInitHelper< TSample >
::PrintSelf(std::ostream & os, Indent indent) const
{
Superclass::PrintSelf(os, indent);
os << indent << "Generalized Init Helper: " << std::endl;
}

} // end of namespace Statistics
} // end of namespace itk

#endif