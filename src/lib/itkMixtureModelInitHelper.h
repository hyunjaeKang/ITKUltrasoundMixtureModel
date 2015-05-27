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
#ifndef __itkMixtureModelInitHelper_h
#define __itkMixtureModelInitHelper_h

#include "itkObject.h"

namespace itk
{
namespace Statistics
{


/** \class MixtureModelInitHelper
 * \brief MixtureModelInitHelper is the base class of the mixture model 
 * initialization strategy (eq. 2) the subclasses should implement the particular
 * strategy according to the membership function.
 *
 * \ingroup ITKStatistics
 */
template< class TSample >
class ITK_EXPORT MixtureModelInitHelper:public Object
  {
  public:
    /** Standard class typedef */
    typedef MixtureModelInitHelper                    Self;
    typedef Object                                       Superclass;
    typedef SmartPointer< Self >                         Pointer;
    typedef SmartPointer< const Self >                   ConstPointer;

    /** Standard macros */
    itkTypeMacro(MixtureModelInitHelper, Object);



    /** TSample template argument related typedefs */
    typedef TSample                                 SampleType;
    typedef typename TSample::MeasurementType       MeasurementType;
    typedef typename TSample::MeasurementVectorType MeasurementVectorType;


    /** Type of the array of the proportion values */
    typedef Array< double > ProportionVectorType;

    /** Type of the array of the proportion values */
    typedef Array< double > ParameterType;
    typedef std::vector< ParameterType > ParameterVectorType;



    /** Sets the target data that will be classified by this */
    void SetSample(TSample *sample) {m_Sample = sample;}

    /** Returns the target data */
    const TSample * GetConstSample() const {return m_Sample;}


    void SetNumberOfClasses(int nClasses){m_NumberOfClasses = nClasses;};
    int GetNumberOfClasses() const{return m_NumberOfClasses;};

    void SetMaximumIteration(int nIter){m_Iterations = nIter;};
    int GetMaximumIteration() const{return m_Iterations;};

    void SetNumberOfParameters(int nParameters){m_NumberOfParameters = nParameters;};
    int GetNumberOfParameters() const{return m_NumberOfParameters;};


    void SetVerbosityOn(){SetVerbosity(true);}
    void SetVerbosityOff(){SetVerbosity(false);}

    void SetVerbosity(bool v){m_Verbosity = v;}
    bool GetVerbosity() const{return m_Verbosity;}


    /** Gets the result proportion values */
    void SetProportions(ProportionVectorType & proportions){m_Proportions = proportions;};
    const ProportionVectorType & GetProportions() const;

    void SetParameters(ParameterVectorType & parameters){m_Parameters = parameters;};
    const ParameterVectorType & GetParameters() const;


    void Update();

  protected:
    MixtureModelInitHelper();
    virtual ~MixtureModelInitHelper() {}
    void PrintSelf(std::ostream & os, Indent indent) const;


    /** Starts the estimation process */
    virtual void GenerateData() = 0;

    TSample * GetSample() const {return m_Sample;}



  private:
  
    SampleType            *m_Sample;
    ProportionVectorType  m_Proportions;
    ParameterVectorType   m_Parameters;
    int                   m_NumberOfClasses;
    int                   m_NumberOfParameters;
    bool                  m_Verbosity;
    int                   m_Iterations;
    
    
  };  // end of class
} // end of namespace Statistics
} // end of namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkMixtureModelInitHelper.hxx"
#endif


#endif
