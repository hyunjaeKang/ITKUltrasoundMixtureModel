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
#ifndef __itkTissueCharacterizationFilter_h
#define __itkTissueCharacterizationFilter_h


#include "itkCommand.h"
#include "itkListSample.h"

// Custom libraries
#include "itkExpectationMaximizationUSMixtureModelEstimator.h"
#include "itkMixtureModelInitHelper.h"





namespace itk
{
/** \class TissueCharacterizationFilter
 * \brief TissueCharacterizationFilter encapsulate the complete procedure
 * to initialize, estimate the mixture model parameters and create the
 * probabilistic tissue characterization as follows:
 * 1) Create the samples.
 * 2) Initialize the mixture model using a k-means strategy.
 * 3) Aplicamos la estimacion utilizando el algoritmo de EM para el modelo de
 * Mesclas de Gammas Generalizadas
 * 4) Sort the component of the mixture model to detect the tissue.
 * 5) Create the likelihood probability map for characterizing the tissue.
 *
 * TissueCharacterizationFilter is a subclass of ImageToImageFilter
 * that implement
 *
 * \ingroup ITKStatistics
 */

  template< class  TInputImage>
  class ITK_EXPORT TissueCharacterizationFilter:public
  ImageToImageFilter< TInputImage, TInputImage >
  {

  public:
  /** Standard class typedefs. */
  typedef TissueCharacterizationFilter Self;
  typedef ImageToImageFilter< TInputImage, TInputImage>     Superclass;

  typedef SmartPointer< Self >       Pointer;
  typedef SmartPointer< const Self > ConstPointer;

  /** Run-time type information (and related methods).   */
  itkTypeMacro(TissueCharacterizationFilter, ImageToImageFilter);

  /** Method for creation through the object factory. */
  itkNewMacro(Self);


  typedef typename Superclass::InputImageType  InputImageType;
  typedef typename Superclass::OutputImageType OutputImageType;
  typedef typename InputImageType::PixelType   InputPixelType;
  typedef typename OutputImageType::PixelType  OutputPixelType;

  typedef itk::Image< double, InputImageType::ImageDimension > TInternalImage;

  /** Image dimension = 3. */
  itkStaticConstMacro(ImageDimension, unsigned int,
                      InputImageType::ImageDimension);


  /** MeasurementVectorSize = 1. */
  itkStaticConstMacro(MeasurementVectorSize, unsigned int, 1);
  typedef itk::Vector< double, MeasurementVectorSize > MeasurementVectorType;
  typedef itk::Statistics::ListSample< MeasurementVectorType > SampleType;


  typedef itk::Statistics::ExpectationMaximizationUSMixtureModelEstimator<
  SampleType > EstimatorType;
  typedef typename EstimatorType::ComponentType ComponentType;
  typedef typename EstimatorType::ComponentVectorType ComponentVectorType;

  typedef typename EstimatorType::MembershipFunctionType MembershipFunctionType;
  typedef typename EstimatorType::MembershipFunctionVectorObjectPointer MembershipFunctionVectorObjectPointer;
  typedef typename EstimatorType::MembershipFunctionVectorObjectType MembershipFunctionVectorObjectType;
  typedef typename EstimatorType::MembershipFunctionVectorType MembershipFunctionVectorType;



  typedef itk::Statistics::MixtureModelInitHelper< SampleType > InitHelperType;


  typedef typename ComponentType::ParametersType ParametersType;

  typedef std::vector<ParametersType> ParameterVectorType;
  typedef typename EstimatorType::ProportionVectorType ProportionVectorType;


  int AddComponent(ComponentType *component);

  unsigned int GetNumberOfComponents() const;





  itkBooleanMacro(Verbosity);
  itkSetMacro(Verbosity, bool);

  itkBooleanMacro(BloodCharacterisation);
  itkSetMacro(BloodCharacterisation, bool);
  itkGetMacro(BloodCharacterisation, bool);

  itkSetMacro(Sampling, unsigned int);
  itkGetMacro(Sampling, unsigned int);

  itkSetMacro(Iterations, unsigned int);
  itkGetMacro(Iterations, unsigned int);

  itkSetMacro(MinimalParametersChange, double);
  itkGetMacro(MinimalParametersChange, double);

  itkSetMacro(MinValue, double);
  itkGetMacro(MinValue, double);




  itkSetMacro(InitHelper, typename InitHelperType::Pointer);




  const MembershipFunctionVectorObjectType * GetOutputMembershipFunctions() const;

  itkGetConstMacro(Parameters, ParameterVectorType);
  itkGetConstMacro(Proportions, ProportionVectorType);



  void SetInitParameters(const ParameterVectorType param){m_InitParameters = param;};
  void SetInitProportions(const ProportionVectorType prop){m_InitProportions = prop;};

  itkGetConstMacro(InitParameters, ParameterVectorType);
  itkGetConstMacro(InitProportions, ProportionVectorType);








  // Funcion auxiliar para hacer el sort
  static bool compare_ascending ( const std::pair<double, uint>& l, const std::pair<double, uint>& r)
  { return l.first < r.first; }

  static bool compare_descending ( const std::pair<double, uint>& l, const std::pair<double, uint>& r)
  { return l.first < r.first; }



    //  The following section of code implements a Command observer
    //  used to monitor the evolution of the registration process.
    //
  template <class TSample>
  class CommandIterationUpdate : public itk::Command
    {
    public:
    typedef  CommandIterationUpdate   Self;
    typedef  itk::Command             Superclass;
    typedef itk::SmartPointer<Self>   Pointer;
    itkNewMacro( Self );

    protected:
    CommandIterationUpdate() {};

    public:

    typedef itk::Statistics::ExpectationMaximizationUSMixtureModelEstimator<
    TSample > EstimatorType;
    typedef const EstimatorType * EstimatorPointerType;

    void Execute(itk::Object *caller, const itk::EventObject & event)
    {
    Execute( (const itk::Object *)caller, event);
    }

    void Execute(const itk::Object * object, const itk::EventObject & event)
    {

    EstimatorPointerType estimator =
    dynamic_cast< EstimatorPointerType >( object );

    if(! itk::IterationEvent().CheckEvent( &event ))
      return;

    if(estimator)
      {
      EstimatorType* est = const_cast<EstimatorType*> (estimator);

      std::cout<< "Iteration: "<<est->GetCurrentIteration();

      const typename EstimatorType::ComponentVectorType components = est->GetComponents();
      std::cout<<"\t[";
      for (unsigned int i=0; i<components.size() -1; ++i)
        {
        std::cout<<components[i]->GetFullParameters()<<",";
        }
      std::cout<<components[components.size() -1]->GetFullParameters()<<"]"<<std::endl;


      }
    }

    };



  protected:
  TissueCharacterizationFilter();

  ~TissueCharacterizationFilter() {}

  /** Generate Data */
  void GenerateData(void);


  private:

  TissueCharacterizationFilter(const Self &); //purposely not
                                                            // implemented
  void operator=(const Self &);                          //purposely not


  double m_MinValue;
  bool m_BloodCharacterisation;
  bool m_Verbosity;
  unsigned int m_Sampling;
  unsigned int m_Iterations;
  ParameterVectorType m_InitParameters;
  ProportionVectorType m_InitProportions;
  double m_MinimalParametersChange;
  typename EstimatorType::Pointer m_Estimator;
  typename InitHelperType::Pointer m_InitHelper;



  ParameterVectorType m_Parameters;
  ProportionVectorType m_Proportions;

  MembershipFunctionVectorObjectPointer  m_MembershipFunctionsObject;

  };

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkTissueCharacterizationFilter.hxx"
#endif

#endif
