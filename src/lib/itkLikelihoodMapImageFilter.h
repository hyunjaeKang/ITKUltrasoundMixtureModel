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
#ifndef __itkLikelihoodMapImageFilter_h
#define __itkLikelihoodMapImageFilter_h

#include "itkSimpleDataObjectDecorator.h"
#include "itkUnaryFunctorImageFilter.h"


#include "itkMixtureModelMembershipFunctionBase.h"



namespace itk
{

  namespace Functor
  {
  template< class TInput, class TOutput >
  class MixPDF
  {
  public:

  typedef itk::Vector< double, 1 > MeasurementVectorType;
  typedef itk::Statistics::MixtureModelMembershipFunctionBase< MeasurementVectorType >
  MembershipFunctionType;
  typedef typename MembershipFunctionType::ConstPointer   MembershipFunctionPointer;
  typedef std::vector< MembershipFunctionPointer >        MembershipFunctionVectorType;
  typedef SimpleDataObjectDecorator<
  MembershipFunctionVectorType >                        MembershipFunctionVectorObjectType;
  typedef typename
  MembershipFunctionVectorObjectType::ConstPointer MembershipFunctionVectorObjectPointer;


  typedef itk::Array< double > ProportionVectorType;


  itkTypeMacro(MixPDF, Object);

  MixPDF(){}

  ~MixPDF(){}


  void SetProportions(ProportionVectorType &prop)
    {
    m_Proportions = prop;
    }

  
  void SetMembershipFunctionsObject(MembershipFunctionVectorObjectPointer membershipFunctions)
    {
    m_MembershipFunctionsObject = membershipFunctions;
    }


  
  bool operator!=(const MixPDF &) const
    {
    return false;
    }

  bool operator==(const MixPDF & other) const
    {
    return !( *this != other );
    }

  inline TOutput operator()(const TInput & value) const
    {
    if(value > 0)
      {
      MembershipFunctionType::MeasurementVectorRealType mv;
      mv[0] = value;



      const MembershipFunctionVectorType & membershipFunctionsVector =
      m_MembershipFunctionsObject->Get();

      double pdf = 0;
      double lmix = 0;
      unsigned int nClasses = membershipFunctionsVector.size();

      for (unsigned int i=0; i < nClasses-1; i++)
      {
        pdf = membershipFunctionsVector[i]->Evaluate(mv);
        lmix += m_Proportions[i] * pdf;
      }


      unsigned int t_id = nClasses-1; // tissue id class
        pdf = membershipFunctionsVector[t_id]->Evaluate(mv);
      double ltissue = m_Proportions[t_id] * pdf;
      lmix += ltissue;

      double likelihood = ltissue / lmix;

      return static_cast< TOutput >(likelihood);
      }
    else
      return NumericTraits< TOutput >::Zero;


    }


  private:
  ProportionVectorType m_Proportions;
  MembershipFunctionVectorObjectPointer  m_MembershipFunctionsObject;

  };
  }


  template< class  TInputImage >
  class ITK_EXPORT LikelihoodMapImageFilter:
  public ImageToImageFilter< TInputImage,
  Image< double, TInputImage::ImageDimension > >
  {
  public:
  /** Standard class typedefs. */
  typedef LikelihoodMapImageFilter Self;
  typedef ImageToImageFilter< TInputImage,
  Image< double, TInputImage::ImageDimension > >     Superclass;

  typedef SmartPointer< Self >       Pointer;
  typedef SmartPointer< const Self > ConstPointer;

  typedef typename Superclass::InputImageType  InputImageType;
  typedef typename Superclass::OutputImageType OutputImageType;
  typedef typename InputImageType::PixelType   InputPixelType;
  typedef typename OutputImageType::PixelType  OutputPixelType;

  /** Image dimension = 3. */
  itkStaticConstMacro(ImageDimension, unsigned int,
                      InputImageType ::ImageDimension);
  itkStaticConstMacro(InputPixelDimension, unsigned int,
                      1);

  /** Run-time type information (and related methods).   */
  itkTypeMacro(LikelihoodMapImageFilter, ImageToImageFilter);

  /** Method for creation through the object factory. */
  itkNewMacro(Self);


  typedef itk::UnaryFunctorImageFilter<InputImageType,OutputImageType,
  Functor::MixPDF<
  typename InputImageType::PixelType,
  typename OutputImageType::PixelType> > UnaryFunctorType;


    // Vector 1-dimensional
  typedef itk::Vector< double, 1 > MeasurementVectorType;
  typedef itk::Statistics::MixtureModelMembershipFunctionBase< MeasurementVectorType >
  MembershipFunctionType;
  typedef typename MembershipFunctionType::ConstPointer   MembershipFunctionPointer;
  typedef std::vector< MembershipFunctionPointer >        MembershipFunctionVectorType;
  typedef SimpleDataObjectDecorator<
  MembershipFunctionVectorType >                        MembershipFunctionVectorObjectType;
  typedef typename
  MembershipFunctionVectorObjectType::ConstPointer MembershipFunctionVectorObjectPointer;


  typedef itk::Array< double > ProportionVectorType;

  itkSetMacro(Proportions, ProportionVectorType);
  itkGetConstMacro(Proportions, ProportionVectorType);

  itkSetMacro(MembershipFunctionsObject, MembershipFunctionVectorObjectPointer);
  itkGetMacro(MembershipFunctionsObject, MembershipFunctionVectorObjectPointer);


  protected:
  LikelihoodMapImageFilter(){}
  ~LikelihoodMapImageFilter() {}
  void PrintSelf(std::ostream & os, Indent indent) const;

  /** Generate Data */
  void GenerateData(void);

  
  private:

  LikelihoodMapImageFilter(const Self &); //purposely not
                                                            // implemented
  void operator=(const Self &);                          //purposely not


  MembershipFunctionVectorObjectPointer  m_MembershipFunctionsObject;
  ProportionVectorType m_Proportions;//size=nClasses
  };

} // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkLikelihoodMapImageFilter.hxx"
#endif

#endif