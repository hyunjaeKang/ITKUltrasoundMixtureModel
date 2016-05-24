
#include <iostream>
#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkCastImageFilter.h>

#include <itkVector.h>
#include <itkListSample.h>


#include <itkTimeProbesCollectorBase.h>
#include <itkMemoryProbesCollectorBase.h>


//--- Custom Libraries

// Tissue Characterization
#include "itkTissueCharacterizationFilter.h"

// Generalized Gamma Mixture Model
#include "itkGeneralizedGammaInitHelper.h"
#include "itkGeneralizedGammaMixtureModelComponent.h"

// Gamma Mixture Model
#include "itkGammaInitHelper.h"
#include "itkGammaMixtureModelComponent.h"

// Gaussian Mixture Model
#include "itkGaussianInitHelper.h"
#include "itkGaussianMixtureModelComponent.h"


// Rayleigh Mixture Model
#include "itkRayleighInitHelper.h"
#include "itkRayleighMixtureModelComponent.h"


int main(int argc, char *argv[])
{


//  DEBUG
//  itk::MultiThreader::SetGlobalDefaultNumberOfThreads(1);



  const unsigned int ImageDimension = 3;
  typedef   unsigned char       OutputPixelType;

  typedef itk::Image< double, ImageDimension > ImageType;

  typedef itk::Image< double, ImageDimension > TInternalImage;
  typedef itk::Image<OutputPixelType,ImageDimension> TOutputImage;
  typedef itk::CastImageFilter<TInternalImage,TOutputImage> TOutputCastFilter;


  typedef itk::CastImageFilter<ImageType,TInternalImage > TInternalCastFilter;

  typedef itk::ImageFileReader<ImageType > TReader;



  std::cout << "ITK Generalized Gamma Mixture Model" << std::endl;


  ImageType::ConstPointer inputImg = NULL;
  TInternalImage::Pointer img = NULL;


  if(argc > 1)
    {
    std::string inputFilename = argv[1];

    std::cout<<"Reading: "<< inputFilename <<std::endl;

    typedef itk::ImageIOBase::IOComponentType ScalarPixelType;

    itk::ImageIOBase::Pointer imageIO =
    itk::ImageIOFactory::CreateImageIO( inputFilename.c_str(),
                                       itk::ImageIOFactory::ReadMode);


    if(!imageIO)
      std::cout<<"Could not read the image information."<<std::endl;


    imageIO->SetFileName(inputFilename);
    imageIO->ReadImageInformation();

    TReader::Pointer reader = TReader::New();

    reader->SetFileName(inputFilename);
    try
      {
      reader->Update();
      inputImg = reader->GetOutput();
      }
    catch( itk::ExceptionObject& err )
      {
      std::cerr << "Could not read the image." << std::endl;
      std::cerr << err << std::endl;
      }

    TInternalCastFilter::Pointer inputCaster = TInternalCastFilter::New();
    inputCaster->SetInput(inputImg);
    inputCaster->Update();

    img = inputCaster->GetOutput();

    }
  else
    {
    std::cout << "Image must be set" << std::endl;
    return 0;
    }

 unsigned int niteration = 30;
  if(argc > 2)
    {
    niteration = atoi(argv[2]);
    }

  unsigned int sampling = 10;
  if(argc > 3)
    {
    sampling = atoi(argv[3]);
    }

  unsigned int th = 10;
  if(argc > 4)
    {
    th = atoi(argv[4]);
    }


  unsigned int numberOfClasses = 2;

  if(argc > 5)
  {
    numberOfClasses = atoi(argv[5]);
  }



  itk::TimeProbesCollectorBase chronometer;
  itk::MemoryProbesCollectorBase memorymeter;

  memorymeter.Start("All");
  chronometer.Start( "All" );

  memorymeter.Start( "Rayleigh Tissue Characterization" );
  chronometer.Start( "Rayleigh Tissue Characterization" );




  typedef itk::TissueCharacterizationFilter<TInternalImage> TissueMapType;
  TissueMapType::Pointer lmap = NULL;

  typedef itk::ImageFileWriter<TInternalImage> WriterType;
  WriterType::Pointer writer = NULL;


  /*
   * Rayleigh Mixter Model
   */

  std::cout << "Rayleigh Mixture Model" << std::endl;

  lmap = TissueMapType::New();
  lmap->SetInput(img);
  lmap->SetIterations(niteration);
  lmap->SetMinimalParametersChange(1e-3);
  lmap->SetSampling(sampling);
  lmap->VerbosityOn();
  lmap->SetMinValue(th);



  typedef itk::Statistics::RayleighMixtureModelComponent< TissueMapType::SampleType >
  RayComponentType;
  typedef itk::Statistics::RayleighInitHelper< TissueMapType::SampleType > RayInitHelperType;

    // Create and Add the  Gamma components
  std::vector<RayComponentType::Pointer> ray_components;
  for ( unsigned int i = 0 ; i < numberOfClasses ; i++ )
    {
    ray_components.push_back(RayComponentType::New());
    }


    // Create and Add the  initializer
  lmap->SetInitHelper((RayInitHelperType::Superclass*) RayInitHelperType::New());


  for ( unsigned int i = 0 ; i < numberOfClasses ; i++ )
    {
    lmap->AddComponent( (RayComponentType::Superclass*) ray_components[i].GetPointer() );
    }


  lmap->Update();

  memorymeter.Stop( "Rayleigh Tissue Characterization" );
  chronometer.Stop( "Rayleigh Tissue Characterization" );


//  DEBUG
//  paramfinal = lmap->GetParameters();
//  propfinal = lmap->GetProportions();
//
//  functions = lmap->GetOutputMembershipFunctions();
//
//  for (unsigned int i=0; i<functions->Get().size(); i++)
//    {
//    functions->Get()[i]->Print(std::cout);
//    }

  writer = WriterType::New();

  writer->SetFileName("tissueRayleigh.mhd");
  writer->SetInput(lmap->GetOutput());
  writer->UseCompressionOn();
  writer->Update();




  memorymeter.Start( "Gamma Tissue Characterization" );
  chronometer.Start( "Gamma Tissue Characterization" );


  /*
   * Gamma Mixter Model
   */

  std::cout << "Gamma Mixture Model" << std::endl;

  lmap = TissueMapType::New();
  lmap->SetInput(img);
  lmap->SetIterations(niteration);
  lmap->SetMinimalParametersChange(1e-3);
  lmap->SetSampling(sampling);
  lmap->VerbosityOn();
  lmap->SetMinValue(th);



  typedef itk::Statistics::GammaMixtureModelComponent< TissueMapType::SampleType >
  GammaComponentType;
  typedef itk::Statistics::GammaInitHelper< TissueMapType::SampleType > GammaInitHelperType;

    // Create and Add the Gamma components
  std::vector<GammaComponentType::Pointer> gamma_components;
  for ( unsigned int i = 0 ; i < numberOfClasses ; i++ )
    {
    gamma_components.push_back(GammaComponentType::New());
    }


    // Create and Add the initializer
  lmap->SetInitHelper((GammaInitHelperType::Superclass*) GammaInitHelperType::New());


  for ( unsigned int i = 0 ; i < numberOfClasses ; i++ )
    {
    lmap->AddComponent( (GammaComponentType::Superclass*) gamma_components[i].GetPointer() );
    }


  lmap->Update();


  memorymeter.Stop( "Gamma Tissue Characterization" );
  chronometer.Stop( "Gamma Tissue Characterization" );


//  DEBUG
//  paramfinal = lmap->GetParameters();
//  propfinal = lmap->GetProportions();
//
//  functions = lmap->GetOutputMembershipFunctions();
//
//  for (unsigned int i=0; i<functions->Get().size(); i++)
//    {
//    functions->Get()[i]->Print(std::cout);
//    }

  writer = WriterType::New();

  writer->SetFileName("tissueGamma.mhd");
  writer->SetInput(lmap->GetOutput());
  writer->UseCompressionOn();
  writer->Update();




  memorymeter.Start( "Generalized Gamma Tissue Characterization" );
  chronometer.Start( "Generalized Gamma Tissue Characterization" );

    /*
     * Generalized Gamma Mixter Model
     */


  std::cout << "Generalized Gamma Mixture Model" << std::endl;

  lmap = TissueMapType::New();
  lmap->SetInput(img);
  lmap->SetIterations(niteration);
  lmap->SetMinimalParametersChange(1e-3);
  lmap->SetSampling(sampling);
  lmap->VerbosityOn();
  lmap->SetMinValue(th);



  typedef itk::Statistics::GeneralizedGammaMixtureModelComponent< TissueMapType::SampleType >
  GGammaComponentType;
  typedef itk::Statistics::GeneralizedGammaInitHelper< TissueMapType::SampleType > GGammaInitHelperType;

    // Create and Add the Generalized Gamma components
  std::vector<GGammaComponentType::Pointer> ggamma_components;
  for ( unsigned int i = 0 ; i < numberOfClasses ; i++ )
    {
    ggamma_components.push_back(GGammaComponentType::New());
    }



  // Create and Add the initializer
  lmap->SetInitHelper((GGammaInitHelperType::Superclass*) GGammaInitHelperType::New());


//   Manual parameters initialization example
//  typename TissueMapType::ParameterVectorType initParameters(2);
//  typename TissueMapType::ParametersType params(3); // [l,a,c]
//  params[0] = 1.4660115657342756;
//  params[1] = 0.703326123878936;
//  params[2] = 3.6398513230768192;
//
//  initParameters[0] = params;
//
//  typename TissueMapType::ParametersType params2(3); // [l,a,c]
//  params2[0] = 1.4690812035118113;
//  params2[1] = 0.4522280978116298;
//  params2[2] = 2.9085593819591686;
//
//  initParameters[1] = params2;
//
//  typename TissueMapType::ProportionVectorType initProportions(2);
//  initProportions[0] = 0.730982;
//  initProportions[1] = 0.269018;
//
//
//  lmap->SetInitParameters(initParameters);
//  lmap->SetInitProportions(initProportions);


  for ( unsigned int i = 0 ; i < numberOfClasses ; i++ )
    {
    lmap->AddComponent( (GGammaComponentType::Superclass*) ggamma_components[i].GetPointer() );
    }


  lmap->Update();

  memorymeter.Stop( "Generalized Gamma Tissue Characterization" );
  chronometer.Stop( "Generalized Gamma Tissue Characterization" );



//     DEBUG
//  typename TissueMapType::ParameterVectorType paramfinal = lmap->GetParameters();
//  typename TissueMapType::ProportionVectorType propfinal = lmap->GetProportions();
//
//  typename TissueMapType::MembershipFunctionVectorObjectType::ConstPointer functions =
//  lmap->GetOutputMembershipFunctions();
//
//  for (unsigned int i=0; i<functions->Get().size(); i++)
//    {
//    functions->Get()[i]->Print(std::cout);
//    }

  writer = WriterType::New();

  writer->SetFileName("tissueGGamma.mhd");
  writer->SetInput(lmap->GetOutput());
  writer->UseCompressionOn();
  writer->Update();




  memorymeter.Start( "Gaussian Tissue Characterization" );
  chronometer.Start( "Gaussian Tissue Characterization" );


  /*
   * Gaussian Mixter Model
   */

  std::cout << "Gaussian Mixture Model" << std::endl;

  lmap = TissueMapType::New();
  lmap->SetInput(img);
  lmap->SetIterations(niteration);
  lmap->SetMinimalParametersChange(1e-3);
  lmap->SetSampling(sampling);
  lmap->VerbosityOn();
  lmap->SetMinValue(th);



  typedef itk::Statistics::GaussianMixtureModelComponent< TissueMapType::SampleType >
  GaussianComponentType;
  typedef itk::Statistics::GaussianInitHelper< TissueMapType::SampleType > GaussianInitHelperType;

    // Create and Add the  Gaussian components
  std::vector<GaussianComponentType::Pointer> gaussian_components;
  for ( unsigned int i = 0 ; i < numberOfClasses ; i++ )
    {
    gaussian_components.push_back(GaussianComponentType::New());
    }


  // Create and Add the  Gaussian initializer
  lmap->SetInitHelper((GaussianInitHelperType::Superclass*) GaussianInitHelperType::New());


//  for ( unsigned int i = 0 ; i < numberOfClasses ; i++ )
//    {
//    lmap->AddComponent( (GaussianComponentType::Superclass*) gaussian_components[i].GetPointer() );
//    }


  lmap->Update();


  memorymeter.Stop( "Gaussian Tissue Characterization" );
  chronometer.Stop( "Gaussian Tissue Characterization" );



//    DEBUG
//  paramfinal = lmap->GetParameters();
//  propfinal = lmap->GetProportions();
//
//  functions = lmap->GetOutputMembershipFunctions();
//
//  for (unsigned int i=0; i<functions->Get().size(); i++)
//    {
//    functions->Get()[i]->Print(std::cout);
//    }

  writer = WriterType::New();

  writer->SetFileName("tissueGaussian.mhd");
  writer->SetInput(lmap->GetOutput());
  writer->UseCompressionOn();
  writer->Update();



  memorymeter.Stop( "All" );
  chronometer.Stop( "All" );

  // Report the time and memory taken by the registration

  chronometer.Report( std::cout<<"Time:"<<std::endl );
  memorymeter.Report( std::cout<<std::endl<<"Memory:"<<std::endl );



  std::cout<<"Bye Bye !"<<std::endl;
  return 0;
}



