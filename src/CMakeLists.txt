#Requerimiento de la versión de cmake
CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

#Change PROJECT_NAME to the name of your project
PROJECT(MixtureModel)

#The following lines are required to use Dart
ENABLE_TESTING()
INCLUDE(Dart)

INCLUDE_DIRECTORIES(./lib/)



FIND_PACKAGE(ITK)
IF (ITK_FOUND)
    INCLUDE_DIRECTORIES(${ITK_INCLUDE_DIR})
    INCLUDE(${ITK_USE_FILE})
ENDIF(ITK_FOUND)


SET(Libraries
    ${ITK_LIBRARIES}
)


FIND_PACKAGE(Boost)
IF (Boost_FOUND)
    INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIR})
    ADD_DEFINITIONS( "-DHAS_BOOST" )
ENDIF(Boost_FOUND)


#the following block of code is an example of how to build an executable in
#cmake.  Unmodified, it will add an executable called "MyExe" to the project.
#MyExe will be built using the files MyClass.h and MyClass.cxx, and it will
#be linked to all the libraries you specified above. 
#You can build more than one executable per project
SET(CurrentExe "MixtureModel")
SET (${CurrentExe}_SOURCE ${CurrentExe}.cxx)

ADD_EXECUTABLE(${CurrentExe} ${${CurrentExe}_SOURCE} 
  ./lib/itkMixtureModelComponentBase.h ./lib/itkMixtureModelComponentBase.hxx
  ./lib/itkMixtureModelInitHelper.h ./lib/itkMixtureModelInitHelper.hxx
  ./lib/itkMixtureModelMembershipFunctionBase.h
  ./lib/itkExpectationMaximizationMixtureModelEstimator.h ./lib/itkExpectationMaximizationMixtureModelEstimator.hxx
  ./lib/itkTissueCharacterizationFilter.h ./lib/itkTissueCharacterizationFilter.hxx
  ./lib/itkLikelihoodMapImageFilter.h ./lib/itkLikelihoodMapImageFilter.hxx
  # Generalized Gamma Mixture Model
  ./lib/itkGeneralizedGammaMixtureModelComponent.h ./lib/itkGeneralizedGammaMixtureModelComponent.hxx
  ./lib/itkGeneralizedGammaMembershipFunction.h ./lib/itkGeneralizedGammaMembershipFunction.hxx
  ./lib/itkGeneralizedGammaInitHelper.h ./lib/itkGeneralizedGammaInitHelper.hxx
  ./lib/itkGeneralizedGammaMLEFunction.h
  # Gamma Mixture Model
  ./lib/itkGammaMixtureModelComponent.h ./lib/itkGammaMixtureModelComponent.hxx 
  ./lib/itkGammaMembershipFunction.h ./lib/itkGammaMembershipFunction.hxx
  ./lib/itkGammaInitHelper.h ./lib/itkGammaInitHelper.hxx
  ./lib/itkGammaMLEFunction.h
  # Rayleigh Mixture Model
  ./lib/itkRayleighMixtureModelComponent.h ./lib/itkRayleighMixtureModelComponent.hxx 
  ./lib/itkRayleighMembershipFunction.h ./lib/itkRayleighMembershipFunction.hxx
  ./lib/itkRayleighInitHelper.h ./lib/itkRayleighInitHelper.hxx
  # Gaussian Mixture Model
  ./lib/itkGaussianMixtureModelComponent.h ./lib/itkGaussianMixtureModelComponent.hxx 
  ./lib/itkGaussianMembershipFunction.h ./lib/itkGaussianMembershipFunction.hxx
  ./lib/itkGaussianInitHelper.h ./lib/itkGaussianInitHelper.hxx
  )


TARGET_LINK_LIBRARIES(${CurrentExe} ${Libraries})


# Auxiliar test code
SET(CompareExe "ImageCompare")
ADD_EXECUTABLE(${CompareExe} ImageCompare.cxx) 
TARGET_LINK_LIBRARIES(${CompareExe} ${Libraries})

#the following line is an example of how to add a test to your project.
#Testname is the title for this particular test.  ExecutableToRun is the
#program which will be running this test.  It can either be a part of this
#project or an external executable.  After that list any args that are needed
#for this test.  Include as many tests as you like.  If your project doesn't have
#any tests you can comment out or delete the following line.

# Synthetic Image
ADD_TEST(CreateTissueMixtureModel-Synthetic
    ${CurrentExe} ${CMAKE_SOURCE_DIR}/../Data/Normal_000.mhd 700 100 10)


ADD_TEST(TestGammaMixtureModel-Synthetic
    ImageCompare ${CMAKE_SOURCE_DIR}/../Data/tissueGamma-Normal_000.mhd
    tissueGamma.mhd)

ADD_TEST(TestGGammaMixtureModel-Synthetic
    ImageCompare ${CMAKE_SOURCE_DIR}/../Data/tissueGGamma-Normal_000.mhd
    tissueGGamma.mhd)

ADD_TEST(TestGaussianMixtureModel-Synthetic
    ImageCompare ${CMAKE_SOURCE_DIR}/../Data/tissueGaussian-Normal_000.mhd 
    tissueGaussian.mhd)


ADD_TEST(TestRayleighMixtureModel-Synthetic
    ImageCompare ${CMAKE_SOURCE_DIR}/../Data/tissueRayleigh-Normal_000.mhd 
    tissueRayleigh.mhd)


# Transesophagic Image
ADD_TEST(CreateTissueMixtureModel-TTE
    ${CurrentExe} ${CMAKE_SOURCE_DIR}/../Data/TEE.mhd 700 100 10)


ADD_TEST(TestGammaMixtureModel-TTE
    ImageCompare ${CMAKE_SOURCE_DIR}/../Data/tissueGamma-TEE.mhd
    tissueGamma.mhd)

ADD_TEST(TestGGammaMixtureModel-TTE
    ImageCompare ${CMAKE_SOURCE_DIR}/../Data/tissueGGamma-TEE.mhd
    tissueGGamma.mhd)

ADD_TEST(TestGaussianMixtureModel-TTE
    ImageCompare ${CMAKE_SOURCE_DIR}/../Data/tissueGaussian-TEE.mhd
    tissueGaussian.mhd)


ADD_TEST(TestRayleighMixtureModel-TTE
    ImageCompare ${CMAKE_SOURCE_DIR}/../Data/tissueRayleigh-TEE.mhd
    tissueRayleigh.mhd)

# Transthoraxic Image
ADD_TEST(CreateTissueMixtureModel-TEE
    ${CurrentExe} ${CMAKE_SOURCE_DIR}/../Data/TTE.mhd 700 100 10)


ADD_TEST(TestGammaMixtureModel-TEE
    ImageCompare ${CMAKE_SOURCE_DIR}/../Data/tissueGamma-TTE.mhd
    tissueGamma.mhd)

ADD_TEST(TestGGammaMixtureModel-TEE
    ImageCompare ${CMAKE_SOURCE_DIR}/../Data/tissueGGamma-TTE.mhd
    tissueGGamma.mhd)

ADD_TEST(TestGaussianMixtureModel-TEE
    ImageCompare ${CMAKE_SOURCE_DIR}/../Data/tissueGaussian-TTE.mhd
    tissueGaussian.mhd)


ADD_TEST(TestRayleighMixtureModel-TEE
    ImageCompare ${CMAKE_SOURCE_DIR}/../Data/tissueRayleigh-TTE.mhd
    tissueRayleigh.mhd)
