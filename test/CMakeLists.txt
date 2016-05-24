FIND_PACKAGE(Boost)
IF (Boost_FOUND)
  message(Boost_INCLUDE_DIR ${Boost_INCLUDE_DIR})
  INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIR})
ENDIF(Boost_FOUND)

itk_module_test()

set(ITKUltrasoundMixtureModelTests
    MixtureModel.cxx
    )

CreateTestDriver(ITKUltrasoundMixtureModel "${ITKUltrasoundMixtureModel_LIBRARIES}" "${ITKUltrasoundMixtureModelTests}")

itk_add_test(NAME CreateTissueMixtureModelSynthetic
  COMMAND ITKUltrasoundMixtureModelTestDriver
  CreateTissueMixtureModelSynthetic DATA{${ITK_DATA_ROOT}/Input/cthead1.png}
  )

