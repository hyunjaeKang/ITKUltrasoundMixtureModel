set(DOCUMENTATION "This project contains the source code for
    Probabilistic Tissue Characterization for Ultrasound Images based on
    the probabilistic nature of speckle.
    The Gaussian Mixture Model, the Rayleigh Mixture Model, the Gamma Mixture Model
    and the Generalized Gamma Mixture Model were implemented using the Insight Toolkit.
    The original version available at http://hdl.handle.net/10380/3517")

itk_module(ITKUltrasoundMixtureModel
  DEPENDS
    ITKCommon
    ITKImageFilterBase
    ITKStatistics
  TEST_DEPENDS
    ITKTestKernel
  EXCLUDE_FROM_DEFAULT
  DESCRIPTION
    "${DOCUMENTATION}"
  )
