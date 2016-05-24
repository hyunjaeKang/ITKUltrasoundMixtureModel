cmake_minimum_required(VERSION 2.9)
project(ITKUltrasoundMixtureModel)

#set(ITKUltrasoundMixtureModel_LIBRARIES ITKUltrasoundMixtureModel)

FIND_PACKAGE(Boost)
IF (Boost_FOUND)
  message(Boost_INCLUDE_DIR ${Boost_INCLUDE_DIR})
  INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIR})
ENDIF(Boost_FOUND)


if(NOT ITK_SOURCE_DIR)
  find_package(ITK REQUIRED)
  list(APPEND CMAKE_MODULE_PATH ${ITK_CMAKE_DIR})
  include(ITKModuleExternal)
else()
  itk_module_impl()
endif()
