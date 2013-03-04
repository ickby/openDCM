project(solvertest)
cmake_minimum_required(VERSION 2.6.0 FATAL_ERROR)

#custom cmake modules
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH}, ${CMAKE_CURRENT_SOURCE_DIR}/cmake )

#options for build-finetuning
option(GENERATE_DOCS "Generate the documentation if doxygen is available" OFF)
option(LOGGING "Log internals to a file. Warning: serious speed tradeoff" OFF)
option(EXTERNALIZE "Use external templates to reduce compile time memory (c++0x feature)" ON)


find_package(Boost COMPONENTS unit_test_framework system filesystem thread REQUIRED)
find_package(Eigen3)

#check for extern templates
include(CheckCXX11Features)
if(HAS_CXX11_EXTERN_TEMPLATES AND EXTERNALIZE)
  add_definitions( -DUSE_EXTERNAL )
endif(HAS_CXX11_EXTERN_TEMPLATES AND EXTERNALIZE)

if(LOGGING)
  add_definitions( -DUSE_LOGGING)
ENDIF(LOGGING)

if(GENERATE_DOCS)
  FIND_PACKAGE(Doxygen)
ENDIF(GENERATE_DOCS)


add_subdirectory(test)
add_subdirectory(doc)