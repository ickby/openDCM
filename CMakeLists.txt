project(newsolver2)
cmake_minimum_required(VERSION 2.6.0 FATAL_ERROR)

#custom cmake modules
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH}, ${CMAKE_CURRENT_SOURCE_DIR}/cmake )

#options for build-finetuning
option(GENERATE_DOCS "Generate the documentation if doxygen is available" OFF)


find_package(Boost COMPONENTS unit_test_framework REQUIRED)
FIND_PACKAGE(Eigen3)

if(GENERATE_DOCS)
  FIND_PACKAGE(Doxygen)
ENDIF(GENERATE_DOCS)


add_subdirectory(test)
add_subdirectory(doc)
