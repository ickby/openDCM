project(newsolver2)
cmake_minimum_required(VERSION 2.6.0 FATAL_ERROR)

#options for build-finetuning
option(GENERATE_DOCS "Generate the documentation if doxygen is available" OFF)


find_package(Boost COMPONENTS unit_test_framework REQUIRED)

if(GENERATE_DOCS)
  FIND_PACKAGE(Doxygen)
ENDIF(GENERATE_DOCS)


add_subdirectory(test)
add_subdirectory(doc)
