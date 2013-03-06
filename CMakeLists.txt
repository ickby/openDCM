project(solvertest)
cmake_minimum_required(VERSION 2.6.0 FATAL_ERROR)

#custom cmake modules
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH}, ${CMAKE_CURRENT_SOURCE_DIR}/cmake )

#options for build-finetuning
option(GENERATE_DOCS "Generate the documentation if doxygen is available" OFF)
option(LOGGING "Log internals to a file. Warning: serious speed tradeoff" OFF)
option(EXTERNALIZE "Explicit instantiation of templates to reduce compile time memory" ON)


find_package(Boost COMPONENTS unit_test_framework system filesystem thread REQUIRED)
find_package(Eigen3)

#check for extern templates
if(EXTERNALIZE)
  add_definitions( -DUSE_EXTERNAL )
endif(EXTERNALIZE)

if(LOGGING)
  add_definitions( -DUSE_LOGGING )
ENDIF(LOGGING)

if(GENERATE_DOCS)
  FIND_PACKAGE(Doxygen)
ENDIF(GENERATE_DOCS)


add_subdirectory(test)
add_subdirectory(doc)