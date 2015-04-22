project(solvertest)
cmake_minimum_required(VERSION 2.6.0 FATAL_ERROR)

#custom cmake modules
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH}, ${CMAKE_CURRENT_SOURCE_DIR}/cmake )

#options for build-finetuning
option(GENERATE_DOCS "Generate the documentation if doxygen is available" OFF)
option(LOGGING "Log internals to a file. Warning: serious speed tradeoff" OFF)
option(EXTERNALIZE "Explicit instantiation of templates to reduce compile time memory" ON)


find_package(Boost COMPONENTS unit_test_framework system filesystem REQUIRED)
find_package(Eigen3)
find_package(TBB)

#we have c++11 code here
add_definitions(-std=c++11)

#check for extern templates
if(EXTERNALIZE)
  add_definitions( -DDCM_EXTERNAL_CORE )
  add_definitions( -DDCM_EXTERNAL_3D )
  add_definitions( -DDCM_EXTERNAL_STATE )
endif(EXTERNALIZE)

if(LOGGING)
  add_definitions( -DDCM_USE_LOGGING )
ENDIF(LOGGING)

if(GENERATE_DOCS)
  FIND_PACKAGE(Doxygen)
ENDIF(GENERATE_DOCS)

IF(CMAKE_BUILD_TYPE MATCHES Debug)
    add_definitions( -DDCM_DEBUG )
ENDIF(CMAKE_BUILD_TYPE MATCHES Debug)

add_definitions(-Wno-deprecated-register)


add_subdirectory(test)
add_subdirectory(doc)