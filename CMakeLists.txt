project(solvertest)
cmake_minimum_required(VERSION 2.6.0 FATAL_ERROR)

#custom cmake modules
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH}, ${CMAKE_CURRENT_SOURCE_DIR}/cmake )

#options for build-finetuning
option(GENERATE_DOCS "Generate the documentation if doxygen is available" OFF)
option(LOGGING "Log internals to a file. Warning: serious speed tradeoff" OFF)

find_package(Boost 1.49.0 COMPONENTS unit_test_framework system filesystem chrono log REQUIRED)
find_package(Eigen3 3.0 REQUIRED)
find_package(TBB 4.2 REQUIRED)

#we have c++11 code here
add_definitions(-std=c++11)

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