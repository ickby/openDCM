project(newsolver2)
cmake_minimum_required(VERSION 2.6.0 FATAL_ERROR)

find_package(Boost COMPONENTS unit_test_framework REQUIRED)
add_subdirectory(test)

