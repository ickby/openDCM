

include_directories(${CMAKE_SOURCE_DIR}/src)
link_directories ( ${Boost_LIBRARY_DIRS} )

set(framework_SRC framework.cpp)
set(test_SRC  clustergraph.cpp
	      module3d.cpp)

add_library(framework ${framework_SRC})

add_executable(solvertest ${test_SRC})
target_link_libraries(solvertest framework)

