

include_directories(${CMAKE_SOURCE_DIR})
link_directories ( ${Boost_LIBRARY_DIRS} )

set(framework_SRC framework.cpp)
set(test_SRC clustergraph.cpp
	      geometry.cpp
	      sys_and_obj.cpp
	      module3d.cpp
	      kernel.cpp
	      dof.cpp
	      modulepart.cpp
)

add_library(framework ${framework_SRC})

add_executable(solvertest ${test_SRC})
target_link_libraries(solvertest framework)

