

include_directories(${CMAKE_SOURCE_DIR})
link_directories ( ${Boost_LIBRARY_DIRS} )

set(framework_SRC framework.cpp)
set(test_SRC  define.cpp
	      #kernel.cpp
	      #clustergraph.cpp
	      #geometry.cpp
	      #sys_and_obj.cpp
	      #clustermath.cpp
	      #module3d.cpp
	      #dof.cpp
	      modulepart.cpp
)

if(LOGGING)
set(log_LIB   boost_log
	      rt
	      ${Boost_SYSTEM_LIBRARY}
	      ${Boost_FILESYSTEM_LIBRARY}
	      ${Boost_THREAD_LIBRARY}
)
endif(LOGGING)

add_library(framework ${framework_SRC})

add_executable(solvertest ${test_SRC})
target_link_libraries(solvertest framework ${log_LIB})

