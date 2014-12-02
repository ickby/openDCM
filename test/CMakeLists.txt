
add_definitions(-DTESTING)

include_directories(${CMAKE_SOURCE_DIR})
include_directories(${EIGEN3_INCLUDE_DIR})
include_directories(${Boost_INCLUDE_DIR})

link_directories ( ${Boost_LIBRARY_DIRS} )

set(framework_SRC framework.cpp)

set(test_SRC  define.cpp
              #scheduler.cpp
	      #kernel.cpp
	      #clustergraph.cpp
	      #system.cpp
	      analyse.cpp
	      #numeric.cpp
	      #clustermath.cpp
	      #constraints3d.cpp
	      #module3d.cpp
	      #dof.cpp
	      #modulepart.cpp
	      #modulepart2.cpp
	      #modulehl3d.cpp
	      #misc.cpp
	      #userbugs.cpp
	      #${state_SRC}	      
)

if(LOGGING)
set(log_LIB   boost_log
	      rt
	      ${Boost_SYSTEM_LIBRARY}
	      ${Boost_FILESYSTEM_LIBRARY}
	      ${Boost_THREAD_LIBRARY}
	      ${Boost_CHRONO_LIBRARY}
	      pthread
	      tbb
)
endif(LOGGING)

add_library(framework ${framework_SRC})

add_executable(solvertest ${test_SRC})

# visual studio has problems which huge cpp files
IF (MSVC)
SET_TARGET_PROPERTIES(solvertest PROPERTIES COMPILE_FLAGS "/bigobj")
ENDIF (MSVC)

target_link_libraries(solvertest framework ${log_LIB})
