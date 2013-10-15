
add_definitions(-DTESTING)

include_directories(${CMAKE_SOURCE_DIR})
include_directories(${EIGEN3_INCLUDE_DIR})
include_directories(${Boost_INCLUDE_DIR})

link_directories ( ${Boost_LIBRARY_DIRS} )

set(framework_SRC framework.cpp)

if(EXTERNALIZE)
  set(state_SRC parser.cpp
		parser_ext4.cpp
		parser_ext2.cpp
		parser_ext3.cpp
		parser_ext1.cpp
	        parser_ext5.cpp
	        parser_ext6.cpp
	        parser_ext7.cpp
	        parser_ext8.cpp
  )
else(EXTERNALIZE)
  set(state_SRC parser.cpp)
endif(EXTERNALIZE)

set(test_SRC  define.cpp
	      kernel.cpp
	      clustergraph.cpp
	      sys_and_obj.cpp
	      geometry.cpp
	      clustermath.cpp
	      constraints3d.cpp
	      module3d.cpp
	      dof.cpp
	      modulepart.cpp
	      modulehl3d.cpp
	      #${state_SRC}	      
)

if(LOGGING)
set(log_LIB   boost_log
	      rt
	      ${Boost_SYSTEM_LIBRARY}
	      ${Boost_FILESYSTEM_LIBRARY}
	      ${Boost_THREAD_LIBRARY}
	      pthread
)
endif(LOGGING)

add_library(framework ${framework_SRC})

add_executable(solvertest ${test_SRC})

# visual studio has problems which huge cpp files
IF (MSVC)
SET_TARGET_PROPERTIES(solvertest PROPERTIES COMPILE_FLAGS "/bigobj")
ENDIF (MSVC)

target_link_libraries(solvertest framework ${log_LIB})

