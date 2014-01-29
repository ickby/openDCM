
add_definitions(-DTESTING)

include_directories(${CMAKE_SOURCE_DIR})
include_directories(${EIGEN3_INCLUDE_DIR})
include_directories(${Boost_INCLUDE_DIR})

link_directories ( ${Boost_LIBRARY_DIRS} )

set(framework_SRC framework.cpp)

if(EXTERNALIZE)

  set(state_SRC parser/parser.cpp
		parser/parser_3d_ext1.cpp	
		parser/parser_3d_ext2.cpp
		parser/parser_3d_ext3.cpp
		parser/parser_ext1.cpp
		parser/parser_ext2.cpp
		parser/parser_ext3.cpp
		parser/parser_ext4.cpp
	        parser/parser_ext5.cpp
	        parser/parser_ext6.cpp
	        parser/parser_ext7.cpp
	        parser/parser_ext8.cpp
)
 set(module3d_SRC module3d/module3d.cpp
		  module3d/m3d_ext1.cpp
		  module3d/m3d_ext2.cpp
		  module3d/m3d_ext3.cpp
		  module3d/m3d_ext4.cpp
 )
 
else(EXTERNALIZE)

  set(state_SRC parser/parser.cpp)
  set(module3d_SRC module3d/module3d.cpp)
  
endif(EXTERNALIZE)

set(test_SRC  define.cpp
	      #kernel.cpp
	      #clustergraph.cpp
	      #sys_and_obj.cpp
	      #geometry.cpp
	      #clustermath.cpp
	      #constraints3d.cpp
	      ${module3d_SRC}
	      #dof.cpp
	      modulepart.cpp
	      modulepart2.cpp
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
