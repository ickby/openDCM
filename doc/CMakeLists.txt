
# Generating documentation with doxygen
#
if(GENERATE_DOCS AND DOXYGEN_FOUND)

  file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/html)
 
  #set the paths in the config file
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/doxygen ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY)
  
  # The initial rm command gets rid of everything previously built by this
  # custom command.
  add_custom_command(
  OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/html/index.html
  COMMAND ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
  DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  )

  add_custom_target(
  documentation ALL
  DEPENDS
  ${CMAKE_CURRENT_BINARY_DIR}/html/index.html
  )
  
  # Install the documentation generated at "make" time.
  #install(DIRECTORY ${CMAKE_SOURCE_DIR}/doc/ DESTINATION ${docdir}/html)

endif(GENERATE_DOCS AND DOXYGEN_FOUND)
