
# Generating documentation with doxygen
#
if(GENERATE_DOCS AND DOXYGEN_FOUND)

  file(MAKE_DIRECTORY ${CMAKE_SOURCE_DIR}/doc/html)
 
  # The initial rm command gets rid of everything previously built by this
  # custom command.
  add_custom_command(
  OUTPUT ${CMAKE_SOURCE_DIR}/doc/html/index.html
  COMMAND rm -rf ${CMAKE_SOURCE_DIR}/doc/html
  COMMAND mkdir ${CMAKE_SOURCE_DIR}/doc/html
  COMMAND ${DOXYGEN_EXECUTABLE} ${CMAKE_SOURCE_DIR}/doc/doxygen
  DEPENDS ${CMAKE_SOURCE_DIR}/doc/doxygen
  WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
  )

  add_custom_target(
  documentation ALL
  DEPENDS
  ${CMAKE_SOURCE_DIR}/doc/html/index.html
  )
  
  # Install the documentation generated at "make" time.
  #install(DIRECTORY ${CMAKE_SOURCE_DIR}/doc/ DESTINATION ${docdir}/html)

endif(GENERATE_DOCS AND DOXYGEN_FOUND)
