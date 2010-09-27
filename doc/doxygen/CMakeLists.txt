#===============================================================================
# @file   CMakeLists.txt
# @author Nicolas Richart <nicolas.richart@epfl.ch>
# @date   Mon Sep 27 17:35:59 2010
#
# @section LICENSE
#
# <insert lisence here>
#
# @section DESCRIPTION
#
#===============================================================================

set(DOXYGEN_INPUT ${CMAKE_CURRENT_BINARY_DIR}/akantu.dox)
set(DOXYGEN_OUTPUT_DIR ${CMAKE_CURRENT_BINARY_DIR}/html)
set(DOXYGEN_OUTPUT ${DOXYGEN_OUTPUT_DIR}/index.html)

make_directory(${DOXYGEN_OUTPUT_DIR})
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/akantu.dox.cmake
  ${CMAKE_CURRENT_BINARY_DIR}/akantu.dox)

add_custom_command(
  OUTPUT ${DOXYGEN_OUTPUT}
  COMMAND ${CMAKE_COMMAND} -E echo "${DOXYGEN_EXECUTABLE} ${DOXYGEN_INPUT}"
  COMMAND ${CMAKE_COMMAND} -E echo_append "Building akantu Documentation..."
  COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_INPUT}
  COMMAND ${CMAKE_COMMAND} -E echo "Done."
  DEPENDS ${DOXYGEN_INPUT}
  )

add_custom_target(akantu-doc ALL DEPENDS ${DOXYGEN_OUTPUT})

add_custom_target(akantu-doc-forced
  COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_INPUT}
  )

install(DIRECTORY ${DOXYGEN_OUTPUT_DIR}
  DESTINATION share/akantu-${akantu_VERSION}/doc)

