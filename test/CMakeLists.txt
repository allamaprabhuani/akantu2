#===============================================================================
# @file   CMakeLists.txt
# @author Nicolas Richart <nicolas.richart@epfl.ch>
# @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
# @date   Fri Jun 11 09:46:59 2010
#
# @section LICENSE
#
# <insert lisence here>
#
# @section DESCRIPTION
#
#===============================================================================

#===============================================================================
# Tests
#===============================================================================

option(AKANTU_BUILD_ALL_TESTS "Build all tests")

#===============================================================================
macro(register_test test_name source_list)
  add_executable(${test_name} ${source_list})
  target_link_libraries(${test_name} akantu ${AKANTU_EXTERNAL_LIBRARIES})
endmacro()

#===============================================================================
macro(add_mesh MESH_TARGET GEO_FILE DIM ORDER)
  set(_geo_file ${CMAKE_CURRENT_SOURCE_DIR}/${GEO_FILE})
  get_filename_component(_msh_file "${GEO_FILE}" NAME_WE)
  set(_msh_file ${CMAKE_CURRENT_BINARY_DIR}/${_msh_file}.msh)
  if(EXISTS ${_geo_file})
    add_custom_command(
      OUTPUT ${_msh_file}
      DEPENDS ${_geo_file}
      COMMAND gmsh
      ARGS -${DIM} -order ${ORDER} -optimize -o ${_msh_file} ${_geo_file} 2>&1 > /dev/null
      )
    add_custom_target(${MESH_TARGET}
      DEPENDS ${_msh_file})
  else(EXISTS ${_geo_file})
    message(FATAL_ERROR "File ${_geo_file} not found")
  endif(EXISTS ${_geo_file})
endmacro()

#===============================================================================
macro(add_example example_name desc)
  string(TOUPPER ${example_name} upper_name)

  option(AKANTU_BUILD_${upper_name} "${desc}")
  mark_as_advanced(AKANTU_BUILD_${upper_name})

  if(AKANTU_BUILD_ALL_TESTS)
    set(AKANTU_BUILD_${upper_name}_OLD
      ${AKANTU_BUILD_${upper_name}}
      CACHE INTERNAL "${desc}" FORCE)

    set(AKANTU_BUILD_${upper_name} ON
      CACHE INTERNAL "${desc}" FORCE)
  else(AKANTU_BUILD_ALL_TESTS)
    if(DEFINED AKANTU_BUILD_${upper_name}_OLD)
      set(AKANTU_BUILD_${upper_name}
	${AKANTU_BUILD_${upper_name}_OLD}
	CACHE BOOL "${desc}" FORCE)

      unset(BUILD_${upper_name}_OLD
	CACHE)
    endif(DEFINED AKANTU_BUILD_${upper_name}_OLD)
  endif(AKANTU_BUILD_ALL_TESTS)

  if(AKANTU_BUILD_${upper_name})
    add_subdirectory(${example_name})
  endif(AKANTU_BUILD_${upper_name})

endmacro()

#===============================================================================
# List of tests
#===============================================================================
add_example(test_facet_extraction "Test mesh utils facet extraction")
add_example(test_fem "Test finite element functionalties")
add_example(test_mesh_io "Test mesh io object")
add_example(test_solid_mechanics_model "Test solid mechanics object")
add_example(test_static_memory "Test static memory")
add_example(test_vector "Test akantu vector")

if(AKANTU_USE_SCOTCH AND SCOTCH_FOUND)
  add_example(test_mesh_partitionate "Test mesh partition creation")
endif(AKANTU_USE_SCOTCH AND SCOTCH_FOUND)
