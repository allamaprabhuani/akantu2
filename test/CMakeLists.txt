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

if(AKANTU_SCOTCH_ON)
  add_example(test_mesh_partitionate "Test mesh partition creation")
endif()

if(AKANTU_MPI_ON)
  add_example(test_synchronizer "Test synchronizers")
endif()
