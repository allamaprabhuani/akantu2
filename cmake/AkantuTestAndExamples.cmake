#===============================================================================
# @file   AkantuTestAndExamples.cmake
# @author Nicolas Richart <nicolas.richart@epfl.ch>
# @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
# @date   Mon Oct 25 09:46:59 2010
#
# @section LICENSE
#
# <insert lisence here>
#
# @section DESCRIPTION
#
#===============================================================================

#===============================================================================
macro(manage_test_and_example et_name desc build_all label)
  string(TOUPPER ${et_name} upper_name)

  option(AKANTU_BUILD${label}${upper_name} "${desc}")
  mark_as_advanced(AKANTU_BUILD_${upper_name})

  if(${build_all})
    set(AKANTU_BUILD${label}${upper_name}_OLD
      ${AKANTU_BUILD${label}${upper_name}}
      CACHE INTERNAL "${desc}" FORCE)

    set(AKANTU_BUILD${label}${upper_name} ON
      CACHE INTERNAL "${desc}" FORCE)
  else(${build_all})
    if(DEFINED AKANTU_BUILD${label}${upper_name}_OLD)
      set(AKANTU_BUILD${label}${upper_name}
	${AKANTU_BUILD${label}${upper_name}_OLD}
	CACHE BOOL "${desc}" FORCE)

      unset(AKANTU_BUILD${label}${upper_name}_OLD
	CACHE)
    endif(DEFINED AKANTU_BUILD${label}${upper_name}_OLD)
  endif(${build_all})

  if(AKANTU_BUILD${label}${upper_name})
    add_subdirectory(${et_name})
  endif(AKANTU_BUILD${label}${upper_name})
endmacro()

#===============================================================================
# Tests
#===============================================================================
if(AKANTU_TESTS)
  option(AKANTU_BUILD_ALL_TESTS "Build all tests")
  mark_as_advanced(AKANTU_BUILD_ALL_TESTS)
endif(AKANTU_TESTS)

#===============================================================================
macro(register_test test_name ${ARGN})
  add_executable(${test_name} ${ARGN})
  target_link_libraries(${test_name} akantu ${AKANTU_EXTERNAL_LIBRARIES})

  if(EXISTS ${CMAKE_CURRENT_BINARY_DIR}/${test_name}.sh)
    add_test(${test_name} ${CMAKE_CURRENT_BINARY_DIR}/${test_name}.sh)
  else(EXISTS ${CMAKE_CURRENT_BINARY_DIR}/${test_name}.sh)
    add_test(${test_name} ${CMAKE_CURRENT_BINARY_DIR}/${test_name})
  endif(EXISTS ${CMAKE_CURRENT_BINARY_DIR}/${test_name}.sh)
endmacro()

#===============================================================================
macro(add_akantu_test test_name desc)
  manage_test_and_example(${test_name} ${desc} AKANTU_BUILD_ALL_TESTS _)
endmacro()


#===============================================================================
# Examples
#===============================================================================
if(AKANTU_EXAMPLES)
  option(AKANTU_BUILD_ALL_EXAMPLES "Build all examples")
  mark_as_advanced(AKANTU_BUILD_ALL_EXAMPLES)
endif(AKANTU_EXAMPLES)

#===============================================================================
macro(register_example example_name source_list)
  add_executable(${example_name} ${source_list})
  target_link_libraries(${example_name} akantu ${AKANTU_EXTERNAL_LIBRARIES})
endmacro()

#===============================================================================
macro(add_example example_name desc)
  manage_test_and_example(${example_name} ${desc} AKANTU_BUILD_ALL_EXAMPLES _EXAMPLE_)
endmacro()
