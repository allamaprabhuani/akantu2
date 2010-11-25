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
if(AKANTU_TESTS)
  option(AKANTU_BUILD_ALL_TESTS "Build all tests")
  mark_as_advanced(AKANTU_BUILD_ALL_TESTS)
endif(AKANTU_TESTS)
#===============================================================================
macro(register_test test_name ${ARGN})
  add_executable(${test_name} ${ARGN})
  target_link_libraries(${test_name} akantu ${AKANTU_EXTERNAL_LIBRARIES})
endmacro()

#===============================================================================
macro(add_akantu_test test_name desc)
  string(TOUPPER ${test_name} upper_name)

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
    add_subdirectory(${test_name})
  endif(AKANTU_BUILD_${upper_name})
endmacro()

#===============================================================================
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
  string(TOUPPER ${example_name} upper_name)

  option(AKANTU_BUILD_${upper_name} "${desc}")
  mark_as_advanced(AKANTU_BUILD_${upper_name})

  if(AKANTU_BUILD_ALL_EXAMPLES)
    set(AKANTU_BUIcLD_${upper_name}_OLD
      ${AKANTU_BUILD_${upper_name}}
      CACHE INTERNAL "${desc}" FORCE)

    set(AKANTU_BUILD_${upper_name} ON
      CACHE INTERNAL "${desc}" FORCE)
  else(AKANTU_BUILD_ALL_EXAMPLES)
    if(DEFINED AKANTU_BUILD_${upper_name}_OLD)
      set(AKANTU_BUILD_${upper_name}
	${AKANTU_BUILD_${upper_name}_OLD}
	CACHE BOOL "${desc}" FORCE)

      unset(BUILD_${upper_name}_OLD
	CACHE)
    endif(DEFINED AKANTU_BUILD_${upper_name}_OLD)
  endif(AKANTU_BUILD_ALL_EXAMPLES)

  if(AKANTU_BUILD_${upper_name})
    add_subdirectory(${example_name})
  endif(AKANTU_BUILD_${upper_name})

endmacro()
