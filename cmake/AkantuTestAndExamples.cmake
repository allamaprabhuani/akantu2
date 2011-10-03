#===============================================================================
# @file   AkantuTestAndExamples.cmake
# @author Nicolas Richart <nicolas.richart@epfl.ch>
# @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
# @date   Mon Oct 25 09:46:59 2010
#
# @section LICENSE
#
# Copyright (©) 2010-2011 EPFL (Ecole Polytechnique Fédérale de Lausanne)
# Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
#
# Akantu is free  software: you can redistribute it and/or  modify it under the
# terms  of the  GNU Lesser  General Public  License as  published by  the Free
# Software Foundation, either version 3 of the License, or (at your option) any
# later version.
#
# Akantu is  distributed in the  hope that it  will be useful, but  WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
# A  PARTICULAR PURPOSE. See  the GNU  Lesser General  Public License  for more
# details.
#
# You should  have received  a copy  of the GNU  Lesser General  Public License
# along with Akantu. If not, see <http://www.gnu.org/licenses/>.
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
#  mark_as_advanced(AKANTU_BUILD_ALL_TESTS)
endif(AKANTU_TESTS)

#===============================================================================
macro(register_test test_name)
  add_executable(${test_name} ${ARGN})
  target_link_libraries(${test_name} akantu ${AKANTU_EXTERNAL_LIBRARIES})

  if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/${test_name}.sh)
    file(COPY ${test_name}.sh DESTINATION .)
  endif()

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
#  mark_as_advanced(AKANTU_BUILD_ALL_EXAMPLES)
endif(AKANTU_EXAMPLES)

#===============================================================================
macro(register_example example_name)
  add_executable(${example_name} ${ARGN})
  target_link_libraries(${example_name} akantu ${AKANTU_EXTERNAL_LIBRARIES})
endmacro()

#===============================================================================
macro(add_example example_name desc)
  manage_test_and_example(${example_name} ${desc} AKANTU_BUILD_ALL_EXAMPLES _EXAMPLE_)
endmacro()
