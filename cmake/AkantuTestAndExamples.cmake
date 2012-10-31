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

set(AKANTU_DIFF_SCRIPT ${AKANTU_CMAKE_DIR}/akantu_diff.sh)

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
  elseif(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/${test_name}.verified)
    add_test(${test_name} ${AKANTU_DIFF_SCRIPT} ${test_name} ${CMAKE_CURRENT_SOURCE_DIR}/${test_name}.verified)
  else()
    add_test(${test_name} ${CMAKE_CURRENT_BINARY_DIR}/${test_name})
  endif()
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

#===============================================================================
macro(register_test_new test_name)
  set(multi_variables
    SOURCES FILES_TO_COPY DEPENDENCIES DIRECTORIES_TO_CREATE COMPILE_OPTIONS
    )

  cmake_parse_arguments(register_test
    ""
    "PACKAGE"
    "${multi_variables}"
    ${ARGN}
    )

  message("register_test_SOURCES :               ${register_test_SOURCES}")
  message("register_test_FILES_TO_COPY :         ${register_test_FILES_TO_COPY}")
  message("register_test_DIRECTORIES_TO_CREATE : ${register_test_DIRECTORIES_TO_CREATE}")
  message("register_test_DEPENDENCIES :          ${register_test_DEPENDENCIES}")
  message("register_test_COMPILE_OPTIONS :       ${register_test_COMPILE_OPTIONS}")
  message("register_test_PACKAGE :               ${register_test_PACKAGE}")

  if(register_test_PACKAGE)
    package_pkg_name(${register_test_PACKAGE} _package_name)
    list(FIND ${_package_name}_TESTS ${test_name} _ret)
    if(_ret EQUAL -1)
      list(APPEND ${_package_name}_TESTS ${test_name})
    endif()
  endif()

  # check if the test should be activated
  set(_activate_test 0)
  foreach(_pkg ${PACKAGE_SYSTEM_PACKAGES_ON})
    package_pkg_name(${_pkg} _package_name)
    list(FIND ${_package_name}_TESTS ${test_name} _ret)
    if(NOT _ret EQUAL -1)
      set(_activate_test 1)
    endif()
  endforeach()

  # check if the package is registered in at least a package
  set(_present_in_packages 0)
  foreach(_pkg ${PACKAGE_SYSTEM_PACKAGES_NAMES_LIST_ALL})
    package_pkg_name(${_pkg} _package_name)
    list(FIND ${_package_name}_TESTS ${test_name} _ret)
    if(NOT _ret EQUAL -1)
      set(_present_in_packages 1)
    endif()
  endforeach()

  if(NOT _present_in_packages)
    message("The test ${test_name} is not registered in any packages")
  endif()

  message("TEST  ${test_name} active: ${_activate_test}")
  if(_activate_test)
    add_executable(${test_name} ${register_test_SOURCES})
    if(register_test_COMPILE_OPTIONS)
      set_target_properties(${test_name}
	PROPERTIES COMPILE_DEFINITIONS ${register_test_COMPILE_OPTIONS})
    endif()
    target_link_libraries(${test_name} akantu ${AKANTU_EXTERNAL_LIBRARIES})

    if(_test_option_FILES_TO_COPY)
      foreach(_file ${register_test_FILES_TO_COPY})
	file(COPY ${_file} DESTINATION .)
      endforeach()
    endif()

    if(register_test_DIRECTORIES_TO_CREATE)
      foreach(_dir ${register_test_DIRECTORIES_TO_CREATE})
	if(IS_ABSOLUTE ${dir})
	  file(MAKE_DIRECTORY ${_dir})
	else()
  	  file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${_dir})
	endif()
      endforeach()
    endif()

    if(register_test_DEPENDENCIES)
      foreach(dep ${register_test_DEPENDENCIES})
	add_dependencies(${test_name} ${dep})
	get_target_property(dep_in_ressources ${dep} RESSOURCES)
      endforeach()
    endif()

    if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/${test_name}.sh)
      file(COPY ${test_name}.sh DESTINATION .)
      add_test(${test_name} ${CMAKE_CURRENT_BINARY_DIR}/${test_name}.sh)
    elseif(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/${test_name}.verified)
      add_test(${test_name} ${AKANTU_DIFF_SCRIPT} ${test_name} ${CMAKE_CURRENT_SOURCE_DIR}/${test_name}.verified)
    else()
      add_test(${test_name} ${CMAKE_CURRENT_BINARY_DIR}/${test_name})
    endif()
  endif()
endmacro(register_test_new)