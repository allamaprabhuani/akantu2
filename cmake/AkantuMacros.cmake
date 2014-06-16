#===============================================================================
# @file   AkantuMacros.cmake
#
# @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
# @author Nicolas Richart <nicolas.richart@epfl.ch>
#
# @date   Thu Feb 17 17:28:21 2011
#
# @brief  Set of macros used by akantu cmake files
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
#===============================================================================

macro(check_for_isnan result)
  include(CheckFunctionExists)
  check_function_exists(std::isnan HAVE_STD_ISNAN)
  if(HAVE_STD_ISNAN)
    set(result "std::isnan(x)")
  else()
    check_function_exists(isnan HAVE_ISNAN)
    if(HAVE_ISNAN)
      set(result "(::isnan(x))")
    else()
      check_function_exists(_isnan HAVE_ISNAN_MATH_H)
      if(HAVE_ISNAN_MATH_H)
        set(result "(_isnan(x))")
      else()
        set(result (x == std::numeric_limits<Real>::quiet_NAN()))
      endif()
    endif()
  endif()
endmacro()

#===============================================================================
macro(copy_files target_depend)
  foreach(_file ${ARGN})
    set(_target ${CMAKE_CURRENT_BINARY_DIR}/${_file})
    set(_source ${CMAKE_CURRENT_SOURCE_DIR}/${_file})
    add_custom_command(
      OUTPUT ${_target}
      COMMAND ${CMAKE_COMMAND} -E copy_if_different ${_source} ${_target}
      DEPENDS ${_source}
      )
    set(_target_name ${target_depend}_${_file})
    add_custom_target(${_target_name} DEPENDS ${_target})
    add_dependencies(${target_depend} ${_target_name})
  endforeach()
endmacro()

#===============================================================================
macro(add_akantu_definitions)
  foreach(_definition ${AKANTU_DEFINITIONS})
    add_definitions(-D${_definition})
  endforeach()
endmacro()


macro(include_boost)
  find_package(Boost REQUIRED)
  mark_as_advanced(Boost_DIR)
  if(Boost_FOUND)
    list(APPEND AKANTU_EXTERNAL_LIB_INCLUDE_DIR ${Boost_INCLUDE_DIRS} ${Boost_INCLUDE_DIR})
  endif()
  message(STATUS "Looking for Boost liraries")
  set(AKANTU_BOOST_COMPONENT_ALREADY_DONE)
  foreach(_comp ${AKANTU_BOOST_COMPONENTS})
    list(FIND AKANTU_BOOST_COMPONENT_ALREADY_DONE ${_comp} _res)
    if(_res EQUAL -1)
      find_package(Boost COMPONENTS ${_comp} QUIET)
      string(TOUPPER ${_comp} _u_comp)
      if(Boost_${_u_comp}_FOUND)
	message(STATUS "   ${_comp}: FOUND")
	set(AKANTU_BOOST_${_u_comp} TRUE CACHE INTERNAL "" FORCE)
	list(APPEND AKANTU_EXTERNAL_LIBRARIES ${Boost_${_u_comp}_LIBRARY})
      else()
	message(STATUS "   ${_comp}: NOT FOUND")
      endif()
      list(APPEND AKANTU_BOOST_COMPONENT_ALREADY_DONE ${_comp})
    endif()
  endforeach()
endmacro()

