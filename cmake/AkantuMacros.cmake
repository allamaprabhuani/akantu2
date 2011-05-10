#===============================================================================
# @file   AkantuMacros.cmake
# @author Nicolas Richart <nicolas.richart@epfl.ch>
# @date   Wed Feb  9 10:59:42 2011
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

#===============================================================================
# macro to include optional packages
macro(add_optional_package PACKAGE DESC DEFAULT)
  string(TOUPPER ${PACKAGE} _u_package)
  set(option_name AKANTU_USE_${_u_package})

  option(${option_name} ${DESC} ${DEFAULT})

  if(${option_name})
    if(${PACKAGE} MATCHES BLAS)
      enable_language(Fortran)
    endif()
    find_package(${PACKAGE} REQUIRED)
    if(${_u_package}_FOUND)
      add_definitions(-DAKANTU_USE_${_u_package})
      set(AKANTU_EXTERNAL_LIB_INCLUDE_PATH ${AKANTU_EXTERNAL_LIB_INCLUDE_PATH}
	${${_u_package}_INCLUDE_PATH}
	)
      set(AKANTU_EXTERNAL_LIBRARIES ${AKANTU_EXTERNAL_LIBRARIES}
	${${_u_package}_LIBRARIES}
	)
      set(AKANTU_${_u_package}_ON ON)
    endif(${_u_package}_FOUND)
  else()
    set(AKANTU_${_u_package}_ON OFF)
  endif(${option_name})
endmacro()

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