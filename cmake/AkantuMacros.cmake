#===============================================================================
# @file   AkantuMacros.cmake
# @author Nicolas Richart <nicolas.richart@epfl.ch>
# @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
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
set(AKANTU_CMAKE_DEBUG TRUE)
macro(akantu_message)
  if(AKANTU_CMAKE_DEBUG)
    message(${ARGN})
  endif()
endmacro()

#===============================================================================
macro(add_cxx_flags flag)
  if(NOT CMAKE_CXX_FLAGS MATCHES "${flag}")
    set(CMAKE_CXX_FLAGS "${flag} ${CMAKE_CXX_FLAGS}"
      CACHE STRING "Flags used by the compiler during all build types." FORCE)
  endif()
endmacro()

#===============================================================================
macro(remove_cxx_flags flag)
  string(REPLACE "${flag} " "" CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}"
    CACHE STRING "Flags used by the compiler during all build types." FORCE)
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

#===============================================================================
macro(add_akantu_definitions)
  foreach(_definition ${AKANTU_DEFINITIONS})
    add_definitions(-D${_definition})
  endforeach()
endmacro()

#===============================================================================
if(__CMAKE_PARSE_ARGUMENTS_INCLUDED)
  return()
endif()
set(__CMAKE_PARSE_ARGUMENTS_INCLUDED TRUE)

function(CMAKE_PARSE_ARGUMENTS prefix _optionNames _singleArgNames _multiArgNames)
  # first set all result variables to empty/FALSE
  foreach(arg_name ${_singleArgNames} ${_multiArgNames})
    set(${prefix}_${arg_name})
  endforeach(arg_name)

  foreach(option ${_optionNames})
    set(${prefix}_${option} FALSE)
  endforeach(option)

  set(${prefix}_UNPARSED_ARGUMENTS)

  set(insideValues FALSE)
  set(currentArgName)

  # now iterate over all arguments and fill the result variables
  foreach(currentArg ${ARGN})
    list(FIND _optionNames "${currentArg}" optionIndex)  # ... then this marks the end of the arguments belonging to this keyword
    list(FIND _singleArgNames "${currentArg}" singleArgIndex)  # ... then this marks the end of the arguments belonging to this keyword
    list(FIND _multiArgNames "${currentArg}" multiArgIndex)  # ... then this marks the end of the arguments belonging to this keyword

    if(${optionIndex} EQUAL -1  AND  ${singleArgIndex} EQUAL -1  AND  ${multiArgIndex} EQUAL -1)
      if(insideValues)
        if("${insideValues}" STREQUAL "SINGLE")
          set(${prefix}_${currentArgName} ${currentArg})
          set(insideValues FALSE)
        elseif("${insideValues}" STREQUAL "MULTI")
          list(APPEND ${prefix}_${currentArgName} ${currentArg})
        endif()
      else(insideValues)
        list(APPEND ${prefix}_UNPARSED_ARGUMENTS ${currentArg})
      endif(insideValues)
    else()
      if(NOT ${optionIndex} EQUAL -1)
        set(${prefix}_${currentArg} TRUE)
        set(insideValues FALSE)
      elseif(NOT ${singleArgIndex} EQUAL -1)
        set(currentArgName ${currentArg})
        set(${prefix}_${currentArgName})
        set(insideValues "SINGLE")
      elseif(NOT ${multiArgIndex} EQUAL -1)
        set(currentArgName ${currentArg})
        set(${prefix}_${currentArgName})
        set(insideValues "MULTI")
      endif()
    endif()

  endforeach(currentArg)

  # propagate the result variables to the caller:
  foreach(arg_name ${_singleArgNames} ${_multiArgNames} ${_optionNames})
    set(${prefix}_${arg_name}  ${${prefix}_${arg_name}} PARENT_SCOPE)
  endforeach(arg_name)
  set(${prefix}_UNPARSED_ARGUMENTS ${${prefix}_UNPARSED_ARGUMENTS} PARENT_SCOPE)
endfunction(CMAKE_PARSE_ARGUMENTS _options _singleArgs _multiArgs)
