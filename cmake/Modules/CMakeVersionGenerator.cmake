#===============================================================================
# @file   CMakeVersionGenerator.cmake
#
# @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
# @author Nicolas Richart <nicolas.richart@epfl.ch>
#
# @date creation: Sun Oct 19 2014
# @date last modification: Mon Jan 18 2016
#
# @brief  Set of macros used by akantu to handle the package system
#
#
# @section LICENSE
#
# Copyright (©) 2015-2021 EPFL (Ecole Polytechnique Fédérale de Lausanne)
# Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
#
# Akantu is free software: you can redistribute it and/or modify it under the
# terms of the GNU Lesser General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option) any
# later version.
# 
# Akantu is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
# A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more
# details.
# 
# You should have received a copy of the GNU Lesser General Public License along
# with Akantu. If not, see <http://www.gnu.org/licenses/>.
#
#===============================================================================


if(__DEFINE_PROJECT_VERSION__)
  return()
endif()
set(__DEFINE_PROJECT_VERSION__ TRUE)

function(_get_version_from_git)
  find_package(Git)
  if(Git_FOUND)
    execute_process(
      COMMAND ${GIT_EXECUTABLE} describe
        --tags
        --dirty
        --always
        --long
        --match v*
      WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
      RESULT_VARIABLE _res
      OUTPUT_VARIABLE _out
      OUTPUT_STRIP_TRAILING_WHITESPACE)

    # git describe to PEP404 version
    set(_version_regex
      "^v([0-9.]+)(-([0-9]+)-g([0-9a-f]+)(-dirty)?)?$")

    if(_out MATCHES ${_version_regex})
      set(_version ${CMAKE_MATCH_1} PARENT_SCOPE)
      if(CMAKE_MATCH_2)
        set(_metadata "${CMAKE_MATCH_3}.${CMAKE_MATCH_4}")
      endif()
      if(CMAKE_MATCH_5)
        set(_metadata "${_metadata}.dirty")
      endif()
    else()
      execute_process(
        COMMAND ${GIT_EXECUTABLE} rev-list HEAD --count
        WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
        RESULT_VARIABLE _res
        OUTPUT_VARIABLE _out_count
        OUTPUT_STRIP_TRAILING_WHITESPACE)

      if(_out MATCHES "^([0-9a-f]+)(-dirty)?$")
        set(_metadata "${CMAKE_MATCH_1}")
        if(_res EQUAL 0)
          set(_metadata "${_out_count}.${_version_metadata}")
        endif()

        if(CMAKE_MATCH_2)
          set(_metadata "${_version_metadata}.dirty")
        endif()
      endif()
    endif()
    set(_version_metadata ${_metadata} PARENT_SCOPE)
  endif()
endfunction()

function(_get_version_from_file)
  if(EXISTS ${PROJECT_SOURCE_DIR}/VERSION)
    file(STRINGS ${PROJECT_SOURCE_DIR}/VERSION _file_version)
    if("_file_version" MATCHES "^([0-9]+(\\.[0-9]+)?(\\.[0-9]+)?)(\\+(.*))?")
      set(_version ${CMAKE_MATCH_1} PARENT_SCOPE)
      if(CMAKE_MATCH_4)
        set(_version_metadata ${CMAKE_MATCH_4} PARENT_SCOPE)
      endif()
    endif()
  endif()
endfunction()

function(define_project_version)
  string(TOUPPER ${PROJECT_NAME} _project)

  _get_version_from_git()
  if(NOT _version)
    _get_version_from_file()
  endif()

  if(_version)
    set(${_project}_VERSION ${_version} PARENT_SCOPE)
    if(_version_metadata)
      set(${_project}_SEMVER "${_version}+${_version_metadata}" PARENT_SCOPE)
      message(STATUS "${PROJECT_NAME} version: ${_version}+${_version_metadata}")
    else()
      message(STATUS "${PROJECT_NAME} version: ${_version}")
    endif()

    if(_version MATCHES "([0-9]+)(\\.([0-9]+))?(\\.([0-9]+))?")
      set(_major_version ${CMAKE_MATCH_1})
      set(${_project}_MAJOR_VERSION ${_major_version} PARENT_SCOPE)
      if(CMAKE_MATCH_2)
        set(_minor_version ${CMAKE_MATCH_3})
        set(${_project}_MINOR_VERSION ${_minor_version} PARENT_SCOPE)
      endif()
      if(CMAKE_MATCH_4)
        set(_patch_version ${CMAKE_MATCH_5})
        set(${_project}_PATCH_VERSION ${_patch_version} PARENT_SCOPE)
      endif()
    endif()
  endif()

  if(NOT ${_project}_NO_LIBRARY_VERSION)
    set(${_project}_LIBRARY_PROPERTIES ${${_project}_LIBRARY_PROPERTIES}
      VERSION "${_version}"
      SOVERSION "${_major_version}.${_minor_version}"
      )
  endif()
endfunction()
