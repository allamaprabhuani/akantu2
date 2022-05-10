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

function(_match_semver _input_semver prefix)
  set(_semver_regexp
    "^([0-9]+(\\.[0-9]+)?(\\.[0-9]+)?)(-([a-zA-Z0-9\-]*))?(\\+(.*))?")

  if(_input_semver MATCHES "^([0-9]+(\\.[0-9]+)?(\\.[0-9]+)?)(-([a-zA-Z0-9-]*))?(\\+(.*))?")
    set(${prefix}_version ${CMAKE_MATCH_1} PARENT_SCOPE)
    if(CMAKE_MATCH_4)
      set(${prefix}_version_prerelease "${CMAKE_MATCH_5}" PARENT_SCOPE)
    endif()

    if(CMAKE_MATCH_6)
      set(${prefix}_version_metadata "${CMAKE_MATCH_7}" PARENT_SCOPE)
    endif()

  endif()
endfunction()

function(_get_version_from_git)
  if(NOT CMAKE_VERSION_GENERATOR_TAG_PREFIX)
    set(CMAKE_VERSION_GENERATOR_TAG_PREFIX "v")
  endif()

  find_package(Git)

  if(NOT Git_FOUND)
    return()
  endif()

  execute_process(
    COMMAND ${GIT_EXECUTABLE} describe
    --tags
    --abbrev=0
    --match ${CMAKE_VERSION_GENERATOR_TAG_PREFIX}*
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
    RESULT_VARIABLE _res
    OUTPUT_VARIABLE _out_tag
    OUTPUT_STRIP_TRAILING_WHITESPACE
    ERROR_VARIABLE _err_tag)

  if(NOT _res EQUAL 0)
    return()
  endif()

  string(REGEX REPLACE "^${CMAKE_VERSION_GENERATOR_TAG_PREFIX}(.*)" "\\1" _tag "${_out_tag}")

  _match_semver("${_tag}" _tag)

  execute_process(
    COMMAND ${GIT_EXECUTABLE} describe
    --tags
    --dirty
    --always
    --long
    --match ${CMAKE_VERSION_GENERATOR_TAG_PREFIX}*
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
    RESULT_VARIABLE _res
    OUTPUT_VARIABLE _out
    OUTPUT_STRIP_TRAILING_WHITESPACE)

  set(_git_version ${_tag_version} PARENT_SCOPE)

  if(_tag_version_prerelease)
    set(_git_version_prerelease ${_tag_version_prerelease} PARENT_SCOPE)
  endif()

  # git describe to PEP404 version
  set(_version_regex
    "^${CMAKE_VERSION_GENERATOR_TAG_PREFIX}${_tag}(-([0-9]+)-(g[0-9a-f]+)(-dirty)?)?$")

  if(_out MATCHES ${_version_regex})
    if(CMAKE_MATCH_1)
      if(_tag_version_metadata)
        set(_metadata "${_tag_version_metadata}.")
      endif()
      set(_metadata "${_metadata}${CMAKE_MATCH_2}.${CMAKE_MATCH_3}")
    endif()
    if(CMAKE_MATCH_4)
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
        set(_metadata "${_out_count}.${_metadata}")
      endif()

      if(CMAKE_MATCH_2)
        set(_metadata "${_metadata}.dirty")
      endif()
    endif()
  endif()
  set(_git_version_metadata ${_metadata} PARENT_SCOPE)

endfunction()

function(_get_version_from_file)
  if(EXISTS ${PROJECT_SOURCE_DIR}/VERSION)
    file(STRINGS ${PROJECT_SOURCE_DIR}/VERSION _file_version)

    _match_semver("${_file_version}" "_file")

    set(_file_version ${_file_version} PARENT_SCOPE)
    if(_file_version_metadata)
      set(_file_version_metadata ${_file_version_metadata} PARENT_SCOPE)
    endif()

    if(_file_version_prerelease)
      set(_file_version_prerelease ${_file_version_prerelease} PARENT_SCOPE)
    endif()
  endif()
endfunction()

function(_get_metadata_from_ci)
  if(NOT DEFINED ENV{CI})
    return()
  endif()

  if(DEFINED ENV{CI_MERGE_REQUEST_ID})
    set(_ci_version_metadata "ci.mr$ENV{CI_MERGE_REQUEST_ID}" PARENT_SCOPE)
  endif()
endfunction()

function(define_project_version)
  string(TOUPPER ${PROJECT_NAME} _project)

  _get_version_from_git()

  if(_git_version)
    set(_version "${_git_version}")
    if(_git_version_metadata)
      set(_version_metadata "${_git_version_metadata}")
    endif()

    if (_git_version_prerelease)
      set(_version_prerelease "${_git_version_prerelease}")
    endif()
  else()
    # we can get metadata if and no version if not tag is properly defined
    if(_git_version_metadata)
      set(git_version_metadata ".${_git_version_metadata}")
    endif()

    _get_version_from_file()

    if(_file_version_metadata)
      set(_version_metadata "${_file_version_metadata}${_git_version_metadata}")
    endif()

    if (_file_version)
      set(_version "${_file_version}")
    endif()

    if (_file_version_prerelease)
      set(_version_prerelease "${_file_version_prerelease}")
    endif()
  endif()

  _get_metadata_from_ci()

  if(_version)
    if(_version_prerelease)
      set(_version_prerelease "-${_version_prerelease}")
    endif()
    if(_version_metadata)
      set(_version_metadata "+${_version_metadata}")
      if(_ci_version_metadata)
        set(_version_metadata "${_version_metadata}.${_ci_version_metadata}")
      endif()
    endif()
    set(${_project}_VERSION ${_version} PARENT_SCOPE)

    set(_semver "${_version}${_version_prerelease}${_version_metadata}")
    set(${_project}_SEMVER "${_semver}" PARENT_SCOPE)
    message(STATUS "${PROJECT_NAME} version: ${_semver}")

    if(_version MATCHES "^([0-9]+)(\\.([0-9]+))?(\\.([0-9]+))?")
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
      if(_version_prerelease)
        set(${_project}_PRERELEASE_VERSION ${_version_prerelease} PARENT_SCOPE)
      endif()
      if(_version_metadata)
        set(${_project}_LOCAL_VERSION ${_version_metadata} PARENT_SCOPE)
      endif()
      if(_version_metadata MATCHES "(\\+([0-9]+))(\\..*)*")
        set(${_project}_TWEAK ${CMAKE_MATCH_1} PARENT_SCOPE)
      endif()
    endif()
  else()
    message(FATAL_ERROR "Could not determine the VERSION for ${PROJECT_NAME}")
  endif()

  if(NOT ${_project}_NO_LIBRARY_VERSION)
    set(${_project}_LIBRARY_PROPERTIES ${${_project}_LIBRARY_PROPERTIES}
      VERSION "${_version}"
      SOVERSION "${_major_version}.${_minor_version}"
      )
  endif()
endfunction()
