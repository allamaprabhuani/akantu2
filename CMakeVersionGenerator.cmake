#===============================================================================
# @file   CMakeVersionGenerator.cmake
#
# @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
# @author Nicolas Richart <nicolas.richart@epfl.ch>
#
# @date   Tue Oct 16 14:05:02 2012
#
# @brief  Set of macros used by akantu to handle the package system
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

if(__DEFINE_PROJECT_VERSION__)
  return()
endif()
set(__DEFINE_PROJECT_VERSION__ TRUE)

macro(define_project_version)
  find_package(Subversion)
  string(TOUPPER ${PROJECT_NAME} _project)

  if(SUBVERSION_FOUND)
    subversion_wc_info(${PROJECT_SOURCE_DIR} ${_project} ERROR_QUIET)
    if(${${_project}_WC_FOUND})
      set(${_project}_BUILD_VERSION ${${_project}_WC_REVISION})
      set(${_project}_VERSION
	"${${_project}_MAJOR_VERSION}.${${_project}_MINOR_VERSION}.${${_project}_BUILD_VERSION}"
	)
      file(WRITE VERSION "${${_project}_VERSION}\n")
    endif()
  endif()

  if(NOT ${_project}_VERSION)
    if(EXISTS ${PROJECT_SOURCE_DIR}/VERSION)
      file(STRINGS VERSION ${_project}_VERSION)
    else()
      set(${_project}_VERSION
        "${${_project}_MAJOR_VERSION}.${${_project}_MINOR_VERSION}"
        )
    endif()
  endif()

  # Append the library version information to the library target properties
  if(NOT ${_project}_NO_LIBRARY_VERSION)
    message(STATUS "${PROJECT_NAME} version: ${${_project}_VERSION}")

    set(${_project}_LIBRARY_PROPERTIES ${${_project}_LIBRARY_PROPERTIES}
      VERSION "${${_project}_VERSION}"
      SOVERSION "${${_project}_MAJOR_VERSION}.${${_project}_MINOR_VERSION}"
      )
  endif()
endmacro()