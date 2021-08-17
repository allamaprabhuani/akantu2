#===============================================================================
# @file   FindParMETIS.cmake
#
# @author Nicolas Richart <nicolas.richart@epfl.ch>
#
# @date creation: Fri Apr 22 2016
# @date last modification: Fri Mar 16 2018
#
# @brief  FindPackage for the parmetis library
#
#
# @section LICENSE
#
# Copyright (©) 2016-2021 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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


find_path(PARMETIS_INCLUDE_DIR parmetis.h
  PATHS "${PARMETIS_DIR}"
  ENV PARMETIS_DIR
  PATH_SUFFIXES include
  )

find_library(PARMETIS_LIBRARY NAMES parmetis
  PATHS "${PARMETIS_DIR}"
  ENV PARMETIS_DIR
  PATH_SUFFIXES lib
  )

mark_as_advanced(PARMETIS_LIBRARY PARMETIS_INCLUDE_DIR)

#===============================================================================
include(FindPackageHandleStandardArgs)
if(CMAKE_VERSION VERSION_GREATER 2.8.12)
  if(PARMETIS_INCLUDE_DIR)
    file(STRINGS ${PARMETIS_INCLUDE_DIR}/parmetis.h _versions
      REGEX "^#define\ +PARMETIS_(MAJOR|MINOR|SUBMINOR)_VERSION .*")
    foreach(_ver ${_versions})
      string(REGEX MATCH "PARMETIS_(MAJOR|MINOR|SUBMINOR)_VERSION *([0-9.]+)" _tmp "${_ver}")
      set(_parmetis_${CMAKE_MATCH_1} ${CMAKE_MATCH_2})
    endforeach()
    set(PARMETIS_VERSION "${_parmetis_MAJOR}.${_parmetis_MINOR}" CACHE INTERNAL "")
  endif()

  find_package_handle_standard_args(ParMETIS
    REQUIRED_VARS
      PARMETIS_LIBRARY
      PARMETIS_INCLUDE_DIR
    VERSION_VAR
      PARMETIS_VERSION)
else()
  find_package_handle_standard_args(ParMETIS DEFAULT_MSG
    PARMETIS_LIBRARY PARMETIS_INCLUDE_DIR)
endif()
