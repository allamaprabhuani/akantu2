#===============================================================================
# @file   FindMETIS.cmake
#
# @author Nicolas Richart <nicolas.richart@epfl.ch>
#
# @date creation: Fri Apr 22 2016
# @date last modification: Fri Mar 16 2018
#
# @brief  FindPackage for the metis library
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


find_path(METIS_INCLUDE_DIR metis.h
  PATHS "${METIS_DIR}"
  ENV METIS_DIR
  PATH_SUFFIXES include
  )

find_library(METIS_LIBRARY NAMES metis
  PATHS "${METIS_DIR}"
  ENV METIS_DIR
  PATH_SUFFIXES lib
  )

mark_as_advanced(METIS_LIBRARY METIS_INCLUDE_DIR)

#===============================================================================
include(FindPackageHandleStandardArgs)
if(CMAKE_VERSION VERSION_GREATER 2.8.12)
  if(METIS_INCLUDE_DIR)
    file(STRINGS ${METIS_INCLUDE_DIR}/metis.h _versions
      REGEX "^#define\ +METIS_VER_(MAJOR|MINOR|SUBMINOR) .*")
    foreach(_ver ${_versions})
      string(REGEX MATCH "METIS_VER_(MAJOR|MINOR|SUBMINOR) *([0-9.]+)" _tmp "${_ver}")
      set(_metis_${CMAKE_MATCH_1} ${CMAKE_MATCH_2})
    endforeach()
    set(METIS_VERSION "${_metis_MAJOR}.${_metis_MINOR}" CACHE INTERNAL "")
  endif()

  find_package_handle_standard_args(METIS
    REQUIRED_VARS
      METIS_LIBRARY
      METIS_INCLUDE_DIR
    VERSION_VAR
      METIS_VERSION)
else()
  find_package_handle_standard_args(METIS DEFAULT_MSG
    METIS_LIBRARY METIS_INCLUDE_DIR)
endif()
