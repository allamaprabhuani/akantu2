#===============================================================================
# @file   AkantuInstall.cmake
# @author Nicolas Richart <nicolas.richart@epfl.ch>
# @date   Wed Oct 17 2012
#
# @brief Create the files that allows users to link with Akantu in an other
# cmake project
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
# Config gen for external packages
#===============================================================================
configure_file(cmake/AkantuBuildTreeSettings.cmake.in  "${CMAKE_BINARY_DIR}/AkantuBuildTreeSettings.cmake" @ONLY)

file(WRITE "${CMAKE_BINARY_DIR}/AkantuConfigInclude.cmake" "
#===============================================================================
# @file   AkantuConfigInclude.cmake
# @author Nicolas Richart <nicolas.richart@epfl.ch>
# @date   Fri Jun 11 09:46:59 2010
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

")

foreach(_option ${PACKAGE_SYSTEM_PACKAGES_NAMES_LIST_ALL})
  list(FIND AKANTU_OPTION_LIST ${_option} _index)
  if (_index EQUAL -1)
    if(NOT "${_option}" STREQUAL "CORE")
      if(NOT AKANTU_${_option})
        set(AKANTU_${_option} OFF)
      endif()
      file(APPEND "${CMAKE_BINARY_DIR}/AkantuConfigInclude.cmake" "
set(AKANTU_HAS_${_option} ${AKANTU_${_option}})")
    endif()
  endif()
endforeach()

file(APPEND "${CMAKE_BINARY_DIR}/AkantuConfigInclude.cmake"
"

set(AKANTU_HAS_PARTITIONER  ${AKANTU_PARTITIONER})
set(AKANTU_HAS_SOLVER       ${AKANTU_SOLVER})

set(AKANTU_OPTION_LIST ${AKANTU_OPTION_LIST})
")

foreach(_option ${AKANTU_OPTION_LIST})
  file(APPEND "${CMAKE_BINARY_DIR}/AkantuConfigInclude.cmake" "
set(AKANTU_USE_${_option} ${AKANTU_${_option}})")
  if(${AKANTU_${_option}_LIBRARIES})
    file(APPEND "${CMAKE_BINARY_DIR}/AkantuConfigInclude.cmake" "
set(AKANTU_${_option}_LIBRARIES ${AKANTU_${_option}_LIBRARIES})
set(AKANTU_${_option}_INCLUDE_DIR ${AKANTU_${_option}_INCLUDE_DIR})
")
  endif()
endforeach()

# Create the AkantuConfig.cmake and AkantuConfigVersion files
get_filename_component(CONF_REL_INCLUDE_DIR "${CMAKE_INSTALL_PREFIX}" ABSOLUTE)
configure_file(cmake/AkantuConfig.cmake.in "${CMAKE_BINARY_DIR}/AkantuConfig.cmake" @ONLY)
configure_file(cmake/AkantuConfigVersion.cmake.in "${CMAKE_BINARY_DIR}/AkantuConfigVersion.cmake" @ONLY)
configure_file(cmake/AkantuUse.cmake "${CMAKE_BINARY_DIR}/AkantuUse.cmake" COPYONLY)

# Install the export set for use with the install-tree
install(FILES ${CMAKE_BINARY_DIR}/AkantuConfig.cmake
  ${CMAKE_BINARY_DIR}/AkantuConfigInclude.cmake
  ${CMAKE_BINARY_DIR}/AkantuConfigVersion.cmake
  ${CMAKE_SOURCE_DIR}/cmake/AkantuUse.cmake
  DESTINATION  lib/akantu
  COMPONENT dev)

install(FILES
  ${CMAKE_SOURCE_DIR}/cmake/FindIOHelper.cmake
  ${CMAKE_SOURCE_DIR}/cmake/FindQVIEW.cmake
  ${CMAKE_SOURCE_DIR}/cmake/FindMumps.cmake
  ${CMAKE_SOURCE_DIR}/cmake/FindScotch.cmake
  ${CMAKE_SOURCE_DIR}/cmake/FindGMSH.cmake
  DESTINATION  lib/akantu/cmake
  COMPONENT dev)
