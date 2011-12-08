#===============================================================================
# @file   CMakeLists.txt
# @author Nicolas Richart <nicolas.richart@epfl.ch>
# @date   Mon Nov 28 15:11:51 2011
#
# @brief  CMake file for the library
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
# Package Management
#===============================================================================

generate_source_list_from_packages(${CMAKE_CURRENT_BINARY_DIR} AKANTU_LIBRARY_SRCS AKANTU_LIBRARY_HDRS AKANTU_LIBRARY_INCLUDE_DIRS)

#===========================================================================
# header precompilation
#===========================================================================

set(AKANTU_COMMON_HDR_TO_PRECOMPILE
  common/aka_vector.hh
  common/aka_math.hh
  common/aka_types.hh
  fem/element_class.hh
  )

set(AKANTU_PRECOMPILE_HDR_INCLUDE_DIRS
  ${CMAKE_CURRENT_BINARY_DIR}/common/
  ${CMAKE_CURRENT_BINARY_DIR}/fem/
  )


list(APPEND AKANTU_INCLUDE_DIRS
  ${AKANTU_LIBRARY_INCLUDE_DIRS}
  ${AKANTU_PRECOMPILE_HDR_INCLUDE_DIRS})

set(AKANTU_INCLUDE_DIRS ${AKANTU_INCLUDE_DIRS} PARENT_SCOPE)

include_directories(${AKANTU_INCLUDE_DIRS} ${AKANTU_EXTERNAL_LIB_INCLUDE_DIR})


if(CMAKE_COMPILER_IS_GNUCXX)
  include(${CMAKE_SOURCE_DIR}/cmake/PCHgcc.cmake)
  foreach(_header ${AKANTU_COMMON_HDR_TO_PRECOMPILE})
    add_pch_rule(${_header} AKANTU_LIBRARY_SRCS)
  endforeach()
endif()

#===============================================================================
# Library generation
#===============================================================================

add_library(akantu ${AKANTU_LIBRARY_SRCS})

add_definitions(-fkeep-inline-functions)
set_target_properties(akantu PROPERTIES ${AKANTU_LIBRARY_PROPERTIES})
set_target_properties(akantu PROPERTIES COMPILE_DEFINITIONS __aka_inline__=inline)
set_target_properties(akantu PROPERTIES PUBLIC_HEADER "${AKANTU_LIBRARY_HDRS}")
#===============================================================================
# Install rules
#===============================================================================

install(TARGETS akantu
  EXPORT AkantuLibraryDepends
  LIBRARY DESTINATION lib COMPONENT dev
  PUBLIC_HEADER DESTINATION include/akantu/ COMPONENT dev
  )
