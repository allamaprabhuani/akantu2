#===============================================================================
# @file   FindScotch.cmake
# @author Nicolas Richart <nicolas.richart@epfl.ch>
# @date   Tue Aug  25 16:53:57 2010
#
# @brief  The find_package file for Scotch
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
#if(SCOTCH_DIR)
#  set(SCOTCH_LIBRARY "NOTFOUND" CACHE INTERNAL "Cleared" FORCE)
#endif(SCOTCH_DIR)

find_library(SCOTCH_LIBRARY scotch
  PATHS ${SCOTCH_DIR}
  PATH_SUFFIXES src/libscotch lib
  )

find_library(SCOTCH_LIBRARY_ERR scotcherr
  PATHS ${SCOTCH_DIR}
  PATH_SUFFIXES src/libscotch lib
  )

find_path(SCOTCH_INCLUDE_PATH scotch.h
  PATHS ${SCOTCH_DIR}
  PATH_SUFFIXES include scotch src/libscotch include/scotch
  )

#===============================================================================
mark_as_advanced(SCOTCH_LIBRARY)
mark_as_advanced(SCOTCH_LIBRARY_ERR)
mark_as_advanced(SCOTCH_INCLUDE_PATH)

set(SCOTCH_LIBRARIES_ALL ${SCOTCH_LIBRARY} ${SCOTCH_LIBRARY_ERR})
set(SCOTCH_LIBRARIES ${SCOTCH_LIBRARIES_ALL} CACHE INTERNAL "Libraries for scotch" FORCE)

#===============================================================================
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(SCOTCH DEFAULT_MSG
  SCOTCH_LIBRARY SCOTCH_LIBRARY_ERR SCOTCH_INCLUDE_PATH)


if(SCOTCH_INCLUDE_PATH)
  file(STRINGS ${SCOTCH_INCLUDE_PATH}/scotch.h SCOTCH_INCLUDE_CONTENT)
  string(REGEX MATCH "_cplusplus" _match ${SCOTCH_INCLUDE_CONTENT})
  if(_match)
    add_definitions(-DAKANTU_SCOTCH_NO_EXTERN)
  endif()
endif()

#===============================================================================
if(NOT SCOTCH_FOUND)
  set(SCOTCH_DIR "" CACHE PATH "Location of Scotch library.")
endif(NOT SCOTCH_FOUND)
