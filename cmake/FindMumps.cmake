#===============================================================================
# @file   FindMumps.cmake
# @author Nicolas Richart <nicolas.richart@epfl.ch>
# @date   Sun Dec 12 18:35:02 2010
#
# @brief  The find_package file for the Mumps solver
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
find_library(MUMPS_LIBRARY_DMUMPS NAME dmumps_scotch
   PATHS ${MUMPS_DIR} /usr
   PATH_SUFFIXES lib
   )

find_library(MUMPS_LIBRARY_COMMON NAME mumps_common_scotch
   PATHS ${MUMPS_DIR}
   PATH_SUFFIXES lib
   )

find_library(MUMPS_LIBRARY_PORD NAME pord_scotch
   PATHS ${MUMPS_DIR}
   PATH_SUFFIXES lib
   )


find_path(MUMPS_INCLUDE_PATH dmumps_c.h
  PATHS ${MUMPS_DIR}
  PATH_SUFFIXES include
  )

set(MUMPS_LIBRARIES
  ${MUMPS_LIB_COMMON}
  ${MUMPS_LIB_DMUMPS}
  ${MUMPS_LIB_PORD}
  CACHE INTERNAL "Libraries for mumps" FORCE
  )

#===============================================================================
mark_as_advanced(MUMPS_LIBRARY_COMMON)
mark_as_advanced(MUMPS_LIBRARY_DMUMPS)
mark_as_advanced(MUMPS_LIBRARY_PROD)
mark_as_advanced(MUMPS_INCLUDE_PATH)

set(MUMPS_LIBRARIES_ALL ${MUMPS_LIBRARY_DMUMPS} ${MUMPS_LIBRARY_COMMON} ${MUMPS_LIBRARY_PROD})
set(MUMPS_LIBRARIES ${MUMPS_LIBRARIES_ALL} CACHE INTERNAL "Libraries for mumps" FORCE)

#===============================================================================
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(MUMPS DEFAULT_MSG
  MUMPS_LIBRARIES MUMPS_INCLUDE_PATH)

#===============================================================================
if(NOT MUMPS_FOUND)
  set(MUMPS_DIR "" CACHE PATH "Prefix of mumps library.")
endif(NOT MUMPS_FOUND)
