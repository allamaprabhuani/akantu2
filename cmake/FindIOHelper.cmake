#===============================================================================
# @file   FindIOHelper.cmake
# @author Nicolas Richart <nicolas.richart@epfl.ch>
# @date   Tue Aug  3 16:29:57 2010
#
# @brief  The find_package file for IOHelper
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
#set(IOHELPER_LIBRARY "NOTFOUND" CACHE INTERNAL "Cleared" FORCE)
find_library(IOHELPER_LIBRARY IOHelper
  PATHS ${IOHELPER_DIR}
  PATH_SUFFIXES src
  )

find_path(IOHELPER_INCLUDE_PATH io_helper.h
  PATHS ${IOHELPER_DIR}
  PATH_SUFFIXES src
  )

#===============================================================================
mark_as_advanced(IOHELPER_LIBRARY)
mark_as_advanced(IOHELPER_INCLUDE_PATH)

#===============================================================================
find_package(ZLIB REQUIRED)

set(IOHELPER_LIBRARIES_ALL ${IOHELPER_LIBRARY} ${ZLIB_LIBRARIES})
set(IOHELPER_LIBRARIES ${IOHELPER_LIBRARIES_ALL} CACHE INTERNAL "Libraries for IOHelper" FORCE)

#===============================================================================
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(IOHELPER DEFAULT_MSG IOHELPER_LIBRARY IOHELPER_INCLUDE_PATH)

#===============================================================================
if(NOT IOHELPER_FOUND)
  set(IOHELPER_DIR "" CACHE PATH "Location of IOHelper source directory.")
endif(NOT IOHELPER_FOUND)

