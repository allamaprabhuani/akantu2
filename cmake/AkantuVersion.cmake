#===============================================================================
# @file   AkantuVersion.cmake
#
# @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
#
# @date   Wed Oct 17 15:19:18 2012
#
# @brief  Handle the generation of version numbers
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

find_package(Subversion)

if(EXISTS ${CMAKE_SOURCE_DIR}/.svn)
  if(SUBVERSION_FOUND)
    subversion_wc_info(${CMAKE_SOURCE_DIR} MY)
    set(AKANTU_BUILD_VERSION ${MY_WC_REVISION})
    set(AKANTU_VERSION
      "${AKANTU_MAJOR_VERSION}.${AKANTU_MINOR_VERSION}.${AKANTU_PATCH_VERSION}.${AKANTU_BUILD_VERSION}"
      )
    file(WRITE VERSION "${AKANTU_VERSION}\n")
  else(SUBVERSION_FOUND)
    message("SVN control files were found but no subversion executable is present... ")
    set(AKANTU_VERSION 0)
  endif(SUBVERSION_FOUND)
else(EXISTS ${CMAKE_SOURCE_DIR}/.svn)
  if(EXISTS ${CMAKE_SOURCE_DIR}/VERSION)
    file(STRINGS VERSION AKANTU_VERSION)
  else(EXISTS ${CMAKE_SOURCE_DIR}/VERSION)
    message("No SVN control file neither VERSION file could be found. How was this release made ?")
  endif(EXISTS ${CMAKE_SOURCE_DIR}/VERSION)
endif(EXISTS ${CMAKE_SOURCE_DIR}/.svn)


# Append the library version information to the library target properties
if(NOT AKANTU_NO_LIBRARY_VERSION)
  set(AKANTU_LIBRARY_PROPERTIES ${AKANTU_LIBRARY_PROPERTIES}
    VERSION "${AKANTU_VERSION}"
    SOVERSION "${AKANTU_MAJOR_VERSION}.${AKANTU_MINOR_VERSION}"
    )
endif(NOT AKANTU_NO_LIBRARY_VERSION)
