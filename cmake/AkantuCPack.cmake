#===============================================================================
# @file   AkantuCPack.cmake
#
# @author Nicolas Richart <nicolas.richart@epfl.ch>
#
# @date   Wed Oct 17 15:19:18 2012
#
# @brief  Configure the packaging system
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

set(PACKAGE_FILE_NAME "akantu" CACHE STRING "Name of package to be generated")
mark_as_advanced(PACKAGE_FILE_NAME)

set(CPACK_GENERATOR "DEB;TGZ;TBZ2;STGZ")

# Debian config package
set(CPACK_DEBIAN_PACKAGE_MAINTAINER "guillaume.anciaux@epfl.ch, nicolas.richart@epfl.ch")
if(CMAKE_SIZEOF_VOID_P EQUAL 8)
  set(CPACK_DEBIAN_PACKAGE_ARCHITECTURE "amd64" CACHE STRING "Architecture of debian package generation")
else()
  set(CPACK_DEBIAN_PACKAGE_ARCHITECTURE "i386" CACHE STRING "Architecture of debian package generation")
endif()
set(CPACK_DEBIAN_PACKAGE_DESCRIPTION "Akantu library")
set(CPACK_DEBIAN_PACKAGE_DEPENDS "${PACKAGE_SYSTEM_DEBIAN_PACKAGE_DEPENDS}")

# General configuration
set(CPACK_PACKAGE_VENDOR "LSMS")
set(CPACK_PACKAGE_FILE_NAME "${PACKAGE_FILE_NAME}-${AKANTU_VERSION}-${CPACK_DEBIAN_PACKAGE_ARCHITECTURE}")
set(CPACK_PACKAGE_VERSION "${AKANTU_VERSION}")

set(CPACK_COMPONENTS_ALL lib dev)
set(CPACK_COMPONENT_LIB_DISPLAY_NAME "Libraries")
set(CPACK_COMPONENT_DEV_DISPLAY_NAME "C++ Headers")
set(CPACK_COMPONENT_DEV_DEPENDS lib)
 set(CPACK_COMPONENT_LIB_DESCRIPTION
   "Akantu libraries")
 set(CPACK_COMPONENT_DEV_DESCRIPTION
   "Akantu C/C++ header files")
set(CPACK_COMPONENT_LIB_GROUP "Akantu Libraries")
set(CPACK_COMPONENT_DEV_GROUP "Development")

set(CPACK_SOURCE_PACKAGE_FILE_NAME "${PACKAGE_FILE_NAME}-${AKANTU_VERSION}-src")
set(CPACK_RESOURCE_FILE_LICENSE "${PROJECT_SOURCE_DIR}/COPYING")

list(APPEND CPACK_SOURCE_IGNORE_FILES ${AKANTU_EXCLUDE_SOURCE_FILE} ${AKANTU_TESTS_EXCLUDE_FILES} ${AKANTU_DOC_EXCLUDE_FILES})
foreach(_pkg ${PACKAGE_SYSTEM_PACKAGES_OFF})
  list(APPEND CPACK_SOURCE_IGNORE_FILES ${CMAKE_SOURCE_DIR}/packages/${_pkg}.cmake)
endforeach()
list(APPEND CPACK_SOURCE_IGNORE_FILES "/doc/manual/;/.*build.*/;/CVS/;/\\\\.svn/;/\\\\.bzr/;/\\\\.hg/;/\\\\.git/;\\\\.swp$;\\\\.#;/#;~")

include(CPack)
