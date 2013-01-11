#===============================================================================
# @file   CMakeLists.txt
#
# @author Nicolas Richart <nicolas.richart@epfl.ch>
#
# @date   Wed Oct 17 15:19:18 2012
#
# @brief  Build the documentation
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

find_package(Doxygen REQUIRED)

if(DOXYGEN_FOUND)
  set(DOXYGEN_WARNINGS NO)
  set(DOXYGEN_QUIET YES)
  if(CMAKE_VERBOSE_MAKEFILE)
    set(DOXYGEN_WARNINGS YES)
    set(DOXYGEN_QUIET NO)
  endif(CMAKE_VERBOSE_MAKEFILE)

  add_subdirectory(doxygen)
endif(DOXYGEN_FOUND)


foreach(_pkg ${PACKAGE_SYSTEM_PACKAGES_OFF})
  package_pkg_name(${_pkg} _package_name)
  foreach(_file ${${_package_name}_DOC})
    list(APPEND _AKANTU_DOC_EXCLUDE_FILES doc/${_file})
  endforeach()
endforeach()

set(AKANTU_DOC_EXCLUDE_FILES ${_AKANTU_DOC_EXCLUDE_FILES} CACHE INTERNAL "Documentation files to excluse from Akantu Package" FORCE)