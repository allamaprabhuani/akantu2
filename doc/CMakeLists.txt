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


file(GLOB_RECURSE __all_files "*.*")
set(_all_files)
foreach(_file ${__all_files})
  if("${_file}" MATCHES  "${PROJECT_SOURCE_DIR}/doc/manual")
    file(RELATIVE_PATH __file ${PROJECT_SOURCE_DIR} ${_file})
    list(APPEND _all_files ${__file})
  endif()
endforeach()

set(_all_manual_files)
foreach(_pkg ${${_project}_PACKAGE_SYSTEM_PACKAGES_ON})
  string(TOUPPER "${_pkg}" _pkg)
  foreach(_f ${AKANTU_${_pkg}_MANUAL_FILES})
    list(APPEND _all_manual_files doc/manual/${_f})
  endforeach()
endforeach()

set(_akantu_doc_exclude_files)
foreach(_file ${_all_files})
  list(FIND _all_manual_files ${_file} _ret)
  if(_ret EQUAL -1)
    list(APPEND _akantu_doc_exclude_files /${_file})
  endif()
endforeach()

set(AKANTU_DOC_EXCLUDE_FILES ${_akantu_doc_exclude_files} CACHE INTERNAL "Documentation files to excluse from Akantu Package" FORCE)
