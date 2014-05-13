#===============================================================================
# @file   cpp_array.cmake
#
# @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
# @author Nicolas Richart <nicolas.richart@epfl.ch>
#
# @date   Mon Nov 21 18:19:15 2011
#
# @brief  package description for cpp_array project
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
option(AKANTU_CPPARRAY "Use cpp-array library" OFF)
find_package(Subversion)
mark_as_advanced(AKANTU_CPPARRAY)

if(SUBVERSION_FOUND)
  if(AKANTU_CPPARRAY)
    if(EXISTS ${PROJECT_SOURCE_DIR}/third-party/cpp-array)
      execute_process(
	COMMAND ${Subversion_SVN_EXECUTABLE} up ${PROJECT_SOURCE_DIR}/third-party/cpp-array
	OUTPUT_VARIABLE _revision)
      string(REGEX REPLACE ".*At revision ([0-9]*)\\..*" "\\1" _rev "${_revision}")
      message(STATUS "Updating Cpp-Array: r${_rev}")
    else()
      message(STATUS "Checking out Cpp-Array")
      execute_process(
	COMMAND ${Subversion_SVN_EXECUTABLE} co http://cpp-array.googlecode.com/svn/trunk ${PROJECT_SOURCE_DIR}/third-party/cpp-array
	OUTPUT_QUIET)
    endif()

    add_subdirectory(${PROJECT_SOURCE_DIR}/third-party/cpp-array/)
    
#    set(cpp-array_TESTS OFF CACHE BOOL "cpparray tests" FORCE)

    list(APPEND AKANTU_EXTERNAL_LIB_INCLUDE_DIR ${cpp-array_INCLUDE_DIRS})

    list(APPEND CPACK_SOURCE_IGNORE_FILES ${PROJECT_SOURCE_DIR}/third-party/cpp-array/)

    set(AKANTU_CPPARRAY_INCLUDE_DIR ${CMAKE_PREFIX}/include)
    list(APPEND AKANTU_OPTION_LIST CPPARRAY)

  endif()
else()
  set(AKANTU_USE_CPPARRAY ${AKANTU_CPPARRAY} CACHE BOOL "Use cpp-array library" FORCE)
  add_optional_external_package(CppArray "Use cpp-array library" OFF)
endif()
