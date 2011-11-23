#===============================================================================
# @file   CMakeLists.txt
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

#===============================================================================
#            _    ,-,    _
#     ,--, /: :\/': :`\/: :\
#    |`;  ' `,'   `.;    `: |                  _               _
#    |    |     |  '  |     |.                | |             | |
#    | :  | F E | A S | T++ ||            __ _| | ____ _ _ __ | |_ _   _
#    | :. |  :  |  :  |  :  | \          / _` | |/ / _` | '_ \| __| | | |
#     \__/: :.. : :.. | :.. |  )        | (_| |   < (_| | | | | |_| |_| |
#          `---',\___/,\___/ /'          \__,_|_|\_\__,_|_| |_|\__|\__,_|
#               `==._ .. . /'
#                    `-::-'
#===============================================================================

#===============================================================================
# CMake Project
#===============================================================================

cmake_minimum_required(VERSION 2.6)

project(AKANTU)

enable_language(CXX)

#===============================================================================
# Misc.
#===============================================================================

set(AKANTU_CMAKE_DIR "${AKANTU_SOURCE_DIR}/cmake")
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake")
set(BUILD_SHARED_LIBS ON CACHE BOOL "Build shared libraries.")
INCLUDE(${AKANTU_CMAKE_DIR}/AkantuMacros.cmake)

#===============================================================================
# Version Number
#===============================================================================

# AKANTU version number.  An even minor number corresponds to releases.
set(AKANTU_MAJOR_VERSION 0)
set(AKANTU_MINOR_VERSION 1)
execute_process(COMMAND ./get_revision.sh WORKING_DIRECTORY ${AKANTU_SOURCE_DIR} OUTPUT_VARIABLE AKANTU_BUILD_VERSION)
string(STRIP ${AKANTU_BUILD_VERSION} AKANTU_BUILD_VERSION)
MESSAGE(${AKANTU_BUILD_VERSION})
set(AKANTU_VERSION
  "${AKANTU_MAJOR_VERSION}.${AKANTU_MINOR_VERSION}.${AKANTU_BUILD_VERSION}"
  )
# Append the library version information to the library target properties
if(NOT AKANTU_NO_LIBRARY_VERSION)
  set(AKANTU_LIBRARY_PROPERTIES ${AKANTU_LIBRARY_PROPERTIES}
    VERSION "${AKANTU_VERSION}"
    SOVERSION "${AKANTU_MAJOR_VERSION}.${AKANTU_MINOR_VERSION}"

    )
endif(NOT AKANTU_NO_LIBRARY_VERSION)
#===============================================================================
# Package builder target module of CPack
#===============================================================================
SET(CPACK_GENERATOR "DEB")
SET(CPACK_DEBIAN_PACKAGE_MAINTAINER "Guillaume Anciaux") 
SET(CPACK_DEBIAN_PACKAGE_DEPENDS "libc6")
SET(CPACK_DEBIAN_PACKAGE_ARCHITECTURE "amd64")
SET(CPACK_DEBIAN_PACKAGE_DESCRIPTION "Akantu library")
SET(CPACK_PACKAGE_VENDOR "LSMS")
SET(CPACK_PACKAGE_VERSION ${AKANTU_VERSION})
INCLUDE(CPack)

#===============================================================================
# Options
#===============================================================================
option(AKANTU_DEBUG "Compiles akantu with the debug messages" ON)
option(AKANTU_DOCUMENTATION "Build source documentation using Doxygen." OFF)

# features
add_optional_package(IOHelper "Add IOHelper support in akantu" OFF)
add_optional_package(QVIEW "Add QView support in akantu" OFF)
add_optional_package(BLAS "Use BLAS for arithmetic operations" OFF)
add_optional_package(EPSN "Use EPSN streering environment" OFF)


# Debug
if(NOT AKANTU_DEBUG)
  add_definitions(-DAKANTU_NDEBUG)
endif(NOT AKANTU_DEBUG)

#===========================================================================
# Dependencies
#===========================================================================
find_package(Boost REQUIRED)

if(Boost_FOUND)
  include_directories(${Boost_INCLUDE_DIRS})
endif()

#check_for_isnan(AKANTU_ISNAN)
#message("ISNAN = "${AKANTU_ISNAN})


#===============================================================================
# File package section
#===============================================================================

file(GLOB AKANTU_PACKAGE_LIST RELATIVE ${AKANTU_SOURCE_DIR} "packages/*.cmake")

foreach(pkg ${AKANTU_PACKAGE_LIST})
  get_filename_component(basename ${pkg} NAME)
  string(REGEX REPLACE "\\.cmake" "" OPTION_NAME ${basename})
  string(TOUPPER "${OPTION_NAME}" OPTION_NAME)
  include(${pkg})


  # differentiate the file types
  foreach(file ${${OPTION_NAME}_FILES})
    set(RELEASE_ALL_FILE ${RELEASE_ALL_FILE} ${file})  
    if (AKANTU_${OPTION_NAME})
      if    (${file} MATCHES ".*inline.*\\.cc")
	set(${OPTION_NAME}_INLINE ${${OPTION_NAME}_INLINE} ${file})
      elseif(${file} MATCHES ".*\\.hh")
	set(${OPTION_NAME}_HEADERS ${${OPTION_NAME}_HEADERS} ${file})
      else  ()
	set(${OPTION_NAME}_SRC ${${OPTION_NAME}_SRC} ${file})
      endif()
    endif()
  endforeach()


  # generates the include directory variable
  foreach(inc ${${OPTION_NAME}_HEADERS})
    get_filename_component(INC_DIR ${inc} PATH)    
    set(${OPTION_NAME}_INCLUDE_DIRS ${${OPTION_NAME}_INCLUDE_DIRS} ${INC_DIR})
    list(REMOVE_DUPLICATES ${OPTION_NAME}_INCLUDE_DIRS)
  endforeach()

  # generate global lists for akantu to know what to build
  if(AKANTU_${OPTION_NAME})
    set(AKANTU_COMMON_SRC
      ${AKANTU_COMMON_SRC}
      ${${OPTION_NAME}_SRC}
      )
    set(AKANTU_COMMON_HEADERS
      ${AKANTU_COMMON_HEADERS}
      ${${OPTION_NAME}_HEADERS}
      )
    set(AKANTU_INCLUDE_DIRS
      ${AKANTU_INCLUDE_DIRS}
      ${${OPTION_NAME}_INCLUDE_DIRS})
  endif()
 
endforeach()

#MESSAGE(${RELEASE_ALL_FILE})
#MESSAGE(${AKANTU_COMMON_SRC})

#check if there are some file in the release that are not registered in a package
file(GLOB_RECURSE all_files RELATIVE ${AKANTU_SOURCE_DIR} "*.cc" "*.hh")
foreach(file ${all_files})
 if(NOT ${file} MATCHES "test.*" AND NOT ${file} MATCHES "example.*")
   list(FIND RELEASE_ALL_FILE ${file} index)
   if (index EQUAL -1)
     MESSAGE("file ${file} is not registered in any package.")
     MESSAGE("Please append the file in one of the files within directory ${AKANTU_SOURCE_DIR}/packages")
   endif()
 endif()
endforeach()

#check if there are some file in the package list that are not in the release file
file(GLOB_RECURSE all_files RELATIVE ${AKANTU_SOURCE_DIR} "*.cc" "*.hh")
foreach(file ${RELEASE_ALL_FILE})
 if(NOT ${file} MATCHES "test.*" AND NOT ${file} MATCHES "example.*")
   list(FIND  all_files ${file} index)
   if (index EQUAL -1)
     MESSAGE("file ${file} is registered in packages but not existing in the release.")
   endif()
 endif()
endforeach()

#check dependencies
foreach(pkg ${AKANTU_PACKAGE_LIST})
  get_filename_component(basename ${pkg} NAME)
  string(REGEX REPLACE "\\.cmake" "" OPTION_NAME ${basename})
  string(TOUPPER "${OPTION_NAME}" OPTION_NAME)
  # differentiate the file types
  if (AKANTU_${OPTION_NAME})
    foreach(dep ${${OPTION_NAME}_DEPENDS})
      if (NOT AKANTU_${dep})
	MESSAGE(FATAL_ERROR "Package ${pkg} depends on package ${dep}. You need to activate it to make it work")
      endif()
    endforeach()
  endif()
endforeach()


#===========================================================================
# header precompilation
#===========================================================================

set(AKANTU_COMMON_HDR_TO_PRECOMPILE
  common/aka_vector.hh
  common/aka_math.hh
  common/aka_types.hh
  fem/by_element_type.hh
  fem/element_class.hh
  )

set(AKANTU_PRECOMPILE_HDR_INCLUDE_DIRS
  ${CMAKE_BINARY_DIR}/common/
  ${CMAKE_BINARY_DIR}/fem/
  )

include_directories(
  ${AKANTU_INCLUDE_DIRS}
  ${AKANTU_PRECOMPILE_HDR_INCLUDE_DIRS}
  ${AKANTU_EXTERNAL_LIB_INCLUDE_PATH}
  )


if(CMAKE_COMPILER_IS_GNUCXX)
  include(cmake/PCHgcc.cmake)
  foreach(_header ${AKANTU_COMMON_HDR_TO_PRECOMPILE})
    add_pch_rule(${_header} AKANTU_COMMON_SRC)
  endforeach()
endif()

add_library(akantu ${AKANTU_COMMON_SRC})

add_definitions(-fkeep-inline-functions)
set_target_properties(akantu PROPERTIES ${AKANTU_LIBRARY_PROPERTIES})
set_target_properties(akantu
  PROPERTIES COMPILE_DEFINITIONS __aka_inline__=inline)

#include(cmake/PCHSupport.cmake)
#if(PCHSupport_FOUND)
#  foreach(_header ${AKANTU_COMMON_HDR_TO_PRECOMPILE})
#    add_precompiled_header(akantu ${_header})
#  endforeach()
#endif()


#===========================================================================
# Documentation
#===========================================================================
if(AKANTU_DOCUMENTATION)
  find_package(Doxygen REQUIRED)

  if(DOXYGEN_FOUND)
    set(DOXYGEN_WARNINGS NO)
    set(DOXYGEN_QUIET YES)
    if(CMAKE_VERBOSE_MAKEFILE)
      set(DOXYGEN_WARNINGS YES)
      set(DOXYGEN_QUIET NO)
    endif(CMAKE_VERBOSE_MAKEFILE)

    add_subdirectory(doc/doxygen)
  endif(DOXYGEN_FOUND)
endif(AKANTU_DOCUMENTATION)

#===============================================================================
# Tests
#===============================================================================
ENABLE_TESTING()
INCLUDE(CTest)

INCLUDE(${AKANTU_CMAKE_DIR}/AkantuTestAndExamples.cmake)

option(AKANTU_TESTS "Activate tests" OFF)
if(AKANTU_TESTS)
  find_package(GMSH REQUIRED)
  add_subdirectory(test)
endif(AKANTU_TESTS)

#===============================================================================
# Examples
#===============================================================================
option(AKANTU_EXAMPLES "Activate examples" OFF)
if(AKANTU_EXAMPLES)
  find_package(GMSH REQUIRED)
  add_subdirectory(examples)
endif(AKANTU_EXAMPLES)


#===============================================================================
# Install rules
#===============================================================================

install(TARGETS akantu 
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib)


install(FILES 
  ${AKANTU_COMMON_HEADERS}
  DESTINATION include/akantu/
)

