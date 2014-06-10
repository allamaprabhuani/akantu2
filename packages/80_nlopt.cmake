#===============================================================================
# @file   cpp_array.cmake
#
# @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
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
option(AKANTU_NLOPT "Use NLOPT library" OFF)
option(AKANTU_NLOPT_AUTO_DOWNLOAD "Automatic download of the NLOPT library" ON)

if(AKANTU_NLOPT)

  if (AKANTU_NLOPT_AUTO_DOWNLOAD)
    find_package(Wget)

    set(NLOPT_VERSION "2.4.2")
    set(NLOPT_ARCHIVE "nlopt-${NLOPT_VERSION}.tar.gz")
    set(NLOPT_SOURCE_DIR "${CMAKE_BINARY_DIR}/nlopt-${NLOPT_VERSION}")
    set(NLOPT_PRODUCED_LIB "${CMAKE_BINARY_DIR}/nlopt-${NLOPT_VERSION}/.libs/libnlopt_cxx.so")
    set(NLOPT_CONFIGURE_COMMAND "./configure --enable-shared --with-cxx")

    if (CMAKE_BUILD_TYPE)
      if (${CMAKE_BUILD_TYPE} STREQUAL "Debug")
	set(NLOPT_CONFIGURE_COMMAND "${NLOPT_CONFIGURE_COMMAND} --enable-debug")
      endif()
    endif()

    set(NLOPT_NEED_REBUILD 1)
    if (NLOPT_LAST_CONFIGURE_COMMAND)
      if (${NLOPT_LAST_CONFIGURE_COMMAND} STREQUAL ${NLOPT_CONFIGURE_COMMAND})
	set(NLOPT_NEED_REBUILD 0)
      endif()
    endif()

    set(NLOPT_LAST_CONFIGURE_COMMAND ${NLOPT_CONFIGURE_COMMAND} CACHE STRING "last configuration command")

    if(NOT EXISTS ${NLOPT_PRODUCED_LIB} OR ${NLOPT_NEED_REBUILD})    
      message(STATUS "need to rebuild NLOPT")
      if(WGET_FOUND)
	if(NOT EXISTS ${CMAKE_BINARY_DIR}/${NLOPT_ARCHIVE})
	  message(STATUS "Downloading NLOPT library")
	  execute_process(COMMAND 
	    ${WGET_EXECUTABLE} http://ab-initio.mit.edu/nlopt/${NLOPT_ARCHIVE} --directory-prefix=${CMAKE_BINARY_DIR}/
	    OUTPUT_QUIET ERROR_QUIET)
	  
	endif()
      endif()
    
      if(NOT EXISTS ${CMAKE_BINARY_DIR}/${NLOPT_ARCHIVE})
	message(FATAL_ERROR "could not download NLOPT archive")
      endif()

      find_program(TAR_EXECUTABLE tar)
      if (NOT TAR_EXECUTABLE)
	message(FATAL_ERROR "need tar executable in order to uncompress NLOPT archive")
      endif()

      message(STATUS "Extract NLOPT library")
      execute_process(COMMAND ${TAR_EXECUTABLE} xfz ${NLOPT_ARCHIVE} OUTPUT_QUIET)
      if(NOT EXISTS ${NLOPT_SOURCE_DIR}/)
	message(FATAL_ERROR "could not untar NLOPT archive")
      endif()
   
      message(STATUS "Configure NLOPT library with: ${NLOPT_CONFIGURE_COMMAND} in ${NLOPT_SOURCE_DIR}")
      string(REPLACE " " ";" NLOPT_CONFIGURE_COMMAND ${NLOPT_CONFIGURE_COMMAND})
      execute_process(COMMAND ${NLOPT_CONFIGURE_COMMAND} WORKING_DIRECTORY ${NLOPT_SOURCE_DIR} OUTPUT_QUIET ERROR_QUIET)
      if(NOT EXISTS ${NLOPT_SOURCE_DIR}/Makefile)
	message(STATUS "Make NLOPT library")
	message(FATAL_ERROR "Could not configure NLOPT")
      endif()
      execute_process(COMMAND make WORKING_DIRECTORY ${NLOPT_SOURCE_DIR} OUTPUT_QUIET ERROR_QUIET)

      if(NOT EXISTS ${NLOPT_PRODUCED_LIB})
	message("NLopt was not correctly compiled")
      endif()

    endif()

    

    set(NLOPT_INTERNAL_DIR ${NLOPT_SOURCE_DIR} CACHE PATH "Location of NLOPT source directory." FORCE)
    mark_as_advanced(NLOPT_INTERNAL_DIR)

  endif()
  add_external_package(NLopt "Use NLOPT library")
else()
  set(AKANTU_USE_NLOPT ${AKANTU_NLOPT} CACHE BOOL "Use NLOPT library" FORCE)
  add_optional_external_package(NLopt "Use NLOPT library" OFF)
endif()



mark_as_advanced(AKANTU_NLOPT)
mark_as_advanced(NLOPT_DIR)
mark_as_advanced(AKANTU_NLOPT_AUTO_DOWNLOAD)
mark_as_advanced(AKANTU_USE_NLOPT)
mark_as_advanced(NLOPT_LAST_CONFIGURE_COMMAND)
mark_as_advanced(TAR_EXECUTABLE)

set(AKANTU_NLOPT_DOCUMENTATION "
When enabled this package download, configure and compiles the \\href{http://ab-initio.mit.edu/wiki/index.php/NLopt}{NLOPT} library (if access to internet is given). 
")
