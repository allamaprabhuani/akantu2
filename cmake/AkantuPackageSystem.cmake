#===============================================================================
# @file   AkantuMacros.cmake
# @author Nicolas Richart <nicolas.richart@epfl.ch>
# @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
# @date   Wed Feb  9 10:59:42 2011
#
# @brief  Set of macros used by akantu to handle the package system
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
# Package Management
#===============================================================================

set(PACKAGE_SYSTEM_PKG_PREFIX AKANTU_)
set(PACKAGE_SYSTEM_OPT_PREFIX AKANTU_USE_)
set(PACKAGE_SYSTEM_DESC_PREFIX AKANTU_DESC_)

macro(package_pkg_name PKG PKG_NAME)
  string(TOUPPER ${PKG} _u_package)
  set(${PKG_NAME} ${PACKAGE_SYSTEM_PKG_PREFIX}${_u_package})
endmacro()

macro(package_opt_name PKG OPT_NAME)
  string(TOUPPER ${PKG} _u_package)
  set(${OPT_NAME} ${PACKAGE_SYSTEM_OPT_PREFIX}${_u_package})
endmacro()

macro(package_desc_name PKG DESC_NAME)
  string(TOUPPER ${PKG} _u_package)
  set(${DESC_NAME} ${PACKAGE_SYSTEM_DESC_PREFIX}${_u_package})
endmacro()

#===============================================================================
macro(add_all_packages package_dir)
  akantu_message("add_all_packages: PKG DIR : ${package_dir}")
  file(GLOB _akantu_package_list "${package_dir}/*.cmake")

  foreach(_pkg ${_akantu_package_list})
    get_filename_component(_basename ${_pkg} NAME)
    string(REGEX REPLACE "\\.cmake" "" _option_name ${_basename})
    string(TOUPPER "${_option_name}" _option_name)
    list(APPEND PACKAGE_SYSTEM_PACKAGES_NAMES_LIST_ALL ${_option_name})
  endforeach()

  akantu_message("add_all_packages: PKG LIST : "${PACKAGE_SYSTEM_PACKAGES_NAMES_LIST_ALL})

  foreach(_pkg ${PACKAGE_SYSTEM_PACKAGES_NAMES_LIST_ALL})
    akantu_message("add_all_packages: including ${_pkg}")
    string(TOLOWER "${_pkg}" _pkg)
    include(${package_dir}/${_pkg}.cmake)
    package_pkg_name(${_pkg} _package_name)
    if (${_package_name})
      list(APPEND PACKAGE_SYSTEM_PACKAGES_ON ${_pkg})
    else (${_package_name})
      list(APPEND PACKAGE_SYSTEM_PACKAGES_OFF ${_pkg})
    endif()

    foreach(_file ${${_pkg_name}_FILES})
      list(APPEND _release_all_file ${_file})
    endforeach()
  endforeach()

  akantu_message("add_all_packages: ON  PKG : ${PACKAGE_SYSTEM_PACKAGES_ON}")
  akantu_message("add_all_packages: ALL FILE LIST : ${_release_all_file}")

  #check if there are some file in the release that are not registered in a package
  file(GLOB_RECURSE _all_files RELATIVE ${CMAKE_SOURCE_DIR}/src "*.cc" "*.hh")
  foreach(_file ${all_files})
    if(NOT ${_file} MATCHES "test.*")
      list(FIND _release_all_file ${_file} _index)
      if (_index EQUAL -1)
        message("The file ${file} is not registered in any package.")
        message("Please append the file in one of the files within directory ${CMAKE_SOURCE_DIR}/packages")
      endif()
    endif()
  endforeach()

  #check if there are some file in the package list that are not on the current directory
  foreach(file ${all_files})
    list(FIND RELEASE_ALL_FILE ${file} index)
    if (index EQUAL -1)
      message("file ${file} is not registered in any package.")
      message("Please append the file in one of the files within directory ${AKANTU_SOURCE_DIR}/packages")
    endif()
  endforeach()

  foreach(_file ${_release_all_file})
    list(FIND _all_files ${_file} _index)
    if (_index EQUAL -1)
      message("The file ${_file} is registered in packages but is not present in the source directory.")
    endif()
  endforeach()

  #construct list of files for unactivated packages
  foreach(_pkg ${PACKAGE_SYSTEM_PACKAGES_OFF})
    package_pkg_name(${_pkg} _pkg_name)
    foreach(_file ${${_pkg_name}_FILES})
      #      string(REGEX REPLACE "\\/" "\\\\\\\\/" __file ${_file})
      akantu_message("add_all_packages: ${_file} ${__file}")
      #      list(APPEND _exclude_source_file "/${__file}/")
      list(APPEND AKANTU_EXCLUDE_SOURCE_FILE ${_file})
    endforeach()
  endforeach()

  #check dependencies
  foreach(_pkg ${PACKAGE_SYSTEM_PACKAGES_OFF})
    # differentiate the file types
    akantu_message("add_all_packages: DEPENDS PKG : ${_pkg}")
    akantu_message("add_all_packages: DEPENDS LST : ${${_pkg}_DEPENDS}")
    package_pkg_name(${_pkg} _pkg_name)
    if (NOT "${${_pkg_name}_DEB_DEPEND}" STREQUAL "")
      set(deb_depend "${deb_depend}, ${${_pkg}_DEB_DEPEND}")
    endif()
  endforeach()
  set(PACKAGE_SYSTEM_DEBIAN_PACKAGE_DEPENDS "${deb_depend}")
endmacro()

#===============================================================================
macro(generate_source_list_from_packages source_dir source_files headers_files include_dirs)
  set(deb_depend "libc6")

  akantu_message("generate_source_list_from_packages: SRC DIR : ${source_dir}")
  foreach(_pkg ${PACKAGE_SYSTEM_PACKAGES_ON})
    # differentiate the file types
    package_pkg_name(${_pkg} _package_name)
    foreach(_file ${${_package_name}_FILES})
      if(${_file} MATCHES ".*inline.*\\.cc")
#        list(APPEND ${_package_name}_inlines ${_file})
        list(APPEND ${_package_name}_headers ${_file})
      elseif(${_file} MATCHES ".*\\.hh")
        list(APPEND ${_package_name}_headers ${_file})
      else()
        list(APPEND ${_package_name}_srcs ${_file})
      endif()
    endforeach()

    # generates the include directory variable
    foreach(_file ${${_package_name}_headers})
      get_filename_component(_absolute_name ${_file} ABSOLUTE)
      get_filename_component(_include_dir ${_absolute_name} PATH)
      list(APPEND ${_package_name}_include_dirs ${_include_dir})
      list(REMOVE_DUPLICATES ${_package_name}_include_dirs)
    endforeach()

    # generate global lists for akantu to know what to build
    list(APPEND ${source_files} ${${_package_name}_srcs})
    list(APPEND ${headers_files} ${${_package_name}_headers})
    list(APPEND ${include_dirs}  ${${_package_name}_include_dirs})

    akantu_message("generate_source_list_from_packages: PKG ${_package_name} SRCS : ${${_package_name}_srcs}")
    akantu_message("generate_source_list_from_packages: PKG ${_package_name} HRDS : ${${_package_name}_headers}")
    akantu_message("generate_source_list_from_packages: PKG ${_package_name} INCS : ${${_package_name}_include_dirs}")
  endforeach()

  akantu_message("generate_source_list_from_packages: SRCS : ${${source_files}}")
  akantu_message("generate_source_list_from_packages: HRDS : ${${headers_files}}")
  akantu_message("generate_source_list_from_packages: INCS : ${${include_dirs}}")
endmacro()

#===============================================================================
# macro to include optional packages
macro(add_optional_package PACKAGE DESC DEFAULT)
  cmake_parse_arguments(_opt_pkg "" "LANGUAGE" "DEPENDS;PREFIX" ${ARGN})

  package_pkg_name (${PACKAGE} _pkg_name)
  package_opt_name (${PACKAGE} _option_name)
  package_desc_name(${PACKAGE} _desc_name)

  akantu_message("add_optional_package: Registering ${PACKAGE} ${DESC} -> ${_option_name}")

  if(_opt_pkg_PREFIX)
    set(_package_prefix ${_opt_pkg_PREFIX})
  else()
    string(TOUPPER ${PACKAGE} _u_package)
    set(_package_prefix ${_u_package})
  endif()

  set(${_desc_name} ${DESC})
  option(${_option_name} ${DESC} ${DEFAULT})

  if(${_option_name})
    if(_opt_pkg_LANGUAGE)
      foreach(_language ${_opt_pkg_LANGUAGE})
	akantu_message("add_optional_package: Package ${PACKAGE} asked for language ${_language}")
	enable_language(${_language})
      endforeach()
    endif()

    foreach(_dep ${_opt_pkg_DEPENDS})
      add_package_dependecies(${PACKAGE} ${_dep})
    endforeach()

    find_package(${PACKAGE} REQUIRED)

    foreach(_prefix ${_package_prefix})
      if(${_prefix}_FOUND)
	list(APPEND AKANTU_DEFINITIONS ${_option_name})
	if(DEFINED ${_prefix}_INCLUDE_DIR)
          list(APPEND AKANTU_EXTERNAL_LIB_INCLUDE_DIR ${${_prefix}_INCLUDE_DIR})
          set(${_pkg_name}_INCLUDE_DIR ${${_prefix}_INCLUDE_DIR})
	else()
          list(APPEND AKANTU_EXTERNAL_LIB_INCLUDE_DIR ${${_prefix}_INCLUDE_PATH})
          set(${_pkg_name}_INCLUDE_DIR ${${_prefix}_INCLUDE_PATH})
	endif()
	list(APPEND AKANTU_EXTERNAL_LIBRARIES ${${_prefix}_LIBRARIES})
	set(${_pkg_name}_LIBRARIES ${${_prefix}_LIBRARIES})
	set(${_pkg_name} ON)
        string(TOUPPER ${PACKAGE} _u_package)
	list(APPEND AKANTU_OPTION_LIST ${_u_package})
	akantu_message("add_optional_package: Package ${PACKAGE} found!")
	akantu_message("add_optional_package: Package ${PACKAGE} includes : ${${_pkg_name}_INCLUDE_DIR}")
	akantu_message("add_optional_package: Package ${PACKAGE} libraries: ${${_pkg_name}_LIBRARIES}")
	akantu_message("add_optional_package: option list: ${AKANTU_OPTION_LIST}")
      else(${_prefix}_FOUND)
	akantu_message("add_optional_package: Package ${PACKAGE} not found!")
	set(${_pkg_name} OFF)
      endif(${_prefix}_FOUND)
    endforeach()
  endif(${_option_name})
endmacro()


#===============================================================================
macro(add_package_dependecies PKG DEP)
  package_pkg_name (${PKG} _opt_name)
  package_opt_name (${DEP} _dep_name)
  package_desc_name(${DEP} _var_dep_desc)

  akantu_message("add_package_dependecies: add dependence between ${_opt_name} and ${_dep_name}")
  set(_dep_desc ${_var_dep_desc})

  akantu_message("add_package_dependecies: ON dependencies of ${_dep_name} are: ${${_dep_name}_DEPS}")
  akantu_message("add_package_dependecies: saved value for ${_dep_name} is: ${${_dep_name}_OLD}")
  if(${_opt_name})
    if("${${_dep_name}_DEPS}" STREQUAL "")
      akantu_message("add_package_dependecies: Save dep state ${_dep_name}:${${_dep_name}}")
      set(${_dep_name}_OLD ${${_dep_name}} CACHE INTERNAL "${_dep_desc}" FORCE)
    endif()

    akantu_message("add_package_dependecies: force value to ON ${_dep_name}")
    set(${_dep_name} ON CACHE BOOL "${_dep_desc}" FORCE)

    list(FIND ${_dep_name}_DEPS ${_opt_name} pos)
    if(pos EQUAL -1)
      list(APPEND ${_dep_name}_DEPS ${_opt_name})
      set(${_dep_name}_DEPS ${${_dep_name}_DEPS} CACHE INTERNAL "Dependencies ON with package ${_dep_name}" FORCE)
    endif()
  else()
    list(LENGTH ${_dep_name}_DEPS len)
    list(FIND ${_dep_name}_DEPS ${_opt_name} pos)
    if((len EQUAL 1) AND (NOT pos EQUAL -1))
      akantu_message("add_package_dependecies: Restore state ${_dep_name}:${${_dep_name}} (${pos})")
      set(${_dep_name} ${${_dep_name}_OLD} CACHE BOOL "${_dep_desc}" FORCE)
      unset(${_dep_name}_OLD CACHE)
    endif()

    if(NOT pos EQUAL -1)
      list(REMOVE_AT ${_dep_name}_DEPS ${pos})
      set(${_dep_name}_DEPS ${${_dep_name}_DEPS} CACHE INTERNAL "Nb dependencies with package ${_dep_name}" FORCE)
    endif()
  endif()
endmacro()

#===============================================================================
# macro to add meta packages
macro(add_meta_package PKG DESC DEFAULT)
  akantu_message("add_meta_package: register meta option ${PKG} ${DESC} ${DEFAULT}")
  package_pkg_name (${PKG} _pkg_name)
  package_desc_name(${PKG} _desc_name)

  set(${_desc_name} ${DESC})
  option(${_pkg_name} ${DESC} ${DEFAULT})

  foreach(_dep ${ARGN})
    package_opt_name (${_dep} _dep_name)
    mark_as_advanced(${_dep_name})
    add_package_dependecies(${PKG} ${_dep})
  endforeach()
endmacro()