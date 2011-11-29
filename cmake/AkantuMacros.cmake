#===============================================================================
# @file   AkantuMacros.cmake
# @author Nicolas Richart <nicolas.richart@epfl.ch>
# @date   Wed Feb  9 10:59:42 2011
#
# @brief  Set of macros used by akantu cmake files
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
# macro to include optional packages
macro(add_optional_package PACKAGE DESC DEFAULT)
  string(TOUPPER ${PACKAGE} _u_package)
  set(_option_name AKANTU_USE_${_u_package})

  option(${_option_name} ${DESC} ${DEFAULT})

  if(${_option_name})
    if(${PACKAGE} MATCHES BLAS)
      enable_language(Fortran)
    endif()
    find_package(${PACKAGE} REQUIRED)
    if(${_u_package}_FOUND)
      add_definitions(-DAKANTU_USE_${_u_package})
      set(AKANTU_EXTERNAL_LIB_INCLUDE_PATH ${AKANTU_EXTERNAL_LIB_INCLUDE_PATH}
        ${${_u_package}_INCLUDE_PATH}
        )
      set(AKANTU_EXTERNAL_LIBRARIES ${AKANTU_EXTERNAL_LIBRARIES}
        ${${_u_package}_LIBRARIES}
        )
      set(AKANTU_${_u_package} ON)
#      MESSAGE("found ${_u_package} by setting AKANTU_${_u_package}" )
    else(${_u_package}_FOUND)
 #     MESSAGE("not found ${_u_package}" )
      set(AKANTU_${_u_package} OFF)
    endif(${_u_package}_FOUND)
  endif(${_option_name})
endmacro()

#===============================================================================
macro(check_for_isnan result)
  include(CheckFunctionExists)
  check_function_exists(std::isnan HAVE_STD_ISNAN)
  if(HAVE_STD_ISNAN)
    set(result "std::isnan(x)")
  else()
    check_function_exists(isnan HAVE_ISNAN)
    if(HAVE_ISNAN)
      set(result "(::isnan(x))")
    else()
      check_function_exists(_isnan HAVE_ISNAN_MATH_H)
      if(HAVE_ISNAN_MATH_H)
        set(result "(_isnan(x))")
      else()
        set(result (x == std::numeric_limits<Real>::quiet_NAN()))
      endif()
    endif()
  endif()
endmacro()

#===============================================================================
macro(copy_files target_depend)
  foreach(_file ${ARGN})
    set(_target ${CMAKE_CURRENT_BINARY_DIR}/${_file})
    set(_source ${CMAKE_CURRENT_SOURCE_DIR}/${_file})
    add_custom_command(
      OUTPUT ${_target}
      COMMAND ${CMAKE_COMMAND} -E copy_if_different ${_source} ${_target}
      DEPENDS ${_source}
      )
    set(_target_name ${target_depend}_${_file})
    add_custom_target(${_target_name} DEPENDS ${_target})
    add_dependencies(${target_depend} ${_target_name})
  endforeach()
endmacro()


#===============================================================================
# Package Management
#===============================================================================

#===============================================================================
macro(add_all_packages package_dir)
#  message("PKG DIR : ${package_dir}")
  file(GLOB _akantu_package_list "${package_dir}/*.cmake")
#  message("PKG LIST : ${_akantu_package_list}")

  foreach(_pkg ${_akantu_package_list})
    include(${_pkg})
    get_filename_component(_basename ${_pkg} NAME)
    string(REGEX REPLACE "\\.cmake" "" _option_name ${_basename})
    string(TOUPPER "${_option_name}" _option_name)
    list(APPEND AKANTU_PACKAGE_NAMES_LIST_ALL ${_option_name})
    if (AKANTU_${_option_name})
      list(APPEND AKANTU_PACKAGE_NAMES_LIST_ON ${_option_name})
    endif()

    foreach(_file ${${_option_name}_FILES})
      list(APPEND _release_all_file ${_file})
    endforeach()
  endforeach()
#  message("ALL PKG : ${AKANTU_PACKAGE_NAMES_LIST_ALL}")
#  message("ON  PKG : ${AKANTU_PACKAGE_NAMES_LIST_ON}")
#  message("ALL FILE LIST : ${_release_all_file}")

  #check if there are some file in the package list that are not on the current directory
  file(GLOB_RECURSE _all_files RELATIVE ${CMAKE_SOURCE_DIR}/src "*.cc" "*.hh")
  foreach(_file ${_release_all_file})
    list(FIND _all_files ${_file} _index)
    if (_index EQUAL -1)
      message("The file ${_file} is registered in packages a package but is not present in the source directory.")
    endif()
  endforeach()

  #check dependencies
  foreach(_pkg ${AKANTU_PACKAGE_NAMES_LIST_ON})
    # differentiate the file types
    foreach(_dep ${${_option_name}_DEPENDS})
      if (NOT AKANTU_${dep})
        message(FATAL_ERROR "Package ${_pkg} depends on package ${_dep}. You need to activate it to make it work")
      endif()
    endforeach()
  endforeach()
endmacro()

#===============================================================================
macro(generate_source_list_from_packages source_dir source_files headers_files include_dirs)
  set(deb_depend "libc6")
  
#  message("SRC DIR : ${source_dir}")
  foreach(_option_name ${AKANTU_PACKAGE_NAMES_LIST_ON})
    # differentiate the file types
    foreach(_file ${${_option_name}_FILES})
      if(${_file} MATCHES ".*inline.*\\.cc")
        list(APPEND ${_option_name}_inlines ${_file})
      elseif(${_file} MATCHES ".*\\.hh")
        list(APPEND ${_option_name}_headers ${_file})
      else()
        list(APPEND ${_option_name}_srcs ${_file})
      endif()
    endforeach()

    # generates the include directory variable
    foreach(_file ${${_option_name}_headers})
      get_filename_component(_absolute_name ${_file} ABSOLUTE)
      get_filename_component(_include_dir ${_absolute_name} PATH)
      list(APPEND ${_option_name}_include_dirs ${_include_dir})
      list(REMOVE_DUPLICATES ${_option_name}_include_dirs)
    endforeach()

    # generate global lists for akantu to know what to build
    list(APPEND ${source_files} ${${_option_name}_srcs})
    list(APPEND ${headers_files} ${${_option_name}_headers})
    list(APPEND ${include_dirs}  ${${_option_name}_include_dirs})
    if (NOT "${${_option_name}_DEB_DEPEND}" STREQUAL "")
      set(deb_depend "${deb_depend}, ${${_option_name}_DEB_DEPEND}")
    endif()

#    message("PKG ${_option_name} SRCS : ${${_option_name}_srcs}")
#    message("PKG ${_option_name} HRDS : ${${_option_name}_headers}")
#    message("PKG ${_option_name} INCS : ${${_option_name}_include_dirs}")
  endforeach()

#  message("SRCS : ${${source_files}}")
#  message("HRDS : ${${headers_files}}")
#  message("INCS : ${${include_dirs}}")

SET(CPACK_DEBIAN_PACKAGE_DEPENDS "${deb_depend}" PARENT_SCOPE)
endmacro()