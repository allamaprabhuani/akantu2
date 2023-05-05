#===============================================================================
# Copyright (©) 2015-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
# Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
#
# This file is part of Akantu
# 
# Akantu is free software: you can redistribute it and/or modify it under the
# terms of the GNU Lesser General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option) any
# later version.
# 
# Akantu is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
# A PARTICULAR PURPOSE. See the GNU Lesser General Public License for more
# details.
# 
# You should have received a copy of the GNU Lesser General Public License along
# with Akantu. If not, see <http://www.gnu.org/licenses/>.
#
#===============================================================================


#===============================================================================
# Material specific
#===============================================================================
#-------------------------------------------------------------------------------
function(package_declare_material_infos pkg)
  cmake_parse_arguments(_opt_pkg
    ""
    ""
    "LIST;INCLUDE"
    ${ARGN})

  package_set_variable(MATERIAL_LIST ${pkg} ${_opt_pkg_LIST})
  package_set_variable(MATERIAL_INCLUDE ${pkg} ${_opt_pkg_INCLUDE})
endfunction()

#-------------------------------------------------------------------------------
function(package_get_all_material_includes includes)
  _package_get_variable_for_activated(MATERIAL_INCLUDE _includes)

  foreach(_mat_inc ${_includes})
    if(DEFINED _mat_includes)
      set(_mat_includes "${_mat_includes}\n#include \"${_mat_inc}\"")
    else()
      set(_mat_includes "#include \"${_mat_inc}\"")
    endif()
  endforeach()

  set(${includes} ${_mat_includes} PARENT_SCOPE)
endfunction()

#-------------------------------------------------------------------------------
function(package_get_all_material_lists lists)
  _package_get_variable_for_activated(MATERIAL_LIST _lists)

  foreach(_mat_list ${_lists})
    if(DEFINED _mat_lists)
      set(_mat_lists "${_mat_lists}\n  ${_mat_list}\t\t\t\\")
    else()
      set(_mat_lists "  ${_mat_list}\t\t\t\\")
    endif()
  endforeach()

  set(${lists} ${_mat_lists} PARENT_SCOPE)
endfunction()

# ------------------------------------------------------------------------------
# Extra files to consider in source package generated by CPack
# ------------------------------------------------------------------------------
function(package_declare_extra_files_to_package pkg)
  set(_types SOURCES MANUAL TESTS PROJECT)
  cmake_parse_arguments(_extra_files
    ""
    ""
    "${_types}"
    ${ARGN})

  set(_files ${_extra_files_UNPARSED_ARGUMENTS})

  package_get_sources_folder(${pkg} _folder_SOURCES)
  package_get_manual_folder(${pkg} _folder_MANUAL)
  package_get_tests_folder(${pkg} _folder_TESTS)
  set(_folder_PROJECT ${PROJECT_SOURCE_DIR})

  foreach(_type ${_types})
    if(_extra_files_${_type})
      foreach(_file ${_extra_files_${_type}})
        list(APPEND _files ${_folder_${_type}}/${_file})
        if(NOT EXISTS ${_folder_${_type}}/${_file})
          message(SEND_ERROR "The package ${pkg} tries to register the file ${_file} (as a ${_type} file).
This file cannot be found.")
        endif()
      endforeach()
    endif()
  endforeach()

  package_set_variable(EXTRA_FILES ${pkg} ${_files})
endfunction()

# ------------------------------------------------------------------------------
function(package_add_files_to_package)
  set(_files)
  foreach(_file ${ARGN})
    list(APPEND _files ${PROJECT_SOURCE_DIR}/${_file})
  endforeach()
  package_add_to_project_variable(EXTRA_FILES ${_files})
endfunction()

function(package_get_files_for_package files)
  package_get_project_variable(EXTRA_FILES _tmp)
  set(${files} ${_tmp} PARENT_SCOPE)
endfunction()


package_add_files_to_package(
  .clang-format
  AUTHORS
  README
  VERSION
  COPYING
  COPYING.lesser
  CTestConfig.cmake
  cmake/akantu_environement.sh.in
  cmake/akantu_environement.csh.in
  cmake/akantu_install_environement.sh.in
  cmake/akantu_install_environement.csh.in
  cmake/Modules/CMakeFlagsHandling.cmake
  cmake/Modules/CMakePackagesSystem.cmake
  cmake/Modules/CMakePackagesSystemGlobalFunctions.cmake
  cmake/Modules/CMakePackagesSystemPrivateFunctions.cmake
  cmake/Modules/CMakeVersionGenerator.cmake
  cmake/Modules/PCHgcc.cmake
  cmake/AkantuBuildTreeSettings.cmake.in
  cmake/AkantuConfig.cmake.in
  cmake/AkantuCPack.cmake
  cmake/AkantuCPackMacros.cmake
  cmake/AkantuInstall.cmake
  cmake/AkantuMacros.cmake
  cmake/AkantuPackagesSystem.cmake
  cmake/AkantuUse.cmake
  cmake/AkantuSimulationMacros.cmake
  cmake/material_lister.cc
  cmake/Modules/FindGMSH.cmake
  )
