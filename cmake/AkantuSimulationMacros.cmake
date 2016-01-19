#===============================================================================
# @file   AkantuSimulationMacros.cmake
#
# @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
# @author Nicolas Richart <nicolas.richart@epfl.ch>
#
# @date creation: Mon Jan 18 2016
#
# @brief  macros for examples
#
# @section LICENSE
#
# Copyright (©) 2015 EPFL (Ecole Polytechnique Fédérale de Lausanne) Laboratory
# (LSMS - Laboratoire de Simulation en Mécanique des Solides)
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

if(__AKANTU_SIMULATION_MACROS)
  return()
endif()
set(__AKANTU_SIMULATION_MACROS TRUE)



#===============================================================================
function(add_akantu_simulation simulation_name)
  _add_akantu_simulation(${simulation_name} ${ARGN})
endfunction()

#===============================================================================
function(_add_akantu_simulation simulation_name)
  set(multi_variables
    SOURCES
    FILES_TO_COPY
    DEPENDS
    DIRECTORIES_TO_CREATE
    COMPILE_OPTIONS
    USE_PACKAGES
    )

  cmake_parse_arguments(_simulation
    ""
    "LIST_FILES"
    "${multi_variables}"
    ${ARGN}
    )

  set(_deps_OK TRUE)
  if(_simulation_USE_PACKAGES)
    foreach(_pkg ${_simulation_USE_PACKAGES})
      package_is_activated(${_pkg} _activated)
      if(_activated)
        package_get_include_dir(${_pkg} _inc_dir)
        list(APPEND _simulation_INCLUDE_DIRS ${_inc_dir})

        package_get_libraries(${_pkg} _libraries)
        list(APPEND _simulation_LIBRARIES ${_libraries})

        package_get_compile_flags(${_pkg} _flags)
        list(APPEND _simulation_COMPILE_FLAGS "${_flags}")
      else()
        message("${simulation_name} use ${_pkg} but Akantu "
          " has been compiled without this package")
        set(_deps_OK FALSE)
      endif()
    endforeach()
  endif()

  if(_deps_OK)
    add_executable(${simulation_name}
      ${_simulation_UNPARSED_ARGUMENTS} ${_simulation_SOURCES})

    target_link_libraries(${simulation_name} akantu ${_simulation_LIBRARIES})
    target_include_directories(${simulation_name} PRIVATE ${AKANTU_INCLUDE_DIRS} ${_simulation_INCLUDE_DIRS})

    if(_simulation_DEPENDS)
      foreach(_deps ${_simulation_DEPENDS})
        get_target_property(_type ${_deps} TYPE)
        if(_type STREQUAL "SHARED_LIBRARY"
            OR _type STREQUAL "STATIC_LIBRARY")
          target_link_libraries(${simulation_name} ${_deps})
        else()
          add_dependencies(${simulation_name} ${_deps})
        endif()
      endforeach()
    endif()

    if(_simulation_COMPILE_OPTIONS OR _simulation_COMPILE_FLAGS)
      set_target_properties(${test_
        PROPERTIES COMPILE_FLAGS "${_simulation_COMPILE_FLAGS} ${_simulation_COMPILE_OPTIONS}")
    endif()

    # copy the needed files to the build folder
    if(_simulation_FILES_TO_COPY)
      file(COPY ${_simulation_FILES_TO_COPY} DESTINATION .)
    endif()

    # create the needed folders in the build folder
    if(_simulation_DIRECTORIES_TO_CREATE)
      foreach(_dir ${_simulation_DIRECTORIES_TO_CREATE})
        if(IS_ABSOLUTE ${dir})
          file(MAKE_DIRECTORY ${_dir})
        else()
          file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${_dir})
        endif()
      endforeach()
    endif()
  endif()


  if(_simulation_LIST_FILES)
    set(_simulation_files)

    foreach(_file ${_simulation_SOURCES} ${_simulation_UNPARSED_ARGUMENTS} ${_simulation_FILES_TO_COPY})
      list(APPEND _simulation_files ${CMAKE_CURRENT_SOURCE_DIR}/${_file})
    endforeach()

    foreach(_dep ${_simulation_DEPENDS})
      get_target_list_of_associated_files(${_dep} _dep_ressources)
      if(_dep_ressources)
        list(APPEND _simulation_files "${_dep_ressources}")
      endif()
    endforeach()

    set(${_simulation_LIST_FILES} ${_simulation_files} PARENT_SCOPE)
  endif()
endfunction()
