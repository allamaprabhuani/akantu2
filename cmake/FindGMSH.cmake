#===============================================================================
# @file   FindGMSH.cmake
# @author Nicolas Richart <nicolas.richart@epfl.ch>
# @date   Thu Oct 14 13:15:47 2010
#
# @brief  Find gmsh and delacre the add_mesh macro
#
# @section LICENSE
#
# Copyright (©) 2010-2011 EPFL (Ecole Polytechnique fédérale de Lausanne)
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

find_program(GMSH gmsh
  DOC "The mesh generetor gmsh")

find_package(PackageHandleStandardArgs)
find_package_handle_standard_args(GMSH DEFAULT_MSG GMSH)

#===============================================================================
macro(PARSE_ARGUMENTS prefix arg_names option_names)
  set(DEFAULT_ARGS)
  foreach(arg_name ${arg_names})
    set(${prefix}_${arg_name})
  endforeach(arg_name)
  foreach(option ${option_names})
    set(${prefix}_${option} FALSE)
  endforeach(option)

  set(current_arg_name DEFAULT_ARGS)
  set(current_arg_list)
  foreach(arg ${ARGN})
    set(larg_names ${arg_names})
    list(FIND larg_names "${arg}" is_arg_name)
    if (is_arg_name GREATER -1)
      set(${prefix}_${current_arg_name} ${current_arg_list})
      set(current_arg_name ${arg})
      set(current_arg_list)
    else (is_arg_name GREATER -1)
      set(loption_names ${option_names})
      list(FIND loption_names "${arg}" is_option)
      if (is_option GREATER -1)
	     set(${prefix}_${arg} TRUE)
      else (is_option GREATER -1)
	     set(current_arg_list ${current_arg_list} ${arg})
      endif (is_option GREATER -1)
    endif (is_arg_name GREATER -1)
  endforeach(arg)
  set(${prefix}_${current_arg_name} ${current_arg_list})
endmacro(PARSE_ARGUMENTS)

#===============================================================================
macro(ADD_MESH MESH_TARGET GEO_FILE DIM ORDER)
  if(GMSH_FOUND)
    set(arguments
      ${MESH_TARGET} ${GEO_FILE} ${DIM} ${ORDER}
      ${ARGN}
      )

    parse_arguments(ADD_MESH
      "OUTPUT"
      ${arguments}
      )

    set(_geo_file ${CMAKE_CURRENT_SOURCE_DIR}/${GEO_FILE})

    if(ADD_MESH_OUTPUT)
      set(_msh_file ${CMAKE_CURRENT_BINARY_DIR}/${ADD_MESH_OUTPUT})
    else(ADD_MESH_OUTPUT)
      get_filename_component(_msh_file "${GEO_FILE}" NAME_WE)
      set(_msh_file ${CMAKE_CURRENT_BINARY_DIR}/${_msh_file}.msh)
    endif(ADD_MESH_OUTPUT)

    if(EXISTS ${_geo_file})
      add_custom_command(
	OUTPUT ${_msh_file}
	DEPENDS ${_geo_file}
	COMMAND ${GMSH}
	ARGS -${DIM} -order ${ORDER} -optimize -o ${_msh_file} ${_geo_file} 2>&1 > /dev/null
	)
      add_custom_target(${MESH_TARGET}
	DEPENDS ${_msh_file})
    else(EXISTS ${_geo_file})
      message(FATAL_ERROR "File ${_geo_file} not found")
    endif(EXISTS ${_geo_file})
  endif(GMSH_FOUND)
endmacro(ADD_MESH)
