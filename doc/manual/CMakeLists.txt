#===============================================================================
# @file   CMakeLists.txt
#
# @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
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

find_program(RUBBER_EXECUTABLE "rubber")
if (NOT RUBBER_EXECUTABLE)
  message(FATAL_ERROR "Manual cannot be built without rubber latex compiler")
endif()

set(_relative_manual_files)
foreach(_f  ${_all_manual_files})
  file(RELATIVE_PATH _rel_f ${CMAKE_CURRENT_SOURCE_DIR} "${PROJECT_SOURCE_DIR}/${_f}")
  list(APPEND _relative_manual_files ${_rel_f})
  message(${_rel_f})
  GET_FILENAME_COMPONENT(_file_path ${_rel_f} DIRECTORY)
  file(COPY ${_rel_f} DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/${_file_path})
endforeach()



set(MANUAL_OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/manual.pdf)

add_custom_command(
  OUTPUT ${MANUAL_OUTPUT}  
  COMMAND ${RUBBER_EXECUTABLE} -f -d manual
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  )

add_custom_target(manual ALL DEPENDS ${MANUAL_OUTPUT})

add_custom_target(manual-forced
  COMMAND ${RUBBER_EXECUTABLE} -f -d manual  
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  )

install(FILES ${MANUAL_OUTPUT} DESTINATION share/akantu/doc)
