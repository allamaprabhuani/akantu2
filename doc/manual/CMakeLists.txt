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

macro(generate_doc_dependency PKG FILENAME LEVEL)

  if (${LEVEL} GREATER 0)
    string(REPLACE "_" "\\_" PKG_TEX ${PKG})
#    file(APPEND ${FILENAME} "\\item ${PKG_TEX}: \\\\")
    file(APPEND ${FILENAME} "\\noindent ${PKG_TEX}:")
  endif()

  if (${PKG}_DOCUMENTATION AND ${LEVEL} EQUAL 0)
    file(APPEND ${FILENAME} ${${PKG}_DOCUMENTATION})
  else()
    string(TOLOWER "${PKG}" _l_pkg)
    string(REPLACE "akantu_" "" _l_pkg_short ${_l_pkg})
    string(REPLACE "_" "\\_" _l_pkg_tex ${_l_pkg_short})

    file(GLOB _file_path "${PROJECT_SOURCE_DIR}/packages/*_${_l_pkg_short}.cmake")
    file(APPEND ${FILENAME} "{\\color{red} TODO: No Documentation in \\href{${_file_path}}{${_l_pkg_tex}.cmake}\\\\")

    file(APPEND ${FILENAME} "looking for the sequence: "
"set(${PKG_TEX}\\_DOCUMENTATION \" documentation text \" ) }" )

  endif()

  list(LENGTH "${${PKG}_dependencies}" _len)
  if (_len GREATER 0)
    if (${LEVEL} EQUAL 0)
#      file(APPEND ${FILENAME} "\\begin{itemize}\n")
      file(APPEND ${FILENAME} "This package depends on the options:\\\\ \\hrule \\vspace{.5cm}\n")
    endif()

    math(EXPR _sub_level "${LEVEL}+1")
    foreach(_pkg ${${PKG}_dependencies})
      generate_doc_dependency(${_pkg} ${FILENAME} ${_sub_level})
    endforeach()

    if (${LEVEL} EQUAL 0)
#      file(APPEND ${FILENAME} "\\end{itemize}\n")
    endif()
  endif()

endmacro()

find_program(RUBBER_EXECUTABLE "rubber")
if (NOT RUBBER_EXECUTABLE)
  message(ERROR "Manual cannot be built without rubber latex compiler")
endif()

set(_relative_manual_files)
foreach(_f  ${_all_manual_files})
  file(RELATIVE_PATH _rel_f ${CMAKE_CURRENT_SOURCE_DIR} "${PROJECT_SOURCE_DIR}/${_f}")
  list(APPEND _relative_manual_files ${_rel_f})
  get_filename_component(_file_path ${_rel_f} PATH)
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


set(DOC_DEPS_TEX_FILENAME "${CMAKE_CURRENT_BINARY_DIR}/manual-package-doc.tex")
file(WRITE ${DOC_DEPS_TEX_FILENAME} "")
file(APPEND ${DOC_DEPS_TEX_FILENAME} "\\section{Package dependencies}\n")
foreach (_pkg ${AKANTU_PACKAGE_LIST})
  string(REPLACE "_" "\\_" PKG_TEX ${_pkg})
  file(APPEND ${DOC_DEPS_TEX_FILENAME} "\\subsection{" ${PKG_TEX} "}\n" )
  generate_doc_dependency(AKANTU_${_pkg} ${DOC_DEPS_TEX_FILENAME} 0)
endforeach()


