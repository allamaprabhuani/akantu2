#===============================================================================
# Copyright (©) 2020-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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


# configured documentation tools and intermediate build results
set(BINARY_BUILD_DIR "${CMAKE_CURRENT_BINARY_DIR}/_build")
# Sphinx cache with pickled ReST documents
set(SPHINX_CACHE_DIR "${CMAKE_CURRENT_BINARY_DIR}/_doctrees")
# HTML output directory
set(SPHINX_HTML_DIR "${CMAKE_CURRENT_BINARY_DIR}/html")
set(SPHINX_OUTPUT "${SPHINX_HTML_DIR}/index.html")
set(SPHINX_INPUT "${CMAKE_CURRENT_BINARY_DIR}/conf.py")

# ---------------------------------------------------------------------------- #
# Sphinx                                                                       #
# ---------------------------------------------------------------------------- #
find_package(Sphinx REQUIRED)

set(SPHINX_VERBOSE_FLAG "-v")
if(CMAKE_VERBOSE_MAKEFILE)
  set(SPHINX_VERBOSE_FLAG -vv)
endif(CMAKE_VERBOSE_MAKEFILE)

configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/conf.py"
  "${SPHINX_INPUT}"
  @ONLY)

configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/manual/manual-bibliography.bib"
  "${CMAKE_CURRENT_BINARY_DIR}/manual/manual-bibliography.bib"
  COPYONLY)

set(SPHINX_PARALLEL_FLAG)
if (SPHINX_VERSION VERSION_GREATER_EQUAL 1.7.0)
  set(SPHINX_PARALLEL_FLAG -j 1)
endif()

set(_sphinx_command ${SPHINX_BUILD_EXECUTABLE}
  ${SPHINX_PARALLEL_FLAG}
  ${SPHINX_VERBOSE_FLAG} -b html
  -c "${CMAKE_CURRENT_BINARY_DIR}"
  -d "${SPHINX_CACHE_DIR}"
  "${CMAKE_CURRENT_SOURCE_DIR}"
  "${SPHINX_HTML_DIR}"
  )

file(GLOB_RECURSE _SPHINX_RST_SRCS
  "*.rst")
file(GLOB_RECURSE _SPHINX_MD_SRCS
  "*.md")

set(ENV(RUNNING_IN_CMAKE) "True")
add_custom_command(
  OUTPUT ${SPHINX_OUTPUT}
  COMMAND ${_sphinx_command}
  DEPENDS ${SPHINX_INPUT} ${_SPHINX_RST_SRCS} ${_SPHINX_MD_SRCS} "${CMAKE_CURRENT_SOURCE_DIR}/akantu.dox.j2"
  COMMENT "Building HTML documentation with Sphinx in ${SPHINX_HTML_DIR}"
  )

add_custom_target(sphinx-doc ALL
  DEPENDS ${SPHINX_OUTPUT})
