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


find_program(SPHINX_BUILD_EXECUTABLE
  NAMES sphinx-build)

if(SPHINX_BUILD_EXECUTABLE)
  execute_process(
    COMMAND ${SPHINX_BUILD_EXECUTABLE} --version
    OUTPUT_VARIABLE SPHINX_VERSION
    )

  string(REPLACE "sphinx-build" "" SPHINX_VERSION "${SPHINX_VERSION}")
  string(STRIP "${SPHINX_VERSION}" SPHINX_VERSION)
endif()

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Sphinx
  REQUIRED_VARS SPHINX_BUILD_EXECUTABLE
  VERSION_VAR SPHINX_VERSION
)

mark_as_advanced(SPHINX_BUILD_EXECUTABLE)
