#===============================================================================
# Copyright (©) 2021-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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


if(DEFINED AKANTU_NEED_PYBIND11_LOADED)
  return()
endif()
set(AKANTU_NEED_PYBIND11_LOADED TRUE)


set(PYBIND11_PYTHON_VERSION ${AKANTU_PREFERRED_PYTHON_VERSION} CACHE INTERNAL "")

find_package(pybind11 QUIET)

if (NOT pybind11_FOUND)
  set(PYBIND11_VERSION "v2.10.3")
  set(PYBIND11_GIT "https://github.com/pybind/pybind11.git")

  include(${PROJECT_SOURCE_DIR}/third-party/cmake/pybind11.cmake)
else()
  message(STATUS "Found pybind11: ${pybind11_INCLUDE_DIRS} (${pybind11_VERSION})")
endif()
