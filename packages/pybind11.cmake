#===============================================================================
# @file   pybind11.cmake
#
# @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
# @author Nicolas Richart <nicolas.richart@epfl.ch>
#
# @date creation: Fri Dec 22 2017
# @date last modification: Wed Dec 04 2019
#
# @brief  package description for the pybind11 binding
#
#
# @section LICENSE
#
# Copyright (©) 2016-2021 EPFL (Ecole Polytechnique Fédérale de Lausanne)
# Laboratory (LSMS - Laboratoire de Simulation en Mécanique des Solides)
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


set(PYBIND11_PYTHON_VERSION ${AKANTU_PREFERRED_PYTHON_VERSION} CACHE INTERNAL "")

package_declare(pybind11 EXTERNAL
  EXTRA_PACKAGE_OPTIONS ARGS "2.4.2;CONFIG" LINK_LIBRARIES pybind11::embed PREFIX pybind11
  DESCRIPTION "Akantu's pybind11 interface"
  SYSTEM AUTO third-party/cmake/pybind11.cmake
  DEPENDS PythonLibsNew
  EXCLUDE_FROM_ALL
  )

package_add_third_party_script_variable(pybind11
  PYBIND11_VERSION "v2.4.2")
package_add_third_party_script_variable(pybind11
  PYBIND11_GIT "https://github.com/pybind/pybind11.git")
