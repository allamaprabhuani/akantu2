#===============================================================================
# @file   documentation_developer_manual.cmake
#
# @author Guillaume Anciaux <guillaume.anciaux@epfl.ch>
# @author Nicolas Richart <nicolas.richart@epfl.ch>
#
# @date creation: Wed Jun 11 2014
# @date last modification: Fri Jan 29 2021
#
# @brief  Doxygen documentation of the code
#
#
# @section LICENSE
#
# Copyright (©) 2015-2021 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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


package_declare(documentation
  DESCRIPTION "Build source documentation using Sphinx/Doxygen.")

package_set_package_system_dependency(documentation deb-src
  python3-sphinx python3-breathe doxygen graphviz)

package_declare_extra_files_to_package(documentation
  PROJECT doc/dev-doc/akantu.dox.in
          doc/dev-doc/conf.py.in
          doc/dev-doc/index.rst
          doc/dev-doc/reference.rst
          cmake/Modules/FindSphinx.cmake
  )
