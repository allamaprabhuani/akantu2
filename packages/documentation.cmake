#===============================================================================
# Copyright (©) 2010-2023 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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


package_declare(documentation
  DESCRIPTION "Build source documentation using Sphinx/Doxygen.")

package_set_package_system_dependency(documentation deb-src
  python3-sphinx python3-breathe doxygen graphviz)

package_declare_extra_files_to_package(documentation
  PROJECT doc/dev-doc/akantu.dox.j2
          doc/dev-doc/conf.py
          doc/dev-doc/index.rst
          doc/dev-doc/reference.rst
          cmake/Modules/FindSphinx.cmake
  )
