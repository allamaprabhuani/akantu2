#===============================================================================
# @file   scalapack.cmake
#
# @author Nicolas Richart <nicolas.richart@epfl.ch>
#
# @date creation: Fri Oct 19 2012
# @date last modification: Fri Jan 22 2016
#
# @brief  package description for mumps support
#
#
# @section LICENSE
#
# Copyright (©) 2010-2021 EPFL (Ecole Polytechnique Fédérale de Lausanne)
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


package_declare(ScaLAPACK EXTERNAL
  DESCRIPTION "Add ScaLAPACK support in akantu"
  SYSTEM OFF third-party/cmake/scalapack.cmake
  DEPENDS MPI
  )


package_add_third_party_script_variable(ScaLAPACK
  SCALAPACK_VERSION "2.0.2")
package_add_third_party_script_variable(ScaLAPACK
  SCALAPACK_ARCHIVE "http://www.netlib.org/scalapack/scalapack-${SCALAPACK_VERSION}.tgz")
package_add_third_party_script_variable(ScaLAPACK
  SCALAPACK_ARCHIVE_HASH_2.0.2 "MD5=2f75e600a2ba155ed9ce974a1c4b536f")

package_set_package_system_dependency(ScaLAPACK deb-src libscalapack-mpi-dev)

package_declare_extra_files_to_package(ScaLAPACK
  PROJECT
    third-party/scalapack_2.0.2.patch
    third-party/cmake/scalapack.cmake
    cmake/Modules/FindScaLAPACK.cmake
  )


package_declare_documentation(ScaLAPACK
  "ScaLAPACK is used to compile mumps as a third-party"
  ""
  "Under Ubuntu (14.04 LTS) the installation can be performed using the commands:"
  "\\begin{command}"
  "  > sudo apt-get install libscalapack-mpi-dev"
  "\\end{command}"
  ""
)
