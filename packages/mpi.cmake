#===============================================================================
# @file   CMakeLists.txt
# @author Richart Nicolas <nicolas.richart@epfl.ch>
# @date   Fri Sep 29 16:46:30 2010 
#
# @brief package description for mpi
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
add_optional_package(MPI "Add MPI support in akantu" OFF PREFIX MPI_C MPI DEPENDS SCOTCH)
set(AKANTU_MPI_FILES
  synchronizer/static_communicator_mpi.cc
  synchronizer/static_communicator_mpi_inline_impl.hh
  synchronizer/static_communicator_mpi.hh
  )
